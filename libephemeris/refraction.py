# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Atmospheric refraction via Gauss-Legendre quadrature of the refraction
integral through the ICAO Standard Atmosphere (ISO 2533:1975).

Instead of relying on empirical curve-fits (Bennett 1982, Saemundsson 1986,
etc.), this module computes refraction from first principles by evaluating
the exact refraction integral through a continuously varying model
atmosphere.

Physical model
--------------
The atmosphere follows the ICAO Standard Atmosphere temperature profile:

    Layer             Altitude (km)    Lapse rate (K/km)    T_base (K)
    ─────────────     ──────────────   ─────────────────    ──────────
    Troposphere       0 – 11           -6.5                 288.15
    Tropopause        11 – 20           0.0                 216.65
    Stratosphere 1    20 – 32          +1.0                 216.65
    Stratosphere 2    32 – 47          +2.8                 228.65
    Upper             47 – 51           0.0                 270.65
    Mesosphere 1      51 – 71          -2.8                 270.65
    Mesosphere 2      71 – 84.852      -2.0                 214.65

Pressure at each altitude follows from the hydrostatic equation
integrated exactly for each layer type (isothermal or linear lapse).

The refractive index of air is computed from:

    n(P, T) = 1 + 8.060e-5 * P / T

where P is in mbar and T in Kelvin.  The coefficient 8.060e-5 is derived
from the Barrell & Sears (1939) group-refractive-index for the visible
band, adjusted to include the average effect of atmospheric humidity
on astronomical refraction (~2% above the dry-air value).

The total refraction is computed by evaluating the integral:

    R = - integral from r_obs to r_top of
            (dn/dr) / (n * sqrt(n^2 * r^2 / C^2 - 1))  dr

where C = n_obs * r_obs * sin(z_obs) is the Bouguer (Snell) invariant.
This integral is evaluated using Gauss-Legendre quadrature with 120
points, giving sub-arcsecond convergence; the near-singular integrand
as z -> 90 deg limits accuracy at the horizon to ~0.2 arcsec (far below
this module's documented deviation envelope from the reference fits).

For the APP_TO_TRUE direction the function inverts the
TRUE_TO_APP computation numerically via Newton-Raphson iteration
(typically 4-5 steps to full float64 convergence).

References
----------
- ISO 2533:1975 "Standard Atmosphere"
- ICAO Doc 7488/3 "Manual of the ICAO Standard Atmosphere"
- Barrell, H. & Sears, J.E. (1939), Phil. Trans. Roy. Soc. A, 238, 1
- Green, R.M. (1985) "Spherical Astronomy", Cambridge Univ. Press, Ch. 4
- Auer, L.H. & Standish, E.M. (2000), AJ, 119, 2472
- Smart, W.M. (1977) "Textbook on Spherical Astronomy", Ch. VI
- Bomford, G. (1980) "Geodesy", 4th ed., Clarendon Press, §2.17-2.20

Deviation from the reference implementation
-------------------------------------------
The reference ephemeris computes refraction from empirical curve fits:
Sinclair's formula (quoted in Bennett 1982) inside its azimuth/altitude
and extended-refraction paths and the Saemundsson/Bennett pair in its
plain refraction path.

Plain ``refrac()`` deliberately keeps the ray-traced model above
(owner decision, 2026-06): it is physically grounded and at least as
accurate against rigorous benchmarks. Measured envelope of the plain
``refrac()`` difference vs. the reference's fits at p=1013.25 hPa, T=10 C:

- above ~15 deg altitude:           < 0.1 arcsec
- 5 .. 15 deg:                      < 1 arcsec
- 0 .. 5 deg (incl. horizon):       up to ~15 arcsec
- below the horizon dip:            up to ~0.18 deg (the fits are
                                    extrapolations there; neither model
                                    is observationally constrained)

``refrac_extended()`` and ``azalt()`` must instead be 1:1 with the
reference, so they use the reference-matched analytic model further down
(see calc_refraction_ref_app_to_true / calc_refraction_ref_true_to_app),
which reproduces the reference to < 30" across observer elevation,
pressure and temperature -- both above and below the dip.

Consequences: apparent altitudes from plain ``refrac()`` can differ from
the reference by the envelope above, which can flip near-horizon
visibility bits (eclipse/occultation *_VISIBLE flags) within a few
seconds of an event's horizon crossing. Rise/set times are NOT affected:
eclipse.rise_trans and everything built on it use the reference's own
Sinclair model (see eclipse._rise_true_to_apparent) so that rise/set,
twilight and eclipse timing agree with the reference to fractions of a
second.
"""

from __future__ import annotations

import math
from typing import Tuple

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
_R_EARTH: float = 6_371_000.0  # Mean Earth radius [m]
_G0: float = 9.80665  # Standard gravity [m/s^2]
_M_AIR: float = 0.028_964_4  # Molar mass of dry air [kg/mol]
_R_GAS: float = 8.314_462_618  # Universal gas constant [J/(mol*K)]
_EXPONENT: float = _G0 * _M_AIR / _R_GAS  # g*M/R ≈ 0.03416 K/m

# Refractive index coefficient.  Derived from Barrell & Sears (1939)
# group-refractive-index for the visible band (n_s - 1 = 2.93e-4 at
# 0 degC, 1013.25 mbar) scaled via the ideal gas law, with a ~2%
# adjustment for the mean effect of atmospheric humidity on optical
# refractivity (Owens 1967).  The value 8.060e-5 is chosen to
# minimise the maximum percentage error across all altitudes when
# compared with standard astronomical refraction tables.
_N_COEFF: float = 8.060e-5

# ICAO atmosphere layers: (base_altitude_m, lapse_rate_K_per_m, T_base_K)
_ICAO_LAYERS: list = [
    (0.0, -0.0065, 288.15),  # Troposphere
    (11_000.0, 0.0, 216.65),  # Tropopause (isothermal)
    (20_000.0, 0.001, 216.65),  # Stratosphere 1
    (32_000.0, 0.0028, 228.65),  # Stratosphere 2
    (47_000.0, 0.0, 270.65),  # Stratopause (isothermal)
    (51_000.0, -0.0028, 270.65),  # Mesosphere 1
    (71_000.0, -0.002, 214.65),  # Mesosphere 2
]
_ATMO_TOP: float = 84_852.0  # Top of modeled atmosphere [m]

# Standard sea-level conditions
_P0: float = 1013.25  # Standard pressure [mbar]
_T0: float = 288.15  # Standard temperature [K]


# ---------------------------------------------------------------------------
# Gauss-Legendre quadrature nodes and weights (120-point)
# ---------------------------------------------------------------------------
# Precomputed once at import time.  The nodes/weights are computed in pure
# Python (Legendre recurrence + Newton-Raphson, no numpy); the runtime
# integration loop is likewise pure Python.


def _gauss_legendre_nodes(n: int) -> Tuple[list, list]:
    """Compute Gauss-Legendre nodes and weights on [-1, 1].

    Uses the Legendre polynomial recurrence and Newton-Raphson
    root finding.  Pure Python, no external dependencies.
    """
    nodes = [0.0] * n
    weights = [0.0] * n
    m = (n + 1) // 2
    for i in range(m):
        # Initial guess via Tricomi approximation
        theta = math.pi * (i + 0.75) / (n + 0.5)
        x = math.cos(theta)

        dp = 1.0  # will be overwritten in the loop
        for _ in range(
            50
        ):  # pragma: no branch - Newton iteration always converges in <50 steps
            p0 = 1.0
            p1 = x
            for j in range(2, n + 1):
                p2 = ((2 * j - 1) * x * p1 - (j - 1) * p0) / j
                p0 = p1
                p1 = p2
            dp = n * (x * p1 - p0) / (x * x - 1.0)
            dx = -p1 / dp
            x += dx
            if abs(dx) < 1e-15:
                break

        w = 2.0 / ((1.0 - x * x) * dp * dp)
        # Place symmetric pair
        nodes[i] = -x
        nodes[n - 1 - i] = x
        weights[i] = w
        weights[n - 1 - i] = w

    return nodes, weights


_N_QUAD: int = 120
_QUAD_NODES, _QUAD_WEIGHTS = _gauss_legendre_nodes(_N_QUAD)


# ---------------------------------------------------------------------------
# ICAO Standard Atmosphere: temperature & pressure profiles
# ---------------------------------------------------------------------------


def _icao_temperature(h: float) -> float:
    """Temperature at altitude *h* [m] per ICAO.  Returns Kelvin."""
    for i in range(len(_ICAO_LAYERS) - 1, -1, -1):
        h_base, lapse, T_base = _ICAO_LAYERS[i]
        if h >= h_base:
            return T_base + lapse * (h - h_base)
    return _ICAO_LAYERS[0][2]


def _icao_pressure(h: float) -> float:
    """Pressure at altitude *h* [m] per ICAO.  Returns mbar."""
    p = _P0
    for i, (h_base, lapse, T_base) in enumerate(_ICAO_LAYERS):
        h_top = _ICAO_LAYERS[i + 1][0] if i + 1 < len(_ICAO_LAYERS) else _ATMO_TOP
        if h <= h_top:
            dh = h - h_base
            if abs(lapse) < 1e-10:
                p *= math.exp(-_EXPONENT * dh / T_base)
            else:
                p *= ((T_base + lapse * dh) / T_base) ** (-_EXPONENT / lapse)
            return p
        dh = h_top - h_base
        if abs(lapse) < 1e-10:
            p *= math.exp(-_EXPONENT * dh / T_base)
        else:
            p *= ((T_base + lapse * dh) / T_base) ** (-_EXPONENT / lapse)
    return p


def _refractive_index(p_mbar: float, T_kelvin: float) -> float:
    """Refractive index of air at pressure *p_mbar* and temperature
    *T_kelvin*.  See module docstring for the coefficient derivation."""
    if T_kelvin <= 0 or p_mbar <= 0:
        return 1.0
    return 1.0 + _N_COEFF * p_mbar / T_kelvin


# ---------------------------------------------------------------------------
# Modified atmosphere: observer-supplied surface conditions
# ---------------------------------------------------------------------------


def _modified_temperature(
    h: float, obs_alt: float, surface_T: float, lapse_rate: float
) -> float:
    """Temperature at altitude *h* using the observer's surface conditions
    in the troposphere and the ICAO profile above 11 km."""
    tropo_top = 11_000.0
    if h < tropo_top:
        return surface_T - lapse_rate * (h - obs_alt)
    T_at_tropo = surface_T - lapse_rate * (tropo_top - obs_alt)
    for i in range(1, len(_ICAO_LAYERS)):
        h_base = _ICAO_LAYERS[i][0]
        lapse_i = _ICAO_LAYERS[i][1]
        h_top = _ICAO_LAYERS[i + 1][0] if i + 1 < len(_ICAO_LAYERS) else _ATMO_TOP
        if h <= h_top:
            if i == 1:
                return T_at_tropo + lapse_i * (h - h_base)
            offset = T_at_tropo - _ICAO_LAYERS[1][2]
            return _ICAO_LAYERS[i][2] + offset + lapse_i * (h - h_base)
    return 200.0


def _modified_pressure(
    h: float, obs_alt: float, obs_P: float, surface_T: float, lapse_rate: float
) -> float:
    """Pressure at altitude *h* using the modified atmosphere.
    Uses exact hydrostatic integrals for each layer type."""
    tropo_top = 11_000.0
    if h <= tropo_top:
        T_obs = surface_T
        T_h = T_obs - lapse_rate * (h - obs_alt)
        if abs(lapse_rate) < 1e-10:
            return obs_P * math.exp(-_EXPONENT * (h - obs_alt) / T_obs)
        return obs_P * (T_h / T_obs) ** (_EXPONENT / lapse_rate)

    T_obs = surface_T
    T_tropo = T_obs - lapse_rate * (tropo_top - obs_alt)
    if abs(lapse_rate) < 1e-10:
        P_tropo = obs_P * math.exp(-_EXPONENT * (tropo_top - obs_alt) / T_obs)
    else:
        P_tropo = obs_P * (T_tropo / T_obs) ** (_EXPONENT / lapse_rate)

    p = P_tropo
    T_at_tropo = T_tropo
    for i in range(1, len(_ICAO_LAYERS)):
        h_base = _ICAO_LAYERS[i][0]
        lapse_i = _ICAO_LAYERS[i][1]
        h_top_i = _ICAO_LAYERS[i + 1][0] if i + 1 < len(_ICAO_LAYERS) else _ATMO_TOP
        if i == 1:
            T_base_i = T_at_tropo
        else:
            T_base_i = _ICAO_LAYERS[i][2] + (T_at_tropo - _ICAO_LAYERS[1][2])
        if h <= h_top_i:
            dh = h - h_base
            if abs(lapse_i) < 1e-10:
                p *= math.exp(-_EXPONENT * dh / T_base_i)
            else:
                p *= ((T_base_i + lapse_i * dh) / T_base_i) ** (-_EXPONENT / lapse_i)
            return p
        dh = h_top_i - h_base
        if abs(lapse_i) < 1e-10:
            p *= math.exp(-_EXPONENT * dh / T_base_i)
        else:
            p *= ((T_base_i + lapse_i * dh) / T_base_i) ** (-_EXPONENT / lapse_i)
    return p


# ---------------------------------------------------------------------------
# Atmospheric profile functions for the integrand
# ---------------------------------------------------------------------------


def _n_at_height(
    h: float, obs_alt: float, obs_P: float, obs_T_K: float, lapse_rate: float
) -> float:
    """Refractive index at altitude *h* in the modified atmosphere."""
    T = _modified_temperature(h, obs_alt, obs_T_K, lapse_rate)
    P = _modified_pressure(h, obs_alt, obs_P, obs_T_K, lapse_rate)
    return _refractive_index(P, T)


def _dn_dr_at_height(
    h: float, obs_alt: float, obs_P: float, obs_T_K: float, lapse_rate: float
) -> float:
    """Derivative dn/dr at altitude *h* via central difference.

    Uses a ±1 m central-difference step (2 m total span) for numerical
    differentiation, which is small enough for accuracy but large enough
    to avoid cancellation errors.
    """
    delta = 1.0  # metres
    n_plus = _n_at_height(h + delta, obs_alt, obs_P, obs_T_K, lapse_rate)
    n_minus = _n_at_height(h - delta, obs_alt, obs_P, obs_T_K, lapse_rate)
    return (n_plus - n_minus) / (2.0 * delta)


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

# Single-slot profile cache: (key, profile) stored atomically as one tuple
_refr_cache_entry: (
    tuple[tuple[float, float, float, float], _AtmosphereProfile] | None
) = None


class _AtmosphereProfile:
    """Precomputed atmosphere profile at Gauss-Legendre quadrature points.

    The change of variables u = sqrt((r - r_obs) / (r_top - r_obs))
    removes the integrable 1/sqrt(r - r_obs) singularity at the
    horizon (z = 90°).  In terms of u on [0, 1]:

        r = r_obs + L * u^2          (L = r_top - r_obs)
        dr = 2 * L * u * du
        sqrt(r - r_obs) = sqrt(L) * u     (cancels with integrand)

    The atmospheric properties (n, dn/dr) at each quadrature u-point
    are precomputed once and cached.
    """

    __slots__ = ("n_obs", "r_obs", "r_top", "L", "_q_r", "_q_n", "_q_dn_dr")

    def __init__(self, obs_alt: float, obs_P: float, obs_T_K: float, lapse_rate: float):
        self.r_obs = _R_EARTH + obs_alt
        self.r_top = _R_EARTH + _ATMO_TOP
        self.L = self.r_top - self.r_obs

        # Map Gauss-Legendre nodes from [-1,1] to u in [0,1]
        # u = (t + 1) / 2,  du = dt / 2
        q_r = []
        q_n = []
        q_dn = []
        for k in range(_N_QUAD):
            u = 0.5 * (_QUAD_NODES[k] + 1.0)
            r = self.r_obs + self.L * u * u
            h = r - _R_EARTH
            n = _n_at_height(h, obs_alt, obs_P, obs_T_K, lapse_rate)
            dn = _dn_dr_at_height(h, obs_alt, obs_P, obs_T_K, lapse_rate)
            q_r.append(r)
            q_n.append(n)
            q_dn.append(dn)

        self._q_r = q_r
        self._q_n = q_n
        self._q_dn_dr = q_dn

        # n at observer (u=0 is the first Gauss point, but compute exactly)
        self.n_obs = _n_at_height(obs_alt, obs_alt, obs_P, obs_T_K, lapse_rate)


def _get_profile(
    obs_alt: float, obs_P: float, obs_T_K: float, lapse_rate: float
) -> _AtmosphereProfile:
    """Return a cached atmosphere profile.

    The single-slot cache stores (key, profile) as one tuple so a
    concurrent reader can never pair the new key with the previous
    profile (the old code assigned the key before the profile).
    """
    global _refr_cache_entry
    key = (obs_alt, obs_P, obs_T_K, lapse_rate)
    entry = _refr_cache_entry
    if entry is not None and entry[0] == key:
        return entry[1]
    prof = _AtmosphereProfile(obs_alt, obs_P, obs_T_K, lapse_rate)
    _refr_cache_entry = (key, prof)
    return prof


# ---------------------------------------------------------------------------
# Core: Gauss-Legendre quadrature of the refraction integral
# ---------------------------------------------------------------------------


def _refraction_integral(
    observed_zenith_deg: float,
    profile: _AtmosphereProfile,
) -> float:
    """Evaluate the refraction integral for a ray arriving at the observer
    at apparent zenith angle *observed_zenith_deg*.

    The refraction integral is (Green 1985, eq. 4.25):

        R = - integral from r_obs to r_top of
                (dn/dr) / (n * sqrt(n^2 * r^2 / C^2 - 1))  dr

    where C = n_obs * r_obs * sin(z_obs) is the Bouguer invariant.

    To handle the integrable singularity at r = r_obs when z -> 90
    (the integrand has a 1/sqrt(r - r_obs) singularity), we use the
    substitution u = sqrt((r - r_obs) / L):

        r = r_obs + L * u^2
        dr = 2 * L * u * du
        sqrt(r - r_obs) = sqrt(L) * u

    In the integrand, sqrt(n^2 * r^2 / C^2 - 1) ~ K * u near the
    observer, so the u from dr cancels the 1/u from the square root,
    leaving a smooth integrand in u.

    Returns refraction in degrees (>= 0).
    """
    z_rad = math.radians(observed_zenith_deg)
    sin_z = math.sin(z_rad)
    C = profile.n_obs * profile.r_obs * sin_z
    C2 = C * C
    L = profile.L

    q_r = profile._q_r
    q_n = profile._q_n
    q_dn = profile._q_dn_dr

    total = 0.0
    for k in range(_N_QUAD):
        u = 0.5 * (_QUAD_NODES[k] + 1.0)
        r = q_r[k]
        n = q_n[k]
        dn_dr = q_dn[k]

        nr = n * r
        arg = nr * nr / C2 - 1.0
        if arg <= 0.0:
            continue  # ray doesn't reach this altitude

        # Integrand in u-coordinates:
        # f(u) = -(dn/dr) / (n * sqrt(arg)) * 2 * L * u
        # The factor u from the substitution dr = 2*L*u*du cancels
        # the 1/u behaviour of 1/sqrt(arg) near u=0.
        f = -dn_dr / (n * math.sqrt(arg)) * 2.0 * L * u

        # Weight includes the [-1,1] -> [0,1] Jacobian (factor 0.5)
        total += 0.5 * _QUAD_WEIGHTS[k] * f

    return max(0.0, math.degrees(total))


def _trace_ray(
    true_zenith_deg: float,
    obs_alt: float = 0.0,
    obs_pressure: float = _P0,
    obs_temperature_C: float = 15.0,
    lapse_rate: float = 0.0065,
) -> float:
    """Compute refraction for a source at true zenith angle
    *true_zenith_deg*.

    Since the Bouguer invariant C depends on the *observed* zenith
    angle z_obs, and we know the *true* zenith angle z_true, we
    iterate:

        z_obs = z_true - R(z_obs)

    Convergence is fast (3-4 iterations) because R is small.

    Returns refraction in degrees (>= 0).
    """
    if true_zenith_deg <= 0.0:
        return 0.0
    if obs_pressure <= 0.0:
        return 0.0

    obs_T_K = obs_temperature_C + 273.15

    # Non-physical temperature: at or below absolute zero (attemp <= -273.15 C)
    # there is no meaningful atmosphere. Report zero refraction, consistent with
    # the pressure<=0 guard above, rather than dividing by a zero absolute
    # temperature (ZeroDivisionError) or silently propagating a negative one.
    if obs_T_K <= 0.0:
        return 0.0

    # Linear extrapolation for objects below the geometric horizon
    if true_zenith_deg > 91.0:
        r_91 = _trace_ray(91.0, obs_alt, obs_pressure, obs_temperature_C, lapse_rate)
        r_90 = _trace_ray(90.0, obs_alt, obs_pressure, obs_temperature_C, lapse_rate)
        slope = r_91 - r_90
        return max(0.0, r_91 + slope * (true_zenith_deg - 91.0))

    profile = _get_profile(obs_alt, obs_pressure, obs_T_K, lapse_rate)

    # Iterate: z_obs = z_true - R(z_obs)
    z_obs = true_zenith_deg
    R = 0.0
    for _ in range(8):
        R = _refraction_integral(z_obs, profile)
        z_obs_new = true_zenith_deg - R
        if abs(z_obs_new - z_obs) < 1e-12:
            break
        z_obs = z_obs_new

    return max(0.0, R)


# ---------------------------------------------------------------------------
# Public functions (used internally by utils.py)
# ---------------------------------------------------------------------------


def calc_refraction_true_to_app(
    true_alt: float,
    pressure: float = _P0,
    temperature_C: float = 15.0,
    obs_alt: float = 0.0,
    lapse_rate: float = 0.0065,
) -> float:
    """Compute the refraction for a given true (geometric) altitude.

    Parameters
    ----------
    true_alt : float
        True altitude above the geometric horizon, in degrees.
    pressure : float
        Atmospheric pressure at the observer, in mbar.
    temperature_C : float
        Temperature at the observer, in degrees Celsius.
    obs_alt : float
        Observer altitude above sea level, in metres.
    lapse_rate : float
        Tropospheric lapse rate in K/m.

    Returns
    -------
    float
        Refraction in degrees (>= 0).
    """
    if pressure <= 0:
        return 0.0
    zenith = 90.0 - true_alt
    return _trace_ray(zenith, obs_alt, pressure, temperature_C, lapse_rate)


def calc_refraction_app_to_true(
    apparent_alt: float,
    pressure: float = _P0,
    temperature_C: float = 15.0,
    obs_alt: float = 0.0,
    lapse_rate: float = 0.0065,
) -> float:
    """Compute the refraction for a given apparent (observed) altitude.

    Uses Newton-Raphson iteration to invert the true-to-apparent
    mapping.

    Parameters
    ----------
    apparent_alt : float
        Apparent (observed) altitude, in degrees.
    pressure : float
        Atmospheric pressure at the observer, in mbar.
    temperature_C : float
        Temperature at the observer, in degrees Celsius.
    obs_alt : float
        Observer altitude above sea level, in metres.
    lapse_rate : float
        Tropospheric lapse rate in K/m.

    Returns
    -------
    float
        Refraction in degrees (>= 0).
    """
    if pressure <= 0:
        return 0.0

    # Initial estimate: true_alt ≈ apparent_alt
    true_est = apparent_alt

    for _ in range(8):
        r = calc_refraction_true_to_app(
            true_est, pressure, temperature_C, obs_alt, lapse_rate
        )
        residual = true_est + r - apparent_alt
        if abs(residual) < 1e-12:
            break

        # Numerical derivative via forward difference
        dt = 0.001
        r_plus = calc_refraction_true_to_app(
            true_est + dt, pressure, temperature_C, obs_alt, lapse_rate
        )
        dapp_dtrue = 1.0 + (r_plus - r) / dt
        if abs(dapp_dtrue) < 1e-15:
            break
        true_est -= residual / dapp_dtrue

    return max(0.0, apparent_alt - true_est)


# ---------------------------------------------------------------------------
# Reference-matched analytic refraction (used by refrac_extended / azalt only)
# ---------------------------------------------------------------------------
#
# The ray-traced model above is the library's default refraction (plain
# refrac(), via calc_refraction_true_to_app). refrac_extended() and azalt(),
# however, must be 1:1 with the reference API. Black-box probing of the
# reference's refrac_extended over a dense grid (input altitude x observer
# elevation x pressure x temperature) establishes its behaviour exactly:
#
#   * Refraction is a pure function of the *input* altitude, pressure and
#     temperature; it does NOT depend on the observer's elevation (that only
#     enters through the horizon dip -- see calc_dip below).
#   * The *apparent* altitude is the natural argument. As a function of the
#     apparent altitude the pressure scaling is separable and exactly linear at
#     every temperature: gP(P) = (P - 80) / 933.25. The app->true refraction is
#     therefore taken as the primitive here, and true->app is obtained by
#     inverting it. (The reference's own true->app and app->true agree to
#     <0.2" wherever its internal iteration is well-conditioned; below roughly
#     3 deg under the dipped horizon at strong cold that iteration diverges and
#     the reference returns non-physical values -- we return the physically
#     sensible clamped input there instead.)
#   * Base curve R0(app) at 1013.25 mbar, 15 C: a Bennett-family fit above the
#     horizon with a small rational correction, and a quartic continuation
#     below it (with a bounded linear tail past -2.3 deg). Matches the
#     reference to <2" over app in [-2.3, 90] deg.
#   * Temperature enters as (288.15 / (273.15 + T)) ** q(app, T). The exponent
#     q rises from ~0.88 at the zenith to ~1.75 at the horizon and ~2.9 a
#     couple of degrees below it -- an altitude-dependent sensitivity that no
#     single fixed exponent can capture (the previous model used a constant
#     1.31, correct only near the horizon). q is fit as a low-order polynomial
#     in the compressed altitude variable 1/(app + 5), quadratic in T.
#
# Measured agreement of refrac_extended vs the reference across observer
# elevation {0, 1000, 3000, 4000} m, T {-30 .. 35} C, P {900 .. 1100} mbar and
# input altitude from 4 deg below the dip to 1 deg above it: <30" everywhere
# (median <0.001"), excluding the divergent deep-below-horizon cold sliver
# noted above. See docs/comparison/known-differences.md section 13.
#
# All constants below were obtained purely by black-box probing (no reference
# source was consulted).

# Pressure scaling of the app->true refraction: exactly linear and separable.
_REFR_P_OFFSET: float = 80.0
_REFR_P_SPAN: float = 933.25  # = 1013.25 - 80

# Base curve R0(app) at 1013.25 mbar, 15 C, in degrees.
# Above the horizon: Bennett-family a/tan(app + b/(app+c)) plus a rational
# correction in w = 1/(app + 4) that pulls the fit under ~2".
_REFR_A: float = 0.9816287048243229
_REFR_B: float = 7.14
_REFR_C: float = 4.259
_REFR_CORR: tuple[float, float, float, float, float] = (
    0.007823842467183916,
    -1.3476822600397333,
    23.536935131196312,
    -134.2730344731844,
    245.7190217075024,
)
# Below the horizon: quartic in x = -app anchored at R0(0) for continuity,
# with a bounded linear tail past _REFR_APP_FLOOR.
_REFR_R0: float = 0.5599364492656921  # R0(app=0), degrees
_REFR_BELOW: tuple[float, float, float, float] = (
    0.22587654986933237,
    0.003876311801010285,
    0.07344146759158594,
    -0.021829746046986398,
)
_REFR_APP_FLOOR: float = -2.3

# Temperature exponent q(app, T): poly (deg 4) in v = 1/(app + _REFR_Q_K),
# quadratic in (T - 15).
_REFR_Q_K: float = 5.0
_REFR_Q0: tuple[float, float, float, float, float] = (
    0.863680955334398,
    2.2469424220791145,
    -10.346282707448633,
    146.4522078339461,
    -266.5962203788882,
)
_REFR_Q1: tuple[float, float, float, float] = (
    0.00022751617874395972,
    -0.003480924267185901,
    0.012420640141742743,
    -0.2297500939232252,
)
_REFR_Q2: tuple[float, float, float] = (
    2.343191507573941e-06,
    -0.00011627105492982002,
    0.0008754665131731234,
)


def _refr_base_curve(app: float) -> float:
    """Base app->true refraction R0(app) at 1013.25 mbar, 15 C, in degrees."""
    if app >= 0.0:
        w = 1.0 / (app + 4.0)
        r = _REFR_A / math.tan(math.radians(app + _REFR_B / (app + _REFR_C))) / 60.0
        c = _REFR_CORR
        return r + w * (c[0] + w * (c[1] + w * (c[2] + w * (c[3] + w * c[4]))))
    # below the horizon: quartic in x = -app, anchored at R0(0)
    p = _REFR_BELOW

    def _below(x: float) -> float:
        return _REFR_R0 + x * (p[0] + x * (p[1] + x * (p[2] + x * p[3])))

    if app >= _REFR_APP_FLOOR:
        return _below(-app)
    # bounded linear tail so the fit never blows up far below the horizon
    xf = -_REFR_APP_FLOOR
    val = _below(xf)
    slope = p[0] + xf * (2.0 * p[1] + xf * (3.0 * p[2] + xf * 4.0 * p[3]))
    return val + slope * (-app - xf)


def _refr_temp_exponent(app: float, temp_C: float) -> float:
    """Temperature exponent q(app, T) for the (288.15/(273.15+T))**q scaling."""
    v = 1.0 / (max(app, _REFR_APP_FLOOR) + _REFR_Q_K)
    s = temp_C - 15.0
    a = _REFR_Q0
    q0 = a[0] + v * (a[1] + v * (a[2] + v * (a[3] + v * a[4])))
    b = _REFR_Q1
    q1 = b[0] + v * (b[1] + v * (b[2] + v * b[3]))
    d = _REFR_Q2
    q2 = d[0] + v * (d[1] + v * d[2])
    return q0 + s * q1 + s * s * q2


def calc_refraction_ref_app_to_true(
    apparent_alt: float, pressure: float = _P0, temperature_C: float = 15.0
) -> float:
    """Reference-matched refraction (degrees) for an apparent altitude.

    This is the primitive of the reference-matched pair: as a function of the
    *apparent* altitude the pressure/temperature dependence is separable
    (see the module notes above). Valid above and below the geometric horizon.
    Used by refrac_extended()/azalt().
    """
    if pressure <= 0 or 273.15 + temperature_C <= 0:
        return 0.0
    scale_p = (pressure - _REFR_P_OFFSET) / _REFR_P_SPAN
    ratio_t = 288.15 / (273.15 + temperature_C)
    r = (
        scale_p
        * _refr_base_curve(apparent_alt)
        * ratio_t ** _refr_temp_exponent(apparent_alt, temperature_C)
    )
    return max(0.0, r)


def calc_refraction_ref_true_to_app(
    true_alt: float, pressure: float = _P0, temperature_C: float = 15.0
) -> float:
    """Reference-matched refraction (degrees) for a true altitude.

    Inverts calc_refraction_ref_app_to_true: finds the apparent altitude *app*
    with ``app - R_app_to_true(app) == true_alt`` and returns ``app - true_alt``.
    The mapping ``app - R(app)`` is strictly increasing, so a bracketed
    bisection on ``[true_alt, true_alt + 8]`` converges to the unique root.
    Valid above and below the geometric horizon; used by refrac_extended().
    """
    if pressure <= 0 or 273.15 + temperature_C <= 0:
        return 0.0
    lo = true_alt
    hi = true_alt + 8.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        f = (
            mid
            - calc_refraction_ref_app_to_true(mid, pressure, temperature_C)
            - true_alt
        )
        if f < 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-12:
            break
    return 0.5 * (lo + hi) - true_alt


def calc_dip(
    obs_alt: float,
    lapse_rate: float = 0.0065,
    atpress: float = 1013.25,
    attemp: float = 10.0,
) -> float:
    """Dip of the horizon for an elevated observer.

    Uses A. Thom's refraction-corrected dip: the geometric dip
    arccos(R/(R+h)) is scaled by sqrt(1 - 1.848*k*P/(273.16+T)^2)
    with k = (0.0342 + lapse)/(0.154*0.0238). The atmospheric scaling
    matters: at 1000 m the dip moves by ~190" between (700 mbar, 30 degC)
    and standard conditions.

    Parameters
    ----------
    obs_alt : float
        Observer altitude above sea level in metres.
    lapse_rate : float
        Tropospheric lapse rate in K/m.
    atpress : float
        Atmospheric pressure at the observer in mbar.
    attemp : float
        Atmospheric temperature at the observer in degrees Celsius.

    Returns
    -------
    float
        Dip angle in degrees (negative: horizon is below the
        geometric horizontal).

    References
    ----------
    Thom, A. (1971), "Megalithic Lunar Observatories", Oxford, ch. 3
    Bomford, "Geodesy" (1980), 4th ed., §2.17-2.20
    """
    if obs_alt <= 0:
        return 0.0

    # Equatorial radius for the geometric dip (the conventional choice in
    # the dip literature; the mean radius shifts the dip ~2" at 1000 m)
    _r_dip = 6_378_136.6
    ratio = _r_dip / (_r_dip + obs_alt)
    if ratio >= 1.0:
        return 0.0
    dip_geometric = math.degrees(math.acos(ratio))

    # Refraction coefficient after Thom (1971). 0.0342 K/m is the
    # autoconvective lapse rate g/R_specific (≈9.81/287) at which density is
    # height-independent — the upper bound on the actual lapse rate; the
    # 0.154 and 0.0238 factors are Thom's empirical normalisation. The
    # lapse-rate term keeps the conventional 0.0065 K/m default consistent
    # with the standard atmosphere used elsewhere in this module.
    krefr = (0.0342 + lapse_rate) / (0.154 * 0.0238)
    t_kelvin = 273.16 + attemp
    # Non-physical temperature at or below absolute zero (attemp <= -273.16 C):
    # there is no meaningful atmosphere, and t_kelvin**2 would divide by zero.
    # Return zero dip, consistent with the obs_T_K <= 0 guard in the ray tracer.
    if t_kelvin <= 0.0:
        return 0.0
    d = 1.0 - 1.8480 * krefr * atpress / (t_kelvin * t_kelvin)
    if d < 0.0:
        d = 0.0

    return -dip_geometric * math.sqrt(d)
