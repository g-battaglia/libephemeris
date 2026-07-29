# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Atmospheric-refraction models for the public and physical APIs.

Plain :func:`libephemeris.refrac` uses independently published closed-form
approximations. :func:`libephemeris.refrac_extended` and
:func:`libephemeris.azalt` use the ICAO Standard Atmosphere ray tracer in this
module.

The physical model does not rely on empirical curve fits. It computes
refraction from first principles by evaluating the exact refraction integral
through a continuously varying model atmosphere.

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

Closed-form public API
----------------------
The compact ``refrac()`` implementation uses the independently published
Saemundsson true-to-apparent and Bennett apparent-to-true approximations, with
the conventional pressure and temperature scaling. These formulas remain
separate from the ray tracer because callers may depend on their historical
closed-form behavior. The true-to-apparent branch is Saemundsson's formula
(Meeus 1998, ch. 16, eq. 16.4, p. 106) below 15 deg with the high-altitude
tangent series above it; the apparent-to-true branch is Bennett's formula
(Meeus 1998, ch. 16, eq. 16.3, p. 106). The two directions are independent
published fits and are therefore not exact numerical inverses of each other.

References
----------
- Bennett, G.G. (1982), *Journal of Navigation*, 35, 255-259
- Saemundsson, Th. (1986), *Sky & Telescope*, 72, 70
- Meeus, J. (1998), *Astronomical Algorithms*, 2nd ed., chapter 16
  (Atmospheric Refraction), pp. 105-108; formulas 16.3 (Bennett,
  apparent->true) and 16.4 (Saemundsson, true->apparent)

Provenance:
    Standard-atmosphere layers, gas laws, refractivity, and closed-form public
    formulas come from the sources cited in their respective sections. The
    spherical ray integral, Gauss-Legendre evaluation, inverse Newton solve,
    quadrature order, coefficient rounding, and guard domains are independently
    documented project engineering; convergence and horizon error are measured.
    Source constants and project-selected numerical parameters are kept in
    separate named tables and are never fitted to compatibility output.
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
_R_GAS: float = 8.314_462_618  # Universal gas constant, CODATA [J/(mol*K)]
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

# Conventional pressure/temperature scaling used with the independently
# published Saemundsson and Bennett closed forms (Meeus 1998, chapter 16).
_COMPAT_PRESSURE_REF: float = 1010.0
_COMPAT_TEMPERATURE_NUMERATOR: float = 283.0
_COMPAT_TEMPERATURE_OFFSET: float = 273.0
_COMPAT_TRUE_LOW_LIMIT: float = -5.0
_COMPAT_TRUE_BRANCH: float = 15.0
_COMPAT_REVERSE_ZENITH_LIMIT: float = 90.0 - 1e-10


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
# Plain refrac() clean-room compatibility model
# ---------------------------------------------------------------------------


def _compat_atmosphere_scale(pressure: float, temperature_C: float) -> float:
    """Return the conventional pressure/temperature refraction scale.

    The explicit zero-denominator handling keeps physically invalid edge
    inputs deterministic instead of raising Python's ``ZeroDivisionError``.
    """
    numerator = pressure / _COMPAT_PRESSURE_REF * _COMPAT_TEMPERATURE_NUMERATOR
    denominator = _COMPAT_TEMPERATURE_OFFSET + temperature_C
    if denominator != 0.0:
        return numerator / denominator
    if numerator == 0.0 or math.isnan(numerator):
        return math.nan
    return math.copysign(math.inf, numerator)


def calc_refrac_compat_true_to_app(
    true_alt: float, pressure: float = _P0, temperature_C: float = 15.0
) -> float:
    """Return plain-API apparent altitude for a true altitude.

    This uses the independently published Saemundsson approximation below 15
    degrees and the standard high-altitude tangent series above it. It is not
    the ICAO ray tracer. Values at or below -5 degrees, and candidates that do
    not rise strictly above zero, are returned unchanged.
    """
    true_alt = float(true_alt)
    pressure = float(pressure)
    temperature_C = float(temperature_C)

    if not math.isfinite(true_alt) or true_alt <= _COMPAT_TRUE_LOW_LIMIT:
        return true_alt

    scale = _compat_atmosphere_scale(pressure, temperature_C)
    if true_alt <= _COMPAT_TRUE_BRANCH:
        angle = true_alt + 10.3 / (true_alt + 5.11)
        correction = 1.02 / math.tan(math.radians(angle)) * scale / 60.0
    else:
        tangent = math.tan(math.radians(90.0 - true_alt))
        correction = (
            (58.276 * tangent - 0.0824 * tangent * tangent * tangent) * scale / 3600.0
        )

    apparent_alt = true_alt + correction
    return float(apparent_alt) if apparent_alt > 0.0 else true_alt


def calc_refrac_compat_app_to_true(
    apparent_alt: float, pressure: float = _P0, temperature_C: float = 15.0
) -> float:
    """Return plain-API true altitude for an apparent altitude.

    The reverse direction uses Bennett's independently published closed form
    rather than numerically inverting :func:`calc_refrac_compat_true_to_app`.
    """
    apparent_alt = float(apparent_alt)
    pressure = float(pressure)
    temperature_C = float(temperature_C)

    if not math.isfinite(apparent_alt):
        return apparent_alt

    try:
        angle = apparent_alt + 7.31 / (apparent_alt + 4.4)
    except ZeroDivisionError:
        return apparent_alt
    if not math.isfinite(angle) or angle >= _COMPAT_REVERSE_ZENITH_LIMIT:
        return apparent_alt

    tangent = math.tan(math.radians(angle))
    if tangent == 0.0:
        return apparent_alt
    correction_arcmin = 1.0 / tangent
    correction = (
        correction_arcmin * _compat_atmosphere_scale(pressure, temperature_C) / 60.0
    )
    true_alt = apparent_alt - correction
    return float(true_alt) if true_alt > 0.0 else apparent_alt


# ---------------------------------------------------------------------------
# Independent ICAO ray-traced functions
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
# Legacy internal aliases
# ---------------------------------------------------------------------------
#
# Older internal callers used ``calc_refraction_ref_*`` for the extended API.
# Keep those names as thin aliases to the independently implemented ICAO ray
# tracer.


def calc_refraction_ref_app_to_true(
    apparent_alt: float,
    pressure: float = _P0,
    temperature_C: float = 15.0,
    obs_alt: float = 0.0,
    lapse_rate: float = 0.0065,
) -> float:
    """Return ICAO ray-traced refraction for an apparent altitude.

    This compatibility alias is retained for internal callers; it delegates
    directly to :func:`calc_refraction_app_to_true`.
    """
    return calc_refraction_app_to_true(
        apparent_alt, pressure, temperature_C, obs_alt, lapse_rate
    )


def calc_refraction_ref_true_to_app(
    true_alt: float,
    pressure: float = _P0,
    temperature_C: float = 15.0,
    obs_alt: float = 0.0,
    lapse_rate: float = 0.0065,
) -> float:
    """Return ICAO ray-traced refraction for a true altitude.

    This compatibility alias is retained for internal callers; it delegates
    directly to :func:`calc_refraction_true_to_app`.
    """
    return calc_refraction_true_to_app(
        true_alt, pressure, temperature_C, obs_alt, lapse_rate
    )


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
