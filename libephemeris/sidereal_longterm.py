# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Long-term sidereal time, precession and mean obliquity (Vondrák 2011).

This module is the single source, across the whole library, of:

* the **of-date mean obliquity** of the ecliptic, and
* the **precession matrix** (J2000 mean frame → mean equator/equinox of date), and
* **Greenwich mean / apparent sidereal time** (the ARMC that house systems need),

all evaluated with a model that stays correct over very long time spans
(≈ ±200 millennia), rather than only a few centuries around J2000.0.

Why a long-term model (and why it beats the usual choice)
---------------------------------------------------------
The Ascendant, Midheaven and all house cusps are obtained from the **ARMC**
(the right ascension of the meridian, i.e. apparent sidereal time at the
observer's longitude) together with the obliquity of the ecliptic. Both of
these quantities are driven by **precession** — the slow (~25 800 yr) wobble of
the Earth's axis that moves the equinox.

The classical IAU 1976/2006 precession is expressed as a **polynomial in time**
fitted near J2000.0. A truncated polynomial is excellent for a few centuries but
**diverges** rapidly outside its fit window: by ±8 000 years the implied sidereal
time and obliquity are wrong by **degrees**, which corrupts every house cusp at
historical or far-future epochs. Many astrology engines inherit exactly this
limitation because they take sidereal time from an IAU-1976/2006 routine.

We instead use the long-term precession of

    Vondrák, J., Capitaine, N. & Wallace, P. (2011),
    "New precession expressions, valid for long time intervals",
    Astronomy & Astrophysics 534, A22.  DOI: 10.1051/0004-6361/201117274
    (with the corrigendum A&A 541, C1, 2012).

It is fitted to a numerical integration and stays accurate over ±200 000 years,
while agreeing with IAU 2006 to sub-milliarcsecond near J2000.0. Adopting it:

* makes ancient/future house cusps **scientifically correct** (no polynomial
  blow-up), and
* leaves modern results **unchanged** (sub-mas agreement near J2000.0).

The precession rotation is built from the Vondrák ecliptic-pole and equator-pole
series (``PQ*`` and ``XY*``), and the of-date mean obliquity is the **angle
between those same two poles**, so precession and obliquity are one
self-consistent realization: a direction lying in the mean ecliptic of date
(e.g. the Sun) reduces to zero mean-ecliptic latitude, and cusps and bodies in
one chart share a single frame. See :func:`mean_obliquity_rad`.

Vondrák 2011 also publishes a *direct* obliquity series (the ``p_A``/``ε_A``
polynomial+periodic terms, transcribed in ``_PEPOL``/``_PEPER`` below and
evaluated by :func:`mean_obliquity_series_rad`). That direct series is a very
good obliquity fit on its own (it tracks the IAU 2006 polynomial to < 1 mas near
J2000, better than the pole angle) and is the of-date obliquity the reference
engine reports. It is **not** used for the of-date ecliptic rotation, because it
is a *separate* fit from the pole series: pairing it with the pole-based
precession matrix tilts the of-date ecliptic away from the pole-defined ecliptic
by up to ~6.5″ at −3000 (~17.6″ at −5000), which shows up as a spurious ecliptic
latitude on the Sun. The pole-angle obliquity removes that inconsistency; the
two realizations agree to < 0.001″ across 1900-2100 (identically 0 at J2000).

How sidereal time is computed (the geometric method)
----------------------------------------------------
Sidereal time is the hour angle of the equinox, i.e. Earth-Rotation-Angle plus
the accumulated precession in right ascension. Evaluating the *precession-in-RA*
as a polynomial diverges at remote epochs for the same reason as above. We avoid
that entirely by a **geometric** construction valid everywhere:

1. take the mean longitude of the Earth from the published secular expression of
   Simon, J.L. et al. (1994), "Numerical expressions for precession formulae and
   mean elements for the Moon and the planets", A&A 282, 663;
2. correct it for Sun-Earth light-time;
3. form the corresponding unit direction on the ecliptic of J2000.0, rotate it to
   the J2000.0 equator, **precess it to the date with the Vondrák matrix**, and
   read off its ecliptic longitude of date;
4. add the equation of the equinoxes (nutation in longitude x cos eps) and the UT1
   hour angle.

Because step 3 transforms a *physical direction* through the long-term precession
matrix (instead of summing a divergent RA polynomial), the result is stable over
the whole supported range. For the modern window 1850-2050 the well-established
IAU 2006 GMST polynomial is used directly (it is most precise there), with tiny
continuity offsets so the two branches join smoothly.

Time scales: precession and obliquity use **TT**; the Earth-rotation hour angle
uses **UT1**. The TT↔UT1 difference (ΔT) is obtained from the library's own ΔT
model (:func:`libephemeris.time_utils.deltat`), so houses and planetary positions
use one consistent ΔT.

Provenance
----------
All coefficients are transcribed from the cited peer-reviewed papers (Vondrák
2011 + corrigendum A&A 541 C1; Simon 1994) and the IAU 2006 GMST expression
(Capitaine, Wallace & Chapront 2003, A&A 412, 567). No third-party source code
is used. The series and the geometric construction are standard published
astronomy.

Comparison vs the reference ephemeris (matched ΔT)
--------------------------------------------------
Inside 1850-2050 the ARMC matches the reference ephemeris to ~0.002" (identical
IAU-2006 GMST branch). Outside that window the two engines use *different*
long-term sidereal-time realizations, so house cusps (ASC/MC) diverge by a
secular, sign-changing amount even at equal ΔT (ASC ≈ +1.5" at 2050, +0.9" at
2100, -1.3" at 2200, -5.6" at 2300, ~0.1-1° at ±3000-5000 yr). This is a benign
*model* difference: obliquity and nutation match the reference to < 0.002" at
every epoch, so the residual is entirely the long-term ARMC model, and
libephemeris's geometric Vondrák construction is the more physically correct of
the two (the reference continues an IAU-2006-style precession-in-RA polynomial
that itself diverges at remote epochs). NB: the ~1.9" step seen near 2050 is a
discontinuity on the reference side (its sidereal time has a one-time ~1.908"
jump exactly at JD 2469807.5); libephemeris's two branches join to 0.000000"
there, i.e. it is self-continuous. See
``docs/methodology/sidereal-time-longterm.md`` and
``docs/comparison/known-differences.md`` for the full analysis.
"""

from __future__ import annotations

import math
from functools import lru_cache

# --------------------------------------------------------------------------
# constants
# --------------------------------------------------------------------------
_J2000 = 2451545.0
_AS2R = (math.pi / 180.0) / 3600.0  # arcseconds -> radians
_D2PI = 2.0 * math.pi
_DEG = math.pi / 180.0

# Obliquity of the ecliptic at J2000.0 (IAU 2006): 84381.406".
_EPS0 = 84381.406 * _AS2R

# Sun-Earth light-time, in days: AU / c. AU = 1.495978707e11 m (IAU 2012),
# c = 299792458 m/s, 86400 s/day  ->  ~0.0057755183 day.
_LIGHT_TIME_DAYS = 1.495978707e11 / 299792458.0 / 86400.0

# Modern window where the IAU 2006 GMST polynomial branch is used directly.
_LTERM_T0 = 2396758.5  # 1 Jan 1850
_LTERM_T1 = 2469807.5  # 1 Jan 2050
# Continuity offsets that join the long-term branch onto the modern polynomial
# branch at the two boundaries are computed at runtime from the live ΔT model
# (see _lterm_offset below), not hardcoded — so they never go stale when the ΔT
# model changes (set_delta_t_userdef / IERS). They are only evaluated on the
# rare long-term branch (dates outside 1850-2050).

# --------------------------------------------------------------------------
# Vondrák 2011 (A&A 534, A22; corrigendum A&A 541, C1) coefficient series.
# Units: arcseconds (polynomials) / arcseconds and years (periodic terms).
# --------------------------------------------------------------------------

# General precession in longitude / obliquity — polynomial part.
_PEPOL = (
    (8134.017132, 84028.206305),
    (5043.0520035, 0.3624445),
    (-0.00710733, -0.00004039),
    (0.000000271, -0.000000110),
)
# General precession / obliquity — periodic part: [period, p_cos, q_cos, p_sin, q_sin].
_PEPER = (
    (409.90, 396.15, 537.22, 402.90, 417.15, 288.92, 4043.00, 306.00, 277.00, 203.00),
    (-6908.287473, -3198.706291, 1453.674527, -857.748557, 1173.231614,
     -156.981465, 371.836550, -216.619040, 193.691479, 11.891524),
    (753.872780, -247.805823, 379.471484, -53.880558, -90.109153,
     -353.600190, -63.115353, -28.248187, 17.703387, 38.911307),
    (-2845.175469, 449.844989, -1255.915323, 886.736783, 418.887514,
     997.912441, -240.979710, 76.541307, -36.788069, -170.964086),
    (-1704.720302, -862.308358, 447.832178, -889.571909, 190.402846,
     -56.564991, -296.222622, -75.859952, 67.473503, 3.014055),
)

# Ecliptic pole P, Q — polynomial part.
_PQPOL = (
    (5851.607687, -1600.886300),
    (-0.1189000, 1.1689818),
    (-0.00028913, -0.00000020),
    (0.000000101, -0.000000437),
)
# Ecliptic pole P, Q — periodic part. Index 6 of the q_cos row carries the
# A&A 541, C1 (2012) corrigendum value 198.296701 (the original printing had a typo).
_PQPER = (
    (708.15, 2309.0, 1620.0, 492.2, 1183.0, 622.0, 882.0, 547.0),
    (-5486.751211, -17.127623, -617.517403, 413.44294, 78.614193,
     -180.732815, -87.676083, 46.140315),
    (-684.66156, 2446.28388, 399.671049, -356.652376, -186.387003,
     -316.80007, 198.296701, 101.135679),
    (667.66673, -2354.886252, -428.152441, 376.202861, 184.778874,
     335.321713, -185.138669, -120.97283),
    (-5523.863691, -549.74745, -310.998056, 421.535876, -36.776172,
     -145.278396, -34.74445, 22.885731),
)

# Equator pole X, Y — polynomial part.
_XYPOL = (
    (5453.282155, -73750.930350),
    (0.4252841, -0.7675452),
    (-0.00037173, -0.00018725),
    (-0.000000152, 0.000000231),
)
# Equator pole X, Y — periodic part.
_XYPER = (
    (256.75, 708.15, 274.2, 241.45, 2309.0, 492.2, 396.1, 288.9, 231.1,
     1610.0, 620.0, 157.87, 220.3, 1200.0),
    (-819.940624, -8444.676815, 2600.009459, 2755.17563, -167.659835,
     871.855056, 44.769698, -512.313065, -819.415595, -538.071099,
     -189.793622, -402.922932, 179.516345, -9.814756),
    (75004.344875, 624.033993, 1251.136893, -1102.212834, -2660.66498,
     699.291817, 153.16722, -950.865637, 499.754645, -145.18821,
     558.116553, -23.923029, -165.405086, 9.344131),
    (81491.287984, 787.163481, 1251.296102, -1257.950837, -2966.79973,
     639.744522, 131.600209, -445.040117, 584.522874, -89.756563,
     524.42963, -13.549067, -210.157124, -44.919798),
    (1558.515853, 7774.939698, -2219.534038, -2523.969396, 247.850422,
     -846.485643, -1393.124055, 368.526116, 749.045012, 444.704518,
     235.934465, 374.049623, -171.33018, -22.899655),
)


# --------------------------------------------------------------------------
# Vondrák series evaluators
# --------------------------------------------------------------------------
@lru_cache(maxsize=4096)
def mean_obliquity_series_rad(jd_tt: float) -> float:
    """Direct Vondrák 2011 obliquity series (``p_A``/``ε_A``) in radians.

    Evaluates the published ``_PEPOL``/``_PEPER`` obliquity polynomial+periodic
    terms directly. This is an independent Vondrák-2011 fit of the of-date mean
    obliquity, and it is the value the reference engine reports for the of-date
    obliquity (it also tracks the IAU 2006 obliquity to < 1 mas near J2000).

    It is retained for provenance and reference-parity comparisons, but is NOT
    used for the of-date ecliptic rotation: because it is a separate fit from the
    equator/ecliptic pole series that build the precession matrix, pairing it
    with that matrix is frame-inconsistent and puts a spurious latitude on the
    Sun at remote epochs. Use :func:`mean_obliquity_rad` (the pole angle) for any
    equator↔ecliptic-of-date transform.

    Args:
        jd_tt: Julian Date in TT.

    Returns:
        The of-date mean obliquity in radians (direct series).
    """
    t = (jd_tt - _J2000) / 36525.0
    q = 0.0
    periods, _pc, q_cos, _ps, q_sin = _PEPER
    for i in range(len(periods)):
        a = _D2PI * t / periods[i]
        q += math.cos(a) * q_cos[i] + math.sin(a) * q_sin[i]
    w = 1.0
    for poly in _PEPOL:
        q += poly[1] * w
        w *= t
    return q * _AS2R


@lru_cache(maxsize=4096)
def mean_obliquity_rad(jd_tt: float) -> float:
    """Of-date mean obliquity of the ecliptic in radians (Vondrák 2011).

    Long-term-valid replacement for the IAU 2006 obliquity polynomial; this is
    the single obliquity realization used by both the houses path and the
    position pipeline.

    It is the **angle between the of-date ecliptic pole and equator pole**
    (:func:`_ecliptic_pole`, :func:`_equator_pole`) — the same two Vondrák pole
    vectors that build :func:`precession_matrix` (and ``erfa.ltp``/``ltpb`` on the
    position pipeline). Deriving the obliquity from those poles guarantees the
    equator-of-date → ecliptic-of-date rotation is frame-consistent with the
    precession: a direction lying in the mean ecliptic of date reduces to zero
    mean-ecliptic latitude. Using the independently-fitted direct series
    (:func:`mean_obliquity_series_rad`) here instead would tilt the of-date
    ecliptic by up to ~6.5″ at −3000 and put a spurious latitude on the Sun.

    Near J2000 the pole angle and the direct series agree to < 0.001″ (identically
    0 at J2000), so modern obliquity/positions/houses are unchanged.

    Args:
        jd_tt: Julian Date in TT.

    Returns:
        The of-date mean obliquity in radians (pole angle).
    """
    ex, ey, ez = _ecliptic_pole(jd_tt)
    qx, qy, qz = _equator_pole(jd_tt)
    d = ex * qx + ey * qy + ez * qz
    if d > 1.0:
        d = 1.0
    elif d < -1.0:
        d = -1.0
    return math.acos(d)


def mean_obliquity_deg(jd_tt: float) -> float:
    """Of-date mean obliquity of the ecliptic in degrees (Vondrák 2011)."""
    return math.degrees(mean_obliquity_rad(jd_tt))


@lru_cache(maxsize=4096)
def _ecliptic_pole(jd_tt: float) -> tuple[float, float, float]:
    """Unit vector of the ecliptic pole at ``jd_tt`` (Vondrák 2011)."""
    t = (jd_tt - _J2000) / 36525.0
    p = 0.0
    q = 0.0
    periods, p_cos, q_cos, p_sin, q_sin = _PQPER
    for i in range(len(periods)):
        a = _D2PI * t / periods[i]
        s = math.sin(a)
        c = math.cos(a)
        p += c * p_cos[i] + s * p_sin[i]
        q += c * q_cos[i] + s * q_sin[i]
    w = 1.0
    for poly in _PQPOL:
        p += poly[0] * w
        q += poly[1] * w
        w *= t
    p *= _AS2R
    q *= _AS2R
    z = 1.0 - p * p - q * q
    z = math.sqrt(z) if z > 0.0 else 0.0
    se = math.sin(_EPS0)
    ce = math.cos(_EPS0)
    return (p, -q * ce - z * se, -q * se + z * ce)


@lru_cache(maxsize=4096)
def _equator_pole(jd_tt: float) -> tuple[float, float, float]:
    """Unit vector of the equator pole (CIP) at ``jd_tt`` (Vondrák 2011)."""
    t = (jd_tt - _J2000) / 36525.0
    x = 0.0
    y = 0.0
    periods, x_cos, y_cos, x_sin, y_sin = _XYPER
    for i in range(len(periods)):
        a = _D2PI * t / periods[i]
        s = math.sin(a)
        c = math.cos(a)
        x += c * x_cos[i] + s * x_sin[i]
        y += c * y_cos[i] + s * y_sin[i]
    w = 1.0
    for poly in _XYPOL:
        x += poly[0] * w
        y += poly[1] * w
        w *= t
    x *= _AS2R
    y *= _AS2R
    w = x * x + y * y
    z = math.sqrt(1.0 - w) if w < 1.0 else 0.0
    return (x, y, z)


def _cross(a: tuple, b: tuple) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


@lru_cache(maxsize=4096)
def precession_matrix(jd_tt: float) -> tuple:
    """Precession matrix: J2000 mean frame → mean equator/equinox of date.

    Built from the Vondrák ecliptic-pole and equator-pole vectors. Returned as a
    9-tuple in row-major order (r00, r01, r02, r10, r11, r12, r20, r21, r22).
    Multiplying a J2000-equator vector ``v`` by the ROWS gives the of-date vector
    (``v_date[i] = sum_j R[i][j] v[j]``).

    Args:
        jd_tt: Julian Date in TT.
    """
    peqr = _equator_pole(jd_tt)
    pecl = _ecliptic_pole(jd_tt)
    v = _cross(peqr, pecl)
    w = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    eqx = (v[0] / w, v[1] / w, v[2] / w)
    mid = _cross(peqr, eqx)
    return (
        eqx[0], eqx[1], eqx[2],
        mid[0], mid[1], mid[2],
        peqr[0], peqr[1], peqr[2],
    )


def _precess_j2000_to_date(vec: tuple, jd_tt: float) -> tuple[float, float, float]:
    """Rotate a J2000-equator vector to the mean equator of date."""
    r = precession_matrix(jd_tt)
    return (
        r[0] * vec[0] + r[1] * vec[1] + r[2] * vec[2],
        r[3] * vec[0] + r[4] * vec[1] + r[5] * vec[2],
        r[6] * vec[0] + r[7] * vec[1] + r[8] * vec[2],
    )


def _rot_x(vec: tuple, eps: float) -> tuple[float, float, float]:
    """Rotate ``vec`` about the x-axis by ``eps`` (radians).

    Convention: ecliptic→equator uses ``-eps``; equator→ecliptic uses ``+eps``.
    """
    c = math.cos(eps)
    s = math.sin(eps)
    return (vec[0], vec[1] * c + vec[2] * s, -vec[1] * s + vec[2] * c)


# --------------------------------------------------------------------------
# Sidereal time
# --------------------------------------------------------------------------
def _deltat_days(jd_ut1: float) -> float:
    # Imported lazily to avoid a module import cycle (time_utils imports constants
    # but not this module; keeping the import local is defensive and cheap).
    from .time_utils import deltat

    return deltat(jd_ut1)


def _gmst_iau2006_deg(jd_ut1: float, jd_tt: float) -> float:
    """Greenwich Mean Sidereal Time (deg) via the IAU 2006 expression.

    GMST = ERA(UT1) + p(T), with the precession-in-RA polynomial p(T) of
    Capitaine, Wallace & Chapront (2003), A&A 412, 567. Used in the modern
    window where it is most precise.
    """
    du = jd_ut1 - _J2000
    era = _D2PI * (0.7790572732640 + 1.00273781191135448 * du)
    t = (jd_tt - _J2000) / 36525.0
    p = (
        0.014506
        + t * (4612.156534
               + t * (1.3915817
                      + t * (-0.00000044
                             + t * (-0.000029956
                                    + t * -0.0000000368))))
    )
    return math.degrees(era) + p / 3600.0


def _mean_sidereal_longterm_deg(jd_ut1: float) -> float:
    """Greenwich mean sidereal time (deg) via the long-term geometric method."""
    jd_tt = jd_ut1 + _deltat_days(jd_ut1)
    t = (jd_tt - _J2000) / 365250.0  # Julian millennia (TT)
    t2 = t * t
    t3 = t * t2
    # Simon et al. 1994 mean longitude of the Earth (degrees).
    dlon = 100.46645683 + (1295977422.83429 * t - 2.04411 * t2 - 0.00523 * t3) / 3600.0
    # Sun-Earth light-time correction.
    dlon = (dlon - _LIGHT_TIME_DAYS * 360.0 / 365.2425) % 360.0
    # Ecliptic-of-J2000 unit direction, then to the J2000 equator.
    rad = dlon * _DEG
    vec = (math.cos(rad), math.sin(rad), 0.0)
    # _J2000 (2451545.0) is already the J2000.0 TT epoch and mean_obliquity_rad
    # takes TT, so it is passed directly — adding ΔT would double-count the
    # UT→TT offset (the eps_date term below correctly uses jd_tt as-is).
    eps_j2000 = mean_obliquity_rad(_J2000)
    vec = _rot_x(vec, -eps_j2000)
    # Precess to the mean equator of date, then back to the ecliptic of date.
    vec = _precess_j2000_to_date(vec, jd_tt)
    eps_date = mean_obliquity_rad(jd_tt)
    vec = _rot_x(vec, eps_date)
    lon = math.degrees(math.atan2(vec[1], vec[0]))
    # UT1 hour angle.
    dhour = math.fmod(jd_ut1 - 0.5, 1.0) * 360.0
    return (lon + dhour) % 360.0


def _lterm_offset(boundary_jd: float) -> float:
    """Continuity offset (deg) joining the long-term branch to the modern one.

    Equals (long-term branch - modern branch) evaluated at ``boundary_jd``, so
    the piecewise GMST is continuous there. Computed from the live delta-T model
    rather than hardcoded, so it can never drift out of sync when the delta-T
    model changes. Only called on the long-term branch (rare), so the extra cost
    is negligible.

    Args:
        boundary_jd: Julian Day (UT1) of the branch boundary (1 Jan 1850 or
            1 Jan 2050).

    Returns:
        Continuity offset in degrees, reduced to (-180, 180] to absorb any
        0/360 wrap between the two branches.
    """
    jd_tt = boundary_jd + _deltat_days(boundary_jd)
    diff = _mean_sidereal_longterm_deg(boundary_jd) - (
        _gmst_iau2006_deg(boundary_jd, jd_tt) % 360.0
    )
    return (diff + 180.0) % 360.0 - 180.0


def mean_sidereal_time_deg(jd_ut1: float) -> float:
    """Greenwich Mean Sidereal Time at ``jd_ut1`` in degrees [0, 360).

    Piecewise: the IAU 2006 polynomial in 1850-2050 (most precise there), the
    long-term geometric method outside it (correct over the whole supported
    range), joined with small continuity offsets.
    """
    if _LTERM_T0 < jd_ut1 < _LTERM_T1:
        jd_tt = jd_ut1 + _deltat_days(jd_ut1)
        return _gmst_iau2006_deg(jd_ut1, jd_tt) % 360.0
    g = _mean_sidereal_longterm_deg(jd_ut1)
    if jd_ut1 <= _LTERM_T0:
        g -= _lterm_offset(_LTERM_T0)
    else:
        g -= _lterm_offset(_LTERM_T1)
    return g % 360.0


def apparent_sidereal_time_deg(
    jd_ut1: float, longitude: float = 0.0, *, dpsi_deg: float, eps_true_deg: float
) -> float:
    """Local Apparent Sidereal Time at ``jd_ut1`` in degrees [0, 360).

    This is the ARMC used by house systems when ``longitude`` is the observer's
    geographic longitude. Apparent = mean + equation of the equinoxes
    (``dpsi · cos ε_true``), with nutation supplied by the caller so that the
    library's tuned nutation series is reused unchanged.

    Args:
        jd_ut1: Julian Date in UT1.
        longitude: Geographic longitude in degrees (East positive). Default 0
            gives Greenwich apparent sidereal time (= ARMC at longitude 0).
        dpsi_deg: Nutation in longitude (degrees).
        eps_true_deg: True obliquity of date (degrees).

    Returns:
        Local apparent sidereal time in degrees [0, 360).
    """
    gmst = mean_sidereal_time_deg(jd_ut1)
    eoe = dpsi_deg * math.cos(math.radians(eps_true_deg))
    return (gmst + eoe + longitude) % 360.0
