# SPDX-License-Identifier: AGPL-3.0-only
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

The model is evaluated by **ERFA** (pyerfa), the IAU SOFA-derived library that
is already a hard dependency of this package, rather than by coefficient tables
kept in this module:

* ``erfa.ltpecl`` — the long-term ecliptic pole of date,
* ``erfa.ltpequ`` — the long-term equator pole (CIP) of date,
* ``erfa.ltp``    — the precession matrix built from those two poles
  (J2000.0 mean equator/equinox → mean equator/equinox of date).

The public of-date mean obliquity is the angle between those same two poles, so
the precession matrix and every equator/ecliptic rotation share one
self-consistent frame. :mod:`libephemeris.precession_vondrak` takes its
of-date mean obliquity from here and its matrices from the same ERFA routines
(``erfa.ltp``/``erfa.ltpb``), so one chart's bodies and angles sit in one frame.

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

Provenance:
    The Vondrák 2011 poles and precession matrix are evaluated by ERFA
    (``erfa.ltpecl``, ``erfa.ltpequ``, ``erfa.ltp``), the IAU SOFA-derived
    library; no coefficient of that model is kept in this module. The mean
    longitude of the Earth is transcribed from Simon et al. (1994) and the
    modern-window GMST from the IAU 2006 expression (Capitaine, Wallace &
    Chapront 2003, A&A 412, 567). The geometric construction is standard
    published astronomy. The modern/long-term branch boundary, continuity
    offsets, the refusal of an infinite epoch, caching, and floating-point
    evaluation order are project choices; their derivation and tests are
    recorded in ``docs/methodology/sidereal-time-longterm.md``.
"""

from __future__ import annotations

import math
from functools import lru_cache

import erfa

# --------------------------------------------------------------------------
# constants
# --------------------------------------------------------------------------
_J2000 = 2451545.0
_D2PI = 2.0 * math.pi
_DEG = math.pi / 180.0

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
# Vondrák 2011 through ERFA: poles, obliquity, precession matrix
# --------------------------------------------------------------------------
def _julian_epoch(jd_tt: float) -> float:
    """Convert a Julian Date (TT) to the Julian epoch ERFA's ``ltp*`` take.

    Args:
        jd_tt: Julian Date in TT.

    Returns:
        The Julian epoch (2000.0 at J2000.0).

    Raises:
        ValueError: if ``jd_tt`` is infinite. An infinite epoch has no place
            on the model's time axis; it is refused the way the transcendental
            functions refuse it (``math.sin(inf)`` raises this same error)
            instead of letting ERFA answer NaN silently. NaN passes through
            and yields NaN, as every other arithmetic path here does.
    """
    if math.isinf(jd_tt):
        raise ValueError("math domain error")
    return 2000.0 + (jd_tt - _J2000) / 365.25


@lru_cache(maxsize=4096)
def mean_obliquity_rad(jd_tt: float) -> float:
    """Return the public of-date mean obliquity in radians.

    The canonical value is the angle between the Vondrák ecliptic pole
    (``erfa.ltpecl``) and equator pole (``erfa.ltpequ``) — the same two poles
    ERFA builds the long-term precession matrix from. Deriving the obliquity
    from those poles makes the precession and every equator↔ecliptic-of-date
    rotation one self-consistent frame: a direction lying in the mean ecliptic
    of date (the Sun, by definition) reduces to ~0 ecliptic latitude at every
    epoch. Pairing the pole-based precession with Vondrák's separately fitted
    direct ``ε_A`` series would instead tilt the of-date ecliptic away from
    its own pole by up to ~6.5″ at −3000, surfacing as a spurious ecliptic
    latitude on the Sun; that series is therefore not used anywhere here.

    Args:
        jd_tt: Julian Date in TT.

    Returns:
        The of-date mean obliquity in radians.
    """
    epj = _julian_epoch(jd_tt)
    ex, ey, ez = (float(c) for c in erfa.ltpecl(epj))
    qx, qy, qz = (float(c) for c in erfa.ltpequ(epj))
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
def _ltp_rows(
    jd_tt: float,
) -> tuple[float, float, float, float, float, float, float, float, float]:
    """Precession matrix J2000.0 mean frame → mean equator/equinox of date.

    ``erfa.ltp`` (Vondrák 2011, no ICRS frame bias), returned as a 9-tuple of
    native floats in row-major order (r00, r01, r02, r10, r11, r12, r20, r21,
    r22). Multiplying a J2000-equator vector ``v`` by the ROWS gives the
    of-date vector (``v_date[i] = sum_j R[i][j] v[j]``).

    Args:
        jd_tt: Julian Date in TT.
    """
    r = erfa.ltp(_julian_epoch(jd_tt))
    return (
        float(r[0][0]),
        float(r[0][1]),
        float(r[0][2]),
        float(r[1][0]),
        float(r[1][1]),
        float(r[1][2]),
        float(r[2][0]),
        float(r[2][1]),
        float(r[2][2]),
    )


def _precess_j2000_to_date(vec: tuple, jd_tt: float) -> tuple[float, float, float]:
    """Rotate a J2000-equator vector to the mean equator of date."""
    r = _ltp_rows(jd_tt)
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
    """Return canonical TT-minus-UT in days for the supplied UT Julian date."""
    # Imported lazily to avoid a module import cycle (time_utils imports constants
    # but not this module; keeping the import local is defensive and cheap).
    from .time_utils import deltat

    return deltat(jd_ut1)


def _gmst_iau2006_deg(jd_ut1: float, jd_tt: float) -> float:
    """Greenwich Mean Sidereal Time (deg) via the IAU 2006 expression.

    GMST = ERA(UT1) + p(T), with the precession-in-RA polynomial p(T) of
    Capitaine, Wallace & Chapront (2003), A&A 412, 567. The Earth Rotation
    Angle ERA(UT1) uses the constants of Capitaine, Guinot & McCarthy (2000),
    A&A 355, 398 (IERS Conventions 2010, eq. 5.15). Used in the modern window
    where it is most precise.
    """
    du = jd_ut1 - _J2000
    # ERA(UT1): Capitaine, Guinot & McCarthy (2000), A&A 355, 398;
    # IERS Conventions 2010, eq. 5.15.
    era = _D2PI * (0.7790572732640 + 1.00273781191135448 * du)
    t = (jd_tt - _J2000) / 36525.0
    p = 0.014506 + t * (
        4612.156534
        + t * (1.3915817 + t * (-0.00000044 + t * (-0.000029956 + t * -0.0000000368)))
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
