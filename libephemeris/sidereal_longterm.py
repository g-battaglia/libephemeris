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

How sidereal time is computed
----------------------------
Greenwich mean sidereal time is the hour angle, at the Greenwich meridian, of
the **mean equinox of date** — the ascending intersection of the mean equator
of date with the mean ecliptic of date. Three things are needed to place it:
where the equinox of date lies, where the Sun's mean place lies, and how far
the Earth has turned. What holds here is:

* **In the era the international expression was fitted for**, the answer *is*
  the IAU 2006 Greenwich mean sidereal time — the Earth Rotation Angle of
  Capitaine, Guinot & McCarthy (2000), a strictly linear function of UT1, plus
  the precession-in-right-ascension polynomial of Capitaine, Wallace &
  Chapront (2003). That polynomial is a truncated series fitted near J2000.0:
  excellent for a few centuries, divergent outside its fit interval, by degrees
  at the ends of the DE441 span.
* **Everywhere**, the answer follows the definition, evaluated from quantities
  that are themselves valid over the whole span: the equinox of date comes from
  the Vondrák 2011 long-term precession as ERFA realizes it, the Sun's mean
  place from the secular mean longitude of the Earth of Simon et al. (1994)
  retarded by the Sun–Earth light time, and the Earth's rotation from UT1. No
  polynomial in right ascension is summed, so nothing diverges.
* **The two descriptions are one curve.** They are tied together at the ends of
  the interval over which the international expression is adopted, and away
  from that interval the answer is carried by the definition alone. Inside the
  interval the answer is the IAU 2006 value itself; across the whole span the
  value advances with the sidereal rotation of every interval, without a step,
  a kink, or a reversal.
* **The scales are the standard ones.** Sidereal time is an angle of the
  rotating Earth, so its argument is UT1; precession and obliquity are
  dynamical, so they are functions of TT. The two are joined by the library's
  own ΔT, read live, so that one chart's angles and bodies share one ΔT. At
  remote epochs the answer is correspondingly sensitive to ΔT: that sensitivity
  is physical.
* **No ephemeris is consulted**, so there is no supported range and no coverage
  refusal: every finite Julian Day is answered with a number in ``[0, 360)``.
  An infinite epoch is refused, and NaN propagates.

Apparent sidereal time is the hour angle of the *true* equinox of date: the
mean value plus the equation of the equinoxes, the classical projection of the
nutation in longitude on the equator, with the nutation supplied by the caller
so that one chart's angles and bodies share one nutation model.

Provenance:
    The Vondrák 2011 poles and precession matrix are evaluated by ERFA
    (``erfa.ltpecl``, ``erfa.ltpequ``, ``erfa.ltp``), the IAU SOFA-derived
    library; no coefficient of that model is kept in this module. The Earth
    Rotation Angle constants are those of Capitaine, Guinot & McCarthy (2000)
    as adopted in the IERS Conventions (2010), eq. 5.15; the precession in
    right ascension is the series of Capitaine, Wallace & Chapront (2003),
    A&A 412, 567; the mean longitude of the Earth is transcribed from Simon
    et al. (1994), A&A 282, 663; the light time for unit distance is the IAU
    (2009) value. The construction from those ingredients is standard published
    astronomy. The interval over which the international expression is adopted,
    the refusal of an infinite epoch, caching, and floating-point evaluation
    order are project choices; their derivation and tests are recorded in
    ``docs/methodology/sidereal-time-longterm.md``.
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


# --------------------------------------------------------------------------
# the international expression, and the era it was fitted for
# --------------------------------------------------------------------------
#: Earth Rotation Angle at J2000.0 and its rate, in turns and turns per day of
#: UT1. Capitaine, Guinot & McCarthy (2000), A&A 355, 398, as adopted in the
#: IERS Conventions (2010), eq. 5.15.
_ERA_TURNS_AT_J2000 = 0.7790572732640
_ERA_TURNS_PER_UT1_DAY = 1.00273781191135448

#: Precession accumulated in right ascension since J2000.0, in arcseconds, as
#: coefficients of a polynomial in Julian centuries of TT. Capitaine, Wallace &
#: Chapront (2003), "Expressions for IAU 2000 precession quantities",
#: A&A 412, 567 — the series the IAU 2006 sidereal time is built on.
_PRECESSION_IN_RA_ARCSEC = (
    0.014506,
    4612.156534,
    1.3915817,
    -0.00000044,
    -0.000029956,
    -0.0000000368,
)

#: Days in a Julian century and in a Julian millennium.
_JULIAN_CENTURY_DAYS = 36525.0
_JULIAN_MILLENNIUM_DAYS = 365250.0

#: Ends of the interval over which the IAU 2006 expression is adopted as the
#: realization of the definition: 1850 Jan 1.0 and 2050 Jan 1.0 as Julian Days.
#: The polynomial of ``_PRECESSION_IN_RA_ARCSEC`` is a truncated fit around
#: J2000.0; two centuries centred on the fit is where it is the best available
#: realization, and the definition of :func:`_mean_sidereal_longterm_deg`
#: carries the answer away from it.
_ADOPTION_FIRST_JD = 2396758.5
_ADOPTION_LAST_JD = 2469807.5


def _earth_rotation_angle_deg(jd_ut1: float) -> float:
    """Earth Rotation Angle at ``jd_ut1``, in degrees, not reduced to a turn.

    The rotation angle of the Earth about the Celestial Intermediate Pole is a
    strictly linear function of UT1; the accumulated angle is returned as it
    stands, because the caller reduces the total.

    Args:
        jd_ut1: Julian Date on the UT1 scale.

    Returns:
        The accumulated Earth Rotation Angle in degrees.
    """
    turns = _ERA_TURNS_AT_J2000 + _ERA_TURNS_PER_UT1_DAY * (jd_ut1 - _J2000)
    return math.degrees(_D2PI * turns)


def _precession_in_right_ascension_deg(jd_tt: float) -> float:
    """Precession accumulated in right ascension since J2000.0, in degrees.

    Args:
        jd_tt: Julian Date in TT (precession is a dynamical quantity).

    Returns:
        The accumulated precession in right ascension, in degrees.
    """
    t = (jd_tt - _J2000) / _JULIAN_CENTURY_DAYS
    a0, a1, a2, a3, a4, a5 = _PRECESSION_IN_RA_ARCSEC
    return (a0 + t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5))))) / 3600.0


def _gmst_iau2006_deg(jd_ut1: float, jd_tt: float) -> float:
    """Greenwich mean sidereal time from the IAU 2006 expression, in degrees.

    The rotation term is a function of UT1 and the precession term a function
    of TT, which is why both scales are arguments.

    Args:
        jd_ut1: Julian Date on the UT1 scale.
        jd_tt: the same instant on the TT scale.

    Returns:
        Greenwich mean sidereal time in degrees, reduced to ``[0, 360)``.
    """
    return (
        _earth_rotation_angle_deg(jd_ut1) + _precession_in_right_ascension_deg(jd_tt)
    ) % 360.0


# --------------------------------------------------------------------------
# the definition, valid over the whole span
# --------------------------------------------------------------------------
#: Mean longitude of the Earth referred to the mean dynamical ecliptic and
#: equinox of J2000.0: the epoch value in degrees and the secular terms in
#: arcseconds per Julian millennium of TT. Simon, J. L., Bretagnon, P.,
#: Chapront, J., Chapront-Touzé, M., Francou, G. & Laskar, J. (1994),
#: "Numerical expressions for precession formulae and mean elements for the
#: Moon and the planets", A&A 282, 663.
_EARTH_MEAN_LONGITUDE_DEG = 100.46645683
_EARTH_MEAN_LONGITUDE_ARCSEC = (1295977422.83429, -2.04411, -0.00523)

#: Light time for unit distance, tau_A = 499.004 783 8061 s, in days. IAU (2009)
#: system of astronomical constants; Luzum et al. (2011), CMDA 110, 293.
_LIGHT_TIME_FOR_UNIT_DISTANCE_DAYS = 499.0047838061 / 86400.0

#: The Sun's mean motion in longitude, at which the light time above is turned
#: into an angle: one revolution per mean year of the Gregorian calendar. Which
#: mean year is a convention, and the convention cannot matter here: reckoning
#: the same light time at the sidereal rate of the mean longitude below, or at
#: the mean tropical year, moves the sidereal time by 1e-7 to 4e-7 arcsec at
#: the ends of the span and by less than that everywhere else -- four orders of
#: magnitude below the 0.001 arcsec the library claims, and below the last
#: representable place of the accumulated longitude it is subtracted from.
_SUN_MEAN_MOTION_DEG_PER_DAY = 360.0 / 365.2425

#: The angle by which light travel time retards the Sun's mean place: about
#: 20.49 arcseconds. A geometric direction would leave the equinox that far
#: away, ten orders of magnitude outside the precision claimed here.
_LIGHT_TIME_RETARDATION_DEG = (
    _SUN_MEAN_MOTION_DEG_PER_DAY * _LIGHT_TIME_FOR_UNIT_DISTANCE_DAYS
)


def _sun_apparent_mean_longitude_deg(jd_tt: float) -> float:
    """Apparent mean longitude of the Sun in the J2000.0 mean ecliptic.

    The Sun's mean place is the Earth's mean longitude turned by half a circle
    and retarded by the Sun–Earth light time; the secular expression is valid
    over the whole span the ephemeris covers.

    Args:
        jd_tt: Julian Date in TT.

    Returns:
        The longitude in degrees, reduced to ``[0, 360)``, referred to the mean
        dynamical ecliptic and equinox of J2000.0.
    """
    t = (jd_tt - _J2000) / _JULIAN_MILLENNIUM_DAYS
    c1, c2, c3 = _EARTH_MEAN_LONGITUDE_ARCSEC
    earth = _EARTH_MEAN_LONGITUDE_DEG + (c1 * t + c2 * t * t + c3 * t**3) / 3600.0
    return (earth - _LIGHT_TIME_RETARDATION_DEG + 180.0) % 360.0


def _ecliptic_longitude_of_date_deg(longitude_j2000_deg: float, jd_tt: float) -> float:
    """Carry a J2000.0 ecliptic longitude to the mean ecliptic of date.

    The direction is taken to the J2000.0 equator with the obliquity of that
    epoch, precessed with the Vondrák long-term matrix, and read back on the
    ecliptic of date with the obliquity of date. Both obliquities are the angle
    between the same pair of poles the matrix is built from, so the equinox,
    the equator and the ecliptic of date are one frame.

    Args:
        longitude_j2000_deg: ecliptic longitude in the J2000.0 mean frame.
        jd_tt: Julian Date in TT of the frame to carry it to.

    Returns:
        The ecliptic longitude of date in degrees, in ``(-180, 180]``.
    """
    longitude = math.radians(longitude_j2000_deg)
    on_j2000_ecliptic = (math.cos(longitude), math.sin(longitude), 0.0)
    on_j2000_equator = _rot_x(on_j2000_ecliptic, -mean_obliquity_rad(_J2000))
    on_equator_of_date = _precess_j2000_to_date(on_j2000_equator, jd_tt)
    on_ecliptic_of_date = _rot_x(on_equator_of_date, mean_obliquity_rad(jd_tt))
    return math.degrees(math.atan2(on_ecliptic_of_date[1], on_ecliptic_of_date[0]))


def _mean_sidereal_longterm_deg(jd_ut1: float, jd_tt: float) -> float:
    """Greenwich hour angle of the mean equinox of date, from the definition.

    Sidereal time is the hour angle of the equinox, which is the hour angle of
    the mean Sun plus the mean Sun's right ascension; the fictitious mean Sun
    runs on the equator with the right ascension the true Sun has as a mean
    longitude, so the second term is the Sun's mean place carried to the
    ecliptic of date. Universal Time gives the first term: it is mean solar
    time at Greenwich, so the mean Sun's hour angle is the UT1 time of day less
    twelve hours. Nothing here is a polynomial fitted near J2000.0, so the
    construction holds over the whole span of the ephemeris.

    Args:
        jd_ut1: Julian Date on the UT1 scale.
        jd_tt: the same instant on the TT scale. It is an argument rather than
            a lookup so that the caller reads ΔT once, and nothing here is
            memoised on a UT1 epoch: a session that changes the ΔT model, or
            supplies its own value, changes the sidereal time with it.

    Returns:
        The hour angle in degrees, reduced to ``[0, 360)``.

    Raises:
        ValueError: if ``jd_ut1`` is infinite; an infinite epoch has no place
            on the time axis of a rotating Earth. NaN propagates instead.
    """
    hour_angle_of_mean_sun = math.fmod(jd_ut1 + 0.5, 1.0) * 360.0 - 180.0
    mean_sun = _sun_apparent_mean_longitude_deg(jd_tt)
    return (
        hour_angle_of_mean_sun + _ecliptic_longitude_of_date_deg(mean_sun, jd_tt)
    ) % 360.0


# --------------------------------------------------------------------------
# one curve over the whole span
# --------------------------------------------------------------------------
def _adoption_epoch(jd_ut1: float) -> float:
    """The epoch at which the international expression is read for ``jd_ut1``.

    Inside the interval where the IAU 2006 expression is adopted this is the
    epoch itself, so the answer is that expression; outside it, it is the end
    of the interval the epoch lies beyond, so the definition carries the answer
    on from the last place the expression is read. NaN passes through.

    Args:
        jd_ut1: Julian Date on the UT1 scale.

    Returns:
        The epoch to read the international expression at.
    """
    if jd_ut1 < _ADOPTION_FIRST_JD:
        return _ADOPTION_FIRST_JD
    if jd_ut1 > _ADOPTION_LAST_JD:
        return _ADOPTION_LAST_JD
    return jd_ut1


def _to_signed_turn_deg(delta_deg: float) -> float:
    """Reduce an angular difference in degrees to ``(-180, 180]``."""
    reduced = delta_deg % 360.0
    return reduced - 360.0 if reduced > 180.0 else reduced


def mean_sidereal_time_deg(jd_ut1: float) -> float:
    """Greenwich mean sidereal time in degrees for a Julian Day in UT1.

    The answer is the IAU 2006 expression read at the adoption epoch, advanced
    by the sidereal rotation the definition puts between that epoch and the
    date. In the era the expression was fitted for the two epochs coincide and
    the answer is the expression itself; away from it the advance comes
    entirely from the long-term construction, which is what keeps the curve
    continuous and free of the polynomial's divergence.

    Args:
        jd_ut1: Julian Date on the UT1 scale. Any finite value is answered,
            inside or outside the range of any ephemeris.

    Returns:
        Greenwich mean sidereal time in degrees, reduced to ``[0, 360)``.

    Raises:
        ValueError: if ``jd_ut1`` is infinite. NaN propagates and yields NaN.
    """
    epoch = _adoption_epoch(jd_ut1)
    epoch_tt = epoch + _deltat_days(epoch)
    at_epoch = _mean_sidereal_longterm_deg(epoch, epoch_tt)
    at_date = _mean_sidereal_longterm_deg(jd_ut1, jd_ut1 + _deltat_days(jd_ut1))
    advance = _to_signed_turn_deg(at_date - at_epoch)
    return (_gmst_iau2006_deg(epoch, epoch_tt) + advance) % 360.0


def apparent_sidereal_time_deg(
    jd_ut1: float, longitude: float = 0.0, *, dpsi_deg: float, eps_true_deg: float
) -> float:
    """Local apparent sidereal time in degrees, the ARMC of the house engine.

    Apparent sidereal time is the hour angle of the *true* equinox of date: the
    mean value plus the equation of the equinoxes, the nutation in longitude
    projected on the equator. The nutation pair is supplied by the caller, so
    one chart's angles and bodies share one nutation model, and the
    complementary terms of the IAU 2000 equation of the equinoxes are not
    applied here. East longitude enters as a plain additive angle, so the
    answer at longitude zero is the Greenwich value.

    Args:
        jd_ut1: Julian Date on the UT1 scale.
        longitude: geographic longitude in degrees, East positive.
        dpsi_deg: nutation in longitude at the same instant, in degrees.
        eps_true_deg: true obliquity of the ecliptic of date, in degrees.

    Returns:
        Local apparent sidereal time in degrees, reduced to ``[0, 360)``.

    Raises:
        ValueError: if ``jd_ut1`` is infinite. NaN propagates and yields NaN.
    """
    equation_of_the_equinoxes = dpsi_deg * math.cos(math.radians(eps_true_deg))
    return (
        mean_sidereal_time_deg(jd_ut1) + equation_of_the_equinoxes + longitude
    ) % 360.0
