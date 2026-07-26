# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Time conversion utilities for libephemeris.

Implements standard astronomical time functions for conversions between:
- Calendar dates and Julian Day numbers
- Gregorian and Julian calendar systems
- UT1 (Universal Time) and TT (Terrestrial Time)

Functions provide reference-API-compatible signatures.
All algorithms follow Meeus "Astronomical Algorithms" (1998).

Provenance:
    Calendar/Julian-day transforms follow the cited published algorithms;
    TT/UT1, UTC, leap-second, and fundamental-argument relations follow IERS
    Conventions (2010) and observed IERS data. Historical Delta-T branches are
    documented in ``docs/methodology/delta-t.md``. Blend dates, cache use,
    compatibility signatures, and user-override precedence are project choices;
    no coefficient is recovered from compatibility output.
"""

from __future__ import annotations

from math import floor as _floor
from typing import Any, Optional

import math

from .constants import (
    GREG_CAL,
    JUL_CAL,
    FLG_EQUATORIAL,
    FLG_MOSEPH,
    FLG_SWIEPH,
    TIDAL_DEFAULT,
    TIDAL_MOSEPH,
)
from .cache import get_cached_time_ut1
from .state import (
    _tid_acc_is_set,
    get_delta_t_model,
    get_delta_t_userdef,
    get_iers_delta_t_enabled,
    get_tid_acc,
    get_timescale,
)

# Julian Day of Gregorian calendar reform: Oct 15, 1582
JD_GREGORIAN_REFORM = 2299161

# --- Tidal-acceleration adjustment of Delta T -------------------------------
# Historical lunar reductions are conventionally referred to 1955.5; see
# Morrison & Stephenson (2004), JHA 35, 327-336, and the LUNAR97 series.
# Expressing that Julian epoch from J2000 avoids an opaque recovered JD.
_TID_ACC_REFERENCE_JD = 2451545.0 - 44.5 * 365.25

# IERS Conventions (2010), TN 36, Eq. 5.43: the linear rates of
# F = L - Omega and Omega in arcsec per Julian century.  Their sum is the
# mean lunar longitude rate.  A change dacc in lunar longitude acceleration
# accumulates 0.5*dacc*T^2 arcsec; dividing by mean lunar motion converts the
# phase displacement to seconds of time.  Thus the familiar ~0.91072 factor
# is derived from IERS Conventions (2010) rather than an empirical fit.
_IERS_F_RATE_ARCSEC_CY = 1739527262.8478
_IERS_OMEGA_RATE_ARCSEC_CY = -6962890.5431

# --- UTC epoch boundary (jd*_to_utc classification) ------------------------
# JD of 1972-01-01 00:00:00, the start of leap-second UTC. Before this instant
# there is no UTC and the reference API returns the UT1 calendar date directly.
# Instants within a guard band of this epoch are ambiguous when classified by
# the UT1 year alone (UT1 leads/lags UTC by up to ~0.9 s), so they are routed
# through Skyfield's UTC to be classified on the reconstructed UTC year.
_UTC_EPOCH_JD = 2441317.5
_UTC_EPOCH_BAND = 2.0 / 86400.0  # ~2 s window treated as "at the UTC epoch"


def _validate_calendar(cal: int, func_name: str) -> None:
    """Reject calendar flags other than GREG_CAL/JUL_CAL.

    The reference binding validates the calendar flag in every calendar
    conversion function (e.g. ``julday: invalid calendar (99)``);
    all five libephemeris counterparts share this guard.

    Raises:
        ValueError: If cal is neither GREG_CAL nor JUL_CAL.
    """
    if cal not in (GREG_CAL, JUL_CAL):
        raise ValueError(f"{func_name}: invalid calendar ({cal})")


def julday(
    year: int, month: int, day: int, hour: float = 12.0, cal: int = GREG_CAL
) -> float:
    """
    Convert calendar date to Julian Day number.

    Args:
        year: Calendar year (negative for BCE)
        month: Month (1-12)
        day: Day of month (1-31)
        hour: Decimal hour (0.0-23.999...)
        cal: GREG_CAL (1) for Gregorian, JUL_CAL (0) for Julian

    Returns:
        float: Julian Day number (days since JD 0.0 = noon Jan 1, 4713 BCE)

    Note:
        Transition date: Oct 15, 1582 (Gregorian) = Oct 5, 1582 (Julian)
        JD 2451545.0 = Jan 1, 2000 12:00 TT (J2000.0 epoch)
    """
    _validate_calendar(cal, "julday")

    if month <= 2:
        year -= 1
        month += 12

    # Meeus INT() is defined as floor(), not truncation toward zero.
    # Python's int() truncates toward zero which gives wrong results
    # for negative years (e.g. int(-30/4)=-7 but floor(-30/4)=-8).
    a = _floor(year / 100)

    if cal == GREG_CAL:
        b = 2 - a + _floor(a / 4)
    else:
        b = 0

    jd = (
        _floor(365.25 * (year + 4716))
        + _floor(30.6001 * (month + 1))
        + day
        + hour / 24.0
        + b
        - 1524.5
    )
    return jd


def revjul(jd: float, cal: int = GREG_CAL) -> tuple[int, int, int, float]:
    """
    Convert Julian Day number to calendar date.

    Args:
        jd: Julian Day number
        cal: GREG_CAL (1) for Gregorian, JUL_CAL (0) for Julian

    Returns:
        tuple: (year, month, day, hour) where:
            - year: Calendar year
            - month: Month (1-12)
            - day: Integer day of month
            - hour: Decimal hour (0.0-23.999...)

    Note:
        Automatic Gregorian calendar used for JD >= 2299161 (Oct 15, 1582)
        unless Julian calendar explicitly requested.
    """
    _validate_calendar(cal, "revjul")

    jd = jd + 0.5
    z = _floor(jd)
    f = jd - z

    # Always respect cal — the reference uses proleptic Gregorian for ancient dates
    # when GREG_CAL is requested, not auto-detection by JD threshold.
    if cal == GREG_CAL:
        alpha = _floor((z - 1867216.25) / 36524.25)
        a = z + 1 + alpha - _floor(alpha / 4)
    else:
        a = z

    b = a + 1524
    c = _floor((b - 122.1) / 365.25)
    d = _floor(365.25 * c)
    e = _floor((b - d) / 30.6001)

    day = b - d - _floor(30.6001 * e) + f
    if e < 14:
        month = e - 1
    else:
        month = e - 13

    if month > 2:
        year = c - 4716
    else:
        year = c - 4715

    d_int = int(day)
    d_frac = day - d_int
    hour = d_frac * 24.0
    day = d_int

    return year, month, day, hour


def _deltat_espenak_meeus(year: float) -> float:
    """Delta T in **seconds** via the Espenak & Meeus (2006) polynomial set.

    This is the classic model used by NASA's eclipse predictions, valid for
    roughly -1999 to +3000 (with a parabolic tail outside). ``year`` is a
    decimal Gregorian year (``year + (month-0.5)/12``).

    Reference: Espenak, F. & Meeus, J. (2006), "Five Millennium Canon of Solar
    Eclipses: -1999 to +3000"; NASA GSFC polynomial expressions for Delta T.
    """
    y = year
    if y < -500:
        u = (y - 1820) / 100.0
        return -20 + 32 * u * u
    if y < 500:
        u = y / 100.0
        return (
            10583.6
            - 1014.41 * u
            + 33.78311 * u**2
            - 5.952053 * u**3
            - 0.1798452 * u**4
            + 0.022174192 * u**5
            + 0.0090316521 * u**6
        )
    if y < 1600:
        u = (y - 1000) / 100.0
        return (
            1574.2
            - 556.01 * u
            + 71.23472 * u**2
            + 0.319781 * u**3
            - 0.8503463 * u**4
            - 0.005050998 * u**5
            + 0.0083572073 * u**6
        )
    if y < 1700:
        t = y - 1600
        return 120 - 0.9808 * t - 0.01532 * t**2 + t**3 / 7129.0
    if y < 1800:
        t = y - 1700
        return (
            8.83 + 0.1603 * t - 0.0059285 * t**2 + 0.00013336 * t**3 - t**4 / 1174000.0
        )
    if y < 1860:
        t = y - 1800
        return (
            13.72
            - 0.332447 * t
            + 0.0068612 * t**2
            + 0.0041116 * t**3
            - 0.00037436 * t**4
            + 0.0000121272 * t**5
            - 0.0000001699 * t**6
            + 0.000000000875 * t**7
        )
    if y < 1900:
        t = y - 1860
        return (
            7.62
            + 0.5737 * t
            - 0.251754 * t**2
            + 0.01680668 * t**3
            - 0.0004473624 * t**4
            + t**5 / 233174.0
        )
    if y < 1920:
        t = y - 1900
        return (
            -2.79 + 1.494119 * t - 0.0598939 * t**2 + 0.0061966 * t**3 - 0.000197 * t**4
        )
    if y < 1941:
        t = y - 1920
        return 21.20 + 0.84493 * t - 0.076100 * t**2 + 0.0020936 * t**3
    if y < 1961:
        t = y - 1950
        return 29.07 + 0.407 * t - t**2 / 233.0 + t**3 / 2547.0
    if y < 1986:
        t = y - 1975
        return 45.45 + 1.067 * t - t**2 / 260.0 - t**3 / 718.0
    if y < 2005:
        t = y - 2000
        return (
            63.86
            + 0.3345 * t
            - 0.060374 * t**2
            + 0.0017275 * t**3
            + 0.000651814 * t**4
            + 0.00002373599 * t**5
        )
    if y < 2050:
        t = y - 2000
        return 62.92 + 0.32217 * t + 0.005589 * t**2
    if y < 2150:
        return -20 + 32 * ((y - 1820) / 100.0) ** 2 - 0.5628 * (2150 - y)
    u = (y - 1820) / 100.0
    return -20 + 32 * u * u


def _tid_acc_adjustment_seconds(tjdut: float, tid_acc: float) -> float:
    """Delta T adjustment (seconds) for a non-default tidal acceleration.

    Uses the standard conversion of a lunar secular-acceleration difference
    into clock seconds, referred to Julian epoch 1955.5.  The conversion is
    derived from the IERS 2010 mean-lunar-longitude rate; it is linear in
    ``tid_acc`` and quadratic in elapsed Julian centuries.  Modern observed
    Delta T is left unchanged from the reference epoch onward.

    Args:
        tjdut: Julian Day number in UT1.
        tid_acc: Tidal acceleration in arcsec/century^2.

    Returns:
        Adjustment to add to Delta T, in seconds.
    """
    dacc = tid_acc - TIDAL_DEFAULT
    if dacc == 0.0 or tjdut >= _TID_ACC_REFERENCE_JD:
        return 0.0
    u = (tjdut - _TID_ACC_REFERENCE_JD) / 36525.0
    lunar_motion_arcsec_day = (
        _IERS_F_RATE_ARCSEC_CY + _IERS_OMEGA_RATE_ARCSEC_CY
    ) / 36525.0
    seconds_per_arcsec = 86400.0 / lunar_motion_arcsec_day
    return -0.5 * seconds_per_arcsec * dacc * u * u


def _deltat_with_tid_acc(tjdut: float, tid_acc: float) -> float:
    """Delta T in days for an explicit tidal acceleration.

    Shared core of :func:`deltat` and :func:`deltat_ex`: resolves the
    userdef / IERS / model priority chain documented in :func:`deltat`,
    then applies the tidal-acceleration adjustment for ``tid_acc`` (zero
    at TIDAL_DEFAULT, and zero for dates from the 1955.5 reference epoch on).
    A user-defined Delta T is a hard override and is never adjusted.
    """
    # Check for user-defined Delta T first
    userdef_dt = get_delta_t_userdef()
    if userdef_dt is not None:
        return userdef_dt

    adjust_seconds = _tid_acc_adjustment_seconds(tjdut, tid_acc)

    # Check for IERS Delta T if enabled
    if get_iers_delta_t_enabled():
        try:
            from . import iers_data

            delta_t_seconds = iers_data.get_delta_t_iers(tjdut)
            if delta_t_seconds is not None:
                # IERS returns seconds, convert to days
                return (delta_t_seconds + adjust_seconds) / 86400.0
        except Exception as exc:  # noqa: BLE001  # robust fallback, any error
            # Fall back to Skyfield if IERS data fails for ANY reason — Delta T
            # feeds every downstream position calculation, so a robust fallback
            # must not let an unexpected exception type (OSError, KeyError, ...)
            # escape and crash the whole pipeline. Log at debug so the silent
            # fallback is still observable when diagnosing Delta T discrepancies.
            from .logging_config import get_logger

            get_logger().debug(
                "IERS Delta T lookup failed (%s); falling back to Skyfield.", exc
            )

    # Selected Delta T model (after userdef / IERS). Default is SMH-2016,
    # obtained from Skyfield (an MIT dependency); ``espenak_meeus`` is a native
    # reimplementation of the classic NASA polynomials. Neither derives from
    # the reference (copyleft) ephemeris. Exact parity with the reference,
    # when needed for validation, is injected externally via
    # set_delta_t_userdef().
    if get_delta_t_model() == "espenak_meeus":
        year, month, _day, _hour = revjul(tjdut)
        decimal_year = year + (month - 0.5) / 12.0
        return (_deltat_espenak_meeus(decimal_year) + adjust_seconds) / 86400.0

    t = get_cached_time_ut1(tjdut)
    delta_t_seconds = float(t.delta_t)
    return (delta_t_seconds + adjust_seconds) / 86400.0


def deltat(tjdut: float) -> float:
    """Calculate Delta T (TT minus UT1) for a Julian Day.

    Args:
        tjdut: Julian Day in UT1.

    Returns:
        Delta T in days.

    Notes:
        A user-defined value has first priority. Observed IERS data is next
        when enabled and available. Otherwise the configured model is used:
        ``smh2016`` (the default) or ``espenak_meeus``. Before 1955, an
        explicitly selected tidal acceleration can adjust the model result.

        Outside the SMH-2016 spline interval, Skyfield's implementation uses
        the long-term parabola ``Delta T = -320 + 32.5 * u**2`` seconds, where
        ``u = (year - 1825) / 100``.

        The selected Delta-T value drives calculations that convert UT to TT
        through this function, including LEB, Horizons, event searches, and
        long-term sidereal time. The direct Skyfield position path uses the
        Delta-T realization in its Skyfield time object.

    References:
        Stephenson, Morrison, and Hohenkerk (2016), *Proceedings of the Royal
        Society A* 472, 20160404.
    """
    # Pure functional call: the effective tidal acceleration is resolved for
    # this computation only. get_tid_acc() changes exclusively through
    # set_tid_acc(), never as a side effect of a Delta-T query.
    return _deltat_with_tid_acc(tjdut, get_tid_acc())


def deltat_ex(tjdut: float, flag: int = FLG_SWIEPH) -> float:
    """Calculate Delta T with an explicit ephemeris-selection flag.

    Args:
        tjdut: Julian Day in UT1.
        flag: Ephemeris-selection flags. ``FLG_SWIEPH`` and
            ``FLG_JPLEPH`` use the standard tidal acceleration;
            ``FLG_MOSEPH`` implies ``TIDAL_MOSEPH`` unless an explicit
            ``set_tid_acc()`` setting takes precedence.

    Returns:
        Delta T in days.

    Notes:
        The user-defined value, optional IERS observations, and configured
        model have the same priority as in ``deltat()``. A user-defined
        Delta T is never tidally adjusted. The default configured model is the
        Stephenson, Morrison, and Hohenkerk (2016) realization supplied by
        Skyfield, including its long-term parabola outside the spline range.

        This function is pure with respect to library state: the
        flag-selected tidal acceleration is applied to this computation
        only, and ``get_tid_acc()`` changes exclusively through
        ``set_tid_acc()``.
    """
    # A hard Delta-T override wins before any ephemeris selection.
    if get_delta_t_userdef() is not None:
        return _deltat_with_tid_acc(tjdut, get_tid_acc())

    # FLG_MOSEPH selects the semi-analytical ephemeris value whenever its bit
    # is present. Explicit set_tid_acc() values remain authoritative. The
    # selection is functional: it never leaks into get_tid_acc().
    if _tid_acc_is_set():
        acc = get_tid_acc()
    else:
        acc = TIDAL_MOSEPH if flag & FLG_MOSEPH else TIDAL_DEFAULT
    return _deltat_with_tid_acc(tjdut, acc)


def date_conversion(
    year: int, month: int, day: int, hour: float, calendar: str
) -> tuple[int, int, int, float]:
    """
    Convert a date between Julian and Gregorian calendars.

    The function automatically detects the input calendar based on the date:
    - Dates before Oct 15, 1582 are assumed to be Julian calendar
    - Dates from Oct 15, 1582 onwards are assumed to be Gregorian calendar

    Args:
        year: Calendar year
        month: Month (1-12)
        day: Day of month (1-31)
        hour: Decimal hour (0.0-23.999...)
        calendar: Target calendar - 'j' for Julian or 'g' for Gregorian

    Returns:
        tuple: (year, month, day, hour) in the requested calendar

    Raises:
        ValueError: If calendar is not 'j' or 'g'

    Note:
        The Gregorian calendar reform occurred on Oct 15, 1582. On this date,
        the Julian calendar was 10 days behind. The function uses Julian Day
        numbers as an intermediate representation to convert between calendars.

    Example:
        >>> # Convert first Gregorian date to Julian
        >>> date_conversion(1582, 10, 15, 12.0, 'j')
        (1582, 10, 5, 12.0)
        >>> # Convert Julian date to Gregorian
        >>> date_conversion(1582, 10, 5, 12.0, 'g')
        (1582, 10, 15, 12.0)
    """
    # Black-box compatibility: the public API accepts bytes calendar codes.
    if isinstance(calendar, bytes):
        calendar = calendar.decode("ascii")
    calendar = calendar.lower()
    if calendar not in ("j", "g"):
        raise ValueError(f"calendar must be 'j' or 'g', got: {calendar!r}")

    # Determine the input calendar from the calendar DATE, not an
    # hour-inclusive JD. JD_GREGORIAN_REFORM is NOON of 1582-10-15, so an
    # instant in the 00:00-11:59 window of the reform day would compare below
    # it and be misclassified as Julian input (a spurious ~10-day shift).
    jd_as_greg = julday(year, month, day, hour, GREG_CAL)
    if (year, month, day) < (1582, 10, 15):
        input_cal = JUL_CAL
        jd = julday(year, month, day, hour, JUL_CAL)
    else:
        input_cal = GREG_CAL
        jd = jd_as_greg

    # Determine target calendar flag
    target_cal = JUL_CAL if calendar == "j" else GREG_CAL

    # If input and target are the same, just return the original values
    if input_cal == target_cal:
        return year, month, day, hour

    # Convert via Julian Day to target calendar
    return revjul(jd, target_cal)


# Dates on which a positive leap second was inserted at 23:59:60 UTC.
# Leap seconds are always inserted at the end of June 30 or December 31.
# Source: IERS Bulletin C / BIPM.  Last entry: 2016-12-31 (TAI-UTC = 37s).
_LEAP_SECOND_DATES: frozenset[tuple[int, int, int]] = frozenset(
    [
        (1972, 6, 30),
        (1972, 12, 31),
        (1973, 12, 31),
        (1974, 12, 31),
        (1975, 12, 31),
        (1976, 12, 31),
        (1977, 12, 31),
        (1978, 12, 31),
        (1979, 12, 31),
        (1981, 6, 30),
        (1982, 6, 30),
        (1983, 6, 30),
        (1985, 6, 30),
        (1987, 12, 31),
        (1989, 12, 31),
        (1990, 12, 31),
        (1992, 6, 30),
        (1993, 6, 30),
        (1994, 6, 30),
        (1995, 12, 31),
        (1997, 6, 30),
        (1998, 12, 31),
        (2005, 12, 31),
        (2008, 12, 31),
        (2012, 6, 30),
        (2015, 6, 30),
        (2016, 12, 31),
    ]
)


def _is_leap_second_date(year: int, month: int, day: int) -> bool:
    """Return True if a positive leap second was inserted at the end of this date."""
    return (year, month, day) in _LEAP_SECOND_DATES


def _days_in_month(year: int, month: int, cal: int) -> int:
    """Return the number of days in a calendar month.

    The leap-year rule is calendar-aware, using astronomical year numbering
    (year 0 exists, and e.g. -100 is divisible by 4): the Julian calendar
    intercalates every 4th year, the proleptic Gregorian calendar skips
    century years not divisible by 400.

    Args:
        year: Astronomical year number.
        month: Month (1-12).
        cal: GREG_CAL or JUL_CAL.

    Returns:
        Number of days (28-31).
    """
    if month in (4, 6, 9, 11):
        return 30
    if month != 2:
        return 31
    if cal == GREG_CAL:
        leap = year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)
    else:
        leap = year % 4 == 0
    return 29 if leap else 28


def utc_to_jd(
    year: int,
    month: int,
    day: int,
    hour: int,
    minute: int,
    second: float,
    calendar: int = GREG_CAL,
) -> tuple[float, float]:
    """
    Convert UTC date/time to Julian Day numbers, properly handling leap seconds.

    Unlike julday() which assumes UT1 input, this function takes UTC input and
    correctly accounts for the difference between UTC and UT1 (DUT1) and between
    UTC and TT (including leap seconds).

    Args:
        year: Calendar year (negative for BCE)
        month: Month (1-12)
        day: Day of month (1-31)
        hour: Hour (0-23)
        minute: Minute (0-59)
        second: Second (0-60, allowing for leap seconds)
        calendar: GREG_CAL (1) for Gregorian, JUL_CAL (0) for Julian

    Returns:
        tuple: (jd_et, jd_ut) where:
            - jd_et: Julian Day in TT (Terrestrial Time / Ephemeris Time)
            - jd_ut: Julian Day in UT1 (Universal Time)

    Raises:
        Error: If the date does not exist in the requested calendar
            ("invalid date: ...") or the clock time is out of range
            ("invalid time: ..."); seconds 60.0-60.999... are accepted
            only at 23:59 on an IERS leap-second date.
        ValueError: If calendar is neither GREG_CAL nor JUL_CAL.

    Note:
        - UTC includes leap seconds while UT1 follows Earth's rotation
        - ``abs(UTC - UT1)`` is always less than 0.9 seconds by definition
        - TT = TAI + 32.184 seconds, where TAI is atomic time
        - For dates before 1972 (when UTC was standardized), the function
          treats the input as UT1 approximation and still provides proper
          TT/UT1 conversion using historical Delta T values

    Example:
        >>> from libephemeris import utc_to_jd, GREG_CAL
        >>> # J2000.0 epoch: Jan 1, 2000 12:00:00 UTC
        >>> jd_tt, jd_ut = utc_to_jd(2000, 1, 1, 12, 0, 0.0, GREG_CAL)
        >>> print(f"JD(TT): {jd_tt:.6f}, JD(UT1): {jd_ut:.6f}")
        JD(TT): 2451545.000743, JD(UT1): 2451545.000004
    """
    # Reference-API parity: the binding validates the calendar flag here
    # too, like julday()/revjul(). Without this guard cal=99 would mean
    # Julian in julday() but Gregorian here (and even flip interpretation
    # across the 1972 boundary, where the pre-UTC branch defers to julday).
    _validate_calendar(calendar, "utc_to_jd")

    # Reference-API parity: validate the calendar date before the clock time
    # (a bad date wins even when the time is also invalid), against the raw
    # input fields in the requested calendar. Rejection, not normalization:
    # month 13 or Feb 30 must not roll over into the next month.
    if (
        month < 1
        or month > 12
        or day < 1
        or day > _days_in_month(year, month, calendar)
    ):
        from .exceptions import Error

        raise Error(f"invalid date: year = {year}, month = {month}, day = {day}")

    # Reference-API parity: reject out-of-range clock fields (hour 24 included).
    # Seconds in [60, 61) are provisionally admitted here; the leap-second
    # check below decides whether such a 23:59:60.x label really exists.
    if hour < 0 or hour > 23 or minute < 0 or minute > 59 or second < 0.0:
        from .exceptions import Error

        raise Error(f"invalid time: {hour}:{minute}:{second:.2f}")

    # Leap-second dates and the beginning of the UTC era are defined in the
    # Gregorian calendar.  Map a Julian-calendar label to that physical
    # Gregorian date once, before either decision.  In particular, Julian
    # 1971-12-19 is Gregorian 1972-01-01 and must use the UTC-era path.
    if calendar == JUL_CAL:
        greg_jd = julday(year, month, day, 12.0, JUL_CAL)
        greg_year, greg_month, greg_day, _ = revjul(greg_jd, GREG_CAL)
    else:
        greg_year, greg_month, greg_day = year, month, day

    # Validate leap second: second=60 is only valid at 23:59:60 on dates
    # where a leap second was actually inserted (end of June 30 or Dec 31).
    # Raises the project Error (reference API raises its Error class too).
    if second >= 60.0:
        from .exceptions import Error

        # The leap-second table holds Gregorian dates, so a Julian calendar
        # input must be mapped to its Gregorian date before the lookup: e.g.
        # Julian 2016-12-18 is the real Gregorian 2016-12-31 leap second,
        # while Julian 2016-12-31 (Gregorian 2017-01-13) is not. The calendar
        # shift only relabels the date, never the clock time, so converting at
        # noon yields the right date for any time within the day.
        # Keep structurally invalid clock fields distinct from a valid 23:59:60
        # label on a date without an IERS leap-second insertion. Pre-1972 dates
        # have no leap-second UTC and therefore fall in the former category.
        if second >= 61.0 or hour != 23 or minute != 59 or greg_year < 1972:
            raise Error(f"invalid time: {hour}:{minute}:{second:.2f}")
        if not _is_leap_second_date(greg_year, greg_month, greg_day):
            raise Error(f"invalid time (no leap second!): {hour}:{minute}:{second:.2f}")

    # Before 1972 UTC (with leap seconds) did not exist; the reference API
    # treats the input as UT1 directly: jd_ut1 is the literal calendar JD
    # and jd_et = jd_ut1 + Delta T. (Routing through Skyfield's ts.utc()
    # would clamp TAI-UTC to 10 s, drifting up to tens of minutes in
    # antiquity — verified +24 s at 1800, -27 min at year 1000.)
    if greg_year < 1972:
        decimal_hour = hour + minute / 60.0 + second / 3600.0
        jd_ut1 = julday(year, month, day, decimal_hour, calendar)
        jd_et = jd_ut1 + deltat(jd_ut1)
        return float(jd_et), float(jd_ut1)

    ts = get_timescale()

    if calendar == JUL_CAL:
        # Convert only the Julian calendar date to Gregorian for Skyfield
        # (which expects the proleptic Gregorian calendar). The calendar shift
        # relabels the date, not the clock, so hour/minute/second — including a
        # :60 leap second — pass straight through. Folding them into a decimal
        # hour instead would push 23:59:60 into the next midnight and drop the
        # leap second entirely.
        t = ts.utc(greg_year, greg_month, greg_day, hour, minute, second)
    else:
        # Gregorian calendar - use directly with Skyfield
        t = ts.utc(year, month, day, hour, minute, second)

    # Explicit float cast to satisfy type checker (Skyfield uses lazy reify decorator)
    return float(t.tt), float(t.ut1)


def _jd_to_calendar_tuple(
    jd: float, calendar: int
) -> tuple[int, int, int, int, int, float]:
    """Split a Julian Day into (y, m, d, hh, mm, ss.s) calendar components."""
    y, m, d, decimal_hour = revjul(jd, calendar)
    hh = int(decimal_hour)
    minute_frac = (decimal_hour - hh) * 60.0
    mm = int(minute_frac)
    ss = (minute_frac - mm) * 60.0
    return y, m, d, hh, mm, ss


def jdet_to_utc(
    jd_et: float, calendar: int = GREG_CAL
) -> tuple[int, int, int, int, int, float]:
    """
    Convert Julian Day in Ephemeris Time (TT/ET) to UTC date/time.

    This is the inverse of utc_to_jd(): it converts a Julian Day number
    in Terrestrial Time (TT, formerly known as Ephemeris Time/ET) back to
    a UTC calendar date and time. The conversion properly accounts for
    Delta-T and leap seconds.

    Args:
        jd_et: Julian Day number in TT (Terrestrial Time / Ephemeris Time)
        calendar: GREG_CAL (1) for Gregorian, JUL_CAL (0) for Julian

    Returns:
        tuple: (year, month, day, hour, minute, second) where:
            - year: Calendar year (negative for BCE)
            - month: Month (1-12)
            - day: Day of month (1-31)
            - hour: Hour (0-23)
            - minute: Minute (0-59)
            - second: Second (0.0-59.999..., or 60.x during leap second)

    Note:
        - TT (Terrestrial Time) is the modern successor to Ephemeris Time (ET)
        - TT = TAI + 32.184 seconds, where TAI is International Atomic Time
        - UTC may include leap seconds (second = 60) on certain dates
        - Delta-T (TT - UT1) is automatically applied using IERS data

    Example:
        >>> from libephemeris import jdet_to_utc, utc_to_jd, GREG_CAL
        >>> # Convert J2000.0 epoch (JD 2451545.0 TT) to UTC
        >>> year, month, day, hour, minute, second = jdet_to_utc(2451545.0)
        >>> print(f"{year}-{month:02d}-{day:02d} {hour:02d}:{minute:02d}:{second:05.2f}")
        2000-01-01 11:58:55.82
    """
    # Reference-API parity: validate the calendar flag like julday()/revjul()
    _validate_calendar(calendar, "jdet_to_utc")

    # Pre-1972 there is no UTC: the reference API returns the UT1 calendar
    # date directly (jd_ut1 = jd_et - Delta T, fixed-point refined).
    jd_ut1_est = jd_et - deltat(jd_et)
    jd_ut1_est = jd_et - deltat(jd_ut1_est)
    # Fast path only when the instant is clearly before the UTC epoch. Within a
    # ~2 s band of 1972-01-01 the UT1 year is unreliable (UT1 leads/lags UTC by
    # up to 0.9 s), so defer to Skyfield's UTC below and classify on the
    # reconstructed UTC year — this keeps the UTC->JD->UTC round-trip exact
    # across the boundary (UTC 1972-01-01 00:00:00.0 has UT1 in year 1971).
    near_utc_epoch = abs(jd_ut1_est - _UTC_EPOCH_JD) < _UTC_EPOCH_BAND
    if revjul(jd_ut1_est, GREG_CAL)[0] < 1972 and not near_utc_epoch:
        return _jd_to_calendar_tuple(jd_ut1_est, calendar)

    ts = get_timescale()

    # Create a Skyfield Time object from the Julian Day. In the boundary band
    # reconstruct from the library-consistent UT1 estimate: Skyfield's UT1->UTC
    # lands on the correct side of the epoch, whereas the TT->UTC chain can
    # underflow it by a sub-microsecond and flip the calendar day back to 1971.
    t = ts.ut1_jd(jd_ut1_est) if near_utc_epoch else ts.tt_jd(jd_et)

    # Get UTC components from Skyfield (handles leap seconds automatically)
    # The .utc attribute returns a tuple: (year, month, day, hour, minute, second)
    # We cast to Any to work around Skyfield's reify decorator type annotation issues
    utc_data: Any = t.utc
    # Guard-band instants that Skyfield still resolves before 1972 have no UTC:
    # return the UT1 calendar date directly like the reference API.
    if int(utc_data[0]) < 1972:
        return _jd_to_calendar_tuple(jd_ut1_est, calendar)
    g_year = int(utc_data[0])
    g_month = int(utc_data[1])
    g_day = int(utc_data[2])
    g_hour = int(utc_data[3])
    g_minute = int(utc_data[4])
    g_second = float(utc_data[5])

    if calendar == JUL_CAL:
        # Only the calendar DATE differs between Gregorian and Julian (a
        # whole number of days), so remap the date at noon and carry the
        # clock components over verbatim. Folding h/m/s into a decimal hour
        # would roll a 23:59:60 leap second past 24h into the next day,
        # dropping the leap second (the reference preserves it).
        jd_noon = julday(g_year, g_month, g_day, 12.0, GREG_CAL)
        j_year, j_month, j_day, _ = revjul(jd_noon, JUL_CAL)
        return j_year, j_month, j_day, g_hour, g_minute, g_second

    return g_year, g_month, g_day, g_hour, g_minute, g_second


def jdut1_to_utc(
    jd_ut1: float, calendar: int = GREG_CAL
) -> tuple[int, int, int, int, int, float]:
    """
    Convert Julian Day in UT1 (Universal Time) to UTC date/time.

    Converts a Julian Day number in UT1 back to a UTC calendar date and time.
    The difference between UT1 and UTC is always less than 0.9 seconds by
    definition (maintained by adding leap seconds to UTC).

    Args:
        jd_ut1: Julian Day number in UT1 (Universal Time)
        calendar: GREG_CAL (1) for Gregorian, JUL_CAL (0) for Julian

    Returns:
        tuple: (year, month, day, hour, minute, second) where:
            - year: Calendar year (negative for BCE)
            - month: Month (1-12)
            - day: Day of month (1-31)
            - hour: Hour (0-23)
            - minute: Minute (0-59)
            - second: Second (0.0-59.999..., or 60.x during leap second)

    Note:
        - UT1 is based on Earth's rotation and is not perfectly uniform
        - UTC is atomic time adjusted to stay within 0.9s of UT1
        - The difference DUT1 = UT1 - UTC is published by IERS
        - For high-precision astronomical work, this difference matters

    Example:
        >>> from libephemeris import jdut1_to_utc, utc_to_jd, GREG_CAL
        >>> # Get JD(UT1) for a date, then convert back
        >>> jd_tt, jd_ut1 = utc_to_jd(2020, 6, 15, 14, 30, 0.0)
        >>> year, month, day, hour, minute, second = jdut1_to_utc(jd_ut1)
        >>> print(f"{year}-{month:02d}-{day:02d} {hour:02d}:{minute:02d}:{second:05.2f}")
        2020-06-15 14:30:00.00
    """
    # Reference-API parity: validate the calendar flag like julday()/revjul()
    _validate_calendar(calendar, "jdut1_to_utc")

    ts = get_timescale()

    # Pre-1972 there is no UTC: return the UT1 calendar date directly. Within a
    # ~2 s band of 1972-01-01 the UT1 year is unreliable (UT1 leads/lags UTC by
    # up to 0.9 s), so defer to Skyfield's UTC below and classify on the
    # reconstructed UTC year — this keeps the UTC->JD->UTC round-trip exact
    # across the boundary (UTC 1972-01-01 00:00:00.0 has UT1 in year 1971).
    near_utc_epoch = abs(jd_ut1 - _UTC_EPOCH_JD) < _UTC_EPOCH_BAND
    if revjul(jd_ut1, GREG_CAL)[0] < 1972 and not near_utc_epoch:
        return _jd_to_calendar_tuple(jd_ut1, calendar)

    # Create a Skyfield Time object from UT1 Julian Day
    t = ts.ut1_jd(jd_ut1)

    # Get UTC components from Skyfield (handles leap seconds automatically)
    # The .utc attribute returns a tuple: (year, month, day, hour, minute, second)
    # We cast to Any to work around Skyfield's reify decorator type annotation issues
    utc_data: Any = t.utc
    # Guard-band instants that Skyfield still resolves before 1972 have no UTC:
    # return the UT1 calendar date directly like the reference API.
    if int(utc_data[0]) < 1972:
        return _jd_to_calendar_tuple(jd_ut1, calendar)
    g_year = int(utc_data[0])
    g_month = int(utc_data[1])
    g_day = int(utc_data[2])
    g_hour = int(utc_data[3])
    g_minute = int(utc_data[4])
    g_second = float(utc_data[5])

    if calendar == JUL_CAL:
        # Only the calendar DATE differs between Gregorian and Julian (a
        # whole number of days), so remap the date at noon and carry the
        # clock components over verbatim. Folding h/m/s into a decimal hour
        # would roll a 23:59:60 leap second past 24h into the next day,
        # dropping the leap second (the reference preserves it).
        jd_noon = julday(g_year, g_month, g_day, 12.0, GREG_CAL)
        j_year, j_month, j_day, _ = revjul(jd_noon, JUL_CAL)
        return j_year, j_month, j_day, g_hour, g_minute, g_second

    return g_year, g_month, g_day, g_hour, g_minute, g_second


def day_of_week(jd: float) -> int:
    """
    Calculate the day of the week for a given Julian Day number.

    Uses the formula: floor(jd + 0.5) % 7 to get 0=Monday convention.
    This matches the reference day_of_week function.

    Args:
        jd: Julian Day number

    Returns:
        int: Day of the week where:
            - 0 = Monday
            - 1 = Tuesday
            - 2 = Wednesday
            - 3 = Thursday
            - 4 = Friday
            - 5 = Saturday
            - 6 = Sunday

    Example:
        >>> from libephemeris import day_of_week, julday
        >>> # J2000.0 epoch: Jan 1, 2000 was a Saturday
        >>> jd = julday(2000, 1, 1, 12.0)
        >>> day_of_week(jd)
        5
    """
    import math

    return int(math.floor(jd + 0.5)) % 7


def time_equ(jd: float) -> float:
    """
    Calculate the Equation of Time.

    The Equation of Time is the difference between apparent solar time
    (sundial time) and mean solar time (clock time). It arises from:
    1. Earth's orbital eccentricity (elliptical orbit)
    2. The obliquity of the ecliptic (Earth's axial tilt)

    Args:
        jd: Julian Day number in UT

    Returns:
        float: Equation of Time in days (multiply by 1440 for minutes)
               Positive values mean the sundial is ahead of the clock.

    Note:
        The equation of time varies from approximately -14 to +16 minutes
        throughout the year. The extremes occur around:
        - Early November: ~+16 minutes (sundial ahead)
        - Early February: ~-14 minutes (sundial behind)
        - Mid-April, mid-June, early September, late December: ~0 minutes

    Example:
        >>> from libephemeris import time_equ, julday
        >>> # Calculate equation of time for Jan 1, 2000
        >>> jd = julday(2000, 1, 1, 12.0)
        >>> eot = time_equ(jd)
        >>> eot_minutes = eot * 1440
        >>> print(f"Equation of Time: {eot_minutes:.2f} minutes")
        Equation of Time: -3.29 minutes
    """
    # The equation of time is derived from the relationship:
    #   E = GAST - RA_sun + 12h - UT
    # where GAST is Greenwich Apparent Sidereal Time (hours),
    # RA_sun is the Sun's apparent right ascension (hours),
    # and UT is the UT time of day (hours).
    # The result is normalized to ±12 hours, then converted to days.

    # Get Sun's apparent right ascension via equatorial coordinates
    # Import here to avoid circular import
    from . import calc_ut

    sun_result = calc_ut(jd, 0, FLG_EQUATORIAL)
    ra_hours = sun_result[0][0] / 15.0  # RA in degrees -> hours

    # Greenwich Apparent Sidereal Time (in hours)
    gast = sidtime(jd)

    # UT time of day in hours
    ut_hours = math.fmod((jd + 0.5), 1.0) * 24.0

    # Equation of time in hours
    e_hours = math.fmod(gast - ra_hours + 12.0 - ut_hours, 24.0)
    if e_hours > 12.0:
        e_hours -= 24.0
    elif e_hours < -12.0:
        e_hours += 24.0

    return e_hours / 24.0


def lat_to_lmt(jd_lat: float, longitude: float) -> float:
    """
    Convert Local Apparent Time (LAT) to Local Mean Time (LMT).

    Local Apparent Time (LAT) is the time as shown by a sundial - it is
    the true solar time based on the actual position of the Sun. Local
    Mean Time (LMT) is the mean solar time for a specific longitude,
    which runs at a uniform rate.

    The difference between LAT and LMT is the Equation of Time, which
    accounts for Earth's elliptical orbit and axial tilt.

    Args:
        jd_lat: Julian Day number in Local Apparent Time (sundial time)
        longitude: Geographic longitude in degrees (positive East, negative West)

    Returns:
        float: Julian Day number in Local Mean Time

    Note:
        - The longitude parameter is used to calculate the proper Julian Day
          for the Equation of Time lookup. Since the input is already in
          Local Apparent Time (which is inherently local to the observer's
          position), the longitude is used to convert to UT for the EoT
          calculation, then the result is converted back to local time.
        - The Equation of Time varies from approximately -14 to +16 minutes
          throughout the year.

    Example:
        >>> from libephemeris import lat_to_lmt, julday
        >>> # Convert LAT to LMT for Rome (12.5°E longitude)
        >>> jd_lat = julday(2000, 6, 15, 12.0)  # Noon LAT
        >>> jd_lmt = lat_to_lmt(jd_lat, 12.5)
        >>> # The difference should be approximately the Equation of Time
    """
    # Convert longitude to time offset in days (360° = 24h = 1 day)
    # Positive longitude (East) means local time is ahead of UT
    longitude_offset_days = longitude / 360.0

    # lat_to_lmt is the exact inverse of lmt_to_lat, which maps
    #     LAT = LMT + E(UT),   with   UT = LMT - longitude_offset
    # Recovering LMT from LAT therefore needs the Equation of Time E evaluated
    # at the UT of the (unknown) LMT instant, not at the UT of LAT. A single
    # evaluation at UT ≈ LAT - offset is off by up to ~0.2s in late November,
    # where dE/dt is steepest. Refine the UT estimate with three fixed-point
    # iterations (E drifts <30s/day, so this converges to well under a
    # millisecond); the results agree with reference-API output to sub-millisecond.
    ut_est = jd_lat - longitude_offset_days
    eot = 0.0
    for _ in range(3):
        eot = time_equ(ut_est)
        ut_est = jd_lat - longitude_offset_days - eot

    # LMT = LAT - EoT, with EoT taken at the converged UT of the LMT instant.
    # (When the sundial is ahead, EoT is positive, so LMT is less than LAT.)
    jd_lmt = jd_lat - eot

    return jd_lmt


def lmt_to_lat(jd_lmt: float, longitude: float) -> float:
    """
    Convert Local Mean Time (LMT) to Local Apparent Time (LAT).

    This is the inverse operation of lat_to_lmt(). Local Mean Time (LMT) is
    the mean solar time for a specific longitude, which runs at a uniform
    rate. Local Apparent Time (LAT) is the time as shown by a sundial - it
    is the true solar time based on the actual position of the Sun.

    The difference between LAT and LMT is the Equation of Time, which
    accounts for Earth's elliptical orbit and axial tilt.

    Args:
        jd_lmt: Julian Day number in Local Mean Time
        longitude: Geographic longitude in degrees (positive East, negative West)

    Returns:
        float: Julian Day number in Local Apparent Time (sundial time)

    Note:
        - The longitude parameter is used to calculate the proper Julian Day
          for the Equation of Time lookup. Since the input is in Local Mean
          Time, the longitude is used to convert to UT for the EoT
          calculation, then the result is converted back to local time.
        - The Equation of Time varies from approximately -14 to +16 minutes
          throughout the year.

    Example:
        >>> from libephemeris import lmt_to_lat, julday
        >>> # Convert LMT to LAT for Rome (12.5°E longitude)
        >>> jd_lmt = julday(2000, 6, 15, 12.0)  # Noon LMT
        >>> jd_lat = lmt_to_lat(jd_lmt, 12.5)
        >>> # The difference should be approximately the Equation of Time
    """
    # Convert longitude to time offset in days (360° = 24h = 1 day)
    # Positive longitude (East) means local time is ahead of UT
    longitude_offset_days = longitude / 360.0

    # Convert from local time to UT for EoT calculation
    # LMT (local) = UT + longitude_offset
    # So UT ≈ LMT - longitude_offset
    jd_ut_approx = jd_lmt - longitude_offset_days

    # Calculate the Equation of Time at this UT
    eot = time_equ(jd_ut_approx)

    # LAT = LMT + EoT
    # (This is the inverse of LMT = LAT - EoT)
    jd_lat = jd_lmt + eot

    return jd_lat


def _sidtime_internal(
    jd: float, longitude: float, obliquity: float, nutation: float
) -> float:
    """
    Internal function to calculate Local Sidereal Time with all parameters.

    This is the full implementation used by both sidtime() and sidtime0().
    For most applications, use sidtime(jd) or sidtime0(jd, eps, nut) instead.

    GMST comes from the library's single long-term (Vondrák 2011) sidereal-time
    source, the same one the house/ARMC engine uses; inside the modern window
    (1850-2050) that source is the IAU 2006 GMST expression (Capitaine et al.
    2003), so modern results are unchanged.

    Args:
        jd: Julian Day number in UT (Universal Time)
        longitude: Geographic longitude in degrees (positive East, negative West)
        obliquity: Obliquity of the ecliptic in degrees (typically ~23.44°)
        nutation: Nutation in longitude in degrees (typically small, ~±0.005°)

    Returns:
        float: Local sidereal time in hours (0.0 to 24.0)
    """
    import math

    # Compute GMST using IAU 2006 formula (requires both UT1 and TT)
    gmst_rad = _gmst06(jd)
    gmst_hours = math.degrees(gmst_rad) / 15.0

    # Apply equation of equinoxes to get Greenwich Apparent Sidereal Time (GAST)
    # Equation of equinoxes = nutation in longitude * cos(obliquity)
    obliquity_rad = math.radians(obliquity)
    # Nutation in longitude (degrees) * cos(obliquity) -> degrees, /15 -> hours
    equation_of_equinoxes = (nutation * math.cos(obliquity_rad)) / 15.0

    gast = gmst_hours + equation_of_equinoxes

    # Convert longitude to hours (15° = 1 hour)
    longitude_hours = longitude / 15.0

    # Local Apparent Sidereal Time
    last = gast + longitude_hours

    # Normalize to 0-24 hours range
    last = last % 24.0
    if last < 0:  # pragma: no cover - modulo by positive 24.0 is always >= 0
        last += 24.0

    return last


def _gmst06(jd_ut1: float) -> float:
    """Compute Greenwich Mean Sidereal Time (radians).

    Routes through the library's single long-term GMST source,
    :func:`libephemeris.sidereal_longterm.mean_sidereal_time_deg` (Vondrák
    2011) — the same realization the house/ARMC engine uses. Inside the modern
    window (1850-2050) that source *is* the IAU 2006 GMST expression (Capitaine
    et al. 2003), so modern results are unchanged; outside it the long-term
    geometric construction keeps public sidereal time consistent with the
    house cusps at remote epochs, where an IAU 2006 precession-in-RA polynomial
    diverges (by ~0.66° at year -3000).

    Args:
        jd_ut1: Julian Day number in UT1.

    Returns:
        GMST in radians in [0, 2*pi).
    """
    from .sidereal_longterm import mean_sidereal_time_deg

    return math.radians(mean_sidereal_time_deg(jd_ut1))


def sidtime(
    jd: float,
    longitude: float = 0.0,
    obliquity: Optional[float] = None,
    nutation: Optional[float] = None,
) -> float:
    """
    Calculate Local Apparent Sidereal Time for a given Julian Day and longitude.

    This function provides a flexible interface for sidereal time calculation.
    It can be called with just the Julian Day (auto-calculating nutation), or
    with explicit obliquity and nutation values.

    Sidereal time is the hour angle of the vernal equinox (First Point of Aries)
    and is used to determine which stars are visible at a given time.

    Args:
        jd: Julian Day number in UT (Universal Time)
        longitude: Geographic longitude in degrees (positive East, negative West).
                   Default is 0.0 (Greenwich).
        obliquity: Obliquity of the ecliptic in degrees (typically ~23.44°).
                   If None, calculated automatically from ephemeris.
        nutation: Nutation in longitude in degrees (typically small, ~±0.005°).
                  If None, calculated automatically from ephemeris.

    Returns:
        float: Local Apparent Sidereal Time in hours (0.0 to 24.0)

    Note:
        - When obliquity/nutation are None, uses IAU 2000B nutation model
        - GMST comes from the long-term (Vondrák 2011) sidereal-time source
          (the IAU 2006 GMST expression inside 1850-2050); this keeps sidtime
          consistent with the house/ARMC engine at remote epochs
        - Result is normalized to the range 0-24 hours
        - For Greenwich sidereal time, use longitude=0.0

    Example:
        >>> from libephemeris import sidtime
        >>> # Calculate GST for J2000.0 (Jan 1, 2000 at noon)
        >>> gst = sidtime(2451545.0)
        >>> print(f"GST: {gst:.4f} hours")
        GST: 18.6974 hours
        >>> # Calculate LST with explicit values
        >>> lst = sidtime(2451545.0, 0.0, 23.44, 0.0)
    """
    # If obliquity or nutation not provided, get from ephemeris
    if obliquity is None or nutation is None:
        # Import here to avoid circular imports
        from .planets import calc_ut
        from .constants import ECL_NUT

        # Get nutation and obliquity
        nut_data, _ = calc_ut(jd, ECL_NUT, 0)
        # nut_data[0] = true obliquity, nut_data[1] = mean obliquity
        # nut_data[2] = nutation in longitude, nut_data[3] = nutation in obliquity
        if obliquity is None:
            obliquity = nut_data[0]  # True obliquity
        if nutation is None:
            nutation = nut_data[2]  # Nutation in longitude

    return _sidtime_internal(jd, longitude, obliquity, nutation)


def sidtime0(jd: float, obliquity: float, nutation: float) -> float:
    """
    Calculate Greenwich Sidereal Time for a given Julian Day.

    This is the sidereal time at Greenwich (longitude 0°) and is the base
    for calculating local sidereal time at any other longitude. Compatible
    with sidtime0().

    The calculation uses the IAU formula for Greenwich Mean Sidereal Time (GMST)
    and applies the equation of equinoxes to get Greenwich Apparent Sidereal
    Time (GAST).

    Args:
        jd: Julian Day number in UT (Universal Time)
        obliquity: Obliquity of the ecliptic in degrees (typically ~23.44°)
        nutation: Nutation in longitude in degrees (typically small, ~±0.005°)

    Returns:
        float: Greenwich sidereal time in hours (0.0 to 24.0)

    Note:
        - This function is equivalent to calling _sidtime_internal(jd, 0.0, obliquity, nutation)
        - GMST comes from the long-term (Vondrák 2011) sidereal-time source
          (the IAU 2006 GMST expression inside 1850-2050)
        - The equation of equinoxes = nutation_in_longitude * cos(obliquity)
        - Result is normalized to the range 0-24 hours

    Example:
        >>> from libephemeris import sidtime0
        >>> # Calculate GST for J2000.0 (Jan 1, 2000 at noon)
        >>> gst = sidtime0(2451545.0, 23.4393, 0.0)
        >>> print(f"GST: {gst:.4f} hours")
        GST: 18.6974 hours
    """
    return _sidtime_internal(jd, 0.0, obliquity, nutation)


# =============================================================================
# TAI (INTERNATIONAL ATOMIC TIME) FUNCTIONS
# =============================================================================

# TT - TAI offset: TT = TAI + 32.184 seconds
TT_TAI_OFFSET_SECONDS = 32.184
TT_TAI_OFFSET_DAYS = TT_TAI_OFFSET_SECONDS / 86400.0


def get_tai_utc_for_jd(jd: float) -> float:
    """
    Get TAI-UTC (leap seconds) for a given Julian Day.

    TAI (International Atomic Time) runs continuously without leap seconds,
    while UTC is adjusted by leap seconds to stay within 0.9 seconds of UT1.
    This function returns the cumulative number of leap seconds at the given
    Julian Day.

    Args:
        jd: Julian Day number (any time scale, typically UT1 or UTC)

    Returns:
        float: TAI-UTC in seconds (the cumulative leap seconds)

    Note:
        - Before 1972, the TAI-UTC relationship was more complex and this
          function returns approximate values based on the historical record.
        - After 1972, TAI-UTC increases by exactly 1 second with each leap
          second, occurring on June 30 or December 31.
        - As of 2017, TAI-UTC = 37 seconds.

    Example:
        >>> from libephemeris import get_tai_utc_for_jd, julday
        >>> # Get TAI-UTC for Jan 1, 2020
        >>> jd = julday(2020, 1, 1, 12.0)
        >>> tai_utc = get_tai_utc_for_jd(jd)
        >>> print(f"TAI-UTC: {tai_utc:.0f} seconds")
        TAI-UTC: 37 seconds
    """
    from .iers_data import get_tai_utc, _jd_to_mjd

    mjd = _jd_to_mjd(jd)
    return get_tai_utc(mjd)


def utc_to_tai_jd(
    year: int,
    month: int,
    day: int,
    hour: int,
    minute: int,
    second: float,
    calendar: int = GREG_CAL,
) -> float:
    """
    Convert UTC date/time to Julian Day number in TAI (International Atomic Time).

    TAI is a continuous atomic time scale that does not have leap seconds.
    It is the basis for other atomic time scales like TT (Terrestrial Time).

    Args:
        year: Calendar year (negative for BCE)
        month: Month (1-12)
        day: Day of month (1-31)
        hour: Hour (0-23)
        minute: Minute (0-59)
        second: Second (0-60, allowing for leap seconds in UTC)
        calendar: GREG_CAL (1) for Gregorian, JUL_CAL (0) for Julian

    Returns:
        float: Julian Day number in TAI

    Note:
        - TAI = UTC + TAI-UTC (leap seconds)
        - TAI runs ahead of UTC by the cumulative number of leap seconds
        - As of 2017, TAI is 37 seconds ahead of UTC

    Example:
        >>> from libephemeris import utc_to_tai_jd
        >>> # Convert Jan 1, 2020 00:00:00 UTC to TAI JD
        >>> jd_tai = utc_to_tai_jd(2020, 1, 1, 0, 0, 0.0)
        >>> print(f"JD(TAI): {jd_tai:.6f}")
    """
    # First get the JD in UTC (using UT1 as approximation since UTC ≈ UT1)
    decimal_hour = hour + minute / 60.0 + second / 3600.0
    jd_utc = julday(year, month, day, decimal_hour, calendar)

    # Get TAI-UTC offset (leap seconds) for this date
    tai_utc_seconds = get_tai_utc_for_jd(jd_utc)

    # TAI = UTC + TAI-UTC
    jd_tai = jd_utc + tai_utc_seconds / 86400.0

    return jd_tai


def tai_jd_to_utc(
    jd_tai: float, calendar: int = GREG_CAL
) -> tuple[int, int, int, int, int, float]:
    """
    Convert Julian Day in TAI (International Atomic Time) to UTC date/time.

    This is the inverse of utc_to_tai_jd(). It converts a Julian Day number
    in TAI back to a UTC calendar date and time.

    Args:
        jd_tai: Julian Day number in TAI
        calendar: GREG_CAL (1) for Gregorian, JUL_CAL (0) for Julian

    Returns:
        tuple: (year, month, day, hour, minute, second) in UTC where:
            - year: Calendar year (negative for BCE)
            - month: Month (1-12)
            - day: Day of month (1-31)
            - hour: Hour (0-23)
            - minute: Minute (0-59)
            - second: Second (0.0-59.999...)

    Note:
        - UTC = TAI - TAI-UTC (leap seconds)
        - The conversion requires knowing the leap seconds at the target time
        - This uses an iterative approach since TAI-UTC depends on the UTC date

    Example:
        >>> from libephemeris import tai_jd_to_utc, utc_to_tai_jd
        >>> # Round-trip test
        >>> jd_tai = utc_to_tai_jd(2020, 1, 1, 0, 0, 0.0)
        >>> year, month, day, hour, minute, second = tai_jd_to_utc(jd_tai)
        >>> print(f"{year}-{month:02d}-{day:02d} {hour:02d}:{minute:02d}:{second:05.2f}")
        2020-01-01 00:00:00.00
    """
    # First approximation: use the TAI JD to estimate the UTC date
    # This is close since TAI-UTC is only tens of seconds
    tai_utc_seconds = get_tai_utc_for_jd(jd_tai)

    # UTC = TAI - TAI-UTC
    jd_utc_approx = jd_tai - tai_utc_seconds / 86400.0

    # Refine: get the leap seconds at the approximate UTC date
    # (In case we crossed a leap second boundary)
    tai_utc_seconds_refined = get_tai_utc_for_jd(jd_utc_approx)
    jd_utc = jd_tai - tai_utc_seconds_refined / 86400.0

    # Convert JD to calendar date
    year, month, day, decimal_hour = revjul(jd_utc, calendar)

    # Extract time components
    hour = int(decimal_hour)
    minute_frac = (decimal_hour - hour) * 60.0
    minute = int(minute_frac)
    second = (minute_frac - minute) * 60.0

    return year, month, day, hour, minute, second


def tt_to_tai_jd(jd_tt: float) -> float:
    """
    Convert Julian Day in TT (Terrestrial Time) to TAI (International Atomic Time).

    The relationship between TT and TAI is fixed and exact:
    TT = TAI + 32.184 seconds

    Args:
        jd_tt: Julian Day number in TT (Terrestrial Time)

    Returns:
        float: Julian Day number in TAI

    Note:
        - TT = TAI + 32.184 seconds (exactly)
        - Therefore TAI = TT - 32.184 seconds
        - This offset was chosen so that TT would be continuous with ET
          (Ephemeris Time) at the moment of the definition

    Example:
        >>> from libephemeris import tt_to_tai_jd
        >>> # J2000.0 epoch in TT
        >>> jd_tt = 2451545.0
        >>> jd_tai = tt_to_tai_jd(jd_tt)
        >>> print(f"JD(TT): {jd_tt:.6f}, JD(TAI): {jd_tai:.9f}")
        >>> print(f"Difference: {(jd_tt - jd_tai) * 86400:.3f} seconds")
    """
    # TAI = TT - 32.184 seconds
    return jd_tt - TT_TAI_OFFSET_DAYS


def tai_to_tt_jd(jd_tai: float) -> float:
    """
    Convert Julian Day in TAI (International Atomic Time) to TT (Terrestrial Time).

    The relationship between TT and TAI is fixed and exact:
    TT = TAI + 32.184 seconds

    Args:
        jd_tai: Julian Day number in TAI

    Returns:
        float: Julian Day number in TT (Terrestrial Time)

    Note:
        - TT = TAI + 32.184 seconds (exactly)
        - This offset was defined when TT replaced ET in 1984
        - TT is used for planetary ephemeris calculations

    Example:
        >>> from libephemeris import tai_to_tt_jd, tt_to_tai_jd
        >>> # Round-trip test
        >>> jd_tt_orig = 2451545.0  # J2000.0
        >>> jd_tai = tt_to_tai_jd(jd_tt_orig)
        >>> jd_tt_back = tai_to_tt_jd(jd_tai)
        >>> print(f"Original: {jd_tt_orig:.10f}")
        >>> print(f"Recovered: {jd_tt_back:.10f}")
    """
    # TT = TAI + 32.184 seconds
    return jd_tai + TT_TAI_OFFSET_DAYS


def tai_to_utc_jd(jd_tai: float) -> float:
    """
    Convert Julian Day in TAI (International Atomic Time) to UTC Julian Day.

    Args:
        jd_tai: Julian Day number in TAI

    Returns:
        float: Julian Day number in UTC

    Note:
        - UTC = TAI - TAI-UTC (leap seconds)
        - This is a simpler version of tai_jd_to_utc() that returns JD directly

    Example:
        >>> from libephemeris import tai_to_utc_jd, utc_to_tai_jd, julday
        >>> jd_utc = julday(2020, 1, 1, 12.0)
        >>> jd_tai = utc_to_tai_jd(2020, 1, 1, 12, 0, 0.0)
        >>> jd_utc_back = tai_to_utc_jd(jd_tai)
        >>> print(f"Difference: {abs(jd_utc - jd_utc_back) * 86400:.6f} seconds")
    """
    # First approximation
    tai_utc_seconds = get_tai_utc_for_jd(jd_tai)
    jd_utc_approx = jd_tai - tai_utc_seconds / 86400.0

    # Refine in case we crossed a leap second boundary
    tai_utc_seconds_refined = get_tai_utc_for_jd(jd_utc_approx)
    jd_utc = jd_tai - tai_utc_seconds_refined / 86400.0

    return jd_utc


def utc_jd_to_tai(jd_utc: float) -> float:
    """
    Convert Julian Day in UTC to TAI (International Atomic Time) Julian Day.

    Args:
        jd_utc: Julian Day number in UTC

    Returns:
        float: Julian Day number in TAI

    Note:
        - TAI = UTC + TAI-UTC (leap seconds)
        - This is a simpler version of utc_to_tai_jd() that takes JD directly

    Example:
        >>> from libephemeris import utc_jd_to_tai, julday
        >>> jd_utc = julday(2020, 1, 1, 12.0)
        >>> jd_tai = utc_jd_to_tai(jd_utc)
        >>> print(f"JD(UTC): {jd_utc:.6f}, JD(TAI): {jd_tai:.6f}")
        >>> print(f"TAI is ahead by: {(jd_tai - jd_utc) * 86400:.0f} seconds")
    """
    tai_utc_seconds = get_tai_utc_for_jd(jd_utc)
    return jd_utc + tai_utc_seconds / 86400.0


def utc_time_zone(
    year: int,
    month: int,
    day: int,
    hour: int,
    minute: int,
    second: float,
    timezone_offset: float,
) -> tuple[int, int, int, int, int, float]:
    """
    Convert local time to UTC by subtracting the timezone offset.

    Given a date/time in a specific timezone, subtract the timezone offset
    to obtain the equivalent UTC date/time. Handles all date boundary
    crossings (day, month, year) correctly.

    This matches the reference convention: the offset is **subtracted**
    from the input time.

    Args:
        year: Calendar year (negative for BCE)
        month: Month (1-12)
        day: Day of month (1-31)
        hour: Hour (0-23)
        minute: Minute (0-59)
        second: Second (0.0-59.999...)
        timezone_offset: Timezone offset in hours from UTC.
            Positive values for timezones east of UTC (e.g., +1 for CET, +9 for JST)
            Negative values for timezones west of UTC (e.g., -5 for EST, -8 for PST)

    Returns:
        tuple: (year, month, day, hour, minute, second) in UTC where:
            - year: Calendar year
            - month: Month (1-12)
            - day: Day of month (1-31)
            - hour: Hour (0-23)
            - minute: Minute (0-59)
            - second: Second (0.0-59.999...)

    Example:
        >>> from libephemeris import utc_time_zone
        >>> # Convert 2024-01-15 12:00:00 CET (UTC+1) to UTC
        >>> utc_time_zone(2024, 1, 15, 12, 0, 0.0, 1)
        (2024, 1, 15, 11, 0, 0.0)
        >>> # Convert 2024-01-15 12:00:00 EST (UTC-5) to UTC
        >>> utc_time_zone(2024, 1, 15, 12, 0, 0.0, -5)
        (2024, 1, 15, 17, 0, 0.0)
    """
    # Component-wise arithmetic (no decimal-hour round trip): the seconds
    # value passes through bit-exact — including a leap-second input
    # (second >= 60), which the reference API preserves — and no
    # millisecond rounding is applied.
    offset_min_f = timezone_offset * 60.0
    offset_min = round(offset_min_f)
    # Sub-minute offset residue (rare, e.g. historical LMT zones)
    offset_sec = (offset_min_f - offset_min) * 60.0

    out_second = second - offset_sec
    minute_carry = 0
    # Keep a leap second (60 <= s < 61) untouched by normalization
    if out_second >= 61.0 or (out_second >= 60.0 and second < 60.0):
        out_second -= 60.0
        minute_carry += 1
    elif out_second < 0.0:
        out_second += 60.0
        minute_carry -= 1

    total_min = hour * 60 + minute - offset_min + minute_carry
    day_shift, total_min = divmod(total_min, 1440)

    local_hour = total_min // 60
    local_minute = total_min % 60

    if day_shift == 0:
        local_year, local_month, local_day = year, month, day
    else:
        # Date arithmetic via JD at noon (immune to time-of-day rounding)
        jd_noon = julday(year, month, day, 12.0, GREG_CAL) + day_shift
        local_year, local_month, local_day, _ = revjul(jd_noon, GREG_CAL)

    return (
        local_year,
        local_month,
        local_day,
        int(local_hour),
        int(local_minute),
        float(out_second),
    )
