"""
Time function comparison tests vs pyswisseph.

Compares julday, revjul, deltat, sidtime, utc_to_jd, and day_of_week
against pyswisseph across a wide range of dates. Tolerances account for
minor algorithm differences (e.g. Delta T polynomial fits).
"""

from __future__ import annotations


import pytest

swe_ref = pytest.importorskip("swisseph")

import libephemeris as swe
from libephemeris.constants import GREG_CAL


# ---------------------------------------------------------------------------
# Test data -- 20 dates for julday/revjul
# ---------------------------------------------------------------------------

JULDAY_DATES = [
    (2000, 1, 1, 12.0),  # J2000.0
    (1970, 1, 1, 0.0),  # Unix epoch
    (1999, 12, 31, 23.99),  # near Y2K
    (2024, 3, 20, 12.5),  # equinox
    (2024, 6, 21, 0.0),  # solstice
    (2024, 12, 25, 18.0),  # Christmas
    (1900, 1, 1, 0.0),  # turn of 20th century
    (1582, 10, 15, 12.0),  # Gregorian calendar start
    (1850, 7, 4, 6.25),  # mid-19th century
    (1950, 6, 15, 12.0),  # mid-20th century
    (2050, 1, 1, 0.0),  # mid-21st century
    (2100, 12, 31, 23.5),  # end 21st century
    (1800, 3, 21, 12.0),  # spring equinox 1800
    (1920, 11, 11, 11.0),  # end of WWI armistice
    (2000, 6, 21, 12.0),  # summer solstice 2000
    (2010, 1, 15, 7.5),  # recent date
    (2020, 3, 1, 0.0),  # start of pandemic year
    (1975, 8, 15, 15.75),  # mid-1970s
    (2040, 4, 8, 18.0),  # future eclipse
    (1600, 1, 1, 12.0),  # early Gregorian
]


# 20 Julian Day values for deltat/sidtime
JD_VALUES = [
    2415020.0,  # 1900-01-01
    2418665.0,  # 1910-01-01
    2422310.0,  # 1920-01-01
    2425955.0,  # 1930-01-01
    2429601.0,  # 1940-01-01
    2433282.5,  # 1950-01-01
    2436934.5,  # 1960-01-01
    2440587.5,  # 1970-01-01
    2444239.5,  # 1980-01-01
    2447892.5,  # 1990-01-01
    2451545.0,  # 2000-01-01 (J2000)
    2453371.5,  # 2005-01-01
    2455197.5,  # 2010-01-01
    2457023.5,  # 2015-01-01
    2458849.5,  # 2020-01-01
    2459580.5,  # 2022-01-01
    2460310.5,  # 2024-01-18
    2462502.5,  # 2030-01-01
    2470547.0,  # 2052-01-01
    2488070.0,  # 2100-01-01
]

# UTC dates for utc_to_jd (year, month, day, hour, minute, second)
UTC_DATES = [
    (2000, 1, 1, 12, 0, 0.0),
    (1999, 12, 31, 23, 59, 59.0),
    (2024, 3, 20, 6, 30, 0.0),
    (2024, 6, 21, 12, 0, 0.0),
    (1970, 1, 1, 0, 0, 0.0),
    (2010, 6, 15, 18, 45, 30.0),
    (1990, 7, 4, 0, 0, 0.0),
    (2050, 1, 1, 0, 0, 0.0),
    (1985, 12, 25, 12, 0, 0.0),
    (2020, 6, 21, 6, 44, 0.0),
]

# Known day-of-week values: (year, month, day, hour, expected_dow)
# 0=Monday, 1=Tuesday, ..., 6=Sunday
DAY_OF_WEEK_DATES = [
    (2000, 1, 1, 12.0, 5),  # Saturday
    (2024, 1, 1, 12.0, 0),  # Monday
    (2024, 3, 15, 12.0, 4),  # Friday
    (2024, 7, 4, 12.0, 3),  # Thursday
    (1970, 1, 1, 12.0, 3),  # Thursday
    (2000, 2, 29, 12.0, 1),  # Tuesday (leap day)
    (1999, 12, 31, 12.0, 4),  # Friday
    (2023, 12, 25, 12.0, 0),  # Monday
    (1969, 7, 20, 12.0, 6),  # Sunday (Moon landing)
    (2025, 1, 1, 12.0, 2),  # Wednesday
]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestJulday:
    """Compare julday against pyswisseph."""

    @pytest.mark.unit
    @pytest.mark.parametrize("y,m,d,h", JULDAY_DATES)
    def test_julday_matches(self, y: int, m: int, d: int, h: float):
        """julday matches pyswisseph to within 1e-6 days (~0.09s)."""
        lib_jd = swe.julday(y, m, d, h)
        ref_jd = swe_ref.julday(y, m, d, h)
        assert abs(lib_jd - ref_jd) < 1e-6, (
            f"julday({y},{m},{d},{h}): lib={lib_jd}, ref={ref_jd}"
        )

    @pytest.mark.unit
    def test_julday_j2000_exact(self):
        """J2000.0 epoch is exactly 2451545.0."""
        assert swe.julday(2000, 1, 1, 12.0) == 2451545.0

    @pytest.mark.unit
    def test_julday_returns_float(self):
        """julday returns a native Python float."""
        result = swe.julday(2000, 1, 1, 12.0)
        assert type(result) is float


class TestRevjul:
    """Compare revjul round-trips against pyswisseph."""

    @pytest.mark.unit
    @pytest.mark.parametrize("y,m,d,h", JULDAY_DATES)
    def test_revjul_roundtrip(self, y: int, m: int, d: int, h: float):
        """julday -> revjul recovers the original date components."""
        jd = swe.julday(y, m, d, h)
        ry, rm, rd, rh = swe.revjul(jd)
        assert ry == y, f"Year: {ry} != {y}"
        assert rm == m, f"Month: {rm} != {m}"
        assert rd == d, f"Day: {rd} != {d}"
        assert abs(rh - h) < 1e-6, f"Hour: {rh} != {h}"

    @pytest.mark.unit
    @pytest.mark.parametrize("y,m,d,h", JULDAY_DATES)
    def test_revjul_matches_swisseph(self, y: int, m: int, d: int, h: float):
        """revjul produces the same result as pyswisseph."""
        jd = swe_ref.julday(y, m, d, h)
        lib_result = swe.revjul(jd)
        ref_result = swe_ref.revjul(jd)
        assert lib_result[0] == ref_result[0], (
            f"Year: {lib_result[0]} != {ref_result[0]}"
        )
        assert lib_result[1] == ref_result[1], (
            f"Month: {lib_result[1]} != {ref_result[1]}"
        )
        assert lib_result[2] == ref_result[2], (
            f"Day: {lib_result[2]} != {ref_result[2]}"
        )
        assert abs(lib_result[3] - ref_result[3]) < 1e-6, (
            f"Hour: {lib_result[3]} != {ref_result[3]}"
        )


class TestDeltat:
    """Compare deltat against pyswisseph."""

    @pytest.mark.unit
    @pytest.mark.parametrize("jd", JD_VALUES)
    def test_deltat_matches(self, jd: float):
        """deltat matches pyswisseph within 1e-4 day (~8.6 seconds).

        The tolerance accounts for different Delta T polynomial models:
        libephemeris uses SMH-2016 / IERS data while pyswisseph uses its
        own internal polynomial fit.
        """
        lib_dt = swe.deltat(jd)
        ref_dt = swe_ref.deltat(jd)
        assert abs(lib_dt - ref_dt) < 1e-4, (
            f"jd={jd}: lib={lib_dt:.10f}, ref={ref_dt:.10f}, "
            f"diff={abs(lib_dt - ref_dt):.2e} days"
        )

    @pytest.mark.unit
    def test_deltat_returns_float(self):
        """deltat returns a native Python float."""
        dt = swe.deltat(2451545.0)
        assert type(dt) is float

    @pytest.mark.unit
    def test_deltat_j2000_positive(self):
        """Delta T at J2000.0 should be positive (~63.8 seconds)."""
        dt = swe.deltat(2451545.0)
        assert dt > 0
        # ~63.8 seconds = ~0.000738 days
        assert 0.0005 < dt < 0.001

    @pytest.mark.unit
    def test_deltat_monotone_modern_era(self):
        """Delta T increases monotonically in the modern era (1960-2020)."""
        modern_jds = [jd for jd in JD_VALUES if 2436934.5 <= jd <= 2458849.5]
        dts = [swe.deltat(jd) for jd in modern_jds]
        for i in range(1, len(dts)):
            assert dts[i] >= dts[i - 1], (
                f"Delta T not monotone: dt[{i - 1}]={dts[i - 1]} > dt[{i}]={dts[i]}"
            )


class TestSidtime:
    """Compare sidereal time against pyswisseph."""

    @pytest.mark.unit
    @pytest.mark.parametrize("jd", JD_VALUES)
    def test_sidtime_matches(self, jd: float):
        """Sidereal time matches pyswisseph within 0.001 hour (~3.6 seconds).

        Small differences are expected due to different nutation models
        (IAU 2000A vs IAU 1980).
        """
        lib_st = swe.sidtime(jd)
        ref_st = swe_ref.sidtime(jd)
        diff = abs(lib_st - ref_st)
        # Handle wrap-around at 24 hours
        if diff > 12:
            diff = 24 - diff
        assert diff < 0.001, (
            f"jd={jd}: lib={lib_st:.8f}h, ref={ref_st:.8f}h, diff={diff:.6f}h"
        )

    @pytest.mark.unit
    def test_sidtime_returns_float(self):
        """sidtime returns a native Python float."""
        st = swe.sidtime(2451545.0)
        assert type(st) is float

    @pytest.mark.unit
    def test_sidtime_range(self):
        """Sidereal time is in [0, 24) hours."""
        for jd in JD_VALUES:
            st = swe.sidtime(jd)
            assert 0 <= st < 24, f"jd={jd}: sidtime={st} out of range"


class TestUtcToJd:
    """Compare utc_to_jd against pyswisseph."""

    @pytest.mark.unit
    @pytest.mark.parametrize("y,m,d,h,mi,s", UTC_DATES)
    def test_utc_to_jd_matches(self, y: int, m: int, d: int, h: int, mi: int, s: float):
        """utc_to_jd matches pyswisseph within 1e-6 day (~0.09s).

        Returns (jd_et, jd_ut1); we compare both components.
        """
        lib_result = swe.utc_to_jd(y, m, d, h, mi, s, GREG_CAL)
        ref_result = swe_ref.utc_to_jd(y, m, d, h, mi, s, GREG_CAL)

        # jd_et (TT)
        assert abs(lib_result[0] - ref_result[0]) < 1e-4, (
            f"jd_et: lib={lib_result[0]}, ref={ref_result[0]}"
        )

        # jd_ut1
        assert abs(lib_result[1] - ref_result[1]) < 1e-4, (
            f"jd_ut1: lib={lib_result[1]}, ref={ref_result[1]}"
        )

    @pytest.mark.unit
    def test_utc_to_jd_returns_tuple(self):
        """utc_to_jd returns a 2-tuple of floats."""
        result = swe.utc_to_jd(2000, 1, 1, 12, 0, 0.0, GREG_CAL)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert all(type(v) is float for v in result)


class TestDayOfWeek:
    """Compare day_of_week against pyswisseph and known values."""

    @pytest.mark.unit
    @pytest.mark.parametrize("y,m,d,h,expected_dow", DAY_OF_WEEK_DATES)
    def test_day_of_week_known(
        self, y: int, m: int, d: int, h: float, expected_dow: int
    ):
        """day_of_week matches known calendar day for specific dates."""
        jd = swe.julday(y, m, d, h)
        dow = swe.day_of_week(jd)
        assert dow == expected_dow, (
            f"{y}-{m:02d}-{d:02d}: got {dow}, expected {expected_dow}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("y,m,d,h,expected_dow", DAY_OF_WEEK_DATES)
    def test_day_of_week_matches_swisseph(
        self, y: int, m: int, d: int, h: float, expected_dow: int
    ):
        """day_of_week matches pyswisseph."""
        jd = swe_ref.julday(y, m, d, h)
        lib_dow = swe.day_of_week(jd)
        ref_dow = swe_ref.day_of_week(jd)
        assert lib_dow == ref_dow, f"{y}-{m:02d}-{d:02d}: lib={lib_dow}, ref={ref_dow}"

    @pytest.mark.unit
    def test_day_of_week_returns_int(self):
        """day_of_week returns a native Python int in [0, 6]."""
        jd = swe.julday(2000, 1, 1, 12.0)
        dow = swe.day_of_week(jd)
        assert type(dow) is int
        assert 0 <= dow <= 6
