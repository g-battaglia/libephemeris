"""
Comprehensive tests for time conversion functions.

Verifies julday, revjul, deltat, sidtime, utc_to_jd, jdut1_to_utc,
and jdet_to_utc work correctly with round-trips and edge cases.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe


class TestJuldayRevjul:
    """Tests for Julian Day conversion round-trips."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "y,m,d,h",
        [
            (2000, 1, 1, 12.0),
            (1999, 12, 31, 0.0),
            (2024, 7, 15, 18.5),
            (1900, 1, 1, 0.0),
            (1582, 10, 15, 12.0),
            (1, 1, 1, 12.0),
            (2100, 6, 21, 6.25),
            (-4712, 1, 1, 12.0),  # Julian Day epoch
        ],
    )
    def test_julday_revjul_roundtrip(self, y: int, m: int, d: int, h: float):
        """julday -> revjul should recover original date."""
        jd = swe.julday(y, m, d, h)
        ry, rm, rd, rh = swe.revjul(jd)
        assert ry == y, f"Year: {ry} != {y}"
        assert rm == m, f"Month: {rm} != {m}"
        assert rd == d, f"Day: {rd} != {d}"
        assert abs(rh - h) < 1e-6, f"Hour: {rh} != {h}"

    @pytest.mark.unit
    def test_j2000_jd(self):
        """J2000.0 should be JD 2451545.0."""
        jd = swe.julday(2000, 1, 1, 12.0)
        assert abs(jd - 2451545.0) < 1e-6

    @pytest.mark.unit
    def test_unix_epoch_jd(self):
        """Unix epoch (1970-01-01 00:00) should be JD 2440587.5."""
        jd = swe.julday(1970, 1, 1, 0.0)
        assert abs(jd - 2440587.5) < 1e-6

    @pytest.mark.unit
    def test_julday_returns_float(self):
        """julday should return a native Python float."""
        jd = swe.julday(2000, 1, 1, 12.0)
        assert type(jd) is float

    @pytest.mark.unit
    def test_revjul_returns_4_values(self):
        """revjul should return (year, month, day, hour)."""
        result = swe.revjul(2451545.0)
        assert len(result) == 4

    @pytest.mark.unit
    def test_julday_consecutive_days(self):
        """Consecutive days should differ by exactly 1.0."""
        jd1 = swe.julday(2024, 3, 15, 12.0)
        jd2 = swe.julday(2024, 3, 16, 12.0)
        assert abs(jd2 - jd1 - 1.0) < 1e-10

    @pytest.mark.unit
    def test_julday_12_hours_is_half_day(self):
        """12 hours should be 0.5 JD."""
        jd1 = swe.julday(2024, 3, 15, 0.0)
        jd2 = swe.julday(2024, 3, 15, 12.0)
        assert abs(jd2 - jd1 - 0.5) < 1e-10


class TestDeltaT:
    """Tests for Delta T (TT - UT)."""

    @pytest.mark.unit
    def test_deltat_returns_float(self):
        """deltat should return a native Python float."""
        dt = swe.deltat(2451545.0)
        assert type(dt) is float

    @pytest.mark.unit
    def test_deltat_finite(self):
        """deltat should return a finite value."""
        dt = swe.deltat(2451545.0)
        assert math.isfinite(dt)

    @pytest.mark.unit
    def test_deltat_j2000_positive(self):
        """Delta T at J2000.0 should be positive (~63.8s)."""
        dt = swe.deltat(2451545.0)
        # dt is in days, ~63.8s = ~0.000739 days
        dt_seconds = dt * 86400
        assert 60 < dt_seconds < 70, (
            f"Delta T at J2000: {dt_seconds:.1f}s (expected ~63.8s)"
        )

    @pytest.mark.unit
    def test_deltat_increases_over_centuries(self):
        """Delta T should generally increase over centuries."""
        jd_2000 = 2451545.0
        jd_2100 = jd_2000 + 36525.0
        dt_2000 = swe.deltat(jd_2000)
        dt_2100 = swe.deltat(jd_2100)
        assert dt_2100 > dt_2000, f"Delta T not increasing: {dt_2000} -> {dt_2100}"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "year,min_s,max_s",
        [
            (1900, -5, 5),  # Delta T near 0 around 1900
            (2000, 60, 70),  # ~63.8s
            (1800, 10, 20),  # ~13.7s
        ],
    )
    def test_deltat_historical_values(self, year: int, min_s: float, max_s: float):
        """Delta T at historical dates should be in expected range."""
        jd = swe.julday(year, 1, 1, 12.0)
        dt = swe.deltat(jd)
        dt_seconds = dt * 86400
        assert min_s < dt_seconds < max_s, (
            f"Year {year}: Delta T = {dt_seconds:.1f}s, expected [{min_s}, {max_s}]"
        )


class TestSidtime:
    """Tests for sidereal time."""

    @pytest.mark.unit
    def test_sidtime_returns_float(self):
        """sidtime should return a native Python float."""
        st = swe.sidtime(2451545.0)
        assert type(st) is float

    @pytest.mark.unit
    def test_sidtime_in_range(self):
        """Sidereal time should be 0-24 hours."""
        st = swe.sidtime(2451545.0)
        assert 0 <= st < 24, f"Sidereal time {st} out of [0,24)"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "jd",
        [2451545.0, 2440587.5, 2460000.0, 2415020.0],
    )
    def test_sidtime_various_dates(self, jd: float):
        """Sidereal time valid at various dates."""
        st = swe.sidtime(jd)
        assert 0 <= st < 24, f"JD {jd}: sidtime={st}"
        assert math.isfinite(st)

    @pytest.mark.unit
    def test_sidtime_advances_with_jd(self):
        """Sidereal time should advance ~4 minutes/day faster than solar."""
        jd1 = 2451545.0
        jd2 = jd1 + 1.0  # Next day

        st1 = swe.sidtime(jd1)
        st2 = swe.sidtime(jd2)

        # Sidereal day is ~23h 56m, so sidtime gains ~3.94 min/day
        # In 1 solar day, sidtime advances by ~24h 3.94m => wraps around
        # The difference modulo 24 should be small (~0.0657 hours)
        diff = (st2 - st1) % 24
        if diff > 12:
            diff = 24 - diff
        # Allow generous range
        assert diff < 1.0 or diff > 23.0, f"Sidtime diff over 1 day: {diff:.4f}h"


class TestUtcToJd:
    """Tests for utc_to_jd conversion."""

    @pytest.mark.unit
    def test_utc_to_jd_returns_tuple(self):
        """utc_to_jd should return a tuple of 2 JD values."""
        result = swe.utc_to_jd(2000, 1, 1, 12, 0, 0.0)
        assert len(result) >= 2

    @pytest.mark.unit
    def test_utc_to_jd_j2000(self):
        """UTC 2000-01-01 12:00:00 should be close to JD 2451545.0."""
        result = swe.utc_to_jd(2000, 1, 1, 12, 0, 0.0)
        # result[0] is JD in ET, result[1] is JD in UT
        jd_ut = result[1]
        assert abs(jd_ut - 2451545.0) < 0.01, f"JD_UT={jd_ut}"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "y,m,d,h,mi,s",
        [
            (2024, 1, 1, 0, 0, 0.0),
            (2024, 6, 21, 12, 30, 45.0),
            (1999, 12, 31, 23, 59, 59.0),
            (1950, 7, 4, 18, 0, 0.0),
        ],
    )
    def test_utc_to_jd_various_dates(
        self, y: int, m: int, d: int, h: int, mi: int, s: float
    ):
        """utc_to_jd works for various dates."""
        result = swe.utc_to_jd(y, m, d, h, mi, s)
        assert len(result) >= 2
        jd_ut = result[1]
        assert math.isfinite(jd_ut)
        assert jd_ut > 0


class TestJdut1ToUtc:
    """Tests for JD UT1 to UTC conversion."""

    @pytest.mark.unit
    def test_jdut1_to_utc_returns_tuple(self):
        """jdut1_to_utc returns a tuple."""
        result = swe.jdut1_to_utc(2451545.0)
        assert len(result) >= 5  # year, month, day, hour, minute, second

    @pytest.mark.unit
    def test_jdut1_to_utc_j2000(self):
        """JD 2451545.0 should convert to ~2000-01-01."""
        result = swe.jdut1_to_utc(2451545.0)
        y, m, d = result[0], result[1], result[2]
        assert y == 2000
        assert m == 1
        assert d == 1


class TestJdetToUtc:
    """Tests for JD ET to UTC conversion."""

    @pytest.mark.unit
    def test_jdet_to_utc_returns_tuple(self):
        """jdet_to_utc returns a tuple."""
        result = swe.jdet_to_utc(2451545.0)
        assert len(result) >= 5

    @pytest.mark.unit
    def test_jdet_to_utc_j2000(self):
        """JD ET 2451545.0 should convert to ~2000-01-01."""
        result = swe.jdet_to_utc(2451545.0)
        y, m = result[0], result[1]
        assert y == 2000
        assert m == 1


class TestUtcRoundTrip:
    """Test UTC conversion round-trips."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "y,m,d,h,mi,s",
        [
            (2000, 1, 1, 12, 0, 0.0),
            (2024, 6, 15, 8, 30, 0.0),
            (1980, 3, 20, 12, 0, 0.0),  # Use noon to avoid midnight boundary issues
        ],
    )
    def test_utc_to_jd_to_utc_roundtrip(
        self, y: int, m: int, d: int, h: int, mi: int, s: float
    ):
        """UTC -> JD_UT -> UTC should roughly recover original date."""
        jd_et, jd_ut = swe.utc_to_jd(y, m, d, h, mi, s)
        result = swe.jdut1_to_utc(jd_ut)
        ry, rm, rd = int(result[0]), int(result[1]), int(result[2])
        assert ry == y, f"Year: {ry} != {y}"
        assert rm == m, f"Month: {rm} != {m}"
        assert rd == d, f"Day: {rd} != {d}"
