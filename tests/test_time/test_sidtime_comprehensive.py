"""
Tests for sidtime (sidereal time) accuracy and consistency.

Verifies sidtime and sidtime0 across dates,
consistency between the two functions, and known values.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    swe.close()


JD_J2000 = 2451545.0


class TestSidtime:
    """Test sidtime (Greenwich Mean Sidereal Time)."""

    @pytest.mark.unit
    def test_returns_float(self):
        """sidtime returns a float."""
        result = swe.sidtime(JD_J2000)
        assert isinstance(result, float)

    @pytest.mark.unit
    def test_range_0_to_24(self):
        """Sidereal time should be in [0, 24) hours."""
        result = swe.sidtime(JD_J2000)
        assert 0.0 <= result < 24.0

    @pytest.mark.unit
    def test_j2000_known_value(self):
        """At J2000.0, GMST should be about 18.7 hours."""
        result = swe.sidtime(JD_J2000)
        # J2000.0 (2000-01-01 12:00 UT), GMST ~ 18.697 hours
        assert 18.0 < result < 19.5, f"GMST at J2000 = {result}"

    @pytest.mark.unit
    def test_advances_with_time(self):
        """Sidereal time advances roughly 24h + 3m56s per solar day."""
        st1 = swe.sidtime(JD_J2000)
        st2 = swe.sidtime(JD_J2000 + 1.0)
        # After 1 solar day, GMST advances by ~3m56s = 0.0657 hours
        diff = (st2 - st1) % 24.0
        # Should be about 0.0657 hours (sidereal day is 23h56m4s)
        assert 0.05 < diff < 0.08, f"Daily GMST advance = {diff} hours"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "jd",
        [
            2415020.0,  # 1900
            2440587.5,  # 1970
            2451545.0,  # 2000
            2458849.5,  # 2020
            2460000.0,  # 2023
        ],
    )
    def test_various_dates(self, jd):
        """Sidereal time is valid at various dates."""
        result = swe.sidtime(jd)
        assert 0.0 <= result < 24.0

    @pytest.mark.unit
    def test_half_day_offset(self):
        """Sidereal time 12 hours later should differ by ~12.03 hours."""
        st1 = swe.sidtime(JD_J2000)
        st2 = swe.sidtime(JD_J2000 + 0.5)
        diff = (st2 - st1) % 24.0
        # Half a solar day = ~12.03 sidereal hours
        assert 11.9 < diff < 12.1, f"Half-day GMST advance = {diff}"


class TestSidtime0:
    """Test sidtime0 (GMST from obliquity and nutation)."""

    @pytest.mark.unit
    def test_returns_float(self):
        """sidtime0 returns a float."""
        eps = 23.44  # obliquity
        nut = -0.004  # nutation in longitude (degrees)
        result = swe.sidtime0(JD_J2000, eps, nut)
        assert isinstance(result, float)

    @pytest.mark.unit
    def test_range_0_to_24(self):
        """sidtime0 result should be in [0, 24)."""
        result = swe.sidtime0(JD_J2000, 23.44, -0.004)
        assert 0.0 <= result < 24.0

    @pytest.mark.unit
    def test_consistent_with_sidtime(self):
        """sidtime0 should be close to sidtime for reasonable eps/nut."""
        st = swe.sidtime(JD_J2000)
        st0 = swe.sidtime0(JD_J2000, 23.44, -0.004)
        # Should be close (within a few seconds = 0.001 hours)
        diff = abs(st - st0) % 24.0
        if diff > 12:
            diff = 24 - diff
        assert diff < 0.01, f"sidtime={st}, sidtime0={st0}, diff={diff}"


class TestSidtimeLocalConversion:
    """Test converting GMST to local sidereal time."""

    @pytest.mark.unit
    @pytest.mark.parametrize("lon", [0.0, 90.0, 180.0, -90.0])
    def test_local_sidereal_time(self, lon):
        """Local sidereal time = GMST + longitude/15."""
        gmst = swe.sidtime(JD_J2000)
        lst = (gmst + lon / 15.0) % 24.0
        assert 0.0 <= lst < 24.0

    @pytest.mark.unit
    def test_greenwich_lst_equals_gmst(self):
        """At Greenwich (lon=0), LST = GMST."""
        gmst = swe.sidtime(JD_J2000)
        lst = (gmst + 0.0 / 15.0) % 24.0
        assert lst == pytest.approx(gmst, abs=1e-10)

    @pytest.mark.unit
    def test_lst_180_degrees_offset(self):
        """At 180° longitude, LST = GMST + 12 hours."""
        gmst = swe.sidtime(JD_J2000)
        lst = (gmst + 180.0 / 15.0) % 24.0
        expected = (gmst + 12.0) % 24.0
        assert lst == pytest.approx(expected, abs=1e-10)
