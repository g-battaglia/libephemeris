"""
Tests for UTC/JD conversion round-trips and degnorm/radnorm edge cases.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.utils import degnorm, radnorm


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    swe.close()


JD_J2000 = 2451545.0


# ============================================================================
# UTC/JD round-trips
# ============================================================================


class TestUtcToJdRoundTrip:
    """Test utc_to_jd -> jdut1_to_utc round-trip."""

    @pytest.mark.unit
    def test_utc_to_jd_returns_two_values(self):
        """utc_to_jd returns (jd_et, jd_ut1) tuple."""
        result = swe.utc_to_jd(2000, 1, 1, 12, 0, 0.0)
        assert len(result) == 2
        jd_et, jd_ut1 = result
        assert isinstance(jd_et, float)
        assert isinstance(jd_ut1, float)

    @pytest.mark.unit
    def test_j2000_value(self):
        """J2000.0 = 2000-01-01 12:00 UT gives known JD."""
        jd_et, jd_ut1 = swe.utc_to_jd(2000, 1, 1, 12, 0, 0.0)
        # JD UT1 should be very close to 2451545.0
        assert jd_ut1 == pytest.approx(JD_J2000, abs=0.01)

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "y,m,d,h,mi,s",
        [
            (2000, 1, 1, 12, 0, 0.0),
            (2020, 6, 15, 6, 30, 0.0),
            (1990, 12, 31, 23, 59, 59.0),
            (2025, 3, 21, 0, 0, 0.0),
        ],
    )
    def test_utc_jd_roundtrip(self, y, m, d, h, mi, s):
        """UTC -> JD -> UTC round-trip."""
        jd_et, jd_ut1 = swe.utc_to_jd(y, m, d, h, mi, s)
        # Convert back
        y2, m2, d2, h2, mi2, s2 = swe.jdut1_to_utc(jd_ut1)
        assert y2 == y
        assert m2 == m
        assert d2 == d
        # Compare total seconds within the day to avoid minute/second boundary issues
        total_secs_orig = h * 3600 + mi * 60 + s
        total_secs_rt = h2 * 3600 + mi2 * 60 + s2
        assert total_secs_rt == pytest.approx(total_secs_orig, abs=1.0)

    @pytest.mark.unit
    def test_jdet_to_utc(self):
        """jdet_to_utc converts ET JD back to UTC components."""
        jd_et, jd_ut1 = swe.utc_to_jd(2020, 6, 15, 12, 0, 0.0)
        y, m, d, h, mi, s = swe.jdet_to_utc(jd_et)
        assert y == 2020
        assert m == 6
        assert d == 15
        assert h == 12
        assert mi == 0
        assert s == pytest.approx(0.0, abs=2.0)

    @pytest.mark.unit
    def test_et_ut_differ_by_deltat(self):
        """JD_ET and JD_UT1 should differ by approximately delta-T."""
        jd_et, jd_ut1 = swe.utc_to_jd(2000, 1, 1, 12, 0, 0.0)
        dt = swe.deltat(jd_ut1)
        # ET = UT1 + delta-T (in days)
        assert jd_et == pytest.approx(jd_ut1 + dt, abs=1e-6)


# ============================================================================
# julday / revjul round-trips
# ============================================================================


class TestJuldayRevjul:
    """Test julday <-> revjul round-trips."""

    @pytest.mark.unit
    def test_julday_j2000(self):
        """julday gives correct JD for J2000."""
        jd = swe.julday(2000, 1, 1, 12.0)
        assert jd == pytest.approx(JD_J2000, abs=0.001)

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "y,m,d,h",
        [
            (2000, 1, 1, 12.0),
            (1900, 1, 1, 0.0),
            (2020, 7, 4, 18.5),
            (1582, 10, 15, 0.0),  # Gregorian calendar start
        ],
    )
    def test_julday_revjul_roundtrip(self, y, m, d, h):
        """julday -> revjul round-trip."""
        jd = swe.julday(y, m, d, h)
        y2, m2, d2, h2 = swe.revjul(jd)
        assert y2 == y
        assert m2 == m
        assert d2 == d
        assert h2 == pytest.approx(h, abs=1e-6)

    @pytest.mark.unit
    def test_revjul_j2000(self):
        """revjul of J2000 gives 2000-01-01 12:00."""
        y, m, d, h = swe.revjul(JD_J2000)
        assert y == 2000
        assert m == 1
        assert d == 1
        assert h == pytest.approx(12.0, abs=1e-6)


# ============================================================================
# degnorm / radnorm edge cases
# ============================================================================


class TestDegnorm:
    """Test degnorm edge cases."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "input_val,expected",
        [
            (0.0, 0.0),
            (360.0, 0.0),
            (720.0, 0.0),
            (180.0, 180.0),
            (359.999, 359.999),
            (-1.0, 359.0),
            (-180.0, 180.0),
            (-360.0, 0.0),
            (-720.0, 0.0),
            (450.0, 90.0),
            (1080.0, 0.0),
        ],
    )
    def test_degnorm_values(self, input_val, expected):
        """degnorm normalizes to [0, 360)."""
        result = degnorm(input_val)
        assert result == pytest.approx(expected, abs=1e-10)

    @pytest.mark.unit
    def test_degnorm_range(self):
        """degnorm always returns [0, 360)."""
        import random

        rng = random.Random(42)
        for _ in range(100):
            val = rng.uniform(-10000, 10000)
            result = degnorm(val)
            assert 0.0 <= result < 360.0, f"degnorm({val}) = {result}"

    @pytest.mark.unit
    def test_degnorm_large_values(self):
        """degnorm handles very large values."""
        assert 0.0 <= degnorm(1e10) < 360.0
        assert 0.0 <= degnorm(-1e10) < 360.0


class TestRadnorm:
    """Test radnorm edge cases."""

    TWO_PI = 2 * math.pi

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "input_val,expected",
        [
            (0.0, 0.0),
            (math.pi, math.pi),
            (2 * math.pi, 0.0),
            (4 * math.pi, 0.0),
            (-math.pi, math.pi),
            (-2 * math.pi, 0.0),
            (3 * math.pi, math.pi),
        ],
    )
    def test_radnorm_values(self, input_val, expected):
        """radnorm normalizes to [0, 2*pi)."""
        result = radnorm(input_val)
        assert result == pytest.approx(expected, abs=1e-10)

    @pytest.mark.unit
    def test_radnorm_range(self):
        """radnorm always returns [0, 2*pi)."""
        import random

        rng = random.Random(42)
        for _ in range(100):
            val = rng.uniform(-100, 100)
            result = radnorm(val)
            assert 0.0 <= result < self.TWO_PI, f"radnorm({val}) = {result}"


# ============================================================================
# cotrans round-trips
# ============================================================================


class TestCotransRoundTrip:
    """Test cotrans ecliptic <-> equatorial round-trips."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "lon,lat",
        [
            (0.0, 0.0),
            (90.0, 0.0),
            (180.0, 0.0),
            (270.0, 0.0),
            (45.0, 30.0),
            (120.0, -15.0),
            (200.0, 60.0),
            (350.0, -80.0),
        ],
    )
    def test_ecliptic_equatorial_roundtrip(self, lon, lat):
        """Ecliptic -> equatorial -> ecliptic round-trip."""
        eps = 23.44  # obliquity
        dist = 1.0

        # Ecliptic to equatorial: negative eps
        ra, dec, d1 = swe.cotrans((lon, lat, dist), -eps)
        # Equatorial to ecliptic: positive eps
        lon2, lat2, d2 = swe.cotrans((ra, dec, d1), eps)

        assert lon2 == pytest.approx(lon, abs=1e-8), f"Longitude: {lon} -> {lon2}"
        assert lat2 == pytest.approx(lat, abs=1e-8), f"Latitude: {lat} -> {lat2}"

    @pytest.mark.unit
    def test_ecliptic_pole(self):
        """Ecliptic pole (lat=90) -> equatorial."""
        eps = 23.44
        ra, dec, _ = swe.cotrans((0.0, 90.0, 1.0), -eps)
        # Ecliptic north pole -> RA=270, Dec=90-eps=66.56
        assert dec == pytest.approx(90.0 - eps, abs=0.01)

    @pytest.mark.unit
    def test_equinox_point(self):
        """Vernal equinox (lon=0, lat=0) -> (RA=0, Dec=0)."""
        ra, dec, _ = swe.cotrans((0.0, 0.0, 1.0), -23.44)
        assert ra == pytest.approx(0.0, abs=0.01)
        assert dec == pytest.approx(0.0, abs=0.01)

    @pytest.mark.unit
    def test_summer_solstice(self):
        """Summer solstice (lon=90, lat=0) -> Dec=+eps."""
        eps = 23.44
        ra, dec, _ = swe.cotrans((90.0, 0.0, 1.0), -eps)
        assert dec == pytest.approx(eps, abs=0.1)
        assert ra == pytest.approx(90.0, abs=1.0)

    @pytest.mark.unit
    def test_cotrans_sp(self):
        """cotrans_sp also works (with speed components)."""
        result = swe.cotrans_sp((90.0, 0.0, 1.0, 1.0, 0.0, 0.0), -23.44)
        assert len(result) == 6


# ============================================================================
# solcross_ut / mooncross_ut precision
# ============================================================================


class TestCrossingPrecision:
    """Test solcross_ut and mooncross_ut precision."""

    @pytest.mark.unit
    def test_solcross_0_degrees(self):
        """Sun crossing 0° (vernal equinox) near J2000."""
        # Sun crosses 0° around March 20
        jd_start = swe.julday(2000, 3, 1, 0.0)
        jd_cross = swe.solcross_ut(0.0, jd_start, 0)
        # Verify Sun is near 0° at crossing
        pos, _ = swe.calc_ut(jd_cross, 0, 256)  # SUN=0, FLG_SPEED=256
        lon = pos[0]
        diff = min(lon, 360 - lon)
        assert diff < 0.001, f"Sun at crossing: {lon}°"

    @pytest.mark.unit
    def test_solcross_90_degrees(self):
        """Sun crossing 90° (summer solstice)."""
        jd_start = swe.julday(2000, 6, 1, 0.0)
        jd_cross = swe.solcross_ut(90.0, jd_start, 0)
        pos, _ = swe.calc_ut(jd_cross, 0, 256)
        diff = abs(pos[0] - 90.0)
        assert diff < 0.001, f"Sun at crossing: {pos[0]}°"

    @pytest.mark.unit
    def test_mooncross_0_degrees(self):
        """Moon crossing 0° near J2000."""
        jd_start = JD_J2000
        jd_cross = swe.mooncross_ut(0.0, jd_start, 0)
        pos, _ = swe.calc_ut(jd_cross, 1, 256)  # MOON=1
        lon = pos[0]
        diff = min(lon, 360 - lon)
        assert diff < 0.01, f"Moon at crossing: {lon}°"

    @pytest.mark.unit
    def test_solcross_returns_future(self):
        """solcross_ut returns a JD after the start date."""
        jd_start = JD_J2000
        jd_cross = swe.solcross_ut(180.0, jd_start, 0)
        assert jd_cross > jd_start

    @pytest.mark.unit
    def test_mooncross_returns_future(self):
        """mooncross_ut returns a JD after the start date."""
        jd_start = JD_J2000
        jd_cross = swe.mooncross_ut(180.0, jd_start, 0)
        assert jd_cross > jd_start
        # Moon crosses any degree within ~27 days
        assert jd_cross < jd_start + 30
