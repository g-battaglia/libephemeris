"""
Tests for degnorm, radnorm, and other utility functions.

Verifies angle normalization, coordinate transforms, and edge cases.
"""

from __future__ import annotations

import math

import pytest

from libephemeris.utils import degnorm, radnorm
import libephemeris as swe


TWO_PI = 2 * math.pi


class TestDegnorm:
    """Test degnorm angle normalization."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "input_val,expected",
        [
            (0.0, 0.0),
            (180.0, 180.0),
            (359.999, 359.999),
            (360.0, 0.0),
            (361.0, 1.0),
            (720.0, 0.0),
            (540.0, 180.0),
            (-1.0, 359.0),
            (-90.0, 270.0),
            (-180.0, 180.0),
            (-360.0, 0.0),
            (-361.0, 359.0),
            (-720.0, 0.0),
        ],
    )
    def test_known_values(self, input_val: float, expected: float):
        """degnorm returns correct values for known inputs."""
        result = degnorm(input_val)
        assert abs(result - expected) < 1e-10, (
            f"degnorm({input_val}) = {result}, expected {expected}"
        )

    @pytest.mark.unit
    def test_output_range(self):
        """Output should always be in [0, 360)."""
        import numpy as np

        rng = np.random.default_rng(42)
        angles = rng.uniform(-3600, 3600, 100)
        for angle in angles:
            result = degnorm(float(angle))
            assert 0 <= result < 360, f"degnorm({angle}) = {result}"

    @pytest.mark.unit
    def test_returns_float(self):
        """degnorm returns a float."""
        assert isinstance(degnorm(45.0), float)

    @pytest.mark.unit
    def test_very_large_angle(self):
        """degnorm handles very large angles."""
        result = degnorm(36000.5)
        assert abs(result - 0.5) < 1e-8

    @pytest.mark.unit
    def test_very_negative_angle(self):
        """degnorm handles very negative angles."""
        result = degnorm(-36000.5)
        assert abs(result - 359.5) < 1e-8


class TestRadnorm:
    """Test radnorm angle normalization."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "input_val,expected",
        [
            (0.0, 0.0),
            (math.pi, math.pi),
            (TWO_PI, 0.0),
            (3 * math.pi, math.pi),
            (-math.pi / 4, TWO_PI - math.pi / 4),
            (-math.pi, math.pi),
            (-TWO_PI, 0.0),
        ],
    )
    def test_known_values(self, input_val: float, expected: float):
        """radnorm returns correct values for known inputs."""
        result = radnorm(input_val)
        assert abs(result - expected) < 1e-10, (
            f"radnorm({input_val}) = {result}, expected {expected}"
        )

    @pytest.mark.unit
    def test_output_range(self):
        """Output should always be in [0, 2*pi)."""
        import numpy as np

        rng = np.random.default_rng(42)
        angles = rng.uniform(-100, 100, 100)
        for angle in angles:
            result = radnorm(float(angle))
            assert 0 <= result < TWO_PI, f"radnorm({angle}) = {result}"

    @pytest.mark.unit
    def test_returns_float(self):
        """radnorm returns a float."""
        assert isinstance(radnorm(1.5), float)


class TestWrapAround:
    """Test 0/360 degree wrap-around edge cases in calc_ut."""

    @pytest.mark.unit
    def test_longitude_always_in_range(self):
        """calc_ut longitude should always be [0, 360)."""
        import numpy as np

        rng = np.random.default_rng(42)
        # Test at 50 random dates
        jd_start = 2451545.0 - 365 * 50  # ~1950
        jd_end = 2451545.0 + 365 * 50  # ~2050
        jds = rng.uniform(jd_start, jd_end, 50)

        for jd in jds:
            result, _ = swe.calc_ut(float(jd), 0, 256)  # Sun + SPEED
            lon = result[0]
            assert 0 <= lon < 360, f"JD {jd}: Sun lon = {lon}"

    @pytest.mark.unit
    def test_moon_longitude_always_in_range(self):
        """Moon longitude should always be [0, 360)."""
        import numpy as np

        rng = np.random.default_rng(123)
        jds = rng.uniform(2451545.0, 2451545.0 + 365, 30)

        for jd in jds:
            result, _ = swe.calc_ut(float(jd), 1, 256)  # Moon + SPEED
            lon = result[0]
            assert 0 <= lon < 360, f"JD {jd}: Moon lon = {lon}"


class TestEarthGeocentric:
    """Test Earth (body 14) with various flags."""

    @pytest.mark.unit
    def test_earth_geocentric_all_zeros(self):
        """Earth geocentric should return all zeros."""
        result, _ = swe.calc_ut(2451545.0, 14, 256)
        for i in range(6):
            assert result[i] == 0.0, f"Earth result[{i}] = {result[i]}"

    @pytest.mark.unit
    def test_earth_heliocentric_nonzero(self):
        """Earth heliocentric should return nonzero values."""
        from libephemeris.constants import FLG_HELCTR, FLG_SPEED

        result, _ = swe.calc_ut(2451545.0, 14, FLG_HELCTR | FLG_SPEED)
        # Heliocentric Earth should have nonzero distance
        assert result[2] > 0, f"Earth helio dist = {result[2]}"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "flags_name,flags_val",
        [
            ("SPEED", 256),
            ("SWIEPH", 2),
            ("SWIEPH+SPEED", 258),
            ("J2000", 32768),  # FLG_J2000
            ("ICRS", 131072),  # FLG_ICRS
        ],
    )
    def test_earth_geocentric_various_flags(self, flags_name: str, flags_val: int):
        """Earth geocentric with various flags should return zeros without error."""
        result, _ = swe.calc_ut(2451545.0, 14, flags_val)
        for i in range(6):
            assert result[i] == 0.0, f"Earth ({flags_name}) result[{i}] = {result[i]}"


class TestCotransComprehensive:
    """Extended tests for cotrans."""

    @pytest.mark.unit
    def test_identity_at_zero_obliquity(self):
        """With eps=0, cotrans should return unchanged coordinates."""
        lon, lat, dist = 120.0, 30.0, 1.5
        # Negative eps = ecliptic to equatorial
        result = swe.cotrans((lon, lat, dist), 0.0)
        assert abs(result[0] - lon) < 1e-8
        assert abs(result[1] - lat) < 1e-8
        assert abs(result[2] - dist) < 1e-8

    @pytest.mark.unit
    def test_roundtrip_ecl_equ_ecl(self):
        """Ecliptic -> equatorial -> ecliptic should recover original."""
        lon, lat, dist = 85.0, 5.0, 1.0
        eps = 23.44

        # Ecliptic to equatorial (negative eps)
        equ = swe.cotrans((lon, lat, dist), -eps)
        # Equatorial to ecliptic (positive eps)
        ecl = swe.cotrans((equ[0], equ[1], equ[2]), eps)

        assert abs(ecl[0] - lon) < 1e-6, f"Lon: {ecl[0]} != {lon}"
        assert abs(ecl[1] - lat) < 1e-6, f"Lat: {ecl[1]} != {lat}"
        assert abs(ecl[2] - dist) < 1e-6, f"Dist: {ecl[2]} != {dist}"

    @pytest.mark.unit
    def test_ecliptic_pole_to_equatorial(self):
        """Ecliptic north pole (lat=90) should map to specific equatorial coords."""
        eps = 23.44
        result = swe.cotrans((0.0, 90.0, 1.0), -eps)
        # Ecliptic pole should have equatorial declination = 90 - eps
        assert abs(result[1] - (90.0 - eps)) < 0.01, (
            f"Ecliptic pole dec: {result[1]} (expected ~{90.0 - eps})"
        )
