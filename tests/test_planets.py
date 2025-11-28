"""
Unit tests for planetary position calculations.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


@pytest.mark.unit
class TestGeocentricPositions:
    """Tests for geocentric planetary positions."""

    def test_sun_position_j2000(self, standard_jd, compare_with_swisseph):
        """Test Sun position at J2000.0."""
        passed, diffs = compare_with_swisseph(standard_jd, SE_SUN, SEFLG_SWIEPH)
        assert passed, f"Sun position diff: {diffs}"

    def test_all_planets_geocentric(
        self, standard_jd, all_planets, compare_with_swisseph
    ):
        """Test all major planets in geocentric mode."""
        for planet_id, planet_name in all_planets:
            passed, diffs = compare_with_swisseph(standard_jd, planet_id, SEFLG_SWIEPH)
            assert passed, f"{planet_name} position diff: {diffs}"

    def test_moon_high_precision(self, standard_jd):
        """Test Moon position with high precision."""
        res_swe, _ = swe.calc_ut(standard_jd, SE_MOON, SEFLG_SWIEPH)
        res_py, _ = ephem.swe_calc_ut(standard_jd, SE_MOON, SEFLG_SWIEPH)

        lon_diff = abs(res_swe[0] - res_py[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        # Moon should be accurate to ~0.001 degrees
        assert lon_diff < 0.001


@pytest.mark.unit
class TestTopocentric:
    """Tests for topocentric (observer-location-specific) positions."""

    def test_topocentric_moon(self, standard_jd):
        """Test topocentric Moon position."""
        lat, lon, alt = 41.9028, 12.4964, 0  # Rome

        # Set topocentric location
        swe.set_topo(lon, lat, alt)
        ephem.swe_set_topo(lon, lat, alt)

        # Calculate topocentric position
        res_swe, _ = swe.calc_ut(standard_jd, SE_MOON, SEFLG_SWIEPH | SEFLG_TOPOCTR)
        res_py, _ = ephem.swe_calc_ut(
            standard_jd, SE_MOON, SEFLG_SWIEPH | SEFLG_TOPOCTR
        )

        lon_diff = abs(res_swe[0] - res_py[0])
        assert lon_diff < 0.001


@pytest.mark.unit
class TestSiderealMode:
    """Tests for sidereal zodiac calculations."""

    def test_sidereal_sun(self, standard_jd):
        """Test sidereal Sun position."""
        # Set Lahiri ayanamsha
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        res_swe, _ = swe.calc_ut(standard_jd, SE_SUN, SEFLG_SWIEPH | SEFLG_SIDEREAL)
        res_py, _ = ephem.swe_calc_ut(
            standard_jd, SE_SUN, SEFLG_SWIEPH | SEFLG_SIDEREAL
        )

        lon_diff = abs(res_swe[0] - res_py[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        assert lon_diff < 0.001

    def test_sidereal_velocity(self, standard_jd):
        """Test sidereal velocity calculation."""
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        res_swe, _ = swe.calc_ut(
            standard_jd, SE_MARS, SEFLG_SWIEPH | SEFLG_SIDEREAL | SEFLG_SPEED
        )
        res_py, _ = ephem.swe_calc_ut(
            standard_jd, SE_MARS, SEFLG_SWIEPH | SEFLG_SIDEREAL | SEFLG_SPEED
        )

        # Velocity should match within 0.01 deg/day
        speed_diff = abs(res_swe[3] - res_py[3])
        assert speed_diff < 0.01, f"Sidereal velocity diff: {speed_diff}"


@pytest.mark.unit
class TestVelocityCalculations:
    """Tests for velocity (speed) calculations."""

    def test_sun_velocity(self, standard_jd):
        """Test Sun velocity calculation."""
        res_swe, _ = swe.calc_ut(standard_jd, SE_SUN, SEFLG_SWIEPH | SEFLG_SPEED)
        res_py, _ = ephem.swe_calc_ut(standard_jd, SE_SUN, SEFLG_SWIEPH | SEFLG_SPEED)

        # Sun longitude speed ~1 deg/day
        assert 0.9 < res_py[3] < 1.1
        assert abs(res_swe[3] - res_py[3]) < 0.01

    def test_all_planets_velocity(self, standard_jd, all_planets):
        """Test velocity for all planets."""
        for planet_id, planet_name in all_planets:
            res_swe, _ = swe.calc_ut(standard_jd, planet_id, SEFLG_SWIEPH | SEFLG_SPEED)
            res_py, _ = ephem.swe_calc_ut(
                standard_jd, planet_id, SEFLG_SWIEPH | SEFLG_SPEED
            )

            speed_diff = abs(res_swe[3] - res_py[3])
            assert speed_diff < 0.01, f"{planet_name} velocity diff: {speed_diff}"


@pytest.mark.integration
class TestCoordinateFrames:
    """Tests for different coordinate frames."""

    def test_j2000_coordinates(self, standard_jd):
        """Test J2000 coordinate frame."""
        res_swe, _ = swe.calc_ut(standard_jd, SE_SUN, SEFLG_SWIEPH | SEFLG_J2000)
        res_py, _ = ephem.swe_calc_ut(standard_jd, SE_SUN, SEFLG_SWIEPH | SEFLG_J2000)

        lon_diff = abs(res_swe[0] - res_py[0])
        assert lon_diff < 0.001

    def test_equatorial_coordinates(self, standard_jd):
        """Test equatorial (RA/Dec) coordinates."""
        res_swe, _ = swe.calc_ut(standard_jd, SE_SUN, SEFLG_SWIEPH | SEFLG_EQUATORIAL)
        res_py, _ = ephem.swe_calc_ut(
            standard_jd, SE_SUN, SEFLG_SWIEPH | SEFLG_EQUATORIAL
        )

        ra_diff = abs(res_swe[0] - res_py[0])
        if ra_diff > 180:
            ra_diff = 360 - ra_diff

        assert ra_diff < 0.001
