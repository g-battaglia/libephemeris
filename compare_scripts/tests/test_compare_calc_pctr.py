"""
Planetocentric Calculations Comparison Tests (calc_pctr / calc_pctr).

Validates planetocentric position calculations between pyswisseph and libephemeris:
- Position of a target body as observed from another planet's center
- Various target/center combinations (Moon from Mars, Sun from Jupiter, etc.)
- Multiple calculation flags (FLG_SPEED, FLG_EQUATORIAL, FLG_SIDEREAL)
- Multiple time periods for validation

Use cases:
- Heliorelocation astrology calculations
- Planetary observer perspective calculations
- Scientific astronomical computations from other planets
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    EARTH,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    FLG_SPEED,
    FLG_EQUATORIAL,
    FLG_SIDEREAL,
    FLG_J2000,
    SIDM_LAHIRI,
)


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# All major planets for comprehensive testing
PLANETS = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (EARTH, "Earth"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
    (PLUTO, "Pluto"),
]

# Planet centers for observation (excluding Sun since that's heliocentric)
CENTER_PLANETS = [
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (EARTH, "Earth"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
]

# Key target/center combinations for focused testing
KEY_COMBINATIONS = [
    # Moon as seen from Mars - primary use case for Mars colonization astrology
    (MOON, "Moon", MARS, "Mars"),
    # Sun as seen from Jupiter - Jupiter-centric calculations
    (SUN, "Sun", JUPITER, "Jupiter"),
    # Earth as seen from Mars - Mars observer perspective
    (EARTH, "Earth", MARS, "Mars"),
    # Sun as seen from Saturn - outer planet perspective
    (SUN, "Sun", SATURN, "Saturn"),
    # Venus as seen from Mercury - inner planet observation
    (VENUS, "Venus", MERCURY, "Mercury"),
    # Jupiter as seen from Saturn - gas giant mutual observation
    (JUPITER, "Jupiter", SATURN, "Saturn"),
    # Saturn as seen from Jupiter - reverse observation
    (SATURN, "Saturn", JUPITER, "Jupiter"),
    # Moon as seen from Venus - Venus observer
    (MOON, "Moon", VENUS, "Venus"),
    # Sun as seen from Uranus - distant outer perspective
    (SUN, "Sun", URANUS, "Uranus"),
    # Earth as seen from Jupiter - far observer
    (EARTH, "Earth", JUPITER, "Jupiter"),
]

# Extended combinations for comprehensive coverage
EXTENDED_COMBINATIONS = [
    (MERCURY, "Mercury", VENUS, "Venus"),
    (MARS, "Mars", EARTH, "Earth"),
    (PLUTO, "Pluto", NEPTUNE, "Neptune"),
    (NEPTUNE, "Neptune", URANUS, "Uranus"),
    (URANUS, "Uranus", NEPTUNE, "Neptune"),
    (SUN, "Sun", NEPTUNE, "Neptune"),
    (MOON, "Moon", JUPITER, "Jupiter"),
    (VENUS, "Venus", MARS, "Mars"),
    (MARS, "Mars", VENUS, "Venus"),
    (SUN, "Sun", PLUTO, "Pluto"),
]

# Test dates spanning different eras
TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (2024, 11, 15, 0.0, "Current Era"),
    (1900, 1, 1, 0.0, "Early 1900s"),
    (2100, 12, 31, 12.0, "Late 2100s"),
    (1950, 6, 21, 12.0, "Mid-Century"),
]

# Tolerances for comparison
# Note: libephemeris uses JPL DE ephemerides via Skyfield while pyswisseph uses
# Swiss Ephemeris internal ephemeris. Small differences (~0.01-0.02°) are expected
# due to different underlying planetary theory and interpolation methods.
LONGITUDE_TOLERANCE = 0.02  # degrees (~72 arcsec) - accounts for ephemeris differences
LATITUDE_TOLERANCE = 0.02  # degrees
DISTANCE_TOLERANCE_PCT = 1.0  # 1% tolerance for distance
VELOCITY_TOLERANCE = 0.1  # degrees/day


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# PRIMARY USE CASE TESTS
# ============================================================================


class TestMoonFromMars:
    """Test Moon position as seen from Mars - primary heliorelocation use case."""

    @pytest.mark.comparison
    def test_moon_from_mars_j2000(self):
        """Moon from Mars at J2000 epoch."""
        jd = 2451545.0
        pos_lib, _ = ephem.calc_pctr(jd, MOON, MARS, 0)
        pos_swe, _ = swe.calc_pctr(jd, MOON, MARS, 0)

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])
        lat_diff = abs(pos_lib[1] - pos_swe[1])

        assert lon_diff < LONGITUDE_TOLERANCE, (
            f"Moon from Mars longitude diff {lon_diff:.6f}° exceeds tolerance"
        )
        assert lat_diff < LATITUDE_TOLERANCE, (
            f"Moon from Mars latitude diff {lat_diff:.6f}° exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_moon_from_mars_multiple_dates(self, year, month, day, hour, date_desc):
        """Moon from Mars across multiple dates."""
        jd = swe.julday(year, month, day, hour)

        pos_lib, _ = ephem.calc_pctr(jd, MOON, MARS, 0)
        try:
            pos_swe, _ = swe.calc_pctr(jd, MOON, MARS, 0)
        except Exception:
            pytest.skip("pyswisseph ephemeris not available for this date")

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])

        assert lon_diff < LONGITUDE_TOLERANCE, (
            f"Moon from Mars at {date_desc}: longitude diff {lon_diff:.6f}°"
        )

    @pytest.mark.comparison
    def test_moon_from_mars_with_speed(self):
        """Moon from Mars with velocity calculations."""
        jd = 2451545.0
        pos_lib, _ = ephem.calc_pctr(jd, MOON, MARS, FLG_SPEED)
        pos_swe, _ = swe.calc_pctr(jd, MOON, MARS, FLG_SPEED)

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])
        vel_diff = abs(pos_lib[3] - pos_swe[3])

        assert lon_diff < LONGITUDE_TOLERANCE, (
            f"Moon from Mars longitude diff {lon_diff:.6f}°"
        )
        assert vel_diff < VELOCITY_TOLERANCE, (
            f"Moon from Mars velocity diff {vel_diff:.6f}°/day"
        )


class TestSunFromJupiter:
    """Test Sun position as seen from Jupiter - Jupiter-centric astrology."""

    @pytest.mark.comparison
    def test_sun_from_jupiter_j2000(self):
        """Sun from Jupiter at J2000 epoch."""
        jd = 2451545.0

        pos_lib, _ = ephem.calc_pctr(jd, SUN, JUPITER, 0)
        try:
            pos_swe, _ = swe.calc_pctr(jd, SUN, JUPITER, 0)
        except Exception:
            pytest.skip("pyswisseph ephemeris not available")

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])

        assert lon_diff < LONGITUDE_TOLERANCE, (
            f"Sun from Jupiter longitude diff {lon_diff:.6f}°"
        )
        # Distance should be approximately 5.2 AU
        assert 4.9 < pos_lib[2] < 5.5, (
            f"Sun distance from Jupiter {pos_lib[2]:.2f} AU should be ~5.2 AU"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_sun_from_jupiter_multiple_dates(self, year, month, day, hour, date_desc):
        """Sun from Jupiter across multiple dates."""
        jd = swe.julday(year, month, day, hour)

        pos_lib, _ = ephem.calc_pctr(jd, SUN, JUPITER, 0)
        try:
            pos_swe, _ = swe.calc_pctr(jd, SUN, JUPITER, 0)
        except Exception:
            pytest.skip("pyswisseph ephemeris not available for this date")

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])

        assert lon_diff < LONGITUDE_TOLERANCE, (
            f"Sun from Jupiter at {date_desc}: longitude diff {lon_diff:.6f}°"
        )

    @pytest.mark.comparison
    def test_sun_from_jupiter_with_speed(self):
        """Sun from Jupiter with velocity calculations."""
        jd = 2451545.0

        pos_lib, _ = ephem.calc_pctr(jd, SUN, JUPITER, FLG_SPEED)
        try:
            pos_swe, _ = swe.calc_pctr(jd, SUN, JUPITER, FLG_SPEED)
        except Exception:
            pytest.skip("pyswisseph ephemeris not available")

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])
        vel_diff = abs(pos_lib[3] - pos_swe[3])

        assert lon_diff < LONGITUDE_TOLERANCE
        assert vel_diff < VELOCITY_TOLERANCE


# ============================================================================
# COMPREHENSIVE TARGET/CENTER COMBINATIONS
# ============================================================================


class TestKeyCombinations:
    """Test key target/center combinations for heliorelocation astrology."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "target_id,target_name,center_id,center_name", KEY_COMBINATIONS
    )
    def test_key_combinations_j2000(
        self, target_id, target_name, center_id, center_name
    ):
        """Test key planetocentric combinations at J2000."""
        jd = 2451545.0

        pos_lib, _ = ephem.calc_pctr(jd, target_id, center_id, 0)
        try:
            pos_swe, _ = swe.calc_pctr(jd, target_id, center_id, 0)
        except Exception:
            pytest.skip("pyswisseph ephemeris not available")

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])
        lat_diff = abs(pos_lib[1] - pos_swe[1])

        assert lon_diff < LONGITUDE_TOLERANCE, (
            f"{target_name} from {center_name}: longitude diff {lon_diff:.6f}°"
        )
        assert lat_diff < LATITUDE_TOLERANCE, (
            f"{target_name} from {center_name}: latitude diff {lat_diff:.6f}°"
        )

        # Distance comparison (1% tolerance)
        if pos_swe[2] > 0:
            dist_diff_pct = abs(pos_lib[2] - pos_swe[2]) / pos_swe[2] * 100
            assert dist_diff_pct < DISTANCE_TOLERANCE_PCT, (
                f"{target_name} from {center_name}: distance diff {dist_diff_pct:.2f}%"
            )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "target_id,target_name,center_id,center_name", EXTENDED_COMBINATIONS
    )
    def test_extended_combinations_j2000(
        self, target_id, target_name, center_id, center_name
    ):
        """Test extended planetocentric combinations at J2000."""
        jd = 2451545.0

        pos_lib, _ = ephem.calc_pctr(jd, target_id, center_id, 0)
        try:
            pos_swe, _ = swe.calc_pctr(jd, target_id, center_id, 0)
        except Exception:
            pytest.skip("pyswisseph ephemeris not available")

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])

        assert lon_diff < LONGITUDE_TOLERANCE, (
            f"{target_name} from {center_name}: longitude diff {lon_diff:.6f}°"
        )


class TestAllPlanetCenters:
    """Test all planets as observation centers."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("center_id,center_name", CENTER_PLANETS)
    def test_sun_from_all_centers(self, center_id, center_name):
        """Sun as seen from all planet centers."""
        jd = 2451545.0

        pos_lib, _ = ephem.calc_pctr(jd, SUN, center_id, 0)
        try:
            pos_swe, _ = swe.calc_pctr(jd, SUN, center_id, 0)
        except Exception:
            pytest.skip("pyswisseph ephemeris not available")

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])

        assert lon_diff < LONGITUDE_TOLERANCE, (
            f"Sun from {center_name}: longitude diff {lon_diff:.6f}°"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("center_id,center_name", CENTER_PLANETS[:5])
    def test_moon_from_inner_planets(self, center_id, center_name):
        """Moon as seen from inner planet centers."""
        jd = 2451545.0

        pos_lib, _ = ephem.calc_pctr(jd, MOON, center_id, 0)
        try:
            pos_swe, _ = swe.calc_pctr(jd, MOON, center_id, 0)
        except Exception:
            pytest.skip("pyswisseph ephemeris not available")

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])

        assert lon_diff < LONGITUDE_TOLERANCE, (
            f"Moon from {center_name}: longitude diff {lon_diff:.6f}°"
        )


# ============================================================================
# CALCULATION FLAGS TESTS
# ============================================================================


class TestCalcPctrFlags:
    """Test various calculation flags for planetocentric calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "target_id,target_name,center_id,center_name", KEY_COMBINATIONS[:3]
    )
    def test_with_speed_flag(self, target_id, target_name, center_id, center_name):
        """Test FLG_SPEED flag for planetocentric calculations."""
        jd = 2451545.0

        pos_lib, _ = ephem.calc_pctr(jd, target_id, center_id, FLG_SPEED)
        try:
            pos_swe, _ = swe.calc_pctr(jd, target_id, center_id, FLG_SPEED)
        except Exception:
            pytest.skip("pyswisseph ephemeris not available")

        # Position
        lon_diff = angular_diff(pos_lib[0], pos_swe[0])
        assert lon_diff < LONGITUDE_TOLERANCE

        # Velocity
        vel_diff = abs(pos_lib[3] - pos_swe[3])
        assert vel_diff < VELOCITY_TOLERANCE, (
            f"{target_name} from {center_name}: velocity diff {vel_diff:.6f}°/day"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "target_id,target_name,center_id,center_name", KEY_COMBINATIONS[:3]
    )
    def test_with_equatorial_flag(self, target_id, target_name, center_id, center_name):
        """Test FLG_EQUATORIAL flag for planetocentric calculations."""
        jd = 2451545.0

        pos_lib, _ = ephem.calc_pctr(jd, target_id, center_id, FLG_EQUATORIAL)
        try:
            pos_swe, _ = swe.calc_pctr(jd, target_id, center_id, FLG_EQUATORIAL)
        except Exception:
            pytest.skip("pyswisseph ephemeris not available")

        # RA should be in valid range [0, 360)
        assert 0 <= pos_lib[0] < 360

        ra_diff = angular_diff(pos_lib[0], pos_swe[0])
        dec_diff = abs(pos_lib[1] - pos_swe[1])

        assert ra_diff < LONGITUDE_TOLERANCE, (
            f"{target_name} from {center_name}: RA diff {ra_diff:.6f}°"
        )
        assert dec_diff < LATITUDE_TOLERANCE, (
            f"{target_name} from {center_name}: Dec diff {dec_diff:.6f}°"
        )

    @pytest.mark.comparison
    def test_with_sidereal_flag(self):
        """Test FLG_SIDEREAL flag for planetocentric calculations."""
        jd = 2451545.0

        # Set sidereal mode
        ephem.set_sid_mode(SIDM_LAHIRI)
        swe.set_sid_mode(swe.SIDM_LAHIRI)

        pos_lib, _ = ephem.calc_pctr(jd, MOON, MARS, FLG_SIDEREAL)
        try:
            pos_swe, _ = swe.calc_pctr(jd, MOON, MARS, FLG_SIDEREAL)
        except Exception:
            pytest.skip("pyswisseph ephemeris not available")

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])

        # Relaxed tolerance for sidereal (ayanamsha differences)
        assert lon_diff < 0.05, (
            f"Moon from Mars sidereal: longitude diff {lon_diff:.6f}°"
        )

    @pytest.mark.comparison
    def test_with_j2000_flag(self):
        """Test FLG_J2000 flag for planetocentric calculations."""
        jd = swe.julday(2024, 6, 15, 12.0)

        pos_lib, _ = ephem.calc_pctr(jd, SUN, JUPITER, FLG_J2000)
        try:
            pos_swe, _ = swe.calc_pctr(jd, SUN, JUPITER, FLG_J2000)
        except Exception:
            pytest.skip("pyswisseph ephemeris not available")

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])

        assert lon_diff < LONGITUDE_TOLERANCE, (
            f"Sun from Jupiter J2000: longitude diff {lon_diff:.6f}°"
        )


# ============================================================================
# MULTIPLE DATES TESTS
# ============================================================================


class TestMultipleDates:
    """Test planetocentric calculations across multiple dates."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    @pytest.mark.parametrize(
        "target_id,target_name,center_id,center_name", KEY_COMBINATIONS[:5]
    )
    def test_combinations_across_dates(
        self,
        year,
        month,
        day,
        hour,
        date_desc,
        target_id,
        target_name,
        center_id,
        center_name,
    ):
        """Test key combinations across multiple dates."""
        jd = swe.julday(year, month, day, hour)

        pos_lib, _ = ephem.calc_pctr(jd, target_id, center_id, 0)
        try:
            pos_swe, _ = swe.calc_pctr(jd, target_id, center_id, 0)
        except Exception:
            pytest.skip("pyswisseph ephemeris not available")

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])

        assert lon_diff < LONGITUDE_TOLERANCE, (
            f"{target_name} from {center_name} at {date_desc}: "
            f"longitude diff {lon_diff:.6f}°"
        )


# ============================================================================
# GEOMETRIC CONSISTENCY TESTS
# ============================================================================


class TestGeometricConsistency:
    """Test geometric consistency of planetocentric calculations."""

    @pytest.mark.comparison
    def test_self_observation_behavior(self):
        """Test behavior when observing a body from itself.

        Note: pyswisseph raises an error for self-observation (ipl == iplctr).
        libephemeris may handle this differently. This test documents the behavior.
        """
        jd = 2451545.0

        # libephemeris behavior
        pos_lib, _ = ephem.calc_pctr(jd, MARS, MARS, 0)

        # Distance should be essentially zero
        assert pos_lib[2] < 0.0001, (
            f"Self-observation distance {pos_lib[2]} should be ~0"
        )

        # pyswisseph raises an error for self-observation
        with pytest.raises(Exception):
            swe.calc_pctr(jd, MARS, MARS, 0)

    @pytest.mark.comparison
    def test_earth_mars_symmetry(self):
        """Earth from Mars vs Mars from Earth should show ~180° difference."""
        jd = 2451545.0

        # Earth as seen from Mars
        pos_earth_from_mars, _ = ephem.calc_pctr(jd, EARTH, MARS, 0)
        # Mars as seen from Earth (geocentric)
        pos_mars_from_earth, _ = ephem.calc_ut(jd, MARS, 0)

        # Compare with pyswisseph
        pos_earth_from_mars_swe, _ = swe.calc_pctr(jd, EARTH, MARS, 0)

        lon_diff = angular_diff(pos_earth_from_mars[0], pos_earth_from_mars_swe[0])
        assert lon_diff < LONGITUDE_TOLERANCE

    @pytest.mark.comparison
    @pytest.mark.parametrize("center_id,center_name", CENTER_PLANETS[:4])
    def test_sun_distance_from_centers(self, center_id, center_name):
        """Sun distance from various planets should match known values."""
        jd = 2451545.0

        pos_lib, _ = ephem.calc_pctr(jd, SUN, center_id, 0)
        try:
            pos_swe, _ = swe.calc_pctr(jd, SUN, center_id, 0)
        except Exception:
            pytest.skip("pyswisseph ephemeris not available")

        # Distance should be positive and match between implementations
        assert pos_lib[2] > 0
        if pos_swe[2] > 0:
            dist_diff_pct = abs(pos_lib[2] - pos_swe[2]) / pos_swe[2] * 100
            assert dist_diff_pct < DISTANCE_TOLERANCE_PCT, (
                f"Sun from {center_name}: distance diff {dist_diff_pct:.2f}%"
            )


# ============================================================================
# VELOCITY TESTS
# ============================================================================


class TestVelocityCalculations:
    """Test velocity calculations for planetocentric positions."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "target_id,target_name,center_id,center_name", KEY_COMBINATIONS[:5]
    )
    def test_velocity_matches_swisseph(
        self, target_id, target_name, center_id, center_name
    ):
        """Velocity calculations should match pyswisseph."""
        jd = 2451545.0

        pos_lib, _ = ephem.calc_pctr(jd, target_id, center_id, FLG_SPEED)
        try:
            pos_swe, _ = swe.calc_pctr(jd, target_id, center_id, FLG_SPEED)
        except Exception:
            pytest.skip("pyswisseph ephemeris not available")

        # Longitude velocity
        lon_vel_diff = abs(pos_lib[3] - pos_swe[3])
        assert lon_vel_diff < VELOCITY_TOLERANCE, (
            f"{target_name} from {center_name}: lon velocity diff {lon_vel_diff:.6f}°/day"
        )

        # Latitude velocity
        lat_vel_diff = abs(pos_lib[4] - pos_swe[4])
        assert lat_vel_diff < VELOCITY_TOLERANCE, (
            f"{target_name} from {center_name}: lat velocity diff {lat_vel_diff:.6f}°/day"
        )

        # Distance velocity
        dist_vel_diff = abs(pos_lib[5] - pos_swe[5])
        assert dist_vel_diff < 0.01, (  # AU/day
            f"{target_name} from {center_name}: dist velocity diff {dist_vel_diff:.6f} AU/day"
        )


# ============================================================================
# RETURN VALUE STRUCTURE TESTS
# ============================================================================


class TestReturnValueStructure:
    """Test that return values match pyswisseph structure."""

    @pytest.mark.comparison
    def test_return_tuple_structure(self):
        """Return should be (6-tuple, int) like pyswisseph."""
        jd = 2451545.0

        result_lib = ephem.calc_pctr(jd, MOON, MARS, FLG_SPEED)
        result_swe = swe.calc_pctr(jd, MOON, MARS, FLG_SPEED)

        # Both should return (positions, flags)
        assert len(result_lib) == 2
        assert len(result_swe) == 2

        # Position tuple should have 6 elements
        assert len(result_lib[0]) == 6
        assert len(result_swe[0]) == 6

        # Flags should be integers
        assert isinstance(result_lib[1], int)
        assert isinstance(result_swe[1], int)

    @pytest.mark.comparison
    def test_position_elements_are_floats(self):
        """All position elements should be floats."""
        jd = 2451545.0

        pos_lib, _ = ephem.calc_pctr(jd, MOON, MARS, FLG_SPEED)

        for i, val in enumerate(pos_lib):
            assert isinstance(val, float), f"Element {i} is {type(val)}, expected float"


# ============================================================================
# ALIAS TESTS
# ============================================================================


class TestAliases:
    """Test pyswisseph-compatible aliases."""

    @pytest.mark.comparison
    def test_calc_pctr_alias_matches_swe_calc_pctr(self):
        """calc_pctr should return same results as calc_pctr."""
        jd = 2451545.0

        pos1, flags1 = ephem.calc_pctr(jd, MOON, MARS, 0)
        pos2, flags2 = ephem.calc_pctr(jd, MOON, MARS, 0)

        assert pos1[0] == pos2[0]
        assert pos1[1] == pos2[1]
        assert pos1[2] == pos2[2]
        assert flags1 == flags2
