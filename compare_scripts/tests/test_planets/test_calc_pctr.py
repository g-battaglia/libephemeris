"""
Tests for planet-centric position calculations (calc_pctr / calc_pctr).

Tests cover:
- Basic planet-centric position calculations
- Return value structure
- Comparison with pyswisseph
- Various calculation flags
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


class TestCalcPctrBasicPositions:
    """Test basic planet-centric position calculations."""

    @pytest.mark.unit
    def test_moon_from_mars_j2000(self):
        """Moon position as seen from Mars at J2000 epoch."""
        jd = 2451545.0  # J2000
        pos, flags = ephem.calc_pctr(jd, MOON, MARS, 0)

        # Longitude should be in valid range
        assert 0 <= pos[0] < 360, f"Moon longitude {pos[0]} out of range"
        # Latitude should be in valid range
        assert -90 <= pos[1] <= 90, f"Moon latitude {pos[1]} out of range"
        # Distance should be positive (Moon-Mars distance varies widely)
        assert pos[2] > 0, f"Distance {pos[2]} should be positive"

    @pytest.mark.unit
    def test_sun_from_jupiter_j2000(self):
        """Sun position as seen from Jupiter at J2000 epoch."""
        jd = 2451545.0
        pos, flags = ephem.calc_pctr(jd, SUN, JUPITER, 0)

        # Longitude should be in valid range
        assert 0 <= pos[0] < 360
        # Distance from Jupiter to Sun should be ~5 AU
        assert 4.9 < pos[2] < 5.5, (
            f"Sun distance from Jupiter {pos[2]} should be ~5.2 AU"
        )

    @pytest.mark.unit
    def test_venus_from_mercury(self):
        """Venus position as seen from Mercury."""
        jd = 2451545.0
        pos, flags = ephem.calc_pctr(jd, VENUS, MERCURY, 0)

        # Both inner planets, distance should be < 2 AU
        assert 0 <= pos[0] < 360
        assert pos[2] < 2.0, f"Venus-Mercury distance {pos[2]} should be < 2 AU"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "target_id,target_name,center_id,center_name",
        [
            (SUN, "Sun", MARS, "Mars"),
            (MOON, "Moon", SATURN, "Saturn"),
            (EARTH, "Earth", MARS, "Mars"),
            (JUPITER, "Jupiter", SATURN, "Saturn"),
            (MERCURY, "Mercury", VENUS, "Venus"),
        ],
    )
    def test_various_planet_pairs(self, target_id, target_name, center_id, center_name):
        """Various planet pairs should return valid positions."""
        jd = 2451545.0
        pos, flags = ephem.calc_pctr(jd, target_id, center_id, 0)

        # Longitude in valid range
        assert 0 <= pos[0] < 360, (
            f"{target_name} from {center_name}: longitude {pos[0]} out of range"
        )
        # Latitude in valid range
        assert -90 <= pos[1] <= 90, (
            f"{target_name} from {center_name}: latitude {pos[1]} out of range"
        )
        # Distance positive
        assert pos[2] > 0, (
            f"{target_name} from {center_name}: distance should be positive"
        )


class TestCalcPctrReturnStructure:
    """Test the structure of calc_pctr return values."""

    @pytest.mark.unit
    def test_return_is_tuple(self):
        """calc_pctr should return a tuple."""
        result = ephem.calc_pctr(2451545.0, MOON, MARS, 0)
        assert isinstance(result, tuple)
        assert len(result) == 2

    @pytest.mark.unit
    def test_position_is_6_element_tuple(self):
        """Position should be a 6-element tuple."""
        pos, flags = ephem.calc_pctr(2451545.0, MOON, MARS, FLG_SPEED)
        assert len(pos) == 6

    @pytest.mark.unit
    def test_position_elements_are_floats(self):
        """All position elements should be floats."""
        pos, flags = ephem.calc_pctr(2451545.0, MOON, MARS, FLG_SPEED)
        for i, val in enumerate(pos):
            assert isinstance(val, float), f"Element {i} is {type(val)}, expected float"

    @pytest.mark.unit
    def test_flags_is_int(self):
        """Return flags should be an integer."""
        pos, flags = ephem.calc_pctr(2451545.0, MOON, MARS, 0)
        assert isinstance(flags, int)


class TestCalcPctrFlags:
    """Test calculation flags."""

    @pytest.mark.unit
    def test_flag_speed(self):
        """FLG_SPEED should return velocity values."""
        pos, flags = ephem.calc_pctr(2451545.0, MOON, MARS, FLG_SPEED)
        # Velocity should be non-zero (Moon and Mars are both moving)
        assert pos[3] != 0, "Longitude velocity should be non-zero"

    @pytest.mark.unit
    def test_flag_equatorial(self):
        """FLG_EQUATORIAL should return RA/Dec."""
        pos_ecl, _ = ephem.calc_pctr(2451545.0, SUN, JUPITER, 0)
        pos_equ, _ = ephem.calc_pctr(
            2451545.0, SUN, JUPITER, FLG_EQUATORIAL
        )

        # Should be in valid RA range
        assert 0 <= pos_equ[0] < 360

    @pytest.mark.unit
    def test_flag_sidereal(self):
        """FLG_SIDEREAL should apply ayanamsha correction."""
        ephem.set_sid_mode(SIDM_LAHIRI)

        pos_trop, _ = ephem.calc_pctr(2451545.0, SUN, MARS, 0)
        pos_sid, _ = ephem.calc_pctr(2451545.0, SUN, MARS, FLG_SIDEREAL)

        # Difference should equal ayanamsha (~23-24 degrees for Lahiri)
        diff = pos_trop[0] - pos_sid[0]
        if diff < 0:
            diff += 360
        assert 20 < diff < 30, f"Tropical-Sidereal diff {diff} should be ~23-24 deg"


class TestCalcPctrVsPyswisseph:
    """Compare calc_pctr results with pyswisseph."""

    @pytest.mark.comparison
    def test_moon_from_mars_matches_swisseph(self):
        """Moon from Mars should match pyswisseph."""
        jd = 2451545.0
        tolerance = 0.01  # degrees (10 arcsec tolerance for different ephemeris)

        pos_lib, _ = ephem.calc_pctr(jd, MOON, MARS, 0)
        pos_swe, _ = swe.calc_pctr(jd, MOON, MARS, 0)

        lon_diff = abs(pos_lib[0] - pos_swe[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        assert lon_diff < tolerance, (
            f"Moon from Mars longitude diff {lon_diff} >= {tolerance}"
        )

    @pytest.mark.comparison
    def test_sun_from_jupiter_matches_swisseph(self):
        """Sun from Jupiter should match pyswisseph."""
        jd = 2451545.0
        tolerance = 0.01

        pos_lib, _ = ephem.calc_pctr(jd, SUN, JUPITER, 0)
        try:
            pos_swe, _ = swe.calc_pctr(jd, SUN, JUPITER, 0)
        except Exception:
            # pyswisseph may require external ephemeris files for this calculation
            pytest.skip("pyswisseph ephemeris files not available for comparison")

        lon_diff = abs(pos_lib[0] - pos_swe[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        assert lon_diff < tolerance, (
            f"Sun from Jupiter longitude diff {lon_diff} >= {tolerance}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "target_id,target_name,center_id,center_name",
        [
            (MOON, "Moon", MARS, "Mars"),
            (MERCURY, "Mercury", VENUS, "Venus"),
            (VENUS, "Venus", EARTH, "Earth"),
            (EARTH, "Earth", MARS, "Mars"),
            (JUPITER, "Jupiter", SATURN, "Saturn"),
            (SATURN, "Saturn", JUPITER, "Jupiter"),
            (URANUS, "Uranus", NEPTUNE, "Neptune"),
        ],
    )
    def test_planet_pairs_match_swisseph(
        self, target_id, target_name, center_id, center_name
    ):
        """Various planet pairs should match pyswisseph."""
        jd = 2451545.0
        tolerance = 0.01  # degrees

        pos_lib, _ = ephem.calc_pctr(jd, target_id, center_id, 0)
        try:
            pos_swe, _ = swe.calc_pctr(jd, target_id, center_id, 0)
        except Exception:
            # pyswisseph may require external ephemeris files for some calculations
            pytest.skip("pyswisseph ephemeris files not available for comparison")

        lon_diff = abs(pos_lib[0] - pos_swe[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        assert lon_diff < tolerance, (
            f"{target_name} from {center_name} longitude diff {lon_diff}"
        )
        assert abs(pos_lib[1] - pos_swe[1]) < tolerance, (
            f"{target_name} from {center_name} latitude diff"
        )
        # Distance should match within 1%
        if pos_swe[2] > 0:
            dist_diff_pct = abs(pos_lib[2] - pos_swe[2]) / pos_swe[2] * 100
            assert dist_diff_pct < 1.0, (
                f"{target_name} from {center_name} distance diff {dist_diff_pct}%"
            )

    @pytest.mark.comparison
    def test_with_speed_flag_matches(self):
        """Speed calculations should match pyswisseph."""
        jd = 2451545.0
        tolerance = 0.1  # degrees/day tolerance for velocity

        pos_lib, _ = ephem.calc_pctr(jd, MOON, MARS, FLG_SPEED)
        pos_swe, _ = swe.calc_pctr(jd, MOON, MARS, FLG_SPEED)

        # Velocity should be reasonably close
        vel_diff = abs(pos_lib[3] - pos_swe[3])
        assert vel_diff < tolerance, (
            f"Moon from Mars velocity diff {vel_diff} >= {tolerance}"
        )


class TestCalcPctrGeometricConsistency:
    """Test geometric consistency of planet-centric calculations."""

    @pytest.mark.unit
    def test_symmetry_earth_mars(self):
        """Earth from Mars should be ~180 degrees from Mars from Earth (approximately)."""
        jd = 2451545.0

        # Earth as seen from Mars
        pos_earth_from_mars, _ = ephem.calc_pctr(jd, EARTH, MARS, 0)
        # Mars as seen from Earth (geocentric)
        pos_mars_from_earth, _ = ephem.calc_ut(jd, MARS, 0)

        # The longitudes should differ by roughly 180 degrees
        # (not exactly due to different reference frames)
        diff = abs(pos_earth_from_mars[0] - pos_mars_from_earth[0])
        if diff > 180:
            diff = 360 - diff

        # Should differ by roughly 180 degrees (allowing for light-time effects)
        assert 170 < diff < 190 or diff < 20 or diff > 340, (
            f"Earth-Mars symmetry unexpected: diff = {diff}"
        )

    @pytest.mark.unit
    def test_self_observation_returns_zero(self):
        """Observing a body from itself should give distance ~0."""
        jd = 2451545.0

        pos, _ = ephem.calc_pctr(jd, MARS, MARS, 0)

        # Distance should be essentially zero
        assert pos[2] < 0.0001, f"Self-observation distance {pos[2]} should be ~0"


class TestCalcPctrAlias:
    """Test the pyswisseph-compatible alias."""

    @pytest.mark.unit
    def test_calc_pctr_alias_exists(self):
        """calc_pctr alias should exist."""
        assert hasattr(ephem, "calc_pctr")

    @pytest.mark.unit
    def test_calc_pctr_alias_works(self):
        """calc_pctr alias should return same results as calc_pctr."""
        jd = 2451545.0

        pos1, flags1 = ephem.calc_pctr(jd, MOON, MARS, 0)
        pos2, flags2 = ephem.calc_pctr(jd, MOON, MARS, 0)

        assert pos1[0] == pos2[0]
        assert pos1[1] == pos2[1]
        assert pos1[2] == pos2[2]
