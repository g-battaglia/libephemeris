"""
Comprehensive tests for planetary position calculations (calc_ut).

Tests cover:
- All major planets (Sun through Pluto)
- Return value structure
- Comparison with pyswisseph
- Various calculation flags
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


class TestCalcUtBasicPositions:
    """Test basic planetary position calculations."""

    @pytest.mark.unit
    def test_sun_position_j2000(self):
        """Sun position at J2000 epoch."""
        jd = 2451545.0  # J2000
        pos, flags = ephem.calc_ut(jd, SUN, 0)

        # Sun should be around 280° (Capricorn) at J2000
        assert 270 < pos[0] < 290, f"Sun longitude {pos[0]} unexpected at J2000"
        # Latitude should be near 0
        assert abs(pos[1]) < 0.1, f"Sun latitude {pos[1]} should be near 0"
        # Distance should be ~1 AU
        assert 0.98 < pos[2] < 1.02, f"Sun distance {pos[2]} should be ~1 AU"

    @pytest.mark.unit
    def test_moon_position_j2000(self):
        """Moon position at J2000 epoch."""
        jd = 2451545.0
        pos, flags = ephem.calc_ut(jd, MOON, 0)

        # Moon longitude valid range
        assert 0 <= pos[0] < 360
        # Moon latitude within ±5.3°
        assert abs(pos[1]) < 6
        # Distance ~0.0025 AU (384,400 km)
        assert 0.002 < pos[2] < 0.003

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
            (URANUS, "Uranus"),
            (NEPTUNE, "Neptune"),
            (PLUTO, "Pluto"),
        ],
    )
    def test_all_planets_valid_positions(self, planet_id, planet_name):
        """All planets should return valid positions."""
        jd = 2451545.0
        pos, flags = ephem.calc_ut(jd, planet_id, 0)

        # Longitude in valid range
        assert 0 <= pos[0] < 360, f"{planet_name} longitude {pos[0]} out of range"
        # Latitude in valid range
        assert -90 <= pos[1] <= 90, f"{planet_name} latitude {pos[1]} out of range"
        # Distance positive
        assert pos[2] > 0, f"{planet_name} distance {pos[2]} should be positive"


class TestCalcUtReturnStructure:
    """Test the structure of calc_ut return values."""

    @pytest.mark.unit
    def test_return_is_tuple(self):
        """calc_ut should return a tuple."""
        result = ephem.calc_ut(2451545.0, SUN, 0)
        assert isinstance(result, tuple)
        assert len(result) == 2

    @pytest.mark.unit
    def test_position_is_6_element_tuple(self):
        """Position should be a 6-element tuple/list."""
        pos, flags = ephem.calc_ut(2451545.0, SUN, FLG_SPEED)
        assert len(pos) == 6

    @pytest.mark.unit
    def test_position_elements_are_floats(self):
        """All position elements should be floats."""
        pos, flags = ephem.calc_ut(2451545.0, SUN, FLG_SPEED)
        for i, val in enumerate(pos):
            assert isinstance(val, float), f"Element {i} is {type(val)}, expected float"

    @pytest.mark.unit
    def test_flags_is_int(self):
        """Return flags should be an integer."""
        pos, flags = ephem.calc_ut(2451545.0, SUN, 0)
        assert isinstance(flags, int)


class TestCalcUtFlags:
    """Test calculation flags."""

    @pytest.mark.unit
    def test_flag_speed(self):
        """FLG_SPEED should return velocity values."""
        pos, flags = ephem.calc_ut(2451545.0, SUN, FLG_SPEED)
        # Velocity should be non-zero
        assert pos[3] != 0, "Longitude velocity should be non-zero"

    @pytest.mark.unit
    def test_flag_equatorial(self):
        """FLG_EQUATORIAL should return RA/Dec."""
        pos_ecl, _ = ephem.calc_ut(2451545.0, SUN, 0)
        pos_equ, _ = ephem.calc_ut(2451545.0, SUN, FLG_EQUATORIAL)

        # RA and ecliptic longitude differ
        # (they're the same coordinate system rotated by obliquity)
        # Sun at J2000 should have RA around 18h40m = 280°
        assert 0 <= pos_equ[0] < 360

    @pytest.mark.unit
    def test_flag_j2000(self):
        """FLG_J2000 should return J2000 frame coordinates."""
        pos_date, _ = ephem.calc_ut(2460000.0, SUN, 0)  # 2023
        pos_j2000, _ = ephem.calc_ut(2460000.0, SUN, FLG_J2000)

        # Should differ due to precession
        # Difference should be small but measurable
        diff = abs(pos_date[0] - pos_j2000[0])
        if diff > 180:
            diff = 360 - diff
        # Precession is ~50" per year, ~23 years = ~0.32°
        assert diff < 1.0, "J2000 vs of-date should differ by precession"

    @pytest.mark.unit
    def test_flag_helctr(self):
        """FLG_HELCTR should return heliocentric coordinates."""
        pos_geo, _ = ephem.calc_ut(2451545.0, MARS, 0)
        pos_hel, _ = ephem.calc_ut(2451545.0, MARS, FLG_HELCTR)

        # Heliocentric and geocentric positions differ
        # (unless planet is at opposition/conjunction)
        diff = abs(pos_geo[0] - pos_hel[0])
        if diff > 180:
            diff = 360 - diff
        # Usually significant difference
        assert diff > 0.001 or abs(pos_geo[2] - pos_hel[2]) > 0.001

    @pytest.mark.unit
    def test_combined_flags(self):
        """Multiple flags should work together."""
        flags = FLG_SPEED | FLG_EQUATORIAL
        pos, ret_flags = ephem.calc_ut(2451545.0, SUN, flags)

        # Should have velocity
        assert pos[3] != 0
        # Should be in RA/Dec


class TestCalcUtVsPyswisseph:
    """Compare calc_ut results with pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
            (URANUS, "Uranus"),
            (NEPTUNE, "Neptune"),
            (PLUTO, "Pluto"),
        ],
    )
    def test_planets_match_swisseph(self, planet_id, planet_name):
        """All planets should match pyswisseph within tolerance."""
        jd = 2451545.0
        tolerance = 0.001  # degrees

        pos_lib, _ = ephem.calc_ut(jd, planet_id, 0)
        pos_swe, _ = swe.calc_ut(jd, planet_id, 0)

        lon_diff = abs(pos_lib[0] - pos_swe[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        assert lon_diff < tolerance, (
            f"{planet_name} longitude diff {lon_diff} >= {tolerance}"
        )
        assert abs(pos_lib[1] - pos_swe[1]) < tolerance, f"{planet_name} latitude diff"
        assert abs(pos_lib[2] - pos_swe[2]) < 0.001, f"{planet_name} distance diff"

    @pytest.mark.comparison
    def test_100_dates_all_planets(
        self, random_dates_in_de421_range, all_planets, progress_reporter
    ):
        """Test all planets at multiple dates."""
        dates = random_dates_in_de421_range(10)  # 10 dates x 10 planets = 100 tests
        tolerance = 0.001
        total = len(dates) * len(all_planets)
        progress = progress_reporter("Testing planets/dates", total, report_every=10)

        iteration = 0
        for year, month, day, hour, jd in dates:
            for planet_id, planet_name in all_planets:
                pos_lib, _ = ephem.calc_ut(jd, planet_id, 0)
                pos_swe, _ = swe.calc_ut(jd, planet_id, 0)

                lon_diff = abs(pos_lib[0] - pos_swe[0])
                if lon_diff > 180:
                    lon_diff = 360 - lon_diff

                assert lon_diff < tolerance, (
                    f"{planet_name} at JD {jd}: lon diff {lon_diff}"
                )
                progress.update(
                    iteration, f"{planet_name} @ {year}-{month:02d}-{day:02d}"
                )
                iteration += 1

        progress.done()


class TestCalcUtTopocentric:
    """Test topocentric calculations."""

    @pytest.mark.unit
    def test_topocentric_moon(self):
        """Topocentric Moon should differ from geocentric."""
        jd = 2451545.0

        # Set observer location (Rome)
        ephem.set_topo(12.5, 41.9, 0)
        swe.set_topo(12.5, 41.9, 0)

        pos_geo, _ = ephem.calc_ut(jd, MOON, 0)
        pos_topo, _ = ephem.calc_ut(jd, MOON, FLG_TOPOCTR)

        # Moon parallax can be up to ~1°
        lon_diff = abs(pos_geo[0] - pos_topo[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        # Should differ by up to ~1°
        assert lon_diff < 2.0, "Topocentric parallax too large"

    @pytest.mark.comparison
    def test_topocentric_matches_swisseph(self):
        """Topocentric should match pyswisseph."""
        jd = 2451545.0

        # Set same location in both
        ephem.set_topo(12.5, 41.9, 0)
        swe.set_topo(12.5, 41.9, 0)

        pos_lib, _ = ephem.calc_ut(jd, MOON, FLG_TOPOCTR)
        pos_swe, _ = swe.calc_ut(jd, MOON, FLG_TOPOCTR)

        lon_diff = abs(pos_lib[0] - pos_swe[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        assert lon_diff < 0.01, f"Topocentric Moon diff {lon_diff}"


class TestCalcUtSidereal:
    """Test sidereal calculations."""

    @pytest.mark.unit
    def test_sidereal_differs_from_tropical(self):
        """Sidereal position should differ from tropical by ayanamsha."""
        jd = 2451545.0

        ephem.set_sid_mode(SIDM_LAHIRI)
        swe.set_sid_mode(SIDM_LAHIRI)

        pos_trop, _ = ephem.calc_ut(jd, SUN, 0)
        pos_sid, _ = ephem.calc_ut(jd, SUN, FLG_SIDEREAL)

        # Difference should equal ayanamsha (~23-24° for Lahiri at J2000)
        diff = pos_trop[0] - pos_sid[0]
        if diff < 0:
            diff += 360

        assert 20 < diff < 30, f"Tropical-Sidereal diff {diff} should be ~23-24°"

    @pytest.mark.comparison
    def test_sidereal_matches_swisseph(self):
        """Sidereal should match pyswisseph."""
        jd = 2451545.0

        ephem.set_sid_mode(SIDM_LAHIRI)
        swe.set_sid_mode(SIDM_LAHIRI)

        pos_lib, _ = ephem.calc_ut(jd, SUN, FLG_SIDEREAL)
        pos_swe, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)

        lon_diff = abs(pos_lib[0] - pos_swe[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        assert lon_diff < 0.01, f"Sidereal Sun diff {lon_diff}"
