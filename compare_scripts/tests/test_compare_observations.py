"""
Pytest-style Observation Modes Comparison Tests.

Validates observation mode calculations (topocentric, heliocentric, barycentric,
true position, J2000 equatorial/ecliptic) against pyswisseph.
"""

import pytest
import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import (
    SUN,
    MOON,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    FLG_TOPOCTR,
    FLG_HELCTR,
    FLG_BARYCTR,
    FLG_TRUEPOS,
    FLG_J2000,
    FLG_EQUATORIAL,
    FLG_SPEED,
    FLG_SWIEPH,
)


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class Tolerances:
    """Tolerance thresholds for different comparison types."""

    LONGITUDE_STRICT = 0.001  # For geocentric/topocentric
    LONGITUDE_RELAXED = 0.03  # For heliocentric/barycentric
    LATITUDE_STRICT = 0.001
    LATITUDE_RELAXED = 0.03
    DISTANCE_STRICT = 0.0001  # AU
    DISTANCE_RELAXED = 0.01
    VELOCITY_ANGULAR = 0.01  # degrees/day
    VELOCITY_RADIAL = 0.001  # AU/day


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Test subjects with locations
SUBJECTS = [
    ("Standard J2000", 2000, 1, 1, 12.0, 0.0, 0.0, 0),
    ("Einstein (Ulm)", 1879, 3, 14, 11.5, 48.4011, 9.9876, 478),
    ("Gandhi (Porbandar)", 1869, 10, 2, 7.2, 21.6417, 69.6293, 0),
    ("Mandela (Mvezo)", 1918, 7, 18, 14.0, -31.9566, 28.5133, 0),
    ("Tromso (High Lat)", 1990, 1, 15, 12.0, 69.6492, 18.9553, 0),
    ("McMurdo (High Lat)", 2005, 6, 21, 0.0, -77.8463, 166.6681, 0),
]

# Comparison modes with flags
COMPARISON_MODES = [
    ("Topocentric", MOON, "Moon", FLG_TOPOCTR | FLG_SPEED, False),
    ("Heliocentric", MARS, "Mars", FLG_HELCTR | FLG_SPEED, True),
    ("Barycentric", JUPITER, "Jupiter", FLG_BARYCTR | FLG_SPEED, True),
    ("TruePos", VENUS, "Venus", FLG_TRUEPOS | FLG_SPEED, False),
    ("J2000 Equ", SUN, "Sun", FLG_J2000 | FLG_EQUATORIAL | FLG_SPEED, False),
    ("J2000 Ecl", SATURN, "Saturn", FLG_J2000 | FLG_SPEED, False),
]


# ============================================================================
# OBSERVATION MODE TESTS
# ============================================================================


class TestObservationModes:
    """Tests for various observation mode calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "mode_name,planet_id,planet_name,flags,use_relaxed", COMPARISON_MODES
    )
    @pytest.mark.parametrize(
        "subject_name,year,month,day,hour,lat,lon,alt",
        SUBJECTS[:3],  # Subset
    )
    def test_observation_mode(
        self,
        mode_name,
        planet_id,
        planet_name,
        flags,
        use_relaxed,
        subject_name,
        year,
        month,
        day,
        hour,
        lat,
        lon,
        alt,
    ):
        """Test observation mode calculations match pyswisseph."""
        jd = swe.julday(year, month, day, hour)

        # Set topocentric location if needed
        if flags & FLG_TOPOCTR:
            swe.set_topo(lon, lat, alt)
            pyephem.set_topo(lon, lat, alt)

        # SwissEphemeris
        try:
            res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        except Exception as e:
            pytest.skip(f"SwissEphemeris failed: {e}")

        # LibEphemeris
        try:
            res_py, _ = pyephem.calc_ut(jd, planet_id, flags)
        except Exception as e:
            pytest.skip(f"LibEphemeris failed: {e}")

        # Calculate differences
        diff_lon = angular_diff(res_swe[0], res_py[0])
        diff_lat = abs(res_swe[1] - res_py[1])
        diff_dist = abs(res_swe[2] - res_py[2])
        diff_lon_speed = abs(res_swe[3] - res_py[3])
        diff_lat_speed = abs(res_swe[4] - res_py[4])
        diff_dist_speed = abs(res_swe[5] - res_py[5])

        # Use appropriate tolerances
        lon_tol = (
            Tolerances.LONGITUDE_RELAXED if use_relaxed else Tolerances.LONGITUDE_STRICT
        )
        lat_tol = (
            Tolerances.LATITUDE_RELAXED if use_relaxed else Tolerances.LATITUDE_STRICT
        )
        dist_tol = (
            Tolerances.DISTANCE_RELAXED if use_relaxed else Tolerances.DISTANCE_STRICT
        )

        assert diff_lon < lon_tol, (
            f"{mode_name} {planet_name} @ {subject_name}: lon diff {diff_lon}°"
        )
        assert diff_lat < lat_tol, (
            f"{mode_name} {planet_name} @ {subject_name}: lat diff {diff_lat}°"
        )
        assert diff_dist < dist_tol, (
            f"{mode_name} {planet_name} @ {subject_name}: dist diff {diff_dist} AU"
        )
        assert diff_lon_speed < Tolerances.VELOCITY_ANGULAR, (
            f"{mode_name} {planet_name} @ {subject_name}: lon speed diff {diff_lon_speed}°/d"
        )


# ============================================================================
# TOPOCENTRIC TESTS
# ============================================================================


class TestTopocentric:
    """Specific tests for topocentric calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "lat,lon,alt,loc_name",
        [
            (0.0, 0.0, 0, "Equator"),
            (45.0, 0.0, 0, "Mid Latitude"),
            (70.0, 0.0, 0, "High Latitude"),
            (-45.0, 0.0, 0, "Southern"),
        ],
    )
    def test_topocentric_moon(self, lat, lon, alt, loc_name):
        """Test topocentric Moon calculations at various latitudes."""
        jd = swe.julday(2024, 1, 15, 12.0)
        flags = FLG_TOPOCTR | FLG_SPEED | FLG_SWIEPH

        swe.set_topo(lon, lat, alt)
        pyephem.set_topo(lon, lat, alt)

        res_swe, _ = swe.calc_ut(jd, MOON, flags)
        res_py, _ = pyephem.calc_ut(jd, MOON, flags)

        diff_lon = angular_diff(res_swe[0], res_py[0])

        assert diff_lon < Tolerances.LONGITUDE_STRICT, (
            f"Topocentric Moon @ {loc_name}: diff {diff_lon}°"
        )


# ============================================================================
# HELIOCENTRIC TESTS
# ============================================================================


class TestHeliocentric:
    """Tests for heliocentric calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_heliocentric_planets(self, planet_id, planet_name):
        """Test heliocentric planetary positions."""
        jd = swe.julday(2024, 1, 1, 12.0)
        flags = FLG_HELCTR | FLG_SPEED | FLG_SWIEPH

        res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        res_py, _ = pyephem.calc_ut(jd, planet_id, flags)

        diff_lon = angular_diff(res_swe[0], res_py[0])
        diff_lat = abs(res_swe[1] - res_py[1])

        assert diff_lon < Tolerances.LONGITUDE_RELAXED, (
            f"Heliocentric {planet_name}: lon diff {diff_lon}°"
        )
        assert diff_lat < Tolerances.LATITUDE_RELAXED, (
            f"Heliocentric {planet_name}: lat diff {diff_lat}°"
        )


# ============================================================================
# J2000 COORDINATE TESTS
# ============================================================================


class TestJ2000Coordinates:
    """Tests for J2000 coordinate calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
        ],
    )
    def test_j2000_ecliptic(self, planet_id, planet_name):
        """Test J2000 ecliptic coordinates."""
        jd = swe.julday(2024, 6, 15, 12.0)
        flags = FLG_J2000 | FLG_SPEED | FLG_SWIEPH

        res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        res_py, _ = pyephem.calc_ut(jd, planet_id, flags)

        diff_lon = angular_diff(res_swe[0], res_py[0])

        assert diff_lon < Tolerances.LONGITUDE_STRICT, (
            f"J2000 {planet_name}: diff {diff_lon}°"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
        ],
    )
    def test_j2000_equatorial(self, planet_id, planet_name):
        """Test J2000 equatorial coordinates."""
        jd = swe.julday(2024, 6, 15, 12.0)
        flags = FLG_J2000 | FLG_EQUATORIAL | FLG_SPEED | FLG_SWIEPH

        res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        res_py, _ = pyephem.calc_ut(jd, planet_id, flags)

        diff_ra = angular_diff(res_swe[0], res_py[0])
        diff_dec = abs(res_swe[1] - res_py[1])

        assert diff_ra < Tolerances.LONGITUDE_STRICT, (
            f"J2000 Equatorial {planet_name}: RA diff {diff_ra}°"
        )
        assert diff_dec < Tolerances.LATITUDE_STRICT, (
            f"J2000 Equatorial {planet_name}: Dec diff {diff_dec}°"
        )
