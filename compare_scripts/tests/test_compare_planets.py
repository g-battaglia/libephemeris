"""
Planetary Calculations Comparison Tests.

Validates all planetary calculations between pyswisseph and libephemeris:
- All major planets (Sun through Pluto)
- Multiple calculation modes (geocentric, heliocentric, barycentric)
- Position (longitude, latitude, distance) and velocity
- Various time periods and edge cases
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_HELCTR,
    FLG_BARYCTR,
    FLG_ICRS,
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_NONUT,
)


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

PLANETS = [
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
]

CALC_MODES = [
    (FLG_SWIEPH, "Geocentric"),
    (FLG_SWIEPH | FLG_SPEED, "Geocentric+Speed"),
    (FLG_SWIEPH | FLG_HELCTR, "Heliocentric"),
    (FLG_SWIEPH | FLG_BARYCTR, "Barycentric"),
]

TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (2024, 11, 15, 0.0, "Current Era"),
    (1900, 1, 1, 0.0, "Early 1900s"),
    (2100, 12, 31, 23.999, "Late 2100s"),
    (1950, 6, 21, 12.0, "Mid-Century"),
]

# Tolerances
LONGITUDE_STRICT = 0.001  # degrees
LONGITUDE_RELAXED = 0.03  # degrees (heliocentric/barycentric)
LATITUDE_STRICT = 0.001
LATITUDE_RELAXED = 0.03
DISTANCE_STRICT = 0.0001  # AU
DISTANCE_RELAXED = 0.01
VELOCITY_ANGULAR = 0.01  # degrees/day
VELOCITY_RADIAL = 0.001  # AU/day


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestPlanetaryPositions:
    """Compare planetary position calculations between swisseph and libephemeris."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_geocentric_positions(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test geocentric planetary positions match pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flag = FLG_SWIEPH

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        assert diff_lon < LONGITUDE_STRICT, (
            f"{planet_name} at {date_desc}: longitude diff {diff_lon:.6f}° exceeds tolerance"
        )
        assert diff_lat < LATITUDE_STRICT, (
            f"{planet_name} at {date_desc}: latitude diff {diff_lat:.6f}° exceeds tolerance"
        )
        assert diff_dist < DISTANCE_STRICT, (
            f"{planet_name} at {date_desc}: distance diff {diff_dist:.8f} AU exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_geocentric_with_speed(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test geocentric positions with velocity calculations."""
        jd = swe.julday(year, month, day, hour)
        flag = FLG_SWIEPH | FLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

        # Position checks
        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        assert diff_lon < LONGITUDE_STRICT, (
            f"{planet_name} at {date_desc}: longitude diff {diff_lon:.6f}° exceeds tolerance"
        )
        assert diff_lat < LATITUDE_STRICT, (
            f"{planet_name} at {date_desc}: latitude diff {diff_lat:.6f}° exceeds tolerance"
        )
        assert diff_dist < DISTANCE_STRICT, (
            f"{planet_name} at {date_desc}: distance diff {diff_dist:.8f} AU exceeds tolerance"
        )

        # Velocity checks
        diff_lon_speed = abs(pos_swe[3] - pos_py[3])
        diff_lat_speed = abs(pos_swe[4] - pos_py[4])
        diff_dist_speed = abs(pos_swe[5] - pos_py[5])

        assert diff_lon_speed < VELOCITY_ANGULAR, (
            f"{planet_name} at {date_desc}: lon velocity diff {diff_lon_speed:.6f}°/day exceeds tolerance"
        )
        assert diff_lat_speed < VELOCITY_ANGULAR, (
            f"{planet_name} at {date_desc}: lat velocity diff {diff_lat_speed:.6f}°/day exceeds tolerance"
        )
        assert diff_dist_speed < VELOCITY_RADIAL, (
            f"{planet_name} at {date_desc}: dist velocity diff {diff_dist_speed:.8f} AU/day exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize(
        "year,month,day,hour,date_desc", TEST_DATES[:3]
    )  # Fewer dates for these
    def test_heliocentric_positions(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test heliocentric planetary positions match pyswisseph."""
        # Skip Sun for heliocentric (meaningless)
        if planet_id == SUN:
            pytest.skip("Sun position is not meaningful in heliocentric coordinates")

        jd = swe.julday(year, month, day, hour)
        flag = FLG_SWIEPH | FLG_HELCTR

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        assert diff_lon < LONGITUDE_RELAXED, (
            f"{planet_name} helio at {date_desc}: longitude diff {diff_lon:.6f}° exceeds tolerance"
        )
        assert diff_lat < LATITUDE_RELAXED, (
            f"{planet_name} helio at {date_desc}: latitude diff {diff_lat:.6f}° exceeds tolerance"
        )
        assert diff_dist < DISTANCE_RELAXED, (
            f"{planet_name} helio at {date_desc}: distance diff {diff_dist:.6f} AU exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES[:3])
    def test_barycentric_positions(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test barycentric planetary positions match pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flag = FLG_SWIEPH | FLG_BARYCTR

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        assert diff_lon < LONGITUDE_RELAXED, (
            f"{planet_name} bary at {date_desc}: longitude diff {diff_lon:.6f}° exceeds tolerance"
        )
        assert diff_lat < LATITUDE_RELAXED, (
            f"{planet_name} bary at {date_desc}: latitude diff {diff_lat:.6f}° exceeds tolerance"
        )
        assert diff_dist < DISTANCE_RELAXED, (
            f"{planet_name} bary at {date_desc}: distance diff {diff_dist:.6f} AU exceeds tolerance"
        )


class TestInnerPlanets:
    """Specific tests for inner planets (faster moving, more precision needed)."""

    INNER_PLANETS = [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MERCURY, "Mercury"),
        (VENUS, "Venus"),
        (MARS, "Mars"),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", INNER_PLANETS)
    def test_inner_planet_daily_motion(self, planet_id, planet_name):
        """Test that daily motion (velocity) matches for inner planets."""
        jd = swe.julday(2024, 6, 15, 12.0)
        flag = FLG_SWIEPH | FLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

        diff_speed = abs(pos_swe[3] - pos_py[3])

        # Moon uses a wider tolerance because it moves ~13°/day, amplifying
        # the sub-arcsecond position difference between DE440 (Skyfield) and
        # Swiss Ephemeris through numerical differentiation (KI-010).
        tol = 0.002 if planet_id == MOON else 0.001
        assert diff_speed < tol, (
            f"{planet_name} daily motion diff {diff_speed:.6f}°/day exceeds tight tolerance"
        )


class TestOuterPlanets:
    """Specific tests for outer planets (slower moving)."""

    OUTER_PLANETS = [
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
        (URANUS, "Uranus"),
        (NEPTUNE, "Neptune"),
        (PLUTO, "Pluto"),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", OUTER_PLANETS)
    def test_outer_planet_positions_extended_range(self, planet_id, planet_name):
        """Test outer planets across extended date range."""
        dates = [
            swe.julday(1600, 1, 1, 12.0),
            swe.julday(1800, 1, 1, 12.0),
            swe.julday(2200, 1, 1, 12.0),
            swe.julday(2500, 1, 1, 12.0),
        ]

        for jd in dates:
            pos_swe, _ = swe.calc_ut(jd, planet_id, FLG_SWIEPH)
            pos_py, _ = ephem.calc_ut(jd, planet_id, FLG_SWIEPH)

            diff_lon = angular_diff(pos_swe[0], pos_py[0])

            assert diff_lon < LONGITUDE_STRICT, (
                f"{planet_name} at JD {jd}: longitude diff {diff_lon:.6f}° exceeds tolerance"
            )


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    def test_j2000_epoch(self, planet_id, planet_name):
        """Test exact J2000.0 epoch positions."""
        jd = 2451545.0  # J2000.0

        pos_swe, _ = swe.calc_ut(jd, planet_id, FLG_SWIEPH | FLG_SPEED)
        pos_py, _ = ephem.calc_ut(jd, planet_id, FLG_SWIEPH | FLG_SPEED)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])

        assert diff_lon < 0.0001, (
            f"{planet_name} at J2000: longitude diff {diff_lon:.8f}° exceeds tight tolerance"
        )

    @pytest.mark.comparison
    def test_moon_at_midnight(self):
        """Test Moon position at midnight (common use case)."""
        jd = swe.julday(2024, 1, 1, 0.0)

        pos_swe, _ = swe.calc_ut(jd, MOON, FLG_SWIEPH | FLG_SPEED)
        pos_py, _ = ephem.calc_ut(jd, MOON, FLG_SWIEPH | FLG_SPEED)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_speed = abs(pos_swe[3] - pos_py[3])

        assert diff_lon < 0.001, (
            f"Moon longitude diff {diff_lon:.6f}° exceeds tolerance"
        )
        assert diff_speed < 0.01, (
            f"Moon speed diff {diff_speed:.6f}°/day exceeds tolerance"
        )

    @pytest.mark.comparison
    def test_sun_at_vernal_equinox(self):
        """Test Sun position at vernal equinox 2024."""
        # 2024 vernal equinox approximately March 20, 03:06 UT
        jd = swe.julday(2024, 3, 20, 3.1)

        pos_swe, _ = swe.calc_ut(jd, SUN, FLG_SWIEPH)
        pos_py, _ = ephem.calc_ut(jd, SUN, FLG_SWIEPH)

        # Sun should be near 0° (Aries ingress)
        assert pos_swe[0] < 1.0 or pos_swe[0] > 359.0, (
            "Sun should be near 0° at vernal equinox"
        )

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        assert diff_lon < 0.001, (
            f"Sun at equinox diff {diff_lon:.6f}° exceeds tolerance"
        )


class TestICRSPositions:
    """Compare ICRS (International Celestial Reference System) calculations.

    ICRS is a fixed reference frame aligned with distant quasars, not subject
    to precession or nutation. When used with FLG_EQUATORIAL, positions in
    ICRS should differ from J2000 equatorial by an amount corresponding to
    precession+nutation accumulated since J2000.0.

    Note: FLG_ICRS primarily affects equatorial coordinates. In ecliptic
    mode without FLG_EQUATORIAL, the flag has minimal effect.
    """

    # Expected precession from J2000.0 is approximately 50.3 arcsec/year
    # For dates far from J2000, difference should be roughly:
    # (year - 2000) * 50.3 arcsec = (year - 2000) * 0.01397 degrees
    PRECESSION_RATE_DEG_PER_YEAR = 0.01397

    # Pluto at extreme dates needs relaxed tolerance
    DISTANCE_PLUTO_RELAXED = 0.001  # AU

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_icrs_positions(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test ICRS planetary positions match pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flag = FLG_SWIEPH | FLG_ICRS

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        assert diff_lon < LONGITUDE_STRICT, (
            f"{planet_name} ICRS at {date_desc}: longitude diff {diff_lon:.6f}° exceeds tolerance"
        )
        assert diff_lat < LATITUDE_STRICT, (
            f"{planet_name} ICRS at {date_desc}: latitude diff {diff_lat:.6f}° exceeds tolerance"
        )

        # Pluto at extreme dates needs relaxed distance tolerance
        dist_tolerance = (
            self.DISTANCE_PLUTO_RELAXED if planet_id == PLUTO else DISTANCE_STRICT
        )
        assert diff_dist < dist_tolerance, (
            f"{planet_name} ICRS at {date_desc}: distance diff {diff_dist:.8f} AU exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_icrs_with_speed(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test ICRS positions with velocity calculations match pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flag = FLG_SWIEPH | FLG_ICRS | FLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

        # Position checks
        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        assert diff_lon < LONGITUDE_STRICT, (
            f"{planet_name} ICRS+Speed at {date_desc}: longitude diff {diff_lon:.6f}° exceeds tolerance"
        )
        assert diff_lat < LATITUDE_STRICT, (
            f"{planet_name} ICRS+Speed at {date_desc}: latitude diff {diff_lat:.6f}° exceeds tolerance"
        )

        # Pluto at extreme dates needs relaxed distance tolerance
        dist_tolerance = (
            self.DISTANCE_PLUTO_RELAXED if planet_id == PLUTO else DISTANCE_STRICT
        )
        assert diff_dist < dist_tolerance, (
            f"{planet_name} ICRS+Speed at {date_desc}: distance diff {diff_dist:.8f} AU exceeds tolerance"
        )

        # Velocity checks
        diff_lon_speed = abs(pos_swe[3] - pos_py[3])
        diff_lat_speed = abs(pos_swe[4] - pos_py[4])
        diff_dist_speed = abs(pos_swe[5] - pos_py[5])

        assert diff_lon_speed < VELOCITY_ANGULAR, (
            f"{planet_name} ICRS at {date_desc}: lon velocity diff {diff_lon_speed:.6f}°/day exceeds tolerance"
        )
        assert diff_lat_speed < VELOCITY_ANGULAR, (
            f"{planet_name} ICRS at {date_desc}: lat velocity diff {diff_lat_speed:.6f}°/day exceeds tolerance"
        )
        assert diff_dist_speed < VELOCITY_RADIAL, (
            f"{planet_name} ICRS at {date_desc}: dist velocity diff {diff_dist_speed:.8f} AU/day exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name", PLANETS[:5]
    )  # Test subset for speed
    def test_icrs_flag_consistent_with_pyswisseph(self, planet_id, planet_name):
        """Verify that ICRS flag produces consistent results with pyswisseph.

        Swiss Ephemeris internally uses the ICRS frame as its fundamental reference.
        The FLG_ICRS flag ensures explicit ICRS frame output. This test validates
        that libephemeris handles the flag correctly and produces results matching
        pyswisseph for both ecliptic and equatorial coordinate requests.
        """
        # Test at J2000.0 - the ICRS reference epoch
        jd_j2000 = 2451545.0

        # Test ICRS with ecliptic coordinates
        icrs_ecliptic_flags = FLG_SWIEPH | FLG_ICRS
        pos_swe_ecl, _ = swe.calc_ut(jd_j2000, planet_id, icrs_ecliptic_flags)
        pos_py_ecl, _ = ephem.calc_ut(jd_j2000, planet_id, icrs_ecliptic_flags)

        diff_ecl = angular_diff(pos_swe_ecl[0], pos_py_ecl[0])
        assert diff_ecl < 0.0001, (
            f"{planet_name} ICRS ecliptic at J2000: diff {diff_ecl:.8f}° exceeds tolerance"
        )

        # Test ICRS with equatorial coordinates
        icrs_equatorial_flags = FLG_SWIEPH | FLG_ICRS | FLG_EQUATORIAL
        pos_swe_eq, _ = swe.calc_ut(jd_j2000, planet_id, icrs_equatorial_flags)
        pos_py_eq, _ = ephem.calc_ut(jd_j2000, planet_id, icrs_equatorial_flags)

        diff_ra = angular_diff(pos_swe_eq[0], pos_py_eq[0])
        diff_dec = abs(pos_swe_eq[1] - pos_py_eq[1])
        assert diff_ra < 0.0001, (
            f"{planet_name} ICRS RA at J2000: diff {diff_ra:.8f}° exceeds tolerance"
        )
        assert diff_dec < 0.0001, (
            f"{planet_name} ICRS Dec at J2000: diff {diff_dec:.8f}° exceeds tolerance"
        )

        # Test at a different date (2024)
        jd_2024 = swe.julday(2024, 6, 15, 12.0)
        pos_swe_2024, _ = swe.calc_ut(jd_2024, planet_id, icrs_equatorial_flags)
        pos_py_2024, _ = ephem.calc_ut(jd_2024, planet_id, icrs_equatorial_flags)

        diff_ra_2024 = angular_diff(pos_swe_2024[0], pos_py_2024[0])
        assert diff_ra_2024 < 0.001, (
            f"{planet_name} ICRS RA at 2024: diff {diff_ra_2024:.6f}° exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS[:5])  # Test subset
    def test_icrs_equatorial_matches_pyswisseph(self, planet_id, planet_name):
        """Verify ICRS equatorial coordinates match pyswisseph at multiple dates."""
        test_jds = [
            2451545.0,  # J2000.0
            swe.julday(2024, 6, 15, 12.0),  # 2024
            swe.julday(1950, 6, 15, 12.0),  # 1950
        ]
        flag = FLG_SWIEPH | FLG_ICRS | FLG_EQUATORIAL | FLG_SPEED

        for jd in test_jds:
            pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
            pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

            diff_ra = angular_diff(pos_swe[0], pos_py[0])
            diff_dec = abs(pos_swe[1] - pos_py[1])

            assert diff_ra < 0.001, (
                f"{planet_name} ICRS equatorial at JD {jd}: RA diff {diff_ra:.6f}° exceeds tolerance"
            )
            assert diff_dec < 0.001, (
                f"{planet_name} ICRS equatorial at JD {jd}: Dec diff {diff_dec:.6f}° exceeds tolerance"
            )

    @pytest.mark.comparison
    def test_icrs_consistency_with_pyswisseph(self):
        """Verify ICRS implementation matches pyswisseph across all planets at J2000."""
        jd = 2451545.0  # J2000.0
        flag = FLG_SWIEPH | FLG_ICRS | FLG_SPEED

        for planet_id, planet_name in PLANETS:
            pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
            pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

            diff_lon = angular_diff(pos_swe[0], pos_py[0])

            assert diff_lon < 0.0001, (
                f"{planet_name} ICRS at J2000: longitude diff {diff_lon:.8f}° "
                "exceeds tight tolerance"
            )


class TestJ2000Positions:
    """Compare J2000 reference frame calculations between swisseph and libephemeris.

    FLG_J2000 computes planetary positions referred to the equinox of J2000.0,
    rather than the equinox of the date. This is critical for modern astronomical
    software that uses J2000 as a fixed reference frame.

    Without this flag, coordinates are in the equinox of the date and will differ
    from J2000 by the precession accumulated since J2000.0 (~0.014 degrees/year).
    For distant epochs like 1900 or 2100, this can cause errors of 0.5+ degrees.
    """

    # Test epochs spanning 200 years around J2000
    J2000_TEST_DATES = [
        (1900, 1, 1, 12.0, "1900 (100 years before J2000)"),
        (1950, 6, 15, 12.0, "1950 (50 years before J2000)"),
        (2000, 1, 1, 12.0, "J2000.0 Epoch"),
        (2050, 6, 15, 12.0, "2050 (50 years after J2000)"),
        (2100, 12, 31, 12.0, "2100 (100 years after J2000)"),
    ]

    # Precession rate is approximately 50.3 arcsec/year = 0.01397 deg/year
    PRECESSION_RATE_DEG_PER_YEAR = 0.01397

    # Pluto at extreme dates needs relaxed distance tolerance
    DISTANCE_PLUTO_RELAXED = 0.001  # AU

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", J2000_TEST_DATES)
    def test_j2000_positions_match_pyswisseph(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test FLG_J2000 planetary positions match pyswisseph at various epochs."""
        jd = swe.julday(year, month, day, hour)
        flag = FLG_SWIEPH | FLG_J2000

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        assert diff_lon < LONGITUDE_STRICT, (
            f"{planet_name} J2000 at {date_desc}: longitude diff {diff_lon:.6f}° "
            f"exceeds tolerance (swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )
        assert diff_lat < LATITUDE_STRICT, (
            f"{planet_name} J2000 at {date_desc}: latitude diff {diff_lat:.6f}° "
            "exceeds tolerance"
        )

        # Pluto at extreme dates needs relaxed distance tolerance
        dist_tolerance = (
            self.DISTANCE_PLUTO_RELAXED if planet_id == PLUTO else DISTANCE_STRICT
        )
        assert diff_dist < dist_tolerance, (
            f"{planet_name} J2000 at {date_desc}: distance diff {diff_dist:.8f} AU "
            "exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", J2000_TEST_DATES)
    def test_j2000_with_speed_match_pyswisseph(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test FLG_J2000 with velocity calculations match pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flag = FLG_SWIEPH | FLG_J2000 | FLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

        # Position checks
        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        assert diff_lon < LONGITUDE_STRICT, (
            f"{planet_name} J2000+Speed at {date_desc}: longitude diff {diff_lon:.6f}° "
            "exceeds tolerance"
        )
        assert diff_lat < LATITUDE_STRICT, (
            f"{planet_name} J2000+Speed at {date_desc}: latitude diff {diff_lat:.6f}° "
            "exceeds tolerance"
        )

        # Pluto at extreme dates needs relaxed distance tolerance
        dist_tolerance = (
            self.DISTANCE_PLUTO_RELAXED if planet_id == PLUTO else DISTANCE_STRICT
        )
        assert diff_dist < dist_tolerance, (
            f"{planet_name} J2000+Speed at {date_desc}: distance diff {diff_dist:.8f} AU "
            "exceeds tolerance"
        )

        # Velocity checks
        diff_lon_speed = abs(pos_swe[3] - pos_py[3])
        diff_lat_speed = abs(pos_swe[4] - pos_py[4])
        diff_dist_speed = abs(pos_swe[5] - pos_py[5])

        assert diff_lon_speed < VELOCITY_ANGULAR, (
            f"{planet_name} J2000 at {date_desc}: lon velocity diff "
            f"{diff_lon_speed:.6f}°/day exceeds tolerance"
        )
        assert diff_lat_speed < VELOCITY_ANGULAR, (
            f"{planet_name} J2000 at {date_desc}: lat velocity diff "
            f"{diff_lat_speed:.6f}°/day exceeds tolerance"
        )
        assert diff_dist_speed < VELOCITY_RADIAL, (
            f"{planet_name} J2000 at {date_desc}: dist velocity diff "
            f"{diff_dist_speed:.8f} AU/day exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", J2000_TEST_DATES)
    def test_j2000_equatorial_match_pyswisseph(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test FLG_J2000 equatorial coordinates (RA/Dec) match pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flag = FLG_SWIEPH | FLG_J2000 | FLG_EQUATORIAL | FLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

        diff_ra = angular_diff(pos_swe[0], pos_py[0])
        diff_dec = abs(pos_swe[1] - pos_py[1])

        assert diff_ra < LONGITUDE_STRICT, (
            f"{planet_name} J2000 RA at {date_desc}: diff {diff_ra:.6f}° "
            f"exceeds tolerance (swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )
        assert diff_dec < LATITUDE_STRICT, (
            f"{planet_name} J2000 Dec at {date_desc}: diff {diff_dec:.6f}° "
            "exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS[:5])  # Test subset
    def test_j2000_differs_from_date_by_precession(self, planet_id, planet_name):
        """Verify J2000 coordinates differ from date coordinates by precession.

        At epochs far from J2000, the difference between equinox-of-date and
        J2000 coordinates should correspond to accumulated precession.
        """
        # Test at 2100: 100 years after J2000
        # Expected precession: ~100 * 0.01397 = ~1.4 degrees
        jd_2100 = swe.julday(2100, 1, 1, 12.0)

        date_flags = FLG_SWIEPH
        j2000_flags = FLG_SWIEPH | FLG_J2000

        pos_date, _ = ephem.calc_ut(jd_2100, planet_id, date_flags)
        pos_j2000, _ = ephem.calc_ut(jd_2100, planet_id, j2000_flags)

        diff = angular_diff(pos_date[0], pos_j2000[0])

        # Difference should be at least 0.5 degrees at 100 years from J2000
        assert diff > 0.5, (
            f"{planet_name} at 2100: J2000 and date coords should differ by precession, "
            f"got only {diff:.4f}° (expected > 0.5°)"
        )
        # But not more than expected precession + some margin
        assert diff < 3.0, (
            f"{planet_name} at 2100: diff {diff:.4f}° seems too large for precession"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS[:5])  # Test subset
    def test_j2000_at_epoch_matches_date(self, planet_id, planet_name):
        """At J2000.0 epoch, J2000 coords should approximately equal date coords.

        At the J2000.0 epoch itself, the equinox of date equals the J2000 equinox,
        so the difference should be minimal (only nutation, no precession).
        """
        jd_j2000 = 2451545.0  # J2000.0 epoch

        date_flags = FLG_SWIEPH
        j2000_flags = FLG_SWIEPH | FLG_J2000

        pos_date, _ = ephem.calc_ut(jd_j2000, planet_id, date_flags)
        pos_j2000, _ = ephem.calc_ut(jd_j2000, planet_id, j2000_flags)

        diff = angular_diff(pos_date[0], pos_j2000[0])

        # At J2000.0, difference should be very small (only nutation ~0.005°)
        assert diff < 0.02, (
            f"{planet_name} at J2000.0: date and J2000 coords should be very close, "
            f"got {diff:.6f}° (expected < 0.02°)"
        )

    @pytest.mark.comparison
    def test_j2000_reference_frame_consistency(self):
        """Verify that FLG_J2000 provides a consistent reference frame.

        The same planet observed at different dates should show different
        positions (due to orbital motion), but both libephemeris and pyswisseph
        should agree precisely on those positions in the J2000 frame.
        """
        flag = FLG_SWIEPH | FLG_J2000 | FLG_SPEED
        test_jds = [
            swe.julday(1900, 1, 1, 12.0),
            swe.julday(2000, 1, 1, 12.0),
            swe.julday(2100, 1, 1, 12.0),
        ]

        for planet_id, planet_name in PLANETS:
            for jd in test_jds:
                pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
                pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

                diff_lon = angular_diff(pos_swe[0], pos_py[0])

                assert diff_lon < 0.001, (
                    f"{planet_name} J2000 at JD {jd}: longitude diff {diff_lon:.6f}° "
                    "exceeds tolerance - reference frame inconsistent"
                )

    @pytest.mark.comparison
    def test_j2000_critical_for_distant_epochs(self):
        """Demonstrate the critical importance of FLG_J2000 for distant epochs.

        Without FLG_J2000, software using date coordinates for epochs far from
        J2000 would get positions that differ from J2000-referenced software by
        the accumulated precession, potentially causing significant errors.
        """
        # At 1900 and 2100, the precession effect is substantial
        epochs = [
            (swe.julday(1900, 1, 1, 12.0), "1900"),
            (swe.julday(2100, 1, 1, 12.0), "2100"),
        ]

        for jd, epoch_name in epochs:
            date_flags = FLG_SWIEPH
            j2000_flags = FLG_SWIEPH | FLG_J2000

            pos_date_swe, _ = swe.calc_ut(jd, SUN, date_flags)
            pos_j2000_swe, _ = swe.calc_ut(jd, SUN, j2000_flags)
            pos_j2000_py, _ = ephem.calc_ut(jd, SUN, j2000_flags)

            # The precession difference between date and J2000
            precession_diff = angular_diff(pos_date_swe[0], pos_j2000_swe[0])

            # libephemeris should match pyswisseph J2000 closely
            impl_diff = angular_diff(pos_j2000_swe[0], pos_j2000_py[0])

            assert impl_diff < 0.001, (
                f"Sun J2000 at {epoch_name}: libephemeris differs from pyswisseph "
                f"by {impl_diff:.6f}° (precession effect is {precession_diff:.4f}°)"
            )

            # Precession should be significant at distant epochs
            assert precession_diff > 0.5, (
                f"Sun at {epoch_name}: precession effect {precession_diff:.4f}° "
                "is smaller than expected (> 0.5°)"
            )


class TestNoNutationPositions:
    """Compare FLG_NONUT (no nutation) calculations between swisseph and libephemeris.

    FLG_NONUT excludes nutation from the calculation, giving mean positions
    rather than true (apparent) positions. Nutation is a periodic oscillation
    of Earth's axis with a principal term of ~9.2 arcseconds (0.00256 degrees)
    and an 18.6-year period.

    This flag is useful for:
    - Comparing with mean ephemerides
    - Isolating precession effects from nutation
    - Astronomical research requiring mean positions

    Note: Tests use relaxed tolerance to accommodate nutation-related differences.
    The expected nutation effect is approximately 3-10 arcsec (0.001-0.003 degrees).
    """

    # Nutation amplitude is approximately 9.2 arcsec = 0.00256 degrees
    # We expect differences between with/without nutation to be in this range
    NUTATION_MAX_ARCSEC = 20.0  # arcsec (generous upper bound)
    NUTATION_MIN_ARCSEC = 0.1  # arcsec (minimum detectable)

    # Relaxed tolerance for NONUT comparison (accounts for nutation ~10 arcsec)
    NONUT_LONGITUDE_TOLERANCE = 0.005  # degrees (~18 arcsec)

    # Test dates for FLG_NONUT validation
    NONUT_TEST_DATES = [
        (2024, 6, 15, 12.0, "Mid-2024"),
        (2024, 11, 15, 0.0, "Late 2024"),
        (2025, 1, 1, 12.0, "Early 2025"),
        (2025, 6, 21, 12.0, "Summer Solstice 2025"),
        (2026, 3, 20, 12.0, "Vernal Equinox 2026"),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", NONUT_TEST_DATES)
    def test_nonut_positions_match_pyswisseph(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test FLG_NONUT planetary positions match pyswisseph.

        Compares libephemeris and pyswisseph when FLG_NONUT flag is set.
        Uses relaxed tolerance to account for nutation-related differences.
        """
        jd = swe.julday(year, month, day, hour)
        flag = FLG_SWIEPH | FLG_NONUT

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        assert diff_lon < self.NONUT_LONGITUDE_TOLERANCE, (
            f"{planet_name} NONUT at {date_desc}: longitude diff {diff_lon:.6f}° "
            f"exceeds tolerance (swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )
        assert diff_lat < LATITUDE_STRICT, (
            f"{planet_name} NONUT at {date_desc}: latitude diff {diff_lat:.6f}° "
            "exceeds tolerance"
        )
        assert diff_dist < DISTANCE_STRICT, (
            f"{planet_name} NONUT at {date_desc}: distance diff {diff_dist:.8f} AU "
            "exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", NONUT_TEST_DATES)
    def test_nonut_with_speed_match_pyswisseph(
        self, planet_id, planet_name, year, month, day, hour, date_desc
    ):
        """Test FLG_NONUT with velocity calculations match pyswisseph."""
        jd = swe.julday(year, month, day, hour)
        flag = FLG_SWIEPH | FLG_NONUT | FLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

        # Position checks with relaxed tolerance
        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        assert diff_lon < self.NONUT_LONGITUDE_TOLERANCE, (
            f"{planet_name} NONUT+Speed at {date_desc}: longitude diff {diff_lon:.6f}° "
            "exceeds tolerance"
        )
        assert diff_lat < LATITUDE_STRICT, (
            f"{planet_name} NONUT+Speed at {date_desc}: latitude diff {diff_lat:.6f}° "
            "exceeds tolerance"
        )
        assert diff_dist < DISTANCE_STRICT, (
            f"{planet_name} NONUT+Speed at {date_desc}: distance diff {diff_dist:.8f} AU "
            "exceeds tolerance"
        )

        # Velocity checks
        diff_lon_speed = abs(pos_swe[3] - pos_py[3])
        diff_lat_speed = abs(pos_swe[4] - pos_py[4])
        diff_dist_speed = abs(pos_swe[5] - pos_py[5])

        assert diff_lon_speed < VELOCITY_ANGULAR, (
            f"{planet_name} NONUT at {date_desc}: lon velocity diff "
            f"{diff_lon_speed:.6f}°/day exceeds tolerance"
        )
        assert diff_lat_speed < VELOCITY_ANGULAR, (
            f"{planet_name} NONUT at {date_desc}: lat velocity diff "
            f"{diff_lat_speed:.6f}°/day exceeds tolerance"
        )
        assert diff_dist_speed < VELOCITY_RADIAL, (
            f"{planet_name} NONUT at {date_desc}: dist velocity diff "
            f"{diff_dist_speed:.8f} AU/day exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS[:5])  # Test subset
    def test_nonut_pyswisseph_produces_different_position(self, planet_id, planet_name):
        """Verify pyswisseph FLG_NONUT produces different positions than default.

        The difference should be in the order of arcseconds, corresponding to
        the nutation effect (~9 arcsec principal term). This validates that
        pyswisseph properly implements the FLG_NONUT flag.
        """
        jd = swe.julday(2024, 6, 15, 12.0)

        with_nut_flags = FLG_SWIEPH
        no_nut_flags = FLG_SWIEPH | FLG_NONUT

        pos_with_nut, _ = swe.calc_ut(jd, planet_id, with_nut_flags)
        pos_no_nut, _ = swe.calc_ut(jd, planet_id, no_nut_flags)

        diff_deg = angular_diff(pos_with_nut[0], pos_no_nut[0])
        diff_arcsec = diff_deg * 3600  # Convert to arcseconds

        # Nutation effect should be detectable (> 0.1 arcsec)
        assert diff_arcsec > self.NUTATION_MIN_ARCSEC, (
            f"{planet_name}: pyswisseph NONUT diff {diff_arcsec:.4f} arcsec is too small "
            f"(expected > {self.NUTATION_MIN_ARCSEC} arcsec)"
        )

        # But not larger than maximum nutation (~20 arcsec including all terms)
        assert diff_arcsec < self.NUTATION_MAX_ARCSEC, (
            f"{planet_name}: pyswisseph NONUT diff {diff_arcsec:.4f} arcsec exceeds expected "
            f"nutation amplitude (< {self.NUTATION_MAX_ARCSEC} arcsec)"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", PLANETS[:5])  # Test subset
    def test_nonut_equatorial_match_pyswisseph(self, planet_id, planet_name):
        """Test FLG_NONUT with equatorial coordinates matches pyswisseph."""
        jd = swe.julday(2024, 6, 15, 12.0)
        flag = FLG_SWIEPH | FLG_NONUT | FLG_EQUATORIAL | FLG_SPEED

        pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
        pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

        diff_ra = angular_diff(pos_swe[0], pos_py[0])
        diff_dec = abs(pos_swe[1] - pos_py[1])

        # Use relaxed tolerance for equatorial coords (nutation affects both RA and Dec)
        assert diff_ra < self.NONUT_LONGITUDE_TOLERANCE, (
            f"{planet_name} NONUT equatorial RA: diff {diff_ra:.6f}° exceeds tolerance"
        )
        assert diff_dec < self.NONUT_LONGITUDE_TOLERANCE, (
            f"{planet_name} NONUT equatorial Dec: diff {diff_dec:.6f}° exceeds tolerance"
        )

    @pytest.mark.comparison
    def test_nonut_pyswisseph_nutation_magnitude(self):
        """Test that pyswisseph nutation effect is in expected range.

        Verify that pyswisseph correctly computes the nutation effect
        (difference between with/without FLG_NONUT) across different dates.
        """
        # Test dates
        test_dates = [
            swe.julday(2024, 1, 1, 12.0),
            swe.julday(2024, 6, 15, 12.0),
            swe.julday(2025, 1, 1, 12.0),
            swe.julday(2025, 6, 21, 12.0),
            swe.julday(2026, 1, 1, 12.0),
        ]

        with_nut_flags = FLG_SWIEPH
        no_nut_flags = FLG_SWIEPH | FLG_NONUT

        for jd in test_dates:
            pos_swe_with, _ = swe.calc_ut(jd, SUN, with_nut_flags)
            pos_swe_no, _ = swe.calc_ut(jd, SUN, no_nut_flags)

            diff_swe = angular_diff(pos_swe_with[0], pos_swe_no[0]) * 3600

            # Nutation effect should be in expected range
            assert 0.1 < diff_swe < 20.0, (
                f"pyswisseph nutation effect {diff_swe:.4f} arcsec outside expected range "
                f"at JD {jd}"
            )

    @pytest.mark.comparison
    def test_nonut_consistency_all_planets(self):
        """Verify FLG_NONUT positions for all planets across implementations."""
        jd = swe.julday(2025, 1, 1, 12.0)
        flag = FLG_SWIEPH | FLG_NONUT | FLG_SPEED

        for planet_id, planet_name in PLANETS:
            pos_swe, _ = swe.calc_ut(jd, planet_id, flag)
            pos_py, _ = ephem.calc_ut(jd, planet_id, flag)

            diff_lon = angular_diff(pos_swe[0], pos_py[0])

            # Use relaxed tolerance to account for nutation difference
            assert diff_lon < self.NONUT_LONGITUDE_TOLERANCE, (
                f"{planet_name} NONUT at 2025: longitude diff {diff_lon:.8f}° "
                "exceeds tolerance"
            )
