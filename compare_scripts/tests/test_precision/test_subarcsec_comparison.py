"""
Sub-arcsecond precision comparison tests against pyswisseph.

Tests planetary positions for Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn
at 1000+ random dates and verifies that the maximum difference is within
per-planet tolerances:
- Sun, Mercury, Venus: < 1 arcsecond (high precision)
- Mars: < 2 arcseconds (slightly relaxed due to perturbation differences)
- Moon: < 5 arcseconds (relaxed for lunar theory differences)
- Jupiter, Saturn: < 5 arcseconds (outer planets)

These tests ensure that libephemeris maintains high-precision compatibility
with the reference pyswisseph implementation.
"""

import pytest
import random
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
)


# =============================================================================
# CONSTANTS
# =============================================================================

# Tolerance in arcseconds
# Note: The Moon requires a relaxed tolerance due to differences in lunar theory
# implementations between libephemeris (Skyfield/JPL DE) and pyswisseph.
# Sun, Mercury, Venus: < 1 arcsecond (high precision)
# Mars: < 2 arcseconds (slightly relaxed due to perturbation differences)
# Moon: < 5 arcseconds (lunar theory differences)
# Jupiter, Saturn: < 5 arcseconds (outer planets)

SUN_TOLERANCE_ARCSEC = 1.0
MOON_TOLERANCE_ARCSEC = 5.0  # Relaxed for lunar theory differences
MERCURY_TOLERANCE_ARCSEC = 1.0
VENUS_TOLERANCE_ARCSEC = 1.0
MARS_TOLERANCE_ARCSEC = 2.0  # Slightly relaxed for Mars perturbations
JUPITER_TOLERANCE_ARCSEC = 5.0
SATURN_TOLERANCE_ARCSEC = 5.0

# Number of random dates for massive comparison tests
NUM_RANDOM_DATES = 1000

# DE421 valid range: 1900-2050
DE421_START_YEAR = 1900
DE421_END_YEAR = 2050

# Planet definitions with individual tolerances
PLANET_TOLERANCES = {
    SUN: SUN_TOLERANCE_ARCSEC,
    MOON: MOON_TOLERANCE_ARCSEC,
    MERCURY: MERCURY_TOLERANCE_ARCSEC,
    VENUS: VENUS_TOLERANCE_ARCSEC,
    MARS: MARS_TOLERANCE_ARCSEC,
    JUPITER: JUPITER_TOLERANCE_ARCSEC,
    SATURN: SATURN_TOLERANCE_ARCSEC,
}

INNER_PLANETS = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
]

OUTER_PLANETS = [
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
]


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def angle_diff(a1: float, a2: float) -> float:
    """
    Calculate the smallest difference between two angles.

    Handles wrap-around at 360 degrees. For example:
    - angle_diff(359.0, 1.0) returns 2.0, not 358.0
    - angle_diff(0.5, 359.5) returns 1.0, not 359.0

    Args:
        a1: First angle in degrees
        a2: Second angle in degrees

    Returns:
        Absolute difference in degrees (0 to 180)
    """
    diff = abs(a1 - a2)
    if diff > 180:
        diff = 360 - diff
    return diff


def generate_random_dates(n: int, seed: int = 42) -> list:
    """
    Generate n random dates within DE421 valid range.

    Args:
        n: Number of dates to generate
        seed: Random seed for reproducibility

    Returns:
        List of tuples: (year, month, day, hour, jd)
    """
    random.seed(seed)
    dates = []
    for _ in range(n):
        year = random.randint(DE421_START_YEAR, DE421_END_YEAR)
        month = random.randint(1, 12)
        day = random.randint(1, 28)  # Safe for all months
        hour = random.uniform(0, 24)
        jd = ephem.julday(year, month, day, hour)
        dates.append((year, month, day, hour, jd))
    return dates


def calculate_planet_position(jd: float, planet_id: int) -> tuple:
    """
    Calculate planet position using both libraries.

    Args:
        jd: Julian Day
        planet_id: Planet identifier

    Returns:
        Tuple of (lib_longitude, swe_longitude, diff_arcsec)
    """
    pos_lib, _ = ephem.calc_ut(jd, planet_id, 0)
    pos_swe, _ = swe.calc_ut(jd, planet_id, 0)

    lon_diff = angle_diff(pos_lib[0], pos_swe[0])
    diff_arcsec = lon_diff * 3600  # Convert to arcseconds

    return pos_lib[0], pos_swe[0], diff_arcsec


# =============================================================================
# TEST CLASSES
# =============================================================================


class TestSubArcsecondPrecisionInnerPlanets:
    """
    Test sub-arcsecond precision for inner planets (Sun, Moon, Mercury, Venus, Mars).

    These planets are fast-moving and require high precision.
    Tolerance is per-planet based on PLANET_TOLERANCES.
    """

    @pytest.mark.precision
    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", INNER_PLANETS)
    def test_inner_planet_1000_random_dates(
        self, planet_id, planet_name, progress_reporter
    ):
        """
        Test inner planet positions at 1000+ random dates.

        Verifies that the maximum longitude difference is within tolerance.
        """
        dates = generate_random_dates(NUM_RANDOM_DATES)
        tolerance_arcsec = PLANET_TOLERANCES[planet_id]

        max_diff_arcsec = 0.0
        worst_case = None
        progress = progress_reporter(
            f"{planet_name} precision", len(dates), report_every=10
        )

        for i, (year, month, day, hour, jd) in enumerate(dates):
            lon_lib, lon_swe, diff_arcsec = calculate_planet_position(jd, planet_id)

            if diff_arcsec > max_diff_arcsec:
                max_diff_arcsec = diff_arcsec
                worst_case = {
                    "jd": jd,
                    "date": f"{year}-{month:02d}-{day:02d} {hour:.2f}h",
                    "lib": lon_lib,
                    "swe": lon_swe,
                    "diff_arcsec": diff_arcsec,
                }

            progress.update(i)

        progress.done(f"max diff: {max_diff_arcsec:.4f} arcsec")

        # Final assertion
        assert max_diff_arcsec < tolerance_arcsec, (
            f"{planet_name} max longitude difference {max_diff_arcsec:.4f} arcsec "
            f">= {tolerance_arcsec} arcsec at {worst_case}"
        )


class TestSubArcsecondPrecisionOuterPlanets:
    """
    Test sub-arcsecond precision for outer planets (Jupiter, Saturn).

    These planets are slower-moving but still require good precision.
    Maximum allowed difference: < 5 arcseconds.
    """

    @pytest.mark.precision
    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", OUTER_PLANETS)
    def test_outer_planet_1000_random_dates(
        self, planet_id, planet_name, progress_reporter
    ):
        """
        Test outer planet positions at 1000+ random dates.

        Verifies that the maximum longitude difference is < 5 arcseconds.
        """
        dates = generate_random_dates(NUM_RANDOM_DATES)
        tolerance_arcsec = PLANET_TOLERANCES[planet_id]

        max_diff_arcsec = 0.0
        worst_case = None
        progress = progress_reporter(
            f"{planet_name} precision", len(dates), report_every=10
        )

        for i, (year, month, day, hour, jd) in enumerate(dates):
            lon_lib, lon_swe, diff_arcsec = calculate_planet_position(jd, planet_id)

            if diff_arcsec > max_diff_arcsec:
                max_diff_arcsec = diff_arcsec
                worst_case = {
                    "jd": jd,
                    "date": f"{year}-{month:02d}-{day:02d} {hour:.2f}h",
                    "lib": lon_lib,
                    "swe": lon_swe,
                    "diff_arcsec": diff_arcsec,
                }

            progress.update(i)

        progress.done(f"max diff: {max_diff_arcsec:.4f} arcsec")

        # Final assertion
        assert max_diff_arcsec < tolerance_arcsec, (
            f"{planet_name} max longitude difference {max_diff_arcsec:.4f} arcsec "
            f">= {tolerance_arcsec} arcsec at {worst_case}"
        )


class TestAllTargetPlanetsCombined:
    """
    Combined test for all 7 target planets at 1000+ dates.

    Tests Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn simultaneously
    and reports the maximum difference for each planet.
    """

    @pytest.mark.precision
    @pytest.mark.comparison
    @pytest.mark.slow
    def test_all_planets_1000_dates_combined(self, progress_reporter):
        """
        Test all 7 planets at 1000+ random dates in a single test.

        Reports maximum differences for each planet and asserts all are within tolerance.
        """
        dates = generate_random_dates(NUM_RANDOM_DATES)

        # All planets with their tolerances
        all_planets = [
            (SUN, "Sun", PLANET_TOLERANCES[SUN]),
            (MOON, "Moon", PLANET_TOLERANCES[MOON]),
            (MERCURY, "Mercury", PLANET_TOLERANCES[MERCURY]),
            (VENUS, "Venus", PLANET_TOLERANCES[VENUS]),
            (MARS, "Mars", PLANET_TOLERANCES[MARS]),
            (JUPITER, "Jupiter", PLANET_TOLERANCES[JUPITER]),
            (SATURN, "Saturn", PLANET_TOLERANCES[SATURN]),
        ]

        # Track maximum differences per planet
        max_diffs = {name: 0.0 for _, name, _ in all_planets}
        worst_cases = {name: None for _, name, _ in all_planets}

        total_iterations = len(dates) * len(all_planets)
        progress = progress_reporter("All planets", total_iterations, report_every=5)

        iteration = 0
        for year, month, day, hour, jd in dates:
            for planet_id, planet_name, _ in all_planets:
                lon_lib, lon_swe, diff_arcsec = calculate_planet_position(jd, planet_id)

                if diff_arcsec > max_diffs[planet_name]:
                    max_diffs[planet_name] = diff_arcsec
                    worst_cases[planet_name] = {
                        "jd": jd,
                        "date": f"{year}-{month:02d}-{day:02d} {hour:.2f}h",
                        "lib": lon_lib,
                        "swe": lon_swe,
                        "diff_arcsec": diff_arcsec,
                    }

                progress.update(iteration)
                iteration += 1

        # Report results
        print("\n  === Sub-arcsecond Precision Results ===")
        for planet_id, planet_name, tolerance in all_planets:
            status = "PASS" if max_diffs[planet_name] < tolerance else "FAIL"
            print(
                f"  {planet_name:10s}: max diff = {max_diffs[planet_name]:8.4f} arcsec "
                f"(limit: {tolerance:.1f} arcsec) [{status}]"
            )
        print("  ========================================")

        progress.done("all planets tested")

        # Assert all planets are within tolerance
        failures = []
        for planet_id, planet_name, tolerance in all_planets:
            if max_diffs[planet_name] >= tolerance:
                failures.append(
                    f"{planet_name}: {max_diffs[planet_name]:.4f} >= {tolerance} arcsec"
                )

        assert not failures, f"Planets exceeding tolerance: {failures}"


class TestStatisticalAnalysis:
    """
    Statistical analysis of precision across all target planets.

    Provides mean, median, and percentile statistics for the differences.
    """

    @pytest.mark.precision
    @pytest.mark.comparison
    @pytest.mark.slow
    def test_statistical_precision_analysis(self, progress_reporter):
        """
        Perform statistical analysis of longitude differences.

        Calculates mean, median, 95th percentile, and max for each planet.
        """
        dates = generate_random_dates(NUM_RANDOM_DATES)

        all_planets = [
            (SUN, "Sun", PLANET_TOLERANCES[SUN]),
            (MOON, "Moon", PLANET_TOLERANCES[MOON]),
            (MERCURY, "Mercury", PLANET_TOLERANCES[MERCURY]),
            (VENUS, "Venus", PLANET_TOLERANCES[VENUS]),
            (MARS, "Mars", PLANET_TOLERANCES[MARS]),
            (JUPITER, "Jupiter", PLANET_TOLERANCES[JUPITER]),
            (SATURN, "Saturn", PLANET_TOLERANCES[SATURN]),
        ]

        # Collect all differences per planet
        diffs_per_planet = {name: [] for _, name, _ in all_planets}

        total_iterations = len(dates) * len(all_planets)
        progress = progress_reporter(
            "Statistical analysis", total_iterations, report_every=5
        )

        iteration = 0
        for year, month, day, hour, jd in dates:
            for planet_id, planet_name, _ in all_planets:
                _, _, diff_arcsec = calculate_planet_position(jd, planet_id)
                diffs_per_planet[planet_name].append(diff_arcsec)
                progress.update(iteration)
                iteration += 1

        # Calculate statistics
        print("\n  === Statistical Precision Analysis ===")
        print(
            f"  {'Planet':10s} {'Mean':>10s} {'Median':>10s} {'P95':>10s} {'Max':>10s}"
        )
        print("  " + "-" * 50)

        for planet_id, planet_name, tolerance in all_planets:
            diffs = sorted(diffs_per_planet[planet_name])
            n = len(diffs)
            mean_diff = sum(diffs) / n
            median_diff = diffs[n // 2]
            p95_diff = diffs[int(n * 0.95)]
            max_diff = diffs[-1]

            print(
                f"  {planet_name:10s} {mean_diff:10.4f} {median_diff:10.4f} "
                f"{p95_diff:10.4f} {max_diff:10.4f}"
            )

        print("  " + "=" * 50)
        progress.done("statistical analysis complete")

        # Assert maximum is within tolerance
        for planet_id, planet_name, tolerance in all_planets:
            max_diff = max(diffs_per_planet[planet_name])
            assert max_diff < tolerance, (
                f"{planet_name} max diff {max_diff:.4f} arcsec >= {tolerance} arcsec"
            )


class TestEdgeCaseDates:
    """
    Test precision at edge case dates (ephemeris boundaries, solstices, etc.).
    """

    EDGE_CASE_DATES = [
        # Ephemeris boundary dates
        (1900, 1, 1, 0.0, "DE421 start"),
        (2050, 12, 31, 23.99, "DE421 end"),
        # Y2K dates
        (1999, 12, 31, 23.99, "Y2K eve"),
        (2000, 1, 1, 0.01, "Y2K"),
        # J2000 epoch
        (2000, 1, 1, 12.0, "J2000 epoch"),
        # Solstices and equinoxes
        (2020, 3, 20, 3.83, "Vernal equinox 2020"),
        (2020, 6, 20, 21.73, "Summer solstice 2020"),
        (2020, 9, 22, 13.31, "Autumn equinox 2020"),
        (2020, 12, 21, 10.02, "Winter solstice 2020"),
        # Leap year dates
        (2000, 2, 29, 12.0, "Leap day 2000"),
        (2024, 2, 29, 12.0, "Leap day 2024"),
    ]

    @pytest.mark.precision
    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", INNER_PLANETS + OUTER_PLANETS)
    @pytest.mark.parametrize(
        "year,month,day,hour,description",
        EDGE_CASE_DATES,
        ids=[d[4] for d in EDGE_CASE_DATES],
    )
    def test_edge_case_precision(
        self, planet_id, planet_name, year, month, day, hour, description
    ):
        """
        Test precision at edge case dates.
        """
        jd = ephem.julday(year, month, day, hour)
        _, _, diff_arcsec = calculate_planet_position(jd, planet_id)

        # Use per-planet tolerance
        tolerance = PLANET_TOLERANCES[planet_id]

        assert diff_arcsec < tolerance, (
            f"{planet_name} at {description}: diff {diff_arcsec:.4f} arcsec "
            f">= {tolerance} arcsec"
        )


class TestLatitudeAndDistancePrecision:
    """
    Test precision of latitude and distance in addition to longitude.
    """

    @pytest.mark.precision
    @pytest.mark.comparison
    def test_full_position_precision_100_dates(self, progress_reporter):
        """
        Test longitude, latitude, and distance precision at 100 dates.

        Latitude and distance use different tolerance thresholds:
        - Latitude: < 1 arcsecond (for all planets)
        - Distance: < 0.0001 AU
        """
        dates = generate_random_dates(100, seed=123)  # Different seed for variety

        all_planets = [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ]

        # Latitude tolerances (Moon needs relaxed tolerance for lunar theory)
        lat_tolerance_arcsec = {
            "Sun": 1.0,
            "Moon": 5.0,  # Relaxed for Moon latitude due to lunar theory differences
            "Mercury": 1.0,
            "Venus": 1.0,
            "Mars": 1.0,
            "Jupiter": 1.0,
            "Saturn": 1.0,
        }
        dist_tolerance_au = 0.0001  # AU

        max_lat_diff = {name: 0.0 for _, name in all_planets}
        max_dist_diff = {name: 0.0 for _, name in all_planets}

        total = len(dates) * len(all_planets)
        progress = progress_reporter("Full position check", total, report_every=20)

        iteration = 0
        for year, month, day, hour, jd in dates:
            for planet_id, planet_name in all_planets:
                pos_lib, _ = ephem.calc_ut(jd, planet_id, 0)
                pos_swe, _ = swe.calc_ut(jd, planet_id, 0)

                lat_diff = angle_diff(pos_lib[1], pos_swe[1]) * 3600  # arcsec
                dist_diff = abs(pos_lib[2] - pos_swe[2])  # AU

                max_lat_diff[planet_name] = max(max_lat_diff[planet_name], lat_diff)
                max_dist_diff[planet_name] = max(max_dist_diff[planet_name], dist_diff)

                progress.update(iteration)
                iteration += 1

        progress.done("position check complete")

        # Print results
        print("\n  === Latitude and Distance Precision ===")
        print(f"  {'Planet':10s} {'Max Lat (arcsec)':>18s} {'Max Dist (AU)':>15s}")
        print("  " + "-" * 45)
        for _, planet_name in all_planets:
            print(
                f"  {planet_name:10s} {max_lat_diff[planet_name]:18.4f} "
                f"{max_dist_diff[planet_name]:15.6f}"
            )
        print("  " + "=" * 45)

        # Assert tolerances
        for _, planet_name in all_planets:
            assert max_lat_diff[planet_name] < lat_tolerance_arcsec[planet_name], (
                f"{planet_name} latitude diff {max_lat_diff[planet_name]:.4f} arcsec "
                f">= {lat_tolerance_arcsec[planet_name]} arcsec"
            )
            assert max_dist_diff[planet_name] < dist_tolerance_au, (
                f"{planet_name} distance diff {max_dist_diff[planet_name]:.6f} AU "
                f">= {dist_tolerance_au} AU"
            )
