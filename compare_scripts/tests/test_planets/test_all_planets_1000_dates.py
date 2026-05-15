"""
Test all planets (Sun through Pluto) at 1000 random dates.

This comprehensive test verifies that libephemeris matches pyswisseph
for all major planet positions across a large random sample of dates
within the DE421 ephemeris range (1900-2050).

Tests cover:
- Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
- 1000 random dates × 10 planets = 10,000 position comparisons
- Longitude, latitude, and distance accuracy
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
)


class TestAllPlanets1000Dates:
    """Test all planets at 1000 random dates against pyswisseph."""

    # Define all planets from Sun through Pluto
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

    @pytest.mark.comparison
    @pytest.mark.slow
    def test_all_planets_1000_dates(
        self, random_dates_in_de421_range, progress_reporter
    ):
        """
        Compare all planets at 1000 random dates with pyswisseph.

        This test generates 1000 random Julian Dates within the DE421
        ephemeris range (1900-2050) and compares the calculated positions
        for all 10 major planets (Sun through Pluto) between libephemeris
        and pyswisseph.

        Total comparisons: 10 planets × 1000 dates = 10,000
        """
        # Generate 1000 random dates
        dates = random_dates_in_de421_range(1000)

        # Tolerances for comparison
        # Note: Distance tolerance relaxed to 0.0002 AU to accommodate
        # outer planets (especially Pluto) which use different ephemeris
        # segments between DE440 (libephemeris) and DE431 (pyswisseph)
        lon_tolerance = 0.001  # degrees (~3.6 arcseconds)
        lat_tolerance = 0.001  # degrees
        dist_tolerance = 0.0002  # AU (relaxed for outer planets)

        total = len(dates) * len(self.PLANETS)
        progress = progress_reporter(
            "All planets × 1000 dates", total, report_every=100
        )

        # Track statistics
        max_lon_diff = 0.0
        max_lat_diff = 0.0
        max_dist_diff = 0.0
        worst_case = None

        iteration = 0
        for year, month, day, hour, jd in dates:
            for planet_id, planet_name in self.PLANETS:
                # Calculate with libephemeris
                pos_lib, _ = ephem.calc_ut(jd, planet_id, 0)

                # Calculate with pyswisseph
                pos_swe, _ = swe.calc_ut(jd, planet_id, 0)

                # Calculate differences
                lon_diff = abs(pos_lib[0] - pos_swe[0])
                if lon_diff > 180:
                    lon_diff = 360 - lon_diff

                lat_diff = abs(pos_lib[1] - pos_swe[1])
                dist_diff = abs(pos_lib[2] - pos_swe[2])

                # Track maximum differences
                if lon_diff > max_lon_diff:
                    max_lon_diff = lon_diff
                    worst_case = (planet_name, jd, "longitude", lon_diff)
                if lat_diff > max_lat_diff:
                    max_lat_diff = lat_diff
                if dist_diff > max_dist_diff:
                    max_dist_diff = dist_diff

                # Assert within tolerance
                assert lon_diff < lon_tolerance, (
                    f"{planet_name} at JD {jd} ({year}-{month:02d}-{day:02d}): "
                    f"longitude diff {lon_diff:.6f}° >= {lon_tolerance}° "
                    f"(lib: {pos_lib[0]:.6f}°, swe: {pos_swe[0]:.6f}°)"
                )
                assert lat_diff < lat_tolerance, (
                    f"{planet_name} at JD {jd} ({year}-{month:02d}-{day:02d}): "
                    f"latitude diff {lat_diff:.6f}° >= {lat_tolerance}° "
                    f"(lib: {pos_lib[1]:.6f}°, swe: {pos_swe[1]:.6f}°)"
                )
                assert dist_diff < dist_tolerance, (
                    f"{planet_name} at JD {jd} ({year}-{month:02d}-{day:02d}): "
                    f"distance diff {dist_diff:.8f} AU >= {dist_tolerance} AU "
                    f"(lib: {pos_lib[2]:.8f}, swe: {pos_swe[2]:.8f})"
                )

                progress.update(
                    iteration, f"{planet_name} @ {year}-{month:02d}-{day:02d}"
                )
                iteration += 1

        progress.done(
            f"max lon diff: {max_lon_diff:.6f}°, "
            f"max lat diff: {max_lat_diff:.6f}°, "
            f"max dist diff: {max_dist_diff:.8f} AU"
        )

    @pytest.mark.comparison
    @pytest.mark.slow
    def test_all_planets_1000_dates_with_speed(
        self, random_dates_in_de421_range, progress_reporter
    ):
        """
        Compare all planets with velocity at 1000 random dates.

        This test additionally validates the velocity (speed) calculations
        by comparing the longitude, latitude, and distance velocities
        between libephemeris and pyswisseph.
        """
        from libephemeris.constants import FLG_SPEED

        # Generate 1000 random dates
        dates = random_dates_in_de421_range(1000)

        # Tolerances
        lon_tolerance = 0.001  # degrees
        lat_tolerance = 0.001  # degrees
        dist_tolerance = 0.0002  # AU (relaxed for outer planets)
        vel_tolerance = 0.01  # degrees/day

        total = len(dates) * len(self.PLANETS)
        progress = progress_reporter(
            "All planets × 1000 dates (with speed)", total, report_every=100
        )

        iteration = 0
        for year, month, day, hour, jd in dates:
            for planet_id, planet_name in self.PLANETS:
                # Calculate with libephemeris (with speed flag)
                pos_lib, _ = ephem.calc_ut(jd, planet_id, FLG_SPEED)

                # Calculate with pyswisseph (with speed flag)
                pos_swe, _ = swe.calc_ut(jd, planet_id, FLG_SPEED)

                # Position differences
                lon_diff = abs(pos_lib[0] - pos_swe[0])
                if lon_diff > 180:
                    lon_diff = 360 - lon_diff
                lat_diff = abs(pos_lib[1] - pos_swe[1])
                dist_diff = abs(pos_lib[2] - pos_swe[2])

                # Velocity differences
                lon_vel_diff = abs(pos_lib[3] - pos_swe[3])
                lat_vel_diff = abs(pos_lib[4] - pos_swe[4])
                dist_vel_diff = abs(pos_lib[5] - pos_swe[5])

                # Assert positions within tolerance
                assert lon_diff < lon_tolerance, (
                    f"{planet_name} at JD {jd}: longitude diff {lon_diff:.6f}°"
                )
                assert lat_diff < lat_tolerance, (
                    f"{planet_name} at JD {jd}: latitude diff {lat_diff:.6f}°"
                )
                assert dist_diff < dist_tolerance, (
                    f"{planet_name} at JD {jd}: distance diff {dist_diff:.8f} AU"
                )

                # Assert velocities within tolerance
                assert lon_vel_diff < vel_tolerance, (
                    f"{planet_name} at JD {jd}: lon velocity diff {lon_vel_diff:.6f}°/day "
                    f"(lib: {pos_lib[3]:.6f}, swe: {pos_swe[3]:.6f})"
                )
                assert lat_vel_diff < vel_tolerance, (
                    f"{planet_name} at JD {jd}: lat velocity diff {lat_vel_diff:.6f}°/day"
                )
                assert dist_vel_diff < vel_tolerance, (
                    f"{planet_name} at JD {jd}: dist velocity diff {dist_vel_diff:.8f} AU/day"
                )

                progress.update(
                    iteration, f"{planet_name} @ {year}-{month:02d}-{day:02d}"
                )
                iteration += 1

        progress.done()


class TestAllPlanets1000DatesStatistics:
    """Statistical analysis of planet position accuracy across 1000 dates."""

    PLANETS = TestAllPlanets1000Dates.PLANETS

    @pytest.mark.comparison
    @pytest.mark.slow
    def test_position_accuracy_statistics(
        self, random_dates_in_de421_range, progress_reporter
    ):
        """
        Compute statistical summary of position accuracy.

        This test calculates the mean, max, and standard deviation of
        position differences for each planet across 1000 dates.
        """
        import statistics

        dates = random_dates_in_de421_range(1000)

        total = len(dates) * len(self.PLANETS)
        progress = progress_reporter("Statistics analysis", total, report_every=200)

        # Collect differences per planet
        planet_stats = {
            name: {"lon": [], "lat": [], "dist": []} for _, name in self.PLANETS
        }

        iteration = 0
        for year, month, day, hour, jd in dates:
            for planet_id, planet_name in self.PLANETS:
                pos_lib, _ = ephem.calc_ut(jd, planet_id, 0)
                pos_swe, _ = swe.calc_ut(jd, planet_id, 0)

                lon_diff = abs(pos_lib[0] - pos_swe[0])
                if lon_diff > 180:
                    lon_diff = 360 - lon_diff

                planet_stats[planet_name]["lon"].append(lon_diff)
                planet_stats[planet_name]["lat"].append(abs(pos_lib[1] - pos_swe[1]))
                planet_stats[planet_name]["dist"].append(abs(pos_lib[2] - pos_swe[2]))

                progress.update(iteration, f"{planet_name}")
                iteration += 1

        progress.done()

        # Print statistics summary
        print("\n" + "=" * 70)
        print("Position Accuracy Statistics (1000 dates × 10 planets)")
        print("=" * 70)
        print(
            f"{'Planet':<10} {'Mean Lon':<12} {'Max Lon':<12} {'Mean Lat':<12} {'Max Dist':<12}"
        )
        print("-" * 70)

        for planet_id, planet_name in self.PLANETS:
            stats = planet_stats[planet_name]
            mean_lon = statistics.mean(stats["lon"])
            max_lon = max(stats["lon"])
            mean_lat = statistics.mean(stats["lat"])
            max_dist = max(stats["dist"])

            print(
                f"{planet_name:<10} {mean_lon:.6f}°    {max_lon:.6f}°    "
                f"{mean_lat:.6f}°    {max_dist:.8f} AU"
            )

            # Verify all differences are within acceptable bounds
            assert max_lon < 0.001, (
                f"{planet_name} max longitude diff {max_lon} exceeds 0.001°"
            )
            assert max(stats["lat"]) < 0.001, (
                f"{planet_name} max latitude diff exceeds 0.001°"
            )
            # Distance tolerance relaxed for outer planets (Pluto uses different ephemeris)
            assert max_dist < 0.0002, (
                f"{planet_name} max distance diff {max_dist} exceeds 0.0002 AU"
            )

        print("=" * 70)
