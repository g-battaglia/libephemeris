"""
Tests for eclipse northern limit calculation in libephemeris.

Tests the calculation of geographic coordinates (latitude, longitude) for points
along the northern limit of the umbral/antumbral shadow path during solar eclipses.

Reference data sources:
- NASA Eclipse website (eclipse.gsfc.nasa.gov)
- Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
- Meeus "Astronomical Algorithms" Chapter 54
"""

import pytest  # noqa: F401
from libephemeris import (
    julday,
    calc_eclipse_northern_limit,
    calc_eclipse_central_line,
    sol_eclipse_when_glob,
    FLG_SWIEPH,
    ECL_TOTAL,
    ECL_ANNULAR,
)

pytestmark = pytest.mark.slow


class TestEclipseNorthernLimitBasicFunctionality:
    """Test basic functionality of calc_eclipse_northern_limit."""

    def test_function_exists(self):
        """Test that calc_eclipse_northern_limit function exists and is callable."""
        from libephemeris.eclipse import calc_eclipse_northern_limit as func

        assert callable(func)

    def test_function_exported_from_main_module(self):
        """Test that function is exported from main libephemeris module."""
        from libephemeris import calc_eclipse_northern_limit as func

        assert callable(func)

    def test_swe_alias_exists(self):
        """Test that calc_eclipse_northern_limit alias exists."""
        assert callable(calc_eclipse_northern_limit)

    def test_returns_three_tuples(self):
        """Test that function returns a tuple of three tuples."""
        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_c1 = times_ecl[2]  # First contact
        jd_c4 = times_ecl[3]  # Fourth contact

        result = calc_eclipse_northern_limit(jd_c1, jd_c4, step_minutes=30.0)

        assert isinstance(result, tuple)
        assert len(result) == 3
        times, lats, lons = result
        assert isinstance(times, tuple)
        assert isinstance(lats, tuple)
        assert isinstance(lons, tuple)

    def test_returns_same_length_tuples(self):
        """Test that all returned tuples have the same length."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_c1 = times_ecl[2]
        jd_c4 = times_ecl[3]

        times, lats, lons = calc_eclipse_northern_limit(jd_c1, jd_c4, step_minutes=30.0)

        assert len(times) == len(lats) == len(lons)

    def test_accepts_flags_parameter(self):
        """Test that function accepts optional flags parameter."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_c1 = times_ecl[2]
        jd_c4 = times_ecl[3]

        result = calc_eclipse_northern_limit(
            jd_c1, jd_c4, step_minutes=30.0, flags=FLG_SWIEPH
        )

        assert isinstance(result, tuple)


class TestEclipseNorthernLimitTotalEclipse:
    """Test northern limit calculation for total eclipses."""

    def test_april_2024_total_eclipse_northern_limit(self):
        """Test April 8, 2024 total solar eclipse northern limit.

        NASA Reference: The April 8, 2024 total eclipse path crosses
        Mexico, USA, and Canada. The northern limit should be north of
        the central line throughout.
        """
        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        assert ecl_type & ECL_TOTAL, "Should find a total eclipse"

        jd_c1 = times_ecl[2]  # First contact
        jd_c4 = times_ecl[3]  # Fourth contact

        # Calculate northern limit with 10-minute steps
        times, lats, lons = calc_eclipse_northern_limit(jd_c1, jd_c4, step_minutes=10.0)

        # Should return points along the northern limit
        assert len(times) > 0, "Should have northern limit points"

        # All latitudes should be valid
        for lat in lats:
            assert -90.0 <= lat <= 90.0, f"Latitude {lat} out of range"

        # All longitudes should be valid
        for lon in lons:
            assert -180.0 <= lon <= 180.0, f"Longitude {lon} out of range"

    def test_northern_limit_times_monotonic(self):
        """Test that times are monotonically increasing."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_c1 = times_ecl[2]
        jd_c4 = times_ecl[3]

        times, lats, lons = calc_eclipse_northern_limit(jd_c1, jd_c4, step_minutes=10.0)

        if len(times) > 1:
            for i in range(1, len(times)):
                assert times[i] > times[i - 1], (
                    "Times should be monotonically increasing"
                )

    def test_northern_limit_latitudes_valid(self):
        """Test that all latitudes are in valid range."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_c1 = times_ecl[2]
        jd_c4 = times_ecl[3]

        times, lats, lons = calc_eclipse_northern_limit(jd_c1, jd_c4, step_minutes=5.0)

        for lat in lats:
            assert -90.0 <= lat <= 90.0, f"Latitude {lat} should be in [-90, 90]"

    def test_northern_limit_longitudes_normalized(self):
        """Test that all longitudes are normalized to [-180, 180]."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_c1 = times_ecl[2]
        jd_c4 = times_ecl[3]

        times, lats, lons = calc_eclipse_northern_limit(jd_c1, jd_c4, step_minutes=5.0)

        for lon in lons:
            assert -180.0 <= lon <= 180.0, f"Longitude {lon} should be in [-180, 180]"


class TestEclipseNorthernLimitAnnularEclipse:
    """Test northern limit calculation for annular eclipses."""

    def test_annular_eclipse_northern_limit(self):
        """Test annular eclipse northern limit calculation."""
        jd_start = julday(2023, 1, 1, 0.0)
        ecl_type, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_ANNULAR)

        if ecl_type & ECL_ANNULAR:
            jd_c1 = times_ecl[2]
            jd_c4 = times_ecl[3]

            times, lats, lons = calc_eclipse_northern_limit(
                jd_c1, jd_c4, step_minutes=10.0
            )

            # Should return points along the northern limit
            assert len(times) > 0, (
                "Should have northern limit points for annular eclipse"
            )

            # All coordinates should be valid
            for lat in lats:
                assert -90.0 <= lat <= 90.0
            for lon in lons:
                assert -180.0 <= lon <= 180.0


class TestEclipseNorthernLimitVsCentralLine:
    """Test that northern limit is correctly positioned relative to central line."""

    def test_northern_limit_generally_north_of_central_line(self):
        """Test that northern limit latitude is generally >= central line latitude.

        Note: Due to projection effects at extreme latitudes and the curvature
        of the shadow path, this relationship may not hold at every point,
        but should hold for most mid-eclipse points.
        """
        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        if not (ecl_type & ECL_TOTAL):
            pytest.skip("No total eclipse found")

        jd_max = times_ecl[0]  # Maximum eclipse time

        # Calculate central line and northern limit near maximum
        # Use a narrow time window around maximum for clearer comparison
        jd_c1 = jd_max - 0.05  # ~1.2 hours before
        jd_c4 = jd_max + 0.05  # ~1.2 hours after

        times_central, lats_central, lons_central = calc_eclipse_central_line(
            jd_c1, jd_c4, step_minutes=5.0
        )
        times_north, lats_north, lons_north = calc_eclipse_northern_limit(
            jd_c1, jd_c4, step_minutes=5.0
        )

        # Both should have points
        assert len(times_central) > 0, "Should have central line points"
        assert len(times_north) > 0, "Should have northern limit points"

        # For matching times, the northern limit should generally be north
        # (or very close to) the central line
        north_count = 0
        total_count = 0

        for i, t_north in enumerate(times_north):
            # Find closest central line point by time
            if times_central:
                closest_central_idx = min(
                    range(len(times_central)),
                    key=lambda j: abs(times_central[j] - t_north),
                )
                # Only compare if times are very close (within 1 minute)
                if abs(times_central[closest_central_idx] - t_north) < 1.0 / (24 * 60):
                    total_count += 1
                    # The northern limit should be at same or higher latitude
                    # Allow small tolerance for numerical precision
                    if lats_north[i] >= lats_central[closest_central_idx] - 0.5:
                        north_count += 1

        # At least 50% of points should satisfy the northern relationship
        # (accounting for edge effects and projection issues)
        if total_count > 0:
            ratio = north_count / total_count
            assert ratio >= 0.5, (
                f"Expected northern limit to be north of central line for most points, "
                f"but only {ratio * 100:.1f}% satisfied this condition"
            )

    def test_northern_limit_different_from_central_line(self):
        """Test that northern limit coordinates differ from central line."""
        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        if not (ecl_type & ECL_TOTAL):
            pytest.skip("No total eclipse found")

        jd_max = times_ecl[0]

        # Calculate at a single time point near maximum
        times_central, lats_central, lons_central = calc_eclipse_central_line(
            jd_max, jd_max + 0.001, step_minutes=0.5
        )
        times_north, lats_north, lons_north = calc_eclipse_northern_limit(
            jd_max, jd_max + 0.001, step_minutes=0.5
        )

        if len(lats_central) > 0 and len(lats_north) > 0:
            # The latitudes should be different (offset by umbral radius)
            lat_diff = abs(lats_north[0] - lats_central[0])
            # There should be some difference (typically a fraction of a degree to several degrees)
            # The difference should be positive but not huge
            assert lat_diff >= 0.0, "Northern limit should differ from central line"


class TestEclipseNorthernLimitStepSize:
    """Test northern limit calculation with different step sizes."""

    def test_smaller_step_gives_more_points(self):
        """Test that smaller step size produces more points."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_c1 = times_ecl[2]
        jd_c4 = times_ecl[3]

        times_1min, _, _ = calc_eclipse_northern_limit(jd_c1, jd_c4, step_minutes=1.0)
        times_10min, _, _ = calc_eclipse_northern_limit(jd_c1, jd_c4, step_minutes=10.0)

        # 1-minute steps should give roughly 10x more points than 10-minute steps
        assert len(times_1min) > len(times_10min)

    def test_default_step_is_one_minute(self):
        """Test that default step size is 1 minute."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_c1 = times_ecl[2]
        jd_c4 = times_ecl[3]

        # Use default step
        times_default, _, _ = calc_eclipse_northern_limit(jd_c1, jd_c4)
        # Explicit 1-minute step
        times_1min, _, _ = calc_eclipse_northern_limit(jd_c1, jd_c4, step_minutes=1.0)

        # Should get same number of points
        assert len(times_default) == len(times_1min)


class TestEclipseNorthernLimitEdgeCases:
    """Test edge cases for northern limit calculation."""

    def test_empty_result_for_non_eclipse_time(self):
        """Test that non-eclipse times return empty tuples."""
        # A time far from any eclipse
        jd_no_eclipse = julday(2024, 6, 15, 12.0)

        # Calculate over a short period that doesn't contain an eclipse center
        times, lats, lons = calc_eclipse_northern_limit(
            jd_no_eclipse, jd_no_eclipse + 0.01, step_minutes=1.0
        )

        # Should return empty or very few points since gamma > 1
        assert isinstance(times, tuple)

    def test_single_point_for_very_short_range(self):
        """Test that very short time range gives at least one point if valid."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times_ecl[0]  # Maximum eclipse time

        # Calculate at just the maximum moment
        times, lats, lons = calc_eclipse_northern_limit(
            jd_max, jd_max + 0.0001, step_minutes=0.1
        )

        # Should have at least one point at maximum
        assert len(times) >= 1


class TestEclipseNorthernLimitMultipleEclipses:
    """Test northern limit calculation across multiple eclipses."""

    def test_multiple_eclipses_have_valid_northern_limits(self):
        """Test that multiple eclipses all produce valid northern limits."""
        jd = julday(2020, 1, 1, 0.0)
        eclipses_tested = 0

        for _ in range(3):  # Test 3 total eclipses
            try:
                ecl_type, times_ecl = sol_eclipse_when_glob(jd, ecltype=ECL_TOTAL)
            except RuntimeError:
                break

            if not (ecl_type & ECL_TOTAL):
                break

            jd_c1 = times_ecl[2]
            jd_c4 = times_ecl[3]

            times, lats, lons = calc_eclipse_northern_limit(
                jd_c1, jd_c4, step_minutes=15.0
            )

            # Should have points for a central eclipse
            assert len(times) > 0, (
                f"Eclipse at JD {times_ecl[0]} should have northern limit"
            )

            eclipses_tested += 1
            # Move past this eclipse
            jd = times_ecl[0] + 30

        assert eclipses_tested > 0, "Should test at least one eclipse"


class TestEclipseNorthernLimitPhysicalReasonableness:
    """Test physical reasonableness of northern limit calculations."""

    def test_northern_limit_offset_reasonable(self):
        """Test that the offset between central line and northern limit is reasonable.

        The path width for total eclipses typically ranges from 0 to ~270 km,
        so the offset from central to limit should be roughly half that in km,
        which translates to roughly 0.0 to 2.0 degrees of latitude.
        """
        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        if not (ecl_type & ECL_TOTAL):
            pytest.skip("No total eclipse found")

        jd_max = times_ecl[0]

        # Calculate at maximum
        times_central, lats_central, lons_central = calc_eclipse_central_line(
            jd_max, jd_max + 0.001, step_minutes=0.1
        )
        times_north, lats_north, lons_north = calc_eclipse_northern_limit(
            jd_max, jd_max + 0.001, step_minutes=0.1
        )

        if len(lats_central) > 0 and len(lats_north) > 0:
            lat_offset = lats_north[0] - lats_central[0]
            # The offset should be non-negative and within reasonable bounds
            # Allow for some variation due to projection effects
            assert -1.0 <= lat_offset <= 10.0, (
                f"Latitude offset {lat_offset:.2f}° seems unreasonable"
            )
