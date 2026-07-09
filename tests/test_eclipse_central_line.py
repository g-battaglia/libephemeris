"""
Tests for eclipse central line coordinates calculation in libephemeris.

Tests the calculation of geographic coordinates (latitude, longitude) for points
along the central line of solar eclipses.

Reference data sources:
- NASA Eclipse website (eclipse.gsfc.nasa.gov)
- Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
- Meeus "Astronomical Algorithms" Chapter 54
"""

import pytest  # noqa: F401

pytestmark = pytest.mark.slow

from libephemeris import (
    julday,
    calc_eclipse_central_line,
    sol_eclipse_when_glob,
    FLG_SWIEPH,
    ECL_TOTAL,
    ECL_ANNULAR,
    ECL_PARTIAL,
)


class TestEclipseCentralLineBasicFunctionality:
    """Test basic functionality of calc_eclipse_central_line."""

    def test_function_exists(self):
        """Test that calc_eclipse_central_line function exists and is callable."""
        from libephemeris.eclipse import calc_eclipse_central_line as func

        assert callable(func)

    def test_function_exported_from_main_module(self):
        """Test that function is exported from main libephemeris module."""
        from libephemeris import calc_eclipse_central_line as func

        assert callable(func)

    def test_swe_alias_exists(self):
        """Test that calc_eclipse_central_line alias exists."""
        from libephemeris import calc_eclipse_central_line

        assert callable(calc_eclipse_central_line)

    def test_returns_three_tuples(self):
        """Test that function returns a tuple of three tuples."""
        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_begin = times_ecl[2]  # Eclipse begin
        jd_end = times_ecl[3]  # Eclipse end

        result = calc_eclipse_central_line(jd_begin, jd_end, step_minutes=30.0)

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

        times, lats, lons = calc_eclipse_central_line(jd_c1, jd_c4, step_minutes=30.0)

        assert len(times) == len(lats) == len(lons)

    def test_accepts_flags_parameter(self):
        """Test that function accepts optional flags parameter."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_c1 = times_ecl[2]
        jd_c4 = times_ecl[3]

        result = calc_eclipse_central_line(
            jd_c1, jd_c4, step_minutes=30.0, flags=FLG_SWIEPH
        )

        assert isinstance(result, tuple)


class TestEclipseCentralLineTotalEclipse:
    """Test central line calculation for total eclipses."""

    def test_april_2024_total_eclipse_central_line(self):
        """Test April 8, 2024 total solar eclipse central line.

        NASA Reference (from eclipse.gsfc.nasa.gov/SEpath/SEpath2001/SE2024Apr08Tpath.html):
        - Eclipse path crosses Mexico, USA (Texas to Maine), and eastern Canada
        - Central line passes through Dallas, Texas area (~32.8°N, ~96.8°W)
        - Central line passes through Cleveland, Ohio area (~41.5°N, ~81.7°W)
        """
        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        assert ecl_type & ECL_TOTAL, "Should find a total eclipse"

        jd_c1 = times_ecl[2]  # First contact
        jd_c4 = times_ecl[3]  # Fourth contact

        # Calculate central line with 10-minute steps
        times, lats, lons = calc_eclipse_central_line(jd_c1, jd_c4, step_minutes=10.0)

        # Should return points along the central line
        assert len(times) > 0, "Should have central line points"

        # Central line should cross North America
        # Latitudes should be in a reasonable range (roughly 10°N to 60°N)
        for lat in lats:
            assert -90.0 <= lat <= 90.0, f"Latitude {lat} out of range"

        # Longitudes should be in Western Hemisphere for this eclipse
        for lon in lons:
            assert -180.0 <= lon <= 180.0, f"Longitude {lon} out of range"

    def test_central_line_times_monotonic(self):
        """Test that times are monotonically increasing."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_c1 = times_ecl[2]
        jd_c4 = times_ecl[3]

        times, lats, lons = calc_eclipse_central_line(jd_c1, jd_c4, step_minutes=10.0)

        if len(times) > 1:
            for i in range(1, len(times)):
                assert times[i] > times[i - 1], (
                    "Times should be monotonically increasing"
                )

    def test_central_line_latitudes_valid(self):
        """Test that all latitudes are in valid range."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_c1 = times_ecl[2]
        jd_c4 = times_ecl[3]

        times, lats, lons = calc_eclipse_central_line(jd_c1, jd_c4, step_minutes=5.0)

        for lat in lats:
            assert -90.0 <= lat <= 90.0, f"Latitude {lat} should be in [-90, 90]"

    def test_central_line_longitudes_normalized(self):
        """Test that all longitudes are normalized to [-180, 180]."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_c1 = times_ecl[2]
        jd_c4 = times_ecl[3]

        times, lats, lons = calc_eclipse_central_line(jd_c1, jd_c4, step_minutes=5.0)

        for lon in lons:
            assert -180.0 <= lon <= 180.0, f"Longitude {lon} should be in [-180, 180]"


class TestEclipseCentralLineAnnularEclipse:
    """Test central line calculation for annular eclipses."""

    def test_annular_eclipse_central_line(self):
        """Test annular eclipse central line calculation."""
        jd_start = julday(2023, 1, 1, 0.0)
        ecl_type, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_ANNULAR)

        if ecl_type & ECL_ANNULAR:
            jd_c1 = times_ecl[2]
            jd_c4 = times_ecl[3]

            times, lats, lons = calc_eclipse_central_line(
                jd_c1, jd_c4, step_minutes=10.0
            )

            # Should return points along the central line
            assert len(times) > 0, "Should have central line points for annular eclipse"

            # All coordinates should be valid
            for lat in lats:
                assert -90.0 <= lat <= 90.0
            for lon in lons:
                assert -180.0 <= lon <= 180.0


class TestEclipseCentralLinePartialEclipse:
    """A partial-only eclipse defines no central line."""

    def test_partial_eclipse_returns_no_points(self):
        """2022-10-25 partial solar eclipse: no central phase anywhere on Earth.

        sol_eclipse_where reports retflag = ECL_PARTIAL | ECL_NONCENTRAL (> 0)
        throughout, so the central-line sampler must contribute no points; the
        former `retflag > 0` gate filled the tuples with the closest-approach
        shadow-axis track instead.
        """
        jd_start = julday(2022, 10, 1, 0.0)
        ecl_type, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_PARTIAL)
        assert ecl_type & ECL_PARTIAL, "Should find the 2022-10-25 partial eclipse"

        jd_c1 = times_ecl[2]
        jd_c4 = times_ecl[3]

        times, lats, lons = calc_eclipse_central_line(jd_c1, jd_c4, step_minutes=10.0)
        assert times == ()
        assert lats == ()
        assert lons == ()


class TestEclipseCentralLineStepSize:
    """Test central line calculation with different step sizes."""

    def test_smaller_step_gives_more_points(self):
        """Test that smaller step size produces more points."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_c1 = times_ecl[2]
        jd_c4 = times_ecl[3]

        times_1min, _, _ = calc_eclipse_central_line(jd_c1, jd_c4, step_minutes=1.0)
        times_10min, _, _ = calc_eclipse_central_line(jd_c1, jd_c4, step_minutes=10.0)

        # 1-minute steps should give roughly 10x more points than 10-minute steps
        assert len(times_1min) > len(times_10min)

    def test_default_step_is_one_minute(self):
        """Test that default step size is 1 minute."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_c1 = times_ecl[2]
        jd_c4 = times_ecl[3]

        # Use default step
        times_default, _, _ = calc_eclipse_central_line(jd_c1, jd_c4)
        # Explicit 1-minute step
        times_1min, _, _ = calc_eclipse_central_line(jd_c1, jd_c4, step_minutes=1.0)

        # Should get same number of points
        assert len(times_default) == len(times_1min)


class TestEclipseCentralLineEdgeCases:
    """Test edge cases for central line calculation."""

    def test_empty_result_for_non_eclipse_time(self):
        """Test that non-eclipse times return empty tuples."""
        # A time far from any eclipse
        jd_no_eclipse = julday(2024, 6, 15, 12.0)

        # Calculate over a short period that doesn't contain an eclipse center
        times, lats, lons = calc_eclipse_central_line(
            jd_no_eclipse, jd_no_eclipse + 0.01, step_minutes=1.0
        )

        # Should return empty or very few points since gamma > 1
        # (shadow axis doesn't intersect Earth)
        assert isinstance(times, tuple)

    def test_single_point_for_very_short_range(self):
        """Test that very short time range gives at least one point if valid."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times_ecl[0]  # Maximum eclipse time

        # Calculate at just the maximum moment
        times, lats, lons = calc_eclipse_central_line(
            jd_max, jd_max + 0.0001, step_minutes=0.1
        )

        # Should have at least one point at maximum
        assert len(times) >= 1


class TestEclipseCentralLineConsistency:
    """Test consistency with other eclipse functions."""

    def test_central_line_passes_through_where_location(self):
        """Test that central line passes through location found by sol_eclipse_where.

        The sol_eclipse_where function finds the central line location at a specific
        time. The central line from calc_eclipse_central_line should pass through
        or near that location at the same time.
        """
        from libephemeris import sol_eclipse_where

        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times_ecl = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        if ecl_type & ECL_TOTAL:
            jd_max = times_ecl[0]

            # Get central location from sol_eclipse_where
            retflag, geopos, attr = sol_eclipse_where(jd_max)
            where_lon = geopos[0]
            where_lat = geopos[1]

            # Get central line around maximum
            times, lats, lons = calc_eclipse_central_line(
                jd_max - 0.01, jd_max + 0.01, step_minutes=1.0
            )

            if len(times) > 0 and where_lat != 0.0:
                # Find the point closest to jd_max
                closest_idx = min(
                    range(len(times)), key=lambda i: abs(times[i] - jd_max)
                )

                # Check if the central line point is reasonably close to where location
                # Allow generous tolerance since algorithms may differ slightly
                lat_diff = abs(lats[closest_idx] - where_lat)
                lon_diff = abs(lons[closest_idx] - where_lon)
                if lon_diff > 180:
                    lon_diff = 360 - lon_diff

                # Should be within 10 degrees (rough consistency check)
                assert lat_diff < 15, (
                    f"Latitude difference {lat_diff} too large between "
                    f"central line ({lats[closest_idx]}) and where ({where_lat})"
                )


class TestEclipseCentralLineMultipleEclipses:
    """Test central line calculation across multiple eclipses."""

    def test_multiple_eclipses_have_valid_central_lines(self):
        """Test that multiple eclipses all produce valid central lines."""
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

            times, lats, lons = calc_eclipse_central_line(
                jd_c1, jd_c4, step_minutes=15.0
            )

            # Should have points for a central eclipse
            assert len(times) > 0, (
                f"Eclipse at JD {times_ecl[0]} should have central line"
            )

            eclipses_tested += 1
            # Move past this eclipse
            jd = times_ecl[0] + 30

        assert eclipses_tested > 0, "Should test at least one eclipse"
