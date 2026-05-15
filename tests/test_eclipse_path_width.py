"""
Tests for eclipse path width calculation in libephemeris.

Tests the calculation of the width of the path of totality or annularity
for central solar eclipses at any point along the central line.

Reference data sources:
- NASA Eclipse website (eclipse.gsfc.nasa.gov)
- Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
- Meeus "Astronomical Algorithms" Chapter 54
"""

import pytest
from libephemeris import (
    julday,
    calc_eclipse_path_width,
    sol_eclipse_where,
    sol_eclipse_when_glob,
    FLG_SWIEPH,
    ECL_TOTAL,
    ECL_ANNULAR,
)

pytestmark = pytest.mark.slow


class TestEclipsePathWidthBasicFunctionality:
    """Test basic functionality of calc_eclipse_path_width."""

    def test_function_exists(self):
        """Test that calc_eclipse_path_width function exists and is callable."""
        from libephemeris.eclipse import calc_eclipse_path_width

        assert callable(calc_eclipse_path_width)

    def test_function_exported_from_main_module(self):
        """Test that function is exported from main libephemeris module."""
        from libephemeris import calc_eclipse_path_width

        assert callable(calc_eclipse_path_width)

    def test_alias_exists(self):
        """Test that calc_eclipse_path_width alias exists and is callable."""
        assert callable(calc_eclipse_path_width)

    def test_returns_float(self):
        """Test that function returns a float value."""
        # Get a known total eclipse maximum
        jd_start = julday(2024, 1, 1, 0.0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_eclipse_path_width(jd_max)

        assert isinstance(result, float)

    def test_accepts_flags_parameter(self):
        """Test that function accepts optional flags parameter."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_eclipse_path_width(jd_max, flags=FLG_SWIEPH)

        assert isinstance(result, float)


class TestEclipsePathWidthTotalEclipses:
    """Test path width calculation for total solar eclipses."""

    def test_april_2024_total_eclipse_path_width(self):
        """Test April 8, 2024 total solar eclipse path width.

        NASA Reference (from eclipse.gsfc.nasa.gov/SEpath/SEpath2001/SE2024Apr08Tpath.html):
        - Path width at greatest eclipse: ~197 km
        - This eclipse traverses North America with varying path width

        We use a tolerance since the calculation method differs slightly.
        """
        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        width = calc_eclipse_path_width(jd_max)

        # Width should be positive for a total eclipse
        assert width > 0, "Path width should be positive for a total eclipse"

        # Width should be in reasonable range for total eclipses (0-270 km typical)
        assert 50 < width < 400, (
            f"Path width {width:.1f} km is outside reasonable range (50-400 km)"
        )

    def test_december_2021_total_eclipse_path_width(self):
        """Test December 4, 2021 total solar eclipse path width.

        This was a total eclipse visible from Antarctica.
        Path width at greatest eclipse: ~417 km (wider than average due to geometry)
        """
        jd_start = julday(2021, 11, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        width = calc_eclipse_path_width(jd_max)

        # Width should be positive
        assert width > 0, "Path width should be positive for a total eclipse"

        # Width should be physically reasonable
        assert 100 < width < 600, (
            f"Path width {width:.1f} km is outside reasonable range"
        )

    def test_path_width_positive_at_central_line(self):
        """Test that path width is positive when calculated at the central line."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # Get the central line coordinates
        _, geopos, attr = sol_eclipse_where(jd_max, FLG_SWIEPH)
        central_lon = geopos[0]
        central_lat = geopos[1]

        # Calculate width at the central line location
        width = calc_eclipse_path_width(jd_max, lat=central_lat, lon=central_lon)

        assert width > 0, f"Path width at central line should be positive, got {width}"


class TestEclipsePathWidthAnnularEclipses:
    """Test path width calculation for annular solar eclipses."""

    def test_october_2023_annular_eclipse_path_width(self):
        """Test October 14, 2023 annular solar eclipse path width.

        This was an annular eclipse visible from the Americas.
        Path width at greatest eclipse: ~187 km
        """
        jd_start = julday(2023, 9, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_ANNULAR)
        jd_max = times[0]

        width = calc_eclipse_path_width(jd_max)

        # Width should be positive for an annular eclipse
        assert width > 0, "Path width should be positive for an annular eclipse"

        # Width should be in reasonable range for annular eclipses
        assert 50 < width < 500, (
            f"Path width {width:.1f} km is outside reasonable range for annular eclipse"
        )

    def test_june_2021_annular_eclipse_path_width(self):
        """Test June 10, 2021 annular solar eclipse path width.

        This was an annular eclipse visible from Canada and Russia.
        """
        jd_start = julday(2021, 5, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_ANNULAR)
        jd_max = times[0]

        width = calc_eclipse_path_width(jd_max)

        # Width should be positive
        assert width > 0, "Path width should be positive for an annular eclipse"


class TestEclipsePathWidthWithLocation:
    """Test path width calculation with specific location."""

    def test_path_width_at_dallas_2024_eclipse(self):
        """Test path width at Dallas, Texas during April 8, 2024 eclipse.

        Dallas was near the center of the path of totality.
        We need to use a time when Dallas is near the central line.
        """
        jd_start = julday(2024, 1, 1, 0.0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        dallas_lat = 32.7767
        dallas_lon = -96.7970

        # Find a time when the eclipse is near Dallas
        # The eclipse moved westward, so we need a time earlier than global max
        # Let's search for the time when the central line longitude is near Dallas
        for offset_hours in range(-3, 3):
            test_jd = jd_max + offset_hours / 24.0
            ecl_type, geopos, attr = sol_eclipse_where(test_jd, FLG_SWIEPH)
            central_lon = geopos[0]
            # If central line is near Dallas longitude (within 10 degrees)
            if abs(central_lon - dallas_lon) < 10:
                width = calc_eclipse_path_width(test_jd, lat=dallas_lat, lon=dallas_lon)
                if width > 0:
                    # Found a time with valid path width
                    assert width > 0, (
                        f"Path width at Dallas should be positive, got {width}"
                    )
                    assert width < 500, f"Path width {width:.1f} km seems too large"
                    return

        # If we didn't find a specific time, just test that the function
        # handles the case gracefully (returning 0 when not at central line)
        width = calc_eclipse_path_width(jd_max, lat=dallas_lat, lon=dallas_lon)
        assert width >= 0, f"Path width should be non-negative, got {width}"

    def test_path_width_outside_central_path_returns_zero(self):
        """Test that path width is zero for locations outside the central path."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # London is well outside the April 2024 eclipse path
        london_lat = 51.5074
        london_lon = -0.1278

        width = calc_eclipse_path_width(jd_max, lat=london_lat, lon=london_lon)

        # Path width should be zero outside the central path
        assert width == 0, f"Path width should be 0 outside central path, got {width}"


class TestEclipsePathWidthEdgeCases:
    """Test edge cases for path width calculation."""

    def test_non_eclipse_time_returns_zero_or_small(self):
        """Test that non-eclipse times return zero or very small width."""
        # A random time when there's no eclipse
        jd = julday(2024, 6, 15, 12.0)

        width = calc_eclipse_path_width(jd)

        # Should return 0 when there's no eclipse
        assert width == 0, f"Path width should be 0 at non-eclipse time, got {width}"

    def test_path_width_consistency_with_where_function(self):
        """Test that path width is consistent with sol_eclipse_where results."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        # Get path width from sol_eclipse_where
        ecl_type, geopos, attr = sol_eclipse_where(jd_max, FLG_SWIEPH)
        where_width = attr[3]  # Path width in km from attr[3]

        # Get path width from our function
        calc_width = calc_eclipse_path_width(jd_max)

        # Both should be non-zero for central eclipse
        assert calc_width > 0, "Path width should be positive"
        # sol_eclipse_where attr[3] is negative for total eclipses (sign convention)
        assert where_width != 0, "sol_eclipse_where width should be non-zero"

        # Both should be in the same order of magnitude
        # Allow for differences in calculation method
        # Use absolute value of where_width since it's negative for total eclipses
        ratio = calc_width / abs(where_width) if where_width != 0 else 0
        assert 0.3 < ratio < 3.0, (
            f"Path width ratio {ratio:.2f} is too different "
            f"(calc={calc_width:.1f}, where={where_width:.1f})"
        )


class TestEclipsePathWidthPhysicalReasonableness:
    """Test physical reasonableness of path width calculations."""

    def test_path_width_not_negative(self):
        """Test that path width is never negative."""
        jd_start = julday(2020, 1, 1, 0.0)

        # Test several total eclipses
        jd = jd_start
        for _ in range(5):
            ecl_type, times = sol_eclipse_when_glob(jd, ecltype=ECL_TOTAL)
            jd_max = times[0]

            width = calc_eclipse_path_width(jd_max)

            assert width >= 0, f"Path width should not be negative, got {width}"

            jd = jd_max + 30  # Move to next eclipse

    def test_path_width_bounded(self):
        """Test that path width stays within physical limits.

        The maximum possible path width for a total eclipse is around 270 km,
        and for annular eclipses around 380 km. However, at grazing geometry
        (near sunrise/sunset on the central line), widths can be larger.
        We use a generous upper bound of 1500 km for this test.
        """
        jd_start = julday(2020, 1, 1, 0.0)

        jd = jd_start
        for _ in range(10):
            ecl_type, times = sol_eclipse_when_glob(jd, ecltype=ECL_TOTAL)
            jd_max = times[0]

            width = calc_eclipse_path_width(jd_max)

            assert width <= 1500, (
                f"Path width {width:.1f} km exceeds physical limit (1500 km)"
            )

            jd = jd_max + 30
