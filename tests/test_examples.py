"""
Tests for the libephemeris example scripts.

These tests verify that all example scripts can be imported and their main
functions execute without errors, producing valid astronomical data.
"""

import subprocess
import sys
from pathlib import Path

import pytest


# Get the examples directory
EXAMPLES_DIR = Path(__file__).parent.parent / "examples"


class TestExampleImports:
    """Tests that all example modules can be imported."""

    def test_basic_usage_import(self):
        """Test that basic_usage can be imported."""
        sys.path.insert(0, str(EXAMPLES_DIR))
        try:
            import basic_usage  # noqa: F401
        finally:
            sys.path.pop(0)

    def test_sidereal_calculations_import(self):
        """Test that sidereal_calculations can be imported."""
        sys.path.insert(0, str(EXAMPLES_DIR))
        try:
            import sidereal_calculations  # noqa: F401
        finally:
            sys.path.pop(0)

    def test_house_systems_import(self):
        """Test that house_systems can be imported."""
        sys.path.insert(0, str(EXAMPLES_DIR))
        try:
            import house_systems  # noqa: F401
        finally:
            sys.path.pop(0)

    def test_fixed_stars_import(self):
        """Test that fixed_stars can be imported."""
        sys.path.insert(0, str(EXAMPLES_DIR))
        try:
            import fixed_stars  # noqa: F401
        finally:
            sys.path.pop(0)

    def test_rise_set_transit_import(self):
        """Test that rise_set_transit can be imported."""
        sys.path.insert(0, str(EXAMPLES_DIR))
        try:
            import rise_set_transit  # noqa: F401
        finally:
            sys.path.pop(0)

    def test_thread_safe_context_import(self):
        """Test that thread_safe_context can be imported."""
        sys.path.insert(0, str(EXAMPLES_DIR))
        try:
            import thread_safe_context  # noqa: F401
        finally:
            sys.path.pop(0)


class TestBasicUsageExample:
    """Tests for basic_usage.py example functions."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Add examples directory to path."""
        sys.path.insert(0, str(EXAMPLES_DIR))
        yield
        sys.path.pop(0)

    def test_julian_day_example(self):
        """Test Julian Day conversion example runs without error."""
        from basic_usage import example_julian_day

        # Should run without raising exceptions
        example_julian_day()

    def test_planetary_positions_example(self):
        """Test planetary positions example runs without error."""
        from basic_usage import example_planetary_positions

        example_planetary_positions()

    def test_house_cusps_example(self):
        """Test house cusps example runs without error."""
        from basic_usage import example_house_cusps

        example_house_cusps()

    def test_zodiac_formatting_example(self):
        """Test zodiac formatting example runs without error."""
        from basic_usage import example_zodiac_formatting

        example_zodiac_formatting()

    def test_ayanamsha_example(self):
        """Test ayanamsha example runs without error."""
        from basic_usage import example_ayanamsha

        example_ayanamsha()

    def test_planet_name_example(self):
        """Test planet name example runs without error."""
        from basic_usage import example_planet_name

        example_planet_name()


class TestSiderealCalculationsExample:
    """Tests for sidereal_calculations.py example functions."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Add examples directory to path."""
        sys.path.insert(0, str(EXAMPLES_DIR))
        yield
        sys.path.pop(0)

    def test_ayanamsha_overview_example(self):
        """Test ayanamsha overview example runs without error."""
        from sidereal_calculations import example_ayanamsha_overview

        example_ayanamsha_overview()

    def test_tropical_vs_sidereal_example(self):
        """Test tropical vs sidereal example runs without error."""
        from sidereal_calculations import example_tropical_vs_sidereal

        example_tropical_vs_sidereal()

    def test_vedic_chart_example(self):
        """Test Vedic chart example runs without error."""
        from sidereal_calculations import example_vedic_chart

        example_vedic_chart()

    def test_ayanamsha_change_over_time_example(self):
        """Test ayanamsha change over time example runs without error."""
        from sidereal_calculations import example_ayanamsha_change_over_time

        example_ayanamsha_change_over_time()

    def test_nakshatra_calculation_example(self):
        """Test nakshatra calculation example runs without error."""
        from sidereal_calculations import example_nakshatra_calculation

        example_nakshatra_calculation()


class TestHouseSystemsExample:
    """Tests for house_systems.py example functions."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Add examples directory to path."""
        sys.path.insert(0, str(EXAMPLES_DIR))
        yield
        sys.path.pop(0)

    def test_house_systems_overview_example(self):
        """Test house systems overview example runs without error."""
        from house_systems import example_house_systems_overview

        example_house_systems_overview()

    def test_compare_house_systems_example(self):
        """Test compare house systems example runs without error."""
        from house_systems import example_compare_house_systems

        example_compare_house_systems()

    def test_angles_and_vertices_example(self):
        """Test angles and vertices example runs without error."""
        from house_systems import example_angles_and_vertices

        example_angles_and_vertices()

    def test_house_position_example(self):
        """Test house position example runs without error."""
        from house_systems import example_house_position

        example_house_position()

    def test_polar_latitudes_example(self):
        """Test polar latitudes example runs without error."""
        from house_systems import example_polar_latitudes

        example_polar_latitudes()

    def test_extended_houses_example(self):
        """Test extended houses example runs without error."""
        from house_systems import example_extended_houses

        example_extended_houses()


class TestFixedStarsExample:
    """Tests for fixed_stars.py example functions."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Add examples directory to path."""
        sys.path.insert(0, str(EXAMPLES_DIR))
        yield
        sys.path.pop(0)

    def test_star_positions_example(self):
        """Test star positions example runs without error."""
        from fixed_stars import example_star_positions

        example_star_positions()

    def test_star_lookup_example(self):
        """Test star lookup example runs without error."""
        from fixed_stars import example_star_lookup

        example_star_lookup()

    def test_star_magnitude_example(self):
        """Test star magnitude example runs without error."""
        from fixed_stars import example_star_magnitude

        example_star_magnitude()

    def test_proper_motion_example(self):
        """Test proper motion example runs without error."""
        from fixed_stars import example_proper_motion

        example_proper_motion()

    def test_conjunction_with_planet_example(self):
        """Test conjunction with planet example runs without error."""
        from fixed_stars import example_conjunction_with_planet

        example_conjunction_with_planet()

    def test_royal_stars_example(self):
        """Test royal stars example runs without error."""
        from fixed_stars import example_royal_stars

        example_royal_stars()


class TestRiseSetTransitExample:
    """Tests for rise_set_transit.py example functions."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Add examples directory to path."""
        sys.path.insert(0, str(EXAMPLES_DIR))
        yield
        sys.path.pop(0)

    def test_sunrise_sunset_example(self):
        """Test sunrise sunset example runs without error."""
        from rise_set_transit import example_sunrise_sunset

        example_sunrise_sunset()

    def test_twilight_times_example(self):
        """Test twilight times example runs without error."""
        from rise_set_transit import example_twilight_times

        example_twilight_times()

    def test_moonrise_moonset_example(self):
        """Test moonrise moonset example runs without error."""
        from rise_set_transit import example_moonrise_moonset

        example_moonrise_moonset()

    def test_planet_visibility_example(self):
        """Test planet visibility example runs without error."""
        from rise_set_transit import example_planet_visibility

        example_planet_visibility()

    def test_week_of_sunrise_example(self):
        """Test week of sunrise example runs without error."""
        from rise_set_transit import example_week_of_sunrise

        example_week_of_sunrise()

    def test_no_refraction_example(self):
        """Test no refraction example runs without error."""
        from rise_set_transit import example_no_refraction

        example_no_refraction()

    def test_disc_center_vs_limb_example(self):
        """Test disc center vs limb example runs without error."""
        from rise_set_transit import example_disc_center_vs_limb

        example_disc_center_vs_limb()


class TestThreadSafeContextExample:
    """Tests for thread_safe_context.py example functions."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Add examples directory to path."""
        sys.path.insert(0, str(EXAMPLES_DIR))
        yield
        sys.path.pop(0)

    def test_basic_context_example(self):
        """Test basic context example runs without error."""
        from thread_safe_context import example_basic_context

        example_basic_context()

    def test_independent_sidereal_modes_example(self):
        """Test independent sidereal modes example runs without error."""
        from thread_safe_context import example_independent_sidereal_modes

        example_independent_sidereal_modes()

    def test_topocentric_settings_example(self):
        """Test topocentric settings example runs without error."""
        from thread_safe_context import example_topocentric_settings

        example_topocentric_settings()

    def test_multithreaded_calculations_example(self):
        """Test multithreaded calculations example runs without error."""
        from thread_safe_context import example_multithreaded_calculations

        example_multithreaded_calculations()

    def test_context_isolation_example(self):
        """Test context isolation example runs without error."""
        from thread_safe_context import example_context_isolation

        example_context_isolation()

    def test_house_calculations_example(self):
        """Test house calculations example runs without error."""
        from thread_safe_context import example_house_calculations

        example_house_calculations()

    def test_context_reuse_example(self):
        """Test context reuse example runs without error."""
        from thread_safe_context import example_context_reuse

        example_context_reuse()


class TestExampleScriptsExecution:
    """
    Integration tests that run example scripts as subprocesses.
    These tests verify that the scripts can be executed from the command line.
    """

    @pytest.mark.slow
    def test_basic_usage_script(self):
        """Test basic_usage.py runs without error."""
        script_path = EXAMPLES_DIR / "basic_usage.py"
        result = subprocess.run(
            [sys.executable, str(script_path)],
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode == 0, f"Script failed: {result.stderr}"

    @pytest.mark.slow
    def test_sidereal_calculations_script(self):
        """Test sidereal_calculations.py runs without error."""
        script_path = EXAMPLES_DIR / "sidereal_calculations.py"
        result = subprocess.run(
            [sys.executable, str(script_path)],
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode == 0, f"Script failed: {result.stderr}"

    @pytest.mark.slow
    def test_house_systems_script(self):
        """Test house_systems.py runs without error."""
        script_path = EXAMPLES_DIR / "house_systems.py"
        result = subprocess.run(
            [sys.executable, str(script_path)],
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode == 0, f"Script failed: {result.stderr}"

    @pytest.mark.slow
    def test_fixed_stars_script(self):
        """Test fixed_stars.py runs without error."""
        script_path = EXAMPLES_DIR / "fixed_stars.py"
        result = subprocess.run(
            [sys.executable, str(script_path)],
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode == 0, f"Script failed: {result.stderr}"

    @pytest.mark.slow
    def test_rise_set_transit_script(self):
        """Test rise_set_transit.py runs without error."""
        script_path = EXAMPLES_DIR / "rise_set_transit.py"
        result = subprocess.run(
            [sys.executable, str(script_path)],
            capture_output=True,
            text=True,
            timeout=120,  # Longer timeout for rise/set calculations
        )
        assert result.returncode == 0, f"Script failed: {result.stderr}"

    @pytest.mark.slow
    def test_thread_safe_context_script(self):
        """Test thread_safe_context.py runs without error."""
        script_path = EXAMPLES_DIR / "thread_safe_context.py"
        result = subprocess.run(
            [sys.executable, str(script_path)],
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode == 0, f"Script failed: {result.stderr}"
