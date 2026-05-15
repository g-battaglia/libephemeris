"""
Tests for get_current_file_data() function.

This tests the function that returns information about the currently loaded
ephemeris file: file path, date range, and ephemeris type.
"""

import os
import pytest
from libephemeris import (
    get_current_file_data,
    calc_ut,
    SUN,
    FLG_SPEED,
    set_ephemeris_file,
    set_ephe_path,
    close,
)


@pytest.fixture(autouse=True)
def reset_state(monkeypatch):
    """Reset ephemeris state before and after each test.

    Disables LEB so calc_ut goes through Skyfield and populates _PLANETS.
    """
    # Disable LEB env vars
    monkeypatch.delenv("LIBEPHEMERIS_LEB", raising=False)
    monkeypatch.delenv("LIBEPHEMERIS_MODE", raising=False)
    # Disable LEB auto-discovery and reader
    from libephemeris import state

    monkeypatch.setattr(state, "_LEB_FILE", None)
    monkeypatch.setattr(state, "_LEB_READER", None)
    monkeypatch.setattr(state, "_discover_leb_file", lambda: None)
    close()
    set_ephemeris_file("de440.bsp")
    set_ephe_path(None)
    yield
    close()
    set_ephemeris_file("de440.bsp")
    set_ephe_path(None)


def test_get_current_file_data_before_loading():
    """Test that get_current_file_data returns empty data before ephemeris is loaded."""
    path, start, end, denum = get_current_file_data(0)

    assert path == ""
    assert start == 0.0
    assert end == 0.0
    assert denum == 0


def test_get_current_file_data_after_loading():
    """Test that get_current_file_data returns valid data after ephemeris is loaded."""
    # Trigger ephemeris loading
    calc_ut(2451545.0, SUN, FLG_SPEED)

    path, start, end, denum = get_current_file_data(0)

    # Path should be non-empty and contain de440.bsp
    assert path != ""
    assert "de440.bsp" in path

    # DE440 covers 1549-12-31 to 2650-01-25
    # JD 2287184.5 = 1549-12-31
    # JD 2688976.5 = 2650-01-25
    assert start == pytest.approx(2287184.5, rel=1e-6)
    assert end == pytest.approx(2688976.5, rel=1e-6)

    # DE number should be 440
    assert denum == 440


def test_get_current_file_data_ifno_0_and_1_same():
    """Test that ifno=0 (planets) and ifno=1 (moon) return same data for JPL files."""
    # Trigger ephemeris loading
    calc_ut(2451545.0, SUN, FLG_SPEED)

    path0, start0, end0, denum0 = get_current_file_data(0)
    path1, start1, end1, denum1 = get_current_file_data(1)

    # Both should return the same data since JPL files contain both
    assert path0 == path1
    assert start0 == start1
    assert end0 == end1
    assert denum0 == denum1


def test_get_current_file_data_ifno_2_not_applicable():
    """Test that ifno=2 (asteroids) returns empty data."""
    # Trigger ephemeris loading
    calc_ut(2451545.0, SUN, FLG_SPEED)

    path, start, end, denum = get_current_file_data(2)

    assert path == ""
    assert start == 0.0
    assert end == 0.0
    assert denum == 0


def test_get_current_file_data_ifno_3_not_applicable():
    """Test that ifno=3 (other asteroids) returns empty data."""
    # Trigger ephemeris loading
    calc_ut(2451545.0, SUN, FLG_SPEED)

    path, start, end, denum = get_current_file_data(3)

    assert path == ""
    assert start == 0.0
    assert end == 0.0
    assert denum == 0


def test_get_current_file_data_ifno_4_not_applicable():
    """Test that ifno=4 (stars) returns empty data."""
    # Trigger ephemeris loading
    calc_ut(2451545.0, SUN, FLG_SPEED)

    path, start, end, denum = get_current_file_data(4)

    assert path == ""
    assert start == 0.0
    assert end == 0.0
    assert denum == 0


def test_get_current_file_data_returns_tuple():
    """Test that get_current_file_data returns a tuple of the correct types."""
    calc_ut(2451545.0, SUN, FLG_SPEED)

    result = get_current_file_data(0)

    assert isinstance(result, tuple)
    assert len(result) == 4

    path, start, end, denum = result
    assert isinstance(path, str)
    assert isinstance(start, float)
    assert isinstance(end, float)
    assert isinstance(denum, int)


def test_swe_get_current_file_data_alias():
    """Test that get_current_file_data is an alias for get_current_file_data."""
    calc_ut(2451545.0, SUN, FLG_SPEED)

    result1 = get_current_file_data(0)
    result2 = get_current_file_data(0)

    assert result1 == result2


def test_get_current_file_data_after_close():
    """Test that get_current_file_data returns empty data after close()."""
    # Trigger ephemeris loading
    calc_ut(2451545.0, SUN, FLG_SPEED)

    # Verify data is available
    path, start, end, denum = get_current_file_data(0)
    assert path != ""

    # Close and reset
    close()

    # Now should return empty data
    path, start, end, denum = get_current_file_data(0)
    assert path == ""
    assert start == 0.0
    assert end == 0.0
    assert denum == 0


def test_get_current_file_data_date_range_valid():
    """Test that the date range returned is sensible."""
    calc_ut(2451545.0, SUN, FLG_SPEED)

    path, start, end, denum = get_current_file_data(0)

    # End should be after start
    assert end > start

    # Date range should span at least 100 years (DE files typically span more)
    days_span = end - start
    years_span = days_span / 365.25
    assert years_span > 100
