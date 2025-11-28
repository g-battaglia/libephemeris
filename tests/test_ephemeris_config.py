"""
Tests for ephemeris file configuration.

Tests the set_ephemeris_file() and set_ephe_path() functions to ensure
users can configure which ephemeris file to use.
"""

import pytest
from libephemeris import (
    swe_calc_ut,
    SE_SUN,
    SEFLG_SPEED,
    swe_julday,
    set_ephemeris_file,
    set_ephe_path,
)
from libephemeris.state import get_planets


def test_set_ephemeris_file_defaults_to_de421():
    """Test that default ephemeris file is de421.bsp"""
    from libephemeris.state import _EPHEMERIS_FILE
    
    assert _EPHEMERIS_FILE == "de421.bsp"


def test_set_ephemeris_file_changes_file():
    """Test that set_ephemeris_file() changes the ephemeris file"""
    from libephemeris.state import _EPHEMERIS_FILE, _PLANETS
    
    # Set a different file
    set_ephemeris_file("de422.bsp")
    
    # Verify the global variable was updated
    from libephemeris import state
    assert state._EPHEMERIS_FILE == "de422.bsp"
    
    # Verify planets cache was cleared
    assert state._PLANETS is None
    
    # Reset to default
    set_ephemeris_file("de421.bsp")


def test_set_ephe_path_clears_cache():
    """Test that set_ephe_path() clears the planets cache"""
    from libephemeris import state
    
    # Load planets first
    _ = get_planets()
    assert state._PLANETS is not None
    
    # Set ephemeris path
    set_ephe_path("/tmp")
    
    # Verify cache was cleared
    assert state._PLANETS is None
    
    # Reset to None
    set_ephe_path(None)


def test_calculation_works_with_default_de421():
    """Test that calculations work with the default de421.bsp file"""
    # Reset to default
    set_ephemeris_file("de421.bsp")
    set_ephe_path(None)
    
    # Date within DE421 range (1900-2050)
    year, month, day = 2025, 11, 29
    tjd_ut = swe_julday(year, month, day, 0.0)
    
    # Should work without error
    pos, retflag = swe_calc_ut(tjd_ut, SE_SUN, SEFLG_SPEED)
    
    # Verify we got a valid result
    assert len(pos) == 6
    assert 0 <= pos[0] < 360  # Longitude should be in valid range
    assert pos[2] > 0  # Distance should be positive


def test_calculation_fails_with_de421_out_of_range():
    """Test that calculations fail with DE421 for dates out of range"""
    # Reset to default
    set_ephemeris_file("de421.bsp")
    set_ephe_path(None)
    
    # Date outside DE421 range
    year, month, day = 3000, 1, 1
    tjd_ut = swe_julday(year, month, day, 0.0)
    
    # Should raise an error about date range
    with pytest.raises(Exception) as exc_info:
        pos, retflag = swe_calc_ut(tjd_ut, SE_SUN, SEFLG_SPEED)
    
    # Verify error message mentions date range
    assert "ephemeris segment only covers dates" in str(exc_info.value)
