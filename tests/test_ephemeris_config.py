"""
Tests for ephemeris file configuration.

Tests the set_ephemeris_file() and set_ephe_path() functions to ensure
users can configure which ephemeris file to use.
"""

import os
import pytest
from libephemeris import (
    calc_ut,
    SUN,
    FLG_SPEED,
    julday,
    set_ephemeris_file,
    set_ephe_path,
    set_jpl_file,
    get_library_path,
)
from libephemeris.state import get_planets


def test_set_ephemeris_file_defaults_to_de440():
    """Test that default ephemeris file is de440.bsp"""
    from libephemeris.state import _EPHEMERIS_FILE

    assert _EPHEMERIS_FILE == "de440.bsp"


def test_set_ephemeris_file_changes_file():
    """Test that set_ephemeris_file() changes the ephemeris file"""

    # Set a different file
    set_ephemeris_file("de422.bsp")

    # Verify the global variable was updated
    from libephemeris import state

    assert state._EPHEMERIS_FILE == "de422.bsp"

    # Verify planets cache was cleared
    assert state._PLANETS is None

    # Reset to default
    set_ephemeris_file("de440.bsp")


@pytest.mark.skipif(
    os.environ.get("LIBEPHEMERIS_MODE", "").lower() == "leb",
    reason="get_planets() opens the JPL/SPICE source, disabled in sealed leb mode",
)
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


def test_calculation_works_with_default_de440():
    """Test that calculations work with the default de440.bsp file"""
    # Reset to default
    set_ephemeris_file("de440.bsp")
    set_ephe_path(None)

    # Date within DE440 range (1550-2650)
    year, month, day = 2025, 11, 29
    tjd_ut = julday(year, month, day, 0.0)

    # Should work without error
    pos, retflag = calc_ut(tjd_ut, SUN, FLG_SPEED)

    # Verify we got a valid result
    assert len(pos) == 6
    assert 0 <= pos[0] < 360  # Longitude should be in valid range
    assert pos[2] > 0  # Distance should be positive


def test_calculation_fails_with_de440_out_of_range():
    """Test that calculations fail with DE440 for dates out of range"""
    # Reset to default
    set_ephemeris_file("de440.bsp")
    set_ephe_path(None)

    # Date outside DE440 range
    year, month, day = 3000, 1, 1
    tjd_ut = julday(year, month, day, 0.0)

    # Should raise an error about date range
    with pytest.raises(Exception) as exc_info:
        pos, retflag = calc_ut(tjd_ut, SUN, FLG_SPEED)

    # Verify error message mentions date range (may use different wording)
    error_msg = str(exc_info.value).lower()
    assert "date" in error_msg or "range" in error_msg or "ephemeris" in error_msg


def test_set_jpl_file_changes_file():
    """Test that set_jpl_file() changes the ephemeris file (alias for set_ephemeris_file)"""

    # Set a different file
    set_jpl_file("de422.bsp")

    # Verify the global variable was updated
    from libephemeris import state

    assert state._EPHEMERIS_FILE == "de422.bsp"

    # Verify planets cache was cleared
    assert state._PLANETS is None

    # Reset to default
    set_jpl_file("de440.bsp")


def test_swe_set_jpl_file_alias_works():
    """Test that set_jpl_file() alias works the same as set_jpl_file()"""

    # Set a different file using the swe_ prefixed version
    set_jpl_file("de430.bsp")

    # Verify the global variable was updated
    from libephemeris import state

    assert state._EPHEMERIS_FILE == "de430.bsp"

    # Verify planets cache was cleared
    assert state._PLANETS is None

    # Reset to default
    set_jpl_file("de440.bsp")


def test_set_jpl_file_with_local_path():
    """Test that set_jpl_file() works with set_ephe_path() for local files"""
    import tempfile

    # Create a temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        # Set the ephemeris path and file
        set_ephe_path(tmpdir)
        set_jpl_file("de441.bsp")

        # Verify settings were applied
        from libephemeris import state

        assert state._EPHEMERIS_PATH == tmpdir
        assert state._EPHEMERIS_FILE == "de441.bsp"

        # Reset to defaults
        set_ephe_path(None)
        set_jpl_file("de440.bsp")


def test_get_library_path_returns_default():
    """Test that get_library_path() returns default path when no custom path set"""
    # Reset to default
    set_ephe_path(None)

    path = get_library_path()

    # Should return an absolute path
    assert os.path.isabs(path)

    # Should be ~/.libephemeris (or LIBEPHEMERIS_DATA_DIR if set)
    expected_base = os.path.join(os.path.expanduser("~"), ".libephemeris")
    assert os.path.normpath(path) == os.path.normpath(expected_base)


def test_get_library_path_returns_custom_path():
    """Test that get_library_path() returns custom path when set"""
    import tempfile

    # Reset first
    set_ephe_path(None)

    with tempfile.TemporaryDirectory() as tmpdir:
        # Set custom ephemeris path
        set_ephe_path(tmpdir)

        path = get_library_path()

        # Should return the custom path
        assert os.path.normpath(path) == os.path.normpath(tmpdir)

        # Reset to default
        set_ephe_path(None)


def test_get_library_path_returns_absolute_path():
    """Test that get_library_path() always returns absolute path"""
    # Test with relative custom path
    set_ephe_path("./relative/path")

    path = get_library_path()

    # Should be absolute
    assert os.path.isabs(path)

    # Reset to default
    set_ephe_path(None)


def test_swe_get_library_path_alias():
    """Test that get_library_path() is an alias for get_library_path()"""
    # Reset to default
    set_ephe_path(None)

    path1 = get_library_path()
    path2 = get_library_path()

    # Both should return the same path
    assert path1 == path2


def test_get_library_path_after_close():
    """Test that get_library_path() returns default after close()"""
    import tempfile
    from libephemeris import close

    with tempfile.TemporaryDirectory() as tmpdir:
        # Set custom path
        set_ephe_path(tmpdir)

        # close() should reset the path
        close()

        path = get_library_path()

        # Should return default path (~/.libephemeris), not the custom one
        expected_base = os.path.join(os.path.expanduser("~"), ".libephemeris")
        assert os.path.normpath(path) == os.path.normpath(expected_base)


def test_set_ephe_path_no_argument_resets():
    """set_ephe_path() with no argument resets to default (reference parity).

    The reference API accepts set_ephe_path() with no argument as a reset to
    its default search location; a required-positional signature previously
    raised TypeError on that call.
    """
    from libephemeris import state

    set_ephe_path("/tmp")
    assert state._EPHEMERIS_PATH == "/tmp"
    # No-argument call must not raise TypeError and resets to default (None).
    set_ephe_path()
    assert state._EPHEMERIS_PATH is None
    set_ephe_path(None)


def test_missing_arbitrary_jpl_file_falls_back_without_download():
    """A missing, unrecognized JPL filename falls back to the tier default.

    The reference API, asked for a JPL file it cannot find (e.g. after
    set_jpl_file('bogus.bsp')), falls back to a valid default ephemeris rather
    than fetching the exact missing name. libephemeris must never turn an
    arbitrary configured filename into a network download of that name.
    """
    from libephemeris.state import (
        _is_recognized_de_kernel,
        _resolve_missing_ephemeris_fallback,
        _get_current_tier,
    )
    from libephemeris.logging_config import get_logger

    logger = get_logger()

    # Recognized JPL DE kernels are still eligible for on-demand download.
    assert _is_recognized_de_kernel("de440.bsp")
    assert _is_recognized_de_kernel("de440s.bsp")
    assert _is_recognized_de_kernel("de441.bsp")
    assert not _is_recognized_de_kernel("de_nonexistent_xyz.bsp")
    assert not _is_recognized_de_kernel("random.bsp")

    tier_default = _get_current_tier().ephemeris_file
    # An unrecognized, absent filename resolves to the tier default kernel.
    assert (
        _resolve_missing_ephemeris_fallback("de_nonexistent_xyz.bsp", logger)
        == tier_default
    )
    # A recognized DE kernel name is preserved (may be downloaded on demand),
    # even when absent locally.
    assert _resolve_missing_ephemeris_fallback("de441.bsp", logger) == "de441.bsp"
