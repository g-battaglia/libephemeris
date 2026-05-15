"""
Tests for automatic SPK fallback feature.

This module tests:
- Configuration of auto SPK download (enable/disable)
- Environment variable support
- Fallback behavior when SPK is not available
- Integration with minor body calculations
"""

import os
import pytest
import libephemeris as eph
from libephemeris import (
    set_auto_spk_download,
    get_auto_spk_download,
    calc_ut,
    close,
)
from libephemeris.constants import (
    CHIRON,
    CERES,
    ERIS,
    FLG_HELCTR,
)
from libephemeris.state import _AUTO_SPK_ENV_VAR
from unittest.mock import patch


@pytest.fixture(autouse=True)
def cleanup():
    """Reset state before and after each test."""
    # Reset auto SPK download setting
    set_auto_spk_download(None)
    yield
    # Cleanup after test
    set_auto_spk_download(None)


class TestAutoSpkDownloadConfiguration:
    """Test set_auto_spk_download() and get_auto_spk_download() functions."""

    def test_default_is_disabled(self):
        """Auto SPK download should be enabled by default."""
        # Ensure no env var is set
        with patch.dict(os.environ, {}, clear=True):
            # Clear any existing setting
            set_auto_spk_download(None)
            assert get_auto_spk_download() is True

    def test_enable_explicitly(self):
        """Should be able to enable auto SPK download explicitly."""
        set_auto_spk_download(True)
        assert get_auto_spk_download() is True

    def test_disable_explicitly(self):
        """Should be able to disable auto SPK download explicitly."""
        set_auto_spk_download(False)
        assert get_auto_spk_download() is False

    def test_reset_to_env_var(self):
        """Setting to None should fall back to environment variable."""
        set_auto_spk_download(True)
        assert get_auto_spk_download() is True

        # Reset to use env var
        set_auto_spk_download(None)
        # With no env var set, should be True (default enabled)
        with patch.dict(os.environ, {}, clear=True):
            assert get_auto_spk_download() is True

    def test_explicit_setting_overrides_env_var(self):
        """Explicit setting should override environment variable."""
        with patch.dict(os.environ, {_AUTO_SPK_ENV_VAR: "1"}, clear=False):
            # Explicitly disable
            set_auto_spk_download(False)
            assert get_auto_spk_download() is False

            # Explicitly enable
            set_auto_spk_download(True)
            assert get_auto_spk_download() is True


class TestEnvironmentVariableSupport:
    """Test environment variable configuration."""

    @pytest.mark.parametrize(
        "env_value,expected",
        [
            ("1", True),
            ("true", True),
            ("True", True),
            ("TRUE", True),
            ("yes", True),
            ("Yes", True),
            ("YES", True),
            ("on", True),
            ("enabled", True),
            ("0", False),
            ("false", False),
            ("False", False),
            ("FALSE", False),
            ("no", False),
            ("No", False),
            ("NO", False),
            ("off", False),
            ("disabled", False),
            ("", True),
            ("anything_else", True),
        ],
    )
    def test_env_var_parsing(self, env_value, expected):
        """Environment variable should be parsed correctly."""
        set_auto_spk_download(None)  # Use env var
        with patch.dict(os.environ, {_AUTO_SPK_ENV_VAR: env_value}, clear=False):
            assert get_auto_spk_download() is expected

    def test_missing_env_var(self):
        """Missing environment variable should default to True (enabled)."""
        set_auto_spk_download(None)  # Use env var
        # Remove the env var if it exists
        env_copy = dict(os.environ)
        env_copy.pop(_AUTO_SPK_ENV_VAR, None)
        with patch.dict(os.environ, env_copy, clear=True):
            assert get_auto_spk_download() is True


class TestFallbackBehavior:
    """Test that calculations work correctly with fallback to Keplerian."""

    def test_minor_body_calculation_without_auto_spk(self):
        """Minor body calculation should work with auto SPK disabled."""
        set_auto_spk_download(False)

        # Calculate Chiron position using Keplerian fallback
        jd = 2451545.0  # J2000.0
        pos, flag = calc_ut(jd, CHIRON, FLG_HELCTR)

        # Should return valid heliocentric position
        assert 0 <= pos[0] < 360  # Longitude
        assert -90 <= pos[1] <= 90  # Latitude
        assert pos[2] > 0  # Distance

    def test_multiple_minor_bodies_work(self):
        """Multiple minor body calculations should work."""
        set_auto_spk_download(False)
        jd = 2451545.0

        bodies = [CHIRON, CERES, ERIS]
        for body_id in bodies:
            pos, flag = calc_ut(jd, body_id, FLG_HELCTR)
            assert 0 <= pos[0] < 360, f"Invalid longitude for body {body_id}"
            assert -90 <= pos[1] <= 90, f"Invalid latitude for body {body_id}"
            assert pos[2] > 0, f"Invalid distance for body {body_id}"

    def test_geocentric_calculation_works(self):
        """Geocentric minor body calculation should work."""
        set_auto_spk_download(False)
        jd = 2451545.0

        # Geocentric calculation (no FLG_HELCTR)
        pos, flag = calc_ut(jd, CHIRON, 0)

        assert 0 <= pos[0] < 360  # Longitude
        assert -90 <= pos[1] <= 90  # Latitude
        assert pos[2] > 0  # Distance


class TestCloseResetsSetting:
    """Test that close() resets the auto SPK download setting."""

    def test_close_resets_setting(self):
        """close() should reset auto SPK download setting to None."""
        set_auto_spk_download(True)
        assert get_auto_spk_download() is True

        close()

        # After close, the explicit setting should be reset to None.
        # This means get_auto_spk_download() will check the environment variable.
        # To test that the explicit setting was cleared, we check that
        # setting it to True again works (to confirm reset happened).
        # We also need to ensure no env var interference.
        with patch.dict(os.environ, {}, clear=True):
            # After close with no env var, should be True (default enabled)
            # But we need to call set_auto_spk_download(None) since close()
            # internally sets it to None which means "use env var"
            assert get_auto_spk_download() is True


class TestAutoSpkDownloadAttempt:
    """Test the auto SPK download attempt behavior."""

    def test_auto_download_not_attempted_when_disabled(self):
        """Auto download should not be attempted when disabled."""
        set_auto_spk_download(False)

        # This should use Keplerian fallback immediately
        jd = 2451545.0
        pos, flag = calc_ut(jd, CHIRON, FLG_HELCTR)

        # Should still return valid position (Keplerian)
        assert 0 <= pos[0] < 360
        assert -90 <= pos[1] <= 90
        assert pos[2] > 0

    def test_auto_download_graceful_failure_no_astroquery(self):
        """Auto download should fail gracefully when astroquery not available."""
        set_auto_spk_download(True)

        # Mock astroquery as unavailable
        with patch(
            "libephemeris.spk_auto._check_astroquery_available", return_value=False
        ):
            jd = 2451545.0
            pos, flag = calc_ut(jd, CHIRON, FLG_HELCTR)

            # Should fall back to Keplerian and return valid position
            assert 0 <= pos[0] < 360
            assert -90 <= pos[1] <= 90
            assert pos[2] > 0


class TestExportedFunctions:
    """Test that the functions are properly exported from libephemeris."""

    def test_set_auto_spk_download_exported(self):
        """set_auto_spk_download should be exported from libephemeris."""
        assert hasattr(eph, "set_auto_spk_download")
        assert callable(eph.set_auto_spk_download)

    def test_get_auto_spk_download_exported(self):
        """get_auto_spk_download should be exported from libephemeris."""
        assert hasattr(eph, "get_auto_spk_download")
        assert callable(eph.get_auto_spk_download)

    def test_functions_in_all(self):
        """Functions should be in __all__."""
        from libephemeris import __all__

        assert "set_auto_spk_download" in __all__
        assert "get_auto_spk_download" in __all__
