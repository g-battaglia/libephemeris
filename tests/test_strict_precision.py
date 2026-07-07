"""
Tests for strict precision mode.

This module tests:
- Strict precision mode configuration (enable/disable)
- Environment variable support
- SPKRequiredError behavior
- Integration with minor body calculations
"""

import os
import pytest
import libephemeris as eph
from libephemeris.constants import (
    CHIRON,
    CERES,
    ERIS,
    PHOLUS,
    FLG_SPEED,
)
from libephemeris import state
from unittest.mock import patch


@pytest.fixture(autouse=True)
def cleanup(monkeypatch):
    """Reset state before and after each test."""
    # Clear any SPK registrations
    state._SPK_BODY_MAP.clear()
    # Disable auto-SPK download so cached SPK files in ~/.libephemeris/spk/
    # don't get re-registered during calc_ut()
    old_auto_spk = state._AUTO_SPK_DOWNLOAD
    state._AUTO_SPK_DOWNLOAD = False
    # Disable LEB so the fast path doesn't intercept calc_ut() before
    # the SPK/strict-precision code path is reached
    monkeypatch.delenv("LIBEPHEMERIS_LEB", raising=False)
    monkeypatch.delenv("LIBEPHEMERIS_MODE", raising=False)
    monkeypatch.setattr(state, "_LEB_FILE", None)
    monkeypatch.setattr(state, "_LEB_READER", None)
    monkeypatch.setattr(state, "_discover_leb_file", lambda: None)
    # Override any --calc-mode leb from root conftest
    eph.set_calc_mode("auto")
    # Reset strict precision setting
    eph.set_strict_precision(None)
    yield
    # Cleanup after test
    state._SPK_BODY_MAP.clear()
    state._AUTO_SPK_DOWNLOAD = old_auto_spk
    eph.set_strict_precision(None)


class TestStrictPrecisionConfiguration:
    """Test set_strict_precision() and get_strict_precision() functions."""

    def test_default_is_enabled(self):
        """Strict precision should be enabled by default."""
        # Ensure no env var is set
        with patch.dict(os.environ, {}, clear=True):
            eph.set_strict_precision(None)
            assert eph.get_strict_precision() is True

    def test_enable_explicitly(self):
        """Should be able to enable strict precision explicitly."""
        eph.set_strict_precision(True)
        assert eph.get_strict_precision() is True

    def test_disable_explicitly(self):
        """Should be able to disable strict precision explicitly."""
        eph.set_strict_precision(False)
        assert eph.get_strict_precision() is False

    def test_reset_to_default(self):
        """Setting None should reset to default/env var."""
        eph.set_strict_precision(True)
        assert eph.get_strict_precision() is True
        eph.set_strict_precision(None)
        # Default is True
        with patch.dict(os.environ, {}, clear=True):
            eph.set_strict_precision(None)
            assert eph.get_strict_precision() is True


class TestStrictPrecisionEnvVar:
    """Test environment variable configuration."""

    def test_env_var_enable(self):
        """LIBEPHEMERIS_STRICT_PRECISION=1 should enable."""
        with patch.dict(os.environ, {"LIBEPHEMERIS_STRICT_PRECISION": "1"}):
            eph.set_strict_precision(None)
            assert eph.get_strict_precision() is True

    def test_env_var_disable(self):
        """LIBEPHEMERIS_STRICT_PRECISION=0 should disable."""
        with patch.dict(os.environ, {"LIBEPHEMERIS_STRICT_PRECISION": "0"}):
            eph.set_strict_precision(None)
            assert eph.get_strict_precision() is False

    def test_env_var_false(self):
        """LIBEPHEMERIS_STRICT_PRECISION=false should disable."""
        with patch.dict(os.environ, {"LIBEPHEMERIS_STRICT_PRECISION": "false"}):
            eph.set_strict_precision(None)
            assert eph.get_strict_precision() is False

    def test_explicit_overrides_env_var(self):
        """Explicit setting should override environment variable."""
        with patch.dict(os.environ, {"LIBEPHEMERIS_STRICT_PRECISION": "0"}):
            eph.set_strict_precision(True)
            assert eph.get_strict_precision() is True


class TestSPKRequiredError:
    """Test SPKRequiredError exception.

    Strict mode refuses the Keplerian fallback only when NO better source
    exists. ASSIST (local N-body data) is such a source, so these tests
    disable it to exercise the genuine "no source at all" raise path
    deterministically -- independent of whether ASSIST data happens to be
    installed on the machine running the suite.
    """

    @pytest.fixture(autouse=True)
    def _no_assist(self):
        with patch(
            "libephemeris.rebound_integration.check_assist_data_available",
            return_value=False,
        ):
            yield

    def test_error_raised_for_chiron_without_spk(self):
        """Chiron without SPK should raise SPKRequiredError in strict mode."""
        eph.set_strict_precision(True)
        with pytest.raises(eph.SPKRequiredError) as exc_info:
            eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)

        assert exc_info.value.body_name == "Chiron"
        assert exc_info.value.body_id == CHIRON
        assert "2060" in str(exc_info.value)  # Horizons ID

    def test_error_raised_for_ceres_without_spk(self):
        """Ceres without SPK should raise SPKRequiredError in strict mode."""
        eph.set_strict_precision(True)
        with pytest.raises(eph.SPKRequiredError) as exc_info:
            eph.calc_ut(2451545.0, CERES, FLG_SPEED)

        assert exc_info.value.body_name == "Ceres"
        assert exc_info.value.body_id == CERES

    def test_no_error_with_strict_disabled(self):
        """Chiron should work with Keplerian fallback when strict disabled."""
        eph.set_strict_precision(False)
        pos, flags = eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)

        # Should return valid position
        assert 0 <= pos[0] < 360  # Longitude in range
        assert -90 <= pos[1] <= 90  # Latitude in range
        assert pos[2] > 0  # Distance positive

    def test_pholus_requires_spk_in_strict_mode(self):
        """Pholus requires SPK in strict mode (it's in SPK_BODY_NAME_MAP)."""
        eph.set_strict_precision(True)
        with pytest.raises(eph.SPKRequiredError):
            eph.calc_ut(2451545.0, PHOLUS, FLG_SPEED)

    def test_eris_requires_spk_in_strict_mode(self):
        """Eris requires SPK in strict mode (it's in SPK_BODY_NAME_MAP)."""
        eph.set_strict_precision(True)
        with pytest.raises(eph.SPKRequiredError):
            eph.calc_ut(2451545.0, ERIS, FLG_SPEED)


class TestExportedFunctions:
    """Test that strict precision functions are properly exported."""

    def test_set_strict_precision_exported(self):
        """set_strict_precision should be exported from libephemeris."""
        assert hasattr(eph, "set_strict_precision")
        assert callable(eph.set_strict_precision)

    def test_get_strict_precision_exported(self):
        """get_strict_precision should be exported from libephemeris."""
        assert hasattr(eph, "get_strict_precision")
        assert callable(eph.get_strict_precision)

    def test_spk_required_error_exported(self):
        """SPKRequiredError should be exported from libephemeris."""
        assert hasattr(eph, "SPKRequiredError")
        assert issubclass(eph.SPKRequiredError, Exception)

    def test_functions_in_all(self):
        """Strict precision functions should be in __all__."""
        from libephemeris import __all__

        assert "set_strict_precision" in __all__
        assert "get_strict_precision" in __all__
        assert "SPKRequiredError" in __all__
