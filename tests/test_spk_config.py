"""
Tests for SPK configuration options in libephemeris.

These tests verify:
- set_spk_cache_dir / get_spk_cache_dir configuration
- set_spk_date_padding / get_spk_date_padding configuration
- Integration of configuration with auto_get_spk and enable_auto_spk
"""

import os
import pytest

import libephemeris as eph
from libephemeris import state, spk_auto


class TestSpkCacheDir:
    """Test SPK cache directory configuration."""

    def setup_method(self):
        """Clear state before each test."""
        state.close()

    def teardown_method(self):
        """Clear state after each test."""
        state.close()

    def test_get_spk_cache_dir_default(self):
        """Default SPK cache dir is None."""
        assert eph.get_spk_cache_dir() is None

    def test_set_spk_cache_dir(self):
        """set_spk_cache_dir sets the cache directory."""
        eph.set_spk_cache_dir("/custom/cache/path")
        assert eph.get_spk_cache_dir() == "/custom/cache/path"

    def test_set_spk_cache_dir_none(self):
        """set_spk_cache_dir(None) clears the setting."""
        eph.set_spk_cache_dir("/custom/path")
        assert eph.get_spk_cache_dir() == "/custom/path"
        eph.set_spk_cache_dir(None)
        assert eph.get_spk_cache_dir() is None

    def test_close_resets_spk_cache_dir(self):
        """close() resets SPK cache dir to None."""
        eph.set_spk_cache_dir("/custom/path")
        eph.close()
        assert eph.get_spk_cache_dir() is None

    def test_spk_cache_dir_exported(self):
        """set_spk_cache_dir and get_spk_cache_dir are exported."""
        assert hasattr(eph, "set_spk_cache_dir")
        assert hasattr(eph, "get_spk_cache_dir")


class TestSpkDatePadding:
    """Test SPK date padding configuration."""

    def setup_method(self):
        """Clear state before each test."""
        state.close()

    def teardown_method(self):
        """Clear state after each test."""
        state.close()

    def test_get_spk_date_padding_default(self):
        """Default SPK date padding is 0."""
        assert eph.get_spk_date_padding() == 0

    def test_set_spk_date_padding(self):
        """set_spk_date_padding sets the padding value."""
        eph.set_spk_date_padding(365)
        assert eph.get_spk_date_padding() == 365

    def test_set_spk_date_padding_zero(self):
        """set_spk_date_padding(0) sets padding to zero."""
        eph.set_spk_date_padding(30)
        eph.set_spk_date_padding(0)
        assert eph.get_spk_date_padding() == 0

    def test_set_spk_date_padding_negative_raises(self):
        """set_spk_date_padding raises ValueError for negative values."""
        with pytest.raises(ValueError, match="non-negative"):
            eph.set_spk_date_padding(-1)

    def test_close_resets_spk_date_padding(self):
        """close() resets SPK date padding to 0."""
        eph.set_spk_date_padding(100)
        eph.close()
        assert eph.get_spk_date_padding() == 0

    def test_spk_date_padding_exported(self):
        """set_spk_date_padding and get_spk_date_padding are exported."""
        assert hasattr(eph, "set_spk_date_padding")
        assert hasattr(eph, "get_spk_date_padding")


class TestAutoSpkDownload:
    """Test auto SPK download configuration (existing functionality)."""

    def setup_method(self):
        """Clear state before each test."""
        state.close()

    def teardown_method(self):
        """Clear state after each test."""
        state.close()

    def test_get_auto_spk_download_default(self):
        """Default auto SPK download is True."""
        eph.set_auto_spk_download(None)
        assert eph.get_auto_spk_download() is True

    def test_set_auto_spk_download_true(self):
        """set_auto_spk_download(True) enables auto download."""
        eph.set_auto_spk_download(True)
        assert eph.get_auto_spk_download() is True

    def test_set_auto_spk_download_false(self):
        """set_auto_spk_download(False) disables auto download."""
        eph.set_auto_spk_download(True)
        eph.set_auto_spk_download(False)
        assert eph.get_auto_spk_download() is False

    def test_set_auto_spk_download_none(self):
        """set_auto_spk_download(None) uses environment variable."""
        eph.set_auto_spk_download(True)
        eph.set_auto_spk_download(None)
        # Should return True since default is now enabled
        assert eph.get_auto_spk_download() is True

    def test_auto_spk_download_exported(self):
        """set_auto_spk_download and get_auto_spk_download are exported."""
        assert hasattr(eph, "set_auto_spk_download")
        assert hasattr(eph, "get_auto_spk_download")


class TestAutoSpkConfigIntegration:
    """Test that AutoSpkConfig uses global configuration."""

    def setup_method(self):
        """Clear state before each test."""
        state.close()

    def teardown_method(self):
        """Clear state after each test."""
        state.close()

    def test_autospkconfig_uses_global_cache_dir(self, tmp_path):
        """AutoSpkConfig.get_cache_dir() uses global cache_dir setting."""
        custom_cache = str(tmp_path / "custom_cache")
        eph.set_spk_cache_dir(custom_cache)

        config = spk_auto.AutoSpkConfig(
            ipl=15,
            body_id="2060",
            start="2000-01-01",
            end="2100-01-01",
            cache_dir=None,  # Not specified - should use global
        )

        cache_dir = config.get_cache_dir()
        assert cache_dir == custom_cache
        assert os.path.exists(custom_cache)

    def test_autospkconfig_explicit_overrides_global(self, tmp_path):
        """Explicit cache_dir in AutoSpkConfig overrides global setting."""
        global_cache = str(tmp_path / "global_cache")
        explicit_cache = str(tmp_path / "explicit_cache")
        eph.set_spk_cache_dir(global_cache)

        config = spk_auto.AutoSpkConfig(
            ipl=15,
            body_id="2060",
            start="2000-01-01",
            end="2100-01-01",
            cache_dir=explicit_cache,  # Explicit - should override global
        )

        cache_dir = config.get_cache_dir()
        assert cache_dir == explicit_cache
        assert os.path.exists(explicit_cache)


class TestConfigurationDocumentation:
    """Test that configuration functions have proper docstrings."""

    def test_set_spk_cache_dir_has_docstring(self):
        """set_spk_cache_dir has a docstring."""
        assert state.set_spk_cache_dir.__doc__ is not None
        assert "cache" in state.set_spk_cache_dir.__doc__.lower()

    def test_get_spk_cache_dir_has_docstring(self):
        """get_spk_cache_dir has a docstring."""
        assert state.get_spk_cache_dir.__doc__ is not None

    def test_set_spk_date_padding_has_docstring(self):
        """set_spk_date_padding has a docstring."""
        assert state.set_spk_date_padding.__doc__ is not None
        assert "padding" in state.set_spk_date_padding.__doc__.lower()

    def test_get_spk_date_padding_has_docstring(self):
        """get_spk_date_padding has a docstring."""
        assert state.get_spk_date_padding.__doc__ is not None


class TestAllExports:
    """Test that all new functions are in __all__."""

    def test_set_spk_cache_dir_in_all(self):
        """set_spk_cache_dir is in __all__."""
        assert "set_spk_cache_dir" in eph.__all__

    def test_get_spk_cache_dir_in_all(self):
        """get_spk_cache_dir is in __all__."""
        assert "get_spk_cache_dir" in eph.__all__

    def test_set_spk_date_padding_in_all(self):
        """set_spk_date_padding is in __all__."""
        assert "set_spk_date_padding" in eph.__all__

    def test_get_spk_date_padding_in_all(self):
        """get_spk_date_padding is in __all__."""
        assert "get_spk_date_padding" in eph.__all__
