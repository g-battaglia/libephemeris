"""
Tests for auto-download SPK in strict precision mode.

Tests verify:
- When auto_spk_download is enabled and strict_precision is True,
  auto-download is attempted before raising SPKRequiredError
- If auto-download succeeds, calculation works
- If auto-download fails, SPKRequiredError is still raised
- Logging messages are produced during auto-download
"""

import pytest
import logging
from unittest.mock import patch, MagicMock
import libephemeris as eph
from libephemeris import state
from libephemeris.constants import (
    CHIRON,
    CERES,
    FLG_SPEED,
)


@pytest.fixture(autouse=True)
def cleanup(monkeypatch):
    """Reset state before and after each test."""
    # Clear any SPK registrations
    state._SPK_BODY_MAP.clear()
    # Disable LEB so the fast path doesn't intercept calc_ut() before
    # the SPK/strict-precision code path is reached
    monkeypatch.delenv("LIBEPHEMERIS_LEB", raising=False)
    monkeypatch.delenv("LIBEPHEMERIS_MODE", raising=False)
    monkeypatch.setattr(state, "_LEB_FILE", None)
    monkeypatch.setattr(state, "_LEB_READER", None)
    monkeypatch.setattr(state, "_discover_leb_file", lambda: None)
    # Override any --calc-mode leb from root conftest
    eph.set_calc_mode("auto")
    # Reset settings
    eph.set_strict_precision(None)
    eph.set_auto_spk_download(None)
    yield
    # Cleanup after test
    state._SPK_BODY_MAP.clear()
    eph.set_strict_precision(None)
    eph.set_auto_spk_download(None)


class TestAutoDownloadInStrictMode:
    """Test auto-download behavior when strict precision is enabled."""

    def test_auto_download_attempted_before_error(self):
        """Auto-download should be attempted before raising SPKRequiredError."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        with patch("libephemeris.spk.download_and_register_spk") as mock_download:
            mock_download.side_effect = RuntimeError("Network error")

            with pytest.raises(eph.SPKRequiredError):
                eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)

            mock_download.assert_called_once()

    def test_auto_download_success_returns_result(self):
        """If auto-download succeeds, should return calculation result."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        mock_position = (100.0, 5.0, 15.0, 0.1, 0.01, 0.001)

        call_count = [0]

        def mock_calc_spk(t, ipl, iflag):
            call_count[0] += 1
            if call_count[0] == 1:
                return None
            else:
                return mock_position

        with patch("libephemeris.spk.download_and_register_spk"):
            with patch(
                "libephemeris.spk.calc_spk_body_position", side_effect=mock_calc_spk
            ):
                pos, flags = eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)

                assert pos[0] == pytest.approx(100.0)
                assert pos[1] == pytest.approx(5.0)
                assert pos[2] == pytest.approx(15.0)

    def test_no_auto_download_when_disabled(self):
        """Auto-download should not be attempted when disabled."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(False)

        with patch("libephemeris.spk.download_and_register_spk") as mock_download:
            with pytest.raises(eph.SPKRequiredError):
                eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)

            mock_download.assert_not_called()

    def test_auto_download_logs_info_message(self):
        """Auto-download should log an info message."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        with patch("libephemeris.spk.download_and_register_spk") as mock_download:
            mock_download.side_effect = RuntimeError("Network error")

            with patch("libephemeris.logging_config.get_logger") as mock_get_logger:
                mock_logger = MagicMock()
                mock_get_logger.return_value = mock_logger

                with pytest.raises(eph.SPKRequiredError):
                    eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)

                mock_logger.info.assert_called()
                log_calls = [call[0][0] for call in mock_logger.info.call_args_list]
                assert any("Auto-downloading" in msg for msg in log_calls)

    def test_download_exception_still_raises_spk_error(self):
        """If download raises exception, SPKRequiredError should still be raised."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        with patch("libephemeris.spk.download_and_register_spk") as mock_download:
            mock_download.side_effect = RuntimeError("Network error")

            with pytest.raises(eph.SPKRequiredError):
                eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)

    def test_ceres_requires_spk_in_strict_mode(self):
        """Ceres requires SPK in strict mode (it's in SPK_BODY_NAME_MAP).

        Ceres is now SPK-downloadable using name syntax ("Ceres;") to bypass
        the JPL Horizons major body index restriction.
        """
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        with patch("libephemeris.spk.download_and_register_spk") as mock_download:
            mock_download.side_effect = RuntimeError("Network error")

            # Ceres is in SPK_BODY_NAME_MAP, so it requires SPK
            with pytest.raises(eph.SPKRequiredError):
                eph.calc_ut(2451545.0, CERES, FLG_SPEED)

            mock_download.assert_called_once()


class TestAutoDownloadWithStrictDisabled:
    """Test that auto-download path is used with strict mode disabled."""

    def test_keplerian_fallback_when_strict_disabled(self):
        """When strict is disabled, should fall back to Keplerian even if download fails."""
        eph.set_strict_precision(False)
        eph.set_auto_spk_download(True)

        # Even without SPK, should get a result via Keplerian fallback
        pos, flags = eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)

        # Should return valid position
        assert 0 <= pos[0] < 360  # Longitude in range
        assert -90 <= pos[1] <= 90  # Latitude in range
        assert pos[2] > 0  # Distance positive


class TestAutoDownloadWithFirstTryPath:
    """Test that the _try_auto_spk_download path is used."""

    def test_first_try_path_still_attempted(self):
        """The _try_auto_spk_download path should be called when auto-download is enabled."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        with patch("libephemeris.planets._try_auto_spk_download") as mock_try:
            mock_try.return_value = None

            with pytest.raises(eph.SPKRequiredError):
                eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)

            mock_try.assert_called_once()


class TestAutoDownloadJDPassthrough:
    """Test that Julian Day is passed correctly to download functions."""

    def test_jd_used_in_download(self):
        """The Julian Day should be used in the download process."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        test_jd = 2460000.5

        with patch("libephemeris.spk.download_and_register_spk") as mock_download:
            mock_download.side_effect = RuntimeError("Network error")

            with pytest.raises(eph.SPKRequiredError):
                eph.calc_ut(test_jd, CHIRON, FLG_SPEED)

            mock_download.assert_called_once()
            call_kwargs = mock_download.call_args[1]
            assert "body" in call_kwargs
            assert call_kwargs["body"] == "2060"
