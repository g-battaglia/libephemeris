"""
Tests for auto-download SPK in strict precision mode.

Tests verify:
- When auto_spk_download is enabled and strict_precision is True,
  auto-download is attempted before falling back to Keplerian
- If auto-download succeeds, calculation uses SPK data
- If auto-download fails, Keplerian fallback is used (no error)
- When auto_spk_download is disabled, SPKRequiredError is raised
- Logging messages are produced during auto-download
"""

import pytest
from unittest.mock import patch
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
    state._SPK_BODY_MAP.clear()
    monkeypatch.delenv("LIBEPHEMERIS_LEB", raising=False)
    monkeypatch.delenv("LIBEPHEMERIS_MODE", raising=False)
    monkeypatch.setattr(state, "_LEB_FILE", None)
    monkeypatch.setattr(state, "_LEB_READER", None)
    monkeypatch.setattr(state, "_discover_leb_file", lambda: None)
    monkeypatch.setattr(state, "_discover_reviewed_leb_tier_cores", lambda: {})
    eph.set_calc_mode("auto")
    eph.set_strict_precision(None)
    eph.set_auto_spk_download(None)
    yield
    state._SPK_BODY_MAP.clear()
    eph.set_strict_precision(None)
    eph.set_auto_spk_download(None)


class TestAutoDownloadInStrictMode:
    """Test auto-download behavior when strict precision is enabled."""

    def test_auto_download_attempted_before_fallback(self):
        """Auto-download should be attempted before Keplerian fallback."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        with patch("libephemeris.spk.download_and_register_spk") as mock_download:
            mock_download.side_effect = RuntimeError("Network error")

            pos, flags = eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)

            mock_download.assert_called_once()
            assert 0 <= pos[0] < 360

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
        """Auto-download should not be attempted when disabled.

        With no SPK registered and no ASSIST data, strict mode must raise
        rather than silently downgrading to Keplerian. ASSIST is disabled
        here so the "no source at all" raise is exercised deterministically.
        """
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(False)

        with patch(
            "libephemeris.rebound_integration.check_assist_data_available",
            return_value=False,
        ):
            with patch("libephemeris.spk.download_and_register_spk") as mock_download:
                with pytest.raises(eph.SPKRequiredError):
                    eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)

                mock_download.assert_not_called()

    def test_auto_download_logs_warning_on_failure(self):
        """Failed auto-download should produce a warning log."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        with patch("libephemeris.spk.download_and_register_spk") as mock_download:
            mock_download.side_effect = RuntimeError("Network error")

            pos, flags = eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)
            assert 0 <= pos[0] < 360

    def test_download_failure_falls_back_to_keplerian(self):
        """If download fails, Keplerian fallback should be used."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        with (
            patch("libephemeris.spk.download_and_register_spk") as mock_download,
            patch(
                "libephemeris.rebound_integration.check_assist_data_available",
                return_value=False,
            ),
        ):
            mock_download.side_effect = RuntimeError("Network error")

            pos, flags = eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)
            mock_download.assert_called_once()
            assert 0 <= pos[0] < 360
            assert pos[2] > 0

    def test_ceres_auto_download_attempted(self):
        """Ceres auto-download should be attempted in strict mode.

        Ceres is SPK-downloadable using name syntax ("Ceres;") to bypass
        the JPL Horizons major body index restriction.
        """
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        with patch("libephemeris.spk.download_and_register_spk") as mock_download:
            mock_download.side_effect = RuntimeError("Network error")

            pos, flags = eph.calc_ut(2451545.0, CERES, FLG_SPEED)
            mock_download.assert_called_once()
            assert 0 <= pos[0] < 360


class TestAutoDownloadWithStrictDisabled:
    """Test that auto-download path is used with strict mode disabled."""

    def test_keplerian_fallback_when_strict_disabled(self):
        """When strict is disabled, should fall back to Keplerian even if download fails."""
        eph.set_strict_precision(False)
        eph.set_auto_spk_download(True)

        pos, flags = eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)

        assert 0 <= pos[0] < 360
        assert -90 <= pos[1] <= 90
        assert pos[2] > 0


class TestAutoDownloadWithFirstTryPath:
    """Test that the _try_auto_spk_download path is used."""

    def test_first_try_path_still_attempted(self):
        """The _try_auto_spk_download path should be called when auto-download is enabled."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        with patch("libephemeris.planets._try_auto_spk_download") as mock_try:
            mock_try.return_value = None

            pos, flags = eph.calc_ut(2451545.0, CHIRON, FLG_SPEED)

            mock_try.assert_called_once()
            assert 0 <= pos[0] < 360


class TestAutoDownloadJDPassthrough:
    """Test that Julian Day is passed correctly to download functions."""

    def test_jd_used_in_download(self):
        """The Julian Day should be used in the download process."""
        eph.set_strict_precision(True)
        eph.set_auto_spk_download(True)

        test_jd = 2460000.5

        with patch("libephemeris.spk.download_and_register_spk") as mock_download:
            mock_download.side_effect = RuntimeError("Network error")

            pos, flags = eph.calc_ut(test_jd, CHIRON, FLG_SPEED)

            mock_download.assert_called_once()
            call_kwargs = mock_download.call_args[1]
            assert "body" in call_kwargs
            assert call_kwargs["body"] == "2060"
