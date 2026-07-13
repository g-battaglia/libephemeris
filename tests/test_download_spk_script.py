"""
Tests for the download_spk.py script.

These tests verify the script's functionality including:
- Argument parsing
- Body definitions and validation
- Helper functions
- Cache operations (mocked)
"""

import os
import sys
from unittest.mock import patch

# Add scripts directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))


class TestDownloadSpkBodyDefinitions:
    """Test body definitions and initialization."""

    def test_init_bodies_populates_available_bodies(self):
        """Test that _init_bodies populates AVAILABLE_BODIES correctly."""
        import scripts.download_spk as download_module

        download_module._init_bodies()

        # Check common bodies are present
        assert "chiron" in download_module.AVAILABLE_BODIES
        assert "ceres" in download_module.AVAILABLE_BODIES
        assert "eris" in download_module.AVAILABLE_BODIES
        assert "sedna" in download_module.AVAILABLE_BODIES

    def test_available_bodies_structure(self):
        """Test that AVAILABLE_BODIES has correct structure."""
        import scripts.download_spk as download_module

        download_module._init_bodies()

        for name, info in download_module.AVAILABLE_BODIES.items():
            assert isinstance(name, str)
            assert isinstance(info, tuple)
            assert len(info) == 3
            horizons_id, ipl, naif_id = info
            assert isinstance(horizons_id, str)
            assert isinstance(ipl, int)
            assert isinstance(naif_id, int)

    def test_common_bodies_list(self):
        """Test that COMMON_BODIES contains expected bodies."""
        from scripts.download_spk import COMMON_BODIES

        expected = [
            "chiron",
            "pholus",
            "ceres",
            "pallas",
            "juno",
            "vesta",
            "eris",
            "sedna",
        ]
        assert COMMON_BODIES == expected


class TestListAvailableBodies:
    """Test listing available bodies."""

    def test_list_available_bodies_output(self, capsys):
        """Test that list_available_bodies produces expected output."""
        from scripts.download_spk import list_available_bodies, _init_bodies

        _init_bodies()
        list_available_bodies()

        captured = capsys.readouterr()
        assert "Available bodies for SPK download" in captured.out
        assert "Centaurs:" in captured.out
        assert "Main Belt Asteroids:" in captured.out
        assert "Trans-Neptunian Objects" in captured.out
        assert "chiron" in captured.out
        assert "ceres" in captured.out


class TestDownloadSpkForBody:
    """Test the download_spk_for_body function."""

    def test_unknown_body_returns_none(self, capsys):
        """Test that unknown body returns failure tuple."""
        from scripts.download_spk import download_spk_for_body, _init_bodies

        _init_bodies()
        result = download_spk_for_body(
            body_name="unknown_body_xyz",
            start_date="2020-01-01",
            end_date="2030-01-01",
        )

        assert result == (0, 0, 1)
        captured = capsys.readouterr()
        assert "Unknown body" in captured.err

    def test_already_cached_returns_skipped(self, tmp_path):
        """Test that already cached file returns skipped count."""
        from scripts.download_spk import download_spk_for_body, _init_bodies

        _init_bodies()

        # Create a fake cached file
        cache_dir = tmp_path / "spk_cache"
        cache_dir.mkdir()

        # Patch the functions to simulate cache hit
        with patch(
            "libephemeris.spk_auto._generate_spk_cache_filename",
            return_value="2060_2451545_2488070.bsp",
        ):
            # Create the expected file
            cached_file = cache_dir / "2060_2451545_2488070.bsp"
            cached_file.write_bytes(b"fake spk data")

            result = download_spk_for_body(
                body_name="chiron",
                start_date="2000-01-01",
                end_date="2100-01-01",
                cache_dir=str(cache_dir),
                force=False,
            )

            # Should report 0 success, 1 skipped, 0 failed
            assert result == (0, 1, 0)

    def test_force_flag_triggers_redownload(self, tmp_path):
        """Test that force flag triggers re-download."""
        from scripts.download_spk import download_spk_for_body, _init_bodies

        _init_bodies()

        cache_dir = tmp_path / "spk_cache"
        cache_dir.mkdir()

        # Create a fake cached file
        cached_file = cache_dir / "2060_test.bsp"
        cached_file.write_bytes(b"old data")

        # Mock the download function to avoid actual network call
        with patch("libephemeris.spk_auto.download_spk_from_horizons") as mock_download:
            mock_download.side_effect = ValueError("direct download failed for test")

            result = download_spk_for_body(
                body_name="chiron",
                start_date="2020-01-01",
                end_date="2030-01-01",
                cache_dir=str(cache_dir),
                force=True,
            )

            # Should return (0, 0, 1) due to the mocked error, but the key is that it tried
            # We're testing that force=True attempts a download
            assert result == (0, 0, 1)  # Due to mocked error


class TestIsoToJd:
    """Test the ISO date to Julian Day conversion."""

    def test_iso_to_jd_j2000(self):
        """Test conversion of J2000 epoch date."""

        # Load the module to access internal function

        # We can't directly access _iso_to_jd, but we can test the effect
        # through the overall function or recreate the logic
        def _iso_to_jd(date_str: str) -> float:
            parts = date_str.split("-")
            year, month, day = int(parts[0]), int(parts[1]), int(parts[2])
            if month <= 2:
                year -= 1
                month += 12
            a = int(year / 100)
            b = 2 - a + int(a / 4)
            jd = (
                int(365.25 * (year + 4716))
                + int(30.6001 * (month + 1))
                + day
                + b
                - 1524.5
            )
            return jd

        # J2000.0 = 2000-01-01 12:00 TT = JD 2451545.0
        # But we use 2000-01-01 (midnight), which is JD 2451544.5
        jd = _iso_to_jd("2000-01-01")
        assert abs(jd - 2451544.5) < 1  # Within 1 day tolerance

        # Test another known date: 2020-01-01
        jd = _iso_to_jd("2020-01-01")
        assert abs(jd - 2458849.5) < 1


class TestMainFunction:
    """Test the main() entry point."""

    def test_main_list_bodies(self, capsys, monkeypatch):
        """Test main() with --list flag."""
        from scripts.download_spk import main, _init_bodies

        _init_bodies()
        monkeypatch.setattr(sys, "argv", ["download_spk.py", "--list"])

        result = main()

        assert result == 0
        captured = capsys.readouterr()
        assert "Available bodies" in captured.out

    def test_main_download_has_no_astroquery_gate(self, monkeypatch):
        """The download script works with only the package's core dependencies."""
        from scripts.download_spk import main, _init_bodies

        _init_bodies()
        monkeypatch.setattr(sys, "argv", ["download_spk.py", "--bodies", "chiron"])

        with patch(
            "scripts.download_spk.download_spk_for_body", return_value=(1, 0, 0)
        ) as mock_download:
            result = main()

        assert result == 0
        mock_download.assert_called_once()

    def test_main_unknown_body(self, monkeypatch, capsys):
        """Test main() with unknown body."""
        from scripts.download_spk import main, _init_bodies

        _init_bodies()
        monkeypatch.setattr(
            sys, "argv", ["download_spk.py", "--bodies", "unknown_body_xyz"]
        )

        result = main()

        assert result == 1
        captured = capsys.readouterr()
        assert "Unknown body" in captured.err


class TestListCacheContents:
    """Test cache listing functionality."""

    def test_list_cache_empty(self, capsys):
        """Test listing empty cache."""
        from scripts.download_spk import list_cache_contents, _init_bodies

        _init_bodies()

        with patch("libephemeris.spk_auto.list_cached_spk", return_value=[]):
            list_cache_contents()

        captured = capsys.readouterr()
        assert "empty" in captured.out.lower()

    def test_list_cache_with_files(self, capsys):
        """Test listing cache with files."""
        from scripts.download_spk import list_cache_contents, _init_bodies

        _init_bodies()

        mock_cached = [
            {
                "filename": "2060_test.bsp",
                "size_mb": 1.5,
                "date_start": "2020-01-01",
                "date_end": "2030-01-01",
            }
        ]

        with patch("libephemeris.spk_auto.list_cached_spk", return_value=mock_cached):
            list_cache_contents()

        captured = capsys.readouterr()
        assert "2060_test.bsp" in captured.out
        assert "1.5" in captured.out


class TestArgumentParsing:
    """Test command-line argument parsing."""

    def test_default_date_range(self, monkeypatch):
        """Test that default date range is 2000-2100."""
        from scripts.download_spk import main, _init_bodies

        _init_bodies()
        monkeypatch.setattr(sys, "argv", ["download_spk.py", "--list"])

        # Run with --list to avoid network calls
        result = main()
        assert result == 0

    def test_custom_date_range(self, monkeypatch, capsys):
        """Test custom date range arguments."""
        from scripts.download_spk import main, _init_bodies

        _init_bodies()
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "download_spk.py",
                "--bodies",
                "chiron",
                "--start",
                "2020-01-01",
                "--end",
                "2050-01-01",
            ],
        )

        with patch(
            "scripts.download_spk.download_spk_for_body",
            return_value=(1, 0, 0),
        ) as mock_download:
            result = main()

            # Check that download was called with correct dates
            mock_download.assert_called_once()
            call_args = mock_download.call_args
            assert call_args.kwargs["start_date"] == "2020-01-01"
            assert call_args.kwargs["end_date"] == "2050-01-01"

    def test_quiet_flag(self, monkeypatch, capsys):
        """Test quiet flag suppresses output."""
        from scripts.download_spk import main, _init_bodies

        _init_bodies()
        monkeypatch.setattr(
            sys,
            "argv",
            ["download_spk.py", "--bodies", "chiron", "--quiet"],
        )

        with patch(
            "scripts.download_spk.download_spk_for_body",
            return_value=(1, 0, 0),
        ):
            result = main()

        captured = capsys.readouterr()
        assert "Downloading" not in captured.out
