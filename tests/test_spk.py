"""
Tests for SPK kernel support in libephemeris.

These tests verify:
- SPK body registration and unregistration
- NAIF ID deduction
- Fallback to Keplerian when SPK not available
- Coverage checking
- Error handling (SPKNotFoundError)
- jplephem version compatibility

Note: Integration tests that download actual SPK files from Horizons
are skipped by default. Set LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1 to run them.
"""

import json
import os
import pytest

pytestmark = pytest.mark.usefixtures("allow_mocked_network")
from unittest.mock import patch, MagicMock
from packaging import version

import jplephem

import libephemeris as eph
from libephemeris import spk, state
from libephemeris.constants import (
    CHIRON,
    CERES,
    ERIS,
    NAIF_CHIRON,
    NAIF_CERES,
    NAIF_ERIS,
    NAIF_ASTEROID_OFFSET,
    FLG_SPEED,
)
from libephemeris.exceptions import SPKNotFoundError


class TestJplephemVersion:
    """Test jplephem version requirements for optimal SPK reading performance."""

    # Minimum required version for all features (from skyfield dependency)
    MINIMUM_VERSION = "2.13"
    # Recommended version with NumPy compatibility fixes
    RECOMMENDED_VERSION = "2.24"

    def test_jplephem_installed(self):
        """Verify jplephem is installed."""
        assert jplephem is not None
        assert hasattr(jplephem, "__version__")

    def test_jplephem_minimum_version(self):
        """Verify jplephem meets minimum version requirement.

        Version 2.13+ is required by skyfield for:
        - OutOfRangeError with array attribute for date errors
        - Essential SPK reading capabilities
        """
        installed = version.parse(jplephem.__version__)
        minimum = version.parse(self.MINIMUM_VERSION)
        assert installed >= minimum, (
            f"jplephem {jplephem.__version__} is below minimum version "
            f"{self.MINIMUM_VERSION} required for SPK support"
        )

    def test_jplephem_recommended_version(self):
        """Verify jplephem meets recommended version for optimal performance.

        Version 2.24+ includes:
        - NumPy deprecation warning fixes (.shape -> .reshape())
        - Better compatibility with modern NumPy versions
        """
        installed = version.parse(jplephem.__version__)
        recommended = version.parse(self.RECOMMENDED_VERSION)
        assert installed >= recommended, (
            f"jplephem {jplephem.__version__} is below recommended version "
            f"{self.RECOMMENDED_VERSION}. Consider upgrading for NumPy compatibility."
        )

    def test_jplephem_spk_module_available(self):
        """Verify jplephem.spk module is available for SPK reading."""
        from jplephem.spk import SPK

        assert SPK is not None
        # Verify SPK has essential methods
        assert hasattr(SPK, "open")

    def test_jplephem_daf_module_available(self):
        """Verify jplephem.daf module is available for DAF file reading."""
        from jplephem.daf import DAF

        assert DAF is not None


class TestNaifIdDeduction:
    """Test NAIF ID deduction from body identifiers."""

    def test_extract_number_from_digits(self):
        """Extract asteroid number from pure digit string."""
        assert spk._extract_asteroid_number("2060") == 2060
        assert spk._extract_asteroid_number("136199") == 136199

    def test_extract_number_with_name(self):
        """Extract asteroid number from 'number name' format."""
        assert spk._extract_asteroid_number("2060 Chiron") == 2060
        assert spk._extract_asteroid_number("136199 Eris") == 136199

    def test_extract_number_with_parentheses(self):
        """Extract asteroid number from '(number)' format."""
        assert spk._extract_asteroid_number("(2060)") == 2060
        assert spk._extract_asteroid_number("(136199) Eris") == 136199

    def test_extract_number_name_only(self):
        """Return None for name-only identifiers."""
        assert spk._extract_asteroid_number("Chiron") is None
        assert spk._extract_asteroid_number("Eris") is None

    def test_extract_number_rejects_provisional_designations(self):
        """Provisional designations (year + survey code) are NOT catalog numbers.

        Regression: the leading discovery year was mistaken for an asteroid
        number ("2003 UB313" -> 2003), fabricating a NAIF ID absent from the
        downloaded kernel.
        """
        assert spk._extract_asteroid_number("2003 UB313") is None
        assert spk._extract_asteroid_number("1998 KY26") is None
        assert spk._extract_asteroid_number("2014 MU69") is None

    def test_extract_number_rejects_cometary_designations(self):
        """Cometary designations (orbit-type letter + slash) are not numbers."""
        assert spk._extract_asteroid_number("1P/Halley") is None
        assert spk._extract_asteroid_number("C/1995 O1") is None
        assert spk._extract_asteroid_number("2P/Encke") is None

    def test_extract_number_year_named_asteroid_still_works(self):
        """A real numbered asteroid whose number reads like a year (2003
        Harding) is still extracted -- the trailing token is a proper name,
        not a two-letter survey code."""
        assert spk._extract_asteroid_number("2003 Harding") == 2003

    def test_deduce_naif_id(self):
        """Deduce NAIF ID from body identifier."""
        assert spk._deduce_naif_id("2060") == 2060 + NAIF_ASTEROID_OFFSET
        assert spk._deduce_naif_id("2060 Chiron") == 2060 + NAIF_ASTEROID_OFFSET
        assert spk._deduce_naif_id("Chiron") is None


class TestBuildHorizonsUrl:
    """Test Horizons API URL building."""

    def test_basic_url(self):
        """Build basic Horizons URL."""
        url = spk._build_horizons_url("Chiron", "2000-01-01", "2100-01-01")
        assert "COMMAND='Chiron'" in url
        assert "EPHEM_TYPE=SPK" in url
        assert "CENTER='500@0'" in url
        assert "START_TIME='2000-01-01'" in url
        assert "STOP_TIME='2100-01-01'" in url

    def test_custom_center(self):
        """Build URL with custom center."""
        url = spk._build_horizons_url(
            "Chiron", "2000-01-01", "2100-01-01", center="@sun"
        )
        assert "CENTER='@sun'" in url


class TestSpkCoverageReuse:
    """Test cached SPK coverage checks."""

    def test_auto_spk_date_padding_crosses_year_boundary(self):
        from libephemeris.planets import _pad_iso_date

        assert _pad_iso_date("1900-01-01", -45) == "1899-11-17"
        assert _pad_iso_date("2100-01-01", 45) == "2100-02-15"

    def test_spk_covers_dates_accepts_covering_kernel(self, monkeypatch):
        start_jd = spk._date_to_jd("1900-01-01")
        end_jd = spk._date_to_jd("2100-01-01")
        monkeypatch.setattr(
            spk, "get_spk_coverage", lambda _path: (start_jd - 1.0, end_jd + 1.0)
        )

        assert spk._spk_covers_dates("cached.bsp", "1900-01-01", "2100-01-01")

    def test_spk_covers_dates_rejects_narrow_kernel(self, monkeypatch):
        start_jd = spk._date_to_jd("1900-01-01")
        end_jd = spk._date_to_jd("2100-01-01")
        monkeypatch.setattr(
            spk, "get_spk_coverage", lambda _path: (start_jd + 2.0, end_jd + 1.0)
        )

        assert not spk._spk_covers_dates("cached.bsp", "1900-01-01", "2100-01-01")

    def test_download_and_register_redownloads_narrow_cache(self, tmp_path):
        cached = str(tmp_path / "cached.bsp")
        overwrites = []

        def fake_download(**kwargs):
            overwrites.append(kwargs["overwrite"])
            return cached

        with (
            patch("libephemeris.spk.download_spk", side_effect=fake_download),
            patch("libephemeris.spk._spk_covers_dates", return_value=False),
            patch("libephemeris.spk.register_spk_body") as register,
        ):
            out = spk.download_and_register_spk(
                "2060",
                CHIRON,
                "1900-01-01",
                "2100-01-01",
                naif_id=NAIF_CHIRON,
            )

        assert out == cached
        assert overwrites == [False, True]
        register.assert_called_once_with(CHIRON, cached, NAIF_CHIRON)


class TestSanitizeFilename:
    """Test filename sanitization."""

    def test_simple_name(self):
        """Sanitize simple name."""
        assert spk._sanitize_filename("Chiron") == "chiron"

    def test_name_with_spaces(self):
        """Sanitize name with spaces."""
        assert spk._sanitize_filename("2060 Chiron") == "2060_chiron"

    def test_name_with_parentheses(self):
        """Sanitize name with special characters."""
        assert spk._sanitize_filename("(136199) Eris") == "136199_eris"


class TestSpkRegistration:
    """Test SPK body registration."""

    def setup_method(self):
        """Clear state before each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def teardown_method(self):
        """Clear state after each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def test_list_spk_bodies_empty(self):
        """List SPK bodies when none registered."""
        assert eph.list_spk_bodies() == {}

    def test_get_spk_body_info_not_registered(self):
        """Get info for unregistered body."""
        assert eph.get_spk_body_info(CHIRON) is None

    def test_unregister_nonexistent(self):
        """Unregister non-existent body (should not error)."""
        eph.unregister_spk_body(CHIRON)  # Should not raise

    def test_register_file_not_found(self):
        """Register with non-existent file raises SPKNotFoundError."""
        with pytest.raises(SPKNotFoundError) as exc_info:
            eph.register_spk_body(CHIRON, "/nonexistent/file.bsp", NAIF_CHIRON)

        # Verify the error contains helpful information
        error = exc_info.value
        assert "/nonexistent/file.bsp" in str(error)
        assert error.filepath == "/nonexistent/file.bsp"

    def test_register_file_not_found_with_body_name(self):
        """SPKNotFoundError includes body name when available."""
        with pytest.raises(SPKNotFoundError) as exc_info:
            eph.register_spk_body(CHIRON, "/path/to/missing.bsp", NAIF_CHIRON)

        error = exc_info.value
        # Should have Chiron as body name
        assert error.body_name == "Chiron"
        # Should include Horizons ID
        assert error.body_id == "2060"
        # Error message should include instructions
        assert "download_spk" in str(error)

    def test_register_file_not_found_with_eris(self):
        """SPKNotFoundError works for ERIS (TNO with high ID offset)."""
        with pytest.raises(SPKNotFoundError) as exc_info:
            eph.register_spk_body(ERIS, "/path/to/missing_eris.bsp", NAIF_ERIS)

        error = exc_info.value
        assert error.body_name == "Eris"
        assert error.body_id == "136199"

    def test_spk_not_found_error_helpful_message(self):
        """SPKNotFoundError message includes multiple options to obtain SPK."""
        with pytest.raises(SPKNotFoundError) as exc_info:
            eph.register_spk_body(CERES, "/missing/ceres.bsp", NAIF_CERES)

        message = str(exc_info.value)
        # Check that key instructions are included
        assert "download_spk" in message
        assert "download_and_register_spk" in message
        assert "set_auto_spk_download" in message
        assert "scripts/download_spk.py" in message


class TestCalcWithoutSpk:
    """Test calculations without SPK (Keplerian fallback)."""

    def setup_method(self):
        """Clear SPK state before each test."""
        state._SPK_BODY_MAP.clear()

    def teardown_method(self):
        """Clear SPK state after each test."""
        state._SPK_BODY_MAP.clear()

    def test_chiron_keplerian(self):
        """Calculate Chiron position using Keplerian model."""
        # J2000.0 epoch
        jd = 2451545.0
        pos, flags = eph.calc_ut(jd, CHIRON, FLG_SPEED)

        # Should return valid position
        assert 0 <= pos[0] < 360  # Longitude in range
        assert -90 <= pos[1] <= 90  # Latitude in range
        assert pos[2] > 0  # Distance positive

    def test_ceres_keplerian(self):
        """Calculate Ceres position using Keplerian model."""
        jd = 2451545.0
        pos, flags = eph.calc_ut(jd, CERES, FLG_SPEED)

        assert 0 <= pos[0] < 360
        assert -90 <= pos[1] <= 90
        assert pos[2] > 0


class TestSpkCalcIntegration:
    """Integration tests for SPK calculations.

    These tests mock the SPK kernel to avoid network dependency.
    """

    def setup_method(self):
        """Clear state before each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def teardown_method(self):
        """Clear state after each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def test_spk_calc_returns_none_when_not_registered(self):
        """calc_spk_body_position returns None for unregistered bodies."""
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        result = spk.calc_spk_body_position(t, CHIRON, 0)
        assert result is None


class TestNaifConstants:
    """Test NAIF ID constants."""

    def test_naif_offset(self):
        """NAIF offset is 2000000."""
        assert NAIF_ASTEROID_OFFSET == 2000000

    def test_naif_chiron(self):
        """Chiron NAIF ID is correct."""
        assert NAIF_CHIRON == 2002060

    def test_naif_ceres(self):
        """Ceres NAIF ID is correct."""
        assert NAIF_CERES == 2000001

    def test_naif_eris(self):
        """Eris NAIF ID is correct."""
        assert NAIF_ERIS == 2136199


class TestEphemerisContextSpk:
    """Test SPK support in EphemerisContext."""

    def setup_method(self):
        """Clear state before each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def teardown_method(self):
        """Clear state after each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def test_context_has_spk_methods(self):
        """EphemerisContext has SPK registration methods."""
        ctx = eph.EphemerisContext()
        assert hasattr(ctx, "register_spk_body")
        assert hasattr(ctx, "unregister_spk_body")
        assert hasattr(ctx, "get_spk_body_info")
        assert hasattr(ctx, "list_spk_bodies")

    def test_context_list_spk_bodies_empty(self):
        """Context lists empty SPK bodies when none registered."""
        ctx = eph.EphemerisContext()
        assert ctx.list_spk_bodies() == {}

    def test_context_get_spk_body_info_not_registered(self):
        """Context returns None for unregistered body."""
        ctx = eph.EphemerisContext()
        assert ctx.get_spk_body_info(CHIRON) is None


# =============================================================================
# NETWORK-DEPENDENT TESTS (Skipped by default)
# =============================================================================


@pytest.mark.skipif(
    os.environ.get("LIBEPHEMERIS_TEST_SPK_DOWNLOAD") != "1",
    reason="Set LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1 to run network tests",
)
class TestSpkDownloadIntegration:
    """Integration tests that download actual SPK files.

    These tests require network access and may be slow.
    Skipped by default; set LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1 to run.
    """

    def setup_method(self):
        """Clear state before each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def teardown_method(self):
        """Clear state after each test."""
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def test_download_chiron_spk(self, tmp_path):
        """Download Chiron SPK from Horizons."""
        path = eph.download_spk(
            body="2060",
            start="2020-01-01",
            end="2025-01-01",
            path=str(tmp_path),
        )

        assert os.path.exists(path)
        assert path.endswith(".bsp")

    def test_download_and_register(self, tmp_path):
        """Download and register Chiron SPK."""
        path = eph.download_and_register_spk(
            body="2060",
            ipl=CHIRON,
            start="2020-01-01",
            end="2025-01-01",
            path=str(tmp_path),
        )

        assert os.path.exists(path)
        assert eph.get_spk_body_info(CHIRON) is not None

        # Calculate position using SPK
        jd = 2459215.5  # 2021-01-01
        pos, _ = eph.calc_ut(jd, CHIRON, FLG_SPEED)

        assert 0 <= pos[0] < 360
        assert pos[2] > 0


# =============================================================================
# LOGGING TESTS
# =============================================================================


class TestDownloadSpkProgressBar:
    """Test progress bar functionality for SPK downloads.

    Note: download_spk() now uses direct Horizons API with base64-encoded
    response, so there's no progress bar during download. These tests
    are skipped but kept for documentation purposes.
    """

    @pytest.mark.skip(reason="download_spk uses direct API, no progress bar")
    def test_progress_bar_shows_for_large_files_on_tty(self, tmp_path):
        pass

    @pytest.mark.skip(reason="download_spk uses direct API, no progress bar")
    def test_no_progress_bar_for_small_files(self, tmp_path):
        pass

    @pytest.mark.skip(reason="download_spk uses direct API, no progress bar")
    def test_no_progress_bar_when_not_tty(self, tmp_path):
        pass

    @pytest.mark.skip(reason="download_spk uses direct API, no progress bar")
    def test_no_progress_bar_when_content_length_unknown(self, tmp_path):
        pass


@pytest.fixture
def capture_libephemeris_logs():
    """Fixture to capture logs from libephemeris logger.

    The libephemeris logger is configured when get_logger() is first called
    by download_spk(). We need to configure the logger ourselves BEFORE
    that happens to ensure DEBUG level and our capture handler are set up.
    """
    import logging
    from libephemeris import logging_config

    captured_records = []

    class RecordCapture(logging.Handler):
        def emit(self, record):
            captured_records.append(record)

    logger = logging.getLogger("libephemeris")

    original_level = logger.level
    original_handlers = list(logger.handlers)
    original_configured = logging_config._logger_configured

    logger.handlers.clear()
    logger.setLevel(logging.DEBUG)

    handler = RecordCapture()
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)

    logging_config._logger_configured = True

    yield captured_records

    logging_config._logger_configured = original_configured
    logger.removeHandler(handler)
    logger.handlers.clear()
    for h in original_handlers:
        logger.addHandler(h)
    logger.setLevel(original_level)


class TestDownloadSpkLogging:
    """Test logging output from download_spk() function."""

    def test_download_logs_start_message_with_naif_id(
        self, tmp_path, capture_libephemeris_logs
    ):
        """download_spk() logs start message with NAIF ID for numbered asteroids."""
        import base64

        mock_spk_data = b"\x00" * 1024
        mock_json_response = {
            "spk_file_id": "test_id",
            "spk": base64.b64encode(mock_spk_data).decode("ascii"),
        }

        with (
            patch("urllib.request.urlopen") as mock_urlopen,
            patch("libephemeris.spk._is_valid_bsp", return_value=True),
        ):
            mock_response = MagicMock()
            mock_response.read.return_value = json.dumps(mock_json_response).encode()
            mock_response.__enter__ = MagicMock(return_value=mock_response)
            mock_response.__exit__ = MagicMock(return_value=False)
            mock_urlopen.return_value = mock_response

            path = eph.download_spk(
                body="2060",
                start="2020-01-01",
                end="2025-01-01",
                path=str(tmp_path),
            )

        records = capture_libephemeris_logs
        messages = [r.getMessage() for r in records]
        assert any("NAIF 2002060" in msg for msg in messages)
        assert any("Downloading SPK" in msg for msg in messages)

    def test_download_logs_start_message_without_naif_id(
        self, tmp_path, capture_libephemeris_logs
    ):
        """download_spk() logs start message without NAIF ID for named bodies."""
        import base64

        mock_spk_data = b"\x00" * 1024
        mock_json_response = {
            "spk_file_id": "test_id",
            "spk": base64.b64encode(mock_spk_data).decode("ascii"),
        }

        with (
            patch("urllib.request.urlopen") as mock_urlopen,
            patch("libephemeris.spk._is_valid_bsp", return_value=True),
        ):
            mock_response = MagicMock()
            mock_response.read.return_value = json.dumps(mock_json_response).encode()
            mock_response.__enter__ = MagicMock(return_value=mock_response)
            mock_response.__exit__ = MagicMock(return_value=False)
            mock_urlopen.return_value = mock_response

            path = eph.download_spk(
                body="Chiron",
                start="2020-01-01",
                end="2025-01-01",
                path=str(tmp_path),
            )

        records = capture_libephemeris_logs
        messages = [r.getMessage() for r in records]
        start_messages = [m for m in messages if "Downloading SPK" in m]
        assert len(start_messages) >= 1
        assert any("Chiron" in msg for msg in start_messages)

    def test_download_logs_completion_with_file_size(
        self, tmp_path, capture_libephemeris_logs
    ):
        """download_spk() logs completion message with file size."""
        import base64

        mock_spk_data = b"\x00" * (2 * 1024 * 1024)
        mock_json_response = {
            "spk_file_id": "test_id",
            "spk": base64.b64encode(mock_spk_data).decode("ascii"),
        }

        with (
            patch("urllib.request.urlopen") as mock_urlopen,
            patch("libephemeris.spk._is_valid_bsp", return_value=True),
        ):
            mock_response = MagicMock()
            mock_response.read.return_value = json.dumps(mock_json_response).encode()
            mock_response.__enter__ = MagicMock(return_value=mock_response)
            mock_response.__exit__ = MagicMock(return_value=False)
            mock_urlopen.return_value = mock_response

            path = eph.download_spk(
                body="2060",
                start="2020-01-01",
                end="2025-01-01",
                path=str(tmp_path),
            )

        records = capture_libephemeris_logs
        messages = [r.getMessage() for r in records]
        saved_messages = [m for m in messages if "SPK saved" in m]
        assert len(saved_messages) == 1
        assert "MB" in saved_messages[0]

    def test_download_logs_warning_on_network_error(
        self, tmp_path, capture_libephemeris_logs
    ):
        """download_spk() logs warning on network error."""
        import urllib.error

        with patch("urllib.request.urlopen") as mock_urlopen:
            mock_urlopen.side_effect = urllib.error.URLError("Connection refused")

            with pytest.raises(ConnectionError):
                eph.download_spk(
                    body="2060",
                    start="2020-01-01",
                    end="2025-01-01",
                    path=str(tmp_path),
                )

        # Verify warning was logged
        records = capture_libephemeris_logs
        warning_records = [r for r in records if r.levelname == "WARNING"]
        assert len(warning_records) >= 1
        messages = [r.getMessage() for r in warning_records]
        assert any("SPK download" in m and "failed" in m for m in messages)

    def test_download_logs_warning_on_http_error(
        self, tmp_path, capture_libephemeris_logs
    ):
        """download_spk() logs warning on HTTP error."""
        import urllib.error

        with patch("urllib.request.urlopen") as mock_urlopen:
            mock_error = urllib.error.HTTPError(
                url="https://test.com",
                code=500,
                msg="Internal Server Error",
                hdrs=None,  # type: ignore[arg-type]
                fp=None,
            )
            mock_urlopen.side_effect = mock_error

            with pytest.raises(ConnectionError):
                eph.download_spk(
                    body="2060",
                    start="2020-01-01",
                    end="2025-01-01",
                    path=str(tmp_path),
                )

        # Verify warning was logged
        records = capture_libephemeris_logs
        warning_records = [r for r in records if r.levelname == "WARNING"]
        assert len(warning_records) >= 1
        messages = [r.getMessage() for r in warning_records]
        assert any("HTTP 500" in m for m in messages)

    def test_download_skips_existing_file_no_log(
        self, tmp_path, capture_libephemeris_logs
    ):
        """download_spk() doesn't log when file already exists."""
        # Create existing file
        existing_file = tmp_path / "2060_202001_202501.bsp"
        existing_file.write_bytes(b"\x00" * 100)

        with patch("libephemeris.spk._is_valid_bsp", return_value=True):
            path = eph.download_spk(
                body="2060",
                start="2020-01-01",
                end="2025-01-01",
                path=str(tmp_path),
                overwrite=False,
            )

        # Should return immediately without logging
        assert path == str(existing_file)
        records = capture_libephemeris_logs
        download_messages = [
            r.getMessage() for r in records if "Downloading SPK" in r.getMessage()
        ]
        assert len(download_messages) == 0
