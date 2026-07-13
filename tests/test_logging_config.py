"""
Tests for libephemeris centralized logging configuration.

This module tests the logging_config module to ensure:
- Logger is properly configured
- Log messages have the correct format
- Logging can be enabled/disabled
- Log levels can be changed
"""

from __future__ import annotations

import io
import logging
import os
import sys
from unittest.mock import patch


class TestGetLogger:
    """Tests for get_logger() function."""

    def test_get_logger_returns_logger(self):
        """get_logger should return a logging.Logger instance."""
        from libephemeris.logging_config import get_logger

        # Reset the logger configured flag for testing
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()
        assert isinstance(logger, logging.Logger)
        assert logger.name == "libephemeris"

    def test_get_logger_has_handler(self):
        """Logger should have at least one handler configured."""
        from libephemeris.logging_config import get_logger
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()
        assert len(logger.handlers) >= 1

    def test_get_logger_handler_is_stderr(self):
        """Logger handler should output to stderr."""
        from libephemeris.logging_config import get_logger
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()
        # At least one handler should be a StreamHandler to stderr
        stderr_handlers = [
            h
            for h in logger.handlers
            if isinstance(h, logging.StreamHandler) and h.stream == sys.stderr
        ]
        assert len(stderr_handlers) >= 1

    def test_get_logger_is_idempotent(self):
        """Calling get_logger multiple times should return the same logger."""
        from libephemeris.logging_config import get_logger
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger1 = get_logger()
        handler_count = len(logger1.handlers)

        logger2 = get_logger()
        # Should not add duplicate handlers
        assert len(logger2.handlers) == handler_count
        assert logger1 is logger2


class TestLogFormat:
    """Tests for log message formatting."""

    def test_log_format_contains_libephemeris_prefix(self):
        """Log messages should be prefixed with [libephemeris]."""
        from libephemeris.logging_config import get_logger
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()

        # Capture stderr output
        captured = io.StringIO()
        with patch.object(sys, "stderr", captured):
            # Create a new handler that writes to our captured stream
            test_handler = logging.StreamHandler(captured)
            test_handler.setLevel(logging.INFO)
            from libephemeris.logging_config import LibephemerisFormatter

            test_handler.setFormatter(LibephemerisFormatter())

            # Temporarily add our test handler and set logger level to INFO
            original_level = logger.level
            logger.setLevel(logging.INFO)
            logger.addHandler(test_handler)
            try:
                logger.info("Test message")
            finally:
                logger.removeHandler(test_handler)
                logger.setLevel(original_level)

        output = captured.getvalue()
        assert "[libephemeris]" in output

    def test_log_format_contains_level(self):
        """Log messages should contain the log level."""
        from libephemeris.logging_config import LibephemerisFormatter

        formatter = LibephemerisFormatter()
        record = logging.LogRecord(
            name="libephemeris",
            level=logging.INFO,
            pathname="",
            lineno=0,
            msg="Test message",
            args=(),
            exc_info=None,
        )

        formatted = formatter.format(record)
        assert "INFO" in formatted
        assert "Test message" in formatted


class TestSetLogLevel:
    """Tests for set_log_level() function."""

    def test_set_log_level_changes_logger_level(self):
        """set_log_level should change the logger level."""
        from libephemeris.logging_config import get_logger, set_log_level
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()

        set_log_level(logging.DEBUG)
        assert logger.level == logging.DEBUG

        set_log_level(logging.WARNING)
        assert logger.level == logging.WARNING

        # Reset to default
        set_log_level(logging.INFO)


class TestDisableLogging:
    """Tests for disable_logging() function."""

    def test_disable_logging_silences_logger(self):
        """disable_logging should prevent all log output."""
        from libephemeris.logging_config import (
            get_logger,
            disable_logging,
            enable_logging,
        )
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()

        # Disable logging
        disable_logging()

        # The level should be above CRITICAL
        assert logger.level > logging.CRITICAL

        # Re-enable for other tests
        enable_logging()


class TestEnableLogging:
    """Tests for enable_logging() function."""

    def test_enable_logging_restores_default_level(self):
        """enable_logging should restore the default INFO level."""
        from libephemeris.logging_config import (
            get_logger,
            disable_logging,
            enable_logging,
            DEFAULT_LEVEL,
        )
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()

        disable_logging()
        enable_logging()

        assert logger.level == DEFAULT_LEVEL

    def test_enable_logging_with_custom_level(self):
        """enable_logging should accept a custom level."""
        from libephemeris.logging_config import get_logger, enable_logging
        import libephemeris.logging_config as lc

        lc._logger_configured = False

        logger = get_logger()

        enable_logging(logging.DEBUG)
        assert logger.level == logging.DEBUG

        # Reset to default
        enable_logging()


class TestFormatFileSize:
    """Tests for format_file_size() function."""

    def test_format_bytes(self):
        """format_file_size should format bytes correctly."""
        from libephemeris.logging_config import format_file_size

        assert format_file_size(0) == "0.0 B"
        assert format_file_size(512) == "512.0 B"

    def test_format_kilobytes(self):
        """format_file_size should format kilobytes correctly."""
        from libephemeris.logging_config import format_file_size

        assert format_file_size(1024) == "1.0 KB"
        assert format_file_size(2048) == "2.0 KB"

    def test_format_megabytes(self):
        """format_file_size should format megabytes correctly."""
        from libephemeris.logging_config import format_file_size

        assert format_file_size(1024 * 1024) == "1.0 MB"
        # DE440 is about 114 MB
        result = format_file_size(119537664)
        assert "MB" in result
        assert float(result.split()[0]) > 100

    def test_format_gigabytes(self):
        """format_file_size should format gigabytes correctly."""
        from libephemeris.logging_config import format_file_size

        assert format_file_size(1024 * 1024 * 1024) == "1.0 GB"


class TestLoggingExportedFromInit:
    """Tests that logging functions are exported from libephemeris.__init__."""

    def test_get_logger_exported(self):
        """get_logger should be available from libephemeris."""
        import libephemeris

        assert hasattr(libephemeris, "get_logger")

    def test_set_log_level_exported(self):
        """set_log_level should be available from libephemeris."""
        import libephemeris

        assert hasattr(libephemeris, "set_log_level")

    def test_disable_logging_exported(self):
        """disable_logging should be available from libephemeris."""
        import libephemeris

        assert hasattr(libephemeris, "disable_logging")

    def test_enable_logging_exported(self):
        """enable_logging should be available from libephemeris."""
        import libephemeris

        assert hasattr(libephemeris, "enable_logging")

    def test_format_file_size_exported(self):
        """format_file_size should be available from libephemeris."""
        import libephemeris

        assert hasattr(libephemeris, "format_file_size")


class TestLoggerNameConstant:
    """Tests for the LOGGER_NAME constant."""

    def test_logger_name_is_libephemeris(self):
        """LOGGER_NAME should be 'libephemeris'."""
        from libephemeris.logging_config import LOGGER_NAME

        assert LOGGER_NAME == "libephemeris"


class TestGetPlanetsLogging:
    """Tests for logging in get_planets() function."""

    def test_logs_info_when_loading_cached_ephemeris(self, caplog, tmp_path):
        """get_planets should log INFO when loading cached ephemeris."""
        import libephemeris.state as state
        from libephemeris.logging_config import get_logger, LOGGER_NAME
        from unittest.mock import MagicMock, patch

        # Create a fake ephemeris file
        fake_bsp = tmp_path / "de440.bsp"
        fake_bsp.touch()

        # Reset state
        state._PLANETS = None
        state._EPHEMERIS_FILE = "de440.bsp"
        original_path = state._EPHEMERIS_PATH
        state._EPHEMERIS_PATH = str(tmp_path)

        # Mock the loader to avoid actually loading the file
        mock_kernel = MagicMock()
        mock_loader = MagicMock(return_value=mock_kernel)

        # Get the logger and temporarily set propagate=True for caplog to capture
        logger = get_logger()
        original_propagate = logger.propagate
        logger.propagate = True

        try:
            with patch.object(state, "get_loader", return_value=mock_loader):
                with caplog.at_level(logging.DEBUG, logger=LOGGER_NAME):
                    state.get_planets()

            # Check INFO log for "Ephemeris loaded"
            info_messages = [
                r.message for r in caplog.records if r.levelno == logging.INFO
            ]
            assert any("Ephemeris loaded:" in msg for msg in info_messages)

            # Check DEBUG log for "Using cached ephemeris"
            debug_messages = [
                r.message for r in caplog.records if r.levelno == logging.DEBUG
            ]
            assert any("Using cached ephemeris:" in msg for msg in debug_messages)
        finally:
            # Restore state
            logger.propagate = original_propagate
            state._EPHEMERIS_PATH = original_path
            state._PLANETS = None

    def test_logs_info_when_downloading_ephemeris(self, caplog):
        """get_planets should log INFO when downloading ephemeris."""
        import libephemeris.state as state
        from libephemeris.logging_config import get_logger, LOGGER_NAME
        from unittest.mock import MagicMock, patch

        # Reset state
        state._PLANETS = None
        state._EPHEMERIS_FILE = "de440.bsp"
        original_path = state._EPHEMERIS_PATH
        state._EPHEMERIS_PATH = None

        # Mock the loader to simulate download
        mock_kernel = MagicMock()
        mock_loader = MagicMock(return_value=mock_kernel)

        # Get the logger and temporarily set propagate=True for caplog to capture
        logger = get_logger()
        original_propagate = logger.propagate

        try:
            logger.propagate = True
            with patch.object(state, "get_loader", return_value=mock_loader):
                # Mock os.path.exists to return False for all paths (force download path)
                # Force non-LEB calc mode so we hit the standard download branch
                with (
                    patch("os.path.exists", return_value=False),
                    patch.object(state, "get_calc_mode", return_value="auto"),
                ):
                    with caplog.at_level(logging.INFO, logger=LOGGER_NAME):
                        state.get_planets()

            # Check that download message was logged
            info_messages = [
                r.message for r in caplog.records if r.levelno == logging.INFO
            ]
            assert any("Downloading JPL ephemeris" in msg for msg in info_messages)
            assert any("Ephemeris downloaded to:" in msg for msg in info_messages)
        finally:
            # Restore state
            logger.propagate = original_propagate
            state._EPHEMERIS_PATH = original_path
            state._PLANETS = None

    def test_log_message_format_matches_expected(self, caplog):
        """Log messages should match the expected format for user visibility."""
        import libephemeris.state as state
        from libephemeris.logging_config import get_logger, LOGGER_NAME
        from unittest.mock import MagicMock, patch

        # Reset state
        state._PLANETS = None
        state._EPHEMERIS_FILE = "de440.bsp"
        original_path = state._EPHEMERIS_PATH
        state._EPHEMERIS_PATH = None

        # Mock the loader
        mock_kernel = MagicMock()
        mock_loader = MagicMock(return_value=mock_kernel)

        # Get the logger and temporarily set propagate=True for caplog to capture
        logger = get_logger()
        original_propagate = logger.propagate

        try:
            logger.propagate = True
            with patch.object(state, "get_loader", return_value=mock_loader):
                # Mock os.path.exists to return False (force download path)
                # Force non-LEB calc mode so we hit the standard download branch
                with (
                    patch("os.path.exists", return_value=False),
                    patch.object(state, "get_calc_mode", return_value="auto"),
                ):
                    with caplog.at_level(logging.INFO, logger=LOGGER_NAME):
                        state.get_planets()

            # Verify the download message format
            download_msgs = [
                r.message
                for r in caplog.records
                if "Downloading JPL ephemeris" in r.message
            ]
            assert len(download_msgs) >= 1
            assert "de440.bsp" in download_msgs[0]
        finally:
            # Restore state
            logger.propagate = original_propagate
            state._EPHEMERIS_PATH = original_path
            state._PLANETS = None


class TestLibephemerisLogLevelEnvVar:
    """Tests for LIBEPHEMERIS_LOG_LEVEL environment variable support."""

    def test_env_var_constant_exported(self):
        """LIBEPHEMERIS_LOG_LEVEL_ENV should be exported from libephemeris."""
        import libephemeris

        assert hasattr(libephemeris, "LIBEPHEMERIS_LOG_LEVEL_ENV")
        assert libephemeris.LIBEPHEMERIS_LOG_LEVEL_ENV == "LIBEPHEMERIS_LOG_LEVEL"

    def test_get_log_level_from_env_default_is_warning(self):
        """Default log level should be WARNING when env var is not set."""
        from libephemeris.logging_config import _get_log_level_from_env

        with patch.dict(os.environ, {}, clear=False):
            # Remove the env var if it exists
            os.environ.pop("LIBEPHEMERIS_LOG_LEVEL", None)
            level = _get_log_level_from_env()
            assert level == logging.WARNING

    def test_get_log_level_from_env_debug(self):
        """LIBEPHEMERIS_LOG_LEVEL=DEBUG should return logging.DEBUG."""
        from libephemeris.logging_config import _get_log_level_from_env

        with patch.dict(os.environ, {"LIBEPHEMERIS_LOG_LEVEL": "DEBUG"}):
            level = _get_log_level_from_env()
            assert level == logging.DEBUG

    def test_get_log_level_from_env_info(self):
        """LIBEPHEMERIS_LOG_LEVEL=INFO should return logging.INFO."""
        from libephemeris.logging_config import _get_log_level_from_env

        with patch.dict(os.environ, {"LIBEPHEMERIS_LOG_LEVEL": "INFO"}):
            level = _get_log_level_from_env()
            assert level == logging.INFO

    def test_get_log_level_from_env_warning(self):
        """LIBEPHEMERIS_LOG_LEVEL=WARNING should return logging.WARNING."""
        from libephemeris.logging_config import _get_log_level_from_env

        with patch.dict(os.environ, {"LIBEPHEMERIS_LOG_LEVEL": "WARNING"}):
            level = _get_log_level_from_env()
            assert level == logging.WARNING

    def test_get_log_level_from_env_error(self):
        """LIBEPHEMERIS_LOG_LEVEL=ERROR should return logging.ERROR."""
        from libephemeris.logging_config import _get_log_level_from_env

        with patch.dict(os.environ, {"LIBEPHEMERIS_LOG_LEVEL": "ERROR"}):
            level = _get_log_level_from_env()
            assert level == logging.ERROR

    def test_get_log_level_from_env_critical(self):
        """LIBEPHEMERIS_LOG_LEVEL=CRITICAL should return logging.CRITICAL."""
        from libephemeris.logging_config import _get_log_level_from_env

        with patch.dict(os.environ, {"LIBEPHEMERIS_LOG_LEVEL": "CRITICAL"}):
            level = _get_log_level_from_env()
            assert level == logging.CRITICAL

    def test_get_log_level_from_env_case_insensitive(self):
        """LIBEPHEMERIS_LOG_LEVEL should be case-insensitive."""
        from libephemeris.logging_config import _get_log_level_from_env

        # Test lowercase
        with patch.dict(os.environ, {"LIBEPHEMERIS_LOG_LEVEL": "debug"}):
            level = _get_log_level_from_env()
            assert level == logging.DEBUG

        # Test mixed case
        with patch.dict(os.environ, {"LIBEPHEMERIS_LOG_LEVEL": "Info"}):
            level = _get_log_level_from_env()
            assert level == logging.INFO

    def test_get_log_level_from_env_strips_whitespace(self):
        """LIBEPHEMERIS_LOG_LEVEL should handle leading/trailing whitespace."""
        from libephemeris.logging_config import _get_log_level_from_env

        with patch.dict(os.environ, {"LIBEPHEMERIS_LOG_LEVEL": "  DEBUG  "}):
            level = _get_log_level_from_env()
            assert level == logging.DEBUG

    def test_get_log_level_from_env_invalid_value_falls_back_to_warning(self):
        """Invalid LIBEPHEMERIS_LOG_LEVEL values should fall back to WARNING."""
        from libephemeris.logging_config import _get_log_level_from_env

        with patch.dict(os.environ, {"LIBEPHEMERIS_LOG_LEVEL": "INVALID"}):
            level = _get_log_level_from_env()
            assert level == logging.WARNING

        with patch.dict(os.environ, {"LIBEPHEMERIS_LOG_LEVEL": "TRACE"}):
            level = _get_log_level_from_env()
            assert level == logging.WARNING

    def test_get_log_level_from_env_empty_string_falls_back_to_warning(self):
        """Empty LIBEPHEMERIS_LOG_LEVEL should fall back to WARNING."""
        from libephemeris.logging_config import _get_log_level_from_env

        with patch.dict(os.environ, {"LIBEPHEMERIS_LOG_LEVEL": ""}):
            level = _get_log_level_from_env()
            assert level == logging.WARNING

    def test_valid_log_levels_constant(self):
        """_VALID_LOG_LEVELS should contain all standard logging levels."""
        from libephemeris.logging_config import _VALID_LOG_LEVELS

        assert "DEBUG" in _VALID_LOG_LEVELS
        assert "INFO" in _VALID_LOG_LEVELS
        assert "WARNING" in _VALID_LOG_LEVELS
        assert "ERROR" in _VALID_LOG_LEVELS
        assert "CRITICAL" in _VALID_LOG_LEVELS
        # NOTSET should not be valid (too permissive)
        assert "NOTSET" not in _VALID_LOG_LEVELS


class TestLoggerConfiguredWithEnvVar:
    """Tests that logger uses env var on configuration."""

    def test_logger_respects_env_var_at_configuration(self):
        """Logger should use LIBEPHEMERIS_LOG_LEVEL at configuration time."""
        import libephemeris.logging_config as lc

        # Store original state
        original_configured = lc._logger_configured
        original_default_level = lc.DEFAULT_LEVEL

        try:
            # Reset logger state
            lc._logger_configured = False
            logger = logging.getLogger(lc.LOGGER_NAME)
            logger.handlers.clear()

            # Override DEFAULT_LEVEL to simulate env var being DEBUG
            lc.DEFAULT_LEVEL = logging.DEBUG

            # Configure the logger
            configured_logger = lc._configure_logger()

            assert configured_logger.level == logging.DEBUG
        finally:
            # Restore original state
            lc._logger_configured = original_configured
            lc.DEFAULT_LEVEL = original_default_level


class TestEnvVarIntegration:
    """Integration tests for env var behavior."""

    def test_debug_level_shows_debug_messages(self, caplog):
        """DEBUG level should show debug-level log messages."""
        from libephemeris.logging_config import get_logger, LOGGER_NAME, set_log_level

        logger = get_logger()
        original_level = logger.level
        original_propagate = logger.propagate

        try:
            logger.propagate = True
            set_log_level(logging.DEBUG)

            with caplog.at_level(logging.DEBUG, logger=LOGGER_NAME):
                logger.debug("Debug message for testing")
                logger.info("Info message for testing")

            debug_messages = [
                r.message for r in caplog.records if r.levelno == logging.DEBUG
            ]
            info_messages = [
                r.message for r in caplog.records if r.levelno == logging.INFO
            ]

            assert any("Debug message" in msg for msg in debug_messages)
            assert any("Info message" in msg for msg in info_messages)
        finally:
            logger.propagate = original_propagate
            set_log_level(original_level)

    def test_warning_level_hides_info_messages(self):
        """WARNING level should filter INFO and DEBUG messages from handlers."""
        from libephemeris.logging_config import get_logger, set_log_level

        logger = get_logger()
        original_level = logger.level

        try:
            set_log_level(logging.WARNING)

            # Verify the logger is configured at WARNING level
            assert logger.level == logging.WARNING

            # Verify that messages below WARNING would be filtered
            # isEnabledFor() checks if a message of the given level would be emitted
            assert not logger.isEnabledFor(logging.DEBUG)
            assert not logger.isEnabledFor(logging.INFO)
            assert logger.isEnabledFor(logging.WARNING)
            assert logger.isEnabledFor(logging.ERROR)
        finally:
            set_log_level(original_level)


class TestLogLevelTomlFallback:
    """Cover the TOML fallback branch in _get_log_level_from_env()."""

    def test_toml_value_used_when_env_absent(self, monkeypatch):
        import libephemeris._config_toml as toml_mod
        from libephemeris import logging_config

        monkeypatch.delenv(logging_config.LIBEPHEMERIS_LOG_LEVEL_ENV, raising=False)
        monkeypatch.setattr(toml_mod, "get_str", lambda _key: "debug")

        assert logging_config._get_log_level_from_env() == logging.DEBUG

    def test_toml_lookup_failure_falls_back_to_warning(self, monkeypatch):
        import libephemeris._config_toml as toml_mod
        from libephemeris import logging_config

        def _boom(_key):
            raise RuntimeError("toml read failed")

        monkeypatch.delenv(logging_config.LIBEPHEMERIS_LOG_LEVEL_ENV, raising=False)
        monkeypatch.setattr(toml_mod, "get_str", _boom)

        assert logging_config._get_log_level_from_env() == logging.WARNING


class TestFormatFileSizePetabyte:
    """Cover format_file_size petabyte tail."""

    def test_petabyte_range(self):
        from libephemeris.logging_config import format_file_size

        result = format_file_size(2 * 1024**5)
        assert result.endswith("PB")
        assert result.startswith("2.0")
