# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Centralized logging configuration for libephemeris.

This module provides a dedicated logger for the libephemeris library,
with a stderr handler for user-visible messages about downloads and
other long-running operations.

Usage:
    from libephemeris.logging_config import get_logger

    logger = get_logger()
    logger.info("Downloading DE440 ephemeris (114 MB)...")

The logger outputs to stderr with the format:
    [libephemeris] INFO: Downloading DE440 ephemeris (114 MB)...

Log Level Configuration:
    The log level can be configured via the LIBEPHEMERIS_LOG_LEVEL environment
    variable. Valid values: DEBUG, INFO, WARNING, ERROR, CRITICAL.
    Default is WARNING for quiet production operation.

    Examples:
        LIBEPHEMERIS_LOG_LEVEL=DEBUG pytest -s    # Show all debug messages
        LIBEPHEMERIS_LOG_LEVEL=INFO python app.py # Show download progress
        LIBEPHEMERIS_LOG_LEVEL=ERROR              # Only errors

    Programmatic configuration:
        import logging
        logging.getLogger("libephemeris").setLevel(logging.DEBUG)

    Or disable logging entirely:
        logging.getLogger("libephemeris").setLevel(logging.CRITICAL + 1)

DEBUG-Level Source Tracing:
    At DEBUG level, libephemeris logs which calculation backend was used for
    every celestial body at every dispatch point. The log format is:

        body=<id> jd=<julian_day> source=<SOURCE>

    Possible source values: LEB, Skyfield, Horizons, SPK, ASSIST (n-body),
    Keplerian (fallback). See docs/development/testing.md for details.
"""

from __future__ import annotations

import logging
import os
import sys

# Module-level logger name
LOGGER_NAME = "libephemeris"

# Environment variable name for log level configuration
LIBEPHEMERIS_LOG_LEVEL_ENV = "LIBEPHEMERIS_LOG_LEVEL"

# Valid log level names
_VALID_LOG_LEVELS = frozenset({"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"})


def _get_log_level_from_env() -> int:
    """Get log level from environment variable.

    Reads the LIBEPHEMERIS_LOG_LEVEL environment variable and returns
    the corresponding logging level. Falls back to WARNING if the
    variable is not set or has an invalid value.

    Returns:
        int: The logging level (e.g., logging.DEBUG, logging.WARNING).

    Example:
        >>> os.environ["LIBEPHEMERIS_LOG_LEVEL"] = "DEBUG"
        >>> _get_log_level_from_env()
        10  # logging.DEBUG
    """
    env_level = os.environ.get(LIBEPHEMERIS_LOG_LEVEL_ENV, "").upper().strip()

    if env_level and env_level in _VALID_LOG_LEVELS:
        return getattr(logging, env_level)

    # TOML config fallback
    try:
        from ._config_toml import get_str as _toml_str

        toml_value = _toml_str("log_level")
        if toml_value and toml_value.upper().strip() in _VALID_LOG_LEVELS:
            return getattr(logging, toml_value.upper().strip())
    except Exception:
        pass

    # Default for production quietness
    return logging.WARNING


# Default logging level (from env var or WARNING)
DEFAULT_LEVEL = _get_log_level_from_env()

# Flag to track if the logger has been configured
_logger_configured = False


class LibephemerisFormatter(logging.Formatter):
    """Custom formatter for libephemeris log messages.

    Formats messages as: [libephemeris] LEVEL: message
    """

    def format(self, record: logging.LogRecord) -> str:
        """Format the log record."""
        return f"[libephemeris] {record.levelname}: {record.getMessage()}"


def _configure_logger() -> logging.Logger:
    """
    Configure and return the libephemeris logger.

    This function sets up the logger with a stderr handler if it hasn't
    been configured yet. It's safe to call multiple times.

    Returns:
        The configured logger instance.
    """
    global _logger_configured

    logger = logging.getLogger(LOGGER_NAME)

    # Only configure once to avoid duplicate handlers
    if not _logger_configured:
        # Set the logger level
        logger.setLevel(DEFAULT_LEVEL)

        # Create stderr handler
        handler = logging.StreamHandler(sys.stderr)
        handler.setLevel(DEFAULT_LEVEL)

        # Set custom formatter
        formatter = LibephemerisFormatter()
        handler.setFormatter(formatter)

        # Add handler to logger
        logger.addHandler(handler)

        # Prevent propagation to root logger to avoid duplicate messages
        logger.propagate = False

        _logger_configured = True

    return logger


def get_logger() -> logging.Logger:
    """
    Get the libephemeris logger instance.

    Returns the configured logger for libephemeris. The logger is
    automatically configured on first call with a stderr handler.

    Returns:
        logging.Logger: The libephemeris logger instance.

    Example:
        >>> from libephemeris.logging_config import get_logger
        >>> logger = get_logger()
        >>> logger.info("Starting download...")
        [libephemeris] INFO: Starting download...
    """
    return _configure_logger()


def set_log_level(level: int) -> None:
    """
    Set the logging level for libephemeris.

    This is a convenience function to change the logging level
    without directly accessing the logger.

    Args:
        level: The logging level (e.g., logging.DEBUG, logging.INFO,
               logging.WARNING, logging.ERROR, logging.CRITICAL)

    Example:
        >>> from libephemeris.logging_config import set_log_level
        >>> import logging
        >>> set_log_level(logging.DEBUG)  # Enable debug messages
        >>> set_log_level(logging.WARNING)  # Only warnings and above
    """
    logger = get_logger()
    logger.setLevel(level)
    for handler in logger.handlers:
        handler.setLevel(level)


def disable_logging() -> None:
    """
    Disable all libephemeris logging.

    This completely silences all log messages from libephemeris.
    Useful for testing or when embedding libephemeris in other
    applications that have their own logging setup.

    Example:
        >>> from libephemeris.logging_config import disable_logging
        >>> disable_logging()  # All logging is now silenced
    """
    logger = get_logger()
    logger.setLevel(logging.CRITICAL + 1)


def enable_logging(level: int = DEFAULT_LEVEL) -> None:
    """
    Enable libephemeris logging at the specified level.

    Re-enables logging after it has been disabled, or changes
    the logging level.

    Args:
        level: The logging level to use (default: INFO)

    Example:
        >>> from libephemeris.logging_config import enable_logging
        >>> import logging
        >>> enable_logging()  # Re-enable at INFO level
        >>> enable_logging(logging.DEBUG)  # Enable at DEBUG level
    """
    set_log_level(level)


def format_file_size(size_bytes: int) -> str:
    """
    Format a file size in bytes to a human-readable string.

    Args:
        size_bytes: Size in bytes

    Returns:
        Human-readable size string (e.g., "114 MB", "1.5 GB")

    Example:
        >>> format_file_size(119537664)
        '114.0 MB'
    """
    size = float(size_bytes)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} PB"
