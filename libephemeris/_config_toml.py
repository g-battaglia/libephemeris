# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""TOML configuration file loader for libephemeris.

Loads configuration from ``libephemeris-config.toml`` files,
providing a structured, version-controllable alternative to ``.env``
files for projects that use libephemeris as a dependency.

Search order (first file found wins):

1. ``LIBEPHEMERIS_CONFIG`` environment variable (if set).
2. ``./libephemeris-config.toml`` -- per-project configuration.
3. ``~/.libephemeris/config.toml`` -- global user configuration.

The file must contain a ``[libephemeris]`` section.  Example::

    [libephemeris]
    precision = "extended"
    mode = "leb"
    auto_spk = true
    strict_precision = true

Resolution order (highest to lowest priority):

1. Programmatic override via ``set_*()`` functions
2. Environment variables (including ``.env`` file)
3. **TOML configuration file** (this module)
4. Hardcoded defaults
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Union

# TOML parser: stdlib on 3.11+, tomli package on 3.9/3.10
if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - Python <3.11 tomli shim; project requires 3.12+
    try:
        import tomli as tomllib  # type: ignore[no-redef]
    except ImportError:
        tomllib = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Environment variable to override the config file path
CONFIG_FILE_VAR = "LIBEPHEMERIS_CONFIG"

# Default config file names
_PROJECT_CONFIG_NAME = "libephemeris-config.toml"
_GLOBAL_CONFIG_NAME = "config.toml"

# Default search locations (resolved lazily, in priority order)
_SEARCH_PATHS = (
    lambda: Path.cwd() / _PROJECT_CONFIG_NAME,
    lambda: Path.home() / ".libephemeris" / _GLOBAL_CONFIG_NAME,
)

# Section name inside the TOML file
_SECTION = "libephemeris"

# ---------------------------------------------------------------------------
# Cached state (populated once by ``load_config()``)
# ---------------------------------------------------------------------------

_CONFIG: Dict[str, Any] = {}
_CONFIG_PATH: Optional[str] = None
_CONFIG_LOADED: bool = False

# Valid keys and their expected Python types (for validation)
_VALID_KEYS: Dict[str, type] = {
    "precision": str,
    "mode": str,
    "ephemeris": str,
    "leb_file": str,
    "data_dir": str,
    "log_level": str,
    "auto_spk": bool,
    "spk_dir": str,
    "strict_precision": bool,
    "iers_auto_download": bool,
    "iers_delta_t": bool,
    "deltat_model": str,
    "mmap_preload": bool,
    "mmap_preload_start": int,
    "mmap_preload_end": int,
}


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


def _find_config_file() -> Optional[Path]:
    """Locate the TOML config file using the default search order.

    Returns the first existing file, or ``None`` if no file is found.
    """
    # 1. Explicit override via environment variable
    env_override = os.environ.get(CONFIG_FILE_VAR, "").strip()
    if env_override:
        p = Path(env_override)
        return p if p.is_file() else None

    # 2. Search default locations
    for path_fn in _SEARCH_PATHS:
        try:
            path = path_fn()
            if path.is_file():
                return path
        except (OSError, RuntimeError):
            continue

    return None


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def load_config(
    path: Optional[Union[str, Path]] = None,
) -> bool:
    """Load configuration from a TOML file.

    Args:
        path: Explicit path to a TOML config file.  If ``None``, the
            default search order is used (see module docstring).

    Returns:
        ``True`` if a config file was found and loaded, ``False`` otherwise.
    """
    global _CONFIG, _CONFIG_PATH, _CONFIG_LOADED

    if tomllib is None:
        # No TOML parser available (Python <3.11 without tomli installed)
        _CONFIG_LOADED = True
        return False

    config_path: Optional[Path]
    if path is not None:
        config_path = Path(path)
        if not config_path.is_file():
            _CONFIG_LOADED = True
            return False
    else:
        config_path = _find_config_file()
        if config_path is None:
            _CONFIG_LOADED = True
            return False

    try:
        with open(config_path, "rb") as f:
            data = tomllib.load(f)
    except OSError:
        _CONFIG_LOADED = True
        return False
    except tomllib.TOMLDecodeError as e:
        # A malformed config silently ignored looks like "my settings
        # don't apply" — say what is wrong and where.
        import warnings

        warnings.warn(
            f"Ignoring malformed libephemeris config {config_path}: {e}",
            stacklevel=2,
        )
        _CONFIG_LOADED = True
        return False

    section = data.get(_SECTION, {})
    if not isinstance(section, dict):
        _CONFIG_LOADED = True
        return False

    # Validate and store only known keys with correct types
    validated: Dict[str, Any] = {}
    for key, expected_type in _VALID_KEYS.items():
        if key in section:
            value = section[key]
            if isinstance(value, expected_type):
                validated[key] = value

    _CONFIG = validated
    _CONFIG_PATH = str(config_path)
    _CONFIG_LOADED = True
    return True


# ---------------------------------------------------------------------------
# Typed accessors (called by getter functions in state.py et al.)
# ---------------------------------------------------------------------------


def _ensure_loaded() -> None:
    """Ensure the config has been loaded at least once."""
    if not _CONFIG_LOADED:
        load_config()


def get_str(key: str) -> Optional[str]:
    """Get a string value from the loaded TOML config.

    Returns ``None`` if the key is not set or is not a string.
    """
    _ensure_loaded()
    value = _CONFIG.get(key)
    return value if isinstance(value, str) else None


def get_bool(key: str) -> Optional[bool]:
    """Get a boolean value from the loaded TOML config.

    Returns ``None`` if the key is not set or is not a boolean.
    """
    _ensure_loaded()
    value = _CONFIG.get(key)
    return value if isinstance(value, bool) else None


def get_int(key: str) -> Optional[int]:
    """Get an integer value from the loaded TOML config.

    Returns ``None`` if the key is not set or is not an integer.
    """
    _ensure_loaded()
    value = _CONFIG.get(key)
    return value if isinstance(value, int) and not isinstance(value, bool) else None


# ---------------------------------------------------------------------------
# Introspection (used by ``libephemeris config`` and ``libephemeris status``)
# ---------------------------------------------------------------------------


def get_config_path() -> Optional[str]:
    """Return the path of the loaded config file, or ``None``."""
    _ensure_loaded()
    return _CONFIG_PATH


def get_all() -> Dict[str, Any]:
    """Return a copy of all loaded config values."""
    _ensure_loaded()
    return dict(_CONFIG)


def is_loaded() -> bool:
    """Return ``True`` if a TOML config file was successfully loaded."""
    _ensure_loaded()
    return bool(_CONFIG)


# ---------------------------------------------------------------------------
# Reset (for testing and ``close()``)
# ---------------------------------------------------------------------------


def reset() -> None:
    """Reset the cached config state.

    After calling this, the next accessor call will re-discover and
    re-load the config file.
    """
    global _CONFIG, _CONFIG_PATH, _CONFIG_LOADED
    _CONFIG = {}
    _CONFIG_PATH = None
    _CONFIG_LOADED = False
