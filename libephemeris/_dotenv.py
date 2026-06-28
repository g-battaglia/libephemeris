# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Minimal .env file loader using only the Python standard library.

Provides automatic loading of environment variables from ``.env`` files
at package import time, with no external dependencies.

Search order (first file found wins):

1. Explicit *path* argument (if provided).
2. ``LIBEPHEMERIS_ENV_FILE`` environment variable (if set).
3. ``./.env`` -- per-project configuration.
4. ``~/.libephemeris/.env`` -- global user configuration.

Syntax supported::

    # comment
    KEY=value
    KEY="value with spaces"
    KEY='value with spaces'
    export KEY=value
    KEY=value  # inline comment (unquoted only)

Existing environment variables are **not** overwritten unless
*override=True* is passed.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional, Tuple, Union


# Environment variable to override the .env file path
_ENV_FILE_VAR = "LIBEPHEMERIS_ENV_FILE"

# Default search locations (resolved lazily, in priority order)
_SEARCH_PATHS = (
    lambda: Path.cwd() / ".env",
    lambda: Path.home() / ".libephemeris" / ".env",
)


def _parse_line(line: str) -> Optional[Tuple[str, str]]:
    """Parse a single .env line into a (key, value) pair.

    Returns ``None`` for blank lines and comments.
    """
    line = line.strip()

    # Skip empty lines and comments
    if not line or line.startswith("#"):
        return None

    # Strip optional 'export ' prefix
    if line.startswith("export "):
        line = line[7:].strip()

    # Must contain '='
    if "=" not in line:
        return None

    key, _, value = line.partition("=")
    key = key.strip()

    if not key:
        return None

    value = value.strip()

    # Remove matching quotes (double or single)
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ('"', "'"):
        value = value[1:-1]
    else:
        # Unquoted: strip inline comments (KEY=value # comment)
        comment_idx = value.find(" #")
        if comment_idx != -1:
            value = value[:comment_idx].rstrip()

    return key, value


def _find_env_file() -> Optional[Path]:
    """Locate the .env file using the default search order.

    Returns the first existing file, or ``None`` if no file is found.
    """
    # 1. Explicit override via environment variable
    env_override = os.environ.get(_ENV_FILE_VAR, "").strip()
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


def load_dotenv(
    path: Optional[Union[str, Path]] = None,
    *,
    override: bool = False,
) -> bool:
    """Load variables from a ``.env`` file into ``os.environ``.

    Args:
        path: Explicit path to a ``.env`` file.  If ``None``, the default
            search order is used (see module docstring).
        override: If ``True``, overwrite existing environment variables.
            Defaults to ``False`` (existing variables take precedence).

    Returns:
        ``True`` if a ``.env`` file was found and loaded, ``False`` otherwise.

    Example::

        >>> from libephemeris._dotenv import load_dotenv
        >>> load_dotenv("/path/to/my/.env")
        True
    """
    env_path: Optional[Path]
    if path is not None:
        env_path = Path(path)
        if not env_path.is_file():
            return False
    else:
        env_path = _find_env_file()
        if env_path is None:
            return False

    try:
        text = env_path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False

    for line in text.splitlines():
        parsed = _parse_line(line)
        if parsed is None:
            continue

        key, value = parsed

        # Only the library's own namespace: a project .env regularly
        # holds unrelated secrets (API keys, DB urls) that must not be
        # exported into the process just because libephemeris imported.
        if not key.startswith("LIBEPHEMERIS_"):
            continue

        if override or key not in os.environ:
            os.environ[key] = value

    return True
