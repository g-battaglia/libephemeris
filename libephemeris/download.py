# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Data file download utilities for libephemeris.

This module provides functionality to download optional data files that enhance
precision for certain calculations. The main file is planet_centers.bsp which
provides precise planet center positions for outer planets (Jupiter-Pluto).

Usage:
    # From command line
    libephemeris download medium       Download data for 'medium' tier
    libephemeris download base         Download data for 'base' tier
    libephemeris download extended     Download data for 'extended' tier

    # From Python
    from libephemeris.download import download_for_tier
    download_for_tier("medium")

Without the data files, libephemeris will fall back to analytical approximations
which are still accurate to ~0.1 arcseconds.
"""

from __future__ import annotations

import hashlib
import http.client
import os
import ssl
import sys
import tempfile
import urllib.request
from pathlib import Path
from typing import Any, Optional

import certifi

from .leb_groups import LEB2_GROUPS as _CANONICAL_LEB2_GROUPS
from .logging_config import get_logger


def _is_valid_bsp(filepath: str) -> bool:
    """Check if a BSP file can be opened and its segments enumerated by jplephem.

    Uses jplephem.spk.SPK.open() which validates:
    - DAF header (magic bytes, endianness)
    - FTP corruption test string
    - All segment summary records are readable

    This catches truncated files where the header is intact but segment
    records are missing deeper in the file.

    Args:
        filepath: Path to the BSP file to validate

    Returns:
        True if file can be opened and all segments enumerated, False otherwise
    """
    try:
        from jplephem.spk import SPK

        with SPK.open(filepath) as kernel:
            # Force iteration through all segments to detect truncation
            for segment in kernel.segments:
                _ = segment.center, segment.target
        return True
    except (OSError, ValueError, KeyError, RuntimeError):
        return False


# GitHub Releases URL for data files (single unified release)
GITHUB_RELEASES = (
    "https://github.com/g-battaglia/libephemeris/releases/download/data-v1"
)

# Data file definitions: (filename, sha256 hash, description)
DATA_FILES: dict[str, dict[str, Any]] = {
    # Legacy single file (kept for backward compatibility)
    "planet_centers.bsp": {
        "url": f"{GITHUB_RELEASES}/planet_centers.bsp",
        "sha256": None,  # Will be set after first release
        "size_mb": 25.4,
        "description": "Precise planet center positions for Jupiter, Saturn, Uranus, Neptune, Pluto (1989-2049)",
    },
    # Tier-specific files
    "planet_centers_base.bsp": {
        "url": f"{GITHUB_RELEASES}/planet_centers_base.bsp",
        "sha256": "a9ec744ff412b095129166587ea0814f81c850faebf92586a738cb5dc103c92a",
        "size_mb": 25.4,
        "description": "Planet centers for 'base' tier (1850-2150)",
    },
    "planet_centers_medium.bsp": {
        "url": f"{GITHUB_RELEASES}/planet_centers_medium.bsp",
        "sha256": "b4fd366f2d00958ee3dd4a8884164d339674d5f7f2ea25c5f6959705c2b66852",
        "size_mb": 72.6,
        "description": "Planet centers for 'medium' tier (1550-2650)",
    },
    "planet_centers_extended.bsp": {
        "url": f"{GITHUB_RELEASES}/planet_centers_extended.bsp",
        "sha256": "a07b046b89a9992fc7fda445b00e656341a3bab66a035adb8108de7d4bd69edc",
        "size_mb": 222.6,
        "description": "Planet centers for 'extended' tier (partial -12000 to +17000)",
    },
    # LEB (LibEphemeris Binary) precomputed ephemeris files
    "ephemeris_base.leb": {
        "url": f"{GITHUB_RELEASES}/ephemeris_base.leb",
        "sha256": "006073f4c1b7926b94fb9137322dc5dd0939d079a12fd81b42e771c5a1c9cb61",
        "size_mb": 53.1,
        "description": "LEB binary ephemeris for 'base' tier (1850-2150, ~14x speedup)",
    },
    "ephemeris_medium.leb": {
        "url": f"{GITHUB_RELEASES}/ephemeris_medium.leb",
        "sha256": "dbd7239ac2aac96ce2dc33c322a5ff9bd9fd7b3d7ec5792cfba3250e6ab3b665",
        "size_mb": 174.6,
        "description": "LEB binary ephemeris for 'medium' tier (1550-2650, ~14x speedup)",
    },
    "ephemeris_extended.leb": {
        "url": f"{GITHUB_RELEASES}/ephemeris_extended.leb",
        "sha256": "8f3d9ca0efac7c8616041c96efa63fa0c64f375e5a4671cb37c70f98f0193a27",
        "size_mb": 1603.6,
        "description": "LEB binary ephemeris for 'extended' tier (-5000 to +5000, ~14x speedup)",
    },
    # LEB2 compressed modular files
    # Core (14 bodies): Sun-Pluto, Earth, Mean/True Node, Mean Apogee
    "base_core.leb2": {
        "url": f"{GITHUB_RELEASES}/base_core.leb2",
        "sha256": "8e81f5d69aabbea6a53bc693fa0949bc546ccd1f52dddc1393a382fdf4a8a3de",
        "size_mb": 10.6,
        "description": "LEB2 core bodies for 'base' tier (1850-2150)",
        "dest_subdir": "leb",
    },
    "base_asteroids.leb2": {
        "url": f"{GITHUB_RELEASES}/base_asteroids.leb2",
        "sha256": "a653443ffe662e54782b404ff736b92e441cac1643658a3b955c032333e7ce6a",
        "size_mb": 8.7,
        "description": "LEB2 asteroids for 'base' tier (Chiron, Ceres, Pallas, Juno, Vesta)",
        "dest_subdir": "leb",
    },
    "base_apogee.leb2": {
        "url": f"{GITHUB_RELEASES}/base_apogee.leb2",
        "sha256": "abd436fd18c21650a95dee38c6e9879e69162620179be5e12062bd172ba0d657",
        "size_mb": 11.4,
        "description": "LEB2 apogee variants for 'base' tier (Oscu/Interp Apogee, Interp Perigee)",
        "dest_subdir": "leb",
    },
    "base_uranians.leb2": {
        "url": f"{GITHUB_RELEASES}/base_uranians.leb2",
        "sha256": "0ea94c9591dc4e90831aec59bf8291075fdd36ad430a01b9cadfd1d6b45f3bab",
        "size_mb": 2.1,
        "description": "LEB2 Uranian hypotheticals for 'base' tier (Cupido-Transpluto)",
        "dest_subdir": "leb",
    },
    "medium_core.leb2": {
        "url": f"{GITHUB_RELEASES}/medium_core.leb2",
        "sha256": "4655d490ed951bdfd214c0a94fc08e8113a724d99b5afb1a026400cc290e37ad",
        "size_mb": 38.3,
        "description": "LEB2 core bodies for 'medium' tier (1550-2650)",
        "dest_subdir": "leb",
    },
    "medium_asteroids.leb2": {
        "url": f"{GITHUB_RELEASES}/medium_asteroids.leb2",
        "sha256": "14502a2710ec5aab42c5bceec3d3f12045b58c9d6d59568172b3760cad53a9ce",
        "size_mb": 29.2,
        "description": "LEB2 asteroids for 'medium' tier",
        "dest_subdir": "leb",
    },
    "medium_apogee.leb2": {
        "url": f"{GITHUB_RELEASES}/medium_apogee.leb2",
        "sha256": "1c4d11ff90dd304a719a673f1a66d25f16529f5f3717154a6a739490b46bd193",
        "size_mb": 42.1,
        "description": "LEB2 apogee variants for 'medium' tier",
        "dest_subdir": "leb",
    },
    "medium_uranians.leb2": {
        "url": f"{GITHUB_RELEASES}/medium_uranians.leb2",
        "sha256": "872873994c5ab904cd86f5aebce37ab4ce53aa727e4afd8a3731b4ac21596798",
        "size_mb": 9.3,
        "description": "LEB2 Uranian hypotheticals for 'medium' tier",
        "dest_subdir": "leb",
    },
    "extended_core.leb2": {
        "url": f"{GITHUB_RELEASES}/extended_core.leb2",
        "sha256": "38e244d2cbcbb216269f5ea97316b543966d368a8f62a507be669aae95003389",
        "size_mb": 334.9,
        "description": "LEB2 core bodies for 'extended' tier (-5000 to +5000)",
        "dest_subdir": "leb",
    },
    "extended_asteroids.leb2": {
        "url": f"{GITHUB_RELEASES}/extended_asteroids.leb2",
        "sha256": "e7b9f5b51ec3321864774c4a622b190bda5a8502b1cb48e39e673da7260d35ba",
        "size_mb": 86.2,
        "description": "LEB2 asteroids for 'extended' tier",
        "dest_subdir": "leb",
    },
    "extended_apogee.leb2": {
        "url": f"{GITHUB_RELEASES}/extended_apogee.leb2",
        "sha256": "e129b31582f0666ec335e0916695bef9e9d5bdf0ccd7735d235bcaa79f74b5af",
        "size_mb": 391.7,
        "description": "LEB2 apogee variants for 'extended' tier",
        "dest_subdir": "leb",
    },
    "extended_uranians.leb2": {
        "url": f"{GITHUB_RELEASES}/extended_uranians.leb2",
        "sha256": "deda2a45ed3df0a7efec6867f9d180007c1a1979f466a552533e1675a6ca02fe",
        "size_mb": 84.0,
        "description": "LEB2 Uranian hypotheticals for 'extended' tier",
        "dest_subdir": "leb",
    },
}


def get_data_dir() -> Path:
    """Get the data directory path.

    Returns the path to ~/.libephemeris (or LIBEPHEMERIS_DATA_DIR env var).
    Creates the directory if it doesn't exist.

    Returns:
        Path to the data directory
    """
    from .state import _get_data_dir

    return Path(_get_data_dir())


def _format_size(size_bytes: float) -> str:
    """Format byte size to human readable string."""
    for unit in ["B", "KB", "MB", "GB"]:
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} TB"


class SimpleProgressBar:
    """Simple progress bar that works without external dependencies."""

    def __init__(
        self,
        total: int,
        description: str = "",
        width: int = 40,
        file=None,
    ):
        self.total = total
        self.description = description
        self.width = width
        self.current = 0
        self.file = file or sys.stdout
        self._last_percent = -1

    def update(self, amount: int = 1):
        """Update progress by amount."""
        self.current += amount
        self._render()

    def _render(self):
        """Render the progress bar."""
        if self.total == 0:
            return

        percent = int(100 * self.current / self.total)

        # Only update if percent changed (avoid too many writes)
        if percent == self._last_percent:
            return
        self._last_percent = percent

        filled = int(self.width * self.current / self.total)
        bar = "=" * filled + "-" * (self.width - filled)

        # Format: Description [====----] 50% (12.5 MB / 25.0 MB)
        current_str = _format_size(self.current)
        total_str = _format_size(self.total)

        line = (
            f"\r{self.description} [{bar}] {percent:3d}% ({current_str} / {total_str})"
        )
        self.file.write(line)
        self.file.flush()

    def close(self):
        """Finish the progress bar."""
        self.file.write("\n")
        self.file.flush()


def _get_progress_bar(total: int, description: str) -> Any:
    """Get the best available progress bar.

    Tries to use rich if available, falls back to simple progress bar.

    Args:
        total: Total size in bytes
        description: Description to show

    Returns:
        Progress bar object with update() and close() methods
    """
    # Try rich first (best experience)
    try:
        from rich.progress import (
            Progress,
            BarColumn,
            DownloadColumn,
            TransferSpeedColumn,
            TimeRemainingColumn,
        )

        class RichProgressWrapper:
            def __init__(self, total: int, description: str):
                self.progress = Progress(
                    "[progress.description]{task.description}",
                    BarColumn(),
                    "[progress.percentage]{task.percentage:>3.0f}%",
                    DownloadColumn(),
                    TransferSpeedColumn(),
                    TimeRemainingColumn(),
                )
                self.progress.start()
                self.task = self.progress.add_task(description, total=total)

            def update(self, amount: int):
                self.progress.update(self.task, advance=amount)

            def close(self):
                self.progress.stop()

        return RichProgressWrapper(total, description)
    except ImportError:
        pass

    # Fall back to simple progress bar
    return SimpleProgressBar(total, description)


def publish_temp_file(temp_path: str | Path, dest_path: str | Path) -> None:
    """Publish an atomically-downloaded temp file at its final path.

    Shared by every mkstemp-based downloader (DE kernels, LEB2 groups, IERS
    tables, SPK kernels, ASSIST data). mkstemp creates files 0600 and
    os.replace preserves the mode, which would lock other users out of
    shared data dirs (LIBEPHEMERIS_DATA_DIR, shared IERS/SPK caches) —
    restore world-readable permissions first. os.chmod is used because
    os.fchmod is Unix-only before Python 3.13.

    The chmod is best-effort: filesystems without POSIX permission support
    (some CIFS/SMB or vfat/exFAT mounts) raise on chmod, and a
    permission-cosmetic step must not discard an already-downloaded,
    hash-verified file.

    Args:
        temp_path: The fully-written temporary file to publish.
        dest_path: Final destination path (atomically replaced).

    Raises:
        OSError: If the atomic replace fails.
    """
    try:
        os.chmod(temp_path, 0o644)
    except OSError as exc:
        get_logger().debug(
            "Could not set 0644 on %s (non-POSIX filesystem?): %s", temp_path, exc
        )
    os.replace(temp_path, dest_path)


def download_file(
    url: str,
    dest_path: Path,
    description: str = "Downloading",
    expected_sha256: Optional[str] = None,
    show_progress: bool = True,
    timeout: int = 30,
) -> bool:
    """Download a file with progress bar.

    Args:
        url: URL to download from
        dest_path: Destination path for the file
        description: Description to show in progress bar
        expected_sha256: Expected SHA256 hash (optional, for verification)
        show_progress: Whether to show progress bar
        timeout: Network timeout in seconds

    Returns:
        True if download successful, False otherwise

    Raises:
        urllib.error.URLError: If download fails
        ValueError: If hash verification fails
    """
    logger = get_logger()

    # Use a temporary file for atomic download
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    temp_fd, temp_path = tempfile.mkstemp(dir=dest_path.parent, suffix=".download")

    try:
        # Open URL and get content length
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "libephemeris-download/1.0"},
        )

        logger.info("Downloading %s...", description)

        _ssl_ctx = ssl.create_default_context(cafile=certifi.where())
        with urllib.request.urlopen(req, timeout=timeout, context=_ssl_ctx) as response:
            total_size = int(response.headers.get("Content-Length", 0))

            if total_size > 0:
                size_mb = total_size / (1024 * 1024)
                logger.info("File size: %.1f MB", size_mb)

            # Setup progress bar
            if show_progress and total_size > 0:
                progress = _get_progress_bar(total_size, description)
            else:
                progress = None

            # Download with progress updates
            sha256 = hashlib.sha256()
            chunk_size = 64 * 1024  # 64KB chunks

            with os.fdopen(temp_fd, "wb") as f:
                # Ownership of the fd has transferred to f, which closes it on
                # exit; mark it so the except clause below does not also close.
                temp_fd = -1
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    sha256.update(chunk)
                    if progress:
                        progress.update(len(chunk))

            if progress:
                progress.close()

        # Verify hash if provided
        if expected_sha256:
            actual_hash = sha256.hexdigest()
            if actual_hash != expected_sha256:
                os.unlink(temp_path)
                raise ValueError(
                    f"Hash mismatch: expected {expected_sha256}, got {actual_hash}"
                )

        # Atomic move to final destination
        publish_temp_file(temp_path, dest_path)
        logger.info("Download complete: %s", dest_path.name)
        return True

    # http.client.HTTPException (e.g. IncompleteRead when the server closes
    # the connection mid-body) is NOT an OSError, so it must be caught
    # explicitly to clean up the temp file instead of orphaning it.
    except (OSError, ValueError, KeyError, RuntimeError, http.client.HTTPException):
        # Close the temp fd if it was never handed to os.fdopen (e.g. urlopen
        # raised first); otherwise it leaks until the process exits.
        if temp_fd != -1:
            try:
                os.close(temp_fd)
            except OSError:
                pass
        # Clean up temp file on error
        if os.path.exists(temp_path):
            os.unlink(temp_path)
        raise


def download_planet_centers(
    force: bool = False,
    show_progress: bool = True,
    quiet: bool = False,
) -> Path:
    """Download the planet_centers.bsp data file.

    This file provides precise planet center positions for Jupiter, Saturn,
    Uranus, Neptune, and Pluto. It enables sub-arcsecond precision for
    outer planet calculations.

    Args:
        force: If True, download even if file already exists
        show_progress: If True, show progress bar during download
        quiet: If True, suppress all output except errors

    Returns:
        Path to the downloaded file

    Raises:
        urllib.error.URLError: If download fails
        ValueError: If hash verification fails
    """
    file_info = DATA_FILES["planet_centers.bsp"]
    dest_path = get_data_dir() / "planet_centers.bsp"
    logger = get_logger()

    # Check if already exists and valid
    if dest_path.exists() and not force:
        if _is_valid_bsp(str(dest_path)):
            if not quiet:
                print(f"planet_centers.bsp already exists at {dest_path}")
                print("Use --force to re-download.")
            return dest_path
        logger.warning("Cached file %s is corrupted, re-downloading", dest_path)
        try:
            os.remove(dest_path)
        except OSError:
            pass

    if not quiet:
        print(f"Downloading planet_centers.bsp (~{file_info['size_mb']:.1f} MB)...")
        print(f"  {file_info['description']}")
        print()

    try:
        download_file(
            url=str(file_info["url"]),
            dest_path=dest_path,
            description="planet_centers.bsp",
            expected_sha256=str(file_info["sha256"])
            if file_info.get("sha256")
            else None,
            show_progress=show_progress,
        )

        if not _is_valid_bsp(str(dest_path)):
            try:
                os.remove(dest_path)
            except OSError:
                pass
            raise ValueError(
                "Downloaded file planet_centers.bsp failed validation - "
                "corrupt or incomplete"
            )

        if not quiet:
            print()
            print(f"Downloaded to: {dest_path}")
            print(
                "Planet center data is now available for high-precision calculations."
            )

        return dest_path

    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(
                "\nError: Data file not found on server.",
                file=sys.stderr,
            )
            print(
                "The data release may not be available yet.",
                file=sys.stderr,
            )
            print(
                "libephemeris will use analytical approximations (~0.1 arcsec precision).",
                file=sys.stderr,
            )
        raise
    except (OSError, ValueError, KeyError, RuntimeError) as e:
        print(f"\nError downloading file: {e}", file=sys.stderr)
        raise


def download_all(
    force: bool = False,
    show_progress: bool = True,
    quiet: bool = False,
) -> list[Path]:
    """Download all optional data files.

    Args:
        force: If True, download even if files already exist
        show_progress: If True, show progress bar during download
        quiet: If True, suppress all output except errors

    Returns:
        List of paths to downloaded files
    """
    paths = []
    paths.append(
        download_planet_centers(force=force, show_progress=show_progress, quiet=quiet)
    )
    return paths


def check_data_status() -> dict[str, dict]:
    """Check the status of all data files.

    Returns:
        Dictionary with file names as keys and status info as values.
        Each status dict contains:
        - exists: bool
        - path: Path or None
        - size: int or None (file size in bytes)
        - description: str
    """
    data_dir = get_data_dir()
    status = {}

    for filename, info in DATA_FILES.items():
        # Resolve path based on file type.  All downloads land under the
        # data dir (download_planet_centers_for_tier and state.py agree
        # on data_dir/planet_centers_<tier>.bsp).
        if filename.endswith((".leb", ".leb2")):
            path = data_dir / "leb" / filename
        else:
            path = data_dir / filename

        exists = path.exists()
        status[filename] = {
            "exists": exists,
            "path": path if exists else None,
            "size": path.stat().st_size if exists else None,
            "description": info["description"],
            "expected_size_mb": info["size_mb"],
        }

    return status


def print_data_status(as_json: bool = False, verbose: int = 0) -> None:
    """Print comprehensive library and data file status.

    Shows version, calc mode, precision tier, LEB file, data directory,
    and status of all data files: DE kernels, planet centers, LEB1, LEB2,
    SPK cache, ASSIST, and IERS Earth orientation data.

    Args:
        as_json: If True, output machine-readable JSON instead of formatted text.
        verbose: Detail level for human-readable output.
    """
    import json as json_mod
    import os

    from .state import get_calc_mode, get_precision_tier

    data_dir = get_data_dir()
    workspace_root = Path(__file__).parent.parent
    current_tier = get_precision_tier()

    # Collect all status data into a dict (used for both JSON and text output)
    result: dict[str, Any] = {}

    # --- Library info ---
    try:
        from . import __version__ as ver
    except Exception:
        ver = "unknown"

    calc_mode = "unknown"
    try:
        calc_mode = get_calc_mode()
    except Exception:
        pass

    leb_path_str = "(none)"
    try:
        leb_path_val = os.environ.get("LIBEPHEMERIS_LEB")
        if not leb_path_val:
            from . import state as _state

            leb_path_val = getattr(_state, "_LEB_FILE", None)
        if leb_path_val:
            leb_path_str = str(leb_path_val)
    except Exception:
        pass

    result["version"] = str(ver)
    result["calc_mode"] = str(calc_mode)
    result["precision_tier"] = current_tier
    result["leb_file"] = leb_path_str
    result["data_directory"] = str(data_dir)

    # --- DE Ephemeris Kernels ---
    de_kernels_info: dict[str, dict[str, Any]] = {}
    de_files = {
        "base": ("de440s.bsp", "DE440s — lightweight (1849-2150)"),
        "medium": ("de440.bsp", "DE440 — general purpose (1549-2650)"),
        "extended": ("de441.bsp", "DE441 — full range (-13198 to +17191)"),
    }
    for tier_name, (filename, desc) in de_files.items():
        # DE kernels can be in the Skyfield data directory or data_dir
        path = Path(data_dir) / filename
        if not path.exists():
            # Check Skyfield's default location
            try:
                from skyfield.api import Loader

                loader = Loader(str(data_dir))
                sf_path = Path(loader.directory) / filename
                if sf_path.exists():
                    path = sf_path
            except Exception:
                pass
        exists = path.exists()
        size = path.stat().st_size if exists else None
        de_kernels_info[filename] = {
            "exists": exists,
            "size": size,
            "description": desc,
            "active": tier_name == current_tier,
            "path": str(path),
        }
    result["de_kernels"] = de_kernels_info

    # --- Planet Center Corrections ---
    planet_centers_info: dict[str, dict[str, Any]] = {}
    for tier_name in ["base", "medium", "extended"]:
        filename = f"planet_centers_{tier_name}.bsp"
        path = Path(data_dir) / filename
        if not path.exists():
            path = workspace_root / filename
        exists = path.exists()
        size = path.stat().st_size if exists else None
        planet_centers_info[filename] = {
            "exists": exists,
            "size": size,
            "active": tier_name == current_tier,
            "path": str(path),
        }
    result["planet_centers"] = planet_centers_info

    # --- LEB1 Binary Ephemeris ---
    leb1_info: dict[str, dict[str, Any]] = {}
    for tier_name in ["base", "medium", "extended"]:
        filename = f"ephemeris_{tier_name}.leb"
        path = Path(data_dir) / "leb" / filename
        exists = path.exists()
        size = path.stat().st_size if exists else None
        active = leb_path_str.endswith(filename) if leb_path_str != "(none)" else False
        leb1_info[filename] = {
            "exists": exists,
            "size": size,
            "active": active,
            "path": str(path),
        }
    result["leb1_files"] = leb1_info

    # --- LEB2 Compressed Files ---
    leb2_info: dict[str, dict[str, Any]] = {}
    leb2_groups = list(_CANONICAL_LEB2_GROUPS)
    for tier_name in ["base", "medium", "extended"]:
        for group in leb2_groups:
            filename = f"{tier_name}_{group}.leb2"
            path = Path(data_dir) / "leb" / filename
            exists = path.exists()
            size = path.stat().st_size if exists else None
            leb2_info[filename] = {
                "exists": exists,
                "size": size,
                "path": str(path),
            }
    result["leb2_files"] = leb2_info

    # --- SPK Asteroid Cache ---
    spk_info: dict[str, Any] = {}
    try:
        from .spk_auto import DEFAULT_AUTO_SPK_DIR
        from .state import get_spk_cache_dir

        spk_dir = get_spk_cache_dir() or DEFAULT_AUTO_SPK_DIR
        spk_path = Path(spk_dir)
        if spk_path.exists():
            spk_files = sorted(spk_path.glob("*.bsp"))
            total_size = sum(f.stat().st_size for f in spk_files)
            spk_info = {
                "directory": str(spk_path),
                "file_count": len(spk_files),
                "total_size": total_size,
                "files": [
                    {
                        "name": f.name,
                        "path": str(f),
                        "size": f.stat().st_size,
                    }
                    for f in spk_files
                ],
            }
        else:
            spk_info = {
                "directory": str(spk_path),
                "file_count": 0,
                "total_size": 0,
                "files": [],
            }
    except Exception:
        spk_info = {
            "directory": "(unknown)",
            "file_count": 0,
            "total_size": 0,
            "files": [],
        }
    result["spk_cache"] = spk_info

    # --- ASSIST N-body Data ---
    assist_info: dict[str, dict[str, Any]] = {}
    try:
        assist_dir = Path(os.path.expanduser("~/.libephemeris/assist"))
        for filename, desc in [
            ("linux_p1550p2650.440", "Planet ephemeris (~98 MB)"),
            ("sb441-n16.bsp", "Asteroid perturbers (~616 MB)"),
        ]:
            path = assist_dir / filename
            exists = path.exists()
            size = path.stat().st_size if exists else None
            assist_info[filename] = {
                "exists": exists,
                "size": size,
                "description": desc,
                "path": str(path),
            }
    except Exception:
        pass
    result["assist_files"] = assist_info

    # --- IERS Earth Orientation Data ---
    iers_info: dict[str, Any] = {}
    try:
        from .iers_data import get_iers_cache_info

        cache_info = get_iers_cache_info()
        cache_dir = Path(str(cache_info.get("cache_dir", "")))

        iers_info = {
            "cache_dir": str(cache_dir),
            "finals": {
                "exists": bool(cache_info.get("finals_exists")),
                "age_days": cache_info.get("finals_age_days"),
                "path": str(cache_dir / "finals2000A.data"),
            },
            "leap_seconds": {
                "exists": bool(cache_info.get("leap_seconds_exists")),
                "age_days": cache_info.get("leap_seconds_age_days"),
                "path": str(cache_dir / "leap_seconds.dat"),
            },
            "delta_t": {
                "exists": bool(cache_info.get("delta_t_exists")),
                "age_days": cache_info.get("delta_t_age_days"),
                "path": str(cache_dir / "deltat.data"),
            },
        }
    except Exception:
        pass
    result["iers_data"] = iers_info

    # --- Output ---
    import click as _click

    if as_json:
        # Convert Path objects and other non-serializable types
        def _serialize(obj: Any) -> Any:
            if isinstance(obj, Path):
                return str(obj)
            return obj

        _click.echo(json_mod.dumps(result, indent=2, default=_serialize))
        return

    # --- Formatted text output ---

    def _ok(text: str) -> str:
        return _click.style("[OK]", fg="green") + " " + text

    def _missing(text: str) -> str:
        return _click.style("[--]", fg="yellow") + " " + text

    def _active_marker(is_active: bool) -> str:
        return _click.style(" *", fg="cyan") if is_active else ""

    def _print_path(info: dict[str, Any]) -> None:
        if verbose < 1:
            return
        path = info.get("path")
        if path:
            _click.echo(f"    Path: {path}")

    _click.echo(_click.style(f"libephemeris {ver}", bold=True))
    _click.echo()

    # Library configuration
    _click.echo(_click.style("Configuration", bold=True))
    _click.echo(f"  Calc mode:       {calc_mode}")
    _click.echo(f"  Precision tier:  {current_tier}")
    _click.echo(f"  LEB file:        {leb_path_str}")
    _click.echo(f"  Data directory:  {data_dir}")
    _click.echo()

    # DE Ephemeris Kernels
    _click.echo(_click.style("Ephemeris Kernels", bold=True))
    for filename, info in de_kernels_info.items():
        marker = _active_marker(info["active"])
        if info["exists"]:
            _click.echo(f"  {_ok(filename)} ({_format_size(info['size'])}){marker}")
        else:
            _click.echo(f"  {_missing(filename)}{marker}")
        _print_path(info)
    _click.echo()

    # Planet Center Corrections
    _click.echo(_click.style("Planet Center Corrections", bold=True))
    for filename, info in planet_centers_info.items():
        marker = _active_marker(info["active"])
        if info["exists"]:
            _click.echo(f"  {_ok(filename)} ({_format_size(info['size'])}){marker}")
        else:
            _click.echo(f"  {_missing(filename)}{marker}")
        _print_path(info)
    _click.echo()

    # LEB1 Binary Ephemeris
    _click.echo(_click.style("LEB1 Binary Ephemeris (~14x speedup)", bold=True))
    for filename, info in leb1_info.items():
        marker = _active_marker(info["active"])
        if info["exists"]:
            _click.echo(f"  {_ok(filename)} ({_format_size(info['size'])}){marker}")
        else:
            _click.echo(f"  {_missing(filename)}{marker}")
        _print_path(info)
    _click.echo()

    # LEB2 Compressed Files
    _click.echo(_click.style("LEB2 Compressed Ephemeris (4-10x smaller)", bold=True))
    for tier_name in ["base", "medium", "extended"]:
        tier_files = [f"{tier_name}_{g}.leb2" for g in leb2_groups]
        present = sum(1 for f in tier_files if leb2_info.get(f, {}).get("exists"))
        total = len(tier_files)
        total_size = sum(
            leb2_info[f]["size"]
            for f in tier_files
            if leb2_info.get(f, {}).get("exists") and leb2_info[f]["size"]
        )
        if present == total:
            _click.echo(
                f"  {_ok(tier_name)}: {present}/{total} groups"
                f" ({_format_size(total_size)})"
            )
        elif present > 0:
            _click.echo(
                _click.style("[~~]", fg="yellow")
                + f" {tier_name}: {present}/{total} groups"
                f" ({_format_size(total_size)})"
            )
        else:
            _click.echo(f"  {_missing(tier_name)}: 0/{total} groups")
        if verbose >= 1:
            for filename in tier_files:
                info = leb2_info[filename]
                if info["exists"]:
                    _click.echo(f"    {_ok(filename)} ({_format_size(info['size'])})")
                else:
                    _click.echo(f"    {_missing(filename)}")
                _print_path(info)
    _click.echo()

    # SPK Asteroid Cache
    _click.echo(_click.style("SPK Asteroid Cache", bold=True))
    _click.echo(f"  Directory:  {spk_info.get('directory', '(unknown)')}")
    fc = spk_info.get("file_count", 0)
    ts = spk_info.get("total_size", 0)
    if fc > 0:
        _click.echo(f"  Files:      {fc} ({_format_size(ts)})")
    else:
        _click.echo("  Files:      (none)")
    if verbose >= 2:
        for info in spk_info.get("files", []):
            _click.echo(f"  {_ok(info['name'])} ({_format_size(info['size'])})")
            _click.echo(f"    Path: {info['path']}")
    _click.echo()

    # ASSIST N-body Data
    _click.echo(_click.style("ASSIST N-body Data", bold=True))
    if assist_info:
        for filename, info in assist_info.items():
            if info["exists"]:
                _click.echo(f"  {_ok(filename)} ({_format_size(info['size'])})")
            else:
                _click.echo(f"  {_missing(filename)}")
            _print_path(info)
    else:
        _click.echo("  (could not check)")
    _click.echo()

    # IERS Earth Orientation Data
    _click.echo(_click.style("IERS Earth Orientation Data", bold=True))
    if iers_info:
        for label, key in [
            ("finals2000A.data", "finals"),
            ("leap_seconds.dat", "leap_seconds"),
            ("deltat.data", "delta_t"),
        ]:
            finfo = iers_info.get(key, {})
            if finfo.get("exists"):
                age = finfo.get("age_days")
                age_str = f", {age:.0f} days old" if age is not None else ""
                _click.echo(f"  {_ok(label)}{age_str}")
            else:
                _click.echo(f"  {_missing(label)}")
            _print_path(finfo)
    else:
        _click.echo("  (could not check)")
    _click.echo()

    # Setup hints
    _click.echo(_click.style("Commands", bold=True))
    _click.echo("  libephemeris download <tier>         Download DE kernel + SPKs")
    _click.echo("  libephemeris download leb-<tier>     Download LEB1 binary ephemeris")
    _click.echo(
        "  libephemeris download leb2-<tier>    Download LEB2 compressed ephemeris"
    )
    _click.echo("  libephemeris download assist         Download ASSIST n-body data")
    _click.echo("  Available tiers: base, medium, extended")


def init_all(
    cache_dir: Optional[str] = None,
    force: bool = False,
    show_progress: bool = True,
    quiet: bool = False,
    start_year: int = 1550,
    end_year: int = 2650,
) -> dict:
    """Initialize libephemeris with all required data files.

    Downloads:
    1. DE440.bsp planetary ephemeris (~128 MB)
    2. planet_centers.bsp precision data (~25 MB)
    3. SPK kernels for all minor bodies defined in SPK_BODY_NAME_MAP
       (default 1550-2650, 20-year chunks)

    Args:
        cache_dir: Custom SPK cache directory (default: ~/.libephemeris/spk/).
            Also configurable via LIBEPHEMERIS_SPK_DIR env var or
            set_spk_cache_dir().
        force: Re-download existing files
        show_progress: Show progress output
        quiet: Suppress non-error output
        start_year: First year of SPK coverage (default: 1550)
        end_year: Last year of SPK coverage (default: 2650)

    Returns:
        dict with summary:
            - de440: bool (success)
            - planet_centers: bool (success)
            - spk_success: int (number of SPK chunks downloaded)
            - spk_failed: int (number of SPK chunks that failed)
            - spk_skipped: int (number of SPK chunks already cached)
            - spk_out_of_range: int (chunks outside body's Horizons SPK limits)
            - failed_details: list of (body_name, chunk_start, chunk_end, error)
    """
    import time

    logger = get_logger()

    result: dict[str, Any] = {
        "de440": False,
        "planet_centers": False,
        "spk_success": 0,
        "spk_failed": 0,
        "spk_skipped": 0,
        "spk_out_of_range": 0,
        "failed_details": [],
    }

    # -------------------------------------------------------------------------
    # Step 1: DE440.bsp
    # -------------------------------------------------------------------------
    if not quiet:
        print("Step 1/3: DE440.bsp planetary ephemeris")

    try:
        from .state import get_planets

        get_planets()  # Downloads if not present
        result["de440"] = True
        if not quiet:
            print("  [OK] DE440.bsp ready")
    except (OSError, ValueError, KeyError, RuntimeError) as e:
        logger.error("Failed to initialize DE440: %s", e)
        if not quiet:
            print(f"  [FAIL] DE440.bsp: {e}", file=sys.stderr)

    # -------------------------------------------------------------------------
    # Step 2: planet_centers.bsp
    # -------------------------------------------------------------------------
    if not quiet:
        print("Step 2/3: planet_centers.bsp precision data")

    try:
        download_planet_centers(
            force=force,
            show_progress=show_progress,
            quiet=quiet,
        )
        result["planet_centers"] = True
        if not quiet:
            print("  [OK] planet_centers.bsp ready")
    except (OSError, ValueError, KeyError, RuntimeError) as e:
        logger.error("Failed to download planet_centers.bsp: %s", e)
        if not quiet:
            print(f"  [FAIL] planet_centers.bsp: {e}", file=sys.stderr)

    # -------------------------------------------------------------------------
    # Step 3: SPK for all minor bodies
    # -------------------------------------------------------------------------
    if not quiet:
        print("Step 3/3: SPK kernels for minor bodies")

    from .constants import SPK_BODY_NAME_MAP
    from .spk_auto import (
        download_spk_from_horizons,
        _generate_spk_cache_filename,
        _iso_to_jd,
        ensure_cache_dir,
    )
    from .spk import _get_body_name

    effective_cache_dir = ensure_cache_dir(cache_dir)

    if not quiet:
        print(f"  Cache directory: {effective_cache_dir}")

    CHUNK_SIZE_YEARS = 20
    START_YEAR = start_year
    END_YEAR = end_year
    chunks_per_body = (END_YEAR - START_YEAR) // CHUNK_SIZE_YEARS
    total_bodies = len(SPK_BODY_NAME_MAP)

    if not quiet:
        print(
            f"  Bodies: {total_bodies} | "
            f"Chunks per body: {chunks_per_body} | "
            f"Total: {chunks_per_body * total_bodies}"
        )
        print()

    from .constants import SPK_AUTO_DOWNLOAD_BLOCKED

    for body_idx, (ipl, (horizons_id, naif_id)) in enumerate(SPK_BODY_NAME_MAP.items()):
        body_name = _get_body_name(ipl) or f"body_{ipl}"

        # Skip bodies where JPL blocks SPK generation
        if ipl in SPK_AUTO_DOWNLOAD_BLOCKED:
            if not quiet:
                print(
                    f"  [{body_idx + 1}/{total_bodies}] {body_name} "
                    f"— skipped ({SPK_AUTO_DOWNLOAD_BLOCKED[ipl].split('.')[0]})"
                )
            result.setdefault("spk_blocked", []).append(body_name)
            continue

        if not quiet:
            print(
                f"  [{body_idx + 1}/{total_bodies}] {body_name} "
                f"(Horizons: {horizons_id})"
            )

        body_success = 0
        body_skipped = 0
        body_failed = 0
        body_out_of_range = 0
        chunk_idx = 0

        for chunk_start in range(START_YEAR, END_YEAR, CHUNK_SIZE_YEARS):
            chunk_end = min(chunk_start + CHUNK_SIZE_YEARS, END_YEAR)
            chunk_idx += 1

            jd_start = _iso_to_jd(f"{chunk_start}-01-01")
            jd_end = _iso_to_jd(f"{chunk_end}-01-01")

            filename = _generate_spk_cache_filename(horizons_id, jd_start, jd_end)
            output_path = os.path.join(effective_cache_dir, filename)

            # Show per-chunk progress
            if not quiet:
                sys.stdout.write(
                    f"\r    chunk {chunk_idx}/{chunks_per_body} "
                    f"({chunk_start}-{chunk_end}) ... "
                )
                sys.stdout.flush()

            if os.path.exists(output_path) and not force:
                body_skipped += 1
                result["spk_skipped"] += 1
                continue

            try:
                # naif_id deliberately omitted: Horizons-generated kernels
                # may use 20000000+N target IDs while SPK_BODY_NAME_MAP holds
                # the classic 2000000+N; registration auto-detects the ID
                # from the downloaded file.
                download_spk_from_horizons(
                    body_id=horizons_id,
                    jd_start=jd_start,
                    jd_end=jd_end,
                    output_path=output_path,
                    ipl=ipl,
                )
                body_success += 1
                result["spk_success"] += 1

                # Rate limiting: small delay between Horizons requests
                time.sleep(1.5)

            except KeyboardInterrupt:
                if not quiet:
                    sys.stdout.write("\n")
                raise
            except (OSError, ValueError, KeyError, RuntimeError) as e:
                error_msg = str(e)

                # Detect Horizons SPK limit errors to skip gracefully
                if "START time outside set SPK limits" in error_msg:
                    # Body data doesn't go back this far; try later chunks
                    body_out_of_range += 1
                    result["spk_out_of_range"] += 1
                    logger.debug(
                        "%s: chunk %d-%d before SPK start, skipping",
                        body_name,
                        chunk_start,
                        chunk_end,
                    )
                    continue

                if "STOP time outside set SPK limits" in error_msg:
                    # Body data ends before this chunk; skip remaining chunks
                    remaining = (END_YEAR - chunk_start) // CHUNK_SIZE_YEARS
                    body_out_of_range += remaining
                    result["spk_out_of_range"] += remaining
                    logger.debug(
                        "%s: chunk %d-%d past SPK end, skipping %d remaining",
                        body_name,
                        chunk_start,
                        chunk_end,
                        remaining,
                    )
                    break

                # Genuine failure
                body_failed += 1
                result["spk_failed"] += 1
                result["failed_details"].append(
                    (body_name, chunk_start, chunk_end, error_msg)
                )
                logger.warning(
                    "Failed to download %s chunk %d-%d: %s",
                    body_name,
                    chunk_start,
                    chunk_end,
                    e,
                )

        if not quiet:
            # Clear the per-chunk line and print summary
            sys.stdout.write("\r" + " " * 60 + "\r")
            parts = []
            if body_success:
                parts.append(f"{body_success} downloaded")
            if body_skipped:
                parts.append(f"{body_skipped} cached")
            if body_out_of_range:
                parts.append(f"{body_out_of_range} out of range")
            if body_failed:
                parts.append(f"{body_failed} failed")
            print(f"    {', '.join(parts)}")

    return result


def download_for_tier(
    tier_name: str,
    force: bool = False,
    show_progress: bool = True,
    quiet: bool = False,
) -> dict:
    """Download all data files required for a specific precision tier.

    Sets the precision tier, then downloads:
    1. The tier's ephemeris file (de440s / de440 / de441)
    2. planet_centers.bsp precision offsets
    3. SPK kernels for all bodies in ``SPK_BODY_NAME_MAP`` (full tier range)

    For the 'extended' tier, SPK files are downloaded as single max-range
    files (1600-2500) rather than chunked, since that's the full Horizons
    range for minor bodies.

    Args:
        tier_name: One of "base", "medium", "extended"
        force: If True, re-download even if files already exist
        show_progress: If True, show progress output
        quiet: If True, suppress non-error output

    Returns:
        Dict with summary (same format as ensure_all_ephemerides)

    Raises:
        ValueError: If tier_name is not valid
    """
    from .state import TIERS, set_precision_tier

    if tier_name not in TIERS:
        valid = ", ".join(sorted(TIERS.keys()))
        raise ValueError(f"Unknown tier '{tier_name}'. Valid tiers: {valid}")

    tier = TIERS[tier_name]

    if not quiet:
        print(f"libephemeris: downloading data for tier '{tier_name}'")
        print(f"  Ephemeris file: {tier.ephemeris_file}")
        print(f"  Description:    {tier.description}")
        print(f"  SPK date range: {tier.spk_date_range[0]} to {tier.spk_date_range[1]}")
        print()

    # Activate the tier so ensure_all_ephemerides uses the right config
    set_precision_tier(tier_name)

    # Step 1: planet_centers file for this tier (workspace root)
    if not quiet:
        print("Step 1/2: planet_centers precision data")

    try:
        _download_planet_centers_for_tier(
            tier_name=tier_name,
            force=force,
            show_progress=show_progress,
            quiet=quiet,
        )
        if not quiet:
            print("  [OK] planet_centers ready")
            print()
    except (OSError, ValueError, KeyError, RuntimeError) as e:
        if not quiet:
            print(f"  [WARN] planet_centers: {e}", file=sys.stderr)
            print("  (non-critical, continuing...)")
            print()

    # Step 2: Ephemeris file + SPK kernels (via ensure_all_ephemerides)
    if not quiet:
        print(f"Step 2/2: {tier.ephemeris_file} + SPK kernels for minor bodies")

    from .spk_auto import ensure_all_ephemerides

    results = ensure_all_ephemerides(
        force_download=force,
        show_progress=show_progress and not quiet,
    )

    # Print final summary
    if not quiet:
        summary = results.get("summary", {})
        print()
        print("Download complete:")
        print(f"  Tier:         {tier_name}")
        print(f"  Ephemeris:    {tier.ephemeris_file}")
        if isinstance(summary, dict):
            print(f"  SPK cached:   {summary.get('cached', 0)}")
            print(f"  SPK downloaded: {summary.get('downloaded', 0)}")
            if summary.get("fallback", 0):
                print(f"  SPK fallback: {summary['fallback']}")
            if summary.get("errors", 0):
                print(f"  SPK errors:   {summary['errors']}")

    return results


def _is_valid_leb(filepath: str) -> bool:
    """Check if a LEB file can be opened and parsed by LEBReader.

    Validates magic bytes, version, section directory, and body index.

    Args:
        filepath: Path to the LEB file to validate

    Returns:
        True if file can be opened and parsed, False otherwise
    """
    try:
        from .leb_reader import open_leb

        reader = open_leb(filepath)
        reader.close()
        return True
    except (OSError, ValueError, KeyError, RuntimeError):
        return False


def get_leb_dir() -> Path:
    """Get the LEB data directory path.

    Returns the path to ~/.libephemeris/leb/ (or under LIBEPHEMERIS_DATA_DIR).
    Creates the directory if it doesn't exist.

    Returns:
        Path to the LEB data directory
    """
    leb_dir = get_data_dir() / "leb"
    leb_dir.mkdir(parents=True, exist_ok=True)
    return leb_dir


def get_leb_path_for_tier(tier_name: str) -> Path:
    """Get the expected LEB file path for a tier.

    Args:
        tier_name: One of "base", "medium", "extended"

    Returns:
        Path where the LEB file should be stored
    """
    return get_leb_dir() / f"ephemeris_{tier_name}.leb"


def download_leb_for_tier(
    tier_name: str,
    force: bool = False,
    show_progress: bool = True,
    quiet: bool = False,
    activate: bool = True,
) -> Path:
    """Download the precomputed LEB binary ephemeris file for a specific tier.

    Downloads the .leb file from GitHub Releases to ~/.libephemeris/leb/.
    After a successful download, optionally activates it via set_leb_file().

    LEB files contain precomputed Chebyshev polynomial approximations for
    all celestial bodies, providing ~14x speedup over the Skyfield pipeline.

    Args:
        tier_name: One of "base", "medium", "extended"
        force: If True, re-download even if file already exists
        show_progress: If True, show progress bar during download
        quiet: If True, suppress all output except errors
        activate: If True, call set_leb_file() after successful download

    Returns:
        Path to the downloaded LEB file

    Raises:
        ValueError: If tier_name is invalid or hash verification fails
        urllib.error.URLError: If download fails
        RuntimeError: If the LEB file has no download hash configured
    """
    logger = get_logger()
    filename = f"ephemeris_{tier_name}.leb"
    file_info = DATA_FILES.get(filename)

    if file_info is None:
        valid_tiers = ["base", "medium", "extended"]
        raise ValueError(
            f"Unknown tier '{tier_name}'. Valid tiers: {', '.join(valid_tiers)}"
        )

    # Guard: file must have a download hash configured
    if file_info.get("sha256") is None:
        raise RuntimeError(
            f"LEB file for '{tier_name}' tier has no download hash configured. "
            f"You can generate it locally with: leph leb generate {tier_name}"
        )

    dest_path = get_leb_path_for_tier(tier_name)

    # Check if already exists and valid
    if dest_path.exists() and not force:
        if _is_valid_leb(str(dest_path)):
            if not quiet:
                print(f"  {filename} already exists at {dest_path}")
                print("  Use --force to re-download.")
            # Activate if requested
            if activate:
                from .state import set_leb_file

                set_leb_file(str(dest_path))
            return dest_path
        else:
            logger.warning("Cached LEB file %s is corrupted, re-downloading", dest_path)
            try:
                os.remove(dest_path)
            except OSError:
                pass

    size_mb = file_info.get("size_mb", 0)
    if not quiet:
        print(f"  Downloading {filename} (~{size_mb:.0f} MB)...")
        print(f"  {file_info['description']}")
        print()

    # download_file streams to a temp file and only replaces dest_path
    # after the hash verifies, so a failure here leaves any existing
    # (force=True) valid file untouched — do not delete dest_path.
    download_file(
        url=str(file_info["url"]),
        dest_path=dest_path,
        description=filename,
        expected_sha256=str(file_info["sha256"]),
        show_progress=show_progress,
    )

    # Validate the downloaded file
    if not _is_valid_leb(str(dest_path)):
        try:
            os.remove(dest_path)
        except OSError:
            pass
        raise ValueError(
            f"Downloaded {filename} failed LEB validation — corrupt or incompatible"
        )

    if not quiet:
        print()
        print(f"  Downloaded to: {dest_path}")
        print(f"  LEB binary ephemeris for '{tier_name}' tier is now available.")
        if activate:
            print("  Activated automatically for this session.")

    # Activate if requested
    if activate:
        from .state import set_leb_file

        set_leb_file(str(dest_path))
        logger.info("Activated LEB file: %s", dest_path)

    return dest_path


# Canonical distribution groups (incl. exotics). A group with no DATA_FILES
# manifest entry yet (e.g. exotics until its artifact is published) is skipped
# gracefully with a visible message rather than silently omitted.
LEB2_GROUPS = list(_CANONICAL_LEB2_GROUPS)


def download_leb2_for_tier(
    tier_name: str,
    groups: Optional[list] = None,
    force: bool = False,
    show_progress: bool = True,
    quiet: bool = False,
    activate: bool = True,
) -> list:
    """Download LEB2 compressed modular ephemeris files for a tier.

    Downloads group files (core, asteroids, exotics, apogee, uranians)
    from GitHub Releases to ~/.libephemeris/leb/.

    Args:
        tier_name: One of "base", "medium", "extended"
        groups: List of groups to download. Default: all groups.
        force: Re-download even if files exist.
        show_progress: Show progress bars.
        quiet: Suppress output.
        activate: Activate the core file via set_leb_file() after download.

    Returns:
        List of downloaded file paths.
    """
    logger = get_logger()
    if groups is None:
        groups = list(LEB2_GROUPS)

    leb_dir = get_leb_dir()
    downloaded = []

    for group in groups:
        filename = f"{tier_name}_{group}.leb2"
        file_info = DATA_FILES.get(filename)
        if file_info is None:
            if not quiet:
                print(f"  [SKIP] {filename}: not in DATA_FILES")
            continue

        dest_path = leb_dir / filename

        if dest_path.exists() and not force:
            if _is_valid_leb(str(dest_path)):
                if not quiet:
                    print(f"  [OK] {filename} already exists")
                downloaded.append(dest_path)
                continue

        size_mb = file_info.get("size_mb", 0)
        if not quiet:
            print(f"  Downloading {filename} (~{size_mb:.0f} MB)...")

        try:
            download_file(
                url=str(file_info["url"]),
                dest_path=dest_path,
                description=filename,
                expected_sha256=file_info.get("sha256"),
                show_progress=show_progress and not quiet,
            )
            # Validate the downloaded file (download_file only checks the hash
            # when one is present; a DATA_FILES entry without a sha256 would
            # otherwise let a truncated/corrupt file through). Mirror
            # download_leb_for_tier / _download_planet_centers_for_tier.
            if not _is_valid_leb(str(dest_path)):
                try:
                    dest_path.unlink()
                except OSError:
                    pass
                if not quiet:
                    print(f"  [FAIL] {filename}: downloaded file is corrupt")
                continue
            downloaded.append(dest_path)
            if not quiet:
                print(f"  [OK] {filename}")
        except (OSError, ValueError, KeyError, RuntimeError) as e:
            if not quiet:
                print(f"  [FAIL] {filename}: {e}")

    # Activate core file
    if activate and any(str(p).endswith("_core.leb2") for p in downloaded):
        core_path = leb_dir / f"{tier_name}_core.leb2"
        if core_path.exists():
            from .state import set_leb_file

            set_leb_file(str(core_path))
            logger.info("Activated LEB2 file: %s", core_path)

    return downloaded


def _download_planet_centers_for_tier(
    tier_name: str,
    force: bool = False,
    show_progress: bool = True,
    quiet: bool = False,
) -> Path:
    """Download planet_centers SPK file for a specific tier.

    Saves to data directory (~/.libephemeris).

    Args:
        tier_name: One of "base", "medium", "extended"
        force: If True, re-download even if file exists
        show_progress: If True, show progress bar
        quiet: If True, suppress output

    Returns:
        Path to the downloaded file

    Raises:
        Exception: If download fails
    """

    logger = get_logger()
    filename = f"planet_centers_{tier_name}.bsp"
    file_info = DATA_FILES.get(filename)

    if file_info is None:
        raise ValueError(f"No planet_centers file for tier '{tier_name}'")

    dest_path = get_data_dir() / filename

    if dest_path.exists() and not force:
        if _is_valid_bsp(str(dest_path)):
            if not quiet:
                print(f"  {filename} already exists at {dest_path}")
            return dest_path
        else:
            logger.warning("Cached file %s is corrupted, re-downloading", dest_path)
            try:
                os.remove(dest_path)
            except OSError:
                pass

    url = file_info["url"]
    expected_size_mb = file_info.get("size_mb", 50)

    if not quiet:
        print(f"  Downloading {filename} (~{expected_size_mb:.0f} MB)...")

    # download_file streams to a temp file, verifies the release sha256,
    # and atomically replaces dest_path — an existing (force=True) valid
    # file survives any failure.
    download_file(
        url=str(url),
        dest_path=dest_path,
        description=filename,
        expected_sha256=str(file_info["sha256"]) if file_info.get("sha256") else None,
        show_progress=show_progress and not quiet,
    )

    if not _is_valid_bsp(str(dest_path)):
        try:
            os.remove(dest_path)
        except OSError:
            pass
        raise ValueError(
            f"Downloaded file {filename} failed validation - corrupt or incomplete"
        )

    if not quiet:
        print(f"  Downloaded to {dest_path}")

    return dest_path
