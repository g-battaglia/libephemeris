# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Automatic SPK download and caching for minor bodies.

Provenance:
    Project-authored cache, locking, coverage, and integrity orchestration over
    the public NASA/JPL Horizons downloader in ``spk.py``. Body identifiers and
    SPK bytes retain their JPL/NAIF provenance; cache policy cannot change a
    state vector and produces no astronomical coefficient.

This module provides automatic SPK file management for minor body calculations:
- Automatic download of SPK kernels from JPL Horizons on demand
- Local caching to avoid redundant downloads
- Seamless integration with calc_ut() for transparent SPK usage

Downloads use :func:`libephemeris.spk.download_spk`, the project's direct JPL
Horizons HTTP client.  No optional download library is required.

Usage:
    >>> from libephemeris import spk_auto
    >>> # Enable auto-downloading for a body
    >>> spk_auto.enable_auto_spk(
    ...     ipl=CHIRON,
    ...     body_id="2060",  # Horizons identifier
    ...     start="2000-01-01",
    ...     end="2100-01-01",
    ... )
    >>> # Now calc_ut automatically downloads and uses SPK if needed
    >>> pos, _ = calc_ut(2451545.0, CHIRON, 0)

References:
    - JPL Horizons: https://ssd.jpl.nasa.gov/horizons/
"""

from __future__ import annotations

import hashlib
import os
import threading
from typing import Dict, Optional, Union

from .logging_config import get_logger
from .download import _is_valid_bsp
from .exceptions import Error


def _detect_naif_id_from_file(
    spk_path: str, preferred: Optional[int] = None
) -> Optional[int]:
    """Read the target NAIF ID from a downloaded SPK file.

    Horizons-generated small-body kernels may number their targets
    20000000+N while the classic convention is 2000000+N; the file itself
    is the source of truth.  Returns ``preferred`` when it is present in
    the file, the single small-body target when the file has exactly one,
    otherwise None.
    """
    from . import spk

    try:
        targets = spk._get_spk_targets(spk_path)
    except (OSError, ValueError, KeyError, ImportError):
        return None
    if not targets:
        return None
    if preferred is not None and preferred in targets:
        return preferred
    small_bodies = [t for t in targets if t > 1000000]
    if len(small_bodies) == 1:
        return small_bodies[0]
    return None


# =============================================================================
# MODULE STATE
# =============================================================================

# Lock for thread-safe operations
_AUTO_SPK_LOCK = threading.RLock()

# Registry of auto-SPK configurations: {ipl: AutoSpkConfig}
_AUTO_SPK_REGISTRY: dict[int, "AutoSpkConfig"] = {}

# Default cache directory for SPK files
DEFAULT_AUTO_SPK_DIR = os.path.join(os.path.expanduser("~"), ".libephemeris", "spk")

# Default date range for auto-downloaded SPK files
DEFAULT_START_DATE = "2000-01-01"
DEFAULT_END_DATE = "2100-01-01"


# =============================================================================
# CACHE DIRECTORY MANAGEMENT
# =============================================================================


def ensure_cache_dir(cache_dir: Optional[str] = None) -> str:
    """
    Ensure the SPK cache directory exists, creating it if necessary.

    This function creates the cache directory if it doesn't exist and returns
    the absolute path to the cache directory.

    Resolution order when cache_dir is None:
        1. Value from set_spk_cache_dir() (programmatic override)
        2. LIBEPHEMERIS_SPK_DIR environment variable
        3. DEFAULT_AUTO_SPK_DIR (~/.libephemeris/spk/)

    Args:
        cache_dir: Optional custom cache directory path. If None, resolves
            via set_spk_cache_dir() / env var / default.

    Returns:
        str: Absolute path to the cache directory.

    Example:
        >>> from libephemeris.spk_auto import ensure_cache_dir
        >>> cache_path = ensure_cache_dir()
        >>> print(cache_path)
        /home/user/.libephemeris/spk
        >>> # Custom cache directory
        >>> cache_path = ensure_cache_dir("/custom/cache/path")
    """
    if cache_dir is None:
        from .state import get_spk_cache_dir

        cache_dir = get_spk_cache_dir()
    if cache_dir is None:
        cache_dir = DEFAULT_AUTO_SPK_DIR

    cache_path = os.path.abspath(cache_dir)
    os.makedirs(cache_path, exist_ok=True)
    return cache_path


def get_cache_path(body_id: Union[int, str], cache_dir: Optional[str] = None) -> str:
    """
    Get the path where an SPK file for a body should be stored.

    This function returns the expected cache file path for a given body ID.
    The path is constructed using the body_id as part of the filename.
    Note: This does not create the file or directory - use ensure_cache_dir()
    to create the cache directory if needed.

    Args:
        body_id: JPL Horizons body identifier. Can be:
            - Asteroid number (int or str): 2060, "2060", "136199"
            - Name (str): "Chiron", "Eris"
        cache_dir: Optional custom cache directory path. If None, uses the
            default cache directory ({library_path}/spk_cache/).

    Returns:
        str: Full path where the SPK file would be stored.
            The filename format is: {sanitized_body_id}.bsp

    Example:
        >>> from libephemeris.spk_auto import get_cache_path
        >>> path = get_cache_path("2060")
        >>> print(path)
        /path/to/libephemeris/spk_cache/2060.bsp
        >>> path = get_cache_path("Chiron")
        >>> print(path)
        /path/to/libephemeris/spk_cache/chiron.bsp
    """
    if cache_dir is not None:
        dir_path = os.path.abspath(cache_dir)
    else:
        dir_path = DEFAULT_AUTO_SPK_DIR

    # Sanitize body_id for filename (shared writer/reader sanitizer)
    from .spk import _sanitize_filename

    safe_body_id = _sanitize_filename(body_id)

    filename = f"{safe_body_id}.bsp"
    return os.path.join(dir_path, filename)


def cache_info(
    cache_dir: Optional[str] = None,
) -> dict[str, Union[int, float, str, list[str]]]:
    """
    Get information about the SPK cache directory.

    This function returns statistics about cached SPK files including the
    number of files, total size, and list of cached filenames.

    Args:
        cache_dir: Optional custom cache directory path. If None, uses the
            default cache directory ({library_path}/spk_cache/).

    Returns:
        dict: Cache information with the following keys:
            - cache_dir (str): Absolute path to the cache directory
            - num_files (int): Number of cached SPK files (.bsp files)
            - total_size_mb (float): Total size of all cached files in megabytes
            - files (list[str]): List of cached SPK filenames

    Example:
        >>> from libephemeris.spk_auto import cache_info
        >>> info = cache_info()
        >>> print(f"Cache at: {info['cache_dir']}")
        >>> print(f"Files cached: {info['num_files']}")
        >>> print(f"Total size: {info['total_size_mb']:.2f} MB")
    """
    if cache_dir is not None:
        dir_path = os.path.abspath(cache_dir)
    else:
        dir_path = DEFAULT_AUTO_SPK_DIR

    if not os.path.exists(dir_path):
        return {
            "cache_dir": dir_path,
            "num_files": 0,
            "total_size_mb": 0.0,
            "files": [],
        }

    files = [f for f in os.listdir(dir_path) if f.endswith(".bsp")]
    total_size = sum(os.path.getsize(os.path.join(dir_path, f)) for f in files)

    return {
        "cache_dir": dir_path,
        "num_files": len(files),
        "total_size_mb": round(total_size / (1024 * 1024), 2),
        "files": files,
    }


# =============================================================================
# CONFIGURATION CLASS
# =============================================================================


class AutoSpkConfig:
    """Configuration for automatic SPK download for a minor body.

    Attributes:
        ipl: libephemeris body ID (e.g., CHIRON)
        body_id: JPL Horizons identifier (e.g., "2060", "Chiron")
        naif_id: NAIF ID for SPK lookup (e.g., 2002060)
        start: Start date for SPK coverage (YYYY-MM-DD)
        end: End date for SPK coverage (YYYY-MM-DD)
        cache_dir: Directory for storing downloaded SPK files
        enabled: Whether auto-download is enabled
        spk_path: Path to downloaded SPK file (None if not yet downloaded)
    """

    def __init__(
        self,
        ipl: int,
        body_id: str,
        naif_id: Optional[int] = None,
        start: str = DEFAULT_START_DATE,
        end: str = DEFAULT_END_DATE,
        cache_dir: Optional[str] = None,
    ):
        self.ipl = ipl
        self.body_id = body_id
        self.naif_id = naif_id
        self.start = start
        self.end = end
        self.cache_dir = cache_dir
        self.enabled = True
        self.spk_path: Optional[str] = None

    def get_cache_dir(self) -> str:
        """Get the cache directory, creating if necessary."""
        from . import state

        cache_path: str
        if self.cache_dir:
            cache_path = self.cache_dir
        elif state.get_spk_cache_dir() is not None:
            cache_path = state.get_spk_cache_dir()  # type: ignore
        else:
            cache_path = DEFAULT_AUTO_SPK_DIR

        if not os.path.exists(cache_path):
            os.makedirs(cache_path, exist_ok=True)

        return cache_path

    def get_cache_filename(self) -> str:
        """Generate a unique cache filename based on configuration."""
        # Create a hash of the configuration for unique identification
        config_str = f"{self.body_id}_{self.start}_{self.end}"
        config_hash = hashlib.md5(config_str.encode()).hexdigest()[:8]

        # Sanitize body_id for filename
        safe_body_id = "".join(
            c if c.isalnum() or c in "-_" else "_" for c in self.body_id.lower()
        )

        return f"{safe_body_id}_{config_hash}.bsp"

    def get_cache_path(self) -> str:
        """Get the full path where the SPK file should be cached."""
        return os.path.join(self.get_cache_dir(), self.get_cache_filename())


# =============================================================================
# DIRECT HORIZONS DOWNLOAD
# =============================================================================


def _download_spk_astroquery(
    body_id: str,
    start: str,
    end: str,
    output_path: str,
    location: str = "@0",
) -> str:
    """Download an SPK file using the direct Horizons HTTP client.

    The private name is retained for compatibility with existing callers and
    tests, but the implementation no longer uses Astroquery.  It wraps
    :func:`libephemeris.spk.download_spk` and normalizes the returned path to
    the cache filename requested by this module.

    Args:
        body_id: JPL Horizons target identifier (e.g., "2060", "Chiron")
        start: Start date (YYYY-MM-DD)
        end: End date (YYYY-MM-DD)
        output_path: Path to save the SPK file
        location: Observer location code (default: "@0" = solar system barycenter)

    Returns:
        str: Path to the downloaded SPK file

    Raises:
        ValueError: If body not found or download fails
    """
    from . import spk

    logger = get_logger()

    logger.info(
        "Downloading SPK for %s (%s to %s) via JPL Horizons API...", body_id, start, end
    )

    # Get the output directory from the path
    output_dir = os.path.dirname(output_path)
    if not output_dir:
        output_dir = "."

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    try:
        # Use the main download_spk function which handles the API correctly
        result_path = spk.download_spk(
            body=body_id,
            start=start,
            end=end,
            path=output_dir,
            center="500@0",
            overwrite=True,
        )

        # download_spk names the file with its own convention; callers rely
        # on the file living at exactly output_path (cache lookups,
        # registration), so move it there.
        if os.path.abspath(result_path) != os.path.abspath(output_path):
            os.replace(result_path, output_path)

        logger.info("SPK download complete: %s", os.path.basename(output_path))
        return output_path
    except (OSError, ValueError, KeyError, ImportError) as e:
        raise ValueError(
            f"Failed to download SPK for '{body_id}' from JPL Horizons: {e}"
        ) from e


def _check_astroquery_available() -> bool:
    """Return whether Astroquery is installed (legacy compatibility probe).

    Auto-SPK downloads no longer consult this helper or require Astroquery.
    Catalog-building and validation scripts may still use the dependency.
    """
    try:
        from astroquery.jplhorizons import Horizons  # noqa: F401 (availability probe)

        return True
    except ImportError:
        return False


# =============================================================================
# PUBLIC API
# =============================================================================


def enable_auto_spk(
    ipl: int,
    body_id: str,
    naif_id: Optional[int] = None,
    start: str = DEFAULT_START_DATE,
    end: str = DEFAULT_END_DATE,
    cache_dir: Optional[str] = None,
) -> None:
    """
    Enable automatic SPK download for a minor body.

    When enabled, the library will automatically download and cache an SPK
    kernel from JPL Horizons when the body's position is calculated.

    Args:
        ipl: libephemeris body ID (e.g., CHIRON, ERIS)
        body_id: JPL Horizons target identifier. Can be:
            - Asteroid number: "2060", "136199"
            - Name: "Chiron", "Eris"
            - Designation: "2003 UB313"
        naif_id: NAIF ID for SPK lookup. If None, will be deduced from body_id
            using convention: asteroid_number + 2000000
        start: Start date for SPK coverage (YYYY-MM-DD). Default: "2000-01-01"
        end: End date for SPK coverage (YYYY-MM-DD). Default: "2100-01-01"
        cache_dir: Directory for storing downloaded SPK files.
            Default: {library_path}/spk_cache/

    Example:
        >>> from libephemeris import spk_auto, CHIRON
        >>> spk_auto.enable_auto_spk(
        ...     ipl=CHIRON,
        ...     body_id="2060",
        ...     start="2000-01-01",
        ...     end="2100-01-01",
        ... )
        >>> # Future calculations for Chiron will use SPK automatically

    Note:
        This function only registers the configuration. The actual download
        occurs lazily when the body's position is first calculated, or
        immediately if download_now() is called.
    """
    # Import state to get configuration settings
    from . import state

    # Use parameter if provided, otherwise fall back to global config
    if cache_dir is None:
        cache_dir = state.get_spk_cache_dir()

    config = AutoSpkConfig(
        ipl=ipl,
        body_id=body_id,
        naif_id=naif_id,
        start=start,
        end=end,
        cache_dir=cache_dir,
    )

    with _AUTO_SPK_LOCK:
        _AUTO_SPK_REGISTRY[ipl] = config


def disable_auto_spk(ipl: int) -> None:
    """
    Disable automatic SPK download for a body.

    This removes the auto-SPK configuration but does not delete any
    already-downloaded SPK files or unregister them from the SPK registry.

    Args:
        ipl: libephemeris body ID (e.g., CHIRON)
    """
    with _AUTO_SPK_LOCK:
        if ipl in _AUTO_SPK_REGISTRY:
            del _AUTO_SPK_REGISTRY[ipl]


def is_auto_spk_enabled(ipl: int) -> bool:
    """
    Check if automatic SPK download is enabled for a body.

    Args:
        ipl: libephemeris body ID

    Returns:
        bool: True if auto-SPK is enabled for this body
    """
    with _AUTO_SPK_LOCK:
        return ipl in _AUTO_SPK_REGISTRY and _AUTO_SPK_REGISTRY[ipl].enabled


def get_auto_spk_config(ipl: int) -> Optional[AutoSpkConfig]:
    """
    Get the auto-SPK configuration for a body.

    Args:
        ipl: libephemeris body ID

    Returns:
        AutoSpkConfig or None if not configured
    """
    with _AUTO_SPK_LOCK:
        return _AUTO_SPK_REGISTRY.get(ipl)


def list_auto_spk_bodies() -> dict[int, AutoSpkConfig]:
    """
    List all bodies with auto-SPK enabled.

    Returns:
        Dict mapping ipl -> AutoSpkConfig
    """
    with _AUTO_SPK_LOCK:
        return dict(_AUTO_SPK_REGISTRY)


def download_now(ipl: int, force: bool = False) -> str:
    """
    Download SPK for a body immediately.

    This triggers immediate download of the SPK file for a body that has
    auto-SPK enabled, rather than waiting for the first calculation.

    Args:
        ipl: libephemeris body ID (must have auto-SPK enabled)
        force: If True, re-download even if cached file exists

    Returns:
        str: Path to the downloaded SPK file

    Raises:
        ValueError: If auto-SPK is not enabled for this body or the direct
            Horizons download fails.
    """
    config = get_auto_spk_config(ipl)
    if config is None:
        raise ValueError(
            f"Auto-SPK not enabled for body {ipl}. Call enable_auto_spk() first."
        )

    return _ensure_spk_downloaded(config, force=force)


def _ensure_spk_downloaded(config: AutoSpkConfig, force: bool = False) -> str:
    """
    Ensure SPK is downloaded and registered for a body.

    Internal function that handles the download and registration logic.

    Args:
        config: Auto-SPK configuration
        force: If True, re-download even if cached

    Returns:
        str: Path to the SPK file
    """
    from . import spk

    logger = get_logger()
    cache_path = config.get_cache_path()

    # Check if already cached and valid
    if os.path.exists(cache_path) and not force:
        if _is_valid_bsp(cache_path):
            config.spk_path = cache_path
        else:
            # Left in place on purpose: the download below publishes
            # atomically over it, so a failed re-download keeps whatever the
            # user already had instead of leaving an empty cache.
            logger.warning(
                "Cached SPK file %s is corrupted, re-downloading", cache_path
            )
            force = True

    if force or not os.path.exists(cache_path):
        with _AUTO_SPK_LOCK:
            # Re-check under the lock; honor force even if another thread
            # (or a previous run) left a file here.
            if force or not os.path.exists(cache_path):
                _download_spk_astroquery(
                    body_id=config.body_id,
                    start=config.start,
                    end=config.end,
                    output_path=cache_path,
                )
        config.spk_path = cache_path

    # Resolve the NAIF ID: explicit > detected from the file > deduced.
    # Horizons-generated kernels may use 20000000+N target IDs while the
    # classic convention is 2000000+N — reading the file is authoritative.
    naif_id = config.naif_id
    detected = _detect_naif_id_from_file(cache_path, naif_id)
    if detected is not None:
        naif_id = detected
    if naif_id is None:
        naif_id = spk._deduce_naif_id(config.body_id)
        if naif_id is None:
            raise ValueError(
                f"Cannot deduce NAIF ID for '{config.body_id}'. "
                "Please provide naif_id explicitly in enable_auto_spk()."
            )
    config.naif_id = naif_id

    # Register the SPK body, or re-register if it now resolves to a different
    # (path, NAIF id). A forced re-download with a wider date range produces a
    # new cache path (the filename hashes the range), and _detect_naif_id_from_file
    # may correct naif_id for the same path; a presence- or path-only check would
    # leave the stale registration in place and calc_ut would keep using it.
    from . import state

    existing = state._SPK_BODY_MAP.get(config.ipl)
    if existing is None or existing != (cache_path, naif_id):
        spk.register_spk_body(config.ipl, cache_path, naif_id)

    return cache_path


def try_auto_download(ipl: int) -> Optional[str]:
    """
    Try to auto-download SPK for a body if configured.

    This function is called internally during position calculations
    to transparently download SPK files when needed.

    Args:
        ipl: libephemeris body ID

    Returns:
        Path to SPK file if downloaded, None if not configured or failed
    """
    config = get_auto_spk_config(ipl)
    if config is None or not config.enabled:
        return None

    try:
        return _ensure_spk_downloaded(config)
    except (OSError, ValueError, KeyError, ImportError, Error):
        # Log error but don't raise - fall back to Keplerian
        # (Error covers SPKNotFoundError raised by registration)
        return None


def clear_cache(ipl: Optional[int] = None) -> int:
    """
    Clear cached SPK files.

    Args:
        ipl: If provided, only clear cache for this body.
             If None, clear all cached SPK files.

    Returns:
        int: Number of files deleted
    """
    deleted = 0

    with _AUTO_SPK_LOCK:
        if ipl is not None:
            config = _AUTO_SPK_REGISTRY.get(ipl)
            if config and config.spk_path and os.path.exists(config.spk_path):
                os.remove(config.spk_path)
                config.spk_path = None
                deleted = 1
        else:
            # Clear all cached files
            for config in _AUTO_SPK_REGISTRY.values():
                if config.spk_path and os.path.exists(config.spk_path):
                    try:
                        os.remove(config.spk_path)
                        config.spk_path = None
                        deleted += 1
                    except OSError:
                        pass

    return deleted


def get_cache_info() -> dict[str, Union[int, float, str, list[str]]]:
    """
    Get information about the SPK cache.

    Returns:
        Dict with cache information:
            - cache_dir: Path to cache directory
            - num_files: Number of cached SPK files
            - total_size_mb: Total size in MB
            - files: List of cached file names
    """
    cache_dir = DEFAULT_AUTO_SPK_DIR

    if not os.path.exists(cache_dir):
        return {
            "cache_dir": cache_dir,
            "num_files": 0,
            "total_size_mb": 0.0,
            "files": [],
        }

    files = [f for f in os.listdir(cache_dir) if f.endswith(".bsp")]
    total_size = sum(os.path.getsize(os.path.join(cache_dir, f)) for f in files)

    return {
        "cache_dir": cache_dir,
        "num_files": len(files),
        "total_size_mb": round(total_size / (1024 * 1024), 2),
        "files": files,
    }


# =============================================================================
# CACHE MANAGEMENT FUNCTIONS
# =============================================================================


def list_cached_spk(
    cache_dir: Optional[str] = None,
) -> list[dict[str, Union[str, float, int, None]]]:
    """
    List all cached SPK files with their date ranges and sizes.

    This function scans the cache directory for SPK files and returns detailed
    information about each file including its size, date range coverage, and
    last access time.

    Args:
        cache_dir: Optional custom cache directory path. If None, checks both:
            - The default cache directory ({library_path}/spk_cache/)
            - The auto SPK directory (~/.libephemeris/spk/)

    Returns:
        list: A list of dictionaries, each containing:
            - filename (str): Name of the SPK file
            - path (str): Full path to the file
            - size_bytes (int): File size in bytes
            - size_mb (float): File size in megabytes
            - jd_start (float or None): Start Julian Day from SPK coverage (if readable)
            - jd_end (float or None): End Julian Day from SPK coverage (if readable)
            - date_start (str or None): Start date in YYYY-MM-DD format (if readable)
            - date_end (str or None): End date in YYYY-MM-DD format (if readable)
            - last_accessed (float): Last access time as Unix timestamp

    Example:
        >>> from libephemeris.spk_auto import list_cached_spk
        >>> cached_files = list_cached_spk()
        >>> for f in cached_files:
        ...     print(f"{f['filename']}: {f['size_mb']:.2f} MB")
        ...     if f['date_start']:
        ...         print(f"  Coverage: {f['date_start']} to {f['date_end']}")
    """
    from .spk import get_spk_coverage

    result = []
    dirs_to_check = []

    if cache_dir is not None:
        dirs_to_check = [os.path.abspath(cache_dir)]
    else:
        # Check both the resolved default and the legacy cache directory
        from .state import get_spk_cache_dir

        resolved_dir = get_spk_cache_dir() or DEFAULT_AUTO_SPK_DIR
        default_cache = DEFAULT_AUTO_SPK_DIR
        # Dedupe: in the default configuration both entries resolve to the
        # same directory, which would double-count sizes and listings.
        dirs_to_check = list(dict.fromkeys([default_cache, resolved_dir]))

    for dir_path in dirs_to_check:
        if not os.path.exists(dir_path):
            continue

        for filename in os.listdir(dir_path):
            if not filename.endswith(".bsp"):
                continue

            filepath = os.path.join(dir_path, filename)
            stat_info = os.stat(filepath)

            file_info: dict[str, Union[str, float, int, None]] = {
                "filename": filename,
                "path": filepath,
                "size_bytes": stat_info.st_size,
                "size_mb": round(stat_info.st_size / (1024 * 1024), 2),
                "jd_start": None,
                "jd_end": None,
                "date_start": None,
                "date_end": None,
                "last_accessed": stat_info.st_atime,
            }

            # Try to get actual coverage from SPK file
            try:
                coverage = get_spk_coverage(filepath)
                if coverage is not None:
                    jd_start, jd_end = coverage
                    file_info["jd_start"] = jd_start
                    file_info["jd_end"] = jd_end
                    file_info["date_start"] = _jd_to_iso_date(jd_start)
                    file_info["date_end"] = _jd_to_iso_date(jd_end)
            except (OSError, ValueError, KeyError, ImportError):
                # If we can't read coverage, try to extract from filename
                parts = filename[:-4].split("_")  # Remove .bsp
                if len(parts) >= 3:
                    try:
                        file_jd_start = float(parts[-2])
                        file_jd_end = float(parts[-1])
                        file_info["jd_start"] = file_jd_start
                        file_info["jd_end"] = file_jd_end
                        file_info["date_start"] = _jd_to_iso_date(file_jd_start)
                        file_info["date_end"] = _jd_to_iso_date(file_jd_end)
                    except ValueError:
                        pass

            result.append(file_info)

    return result


def clear_spk_cache(cache_dir: Optional[str] = None) -> int:
    """
    Delete all cached SPK files.

    This function removes all .bsp files from the cache directory. It does not
    delete the directory itself, only the SPK files within it.

    Warning: This operation cannot be undone. All cached SPK files will need
    to be re-downloaded if needed.

    Args:
        cache_dir: Optional custom cache directory path. If None, clears both:
            - The default cache directory ({library_path}/spk_cache/)
            - The auto SPK directory (~/.libephemeris/spk/)

    Returns:
        int: Number of files deleted

    Example:
        >>> from libephemeris.spk_auto import clear_spk_cache, get_cache_size
        >>> print(f"Cache size before: {get_cache_size():.2f} MB")
        Cache size before: 15.32 MB
        >>> deleted = clear_spk_cache()
        >>> print(f"Deleted {deleted} files")
        Deleted 5 files
        >>> print(f"Cache size after: {get_cache_size():.2f} MB")
        Cache size after: 0.00 MB
    """
    deleted_count = 0
    dirs_to_clear = []

    if cache_dir is not None:
        dirs_to_clear = [os.path.abspath(cache_dir)]
    else:
        # Clear both the resolved default and the legacy cache directory
        from .state import get_spk_cache_dir

        resolved_dir = get_spk_cache_dir() or DEFAULT_AUTO_SPK_DIR
        default_cache = DEFAULT_AUTO_SPK_DIR
        dirs_to_clear = [default_cache, resolved_dir]

    for dir_path in dirs_to_clear:
        if not os.path.exists(dir_path):
            continue

        for filename in os.listdir(dir_path):
            if not filename.endswith(".bsp"):
                continue

            filepath = os.path.join(dir_path, filename)
            try:
                os.remove(filepath)
                deleted_count += 1
            except OSError:
                # Skip files that can't be deleted (e.g., permission issues)
                pass

    # Also clear the spk_path from any registered configurations
    with _AUTO_SPK_LOCK:
        for config in _AUTO_SPK_REGISTRY.values():
            if config.spk_path and not os.path.exists(config.spk_path):
                config.spk_path = None

    return deleted_count


def get_cache_size(cache_dir: Optional[str] = None) -> float:
    """
    Get the total disk usage of the SPK cache in megabytes.

    This function calculates the combined size of all .bsp files in the
    cache directory.

    Args:
        cache_dir: Optional custom cache directory path. If None, calculates
            the combined size of both:
            - The default cache directory ({library_path}/spk_cache/)
            - The auto SPK directory (~/.libephemeris/spk/)

    Returns:
        float: Total size of cached SPK files in megabytes

    Example:
        >>> from libephemeris.spk_auto import get_cache_size
        >>> size = get_cache_size()
        >>> print(f"Total cache size: {size:.2f} MB")
        Total cache size: 15.32 MB
        >>> if size > 100:
        ...     print("Consider pruning old cache files")
    """
    total_size = 0
    dirs_to_check = []

    if cache_dir is not None:
        dirs_to_check = [os.path.abspath(cache_dir)]
    else:
        # Check both the resolved default and the legacy cache directory
        from .state import get_spk_cache_dir

        resolved_dir = get_spk_cache_dir() or DEFAULT_AUTO_SPK_DIR
        default_cache = DEFAULT_AUTO_SPK_DIR
        # Dedupe: in the default configuration both entries resolve to the
        # same directory, which would double-count sizes and listings.
        dirs_to_check = list(dict.fromkeys([default_cache, resolved_dir]))

    for dir_path in dirs_to_check:
        if not os.path.exists(dir_path):
            continue

        for filename in os.listdir(dir_path):
            if not filename.endswith(".bsp"):
                continue

            filepath = os.path.join(dir_path, filename)
            try:
                total_size += os.path.getsize(filepath)
            except OSError:
                pass

    return round(total_size / (1024 * 1024), 2)


def prune_old_cache(
    max_age_days: int,
    cache_dir: Optional[str] = None,
) -> int:
    """
    Remove SPK files that have not been accessed recently.

    This function removes cached SPK files whose last access time exceeds
    the specified age threshold. This helps manage disk space by removing
    files that are no longer in active use.

    Args:
        max_age_days: Maximum age in days since last access. Files not accessed
            within this period will be deleted. Must be a positive integer.
        cache_dir: Optional custom cache directory path. If None, prunes both:
            - The default cache directory ({library_path}/spk_cache/)
            - The auto SPK directory (~/.libephemeris/spk/)

    Returns:
        int: Number of files deleted

    Raises:
        ValueError: If max_age_days is not a positive integer

    Example:
        >>> from libephemeris.spk_auto import prune_old_cache, list_cached_spk
        >>> # Remove files not accessed in the last 30 days
        >>> pruned = prune_old_cache(max_age_days=30)
        >>> print(f"Removed {pruned} old cache files")
        Removed 3 old cache files
        >>> # Check remaining files
        >>> remaining = list_cached_spk()
        >>> print(f"Remaining cached files: {len(remaining)}")

    Note:
        The access time is based on the file's last access timestamp (atime),
        which may not be updated on all file systems. Some file systems mount
        with 'noatime' for performance, which would prevent this from working
        correctly.
    """
    import time

    if not isinstance(max_age_days, int) or max_age_days <= 0:
        raise ValueError("max_age_days must be a positive integer")

    deleted_count = 0
    current_time = time.time()
    max_age_seconds = max_age_days * 24 * 60 * 60
    dirs_to_prune = []

    if cache_dir is not None:
        dirs_to_prune = [os.path.abspath(cache_dir)]
    else:
        # Prune both the resolved default and the legacy cache directory
        from .state import get_spk_cache_dir

        resolved_dir = get_spk_cache_dir() or DEFAULT_AUTO_SPK_DIR
        default_cache = DEFAULT_AUTO_SPK_DIR
        dirs_to_prune = [default_cache, resolved_dir]

    for dir_path in dirs_to_prune:
        if not os.path.exists(dir_path):
            continue

        for filename in os.listdir(dir_path):
            if not filename.endswith(".bsp"):
                continue

            filepath = os.path.join(dir_path, filename)
            try:
                stat_info = os.stat(filepath)
                last_access = stat_info.st_atime
                age_seconds = current_time - last_access

                if age_seconds > max_age_seconds:
                    os.remove(filepath)
                    deleted_count += 1
            except OSError:
                # Skip files that can't be accessed or deleted
                pass

    # Also clear the spk_path from any registered configurations for deleted files
    with _AUTO_SPK_LOCK:
        for config in _AUTO_SPK_REGISTRY.values():
            if config.spk_path and not os.path.exists(config.spk_path):
                config.spk_path = None

    return deleted_count


# =============================================================================
# PRESET CONFIGURATIONS
# =============================================================================


def enable_common_bodies(
    start: str = DEFAULT_START_DATE,
    end: str = DEFAULT_END_DATE,
    cache_dir: Optional[str] = None,
) -> None:
    """
    Enable auto-SPK for commonly used minor bodies.

    This is a convenience function that enables automatic SPK downloads
    for popular asteroids and TNOs used in astrological and astronomical
    calculations.

    Bodies enabled:
        - Chiron (2060)
        - Pholus (5145)
        - Ceres (1)
        - Pallas (2)
        - Juno (3)
        - Vesta (4)
        - Eris (136199)
        - Sedna (90377)

    Args:
        start: Start date for SPK coverage (YYYY-MM-DD)
        end: End date for SPK coverage (YYYY-MM-DD)
        cache_dir: Directory for storing downloaded SPK files

    Example:
        >>> from libephemeris import spk_auto
        >>> spk_auto.enable_common_bodies()
        >>> # Now all common minor bodies use SPK automatically
    """
    from .constants import (
        CHIRON,
        PHOLUS,
        CERES,
        PALLAS,
        JUNO,
        VESTA,
        ERIS,
        SEDNA,
        SPK_BODY_NAME_MAP,
    )

    # Use the canonical Horizons command strings from SPK_BODY_NAME_MAP:
    # bare numbers like "1".."4" collide with planet-barycenter indices on
    # Horizons ("SPK creation is not available for pre-computed objects"),
    # so Ceres..Vesta need the "Name;" small-body syntax.
    for ipl in (CHIRON, PHOLUS, CERES, PALLAS, JUNO, VESTA, ERIS, SEDNA):
        body_id, naif_id = SPK_BODY_NAME_MAP[ipl]
        enable_auto_spk(
            ipl=ipl,
            body_id=body_id,
            naif_id=naif_id,
            start=start,
            end=end,
            cache_dir=cache_dir,
        )


def disable_all() -> None:
    """
    Disable auto-SPK for all bodies.

    This clears the auto-SPK registry but does not delete cached files
    or unregister already-registered SPK bodies.
    """
    with _AUTO_SPK_LOCK:
        _AUTO_SPK_REGISTRY.clear()


def _jd_to_iso_date(jd: float) -> str:
    """
    Convert Julian Day to ISO date string (YYYY-MM-DD).

    Args:
        jd: Julian Day number

    Returns:
        str: ISO date string (YYYY-MM-DD)
    """
    # Julian Day to calendar date conversion
    # Algorithm from Meeus, Astronomical Algorithms.
    # Use the proleptic Gregorian calendar for ALL dates so this is the exact
    # inverse of _iso_to_jd (which always applies the Gregorian correction
    # b = 2 - a + a//4). A z < 2299161 Julian branch would drift these two
    # apart by ~10 days for pre-1582 dates (e.g. medium/extended tier starts),
    # corrupting the SPK request window and cache filenames.
    # Meeus's INT() is a FLOOR: int() truncates toward zero, which breaks
    # every negative intermediate (JD < 0, i.e. years before -4712) and
    # silently returned dates hundreds of days off there. math.floor keeps
    # this the exact inverse of _iso_to_jd over the whole range.
    import math

    jd = jd + 0.5
    z = math.floor(jd)
    # The sub-day fraction (jd - z) is intentionally discarded: coverage
    # windows and cache filenames are day-granular, so the result is the
    # civil date containing the instant.

    alpha = math.floor((z - 1867216.25) / 36524.25)
    a = z + 1 + alpha - math.floor(alpha / 4)

    b = a + 1524
    c = math.floor((b - 122.1) / 365.25)
    d = math.floor(365.25 * c)
    e = math.floor((b - d) / 30.6001)

    day = b - d - math.floor(30.6001 * e)

    if e < 14:
        month = e - 1
    else:
        month = e - 13

    if month > 2:
        year = c - 4716
    else:
        year = c - 4715

    return f"{year:04d}-{month:02d}-{day:02d}"


def _iso_to_jd(date_str: str) -> float:
    """
    Convert ISO date string (YYYY-MM-DD) to Julian Day.

    Uses the Meeus algorithm (Astronomical Algorithms) for the conversion.
    This is the inverse of _jd_to_iso_date().

    Args:
        date_str: ISO date string in YYYY-MM-DD format

    Returns:
        float: Julian Day number

    Example:
        >>> _iso_to_jd("2000-01-01")
        2451544.5
        >>> _iso_to_jd("1550-01-01")
        2287185.5
    """
    # A BCE date from _jd_to_iso_date renders the year with a leading '-'
    # (e.g. "-1975-11-07"); strip it before splitting so the year sign is not
    # mistaken for a field separator (which would make int("") raise).
    neg_year = date_str.startswith("-")
    core = date_str[1:] if neg_year else date_str
    parts = core.split("-")
    year, month, day = int(parts[0]), int(parts[1]), int(parts[2])
    if neg_year:
        year = -year

    # Algorithm from Meeus, Astronomical Algorithms. Meeus's INT() is a
    # FLOOR, not a truncation: Python's int() rounds toward zero, which for
    # negative intermediates (BCE years, and year < -4716 in the 365.25
    # term) yields a JD one day too large. math.floor matches erfa.cal2jd
    # across the whole proleptic-Gregorian range.
    import math

    if month <= 2:
        year -= 1
        month += 12

    a = math.floor(year / 100)
    b = 2 - a + math.floor(a / 4)

    jd = (
        math.floor(365.25 * (year + 4716))
        + math.floor(30.6001 * (month + 1))
        + day
        + b
        - 1524.5
    )
    return jd


def _generate_spk_cache_filename(
    body_id: Union[int, str], jd_start: float, jd_end: float
) -> str:
    """
    Generate a cache filename for an SPK file.

    The filename includes body identifier and truncated Julian Day range
    to enable cache lookup.

    Args:
        body_id: Body identifier
        jd_start: Start Julian Day
        jd_end: End Julian Day

    Returns:
        str: Filename like "2060_2458849_2462502.bsp"
    """
    # Sanitize body_id for filename (shared writer/reader sanitizer)
    from .spk import _sanitize_filename

    safe_body_id = _sanitize_filename(body_id)

    # Truncate JD to integer for filename (sufficient precision for cache lookup)
    jd_start_int = int(jd_start)
    jd_end_int = int(jd_end)

    return f"{safe_body_id}_{jd_start_int}_{jd_end_int}.bsp"


def _find_covering_spk(
    body_id: Union[int, str],
    jd_start: float,
    jd_end: float,
    cache_dir: str,
) -> Optional[str]:
    """
    Find an existing SPK file in the cache that covers the requested date range.

    Searches the cache directory for any SPK file for the given body that
    covers the requested Julian Day range.

    Args:
        body_id: Body identifier
        jd_start: Start Julian Day
        jd_end: End Julian Day
        cache_dir: Cache directory path

    Returns:
        str: Path to covering SPK file if found, None otherwise
    """
    if not os.path.exists(cache_dir):
        return None

    # Sanitize body_id to match filename pattern (shared writer/reader sanitizer)
    from .spk import _sanitize_filename

    safe_body_id = _sanitize_filename(body_id)

    # Read the ACTUAL coverage from each candidate kernel's segment metadata
    # rather than parsing the filename. Cached files appear under several naming
    # schemes and only one encodes the JD range:
    #   - body_<jdstart>_<jdend>.bsp  (_generate_spk_cache_filename)
    #   - body_<YYYYMM>_<YYYYMM>.bsp  (download_spk)
    #   - body_<md5hash>.bsp          (AutoSpkConfig.get_cache_filename)
    # The old filename-parsing path recognised only the first, so YYYYMM- and
    # hash-named kernels were never detected as covering and triggered
    # redundant re-downloads. get_spk_coverage() reads the real time span from
    # the kernel, so all three schemes are handled correctly.
    from .spk import (
        _SPK_COVERAGE_END_TOLERANCE_DAYS,
        _SPK_COVERAGE_START_TOLERANCE_DAYS,
        get_spk_coverage,
    )

    for filename in os.listdir(cache_dir):
        if not filename.endswith(".bsp"):
            continue

        # Every scheme prefixes the sanitized body id followed by "_".
        if not filename.startswith(safe_body_id + "_"):
            continue

        filepath = os.path.abspath(os.path.join(cache_dir, filename))

        try:
            coverage = get_spk_coverage(filepath)
        except (OSError, ValueError, KeyError, IndexError, ImportError):
            continue
        if coverage is None:
            continue

        file_jd_start, file_jd_end = coverage
        # Does this kernel actually cover the requested range? The coverage
        # is the usable span (its start sits one light-time inside the
        # stored data), so compare with the documented edge tolerances.
        if (
            file_jd_start <= jd_start + _SPK_COVERAGE_START_TOLERANCE_DAYS
            and file_jd_end >= jd_end - _SPK_COVERAGE_END_TOLERANCE_DAYS
        ):
            if _is_valid_bsp(filepath):
                return filepath
            get_logger().warning("Cached SPK file %s is corrupted, skipping", filepath)

    return None


def is_spk_cached(
    body_id: Union[int, str],
    jd_start: float,
    jd_end: float,
    cache_dir: Optional[str] = None,
) -> bool:
    """
    Check if a cached SPK file exists and covers the requested date range.

    This function searches the cache directory for SPK files matching the given
    body_id and verifies that the file covers the requested Julian Day range by
    parsing the SPK file header (not just relying on filename patterns).

    Args:
        body_id: JPL Horizons body identifier. Can be:
            - Asteroid number (int): 2060, 136199
            - Name (str): "Chiron", "Eris"
            - Combined: "2060 Chiron"
        jd_start: Start of the date range (Julian Day)
        jd_end: End of the date range (Julian Day)
        cache_dir: Directory to search for SPK files.
            Default: ~/.libephemeris/spk/

    Returns:
        bool: True if a cached SPK file exists that covers the requested range,
              False otherwise.

    Example:
        >>> from libephemeris.spk_auto import is_spk_cached
        >>> # Check if Chiron SPK is cached for 2020-2030
        >>> jd_start = 2458849.5  # 2020-01-01
        >>> jd_end = 2462502.5    # 2030-01-01
        >>> if is_spk_cached("2060", jd_start, jd_end):
        ...     print("Chiron SPK is cached")
        ... else:
        ...     print("Need to download Chiron SPK")
    """
    # Determine cache directory (resolve via state/env/default)
    if cache_dir is None:
        from .state import get_spk_cache_dir

        cache_dir = get_spk_cache_dir()
    if cache_dir is None:
        cache_dir = DEFAULT_AUTO_SPK_DIR

    if not os.path.exists(cache_dir):
        return False

    # Sanitize body_id to match filename pattern (shared writer/reader sanitizer)
    from .spk import (
        _SPK_COVERAGE_END_TOLERANCE_DAYS,
        _SPK_COVERAGE_START_TOLERANCE_DAYS,
        _sanitize_filename,
        get_spk_coverage,
    )

    safe_body_id = _sanitize_filename(body_id)

    # Look for files matching this body
    for filename in os.listdir(cache_dir):
        if not filename.endswith(".bsp"):
            continue

        # Check if filename matches our body
        if not filename.startswith(safe_body_id + "_"):
            continue

        # Found a potential match - parse the actual SPK file to get coverage
        spk_path = os.path.join(cache_dir, filename)

        try:
            coverage = get_spk_coverage(spk_path)
            if coverage is not None:
                file_jd_start, file_jd_end = coverage
                # Check if this file covers our requested range, with the
                # same edge tolerances as _find_covering_spk (the coverage
                # is the usable span, not the stored one).
                if (
                    file_jd_start <= jd_start + _SPK_COVERAGE_START_TOLERANCE_DAYS
                    and file_jd_end >= jd_end - _SPK_COVERAGE_END_TOLERANCE_DAYS
                ):
                    return True
        except (OSError, ValueError, KeyError, ImportError):
            # If we can't read the file, skip it and try others
            continue

    return False


def auto_get_spk(
    body_id: Union[int, str],
    jd_start: float,
    jd_end: float,
    cache_dir: Optional[str] = None,
    ipl: Optional[int] = None,
    naif_id: Optional[int] = None,
) -> str:
    """
    Automatically get SPK file for a body and date range, downloading if needed.

    This function checks if an SPK file for the requested body and date range
    exists in the local cache directory. If not, it automatically downloads
    the file with the direct JPL Horizons HTTP client.

    If ``ipl`` is provided, the SPK body is automatically registered with
    libephemeris after download, making it immediately usable by ``calc_ut()``.

    Args:
        body_id: JPL Horizons body identifier. Can be:
            - Asteroid number (int): 2060, 136199
            - Name (str): "Chiron", "Eris"
            - Combined: "2060 Chiron"
        jd_start: Start of the date range (Julian Day)
        jd_end: End of the date range (Julian Day)
        cache_dir: Directory for storing/finding SPK files.
            Default: Uses set_spk_cache_dir() value if set, otherwise
            ~/.libephemeris/spk/
        ipl: Optional libephemeris body ID (e.g., CHIRON). If provided,
            the SPK body is automatically registered after download.
        naif_id: Optional NAIF ID for SPK lookup. If not provided but ipl is,
            the NAIF ID is deduced from body_id using the convention:
            naif_id = asteroid_number + 2000000

    Returns:
        str: Path to the SPK file (either existing cached or newly downloaded)

    Raises:
        ValueError: If the body is not found, download fails, or naif_id cannot be
            deduced when ipl is provided

    Note:
        The date range may be extended by the value set via set_spk_date_padding()
        to provide a buffer around the requested range.

    Example:
        >>> from libephemeris.spk_auto import auto_get_spk
        >>> from libephemeris.constants import CHIRON
        >>> # Get SPK for Chiron covering 2020-2030 and register it
        >>> jd_start = 2458849.5  # 2020-01-01
        >>> jd_end = 2462502.5    # 2030-01-01
        >>> spk_path = auto_get_spk("2060", jd_start, jd_end, ipl=CHIRON)
        >>> print(spk_path)
        /home/user/.libephemeris/spk/2060_2458849_2462502.bsp
        >>> # Now calc_ut(jd, CHIRON, ...) uses SPK data automatically
    """
    # Validate inputs
    if jd_end <= jd_start:
        raise ValueError(
            f"jd_end ({jd_end}) must be greater than jd_start ({jd_start})"
        )

    # Import state to get configuration settings
    from . import state

    # Determine cache directory (parameter > global config > default)
    if cache_dir is None:
        cache_dir = state.get_spk_cache_dir()
    if cache_dir is None:
        cache_dir = DEFAULT_AUTO_SPK_DIR

    # Apply date padding from global configuration
    date_padding = state.get_spk_date_padding()
    if date_padding > 0:
        jd_start = jd_start - date_padding
        jd_end = jd_end + date_padding

    # Create cache directory if it doesn't exist
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir, exist_ok=True)

    # Check if we already have a cached SPK that covers this range
    existing_spk = _find_covering_spk(body_id, jd_start, jd_end, cache_dir)
    if existing_spk is not None:
        # Even if cached, register if ipl is provided
        if ipl is not None:
            _register_spk_after_download(existing_spk, body_id, ipl, naif_id)
        return existing_spk

    # Need to download - convert Julian Days to date strings
    start_date = _jd_to_iso_date(jd_start)
    end_date = _jd_to_iso_date(jd_end)

    # Generate cache filename
    filename = _generate_spk_cache_filename(body_id, jd_start, jd_end)
    output_path = os.path.join(cache_dir, filename)

    # Download through the direct JPL Horizons HTTP client.  The helper keeps
    # its historical private name for compatibility.
    with _AUTO_SPK_LOCK:
        # Double-check after acquiring lock (another thread may have downloaded)
        if os.path.exists(output_path) and _is_valid_bsp(output_path):
            # Register if ipl is provided
            if ipl is not None:
                _register_spk_after_download(output_path, body_id, ipl, naif_id)
            return output_path

        # Check again for covering SPK (may have been downloaded with different range)
        existing_spk = _find_covering_spk(body_id, jd_start, jd_end, cache_dir)
        if existing_spk is not None:
            # Register if ipl is provided
            if ipl is not None:
                _register_spk_after_download(existing_spk, body_id, ipl, naif_id)
            return existing_spk

        # Download SPK file
        _download_spk_astroquery(
            body_id=str(body_id),
            start=start_date,
            end=end_date,
            output_path=output_path,
        )

    # Register the SPK body if ipl is provided
    if ipl is not None:
        _register_spk_after_download(output_path, body_id, ipl, naif_id)

    return output_path


def _register_spk_after_download(
    spk_path: str,
    body_id: Union[int, str],
    ipl: int,
    naif_id: Optional[int] = None,
) -> None:
    """
    Register an SPK body after download.

    Internal helper function that registers the SPK body with libephemeris
    so it can be used by calc_ut().

    Args:
        spk_path: Path to the downloaded SPK file
        body_id: JPL Horizons body identifier (used for NAIF ID deduction)
        ipl: libephemeris body ID (e.g., CHIRON)
        naif_id: Optional NAIF ID. If not provided, deduced from body_id.

    Raises:
        ValueError: If naif_id cannot be deduced and is not provided
    """
    from . import spk

    # The file's own target IDs are authoritative (Horizons kernels may use
    # 20000000+N where the classic convention is 2000000+N).
    detected = _detect_naif_id_from_file(spk_path, naif_id)
    if detected is not None:
        naif_id = detected

    # Deduce NAIF ID if still unknown
    if naif_id is None:
        naif_id = spk._deduce_naif_id(str(body_id))
        if naif_id is None:
            raise ValueError(
                f"Cannot deduce NAIF ID for '{body_id}'. "
                f"Please provide naif_id explicitly. "
                f"For numbered asteroids: naif_id = asteroid_number + 2000000"
            )

    # Register the SPK body
    spk.register_spk_body(ipl, spk_path, naif_id)


# =============================================================================
# JPL HORIZONS SPK DOWNLOAD (PUBLIC API)
# =============================================================================


def download_spk_from_horizons(
    body_id: Union[int, str],
    jd_start: float,
    jd_end: float,
    output_path: str,
    location: str = "@0",
    ipl: Optional[int] = None,
    naif_id: Optional[int] = None,
) -> str:
    """
    Download an SPK file from JPL Horizons for a specified body and date range.

    This function uses the JPL Horizons HTTP API to download an SPK (SPICE kernel)
    file for the specified body and Julian Day range. The file can then be used
    for high-precision ephemeris calculations.

    If ``ipl`` is provided, the SPK body is automatically registered with
    libephemeris after download, making it immediately usable by ``calc_ut()``.

    Args:
        body_id: JPL Horizons body identifier. Can be:
            - Asteroid number (int or str): 2060, "2060", "136199"
            - Name (str): "Chiron", "Eris"
            - Designation (str): "2003 UB313"
        jd_start: Start of the date range (Julian Day)
        jd_end: End of the date range (Julian Day)
        output_path: Full path where the SPK file should be saved (including filename)
        location: Observer location code for Horizons (default: "@0" = solar system barycenter)
        ipl: Optional libephemeris body ID (e.g., CHIRON). If provided,
            the SPK body is automatically registered after download.
        naif_id: Optional NAIF ID for SPK lookup. If not provided but ipl is,
            the NAIF ID is deduced from body_id using the convention:
            naif_id = asteroid_number + 2000000

    Returns:
        str: Path to the downloaded SPK file (same as output_path)

    Raises:
        ValueError: If body is not found on Horizons, date range is invalid,
            the date range is too large for Horizons to process, or naif_id
            cannot be deduced when ipl is provided
        ConnectionError: If network request fails

    Example:
        >>> from libephemeris.spk_auto import download_spk_from_horizons
        >>> from libephemeris.constants import CHIRON
        >>> # Download Chiron SPK for 2020-2030 and register it
        >>> jd_start = 2458849.5  # 2020-01-01
        >>> jd_end = 2462502.5    # 2030-01-01
        >>> path = download_spk_from_horizons(
        ...     body_id="2060",
        ...     jd_start=jd_start,
        ...     jd_end=jd_end,
        ...     output_path="/path/to/chiron.bsp",
        ...     ipl=CHIRON,  # Automatically registers the SPK
        ... )
        >>> print(path)
        /path/to/chiron.bsp
        >>> # Now calc_ut(jd, CHIRON, ...) uses SPK data automatically

    Notes:
        - Horizons has limits on the date range for SPK generation. For most
          bodies, a 100-year range is typically acceptable.
        - For numbered asteroids, use the asteroid number (e.g., "2060" for Chiron).
        - For major planets, use the body name (e.g., "Mars", "Jupiter").
    """
    # Validate inputs
    if jd_end <= jd_start:
        raise ValueError(
            f"jd_end ({jd_end}) must be greater than jd_start ({jd_start})"
        )

    # Convert Julian Days to ISO date strings
    start_date = _jd_to_iso_date(jd_start)
    end_date = _jd_to_iso_date(jd_end)

    # Ensure output directory exists
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    from . import spk as _spk_mod

    logger = get_logger()
    body_id_str = str(body_id)

    # Map location format for spk.download_spk(): "@0" -> "500@0"
    if location.startswith("@"):
        center = f"500{location}"
    else:
        center = location

    try:
        # Download SPK via direct JPL Horizons HTTP API
        result_path = _spk_mod.download_spk(
            body=body_id_str,
            start=start_date,
            end=end_date,
            path=output_dir or ".",
            center=center,
            overwrite=True,
        )

        # Rename to the requested output_path if the API generated a
        # different filename (download_spk uses its own naming convention)
        if os.path.abspath(result_path) != os.path.abspath(output_path):
            os.replace(result_path, output_path)

        logger.info("SPK saved: %s", output_path)

    except ValueError as e:
        # Body not found or invalid parameters from the Horizons API
        raise ValueError(
            f"Failed to download SPK for '{body_id}' from JPL Horizons: {e}"
        ) from e
    except ConnectionError as e:
        raise ConnectionError(
            f"Network error downloading SPK for '{body_id}': {e}"
        ) from e
    except (OSError, ValueError, KeyError, ImportError) as e:
        raise ValueError(
            f"Failed to download SPK for '{body_id}' from JPL Horizons: {e}"
        ) from e

    # Register the SPK body if ipl is provided
    if ipl is not None:
        _register_spk_after_download(output_path, body_id_str, ipl, naif_id)

    return output_path


# =============================================================================
# LOCAL SPK DISCOVERY
# =============================================================================


def discover_local_spks(path: str) -> Dict[str, str]:
    """
    Scan a directory for SPK files and register known minor bodies.

    This function reads the NAIF IDs contained in each .bsp file and
    matches them against the bodies in SPK_BODY_NAME_MAP. Matched bodies
    are automatically registered for use in calculations.

    No network requests are made. This is a purely local file operation.

    Args:
        path: Directory path to scan for .bsp files.

    Returns:
        Dict mapping body name -> status string for each processed body.
        Possible status values: "registered", "already_registered",
        "error: <message>".

    Example:
        >>> from libephemeris.spk_auto import discover_local_spks
        >>> results = discover_local_spks("/path/to/ephemeris/")
        >>> for body, status in results.items():
        ...     print(f"{body}: {status}")
        Chiron: registered
        Ceres: registered
    """
    from . import spk, state
    from .constants import SPK_BODY_NAME_MAP

    logger = get_logger()
    results: Dict[str, str] = {}

    if not os.path.isdir(path):
        logger.debug("SPK discovery path does not exist: %s", path)
        return results

    # Collect all .bsp files
    bsp_files = [f for f in os.listdir(path) if f.lower().endswith(".bsp")]
    if not bsp_files:
        logger.debug("No .bsp files found in %s", path)
        return results

    # Build a reverse map: naif_id -> (ipl, horizons_id, body_name)
    naif_to_body: Dict[int, tuple] = {}
    for ipl, (horizons_id, naif_id) in SPK_BODY_NAME_MAP.items():
        body_name = spk._get_body_name(ipl) or horizons_id
        naif_to_body[naif_id] = (ipl, horizons_id, body_name)
        # Also map the Horizons convention (20000000+N)
        # in case the file uses that instead of the 2000000+N convention
        asteroid_number = naif_id - 2000000
        if asteroid_number > 0:
            horizons_naif = asteroid_number + 20000000
            if horizons_naif not in naif_to_body:
                naif_to_body[horizons_naif] = (ipl, horizons_id, body_name)

    logger.debug(
        "Scanning %d .bsp files in %s for known bodies...",
        len(bsp_files),
        path,
    )

    for bsp_file in bsp_files:
        filepath = os.path.join(path, bsp_file)

        # Read NAIF target IDs from the SPK file
        try:
            targets = spk._get_spk_targets(filepath)
        except (OSError, ValueError, KeyError, ImportError):
            continue

        if not targets:
            continue

        for target_naif in targets:
            if target_naif not in naif_to_body:
                continue

            ipl, horizons_id, body_name = naif_to_body[target_naif]

            # If already registered, check if the new file has wider coverage
            if ipl in state._SPK_BODY_MAP:
                existing_file, _ = state._SPK_BODY_MAP[ipl]
                try:
                    existing_coverage = spk.get_spk_coverage(existing_file)
                    new_coverage = spk.get_spk_coverage(filepath)
                    if existing_coverage and new_coverage:
                        existing_span = existing_coverage[1] - existing_coverage[0]
                        new_span = new_coverage[1] - new_coverage[0]
                        if new_span > existing_span:
                            # Re-register with wider-range file
                            spk.register_spk_body(ipl, filepath, target_naif)
                            results[body_name] = "registered"
                            logger.debug(
                                "Re-registered %s with wider SPK %s "
                                "(span %.0f > %.0f days)",
                                body_name,
                                bsp_file,
                                new_span,
                                existing_span,
                            )
                            continue
                except (OSError, ValueError, KeyError, ImportError):
                    pass
                if body_name not in results:
                    results[body_name] = "already_registered"
                continue

            try:
                spk.register_spk_body(ipl, filepath, target_naif)
                results[body_name] = "registered"
                logger.debug(
                    "Discovered and registered %s (NAIF %d) from %s",
                    body_name,
                    target_naif,
                    bsp_file,
                )
            except (OSError, ValueError, KeyError, ImportError) as e:
                results[body_name] = f"error: {e}"
                logger.debug(
                    "Failed to register %s from %s: %s",
                    body_name,
                    bsp_file,
                    e,
                )

    registered_count = sum(1 for v in results.values() if v == "registered")
    if registered_count > 0:
        logger.info(
            "SPK discovery: registered %d bodies from %s",
            registered_count,
            path,
        )

    return results


# =============================================================================
# ENSURE ALL EPHEMERIDES
# =============================================================================


def ensure_all_ephemerides(
    force_download: bool = False,
    show_progress: bool = True,
) -> Dict[str, Union[str, int, dict]]:
    """
    Ensure all ephemeris files are present for the current precision tier.

    Downloads the main ephemeris file and all minor body SPK files needed
    for the current tier if they are not already present/registered.

    This function must be called explicitly by the user. It is NOT triggered
    automatically during library initialization or set_ephe_path().

    Args:
        force_download: If True, re-download even if files exist.
        show_progress: If True, print progress to stdout (if interactive).

    Returns:
        Dict with summary:
            - tier: Current tier name
            - ephemeris: Status of the main ephemeris file
            - spks: Dict mapping body name -> status
            - summary: Dict with counts (cached, downloaded, fallback, errors)

    Example:
        >>> import libephemeris as eph
        >>> eph.set_precision_tier("base")
        >>> results = eph.ensure_all_ephemerides()
        >>> print(f"Tier: {results['tier']}")
        >>> print(f"Summary: {results['summary']}")
    """
    import sys

    from . import spk, state
    from .constants import SPK_BODY_NAME_MAP
    from .state import _get_current_tier, get_spk_date_range_for_tier

    logger = get_logger()
    tier = _get_current_tier()

    logger.info(
        "Ensuring all ephemerides for tier '%s' (%s)",
        tier.name,
        tier.description,
    )

    results: Dict[str, Union[str, int, dict]] = {
        "tier": tier.name,
        "ephemeris": {},
        "spks": {},
        "summary": {
            "total": len(SPK_BODY_NAME_MAP) + 1,
            "cached": 0,
            "downloaded": 0,
            "fallback": 0,
            "errors": 0,
        },
    }
    summary = results["summary"]
    assert isinstance(summary, dict)

    # 1. Ensure main ephemeris file is available
    logger.info("Checking main ephemeris: %s", tier.ephemeris_file)
    eph_file = tier.ephemeris_file

    # Check if the file exists in known locations
    search_paths = []
    if state._EPHEMERIS_PATH:
        search_paths.append(state._EPHEMERIS_PATH)
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    search_paths.append(base_dir)

    eph_found = False
    for sp in search_paths:
        full_path = os.path.join(sp, eph_file)
        if os.path.exists(full_path):
            results["ephemeris"] = {
                "file": eph_file,
                "status": "cached",
                "path": full_path,
            }
            summary["cached"] += 1
            eph_found = True
            break

    if not eph_found or force_download:
        try:
            load = state.get_loader()
            load(eph_file)
            results["ephemeris"] = {"file": eph_file, "status": "downloaded"}
            summary["downloaded"] += 1
        except (OSError, ValueError, KeyError, ImportError) as e:
            results["ephemeris"] = {"file": eph_file, "status": f"error: {e}"}
            summary["errors"] += 1

    # 2. Ensure all minor body SPKs
    start_date, end_date = get_spk_date_range_for_tier()
    total_spks = len(SPK_BODY_NAME_MAP)
    spk_results: Dict[str, str] = {}
    processed = 0

    logger.info(
        "Checking %d minor body SPKs (range %s to %s)...",
        total_spks,
        start_date,
        end_date,
    )

    from .constants import SPK_AUTO_DOWNLOAD_BLOCKED

    for ipl, (horizons_id, naif_id) in SPK_BODY_NAME_MAP.items():
        body_name = spk._get_body_name(ipl) or horizons_id
        processed += 1

        # Show progress
        if show_progress and sys.stdout.isatty():
            pct = processed * 100 // total_spks
            sys.stdout.write(
                f"\r  SPK files: [{processed}/{total_spks}] {pct}% - {body_name:<20s}"
            )
            sys.stdout.flush()

        # Skip bodies where JPL blocks SPK generation
        if ipl in SPK_AUTO_DOWNLOAD_BLOCKED:
            spk_results[body_name] = f"fallback ({SPK_AUTO_DOWNLOAD_BLOCKED[ipl]})"
            summary.setdefault("fallback", 0)
            summary["fallback"] += 1
            continue

        # Check if already registered
        if ipl in state._SPK_BODY_MAP and not force_download:
            spk_results[body_name] = "cached"
            summary["cached"] += 1
            continue

        # Try to download and register
        try:
            spk.download_and_register_spk(
                body=horizons_id,
                ipl=ipl,
                start=start_date,
                end=end_date,
                overwrite=force_download,
            )
            spk_results[body_name] = "downloaded"
            summary["downloaded"] += 1
        except (OSError, ValueError, KeyError, ImportError) as e:
            # Download failed - body will use Keplerian fallback
            spk_results[body_name] = f"fallback ({e})"
            summary["fallback"] += 1
            logger.debug(
                "SPK download failed for %s, will use Keplerian fallback: %s",
                body_name,
                e,
            )

    if show_progress and sys.stdout.isatty():
        sys.stdout.write(
            f"\r  SPK files: [{total_spks}/{total_spks}] 100% - Done{'':20s}\n"
        )
        sys.stdout.flush()

    results["spks"] = spk_results

    logger.info(
        "Ephemerides ready: %d cached, %d downloaded, %d fallback, %d errors",
        summary["cached"],
        summary["downloaded"],
        summary["fallback"],
        summary["errors"],
    )

    return results
