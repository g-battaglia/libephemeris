"""
Global state management for libephemeris.

This module maintains the library's singleton state including:
- Ephemeris data loader (Skyfield Loader)
- Planetary ephemeris (DE440 or DE441 via env var)
- Timescale (for UTC/TT conversions)
- Observer topocentric location
- Sidereal mode configuration
- Cached angles for Arabic parts calculations
- SPK kernel registry for minor body calculations

All state is stored in module-level globals to provide a stateful module-level API
compatible with pyswisseph's threading model (thread-unsafe by design).
"""

from __future__ import annotations

import os
import threading
from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, List, Literal, Optional, Tuple, Union, overload
from skyfield.api import Loader, Topos
from skyfield.timelib import Timescale
from skyfield.jpllib import SpiceKernel

from .logging_config import get_logger


# =============================================================================
# DATA DIRECTORY CONFIGURATION
# =============================================================================

DEFAULT_DATA_DIR = os.path.join(os.path.expanduser("~"), ".libephemeris")
_DATA_DIR_ENV_VAR = "LIBEPHEMERIS_DATA_DIR"


def _resolve_data_dir() -> str:
    """Resolve the data directory path without creating it.

    Resolution order:
        1. LIBEPHEMERIS_DATA_DIR environment variable
        2. TOML config ``data_dir``
        3. DEFAULT_DATA_DIR (~/.libephemeris)

    Returns:
        str: Absolute path to the data directory.
    """
    env_value = os.environ.get(_DATA_DIR_ENV_VAR, "").strip()
    if env_value:
        return os.path.abspath(env_value)
    from ._config_toml import get_str as _toml_str

    toml_value = _toml_str("data_dir")
    return os.path.abspath(toml_value if toml_value else DEFAULT_DATA_DIR)


def _get_data_dir() -> str:
    """Get the base data directory, creating it if needed.

    Returns:
        str: Absolute path to the data directory.
    """
    data_dir = _resolve_data_dir()
    os.makedirs(data_dir, exist_ok=True)
    return data_dir


# =============================================================================
# PRECISION TIER SYSTEM
# =============================================================================

# Three tiers with different date ranges and file sizes:
# - base: de440s.bsp (1849-2150, ~31 MB) - lightweight, modern usage
# - medium: de440.bsp (1549-2650, ~114 MB) - general purpose (DEFAULT)
# - extended: de441.bsp (-13198 to +17191, ~3.1 GB) - historical research


@dataclass(frozen=True)
class PrecisionTier:
    """Configuration for a precision tier."""

    name: str
    ephemeris_file: str
    spk_date_range: Tuple[str, str]
    description: str


TIERS: Dict[str, PrecisionTier] = {
    "base": PrecisionTier(
        name="base",
        ephemeris_file="de440s.bsp",
        spk_date_range=("1850-01-01", "2150-01-01"),
        description="Modern usage (1850-2150), ~31 MB",
    ),
    "medium": PrecisionTier(
        name="medium",
        ephemeris_file="de440.bsp",
        spk_date_range=("1900-01-01", "2100-01-01"),
        description="General purpose (1550-2650), ~114 MB",
    ),
    "extended": PrecisionTier(
        name="extended",
        ephemeris_file="de441.bsp",
        spk_date_range=("1600-01-01", "2500-01-01"),
        description="Extended range (-13200 to +17191), ~3.1 GB",
    ),
}

_DEFAULT_TIER: str = "medium"
_PRECISION_TIER_ENV_VAR: str = "LIBEPHEMERIS_PRECISION"

# =============================================================================
# GLOBAL STATE VARIABLES
# =============================================================================

# Lock to protect global state during context-swap operations.
# This ensures thread-safety when EphemerisContext temporarily swaps
# global state for calculations. Without this lock, concurrent threads
# could interfere with each other's state during the save-set-restore cycle.
# We use RLock (reentrant lock) to allow nested locking in case any
# calculation function internally calls other state-swapping functions.
_CONTEXT_SWAP_LOCK = threading.RLock()

# Lock to protect lazy initialization of global singletons (_LOADER, _TS, _PLANETS).
# Prevents race conditions when multiple threads call get_loader/get_timescale/get_planets
# concurrently before initialization has completed.
_INIT_LOCK = threading.RLock()

# Lock to protect mutable global state (setter/getter pairs).
# This is the SAME lock as _CONTEXT_SWAP_LOCK to prevent deadlocks between
# context-based calculations (which swap global state) and direct setters.
_STATE_LOCK = _CONTEXT_SWAP_LOCK

_EPHEMERIS_PATH: Optional[str] = None  # Custom ephemeris directory
_EPHEMERIS_FILE: str = "de440.bsp"  # Ephemeris file to use (default: DE440)
_EPHEMERIS_FILE_EXPLICIT: bool = False  # True if set_ephemeris_file() was called
_EPHEMERIS_ENV_VAR = "LIBEPHEMERIS_EPHEMERIS"  # Env var for ephemeris file selection
_PRECISION_TIER: Optional[str] = None  # Programmatic tier override
_LOADER: Optional[Loader] = None  # Skyfield data loader
_PLANETS: Optional[SpiceKernel] = None  # Loaded planetary ephemeris
_PLANET_CENTERS: Optional[SpiceKernel] = None  # Planet center offsets (599, 699, etc.)
_PLANET_CENTERS_TIER: Optional[str] = None  # Tier of loaded planet_centers file
_TS: Optional[Timescale] = None  # Timescale object
_TOPO: Optional[Topos] = None  # Observer location
_SIDEREAL_MODE: Optional[int] = None  # Active sidereal mode ID
_SIDEREAL_AYAN_T0: float = 0.0  # Ayanamsha value at reference epoch
_SIDEREAL_T0: float = 0.0  # Reference epoch (JD) for ayanamsha
_ANGLES_CACHE: dict[str, float] = {}  # Pre-calculated angles {name: longitude}
_TIDAL_ACCELERATION: Optional[float] = None  # Tidal acceleration for Delta T
_DELTA_T_USERDEF: Optional[float] = None  # User-defined Delta T value
_LAPSE_RATE: Optional[float] = None  # Atmospheric lapse rate for refraction

# SPK kernel state for minor body calculations
_SPK_KERNELS: dict[str, SpiceKernel] = {}  # Cached SPK kernels {filepath: kernel}
_SPK_TYPE21_KERNELS: dict[
    str, object
] = {}  # Cached SPK type 21 kernels {filepath: SPKType21}
_SPK_BODY_MAP: dict[
    int, tuple[str, int]
] = {}  # Body mappings {ipl: (spk_file, naif_id)}

# Auto SPK download configuration
# When enabled, attempts to download SPK from JPL Horizons for minor bodies
# before falling back to Keplerian propagation
_AUTO_SPK_DOWNLOAD: Optional[bool] = None  # None = check env var, True/False = explicit

# SPK cache directory override
# When set, SPK files are cached in this directory instead of the default
_SPK_CACHE_DIR: Optional[str] = None

# SPK date padding in days
# Extra time buffer added to date ranges when downloading SPK files
# e.g., if requesting 2020-2030, with padding=365 will download 2019-2031
_SPK_DATE_PADDING: int = 0

# IERS Delta T configuration
# When enabled, uses observed Delta T values from IERS for recent dates
_IERS_DELTA_T_ENABLED: Optional[bool] = None  # None = check env var

# LEB (binary ephemeris) configuration
_LEB_FILE: Optional[str] = None  # Path to .leb file
_LEB_READER: Optional["LEBReader"] = None  # Cached LEBReader instance

# Horizons API client
_HORIZONS_CLIENT: Optional[object] = None
_HORIZONS_WARNED: bool = False

# Calculation mode: "auto" (default), "skyfield", "leb", or "horizons"
_CALC_MODE: Optional[str] = None  # None = check env var
_CALC_MODE_ENV_VAR = "LIBEPHEMERIS_MODE"
_VALID_CALC_MODES = ("auto", "skyfield", "leb", "horizons")

if TYPE_CHECKING:
    from .leb_reader import LEBReader


def set_calc_mode(mode: Optional[str]) -> None:
    """Set the calculation mode for the library.

    Controls how calc_ut() and calc() resolve positions:

    - ``"auto"`` (default): Try LEB first, then Horizons API (if no
      local DE440), then Skyfield. This is the standard behavior.
    - ``"skyfield"``: Always use Skyfield/DE440. Fails if DE440 is not
      available locally.
    - ``"leb"``: Require a valid LEB file (configured, auto-discovered,
      or auto-downloaded). Raises RuntimeError if none can be resolved.
      Bodies not in the .leb file fall through to Skyfield.
    - ``"horizons"``: Always use NASA JPL Horizons API (requires internet).
      Bodies not supported by Horizons fall through to Skyfield.

    Args:
        mode: One of ``"auto"``, ``"skyfield"``, ``"leb"``, ``"horizons"``,
              or None to reset to environment variable / default.

    Raises:
        ValueError: If mode is not a valid mode string.

    Environment Variable:
        LIBEPHEMERIS_MODE: Same values as the ``mode`` argument.
        Case-insensitive. Default is ``"auto"``.

    Example:
        >>> from libephemeris import set_calc_mode, get_calc_mode
        >>> set_calc_mode("skyfield")  # Force Skyfield path
        >>> get_calc_mode()
        'skyfield'
        >>> set_calc_mode(None)  # Reset to env var / default
    """
    global _CALC_MODE
    if mode is not None:
        mode = mode.lower().strip()
        if mode not in _VALID_CALC_MODES:
            raise ValueError(
                f"Invalid mode: {mode!r}. Must be one of: {list(_VALID_CALC_MODES)}"
            )
    with _STATE_LOCK:
        _CALC_MODE = mode


def get_calc_mode() -> str:
    """Get the current calculation mode.

    Resolution order:
        1. Programmatic override via set_calc_mode()
        2. LIBEPHEMERIS_MODE environment variable
        3. TOML config ``mode`` key
        4. Default: ``"auto"``

    Returns:
        str: The active calculation mode (``"auto"``, ``"leb"``,
             ``"skyfield"``, or ``"horizons"``).

    Example:
        >>> from libephemeris import get_calc_mode
        >>> get_calc_mode()  # Default
        'auto'
    """
    with _STATE_LOCK:
        if _CALC_MODE is not None:
            return _CALC_MODE
    env_value = os.environ.get(_CALC_MODE_ENV_VAR, "").lower().strip()
    if env_value in _VALID_CALC_MODES:
        return env_value
    # TOML config fallback
    from ._config_toml import get_str as _toml_str

    toml_value = _toml_str("mode")
    if toml_value and toml_value.lower().strip() in _VALID_CALC_MODES:
        return toml_value.lower().strip()
    return "auto"


def set_leb_file(filepath: Optional[str]) -> None:
    """Set the .leb file path for binary ephemeris mode.

    When a .leb file is configured, calc_ut() and calc() will
    attempt to use precomputed Chebyshev polynomials for fast evaluation
    before falling back to the Skyfield pipeline.

    Args:
        filepath: Path to a .leb file, or None to disable binary mode.

    Example:
        >>> from libephemeris import set_leb_file, calc_ut, SUN
        >>> set_leb_file("/path/to/ephemeris.leb")
        >>> pos, _ = calc_ut(2451545.0, SUN, 0)  # uses .leb fast path
        >>> set_leb_file(None)  # disable binary mode
    """
    global _LEB_FILE, _LEB_READER
    with _STATE_LOCK:
        if _LEB_READER is not None:
            try:
                _LEB_READER.close()
            except (AttributeError, OSError):
                get_logger().debug("Failed to close LEB reader: %s", _LEB_FILE)
        _LEB_FILE = filepath
        _LEB_READER = None  # force re-creation on next access

        # Same cleanup as _close_inner(): stale _active_reader and caches
        # would otherwise serve data from the old reader.
        from . import fast_calc as _fast_calc

        _fast_calc._reset_active_reader()
        _fast_calc._leb_frame_cache.clear()


def _discover_leb_file() -> Optional[str]:
    """Auto-discover a LEB file for the active precision tier.

    Checks in order:
    1. LEB2 modular files in user data dir: ``~/.libephemeris/leb/{tier}_core.leb2``
    2. LEB2 legacy extension: ``~/.libephemeris/leb/{tier}_core.leb``
    3. LEB1 merged file in user data dir: ``~/.libephemeris/leb/ephemeris_{tier}.leb``
    4. Bundled LEB2 core file shipped with the package (base tier only)
    5. Auto-download LEB2 from GitHub Releases (in ``"auto"`` mode only)

    Returns:
        Path to the discovered LEB file, or None if not found.
    """
    tier = get_precision_tier()
    leb_dir = os.path.join(_get_data_dir(), "leb")

    # Check for LEB2 modular files first (core is the primary)
    leb2_core = os.path.join(leb_dir, f"{tier}_core.leb2")
    if os.path.isfile(leb2_core):
        return leb2_core

    # Legacy: LEB2 files with .leb extension (backward compatibility)
    leb2_core_legacy = os.path.join(leb_dir, f"{tier}_core.leb")
    if os.path.isfile(leb2_core_legacy):
        return leb2_core_legacy

    # Fall back to LEB1 merged file
    candidate = os.path.join(leb_dir, f"ephemeris_{tier}.leb")
    if os.path.isfile(candidate):
        return candidate

    # Fall back to bundled LEB2 core (ships with the PyPI wheel, base tier only)
    if tier == "base":
        bundled = os.path.join(
            os.path.dirname(__file__), "data", "leb2", "base_core.leb2"
        )
        if os.path.isfile(bundled):
            return bundled
        # Legacy bundled path
        bundled_legacy = os.path.join(
            os.path.dirname(__file__), "data", "leb2", "base_core.leb"
        )
        if os.path.isfile(bundled_legacy):
            return bundled_legacy

    # Auto-download LEB2 from GitHub Releases when no local LEB exists.
    # In "auto" mode: gives new users fast LEB2 calculations (~33 MB)
    # instead of falling through to Skyfield which downloads DE440 (~114 MB).
    # In "leb" mode: download what the user asked for instead of crashing.
    mode = get_calc_mode()
    if mode in ("auto", "leb"):
        logger = get_logger()
        logger.info(
            "No LEB file found for '%s' tier. Downloading LEB2 from GitHub Releases...",
            tier,
        )
        try:
            from .download import download_leb2_for_tier

            downloaded = download_leb2_for_tier(
                tier, quiet=False, show_progress=True, activate=False
            )
            if downloaded:
                # Re-check — the download should have placed the file
                leb2_core = os.path.join(leb_dir, f"{tier}_core.leb2")
                if os.path.isfile(leb2_core):
                    logger.info("LEB2 downloaded and ready: %s", leb2_core)
                    return leb2_core
        except Exception as e:
            logger.warning("LEB2 auto-download failed: %s", e)

    return None


def _year_to_jd(year: int) -> float:
    """Convert a Gregorian year to a Julian Day (Jan 1, 12:00 TT).

    Uses the proleptic Gregorian calendar algorithm.

    Args:
        year: Gregorian calendar year.

    Returns:
        Julian Day number for January 1 of the given year at noon TT.
    """
    a = (14 - 1) // 12
    y = year + 4800 - a
    m = 1 + 12 * a - 3
    return 1 + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045.0


def _maybe_warm_reader(reader: "LEBReader") -> None:
    """Conditionally warm the reader based on TOML ``mmap_preload`` config.

    When ``mmap_preload = true`` in the TOML configuration, converts the
    configured year range to Julian Days and calls ``reader.warm()``.

    Args:
        reader: An open LEBReader, LEB2Reader, or CompositeLEBReader.
    """
    from ._config_toml import get_bool, get_int

    if not get_bool("mmap_preload"):
        return

    start_year = get_int("mmap_preload_start") or 1800
    end_year = get_int("mmap_preload_end") or 2200
    jd_start = _year_to_jd(start_year)
    jd_end = _year_to_jd(end_year)
    try:
        reader.warm(jd_start, jd_end)
        logger = get_logger()
        logger.debug("mmap preload: warmed JD range %s-%s (%d-%d)", jd_start, jd_end, start_year, end_year)
    except (AttributeError, OSError, ValueError) as exc:
        logger = get_logger()
        logger.warning("mmap preload failed: %s", exc)


def release_data_cache() -> None:
    """Advise the kernel to reclaim page cache for ephemeris data files.

    Calls ``posix_fadvise(FADV_DONTNEED)`` on all regular files in the
    data directory.  In containerised environments (cgroup v2), this
    reduces the memory counter without affecting correctness — pages
    are reloaded transparently on next access.

    Note:
        Safe to call on any platform.  No-op where ``posix_fadvise``
        is unavailable (macOS, Windows).  Only processes regular files;
        FIFOs, devices, and symlinks are skipped.
    """
    _fadvise = getattr(os, "posix_fadvise", None)
    _advice = getattr(os, "POSIX_FADV_DONTNEED", None)
    if _fadvise is None or _advice is None:
        return
    data_dir = _resolve_data_dir()
    if not os.path.isdir(data_dir):
        return
    import stat

    for root, _, files in os.walk(data_dir):
        for name in files:
            path = os.path.join(root, name)
            try:
                st = os.stat(path)
                if not stat.S_ISREG(st.st_mode):
                    continue
            except OSError:
                continue
            fd = -1
            try:
                fd = os.open(path, os.O_RDONLY | os.O_NONBLOCK)
                _fadvise(fd, 0, st.st_size, _advice)
            except OSError:
                pass
            finally:
                if fd >= 0:
                    try:
                        os.close(fd)
                    except OSError:
                        pass


def get_leb_reader() -> Optional["LEBReader"]:
    """Get the active LEBReader, if any.

    Respects the calculation mode set via set_calc_mode() or the
    LIBEPHEMERIS_MODE environment variable:

    - ``"skyfield"``: Always returns None (LEB disabled).
    - ``"leb"``: Returns LEBReader or raises RuntimeError if unavailable.
    - ``"auto"`` (default): Returns LEBReader if configured or
      auto-discovered, else None.

    Resolution order for the .leb file path:

    1. Explicit path via ``set_leb_file()``
    2. ``LIBEPHEMERIS_LEB`` environment variable
    3. Auto-discovery: ``~/.libephemeris/leb/ephemeris_{tier}.leb``

    If the .leb file path is invalid or the file is corrupted,
    logs a warning and returns None (silent fallback to Skyfield),
    unless mode is ``"leb"`` in which case RuntimeError is raised.

    Returns:
        LEBReader instance if a .leb file is configured and valid,
        None otherwise.

    Raises:
        RuntimeError: If mode is ``"leb"`` and no valid .leb file is
                      available.
    """
    global _LEB_READER
    mode = get_calc_mode()

    # In skyfield mode, never use LEB
    if mode == "skyfield":
        return None

    if _LEB_READER is None:
        path = _LEB_FILE or os.environ.get("LIBEPHEMERIS_LEB")

        # TOML config fallback
        if path is None:
            from ._config_toml import get_str as _toml_str

            path = _toml_str("leb_file")

        # Auto-discover if no explicit path configured
        if path is None:
            path = _discover_leb_file()
            if path is not None:
                logger = get_logger()
                logger.debug("Auto-discovered LEB file: %s", path)

        if path is not None:
            try:
                basename = os.path.splitext(os.path.basename(path))[0]
                if "_" in basename:
                    from .leb_composite import CompositeLEBReader

                    _LEB_READER = CompositeLEBReader.from_file_with_companions(path)
                else:
                    from .leb_reader import open_leb

                    _LEB_READER = open_leb(path)
                _maybe_warm_reader(_LEB_READER)
            except (FileNotFoundError, ValueError, OSError) as e:
                if mode == "leb":
                    raise RuntimeError(
                        f"LIBEPHEMERIS_MODE=leb but failed to open LEB file {path}: {e}"
                    ) from e
                logger = get_logger()
                logger.warning("Failed to open LEB file %s: %s", path, e)
                return None
        elif mode == "leb":
            raise RuntimeError(
                "LIBEPHEMERIS_MODE=leb but no .leb file configured. "
                "Use set_leb_file() or set LIBEPHEMERIS_LEB env var."
            )
    return _LEB_READER


def get_horizons_client():
    """Get the Horizons API client, if the current mode allows it.

    Returns a HorizonsClient when:
    - mode="horizons" — always
    - mode="auto" — only when no DE440 is available locally

    Returns None when:
    - mode="skyfield" or mode="leb"
    """
    global _HORIZONS_CLIENT, _HORIZONS_WARNED
    mode = get_calc_mode()

    if mode in ("skyfield", "leb"):
        return None

    if mode == "horizons":
        return _get_or_create_horizons_client()

    # mode == "auto": use Horizons only when no local ephemeris
    if mode == "auto":
        # Check if DE440/441 is locally available
        data_dir = _get_data_dir()
        eph_file = _get_effective_ephemeris_file()
        local_paths = [
            os.path.join(data_dir, eph_file),
        ]
        if _EPHEMERIS_PATH:
            local_paths.append(os.path.join(_EPHEMERIS_PATH, eph_file))

        for p in local_paths:
            if os.path.isfile(p):
                return None  # local ephemeris exists, skip Horizons

        return _get_or_create_horizons_client()

    return None


def _get_or_create_horizons_client():
    """Create or return the singleton HorizonsClient."""
    global _HORIZONS_CLIENT, _HORIZONS_WARNED
    if _HORIZONS_CLIENT is None:
        from .horizons_backend import HorizonsClient

        _HORIZONS_CLIENT = HorizonsClient()
        if not _HORIZONS_WARNED:
            logger = get_logger()
            logger.info(
                "No local ephemeris file found. Using NASA JPL Horizons API "
                "(requires internet, slower). For offline/faster calculations:\n"
                "  python -m libephemeris download"
            )
            _HORIZONS_WARNED = True
    return _HORIZONS_CLIENT


def get_loader() -> Loader:
    """
    Get or create the Skyfield data loader.

    Returns:
        Loader: Skyfield Loader instance for downloading/caching ephemeris files

    Note:
        Data files are cached in ~/.libephemeris by default (or
        LIBEPHEMERIS_DATA_DIR environment variable if set).
    """
    global _LOADER
    if _LOADER is None:
        with _INIT_LOCK:
            if _LOADER is None:
                _LOADER = Loader(_get_data_dir())
    return _LOADER


def _build_enhanced_timescale() -> Timescale:
    """Build a Timescale with extended historical Delta T coverage.

    Skyfield's built-in ``iers.npz`` provides daily observed Delta T from
    ~1973 to ~2027.  For earlier dates it falls back to the Stephenson,
    Morrison & Hohenkerk (2016) cubic spline which is designed for
    millennial-scale modelling and can deviate from 20th-century
    observations by up to ~0.7 s (worst case around 1955).

    Skyfield also ships ``historic_deltat.npy`` — semi-annual observed
    values from 1657 to 1984 compiled from the *Astronomical Almanac*
    (McCarthy & Babcock 1986) — but the current ``build_delta_t()``
    function does not incorporate them.

    This function merges the two datasets into a single table that is
    then passed to ``build_delta_t()``, giving continuous observed
    coverage from **1657 to ~2027** with a smooth hand-off at the
    junction (~1973).
    """
    import numpy as np
    from pathlib import Path
    from skyfield.timelib import build_delta_t

    # Locate Skyfield's bundled data directory
    import skyfield.data

    data_dir = Path(skyfield.data.__file__).parent

    # --- Load IERS daily table (1973-~2027) ----------------------------------
    iers_path = data_dir / "iers.npz"
    iers = dict(np.load(str(iers_path)))

    tt_offset = iers["tt_jd_minus_arange"]
    n = len(tt_offset)
    iers_tt = tt_offset + np.arange(n, dtype=np.float64)
    iers_dt = iers["delta_t_1e7"] / 1e7

    leap_dates = iers["leap_dates"]
    leap_offsets = iers["leap_offsets"]

    # --- Load historical observed table (1657-1984) --------------------------
    hist_path = data_dir / "historic_deltat.npy"
    if not hist_path.exists():
        # Graceful fallback: if the file is missing (unlikely but
        # possible with a stripped Skyfield install), just use the
        # standard IERS-only timescale.
        load = get_loader()
        return load.timescale()

    hist = np.load(str(hist_path))  # shape (2, 656)

    # --- Merge: keep historic entries strictly before IERS start --------------
    cutoff_jd = iers_tt[0]
    mask = hist[0] < cutoff_jd
    merged_tt = np.concatenate([hist[0][mask], iers_tt])
    merged_dt = np.concatenate([hist[1][mask], iers_dt])

    # --- Build the composite Delta T function and Timescale ------------------
    delta_t_func = build_delta_t((merged_tt, merged_dt))
    return Timescale(delta_t_func, leap_dates, leap_offsets)


def get_timescale() -> Timescale:
    """
    Get or create the Skyfield timescale object.

    Returns:
        Timescale: Skyfield timescale for time conversions (UTC, TT, etc.)

    Note:
        Uses an enhanced timescale that merges Skyfield's bundled
        historical Delta T observations (1657-1984, from the
        *Astronomical Almanac*) with the modern IERS daily table
        (1973-~2027), providing continuous observed Delta T coverage
        from 1657 onwards.  For dates outside this range, Skyfield's
        standard Stephenson, Morrison & Hohenkerk (2016) long-term
        model is used automatically.
    """
    global _TS
    if _TS is None:
        with _INIT_LOCK:
            if _TS is None:
                try:
                    _TS = _build_enhanced_timescale()
                except (ImportError, FileNotFoundError, RuntimeError, ValueError):
                    # If anything goes wrong with the merge, fall back to the
                    # standard Skyfield timescale so the library never fails to
                    # initialise.
                    get_logger().warning(
                        "Failed to build enhanced timescale with historical "
                        "Delta T data; falling back to standard Skyfield timescale."
                    )
                    load = get_loader()
                    _TS = load.timescale()
    return _TS


# =============================================================================
# PRECISION TIER ACCESSORS
# =============================================================================


def _get_current_tier() -> PrecisionTier:
    """
    Get the current precision tier configuration.

    Priority:
        1. Programmatic override via set_precision_tier()
        2. LIBEPHEMERIS_PRECISION environment variable
        3. Default: "medium" (de440.bsp)

    Returns:
        PrecisionTier: The active precision tier configuration.
    """
    # Priority 1: programmatic override
    if _PRECISION_TIER is not None and _PRECISION_TIER in TIERS:
        return TIERS[_PRECISION_TIER]

    # Priority 2: env var
    env_value = os.environ.get(_PRECISION_TIER_ENV_VAR, "").lower().strip()
    if env_value in TIERS:
        return TIERS[env_value]

    # Priority 3: TOML config
    from ._config_toml import get_str as _toml_str

    toml_value = _toml_str("precision")
    if toml_value and toml_value.lower().strip() in TIERS:
        return TIERS[toml_value.lower().strip()]

    # Priority 4: default
    return TIERS[_DEFAULT_TIER]


def get_precision_tier() -> str:
    """
    Get the name of the current precision tier.

    Returns:
        str: Current tier name ("base", "medium", or "extended").

    Example:
        >>> from libephemeris.state import get_precision_tier
        >>> get_precision_tier()  # Default
        'medium'
    """
    return _get_current_tier().name


def set_precision_tier(tier: str) -> None:
    """
    Set the precision tier programmatically.

    This controls both the main ephemeris file and the date range used
    for auto-downloading SPK files for minor bodies.

    Args:
        tier: Tier name ("base", "medium", or "extended").

    Raises:
        ValueError: If tier name is not valid.

    Note:
        This overrides the LIBEPHEMERIS_PRECISION environment variable.
        Changing the tier clears cached ephemeris data to force a reload.

    Example:
        >>> from libephemeris.state import set_precision_tier, get_precision_tier
        >>> set_precision_tier("extended")
        >>> get_precision_tier()
        'extended'
    """
    global _PRECISION_TIER, _PLANETS, _EPHEMERIS_FILE_EXPLICIT
    if tier not in TIERS:
        raise ValueError(
            f"Invalid tier: {tier!r}. Must be one of: {list(TIERS.keys())}"
        )
    _PRECISION_TIER = tier
    _EPHEMERIS_FILE_EXPLICIT = False
    # Close old kernel and clear caches before reloading
    if _PLANETS is not None:
        try:
            _PLANETS.close()
        except (AttributeError, Exception):
            pass
    _PLANETS = None
    from .cache import clear_caches

    clear_caches()


def list_tiers() -> List[PrecisionTier]:
    """
    List all available precision tiers.

    Returns:
        List of PrecisionTier objects with name, ephemeris_file, etc.
    """
    return list(TIERS.values())


def get_spk_date_range_for_tier(tier_name: Optional[str] = None) -> Tuple[str, str]:
    """
    Get the SPK date range for a tier.

    This determines the date range used when auto-downloading SPK files
    from JPL Horizons for minor bodies.

    Args:
        tier_name: Tier name, or None to use the current tier.

    Returns:
        Tuple of (start_date, end_date) in "YYYY-MM-DD" format.

    Example:
        >>> from libephemeris.state import get_spk_date_range_for_tier
        >>> get_spk_date_range_for_tier()
        ('1900-01-01', '2100-01-01')
        >>> get_spk_date_range_for_tier("extended")
        ('1550-01-01', '2650-01-01')
    """
    if tier_name is not None:
        if tier_name not in TIERS:
            raise ValueError(
                f"Invalid tier: {tier_name!r}. Must be one of: {list(TIERS.keys())}"
            )
        return TIERS[tier_name].spk_date_range
    return _get_current_tier().spk_date_range


def _get_effective_ephemeris_file() -> str:
    """
    Determine the ephemeris file to use.

    Priority:
        1. LIBEPHEMERIS_EPHEMERIS environment variable (if set)
        2. set_ephemeris_file() programmatic override (explicit call)
        3. Precision tier setting (base/medium/extended)
        4. Default: "de440.bsp" (medium tier)

    Returns:
        str: The ephemeris filename to use (e.g., "de440.bsp", "de441.bsp").
    """
    # Priority 1: env var override
    env_value = os.environ.get(_EPHEMERIS_ENV_VAR, "").strip()
    if env_value:
        return env_value

    # Priority 2: explicit programmatic override via set_ephemeris_file()
    if _EPHEMERIS_FILE_EXPLICIT:
        return _EPHEMERIS_FILE

    # Priority 3: TOML config
    from ._config_toml import get_str as _toml_str

    toml_value = _toml_str("ephemeris")
    if toml_value:
        return toml_value

    # Priority 4: precision tier
    tier = _get_current_tier()
    return tier.ephemeris_file


def get_planets() -> SpiceKernel:
    """
    Get or load the planetary ephemeris (DE440 by default).

    Returns:
        SpiceKernel: Loaded JPL ephemeris kernel containing planetary positions

    Raises:
        FileNotFoundError: If ephemeris file cannot be found or downloaded

    Note:
        Uses the ephemeris file determined by (in priority order):
        1. set_ephemeris_file() if called explicitly
        2. LIBEPHEMERIS_EPHEMERIS environment variable (e.g., "de441.bsp")
        3. Default: "de440.bsp"

        Searches in _EPHEMERIS_PATH if set, then ~/.libephemeris (downloads if missing).
    """
    global _PLANETS, _EPHEMERIS_FILE
    logger = get_logger()
    if _PLANETS is None:
        with _INIT_LOCK:
            if _PLANETS is None:
                load = get_loader()
                effective_file = _get_effective_ephemeris_file()
                _EPHEMERIS_FILE = effective_file

                # Try custom ephemeris path first if set
                if _EPHEMERIS_PATH:
                    bsp_path = os.path.join(_EPHEMERIS_PATH, _EPHEMERIS_FILE)
                    if os.path.exists(bsp_path):
                        logger.debug("Using cached ephemeris: %s", bsp_path)
                        _PLANETS = load(bsp_path)
                        logger.info("Ephemeris loaded: %s", bsp_path)
                        return _PLANETS

                # Load from data dir (downloads automatically if missing)
                data_dir = _get_data_dir()
                bsp_path = os.path.join(data_dir, _EPHEMERIS_FILE)
                if os.path.exists(bsp_path):
                    logger.debug("Using cached ephemeris: %s", bsp_path)
                    _PLANETS = load(bsp_path)
                    logger.info("Ephemeris loaded: %s", bsp_path)
                elif get_calc_mode() == "leb":
                    # In LEB mode, the DE kernel is only needed for
                    # internal Skyfield calls (True Node, houses, etc.).
                    # Try lighter alternatives before downloading the
                    # full tier kernel (e.g., de441.bsp at 3.1 GB).
                    _fallback_files = ["de440s.bsp", "de440.bsp", "de441.bsp"]
                    loaded = False
                    for fb in _fallback_files:
                        fb_path = os.path.join(data_dir, fb)
                        if os.path.exists(fb_path):
                            logger.info(
                                "LEB mode: using available %s instead of %s",
                                fb,
                                _EPHEMERIS_FILE,
                            )
                            _PLANETS = load(fb_path)
                            loaded = True
                            break
                    if not loaded:
                        # No DE kernel at all — download the smallest one
                        logger.info(
                            "LEB mode: downloading de440s.bsp (31 MB) for "
                            "internal calculations..."
                        )
                        _PLANETS = load("de440s.bsp")
                else:
                    logger.info("Downloading JPL ephemeris %s...", _EPHEMERIS_FILE)
                    _PLANETS = load(_EPHEMERIS_FILE)
                    logger.info("Ephemeris downloaded to: %s", data_dir)
    return _PLANETS


def get_planet_centers() -> Optional[SpiceKernel]:
    """
    Get or load the planet centers SPK file for the active tier.

    This file contains precise offsets from system barycenters to planet centers
    for Jupiter (599), Saturn (699), Uranus (799), Neptune (899), and Pluto (999).

    The function loads the appropriate file for the current precision tier:
    - base: planet_centers_base.bsp (1850-2150)
    - medium: planet_centers_medium.bsp (1550-2650)
    - extended: planet_centers_extended.bsp (extended range, partial coverage)

    Falls back to legacy planet_centers.bsp if tier-specific file not found.

    Returns:
        SpiceKernel containing planet center segments, or None if not available.

    Note:
        The planet_centers files are generated using the
        scripts/generate_planet_centers_spk.py script and provide <0.001 arcsec
        precision for planet center positions.
    """
    global _PLANET_CENTERS, _PLANET_CENTERS_TIER

    current_tier = get_precision_tier()

    # Reload if tier changed or not loaded
    if _PLANET_CENTERS is None or _PLANET_CENTERS_TIER != current_tier:
        _PLANET_CENTERS = None
        _PLANET_CENTERS_TIER = None
        _PLANET_CENTER_MISSING.clear()

        data_dir = _get_data_dir()

        # Search order: tier-specific → legacy
        candidates = [
            os.path.join(data_dir, f"planet_centers_{current_tier}.bsp"),
            os.path.join(data_dir, "planet_centers.bsp"),
        ]

        load = get_loader()

        for path in candidates:
            if os.path.exists(path):
                try:
                    _PLANET_CENTERS = load(path)
                    _PLANET_CENTERS_TIER = current_tier
                    break
                except (FileNotFoundError, ValueError, KeyError, OSError) as e:
                    get_logger().warning(
                        "Failed to load planet_centers file %s: %s", path, e
                    )

    return _PLANET_CENTERS


# Negative cache: NAIF IDs confirmed absent from the planet_centers SPK file.
# Avoids redundant segment iteration when a requested body is not available.
# Cleared on close() / tier change when _PLANET_CENTERS is reloaded.
_PLANET_CENTER_MISSING: set[int] = set()


def get_planet_center_segment(naif_id: int, jd: Optional[float] = None):
    """
    Get a planet center segment from the planet_centers SPK file.

    This returns a Skyfield segment that can be used with .at(t) to get
    the position of the planet center relative to its system barycenter.

    Args:
        naif_id: NAIF ID of the planet center (599, 699, 799, 899, or 999)
        jd: Optional Julian Date to check coverage. If provided and the date
            is outside the segment's coverage, returns None (triggers fallback).

    Returns:
        Skyfield ChebyshevPosition segment, or None if not available or
        if jd is outside the segment's coverage.

    Example:
        >>> seg = get_planet_center_segment(599)  # Jupiter center
        >>> if seg:
        ...     offset = seg.at(t)  # Position relative to Jupiter barycenter
    """
    # Fast negative-cache check: skip segment iteration for IDs known to be absent
    if naif_id in _PLANET_CENTER_MISSING:
        return None

    centers = get_planet_centers()
    if centers is None:
        return None

    for seg in centers.segments:
        if seg.target == naif_id:
            # If jd provided, check coverage
            if jd is not None:
                try:
                    start_jd = getattr(seg, "start_jd", None)
                    end_jd = getattr(seg, "end_jd", None)
                    if start_jd is not None and end_jd is not None:
                        if not (start_jd <= jd <= end_jd):
                            return None  # Outside coverage, trigger fallback
                except (AttributeError, TypeError):
                    # If we can't check coverage, return segment anyway
                    get_logger().debug(
                        "Could not check SPK segment coverage for jd=%.1f", jd
                    )
            return seg

    # NAIF ID not found in any segment — remember for future calls
    _PLANET_CENTER_MISSING.add(naif_id)
    return None


def set_topo(lon: float, lat: float, alt: float = 0.0) -> None:
    """
    Set observer's topocentric location for planet calculations.

    Args:
        lon: Geographic longitude in degrees (East positive, range -180 to 180)
        lat: Geographic latitude in degrees (North positive, range -90 to 90)
        alt: Elevation above sea level in meters

    Raises:
        CoordinateError: If latitude is outside [-90, 90] or
                        longitude is outside [-180, 180]

    Note:
        Required for topocentric calculations (FLG_TOPOCTR),
        angles (Ascendant, MC), and Arabic parts.

    Example:
        >>> import libephemeris as ephem
        >>> ephem.set_topo(12.5, 41.9, 0)  # Rome
        >>> ephem.set_topo(-74.0, 40.7, 10)  # New York
    """
    from .exceptions import validate_coordinates

    validate_coordinates(lat, lon, "set_topo")
    global _TOPO
    with _STATE_LOCK:
        _TOPO = Topos(latitude_degrees=lat, longitude_degrees=lon, elevation_m=alt)
    # Invalidate the observer-at-time cache.  The cache keys include
    # ``id(observer)`` which is the memory address of the ``earth + topo``
    # VectorSum object.  When ``set_topo`` is called with a new location the
    # old VectorSum is deallocated and Python may reuse the same address for
    # the new one, causing cache hits that return positions computed for the
    # *previous* observer location.
    from .cache import clear_observer_cache

    clear_observer_cache()


def get_topo() -> Optional[Topos]:
    """
    Get the observer's topocentric location.

    Returns:
        Optional[Topos]: Current observer location or None if not set
    """
    with _STATE_LOCK:
        return _TOPO


def set_sid_mode(mode: int, t0: float = 0.0, ayan_t0: float = 0.0) -> None:
    """
    Set the sidereal mode (ayanamsha system) for calculations.

    Args:
        mode: Sidereal mode ID (SIDM_*) or 255 for custom
        t0: Reference epoch (Julian Day) for custom ayanamsha (default: J2000.0)
        ayan_t0: Ayanamsha value at t0 in degrees (for custom mode)

    Note:
        Affects all position calculations when FLG_SIDEREAL is set.
        Default is Lahiri (SIDM_LAHIRI) if never set.
    """
    global _SIDEREAL_MODE, _SIDEREAL_T0, _SIDEREAL_AYAN_T0
    with _STATE_LOCK:
        _SIDEREAL_MODE = mode
        _SIDEREAL_T0 = t0 if t0 != 0.0 else 2451545.0
        _SIDEREAL_AYAN_T0 = ayan_t0


@overload
def get_sid_mode(full: Literal[True]) -> tuple[int, float, float]: ...


@overload
def get_sid_mode(full: Literal[False] = ...) -> int: ...


def get_sid_mode(full: bool = False) -> Union[int, tuple[int, float, float]]:
    """
    Get the current sidereal mode configuration.

    Args:
        full: If True, return (mode, t0, ayan_t0); if False, return only mode ID

    Returns:
        int or tuple: Sidereal mode ID, or full configuration tuple

    Note:
        Returns SIDM_LAHIRI (1) by default if never set.
    """
    with _STATE_LOCK:
        mode = _SIDEREAL_MODE if _SIDEREAL_MODE is not None else 1
        if full:
            return mode, _SIDEREAL_T0, _SIDEREAL_AYAN_T0
        return mode


def get_angles_cache() -> dict[str, float]:
    """
    Get cached astrological angles for the current calculation context.

    Returns:
        dict: Cached angles {name: longitude_degrees}

    Note:
        Used by Arabic parts calculations which require pre-calculated
        planetary positions and angles.
    """
    return _ANGLES_CACHE


def set_angles_cache(angles: dict[str, float]) -> None:
    """
    Cache pre-calculated angles for use in Arabic parts and other virtual points.

    Args:
        angles: Dictionary of angles {name: longitude_degrees}
                e.g., {"Sun": 120.5, "Moon": 240.3, "Asc": 15.7}

    Note:
        Creates a copy to prevent external mutation of cache.
    """
    global _ANGLES_CACHE
    _ANGLES_CACHE = angles.copy()


def clear_angles_cache() -> None:
    """
    Clear the angles cache.

    Use this between unrelated calculation contexts to prevent stale data.
    """
    global _ANGLES_CACHE
    _ANGLES_CACHE = {}


def reset_session() -> None:
    """Reset per-calculation state without closing files or clearing caches.

    Resets: topo, sidereal mode, angles cache, observer cache.
    Keeps alive: LEB reader, Skyfield timescale, SPK kernels, LRU caches.

    Use this between independent calculations that may use different
    topo/sidereal settings. Use close() only for full teardown (e.g.
    switching ephemeris files, test cleanup, or application shutdown).
    """
    global _TOPO, _SIDEREAL_MODE, _SIDEREAL_AYAN_T0, _SIDEREAL_T0
    global _ANGLES_CACHE
    with _STATE_LOCK:
        _TOPO = None
        _SIDEREAL_MODE = None
        _SIDEREAL_AYAN_T0 = 0.0
        _SIDEREAL_T0 = 0.0
        _ANGLES_CACHE = {}
    from .cache import clear_observer_cache

    clear_observer_cache()


def set_ephe_path(path: Optional[str]) -> None:
    """
    Set the path for ephemeris files.

    Args:
        path: Path to directory containing ephemeris files.

    Note:
        This sets the directory where get_planets() will look for the ephemeris
        file specified by set_ephemeris_file(). If the file is not found there,
        it will fall back to the workspace root and then download if needed.

        When a valid directory is provided, any SPK files (.bsp) present in
        it are automatically scanned and registered for known minor bodies.
        This allows pre-downloaded SPK files to be used without manual
        registration. No network requests are made during this scan.

        Idempotent: if the path hasn't changed, no resources are released.
    """
    global _EPHEMERIS_PATH, _PLANETS
    if path == _EPHEMERIS_PATH:
        return  # No-op: path unchanged, keep loaded kernels and caches
    _EPHEMERIS_PATH = path
    # Close old kernel and clear caches before reloading
    if _PLANETS is not None:
        try:
            _PLANETS.close()
        except (AttributeError, Exception):
            pass
    _PLANETS = None
    from .cache import clear_caches

    clear_caches()

    # Auto-discover local SPK files (no network, just file scanning)
    if path and os.path.isdir(path):
        try:
            from .spk_auto import discover_local_spks

            discover_local_spks(path)
        except (ImportError, OSError, ValueError):
            # Discovery is best-effort; don't fail set_ephe_path()
            get_logger().debug("Local SPK discovery failed for path: %s", path)


def set_ephemeris_file(filename: str) -> None:
    """
    Set the ephemeris file to use for planetary calculations.

    Args:
        filename: Name of the JPL ephemeris file (e.g., "de440.bsp", "de441.bsp")

    Note:
        Supported ephemeris files:
        - de440s.bsp: 1849-2150 (31 MB) - lightweight subset of DE440
        - de440.bsp: 1550-2650 (default, 128 MB) - recommended for most uses
        - de441.bsp: -13200-17191 (3.4 GB) - for extended historical work

        The file will be searched in:
        1. The path set by set_ephe_path() if configured
        2. The workspace root directory
        3. Downloaded from JPL if not found locally

        Changing the ephemeris file clears the cached planets and forces a reload.

        This takes priority over the precision tier setting. To revert to
        tier-based ephemeris selection, call close() to reset state.
    """
    global _EPHEMERIS_FILE, _EPHEMERIS_FILE_EXPLICIT, _PLANETS
    _EPHEMERIS_FILE = filename
    _EPHEMERIS_FILE_EXPLICIT = True
    # Close old kernel and clear caches before reloading
    if _PLANETS is not None:
        try:
            _PLANETS.close()
        except (AttributeError, Exception):
            pass
    _PLANETS = None
    from .cache import clear_caches

    clear_caches()


def set_jpl_file(filename: str) -> None:
    """
    Set the JPL ephemeris file to use for planetary calculations.

    This is an alias for set_ephemeris_file() with a more descriptive name
    emphasizing JPL (Jet Propulsion Laboratory) ephemeris files.

    Args:
        filename: Name of the JPL ephemeris file (e.g., "de440.bsp", "de441.bsp")

    Note:
        Supported ephemeris files:
        - de440s.bsp: 1849-2150 (31 MB) - lightweight subset of DE440
        - de440.bsp: 1550-2650 (default, 128 MB) - recommended for most uses
        - de441.bsp: -13200-17191 (3.4 GB) - for extended historical work

        If the file is not found locally, Skyfield will attempt to download it.
        Use set_ephe_path() to specify a local directory containing the file.

        Alternatively, use set_precision_tier() to select a tier (base/medium/
        extended), which automatically selects the appropriate ephemeris file.
        Note: set_ephemeris_file() takes priority over the tier setting.

    Example:
        >>> from libephemeris import set_jpl_file, set_ephe_path
        >>> set_ephe_path("/path/to/ephemeris/files")
        >>> set_jpl_file("de441.bsp")  # Use local de441.bsp file
    """
    set_ephemeris_file(filename)


def set_tid_acc(acc: float) -> None:
    """
    Set the tidal acceleration used in Delta T calculations.

    The tidal acceleration of the Moon affects the long-term extrapolation
    of Delta T (TT - UT1), which is important for historical astronomical
    calculations. Different JPL ephemeris files assume different values.

    Args:
        acc: Tidal acceleration in arcsec/century^2.
             Use TIDAL_* constants for standard ephemeris values,
             or TIDAL_AUTOMATIC (999999) to use the default.

    Note:
        - The default value is based on DE440 (-25.936 arcsec/cy^2)
        - This setting primarily affects Delta T calculations for dates
          far from the present (historical studies, ancient astronomy)
        - For modern dates (1900-2100), the effect is negligible
        - Common values:
          - DE421: -25.85 arcsec/cy^2
          - DE431: -25.80 arcsec/cy^2
          - DE440/DE441: -25.936 arcsec/cy^2 (default)

    Example:
        >>> from libephemeris import set_tid_acc, get_tid_acc
        >>> from libephemeris import TIDAL_DE421, TIDAL_DE440
        >>> set_tid_acc(TIDAL_DE421)  # Use DE421 tidal acceleration
        >>> get_tid_acc()
        -25.85
        >>> set_tid_acc(TIDAL_DE440)  # Use DE440 tidal acceleration
        >>> get_tid_acc()
        -25.936
    """
    global _TIDAL_ACCELERATION
    from .constants import TIDAL_AUTOMATIC

    _TIDAL_ACCELERATION = acc if acc != TIDAL_AUTOMATIC else None


def get_tid_acc() -> float:
    """
    Get the current tidal acceleration used in Delta T calculations.

    Returns:
        float: The tidal acceleration in arcsec/century^2.
               Returns TIDAL_DEFAULT (-25.8) if not explicitly set.

    Note:
        The tidal acceleration affects how Delta T is extrapolated for dates
        outside the range of direct observations. This is particularly
        important for historical astronomical calculations.

        If set_tid_acc() was called with TIDAL_AUTOMATIC (999999),
        this returns the default value.

    Example:
        >>> from libephemeris import get_tid_acc, set_tid_acc, TIDAL_DE441
        >>> get_tid_acc()  # Default value
        -25.936
        >>> set_tid_acc(TIDAL_DE441)
        >>> get_tid_acc()
        -25.936
    """
    from .constants import TIDAL_DEFAULT

    if _TIDAL_ACCELERATION is None:
        return TIDAL_DEFAULT
    return _TIDAL_ACCELERATION


def set_delta_t_userdef(acc: Optional[float]) -> None:
    """
    Set a user-defined Delta T value to use instead of computed values.

    When set, deltat() and deltat_ex() will return this fixed value
    instead of computing Delta T from IERS data. This is useful for:
    - Testing and reproducibility
    - Very ancient dates where Delta T is highly uncertain
    - Very future dates where Delta T cannot be predicted
    - Experimentation with different Delta T assumptions

    Args:
        acc: Delta T value in days (TT - UT1), or None to clear and resume
             using computed values. The value should be in the same units
             as returned by deltat() (days, not seconds).

    Note:
        - To convert from seconds to days, divide by 86400
        - Modern Delta T (year 2000) is approximately 64 seconds (~0.00074 days)
        - For ancient dates, Delta T can be hours or even days
        - Use get_delta_t_userdef() to check if a user-defined value is active

    Example:
        >>> from libephemeris import set_delta_t_userdef, get_delta_t_userdef, deltat
        >>> # Set a fixed Delta T of 65 seconds (in days)
        >>> set_delta_t_userdef(65.0 / 86400.0)
        >>> get_delta_t_userdef()
        7.523148148148148e-04
        >>> deltat(2451545.0)  # Now returns fixed value
        7.523148148148148e-04
        >>> # Clear to resume computed values
        >>> set_delta_t_userdef(None)
        >>> get_delta_t_userdef() is None
        True
    """
    global _DELTA_T_USERDEF
    _DELTA_T_USERDEF = acc


def get_delta_t_userdef() -> Optional[float]:
    """
    Get the current user-defined Delta T value.

    Returns:
        Optional[float]: The user-defined Delta T in days if set,
                        or None if using computed values.

    Note:
        When this returns None, deltat() computes Delta T from IERS data.
        When this returns a float, that value is used directly by deltat().

    Example:
        >>> from libephemeris import set_delta_t_userdef, get_delta_t_userdef
        >>> get_delta_t_userdef()  # No user-defined value
        None
        >>> set_delta_t_userdef(0.001)  # Set ~86 seconds
        >>> get_delta_t_userdef()
        0.001
        >>> set_delta_t_userdef(None)  # Clear
        >>> get_delta_t_userdef()
        None
    """
    return _DELTA_T_USERDEF


# Default atmospheric lapse rate in K/m (standard atmosphere)
LAPSE_RATE_DEFAULT: float = 0.0065


def set_lapse_rate(lrate: Optional[float]) -> None:
    """
    Set the atmospheric temperature lapse rate for refraction calculations.

    The lapse rate is the rate at which temperature decreases with altitude
    in the atmosphere. It affects the calculation of atmospheric refraction,
    particularly the dip of the horizon for elevated observers.

    Args:
        lrate: Temperature lapse rate dT/dh in degrees Kelvin per meter.
               Use None to reset to the default value (0.0065 K/m).
               Typical values range from 0.0034 to 0.010 K/m.

    Note:
        - Standard atmospheric lapse rate is 0.0065 K/m (6.5°C per 1000m)
        - This setting affects refrac_extended() calculations
        - Higher lapse rates result in less refraction of the horizon
        - The lapse rate is used in calculating the dip of the horizon
          for elevated observers

    Example:
        >>> from libephemeris import set_lapse_rate, get_lapse_rate
        >>> get_lapse_rate()  # Default value
        0.0065
        >>> set_lapse_rate(0.005)  # Set custom lapse rate
        >>> get_lapse_rate()
        0.005
        >>> set_lapse_rate(None)  # Reset to default
        >>> get_lapse_rate()
        0.0065
    """
    global _LAPSE_RATE
    _LAPSE_RATE = lrate


def get_lapse_rate() -> float:
    """
    Get the current atmospheric temperature lapse rate.

    Returns:
        float: The lapse rate in degrees Kelvin per meter.
               Returns 0.0065 K/m (standard atmosphere) if not explicitly set.

    Note:
        The lapse rate affects refrac_extended() calculations, particularly
        the dip of the horizon for elevated observers.

    Example:
        >>> from libephemeris import get_lapse_rate, set_lapse_rate
        >>> get_lapse_rate()  # Default value
        0.0065
        >>> set_lapse_rate(0.008)
        >>> get_lapse_rate()
        0.008
    """
    if _LAPSE_RATE is None:
        return LAPSE_RATE_DEFAULT
    return _LAPSE_RATE


def get_library_path() -> str:
    """
    Get the path where ephemeris files are stored.

    For libephemeris, this returns the directory containing the JPL ephemeris
    files (.bsp files used by Skyfield). This is analogous to the reference API's
    get_library_path(), which returns the ephemeris data path.

    Returns:
        str: The absolute path to the directory containing ephemeris files.
             If set_ephe_path() was called, returns that custom path.
             Otherwise returns the default data directory (~/.libephemeris).

    Note:
        - If a custom ephemeris path was set via set_ephe_path(), that path
          is returned regardless of whether files actually exist there.
        - Otherwise, returns ~/.libephemeris (or LIBEPHEMERIS_DATA_DIR env var).

    Example:
        >>> from libephemeris import get_library_path, set_ephe_path
        >>> get_library_path()  # Returns default path
        '/home/user/.libephemeris'
        >>> set_ephe_path('/custom/ephemeris/path')
        >>> get_library_path()  # Returns custom path
        '/custom/ephemeris/path'
    """
    if _EPHEMERIS_PATH is not None:
        return os.path.abspath(_EPHEMERIS_PATH)
    return _get_data_dir()


def close() -> None:
    """
    Close all opened ephemeris files and release resources.

    This function closes the SPK kernel file handles and resets all global
    state to its initial values. Call this when you're done using the library
    or want to reload ephemeris files with different settings.

    Note:
        - This closes the underlying file handles in the JPL ephemeris kernel
        - After calling close(), the next ephemeris calculation will
          automatically reload the ephemeris files as needed
        - This is useful for:
          - Freeing memory and file handles in long-running applications
          - Switching to a different ephemeris file
          - Ensuring clean state in test suites

    Example:
        >>> from libephemeris import calc_ut, close, SUN
        >>> pos, _ = calc_ut(2451545.0, SUN, 0)  # Loads ephemeris
        >>> close()  # Close files and reset state
        >>> pos, _ = calc_ut(2451545.0, SUN, 0)  # Reloads ephemeris
    """
    global _EPHEMERIS_PATH, _EPHEMERIS_FILE, _EPHEMERIS_FILE_EXPLICIT
    global _LOADER, _PLANETS, _PLANET_CENTERS, _TS
    global _TOPO, _SIDEREAL_MODE, _SIDEREAL_AYAN_T0, _SIDEREAL_T0
    global _ANGLES_CACHE, _TIDAL_ACCELERATION, _DELTA_T_USERDEF, _LAPSE_RATE
    global _SPK_KERNELS, _SPK_BODY_MAP, _SPK_TYPE21_KERNELS, _AUTO_SPK_DOWNLOAD
    global _SPK_CACHE_DIR, _SPK_DATE_PADDING, _IERS_DELTA_T_ENABLED
    global _PRECISION_TIER
    global _LEB_FILE, _LEB_READER, _CALC_MODE
    global _HORIZONS_CLIENT, _HORIZONS_WARNED

    with _STATE_LOCK:
        _close_inner()


def _close_inner() -> None:
    """Inner close implementation (must be called with _STATE_LOCK held)."""
    global _EPHEMERIS_PATH, _EPHEMERIS_FILE, _EPHEMERIS_FILE_EXPLICIT
    global _LOADER, _PLANETS, _PLANET_CENTERS, _TS
    global _TOPO, _SIDEREAL_MODE, _SIDEREAL_AYAN_T0, _SIDEREAL_T0
    global _ANGLES_CACHE, _TIDAL_ACCELERATION, _DELTA_T_USERDEF, _LAPSE_RATE
    global _SPK_KERNELS, _SPK_BODY_MAP, _SPK_TYPE21_KERNELS, _AUTO_SPK_DOWNLOAD
    global _SPK_CACHE_DIR, _SPK_DATE_PADDING, _IERS_DELTA_T_ENABLED
    global _PRECISION_TIER
    global _LEB_FILE, _LEB_READER, _CALC_MODE
    global _HORIZONS_CLIENT, _HORIZONS_WARNED

    # Close the Horizons client if loaded
    if _HORIZONS_CLIENT is not None:
        try:
            _HORIZONS_CLIENT.shutdown()
        except (AttributeError, OSError):
            pass
    _HORIZONS_CLIENT = None
    _HORIZONS_WARNED = False

    # Close the LEB reader if loaded
    if _LEB_READER is not None:
        try:
            _LEB_READER.close()
        except (AttributeError, OSError):
            pass
    _LEB_FILE = None
    _LEB_READER = None
    _CALC_MODE = None

    # Clear fast_calc's cached reference to the closed LEB reader. Without
    # this, helpers like _frame_data() would still dispatch through a stale
    # _active_reader whose mmap is now None.
    from . import fast_calc as _fast_calc

    _fast_calc._reset_active_reader()

    # Close the SPK kernel file handles if loaded
    if _PLANETS is not None:
        try:
            _PLANETS.close()
        except (AttributeError, Exception):
            # SpiceKernel may not have close() in all versions,
            # or may already be closed
            pass

    # Close planet centers kernel if loaded
    if _PLANET_CENTERS is not None:
        try:
            _PLANET_CENTERS.close()
        except (AttributeError, Exception):
            pass

    # Clear computation caches for hot path optimization
    from .cache import clear_caches

    clear_caches()

    # Reset all global state to initial values
    _EPHEMERIS_PATH = None
    _EPHEMERIS_FILE = "de440.bsp"
    _EPHEMERIS_FILE_EXPLICIT = False
    _LOADER = None
    _PLANETS = None
    _PLANET_CENTERS = None
    _PLANET_CENTER_MISSING.clear()
    _TS = None
    _TOPO = None
    _SIDEREAL_MODE = None
    _SIDEREAL_AYAN_T0 = 0.0
    _SIDEREAL_T0 = 0.0
    _ANGLES_CACHE = {}
    _TIDAL_ACCELERATION = None
    _DELTA_T_USERDEF = None
    _LAPSE_RATE = None
    _AUTO_SPK_DOWNLOAD = None
    _SPK_CACHE_DIR = None
    _SPK_DATE_PADDING = 0
    _IERS_DELTA_T_ENABLED = None
    _PRECISION_TIER = None

    # Clear IERS cache
    try:
        from . import iers_data

        iers_data.clear_iers_cache()
    except ImportError:
        pass

    # Reset TOML config cache so it re-discovers on next access
    try:
        from ._config_toml import reset as _reset_toml_config

        _reset_toml_config()
    except ImportError:
        pass

    # Close and clear SPK kernels
    for kernel in _SPK_KERNELS.values():
        try:
            kernel.close()
        except (AttributeError, Exception):
            pass
    _SPK_KERNELS = {}
    _SPK_BODY_MAP = {}

    # Close and clear SPK type 21 kernels
    for kernel in _SPK_TYPE21_KERNELS.values():
        try:
            kernel.close()
        except (AttributeError, Exception):
            pass
    _SPK_TYPE21_KERNELS = {}

    # Close planetary moon kernels
    from . import planetary_moons

    planetary_moons.close_moon_kernels()

    # Reset ASSIST data availability cache
    try:
        from .rebound_integration import reset_assist_data_cache

        reset_assist_data_cache()
    except ImportError:
        pass


def get_current_file_data(ifno: int = 0) -> tuple[str, float, float, int]:
    """
    Get information about the currently loaded ephemeris file.

    This function returns details about the JPL ephemeris file currently in use,
    including the file path and the date range covered by the ephemeris.

    Args:
        ifno: File number indicator (for compatibility with reference API):
              - 0: planet file (main ephemeris - DE421, DE431, etc.)
              - 1: moon file (same as planets in JPL ephemeris)
              - 2: main asteroid file (not applicable for JPL ephemeris)
              - 3: other asteroid file (not applicable for JPL ephemeris)
              - 4: star file (not applicable for JPL ephemeris)

    Returns:
        A tuple containing:
            - path (str): Full path to the ephemeris file, or empty string if
                         no ephemeris is loaded or ifno is not applicable
            - start (float): Start date of the ephemeris file (Julian Day)
            - end (float): End date of the ephemeris file (Julian Day)
            - denum (int): JPL Development Ephemeris number (e.g., 421, 431, 441)
                          or 0 if not determinable

    Note:
        - For libephemeris using JPL/Skyfield ephemerides, ifno values 0 and 1
          return the same data since planets and Moon are in the same file.
        - ifno values 2, 3, 4 return empty data as those file types are not
          used by libephemeris.
        - The ephemeris file must have been loaded (by calling calc_ut or
          get_planets) before this function returns meaningful data.

    Example:
        >>> from libephemeris import get_current_file_data, calc_ut, SUN
        >>> calc_ut(2451545.0, SUN, 0)  # Triggers ephemeris loading
        >>> path, start, end, denum = get_current_file_data(0)
        >>> print(f"Using DE{denum}, covering JD {start:.1f} to {end:.1f}")
        Using DE421, covering JD 2414864.5 to 2471184.5
    """
    # Only file types 0 (planets) and 1 (moon) are meaningful for JPL ephemeris
    # Types 2, 3, 4 (asteroids, stars) are not applicable
    if ifno not in (0, 1):
        return ("", 0.0, 0.0, 0)

    # If no ephemeris loaded, return empty data
    if _PLANETS is None:
        return ("", 0.0, 0.0, 0)

    try:
        # Get the file path
        path = str(_PLANETS.path) if hasattr(_PLANETS, "path") else ""
        if not path and hasattr(_PLANETS, "filename"):
            path = str(_PLANETS.filename)

        # Get start/end dates from the SPK segments
        # All segments in a JPL DE file typically have the same date range
        start_jd = 0.0
        end_jd = 0.0

        if hasattr(_PLANETS, "spk") and hasattr(_PLANETS.spk, "segments"):
            segments = list(_PLANETS.spk.segments)
            if segments:
                # DE441 has segments split at 1969, so we need to find the
                # overall min start_jd and max end_jd across ALL segments
                start_jd = float("inf")
                end_jd = float("-inf")
                for seg in segments:
                    if hasattr(seg, "start_jd"):
                        seg_start = float(seg.start_jd)
                        if seg_start < start_jd:
                            start_jd = seg_start
                    if hasattr(seg, "end_jd"):
                        seg_end = float(seg.end_jd)
                        if seg_end > end_jd:
                            end_jd = seg_end
                # Fallback if no valid dates found
                if start_jd == float("inf"):
                    start_jd = 0.0
                if end_jd == float("-inf"):
                    end_jd = 0.0

        # Extract DE number from filename (e.g., "de421.bsp" -> 421)
        denum = 0
        filename = _EPHEMERIS_FILE.lower()
        if filename.startswith("de") and ".bsp" in filename:
            try:
                # Extract digits between "de" and ".bsp"
                num_str = filename[2:].split(".")[0]
                denum = int(num_str)
            except (ValueError, IndexError):
                denum = 0

        return (path, start_jd, end_jd, denum)

    except (ValueError, IndexError, KeyError, OSError):
        # Return empty data on any error
        return ("", 0.0, 0.0, 0)


# =============================================================================
# AUTO SPK DOWNLOAD CONFIGURATION
# =============================================================================

# Environment variable name for auto SPK download
_AUTO_SPK_ENV_VAR = "LIBEPHEMERIS_AUTO_SPK"


def set_auto_spk_download(enabled: Optional[bool]) -> None:
    """
    Enable or disable automatic SPK download for minor bodies.

    When enabled, the library will attempt to download high-precision SPK files
    from JPL Horizons for minor bodies before falling back to Keplerian propagation.
    This provides significantly better accuracy (~arcseconds vs ~arcminutes).

    Args:
        enabled: True to enable automatic SPK download, False to disable,
                 or None to use the environment variable (LIBEPHEMERIS_AUTO_SPK).

    Note:
        - Requires the 'astroquery' package to be installed for downloads
        - Downloads are cached in ~/.libephemeris/spk/ directory
        - If astroquery is not available and auto-download is enabled,
          the library will silently fall back to Keplerian propagation
        - Set to False for offline use or to ensure consistent behavior
        - Default is True (enabled) when no env var is set

    Environment Variable:
        LIBEPHEMERIS_AUTO_SPK: Set to "0", "false", or "no" to disable,
                               "1", "true", or "yes" to enable.
                               Case-insensitive. Default is enabled.

    Example:
        >>> from libephemeris import set_auto_spk_download, get_auto_spk_download
        >>> set_auto_spk_download(True)  # Enable automatic SPK downloads
        >>> get_auto_spk_download()
        True
        >>> set_auto_spk_download(None)  # Use environment variable
    """
    global _AUTO_SPK_DOWNLOAD
    _AUTO_SPK_DOWNLOAD = enabled


def get_auto_spk_download() -> bool:
    """
    Get the current auto SPK download setting.

    Returns:
        bool: True if automatic SPK download is enabled, False otherwise.

    Note:
        - If set_auto_spk_download() was called with an explicit value, returns that.
        - Otherwise, checks the LIBEPHEMERIS_AUTO_SPK environment variable.
        - If the environment variable is not set, returns True (enabled by default).
        - Set LIBEPHEMERIS_AUTO_SPK=0 to disable auto-download via env var.

    Example:
        >>> from libephemeris import get_auto_spk_download
        >>> get_auto_spk_download()  # Default is True
        True
        >>> import os
        >>> os.environ['LIBEPHEMERIS_AUTO_SPK'] = '0'
        >>> get_auto_spk_download()  # Now disabled via env var
        False
    """
    # If explicitly set via function, use that value
    if _AUTO_SPK_DOWNLOAD is not None:
        return _AUTO_SPK_DOWNLOAD

    # Otherwise check environment variable
    env_value = os.environ.get(_AUTO_SPK_ENV_VAR, "").lower().strip()

    # If env var explicitly disables, return False
    if env_value in ("0", "false", "no", "off", "disabled"):
        return False
    # If env var explicitly enables, return True
    if env_value in ("1", "true", "yes", "on", "enabled"):
        return True

    # TOML config fallback
    from ._config_toml import get_bool as _toml_bool

    toml_value = _toml_bool("auto_spk")
    if toml_value is not None:
        return toml_value

    # Default to True (auto-download enabled) because:
    # - strict_precision defaults to True, so SPK kernels are required
    # - Without auto-download, minor bodies like Chiron raise SPKRequiredError
    # - The library already requires network for DE440 via Skyfield
    # - Downloads are cached in ~/.libephemeris/spk/ and only happen once
    # Users can explicitly disable via env var or set_auto_spk_download(False)
    return True


def set_spk_cache_dir(path: Optional[str]) -> None:
    """
    Set the directory for caching SPK files.

    When set, SPK files will be stored in this directory instead of the
    default cache location (~/.libephemeris/spk/).

    Args:
        path: Directory path for SPK cache, or None to use the default.
              The directory will be created if it doesn't exist.

    Note:
        - This affects where auto_get_spk() and enable_auto_spk() store files
        - Already-registered SPK bodies are not affected
        - Set to None to revert to the default cache location

    Example:
        >>> from libephemeris import set_spk_cache_dir, get_spk_cache_dir
        >>> set_spk_cache_dir("/data/ephemeris/spk_cache")
        >>> get_spk_cache_dir()
        '/data/ephemeris/spk_cache'
        >>> set_spk_cache_dir(None)  # Revert to default
        >>> get_spk_cache_dir() is None
        True
    """
    global _SPK_CACHE_DIR
    _SPK_CACHE_DIR = path


_SPK_CACHE_DIR_ENV_VAR = "LIBEPHEMERIS_SPK_DIR"


def get_spk_cache_dir() -> Optional[str]:
    """
    Get the effective SPK cache directory.

    Priority:
        1. Value set via set_spk_cache_dir() (programmatic override)
        2. LIBEPHEMERIS_SPK_DIR environment variable
        3. None (caller uses DEFAULT_AUTO_SPK_DIR as fallback)

    Returns:
        Optional[str]: The effective cache directory path if set via
                      programmatic override or environment variable,
                      or None if using the default location.

    Note:
        When this returns None, SPK files are cached in the default
        location (~/.libephemeris/spk/).

    Example:
        >>> from libephemeris import get_spk_cache_dir, set_spk_cache_dir
        >>> get_spk_cache_dir()  # Default is None (uses default location)
        >>> set_spk_cache_dir("/custom/cache")
        >>> get_spk_cache_dir()
        '/custom/cache'
        >>> set_spk_cache_dir(None)  # Revert to default
        >>> import os; os.environ["LIBEPHEMERIS_SPK_DIR"] = "/env/cache"
        >>> get_spk_cache_dir()
        '/env/cache'
    """
    if _SPK_CACHE_DIR is not None:
        return _SPK_CACHE_DIR
    env_dir = os.environ.get(_SPK_CACHE_DIR_ENV_VAR)
    if env_dir:
        return env_dir
    # TOML config fallback
    from ._config_toml import get_str as _toml_str

    toml_value = _toml_str("spk_dir")
    if toml_value:
        return toml_value
    return None


def set_spk_date_padding(days: int) -> None:
    """
    Set the date padding for SPK downloads.

    When downloading SPK files, this number of days is added before and
    after the requested date range. This provides a buffer to ensure
    coverage even when calculations occur near the edges of the range.

    Args:
        days: Number of days to add as padding on each side.
              Must be non-negative (0 means no padding).

    Note:
        - A padding of 365 days means downloading 1 year extra on each side
        - This is useful when you need to calculate positions for dates
          slightly outside your main date range
        - Default is 0 (no padding)

    Example:
        >>> from libephemeris import set_spk_date_padding, get_spk_date_padding
        >>> set_spk_date_padding(365)  # Add 1 year buffer on each side
        >>> get_spk_date_padding()
        365
        >>> set_spk_date_padding(0)  # No padding
    """
    global _SPK_DATE_PADDING
    if days < 0:
        raise ValueError("spk_date_padding must be non-negative")
    _SPK_DATE_PADDING = days


def get_spk_date_padding() -> int:
    """
    Get the current SPK date padding in days.

    Returns:
        int: Number of days added as padding on each side when downloading
             SPK files. Default is 0.

    Note:
        This value is used by auto_get_spk() to expand the requested
        date range when downloading new SPK files.

    Example:
        >>> from libephemeris import get_spk_date_padding, set_spk_date_padding
        >>> get_spk_date_padding()  # Default
        0
        >>> set_spk_date_padding(30)
        >>> get_spk_date_padding()
        30
    """
    return _SPK_DATE_PADDING


# =============================================================================
# IERS DELTA T CONFIGURATION
# =============================================================================

# Environment variable name for IERS Delta T
_IERS_DELTA_T_ENV_VAR = "LIBEPHEMERIS_IERS_DELTA_T"


def set_iers_delta_t_enabled(enabled: Optional[bool]) -> None:
    """
    Enable or disable using IERS observed Delta T values for recent dates.

    When enabled, the library will use high-precision observed Delta T values
    from IERS (International Earth Rotation and Reference Systems Service)
    for dates where IERS data is available (typically 1973-present).

    Args:
        enabled: True to enable IERS Delta T, False to disable,
                 or None to use the environment variable.

    Note:
        - IERS data must be downloaded first using `iers_data.download_iers_finals()`
        - For dates outside the IERS data range, the standard Delta T model is used
        - Enable `iers_data.set_iers_auto_download(True)` for automatic download

    Environment Variable:
        LIBEPHEMERIS_IERS_DELTA_T: Set to "1", "true", or "yes" to enable.

    Example:
        >>> from libephemeris import set_iers_delta_t_enabled, get_iers_delta_t_enabled
        >>> set_iers_delta_t_enabled(True)  # Enable IERS Delta T
        >>> get_iers_delta_t_enabled()
        True
    """
    global _IERS_DELTA_T_ENABLED
    _IERS_DELTA_T_ENABLED = enabled


def get_iers_delta_t_enabled() -> bool:
    """
    Get the current IERS Delta T setting.

    Returns:
        True if IERS Delta T is enabled, False otherwise.

    Note:
        - If set_iers_delta_t_enabled() was called with an explicit value, returns that.
        - Otherwise, checks the LIBEPHEMERIS_IERS_DELTA_T environment variable.
        - If the environment variable is not set, returns False (disabled by default).

    Example:
        >>> from libephemeris import get_iers_delta_t_enabled
        >>> get_iers_delta_t_enabled()  # Default is False
        False
        >>> import os
        >>> os.environ['LIBEPHEMERIS_IERS_DELTA_T'] = '1'
        >>> get_iers_delta_t_enabled()  # Now enabled via env var
        True
    """
    # If explicitly set via function, use that value
    if _IERS_DELTA_T_ENABLED is not None:
        return _IERS_DELTA_T_ENABLED

    # Otherwise check environment variable
    env_value = os.environ.get(_IERS_DELTA_T_ENV_VAR, "").lower().strip()
    if env_value in ("1", "true", "yes", "on", "enabled"):
        return True
    if env_value in ("0", "false", "no", "off", "disabled"):
        return False

    # TOML config fallback
    from ._config_toml import get_bool as _toml_bool

    toml_value = _toml_bool("iers_delta_t")
    if toml_value is not None:
        return toml_value

    return False


# =============================================================================
# STRICT PRECISION MODE
# =============================================================================

# Global state for strict precision mode
_STRICT_PRECISION: Optional[bool] = None

# Environment variable name
_STRICT_PRECISION_ENV_VAR = "LIBEPHEMERIS_STRICT_PRECISION"


def set_strict_precision(enabled: Optional[bool]) -> None:
    """Enable or disable strict precision mode for major asteroids.

    When strict precision mode is enabled (the default), libephemeris will
    raise SPKRequiredError when calculating positions for major asteroids
    (Ceres, Pallas, Juno, Vesta, Chiron) without an SPK kernel registered.

    These bodies have strongly perturbed orbits that cannot be accurately
    modeled with Keplerian elements. The fallback Keplerian calculation can
    have errors of 1-10 degrees, which is unacceptable for most applications.

    Args:
        enabled: True to enable strict mode, False to disable,
                 None to reset to default (check environment variable).

    Example:
        >>> import libephemeris as eph
        >>> eph.set_strict_precision(True)  # Require SPK for major asteroids
        >>> eph.set_strict_precision(False)  # Allow Keplerian fallback
        >>> eph.set_strict_precision(None)  # Reset to default/env var

    See Also:
        get_strict_precision: Get current strict precision setting
        set_auto_spk_download: Enable automatic SPK downloading
        download_and_register_spk: Manually download and register SPK
    """
    global _STRICT_PRECISION
    _STRICT_PRECISION = enabled


def get_strict_precision() -> bool:
    """Get whether strict precision mode is enabled.

    In strict precision mode, SPK kernels are required for major asteroids
    (Ceres, Pallas, Juno, Vesta, Chiron). Without SPK, SPKRequiredError is
    raised instead of falling back to imprecise Keplerian calculations.

    Returns:
        True if strict precision mode is enabled, False otherwise.

    Note:
        - If set_strict_precision() was called with an explicit value, returns that.
        - Otherwise, checks the LIBEPHEMERIS_STRICT_PRECISION environment variable.
        - If the environment variable is not set, returns True (enabled by default).
        - Set LIBEPHEMERIS_STRICT_PRECISION=0 to disable strict mode via env var.

    Example:
        >>> from libephemeris import get_strict_precision
        >>> get_strict_precision()  # Default is True
        True
        >>> import os
        >>> os.environ['LIBEPHEMERIS_STRICT_PRECISION'] = '0'
        >>> get_strict_precision()  # Now disabled via env var
        False
    """
    # If explicitly set via function, use that value
    if _STRICT_PRECISION is not None:
        return _STRICT_PRECISION

    # Otherwise check environment variable
    env_value = os.environ.get(_STRICT_PRECISION_ENV_VAR, "").lower().strip()

    # If env var explicitly disables, return False
    if env_value in ("0", "false", "no", "off", "disabled"):
        return False
    # If env var explicitly enables, return True
    if env_value in ("1", "true", "yes", "on", "enabled"):
        return True

    # TOML config fallback
    from ._config_toml import get_bool as _toml_bool

    toml_value = _toml_bool("strict_precision")
    if toml_value is not None:
        return toml_value

    # Default to True (strict mode enabled) because:
    # - SPK type 21 is now fully supported via spktype21 library
    # - Auto-download works with JPL Horizons API
    # - Keplerian fallback has poor accuracy (1-10 degrees for asteroids)
    # Users can explicitly disable strict mode if needed
    return True


# =============================================================================
# SPK KERNEL MANAGEMENT
# =============================================================================


def _load_spk_kernel(filepath: str) -> SpiceKernel:
    """
    Load and cache an SPK kernel file.

    Internal function used by spk.py to load SPK kernels for minor bodies.

    Args:
        filepath: Full path to the SPK file

    Returns:
        SpiceKernel: Loaded Skyfield SpiceKernel object

    Raises:
        FileNotFoundError: If file does not exist
        ValueError: If file cannot be loaded as SPK
    """
    global _SPK_KERNELS

    # Return cached kernel if available
    if filepath in _SPK_KERNELS:
        return _SPK_KERNELS[filepath]

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"SPK file not found: {filepath}")

    # Load using Skyfield loader
    load = get_loader()
    try:
        kernel = load(filepath)
        _SPK_KERNELS[filepath] = kernel
        return kernel
    except (FileNotFoundError, ValueError, OSError) as e:
        raise ValueError(f"Failed to load SPK file {filepath}: {e}") from e


def _get_spk_target(ipl: int):
    """
    Get Skyfield target object from registered SPK for a body.

    Internal function used by planets.py for SPK-based calculations.

    Args:
        ipl: libephemeris body ID (e.g., CHIRON)

    Returns:
        Skyfield ephemeris target object, or None if not registered.
    """
    if ipl not in _SPK_BODY_MAP:
        return None

    spk_file, naif_id = _SPK_BODY_MAP[ipl]

    kernel = _SPK_KERNELS.get(spk_file)
    if kernel is None:
        return None

    try:
        return kernel[naif_id]
    except KeyError:
        return None
