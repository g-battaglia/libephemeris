"""
Global state management for libephemeris.

This module maintains the library's singleton state including:
- Ephemeris data loader (Skyfield Loader)
- Planetary ephemeris (DE421 or other JPL files)
- Timescale (for UTC/TT conversions)
- Observer topocentric location
- Sidereal mode configuration
- Cached angles for Arabic parts calculations

All state is stored in module-level globals to provide SwissEphemeris-compatible
stateful API behavior. This is thread-unsafe by design, matching SwissEph behavior.
"""

import os
from typing import Optional, Union
from skyfield.api import Loader, Topos
from skyfield.timelib import Timescale
from skyfield.jpllib import SpiceKernel

# =============================================================================
# GLOBAL STATE VARIABLES
# =============================================================================

_EPHEMERIS_PATH: Optional[str] = None  # Custom ephemeris directory
_EPHEMERIS_FILE: str = "de421.bsp"  # Ephemeris file to use (default: DE421)
_LOADER: Optional[Loader] = None  # Skyfield data loader
_PLANETS: Optional[SpiceKernel]  = None  # Loaded planetary ephemeris
_TS: Optional[Timescale] = None  # Timescale object
_TOPO: Optional[Topos] = None  # Observer location
_SIDEREAL_MODE: Optional[int] = None  # Active sidereal mode ID
_SIDEREAL_AYAN_T0: float = 0.0  # Ayanamsha value at reference epoch
_SIDEREAL_T0: float = 0.0  # Reference epoch (JD) for ayanamsha
_ANGLES_CACHE: dict[str, float] = {}  # Pre-calculated angles {name: longitude}


def get_loader() -> Loader:
    """
    Get or create the Skyfield data loader.
    
    Returns:
        Loader: Skyfield Loader instance for downloading/caching ephemeris files
        
    Note:
        Data files are cached in the parent directory of this module by default.
    """
    global _LOADER
    if _LOADER is None:
        data_dir = os.path.join(os.path.dirname(__file__), "..")
        _LOADER = Loader(data_dir)
    return _LOADER


def get_timescale() -> Timescale:
    """
    Get or create the Skyfield timescale object.
    
    Returns:
        Timescale: Skyfield timescale for time conversions (UTC, TT, etc.)
        
    Note:
        Automatically downloads IERS data for Delta T calculations if needed.
    """
    global _TS
    if _TS is None:
        load = get_loader()
        _TS = load.timescale()
    return _TS


def get_planets() -> SpiceKernel:
    """
    Get or load the planetary ephemeris (DE421 by default).
    
    Returns:
        SpiceKernel: Loaded JPL ephemeris kernel containing planetary positions
        
    Raises:
        FileNotFoundError: If ephemeris file cannot be found or downloaded
        
    Note:
        Uses the ephemeris file set via set_ephemeris_file() (default: de421.bsp).
        Searches in _EPHEMERIS_PATH if set, then workspace root, then downloads.
    """
    global _PLANETS
    if _PLANETS is None:
        load = get_loader()
        
        # Try custom ephemeris path first if set
        if _EPHEMERIS_PATH:
            bsp_path = os.path.join(_EPHEMERIS_PATH, _EPHEMERIS_FILE)
            if os.path.exists(bsp_path):
                _PLANETS = load(bsp_path)
                return _PLANETS
        
        # Try workspace root
        base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        bsp_path = os.path.join(base_dir, _EPHEMERIS_FILE)
        if os.path.exists(bsp_path):
            _PLANETS = load(bsp_path)
        else:
            # Download from internet
            _PLANETS = load(_EPHEMERIS_FILE)
    return _PLANETS


def set_topo(lon: float, lat: float, alt: float) -> None:
    """
    Set observer's topocentric location for planet calculations.
    
    Args:
        lon: Geographic longitude in degrees (East positive)
        lat: Geographic latitude in degrees (North positive)
        alt: Elevation above sea level in meters
        
    Note:
        Required for topocentric calculations (SEFLG_TOPOCTR),
        angles (Ascendant, MC), and Arabic parts.
    """
    global _TOPO
    _TOPO = Topos(latitude_degrees=lat, longitude_degrees=lon, elevation_m=alt)


def get_topo() -> Optional[Topos]:
    """
    Get the observer's topocentric location.
    
    Returns:
        Optional[Topos]: Current observer location or None if not set
    """
    return _TOPO


def set_sid_mode(mode: int, t0: float = 0.0, ayan_t0: float = 0.0) -> None:
    """
    Set the sidereal mode (ayanamsha system) for calculations.
    
    Args:
        mode: Sidereal mode ID (SE_SIDM_*) or 255 for custom
        t0: Reference epoch (Julian Day) for custom ayanamsha (default: J2000.0)
        ayan_t0: Ayanamsha value at t0 in degrees (for custom mode)
        
    Note:
        Affects all position calculations when SEFLG_SIDEREAL is set.
        Default is Lahiri (SE_SIDM_LAHIRI) if never set.
    """
    global _SIDEREAL_MODE, _SIDEREAL_T0, _SIDEREAL_AYAN_T0
    _SIDEREAL_MODE = mode
    _SIDEREAL_T0 = t0 if t0 != 0.0 else 2451545.0
    _SIDEREAL_AYAN_T0 = ayan_t0


def get_sid_mode(full: bool = False) -> Union[int, tuple[int, float, float]]:
    """
    Get the current sidereal mode configuration.
    
    Args:
        full: If True, return (mode, t0, ayan_t0); if False, return only mode ID
        
    Returns:
        int or tuple: Sidereal mode ID, or full configuration tuple
        
    Note:
        Returns SE_SIDM_LAHIRI (1) by default if never set.
    """
    if full:
        return _SIDEREAL_MODE, _SIDEREAL_T0, _SIDEREAL_AYAN_T0
    return _SIDEREAL_MODE if _SIDEREAL_MODE is not None else 1


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

    _ANGLES_CACHE = {}


def set_ephe_path(path: Optional[str]) -> None:
    """
    Set the path for ephemeris files.
    
    Args:
        path: Path to directory containing ephemeris files.
        
    Note:
        This sets the directory where get_planets() will look for the ephemeris
        file specified by set_ephemeris_file(). If the file is not found there,
        it will fall back to the workspace root and then download if needed.
    """
    global _EPHEMERIS_PATH, _PLANETS
    _EPHEMERIS_PATH = path
    # Clear cached planets to force reload from new path
    _PLANETS = None


def set_ephemeris_file(filename: str) -> None:
    """
    Set the ephemeris file to use for planetary calculations.
    
    Args:
        filename: Name of the JPL ephemeris file (e.g., "de421.bsp", "de422.bsp", "de431.bsp")
        
    Note:
        This allows using different ephemeris files with varying date ranges:
        - de421.bsp: 1900-2050 (default, 16 MB)
        - de422.bsp: -3000-3000 (623 MB)
        - de430.bsp: 1550-2650 (128 MB)
        - de431.bsp: -13200-17191 (3.4 GB)
        
        The file will be searched in:
        1. The path set by set_ephe_path() if configured
        2. The workspace root directory
        3. Downloaded from JPL if not found locally
        
        Changing the ephemeris file clears the cached planets and forces a reload.
    """
    global _EPHEMERIS_FILE, _PLANETS
    _EPHEMERIS_FILE = filename
    # Clear cached planets to force reload with new file
    _PLANETS = None
