"""
Thread-safe ephemeris context for libephemeris.

This module provides the EphemerisContext class which encapsulates all state
needed for ephemeris calculations in a thread-safe manner.

Each EphemerisContext instance maintains its own calculation state (observer
location, sidereal mode, angles cache) while sharing expensive resources
(ephemeris files, timescale data) globally in a thread-safe way.

Usage:
    >>> from libephemeris import EphemerisContext, SE_SUN
    >>> ctx = EphemerisContext()
    >>> ctx.set_topo(12.5, 41.9, 0)  # Rome
    >>> pos, _ = ctx.calc_ut(2451545.0, SE_SUN, 0)

Thread Safety:
    Multiple EphemerisContext instances can be used concurrently in different
    threads without interference. Each context has isolated state.
"""

import os
import threading
from typing import Optional, Tuple, Union
from skyfield.api import Loader, Topos
from skyfield.timelib import Timescale
from skyfield.jpllib import SpiceKernel


# =============================================================================
# SHARED RESOURCES (Thread-Safe Singleton Pattern)
# =============================================================================
# These resources are expensive to load (~16MB+ for ephemeris files)
# and thread-safe internally, so we share them across all contexts.

_SHARED_LOADER: Optional[Loader] = None
_SHARED_PLANETS: Optional[SpiceKernel] = None
_SHARED_TS: Optional[Timescale] = None
_SHARED_EPHE_PATH: Optional[str] = None
_SHARED_EPHE_FILE: str = "de421.bsp"
_SHARED_LOCK = threading.Lock()


class EphemerisContext:
    """
    Thread-safe context for ephemeris calculations.
    
    Each instance maintains independent calculation state while sharing
    expensive resources (ephemeris files, timescale) across instances.
    
    Attributes:
        topo: Observer topocentric location (or None for geocentric)
        sidereal_mode: Active sidereal mode ID (1 = Lahiri by default)
        sidereal_t0: Reference epoch for custom ayanamsha (JD)
        sidereal_ayan_t0: Ayanamsha value at reference epoch (degrees)
        
    Example:
        >>> ctx = EphemerisContext()
        >>> ctx.set_topo(12.5, 41.9, 0)  # Rome coordinates
        >>> ctx.set_sid_mode(1)  # Lahiri ayanamsha
        >>> pos, flags = ctx.calc_ut(2451545.0, SE_MARS, SEFLG_SIDEREAL)
        >>> print(f"Mars at {pos[0]:.2f}° sidereal longitude")
        
    Thread Safety:
        Each context instance is independent. Multiple threads can create
        and use their own contexts without interference.
    """
    
    def __init__(self, ephe_path: Optional[str] = None, 
                 ephe_file: str = "de421.bsp"):
        """
        Initialize a new ephemeris context.
        
        Args:
            ephe_path: Optional path to directory containing ephemeris files.
                      If None, uses default workspace directory.
            ephe_file: Ephemeris file to use (default: "de421.bsp")
        """
        # Instance-specific state (NOT shared between contexts)
        self.topo: Optional[Topos] = None
        self.sidereal_mode: int = 1  # Default: Lahiri (SE_SIDM_LAHIRI)
        self.sidereal_t0: float = 2451545.0  # J2000.0
        self.sidereal_ayan_t0: float = 0.0
        self._angles_cache: dict[str, float] = {}
        
        # Ephemeris configuration (for this context)
        self._ephe_path = ephe_path
        self._ephe_file = ephe_file
        
        # Update shared config if needed
        global _SHARED_EPHE_PATH, _SHARED_EPHE_FILE
        if ephe_path is not None:
            _SHARED_EPHE_PATH = ephe_path
        if ephe_file != "de421.bsp":
            _SHARED_EPHE_FILE = ephe_file
    
    def get_loader(self) -> Loader:
        """
        Get shared Skyfield data loader.
        
        Thread-safe lazy initialization. All contexts share the same loader.
        
        Returns:
            Loader: Skyfield Loader instance for downloading/caching files
        """
        global _SHARED_LOADER
        if _SHARED_LOADER is None:
            with _SHARED_LOCK:
                if _SHARED_LOADER is None:  # Double-checked locking
                    data_dir = os.path.join(os.path.dirname(__file__), "..")
                    _SHARED_LOADER = Loader(data_dir)
        return _SHARED_LOADER
    
    def get_timescale(self) -> Timescale:
        """
        Get shared Skyfield timescale object.
        
        Thread-safe lazy initialization. All contexts share the same timescale.
        
        Returns:
            Timescale: Skyfield timescale for time conversions (UTC, TT, etc.)
        """
        global _SHARED_TS
        if _SHARED_TS is None:
            with _SHARED_LOCK:
                if _SHARED_TS is None:  # Double-checked locking
                    load = self.get_loader()
                    _SHARED_TS = load.timescale()
        return _SHARED_TS
    
    def get_planets(self) -> SpiceKernel:
        """Get shared planetary ephemeris.
        
        Thread-safe lazy loading of JPL ephemeris file. All contexts share
        the same ephemeris data to save memory.
        
        Returns:
            SpiceKernel: Loaded JPL ephemeris kernel (DE421 by default)
            
        Raises:
            FileNotFoundError: If ephemeris file cannot be found or downloaded
        """
        global _SHARED_PLANETS
        if _SHARED_PLANETS is None:
            with _SHARED_LOCK:
                if _SHARED_PLANETS is None:  # Double-checked locking
                    load = self.get_loader()
                    
                    # Try custom ephemeris path first if set
                    if _SHARED_EPHE_PATH:
                        bsp_path = os.path.join(_SHARED_EPHE_PATH, _SHARED_EPHE_FILE)
                        if os.path.exists(bsp_path):
                            _SHARED_PLANETS = load(bsp_path)
                            return _SHARED_PLANETS
                    
                    # Try workspace root
                    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
                    bsp_path = os.path.join(base_dir, _SHARED_EPHE_FILE)
                    if os.path.exists(bsp_path):
                        _SHARED_PLANETS = load(bsp_path)
                    else:
                        # Download from internet
                        _SHARED_PLANETS = load(_SHARED_EPHE_FILE)
        return _SHARED_PLANETS
    
    def set_topo(self, lon: float, lat: float, alt: float) -> None:
        """
        Set observer's topocentric location for calculations.
        
        Args:
            lon: Geographic longitude in degrees (East positive)
            lat: Geographic latitude in degrees (North positive)
            alt: Elevation above sea level in meters
            
        Note:
            Required for topocentric calculations (SEFLG_TOPOCTR),
            angles (Ascendant, MC), and Arabic parts.
        """
        self.topo = Topos(latitude_degrees=lat, 
                         longitude_degrees=lon, 
                         elevation_m=alt)
    
    def get_topo(self) -> Optional[Topos]:
        """
        Get the observer's topocentric location.
        
        Returns:
            Optional[Topos]: Current observer location or None if not set
        """
        return self.topo
    
    def set_sid_mode(self, mode: int, t0: float = 0.0, ayan_t0: float = 0.0) -> None:
        """
        Set the sidereal mode (ayanamsha system) for calculations.
        
        Args:
            mode: Sidereal mode ID (SE_SIDM_*) or 255 for custom
            t0: Reference epoch (Julian Day) for custom ayanamsha (default: J2000.0)
            ayan_t0: Ayanamsha value at t0 in degrees (for custom mode)
            
        Note:
            Affects all position calculations when SEFLG_SIDEREAL is set.
            Default is Lahiri (SE_SIDM_LAHIRI = 1) if never set.
        """
        self.sidereal_mode = mode
        self.sidereal_t0 = t0 if t0 != 0.0 else 2451545.0
        self.sidereal_ayan_t0 = ayan_t0
    
    def get_sid_mode(self, full: bool = False) -> Union[int, Tuple[int, float, float]]:
        """
        Get the current sidereal mode configuration.
        
        Args:
            full: If True, return (mode, t0, ayan_t0); if False, return only mode ID
            
        Returns:
            int or tuple: Sidereal mode ID, or full configuration tuple
        """
        if full:
            return self.sidereal_mode, self.sidereal_t0, self.sidereal_ayan_t0
        return self.sidereal_mode
    
    def get_angles_cache(self) -> dict[str, float]:
        """
        Get cached astrological angles for the current calculation context.
        
        Returns:
            dict: Cached angles {name: longitude_degrees}
            
        Note:
            Used by Arabic parts calculations which require pre-calculated
            planetary positions and angles.
        """
        return self._angles_cache
    
    def set_angles_cache(self, angles: dict[str, float]) -> None:
        """
        Cache pre-calculated angles for use in Arabic parts.
        
        Args:
            angles: Dictionary of angles {name: longitude_degrees}
                   e.g., {"Sun": 120.5, "Moon": 240.3, "Asc": 15.7}
                   
        Note:
            Creates a copy to prevent external mutation of cache.
        """
        self._angles_cache = angles.copy()
    
    def clear_angles_cache(self) -> None:
        """
        Clear the angles cache.
        
        Use this between unrelated calculation contexts to prevent stale data.
        """
        self._angles_cache = {}
    
    # =========================================================================
    # Calculation Methods (delegate to planets.py with context)
    # =========================================================================
    
    def calc_ut(self, tjd_ut: float, ipl: int, iflag: int) -> Tuple[Tuple[float, float, float, float, float, float], int]:
        """
        Calculate planetary position for Universal Time.
        
        Thread-safe calculation using this context's state.
        
        Args:
            tjd_ut: Julian Day in Universal Time (UT1)
            ipl: Planet/body ID (SE_SUN, SE_MOON, etc.)
            iflag: Calculation flags (SEFLG_SPEED, SEFLG_HELCTR, etc.)
            
        Returns:
            Tuple containing:
                - Position tuple: (longitude, latitude, distance, 
                                  speed_lon, speed_lat, speed_dist)
                - Return flag: iflag value on success
                
        Example:
            >>> pos, retflag = ctx.calc_ut(2451545.0, SE_MARS, SEFLG_SPEED)
            >>> lon, lat, dist = pos[0], pos[1], pos[2]
        """
        from .planets import _calc_body_with_context
        ts = self.get_timescale()
        t = ts.ut1_jd(tjd_ut)
        return _calc_body_with_context(t, ipl, iflag, self)
    
    def calc(self, tjd: float, ipl: int, iflag: int) -> Tuple[Tuple[float, float, float, float, float, float], int]:
        """
        Calculate planetary position for Ephemeris Time (ET/TT).
        
        Thread-safe calculation using this context's state.
        
        Args:
            tjd: Julian Day in Terrestrial Time (TT/ET)
            ipl: Planet/body ID (SE_SUN, SE_MOON, etc.)
            iflag: Calculation flags (SEFLG_SPEED, SEFLG_HELCTR, etc.)
            
        Returns:
            Tuple containing:
                - Position tuple: (longitude, latitude, distance,
                                  speed_lon, speed_lat, speed_dist)
                - Return flag: iflag value on success
                
        Note:
            TT differs from UT by Delta T (~32s for year 2000).
            For most astrological applications, use calc_ut() instead.
        """
        from .planets import _calc_body_with_context
        ts = self.get_timescale()
        t = ts.tt_jd(tjd)
        return _calc_body_with_context(t, ipl, iflag, self)
    
    def houses(self, tjd_ut: float, lat: float, lon: float, hsys: int) -> Tuple[list[float], list[float]]:
        """
        Calculate house cusps and angles.
        
        Thread-safe house calculation using this context.
        
        Args:
            tjd_ut: Julian Day in Universal Time
            lat: Geographic latitude in degrees
            lon: Geographic longitude in degrees
            hsys: House system character code (ord('P') for Placidus, etc.)
            
        Returns:
            Tuple of (cusps, ascmc) where:
                - cusps: List of 12 house cusp longitudes [1-12]
                - ascmc: List of angles [ASC, MC, ARMC, Vertex, ...]
        """
        from .houses import _swe_houses_with_context
        return _swe_houses_with_context(tjd_ut, lat, lon, hsys, self)
