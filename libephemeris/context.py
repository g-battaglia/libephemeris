# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Thread-safe ephemeris context for libephemeris.

This module provides the EphemerisContext class which encapsulates all state
needed for ephemeris calculations in a thread-safe manner.

Each EphemerisContext instance maintains its own calculation state (observer
location, sidereal mode, angles cache) while sharing expensive resources
(ephemeris files, timescale data) globally in a thread-safe way.

Usage:
    >>> from libephemeris import EphemerisContext, SUN
    >>> ctx = EphemerisContext()
    >>> ctx.set_topo(12.5, 41.9, 0)  # Rome
    >>> pos, _ = ctx.calc_ut(2451545.0, SUN, 0)

Thread Safety:
    Multiple EphemerisContext instances can be used concurrently in different
    threads without interference. Each context has isolated state.
"""

from __future__ import annotations

import os
import threading
from typing import TYPE_CHECKING, Literal, Optional, Tuple, Union, overload

from .constants import MEAN_NODE, TRUE_NODE
from .tracing import _record
from skyfield.api import Loader, Topos
from skyfield.timelib import Timescale
from skyfield.jpllib import SpiceKernel

if TYPE_CHECKING:
    from .leb_composite import CompositeLEBReader
    from .leb2_reader import LEB2Reader
    from .leb_reader import LEBReader


# =============================================================================
# SHARED RESOURCES (Thread-Safe Singleton Pattern)
# =============================================================================
# These resources are expensive to load (~16MB+ for ephemeris files)
# and thread-safe internally, so we share them across all contexts.

_SHARED_LOADER: Optional[Loader] = None
_SHARED_PLANETS: Optional[SpiceKernel] = None
_SHARED_TS: Optional[Timescale] = None
_SHARED_EPHE_PATH: Optional[str] = None
_SHARED_EPHE_FILE: str = "de440.bsp"
_SHARED_LOCK = threading.RLock()


class EphemerisContext:
    """
    Thread-safe context for ephemeris calculations.

    Each instance maintains independent calculation state while sharing
    expensive resources (ephemeris files, timescale) across instances.

    Attributes:
        topo: Observer topocentric location (or None for geocentric)
        sidereal_mode: Active sidereal mode ID (0 = Fagan/Bradley by default,
            matching the module-level API's default when set_sid_mode is never
            called)
        sidereal_t0: Reference epoch for custom ayanamsha (JD)
        sidereal_ayan_t0: Ayanamsha value at reference epoch (degrees)

    Example:
        >>> ctx = EphemerisContext()
        >>> ctx.set_topo(12.5, 41.9, 0)  # Rome coordinates
        >>> ctx.set_sid_mode(1)  # Lahiri ayanamsha
        >>> pos, flags = ctx.calc_ut(2451545.0, MARS, FLG_SIDEREAL)
        >>> print(f"Mars at {pos[0]:.2f}° sidereal longitude")

    Thread Safety:
        Each context instance is independent. Multiple threads can create
        and use their own contexts without interference.
    """

    def __init__(self, ephe_path: Optional[str] = None, ephe_file: str = "de440.bsp"):
        """
        Initialize a new ephemeris context.

        Args:
            ephe_path: Optional path to directory containing ephemeris files.
                      If None, uses default workspace directory.
            ephe_file: Ephemeris file to use (default: "de440.bsp").
                      Supported files: "de440s.bsp" (1849-2150, lightweight),
                      "de440.bsp" (1550-2650, default), "de441.bsp" (-13200 to
                      +17191, extended range).

        Note:
            ``ephe_path``/``ephe_file`` select the *process-wide* shared kernel,
            not a per-context one: all contexts deliberately share a single
            loaded JPL kernel to save memory (see ``get_planets``). The first
            context whose calculation triggers the load wins — a later context
            constructed with a different ``ephe_file`` updates the shared
            selection only if the kernel has not been loaded yet, otherwise the
            already-loaded kernel is reused and the new ``ephe_file`` is ignored.
            Calculation state (topo, sidereal mode, SPK map) IS per-context.
        """
        # Instance-specific state (NOT shared between contexts)
        self.topo: Optional[Topos] = None
        # Default Fagan/Bradley (0) to match the module-level API: get_sid_mode()
        # returns 0 when set_sid_mode was never called, so a fresh context and the
        # module API must agree for an unset FLG_SIDEREAL request.
        self.sidereal_mode: int = 0  # Default: Fagan/Bradley (SIDM_FAGAN_BRADLEY)
        self.sidereal_t0: float = 2451545.0  # J2000.0
        self.sidereal_ayan_t0: float = 0.0
        self._angles_cache: dict[str, float] = {}
        self._spk_body_map: dict[
            int, tuple[str, int]
        ] = {}  # Context-local SPK mappings

        # LEB binary ephemeris configuration (per-context)
        self._leb_file: Optional[str] = None
        self._leb_reader: Optional["LEBReader | LEB2Reader | CompositeLEBReader"] = None

        # Ephemeris configuration (for this context)
        self._ephe_path = ephe_path
        self._ephe_file = ephe_file

        # Update shared config if needed
        global _SHARED_EPHE_PATH, _SHARED_EPHE_FILE
        if ephe_path is not None:
            _SHARED_EPHE_PATH = ephe_path
        if ephe_file != "de440.bsp":
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
                if _SHARED_LOADER is None:  # pragma: no branch - double-checked locking; the race arc only triggers under concurrent first-init
                    from .state import _get_data_dir

                    _SHARED_LOADER = Loader(_get_data_dir())
        return _SHARED_LOADER

    def get_timescale(self) -> Timescale:
        """
        Get shared Skyfield timescale object.

        Thread-safe lazy initialization. All contexts share the same
        timescale — and crucially the SAME one as the module-level API:
        state.get_timescale() builds an enhanced timescale with merged
        historical Delta-T (1657-1973). A plain load.timescale() here
        diverged from module calc_ut by up to ~0.7 s of Delta-T (~0.4"
        of Moon longitude) for pre-1973 dates.

        Returns:
            Timescale: Skyfield timescale for time conversions (UTC, TT, etc.)
        """
        global _SHARED_TS
        if _SHARED_TS is None:
            with _SHARED_LOCK:
                if _SHARED_TS is None:  # Double-checked locking
                    from .state import get_timescale as _state_get_timescale

                    _SHARED_TS = _state_get_timescale()
        return _SHARED_TS

    def get_planets(self) -> SpiceKernel:
        """Get shared planetary ephemeris.

        Thread-safe lazy loading of JPL ephemeris file. All contexts share
        the same ephemeris data to save memory.

        Returns:
            SpiceKernel: Loaded JPL ephemeris kernel (DE440 by default)

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

                    # Load from data dir (downloads automatically if missing)
                    from .state import _get_data_dir

                    data_dir = _get_data_dir()
                    bsp_path = os.path.join(data_dir, _SHARED_EPHE_FILE)
                    if os.path.exists(bsp_path):
                        _SHARED_PLANETS = load(bsp_path)
                    else:
                        _SHARED_PLANETS = load(_SHARED_EPHE_FILE)
        return _SHARED_PLANETS

    def set_leb_file(self, filepath: Optional[str]) -> None:
        """Set the .leb file path for binary ephemeris mode.

        When a .leb file is configured, calc_ut() and calc() will attempt
        to use precomputed Chebyshev polynomials for fast evaluation before
        falling back to the Skyfield pipeline.

        Args:
            filepath: Path to a .leb file, or None to disable binary mode.
        """
        if self._leb_reader is not None:
            try:
                self._leb_reader.close()
            except (AttributeError, OSError):
                pass
        self._leb_file = filepath
        self._leb_reader = None

    def get_leb_reader(self) -> Optional["LEBReader | LEB2Reader | CompositeLEBReader"]:
        """Get the active LEBReader for this context, if any.

        If the .leb file path is invalid or the file is corrupted,
        logs a warning and returns None (silent fallback to Skyfield).

        Returns:
            LEBReader instance if a .leb file is configured and valid,
            None otherwise.
        """
        if self._leb_reader is None and self._leb_file is not None:
            try:
                from .leb_reader import open_leb

                self._leb_reader = open_leb(self._leb_file)
            except (FileNotFoundError, ValueError, OSError) as e:
                from .logging_config import get_logger

                get_logger().warning(
                    "Failed to open LEB file %s: %s", self._leb_file, e
                )
                return None
        return self._leb_reader

    def set_topo(self, lon: float, lat: float, alt: float) -> None:
        """
        Set observer's topocentric location for calculations.

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
            >>> ctx = EphemerisContext()
            >>> ctx.set_topo(12.5, 41.9, 0)  # Rome
        """
        from .exceptions import validate_coordinates

        validate_coordinates(lat, lon, "set_topo")
        self.topo = Topos(latitude_degrees=lat, longitude_degrees=lon, elevation_m=alt)

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
            mode: Sidereal mode ID (SIDM_*) or 255 for custom
            t0: Reference epoch (Julian Day) for custom ayanamsha (default: J2000.0)
            ayan_t0: Ayanamsha value at t0 in degrees (for custom mode)

        Note:
            Affects all position calculations when FLG_SIDEREAL is set.
            Default is Fagan/Bradley (SIDM_FAGAN_BRADLEY = 0) if never set.
        """
        self.sidereal_mode = mode
        self.sidereal_t0 = t0 if t0 != 0.0 else 2451545.0
        self.sidereal_ayan_t0 = ayan_t0

    @overload
    def get_sid_mode(self, full: Literal[True]) -> Tuple[int, float, float]: ...

    @overload
    def get_sid_mode(self, full: Literal[False] = ...) -> int: ...

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
    # SPK Body Registration (Context-local)
    # =========================================================================

    def register_spk_body(
        self, ipl: int, spk_file: str, naif_id: Union[int, str]
    ) -> None:
        """
        Register a mapping between a body ID and an SPK kernel target.

        This is a context-local registration that only affects calculations
        performed through this context instance.

        Args:
            ipl: libephemeris body ID (e.g., CHIRON, ERIS)
            spk_file: Path to the SPK file, or filename if in library path
            naif_id: NAIF ID of the body in the SPK kernel

        Raises:
            FileNotFoundError: If SPK file not found
            ValueError: If naif_id not found in SPK kernel

        Example:
            >>> ctx = EphemerisContext()
            >>> ctx.register_spk_body(CHIRON, "chiron.bsp", 2002060)
        """
        from . import state

        # Convert naif_id to int if string
        if isinstance(naif_id, str):
            naif_id = int(naif_id)

        # Resolve file path
        if not os.path.isabs(spk_file):
            lib_path = state.get_library_path()
            full_path = os.path.join(lib_path, spk_file)
            if os.path.exists(full_path):
                spk_file = full_path
            elif not os.path.exists(spk_file):
                raise FileNotFoundError(f"SPK file not found: {spk_file}")

        if not os.path.exists(spk_file):
            raise FileNotFoundError(f"SPK file not found: {spk_file}")

        # Load kernel using shared loader
        state._load_spk_kernel(spk_file)

        # Validate that naif_id exists in kernel
        kernel = state._SPK_KERNELS.get(spk_file)
        if kernel is not None:
            try:
                _ = kernel[naif_id]
            except KeyError:
                try:
                    _ = kernel[str(naif_id)]
                except KeyError:
                    available = []
                    if hasattr(kernel, "names"):
                        available = list(kernel.names())[:10]
                    raise ValueError(
                        f"NAIF ID {naif_id} not found in SPK kernel {spk_file}. "
                        f"Available targets (first 10): {available}"
                    )

        # Register in context-local map
        self._spk_body_map[ipl] = (spk_file, naif_id)

    def unregister_spk_body(self, ipl: int) -> None:
        """
        Remove SPK registration for a body from this context.

        Args:
            ipl: libephemeris body ID (e.g., CHIRON)
        """
        if ipl in self._spk_body_map:
            del self._spk_body_map[ipl]

    def get_spk_body_info(self, ipl: int) -> Optional[Tuple[str, int]]:
        """
        Get SPK registration info for a body in this context.

        First checks context-local registrations, then falls back to global.

        Args:
            ipl: libephemeris body ID

        Returns:
            Tuple of (spk_file, naif_id) if registered, None otherwise.
        """
        # Check context-local first
        if ipl in self._spk_body_map:
            return self._spk_body_map[ipl]
        # Fall back to global
        from . import state

        return state._SPK_BODY_MAP.get(ipl)

    def list_spk_bodies(self) -> dict[int, Tuple[str, int]]:
        """
        List all SPK body mappings visible to this context.

        Merges global registrations with context-local ones (local takes precedence).

        Returns:
            Dict mapping ipl -> (spk_file, naif_id)
        """
        from . import state

        # Start with global, overlay with local
        result = dict(state._SPK_BODY_MAP)
        result.update(self._spk_body_map)
        return result

    # =========================================================================
    # Calculation Methods (delegate to planets.py with context)
    # =========================================================================

    def calc_ut(
        self, tjd_ut: float, ipl: int, iflag: int
    ) -> Tuple[Tuple[float, float, float, float, float, float], int]:
        """
        Calculate planetary position for Universal Time.

        Thread-safe calculation using this context's state.

        If a .leb file is configured (via set_leb_file()), uses the fast
        binary ephemeris path first, falling back to Skyfield when the body
        is not in the .leb file or the JD is out of range.

        Args:
            tjd_ut: Julian Day in Universal Time (UT1)
            ipl: Planet/body ID (SUN, MOON, etc.)
            iflag: Calculation flags (FLG_SPEED, FLG_HELCTR, etc.)

        Returns:
            Tuple containing:
                - Position tuple: (longitude, latitude, distance,
                                  speed_lon, speed_lat, speed_dist)
                - Return flag: iflag value on success

        Example:
            >>> pos, retflag = ctx.calc_ut(2451545.0, MARS, FLG_SPEED)
            >>> lon, lat, dist = pos[0], pos[1], pos[2]
        """
        # South nodes: derive from the north node via this same context path,
        # mirroring the module-level calc_ut(). The descending node must be
        # intercepted here because no downstream path derives it (_calc_body has
        # no antipode branch), and the antipode is representation-dependent
        # (FLG_XYZ / FLG_RADIANS) — _south_node_from_north handles all three.
        if ipl in (-MEAN_NODE, -TRUE_NODE):
            from .planets import _south_node_from_north

            north_result, retflag = self.calc_ut(tjd_ut, abs(ipl), iflag)
            return _south_node_from_north(north_result, iflag), retflag

        # --- LEB fast path: try binary ephemeris first ---
        reader = self.get_leb_reader()
        if reader is None:
            # Fall back to global reader if context has no .leb
            from . import state

            reader = state.get_leb_reader()

        if reader is not None:
            try:
                from . import fast_calc

                # Pass sidereal params explicitly (thread-safe, no global swap)
                _ctx_topo = None
                if self.topo is not None:
                    _ctx_topo = (
                        float(self.topo.longitude.degrees),
                        float(self.topo.latitude.degrees),
                        float(self.topo.elevation.m),
                    )
                result = fast_calc.fast_calc_ut(
                    reader,
                    tjd_ut,
                    ipl,
                    iflag,
                    sid_mode=self.sidereal_mode,
                    sid_t0=self.sidereal_t0,
                    sid_ayan_t0=self.sidereal_ayan_t0,
                    topo=_ctx_topo,
                )
                from .logging_config import get_logger

                get_logger().debug("body=%d jd=%.1f source=LEB (context)", ipl, tjd_ut)
                _record(ipl, "LEB")
                return result
            except (KeyError, ValueError):
                from .logging_config import get_logger

                get_logger().debug(
                    "body=%d jd=%.1f source=LEB->fallback (context)",
                    ipl,
                    tjd_ut,
                )
        # --- END LEB fast path ---

        from .planets import _calc_body_with_context

        ts = self.get_timescale()
        t = ts.ut1_jd(tjd_ut)
        return _calc_body_with_context(t, ipl, iflag, self)

    def calc(
        self, tjd: float, ipl: int, iflag: int
    ) -> Tuple[Tuple[float, float, float, float, float, float], int]:
        """
        Calculate planetary position for Ephemeris Time (ET/TT).

        Thread-safe calculation using this context's state.

        If a .leb file is configured (via set_leb_file()), uses the fast
        binary ephemeris path first, falling back to Skyfield when the body
        is not in the .leb file or the JD is out of range.

        Args:
            tjd: Julian Day in Terrestrial Time (TT/ET)
            ipl: Planet/body ID (SUN, MOON, etc.)
            iflag: Calculation flags (FLG_SPEED, FLG_HELCTR, etc.)

        Returns:
            Tuple containing:
                - Position tuple: (longitude, latitude, distance,
                                  speed_lon, speed_lat, speed_dist)
                - Return flag: iflag value on success

        Note:
            TT differs from UT by Delta T (~32s for year 2000).
            For most astrological applications, use calc_ut() instead.
        """
        # South nodes: derive from the north node via this same context path,
        # mirroring the module-level calc().
        if ipl in (-MEAN_NODE, -TRUE_NODE):
            from .planets import _south_node_from_north

            north_result, retflag = self.calc(tjd, abs(ipl), iflag)
            return _south_node_from_north(north_result, iflag), retflag

        # --- LEB fast path: try binary ephemeris first ---
        reader = self.get_leb_reader()
        if reader is None:
            from . import state

            reader = state.get_leb_reader()

        if reader is not None:
            try:
                from . import fast_calc

                # Pass sidereal params explicitly (thread-safe, no global swap)
                _ctx_topo_tt = None
                if self.topo is not None:
                    _ctx_topo_tt = (
                        float(self.topo.longitude.degrees),
                        float(self.topo.latitude.degrees),
                        float(self.topo.elevation.m),
                    )
                result = fast_calc.fast_calc_tt(
                    reader,
                    tjd,
                    ipl,
                    iflag,
                    sid_mode=self.sidereal_mode,
                    sid_t0=self.sidereal_t0,
                    sid_ayan_t0=self.sidereal_ayan_t0,
                    topo=_ctx_topo_tt,
                )
                from .logging_config import get_logger

                get_logger().debug("body=%d jd=%.1f source=LEB (context)", ipl, tjd)
                _record(ipl, "LEB")
                return result
            except (KeyError, ValueError):
                from .logging_config import get_logger

                get_logger().debug(
                    "body=%d jd=%.1f source=LEB->fallback (context)",
                    ipl,
                    tjd,
                )
        # --- END LEB fast path ---

        from .planets import _calc_body_with_context

        ts = self.get_timescale()
        t = ts.tt_jd(tjd)
        return _calc_body_with_context(t, ipl, iflag, self)

    def houses(
        self, tjd_ut: float, lat: float, lon: float, hsys: int
    ) -> Tuple[Tuple[float, ...], Tuple[float, ...]]:
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
        from .houses import _houses_with_context

        return _houses_with_context(tjd_ut, lat, lon, hsys, self)

    def calc_pctr(
        self, tjd_ut: float, ipl: int, iplctr: int, iflag: int
    ) -> Tuple[Tuple[float, float, float, float, float, float], int]:
        """
        Calculate planet-centric position (target as seen from another planet).

        Thread-safe calculation using this context's state.

        Args:
            tjd_ut: Julian Day in Universal Time (UT1)
            ipl: Target planet/body ID (SUN, MOON, etc.)
            iplctr: Observer/center planet ID (body from which to observe)
            iflag: Calculation flags (FLG_SPEED, etc.)

        Returns:
            Tuple containing:
                - Position tuple: (longitude, latitude, distance,
                                  speed_lon, speed_lat, speed_dist)
                - Return flag: iflag value on success

        Example:
            >>> # Position of Moon as seen from Mars
            >>> pos, retflag = ctx.calc_pctr(2451545.0, MOON, MARS, FLG_SPEED)
        """
        from .planets import _calc_body_pctr_with_context

        ts = self.get_timescale()
        t = ts.ut1_jd(tjd_ut)
        return _calc_body_pctr_with_context(t, ipl, iplctr, iflag, self)

    @classmethod
    def close(cls) -> None:
        """
        Close all shared ephemeris resources and release file handles.

        This class method closes the shared SPK kernel file handles and resets
        all shared resources to their initial values. Call this when you want to:
        - Free memory and file handles in long-running applications
        - Switch to a different ephemeris file
        - Ensure clean state in test suites

        Note:
            - This affects all EphemerisContext instances since they share
              resources (ephemeris files, timescale, loader)
            - After calling close(), the next calculation on any context
              will automatically reload resources as needed
            - The timescale cache is reset at the state level too, so the
              reload picks up fresh time data (IERS table, leap seconds)
              instead of the previously cached singleton
            - Instance-specific state (topo, sidereal_mode, angles_cache)
              is NOT affected - only shared resources are reset

        Example:
            >>> ctx = EphemerisContext()
            >>> pos, _ = ctx.calc_ut(2451545.0, SUN, 0)  # Loads ephemeris
            >>> EphemerisContext.close()  # Close shared files
            >>> pos, _ = ctx.calc_ut(2451545.0, SUN, 0)  # Reloads ephemeris
        """
        global _SHARED_LOADER, _SHARED_PLANETS, _SHARED_TS
        global _SHARED_EPHE_PATH, _SHARED_EPHE_FILE

        with _SHARED_LOCK:
            # Close the SPK kernel file handles if loaded
            if _SHARED_PLANETS is not None:
                try:
                    _SHARED_PLANETS.close()
                except (AttributeError, OSError, ValueError):
                    # SpiceKernel may not have close() in all versions,
                    # or may already be closed
                    pass

            # Reset all shared state to initial values
            _SHARED_LOADER = None
            _SHARED_PLANETS = None
            _SHARED_TS = None
            _SHARED_EPHE_PATH = None
            _SHARED_EPHE_FILE = "de440.bsp"

            # _SHARED_TS caches the singleton built by state.get_timescale();
            # clearing the reference alone would re-cache the same live
            # object on the next call. Reset the state-level cache too so
            # time data (IERS table, leap seconds) is actually rebuilt.
            # Lock order _SHARED_LOCK -> state._INIT_LOCK matches
            # get_timescale() above; state never takes locks in the
            # reverse order.
            from .state import _reset_timescale

            _reset_timescale()
