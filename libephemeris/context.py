# SPDX-License-Identifier: Apache-2.0
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
import warnings
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
                if (
                    _SHARED_LOADER is None
                ):  # pragma: no branch - double-checked locking; the race arc only triggers under concurrent first-init
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

    def set_topo(self, lon: float, lat: float, alt: float = 0.0) -> None:
        """
        Set observer's topocentric location for calculations.

        Args:
            lon: Geographic longitude in degrees (East positive, range -180 to 360)
            lat: Geographic latitude in degrees (North positive, range -90 to 90)
            alt: Elevation above sea level in meters

        Raises:
            CoordinateError: If latitude is outside [-90, 90] or
                            longitude is outside [-180, 360]

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
        # Apply the same SIDBIT guard as the global setter (state.set_sid_mode):
        # strip the SIDBIT_* projection flags (>= 256: ECL_T0, SSY_PLANE,
        # USER_UT, ECL_DATE, NO_PREC_OFFSET, PREC_ORIG), which select
        # alternative ecliptic/equinox projections not implemented here.
        # Keeping only the base ayanamsha mode avoids the silent wrong
        # fallback the composite value would otherwise trigger downstream
        # (and the LEB fast-path miss it causes). Warn so the caller knows
        # the projection was dropped rather than applied.
        sidbits = mode & ~0xFF
        if sidbits:
            warnings.warn(
                f"SIDBIT projection flags 0x{sidbits:04X} are not supported in "
                f"this version and were ignored; using base ayanamsha mode "
                f"{mode & 0xFF}.",
                UserWarning,
                stacklevel=2,
            )
            mode = mode & 0xFF

        self.sidereal_mode = mode
        # t0 is stored literally (no J2000 sentinel — JD 0.0 is a valid
        # ~4713 BCE epoch), matching state.set_sid_mode and the downstream
        # SIDM_USER consumers (planets._calc_ayanamsa, fast_calc), which all
        # use the raw t0 with method_b_accumulated_precession.
        self.sidereal_t0 = t0
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
        return self._calc_impl(tjd_ut, ipl, iflag, ut=True)

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
        return self._calc_impl(tjd, ipl, iflag, ut=False)

    def _calc_impl(
        self, tjd: float, ipl: int, iflag: int, *, ut: bool
    ) -> Tuple[Tuple[float, float, float, float, float, float], int]:
        """Shared implementation of calc_ut() and calc().

        calc_ut() and calc() differ only in the time scale of ``tjd`` —
        UT1 (``ut=True``) vs TT (``ut=False``) — which selects the
        fast_calc variant, the Skyfield time constructor, and the ECL_NUT
        helper. Keeping the two in one body makes their parity structural
        rather than maintained by hand across near-identical copies.

        Args:
            tjd: Julian Day, in UT1 if ``ut`` else TT.
            ipl: Planet/body ID (SUN, MOON, ...).
            iflag: Calculation flags.
            ut: True for the calc_ut() (UT1) entry point, False for calc()
                (TT).

        Returns:
            ``((lon, lat, dist, dlon, dlat, ddist), retflag)``.
        """
        from .constants import ECL_NUT
        from .planets import _normalize_calc_flags, _remap_ast_offset

        # Handle ECL_NUT (-1), like the module-level entry points: the
        # reference uses plaus_iflag-style flag handling here (MOSEPH kept,
        # SPEED3 not remapped), so it runs on the raw flags.
        if ipl == ECL_NUT:
            from .planets import (
                _calc_nutation_obliquity,
                _calc_nutation_obliquity_tt,
                _plaus_ephemeris_flags,
            )

            iflag = _plaus_ephemeris_flags(iflag)
            if ut:
                return _calc_nutation_obliquity(tjd, iflag)
            return _calc_nutation_obliquity_tt(tjd, iflag)

        # Same normalization preamble as the module-level calc_ut()/calc(),
        # so the context entry points stay 1:1 with the reference API
        # (FLG_MOSEPH strip, ephemeris-bit echo, FLG_SPEED3 mapping) on
        # every path, including the LEB fast path.
        raw_iflag = iflag
        iflag = _normalize_calc_flags(iflag)

        # Built-in asteroids by AST_OFFSET number (see module calc_ut)
        ipl = _remap_ast_offset(ipl)

        # South nodes: derive from the north node via this same context path,
        # mirroring the module-level entry points. The descending node must be
        # intercepted here because no downstream path derives it (_calc_body has
        # no antipode branch), and the antipode is representation-dependent
        # (FLG_XYZ / FLG_RADIANS) — _south_node_from_north handles all three.
        if ipl in (-MEAN_NODE, -TRUE_NODE):
            from .planets import _south_node_from_north

            north_result, retflag = self._calc_impl(tjd, abs(ipl), iflag, ut=ut)
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
                fast = fast_calc.fast_calc_ut if ut else fast_calc.fast_calc_tt
                result = fast(
                    reader,
                    tjd,
                    ipl,
                    iflag,
                    sid_mode=self.sidereal_mode,
                    sid_t0=self.sidereal_t0,
                    sid_ayan_t0=self.sidereal_ayan_t0,
                    topo=_ctx_topo,
                )
                from .logging_config import get_logger

                get_logger().debug("body=%d jd=%.1f source=LEB (context)", ipl, tjd)
                _record(ipl, "LEB")
                from .planets import _implied_retflag_bits, _echo_speed_bit

                from .planets import _to_native_floats

                return _to_native_floats(result[0]), _echo_speed_bit(
                    result[1] | _implied_retflag_bits(iflag), raw_iflag
                )
            except (KeyError, ValueError) as _leb_err:
                # missing body / out-of-range -> DEBUG, corruption -> WARNING
                from .leb_reader import log_leb_fallback

                log_leb_fallback(f"body={ipl} jd={tjd:.1f} (context)", _leb_err)
        # --- END LEB fast path ---

        # --- Horizons path: mirror the module-level dispatch. Without this,
        # set_calc_mode("horizons") silently computed via Skyfield on the
        # context API (or failed without a local DE440). Runs under the
        # context-state swap so per-context sidereal settings are honored.
        from .state import get_horizons_client

        h_client = get_horizons_client()
        if h_client is not None:
            from .logging_config import get_logger
            from . import horizons_backend

            from .planets import _swapped_context_state

            try:
                if ut:
                    jd_ut_approx = tjd
                else:
                    from .time_utils import deltat

                    # calc uses TT, convert to UT for horizons_calc_ut
                    jd_ut_approx = tjd - deltat(tjd)
                # horizons_calc_ut passes sid_mode explicitly, but the SIDM_USER
                # ayanamsha (reference epoch t0 and ayan_t0) is still read from
                # the module-level sidereal globals by get_ayanamsa downstream.
                # Run the call under _swapped_context_state — the same swap the
                # Skyfield path uses — so this context's per-instance t0/ayan_t0
                # (and mode) are honored, and the read happens under
                # _CONTEXT_SWAP_LOCK rather than unlocked against global state.
                with _swapped_context_state(self):
                    result = horizons_backend.horizons_calc_ut(
                        h_client, jd_ut_approx, ipl, iflag, self.sidereal_mode
                    )
                get_logger().debug(
                    "body=%d jd=%.1f source=Horizons (context)", ipl, tjd
                )
                _record(ipl, "Horizons")
                from .planets import _implied_retflag_bits, _echo_speed_bit

                from .planets import _to_native_floats

                return _to_native_floats(result[0]), _echo_speed_bit(
                    result[1] | _implied_retflag_bits(iflag), raw_iflag
                )
            except KeyError as _hz_err:
                get_logger().debug(
                    "body=%d jd=%.1f source=Horizons->fallback (context) reason=%s",
                    ipl,
                    tjd,
                    _hz_err,
                )
            except (ConnectionError, ValueError, OSError) as _hz_err:
                # Transient network/API failure, or a ValueError/OSError raised
                # from the HELCTR / velocity path — fall through to Skyfield
                # rather than let it escape the context calc.
                get_logger().warning(
                    "body=%d jd=%.1f source=Horizons->fallback (context, network): %s",
                    ipl,
                    tjd,
                    _hz_err,
                )
        # --- END Horizons path ---

        from skyfield.errors import EphemerisRangeError as SkyfieldRangeError

        from .exceptions import validate_jd_range
        from .planets import (
            _body_uses_jpl_ephemeris,
            _calc_body_with_context,
            _finalize_output_flags,
            _strip_output_flags,
            _wrap_ephemeris_range_error,
        )

        # Same out-of-range handling as the module-level entry points:
        # validate the JD first and wrap Skyfield's range error in the
        # library EphemerisRangeError so `except EphemerisRangeError`
        # works identically on the context API.
        if _body_uses_jpl_ephemeris(ipl):
            validate_jd_range(tjd, ipl, "calc_ut" if ut else "calc")

        ts = self.get_timescale()
        t = ts.ut1_jd(tjd) if ut else ts.tt_jd(tjd)
        # Output-format flags (FLG_XYZ / FLG_RADIANS) are handled here, like
        # the module-level entry points; _calc_body does not apply them.
        try:
            pos, retflag = _calc_body_with_context(
                t, ipl, _strip_output_flags(iflag), self
            )
        except SkyfieldRangeError as e:
            raise _wrap_ephemeris_range_error(e, tjd, ipl) from e
        except ValueError as e:
            # SPK type 21 kernels raise ValueError for out-of-coverage dates
            # (same conversion as the module-level entry points).
            if "Invalid Time" in str(e) or "time" in str(e).lower():
                raise _wrap_ephemeris_range_error(e, tjd, ipl) from e
            raise
        # Record dispatch source, like the LEB/Horizons branches above and
        # the module-level entry points (the context LEB/Horizons paths call
        # _record; the Skyfield fallback must too for telemetry parity).
        _record(ipl, "Skyfield")
        from .planets import _echo_speed_bit

        pos_out, rf_out = _finalize_output_flags(pos, retflag, iflag)
        return pos_out, _echo_speed_bit(rf_out, raw_iflag)

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

        Note:
            **Time scale — this method takes UT1.** The 1:1 reference-API
            parity is defined against the module-level ``calc_pctr()``,
            which takes TT (the reference API exposes only a TT
            planet-centric function). This context method is the
            deliberate planet-centric analogue of ``calc_ut()`` (UT1 in,
            like ``calc_ut``/``houses``), kept UT1 for backward
            compatibility; for TT input and strict reference-API parity,
            use the module-level ``libephemeris.calc_pctr()``.

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
        from skyfield.errors import EphemerisRangeError as SkyfieldRangeError

        from .exceptions import validate_jd_range
        from .planets import (
            _body_uses_jpl_ephemeris,
            _calc_body_pctr_with_context,
            _run_pctr_pipeline,
            _wrap_ephemeris_range_error,
        )

        # Same out-of-range handling as the module-level calc_pctr()
        if _body_uses_jpl_ephemeris(ipl) or _body_uses_jpl_ephemeris(iplctr):
            validate_jd_range(tjd_ut, ipl, "calc_pctr")

        ts = self.get_timescale()
        t = ts.ut1_jd(tjd_ut)
        # Same plaus_iflag/output-flag pipeline as the module-level
        # calc_pctr(), shared so the semantics cannot drift.
        try:
            return _run_pctr_pipeline(
                lambda calc_iflag: _calc_body_pctr_with_context(
                    t, ipl, iplctr, calc_iflag, self
                ),
                iflag,
            )
        except SkyfieldRangeError as e:
            raise _wrap_ephemeris_range_error(e, tjd_ut, ipl) from e

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

            # Context calculations route through the MODULE-level kernel
            # (state.get_planets() via _calc_body), not _SHARED_PLANETS —
            # which is only loaded by an explicit ctx.get_planets() call.
            # Honoring the documented contract (free file handles, pick up a
            # new ephemeris file on the next calculation) therefore requires
            # releasing the state-level kernel resources too. The narrow
            # helper closes only kernel handles/loader/caches and preserves
            # every user-facing setting. Lock order _SHARED_LOCK ->
            # state._INIT_LOCK matches get_timescale() above; state never
            # takes locks in the reverse order.
            from .state import _close_kernel_resources, _reset_timescale

            _close_kernel_resources()

            # _SHARED_TS caches the singleton built by state.get_timescale();
            # clearing the reference alone would re-cache the same live
            # object on the next call. Reset the state-level cache too so
            # time data (IERS table, leap seconds) is actually rebuilt.
            _reset_timescale()
