# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Planetary moon position calculations for libephemeris.

This module provides support for calculating positions of planetary moons
(satellites) using JPL satellite SPK files (jup365.bsp, sat441.bsp, etc.).

Supported Moons:
- Jupiter's Galilean moons: Io, Europa, Ganymede, Callisto
- Saturn's major moons: Titan, Rhea, Dione, Tethys, Iapetus, Enceladus, Mimas
- Uranus' major moons: Miranda, Ariel, Umbriel, Titania, Oberon
- Neptune's major moon: Triton
- Mars' moons: Phobos, Deimos
- Pluto's moon: Charon

Usage:
    >>> import libephemeris as eph
    >>> # Register satellite SPK file
    >>> eph.register_moon_spk("jup365.bsp")
    >>> # Calculate Io's position
    >>> pos, _ = eph.calc_ut(2451545.0, eph.MOON_IO, eph.FLG_SPEED)
    >>> print(f"Io longitude: {pos[0]:.4f}")

SPK File Sources:
    - Jupiter satellites: jup365.bsp from JPL Horizons
    - Saturn satellites: sat441.bsp from JPL Horizons
    - Uranus satellites: ura116.bsp from JPL Horizons
    - Neptune satellites: nep097.bsp from JPL Horizons
    - Mars satellites: mar097.bsp from JPL Horizons
    - Pluto satellites: plu058.bsp from JPL Horizons

References:
    - JPL NAIF: https://naif.jpl.nasa.gov/naif/
    - Reference API 2.10+ moon support
"""

import math
import os
from typing import Any, Optional, Tuple

from skyfield.framelib import ecliptic_frame
from skyfield.positionlib import ICRF

from .constants import (
    FLG_BARYCTR,
    FLG_HELCTR,
    FLG_NOABERR,
    FLG_SIDEREAL,
    FLG_SPEED,
    FLG_TRUEPOS,
    # Legacy moon ids (deprecated aliases) and NAIF ids: constants.py is
    # the single source of truth; re-exported here for backward compat.
    MOON_OFFSET,
    MOON_IO,
    MOON_EUROPA,
    MOON_GANYMEDE,
    MOON_CALLISTO,
    MOON_MIMAS,
    MOON_ENCELADUS,
    MOON_TETHYS,
    MOON_DIONE,
    MOON_RHEA,
    MOON_TITAN,
    MOON_HYPERION,
    MOON_IAPETUS,
    MOON_MIRANDA,
    MOON_ARIEL,
    MOON_UMBRIEL,
    MOON_TITANIA,
    MOON_OBERON,
    MOON_TRITON,
    MOON_PHOBOS,
    MOON_DEIMOS,
    MOON_CHARON,
    NAIF_IO,
    NAIF_EUROPA,
    NAIF_GANYMEDE,
    NAIF_CALLISTO,
    NAIF_MIMAS,
    NAIF_ENCELADUS,
    NAIF_TETHYS,
    NAIF_DIONE,
    NAIF_RHEA,
    NAIF_TITAN,
    NAIF_HYPERION,
    NAIF_IAPETUS,
    NAIF_MIRANDA,
    NAIF_ARIEL,
    NAIF_UMBRIEL,
    NAIF_TITANIA,
    NAIF_OBERON,
    NAIF_TRITON,
    NAIF_PHOBOS,
    NAIF_DEIMOS,
    NAIF_CHARON,
    PLMOON_OFFSET,
)
from .exceptions import EphemerisRangeError
from .state import get_loader, get_planets, get_timescale


# Planet barycenters (for ephemeris lookup)
NAIF_MARS_BARYCENTER: int = 4
NAIF_JUPITER_BARYCENTER: int = 5
NAIF_SATURN_BARYCENTER: int = 6
NAIF_URANUS_BARYCENTER: int = 7
NAIF_NEPTUNE_BARYCENTER: int = 8
NAIF_PLUTO_BARYCENTER: int = 9


# =============================================================================
# MOON ID TO NAIF ID MAPPING
# =============================================================================

MOON_NAIF_MAP: dict[int, int] = {
    # Jupiter's Galilean moons
    MOON_IO: NAIF_IO,
    MOON_EUROPA: NAIF_EUROPA,
    MOON_GANYMEDE: NAIF_GANYMEDE,
    MOON_CALLISTO: NAIF_CALLISTO,
    # Saturn's major moons
    MOON_MIMAS: NAIF_MIMAS,
    MOON_ENCELADUS: NAIF_ENCELADUS,
    MOON_TETHYS: NAIF_TETHYS,
    MOON_DIONE: NAIF_DIONE,
    MOON_RHEA: NAIF_RHEA,
    MOON_TITAN: NAIF_TITAN,
    MOON_HYPERION: NAIF_HYPERION,
    MOON_IAPETUS: NAIF_IAPETUS,
    # Uranus' major moons
    MOON_MIRANDA: NAIF_MIRANDA,
    MOON_ARIEL: NAIF_ARIEL,
    MOON_UMBRIEL: NAIF_UMBRIEL,
    MOON_TITANIA: NAIF_TITANIA,
    MOON_OBERON: NAIF_OBERON,
    # Neptune's major moon
    MOON_TRITON: NAIF_TRITON,
    # Mars' moons
    MOON_PHOBOS: NAIF_PHOBOS,
    MOON_DEIMOS: NAIF_DEIMOS,
    # Pluto's moon
    MOON_CHARON: NAIF_CHARON,
}

# Moon ID to human-readable name mapping
MOON_NAMES: dict[int, str] = {
    MOON_IO: "Io",
    MOON_EUROPA: "Europa",
    MOON_GANYMEDE: "Ganymede",
    MOON_CALLISTO: "Callisto",
    MOON_MIMAS: "Mimas",
    MOON_ENCELADUS: "Enceladus",
    MOON_TETHYS: "Tethys",
    MOON_DIONE: "Dione",
    MOON_RHEA: "Rhea",
    MOON_TITAN: "Titan",
    MOON_HYPERION: "Hyperion",
    MOON_IAPETUS: "Iapetus",
    MOON_MIRANDA: "Miranda",
    MOON_ARIEL: "Ariel",
    MOON_UMBRIEL: "Umbriel",
    MOON_TITANIA: "Titania",
    MOON_OBERON: "Oberon",
    MOON_TRITON: "Triton",
    MOON_PHOBOS: "Phobos",
    MOON_DEIMOS: "Deimos",
    MOON_CHARON: "Charon",
}

# Moon ID to parent planet barycenter NAIF ID mapping
MOON_PARENT_MAP: dict[int, int] = {
    # Jupiter's moons
    MOON_IO: NAIF_JUPITER_BARYCENTER,
    MOON_EUROPA: NAIF_JUPITER_BARYCENTER,
    MOON_GANYMEDE: NAIF_JUPITER_BARYCENTER,
    MOON_CALLISTO: NAIF_JUPITER_BARYCENTER,
    # Saturn's moons
    MOON_MIMAS: NAIF_SATURN_BARYCENTER,
    MOON_ENCELADUS: NAIF_SATURN_BARYCENTER,
    MOON_TETHYS: NAIF_SATURN_BARYCENTER,
    MOON_DIONE: NAIF_SATURN_BARYCENTER,
    MOON_RHEA: NAIF_SATURN_BARYCENTER,
    MOON_TITAN: NAIF_SATURN_BARYCENTER,
    MOON_HYPERION: NAIF_SATURN_BARYCENTER,
    MOON_IAPETUS: NAIF_SATURN_BARYCENTER,
    # Uranus' moons
    MOON_MIRANDA: NAIF_URANUS_BARYCENTER,
    MOON_ARIEL: NAIF_URANUS_BARYCENTER,
    MOON_UMBRIEL: NAIF_URANUS_BARYCENTER,
    MOON_TITANIA: NAIF_URANUS_BARYCENTER,
    MOON_OBERON: NAIF_URANUS_BARYCENTER,
    # Neptune's moons
    MOON_TRITON: NAIF_NEPTUNE_BARYCENTER,
    # Mars' moons
    MOON_PHOBOS: NAIF_MARS_BARYCENTER,
    MOON_DEIMOS: NAIF_MARS_BARYCENTER,
    # Pluto's moons
    MOON_CHARON: NAIF_PLUTO_BARYCENTER,
}


# =============================================================================
# MODULE STATE
# =============================================================================
# Registered satellite SPK kernels

_MOON_SPK_KERNELS: dict[str, Any] = {}  # {filepath: SpiceKernel}
_MOON_SPK_BY_BODY: dict[int, str] = {}  # {moon_id: filepath}


# =============================================================================
# SATELLITE SPK REGISTRATION
# =============================================================================


def register_moon_spk(
    spk_file: str,
    moons: Optional[list[int]] = None,
) -> None:
    """
    Register a satellite SPK kernel file for planetary moon calculations.

    After registration, calc_ut() will automatically use the SPK kernel
    for the moons contained in the file.

    Args:
        spk_file: Path to the satellite SPK file, or filename if in library path.
                  Common files: jup365.bsp, sat441.bsp, ura116.bsp, nep097.bsp
        moons: Optional list of specific moon IDs to register from this file.
               If None, auto-detects available moons based on NAIF IDs in the kernel.

    Raises:
        FileNotFoundError: If SPK file not found
        ValueError: If no valid moons found in the SPK kernel

    Example:
        >>> register_moon_spk("jup365.bsp")  # Auto-detect Jupiter moons
        >>> register_moon_spk("sat441.bsp", [MOON_TITAN, MOON_ENCELADUS])
    """
    from . import state

    # Resolve file path
    if not os.path.isabs(spk_file):
        data_dir = state._get_data_dir()
        full_path = os.path.join(data_dir, spk_file)
        if os.path.exists(full_path):
            spk_file = full_path
        elif not os.path.exists(spk_file):
            raise FileNotFoundError(f"Satellite SPK file not found: {spk_file}")

    if not os.path.exists(spk_file):
        raise FileNotFoundError(f"Satellite SPK file not found: {spk_file}")

    # Load kernel
    load = get_loader()
    kernel = load(spk_file)
    _MOON_SPK_KERNELS[spk_file] = kernel

    # Determine which moons to register
    if moons is not None:
        # Use specified moons (legacy or canonical ids)
        for moon_id in moons:
            naif_id = _moon_naif_id(moon_id)
            if naif_id is None:
                continue
            # Verify NAIF ID exists in kernel
            try:
                _ = kernel[naif_id]
            except KeyError:
                continue  # Skip moons not in this kernel
            _MOON_SPK_BY_BODY[moon_id] = spk_file
            _MOON_SPK_BY_BODY[PLMOON_OFFSET + naif_id] = spk_file
    else:
        # Auto-detect: register every satellite target the kernel can
        # chain to the SSB, under its canonical id (PLMOON_OFFSET + NAIF);
        # the 21 named moons get their legacy alias ids as well.
        naif_targets: set[int] = set()
        try:
            for seg in kernel.spk.segments:
                tgt = int(seg.target)
                if 400 < tgt < 999:
                    naif_targets.add(tgt)
        except (AttributeError, TypeError):
            naif_targets = set(MOON_NAIF_MAP.values())
        legacy_by_naif = {n: m for m, n in MOON_NAIF_MAP.items()}
        for naif_id in sorted(naif_targets):
            try:
                _ = kernel[naif_id]
            except KeyError:
                continue  # Cannot chain to SSB with loaded kernels
            _MOON_SPK_BY_BODY[PLMOON_OFFSET + naif_id] = spk_file
            if naif_id in legacy_by_naif:
                _MOON_SPK_BY_BODY[legacy_by_naif[naif_id]] = spk_file

    if not any(spk_file == path for path in _MOON_SPK_BY_BODY.values()):
        raise ValueError(f"No valid planetary moons found in SPK kernel: {spk_file}")


def unregister_moon_spk(spk_file: str) -> None:
    """
    Unregister a satellite SPK kernel and its associated moons.

    Args:
        spk_file: Path to the SPK file to unregister

    Example:
        >>> unregister_moon_spk("jup365.bsp")
    """
    # Remove body registrations for this file
    moons_to_remove = [
        moon_id
        for moon_id, path in _MOON_SPK_BY_BODY.items()
        if path == spk_file or path.endswith(os.path.basename(spk_file))
    ]
    for moon_id in moons_to_remove:
        del _MOON_SPK_BY_BODY[moon_id]

    # Remove kernel from cache
    if spk_file in _MOON_SPK_KERNELS:
        try:
            _MOON_SPK_KERNELS[spk_file].close()
        except (AttributeError, OSError, ValueError):
            pass
        del _MOON_SPK_KERNELS[spk_file]


def list_registered_moons() -> dict[int, str]:
    """
    List all registered planetary moons.

    Returns:
        Dict mapping moon_id -> spk_file path for all registered moons.

    Example:
        >>> register_moon_spk("jup365.bsp")
        >>> moons = list_registered_moons()
        >>> for moon_id, path in moons.items():
        ...     print(f"{get_moon_name(moon_id)}: {path}")
    """
    return dict(_MOON_SPK_BY_BODY)


# Reverse map for naming canonical ids: NAIF id -> moon name
_NAIF_NAMES: dict[int, str] = {
    MOON_NAIF_MAP[mid]: name for mid, name in MOON_NAMES.items()
}


def _moon_naif_id(ipl: int) -> Optional[int]:
    """NAIF satellite id for a moon body id, or None if not a moon id.

    Accepts both the canonical scheme (PLMOON_OFFSET + NAIF id, e.g.
    Io = 9501) and the legacy deprecated aliases (MOON_IO = 9001 etc.).
    The canonical range covers every NAIF satellite id 401-998 plus the
    9n99 planet-center ids, so moons without a named constant (e.g.
    Amalthea = 9505) work as long as a registered kernel contains them.
    """
    legacy = MOON_NAIF_MAP.get(ipl)
    if legacy is not None:
        return legacy
    if PLMOON_OFFSET + 400 < ipl < PLMOON_OFFSET + 1000:
        return ipl - PLMOON_OFFSET
    return None


def get_moon_name(moon_id: int) -> str:
    """
    Get the human-readable name of a planetary moon.

    Args:
        moon_id: Moon ID (canonical PLMOON_OFFSET + NAIF id, or legacy MOON_*)

    Returns:
        Moon name as string, or "Unknown Moon (ID)" if not recognized.

    Example:
        >>> get_moon_name(MOON_IO)
        'Io'
        >>> get_moon_name(9501)
        'Io'
    """
    name = MOON_NAMES.get(moon_id)
    if name is not None:
        return name
    naif = _moon_naif_id(moon_id)
    if naif is not None and naif in _NAIF_NAMES:
        return _NAIF_NAMES[naif]
    return f"Unknown Moon ({moon_id})"


def is_planetary_moon(ipl: int) -> bool:
    """
    Check if a body ID corresponds to a planetary moon.

    Args:
        ipl: Body ID to check

    Returns:
        True for canonical moon ids (PLMOON_OFFSET + NAIF satellite id,
        e.g. 9501 for Io) and for the legacy MOON_* aliases.

    Example:
        >>> is_planetary_moon(MOON_TITAN)
        True
        >>> is_planetary_moon(9606)
        True
        >>> is_planetary_moon(SATURN)
        False
    """
    return _moon_naif_id(ipl) is not None


# =============================================================================
# POSITION CALCULATION
# =============================================================================


def calc_moon_position(
    t,
    moon_id: int,
    iflag: int,
) -> Optional[Tuple[float, float, float, float, float, float]]:
    """
    Calculate planetary moon position using satellite SPK kernel.

    Internal function called by planets._calc_body() when a moon is requested.

    Geocentric output is apparent: light-time iteration on the satellite
    (unless FLG_TRUEPOS) and annual aberration (unless FLG_NOABERR).
    Heliocentric and barycentric output stays geometric, matching the
    library's other helio/bary paths.

    Args:
        t: Skyfield Time object
        moon_id: Moon ID (canonical PLMOON_OFFSET + NAIF id, or legacy MOON_*)
        iflag: Calculation flags (FLG_SPEED, FLG_HELCTR, FLG_BARYCTR, etc.)

    Returns:
        Position tuple (lon, lat, dist, speed_lon, speed_lat, speed_dist)
        or None if moon not registered.

    Raises:
        EphemerisRangeError: If JD is outside the SPK kernel's coverage
    """
    from .planets import get_ayanamsa_ut

    naif_id = _moon_naif_id(moon_id)
    if naif_id is None:
        return None

    # Find the registered kernel: by body id first, then by NAIF target
    # scan (covers canonical ids the registration didn't key explicitly).
    spk_file = _MOON_SPK_BY_BODY.get(moon_id) or _MOON_SPK_BY_BODY.get(
        PLMOON_OFFSET + naif_id
    )
    kernel = _MOON_SPK_KERNELS.get(spk_file) if spk_file else None
    moon_target = None
    if kernel is not None:
        try:
            moon_target = kernel[naif_id]
        except KeyError:
            moon_target = None
    if moon_target is None:
        for cand in _MOON_SPK_KERNELS.values():
            try:
                moon_target = cand[naif_id]
                break
            except KeyError:
                continue
    if moon_target is None:
        return None

    # Get main ephemeris for Earth/Sun
    planets = get_planets()

    # Determine observer
    is_heliocentric = bool(iflag & FLG_HELCTR)
    is_barycentric = bool(iflag & FLG_BARYCTR)

    if is_barycentric:
        observer = None  # solar-system barycenter: zero position/velocity
    elif is_heliocentric:
        observer = planets["sun"]
    else:
        observer = planets["earth"]

    # Skyfield chains SPK segments to the solar-system barycenter, so
    # kernel[naif_id].at(t) is already SSB -> moon (satellite kernels like
    # jup365/sat441 ship the SSB -> planet-barycenter segment; if a kernel
    # cannot chain, kernel[naif_id] above raises KeyError).  Adding the
    # parent barycenter again would double-count it (~5-30 AU error).
    from .astrometry import apply_aberration_to_position

    C_AU_DAY = 173.144632674240
    apply_light_time = (
        not (iflag & FLG_TRUEPOS) and not is_heliocentric and not is_barycentric
    )
    apply_aberr = (
        not (iflag & FLG_NOABERR) and not is_heliocentric and not is_barycentric
    )
    ts = get_timescale()

    def _rel_at(t_obs):
        """Light-time corrected satellite-minus-observer vector (AU)."""
        try:
            if observer is None:
                obs_at = None
                obs_xyz = (0.0, 0.0, 0.0)
            else:
                obs_at = observer.at(t_obs)
                obs_xyz = obs_at.position.au
            t_emit = t_obs
            rel = None
            for _ in range(3):
                tgt_xyz = moon_target.at(t_emit).position.au
                rel = [tgt_xyz[i] - obs_xyz[i] for i in range(3)]
                if not apply_light_time:
                    break
                dist_au = math.sqrt(rel[0] ** 2 + rel[1] ** 2 + rel[2] ** 2)
                t_emit = ts.tt_jd(t_obs.tt - dist_au / C_AU_DAY)
        except ValueError as e:
            raise EphemerisRangeError(
                f"JD {t_obs.tt:.1f} outside satellite SPK coverage for "
                f"{get_moon_name(moon_id)}: {e}"
            ) from e
        if apply_aberr and obs_at is not None:
            obs_vel = obs_at.velocity.au_per_d
            rel = list(
                apply_aberration_to_position(
                    (rel[0], rel[1], rel[2]),
                    (obs_vel[0], obs_vel[1], obs_vel[2]),
                )
            )
        return rel

    def _ecl_at(t_obs):
        rel = _rel_at(t_obs)
        rel_pos = ICRF(rel, t=t_obs, center=399)
        ecl = rel_pos.frame_latlon(ecliptic_frame)
        return (
            float(ecl[1].degrees) % 360.0,
            float(ecl[0].degrees),
            float(ecl[2].au),
        )

    lon, lat, dist = _ecl_at(t)

    # Calculate speeds if requested
    speed_lon, speed_lat, speed_dist = 0.0, 0.0, 0.0

    if iflag & FLG_SPEED:
        # Central difference numerical differentiation (1 second timestep)
        dt = 1.0 / 86400.0  # 1 second in days

        lon_prev, lat_prev, dist_prev = _ecl_at(ts.tt_jd(t.tt - dt))
        lon_next, lat_next, dist_next = _ecl_at(ts.tt_jd(t.tt + dt))

        # Compute rates using central difference (per day)
        speed_lon = (lon_next - lon_prev) / (2.0 * dt)
        speed_lat = (lat_next - lat_prev) / (2.0 * dt)
        speed_dist = (dist_next - dist_prev) / (2.0 * dt)

        # Handle 360 wrap
        if speed_lon > 180.0 / (2.0 * dt):
            speed_lon -= 360.0 / (2.0 * dt)
        if speed_lon < -180.0 / (2.0 * dt):
            speed_lon += 360.0 / (2.0 * dt)

    # Apply sidereal correction if requested
    if iflag & FLG_SIDEREAL:
        ayanamsa = get_ayanamsa_ut(t.ut1)
        lon = (lon - ayanamsa) % 360.0

    return (lon, lat, dist, speed_lon, speed_lat, speed_dist)


def get_moon_coverage(moon_id: int) -> Optional[Tuple[float, float]]:
    """
    Get the time coverage for a registered planetary moon.

    Args:
        moon_id: Planetary moon ID (MOON_*)

    Returns:
        Tuple of (start_jd, end_jd) Julian Day coverage, or None if not registered.

    Example:
        >>> register_moon_spk("jup365.bsp")
        >>> start, end = get_moon_coverage(MOON_IO)
        >>> print(f"Io coverage: JD {start:.1f} to {end:.1f}")
    """
    if moon_id not in _MOON_SPK_BY_BODY:
        return None

    spk_file = _MOON_SPK_BY_BODY[moon_id]
    kernel = _MOON_SPK_KERNELS.get(spk_file)
    if kernel is None:
        return None

    try:
        if hasattr(kernel, "spk") and hasattr(kernel.spk, "segments"):
            segments = list(kernel.spk.segments)
            if segments:
                start_jd = min(
                    float(s.start_jd) for s in segments if hasattr(s, "start_jd")
                )
                end_jd = max(float(s.end_jd) for s in segments if hasattr(s, "end_jd"))
                return (start_jd, end_jd)
    except (ValueError, KeyError, AttributeError):
        pass

    return None


# =============================================================================
# CLEANUP
# =============================================================================


def close_moon_kernels() -> None:
    """
    Close all loaded satellite SPK kernels and clear registrations.

    Called automatically by state.close(), but can be called manually
    to free resources.
    """
    global _MOON_SPK_KERNELS, _MOON_SPK_BY_BODY

    for kernel in _MOON_SPK_KERNELS.values():
        try:
            kernel.close()
        except (AttributeError, OSError, ValueError):
            pass

    _MOON_SPK_KERNELS = {}
    _MOON_SPK_BY_BODY = {}
