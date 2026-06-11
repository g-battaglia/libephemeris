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

from .constants import FLG_HELCTR, FLG_SPEED, FLG_SIDEREAL
from .state import get_loader, get_planets, get_timescale


# =============================================================================
# PLANETARY MOON CONSTANTS
# =============================================================================
# Body IDs following reference API 2.10+ convention
# Moon IDs start at MOON_OFFSET (9000) to avoid collision with other bodies

MOON_OFFSET: int = 9000

# Jupiter's Galilean Moons (discovered by Galileo in 1610)
MOON_IO: int = MOON_OFFSET + 1  # Jupiter I - innermost Galilean moon
MOON_EUROPA: int = MOON_OFFSET + 2  # Jupiter II - potential for life
MOON_GANYMEDE: int = MOON_OFFSET + 3  # Jupiter III - largest moon in solar system
MOON_CALLISTO: int = MOON_OFFSET + 4  # Jupiter IV - heavily cratered

# Saturn's Major Moons
MOON_MIMAS: int = MOON_OFFSET + 11  # Saturn I - "Death Star" moon
MOON_ENCELADUS: int = MOON_OFFSET + 12  # Saturn II - geysers, potential life
MOON_TETHYS: int = MOON_OFFSET + 13  # Saturn III - icy moon
MOON_DIONE: int = MOON_OFFSET + 14  # Saturn IV - icy moon
MOON_RHEA: int = MOON_OFFSET + 15  # Saturn V - second largest Saturn moon
MOON_TITAN: int = (
    MOON_OFFSET + 16
)  # Saturn VI - largest Saturn moon, thick atmosphere
MOON_HYPERION: int = MOON_OFFSET + 17  # Saturn VII - irregularly shaped
MOON_IAPETUS: int = MOON_OFFSET + 18  # Saturn VIII - two-toned coloring

# Uranus' Major Moons
MOON_MIRANDA: int = MOON_OFFSET + 21  # Uranus V - extreme geological features
MOON_ARIEL: int = MOON_OFFSET + 22  # Uranus I - brightest Uranian moon
MOON_UMBRIEL: int = MOON_OFFSET + 23  # Uranus II - darkest Uranian moon
MOON_TITANIA: int = MOON_OFFSET + 24  # Uranus III - largest Uranian moon
MOON_OBERON: int = MOON_OFFSET + 25  # Uranus IV - outermost major Uranian moon

# Neptune's Major Moon
MOON_TRITON: int = MOON_OFFSET + 31  # Neptune I - retrograde orbit, captured KBO

# Mars' Moons
MOON_PHOBOS: int = MOON_OFFSET + 41  # Mars I - larger, closer moon
MOON_DEIMOS: int = MOON_OFFSET + 42  # Mars II - smaller, farther moon

# Pluto's Moon
MOON_CHARON: int = (
    MOON_OFFSET + 51
)  # Pluto I - Pluto's largest moon (binary system)


# =============================================================================
# NAIF IDS FOR PLANETARY MOONS
# =============================================================================
# Standard NAIF SPICE IDs for planetary satellites

# Jupiter system (5xx)
NAIF_IO: int = 501
NAIF_EUROPA: int = 502
NAIF_GANYMEDE: int = 503
NAIF_CALLISTO: int = 504

# Saturn system (6xx)
NAIF_MIMAS: int = 601
NAIF_ENCELADUS: int = 602
NAIF_TETHYS: int = 603
NAIF_DIONE: int = 604
NAIF_RHEA: int = 605
NAIF_TITAN: int = 606
NAIF_HYPERION: int = 607
NAIF_IAPETUS: int = 608

# Uranus system (7xx)
NAIF_MIRANDA: int = 705
NAIF_ARIEL: int = 701
NAIF_UMBRIEL: int = 702
NAIF_TITANIA: int = 703
NAIF_OBERON: int = 704

# Neptune system (8xx)
NAIF_TRITON: int = 801

# Mars system (4xx)
NAIF_PHOBOS: int = 401
NAIF_DEIMOS: int = 402

# Pluto system (9xx)
NAIF_CHARON: int = 901

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
        # Use specified moons
        for moon_id in moons:
            if moon_id in MOON_NAIF_MAP:
                naif_id = MOON_NAIF_MAP[moon_id]
                # Verify NAIF ID exists in kernel
                try:
                    _ = kernel[naif_id]
                    _MOON_SPK_BY_BODY[moon_id] = spk_file
                except KeyError:
                    pass  # Skip moons not in this kernel
    else:
        # Auto-detect moons in kernel
        for moon_id, naif_id in MOON_NAIF_MAP.items():
            try:
                _ = kernel[naif_id]
                _MOON_SPK_BY_BODY[moon_id] = spk_file
            except KeyError:
                pass  # Moon not in this kernel

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


def get_moon_name(moon_id: int) -> str:
    """
    Get the human-readable name of a planetary moon.

    Args:
        moon_id: Moon ID (MOON_*)

    Returns:
        Moon name as string, or "Unknown Moon (ID)" if not recognized.

    Example:
        >>> get_moon_name(MOON_IO)
        'Io'
    """
    return MOON_NAMES.get(moon_id, f"Unknown Moon ({moon_id})")


def is_planetary_moon(ipl: int) -> bool:
    """
    Check if a body ID corresponds to a planetary moon.

    Args:
        ipl: Body ID to check

    Returns:
        True if the ID is a planetary moon ID (MOON_*), False otherwise.

    Example:
        >>> is_planetary_moon(MOON_TITAN)
        True
        >>> is_planetary_moon(SATURN)
        False
    """
    return ipl in MOON_NAIF_MAP


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

    Args:
        t: Skyfield Time object
        moon_id: Planetary moon ID (MOON_*)
        iflag: Calculation flags (FLG_SPEED, FLG_HELCTR, etc.)

    Returns:
        Position tuple (lon, lat, dist, speed_lon, speed_lat, speed_dist)
        or None if moon not registered.

    Raises:
        ValueError: If JD is outside SPK coverage
    """
    from .planets import get_ayanamsa_ut

    # Check if moon is registered
    if moon_id not in _MOON_SPK_BY_BODY:
        return None

    spk_file = _MOON_SPK_BY_BODY[moon_id]
    kernel = _MOON_SPK_KERNELS.get(spk_file)
    if kernel is None:
        return None

    naif_id = MOON_NAIF_MAP[moon_id]

    # Get moon target from kernel
    try:
        moon_target = kernel[naif_id]
    except KeyError:
        return None

    # Get main ephemeris for Earth/Sun
    planets = get_planets()

    # Determine observer
    is_heliocentric = bool(iflag & FLG_HELCTR)

    if is_heliocentric:
        observer = planets["sun"]
    else:
        observer = planets["earth"]

    # Skyfield chains SPK segments to the solar-system barycenter, so
    # kernel[naif_id].at(t) is already SSB -> moon (satellite kernels like
    # jup365/sat441 ship the SSB -> planet-barycenter segment; if a kernel
    # cannot chain, kernel[naif_id] above raises KeyError and we returned
    # None).  Adding the parent barycenter again would double-count it
    # (~5-30 AU error).
    moon_at_ssb = moon_target.at(t)
    observer_at_ssb = observer.at(t)

    moon_ssb = moon_at_ssb.position.au
    observer_xyz = observer_at_ssb.position.au

    # Moon position relative to observer
    rel_xyz = [moon_ssb[i] - observer_xyz[i] for i in range(3)]

    # Convert to ecliptic coordinates
    rel_pos = ICRF(rel_xyz, t=t, center=399)

    # Get ecliptic lat/lon
    ecl_pos = rel_pos.frame_latlon(ecliptic_frame)
    lat = ecl_pos[0].degrees
    lon = ecl_pos[1].degrees
    dist = ecl_pos[2].au

    # Normalize longitude to 0-360
    lon = lon % 360.0

    # Calculate speeds if requested
    speed_lon, speed_lat, speed_dist = 0.0, 0.0, 0.0

    if iflag & FLG_SPEED:
        # Central difference numerical differentiation (1 second timestep)
        ts = get_timescale()
        dt = 1.0 / 86400.0  # 1 second in days

        t_prev = ts.tt_jd(t.tt - dt)
        t_next = ts.tt_jd(t.tt + dt)

        # Position at t - dt
        moon_at_ssb_prev = moon_target.at(t_prev)
        observer_at_ssb_prev = observer.at(t_prev)

        moon_ssb_prev = moon_at_ssb_prev.position.au
        observer_xyz_prev = observer_at_ssb_prev.position.au

        rel_xyz_prev = [moon_ssb_prev[i] - observer_xyz_prev[i] for i in range(3)]

        rel_pos_prev = ICRF(rel_xyz_prev, t=t_prev, center=399)
        ecl_pos_prev = rel_pos_prev.frame_latlon(ecliptic_frame)

        lat_prev = ecl_pos_prev[0].degrees
        lon_prev = ecl_pos_prev[1].degrees % 360.0
        dist_prev = ecl_pos_prev[2].au

        # Position at t + dt
        moon_at_ssb_next = moon_target.at(t_next)
        observer_at_ssb_next = observer.at(t_next)

        moon_ssb_next = moon_at_ssb_next.position.au
        observer_xyz_next = observer_at_ssb_next.position.au

        rel_xyz_next = [moon_ssb_next[i] - observer_xyz_next[i] for i in range(3)]

        rel_pos_next = ICRF(rel_xyz_next, t=t_next, center=399)
        ecl_pos_next = rel_pos_next.frame_latlon(ecliptic_frame)

        lat_next = ecl_pos_next[0].degrees
        lon_next = ecl_pos_next[1].degrees % 360.0
        dist_next = ecl_pos_next[2].au

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
