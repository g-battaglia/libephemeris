# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Astrological house system calculations for libephemeris.

Implements 25 distinct house systems across 26 character codes (A as alias of E):
- Placidus (P): Most common, time-based, fails at polar latitudes
- Koch (K): Birthplace system, similar to Placidus
- Porphyry (O): Space-based trisection
- Regiomontanus (R): Medieval rational system
- Campanus (C): Prime vertical system
- Equal (A/E): Equal 30° divisions from Ascendant
- Whole Sign (W): Whole zodiac signs from Ascendant sign
- Meridian (X): Equatorial meridian divisions
- Azimuthal/Horizontal (H): Based on horizon
- Polich-Page/Topocentric (T)
- Alcabitius (B): Ancient Arabic system
- Morinus (M): Equatorial divisions
- Krusinski-Pisa-Goelzer (U)
- Gauquelin (G): 36 sectors
- Vehlow (V): Equal from midpoint
- APC (Y): Astronomical Planetary Cusps
- Carter Poli-Equatorial (F)
- Pullen SD (L): Sinusoidal Delta / Neo-Porphyry
- Pullen SR (Q): Sinusoidal Ratio
- Sripati (S): Divide quadrants equally
- Natural Gradient (N): Equal from 1=Aries
- Equal from MC (D)
- Sunshine/Treindl (I)
- Sunshine/Makransky (i)
- Savard-A (J)

Main Functions:
- houses(): Calculate house cusps and angles (ASCMC)
- houses_ex(): Extended version with sidereal support
- house_pos(): Find which house a point is in
- house_name(): Get house system name

Precision Notes:
Core Calculations (Phase 1 improved):
- GAST (sidereal time): ~0.001 sec = ~0.015 arcsec (Skyfield IAU SOFA)
- Obliquity: Vondrák 2011 long-term mean obliquity + IAU 2006/2000A nutation
- Vertex: Exact at equator, ~0.001° elsewhere (rigorous limiting formula)
- Iterative systems (Placidus/Koch): ~0.00036 arcsec convergence (1e-7°)
- Non-iterative systems: Limited by obliquity/GAST precision (~0.01 arcsec)

Expected Accuracy by System:
- Simple (Equal, Whole Sign, Vehlow): Exact (no astronomical calculations)
- Geometric (Porphyry, Morinus, Meridian): ~0.01-0.1 arcsec
- Iterative (Placidus, Koch): ~0.001-0.01° (a few arcseconds)
- Complex (Campanus, Regiomontanus, Topocentric): ~0.01°
- Horizontal: ~0.01° (with convergence fallback to Porphyry)

Polar Latitudes:
- Placidus, Koch undefined > ~66° latitude (circumpolar ecliptic points)
- Automatic fallback to Porphyry when iteration fails
- Equal/Whole Sign work at all latitudes

Source boundary:
- Meeus, chapter 13, supplies general coordinate transformations, not the
  definitions of the individual house systems.
- Each system's defining construction and bibliographic status are recorded
  separately in ``docs/reference/house-systems.md``.
- IERS/ERFA and Vondrak supply Earth orientation and the equator/ecliptic of
  date; they do not define astrological house divisions.

Provenance:
    Per-system definitions and source status are enumerated in
    ``docs/reference/house-systems.md``. This module independently derives the
    cusp geometry using spherical trigonometry or Cartesian great-circle
    intersections; ARMC, obliquity, precession, and nutation come from the
    registered Vondrak/IAU/IERS chain. Iteration tolerances, centered speed
    differences, polar-domain checks, and fallback policy are project-authored
    numerical/API choices and are documented beside their use. No cusp table or
    fitted output from another implementation is an input.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING, Any, Callable, List, Optional, Tuple, Union, overload
from .constants import *
from .constants import (
    FLG_RADIANS,
    FLG_SIDEREAL,
    FLG_SPEED,
    FLG_EQUATORIAL,
    SUN,
    FLG_JPLEPH,
    FLG_SWIEPH,
    FLG_TOPOCTR,
    SIDBIT_ECL_T0,
    SIDBIT_SSY_PLANE,
    BIT_DISC_CENTER,
    BIT_NO_REFRACTION,
    CALC_RISE,
    CALC_SET,
)
from .compat_names import HOUSE_SYSTEM_NAMES
from .planets import calc_ut
from .cache import get_cached_nutation
from .exceptions import (
    CalculationError,
    EphemerisRangeError,
    Error,
    PolarCircleError,
    validate_coordinates,
    validate_latitude,
)
from .house_constructions import (
    apc_cusp as _apc_cusp,
    houses_krusinski as _houses_krusinski,
    houses_savard_a as _houses_savard_a,
    houses_sunshine_makransky as _houses_sunshine_makransky,
)
from .utils import cotrans, degnorm, difdeg2n
from . import sidereal_longterm as _sidlt
from .time_utils import deltat as _deltat

if TYPE_CHECKING:
    import numpy as np


# Systems whose cusps are pinned to exact 30-degree sign boundaries, so the
# cusp function is piecewise constant in time: Whole Sign ('W') and Aries
# ('N'). Their true cusp derivative is zero, and the only variation across a
# finite-difference stencil is the 30-degree step taken when the Ascendant
# changes sign.
_SIGN_PINNED_HSYS = frozenset({"W", "N"})


# An int outside the character range is an unknown house-system selector like
# any other. Letting chr() raise made -1 fail with an untyped ValueError while
# 999 - equally unmapped - fell through to the documented default. U+FFFD is
# never a house-system code, so it reaches that same default.
_UNKNOWN_HSYS_CHAR = "\ufffd"


def _int_hsys_to_char(code: int) -> str:
    """``chr()`` for a house-system selector given as an integer code.

    An int outside the range ``chr()`` accepts is an unknown selector like
    any other and folds to ``_UNKNOWN_HSYS_CHAR``, so it reaches the same
    default as every other unmapped value instead of leaking ``chr()``'s
    ValueError. Every int-to-char conversion of a selector goes through
    here; the reference answers an unmapped selector with the default.
    """
    try:
        return chr(code)
    except (ValueError, OverflowError):
        return _UNKNOWN_HSYS_CHAR


def _fold_hsys_case(hsys_char: str) -> str:
    """Fold a house-system selector to its canonical case.

    Compatibility contract: house-system letters are accepted
    case-insensitively (``'k'`` selects Koch exactly like ``'K'``), with two
    deliberate exceptions: lowercase ``'i'`` is a distinct system of its own
    (Sunshine/Makransky alternative), not an alias of ``'I'`` (Sunshine),
    and lowercase ``'g'`` is NOT an alias of Gauquelin ``'G'`` — it stays an
    unrecognized selector served by the default 12-cusp fallback (folding it
    would change the return shape from 12 to 36 cusps).

    Only ASCII ``a``-``z`` fold. ``str.upper()`` is Unicode-aware and can
    both change length and cross into the ASCII selector space: byte 0xDF is
    U+00DF, whose uppercase is the two-character ``"SS"``, so an unknown
    one-byte selector was silently resolving to Sripati (``'S'``) through
    the later first-character lookup instead of reaching the
    unknown-selector fallback.

    The selector is one byte, so the ``bytes`` entry points decode it as
    latin-1 (byte value = code point) in :func:`_hsys_to_char`. Decoding as
    UTF-8 instead made every high byte raise ``UnicodeDecodeError`` on calls
    the reference answers with the ordinary unknown-selector fallback.
    """
    # A selector is ONE character. The bytes/str entry points can hand a
    # longer value straight here; a multi-character selector raises
    # (compatibility contract: TypeError, argument 4 must be a byte string of
    # length 1), and the table lookup below is only meaningful for one
    # character anyway.
    if len(hsys_char) != 1:
        raise TypeError(
            f"house-system identifier must be a single character, not {hsys_char!r}"
        )
    return _HSYS_CASE_FOLD.get(hsys_char, hsys_char)


# The case-folding table behind _fold_hsys_case: the ASCII lowercase letters
# that select like their uppercase form. 'i' and 'g' are absent on purpose
# (see the docstring); anything not listed is returned unchanged.
_HSYS_CASE_FOLD: dict[str, str] = {
    lower: lower.upper()
    for lower in "abcdefghijklmnopqrstuvwxyz"
    if lower not in ("i", "g")
}


def _hsys_to_char(hsys: int | bytes | str) -> str:
    """Normalize a house-system selector (int code, bytes, or str) to a char.

    This is the single entry for the three selector forms the public house
    API accepts: an ``int`` code (``ord('P')``), a one-byte ``bytes``
    (``b'P'``, the reference-API default form) or a one-character ``str``.
    An integer outside the Unicode range is an unknown selector like any
    other and folds to ``_UNKNOWN_HSYS_CHAR`` (see ``_int_hsys_to_char``).
    Letting ``chr()`` raise made one class of unknown selector fail with an
    untyped ValueError while every other unmapped value reached the
    documented default.
    """
    if isinstance(hsys, int):
        return _fold_hsys_case(_int_hsys_to_char(hsys))
    if isinstance(hsys, bytes):
        return _fold_hsys_case(hsys.decode("latin-1"))
    return _fold_hsys_case(str(hsys))


@dataclass(frozen=True)
class _HouseFrame:
    """The quantities a cusp construction reads, all in degrees.

    ``armc`` and ``lat`` describe the frame the construction builds on: for
    the systems that follow the midheaven above the horizon
    (:attr:`_McBranch.ROTATE_FRAME`) they are the rotated ARMC and the
    negated latitude, otherwise the original ones. ``asc`` and ``mc`` are
    the Ascendant and the midheaven already reported in ``ascmc``;
    ``sun_dec`` is the Sun's declination, read by the Sunshine systems only.
    """

    armc: float
    lat: float
    eps: float
    asc: float
    mc: float
    sun_dec: float


class _McBranch(Enum):
    """How a house system chooses the branch of the meridian it reports as MC.

    ``UPPER``: the midheaven is the upper culmination of the given ARMC.
    ``ROTATE_FRAME``: when that midheaven is below the horizon the whole
    frame is rotated half a turn (ARMC + 180, latitude negated) before the
    midheaven and the intermediate cusps are built (Regiomontanus,
    Campanus, Polich/Page, Savard-A). ``ROTATE_REPORTED``: the cusps are
    built on the given frame and only the reported MC, with cusps 10 and 4,
    is moved to the branch above the horizon (APC). ``FROM_CUSPS``: the
    construction selects the branch itself and the reported MC is its cusp
    10 (the Sunshine systems).
    """

    UPPER = "upper"
    ROTATE_FRAME = "rotate-frame"
    ROTATE_REPORTED = "rotate-reported"
    FROM_CUSPS = "from-cusps"


@dataclass(frozen=True)
class HouseSystem:
    """One row of the house-system registry consulted by the dispatcher.

    Attributes:
        code: The canonical selector character.
        build_cusps: The construction; it receives a :class:`_HouseFrame`
            and returns the 13-slot (37 for Gauquelin) cusp list with slot 0
            unused.
        mc_branch: Which branch of the meridian is reported as MC.
        raises_at_pole: Whether the system is undefined inside the polar
            circle and raises :class:`PolarCircleError` there (Placidus,
            Koch, Gauquelin, and every selector the registry does not list,
            which is served by the Placidus default).
        sectors: How many cusps the public tuple carries: 12, or 36 for the
            uppercase Gauquelin selector only.
        reads_sun_declination: Whether the construction reads the Sun's
            declination (the Sunshine systems), so :func:`houses` has to
            compute it.
    """

    code: str
    build_cusps: Callable[[_HouseFrame], List[float]]
    mc_branch: _McBranch = _McBranch.UPPER
    raises_at_pole: bool = False
    sectors: int = 12
    reads_sun_declination: bool = False


def _house_armc_obliquity(tjdut: float) -> tuple[float, float]:
    """Sidereal time at longitude 0 (= ARMC) and true obliquity, for houses.

    Both quantities come from the long-term model in :mod:`sidereal_longterm`
    (Vondrák 2011 precession/obliquity + the geometric sidereal-time method),
    so house cusps stay correct over the full supported date range instead of
    diverging at remote epochs the way an IAU-2006 sidereal-time polynomial
    would. The of-date mean obliquity is the same realization used by the
    position pipeline, keeping cusps and bodies in one self-consistent frame.

    Args:
        tjdut: Julian Day in UT1.

    Returns:
        (armc0_deg, eps_true_deg): apparent sidereal time at longitude 0 and
        the true obliquity of the ecliptic, both in degrees.
    """
    jd_tt = tjdut + _deltat(tjdut)
    dpsi_rad, deps_rad = get_cached_nutation(jd_tt)
    eps_mean = _sidlt.mean_obliquity_deg(jd_tt)
    eps_true = eps_mean + math.degrees(deps_rad)
    armc0 = _sidlt.apparent_sidereal_time_deg(
        tjdut, 0.0, dpsi_deg=math.degrees(dpsi_rad), eps_true_deg=eps_true
    )
    return armc0, eps_true


def _init_cardinal_cusps(asc: float, mc: float) -> list:
    """Create 13-element cusp list with cardinal points and their opposites."""
    cusps = [0.0] * 13
    cusps[1] = asc
    cusps[10] = mc
    cusps[7] = (asc + 180) % 360.0
    cusps[4] = (mc + 180) % 360.0
    return cusps


# Intermediate cusps and their diametrically opposite houses: a quadrant
# system constructs 11, 12, 2 and 3 and mirrors them into 5, 6, 8 and 9.
_OPPOSITE_INTERMEDIATE_CUSPS: tuple[tuple[int, int], ...] = (
    (11, 5),
    (12, 6),
    (2, 8),
    (3, 9),
)


def _set_opposite_cusps(cusps: list) -> None:
    """Fill cusps 5, 6, 8 and 9 as the opposites of 11, 12, 2 and 3."""
    for built, opposite in _OPPOSITE_INTERMEDIATE_CUSPS:
        cusps[opposite] = (cusps[built] + 180.0) % 360.0


def _is_polar_circle(lat: float, eps: float) -> bool:
    """
    Check if latitude is within the polar circle for house calculations.

    At polar latitudes (approximately >66.5°), some house systems like Placidus
    and Koch cannot be calculated because the ecliptic does not properly
    intersect the horizon. This occurs when abs(lat) + eps > 90°.

    Args:
        lat: Geographic latitude in degrees
        eps: True obliquity of the ecliptic in degrees

    Returns:
        True if within polar circle (house calculation will fail)
    """
    # The polar circle is where the sun can be circumpolar
    # This happens when |lat| + obliquity > 90°
    # For typical obliquity of ~23.44°, this is lat > ~66.56°
    return abs(lat) + eps > 90.0


def get_polar_latitude_threshold(obliquity: float = 23.44) -> float:
    """
    Calculate the polar latitude threshold for a given obliquity.

    At latitudes beyond this threshold, some house systems (Placidus, Koch,
    Gauquelin) cannot be calculated because ecliptic points can be circumpolar
    (never rising or setting).

    The threshold is calculated as: 90° - obliquity

    For the current epoch (J2000), obliquity is approximately 23.44°,
    giving a threshold of approximately 66.56°.

    Args:
        obliquity: True obliquity of the ecliptic in degrees.
                   Default is 23.44° (approximate J2000 value).

    Returns:
        The polar latitude threshold in degrees. Latitudes with
        abs(lat) > threshold will trigger polar circle errors.

    Example:
        >>> get_polar_latitude_threshold()
        66.56
        >>> get_polar_latitude_threshold(23.5)
        66.5
    """
    return 90.0 - obliquity


def _get_polar_circle_info(lat: float, eps: float, house_system: str) -> dict:
    """
    Get detailed information about a polar circle condition.

    Args:
        lat: Geographic latitude in degrees
        eps: True obliquity of the ecliptic in degrees
        house_system: House system character (e.g., 'P', 'K', 'G')

    Returns:
        Dictionary with polar circle information:
        - is_polar: Whether the location is within the polar circle
        - latitude: The input latitude
        - threshold: The polar circle threshold for the given obliquity
        - obliquity: The obliquity used
        - excess: How many degrees beyond the threshold (if polar)
        - house_system: The house system character
        - hemisphere: 'N' for north, 'S' for south
    """
    threshold = get_polar_latitude_threshold(eps)
    is_polar = abs(lat) > threshold
    excess = abs(lat) - threshold if is_polar else 0.0
    hemisphere = "N" if lat >= 0 else "S"

    return {
        "is_polar": is_polar,
        "latitude": lat,
        "threshold": threshold,
        "obliquity": eps,
        "excess": excess,
        "house_system": house_system,
        "hemisphere": hemisphere,
    }


# Default threshold for extreme latitude (80°). Project convention, not from any
# published standard: the |latitude| beyond which several quadrant house systems
# (Campanus, Regiomontanus, Topocentric) become numerically fragile. Chosen for
# this library and used only to flag reduced-accuracy regions, never to alter a
# computation.
EXTREME_LATITUDE_THRESHOLD = 80.0


def _is_extreme_latitude(
    lat: float, threshold: float = EXTREME_LATITUDE_THRESHOLD
) -> bool:
    """
    Check if latitude is at an extreme value where house calculations may be numerically unstable.

    At latitudes beyond 80°, many house systems (especially Campanus, Regiomontanus,
    and Topocentric) may produce results with reduced accuracy or numerical instability,
    even if they don't technically fail like Placidus/Koch at polar latitudes.

    This is distinct from _is_polar_circle() which checks the mathematical limit
    where Placidus/Koch become undefined (~66.5°).

    Args:
        lat: Geographic latitude in degrees
        threshold: Latitude threshold for "extreme" (default 80°)

    Returns:
        True if latitude is at an extreme value (abs(lat) >= threshold)

    Example:
        >>> _is_extreme_latitude(85.0)
        True
        >>> _is_extreme_latitude(70.0)
        False
    """
    return abs(lat) >= threshold


def get_extreme_latitude_info(lat: float, obliquity: float = 23.44) -> dict:
    """
    Get detailed information about a location's latitude classification for house calculations.

    This function provides comprehensive information about how latitude affects
    house calculations, including whether the location is at extreme latitude,
    within the polar circle, and which house systems are affected.

    Args:
        lat: Geographic latitude in degrees
        obliquity: True obliquity of the ecliptic in degrees.
                   Default is 23.44° (approximate J2000 value).

    Returns:
        Dictionary with latitude information:
        - latitude: The input latitude
        - is_extreme: True if abs(lat) >= 80° (may have numerical instability)
        - is_polar_circle: True if abs(lat) > 90° - obliquity (Placidus/Koch/Gauquelin fail)
        - polar_threshold: The polar circle threshold latitude
        - extreme_threshold: The extreme latitude threshold (80°)
        - hemisphere: 'N' for north, 'S' for south
        - affected_systems: List of house system codes that fail at this latitude
        - unstable_systems: List of house system codes that may be numerically unstable
        - stable_systems: List of house system codes that work reliably at all latitudes

    Example:
        >>> info = get_extreme_latitude_info(85.0)
        >>> info['is_extreme']
        True
        >>> info['is_polar_circle']
        True
        >>> info['affected_systems']
        ['P', 'K', 'G']
        >>> info['stable_systems']
        ['E', 'W', 'O', 'M', 'X', 'V', 'N']
    """
    polar_threshold = get_polar_latitude_threshold(obliquity)
    is_polar = abs(lat) > polar_threshold
    is_extreme = _is_extreme_latitude(lat)
    hemisphere = "N" if lat >= 0 else "S"

    # Systems that completely fail at polar latitudes
    affected_systems = ["P", "K", "G"] if is_polar else []

    # Systems that may be numerically unstable at extreme latitudes (>80°)
    # but don't technically fail
    unstable_systems = []
    if is_extreme:
        unstable_systems = ["C", "R", "T", "B", "H", "U", "Y", "F"]

    # Systems that work reliably at all latitudes
    stable_systems = ["E", "W", "O", "M", "X", "V", "N"]

    return {
        "latitude": lat,
        "is_extreme": is_extreme,
        "is_polar_circle": is_polar,
        "polar_threshold": polar_threshold,
        "extreme_threshold": EXTREME_LATITUDE_THRESHOLD,
        "hemisphere": hemisphere,
        "affected_systems": affected_systems,
        "unstable_systems": unstable_systems,
        "stable_systems": stable_systems,
    }


def _validate_cusps(cusps: list | tuple) -> tuple[bool, str | None]:
    """
    Validate house cusps for numerical sanity.

    Checks that all cusp values are within valid range [0, 360) and are finite.
    This helps detect numerical instability at extreme latitudes.

    Args:
        cusps: List or tuple of house cusp longitudes

    Returns:
        Tuple of (is_valid, error_message or None)

    Example:
        >>> _validate_cusps([0.0, 30.0, 60.0, 90.0])
        (True, None)
        >>> _validate_cusps([float('nan'), 30.0])
        (False, 'House cusp contains NaN value')
    """
    for i, cusp in enumerate(cusps):
        # Check for NaN
        if cusp != cusp:  # NaN check (NaN != NaN is True)
            return False, f"House cusp {i + 1} contains NaN value"

        # Check for infinity
        if cusp == float("inf") or cusp == float("-inf"):
            return False, f"House cusp {i + 1} contains infinite value"

        # Check range (should be in [0, 360))
        if cusp < 0 or cusp >= 360:
            return False, f"House cusp {i + 1} out of range [0, 360): {cusp}"

    return True, None


def _house_system(hsys_char: str) -> HouseSystem:
    """The registry row for a folded selector character.

    A selector the registry does not list is served by the Placidus default,
    which raises inside the polar circle like Placidus itself (the reference
    errors for e.g. hsys 'Z' at latitude 80 instead of quietly computing a
    fallback). The registry, ``_HOUSE_SYSTEMS``, closes the constructions
    section of this module.
    """
    return _HOUSE_SYSTEMS.get(hsys_char, _DEFAULT_HOUSE_SYSTEM)


def _raise_polar_circle_error(
    lat: float, eps: float, house_system: str, func_name: str
) -> None:
    """
    Raise a PolarCircleError with detailed information.

    Args:
        lat: Geographic latitude in degrees
        eps: True obliquity of the ecliptic in degrees
        house_system: House system character (e.g., 'P', 'K', 'G')
        func_name: Name of the function raising the error

    Raises:
        PolarCircleError: Always raised with detailed polar circle information
    """
    info = _get_polar_circle_info(lat, eps, house_system)
    threshold = info["threshold"]
    hemisphere = "Northern" if info["hemisphere"] == "N" else "Southern"

    # Map house system codes to names
    system_names = {
        "P": "Placidus",
        "K": "Koch",
        "G": "Gauquelin",
    }
    system_name = system_names.get(house_system, house_system)

    message = (
        f"{func_name}: {system_name} house system cannot be calculated at "
        f"latitude {abs(lat):.2f}°{info['hemisphere']} (within {hemisphere} polar circle). "
        f"Polar threshold for obliquity {eps:.2f}° is ±{threshold:.2f}°. "
        f"Consider using Porphyry ('O'), Equal ('E'), or Whole Sign ('W') house systems "
        f"which work at all latitudes, or use houses_with_fallback() for automatic fallback."
    )

    raise PolarCircleError(
        message=message,
        latitude=lat,
        threshold=threshold,
        obliquity=eps,
        house_system=house_system,
    )


def angular_diff(a: float, b: float) -> float:
    """
    Calculate signed angular difference (a - b) handling 360° wrapping.

    Used by horizontal house system to find ecliptic longitude for given azimuth.

    Args:
        a: First angle in degrees
        b: Second angle in degrees

    Returns:
        Signed difference in range [-180, 180]
    """
    diff = (a - b) % 360.0
    if diff > 180.0:
        diff -= 360.0
    return diff


def _calc_vertex(armc_deg: float, eps: float, lat: float) -> float:
    """
    Calculate the Vertex (intersection of Prime Vertical and Ecliptic in Western hemisphere).

    The Vertex is an auxiliary angle representing where the Prime Vertical (great circle
    through zenith perpendicular to meridian) intersects the ecliptic in the western sky.
    Often used in astrology for fateful encounters or significant relationships.

    At the equator (lat=0), the formula has a 1/tan(lat) singularity. The
    positive-latitude limit is used as a deterministic convention.

    Args:
        armc_deg: Right Ascension of Midheaven (sidereal time) in degrees
        eps: True obliquity of ecliptic in degrees
        lat: Geographic latitude in degrees

    Returns:
        Vertex longitude in degrees (western hemisphere)

    Precision: ~0.001° for non-equatorial latitudes
    """
    eps_rad = math.radians(eps)

    # Preserve real nonzero latitudes. Only the exact singularity needs a
    # deterministic one-sided limit.
    if lat == 0.0:
        lat = 1e-12

    # Evaluate the prime-vertical/ecliptic intersection derived in the
    # docstring, then choose its western antipode geometrically below.
    armc_rad = math.radians(armc_deg)
    lat_rad = math.radians(lat)

    num = -math.cos(armc_rad)
    den = math.sin(armc_rad) * math.cos(eps_rad) - math.sin(eps_rad) / math.tan(lat_rad)

    # 0/0 singularity at |lat| == eps with armc 90/270 (den collapses to
    # cos(eps) - cos(eps)): atan2 of the rounding noise is meaningless.
    # At this geometric degeneracy choose the western cardinal intersection.
    if abs(num) < 1e-12 and abs(den) < 1e-12:
        return (armc_deg + 270.0) % 360.0

    vtx_rad = math.atan2(num, den)
    vtx = math.degrees(vtx_rad) % 360.0

    # The atan2 above yields one of the two antipodal intersections of the
    # prime vertical with the ecliptic; the Vertex is the WESTERN one.
    # Disambiguate in the equatorial frame: convert the candidate to right
    # ascension and require an hour angle in (0, 180), i.e. west of the
    # upper meridian.  (A previous version compared the candidate's
    # *ecliptic* longitude against the ARMC, an *equatorial* longitude;
    # that mixed-frame test is off by up to ~2.5 deg and flipped the
    # Vertex by 180 deg whenever it fell within that margin of
    # ARMC + 180, which is reachable at |lat| < eps.
    vtx_r = math.radians(vtx)
    ra_vtx = (
        math.degrees(math.atan2(math.cos(eps_rad) * math.sin(vtx_r), math.cos(vtx_r)))
        % 360.0
    )
    hour_angle = (armc_deg - ra_vtx) % 360.0
    if not (0.0 < hour_angle < 180.0):
        vtx = (vtx + 180.0) % 360.0

    return vtx


def _ra_to_ecliptic_longitude(
    ra_deg: float, pole_height_deg: float, sin_obliquity: float, cos_obliquity: float
) -> float:
    """
    Convert right ascension to ecliptic longitude via spherical trigonometry.

    Given a right ascension and pole height, computes the corresponding ecliptic
    longitude using the standard spherical trigonometric relation:

        tan(λ) = sin(α) / (cos(ε)·cos(α) - sin(ε)·tan(φ))

    where α is the right ascension, ε is the obliquity, and φ is the pole height.
    Uses atan2 for unambiguous quadrant determination.

    Reference: Smart, "Textbook on Spherical Astronomy", Ch. 3;
               Meeus, "Astronomical Algorithms", Ch. 13

    Args:
        ra_deg: Right ascension in degrees
        pole_height_deg: Pole height (geographic latitude for ascendant) in degrees
        sin_obliquity: Pre-computed sin(obliquity)
        cos_obliquity: Pre-computed cos(obliquity)

    Returns:
        Ecliptic longitude in degrees [0, 360)
    """
    _NEAR_ZERO = 1e-10

    # Polar degenerate cases
    if abs(90.0 - pole_height_deg) < _NEAR_ZERO:
        return 180.0
    if abs(90.0 + pole_height_deg) < _NEAR_ZERO:
        return 0.0

    ra_rad = math.radians(ra_deg % 360.0)
    tan_pole = math.tan(math.radians(pole_height_deg))

    sin_ra = math.sin(ra_rad)
    cos_ra = math.cos(ra_rad)

    # Standard spherical trigonometry formula
    numerator = sin_ra
    denominator = cos_obliquity * cos_ra - sin_obliquity * tan_pole

    # 0/0 singularity (e.g. the CoAsc Munkasey path at |lat| == eps with
    # ra == 180): atan2 of the rounding noise is meaningless. Use the
    # continuous right-ascension limit.
    if abs(numerator) < 1e-12 and abs(denominator) < 1e-12:
        return ra_deg % 360.0

    # atan2 handles all quadrants correctly
    longitude = math.degrees(math.atan2(numerator, denominator)) % 360.0

    # Snap cardinal points to exact values (avoid floating-point drift)
    for cardinal in (90.0, 180.0, 270.0):
        if abs(longitude - cardinal) < _NEAR_ZERO:
            return cardinal
    if abs(longitude - 360.0) < _NEAR_ZERO or abs(longitude) < _NEAR_ZERO:
        return 0.0

    return longitude


def _calc_ascendant(
    armc_deg: float, eps: float, lat: float, pole_height: float
) -> float:
    """
    Calculate ecliptic longitude from ARMC and pole height.

    Wrapper around _ra_to_ecliptic_longitude that accepts obliquity in degrees
    and computes the trigonometric values internally.

    Based on standard spherical trigonometry:
        tan(λ) = sin(α) / (cos(ε)·cos(α) - sin(ε)·tan(φ))

    Reference: Smart, "Textbook on Spherical Astronomy";
               Meeus, "Astronomical Algorithms"

    Args:
        armc_deg: Right Ascension of MC in degrees
        eps: True obliquity of ecliptic in degrees
        lat: Geographic latitude (unused, kept for call-site compatibility)
        pole_height: Pole height (latitude parameter) in degrees

    Returns:
        Ecliptic longitude in degrees (0-360)
    """
    sin_obliquity = math.sin(math.radians(eps))
    cos_obliquity = math.cos(math.radians(eps))
    return _ra_to_ecliptic_longitude(
        armc_deg, pole_height, sin_obliquity, cos_obliquity
    )


def _sun_declination_analytic(tjdut: float) -> float:
    """Low-precision apparent solar declination (Meeus 1998, ch. 25).

    Fallback for the Sunshine ('I'/'i') Sun-declination fetch when the
    loaded ephemeris does not cover the date. An analytic solar position avoids
    the severe cusp distortion that a 0.0 substitute would cause at high latitude.
    Accuracy is ~0.01 deg, ample for cusp geometry.
    """
    t = (tjdut + _deltat(tjdut) - 2451545.0) / 36525.0
    l0 = (280.46646 + 36000.76983 * t + 0.0003032 * t * t) % 360.0
    m = math.radians((357.52911 + 35999.05029 * t - 0.0001537 * t * t) % 360.0)
    c = (
        (1.914602 - 0.004817 * t - 0.000014 * t * t) * math.sin(m)
        + (0.019993 - 0.000101 * t) * math.sin(2.0 * m)
        + 0.000289 * math.sin(3.0 * m)
    )
    omega = math.radians(125.04 - 1934.136 * t)
    lon_app = math.radians(l0 + c - 0.00569 - 0.00478 * math.sin(omega))
    eps = math.radians(
        23.4392911111
        - 0.0130041667 * t
        - 1.639e-7 * t * t
        + 5.036e-7 * t**3
        + 0.00256 * math.cos(omega)
    )
    return math.degrees(math.asin(math.sin(eps) * math.sin(lon_app)))


def _mc_below_horizon(armc_deg: float, eps: float, lat: float) -> bool:
    """Whether the midheaven of ``armc_deg`` culminates below the horizon.

    The midheaven's declination follows from its right ascension (the ARMC)
    and the obliquity; on the meridian its altitude satisfies
    ``sin(alt) = sin(lat) sin(dec) + cos(lat) cos(dec)``. Inside the polar
    circle the ecliptic point on the upper meridian can be below the horizon,
    and the systems flagged in the registry then report the other branch.
    """
    mc_dec_rad = math.atan(
        math.sin(math.radians(armc_deg)) * math.tan(math.radians(eps))
    )
    lat_rad = math.radians(lat)
    sin_alt = math.sin(lat_rad) * math.sin(mc_dec_rad) + math.cos(lat_rad) * math.cos(
        mc_dec_rad
    )
    return sin_alt < 0


def _rising_longitude(armc_deg: float, eps: float, lat: float) -> float:
    """The Ascendant: the ecliptic longitude rising on the eastern horizon.

    Solves the horizon/ecliptic intersection for the given ARMC, obliquity
    and latitude, then keeps the intersection whose hour angle puts it east
    of the meridian. When the ecliptic is tangent to the horizon (``|lat| ==
    90 - eps`` with the equinoctial colure on the meridian) the rising point
    degenerates and ``atan2(0, 0)`` is float noise: the Ascendant is then
    pinned to the east point, ARMC + 90, as the reference API does.
    """
    num = math.cos(math.radians(armc_deg))
    den = -(
        math.sin(math.radians(armc_deg)) * math.cos(math.radians(eps))
        + math.tan(math.radians(lat)) * math.sin(math.radians(eps))
    )
    if abs(num) < 1e-12 and abs(den) < 1e-12:
        asc = (armc_deg + 90.0) % 360.0
    else:
        asc = math.degrees(math.atan2(num, den)) % 360.0
        # Right ascension of that intersection and its hour angle from the
        # TRUE ARMC: an hour angle in (0, 180) puts the point west of the
        # meridian, i.e. setting, so the opposite intersection is rising.
        asc_r = math.radians(asc)
        eps_r = math.radians(eps)
        ra = (
            math.degrees(math.atan2(math.cos(eps_r) * math.sin(asc_r), math.cos(asc_r)))
            % 360.0
        )
        h_deg = (armc_deg - ra + 360.0) % 360.0
        if 0.0 < h_deg < 180.0:
            asc = (asc + 180.0) % 360.0

    # Ensure ASC is in [0, 360) range (handle floating-point near-360 values)
    if asc >= 360.0 or (360.0 - asc) < 1e-10:
        asc = 0.0
    return asc


def _sunshine_sun_declination(tjdut: float, iflag: int) -> float:
    """The Sun's apparent declination, for the Sunshine house systems.

    The ephemeris flags in ``iflag`` are forwarded to the Sun computation. A
    date outside the loaded ephemeris falls back to the low-precision
    analytic declination; any other failure of the Sun computation falls
    back to zero declination, the equinox behaviour.
    """
    try:
        # Extract ephemeris flags: FLG_JPLEPH=1, FLG_SWIEPH=2
        eph_flags = iflag & (FLG_JPLEPH | FLG_SWIEPH)
        sun_pos, _ = calc_ut(tjdut, SUN, FLG_EQUATORIAL | eph_flags)
        return sun_pos[1]  # Declination is second element in equatorial coords
    except EphemerisRangeError:
        return _sun_declination_analytic(tjdut)
    except (IndexError, TypeError, ValueError, CalculationError):
        return 0.0


def _houses_from_armc(
    armc_deg: float,
    lat: float,
    eps: float,
    hsys_char: str,
    sun_dec: float,
    func_name: str,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Cusps and angles of one house system from ARMC, latitude and obliquity.

    The single dispatcher behind :func:`houses` and :func:`houses_armc`. The
    system is looked up in the registry by its folded selector character;
    its row says which branch of the meridian is reported as MC, whether the
    system raises inside the polar circle and how many cusps the public
    tuple carries.

    Args:
        armc_deg: Right ascension of the midheaven in degrees, in [0, 360).
        lat: Geographic latitude in degrees.
        eps: True obliquity of the ecliptic in degrees.
        hsys_char: House-system selector, already folded by ``_hsys_to_char``.
        sun_dec: The Sun's declination in degrees (Sunshine systems only).
        func_name: Public entry point named in a ``PolarCircleError``.

    Returns:
        ``(cusps, ascmc)``: the 12 (36 for Gauquelin) cusp longitudes and the
        8 angles ``[Asc, MC, ARMC, Vertex, EquAsc, CoAsc (Koch), CoAsc
        (Munkasey), PolarAsc]``, every value in [0, 360).

    Raises:
        PolarCircleError: for a system that is undefined inside the polar
            circle (``raises_at_pole`` in the registry).
    """
    system = _house_system(hsys_char)

    # Frame the midheaven and the intermediate cusps are built on. The
    # systems that follow the MC above the horizon rotate it half a turn when
    # the upper culmination is below the horizon; the Ascendant and the other
    # angles always use the given ARMC and latitude.
    armc_active = armc_deg
    calc_lat = lat
    if system.mc_branch is _McBranch.ROTATE_FRAME and _mc_below_horizon(
        armc_deg, eps, lat
    ):
        armc_active = (armc_deg + 180.0) % 360.0
        calc_lat = -lat

    mc = _armc_to_mc(armc_active, eps)
    asc = _rising_longitude(armc_deg, eps, lat)

    # Vertex uses the original ARMC and the geometric equator convention in
    # _calc_vertex.
    vertex = _calc_vertex(armc_deg, eps, lat)

    # Equatorial Ascendant (East Point)
    # This is the intersection of the ecliptic with the celestial equator in the east
    # It's the ecliptic longitude where RA = ARMC + 90°
    equ_asc_ra = (armc_deg + 90.0) % 360.0
    # Convert RA to ecliptic longitude
    # tan(Lon) = tan(RA) / cos(eps)
    # y = sin(RA)
    # x = cos(RA) * cos(eps)
    equ_asc_ra_r = math.radians(equ_asc_ra)
    eps_r = math.radians(eps)
    y = math.sin(equ_asc_ra_r)
    x = math.cos(equ_asc_ra_r) * math.cos(eps_r)
    equ_asc = math.degrees(math.atan2(y, x)) % 360.0

    # Co-Ascendant W. Koch (coasc1)
    # Co-Ascendant (Koch) formula:
    # coasc1 = Asc(ARMC - 90°, latitude) + 180°
    # This is the Ascendant calculated 90° westward on the equator, then opposite point
    coasc_armc = (armc_deg - 90.0) % 360.0
    co_asc_koch = _calc_ascendant(coasc_armc, eps, lat, lat)

    # Add 180° to get opposite point
    co_asc_koch = (co_asc_koch + 180.0) % 360.0

    # Co-Ascendant M. Munkasey (coasc2)
    # Co-Ascendant (Munkasey) formula:
    # If lat >= 0: coasc2 = Asc(ARMC + 90°, 90° - lat)
    # If lat < 0:  coasc2 = Asc(ARMC + 90°, -90° - lat)
    # At the exact equator the northern (+90-lat) and southern (-90-lat) pole
    # heights give antipodal one-sided limits for coasc2 (180 vs 0).
    # Compatibility contract: the degeneracy is resolved per house system: the horizon
    # system 'H' takes the southern branch (coasc2 -> 0), while every other
    # system takes the northern branch (coasc2 -> 180). Away from lat == 0 the
    # sign of the latitude selects the branch unambiguously and the two
    # implementations already agree.
    coasc2_armc = (armc_deg + 90.0) % 360.0
    if lat > 0.0 or (lat == 0.0 and hsys_char != "H"):
        coasc2_lat = 90.0 - lat
    else:
        coasc2_lat = -90.0 - lat
    co_asc = _calc_ascendant(coasc2_armc, eps, coasc2_lat, coasc2_lat)

    # Polar Ascendant M. Munkasey (polasc)
    # Polar Ascendant formula:
    # polasc = Asc(ARMC - 90°, latitude)
    # Note: This is the same as coasc1 but WITHOUT the +180°
    polar_asc = _calc_ascendant(coasc_armc, eps, lat, lat)

    # ASCMC array with 8 elements (reference API compatible): Asc, MC, ARMC,
    # Vertex, EquAsc, coasc1 (W. Koch), coasc2 (M. Munkasey), polasc.
    ascmc = [asc, mc, armc_deg, vertex, equ_asc, co_asc_koch, co_asc, polar_asc]

    # Placidus, Koch and Gauquelin - and the default an unknown selector
    # reaches - cannot be calculated when abs(lat) + eps > 90.
    if system.raises_at_pole and _is_polar_circle(lat, eps):
        _raise_polar_circle_error(lat, eps, hsys_char, func_name)

    cusps = system.build_cusps(
        _HouseFrame(armc_active, calc_lat, eps, asc, mc, sun_dec)
    )

    if system.mc_branch is _McBranch.FROM_CUSPS:
        # The Sunshine constructions choose the meridian branch themselves;
        # the returned MC follows their cusp 10.
        ascmc[1] = cusps[10]
    elif system.mc_branch is _McBranch.ROTATE_REPORTED and _mc_below_horizon(
        armc_deg, eps, lat
    ):
        # APC: the parallel construction has already selected the eastern
        # Ascendant branch, so only the meridian pair follows the MC above the
        # horizon; rotating every cusp would undo that geometrical selection
        # at polar latitudes.
        ascmc[1] = (ascmc[1] + 180.0) % 360.0
        cusps[10] = ascmc[1]
        cusps[4] = (ascmc[1] + 180.0) % 360.0

    # Slot 0 of the cusp list is unused (reference API compatible: no padding
    # in the returned tuple). degnorm snaps the bare-%360 artifact (a
    # tiny-negative angle wrapping to exactly 360.0) back to 0.0, keeping
    # every output in [0, 360) like the reference API.
    return (
        tuple(degnorm(c) for c in cusps[1 : system.sectors + 1]),
        tuple(degnorm(a) for a in ascmc),
    )


def houses(
    tjdut: float, lat: float, lon: float, hsys: int = ord("P"), iflag: int = 0
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """
    Calculate astrological house cusps and angles for a given time and location.

    Reference API compatible function. Computes house divisions according to
    the specified house system and returns both house cusps and major angles (ASCMC).

    Args:
        tjdut: Julian Day in Universal Time (UT1)
        lat: Geographic latitude in degrees (positive North, negative South)
        lon: Geographic longitude in degrees (positive East, negative West)
        hsys: House system identifier (e.g., ord('P') for Placidus, ord('K') for Koch)
        iflag: Calculation flags (optional, default 0). Use FLG_MOSEPH for extended
               date range (-3000 to +3000 CE). The ephemeris flag is propagated to
               internal calculations (e.g., Sun position for Sunshine houses).

    Returns:
        Tuple containing:
            - cusps: Tuple of 12 house cusp longitudes (houses 1-12) in degrees
            - ascmc: Tuple of 8 angles: [Asc, MC, ARMC, Vertex, EquAsc, CoAsc, CoAscKoch, PolarAsc]

    House Systems Supported:
        'P' = Placidus, 'K' = Koch, 'R' = Regiomontanus, 'C' = Campanus,
        'E'/'A' = Equal, 'W' = Whole Sign, 'O' = Porphyry, 'B' = Alcabitius,
        'T' = Topocentric, 'M' = Morinus, 'X' = Meridian, 'H' = Horizontal,
        'V' = Vehlow, 'G' = Gauquelin, 'U' = Krusinski, 'F' = Carter,
        'Y' = APC, 'N' = Natural Gradient

    Example:
        >>> cusps, ascmc = houses(2451545.0, 51.5, -0.12, ord('P'))  # London, Placidus
        >>> asc, mc = ascmc[0], ascmc[1]
        >>> house_1_start = cusps[0]  # First house cusp
    """
    # Validate latitude and longitude ranges
    validate_coordinates(lat, lon, "houses")

    # Sidereal time (ARMC) and true obliquity from the long-term model.
    # ARMC = apparent sidereal time at longitude 0 + lon. Both the sidereal time
    # and the obliquity use the Vondrák 2011 long-term realization (see
    # sidereal_longterm), so cusps remain correct across the whole supported date
    # range rather than diverging at remote epochs.
    armc0_deg, eps = _house_armc_obliquity(tjdut)
    armc_deg = (armc0_deg + lon) % 360.0

    hsys_char = _hsys_to_char(hsys)

    # The Sun's declination is read by the Sunshine systems only, and this is
    # the one entry point that has the date to compute it.
    sun_dec = 0.0
    if _house_system(hsys_char).reads_sun_declination:
        sun_dec = _sunshine_sun_declination(tjdut, iflag)

    return _houses_from_armc(armc_deg, lat, eps, hsys_char, sun_dec, "houses")


# Names used in the messages of the two *_with_fallback entry points.
_FALLBACK_SYSTEM_NAMES: dict[str, str] = {
    "P": "Placidus",
    "K": "Koch",
    "G": "Gauquelin",
    "O": "Porphyry",
    "E": "Equal",
    "W": "Whole Sign",
    "R": "Regiomontanus",
    "C": "Campanus",
    "T": "Topocentric",
    "B": "Alcabitius",
    "M": "Morinus",
    "X": "Meridian",
    "H": "Horizontal",
    "V": "Vehlow",
    "U": "Krusinski",
    "F": "Carter",
    "Y": "APC",
    "N": "Natural Gradient",
}


def _houses_with_fallback_impl(
    compute: Callable[[Any], tuple[tuple[float, ...], tuple[float, ...]]],
    classify_latitude: Callable[[], dict],
    lat: float,
    hsys: int,
    fallback_hsys: int,
    validate_cusps: bool,
) -> tuple[tuple[float, ...], tuple[float, ...], bool, str | None]:
    """Shared body of :func:`houses_with_fallback` and its ARMC twin.

    ``compute(selector)`` runs the underlying entry point for one house
    system; ``classify_latitude()`` is the extreme-latitude classification
    that entry point consults for its warning (:func:`get_extreme_latitude_info`
    with the obliquity it knows). The selectors are normalised here, the
    primary system is tried, and the fallback is used when the primary raises
    :class:`PolarCircleError` or, if ``validate_cusps``, returns cusps that
    fail :func:`_validate_cusps`.
    """
    hsys_char = _hsys_to_char(hsys)
    fallback_char = (
        _int_hsys_to_char(fallback_hsys)
        if isinstance(fallback_hsys, int)
        else fallback_hsys
    )
    primary_name = _FALLBACK_SYSTEM_NAMES.get(hsys_char, hsys_char)
    fallback_name = _FALLBACK_SYSTEM_NAMES.get(fallback_char, fallback_char)

    # Check for extreme latitude before calculation
    is_extreme = _is_extreme_latitude(lat)

    try:
        cusps, ascmc = compute(hsys)

        # Validate cusps if requested
        if validate_cusps:
            is_valid, validation_error = _validate_cusps(cusps)
            if not is_valid:
                # Invalid cusps detected - fall back to stable system
                cusps, ascmc = compute(fallback_hsys)
                warning = (
                    f"{primary_name} house system produced invalid cusps at latitude "
                    f"{abs(lat):.2f}° ({validation_error}). Using {fallback_name} as fallback."
                )
                return cusps, ascmc, True, warning

        # Generate warning for extreme latitudes even if calculation succeeded
        if is_extreme:
            info = classify_latitude()
            if hsys_char in info["unstable_systems"]:
                warning = (
                    f"{primary_name} house system may have reduced accuracy at "
                    f"extreme latitude {abs(lat):.2f}°{info['hemisphere']}. "
                    f"Consider using a stable system like Porphyry, Equal, or Whole Sign."
                )
                return cusps, ascmc, False, warning

        return cusps, ascmc, False, None
    except PolarCircleError as e:
        # Use fallback house system
        cusps, ascmc = compute(fallback_hsys)

        warning = (
            f"{primary_name} house system unavailable at latitude {abs(lat):.2f}° "
            f"(polar circle threshold: {e.threshold:.2f}°). "
            f"Using {fallback_name} as fallback."
        )
        return cusps, ascmc, True, warning


def houses_with_fallback(
    jd_ut: float,
    lat: float,
    lon: float,
    hsys: int,
    fallback_hsys: int = ord("O"),
    validate_cusps: bool = True,
) -> tuple[tuple[float, ...], tuple[float, ...], bool, str | None]:
    """
    Calculate house cusps with automatic fallback for polar latitudes.

    This convenience function attempts to calculate houses using the requested
    system, and automatically falls back to a polar-safe system (default: Porphyry)
    if the location is within the polar circle.

    For extreme latitudes (>80°), this function also:
    - Validates cusp values for numerical sanity (NaN, infinity, range)
    - Includes warnings about potential numerical instability
    - Falls back to a stable system if cusps are invalid

    This is useful for applications that need to handle polar latitudes gracefully
    without explicit error handling for each calculation.

    Args:
        jd_ut: Julian Day in Universal Time
        lat: Geographic latitude in degrees (positive North, negative South)
        lon: Geographic longitude in degrees (positive East, negative West)
        hsys: Primary house system identifier (e.g., ord('P') for Placidus)
        fallback_hsys: Fallback house system for polar latitudes.
                       Default: ord('O') (Porphyry).
                       Other good choices: ord('E') (Equal), ord('W') (Whole Sign)
        validate_cusps: If True (default), validate cusp values for numerical sanity
                        and fall back if invalid cusps are detected.

    Returns:
        Tuple containing:
            - cusps: Tuple of 12 house cusp longitudes (houses 1-12) in degrees
            - ascmc: Tuple of 8 angles
            - used_fallback: True if fallback system was used
            - warning_message: Informational message if fallback was used or
                               extreme latitude warning, else None

    Example:
        >>> import libephemeris as ephem
        >>> jd = 2451545.0
        >>> # This would fail with Placidus at 70°N, but uses fallback
        >>> cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
        ...     jd, 70.0, 0.0, ord('P')
        ... )
        >>> if used_fallback:
        ...     print(f"Used fallback: {warning}")
        >>>
        >>> # At extreme latitudes (>80°), you may get a warning even for non-failing systems
        >>> cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
        ...     jd, 85.0, 0.0, ord('C')  # Campanus at 85°N
        ... )
        >>> if warning:
        ...     print(f"Warning: {warning}")

    See Also:
        houses: Standard function that raises PolarCircleError at polar latitudes
        get_polar_latitude_threshold: Returns the threshold for a given obliquity
        get_extreme_latitude_info: Returns detailed latitude classification information
    """
    return _houses_with_fallback_impl(
        lambda selector: houses(jd_ut, lat, lon, selector),
        lambda: get_extreme_latitude_info(lat),
        lat,
        hsys,
        fallback_hsys,
        validate_cusps,
    )


def houses_armc_with_fallback(
    armc: float,
    lat: float,
    eps: float,
    hsys: int,
    fallback_hsys: int = ord("O"),
    validate_cusps: bool = True,
    ascmc9: float = 0.0,
) -> tuple[tuple[float, ...], tuple[float, ...], bool, str | None]:
    """
    Calculate house cusps from ARMC with automatic fallback for polar latitudes.

    Similar to houses_with_fallback, but calculates from ARMC instead of
    Julian Day. Includes the same extreme latitude handling and cusp validation.

    Args:
        armc: Right Ascension of Medium Coeli in degrees (0-360)
        lat: Geographic latitude in degrees (positive North, negative South)
        eps: True obliquity of the ecliptic in degrees
        hsys: Primary house system identifier (e.g., ord('P') for Placidus)
        fallback_hsys: Fallback house system for polar latitudes.
                       Default: ord('O') (Porphyry).
        validate_cusps: If True (default), validate cusp values for numerical sanity
                        and fall back if invalid cusps are detected.
        ascmc9: The Sun's declination, used only by the Sunshine house system
                ('I'/'i') and forwarded to the underlying cusp solution; it is
                ignored by every other house system.

    Returns:
        Tuple containing:
            - cusps: Tuple of 12 house cusp longitudes (houses 1-12) in degrees
            - ascmc: Tuple of 8 angles
            - used_fallback: True if fallback system was used
            - warning_message: Informational message if fallback was used or
                               extreme latitude warning, else None

    See Also:
        houses_armc: Standard function that raises PolarCircleError at polar latitudes
        get_extreme_latitude_info: Returns detailed latitude classification information
    """
    # ascmc9 (the Sun's declination) is forwarded so the Sunshine system
    # ('I'/'i') gets it; houses_armc ignores it for every other system.
    return _houses_with_fallback_impl(
        lambda selector: houses_armc(armc, lat, eps, selector, ascmc9),
        lambda: get_extreme_latitude_info(lat, eps),
        lat,
        hsys,
        fallback_hsys,
        validate_cusps,
    )


def houses_armc(
    armc: float, lat: float, eps: float, hsys: int = ord("P"), ascmc9: float = 0.0
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """
    Calculate house cusps and angles from ARMC (Right Ascension of Medium Coeli).

    This function calculates house cusps directly from the ARMC value instead of
    from a Julian Day. This is useful when you have a pre-calculated ARMC or when
    working with house systems that depend only on ARMC, latitude, and obliquity.

    Reference API compatible function (houses_armc equivalent).

    Args:
        armc: Right Ascension of Medium Coeli in degrees (0-360)
        lat: Geographic latitude in degrees (positive North, negative South)
        eps: True obliquity of the ecliptic in degrees
        hsys: House system identifier (e.g., ord('P') for Placidus, ord('K') for Koch)

    Returns:
        Tuple containing:
            - cusps: Tuple of 12 house cusp longitudes (houses 1-12) in degrees
            - ascmc: Tuple of 8 angles: [Asc, MC, ARMC, Vertex, EquAsc, CoAsc, CoAscKoch, PolarAsc]

    House Systems Supported:
        'P' = Placidus, 'K' = Koch, 'R' = Regiomontanus, 'C' = Campanus,
        'E'/'A' = Equal, 'W' = Whole Sign, 'O' = Porphyry, 'B' = Alcabitius,
        'T' = Topocentric, 'M' = Morinus, 'X' = Meridian, 'H' = Horizontal,
        'V' = Vehlow, 'G' = Gauquelin, 'U' = Krusinski, 'F' = Carter,
        'Y' = APC, 'N' = Natural Gradient

    Example:
        >>> # Calculate obliquity for J2000.0
        >>> eps = 23.4393  # approximate true obliquity
        >>> armc = 292.957  # ARMC in degrees
        >>> cusps, ascmc = houses_armc(armc, 41.9, eps, ord('P'))
        >>> asc, mc = ascmc[0], ascmc[1]
    """
    # Validate latitude range (must be in [-90, 90])
    validate_latitude(lat, "houses_armc")

    # Normalize ARMC to 0-360
    armc_deg = armc % 360.0

    hsys_char = _hsys_to_char(hsys)

    # ascmc9 carries the Sun's declination (reference convention); only the
    # Sunshine systems read it.
    return _houses_from_armc(armc_deg, lat, eps, hsys_char, ascmc9, "houses_armc")


# House systems whose reference cusp SPEEDS follow a closed-form analytic rule
# instead of the finite-difference derivative of the reported cusp positions.
def houses_armc_ex2(
    armc: float,
    lat: float,
    eps: float,
    hsys: int = ord("P"),
    ascmc9: float = 0.0,
) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    """
    Extended house calculation from ARMC returning cusps, angles, and their velocities.

    This function combines houses_armc() with velocity calculations similar to
    houses_ex2(). It calculates house cusps directly from the ARMC value and
    also returns the velocities (derivatives) of house cusps and angles.

    Velocities are always calculated (compatibility contract).
    The ``ascmc9`` parameter carries the Sun's declination for the Sunshine
    house system ('I'/'i') and is forwarded to the underlying cusp solution
    (matching the reference, whose houses_armc_ex2 also uses it); it is
    ignored by every other house system.

    Velocities are calculated using centered finite differences, with ARMC
    shifted by ±1 second (Koch/Placidus) or ±1 minute (other systems).

    Note: The sidereal flag (FLG_SIDEREAL) has no effect since sidereal mode
    requires a Julian Day for ayanamsa calculation, which is not available when
    using ARMC directly.

    Args:
        armc: Right Ascension of Medium Coeli in degrees (0-360)
        lat: Geographic latitude in degrees (positive North, negative South)
        eps: True obliquity of the ecliptic in degrees
        hsys: House system identifier (e.g., ord('P') for Placidus, ord('K') for Koch)
        ascmc9: Optional parameter for Sunshine house system (default 0.0)

    Returns:
        Tuple containing:
            - cusps: Tuple of 12 house cusp longitudes in degrees
            - ascmc: Tuple of 8 angles (Asc, MC, ARMC, Vertex, EquAsc, CoAsc, CoAscKoch, PolarAsc)
            - cusps_speed: Tuple of 12 house cusp velocities in degrees/day
            - ascmc_speed: Tuple of 8 angle velocities in degrees/day

    Example:
        >>> eps = 23.4393  # true obliquity
        >>> armc = 292.957  # ARMC in degrees
        >>> cusps, ascmc, cusps_speed, ascmc_speed = houses_armc_ex2(
        ...     armc, 41.9, eps, ord('P')
        ... )
        >>> # cusps_speed[0] is the velocity of the 1st house cusp (same as ASC)
    """
    # Calculate positions at current ARMC. ascmc9 (the Sun's declination) is
    # forwarded so the Sunshine system ('I'/'i') gets it; houses_armc ignores
    # it for every other system.
    cusps, ascmc = houses_armc(armc, lat, eps, hsys, ascmc9)

    # Compute d(cusp)/d(ARMC) via centered finite differences, then
    # scale by the sidereal rotation rate to obtain deg/day.
    #
    # The Earth completes 360.98564736629° of ARMC per mean solar day
    # (one extra rotation relative to the Sun).  Previous code used
    # 360°/day (solar rate), which under-estimated all speeds by the
    # ratio 360/360.9856 ≈ 0.27 %.
    _SIDEREAL_RATE = 360.98564736629  # ARMC deg per mean solar day (Meeus ch. 12)

    # A one-second rotation step balances centered-difference truncation and
    # binary64 roundoff for the smooth spherical functions used here.
    d_armc = _SIDEREAL_RATE / 86400.0  # sidereal degrees per 1 second

    # Calculate positions at ARMC ± d_armc
    cusps_before, ascmc_before = houses_armc(armc - d_armc, lat, eps, hsys, ascmc9)
    cusps_after, ascmc_after = houses_armc(armc + d_armc, lat, eps, hsys, ascmc9)

    def angular_diff_local(pos2: float, pos1: float) -> float:
        """Calculate angular difference handling 360° wraparound."""
        diff = pos2 - pos1
        if diff > 180:
            diff -= 360
        elif diff < -180:
            diff += 360
        return diff

    def fd_speed(
        after: float, before: float, *, sign_pinned_cusp: bool = False
    ) -> float:
        """Centered-difference rate in deg/day; 0.0 across a discontinuity.

        A change larger than 90 degrees across the stencil is treated as a
        branch discontinuity, where no finite derivative exists.
        """
        diff = angular_diff_local(after, before)
        if abs(diff) > 90.0:
            return 0.0
        # Sign-pinned wheels step by exactly 30 degrees when the Ascendant
        # changes sign; that is a discontinuity, not a rate (see houses_ex2).
        # CUSPS ONLY: near the polar circle the Ascendant itself genuinely
        # sweeps ~30 degrees across the stencil, so applying this test to the
        # angles zeroed a real rate — and made the same angle read 0.0 for
        # 'W'/'N' and 1296000 deg/day for every other system.
        if (
            sign_pinned_cusp
            and _hsys_to_char(hsys) in _SIGN_PINNED_HSYS
            and abs(abs(diff) - 30.0) < 1e-6
        ):
            return 0.0
        return diff / (2 * d_armc) * _SIDEREAL_RATE

    # Velocities in deg/day = d(cusp)/d(ARMC) * (ARMC deg/day)
    # Since d_armc already equals _SIDEREAL_RATE * dt, dividing by
    # (2*d_armc) and then multiplying by _SIDEREAL_RATE is equivalent to
    # dividing by (2*dt).  Using d_armc directly avoids a separate dt
    # variable: speed = Δcusp / (2*d_armc) * _SIDEREAL_RATE.
    cusps_speed = tuple(
        fd_speed(cusps_after[i], cusps_before[i], sign_pinned_cusp=True)
        for i in range(len(cusps))
    )

    ascmc_speed = tuple(
        fd_speed(ascmc_after[i], ascmc_before[i]) for i in range(len(ascmc))
    )

    # The speeds are the derivative of the reported cusp functions with respect
    # to ARMC, scaled by the sidereal rate. This is the most accurate speed
    # obtainable from an ARMC input alone: with only the ARMC (and no Julian
    # Day) the small obliquity-rate term dε/dt — worth ~0.01 deg/day — cannot be
    # included. Callers that have the time should prefer houses_ex2(), which
    # finite-differences the full house solution in time and so captures every
    # time-dependent term exactly.
    #
    # degnorm on the angle outputs (not the speeds) snaps the bare-%360
    # artifact (exactly 360.0 from a tiny-negative angle) back to 0.0.
    return (
        tuple(degnorm(c) for c in cusps),
        tuple(degnorm(a) for a in ascmc),
        cusps_speed,
        ascmc_speed,
    )


def _houses_on_reference_plane(
    tjdut: float,
    lat: float,
    hsys_char: str,
    trop_ascmc: tuple,
    flags: int,
    plane_to_equator: Callable[[np.ndarray], np.ndarray],
    zero_point: float,
    node_reflection_epoch: float | None = None,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Build the houses on a reference plane other than the ecliptic of date.

    The fixed-epoch sidereal modes and the ``SIDBIT_ECL_T0`` /
    ``SIDBIT_SSY_PLANE`` projections measure the zodiac on another plane
    (the mean ecliptic of an epoch, or the invariable plane). The house
    construction is then carried out on that plane instead of the ecliptic
    of date, with the plane's node and inclination taking the place of the
    equinox and the obliquity:

    1. the inclination between the plane and the true equator of date is the
       construction obliquity ``eps_p``;
    2. the ascending node of the plane on the equator of date is the origin
       the ARMC is re-based to (``armc_p``), so :func:`houses_armc` returns
       cusp arcs measured from that node along the plane;
    3. adding the node's longitude on the plane and subtracting the zero
       point turns every arc into a sidereal longitude.

    The reported ARMC slot carries the node-rebased ``armc_p`` the
    construction runs on. Aries (``N``) cusps stay anchored at zero degrees
    of the zodiac in use, like the general sidereal path. Both callers keep
    ``i`` as the distinct Makransky construction: the collapse of ``i`` onto
    ``I`` belongs to the plain ayanamsha path of :func:`houses_ex` only.

    Args:
        tjdut: Julian Day in UT1.
        lat: Geographic latitude in degrees.
        hsys_char: Canonical house-system character.
        trop_ascmc: The tropical ASCMC tuple; only ``trop_ascmc[2]`` (the
            plain apparent ARMC) is consumed, the cusps are rebuilt on the
            plane.
        flags: Calculation flags (ephemeris bits forwarded to the Sunshine
            Sun fetch; FLG_RADIANS is applied by the caller).
        plane_to_equator: Maps the precession-nutation matrix (true equator
            of date <- ICRS) to the rotation taking plane coordinates onto
            the true equator of date.
        zero_point: Longitude on the plane, from the plane's own origin, of
            the sidereal zero point.
        node_reflection_epoch: When given (the fixed-epoch modes), the node
            longitude is folded to (-180, 180] and its sign is taken from
            ``t - epoch`` rather than from the geometry: inside the interval
            where nutation puts the node on the other side of the plane's
            origin, the anchor reflects while the ARMC does not. ``None``
            keeps the geometric node.

    Returns:
        (cusps, ascmc) in degrees; the caller applies degnorm / FLG_RADIANS.
    """
    import numpy as np

    from .precession_vondrak import vondrak_pn_matrix

    # True equator of date <- ICRS, using the SAME Vondrak precession + IAU
    # nutation chain as the tropical house frame (_house_armc_obliquity), so
    # the rebased ARMC stays consistent with trop_ascmc[2].
    jd_tt = tjdut + _deltat(tjdut)
    dpsi, deps = get_cached_nutation(jd_tt)
    pn_tuples, _eps_true = vondrak_pn_matrix(jd_tt, dpsi, deps)
    m = plane_to_equator(np.array(pn_tuples))

    # Inclination between the plane and the equator of date (the construction
    # obliquity), and the ascending node of the equator of date on the plane:
    # its longitude on the plane (lon_node) and its right ascension (ang_n).
    z = m.T @ np.array([0.0, 0.0, 1.0])
    eps_p = math.degrees(math.acos(max(-1.0, min(1.0, float(z[2])))))
    n_p = np.cross(z, [0.0, 0.0, 1.0])
    n_p /= np.linalg.norm(n_p)
    lon_node = math.degrees(math.atan2(float(n_p[1]), float(n_p[0]))) % 360.0
    n_eq = m @ n_p
    ang_n = math.degrees(math.atan2(float(n_eq[1]), float(n_eq[0]))) % 360.0
    armc_p = (trop_ascmc[2] - ang_n) % 360.0
    if node_reflection_epoch is not None:
        d_node = lon_node if lon_node <= 180.0 else lon_node - 360.0
        lon_node = -abs(d_node) if jd_tt >= node_reflection_epoch else abs(d_node)

    sun_dec = 0.0
    if hsys_char in ("I", "i"):
        sun_dec = _sunshine_sun_declination(tjdut, flags)

    eng_cusps, eng_ascmc = houses_armc(
        armc_p, lat, eps_p, ord(hsys_char), ascmc9=sun_dec
    )

    # Cusp arcs (measured from the node by houses_armc) -> plane longitude
    # (add the node's plane longitude) -> sidereal longitude (subtract the
    # zero point). Aries ('N') stays anchored at 0 deg of the sidereal zodiac.
    off = lon_node - zero_point
    if hsys_char == "N":
        cusps = tuple(float(i * 30) for i in range(12))
    else:
        cusps = tuple((off + c) % 360.0 for c in eng_cusps)

    ascmc_list = list(trop_ascmc)
    for i in (0, 1, 3, 4, 5, 6, 7):
        ascmc_list[i] = (off + eng_ascmc[i]) % 360.0
    # ascmc[2] (ARMC) reports the node-rebased armc_p the construction runs on.
    ascmc_list[2] = armc_p
    return cusps, tuple(ascmc_list)


def _houses_fixed_epoch_sidereal(
    tjdut: float,
    lat: float,
    hsys_char: str,
    trop_cusps: tuple,
    trop_ascmc: tuple,
    flags: int,
    sid_mode: int,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Sidereal houses for the fixed-epoch modes (J2000/J1900/B1950).

    These modes are evaluated on the mean ecliptic of their epoch ``t0``. The
    reference plane handed to :func:`_houses_on_reference_plane` is that
    ecliptic, reached from the true equator of date through the transposed
    Vondrák precession matrix of ``t0`` and the mean obliquity of ``t0``;
    the cusp arcs are measured from the mean equinox of ``t0`` along it.
    Modes like ``SIDM_GALALIGN_MARDYKS`` carry a constant longitude offset of
    the ``t0`` frame (their defining ayanamsha at ``t0``): that offset is the
    zero point. The node anchor is reflected about ``t0``; the ARMC is not.
    """
    import numpy as np

    import erfa

    from .precession_vondrak import vondrak_precession_matrix
    from .sidereal_epoch import FIXED_EPOCH_LON_OFFSET, FIXED_EPOCH_T0, _rot_x

    t0 = FIXED_EPOCH_T0[sid_mode]
    v_t0 = np.array(vondrak_precession_matrix(t0))
    eps_t0 = erfa.obl06(t0, 0.0)

    def mean_ecliptic_of_t0_to_equator(pn: np.ndarray) -> np.ndarray:
        # mean ecliptic of t0 -> mean equator of t0 (obliquity of t0) -> ICRS
        # (transposed precession of t0) -> true equator of date (pn).
        return pn @ v_t0.T @ _rot_x(-eps_t0)

    return _houses_on_reference_plane(
        tjdut,
        lat,
        hsys_char,
        trop_ascmc,
        flags,
        mean_ecliptic_of_t0_to_equator,
        FIXED_EPOCH_LON_OFFSET.get(sid_mode, 0.0),
        node_reflection_epoch=t0,
    )


def _houses_sidbit_projection(
    tjdut: float,
    lat: float,
    hsys_char: str,
    trop_ascmc: tuple,
    flags: int,
    sid_mode: int,
    sid_bits: int,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Sidereal houses for the SIDBIT_ECL_T0 / SIDBIT_SSY_PLANE projections.

    These two SIDBIT flags request the sidereal zodiac measured on a reference
    plane other than the ecliptic of date:

      * ``SIDBIT_ECL_T0``  -> the mean ecliptic and equinox of the mode's
        reference epoch ``t0`` (``AYANAMSHA_DEFINING``); the zero point is the
        mode's ayanamsha at ``t0``.
      * ``SIDBIT_SSY_PLANE`` -> the solar-system invariable plane (Souami &
        Souchay 2012); the zero point is the J2000 ayanamsha direction projected
        onto that plane.

    The whole house construction is built ON that plane by
    :func:`_houses_on_reference_plane` rather than shifting the
    ecliptic-of-date cusps by a scalar ayanamsha: the ARMC is re-based to the
    ascending node of the true equator of date on the plane (so the reported
    ARMC itself moves -- up to ~3.9 deg for the tilted invariable plane), the
    obliquity becomes the inclination between the two planes, and the cusp
    arcs are measured from the mode's sidereal zero point along the plane.
    Projecting each ecliptic-of-date cusp point by point is not equivalent:
    it leaves the equatorial ARMC unchanged and misplaces the
    invariable-plane cusps by up to several degrees.

    The projection-plane definitions and the sidereal zero point are taken from
    the SAME frame model the position path uses
    (:mod:`libephemeris.sidereal_epoch`, ``_sidbit_projection_calc``), so cusps
    and bodies stay in one self-consistent frame.

    Measured residual envelope against the reference:

      * ``SIDBIT_SSY_PLANE`` inherits the invariable-plane zero point of the
        calc path, whose Souami & Souchay (2012) orientation differs from the
        reference's by a near-constant ~30 arcsec. Houses reproduce the
        reference to that same shared envelope, not more tightly, by design --
        matching the calc positions is required for a body and its house to
        agree. The construction itself (ARMC, obliquity, cusp arcs on the
        plane) agrees to <0.2 arcsec once that shared offset is removed.
      * ``SIDBIT_ECL_T0`` carries no such calc offset (the position path is at
        parity), and at ordinary modern epochs it reproduces the reference
        closely -- e.g. the Ascendant to a few arcsec at |lat| ~ 40 deg. Two
        few-arcsec-class differences remain and are NOT fitted here: (1) the
        reference measures this near-ecliptic-of-date plane in a mean-equinox
        frame whose nutation treatment near the mode's reference epoch t0 the
        node-rebased ``houses_armc`` reduction does not fully reproduce (at the
        exact t0, where the plane coincides with the ecliptic of date, the
        residual is dominated by the of-date nutation, ~20-30 arcsec on the
        Ascendant); (2) the reported ARMC differs from the geometric
        node-rebased value by a mode-specific ~4-6 arcsec. Through the quadrant
        house systems these differences are amplified at high geographic
        latitude, reaching tens of arcsec near |lat| ~ 55-60 deg. That is the
        documented ECL_T0 envelope; the strongly tilted SSY_PLANE construction
        is unaffected because its large reconstruction dominates the of-date
        nutation.

    Args:
        tjdut: Julian Day in UT1.
        lat: Geographic latitude in degrees.
        hsys_char: Canonical house-system character.
        trop_ascmc: The tropical ASCMC tuple; only ``trop_ascmc[2]`` (the plain
            apparent ARMC) is consumed -- the cusps are rebuilt on the plane.
        flags: Calculation flags (ephemeris bits forwarded to the Sunshine Sun
            fetch; FLG_RADIANS applied by the caller).
        sid_mode: The base sidereal mode.
        sid_bits: The active SIDBIT projection flags.

    Returns:
        (cusps, ascmc) in degrees; the caller applies degnorm / FLG_RADIANS.
    """
    import numpy as np

    import erfa

    from .precession_vondrak import vondrak_precession_matrix
    from .sidereal_epoch import (
        _ecliptic_of_t0_matrix,
        _invariable_plane_matrix,
        _rot_x,
        ssy_plane_zero_point_deg,
    )
    from .planets import _calc_ayanamsa

    _J2000 = 2451545.0
    v_j2k = np.array(vondrak_precession_matrix(_J2000))
    eps_j2000 = erfa.obl06(_J2000, 0.0)

    # Projection plane, expressed as J2000-mean-ecliptic -> plane frame, and the
    # sidereal zero point in that frame -- identical to the calc SIDBIT path.
    # Flag precedence: with BOTH projection bits set,
    # SIDBIT_ECL_T0 takes precedence over SIDBIT_SSY_PLANE.
    if sid_bits & SIDBIT_ECL_T0:
        from .planets import _ecl_t0_epoch_jd

        t0_jd = _ecl_t0_epoch_jd(sid_mode)
        base = _ecliptic_of_t0_matrix(t0_jd)
        zero_point = _calc_ayanamsa(t0_jd, sid_mode)
    else:  # SIDBIT_SSY_PLANE
        base = _invariable_plane_matrix()
        zero_point = ssy_plane_zero_point_deg(_calc_ayanamsa(_J2000, sid_mode))

    def plane_to_equator(pn: np.ndarray) -> np.ndarray:
        # plane -> mean ecliptic J2000 (transposed base) -> mean equator J2000
        # (obliquity of J2000) -> ICRS (transposed precession of J2000) -> true
        # equator of date (pn).
        return pn @ v_j2k.T @ _rot_x(-eps_j2000) @ base.T

    # 'i' keeps the Makransky construction here: it differs from 'I' above
    # about 59 degrees of latitude (78.0 deg on cusp 2 at lat 60, 100.7 deg
    # at lat 65) and is identical below it.
    return _houses_on_reference_plane(
        tjdut, lat, hsys_char, trop_ascmc, flags, plane_to_equator, zero_point
    )


def houses_ex(
    tjdut: float, lat: float, lon: float, hsys: int = ord("P"), flags: int = 0
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """
    Extended house calculation with sidereal zodiac support.

    Similar to houses() but applies ayanamsa correction when FLG_SIDEREAL
    flag is set, converting tropical positions to sidereal.

    For Ascendant-based house systems (Whole Sign, Equal, Vehlow), the cusps
    are recalculated using the sidereal Ascendant to ensure proper sign alignment.

    Args:
        tjdut: Julian Day in Universal Time (UT1)
        lat: Geographic latitude in degrees
        lon: Geographic longitude in degrees
        hsys: House system identifier (int or bytes)
        flags: Calculation flags bitmask (e.g., FLG_SIDEREAL, FLG_MOSEPH).
               Use FLG_MOSEPH for extended date range (-3000 to +3000 CE).
               The ephemeris flag is propagated to all internal calculations.

    Returns:
        Tuple containing:
            - cusps: Tuple of 12 house cusp longitudes in degrees
            - ascmc: Tuple of 8 angles (Asc, MC corrected if sidereal)

    Example:
        >>> from libephemeris import set_sid_mode, SIDM_LAHIRI
        >>> set_sid_mode(SIDM_LAHIRI)
        >>> cusps, ascmc = houses_ex(2451545.0, 51.5, -0.12, ord('P'), FLG_SIDEREAL)
    """
    # House-system character, needed up front to seed the tropical computation
    # for the ayanamsha-mode sidereal Sunshine 'i'.
    hsys_char = _hsys_to_char(hsys)

    # In an ayanamsha (non-fixed-epoch) sidereal zodiac the reference computes
    # Sunshine 'i' (Makransky) identically to 'I' (Treindl): the Makransky
    # upper/lower-meridian construction — and its circumpolar guard — is not
    # applied on that path. Seed the tropical computation from 'I' so a
    # circumpolar Sun returns cusps (matching the reference) instead of raising.
    # The fixed-epoch modes keep 'i' on the Makransky path, where they match the
    # reference by raising, so they are deliberately excluded here.
    _tropical_hsys = hsys
    if (flags & FLG_SIDEREAL) and hsys_char == "i":
        from .state import get_sid_mode as _get_sid_mode_seed
        from .sidereal_epoch import FIXED_EPOCH_T0 as _FIXED_EPOCH_T0_seed

        if _get_sid_mode_seed() not in _FIXED_EPOCH_T0_seed:
            _tropical_hsys = ord("I")

    # Propagate flags to houses() so ephemeris flags (FLG_MOSEPH etc.) are used
    cusps, ascmc = houses(tjdut, lat, lon, _tropical_hsys, flags)

    if flags & FLG_SIDEREAL:
        # Fixed-epoch modes (SIDM_J2000/J1900/B1950) are frame requests:
        # houses are computed on the mean ecliptic of the mode's epoch t0,
        # not by an ayanamsha shift (see _houses_fixed_epoch_sidereal).
        from .state import get_sid_mode as _get_sid_mode_hx
        from .sidereal_epoch import FIXED_EPOCH_T0

        _sidm_hx = _get_sid_mode_hx()
        if _sidm_hx in FIXED_EPOCH_T0:
            _hch = _hsys_to_char(hsys)
            cusps, ascmc = _houses_fixed_epoch_sidereal(
                tjdut, lat, _hch, cusps, ascmc, flags, _sidm_hx
            )
            cusps = tuple(degnorm(c) for c in cusps)
            ascmc = tuple(degnorm(a) for a in ascmc)
            if flags & FLG_RADIANS:
                cusps = tuple(math.radians(c) for c in cusps)
                ascmc = tuple(math.radians(a) for a in ascmc)
            return cusps, ascmc

        # SIDBIT_ECL_T0 / SIDBIT_SSY_PLANE frame projections. The sidereal zodiac
        # is measured on the mean ecliptic of the mode's reference epoch t0, or on
        # the solar-system invariable plane, instead of the ecliptic of date.
        # Compatibility contract: the whole house construction is rebuilt on that
        # plane (the reported ARMC itself shifts to the plane's node), so route
        # to the projected-plane reconstruction. The star/galactic "true"
        # modes define their zero point on the ecliptic of date, where the
        # projection is inert (the same suppression the calc path applies): those
        # fall through to the ayanamsha path below, which reproduces the base
        # longitude.
        from .planets import (
            _SIDBIT_PROJECTION_SUPPRESS_MODES as _sidbit_suppress,
            _get_sidereal_bits,
        )

        _bits_hx = _get_sidereal_bits()
        if (_bits_hx & (SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE)) and (
            _sidm_hx not in _sidbit_suppress
        ):
            cusps, ascmc = _houses_sidbit_projection(
                tjdut, lat, hsys_char, ascmc, flags, _sidm_hx, _bits_hx
            )
            cusps = tuple(degnorm(c) for c in cusps)
            ascmc = tuple(degnorm(a) for a in ascmc)
            if flags & FLG_RADIANS:
                cusps = tuple(math.radians(c) for c in cusps)
                ascmc = tuple(math.radians(a) for a in ascmc)
            return cusps, ascmc

        # Sidereal house cusps use the MEAN-equinox ayanamsha (no nutation
        # term — houses are geometric ARMC-frame quantities). FLG_TRUEPOS /
        # FLG_NOABERR / FLG_NOGDEFL are honoured: for the star/galactic-center
        # anchored modes the anchor's annual aberration and solar light
        # deflection are removed from the subtracted value (compatibility
        # contract: every cusp and angle shifts by the same amount as the
        # calc sidereal path under those bits). Dropping FLG_NOGDEFL here
        # left the cusps completely unchanged by that flag, when they should
        # all shift by the anchor's deflection delta (-0.045" for a True
        # Citra solar-conjunction case).
        from .constants import FLG_NOABERR, FLG_NOGDEFL, FLG_TRUEPOS
        from .planets import get_ayanamsa_ex_ut

        ayanamsa = get_ayanamsa_ex_ut(
            tjdut, flags & (FLG_TRUEPOS | FLG_NOABERR | FLG_NOGDEFL)
        )[1]

        # Compute sidereal angles
        # All ecliptic longitudes in ascmc get ayanamsa correction EXCEPT
        # ascmc[2] (ARMC) which is an equatorial/geometric quantity.
        ascmc_list = list(ascmc)
        sid_asc = (ascmc_list[0] - ayanamsa) % 360.0
        sid_mc = (ascmc_list[1] - ayanamsa) % 360.0
        ascmc_list[0] = sid_asc
        ascmc_list[1] = sid_mc
        # ascmc[2] = ARMC — geometric, NOT corrected
        ascmc_list[3] = (ascmc_list[3] - ayanamsa) % 360.0  # Vertex
        ascmc_list[4] = (ascmc_list[4] - ayanamsa) % 360.0  # EquAsc
        ascmc_list[5] = (ascmc_list[5] - ayanamsa) % 360.0  # CoAsc Koch
        ascmc_list[6] = (ascmc_list[6] - ayanamsa) % 360.0  # CoAsc Munkasey
        ascmc_list[7] = (ascmc_list[7] - ayanamsa) % 360.0  # PolarAsc
        ascmc = tuple(ascmc_list)

        # For Ascendant-based house systems, recalculate using sidereal Ascendant
        # Whole Sign (W), Equal (A/E), Vehlow (V)
        if hsys_char == "W":
            # Whole Sign: each house is 0° of consecutive signs from Asc sign
            # cusps[0] = House 1, cusps[1] = House 2, etc.
            start = math.floor(sid_asc / 30.0) * 30.0
            cusps = tuple([(start + i * 30.0) % 360.0 for i in range(12)])
        elif hsys_char in ("A", "E"):
            # Equal: 30° divisions starting from sidereal Asc
            cusps = tuple([(sid_asc + i * 30.0) % 360.0 for i in range(12)])
        elif hsys_char == "V":
            # Vehlow Equal: sidereal Asc at middle of 1st house
            start = (sid_asc - 15.0) % 360.0
            cusps = tuple([(start + i * 30.0) % 360.0 for i in range(12)])
        elif hsys_char == "N":
            # Aries houses: the wheel is anchored at 0 deg of the
            # zodiac in use — in the sidereal zodiac the cusps stay at
            # 0, 30, ... (the reference does not shift them by the
            # ayanamsha).
            cusps = tuple(float(i * 30) for i in range(12))
        else:
            # For every other system — including Sunshine 'I'/'i', whose
            # tropical Treindl cusps already carry the single MC-below-horizon
            # reflection — subtract the ayanamsha uniformly from all twelve
            # cusps. This mirrors the reference (compute tropical, shift by the
            # ayanamsha) and avoids re-running the Sunshine construction on an
            # already-reflected sidereal MC, which would flip it a second time.
            cusps = tuple([(c - ayanamsa) % 360.0 for c in cusps])

    # degnorm snaps the bare-%360 artifact (exactly 360.0 from a
    # tiny-negative angle, e.g. after the ayanamsha subtraction) to 0.0.
    cusps = tuple(degnorm(c) for c in cusps)
    ascmc = tuple(degnorm(a) for a in ascmc)
    # FLG_RADIANS converts every angular position output, including Gauquelin
    # sectors and the ASCMC slots.
    if flags & FLG_RADIANS:
        cusps = tuple(math.radians(c) for c in cusps)
        ascmc = tuple(math.radians(a) for a in ascmc)
    return cusps, ascmc


def houses_ex2(
    tjdut: float, lat: float, lon: float, hsys: int = ord("P"), flags: int = 0
) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    """
    Extended house calculation returning cusps, angles, and their velocities.

    This function is an extended version of houses_ex() that also returns
    the velocities (derivatives) of house cusps and angles.

    Velocities are always calculated. This is useful for progressed chart applications where
    the rate of change of house cusps is needed.

    Velocities are computed via centered finite differences of houses_ex(),
    so the reported rates include the same flag-dependent frame corrections
    (e.g. the FLG_SIDEREAL ayanamsa) and time-dependent terms (ARMC rate,
    obliquity drift, nutation) as the returned cusp and angle positions.

    Args:
        tjdut: Julian Day in Universal Time (UT1)
        lat: Geographic latitude in degrees
        lon: Geographic longitude in degrees
        hsys: House system identifier (int or bytes)
        flags: Calculation flags bitmask. FLG_SIDEREAL can be used for
               sidereal calculations. Velocities are always computed
               regardless of FLG_SPEED.

    Returns:
        Tuple containing:
            - cusps: Tuple of 12 house cusp longitudes in degrees
            - ascmc: Tuple of 8 angles (Asc, MC, etc.)
            - cusps_speed: Tuple of 12 house cusp velocities in degrees/day
            - ascmc_speed: Tuple of 8 angle velocities in degrees/day

    Example:
        >>> cusps, ascmc, cusps_speed, ascmc_speed = houses_ex2(
        ...     2451545.0, 41.9, 12.5, ord('P'), FLG_SPEED
        ... )
        >>> # cusps_speed[0] is the velocity of the 1st house cusp (ASC)
    """
    # Calculate positions at current time
    # Work in degrees internally: FLG_RADIANS is stripped here and applied
    # once to the final position tuples below, so the finite-difference
    # speed stencil stays in degrees. The speed tuples are documented in
    # degrees per day regardless of the position-output unit flag.
    _flags_deg = flags & ~FLG_RADIANS
    cusps, ascmc = houses_ex(tjdut, lat, lon, hsys, _flags_deg)

    # Velocities (daily motion). A cusp or angle longitude is a function of time
    # through the sidereal time (ARMC), the obliquity of the ecliptic and (via
    # the ecliptic frame) nutation — plus, under FLG_SIDEREAL, the ayanamsa that
    # houses() does not apply. The true daily motion is the TOTAL time
    # derivative dλ/dt, obtained by finite-differencing the full houses_ex()
    # solution in time so every time-dependent term is captured automatically.
    # This is the genuine derivative of the reported positions for every system,
    # including the iteratively-solved ones (Placidus, Koch): integrating it
    # reproduces the cusp motion, which an analytic speed approximation of those
    # systems does not.
    #
    # Step choice is a truncation-vs-roundoff trade-off. A centered difference
    # carries O(h²) truncation ≈ f'''·h²/6 and O(ε/h) roundoff, where ε is the
    # position noise floor. The coarse half-step 1/4096 day (≈ 21 s) is an exact
    # binary fraction, so tjdut ± dt is exactly representable and the stencil
    # stays symmetric; it sits near the roundoff floor for the slowly-curving
    # cusp longitudes and for MC/ARMC. (A one-second half-step is unusable: at
    # JD ~2.45e6 the ULP is ~46 µs, so it loses ~5 significant digits to
    # roundoff and biases every rate low by several arcsec/day.)
    #
    # The angles carry a large f''' on the latitude-amplified Ascendant, Vertex
    # and co-ascendant slots. There the fixed-1/4096 truncation is not
    # negligible: it reaches a few arcsec/day at mid-latitude and hundreds of
    # arcsec/day near the polar circle, while MC/ARMC stay at the roundoff floor.
    # Refining the step alone cannot fix this — a step short enough to tame the
    # fast angles pushes the slow ones into roundoff. The ascmc rates therefore
    # use two-step Richardson extrapolation, combining centered differences at
    # h = 1/4096 and h/2 = 1/8192 day (both exact binary fractions, both above
    # the roundoff floor) to cancel the leading f'''·h² term, O(h²) → O(h⁴).
    # The cusp rates keep the plain centered difference at the coarse step.
    _DT_DAYS = 1.0 / 4096.0
    _DT_HALF = 1.0 / 8192.0
    cusps_minus, ascmc_minus = houses_ex(tjdut - _DT_DAYS, lat, lon, hsys, _flags_deg)
    cusps_plus, ascmc_plus = houses_ex(tjdut + _DT_DAYS, lat, lon, hsys, _flags_deg)
    # Fine ± h/2 samples feed only the Richardson step for the angles.
    _, ascmc_hminus = houses_ex(tjdut - _DT_HALF, lat, lon, hsys, _flags_deg)
    _, ascmc_hplus = houses_ex(tjdut + _DT_HALF, lat, lon, hsys, _flags_deg)

    _is_sign_pinned = _hsys_to_char(hsys) in _SIGN_PINNED_HSYS

    def _rate(after: float, before: float) -> float:
        d = after - before
        if d > 180.0:
            d -= 360.0
        elif d < -180.0:
            d += 360.0
        # A change larger than 90 degrees across the stencil is a branch
        # discontinuity, where no finite derivative exists.
        if abs(d) > 90.0:
            return 0.0
        # Sign-pinned systems (Whole Sign, Aries) hold their cusps on exact
        # 30-degree boundaries, so the cusp function is piecewise constant and
        # its derivative is zero — except across the instant the Ascendant
        # changes sign, where the whole wheel steps by exactly 30 degrees.
        # That step is far below the 90-degree guard above, so the stencil
        # reported it as a rate: 30 / (2 * 1/4096) = 61440 deg/day on all
        # twelve cusps, hit by ~0.5% of random (jd, lat, lon) triples.
        if _is_sign_pinned and abs(abs(d) - 30.0) < 1e-6:
            return 0.0
        return d / (2.0 * _DT_DAYS)

    def _rate_richardson(ph: float, mh: float, pq: float, mq: float) -> float:
        # Wrapped central differences at the coarse (± h) and fine (± h/2)
        # symmetric steps, then Richardson-extrapolated to cancel the O(h²) term.
        d_h = ph - mh
        if d_h > 180.0:
            d_h -= 360.0
        elif d_h < -180.0:
            d_h += 360.0
        d_q = pq - mq
        if d_q > 180.0:
            d_q -= 360.0
        elif d_q < -180.0:
            d_q += 360.0
        # A change larger than 90 degrees across either symmetric step is a
        # branch discontinuity, where no finite derivative exists.
        if abs(d_h) > 90.0 or abs(d_q) > 90.0:
            return 0.0
        d_coarse = d_h / (2.0 * _DT_DAYS)
        d_fine = d_q / (2.0 * _DT_HALF)
        # Cancel the leading f'''·h²/6 term: D = [4·D(h/2) − D(h)] / 3.
        return (4.0 * d_fine - d_coarse) / 3.0

    cusps_speed = tuple(_rate(cusps_plus[i], cusps_minus[i]) for i in range(len(cusps)))
    ascmc_speed = tuple(
        _rate_richardson(ascmc_plus[i], ascmc_minus[i], ascmc_hplus[i], ascmc_hminus[i])
        for i in range(len(ascmc))
    )

    # Cusp speeds are ALWAYS the finite-difference time derivative of the
    # reported cusp positions, for every house system: the speed contract is
    # "derivative of what this function returns". An external implementation
    # reports, for Porphyry/Whole-Sign/Aries, analytic values that are not the
    # derivative of its own cusps (see docs/comparison/intentional-divergences
    # "Cusp and angle speeds"); that behavioral rule is not reproducible from
    # any published definition and is intentionally not replicated.

    # degnorm on the angle outputs (not the speeds) snaps the bare-%360
    # artifact (exactly 360.0 from a tiny-negative angle) back to 0.0.
    cusps = tuple(degnorm(c) for c in cusps)
    ascmc = tuple(degnorm(a) for a in ascmc)
    # FLG_RADIANS converts position tuples only; speed tuples retain their
    # documented degree-per-day unit.
    if flags & FLG_RADIANS:
        cusps = tuple(math.radians(c) for c in cusps)
        ascmc = tuple(math.radians(a) for a in ascmc)
    return (cusps, ascmc, cusps_speed, ascmc_speed)


def _houses_with_context(
    tjdut: float, lat: float, lon: float, hsys: int, ctx
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """
    Calculate houses using an explicit EphemerisContext (thread-safe).

    Thread-safe wrapper around houses that uses context state.

    Args:
        tjdut: Julian Day in Universal Time
        lat: Geographic latitude in degrees
        lon: Geographic longitude in degrees
        hsys: House system identifier
        ctx: EphemerisContext instance

    Returns:
        Same as houses: (cusps, ascmc)

    Thread Safety:
        This function acquires state._CONTEXT_SWAP_LOCK to ensure that the
        save-set-restore cycle is atomic across threads.
    """
    from . import state

    with state._CONTEXT_SWAP_LOCK:
        # Save current global state
        old_sid_mode = state._SIDEREAL_MODE
        old_sid_t0 = state._SIDEREAL_T0
        old_sid_ayan_t0 = state._SIDEREAL_AYAN_T0

        try:
            # Temporarily set global state from context
            state._SIDEREAL_MODE = ctx.sidereal_mode
            state._SIDEREAL_T0 = ctx.sidereal_t0
            state._SIDEREAL_AYAN_T0 = ctx.sidereal_ayan_t0

            # Use existing house calculation logic
            return houses(tjdut, lat, lon, hsys)
        finally:
            # Restore global state
            state._SIDEREAL_MODE = old_sid_mode
            state._SIDEREAL_T0 = old_sid_t0
            state._SIDEREAL_AYAN_T0 = old_sid_ayan_t0


def house_name(hsys: int) -> str:
    """
    Get the name of a house system.

    The labels are the :data:`~libephemeris.compat_names.HOUSE_SYSTEM_NAMES`
    table; this function only normalises the selector.

    Args:
        hsys: House-system selector as an int code (``ord("P")``), a
            one-byte ``bytes`` or a one-character ``str``. Lowercase letters
            are accepted like their uppercase form, except ``'i'``, which
            is a system of its own.

    Returns:
        The house-system label, or the empty string for a selector the
        table does not list (an int outside the ``chr()`` range included).
    """
    hsys_char = _hsys_to_char(hsys)
    # Compatibility contract: an unknown selector yields an empty
    # string, and the name lookup folds 'g' to Gauquelin (unlike the
    # houses() tuple shape, which stays 12 for the lowercase byte).
    if hsys_char == "g":
        hsys_char = "G"
    return HOUSE_SYSTEM_NAMES.get(hsys_char, "")


def _houses_placidus(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Placidus house system (time-based divisions of diurnal/nocturnal arcs).

    Most popular house system in modern Western astrology. Divides the time a point
    takes to travel from horizon to meridian (and meridian to horizon) into thirds.

    Historical basis: Placidus de Titis, "Primum Mobile" (1657); modern
    description in Holden, "The Elements of House Division" (1977).

    Algorithm:
        1. Trisect semi-diurnal arc (rising to culmination) for houses 11, 12
        2. Trisect semi-nocturnal arc (setting to anti-culmination) for houses 2, 3
        3. Use iterative solution to find ecliptic longitude at each time division
        4. Calculate opposite houses by adding 180°

    Mathematical Formulas:
        Ascensional Difference: AD = arcsin(tan(φ) · tan(δ))
        Semi-diurnal Arc: SA = 90° + AD (above horizon)
        Semi-nocturnal Arc: SA = 90° - AD (below horizon)

        House 11: H₁₁ = SA/3 = (90° + AD)/3; RA₁₁ = ARMC + H₁₁
        House 12: H₁₂ = 2·SA/3 = 2·(90° + AD)/3; RA₁₂ = ARMC + H₁₂
        House 2:  H₂ = 2·(90° - AD)/3; RA₂ = ARMC + 180° - H₂
        House 3:  H₃ = (90° - AD)/3; RA₃ = ARMC + 180° - H₃

        Iterative convergence (threshold 1e-7°):
            1. tan(δ) = sin(RA) · tan(ε)
            2. AD = arcsin(tan(φ) · tan(δ))
            3. RA_new = ARMC ± H(AD)
            4. Repeat until |RA_new - RA| < 1e-7°

        RA to ecliptic: λ = atan2(sin(RA), cos(RA) · cos(ε))

    Note: Polar latitude handling
        Placidus is undefined at latitudes > ~66° where some ecliptic points never
        rise/set (circumpolar behavior). At polar latitudes, houses() raises
        PolarCircleError with detailed information about the threshold and
        recommended alternatives. Use houses_with_fallback() for automatic
        fallback to Porphyry. The polar threshold is approximately 90° - obliquity.
        Internal fallback to Porphyry is retained as a safety net.

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees (house 1)
        mc: Midheaven longitude in degrees (house 10)

    Returns:
        List of 13 house cusp longitudes (index 0 is 0.0, indices 1-12 are cusps)
    """
    cusps = _init_cardinal_cusps(asc, mc)

    rad_lat = math.radians(lat)
    rad_eps = math.radians(eps)

    def iterate_placidus(offset_deg, is_below_horizon):
        """Solve one non-cardinal Placidus cusp by fixed-point iteration.

        ``offset_deg`` is only a stable internal selector and initial RA
        displacement: 30, 60, 120, and 150 identify cusps 11, 12, 2, and 3.
        ``is_below_horizon`` records whether the selected cusp belongs to the
        nocturnal half; the explicit selector table below determines the exact
        fraction and sign.  The defining equations are those in the parent
        docstring, not empirical offsets inferred from compatibility output.
        """

        ra = (armc + offset_deg) % 360.0
        new_ra = ra  # Initialize for type safety
        prev_diff = float("inf")  # Track convergence/divergence
        converged = False

        # Convergence threshold: 1e-7° = 0.00036 arcsec
        CONVERGENCE_THRESHOLD = 1e-7
        # Maximum iterations - if not converged by then, likely diverging
        MAX_ITERATIONS = 50

        # Iterate to convergence (typically 15-20 iterations needed for 1e-7° threshold)
        for iteration in range(MAX_ITERATIONS):
            # Calculate declination for point at this RA on ecliptic
            # Using spherical astronomy formula: tan(dec) = sin(RA) * tan(eps)

            sin_ra = math.sin(math.radians(ra))
            tan_dec = sin_ra * math.tan(rad_eps)

            # Calculate semi-arc (or part of it)
            # tan(lat) * tan(dec)
            # Within one float-noise step of the polar circle the
            # candidate RA can sweep ecliptic degrees whose declination
            # makes |tan(lat)*tan(dec)| exceed 1 mid-iteration even
            # though the converged cusp itself is well-defined.  Clamp
            # the ascensional-difference argument (AD saturates at
            # +/-90 deg) and keep iterating — the latitude pre-check
            # has already raised for genuinely polar latitudes.
            # Bailing out here sent these cases to the fallback cusps,
            # up to 20 deg of divergence from one external implementation
            # at lat +/-66.55.
            prod = math.tan(rad_lat) * tan_dec
            if prod > 1.0:
                prod = 1.0
            elif prod < -1.0:
                prod = -1.0

            ad = math.asin(prod)
            ad_deg = math.degrees(ad)

            if offset_deg == 30:  # House 11
                # One third of the semi-diurnal arc, east of the upper
                # meridian: alpha_11 = ARMC + (90 deg + AD) / 3.
                h_deg = (90.0 + ad_deg) / 3.0
                new_ra = (armc + h_deg) % 360.0

            elif offset_deg == 60:  # House 12
                # Two thirds of the same semi-diurnal arc.
                h_deg = 2.0 * (90.0 + ad_deg) / 3.0
                new_ra = (armc + h_deg) % 360.0

            elif offset_deg == 120:  # House 2
                # Two thirds of the semi-nocturnal arc measured back from
                # the lower meridian (IC).
                h_deg = 2.0 * (90.0 - ad_deg) / 3.0
                new_ra = (armc + 180.0 - h_deg) % 360.0

            elif offset_deg == 150:  # House 3
                # One third of the same semi-nocturnal arc from the IC.
                # Keep the historical operation order literally unchanged:
                # multiplying by 1.0 is algebraically redundant, but retaining
                # it makes this documentation-only edit byte-for-byte neutral
                # at the floating-point expression level.
                h_deg = 1.0 * (90.0 - ad_deg) / 3.0
                new_ra = (armc + 180.0 - h_deg) % 360.0

            else:
                # Should never reach here - all valid offsets are handled above
                new_ra = ra  # Fallback to prevent unbound variable error

            # Update RA
            diff = abs(new_ra - ra)
            if diff > 180:
                diff = 360 - diff
            ra = new_ra

            # Check for convergence
            if diff < CONVERGENCE_THRESHOLD:
                converged = True
                break

            # Check for divergence: if diff is increasing, iteration is unstable
            # Allow some oscillation (diff can increase slightly) but detect runaway
            if iteration > 5 and diff > prev_diff * 1.5:
                # Diverging - this happens near polar latitudes
                return None

            prev_diff = diff

        # If we exhausted iterations without converging, return None for fallback
        if not converged:
            return None

        # Converged RA. Find Ecliptic Longitude.
        # tan(lon) = tan(ra) / cos(eps)
        # Use atan2 to handle quadrants correctly
        # sin(lon) ~ sin(ra)
        # cos(lon) ~ cos(ra) * cos(eps)
        y = math.sin(math.radians(ra))
        x = math.cos(math.radians(ra)) * math.cos(rad_eps)
        lon = math.degrees(math.atan2(y, x))

        return lon % 360.0

    def _ad_deg_at(ra_deg: float) -> float:
        """Clamped ascensional difference (deg) of the ecliptic point at RA."""
        tan_dec = math.sin(math.radians(ra_deg)) * math.tan(rad_eps)
        prod = math.tan(rad_lat) * tan_dec
        if prod > 1.0:
            prod = 1.0
        elif prod < -1.0:
            prod = -1.0
        return math.degrees(math.asin(prod))

    def bisect_placidus(offset_deg: float, is_below_horizon: bool):
        """Solve the Placidus hour-angle condition by bisection.

        The fixed-point iteration above oscillates within float-noise of
        the polar circle (lat +/-66.55 at eps 23.4393), where the
        reference still converges.  The cusp condition is a continuous
        root problem with a guaranteed sign change on its bracket, so
        bisection always lands it:

          house 11: x = (90 + AD(armc + x)) / 3          on x in [0, 60]
          house 12: x = 2 (90 + AD(armc + x)) / 3        on x in [0, 120]
          house 2:  x = 2 (90 - AD(armc + 180 - x)) / 3  on x in [0, 120]
          house 3:  x = (90 - AD(armc + 180 - x)) / 3    on x in [0, 60]

        (x is the hour-angle fraction from the relevant meridian; AD is
        clamped so the bracket endpoints always satisfy g(lo) <= 0 and
        g(hi) >= 0.)
        """
        if offset_deg == 30:
            frac, hi, from_ic = 1.0 / 3.0, 60.0, False
        elif offset_deg == 60:
            frac, hi, from_ic = 2.0 / 3.0, 120.0, False
        elif offset_deg == 120:
            frac, hi, from_ic = 2.0 / 3.0, 120.0, True
        else:  # 150
            frac, hi, from_ic = 1.0 / 3.0, 60.0, True

        def g(x: float) -> float:
            if from_ic:
                ad = _ad_deg_at(armc + 180.0 - x)
                return x - frac * (90.0 - ad)
            ad = _ad_deg_at(armc + x)
            return x - frac * (90.0 + ad)

        lo_x, hi_x = 0.0, hi
        g_lo = g(lo_x)
        if g_lo > 0.0:
            return None  # bracket failed (cannot happen with clamped AD)
        for _ in range(80):
            mid = 0.5 * (lo_x + hi_x)
            if g(mid) <= 0.0:
                lo_x = mid
            else:
                hi_x = mid
        x = 0.5 * (lo_x + hi_x)
        ra = (armc + 180.0 - x) % 360.0 if from_ic else (armc + x) % 360.0
        y = math.sin(math.radians(ra))
        xx = math.cos(math.radians(ra)) * math.cos(rad_eps)
        return math.degrees(math.atan2(y, xx)) % 360.0

    # Calculate cusps: fast fixed point first, bisection where it fails
    c11 = iterate_placidus(30, False)
    c12 = iterate_placidus(60, False)
    c2 = iterate_placidus(120, True)
    c3 = iterate_placidus(150, True)
    if c11 is None:
        c11 = bisect_placidus(30, False)
    if c12 is None:
        c12 = bisect_placidus(60, False)
    if c2 is None:
        c2 = bisect_placidus(120, True)
    if c3 is None:
        c3 = bisect_placidus(150, True)

    if c11 is None or c12 is None or c2 is None or c3 is None:
        # Safety net only: the bisection bracket is mathematically
        # guaranteed with clamped AD.
        return _houses_porphyry(asc, mc)

    cusps[11] = c11
    cusps[12] = c12
    cusps[2] = c2
    cusps[3] = c3

    # Opposites
    cusps[5] = (c11 + 180) % 360.0
    cusps[6] = (c12 + 180) % 360.0
    cusps[8] = (c2 + 180) % 360.0
    cusps[9] = (c3 + 180) % 360.0

    return cusps


def _houses_koch(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Koch (Birthplace/GOH) house system (code 'K').

    The Koch system (Geburtsort-Häusersystem = Birthplace House System) was
    developed by Walter Koch. It divides houses based on the ascensional
    difference of the MC's declination.

    Primary basis: Walter Koch & Elisabeth Knappich, "Geburtsort-
    Häusertabellen" (Birthplace House Tables, 1960s); modern derivation in
    Holden, 'The Elements of House Division' (1977), Koch section. The Koch arc
    (one-third of the MC's ascensional difference) comes from that
    construction, not from Meeus.

    Derivation from semi-arc concept:
        1. Compute sin(dec_MC) = sin(MC) * sin(eps) (MC declination)
        2. Compute elevation angle c = atan(tan(lat) / cos(dec_MC))
        3. Koch arc = asin(sin(c) * sin(dec_MC)) / 3
           (one-third of the ascensional difference)
        4. Cusps are the ecliptic rising points at ARMC offsets of
           30, 60, 120, 150 degrees, adjusted by multiples of the Koch arc.

    Reference for the rising-point coordinate transform only: Meeus,
    'Astronomical Algorithms' 2nd ed., Ch. 13 (spherical trigonometry). The
    Koch quadrant division itself is not from Meeus.

    Note: Polar latitude handling
        Koch is undefined at latitudes > ~66 deg where some ecliptic points never
        rise/set (circumpolar behavior). At polar latitudes, houses() raises
        PolarCircleError with detailed information about the threshold and
        recommended alternatives. Use houses_with_fallback() for automatic
        fallback to Porphyry. The polar threshold is approximately 90 deg - obliquity.
        Internal fallback to Porphyry is retained as a safety net.

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = _init_cardinal_cusps(asc, mc)

    # Check for polar circle - Koch undefined at |lat| >= 90 - eps
    if abs(lat) >= 90 - eps:
        return _houses_porphyry(asc, mc)

    # Precompute trigonometric values for the obliquity and latitude
    sin_obliquity = math.sin(math.radians(eps))
    tan_lat = math.tan(math.radians(lat))
    cos_lat = math.cos(math.radians(lat))

    # Koch auxiliary: sin(dec_MC) = sin(MC) * sin(ε) / cos(φ). This belongs to
    # the Koch/Knappich house construction (see Holden 1977), NOT to Meeus;
    # Meeus Ch. 13 only supplies the rising-point transform used by
    # _calc_ascendant below.
    mc_dec_sin = max(
        -1.0, min(1.0, math.sin(math.radians(mc)) * sin_obliquity / cos_lat)
    )
    mc_dec_cos = math.sqrt(1.0 - mc_dec_sin * mc_dec_sin)  # always >= 0

    # Auxiliary angle c: elevation of the MC above the prime vertical
    # atan2 handles the degenerate case mc_dec_cos == 0 natively
    zenith_angle = math.degrees(math.atan2(tan_lat, mc_dec_cos))

    # Ascensional difference divided by 3 (used to offset each intermediate cusp)
    # asc_diff_third = asin(sin(zenith_angle) * sin(dec_MC)) / 3
    asc_diff_arg = max(
        -1.0, min(1.0, math.sin(math.radians(zenith_angle)) * mc_dec_sin)
    )
    asc_diff_third = math.degrees(math.asin(asc_diff_arg)) / 3.0

    # Calculate cusps using _calc_ascendant (ecliptic rising-point formula)
    cusps[11] = _calc_ascendant(armc + 30 - 2 * asc_diff_third, eps, lat, lat)
    cusps[12] = _calc_ascendant(armc + 60 - asc_diff_third, eps, lat, lat)
    cusps[2] = _calc_ascendant(armc + 120 + asc_diff_third, eps, lat, lat)
    cusps[3] = _calc_ascendant(armc + 150 + 2 * asc_diff_third, eps, lat, lat)

    # Opposite houses
    _set_opposite_cusps(cusps)

    return cusps


#: The equator points of the intermediate house circles, as right-ascension
#: offsets east of the meridian, keyed by the cusp they define. Cusp 10 is the
#: meridian itself (offset 0) and cusp 1 the horizon (offset 90); cusps 5, 6,
#: 8 and 9 are the antipodes of 11, 12, 2 and 3 (:func:`_set_opposite_cusps`).
_HOUSE_CIRCLE_OFFSETS: tuple[tuple[int, float], ...] = (
    (11, 30.0),
    (12, 60.0),
    (2, 120.0),
    (3, 150.0),
)


def _house_circle_rising_longitude(
    armc: float, offset_deg: float, tan_pole: float, cos_eps: float, sin_eps: float
) -> float:
    """Ecliptic longitude cut by a house circle of pole height ``P``.

    A great circle that crosses the equator ``offset_deg`` east of the
    meridian and whose pole stands at height ``P`` (``tan_pole`` = tan P) is
    the horizon of a fictitious observer: one at latitude ``P`` whose meridian
    is turned by ``offset_deg - 90`` in right ascension. The longitude the
    circle cuts on the ecliptic is that observer's Ascendant, the
    horizon/ecliptic intersection formula evaluated with ``R = ARMC + offset -
    90`` in place of the ARMC and ``P`` in place of the latitude::

        tan(lambda) = cos(R) / -(sin(R) cos(eps) + tan(P) sin(eps))

    Equivalently, the cusp is the ecliptic point whose oblique ascension under
    the pole ``P`` equals ``ARMC + offset``, which is how Holden ("The
    Elements of House Division", 1977) and Munkasey ("An Astrological House
    Formulary") state the Regiomontanus and Topocentric cusps; for
    ``offset_deg = 90`` and ``P = lat`` it is the Ascendant itself. The two
    systems differ only in how they choose ``P``.
    """
    r_rad = math.radians((armc + offset_deg - 90.0) % 360.0)
    num = math.cos(r_rad)
    den = -(math.sin(r_rad) * cos_eps + tan_pole * sin_eps)
    return math.degrees(math.atan2(num, den)) % 360.0


def _houses_regiomontanus(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Regiomontanus (Medieval rational) house system.

    Divides the celestial equator into 12 equal 30° arcs from the meridian and
    projects the divisions onto the ecliptic along the great circles through
    the north and south points of the horizon.

    Historical basis: Regiomontanus (Johannes Müller von Königsberg), "Tabulae
    directionum profectionumque" (compiled 1467; printed 1490); modern
    descriptions in Holden, "The Elements of House Division" (1977) and
    Munkasey, "An Astrological House Formulary".

    Geometry:
        The division point of cusp ``k`` sits on the equator ``H`` degrees
        east of the meridian (:data:`_HOUSE_CIRCLE_OFFSETS`). The house
        circle through it and the north and south points of the horizon has
        its pole on the prime vertical, at height ``P`` with
        ``tan(P) = tan(lat) · sin(H)``. The cusp is the longitude that
        circle cuts on the ecliptic, read as the Ascendant of the fictitious
        observer at latitude ``P`` with meridian ``ARMC + H - 90``
        (:func:`_house_circle_rising_longitude`); the opposite cusps are the
        antipodes.

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = _init_cardinal_cusps(asc, mc)

    tan_lat = math.tan(math.radians(lat))
    cos_eps = math.cos(math.radians(eps))
    sin_eps = math.sin(math.radians(eps))

    for cusp, offset_deg in _HOUSE_CIRCLE_OFFSETS:
        tan_pole = tan_lat * math.sin(math.radians(offset_deg))
        cusps[cusp] = _house_circle_rising_longitude(
            armc, offset_deg, tan_pole, cos_eps, sin_eps
        )

    _set_opposite_cusps(cusps)

    return cusps


def _houses_campanus(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Campanus (Prime Vertical) house system.

    Divides the prime vertical (great circle through Zenith, East, Nadir, West)
    into 12 equal 30° arcs, then projects onto the ecliptic via great circles
    through the North and South points of the horizon.

    Historical basis: Campanus of Novara (13th century); modern descriptions in
    Holden, "The Elements of House Division" (1977) and Munkasey, "An
    Astrological House Formulary". As in Regiomontanus the house circles run
    through the north and south points of the horizon, but the division
    points lie on the prime vertical rather than on the equator, so the pole
    of each circle is derived below instead of from ``tan(lat) · sin(H)``.

    Derivation from spherical trigonometry
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The North point of the horizon (N) is the **pole** of the prime vertical
    (PV), since it is exactly 90° from every point on the PV.

    For a division point D at angular offset *h* from the East point along the
    PV (positive towards zenith, negative towards nadir), the house circle is
    the great circle through N, D, and S.  Because N is the pole of the PV,
    the distance ND = 90° and the house circle is perpendicular to the PV at D.

    The pole of this house circle lies on the PV at 90° from D.  Converting
    the pole's horizontal coordinates (azimuth 270°, altitude 90° − h) to
    equatorial gives:

        pole_height  = arcsin(sin φ · cos h)
        H_pole       = atan2(sin h,  cos h · cos φ)

    The equator crossing of the house circle is at:

        RA_crossing  = ARMC + 90° − H_pole

    Ecliptic longitude follows from _ra_to_ecliptic_longitude().

    Ref: Smart, "Textbook on Spherical Astronomy", Ch. 3 (pole/great-circle
         geometry); Meeus, "Astronomical Algorithms", Ch. 13 (coordinate
         transforms).

    Args:
        armc: Right Ascension of the Midheaven (RAMC) in degrees.
        lat: Geographic latitude in degrees.
        eps: True obliquity of the ecliptic in degrees.
        asc: Ascendant longitude in degrees.
        mc: Midheaven longitude in degrees.

    Returns:
        List of 13 house cusp longitudes (index 0 unused).
    """
    cusps = _init_cardinal_cusps(asc, mc)

    sin_lat = math.sin(math.radians(lat))
    cos_lat = math.cos(math.radians(lat))
    sin_eps = math.sin(math.radians(eps))
    cos_eps = math.cos(math.radians(eps))

    def _cusp_from_pv_offset(h_deg: float) -> float:
        """Ecliptic longitude for a Campanus house circle at PV offset *h_deg*."""
        h = math.radians(h_deg)
        sin_h = math.sin(h)
        cos_h = math.cos(h)

        # Pole height of the house circle (declination of the circle's pole)
        pole_height = math.degrees(math.asin(max(-1.0, min(1.0, sin_lat * cos_h))))

        # Hour angle of the pole → RA of equator crossing
        h_pole = math.degrees(math.atan2(sin_h, cos_h * cos_lat))
        ra_crossing = armc + 90.0 - h_pole

        return _ra_to_ecliptic_longitude(ra_crossing, pole_height, sin_eps, cos_eps)

    # Cusp mapping: angular offset from East point along prime vertical
    #   +30° / +60° = above horizon (cusps 12, 11 — between ASC and MC)
    #   −30° / −60° = below horizon (cusps  2,  3 — between ASC and IC)
    cusps[12] = _cusp_from_pv_offset(30.0)
    cusps[11] = _cusp_from_pv_offset(60.0)
    cusps[2] = _cusp_from_pv_offset(-30.0)
    cusps[3] = _cusp_from_pv_offset(-60.0)

    _set_opposite_cusps(cusps)

    return cusps


def _houses_equal(asc: float) -> List[float]:
    """
    Equal house system (30° divisions from Ascendant).

    Simplest house system: each house is exactly 30° of ecliptic longitude.
    House 1 starts at Ascendant, each subsequent house adds 30°.

    Mathematical Formula:
        λᵢ = (λ_Asc + (i-1) × 30°) mod 360°
        where i = 1, 2, ..., 12

    Properties:
        - Works at all latitudes including polar regions
        - MC may not coincide with 10th house cusp
        - Computationally exact (no astronomical calculations)

    Args:
        asc: Ascendant longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = [0.0] * 13
    for i in range(1, 13):
        cusps[i] = (asc + (i - 1) * 30.0) % 360.0
    return cusps


def _houses_whole_sign(asc: float) -> List[float]:
    """
    Whole Sign house system (ancient Hellenistic method).

    Each house occupies one complete zodiac sign. House 1 starts at 0° of the
    sign containing the Ascendant. Used extensively in ancient astrology.

    Historical basis: Hellenistic and early Indian traditions; modern
    description in Holden, "The Elements of House Division" (1977).

    Algorithm:
        1. Find zodiac sign of Ascendant (floor(asc / 30) * 30)
        2. Each house = one complete sign (30° intervals from sign 0°)

    Mathematical Formula:
        start = floor(λ_Asc / 30°) × 30°
        λᵢ = (start + (i-1) × 30°) mod 360°
        where i = 1, 2, ..., 12

    Example:
        If Asc = 15° Taurus (45°), then:
        - House 1 starts at 30° (0° Taurus)
        - House 2 starts at 60° (0° Gemini)
        - etc.

    Properties:
        - Works at all latitudes including polar regions
        - Computationally exact (no astronomical calculations)
        - House cusps always at 0° of each sign

    Args:
        asc: Ascendant longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = [0.0] * 13
    # Start of sign containing Asc.  The 1e-9 deg nudge resolves the
    # boundary case: at cardinal ARMC the ascendant lands within one
    # float ULP of an exact sign cusp (e.g. 180 - 1.4e-13 at ARMC 90,
    # lat -66) and the reference assigns the NEW sign there.
    start = math.floor((asc + 1e-9) / 30.0) * 30.0
    for i in range(1, 13):
        cusps[i] = (start + (i - 1) * 30.0) % 360.0
    return cusps


def _houses_porphyry(asc: float, mc: float) -> List[float]:
    """
    Porphyry house system (space-based trisection).

    Divides each quadrant (Asc-MC, MC-Desc, Desc-IC, IC-Asc) into three equal
    30° sections along the ecliptic. Simple and well-defined at all latitudes.

    Algorithm:
        1. Calculate arc from Asc to MC, divide by 3
        2. Calculate arc from MC to Desc (MC+180), divide by 3
        3. Opposite houses are 180° from each other

    Mathematical Formulas:
        Quadrant MC→Asc (houses 11, 12):
            step = (λ_Asc - λ_MC) mod 360° / 3
            λ₁₁ = λ_MC + step
            λ₁₂ = λ_MC + 2·step

        Quadrant Asc→IC (houses 2, 3):
            step = (λ_IC - λ_Asc) mod 360° / 3
            λ₂ = λ_Asc + step
            λ₃ = λ_Asc + 2·step

        Opposite houses:
            λ₅ = (λ₁₁ + 180°) mod 360°
            λ₆ = (λ₁₂ + 180°) mod 360°
            λ₈ = (λ₂ + 180°) mod 360°
            λ₉ = (λ₃ + 180°) mod 360°

    Properties:
        - Works at all latitudes including polar regions
        - Used as fallback when Placidus/Koch fail
        - Computationally simple (no trigonometry)

    Args:
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = _init_cardinal_cusps(asc, mc)

    # Trisect the ecliptic arc between angles
    # Arc 10-1
    diff = (asc - mc) % 360.0
    step = diff / 3.0
    cusps[11] = (mc + step) % 360.0
    cusps[12] = (mc + 2 * step) % 360.0

    # Arc 1-4
    ic = cusps[4]
    diff = (ic - asc) % 360.0
    step = diff / 3.0
    cusps[2] = (asc + step) % 360.0
    cusps[3] = (asc + 2 * step) % 360.0

    # Opposites
    _set_opposite_cusps(cusps)

    return cusps


def _houses_sripati(asc: float, mc: float) -> List[float]:
    """
    Sripati house system (traditional Indian method).

    Creates house cusps at the midpoints of Porphyry house cusps. Each Sripati
    cusp is the midpoint between the previous and current Porphyry cusp.

    Historical basis: Jyotiṣa (Indian) bhāva-cakra tradition, associated with
    Śrīpati; modern description in Holden, "The Elements of House Division"
    (1977).

    Algorithm:
        1. Calculate Porphyry house cusps
        2. For each house i, the Sripati cusp is:
           H'[i] = (H[i-1] + H[i]) / 2  (midpoint on the ecliptic circle)

    Mathematical Formula:
        λ'ᵢ = λᵢ₋₁ + ((λᵢ - λᵢ₋₁) mod 360°) / 2

    This effectively shifts the house boundaries so that the "core" of each
    Porphyry house becomes the cusp, making the Porphyry cusp the center of
    each Sripati house.

    Properties:
        - Works at all latitudes (inherits Porphyry's polar stability)
        - Used in traditional Indian astrology (Jyotish)
        - Computationally simple (no trigonometry beyond Porphyry)

    Args:
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes (index 0 unused, 1-12 are house cusps)
    """
    # Get Porphyry cusps as base
    porphyry = _houses_porphyry(asc, mc)

    cusps = [0.0] * 13
    for i in range(1, 13):
        prev_i = 12 if i == 1 else i - 1
        # Calculate midpoint between previous and current Porphyry cusp
        # Handle wrap-around at 360 degrees
        diff = (porphyry[i] - porphyry[prev_i]) % 360.0
        cusps[i] = (porphyry[prev_i] + diff / 2.0) % 360.0

    return cusps


#: The four quadrants of a quadrant-division wheel: the angle each starts at,
#: the angle it ends at, and the two cusps built inside it, in order.
_PULLEN_QUADRANTS: tuple[tuple[int, int, int, int], ...] = (
    (10, 1, 11, 12),  # MC -> Asc: houses 11, 12
    (1, 4, 2, 3),  # Asc -> IC: houses 2, 3
    (4, 7, 5, 6),  # IC -> Desc: houses 5, 6
    (7, 10, 8, 9),  # Desc -> MC: houses 8, 9
)


def _houses_pullen_sd(asc: float, mc: float) -> List[float]:
    """
    Pullen SD (Sinusoidal Delta) house system, also known as Neo-Porphyry.

    Like Porphyry, an ecliptic division of the four quadrants between the
    angles; but the three houses of a quadrant differ by a constant delta
    instead of being equal, so that the house widths around the wheel follow
    a sine wave.

    Author definition: Walter D. Pullen, "Sinusoidal House Systems"
    (astrolog.org/astrolog/astsine.htm), Sinusoidal Delta houses, first
    implemented on 13 February 1994. For a quadrant of size ``q`` the three
    widths are ``x+n, x, x+n`` and those of the opposite quadrant
    ``x+3n, x+4n, x+3n``; the two sum constraints give ``n = (90 - q)/4``
    and ``x = (q - 30)/2``. Written with the deviation ``d = q - 90``, the
    side houses are ``30 + d/4`` and the middle house ``30 + d/2``, the same
    algebra for every quadrant.

    Narrow quadrants: Pullen notes that a quadrant of 30 degrees pinches the
    middle house to zero size and that below 30 degrees the pattern would
    need a house of negative size; his rule is to bisect such a quadrant
    between the two side houses and keep the middle house at zero. The rule
    is reached in practice: inside the polar circle the Ascendant can lie
    within 30 degrees of the Midheaven (the golden ``houses`` and
    ``housepos`` grids reach it at high latitude, and
    ``tests/test_houses/test_cov100_houses.py::test_pullen_sd_negative_middle``
    exercises it directly).

    Args:
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes (index 0 unused, 1-12 are house cusps)
    """
    cusps = _init_cardinal_cusps(asc, mc)

    def house_widths(quadrant_size: float) -> tuple[float, float]:
        """(side, middle) house widths of a quadrant under Pullen's delta rule."""
        d = quadrant_size - 90.0
        middle = 30.0 + d / 2.0
        if middle < 0:
            # Pullen's narrow-quadrant rule: the middle house is pinched to
            # zero and the quadrant bisected between the two side houses.
            return quadrant_size / 2.0, 0.0
        return 30.0 + d / 4.0, middle

    for start, end, first, second in _PULLEN_QUADRANTS:
        side, middle = house_widths((cusps[end] - cusps[start]) % 360.0)
        cusps[first] = (cusps[start] + side) % 360.0
        cusps[second] = (cusps[first] + middle) % 360.0

    return cusps


def _houses_pullen_sr(asc: float, mc: float) -> List[float]:
    """
    Pullen SR (Sinusoidal Ratio) house system.

    Proposed by Walter Pullen in 2016 as an improvement over Pullen SD. Uses
    ratio multipliers instead of additive offsets, avoiding negative house
    sizes for small quadrants.

    Reference: Walter D. Pullen, Astrolog (open-source software) documentation,
    "House Systems" (astrolog.org); SR introduced 2016.

    Algorithm:
        For quadrant size q:
        - Small quadrant houses: rx, x, rx  (sizes sum to q)
        - Large quadrant houses: r³x, r⁴x, r³x  (sizes sum to 180-q)

        Solving: rx + x + rx = q  →  x(2r + 1) = q
                 r³x + r⁴x + r³x = 180 - q  →  x·r³(2 + r) = 180 - q

        Dividing: r³(2 + r) / (2r + 1) = (180 - q) / q
        This is solved numerically for r.

    Properties:
        - Works at all latitudes
        - No negative house sizes even for extreme quadrants
        - More elegant mathematical formulation

    Args:
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes (index 0 unused, 1-12 are house cusps)
    """
    cusps = _init_cardinal_cusps(asc, mc)

    # Calculate quadrant sizes
    q1 = (asc - mc) % 360.0  # MC to Asc
    ic = cusps[4]
    q2 = (ic - asc) % 360.0  # Asc to IC
    desc = cusps[7]
    q3 = (desc - ic) % 360.0  # IC to Desc
    q4 = (mc - desc) % 360.0  # Desc to MC

    def calc_house_sizes_ratio(quadrant_size: float) -> tuple:
        """Calculate house sizes using sinusoidal ratio method."""
        q = quadrant_size
        q_opp = 180.0 - q  # opposite quadrant

        if abs(q - 90.0) < 0.0001:
            # Quadrant is exactly 90°, r = 1, equal houses of 30°
            return (30.0, 30.0, 30.0)

        # Target ratio: r³(2 + r) / (2r + 1) = q_opp / q
        target = q_opp / q if q > 0.0001 else 1000.0

        # Solve for r using Newton-Raphson
        # f(r) = r³(2 + r) / (2r + 1) - target = 0
        r = 1.0  # initial guess
        for _ in range(20):  # max iterations
            numerator = r * r * r * (2.0 + r)
            denominator = 2.0 * r + 1.0
            f = numerator / denominator - target

            # Derivative: d/dr[r³(2+r)/(2r+1)]
            # Using quotient rule
            dn = 3.0 * r * r * (2.0 + r) + r * r * r  # derivative of numerator
            dd = 2.0  # derivative of denominator
            df = (dn * denominator - numerator * dd) / (denominator * denominator)

            if abs(df) < 1e-12:
                break
            r_new = r - f / df
            if r_new < 0.01:
                r_new = 0.01  # keep r positive
            if abs(r_new - r) < 1e-10:
                break
            r = r_new

        # Calculate x from: x(2r + 1) = q
        x = q / (2.0 * r + 1.0) if (2.0 * r + 1.0) > 0.0001 else q / 3.0

        # House sizes: rx, x, rx
        side = r * x
        middle = x

        return (side, middle, side)

    # Quadrant 1: MC to Asc (houses 11, 12)
    sizes1 = calc_house_sizes_ratio(q1)
    cusps[11] = (mc + sizes1[0]) % 360.0
    cusps[12] = (cusps[11] + sizes1[1]) % 360.0

    # Quadrant 2: Asc to IC (houses 2, 3)
    sizes2 = calc_house_sizes_ratio(q2)
    cusps[2] = (asc + sizes2[0]) % 360.0
    cusps[3] = (cusps[2] + sizes2[1]) % 360.0

    # Quadrant 3: IC to Desc (houses 5, 6)
    sizes3 = calc_house_sizes_ratio(q3)
    cusps[5] = (ic + sizes3[0]) % 360.0
    cusps[6] = (cusps[5] + sizes3[1]) % 360.0

    # Quadrant 4: Desc to MC (houses 8, 9)
    sizes4 = calc_house_sizes_ratio(q4)
    cusps[8] = (desc + sizes4[0]) % 360.0
    cusps[9] = (cusps[8] + sizes4[1]) % 360.0

    return cusps


def _houses_alcabitius(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Alcabitius (Alchabitius) house system (ancient Arabic method).

    Medieval Arabic system that divides the diurnal and nocturnal arcs differently
    than Placidus, using a simpler geometric approach.

    Historical basis: al-Qabīṣī (Alcabitius, 10th century). Text: al-Qabīṣī,
    "The Introduction to Astrology", ed. Burnett, Yamamoto & Yano (2004).

    Algorithm:
        1. Calculate RA of Ascendant and MC
        2. Divide RA intervals between angles
        3. Convert RA divisions back to ecliptic longitude

    Mathematical Formulas:
        Right Ascension of Ascendant:
            RA_Asc = atan2(sin(λ_Asc)·cos(ε), cos(λ_Asc))

        Quadrant MC→Asc (houses 11, 12):
            arc = (RA_Asc - ARMC) mod 360°
            step = arc / 3
            RA₁₁ = ARMC + step
            RA₁₂ = ARMC + 2·step

        Quadrant Asc→IC (houses 2, 3):
            RA_IC = (ARMC + 180°) mod 360°
            arc = (RA_IC - RA_Asc) mod 360°
            step = arc / 3
            RA₂ = RA_Asc + step
            RA₃ = RA_Asc + 2·step

        RA to ecliptic conversion:
            λ = atan2(sin(RA), cos(RA)·cos(ε))

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    # Alcabitius
    # Time trisection of Ascendant's diurnal arc, projected by Hour Circles.
    # RA_11 = RAMC + SA/3.
    # SA = 90 + AD_asc.

    cusps = _init_cardinal_cusps(asc, mc)

    math.radians(lat)
    rad_eps = math.radians(eps)

    # RA of Ascendant
    y = math.sin(math.radians(asc)) * math.cos(rad_eps)
    x = math.cos(math.radians(asc))
    ra_asc = math.degrees(math.atan2(y, x)) % 360.0

    # Arc from MC to Asc
    arc = (ra_asc - armc) % 360.0
    step = arc / 3.0

    def get_lon_from_ra(ra):
        y = math.sin(math.radians(ra))
        x = math.cos(math.radians(ra)) * math.cos(rad_eps)
        lon = math.degrees(math.atan2(y, x))
        return lon % 360.0

    cusps[11] = get_lon_from_ra(armc + step)
    cusps[12] = get_lon_from_ra(armc + 2 * step)

    # Sector 2: Asc to IC
    ra_ic = (armc + 180.0) % 360.0
    arc2 = (ra_ic - ra_asc) % 360.0
    step2 = arc2 / 3.0

    cusps[2] = get_lon_from_ra(ra_asc + step2)
    cusps[3] = get_lon_from_ra(ra_asc + 2 * step2)

    _set_opposite_cusps(cusps)

    return cusps


#: Polich and Page choose the poles of the house circles by trisecting the
#: tangent of the latitude: the circles one house from the meridian (cusps 11
#: and 3) run under tan(P) = tan(lat)/3, those two houses away (12 and 2)
#: under 2·tan(lat)/3, and the horizon (cusps 1 and 7) under the full tan(lat).
_TOPOCENTRIC_POLE_FRACTIONS: dict[int, float] = {
    11: 1.0 / 3.0,
    12: 2.0 / 3.0,
    2: 2.0 / 3.0,
    3: 1.0 / 3.0,
}


def _houses_polich_page(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Polich-Page (Topocentric) house system.

    Published by Wendel Polich and A. P. Nelson Page, "The Topocentric System
    of Houses", Spica 3(3), 1964 (the system dates from 1961); described in
    Holden, "The Elements of House Division" (1977) and Munkasey, "An
    Astrological House Formulary".

    Geometry:
        The house circles cross the equator ``H = 30, 60, 120, 150`` degrees
        east of the meridian like the Regiomontanus ones
        (:data:`_HOUSE_CIRCLE_OFFSETS`), but their poles are not derived
        from ``H``: the authors trisect the tangent of the latitude,
        ``tan(P) = tan(lat) · k/3`` with ``k = 1`` one house from the
        meridian and ``k = 2`` two houses away
        (:data:`_TOPOCENTRIC_POLE_FRACTIONS`), so the circles no longer pass
        through the north and south points of the horizon. Each cusp is the
        longitude its circle cuts on the ecliptic, read as the Ascendant of
        the fictitious observer at latitude ``P`` with meridian
        ``ARMC + H - 90`` (:func:`_house_circle_rising_longitude`); the
        opposite cusps are the antipodes.

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = _init_cardinal_cusps(asc, mc)

    tan_lat = math.tan(math.radians(lat))
    cos_eps = math.cos(math.radians(eps))
    sin_eps = math.sin(math.radians(eps))

    for cusp, offset_deg in _HOUSE_CIRCLE_OFFSETS:
        tan_pole = tan_lat * _TOPOCENTRIC_POLE_FRACTIONS[cusp]
        cusps[cusp] = _house_circle_rising_longitude(
            armc, offset_deg, tan_pole, cos_eps, sin_eps
        )

    _set_opposite_cusps(cusps)

    return cusps


def _houses_morinus(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Morinus house system (equatorial divisions).

    Divides the celestial equator into 12 equal 30° sections starting from 0° Aries,
    then projects to ecliptic. Independent of observer location.

    Historical basis: Jean-Baptiste Morin (Morinus), "Astrologia Gallica"
    (1661), Book XVII.

    Algorithm:
        1. Divide equator into 30° RA sections from 0h RA
        2. Convert each RA to ecliptic longitude using obliquity

    Mathematical Formulas:
        For each house at RA = ARMC + (i-10)×30° where i = 10, 11, 12, 1, 2, 3:

        Projection via ecliptic poles:
            tan(λ) = tan(RA) × cos(ε)

        Using atan2 for proper quadrant:
            λ = atan2(sin(RA)·cos(ε), cos(RA))

    Properties:
        - Location-independent (ignores latitude)
        - Morinus 1st cusp differs from standard Ascendant
        - Morinus 10th cusp differs from standard MC

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    # Morinus
    # Projects Equator points (RAMC + 30) to Ecliptic via Ecliptic Poles.
    # Lon = Lon of point on Equator.
    # tan(lon) = tan(ra) * cos(eps).

    cusps = [0.0] * 13
    # Morinus Ascendant? Standard houses returns standard Asc.
    # But cusps[1] should be Morinus Ascendant (RAMC+90 projected).
    # Let's follow standard behavior: cusps array contains system cusps.

    rad_eps = math.radians(eps)

    def get_lon(ra):
        # tan(lon) = tan(ra) * cos(eps)
        y = math.sin(math.radians(ra)) * math.cos(rad_eps)
        x = math.cos(math.radians(ra))
        lon = math.degrees(math.atan2(y, x))
        return lon % 360.0

    cusps[10] = get_lon(armc)
    cusps[11] = get_lon(armc + 30)
    cusps[12] = get_lon(armc + 60)
    cusps[1] = get_lon(armc + 90)
    cusps[2] = get_lon(armc + 120)
    cusps[3] = get_lon(armc + 150)

    cusps[4] = (cusps[10] + 180) % 360.0
    cusps[5] = (cusps[11] + 180) % 360.0
    cusps[6] = (cusps[12] + 180) % 360.0
    cusps[7] = (cusps[1] + 180) % 360.0
    cusps[8] = (cusps[2] + 180) % 360.0
    cusps[9] = (cusps[3] + 180) % 360.0

    return cusps


def _houses_meridian(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Meridian (Zariel/Axial Rotation) house system.

    Based on meridian passages, divides RA from MC in equal 30° intervals.
    Related to Morinus but starts from MC instead of 0° Aries.

    Historical basis: proposed by "Zariel" (David Cope), c. 1910; modern
    description in Holden, "The Elements of House Division" (1977).

    Algorithm:
        Projects equator points (ARMC + n×30°) to ecliptic via celestial poles.

    Mathematical Formulas:
        For each house at RA = ARMC + (i-10)×30° where i = 10, 11, 12, 1, 2, 3:

        Projection via celestial poles:
            tan(λ) = tan(RA) / cos(ε)

        Using atan2 for proper quadrant:
            λ = atan2(sin(RA), cos(RA)·cos(ε))

    Difference from Morinus:
        Morinus:  λ = atan2(sin(RA)·cos(ε), cos(RA))
        Meridian: λ = atan2(sin(RA), cos(RA)·cos(ε))

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    # Meridian (Axial)
    # Projects Equator points to Ecliptic via Celestial Poles.
    # tan(lon) = tan(ra) / cos(eps).

    cusps = [0.0] * 13

    rad_eps = math.radians(eps)

    def get_lon_from_ra(ra):
        y = math.sin(math.radians(ra))
        x = math.cos(math.radians(ra)) * math.cos(rad_eps)
        lon = math.degrees(math.atan2(y, x))
        return lon % 360.0

    cusps[10] = get_lon_from_ra(armc)
    cusps[11] = get_lon_from_ra(armc + 30)
    cusps[12] = get_lon_from_ra(armc + 60)
    cusps[1] = get_lon_from_ra(armc + 90)
    cusps[2] = get_lon_from_ra(armc + 120)
    cusps[3] = get_lon_from_ra(armc + 150)

    cusps[4] = (cusps[10] + 180) % 360.0
    cusps[5] = (cusps[11] + 180) % 360.0
    cusps[6] = (cusps[12] + 180) % 360.0
    cusps[7] = (cusps[1] + 180) % 360.0
    cusps[8] = (cusps[2] + 180) % 360.0
    cusps[9] = (cusps[3] + 180) % 360.0

    return cusps


def _houses_vehlow(asc: float) -> List[float]:
    """
    Vehlow house system (Equal with Asc in middle of House 1).

    Variant of Equal houses where the Ascendant falls at 15° into House 1
    rather than at the cusp. Each house is still 30°.

    Historical basis: Johannes Vehlow (German astrologer), "Erlebte
    Sternenwelt"; the Ascendant is centered in House 1.

    Mathematical Formula:
        start = (λ_Asc - 15°) mod 360°
        λᵢ = (start + (i-1) × 30°) mod 360°
        where i = 1, 2, ..., 12

    Properties:
        - Ascendant falls exactly at 15° into House 1
        - Works at all latitudes including polar regions
        - Computationally exact (no astronomical calculations)

    Args:
        asc: Ascendant longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = [0.0] * 13
    start = (asc - 15.0) % 360.0
    for i in range(1, 13):
        cusps[i] = (start + (i - 1) * 30.0) % 360.0
    return cusps


def _houses_carter(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Carter Poli-Equatorial house system.

    Equal 30° divisions on the celestial equator starting from RA of Ascendant,
    projected to ecliptic via hour circles.

    Historical basis: Charles E. O. Carter, "Essays on the Foundations of
    Astrology" (Theosophical Publishing House, London).

    Algorithm:
        1. Calculate RA of Ascendant
        2. Add 30° RA increments for each house
        3. Convert each RA to ecliptic longitude

    Mathematical Formulas:
        Right Ascension of Ascendant:
            RA_Asc = atan2(sin(λ_Asc)·cos(ε), cos(λ_Asc))

        For each house i (1-12):
            RA = RA_Asc + (i-1) × 30°

        RA to ecliptic conversion via celestial poles:
            λ = atan2(sin(RA), cos(RA)·cos(ε))

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """

    cusps = [0.0] * 13

    rad_eps = math.radians(eps)

    # Get RA of Ascendant
    y = math.sin(math.radians(asc)) * math.cos(rad_eps)
    x = math.cos(math.radians(asc))
    ra_asc = math.degrees(math.atan2(y, x)) % 360.0

    def get_lon_from_ra(ra):
        y = math.sin(math.radians(ra))
        x = math.cos(math.radians(ra)) * math.cos(rad_eps)
        lon = math.degrees(math.atan2(y, x))
        return lon % 360.0

    # Equal 30-degree divisions from RA of Asc
    for i in range(1, 13):
        ra = (ra_asc + (i - 1) * 30.0) % 360.0
        cusps[i] = get_lon_from_ra(ra)

    return cusps


def _gauquelin_cusp_for_sector(
    sector_offset: int,
    base_ramc: float,
    ramc_sign: float,
    ascensional_diff: float,
    tan_obliquity: float,
    sin_obliquity: float,
    tan_lat: float,
    eps: float,
    lat: float,
    niter_max: int,
    convergence_threshold: float,
    near_zero: float,
) -> float:
    """
    Compute one intermediate Gauquelin sector cusp via iterative pole-height refinement.

    The algorithm projects a horizon-fraction (sector_offset/9 of the ascensional
    difference arc) onto the ecliptic by iterating the pole height of the house
    circle until the cusp longitude converges.

    This is the spherical-trigonometry core shared by both quadrant loops of the
    36-sector system (Ref: Gauquelin, "L'influence des astres", 1955; geometric
    derivation follows Meeus, "Astronomical Algorithms", Ch. 13).

    Args:
        sector_offset: Integer fraction index (1-8) of the quadrant arc.
        base_ramc: Base RAMC value for the quadrant (armc or armc+180).
        ramc_sign: +1.0 or -1.0 controlling the direction of the RAMC offset.
        ascensional_diff: Ascensional difference angle (degrees) for the latitude.
        tan_obliquity: tan(ε) — tangent of the ecliptic obliquity.
        sin_obliquity: sin(ε) — sine of the ecliptic obliquity.
        tan_lat: tan(φ) — tangent of the geographic latitude.
        eps: True obliquity of ecliptic (degrees).
        lat: Geographic latitude (degrees).
        niter_max: Maximum number of refinement iterations.
        convergence_threshold: Angular convergence tolerance (degrees).
        near_zero: Threshold for treating a tangent as effectively zero.

    Returns:
        Ecliptic longitude of the sector cusp (degrees, 0-360).
    """
    # Intermediate RAMC for this sector
    intermediate_ramc = (base_ramc + ramc_sign * (90.0 / 9.0) * sector_offset) % 360.0

    # Degenerate obliquity (eps -> 0): the ecliptic coincides with the equator,
    # so every sector cusp lies on the equator at the intermediate RAMC — the
    # same value the equatorial-declination branch below returns. Guard the
    # division by tan(eps) that would otherwise raise ZeroDivisionError.
    if abs(tan_obliquity) < near_zero:
        return intermediate_ramc

    # Initial pole height: latitude of the great circle pole for this sector fraction
    arc_fraction_deg = ascensional_diff * sector_offset / 9.0
    initial_pole_height = math.degrees(
        math.atan(math.sin(math.radians(arc_fraction_deg)) / tan_obliquity)
    )

    # Seed estimate: ecliptic point rising at the intermediate RAMC with the initial pole
    cusp = _calc_ascendant(intermediate_ramc, eps, lat, initial_pole_height)

    # Declination of the seed cusp point
    cusp_dec = math.asin(sin_obliquity * math.sin(math.radians(cusp)))
    tan_cusp_dec = math.tan(cusp_dec)

    # Degenerate case: cusp on the equator, no pole-height correction possible
    if abs(tan_cusp_dec) < near_zero:
        return intermediate_ramc

    # First pole-height refinement from seed
    pole_arg = max(-1.0, min(1.0, tan_lat * tan_cusp_dec))
    pole_height = math.degrees(
        math.atan(
            math.sin(
                math.radians(math.degrees(math.asin(pole_arg)) * sector_offset / 9.0)
            )
            / tan_cusp_dec
        )
    )
    cusp = _calc_ascendant(intermediate_ramc, eps, lat, pole_height)

    # Iterative refinement until convergence
    prev_cusp = 0.0
    for iteration in range(1, niter_max + 1):
        updated_dec = math.asin(sin_obliquity * math.sin(math.radians(cusp)))
        tan_cusp_dec_updated = math.tan(updated_dec)

        if abs(tan_cusp_dec_updated) < near_zero:
            cusp = intermediate_ramc
            break

        pole_arg = max(-1.0, min(1.0, tan_lat * tan_cusp_dec_updated))
        pole_height = math.degrees(
            math.atan(
                math.sin(
                    math.radians(
                        math.degrees(math.asin(pole_arg)) * sector_offset / 9.0
                    )
                )
                / tan_cusp_dec_updated
            )
        )
        cusp = _calc_ascendant(intermediate_ramc, eps, lat, pole_height)

        if iteration > 1:
            angular_diff = abs(cusp - prev_cusp)
            if angular_diff > 180.0:
                angular_diff = 360.0 - angular_diff
            if angular_diff < convergence_threshold:
                break
        prev_cusp = cusp

    return cusp


def _houses_gauquelin(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Gauquelin 36-Sector house system.

    Divides the celestial sphere into 36 sectors of 10° each using a
    Placidus-inspired construction: the diurnal and nocturnal semi-arcs are
    divided into ninths, and the resulting great circles are projected onto
    the ecliptic via iterative pole-height refinement.

    Reference: Michel Gauquelin, "L'influence des astres" (1955).
    Mathematical derivation: Meeus, "Astronomical Algorithms" 2nd ed., Ch. 13.

    Algorithm:
        1. Compute the ascensional difference:
               a = arcsin(tan(φ) · tan(ε))
           where φ = geographic latitude, ε = ecliptic obliquity.
        2. Divide each quadrant arc into ninths.
        3. For each of the 8 intermediate sectors per quadrant, call
           _gauquelin_cusp_for_sector() to find the ecliptic intersection via
           iterative pole-height refinement.
        4. Opposite sectors are set by adding 180°.
        5. Cardinals: sector 1 = Asc, 10 = MC, 19 = Desc, 28 = IC.

    Polar Limitation:
        Fails within polar circle (|φ| >= 90° - ε); falls back to equal
        10°-divisions from the Ascendant.

    Args:
        armc: Right Ascension of the Midheaven (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 37 sector longitudes (index 0 unused, sectors 1-36)
    """
    NEAR_ZERO = 1e-10
    CONVERGENCE_THRESHOLD = 1e-8
    niter_max = 100

    cusps = [0.0] * 37  # index 0 unused; sectors 1-36

    # Polar circle fallback: equal 10° divisions from the Ascendant
    if abs(lat) >= 90.0 - eps:
        for i in range(1, 37):
            cusps[i] = (asc - (i - 1) * 10.0) % 360.0
        return cusps

    # Precompute obliquity and latitude trig quantities
    sin_obliquity = math.sin(math.radians(eps))
    tan_obliquity = math.tan(math.radians(eps))
    tan_lat = math.tan(math.radians(lat))

    # Ascensional difference: a = arcsin(tan(φ) · tan(ε))
    asc_diff_arg = max(-1.0, min(1.0, tan_lat * tan_obliquity))
    ascensional_diff = math.degrees(math.asin(asc_diff_arg))

    # Shared keyword arguments for the cusp helper
    _cusp_kwargs: dict[str, Any] = dict(
        ascensional_diff=ascensional_diff,
        tan_obliquity=tan_obliquity,
        sin_obliquity=sin_obliquity,
        tan_lat=tan_lat,
        eps=eps,
        lat=lat,
        niter_max=niter_max,
        convergence_threshold=CONVERGENCE_THRESHOLD,
        near_zero=NEAR_ZERO,
    )

    # --- Fourth/Second quadrant: sectors 2-9 (RAMC offsets +10° … +80°) ---
    for ih in range(2, 10):
        sector_offset = 10 - ih  # counts down 8..1 as ih goes 2..9
        cusps[ih] = _gauquelin_cusp_for_sector(
            sector_offset=sector_offset,
            base_ramc=armc,
            ramc_sign=+1.0,
            **_cusp_kwargs,
        )
        cusps[ih + 18] = (cusps[ih] + 180.0) % 360.0  # opposite sector

    # --- First/Third quadrant: sectors 29-36 (RAMC offsets -10° … -80°) ---
    for ih in range(29, 37):
        sector_offset = ih - 28  # counts up 1..8 as ih goes 29..36
        cusps[ih] = _gauquelin_cusp_for_sector(
            sector_offset=sector_offset,
            base_ramc=(armc + 180.0) % 360.0,
            ramc_sign=-1.0,
            **_cusp_kwargs,
        )
        cusps[ih - 18] = (cusps[ih] + 180.0) % 360.0  # opposite sector

    # Cardinal sectors
    cusps[1] = asc
    cusps[10] = mc
    cusps[19] = (asc + 180.0) % 360.0
    cusps[28] = (mc + 180.0) % 360.0

    return cusps


def _ecliptic_to_ra_simple(lon: float, eps: float) -> float:
    """Convert ecliptic longitude to right ascension."""
    rad_lon = math.radians(lon)
    rad_eps = math.radians(eps)

    y = math.sin(rad_lon) * math.cos(rad_eps)
    x = math.cos(rad_lon)
    ra = math.degrees(math.atan2(y, x))

    return ra % 360.0


def _ra_to_ecliptic_simple(ra: float, eps: float) -> float:
    """Convert right ascension to ecliptic longitude."""
    rad_ra = math.radians(ra)
    rad_eps = math.radians(eps)

    y = math.sin(rad_ra)
    x = math.cos(rad_ra) * math.cos(rad_eps)
    lon = math.degrees(math.atan2(y, x))

    return lon % 360.0


def _rotate_spherical_x_axis(x: List[float], eps: float) -> List[float]:
    """
    Standard rotation of spherical coordinates around the x-axis.

    Rotates spherical coordinates [lon, lat, r] by angle eps.
    This is the standard x-axis rotation matrix applied to spherical
    coordinates, used for ecliptic/equatorial conversion and similar
    coordinate frame rotations.

    Reference: Montenbruck & Pfleger, 'Astronomy on the Personal Computer',
    par. 1.3; Meeus, 'Astronomical Algorithms' 2nd ed., ch. 13.

    Args:
        x: [longitude, latitude, radius] in degrees
        eps: Rotation angle in degrees

    Returns:
        Transformed [longitude, latitude, radius]
    """
    eps_rad = math.radians(eps)
    lon_rad = math.radians(x[0])
    lat_rad = math.radians(x[1])
    r = x[2]

    # Convert to Cartesian
    cos_lat = math.cos(lat_rad)
    x_cart = r * cos_lat * math.cos(lon_rad)
    y_cart = r * cos_lat * math.sin(lon_rad)
    z_cart = r * math.sin(lat_rad)

    # Rotate around x-axis by eps
    cos_eps = math.cos(eps_rad)
    sin_eps = math.sin(eps_rad)

    x_new = x_cart
    y_new = y_cart * cos_eps + z_cart * sin_eps
    z_new = -y_cart * sin_eps + z_cart * cos_eps

    # Convert back to spherical
    r_new = math.sqrt(x_new**2 + y_new**2 + z_new**2)

    if r_new > 0:
        lon_new = math.degrees(math.atan2(y_new, x_new))
        # r_new >= |z_new| mathematically, but rounding can push the ratio a
        # hair past +-1; clamp so asin never raises a domain error.
        lat_new = math.degrees(math.asin(max(-1.0, min(1.0, z_new / r_new))))
    else:
        lon_new = 0.0
        lat_new = 0.0

    return [(lon_new % 360.0), lat_new, r_new]


def _houses_equal_mc(asc: float, mc: float) -> List[float]:
    """
    Equal houses from MC (code 'D').

    All houses are exactly 30° each, with the MC as the 10th house cusp.
    Unlike Equal from Ascendant, the Ascendant is NOT a house cusp in this system.

    Algorithm:
        H10 = MC
        H11 = MC + 30°
        H12 = MC + 60°
        H1 = MC + 90° (= MC - 270°)
        ... and so on

    Properties:
        - Works at all latitudes
        - MC is exactly on the 10th house cusp
        - Ascendant does NOT coincide with the 1st house cusp (stored in ascmc)

    Args:
        asc: Ascendant longitude in degrees (unused, kept for API compatibility)
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = [0.0] * 13
    # H10 = MC, then each house is 30° apart
    # H1 = MC - 270° = MC + 90°
    for i in range(1, 13):
        # House 10 is at MC, so offset = (i - 10) * 30
        cusps[i] = (mc + (i - 10) * 30.0) % 360.0
    return cusps


def _sunshine_arc_to_ecliptic(
    arc_offset_ra: float,
    is_diurnal: bool,
    armc: float,
    lat: float,
    sun_dec: float,
    sin_lat: float,
    cos_lat: float,
    cos_sun_dec: float,
    tan_sun_dec: float,
    sin_obliquity: float,
    cos_obliquity: float,
) -> float:
    """
    Map one trisected semi-arc offset to an ecliptic longitude.

    The Sunshine system divides the Sun's diurnal and nocturnal semi-arcs into
    thirds.  For each trisection point, this function uses spherical triangle
    geometry to find the corresponding ecliptic cusp (Makransky 1988,
    "Primary Directions").

    Steps:
        1. Convert the RA-arc offset to a great-circle arc length (xhs).
        2. Solve the spherical triangle {pole, arc midpoint, horizon point}
           via the spherical law of cosines to find side c.
        3. Apply the spherical law of sines to get the zenith distance (zd).
        4. Derive the equatorial RA intersection and the pole height of the
           house circle.
        5. Project to the ecliptic via the rising-point formula.

    Args:
        arc_offset_ra: RA offset along the semi-arc (degrees).
        is_diurnal: True for diurnal semi-arc cusps (houses 8-12 side),
                    False for nocturnal semi-arc cusps (houses 2-6 side).
        armc: Right Ascension of the Midheaven (degrees).
        lat: Geographic latitude (degrees).
        sun_dec: Sun's declination (degrees).
        sin_lat: sin(φ) — pre-computed for speed.
        cos_lat: cos(φ) — pre-computed for speed.
        cos_sun_dec: cos(δ☉) — pre-computed for speed.
        tan_sun_dec: tan(δ☉) — pre-computed for speed.
        sin_obliquity: sin(ε) — pre-computed for speed.
        cos_obliquity: cos(ε) — pre-computed for speed.

    Returns:
        Ecliptic longitude of the cusp (degrees, 0-360).
    """
    # Step 1: great-circle arc subtended by the RA offset on the Sun's parallel
    # xhs = 2 · arcsin(cos(δ☉) · sin(arc_offset / 2))
    great_circle_arc = 2.0 * math.degrees(
        math.asin(cos_sun_dec * math.sin(math.radians(arc_offset_ra / 2.0)))
    )

    # Step 2a: interior angle α at the pole of the spherical triangle
    # cos(α) = tan(δ☉) · tan(great_circle_arc / 2)
    cos_alpha = max(
        -1.0, min(1.0, tan_sun_dec * math.tan(math.radians(great_circle_arc / 2.0)))
    )
    alpha = math.degrees(math.acos(cos_alpha))

    # Step 2b: triangle side and opening angle differ for diurnal vs nocturnal arc
    if is_diurnal:
        # Diurnal arc: triangle opening angle is supplement of α
        triangle_angle = 180.0 - alpha
        triangle_side = 90.0 - lat + sun_dec
    else:
        # Nocturnal arc: triangle opening angle equals α
        triangle_angle = alpha
        triangle_side = 90.0 - lat - sun_dec

    # Step 2c: third side c via spherical law of cosines
    cos_c = math.cos(math.radians(great_circle_arc)) * math.cos(
        math.radians(triangle_side)
    ) + math.sin(math.radians(great_circle_arc)) * math.sin(
        math.radians(triangle_side)
    ) * math.cos(math.radians(triangle_angle))
    cos_c = max(-1.0, min(1.0, cos_c))
    side_c = math.degrees(math.acos(cos_c))

    # Step 3: zenith distance via spherical law of sines
    if side_c < 1e-6:
        zenith_dist = 0.0
    else:
        sin_zd = max(
            -1.0,
            min(
                1.0,
                math.sin(math.radians(great_circle_arc))
                * math.sin(math.radians(triangle_angle))
                / math.sin(math.radians(side_c)),
            ),
        )
        zenith_dist = math.degrees(math.asin(sin_zd))

    # Step 4: equatorial RA of the house-circle / equator intersection
    equator_ra_offset = math.degrees(
        math.atan(cos_lat * math.tan(math.radians(zenith_dist)))
    )

    # Pole height of the house circle
    pole_height = math.degrees(math.asin(math.sin(math.radians(zenith_dist)) * sin_lat))

    # Step 5: sign/direction convention and final RA
    if is_diurnal:
        ra_intersection = (armc + equator_ra_offset) % 360.0
    else:
        pole_height = -pole_height
        ra_intersection = (equator_ra_offset + armc + 180.0) % 360.0

    # Project equatorial intersection onto the ecliptic
    return _ra_to_ecliptic_longitude(
        ra_intersection, pole_height, sin_obliquity, cos_obliquity
    )


def _houses_sunshine(
    armc: float, lat: float, eps: float, asc: float, mc: float, sun_dec: float
) -> List[float]:
    """
    Sunshine/Treindl numerical solution (code ``I``) of Makransky's system.

    Bob Makransky's public Sunshine definition divides the diurnal and
    nocturnal arcs of the Sun, then carries the division anchors to the
    ecliptic through the horizon north/south house circles. Code ``I`` denotes
    this project's documented Treindl-style numerical realization; code ``i``
    is the separately implemented Makransky upper/lower-meridian construction
    in ``house_constructions.py``. The exact source/status boundary is recorded
    in ``docs/reference/house-systems.md``.

    Note: Cusps 11, 12, 2, 3 are NOT in exact opposition to cusps 5, 6, 8, 9.

    Algorithm:
        1. Compute the Sun's ascensional difference (AD) from its declination
           and the geographic latitude.
        2. Derive the diurnal semi-arc (DSA = 90° + AD) and the nocturnal
           semi-arc (NSA = 90° − AD); trisect each to get RA offsets.
        3. For nocturnal cusps (2, 3, 5, 6) call _sunshine_arc_to_ecliptic()
           with is_diurnal=False.
        4. For diurnal cusps (8, 9, 11, 12) call _sunshine_arc_to_ecliptic()
           with is_diurnal=True.
        5. If the MC is below the horizon, reflect all intermediate cusps by
           180°.

    Args:
        armc: Right Ascension of the Midheaven (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees
        sun_dec: Sun's declination in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    NEAR_ZERO = 1e-10

    # Detect MC-below-horizon condition.
    # mcdec = declination of the MC, derived from ARMC and obliquity.
    # When |lat - mcdec| > 90°, the MC is below the horizon.
    mcdec = math.degrees(
        math.atan(math.sin(math.radians(armc)) * math.tan(math.radians(eps)))
    )
    mc_under_horizon = abs(lat - mcdec) > 90

    cusps = [0.0] * 13

    # When MC is below the horizon, flip MC and IC by 180° to keep MC
    # above the horizon (analogous to Regiomontanus polar handling).
    # The ASC is already hemisphere-corrected by the dispatch.
    if mc_under_horizon:
        mc = (mc + 180.0) % 360.0

    cusps[1] = asc
    cusps[10] = mc
    cusps[4] = (mc + 180.0) % 360.0
    cusps[7] = (asc + 180.0) % 360.0

    # Ascensional difference: sin(AD) = tan(δ☉) · tan(φ)
    ad_arg = math.tan(math.radians(sun_dec)) * math.tan(math.radians(lat))
    if ad_arg >= 1.0:
        ascensional_diff = 90.0 - NEAR_ZERO
    elif ad_arg <= -1.0:
        ascensional_diff = -90.0 + NEAR_ZERO
    else:
        ascensional_diff = math.degrees(math.asin(ad_arg))

    nocturnal_semi_arc = 90.0 - ascensional_diff  # NSA
    diurnal_semi_arc = 90.0 + ascensional_diff  # DSA

    # Pre-compute repeated trig quantities
    sin_lat = math.sin(math.radians(lat))
    cos_lat = math.cos(math.radians(lat))
    cos_sun_dec = math.cos(math.radians(sun_dec))
    tan_sun_dec = math.tan(math.radians(sun_dec))
    sin_obliquity = math.sin(math.radians(eps))
    cos_obliquity = math.cos(math.radians(eps))

    # RA offsets along each semi-arc (positive = toward MC, negative = toward IC)
    # NSA offsets (nocturnal cusps 2, 3, 5, 6)
    nsa_third = nocturnal_semi_arc / 3.0
    nsa_offsets = {
        2: -2.0 * nsa_third,
        3: -1.0 * nsa_third,
        5: 1.0 * nsa_third,
        6: 2.0 * nsa_third,
    }
    # DSA offsets (diurnal cusps 8, 9, 11, 12)
    dsa_third = diurnal_semi_arc / 3.0
    dsa_offsets = {
        8: -2.0 * dsa_third,
        9: -1.0 * dsa_third,
        11: 1.0 * dsa_third,
        12: 2.0 * dsa_third,
    }

    # Shared kwargs for the geometric helper
    _common = dict(
        armc=armc,
        lat=lat,
        sun_dec=sun_dec,
        sin_lat=sin_lat,
        cos_lat=cos_lat,
        cos_sun_dec=cos_sun_dec,
        tan_sun_dec=tan_sun_dec,
        sin_obliquity=sin_obliquity,
        cos_obliquity=cos_obliquity,
    )

    # --- Nocturnal semi-arc cusps (below the horizon arc) ---
    for house_index, ra_offset in nsa_offsets.items():
        cusps[house_index] = _sunshine_arc_to_ecliptic(
            arc_offset_ra=ra_offset,
            is_diurnal=False,
            **_common,
        )

    # --- Diurnal semi-arc cusps (above the horizon arc) ---
    for house_index, ra_offset in dsa_offsets.items():
        cusps[house_index] = _sunshine_arc_to_ecliptic(
            arc_offset_ra=ra_offset,
            is_diurnal=True,
            **_common,
        )

    # When MC is below the horizon, shift all intermediate cusps by 180°.
    # This keeps the MC above the horizon, consistent with the polar-circle
    # handling used for Regiomontanus and Campanus systems.
    if mc_under_horizon:
        for i in (2, 3, 5, 6, 8, 9, 11, 12):
            cusps[i] = (cusps[i] + 180.0) % 360.0

    return cusps


def _houses_horizontal(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Horizontal (Azimuthal) house system (code 'H').

    The horizon is divided into 12 equal arcs of 30 degrees each, starting
    from the East Point. For each division point, a vertical circle (great
    circle through zenith and nadir) passes through that point and intersects
    the ecliptic. Those intersection points define the house cusps.

    Algorithm derived from the geometric definition described in:
    Ralph William Holden, 'The Elements of House Division' (1977), ch. 10.

    Geometric derivation (standard spherical trigonometry):
        Let phi' = co-latitude (90 deg - |phi|, signed).
        For each horizon division at angle theta from the meridian
        (theta = k * 30 deg, k = 1..5), the vertical circle's parameters are:

            pole_height = arcsin(sin(theta) * sin(phi'))
            ra_offset   = atan2(cos(theta), sin(theta) * cos(phi'))

        These follow from solving the spherical triangle (zenith, celestial
        pole, horizon division point) via the spherical law of sines/cosines.

        The ecliptic intersection is then computed via the rising-point formula
        (Smart, 'Textbook on Spherical Astronomy', ch. 3).

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes (index 0 unused, 1-12 are cusps)
    """
    cusps = [0.0] * 13
    _NEAR_ZERO = 1e-10

    # Co-latitude: complement of geographic latitude
    if lat > 0:
        co_lat = 90.0 - lat
    else:
        co_lat = -90.0 - lat

    # Clamp to avoid singularity at the equator (|co_lat| = 90 deg)
    if abs(abs(co_lat) - 90.0) < _NEAR_ZERO:
        co_lat = (90.0 - _NEAR_ZERO) if co_lat > 0 else (-90.0 + _NEAR_ZERO)

    # ARMC rotated by 180 deg (base reference for house calculations)
    armc_base = (armc + 180.0) % 360.0

    sin_colat = math.sin(math.radians(co_lat))
    cos_colat = math.cos(math.radians(co_lat))

    # Compute 5 primary cusps by looping over horizon division angles.
    # The horizon is divided at theta = k * 30 deg (k = 1..5) measured from
    # the meridian toward the East Point along the horizon.
    #
    # Cusp mapping:  theta=30 -> cusp 3,  theta=60 -> cusp 2,
    #                theta=90 -> cusp 1 (East Point),
    #                theta=120 -> cusp 12, theta=150 -> cusp 11
    _CUSP_MAP = [(3, 30.0), (2, 60.0), (1, 90.0), (12, 120.0), (11, 150.0)]

    for cusp_idx, theta_deg in _CUSP_MAP:
        theta_rad = math.radians(theta_deg)
        sin_theta = math.sin(theta_rad)
        cos_theta = math.cos(theta_rad)

        # Pole height of the vertical circle at this division angle
        pole_height = math.degrees(math.asin(sin_theta * sin_colat))

        # RA offset from the base RA position (ARMC + 270 deg)
        ra_offset = math.degrees(math.atan2(cos_theta, sin_theta * cos_colat))

        cusps[cusp_idx] = _calc_ascendant(
            armc_base + 90.0 + ra_offset, eps, lat, pole_height
        )

    # Polar circle handling: check Asc-MC orientation
    if abs(co_lat) >= 90.0 - eps:
        asc_mc_diff = (asc - mc + 540.0) % 360.0 - 180.0
        if asc_mc_diff < 0:
            asc = (asc + 180.0) % 360.0
            mc = (mc + 180.0) % 360.0
            for i in range(1, 13):
                if 4 <= i < 10:
                    continue
                cusps[i] = (cusps[i] + 180.0) % 360.0

    # Flip primary cusps by 180 deg (convention: cusps represent the opposite
    # hemisphere from the direct calculation)
    for i in [1, 2, 3, 11, 12]:
        cusps[i] = (cusps[i] + 180.0) % 360.0

    # Check Asc/DC orientation
    asc_mc_diff = (asc - mc + 540.0) % 360.0 - 180.0
    if asc_mc_diff < 0:
        asc = (asc + 180.0) % 360.0

    # Set MC/IC and derive opposite cusps (4-9 from 10-3)
    cusps[10] = mc
    cusps[4] = (mc + 180.0) % 360.0
    cusps[7] = (cusps[1] + 180.0) % 360.0
    cusps[8] = (cusps[2] + 180.0) % 360.0
    cusps[9] = (cusps[3] + 180.0) % 360.0
    cusps[5] = (cusps[11] + 180.0) % 360.0
    cusps[6] = (cusps[12] + 180.0) % 360.0

    return cusps


def _houses_natural_gradient(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Natural Gradient house system ('N').

    'N' maps to "Equal houses with 0° Aries as cusp 1".
    This is effectively a Whole Sign system starting from 0° Aries.

    Mathematical Formula:
        λᵢ = (i-1) × 30°
        where i = 1, 2, ..., 12

    Result:
        House 1 = 0° Aries (0°)
        House 2 = 0° Taurus (30°)
        House 3 = 0° Gemini (60°)
        ... etc.

    Properties:
        - Independent of time and location
        - Works at all latitudes
        - Computationally exact

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees (unused)
        lat: Geographic latitude in degrees (unused)
        eps: True obliquity of ecliptic in degrees (unused)
        asc: Ascendant longitude in degrees (unused)
        mc: Midheaven longitude in degrees (unused)

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = [0.0] * 13
    for i in range(1, 13):
        cusps[i] = ((i - 1) * 30.0) % 360.0
    return cusps


def _houses_apc(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    APC (Ascendant-Parallel Circle) house system (code 'Y').

    The parallel of declination passing through the Ascendant is divided
    into 6 equal arcs above and 6 below the horizon. Position circles
    through the North and South points of the horizon pass through each
    division point; their intersections with the ecliptic define the cusps.

    Algorithm derived from the geometric definition by L. Knegt
    (WvA/Ram school, Netherlands).

    Reference: Knegt, 'Handleiding voor het berekenen van huizentabellen'
    (Manual for computing house tables).

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes (index 0 unused, 1-12 are cusps)
    """
    cusps = [0.0] * 13

    # Per-cusp computation works in radians
    lat_rad = math.radians(lat)
    eps_rad = math.radians(eps)
    armc_rad = math.radians(armc)

    # Calculate all 12 house cusps
    for i in range(1, 13):
        cusps[i] = _apc_cusp(i, lat_rad, eps_rad, armc_rad)

    # The per-cusp formula is ill-conditioned for the MC near |lat| = 90;
    # the directly computed MC is exact, so always substitute it.
    cusps[10] = mc
    cusps[4] = (mc + 180.0) % 360.0

    # Within polar circle, handle horizon crossing: when the Ascendant trails
    # the MC the whole wheel is turned half a round (clockwise direction).
    if abs(lat) >= 90.0 - eps:
        if difdeg2n(asc, mc) < 0:
            for i in range(1, 13):
                cusps[i] = (cusps[i] + 180.0) % 360.0

    return cusps


# ---------------------------------------------------------------------------
# House-system registry
# ---------------------------------------------------------------------------
# One row per selector character, consulted by _houses_from_armc through
# _house_system. The adapters below bind each construction's own argument
# list to the _HouseFrame the dispatcher hands out.

_FrameConstruction = Callable[[float, float, float, float, float], List[float]]


def _builds_on_frame(
    construction: _FrameConstruction,
) -> Callable[[_HouseFrame], List[float]]:
    """Adapt a construction taking ``(armc, lat, eps, asc, mc)``."""
    return lambda frame: construction(
        frame.armc, frame.lat, frame.eps, frame.asc, frame.mc
    )


def _builds_on_ascendant(
    construction: Callable[[float], List[float]],
) -> Callable[[_HouseFrame], List[float]]:
    """Adapt a construction taking the Ascendant only."""
    return lambda frame: construction(frame.asc)


def _builds_on_angles(
    construction: Callable[[float, float], List[float]],
) -> Callable[[_HouseFrame], List[float]]:
    """Adapt a construction taking ``(asc, mc)``."""
    return lambda frame: construction(frame.asc, frame.mc)


def _builds_on_frame_and_sun(
    construction: Callable[[float, float, float, float, float, float], List[float]],
) -> Callable[[_HouseFrame], List[float]]:
    """Adapt a construction taking ``(armc, lat, eps, asc, mc, sun_dec)``."""
    return lambda frame: construction(
        frame.armc, frame.lat, frame.eps, frame.asc, frame.mc, frame.sun_dec
    )


#: Placidus, also the system every selector the registry does not list is
#: served with; like Placidus it raises inside the polar circle.
_DEFAULT_HOUSE_SYSTEM = HouseSystem(
    "P", _builds_on_frame(_houses_placidus), raises_at_pole=True
)

_HOUSE_SYSTEMS: dict[str, HouseSystem] = {
    system.code: system
    for system in (
        _DEFAULT_HOUSE_SYSTEM,
        HouseSystem("K", _builds_on_frame(_houses_koch), raises_at_pole=True),
        HouseSystem(
            "R",
            _builds_on_frame(_houses_regiomontanus),
            mc_branch=_McBranch.ROTATE_FRAME,
        ),
        HouseSystem(
            "C", _builds_on_frame(_houses_campanus), mc_branch=_McBranch.ROTATE_FRAME
        ),
        # 'A' is served by the same Ascendant-anchored equal wheel as 'E'.
        HouseSystem("E", _builds_on_ascendant(_houses_equal)),
        HouseSystem("A", _builds_on_ascendant(_houses_equal)),
        HouseSystem("W", _builds_on_ascendant(_houses_whole_sign)),
        HouseSystem("O", _builds_on_angles(_houses_porphyry)),
        HouseSystem("B", _builds_on_frame(_houses_alcabitius)),
        HouseSystem(
            "T",
            _builds_on_frame(_houses_polich_page),
            mc_branch=_McBranch.ROTATE_FRAME,
        ),
        HouseSystem("M", _builds_on_frame(_houses_morinus)),
        HouseSystem("X", _builds_on_frame(_houses_meridian)),
        HouseSystem("V", _builds_on_ascendant(_houses_vehlow)),
        HouseSystem("H", _builds_on_frame(_houses_horizontal)),
        HouseSystem(
            "Y", _builds_on_frame(_houses_apc), mc_branch=_McBranch.ROTATE_REPORTED
        ),
        HouseSystem("F", _builds_on_frame(_houses_carter)),
        HouseSystem("U", _builds_on_frame(_houses_krusinski)),
        HouseSystem("N", _builds_on_frame(_houses_natural_gradient)),
        HouseSystem(
            "G", _builds_on_frame(_houses_gauquelin), raises_at_pole=True, sectors=36
        ),
        # Lowercase 'g' computes the same 36 Gauquelin sectors, but the
        # RETURN SHAPE is keyed on the uppercase byte only (compatibility
        # contract): it yields the first 12 sectors in the ordinary 12-cusp
        # tuple. house_pos and house_name fold it to 'G' instead.
        HouseSystem("g", _builds_on_frame(_houses_gauquelin), raises_at_pole=True),
        HouseSystem("S", _builds_on_angles(_houses_sripati)),
        HouseSystem("L", _builds_on_angles(_houses_pullen_sd)),
        HouseSystem("Q", _builds_on_angles(_houses_pullen_sr)),
        HouseSystem("D", _builds_on_angles(_houses_equal_mc)),
        HouseSystem(
            "I",
            _builds_on_frame_and_sun(_houses_sunshine),
            mc_branch=_McBranch.FROM_CUSPS,
            reads_sun_declination=True,
        ),
        HouseSystem(
            "i",
            _builds_on_frame_and_sun(_houses_sunshine_makransky),
            mc_branch=_McBranch.FROM_CUSPS,
            reads_sun_declination=True,
        ),
        HouseSystem(
            "J", _builds_on_frame(_houses_savard_a), mc_branch=_McBranch.ROTATE_FRAME
        ),
    )
}


# Shared placement constants. Degenerate spherical arguments use a small
# binary64 guard. Cusp placement is unbiased; exact cusps belong to the
# interval selected by the geometric comparisons below.
_NEAR_ZERO = 1e-10
CUSP_BOUNDARY_OFFSET = 0.0


# --- Degree-argument trigonometry ------------------------------------------
# Spherical-astronomy formulas for house placement are stated in degrees in
# every published source (semi-arcs, meridian distances, poles). These
# one-line wrappers keep the formulas below readable in the same units the
# literature uses; the inverse functions clamp their argument to the
# principal domain, the standard guard against |x| creeping past 1 by a few
# ULPs in intermediate float arithmetic.


def _sin_deg(x: float) -> float:
    """Return the sine of an angle supplied in degrees."""
    return math.sin(math.radians(x))


def _cos_deg(x: float) -> float:
    """Return the cosine of an angle supplied in degrees."""
    return math.cos(math.radians(x))


def _tan_deg(x: float) -> float:
    """Return the tangent of an angle supplied in degrees."""
    return math.tan(math.radians(x))


def _atan_deg(x: float) -> float:
    """Return the principal arctangent in degrees."""
    return math.degrees(math.atan(x))


def _asin_deg(x: float) -> float:
    """Return arcsine in degrees after a roundoff-only unit-domain clamp."""
    return math.degrees(math.asin(max(-1.0, min(1.0, x))))


def _acos_deg(x: float) -> float:
    """Return arccosine in degrees after a roundoff-only unit-domain clamp."""
    return math.degrees(math.acos(max(-1.0, min(1.0, x))))


def _rotate_frame(lon: float, lat: float, tilt_deg: float) -> Tuple[float, float]:
    """Tilt a spherical position into a frame rotated about the x-axis.

    Thin wrapper over :func:`libephemeris.utils.cotrans` (this project's
    public coordinate-rotation primitive) for callers that carry no radial
    component. A positive ``tilt_deg`` rotates by the textbook x-axis
    rotation used throughout spherical astronomy for plane changes —
    ecliptic/equator, equator/horizon, equator/prime-vertical (Meeus,
    'Astronomical Algorithms', 2nd ed., Ch. 13; Montenbruck & Pfleger,
    'Astronomy on the Personal Computer', par. 1.3).

    Returns:
        ``(lon_out, lat_out)`` in degrees, longitude normalized to
        ``[0, 360)``.
    """
    lon_out, lat_out, _ = cotrans((lon, lat, 1.0), tilt_deg)
    return lon_out, lat_out


def _ascendant_on_eastern_horizon(
    asc: float, armc: float, eps: float, geolat: float
) -> float:
    """Return the ascendant longitude anchored to the rising intersection.

    The horizon and the ecliptic intersect in two opposite points; the
    ascendant is by definition the *eastern* (rising) one. The closed-form
    ecliptic-longitude expression for the ascendant silently picks the
    intersection nearer the meridian, which inside the polar circles can be
    the setting point — a six-house (180 deg) error for every equal-style
    system anchored to the ascendant.

    Geometric test: the MC degree of the ecliptic has declination
    ``delta_mc = atan(sin ARMC * tan eps)`` (the standard declination of an
    ecliptic point expressed through its right ascension). Its meridian
    altitude is ``90 - phi + delta_mc`` in the northern hemisphere and
    ``-(-90 - phi + delta_mc)`` in the southern one; when that altitude
    goes negative the MC has sunk below the horizon and the two
    horizon-ecliptic intersections have swapped roles, so the raw
    ascendant must be flipped to the opposite point.

    Args:
        asc: Ascendant longitude from the closed-form expression (degrees).
        armc: Right ascension of the MC (degrees).
        eps: True obliquity of the ecliptic (degrees).
        geolat: Geographic latitude (degrees).

    Returns:
        The ascendant longitude on the eastern horizon (degrees).
    """
    delta_mc = _atan_deg(_sin_deg(armc) * _tan_deg(eps))
    if geolat >= 0 and 90.0 - geolat + delta_mc < 0.0:
        asc = (asc + 180.0) % 360.0
    if geolat < 0 and -90.0 - geolat + delta_mc > 0.0:
        asc = (asc + 180.0) % 360.0
    return asc


def _above_horizon_test(dec: float, geolat: float, merid_dist: float) -> float:
    """Signed above-horizon indicator for a body at hour-angle distance.

    From the standard altitude formula ``sin h = sin(phi) sin(delta) +
    cos(phi) cos(delta) cos(H)``: dividing through by ``cos(phi) cos(delta)``
    (positive for |phi|, |delta| < 90) leaves the sign of the altitude in
    ``tan(phi) tan(delta) + cos(H)``. Non-negative means the body sits on or
    above the horizon.
    """
    return _tan_deg(dec) * _tan_deg(geolat) + _cos_deg(merid_dist)


def _asc_diff_saturated(dec: float, geolat: float) -> float:
    """Ascensional difference, saturated at the circumpolar limit.

    ``AD = asin(tan(phi) tan(delta))`` (the classic rising-time correction,
    Meeus Ch. 15). Where the argument exceeds unit magnitude the body never
    crosses the horizon at this latitude; the difference saturates at
    +/-90 deg, which turns the corresponding semi-arc into 180 or 0 deg.
    """
    s = _tan_deg(dec) * _tan_deg(geolat)
    if s >= 1.0:
        return 90.0
    if s <= -1.0:
        return -90.0
    return _asin_deg(s)


#: Latitude at which the Topocentric family is read for a chart placed on a
#: geographic pole. The construction has no answer exactly there: ``tan(lat)``
#: is unbounded, every position circle whose pole share is non-zero collapses
#: onto the celestial equator, and one circle no longer divides the sky. A
#: thousandth of a degree inside the pole the circles are still distinct, the
#: latitude axis stays continuous, and the answers pile up on the meridian --
#: the limit that axis approaches.
_TOPOCENTRIC_POLE_LIMIT = 89.999


def _topocentric_pole_share(house_angle: float) -> float:
    """Return Polich and Page's share of ``tan(latitude)`` at a house angle.

    The pole height ``P`` of the position circle at house angle ``H`` is
    fixed by ``tan(P) = tan(lat) * share(H)``, and the share is the
    triangular wave through the nodes 0 on the meridian, +1 on the eastern
    horizon, 0 on the lower meridian and -1 on the western horizon. The
    linear interpolation puts the circles one and two houses from the
    meridian at one third and two thirds of the tangent -- the trisected
    tangent the system is named for (Polich & Page, "The Topocentric System
    of Houses", 1961/1964; Holden, "The Elements of House Division", 1977).
    Regiomontanus runs ``sin(H)`` through the same three nodes, which is why
    the two families agree at the four angles and nowhere else.

    Args:
        house_angle: Right ascension of the circle's equator crossing,
            measured eastwards from the ARMC, in degrees.

    Returns:
        The share of the tangent, in ``[-1, 1]``.
    """
    angle = house_angle % 360.0
    share_sign = 1.0
    if angle > 180.0:
        # Odd about the meridian: over the western half of the wheel the
        # fictitious observer stands in the opposite hemisphere.
        angle -= 180.0
        share_sign = -1.0
    if angle > 90.0:
        # Even about the horizon: the share climbs to the whole tangent at
        # the horizon and falls back to zero at the next meridian.
        angle = 180.0 - angle
    return share_sign * angle / 90.0


def _hpos_topocentric(
    armc: float, geolat: float, ra: float, dec: float, md_upper: float
) -> float:
    """Topocentric (Polich-Page) house position.

    The system divides the sky with a family of position circles: the circle
    at house angle ``H`` is the horizon of a fictitious observer standing at
    the pole height ``P`` of :func:`_topocentric_pole_share`, whose meridian
    is turned to ``ARMC + H - 90``. A body's position is the house angle of
    the single circle of that family that passes through it, so it is the
    root of the perpendicularity condition between the body and the pole of
    the circle::

        sin(dec) sin(P) - cos(dec) cos(P) sin(A - H) = 0

    with ``A`` the body's meridian distance measured eastwards. Equivalently
    the body's oblique ascension under the pole ``P(H)`` equals ``ARMC + H``.
    The equation is transcendental -- a straight line in ``H`` set equal to a
    sinusoid, so unlike Regiomontanus it collapses to no closed form -- but
    the two circles that bound each quadrant belong to the family themselves,
    which brackets the root in the quadrant the body occupies (Polich & Page,
    "The Topocentric System of Houses", Spica 3(3), 1964, and their 1961
    monograph; Holden, "The Elements of House Division", 1977, for the
    position circles and the oblique ascension; Smart, "Textbook on Spherical
    Astronomy", 6th ed., ch. II and VI).

    Args:
        armc: Right ascension of the medium coeli in degrees. Only its
            distance to the body enters the construction, and the caller has
            already reduced that into ``md_upper``.
        geolat: Geographic latitude in degrees, positive north. Its sign is
            carried through the tangent, so the two hemispheres are not
            mirror images of one another.
        ra: Right ascension of the body in degrees, likewise already reduced
            into ``md_upper``.
        dec: Declination of the body in degrees.
        md_upper: Meridian distance of the body -- its right ascension minus
            the ARMC -- in ``[-180, 180)`` degrees.

    Returns:
        The house position: the house number plus the fraction of that house
        the body has already crossed, in ``[1, 13)``.
    """
    # The body's meridian distance on the same wheel as the house angle:
    # 0 on the upper meridian, 90 on the eastern horizon, 180 on the lower.
    east_dist = md_upper % 360.0

    lat = max(-_TOPOCENTRIC_POLE_LIMIT, min(_TOPOCENTRIC_POLE_LIMIT, geolat))
    tan_lat = _tan_deg(lat)
    sin_dec = _sin_deg(dec)
    cos_dec = _cos_deg(dec)

    def circle_residual(house_angle: float) -> float:
        """Perpendicularity of the body to the pole of one position circle.

        Zero exactly when the body lies on the circle, and its sign tells
        which side of the circle the body is on. Taking ``sin(P) = share / hypot(1,
        share)`` and ``cos(P) = 1 / hypot(1, share)`` keeps the residual
        bounded and exact for every pole height, where the classical
        oblique-ascension form saturates its arcsine as soon as the body
        stops rising for the fictitious observer.
        """
        share = tan_lat * _topocentric_pole_share(house_angle)
        return (
            sin_dec * share - cos_dec * _sin_deg(east_dist - house_angle)
        ) / math.hypot(1.0, share)

    # The quadrant bounds are themselves members of the family: at H = 0 and
    # 180 the fictitious observer stands on the equator, so his horizon is
    # the meridian, and at H = 90 and 270 he is the real observer, so it is
    # the horizon. No body crosses one of those two circles without its house
    # angle crossing the matching bound, so the answer lies in the quadrant
    # the body itself occupies, and exactly one root lies there. The two
    # residuals that decide the quadrant are the classical tests in disguise:
    # at H = 0 the residual is -cos(dec) sin(A), negative east of the
    # meridian, and at H = 90 it is the sine of the body's altitude,
    # non-negative above the horizon.
    east_of_meridian = circle_residual(0.0) <= 0.0
    above_horizon = circle_residual(90.0) >= 0.0
    if east_of_meridian:
        lo = 0.0 if above_horizon else 90.0
    else:
        lo = 270.0 if above_horizon else 180.0
    hi = lo + 90.0

    f_lo = circle_residual(lo)
    f_hi = circle_residual(hi)
    if f_lo == 0.0:
        # A body sitting on one of the two bounding circles -- on the
        # meridian or on the horizon -- is answered by that angle itself.
        root = lo
    elif f_hi == 0.0:
        root = hi
    else:
        # The residual runs from negative at the lower bound to positive at
        # the upper one. Halve the bracket until its ends are neighbouring
        # binary64 angles and keep the end with the smaller residual: the
        # answer is then the root of the equation to the precision of the
        # arithmetic, not an iterate of the search.
        while True:
            mid = 0.5 * (lo + hi)
            if mid <= lo or mid >= hi:
                break
            f_mid = circle_residual(mid)
            if f_mid <= 0.0:
                lo, f_lo = mid, f_mid
            else:
                hi, f_hi = mid, f_mid
        root = lo if abs(f_lo) <= abs(f_hi) else hi

    # Read the house angle off the wheel: H = 0 is the tenth cusp, H = 90 the
    # first, and one house is 30 degrees of house angle.
    pos_deg = (root - 90.0) % 360.0
    if pos_deg >= 360.0:
        # A house angle a hair under 90 can round its way up to a whole turn,
        # which would encode the non-existent house 13.0 rather than the
        # first cusp; the interval is half-open.
        pos_deg -= 360.0
    return pos_deg / 30.0 + 1.0


def _hpos_horizon(md_upper: float, dec: float, geolat: float) -> float:
    """Horizon / azimuthal house position.

    The azimuthal system divides the horizon itself into twelve 30-degree
    houses, so a body is placed by its true azimuth, not by any ecliptic
    coordinate. Tilting the equatorial position (hour angle measured from
    the east point, declination) onto the horizon plane by the observer's
    co-latitude yields that azimuthal longitude directly — the same plane
    change as every other frame rotation in this module.
    """
    az_lon, _ = _rotate_frame((md_upper - 90.0) % 360.0, dec, 90.0 - geolat)
    pos_deg = (az_lon + CUSP_BOUNDARY_OFFSET) % 360.0
    return pos_deg / 30.0 + 1.0


def _hpos_carter(ra: float, armc: float, eps: float, geolat: float) -> float:
    """Carter "Poli-Equatorial" house position.

    Carter's system divides the celestial equator into twelve equal arcs of
    right ascension anchored at the right ascension of the ascendant
    (C.E.O. Carter, "Essays on the Foundations of Astrology"). The body's
    placement is simply its RA distance from the ascendant's RA.
    """
    asc = _calc_ascendant((armc + 90.0) % 360.0, eps, geolat, geolat)
    asc = _ascendant_on_eastern_horizon(asc, armc, eps, geolat)
    asc_ra, _ = _rotate_frame(asc, 0.0, -eps)
    pos_deg = (ra - asc_ra) % 360.0
    return pos_deg / 30.0 + 1.0


_Vec3 = tuple[float, float, float]

#: Half-width, in degrees of house-circle arc, of the rounding window
#: around the ascendant that the Krusinski scale treats as arc zero.
_HOUSE_ARC_SEAM = 1e-9


def _cross3(a: _Vec3, b: _Vec3) -> _Vec3:
    """Return the cross product of two equatorial direction vectors."""
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot3(a: _Vec3, b: _Vec3) -> float:
    """Return the dot product of two equatorial direction vectors."""
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _hpos_krusinski(ra: float, armc: float, eps: float, geolat: float) -> float:
    """Krusinski-Pisa-Goelzer house position.

    The system divides the great circle through the ascendant and the zenith
    into twelve arcs of 30 degrees, starting at the ascendant and running
    towards the nadir, and carries those divisions to the ecliptic along hour
    circles (B. Krusinski's 1995 description of his system, credited
    independently to M. Pisa and to G. Goelzer). A body is therefore placed by
    the arc, measured from the ascendant, of the point of that circle sharing
    its right ascension. That is the anatomy R. W. Holden, "The Elements of
    House Division" (1977), gives for the great-circle family — which circle is
    divided, where the division starts, through which pencil it reaches the
    ecliptic — with a fourth circle in place of the prime vertical of Campanus
    or the equator of Regiomontanus, and the ordinary pencil of hour circles.

    Because the projection runs along hour circles, the answer depends on the
    body only through its right ascension; its declination cannot enter.

    Args:
        ra: Right ascension of the body in degrees.
        armc: Right ascension of the MC of the chart in degrees.
        eps: True obliquity of the ecliptic in degrees.
        geolat: Geographic latitude in degrees, north positive.

    Returns:
        House position in ``[1, 13)``: the integer part names the house, the
        fraction is the part of that house already travelled.
    """
    # The zenith is the direction whose right ascension is the ARMC and whose
    # declination is the geographic latitude — that is what the ARMC is.
    cos_geolat = _cos_deg(geolat)
    zenith: _Vec3 = (
        cos_geolat * _cos_deg(armc),
        cos_geolat * _sin_deg(armc),
        _sin_deg(geolat),
    )
    # The scale runs from the ascendant towards the nadir, so the nadir — not
    # the zenith — is the quarter-turn mark of the circle and the second basis
    # direction below. It puts the IC on cusp 4 and the MC on cusp 10.
    nadir: _Vec3 = (-zenith[0], -zenith[1], -zenith[2])

    # The ascendant is the ecliptic point on the *eastern* horizon. The module
    # owns both the closed-form longitude and the rising-intersection test,
    # which inside the polar circles is what keeps the anchor off the setting
    # point; sharing them keeps this system aligned with every other one.
    asc_lon = _calc_ascendant((armc + 90.0) % 360.0, eps, geolat, geolat)
    asc_lon = _ascendant_on_eastern_horizon(asc_lon, armc, eps, geolat)
    sin_asc_lon = _sin_deg(asc_lon)
    ascendant: _Vec3 = (
        _cos_deg(asc_lon),
        sin_asc_lon * _cos_deg(eps),
        sin_asc_lon * _sin_deg(eps),
    )

    # The house circle is the great circle through the ascendant and the
    # zenith. The ascendant lies on the horizon, so it is already orthogonal
    # to the zenith: the pair is an orthonormal basis of the circle's plane
    # and their cross product is the pole the circle turns about.
    house_pole = _cross3(ascendant, nadir)

    # The body's hour circle is the great circle through the celestial poles
    # and the body; its own pole sits on the equator a quarter turn of right
    # ascension ahead of it. Two great circles meet in a pair of opposite
    # directions, and the wanted one is the member of the pair that shares the
    # body's right ascension rather than the opposite one — the single branch
    # choice of the construction, made once here and nowhere else.
    hour_pole: _Vec3 = (-_sin_deg(ra), _cos_deg(ra), 0.0)
    meeting = _cross3(house_pole, hour_pole)
    if meeting[0] * _cos_deg(ra) + meeting[1] * _sin_deg(ra) < 0.0:
        meeting = (-meeting[0], -meeting[1], -meeting[2])

    # Resolving that direction on the two basis vectors gives the arc from the
    # ascendant over the whole turn: ascendant 0, nadir 90, descendant 180,
    # zenith 270. At a geographic pole the house circle is itself an hour
    # circle, every other hour circle cuts it only at the celestial poles, and
    # this arc collapses to exactly 90 or 270 — the 4th or the 10th cusp.
    arc = math.degrees(math.atan2(_dot3(meeting, nadir), _dot3(meeting, ascendant)))
    # Guard the seam of the scale. The ascendant is arc zero and cusp one, and
    # its own component along the nadir cancels to a rounding residue whose
    # sign is arbitrary; left alone, a negative residue would wrap a body
    # sitting on the ascendant round to just short of a full turn and report it
    # in house 12 instead of house 1. The window is five orders of magnitude
    # above that residue and far below any arc a chart can resolve.
    if -_HOUSE_ARC_SEAM < arc < 0.0:
        arc = 0.0
    return (arc % 360.0) / 30.0 + 1.0


def _savard_prime_vertical_arc(geolat: float, thirds: int) -> float:
    """Arc along the prime vertical from a horizon point to a parallel.

    Parametrize the prime vertical by the arc ``t`` measured from the east
    point of the horizon toward the zenith. That circle crosses the celestial
    equator at the east point under the angle of the geographic latitude and
    reaches the zenith, whose declination is the latitude itself, so a point
    at arc ``t`` has ``sin(dec) = sin(t) * sin(phi)``. The parallel of
    declination that cuts off ``thirds`` thirds of the declination gained
    between the east point and the zenith is therefore met where

        sin(t) = sin(thirds * phi / 3) / sin(phi).

    On the equator numerator and denominator vanish together and the ratio
    has the finite limit ``thirds / 3``, the value taken there. The arc
    depends on the size of the latitude only, since both sines change sign
    with it.

    Args:
        geolat: Geographic latitude, degrees.
        thirds: Which parallel, 1 for one third or 2 for two thirds.

    Returns:
        The arc in degrees, strictly between 0 and 90.
    """
    sin_lat = _sin_deg(geolat)
    if sin_lat == 0.0:
        sin_arc = thirds / 3.0
    else:
        sin_arc = _sin_deg(thirds * geolat / 3.0) / sin_lat
    return _asin_deg(sin_arc)


def _hpos_savard(md_upper: float, dec: float, geolat: float) -> float:
    """Savard-A house position.

    Savard-A divides the prime vertical — the great circle through the east
    point of the horizon, the zenith, the west point and the nadir — by
    parallels of declination instead of into equal arcs. Between the east
    point and the zenith the declination of that circle runs from zero to the
    geographic latitude; the parallels at one third and two thirds of that
    span cut it, and the great circles through those crossings and the north
    and south points of the horizon are the intermediate house circles
    (John Savard, "Astrological House Systems",
    http://www.quadibloc.com/other/as01.htm, the Albategnus entry; the family
    of house circles through the horizon points is described in J. H. Holden,
    "The Elements of House Division", L. N. Fowler & Co., 1977). The cusp side
    of the system builds the same division.

    Every house circle of the family is perpendicular to the prime vertical
    and meets it in a single point, so a body is placed by where its own
    circle crosses that arc: its distance from the prime vertical counts as
    little here as ecliptic latitude counts for a whole-sign house. Referring
    the body to the east point of the equator and tilting the frame by the
    latitude yields that crossing, and the answer is the share of the arc
    between the two bracketing cusps the body has already covered (the usual
    interpolation rule for a house position, Holden 1977).

    Args:
        md_upper: Right ascension of the body less the ARMC of the chart, as
            a signed angle in degrees; it grows eastward.
        dec: Declination of the body, degrees.
        geolat: Geographic latitude, degrees.

    Returns:
        The house position in ``[1, 13)``: the house number plus the fraction
        of that house already behind the body.
    """
    # Longitude on the prime vertical, counted from the east point in the
    # sense east point -> nadir -> west point -> zenith, which is the order
    # the houses are numbered in. The east point lies a quarter turn east of
    # the meridian on the equator, and tilting the equator about that point
    # by the latitude carries it onto the prime vertical. A body exactly at
    # one of the horizon's north/south points lies on the pole of the prime
    # vertical and has no longitude on it; the rotation answers zero there,
    # which puts it on the first cusp.
    pv_lon, _ = _rotate_frame((md_upper - 90.0) % 360.0, dec, -geolat)
    # A longitude infinitesimally below zero can round up to exactly 360.0,
    # which is the first cusp again and not a thirteenth house.
    pv_lon = pv_lon % 360.0

    # The cusps on that arc. The four angles are fixed — ascendant at the east
    # point, IC at the nadir, descendant at the west point, MC at the zenith —
    # and in each quadrant the two intermediate cusps stand at the arcs above,
    # measured from whichever horizon point bounds the quadrant. Opposite
    # cusps share a house circle, so the twelve come from two arcs.
    first = _savard_prime_vertical_arc(geolat, 1)
    second = _savard_prime_vertical_arc(geolat, 2)
    cusps = (
        0.0,
        first,
        second,
        90.0,
        180.0 - second,
        180.0 - first,
        180.0,
        180.0 + first,
        180.0 + second,
        270.0,
        360.0 - second,
        360.0 - first,
        360.0,
    )

    # The twelve arcs are unequal away from the poles but always in
    # increasing order, so the body's house is the one interval that brackets
    # its longitude and the fraction is its progress across that interval.
    # A longitude past the last intermediate cusp belongs to the twelfth.
    house = 12
    for candidate in range(1, 12):
        if pv_lon < cusps[candidate]:
            house = candidate
            break
    lower = cusps[house - 1]
    upper = cusps[house]
    return float(house) + (pv_lon - lower) / (upper - lower)


def _hpos_sripati(lon: float, armc: float, eps: float, geolat: float) -> float:
    """Sripati house position.

    Sripati houses are the classical Porphyry trisection of the ecliptic
    quadrants with every cusp pulled back so it lies at the *midpoint* of
    the Porphyry house — equivalently, the Porphyry placement advanced by
    half a house (standard Hindu-astrology construction; see e.g.
    B.V. Raman, "A Manual of Hindu Astrology").

    The final half-house shift is wrapped continuously across house 12.
    """
    asc = _calc_ascendant((armc + 90.0) % 360.0, eps, geolat, geolat)
    mc = _armc_to_mc(armc, eps)
    asc = _ascendant_on_eastern_horizon(asc, armc, eps, geolat)

    arc_from_asc = (lon - asc) % 360.0
    arc_from_asc = (arc_from_asc + CUSP_BOUNDARY_OFFSET) % 360.0
    if arc_from_asc < 180.0:
        hpos = 1.0
    else:
        hpos = 7.0
        arc_from_asc -= 180.0
    quadrant = difdeg2n(asc, mc)
    if arc_from_asc < 180.0 - quadrant:
        hpos += arc_from_asc * 3.0 / (180.0 - quadrant)
    else:
        hpos += 3.0 + (arc_from_asc - 180.0 + quadrant) * 3.0 / quadrant
    hpos += 0.5

    return (hpos - 1.0) % 12.0 + 1.0


def _hpos_sunshine_apc(
    hsys_char: str,
    armc: float,
    geolat: float,
    dec: float,
    md_upper: float,
    eps: float,
) -> float:
    """Sunshine ('I'/'i') and APC ('Y') house positions.

    Both systems divide the diurnal and nocturnal arcs of a *reference
    declination* into six equal parts along the equator and place the body
    by proportional membership of those arcs:

      * Sunshine (Makransky, "Primary Directions"): the reference
        declination is the Sun's. ``house_pos()`` receives no Sun context,
        so the reference declination falls back to 0 by definition of this
        context-free entry point.
      * APC (Knegt): the reference declination is the ascendant's.

    Construction: start from the body's Regiomontanus position (the pole
    of its position circle on the celestial equator), split the sphere at
    the horizon of the reference declination, and rescale the arc between
    meridian and horizon by the reference semi-arc (diurnal above the
    horizon, nocturnal below). The spherical triangle between the meridian,
    the position circle and the reference-declination parallel supplies
    the offset that maps the body's position-circle angle onto the
    reference parallel before rescaling.
    """
    gl = geolat
    if gl > 90.0 - CUSP_BOUNDARY_OFFSET:
        gl = 90.0 - CUSP_BOUNDARY_OFFSET
    if gl < -90.0 + CUSP_BOUNDARY_OFFSET:
        gl = -90.0 + CUSP_BOUNDARY_OFFSET

    if hsys_char == "Y":
        asc = _calc_ascendant((armc + 90.0) % 360.0, eps, geolat, geolat)
        asc = _ascendant_on_eastern_horizon(asc, armc, eps, geolat)
        _, ref_dec = _rotate_frame(asc, 0.0, -eps)
    else:
        ref_dec = 0.0

    dec_b = dec
    if 90.0 - abs(dec_b) < _NEAR_ZERO:
        dec_b = 90.0 - _NEAR_ZERO if dec_b > 0 else -90.0 + _NEAR_ZERO

    merid_dist = md_upper
    sin_md = _sin_deg(merid_dist)
    if sin_md == 0.0:
        sin_md = 1e-300

    # Regiomontanus position of the body (pole of its position circle).
    alt_sign = _tan_deg(gl) * _tan_deg(dec_b) + _cos_deg(merid_dist)
    circle_arc = _atan_deg(-alt_sign / sin_md) % 360.0
    if merid_dist < 0:
        circle_arc += 180.0
    circle_arc = circle_arc % 360.0

    is_above_horizon = _above_horizon_test(dec_b, gl, merid_dist) >= 0

    # Angle between the celestial pole and the horizon along the meridian,
    # then the body's position-circle angle from the lower meridian.
    pole_arc = 90.0 - gl if gl >= 0 else 90.0 + gl
    from_lower_meridian = (circle_arc - 270.0) % 360.0
    is_western_half = False
    if from_lower_meridian > 180.0:
        is_western_half = True
        from_lower_meridian = 360.0 - from_lower_meridian

    # Semi-arcs of the reference declination.
    ref_ad = _asc_diff_saturated(ref_dec, gl)
    diurnal_arc = 90.0 + ref_ad
    nocturnal_arc = 90.0 - ref_ad

    if diurnal_arc == 0 and is_above_horizon:
        circle_arc = 270.0
    elif nocturnal_arc == 0 and not is_above_horizon:
        circle_arc = 90.0
    else:
        semi_arc = diurnal_arc
        ref_dec_side = ref_dec
        if not is_above_horizon:
            # Nocturnal side: mirror the reference parallel and the angle.
            ref_dec_side = -ref_dec_side
            semi_arc = nocturnal_arc
            from_lower_meridian = 180.0 - from_lower_meridian
            is_western_half = not is_western_half
        # Spherical triangle meridian / position circle / reference
        # parallel: side between pole and body direction, angle at the
        # body, and the offset carrying the position-circle angle onto the
        # reference parallel.
        side_a = _acos_deg(_cos_deg(pole_arc) * _cos_deg(from_lower_meridian))
        if side_a < _NEAR_ZERO:
            side_a = _NEAR_ZERO
        sin_psi = _sin_deg(pole_arc) / _sin_deg(side_a)
        sin_psi = max(-1.0, min(1.0, sin_psi))
        ratio = _sin_deg(ref_dec_side) / sin_psi
        if ratio > 1:
            ratio = 90.0 - _NEAR_ZERO
        elif ratio < -1:
            ratio = -(90.0 - _NEAR_ZERO)
        else:
            ratio = _asin_deg(ratio)
        offset = _acos_deg(_cos_deg(ratio) / _cos_deg(ref_dec_side))
        if ref_dec_side < 0:
            offset = -offset
        if gl < 0:
            offset = -offset
        from_lower_meridian += offset
        if is_western_half:
            circle_arc = 270.0 - (from_lower_meridian / semi_arc) * 90.0
        else:
            circle_arc = 270.0 + (from_lower_meridian / semi_arc) * 90.0
        if not is_above_horizon:
            circle_arc = (circle_arc + 180.0) % 360.0

    pos_deg = (circle_arc + CUSP_BOUNDARY_OFFSET) % 360.0
    return pos_deg / 30.0 + 1.0


@overload
def _house_pos_pythonic(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Tuple[float, float],
    lon_or_hsys: Union[bytes, str],
    lat_body: float = ...,
) -> float: ...


@overload
def _house_pos_pythonic(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Union[int, bytes, str],
    lon_or_hsys: float,
    lat_body: float = ...,
) -> float: ...


@overload
def _house_pos_pythonic(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Tuple[float, float],
    lon_or_hsys: None = ...,
    lat_body: float = ...,
) -> float: ...


@overload
def _house_pos_pythonic(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Union[int, bytes, str],
    lon_or_hsys: None = ...,
    lat_body: float = ...,
) -> float: ...


def _house_pos_pythonic(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Union[int, bytes, str, Tuple[float, float]],
    lon_or_hsys: Optional[Union[float, bytes, str]] = None,
    lat_body: float = 0.0,
) -> float:
    """
    Determine in which house a celestial body is located.

    Returns a decimal value where the integer part is the house number (1-12)
    and the decimal part indicates the position within the house
    (0.0 = start of cusp, 0.999... = end of house, just before next cusp).

    This function supports two calling conventions:

    1. Reference API compatible (5 args):
       _house_pos_pythonic(armc, lat, obliquity, objcoord, hsys)
       where objcoord is a tuple (lon, lat_body) and hsys is bytes/str

    2. Extended signature (6 args):
       _house_pos_pythonic(armc, lat, obliquity, hsys, lon, lat_body)
       where hsys is int/bytes/str and lon/lat_body are separate floats

    Args:
        armc: Right Ascension of Medium Coeli (ARMC) in degrees (0-360)
        lat: Geographic latitude in degrees (positive North, negative South)
        obliquity: True obliquity of the ecliptic in degrees
        hsys_or_objcoord: Either house system (int/bytes/str) or objcoord tuple (lon, lat)
        lon_or_hsys: Either body longitude (float) or house system (bytes/str)
        lat_body: Ecliptic latitude of the body in degrees (only for 6-arg form)

    Returns:
        Decimal value where:
            - Integer part (1-12): House number
            - Decimal part (0.0-0.999...): Position within house

    Example:
        >>> # Sun at 15° Aries, Placidus houses, Rome (5-arg reference API form)
        >>> pos = _house_pos_pythonic(292.957, 41.9, 23.4393, (15.0, 0.0), b'P')
        >>> # Or 6-arg extended form:
        >>> pos = _house_pos_pythonic(292.957, 41.9, 23.4393, ord('P'), 15.0, 0.0)
        >>> house = int(pos)  # House number (e.g., 10)
        >>> position = pos - house  # Position within house (e.g., 0.5 = halfway)
    """
    _NEAR_ZERO = 1e-10
    CUSP_BOUNDARY_OFFSET = 0.0

    # Declare typed local variables for proper type narrowing
    lon: float
    hsys_char: str
    hsys_int: int

    # Detect which calling convention is used.  objcoord may be any
    # sequence (lists are accepted as well as tuples).
    if isinstance(hsys_or_objcoord, (tuple, list)):
        # 5-arg reference API form: (armc, lat, obliquity, objcoord, hsys)
        objcoord = hsys_or_objcoord
        lon = objcoord[0]
        lat_body = objcoord[1] if len(objcoord) > 1 else 0.0
        # hsys comes from lon_or_hsys (bytes/str/int), default to b"P".
        # An int is a character code (the same convention as
        # houses(..., hsys=ord('P')) and the 6-arg form below); silently
        # falling back to Placidus here used to discard the requested
        # system for every int caller.
        if isinstance(lon_or_hsys, bytes):
            hsys_char = chr(lon_or_hsys[0])
            hsys_int = lon_or_hsys[0]
        elif isinstance(lon_or_hsys, str):
            hsys_char = lon_or_hsys[0]
            hsys_int = ord(lon_or_hsys[0])
        elif isinstance(lon_or_hsys, int) and not isinstance(lon_or_hsys, bool):
            hsys_char = _int_hsys_to_char(lon_or_hsys)
            hsys_int = lon_or_hsys
        else:
            # Default to Placidus (hsys omitted in the 4-arg form)
            hsys_char = "P"
            hsys_int = ord("P")
    else:
        # 6-arg form: (armc, lat, obliquity, hsys, lon, lat_body)
        # hsys comes from hsys_or_objcoord (int/bytes/str)
        if isinstance(hsys_or_objcoord, bytes):
            hsys_char = chr(hsys_or_objcoord[0])
            hsys_int = hsys_or_objcoord[0]
        elif isinstance(hsys_or_objcoord, str):
            hsys_char = hsys_or_objcoord[0]
            hsys_int = ord(hsys_or_objcoord[0])
        else:
            # int case
            hsys_char = _int_hsys_to_char(hsys_or_objcoord)
            hsys_int = hsys_or_objcoord
        # lon comes from lon_or_hsys: a float, or an objcoord sequence
        # in the hsys-first calling order
        if isinstance(lon_or_hsys, (tuple, list)):
            lon = float(lon_or_hsys[0])
            if len(lon_or_hsys) > 1:
                lat_body = float(lon_or_hsys[1])
        else:
            lon = float(lon_or_hsys) if lon_or_hsys is not None else 0.0

    # Case-fold the selector like every other entry point ('k' == 'K',
    # lowercase 'i' stays the distinct Sunshine-alternative system).
    # house_pos folds 'g' fully (compatibility contract): Gauquelin
    # SECTOR positions (1-36) are returned for the lowercase selector too —
    # only the houses() tuple shape is keyed on the uppercase byte.
    hsys_char = _fold_hsys_case(hsys_char)
    if hsys_char == "g":
        hsys_char = "G"
    hsys_int = ord(hsys_char[0]) if hsys_char else hsys_int

    # Normalize angular inputs
    armc %= 360.0
    lon = lon % 360.0
    eps = obliquity
    geolat = lat

    # Convert ecliptic coordinates to equatorial via _rotate_spherical_x_axis rotation
    # (rotation around x-axis by obliquity angle)
    # This is crucial for proper house position when body has non-zero latitude
    lon_rad = math.radians(lon)
    lat_body_rad = math.radians(lat_body)
    eps_rad = math.radians(eps)

    # Ecliptic to equatorial transformation
    # RA: tan(RA) = (sin(lon)*cos(eps) - tan(lat)*sin(eps)) / cos(lon)
    # Dec: sin(Dec) = sin(lat)*cos(eps) + cos(lat)*sin(lon)*sin(eps)
    sin_lon = math.sin(lon_rad)
    cos_lon = math.cos(lon_rad)
    sin_lat = math.sin(lat_body_rad)
    cos_lat = math.cos(lat_body_rad)
    sin_eps = math.sin(eps_rad)
    cos_eps = math.cos(eps_rad)

    # Right Ascension
    y_ra = sin_lon * cos_eps - math.tan(lat_body_rad) * sin_eps
    x_ra = cos_lon
    ra = math.degrees(math.atan2(y_ra, x_ra)) % 360.0

    # Declination
    sin_dec = sin_lat * cos_eps + cos_lat * sin_lon * sin_eps
    if sin_dec > 1.0:
        sin_dec = 1.0
    if sin_dec < -1.0:
        sin_dec = -1.0
    dec = math.degrees(math.asin(sin_dec))

    # Signed meridian distances of the body, in (-180, 180]:
    #   md_upper — RA distance from the upper meridian (culmination, RA = ARMC)
    #   md_lower — RA distance from the lower meridian (anti-culmination)
    md_upper = (ra - armc) % 360.0
    if md_upper >= 180:
        md_upper -= 360.0
    md_lower = (md_upper + 180.0) % 360.0
    if md_lower >= 180:
        md_lower -= 360.0

    # Handle different house systems
    if hsys_char == "P" or hsys_char == "G":
        # Placidus / Gauquelin house position.
        # The Placidus position of a body is the fraction of its own
        # semi-arc already traversed, mapped onto a 360° house circle
        # (Holden, 'The Elements of House Division', 1977):
        #
        #   above horizon: position = (md_upper / day_semi_arc + 3) * 90°
        #   below horizon: position = (md_lower / night_semi_arc + 1) * 90°

        if 90.0 - abs(dec) <= abs(geolat):
            # Circumpolar body: it never rises or sets, so no semi-arcs
            # exist.  A never-rising body is placed in the nocturnal half
            # by its half-weighted distance from lower culmination; a
            # never-setting body in the diurnal half from upper culmination.
            if dec * geolat < 0:
                pos_deg = (90.0 + md_lower / 2.0) % 360.0
            else:
                pos_deg = (270.0 + md_upper / 2.0) % 360.0
        else:
            # Ascensional difference of the body: sin(AD) = tan δ · tan φ
            sin_ad = math.tan(math.radians(dec)) * math.tan(math.radians(geolat))
            if abs(sin_ad) > 1.0:
                sin_ad = 1.0 if sin_ad > 0 else -1.0
            asc_diff = math.degrees(math.asin(sin_ad))

            # Above-horizon test.  The altitude h of the body satisfies
            #   sin h = sin φ sin δ + cos φ cos δ cos H;
            # dividing by cos φ cos δ (> 0) gives the equivalent sign test
            #   tan φ tan δ + cos H  =  sin(AD) + cos H  >=  0.
            horizon_test = sin_ad + math.cos(math.radians(md_upper))
            above_horizon = horizon_test >= 0

            day_semi_arc = 90.0 + asc_diff
            night_semi_arc = 90.0 - asc_diff

            if above_horizon:
                pos_deg = (md_upper / day_semi_arc + 3.0) * 90.0
            else:
                pos_deg = (md_lower / night_semi_arc + 1.0) * 90.0

            # Add small offset for cusp precision
            pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET) % 360.0

        if hsys_char == "G":
            # Gauquelin sectors run clockwise. Wrap after reversing so a
            # pos_deg of 0 gives sector 1.0 instead of 37.0 (out of domain).
            pos_deg = (360.0 - pos_deg) % 360.0
            hpos = pos_deg / 10.0 + 1.0
        else:
            hpos = pos_deg / 30.0 + 1.0

        return hpos

    elif hsys_char == "R":
        # Regiomontanus house position.
        #
        # Regiomontanus divides the celestial equator into equal arcs by
        # position circles through the North and South points of the
        # horizon; a body's position is the equator point cut by its own
        # position circle.  Intersecting the plane of the circle through
        # the horizon N/S points and the body (meridian distance M,
        # declination δ) with the equator gives the closed form
        #
        #   tan x = -(tan φ tan δ + cos M) / sin M
        #
        # measured so that x = 270° on the upper meridian (M = 0) and
        # x = 90° on the lower (|M| = 180°) — those two degenerate cases
        # are handled directly.
        # Degenerate meridian cases: the equator intersection is 270
        # (house 10 side) or 90 (house 4 side) according to whether the
        # meridian point is above or below the horizon — at polar
        # latitudes the upper-meridian point can be below it (altitude
        # 90 - |lat - dec| < 0 ⟺ tan(lat) tan(dec) < -1) and the limit
        # of the closed form lands on the opposite branch.
        _tt = math.tan(math.radians(geolat)) * math.tan(math.radians(dec))
        if abs(md_upper) < _NEAR_ZERO:
            pos_deg = 270.0 if _tt > -1.0 else 90.0
        elif 180.0 - abs(md_upper) < _NEAR_ZERO:
            pos_deg = 90.0 if _tt < 1.0 else 270.0
        else:
            if 90.0 - abs(geolat) < _NEAR_ZERO:
                geolat = 90.0 - _NEAR_ZERO if geolat > 0 else -90.0 + _NEAR_ZERO
            if 90.0 - abs(dec) < _NEAR_ZERO:
                dec = 90.0 - _NEAR_ZERO if dec > 0 else -90.0 + _NEAR_ZERO

            pole_term = math.tan(math.radians(geolat)) * math.tan(
                math.radians(dec)
            ) + math.cos(math.radians(md_upper))
            pos_deg = math.degrees(
                math.atan(-pole_term / math.sin(math.radians(md_upper)))
            )
            if md_upper < 0:
                pos_deg += 180.0
            pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET) % 360.0

        # Python's float modulo can round a value infinitesimally below zero
        # to exactly 360.0.  The angles are identical geometrically, but the
        # latter would encode the nonexistent decimal house 13.0 rather than
        # the first cusp.  Keep the public house-position interval half-open.
        if pos_deg >= 360.0:
            pos_deg -= 360.0

        hpos = pos_deg / 30.0 + 1.0
        return hpos

    elif hsys_char == "C":
        # Campanus house position.
        # The Campanus system divides the prime vertical into 12 equal
        # 30° arcs.  A body's position follows from rotating its
        # equatorial coordinates into the frame whose fundamental plane
        # is the prime vertical (pole = horizon North point):
        #
        # 1. Express the body as (meridian distance - 90°, declination),
        #    i.e. longitude measured from the east point of the equator.
        # 2. Rotate about the shared x-axis (the east point) by -latitude,
        #    which carries the equator onto the prime vertical.
        # 3. The longitude in the rotated frame is the Campanus position.
        #
        # Reference: Meeus, 'Astronomical Algorithms' 2nd ed., Ch. 13;
        # Montenbruck & Pfleger, par. 1.3 (x-axis rotation matrix).

        ha_east = md_upper - 90.0

        ha_east_rad = math.radians(ha_east)
        dec_rad = math.radians(dec)
        rot_angle = math.radians(-geolat)

        cos_rot = math.cos(rot_angle)
        sin_rot = math.sin(rot_angle)
        cos_ha = math.cos(ha_east_rad)
        sin_ha = math.sin(ha_east_rad)
        cos_dec = math.cos(dec_rad)
        sin_dec_b = math.sin(dec_rad)

        # Convert to Cartesian
        vec_x = cos_dec * cos_ha
        vec_y = cos_dec * sin_ha
        vec_z = sin_dec_b

        # Rotation around x-axis:
        #   x' = x;  y' = y cos θ + z sin θ;  z' = -y sin θ + z cos θ
        rot_x = vec_x
        rot_y = vec_y * cos_rot + vec_z * sin_rot

        # Convert back to spherical (longitude only needed)
        pos_deg = math.degrees(math.atan2(rot_y, rot_x))
        pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET) % 360.0
        hpos = pos_deg / 30.0 + 1.0
        return hpos

    elif hsys_char in ["E", "A", "W", "V", "D", "N"]:
        # Equal-based systems - simple longitude-based calculation
        asc = _calc_ascendant((armc + 90.0) % 360.0, eps, geolat, geolat)
        mc = _armc_to_mc(armc, eps)
        # Keep the ascendant on the eastern horizon. Without this correction
        # the raw ecliptic-longitude ascendant lands on the western half
        # beyond the polar circle, giving a six-house error for E/A/W/V.
        asc = _ascendant_on_eastern_horizon(asc, armc, eps, geolat)

        if hsys_char == "D":
            pos_deg = (lon - mc - 90.0) % 360.0
        elif hsys_char == "V":
            pos_deg = (lon - asc + 15.0) % 360.0
        elif hsys_char == "W":
            pos_deg = (lon - asc + (asc % 30.0)) % 360.0
        elif hsys_char == "N":
            pos_deg = lon
        else:
            pos_deg = (lon - asc) % 360.0

        pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET) % 360.0
        hpos = pos_deg / 30.0 + 1.0
        return hpos

    elif hsys_char == "X":
        # Meridian (axial rotation)
        hpos = ((md_upper - 90.0) % 360.0) / 30.0 + 1.0
        return hpos

    elif hsys_char == "M":
        # Morinus
        # Project the ecliptic longitude onto the equatorial frame via:
        # tan(ra_equiv) = tan(λ) / cos(ε)
        # atan2(sin(λ), cos(λ)·cos(ε)) gives the RA-equivalent in [-180, 180],
        # then % 360 maps to [0, 360).
        a = lon
        if abs(a - 90.0) > _NEAR_ZERO and abs(a - 270.0) > _NEAR_ZERO:
            hpos_deg = (
                math.degrees(
                    math.atan2(
                        math.sin(math.radians(a)), math.cos(math.radians(a)) * cos_eps
                    )
                )
                % 360.0
            )
        else:
            hpos_deg = 90.0 if abs(a - 90.0) <= _NEAR_ZERO else 270.0
        hpos_deg = (hpos_deg - armc - 90.0) % 360.0
        hpos = hpos_deg / 30.0 + 1.0
        return hpos

    elif hsys_char == "K":
        # Koch (Birthplace / GOH) house position.
        # The Koch system defines every house circle by the *birthplace
        # horizon*: each cusp rises after the MC by an equal share of the
        # MC degree's own semi-arc.  A body's position is therefore the
        # share of the MC's semi-arc that has elapsed between the body's
        # rising (its meridian distance corrected by its own ascensional
        # difference) and the MC's rising (corrected by the MC's AD):
        #
        #   ad_body      = arcsin(tan φ · tan δ_body)
        #   ad_mc        = arcsin(tan ε · tan φ · sin ARMC)
        #                  (δ_MC from sin δ = sin ε · sin λ_MC, folded into
        #                   one expression via the MC's RA = ARMC)
        #   mc_semi_arc  = 90° + ad_mc
        #
        #   east  (md_upper >= 0): fraction = (md_upper - ad_body + ad_mc)
        #                                       / mc_semi_arc
        #                          position = (fraction - 1) · 90°
        #   west  (md_upper <  0): fraction = (md_upper + 180 + ad_body
        #                                       + ad_mc) / mc_semi_arc
        #                          position = (fraction + 1) · 90°
        #
        # Bodies whose fraction falls outside [0, 2] never cross the
        # birthplace horizon inside this quadrant pair — the Koch position
        # is undefined and 0.0 is returned (matching the reference API).

        tan_eps = math.tan(math.radians(eps))

        invalid = False

        # Body's ascensional difference, saturated at ±90° for declinations
        # that never rise / never set at this latitude.
        if 90.0 - geolat < dec or -90.0 - geolat > dec:
            ad_body = 90.0
        elif geolat - 90.0 > dec or geolat + 90.0 < dec:
            ad_body = -90.0
        else:
            ad_arg = math.tan(math.radians(geolat)) * math.tan(math.radians(dec))
            ad_arg = max(-1.0, min(1.0, ad_arg))
            ad_body = math.degrees(math.asin(ad_arg))
            if abs(ad_body) < 1e-12:
                # An ascensional difference of ~1e-15 deg is float noise
                # from a declination that is physically zero (body on
                # the equinoctial points); its sign must not decide the
                # validity boundary below.
                ad_body = 0.0

        # MC's ascensional difference
        ad_mc_arg = (
            tan_eps * math.tan(math.radians(geolat)) * math.sin(math.radians(armc))
        )
        if abs(ad_mc_arg) > 1.0:
            ad_mc_arg = 1.0 if ad_mc_arg > 0 else -1.0
        ad_mc = math.degrees(math.asin(ad_mc_arg))
        if abs(ad_mc) < 1e-12:
            # Same noise floor as ad_body: sin(180 deg) = 1.2e-16 makes
            # the MC's AD ~1e-15 deg at cardinal ARMC; physically zero.
            ad_mc = 0.0

        # MC's semi-arc
        mc_semi_arc = 90.0 + ad_mc
        if mc_semi_arc == 0.0:
            invalid = True

        pos_deg = 0.0
        if not invalid and abs(mc_semi_arc) > 0:
            # The valid fraction range is the closed interval [0, 2]: both
            # endpoints are reached by real bodies.  A body exactly on the
            # meridian can land infinitesimally outside an endpoint because of
            # binary64 roundoff. Accept a small numerical tolerance and clamp
            # back into the geometrically valid interval; fractions genuinely
            # outside it remain undefined.
            _FRACTION_TOL = 1e-9
            if md_upper >= 0:  # east
                arc_fraction = (md_upper - ad_body + ad_mc) / mc_semi_arc
            else:  # west
                md_from_lower = md_upper + 180.0
                arc_fraction = (md_from_lower + ad_body + ad_mc) / mc_semi_arc
            if arc_fraction < -_FRACTION_TOL or arc_fraction > 2.0 + _FRACTION_TOL:
                invalid = True
            else:
                arc_fraction = min(max(arc_fraction, 0.0), 2.0)
                if md_upper >= 0:  # east
                    pos_deg = (arc_fraction - 1.0) * 90.0
                else:  # west
                    pos_deg = (arc_fraction + 1.0) * 90.0

        if invalid:
            # Koch position undefined in the circumpolar area — return 0.0
            # matching the reference API behavior
            return 0.0
        else:
            pos_deg = pos_deg % 360.0
            pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET) % 360.0
            hpos = pos_deg / 30.0 + 1.0
            return hpos

    elif hsys_char == "T":
        return _hpos_topocentric(armc, geolat, ra, dec, md_upper)

    elif hsys_char == "B":
        # Alcabitius house position.
        # Alcabitius cusps are equally spaced in RA within each quadrant.
        # For _house_pos_pythonic, we interpolate the body's RA between the cusp RAs.
        cusps, _ = houses_armc(armc, lat, obliquity, hsys_int)

        # Convert cusp ecliptic longitudes to RA
        cusp_ras = []
        for c in cusps:
            c_r = math.radians(c)
            y_c = math.sin(c_r) * cos_eps
            x_c = math.cos(c_r)
            cusp_ras.append(math.degrees(math.atan2(y_c, x_c)) % 360.0)

        # Find which RA interval contains the body's RA
        for i in range(12):
            ra_start = cusp_ras[i]
            ra_end = cusp_ras[(i + 1) % 12]

            diff_to_body = (ra - ra_start + 360.0) % 360.0
            interval_size = (ra_end - ra_start + 360.0) % 360.0

            if interval_size < 0.0001:
                interval_size = 30.0

            if diff_to_body < interval_size:
                house_num = i + 1
                fraction = diff_to_body / interval_size
                fraction = max(0.0, min(fraction, 0.9999999999))
                return float(house_num) + fraction

        return 1.0

    elif hsys_char == "H":
        return _hpos_horizon(md_upper, dec, geolat)

    elif hsys_char == "F":
        return _hpos_carter(ra, armc, eps, geolat)

    elif hsys_char == "U":
        return _hpos_krusinski(ra, armc, eps, geolat)

    elif hsys_char == "J":
        return _hpos_savard(md_upper, dec, geolat)

    elif hsys_char == "S":
        return _hpos_sripati(lon, armc, eps, geolat)

    elif hsys_char in ("I", "i", "Y"):
        return _hpos_sunshine_apc(hsys_char, armc, geolat, dec, md_upper, eps)

    # Default fallback: use cusp-based method (ecliptic longitude interpolation)
    # This applies to O, L, Q, and any other system whose reference house
    # position is a simple distance-from-cusp / house-size interpolation.
    cusps, ascmc = houses_armc(armc, lat, obliquity, hsys_int)

    # Find the house containing this longitude
    for i in range(12):
        cusp_start = cusps[i]
        cusp_end = cusps[(i + 1) % 12]

        diff_to_body = (lon - cusp_start + 360.0) % 360.0
        house_size = (cusp_end - cusp_start + 360.0) % 360.0

        if house_size < 1e-9:
            # Zero-width house (coincident cusps happen at extreme
            # latitudes, e.g. Pullen SD at lat 89 collapses houses 2-3).
            # A body can only be "in" it when it sits exactly on the
            # cusp; treating the empty house as 30 deg wide used to
            # swallow a whole sign of sky into it.
            if diff_to_body < 1e-9:
                return float(i + 1)
            continue

        if diff_to_body < house_size:
            house_num = i + 1
            fraction = diff_to_body / house_size
            fraction = max(0.0, min(fraction, 0.9999999999))
            return float(house_num) + fraction

    return 1.0


def _armc_to_mc(armc: float, eps: float) -> float:
    """Convert ARMC to MC longitude."""
    mc_rad = math.atan2(math.tan(math.radians(armc)), math.cos(math.radians(eps)))
    mc = math.degrees(mc_rad)
    if mc < 0:
        mc += 360.0
    if 90.0 < armc <= 270.0:
        # Match the canonical houses_armc() quadrant correction (<= boundary).
        if mc <= 90.0 or mc > 270.0:
            mc += 180.0
    elif armc > 270.0:
        if mc <= 270.0:
            mc += 180.0
    elif armc <= 90.0:
        if mc > 90.0:
            mc += 180.0
    return mc % 360.0


@overload
def house_pos(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Tuple[float, float],
    lon_or_hsys: Union[bytes, str],
    lat_body: float = ...,
) -> float: ...


@overload
def house_pos(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Union[int, bytes, str],
    lon_or_hsys: float,
    lat_body: float = ...,
) -> float: ...


@overload
def house_pos(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Tuple[float, float],
    lon_or_hsys: None = ...,
    lat_body: float = ...,
) -> float: ...


@overload
def house_pos(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Union[int, bytes, str],
    lon_or_hsys: None = ...,
    lat_body: float = ...,
) -> float: ...


def house_pos(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Union[int, bytes, str, Tuple[float, float]],
    lon_or_hsys: Optional[Union[float, bytes, str]] = None,
    lat_body: float = 0.0,
) -> float:
    """
    Determine in which house a celestial body is located.

    This function supports two calling conventions:

    1. Reference API compatible (5 args):
       house_pos(armc, lat, obliquity, objcoord, hsys)
       where objcoord is a tuple (lon, lat_body) and hsys is bytes/str

    2. Extended signature (6 args):
       house_pos(armc, lat, obliquity, hsys, lon, lat_body)
       where hsys is int/bytes/str and lon/lat_body are separate floats

    Returns a decimal value where the integer part is the house number (1-12)
    and the decimal part indicates the position within the house.

    Args:
        armc: Right Ascension of Medium Coeli (ARMC) in degrees
        lat: Geographic latitude in degrees
        obliquity: True obliquity of the ecliptic in degrees
        hsys_or_objcoord: Either house system (int/bytes/str) or objcoord tuple
        lon_or_hsys: Either body longitude (float) or house system (bytes/str)
        lat_body: Body's ecliptic latitude in degrees (only for 6-arg form)

    Returns:
        Decimal value: integer part = house number, decimal part = position within house
    """
    # Delegate to _house_pos_pythonic which now supports both calling conventions
    # Use type: ignore since the overloads cover all valid combinations
    return _house_pos_pythonic(
        armc,
        lat,
        obliquity,
        hsys_or_objcoord,  # type: ignore[arg-type]
        lon_or_hsys,  # type: ignore[arg-type]
        lat_body,
    )


#: ``gauquelin_sector`` methods 2-5 time the sector on the body's actual rising
#: and setting; the value is the ``rise_trans`` flag set selecting which point
#: of the disc is timed against which horizon. ``rise_trans`` defaults to the
#: upper limb on the refracted horizon, so method 5 adds no flag.
_GAUQUELIN_RISE_SET_FLAGS: dict[int, int] = {
    2: BIT_DISC_CENTER | BIT_NO_REFRACTION,  # disc centre, geometric horizon
    3: BIT_DISC_CENTER,  # disc centre, refracted horizon
    4: BIT_NO_REFRACTION,  # upper limb, geometric horizon
    5: 0,  # upper limb, refracted horizon
}

#: Methods 0 and 1 place the body on its semi-arcs by hour angle and
#: declination (the 'G' house-position construction); method 1 first drops the
#: body's ecliptic latitude, i.e. it sectors the ecliptic point of the body.
_GAUQUELIN_DROPS_LATITUDE: dict[int, bool] = {0: False, 1: True}


def _gauquelin_sector_from_rise_set(
    jd: float,
    planet: "int | str",
    lat: float,
    lon: float,
    altitude: float,
    pressure: float,
    temperature: float,
    flags: int,
    method: int,
) -> float:
    """
    Calculate Gauquelin sector using actual rise/set times.

    This implements methods 2-5 which use real rise/set event times
    rather than the hour angle approximation.

    Args:
        jd: Julian Day in UT
        planet: Planet ID (int) or fixed-star name (str); rise_trans accepts
            both.
        lat, lon: Observer location
        altitude, pressure, temperature: Atmospheric parameters
        flags: Calculation flags
        method: 2=disc center no refraction, 3=disc center with refraction,
                4=disc edge no refraction, 5=disc edge with refraction

    Returns:
        Sector position in range [1, 37)

    Raises:
        Error: If the body never rises or sets at the given latitude
            (circumpolar); the reference refuses these cases as well.
    """
    from .eclipse import rise_trans

    def _rise_set_not_found() -> Error:
        """No rise/set event exists (circumpolar / never-rises body).

        Methods 2-5 are defined only from real rise/set times, so — like the
        reference API, which refuses here rather than silently substituting
        a different method — signal the caller with an error instead of
        falling back to the method-0 hour-angle approximation.
        """
        return Error(
            f"gauquelin_sector: body {planet} has no rise or set at this "
            "location, so methods 2-5 cannot be evaluated"
        )

    # Determine rise_trans flags based on method
    rsmi_flags = _GAUQUELIN_RISE_SET_FLAGS[method]

    def _event_after(jd_from: float, calc_flag: int) -> float:
        """Time of the first rise (CALC_RISE) or set (CALC_SET) after jd_from."""
        retflag, tret = rise_trans(
            jd_from,
            planet,
            calc_flag | rsmi_flags,
            [lon, lat, altitude],
            pressure,
            temperature,
            flags,
        )
        if retflag == -2:
            raise _rise_set_not_found()  # circumpolar: no such event
        return tret[0]

    def _event_before(jd_next: float, calc_flag: int) -> float:
        """The last event of that kind before ``jd``, found from the next one.

        A body rises and sets about once a day, so the previous event is
        searched from 1.5 days before the next one, and from 2.5 days before
        it if that search still lands on or after ``jd``.
        """
        jd_prev = _event_after(jd_next - 1.5, calc_flag)
        if jd_prev >= jd:
            jd_prev = _event_after(jd_next - 2.5, calc_flag)
            if jd_prev >= jd:
                raise _rise_set_not_found()
        return jd_prev

    try:
        jd_next_rise = _event_after(jd, CALC_RISE)
        jd_next_set = _event_after(jd, CALC_SET)
        jd_prev_rise = _event_before(jd_next_rise, CALC_RISE)
        jd_prev_set = _event_before(jd_next_set, CALC_SET)
    except (IndexError, TypeError, ValueError, ArithmeticError):
        raise _rise_set_not_found()

    def _sector_on_arc(arc_start: float, arc_end: float, first_sector: float) -> float:
        """Position of ``jd`` on an arc divided into 18 equal parts of time."""
        arc = arc_end - arc_start
        if arc <= 0:
            raise _rise_set_not_found()
        fraction = (jd - arc_start) / arc
        fraction = max(0.0, min(1.0, fraction))  # clamp for safety
        return first_sector + fraction * 18.0

    # Gauquelin's sectors (Gauquelin & Gauquelin 1957, "Méthodes pour étudier
    # la répartition des astres dans le mouvement diurne"): the diurnal arc
    # from rising to setting and the nocturnal arc from setting to rising are
    # each divided into 18 equal parts, numbered from 1 at rising (10 at the
    # upper culmination) and from 19 at setting (28 at the lower culmination).
    # The more recent of the two last events tells which arc the body is on.
    if jd_prev_rise > jd_prev_set:
        sector = _sector_on_arc(jd_prev_rise, jd_next_set, 1.0)
    else:
        sector = _sector_on_arc(jd_prev_set, jd_next_rise, 19.0)

    # Normalize to range [1, 37)
    if sector >= 37.0:
        sector -= 36.0
    elif sector < 1.0:
        sector += 36.0

    return sector


def _gauquelin_sector_pythonic(
    jd: float,
    planet: "int | str",
    lat: float,
    lon: float,
    altitude: float = 0.0,
    pressure: float = 1013.25,
    temperature: float = 15.0,
    flags: int = 0,
    method: int = 0,
) -> float:
    """
    Calculate the Gauquelin sector (1-36) in which a planet is located.

    Gauquelin sectors are a 36-fold division of the diurnal and nocturnal arcs
    used in statistical astrology research by Michel Gauquelin. Sectors are
    numbered clockwise from the Ascendant:
    - Sector 1: Rising (Ascendant)
    - Sector 10: Upper culmination (MC)
    - Sector 19: Setting (Descendant)
    - Sector 28: Lower culmination (IC)

    Sector definition: Gauquelin & Gauquelin (1957), "Méthodes pour étudier la
    répartition des astres dans le mouvement diurne"; the diurnal and the
    nocturnal arc are each divided into 18 equal parts. The ``method``
    argument selects how the arcs are measured: :data:`_GAUQUELIN_RISE_SET_FLAGS`
    lists the methods timed on actual rise and set events and the disc/horizon
    they use, :data:`_GAUQUELIN_DROPS_LATITUDE` the two hour-angle methods.

    Reference API compatible function (gauquelin_sector equivalent).

    Args:
        jd: Julian Day in Universal Time (UT)
        planet: Planet/body ID (int, e.g. SUN, MOON) or star name (str)
        lat: Geographic latitude in degrees (positive North)
        lon: Geographic longitude in degrees (positive East)
        altitude: Geographic altitude in meters above sea level (default 0.0)
        pressure: Atmospheric pressure in mbar (default 1013.25)
        temperature: Atmospheric temperature in degrees Celsius (default 15.0)
        flags: Calculation flags (FLG_SWIEPH, FLG_TOPOCTR, etc.)
            The public wrapper defaults to FLG_SWIEPH | FLG_TOPOCTR.
        method: Calculation method:
            - 0: with latitude
            - 1: without latitude
            - 2: from rising/setting times of disc center
            - 3: from rising/setting times of disc center, with refraction
            - 4: from rising/setting times of disc edge
            - 5: from rising/setting times of disc edge, with refraction

    Returns:
        Sector position as float in range [1, 37).
        Integer part is sector number (1-36), decimal part is position within sector.

    Note:
        Methods 2-5 use actual rise/set times calculated via rise_trans().
        These methods may be slower than methods 0-1 but provide more
        accurate sector positions based on real rise/set events. A star name
        (str) is accepted for these methods as well as a planet id.

        For circumpolar objects (that never rise or set at the given
        latitude), methods 2-5 raise an Error, matching the reference API,
        rather than silently falling back to method 0.

    Example:
        >>> sector = _gauquelin_sector_pythonic(2451545.0, MARS, 48.85, 2.35)
        >>> print(f"Mars is in sector {int(sector)}")
    """
    # Methods 2-5: actual rise/set times. rise_trans accepts a star name (str)
    # as well as a planet id, so _gauquelin_sector_from_rise_set handles both.
    if method in _GAUQUELIN_RISE_SET_FLAGS:
        return _gauquelin_sector_from_rise_set(
            jd, planet, lat, lon, altitude, pressure, temperature, flags, method
        )

    # Methods 0-1: Use _house_pos_pythonic with Gauquelin house system ('G')
    # This computes Gauquelin sectors from house position
    from .planets import calc_ut
    from .fixed_stars import fixstar_ut

    # ARMC (sidereal time at location) and true obliquity from the long-term model.
    armc0_deg, eps = _house_armc_obliquity(jd)
    armc_deg = (armc0_deg + lon) % 360.0

    # Get body position - planet (int) or star (str).
    # The reference API computes the body topocentrically from the geopos
    # argument (default flags include FLG_TOPOCTR): set the observer for
    # the sub-call and restore the caller's topo afterwards.
    from . import state as _state
    from .state import get_topo, set_topo

    saved_topo = get_topo()
    try:
        if flags & FLG_TOPOCTR:
            set_topo(lon, lat, altitude)
        if isinstance(planet, str):
            pos, _star_name, retflag = fixstar_ut(planet, jd, flags | FLG_SPEED)
            planet_lon = float(pos[0])
            planet_lat = float(pos[1])
        else:
            pos, retflag = calc_ut(jd, planet, flags | FLG_SPEED)
            planet_lon = pos[0]
            planet_lat = pos[1]
    finally:
        if flags & FLG_TOPOCTR:
            # Restore the exact previous observer object (or None) so the
            # caller's topo state and observer-cache identity are untouched.
            with _state._STATE_LOCK:
                _state._TOPO = saved_topo

    if _GAUQUELIN_DROPS_LATITUDE.get(method, False):
        planet_lat = 0.0

    # Use _house_pos_pythonic with 'G' (Gauquelin sectors) to get the sector position
    # _house_pos_pythonic returns values in [1, 37) for Gauquelin
    sector = _house_pos_pythonic(armc_deg, lat, eps, (planet_lon, planet_lat), "G")

    return sector


def gauquelin_sector(
    tjdut: float,
    body: "int | str",
    method: int,
    geopos: tuple[float, float, float],
    atpress: float = 0.0,
    attemp: float = 0.0,
    flags: int = FLG_SWIEPH | FLG_TOPOCTR,
) -> float:
    """
    Calculate the Gauquelin sector (1-36) in which a planet is located.

    Gauquelin sectors are a 36-fold division of the diurnal and nocturnal arcs
    used in statistical astrology research by Michel Gauquelin. Sectors are
    numbered clockwise from the Ascendant:
    - Sector 1: Rising (Ascendant)
    - Sector 10: Upper culmination (MC)
    - Sector 19: Setting (Descendant)
    - Sector 28: Lower culmination (IC)

    Reference API compatible function (gauquelin_sector equivalent).

    Args:
        jd: Julian Day in Universal Time (UT)
        planet: Planet/body ID (SUN, MOON, etc.)
        lat: Geographic latitude in degrees (positive North)
        lon: Geographic longitude in degrees (positive East)
        altitude: Geographic altitude in meters above sea level (default 0.0)
        pressure: Atmospheric pressure in mbar (default 1013.25)
        temperature: Atmospheric temperature in degrees Celsius (default 15.0)
        flags: Calculation flags (FLG_SWIEPH, FLG_TOPOCTR, etc.)
            The public wrapper defaults to FLG_SWIEPH | FLG_TOPOCTR.
        method: Calculation method:
            - 0: with latitude
            - 1: without latitude
            - 2: from rising/setting times of disc center of planet
            - 3: from rising/setting times of disc center, incl. refraction
            - 4: from rising/setting times of disk edge of planet
            - 5: from rising/setting times of disk edge, incl. refraction
        geopos: Tuple of (longitude, latitude, altitude):
            - longitude: geographic longitude in degrees (eastern positive)
            - latitude: geographic latitude in degrees (northern positive)
            - altitude: geographic altitude in meters above sea level
        atpress: Atmospheric pressure in mbar (if 0, default 1013.25 is used)
        attemp: Atmospheric temperature in degrees Celsius (if 0, default 15 is used)
        flags: Bit flags for ephemeris (FLG_SWIEPH, etc.)

    Returns:
        Sector position as float in range [1, 37).
        Integer part is sector number (1-36), decimal part is position within sector.

    Note:
        This function matches the gauquelin_sector() API signature.
        Sectors are numbered clockwise from the Ascendant:
        - Sector 1: Rising (Ascendant)
        - Sector 10: Upper culmination (MC)
        - Sector 19: Setting (Descendant)
        - Sector 28: Lower culmination (IC)

    Example:
        >>> geopos = (2.35, 48.85, 0.0)  # Paris: (lon, lat, alt)
        >>> sector = gauquelin_sector(2451545.0, MARS, 0, geopos)
    """
    lon, lat, altitude = geopos
    # Use defaults if 0 is passed (standard API behavior)
    pressure = atpress if atpress != 0.0 else 1013.25
    temperature = attemp if attemp != 0.0 else 15.0

    return _gauquelin_sector_pythonic(
        tjdut, body, lat, lon, altitude, pressure, temperature, flags, method
    )
