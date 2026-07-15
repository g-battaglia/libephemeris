# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Arabic Parts (Lots) calculations for libephemeris.

Arabic Parts are classical astrological points calculated from combinations
of planetary positions and angles. They represent specific life areas and
are used in Hellenistic, Medieval, and modern astrology.

The reversible Fortune/Spirit rules map directly to al-Biruni,
*The Book of Instruction in the Elements of the Art of Astrology*, R. Ramsay
Wright translation (1934), section 476. Every result evaluates the explicit
modular arc ``ASC + A - B``. Fortune and Spirit reverse that arc by sect; Love
and Faith retain documented project formulas and are not attributed to that
section. The exact source/status matrix is
``docs/methodology/classical-astrology-helpers.md``.

Supported Parts:
- Part of Fortune (Pars Fortunae): Body, health, material success
- Part of Spirit (Pars Spiritus): Soul, intellect, spiritual matters
- Part of Love (Pars Amoris/Eros): Romance, desire, relationships
- Part of Faith (Pars Fidei): Belief, religion, trust

Provenance:
    Fortune/Spirit and their sect reversal map exactly to al-Biruni section 476;
    the public facsimile, formula table, and field mapping are recorded in
    ``docs/methodology/classical-astrology-helpers.md``. Love/Faith, the
    60-degree high-latitude switch, boundary handling, validation, and API
    shapes are labelled project conventions. Geometric sect delegates to the
    registered unrefracted horizontal-coordinate pipeline. No compatibility
    output supplies a formula or coefficient.
"""

from __future__ import annotations

from typing import Dict, Optional

# Latitude threshold above which 3D horizontal calculation is used
# for more accurate day/night determination. At latitudes beyond this,
# the ecliptic-horizon geometry becomes complex enough that 2D longitude
# comparison can give wrong results.
_EXTREME_LATITUDE_THRESHOLD = 60.0


def _is_sun_above_horizon_3d(
    jd: float,
    sun_lon: float,
    sun_lat: float,
    geo_lat: float,
    geo_lon: float,
) -> bool:
    """
    Determine if Sun is above horizon using 3D horizontal coordinates.

    Uses coordinate transformation from ecliptic to horizontal system
    to calculate Sun's true altitude above the local horizon.

    Args:
        jd: Julian Day in Universal Time
        sun_lon: Sun ecliptic longitude in degrees
        sun_lat: Sun ecliptic latitude in degrees (typically ~0)
        geo_lat: Observer geographic latitude in degrees (North positive)
        geo_lon: Observer geographic longitude in degrees (East positive)

    Returns:
        bool: True if Sun altitude >= 0 (above or on horizon)
    """
    from .utils import azalt, ECL2HOR

    # geopos: (longitude, latitude, altitude_meters)
    geopos = (geo_lon, geo_lat, 0.0)

    # No atmospheric refraction for geometric calculation
    # (refraction would slightly raise apparent position)
    atpress = 0.0
    attemp = 0.0

    # xin: (ecliptic_lon, ecliptic_lat, distance)
    # Distance doesn't affect altitude calculation, use 1.0 AU
    xin = (sun_lon, sun_lat, 1.0)

    # azalt returns (azimuth, true_altitude, apparent_altitude)
    _, true_altitude, _ = azalt(jd, ECL2HOR, geopos, atpress, attemp, xin)

    return true_altitude >= 0.0


def calc_arabic_part_of_fortune(
    asc: float, sun: float, moon: float, is_diurnal: bool = True
) -> float:
    """
    Calculate Part of Fortune (Pars Fortunae / Lot of Fortune).

    The most important Arabic Part, representing body, health, and material fortune.

    Formula (Vettius Valens, Al-Biruni):
        - Day chart (Sun above horizon): ASC + Moon - Sun
        - Night chart (Sun below horizon): ASC + Sun - Moon

    Args:
        asc: Ascendant (ecliptic) longitude in degrees
        sun: Sun ecliptic longitude in degrees
        moon: Moon ecliptic longitude in degrees
        is_diurnal: True if day chart (Sun above horizon), False if night

    Returns:
        float: Part of Fortune ecliptic longitude in degrees (0-360)

    Note:
        The formula reflects sect: in day charts, emphasize the Moon (nocturnal);
        in night charts, emphasize the Sun (diurnal). This balances opposites.

        Some modern astrologers use only the day formula regardless of sect.
        This implementation follows classical tradition.
    """
    if is_diurnal:
        lot = (asc + moon - sun) % 360.0
    else:
        lot = (asc + sun - moon) % 360.0

    return lot


def calc_arabic_part_of_spirit(
    asc: float, sun: float, moon: float, is_diurnal: bool = True
) -> float:
    """
    Calculate Part of Spirit (Pars Spiritus / Lot of Spirit / Daimon).

    Represents the soul, intellect, character, and spiritual development.

    Formula (opposite of Part of Fortune by sect):
        - Day chart: ASC + Sun - Moon
        - Night chart: ASC + Moon - Sun

    Args:
        asc: Ascendant longitude in degrees
        sun: Sun longitude in degrees
        moon: Moon longitude in degrees
        is_diurnal: True if day chart

    Returns:
        float: Part of Spirit longitude in degrees (0-360)

    Note:
        In Hellenistic astrology (Valens), Part of Spirit was called "Daimon"
        and considered complementary to Fortune. Together they represent
        body (Fortune) and soul (Spirit).
    """
    if is_diurnal:
        lot = (asc + sun - moon) % 360.0
    else:
        lot = (asc + moon - sun) % 360.0

    return lot


def calc_arabic_part_of_love(asc: float, venus: float, sun: float) -> float:
    """
    Calculate Part of Love (Pars Amoris / Lot of Eros).

    Represents romantic love, desire, sexual attraction, and relationships.

    Formula: ASC + Venus - Sun

    Args:
        asc: Ascendant longitude in degrees
        venus: Venus longitude in degrees
        sun: Sun longitude in degrees

    Returns:
        float: Part of Love longitude in degrees (0-360)

    Note:
        This sect-independent expression is an explicitly documented project
        convention retained by the API. It is not al-Biruni section 476's
        "Friendship and love" lot, which uses Fortune and Spirit and reverses
        them by sect.
    """
    return (asc + venus - sun) % 360.0


def calc_arabic_part_of_faith(asc: float, mercury: float, moon: float) -> float:
    """
    Calculate Part of Faith (Pars Fidei / Lot of Faith).

    Represents religious belief, trust, philosophical convictions.

    Formula: ASC + Mercury - Moon

    Args:
        asc: Ascendant longitude in degrees
        mercury: Mercury longitude in degrees
        moon: Moon longitude in degrees

    Returns:
        float: Part of Faith longitude in degrees (0-360)

    Note:
        This expression is an explicitly documented project convention. No
        page-level primary locator has been established for this exact
        name/formula pairing, and it is not presented as al-Biruni's Spirit or
        religion lot.
    """
    return (asc + mercury - moon) % 360.0


def is_day_chart(
    sun_lon: float,
    asc: float,
    *,
    jd: Optional[float] = None,
    geo_lat: Optional[float] = None,
    geo_lon: Optional[float] = None,
    sun_lat: float = 0.0,
) -> bool:
    """
    Determine if chart is diurnal (day) or nocturnal (night) based on sect.

    A chart is diurnal if the Sun is above the horizon (in houses 7-12).

    Args:
        sun_lon: Sun ecliptic longitude in degrees (0-360)
        asc: Ascendant ecliptic longitude in degrees (0-360)
        jd: Julian Day in UT (optional, enables 3D calculation)
        geo_lat: Observer geographic latitude in degrees (optional)
        geo_lon: Observer geographic longitude in degrees (optional)
        sun_lat: Sun ecliptic latitude in degrees (default 0.0)

    Returns:
        bool: True if day chart (Sun above horizon), False if night chart

    Algorithm:
        For moderate latitudes (|lat| <= 60°) or when location is not provided,
        uses the traditional 2D method based on ecliptic longitude comparison.

        For extreme latitudes (|lat| > 60°), uses 3D horizontal coordinate
        calculation to accurately determine Sun's altitude above the local
        horizon. This accounts for the complex geometry where the ecliptic
        and horizon planes intersect at steep angles.

    Note:
        At extreme latitudes (polar/subpolar regions), the 2D ecliptic method
        can give incorrect results because the horizon-ecliptic geometry
        differs significantly from the simplified model. The 3D method uses
        actual coordinate transformation to determine if Sun is above horizon.

        To enable 3D calculation, provide jd, geo_lat, and geo_lon parameters.
        If any of these are missing, the traditional 2D method is used.
    """
    # Use 3D calculation at extreme latitudes if location provided
    if (
        jd is not None
        and geo_lat is not None
        and geo_lon is not None
        and abs(geo_lat) > _EXTREME_LATITUDE_THRESHOLD
    ):
        return _is_sun_above_horizon_3d(jd, sun_lon, sun_lat, geo_lat, geo_lon)

    # Traditional 2D calculation using ecliptic longitude.
    # The lower hemisphere (below horizon, houses 1-6) is the open arc
    # (ASC, DSC) going counter-clockwise (increasing longitude through IC).
    # Boundaries: Sun on ASC (sunrise) or DSC (sunset) count as day chart.
    sun_rel = (sun_lon - asc) % 360.0
    return not (0.0 < sun_rel < 180.0)


def calc_all_arabic_parts(
    positions: Dict[str, float],
    *,
    jd: Optional[float] = None,
    geo_lat: Optional[float] = None,
    geo_lon: Optional[float] = None,
) -> Dict[str, float]:
    """
    Calculate all standard Arabic parts from a position dictionary.

    Args:
        positions: Dictionary of celestial positions in ecliptic longitude degrees.
                  Required keys: 'Asc', 'Sun', 'Moon', 'Mercury', 'Venus'
                  Optional key: 'Sun_lat' for Sun ecliptic latitude (default 0.0)
                  All values should be in range 0-360 degrees.
        jd: Julian Day in UT (optional, enables 3D day/night calculation)
        geo_lat: Observer geographic latitude in degrees (optional)
        geo_lon: Observer geographic longitude in degrees (optional)

    Returns:
        Dict[str, float]: Dictionary mapping part names to longitudes:
            - "Pars_Fortunae": Part of Fortune
            - "Pars_Spiritus": Part of Spirit
            - "Pars_Amoris": Part of Love
            - "Pars_Fidei": Part of Faith

    Example:
        >>> positions = {
        ...     'Asc': 15.5, 'Sun': 120.0, 'Moon': 240.0,
        ...     'Mercury': 130.0, 'Venus': 110.0
        ... }
        >>> parts = calc_all_arabic_parts(positions)
        >>> print(parts['Pars_Fortunae'])
        255.5

    Raises:
        KeyError: If any required position ('Asc', 'Sun', 'Moon', 'Mercury',
            'Venus') is missing. The parts are meaningless when computed from a
            placeholder longitude, so a missing input fails fast rather than
            silently producing a plausible-but-wrong result.

    Note:
        Only 'Sun_lat' is optional (defaults to 0.0); the five longitude keys
        above are required.

        For accurate day/night determination at extreme latitudes
        (``abs(lat) > 60°``),
        provide jd, geo_lat, and geo_lon parameters. This enables 3D horizontal
        coordinate calculation instead of the simplified 2D ecliptic method.
    """
    required = ("Asc", "Sun", "Moon", "Mercury", "Venus")
    missing = [k for k in required if k not in positions]
    if missing:
        raise KeyError(
            "calc_all_arabic_parts missing required position(s): " + ", ".join(missing)
        )

    asc = positions["Asc"]
    sun = positions["Sun"]
    sun_lat = positions.get("Sun_lat", 0.0)
    moon = positions["Moon"]
    mercury = positions["Mercury"]
    venus = positions["Venus"]

    is_diurnal = is_day_chart(
        sun, asc, jd=jd, geo_lat=geo_lat, geo_lon=geo_lon, sun_lat=sun_lat
    )

    return {
        "Pars_Fortunae": calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal),
        "Pars_Spiritus": calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal),
        "Pars_Amoris": calc_arabic_part_of_love(asc, venus, sun),
        "Pars_Fidei": calc_arabic_part_of_faith(asc, mercury, moon),
    }
