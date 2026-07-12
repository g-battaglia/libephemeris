#!/usr/bin/env python3
"""
LibEphemeris Cookbook
=====================

A collection of practical examples demonstrating common astrological calculations
using libephemeris. Each example is self-contained, well-commented, and ready to run.

Examples included:
1. Complete Natal Chart Calculation
2. Finding the Next Transit
3. Synastry Calculation (Chart Comparison)
4. Eclipse Search
5. Monthly Ephemeris Calculation

Requirements:
    pip install libephemeris

Usage:
    python cookbook.py

Author: LibEphemeris Contributors
License: Apache-2.0
"""

from __future__ import annotations

import libephemeris as ephem
from libephemeris.constants import (
    # Planet IDs
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    TRUE_NODE,
    CHIRON,
    # Calculation flags
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_SIDEREAL,
    # Sidereal modes
    SIDM_LAHIRI,
    ECL_TOTAL,
    ECL_ANNULAR,
    ECL_PARTIAL,
    ECL_PENUMBRAL,
)

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================


def format_longitude(lon: float) -> str:
    """
    Format a longitude in degrees to zodiacal notation (e.g., "15° 23' 45\" Aries").

    Args:
        lon: Longitude in degrees (0-360)

    Returns:
        Formatted string with sign, degrees, minutes, seconds
    """
    # Zodiac signs in order
    signs = [
        "Aries",
        "Taurus",
        "Gemini",
        "Cancer",
        "Leo",
        "Virgo",
        "Libra",
        "Scorpio",
        "Sagittarius",
        "Capricorn",
        "Aquarius",
        "Pisces",
    ]

    # Normalize longitude to 0-360
    lon = lon % 360

    # Determine sign (each sign is 30 degrees)
    sign_num = int(lon / 30)
    sign_name = signs[sign_num]

    # Calculate position within sign
    pos_in_sign = lon - (sign_num * 30)

    # Split into degrees, minutes, seconds
    degrees = int(pos_in_sign)
    minutes = int((pos_in_sign - degrees) * 60)
    seconds = int(((pos_in_sign - degrees) * 60 - minutes) * 60)

    return f"{degrees:02d}° {minutes:02d}' {seconds:02d}\" {sign_name}"


def format_aspect(angle: float) -> str | None:
    """
    Determine the aspect from an angular separation.

    Args:
        angle: Angular separation in degrees (-180 to 180)

    Returns:
        Aspect name and orb, or None if no major aspect
    """
    # Major aspects with their exact angles and typical orbs
    aspects = [
        (0, "Conjunction", 8),
        (60, "Sextile", 6),
        (90, "Square", 8),
        (120, "Trine", 8),
        (180, "Opposition", 8),
    ]

    # Normalize angle to positive
    angle = abs(angle)

    for exact_angle, name, orb in aspects:
        diff = abs(angle - exact_angle)
        if diff <= orb:
            applying = "Applying" if angle < exact_angle else "Separating"
            return f"{name} ({applying}, orb: {diff:.2f}°)"

    return None


# =============================================================================
# EXAMPLE 1: COMPLETE NATAL CHART CALCULATION
# =============================================================================


def calculate_natal_chart(
    year: int,
    month: int,
    day: int,
    hour: float,
    latitude: float,
    longitude: float,
    house_system: bytes = b"P",
    sidereal: bool = False,
) -> dict:
    """
    Calculate a complete natal chart including planets, houses, and angles.

    This example demonstrates:
    - Converting date/time to Julian Day
    - Calculating planetary positions with velocities
    - Calculating house cusps and angles (ASC, MC, etc.)
    - Optionally using sidereal zodiac

    Args:
        year: Birth year
        month: Birth month (1-12)
        day: Birth day (1-31)
        hour: Birth hour in decimal format (e.g., 14.5 for 2:30 PM)
        latitude: Birth place latitude (positive = North)
        longitude: Birth place longitude (positive = East)
        house_system: House system code:
            b"P" = Placidus (default)
            b"K" = Koch
            b"W" = Whole Sign
            b"E" = Equal
            b"R" = Regiomontanus
            b"C" = Campanus
        sidereal: If True, use Lahiri ayanamsha for sidereal calculations

    Returns:
        Dictionary containing:
        - planets: dict of planet positions with longitude, latitude, speed, sign
        - houses: list of 12 house cusp longitudes
        - angles: dict with ASC, MC, DESC, IC, Vertex
        - meta: calculation metadata

    Example:
        >>> # Calculate chart for Rome, Italy on Jan 1, 2000 at 12:00 noon
        >>> chart = calculate_natal_chart(
        ...     year=2000, month=1, day=1, hour=12.0,
        ...     latitude=41.9028, longitude=12.4964
        ... )
        >>> print(f"Sun: {chart['planets']['Sun']['formatted']}")
        Sun: 10° 37' 26" Capricorn
    """
    # Step 1: Convert date/time to Julian Day
    # The Julian Day is the standard time format for astronomical calculations
    jd = ephem.julday(year, month, day, hour)

    # Step 2: Set up sidereal mode if requested
    calc_flags = FLG_SWIEPH | FLG_SPEED
    if sidereal:
        ephem.set_sid_mode(SIDM_LAHIRI)
        calc_flags |= FLG_SIDEREAL
        ayanamsha = ephem.get_ayanamsa_ut(jd)
    else:
        # Reset to tropical (Fagan-Bradley is default, but we use 0 offset)
        ayanamsha = 0.0

    # Step 3: Calculate planetary positions
    # List of planets to calculate
    planets = [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MERCURY, "Mercury"),
        (VENUS, "Venus"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
        (URANUS, "Uranus"),
        (NEPTUNE, "Neptune"),
        (PLUTO, "Pluto"),
        (TRUE_NODE, "North Node"),
        (CHIRON, "Chiron"),
    ]

    planet_positions = {}
    for planet_id, planet_name in planets:
        # calc_ut returns (position_tuple, flags_used)
        # position_tuple = (longitude, latitude, distance, lon_speed, lat_speed, dist_speed)
        pos, _ = ephem.calc_ut(jd, planet_id, calc_flags)

        # Determine if planet is retrograde (negative longitude speed)
        is_retrograde = pos[3] < 0

        planet_positions[planet_name] = {
            "longitude": pos[0],  # Ecliptic longitude in degrees
            "latitude": pos[1],  # Ecliptic latitude in degrees
            "distance": pos[2],  # Distance in AU
            "speed": pos[3],  # Daily motion in degrees/day
            "retrograde": is_retrograde,
            "formatted": format_longitude(pos[0]),
            "retrograde_symbol": "R" if is_retrograde else "",
        }

    # Step 4: Calculate house cusps and angles
    # houses returns (cusps, ascmc)
    # cusps = tuple of 12 house cusp longitudes (index 0-11 for houses 1-12)
    # ascmc = tuple of important points: [0]=ASC, [1]=MC, [2]=ARMC, [3]=Vertex, etc.
    cusps, ascmc = ephem.houses(jd, latitude, longitude, house_system)

    # Note: Some house systems may have issues at extreme latitudes
    # (> 66.5° North/South). In such cases, consider using Whole Sign houses.

    # Step 5: Extract angles
    angles = {
        "ASC": {"longitude": ascmc[0], "formatted": format_longitude(ascmc[0])},
        "MC": {"longitude": ascmc[1], "formatted": format_longitude(ascmc[1])},
        "DESC": {
            "longitude": (ascmc[0] + 180) % 360,
            "formatted": format_longitude((ascmc[0] + 180) % 360),
        },
        "IC": {
            "longitude": (ascmc[1] + 180) % 360,
            "formatted": format_longitude((ascmc[1] + 180) % 360),
        },
        "Vertex": {"longitude": ascmc[3], "formatted": format_longitude(ascmc[3])},
    }

    # Step 6: Format house cusps
    houses = []
    for i in range(12):
        houses.append(
            {
                "house": i + 1,
                "longitude": cusps[i],
                "formatted": format_longitude(cusps[i]),
            }
        )

    # Return complete chart data
    return {
        "planets": planet_positions,
        "houses": houses,
        "angles": angles,
        "meta": {
            "julian_day": jd,
            "date": f"{year}-{month:02d}-{day:02d}",
            "time": f"{int(hour):02d}:{int((hour % 1) * 60):02d}",
            "location": f"{latitude:.4f}, {longitude:.4f}",
            "house_system": house_system.decode(),
            "zodiac": "Sidereal (Lahiri)" if sidereal else "Tropical",
            "ayanamsha": ayanamsha if sidereal else None,
        },
    }


def print_natal_chart(chart: dict) -> None:
    """Pretty print a natal chart."""
    print("=" * 60)
    print("NATAL CHART")
    print("=" * 60)
    print(f"Date: {chart['meta']['date']} at {chart['meta']['time']}")
    print(f"Location: {chart['meta']['location']}")
    print(f"House System: {chart['meta']['house_system']}")
    print(f"Zodiac: {chart['meta']['zodiac']}")
    print()

    print("PLANETS:")
    print("-" * 40)
    for name, data in chart["planets"].items():
        retro = " (R)" if data["retrograde"] else ""
        print(f"  {name:12} {data['formatted']}{retro}")

    print()
    print("ANGLES:")
    print("-" * 40)
    for name, data in chart["angles"].items():
        print(f"  {name:12} {data['formatted']}")

    print()
    print("HOUSE CUSPS:")
    print("-" * 40)
    for house in chart["houses"]:
        print(f"  House {house['house']:2}:    {house['formatted']}")


# =============================================================================
# EXAMPLE 2: FINDING THE NEXT TRANSIT
# =============================================================================


def find_next_transit(
    planet_id: int,
    target_longitude: float,
    start_year: int,
    start_month: int,
    start_day: int,
    start_hour: float = 0.0,
) -> dict:
    """
    Find when a planet next crosses a specific longitude.

    This is useful for finding:
    - Sign ingresses (when a planet enters a new sign)
    - Planetary returns (when a planet returns to its natal position)
    - Specific degree transits

    This example demonstrates:
    - Using cross_ut for generic planet crossings
    - Using solcross_ut for Sun crossings (optimized)
    - Using mooncross_ut for Moon crossings (optimized)
    - Converting Julian Day back to calendar date

    Args:
        planet_id: Planet ID constant (e.g., SUN, MOON, MARS)
        target_longitude: Target longitude in degrees (0-360)
        start_year: Year to start searching from
        start_month: Month to start searching from
        start_day: Day to start searching from
        start_hour: Hour to start searching from

    Returns:
        Dictionary containing:
        - crossing_jd: Julian Day of the crossing
        - date: Calendar date/time of crossing
        - planet_position: Position at crossing time
        - target: Target longitude searched for

    Example:
        >>> # Find next Sun ingress into Aries (0° Aries = vernal equinox)
        >>> result = find_next_transit(SUN, 0.0, 2024, 1, 1)
        >>> print(f"Vernal Equinox: {result['date']}")
        Vernal Equinox: 2024-03-20 03:06

        >>> # Find next Mars ingress into Leo (120° = Leo)
        >>> result = find_next_transit(MARS, 120.0, 2024, 1, 1)
    """
    # Convert start date to Julian Day
    start_jd = ephem.julday(start_year, start_month, start_day, start_hour)

    # Normalize target longitude to 0-360
    target_longitude = target_longitude % 360

    # Use optimized function for Sun and Moon, generic function for others
    if planet_id == SUN:
        # solcross_ut is optimized for Sun crossings
        crossing_jd = ephem.solcross_ut(target_longitude, start_jd, FLG_SWIEPH)
    elif planet_id == MOON:
        # mooncross_ut is optimized for Moon crossings
        crossing_jd = ephem.mooncross_ut(target_longitude, start_jd, FLG_SWIEPH)
    else:
        # cross_ut works for any planet
        crossing_jd = ephem.cross_ut(planet_id, target_longitude, start_jd, FLG_SWIEPH)

    # Convert Julian Day back to calendar date
    year, month, day, hour = ephem.revjul(crossing_jd)

    # Format time
    hour_int = int(hour)
    minute = int((hour - hour_int) * 60)

    # Get planet position at crossing to verify
    pos, _ = ephem.calc_ut(crossing_jd, planet_id, FLG_SWIEPH | FLG_SPEED)

    # Planet names for display
    planet_names = {
        SUN: "Sun",
        MOON: "Moon",
        MERCURY: "Mercury",
        VENUS: "Venus",
        MARS: "Mars",
        JUPITER: "Jupiter",
        SATURN: "Saturn",
        URANUS: "Uranus",
        NEPTUNE: "Neptune",
        PLUTO: "Pluto",
    }
    planet_name = planet_names.get(planet_id, f"Planet {planet_id}")

    return {
        "planet": planet_name,
        "target": target_longitude,
        "target_formatted": format_longitude(target_longitude),
        "crossing_jd": crossing_jd,
        "date": f"{year}-{month:02d}-{day:02d} {hour_int:02d}:{minute:02d}",
        "planet_position": {
            "longitude": pos[0],
            "formatted": format_longitude(pos[0]),
            "speed": pos[3],
        },
    }


def find_planetary_return(
    planet_id: int,
    natal_longitude: float,
    search_year: int,
) -> dict:
    """
    Find the next planetary return (when planet returns to natal position).

    Useful for:
    - Solar returns (birthday charts)
    - Lunar returns (monthly cycle)
    - Saturn returns, Jupiter returns, etc.

    Args:
        planet_id: Planet ID constant
        natal_longitude: Natal position of the planet in degrees
        search_year: Year to search in

    Returns:
        Transit information dictionary

    Example:
        >>> # Find 2024 Solar Return for someone born with Sun at 15° Leo
        >>> result = find_planetary_return(SUN, 135.0, 2024)
        >>> print(f"Solar Return: {result['date']}")
    """
    return find_next_transit(planet_id, natal_longitude, search_year, 1, 1, 0.0)


def find_sign_ingresses(planet_id: int, year: int) -> list[dict]:
    """
    Find all sign ingresses for a planet in a given year.

    Args:
        planet_id: Planet ID constant
        year: Year to search

    Returns:
        List of transit dictionaries for each ingress

    Example:
        >>> # Find all Sun ingresses (zodiac sign changes) in 2024
        >>> ingresses = find_sign_ingresses(SUN, 2024)
        >>> for ing in ingresses:
        ...     print(f"{ing['target_formatted']}: {ing['date']}")
    """
    ingresses = []
    start_jd = ephem.julday(year, 1, 1, 0.0)
    end_jd = ephem.julday(year + 1, 1, 1, 0.0)

    # Sign boundaries (0° of each sign)
    for sign_num in range(12):
        target = sign_num * 30  # 0, 30, 60, ... 330

        # Use the optimized crossing function
        if planet_id == SUN:
            crossing_jd = ephem.solcross_ut(target, start_jd, FLG_SWIEPH)
        elif planet_id == MOON:
            crossing_jd = ephem.mooncross_ut(target, start_jd, FLG_SWIEPH)
        else:
            crossing_jd = ephem.cross_ut(planet_id, target, start_jd, FLG_SWIEPH)

        # Check if crossing is within the year
        if start_jd <= crossing_jd < end_jd:
            year_c, month, day, hour = ephem.revjul(crossing_jd)
            hour_int = int(hour)
            minute = int((hour - hour_int) * 60)

            ingresses.append(
                {
                    "sign_num": sign_num,
                    "target": target,
                    "target_formatted": format_longitude(target),
                    "crossing_jd": crossing_jd,
                    "date": f"{year_c}-{month:02d}-{day:02d} {hour_int:02d}:{minute:02d}",
                }
            )

    # Sort by crossing time
    ingresses.sort(key=lambda x: x["crossing_jd"])
    return ingresses


# =============================================================================
# EXAMPLE 3: SYNASTRY CALCULATION (CHART COMPARISON)
# =============================================================================


def calculate_synastry(chart1: dict, chart2: dict) -> dict:
    """
    Calculate synastry aspects between two natal charts.

    Synastry compares the planetary positions of two people to analyze
    relationship compatibility and dynamics.

    This example demonstrates:
    - Comparing planetary positions between charts
    - Calculating angular separations (aspects)
    - Using difdeg2n for normalized angle differences

    Args:
        chart1: First person's natal chart (from calculate_natal_chart)
        chart2: Second person's natal chart (from calculate_natal_chart)

    Returns:
        Dictionary containing:
        - aspects: List of aspects between planets
        - summary: Count of harmonious vs challenging aspects

    Example:
        >>> chart_person1 = calculate_natal_chart(1990, 5, 15, 10.0, 40.7, -74.0)
        >>> chart_person2 = calculate_natal_chart(1992, 8, 22, 14.0, 34.0, -118.2)
        >>> synastry = calculate_synastry(chart_person1, chart_person2)
        >>> for aspect in synastry['aspects'][:5]:
        ...     print(f"{aspect['planet1']} - {aspect['planet2']}: {aspect['aspect']}")
    """
    aspects_found = []
    harmonious_count = 0
    challenging_count = 0

    # Planets to compare
    planets_to_check = [
        "Sun",
        "Moon",
        "Mercury",
        "Venus",
        "Mars",
        "Jupiter",
        "Saturn",
        "Uranus",
        "Neptune",
        "Pluto",
    ]

    # Compare each planet in chart1 with each planet in chart2
    for p1_name in planets_to_check:
        if p1_name not in chart1["planets"]:
            continue
        p1_lon = chart1["planets"][p1_name]["longitude"]

        for p2_name in planets_to_check:
            if p2_name not in chart2["planets"]:
                continue
            p2_lon = chart2["planets"][p2_name]["longitude"]

            # Calculate the angular separation
            # difdeg2n normalizes to -180 to +180 range
            separation = ephem.difdeg2n(p1_lon, p2_lon)

            # Check if there's an aspect
            aspect = format_aspect(separation)
            if aspect:
                aspect_entry = {
                    "planet1": f"Person 1's {p1_name}",
                    "planet1_position": format_longitude(p1_lon),
                    "planet2": f"Person 2's {p2_name}",
                    "planet2_position": format_longitude(p2_lon),
                    "separation": abs(separation),
                    "aspect": aspect,
                }
                aspects_found.append(aspect_entry)

                # Count harmonious vs challenging
                if "Trine" in aspect or "Sextile" in aspect:
                    harmonious_count += 1
                elif "Square" in aspect or "Opposition" in aspect:
                    challenging_count += 1
                # Conjunctions can be either depending on planets

    # Also compare angles
    for angle1_name in ["ASC", "MC"]:
        if angle1_name not in chart1["angles"]:
            continue
        angle1_lon = chart1["angles"][angle1_name]["longitude"]

        for p2_name in planets_to_check:
            if p2_name not in chart2["planets"]:
                continue
            p2_lon = chart2["planets"][p2_name]["longitude"]

            separation = ephem.difdeg2n(angle1_lon, p2_lon)
            aspect = format_aspect(separation)

            if aspect:
                aspects_found.append(
                    {
                        "planet1": f"Person 1's {angle1_name}",
                        "planet1_position": format_longitude(angle1_lon),
                        "planet2": f"Person 2's {p2_name}",
                        "planet2_position": format_longitude(p2_lon),
                        "separation": abs(separation),
                        "aspect": aspect,
                    }
                )

    return {
        "aspects": aspects_found,
        "summary": {
            "total_aspects": len(aspects_found),
            "harmonious": harmonious_count,
            "challenging": challenging_count,
        },
    }


def calculate_composite_midpoints(chart1: dict, chart2: dict) -> dict:
    """
    Calculate composite chart midpoints between two natal charts.

    The composite chart uses the midpoint between each pair of planets
    to create a chart representing the relationship itself.

    Args:
        chart1: First person's natal chart
        chart2: Second person's natal chart

    Returns:
        Dictionary with composite planetary positions

    Example:
        >>> composite = calculate_composite_midpoints(chart1, chart2)
        >>> print(f"Composite Sun: {composite['planets']['Sun']['formatted']}")
    """
    composite = {"planets": {}, "angles": {}}

    # Calculate midpoints for planets
    for planet_name in chart1["planets"]:
        if planet_name not in chart2["planets"]:
            continue

        lon1 = chart1["planets"][planet_name]["longitude"]
        lon2 = chart2["planets"][planet_name]["longitude"]

        # Use deg_midp for proper midpoint calculation
        # This handles the 0/360 boundary correctly
        midpoint = ephem.deg_midp(lon1, lon2)

        composite["planets"][planet_name] = {
            "longitude": midpoint,
            "formatted": format_longitude(midpoint),
        }

    # Calculate midpoints for angles
    for angle_name in ["ASC", "MC"]:
        if angle_name in chart1["angles"] and angle_name in chart2["angles"]:
            lon1 = chart1["angles"][angle_name]["longitude"]
            lon2 = chart2["angles"][angle_name]["longitude"]
            midpoint = ephem.deg_midp(lon1, lon2)

            composite["angles"][angle_name] = {
                "longitude": midpoint,
                "formatted": format_longitude(midpoint),
            }

    return composite


# =============================================================================
# EXAMPLE 4: ECLIPSE SEARCH
# =============================================================================


def find_next_solar_eclipse(
    start_year: int,
    start_month: int,
    start_day: int,
    eclipse_type: int | None = None,
) -> dict:
    """
    Find the next solar eclipse after a given date.

    This example demonstrates:
    - Using sol_eclipse_when_glob to find global eclipses
    - Interpreting eclipse type flags
    - Getting eclipse timing details

    Args:
        start_year: Year to start searching from
        start_month: Month to start searching from
        start_day: Day to start searching from
        eclipse_type: Optional filter for eclipse type:
            - ECL_TOTAL: Total eclipses only
            - ECL_ANNULAR: Annular eclipses only
            - ECL_PARTIAL: Partial eclipses only
            - None: Any eclipse type

    Returns:
        Dictionary containing:
        - type: Eclipse type description
        - maximum: Date/time of maximum eclipse
        - times: All eclipse phase times
        - duration: Duration of totality/annularity (if applicable)

    Example:
        >>> # Find next total solar eclipse
        >>> eclipse = find_next_solar_eclipse(2024, 1, 1, ECL_TOTAL)
        >>> print(f"Next total solar eclipse: {eclipse['maximum']}")
        Next total solar eclipse: 2024-04-08 18:17
    """
    # Convert start date to Julian Day
    start_jd = ephem.julday(start_year, start_month, start_day, 0.0)

    # Search for eclipse
    # Returns (times_tuple, eclipse_type_flags)
    # times_tuple indices:
    #   [0] = time of maximum eclipse
    #   [1] = time of first contact (partial eclipse begins)
    #   [2] = time of second contact (total/annular begins)
    #   [3] = time of third contact (total/annular ends)
    #   [4] = time of fourth contact (partial eclipse ends)
    #   [5] = time of sunrise during eclipse (if applicable)
    #   [6] = time of sunset during eclipse (if applicable)
    #   [7] = unused
    if eclipse_type:
        ecl_flags, times = ephem.sol_eclipse_when_glob(
            start_jd, flags=FLG_SWIEPH, ecltype=eclipse_type
        )
    else:
        ecl_flags, times = ephem.sol_eclipse_when_glob(start_jd, flags=FLG_SWIEPH)

    # Interpret eclipse type
    if ecl_flags & ECL_TOTAL:
        type_str = "Total Solar Eclipse"
    elif ecl_flags & ECL_ANNULAR:
        type_str = "Annular Solar Eclipse"
    elif ecl_flags & ECL_PARTIAL:
        type_str = "Partial Solar Eclipse"
    else:
        type_str = "Solar Eclipse (type unknown)"

    # Convert maximum time to date
    year, month, day, hour = ephem.revjul(times[0])
    hour_int = int(hour)
    minute = int((hour - hour_int) * 60)

    # Calculate duration of totality/annularity if applicable
    duration_str = None
    if times[2] > 0 and times[3] > 0:  # Second and third contact times exist
        duration_days = times[3] - times[2]
        duration_minutes = duration_days * 24 * 60
        duration_seconds = duration_minutes * 60
        duration_str = f"{int(duration_minutes)}m {int(duration_seconds % 60)}s"

    return {
        "type": type_str,
        "maximum": f"{year}-{month:02d}-{day:02d} {hour_int:02d}:{minute:02d}",
        "maximum_jd": times[0],
        "times": {
            "maximum": times[0],
            "first_contact": times[1] if times[1] > 0 else None,
            "second_contact": times[2] if times[2] > 0 else None,
            "third_contact": times[3] if times[3] > 0 else None,
            "fourth_contact": times[4] if times[4] > 0 else None,
        },
        "duration": duration_str,
        "eclipse_flags": ecl_flags,
    }


def find_next_lunar_eclipse(
    start_year: int,
    start_month: int,
    start_day: int,
) -> dict:
    """
    Find the next lunar eclipse after a given date.

    Args:
        start_year: Year to start searching from
        start_month: Month to start searching from
        start_day: Day to start searching from

    Returns:
        Dictionary with eclipse information

    Example:
        >>> eclipse = find_next_lunar_eclipse(2024, 1, 1)
        >>> print(f"Next lunar eclipse: {eclipse['maximum']}")
    """
    start_jd = ephem.julday(start_year, start_month, start_day, 0.0)

    # lun_eclipse_when returns (eclipse_type_flags, times)
    ecl_flags, times = ephem.lun_eclipse_when(start_jd)

    # Interpret eclipse type
    if ecl_flags & ECL_TOTAL:
        type_str = "Total Lunar Eclipse"
    elif ecl_flags & ECL_PARTIAL:
        type_str = "Partial Lunar Eclipse"
    elif ecl_flags & ECL_PENUMBRAL:
        type_str = "Penumbral Lunar Eclipse"
    else:
        type_str = "Lunar Eclipse (type unknown)"

    # Convert maximum time to date
    year, month, day, hour = ephem.revjul(times[0])
    hour_int = int(hour)
    minute = int((hour - hour_int) * 60)

    return {
        "type": type_str,
        "maximum": f"{year}-{month:02d}-{day:02d} {hour_int:02d}:{minute:02d}",
        "maximum_jd": times[0],
        "eclipse_flags": ecl_flags,
    }


def find_eclipse_visibility(
    eclipse_jd: float,
    latitude: float,
    longitude: float,
    altitude: float = 0.0,
) -> dict:
    """
    Check if a solar eclipse is visible from a specific location.

    Args:
        eclipse_jd: Julian Day of eclipse maximum
        latitude: Observer latitude
        longitude: Observer longitude
        altitude: Observer altitude in meters

    Returns:
        Dictionary with visibility information including magnitude and obscuration

    Example:
        >>> eclipse = find_next_solar_eclipse(2024, 1, 1, ECL_TOTAL)
        >>> visibility = find_eclipse_visibility(
        ...     eclipse['maximum_jd'], 32.7767, -96.7970  # Dallas
        ... )
        >>> print(f"Magnitude: {visibility['magnitude']:.3f}")
    """
    # Use sol_eclipse_how to get eclipse circumstances at location
    ecl_flags, attr = ephem.sol_eclipse_how(eclipse_jd, (longitude, latitude, altitude))

    # attr indices:
    #   [0] = eclipse magnitude
    #   [1] = ratio of diameters
    #   [2] = fraction of solar disc covered (obscuration)
    #   [3] = azimuth of sun at maximum
    #   [4] = altitude of sun at maximum
    #   [5] = apparent diameter of moon
    #   [6] = apparent diameter of sun
    #   [7] = difference of diameters
    #   [8] = eclipse magnitude for umbral/penumbral phases
    #   [9] = saros series number
    #   [10] = inex series number

    from libephemeris.constants import ECL_VISIBLE

    is_visible = bool(ecl_flags & ECL_VISIBLE)

    return {
        "visible": is_visible,
        "magnitude": attr[0] if is_visible else 0.0,
        "obscuration": attr[2] if is_visible else 0.0,
        "sun_altitude": attr[4],
        "sun_azimuth": attr[3],
        "eclipse_type_local": (
            "Total"
            if ecl_flags & ECL_TOTAL
            else "Annular"
            if ecl_flags & ECL_ANNULAR
            else "Partial"
        ),
    }


# =============================================================================
# EXAMPLE 5: MONTHLY EPHEMERIS CALCULATION
# =============================================================================


def calculate_monthly_ephemeris(
    year: int,
    month: int,
    include_asteroids: bool = False,
) -> list[dict]:
    """
    Calculate planetary positions for each day of a month.

    This is useful for:
    - Creating ephemeris tables
    - Tracking planetary movements
    - Finding retrograde periods

    This example demonstrates:
    - Iterating through dates
    - Calculating multiple planet positions
    - Tracking retrograde status

    Args:
        year: Year
        month: Month (1-12)
        include_asteroids: If True, include Chiron

    Returns:
        List of dictionaries, one per day, containing all planetary positions

    Example:
        >>> ephemeris = calculate_monthly_ephemeris(2024, 3)
        >>> for day_data in ephemeris[:3]:
        ...     print(f"{day_data['date']}: Sun at {day_data['Sun']['formatted']}")
    """
    # Determine number of days in month
    if month == 12:
        next_month_jd = ephem.julday(year + 1, 1, 1, 0.0)
    else:
        next_month_jd = ephem.julday(year, month + 1, 1, 0.0)
    first_day_jd = ephem.julday(year, month, 1, 0.0)
    days_in_month = int(next_month_jd - first_day_jd)

    # Planets to calculate
    planets = [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MERCURY, "Mercury"),
        (VENUS, "Venus"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
        (URANUS, "Uranus"),
        (NEPTUNE, "Neptune"),
        (PLUTO, "Pluto"),
    ]

    if include_asteroids:
        planets.append((CHIRON, "Chiron"))

    # Calculate for each day at noon (12:00)
    ephemeris = []

    for day in range(1, days_in_month + 1):
        jd = ephem.julday(year, month, day, 12.0)  # Noon UTC

        day_data = {
            "date": f"{year}-{month:02d}-{day:02d}",
            "julian_day": jd,
        }

        for planet_id, planet_name in planets:
            pos, _ = ephem.calc_ut(jd, planet_id, FLG_SWIEPH | FLG_SPEED)

            day_data[planet_name] = {
                "longitude": pos[0],
                "formatted": format_longitude(pos[0]),
                "latitude": pos[1],
                "distance": pos[2],
                "speed": pos[3],
                "retrograde": pos[3] < 0,
            }

        ephemeris.append(day_data)

    return ephemeris


def find_retrograde_periods(year: int, planet_id: int) -> list[dict]:
    """
    Find all retrograde periods for a planet in a given year.

    Args:
        year: Year to search
        planet_id: Planet ID constant (outer planets have more pronounced retrogrades)

    Returns:
        List of retrograde periods with start and end dates

    Example:
        >>> # Find Mercury retrograde periods in 2024
        >>> retrogrades = find_retrograde_periods(2024, MERCURY)
        >>> for period in retrogrades:
        ...     print(f"{period['start']} to {period['end']}")
    """
    # Calculate positions for each day of the year
    start_jd = ephem.julday(year, 1, 1, 12.0)
    end_jd = ephem.julday(year + 1, 1, 1, 12.0)

    periods = []
    in_retrograde = False
    retrograde_start = None

    jd = start_jd
    while jd < end_jd:
        pos, _ = ephem.calc_ut(jd, planet_id, FLG_SWIEPH | FLG_SPEED)
        is_retrograde = pos[3] < 0

        if is_retrograde and not in_retrograde:
            # Retrograde begins
            in_retrograde = True
            retrograde_start = jd
        elif not is_retrograde and in_retrograde:
            # Retrograde ends
            in_retrograde = False

            # Get start date
            y1, m1, d1, _ = ephem.revjul(retrograde_start)
            # Get end date
            y2, m2, d2, _ = ephem.revjul(jd)

            # Get station positions
            pos_start, _ = ephem.calc_ut(
                retrograde_start, planet_id, FLG_SWIEPH | FLG_SPEED
            )
            pos_end, _ = ephem.calc_ut(jd, planet_id, FLG_SWIEPH | FLG_SPEED)

            periods.append(
                {
                    "start": f"{y1}-{m1:02d}-{d1:02d}",
                    "end": f"{y2}-{m2:02d}-{d2:02d}",
                    "start_position": format_longitude(pos_start[0]),
                    "end_position": format_longitude(pos_end[0]),
                    "duration_days": int(jd - retrograde_start),
                }
            )

        jd += 1  # Check daily

    return periods


def print_monthly_ephemeris(ephemeris: list[dict]) -> None:
    """Pretty print a monthly ephemeris table."""
    if not ephemeris:
        return

    # Get planet names (excluding metadata fields)
    planet_names = [k for k in ephemeris[0].keys() if k not in ("date", "julian_day")]

    # Print header
    header = f"{'Date':<12}"
    for name in planet_names[:6]:  # Limit to 6 planets for readability
        header += f"{name:<12}"
    print(header)
    print("-" * len(header))

    # Print data
    for day in ephemeris:
        line = f"{day['date']:<12}"
        for name in planet_names[:6]:
            lon = day[name]["longitude"]
            retro = "R" if day[name]["retrograde"] else " "
            # Show just degrees in sign
            sign_num = int(lon / 30)
            deg_in_sign = int(lon % 30)
            signs_abbr = [
                "Ar",
                "Ta",
                "Ge",
                "Cn",
                "Le",
                "Vi",
                "Li",
                "Sc",
                "Sg",
                "Cp",
                "Aq",
                "Pi",
            ]
            line += f"{deg_in_sign:2d}{signs_abbr[sign_num]}{retro:<7}"
        print(line)


# =============================================================================
# MAIN: RUN ALL EXAMPLES
# =============================================================================


def main():
    """Run all cookbook examples."""
    print("=" * 70)
    print("LIBEPHEMERIS COOKBOOK - PRACTICAL EXAMPLES")
    print("=" * 70)
    print()

    # -------------------------------------------------------------------------
    # Example 1: Natal Chart
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("EXAMPLE 1: NATAL CHART CALCULATION")
    print("=" * 70)
    print()
    print("Calculating natal chart for Rome, Italy")
    print("Date: January 1, 2000 at 12:00 noon")
    print()

    chart = calculate_natal_chart(
        year=2000,
        month=1,
        day=1,
        hour=12.0,
        latitude=41.9028,
        longitude=12.4964,
        house_system=b"P",  # Placidus
    )
    print_natal_chart(chart)

    # -------------------------------------------------------------------------
    # Example 2: Finding Transits
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("EXAMPLE 2: FINDING TRANSITS")
    print("=" * 70)
    print()

    # Find next vernal equinox
    print("Finding next vernal equinox (Sun at 0 Aries)...")
    transit = find_next_transit(SUN, 0.0, 2024, 1, 1)
    print(f"  {transit['planet']} enters Aries: {transit['date']}")
    print()

    # Find Sun sign ingresses for 2024
    print("Sun sign ingresses for 2024:")
    ingresses = find_sign_ingresses(SUN, 2024)
    for ing in ingresses:
        print(f"  {ing['target_formatted']}: {ing['date']}")

    # -------------------------------------------------------------------------
    # Example 3: Synastry
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("EXAMPLE 3: SYNASTRY (CHART COMPARISON)")
    print("=" * 70)
    print()

    # Create two sample charts
    print("Comparing two charts...")
    print("  Person 1: May 15, 1990, 10:00 AM, New York")
    print("  Person 2: August 22, 1992, 2:00 PM, Los Angeles")
    print()

    chart1 = calculate_natal_chart(1990, 5, 15, 10.0, 40.7128, -74.0060)
    chart2 = calculate_natal_chart(1992, 8, 22, 14.0, 34.0522, -118.2437)

    synastry = calculate_synastry(chart1, chart2)

    print(f"Total aspects found: {synastry['summary']['total_aspects']}")
    print(f"  Harmonious (trines, sextiles): {synastry['summary']['harmonious']}")
    print(f"  Challenging (squares, oppositions): {synastry['summary']['challenging']}")
    print()
    print("Sample aspects:")
    for aspect in synastry["aspects"][:5]:
        print(f"  {aspect['planet1']} - {aspect['planet2']}")
        print(f"    {aspect['aspect']}")

    # -------------------------------------------------------------------------
    # Example 4: Eclipse Search
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("EXAMPLE 4: ECLIPSE SEARCH")
    print("=" * 70)
    print()

    # Find next total solar eclipse
    print("Searching for next total solar eclipse from 2024...")
    eclipse = find_next_solar_eclipse(2024, 1, 1, ECL_TOTAL)
    print(f"  Type: {eclipse['type']}")
    print(f"  Maximum: {eclipse['maximum']} UTC")
    if eclipse["duration"]:
        print(f"  Duration of totality: {eclipse['duration']}")
    print()

    # Check visibility from Dallas (in path of totality for April 2024)
    print("Checking visibility from Dallas, TX...")
    visibility = find_eclipse_visibility(eclipse["maximum_jd"], 32.7767, -96.7970)
    if visibility["visible"]:
        print("  Visible: Yes")
        print(f"  Magnitude: {visibility['magnitude']:.3f}")
        print(f"  Obscuration: {visibility['obscuration'] * 100:.1f}%")
        print(f"  Local type: {visibility['eclipse_type_local']}")
    else:
        print("  Visible: No")
    print()

    # Find next lunar eclipse
    print("Searching for next lunar eclipse from 2024...")
    lunar_eclipse = find_next_lunar_eclipse(2024, 1, 1)
    print(f"  Type: {lunar_eclipse['type']}")
    print(f"  Maximum: {lunar_eclipse['maximum']} UTC")

    # -------------------------------------------------------------------------
    # Example 5: Monthly Ephemeris
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("EXAMPLE 5: MONTHLY EPHEMERIS")
    print("=" * 70)
    print()

    print("Ephemeris for March 2024 (first 7 days):")
    print()
    ephemeris = calculate_monthly_ephemeris(2024, 3)
    print_monthly_ephemeris(ephemeris[:7])
    print()

    # Find Mercury retrograde periods
    print("Mercury retrograde periods in 2024:")
    retrogrades = find_retrograde_periods(2024, MERCURY)
    for period in retrogrades:
        print(
            f"  {period['start']} to {period['end']} ({period['duration_days']} days)"
        )
        print(f"    From {period['start_position']} to {period['end_position']}")

    print("\n" + "=" * 70)
    print("COOKBOOK EXAMPLES COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
