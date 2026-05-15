#!/usr/bin/env python3
"""
Basic Usage Example for libephemeris.

This script demonstrates the core functionality of libephemeris:
1. Converting dates to Julian Day
2. Calculating planetary positions
3. Calculating house cusps and angles
4. Getting ayanamsha values

This is a good starting point for new users.

Requirements:
    pip install libephemeris

Usage:
    python examples/basic_usage.py
"""

from __future__ import annotations

import libephemeris as eph
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
)


def example_julian_day() -> None:
    """Example 1: Julian Day conversion."""
    print("=" * 60)
    print("Example 1: Julian Day Conversion")
    print("=" * 60)

    # Convert a date to Julian Day
    # Format: year, month, day, hour (decimal)
    jd = eph.julday(2000, 1, 1, 12.0)  # J2000.0 epoch
    print(f"January 1, 2000 at noon = JD {jd}")

    # Convert back to calendar date
    year, month, day, hour = eph.revjul(jd)
    print(
        f"JD {jd} = {year}-{month:02d}-{day:02d} {int(hour):02d}:{int((hour % 1) * 60):02d}"
    )

    # Current moment (example with specific date)
    jd_now = eph.julday(2024, 6, 21, 12.0)  # Summer solstice 2024
    print(f"June 21, 2024 at noon = JD {jd_now}")

    # Get Delta T (difference between TT and UT)
    delta_t = eph.deltat(jd_now)
    print(f"Delta T at JD {jd_now}: {delta_t * 86400:.2f} seconds")


def example_planetary_positions() -> None:
    """Example 2: Calculate planetary positions."""
    print("\n" + "=" * 60)
    print("Example 2: Planetary Positions")
    print("=" * 60)

    # Julian Day for calculation
    jd = eph.julday(2024, 1, 1, 12.0)
    print(f"\nPlanetary positions at JD {jd} (January 1, 2024, noon UT):\n")

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
        (TRUE_NODE, "True Node"),
        (CHIRON, "Chiron"),
    ]

    # Calculate with speed
    flags = FLG_SWIEPH | FLG_SPEED

    print(
        f"{'Planet':<12} {'Longitude':>12} {'Latitude':>10} {'Distance':>10} {'Speed':>10}"
    )
    print("-" * 56)

    for planet_id, name in planets:
        # calc_ut returns (position_tuple, flags_used)
        # position_tuple = (longitude, latitude, distance, lon_speed, lat_speed, dist_speed)
        pos, _ = eph.calc_ut(jd, planet_id, flags)

        longitude = pos[0]
        latitude = pos[1]
        distance = pos[2]
        speed = pos[3]

        # Mark retrograde planets
        retro = " R" if speed < 0 else ""

        print(
            f"{name:<12} {longitude:>10.4f}° {latitude:>9.4f}° {distance:>9.4f} {speed:>9.4f}{retro}"
        )


def example_house_cusps() -> None:
    """Example 3: Calculate house cusps and angles."""
    print("\n" + "=" * 60)
    print("Example 3: House Cusps and Angles")
    print("=" * 60)

    # Location: Rome, Italy
    latitude = 41.9028
    longitude = 12.4964

    # Time: January 1, 2024 at noon
    jd = eph.julday(2024, 1, 1, 12.0)

    print(f"\nHouse cusps for Rome ({latitude}°N, {longitude}°E)")
    print(f"at JD {jd}\n")

    # Calculate houses using Placidus system
    # houses returns (cusps, ascmc)
    # cusps = tuple of 12 house cusp longitudes
    # ascmc = (Asc, MC, ARMC, Vertex, EquatorialAsc, Co-Asc Koch, Co-Asc Munkasey, Polar Asc)
    cusps, ascmc = eph.houses(jd, latitude, longitude, ord("P"))  # P = Placidus

    print("Placidus House Cusps:")
    for i in range(12):
        print(f"  House {i + 1:2}: {cusps[i]:>10.4f}°")

    print("\nAngles:")
    print(f"  Ascendant:     {ascmc[0]:>10.4f}°")
    print(f"  Midheaven:     {ascmc[1]:>10.4f}°")
    print(f"  ARMC:          {ascmc[2]:>10.4f}°")
    print(f"  Vertex:        {ascmc[3]:>10.4f}°")

    # Try different house systems
    print("\nComparing Ascendants in different house systems:")
    systems = [
        ("P", "Placidus"),
        ("K", "Koch"),
        ("W", "Whole Sign"),
        ("E", "Equal"),
        ("C", "Campanus"),
        ("R", "Regiomontanus"),
    ]

    for code, name in systems:
        cusps_sys, ascmc_sys = eph.houses(jd, latitude, longitude, ord(code))
        print(f"  {name:15} ASC: {ascmc_sys[0]:>10.4f}°  MC: {ascmc_sys[1]:>10.4f}°")


def example_zodiac_formatting() -> None:
    """Example 4: Format positions in zodiac notation."""
    print("\n" + "=" * 60)
    print("Example 4: Zodiac Notation")
    print("=" * 60)

    # Helper function to convert longitude to zodiac sign
    def format_zodiac(longitude: float) -> str:
        """Convert ecliptic longitude to zodiac notation."""
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
        sign_num = int(longitude / 30)
        pos_in_sign = longitude - (sign_num * 30)
        degrees = int(pos_in_sign)
        minutes = int((pos_in_sign - degrees) * 60)
        seconds = int(((pos_in_sign - degrees) * 60 - minutes) * 60)
        return f"{degrees:02d}° {minutes:02d}' {seconds:02d}\" {signs[sign_num]}"

    jd = eph.julday(2024, 1, 1, 12.0)

    print("\nPlanetary positions in zodiac notation:\n")

    planets = [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MERCURY, "Mercury"),
        (VENUS, "Venus"),
        (MARS, "Mars"),
    ]

    for planet_id, name in planets:
        pos, _ = eph.calc_ut(jd, planet_id, FLG_SWIEPH | FLG_SPEED)
        formatted = format_zodiac(pos[0])
        retro = " (Retrograde)" if pos[3] < 0 else ""
        print(f"  {name:10} {formatted}{retro}")


def example_ayanamsha() -> None:
    """Example 5: Ayanamsha (precession of equinoxes)."""
    print("\n" + "=" * 60)
    print("Example 5: Ayanamsha Values")
    print("=" * 60)

    jd = eph.julday(2024, 1, 1, 12.0)

    print(f"\nAyanamsha values at JD {jd}:\n")

    # Common sidereal modes
    from libephemeris.constants import (
        SIDM_FAGAN_BRADLEY,
        SIDM_LAHIRI,
        SIDM_RAMAN,
        SIDM_TRUE_CITRA,
        SIDM_KRISHNAMURTI,
    )

    modes = [
        (SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
        (SIDM_LAHIRI, "Lahiri"),
        (SIDM_RAMAN, "Raman"),
        (SIDM_TRUE_CITRA, "True Citra"),
        (SIDM_KRISHNAMURTI, "Krishnamurti"),
    ]

    for mode, name in modes:
        eph.set_sid_mode(mode)
        ayanamsha = eph.get_ayanamsa_ut(jd)
        print(f"  {name:20} {ayanamsha:>10.6f}°")


def example_planet_name() -> None:
    """Example 6: Get planet names."""
    print("\n" + "=" * 60)
    print("Example 6: Planet Names")
    print("=" * 60)

    print("\nPlanet IDs and their names:\n")

    planet_ids = [
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
    ]

    for pid in planet_ids:
        name = eph.get_planet_name(pid)
        print(f"  ID {pid:2}: {name}")


def main() -> None:
    """Run all basic usage examples."""
    print("\n" + "#" * 60)
    print("# libephemeris Basic Usage Examples")
    print("#" * 60)

    example_julian_day()
    example_planetary_positions()
    example_house_cusps()
    example_zodiac_formatting()
    example_ayanamsha()
    example_planet_name()

    print("\n" + "#" * 60)
    print("# Basic examples completed!")
    print("#" * 60 + "\n")


if __name__ == "__main__":
    main()
