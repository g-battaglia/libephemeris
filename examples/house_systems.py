#!/usr/bin/env python3
"""
House Systems Example for libephemeris.

This script demonstrates how to work with different house systems:
1. Overview of available house systems
2. Calculating house cusps with different systems
3. Finding house position of a planet
4. Handling polar latitudes
5. Using Gauquelin sectors

Requirements:
    pip install libephemeris

Usage:
    python examples/house_systems.py
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
    # Calculation flags
    FLG_SWIEPH,
)


# House system codes and their descriptions
HOUSE_SYSTEMS = {
    "P": (
        "Placidus",
        "Most popular Western system. Time-based divisions using semi-arcs.",
    ),
    "K": ("Koch", "Similar to Placidus but uses birthplace latitude more directly."),
    "O": ("Porphyry", "Trisects the semi-arcs between angles. Simple and ancient."),
    "R": ("Regiomontanus", "Space-based divisions using celestial equator."),
    "C": ("Campanus", "Space-based using prime vertical great circles."),
    "E": ("Equal", "Each house is exactly 30°. ASC determines house 1."),
    "W": ("Whole Sign", "Each sign = one house. Oldest house system."),
    "B": ("Alcabitius", "Medieval Arabic system. Semi-arc divisions."),
    "M": ("Morinus", "Uses equatorial divisions. No house cusps at poles."),
    "T": ("Topocentric", "Modern system by Polich/Page. Good at high latitudes."),
    "X": ("Meridian/Axial", "MC-based equal divisions on the ecliptic."),
    "H": ("Azimuthal/Horizontal", "Based on horizon and prime vertical."),
    "V": ("Vehlow Equal", "Equal houses with ASC in middle of 1st house."),
    "U": ("Krusinski", "Based on polar axis. Works at all latitudes."),
    "G": ("Gauquelin Sectors", "36 sectors for statistical research."),
    "Y": ("APC", "Astro-Psychological-Combination. Equal from ASC."),
    "N": ("Natural Graduation", "Hellenistic proportional houses."),
}


def format_zodiac(longitude: float) -> str:
    """Convert ecliptic longitude to short zodiac notation."""
    signs = ["Ar", "Ta", "Ge", "Cn", "Le", "Vi", "Li", "Sc", "Sg", "Cp", "Aq", "Pi"]
    sign_num = int(longitude / 30) % 12
    degrees = longitude % 30
    return f"{degrees:05.2f}° {signs[sign_num]}"


def example_house_systems_overview() -> None:
    """Example 1: Overview of available house systems."""
    print("=" * 70)
    print("Example 1: House Systems Overview")
    print("=" * 70)

    print("\nAvailable house systems in libephemeris:\n")
    print(f"{'Code':<6}{'Name':<20}Description")
    print("-" * 70)

    for code, (name, desc) in HOUSE_SYSTEMS.items():
        print(f"  {code:<4} {name:<18} {desc}")


def example_compare_house_systems() -> None:
    """Example 2: Compare house cusps across different systems."""
    print("\n" + "=" * 70)
    print("Example 2: Compare House Systems")
    print("=" * 70)

    # Location: London, UK
    latitude = 51.5074
    longitude = -0.1278

    # Time: January 15, 2024 at 14:30
    jd = eph.julday(2024, 1, 15, 14.5)

    print(f"\nLocation: London ({latitude}°N, {longitude}°W)")
    print("Time: 2024-01-15 14:30 UT\n")

    # Systems to compare
    systems = ["P", "K", "O", "R", "C", "E", "W", "B"]

    # Print header
    print(f"{'System':<14}", end="")
    for i in range(1, 7):
        print(f"{'House ' + str(i):>10}", end="")
    print()
    print("-" * 74)

    for sys_code in systems:
        sys_name = HOUSE_SYSTEMS[sys_code][0]
        cusps, ascmc = eph.houses(jd, latitude, longitude, ord(sys_code))

        print(f"{sys_name:<14}", end="")
        for i in range(6):
            print(f"{cusps[i]:>10.2f}", end="")
        print()

    print("\nNote: All values are in degrees of the ecliptic (0-360)")


def example_angles_and_vertices() -> None:
    """Example 3: Calculate angles (ASC, MC, Vertex, etc.)."""
    print("\n" + "=" * 70)
    print("Example 3: Angles and Special Points")
    print("=" * 70)

    # Location: New York City
    latitude = 40.7128
    longitude = -74.0060

    jd = eph.julday(2024, 3, 20, 12.0)

    print(f"\nLocation: New York City ({latitude}°N, {longitude}°W)")
    print("Time: 2024-03-20 12:00 UT (Spring Equinox)\n")

    # Calculate houses (system doesn't affect angles, but we need to call it)
    cusps, ascmc = eph.houses(jd, latitude, longitude, ord("P"))

    print("Angles (ASCMC array):")
    print("-" * 40)
    print(f"  Ascendant (ASC):           {format_zodiac(ascmc[0])}")
    print(f"  Midheaven (MC):            {format_zodiac(ascmc[1])}")
    print(f"  ARMC (Sidereal Time):      {ascmc[2]:.4f}°")
    print(f"  Vertex:                    {format_zodiac(ascmc[3])}")
    print(f"  Equatorial Ascendant:      {format_zodiac(ascmc[4])}")
    print(f"  Co-Ascendant (Koch):       {format_zodiac(ascmc[5])}")
    print(f"  Co-Ascendant (Munkasey):   {format_zodiac(ascmc[6])}")
    print(f"  Polar Ascendant:           {format_zodiac(ascmc[7])}")

    # Derived points
    desc = (ascmc[0] + 180) % 360
    ic = (ascmc[1] + 180) % 360
    anti_vertex = (ascmc[3] + 180) % 360

    print("\nDerived Points:")
    print("-" * 40)
    print(f"  Descendant (DSC):          {format_zodiac(desc)}")
    print(f"  Imum Coeli (IC):           {format_zodiac(ic)}")
    print(f"  Anti-Vertex:               {format_zodiac(anti_vertex)}")


def example_house_position() -> None:
    """Example 4: Find which house a planet is in."""
    print("\n" + "=" * 70)
    print("Example 4: House Position of Planets")
    print("=" * 70)

    # Location: Los Angeles
    latitude = 34.0522
    longitude = -118.2437

    jd = eph.julday(2024, 6, 21, 12.0)  # Summer solstice

    print(f"\nLocation: Los Angeles ({latitude}°N, {longitude}°W)")
    print("Time: 2024-06-21 12:00 UT (Summer Solstice)\n")

    # Get house cusps for Placidus
    cusps, ascmc = eph.houses(jd, latitude, longitude, ord("P"))

    # Get planets
    planets = [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MERCURY, "Mercury"),
        (VENUS, "Venus"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
    ]

    print(f"{'Planet':<12}{'Position':>15}{'House':>8}")
    print("-" * 38)

    # Get true obliquity for house_pos calculation
    # Standard value is approximately 23.4 degrees
    true_obliquity = 23.4393  # Approximate value for modern epoch

    for planet_id, name in planets:
        pos, _ = eph.calc_ut(jd, planet_id, FLG_SWIEPH)

        # Use house_pos to get the house position
        # Returns a float: whole part = house number, decimal = position within house
        house_pos = eph.house_pos(
            ascmc[2],  # ARMC in degrees
            latitude,
            true_obliquity,  # True obliquity of ecliptic
            ord("P"),  # House system
            pos[0],  # Planet longitude
            pos[1],  # Planet latitude
        )

        house_num = int(house_pos)
        pos_in_house = (house_pos - house_num) * 100

        print(
            f"{name:<12}{format_zodiac(pos[0]):>15}   {house_num:>3} ({pos_in_house:.0f}%)"
        )


def example_polar_latitudes() -> None:
    """Example 5: Handling polar latitudes."""
    print("\n" + "=" * 70)
    print("Example 5: Polar Latitude Handling")
    print("=" * 70)

    print("\nAt extreme latitudes (>66.5°), some house systems fail because")
    print("certain points (like the MC) may not intersect the ecliptic properly.\n")

    # Polar location: Tromsø, Norway
    latitude = 69.6496
    longitude = 19.0134

    jd = eph.julday(2024, 6, 21, 12.0)  # Summer solstice

    print(f"Location: Tromsø, Norway ({latitude}°N, {longitude}°E)")
    print("Time: 2024-06-21 12:00 UT\n")

    # Test different house systems at polar latitude
    systems_to_test = ["P", "K", "O", "R", "C", "E", "W", "T"]

    print(f"{'System':<16}{'Status':<12}{'ASC':>12}{'MC':>12}")
    print("-" * 54)

    for sys_code in systems_to_test:
        sys_name = HOUSE_SYSTEMS[sys_code][0]

        try:
            cusps, ascmc = eph.houses(jd, latitude, longitude, ord(sys_code))
            print(f"{sys_name:<16}{'OK':<12}{ascmc[0]:>11.2f}°{ascmc[1]:>11.2f}°")
        except Exception:
            print(f"{sys_name:<16}{'FAILED':<12} (Polar circle issue)")

    print("\nRecommendations for polar latitudes:")
    print("  - Use Whole Sign (W) - always works")
    print("  - Use Equal (E) - always works")
    print("  - Use Topocentric (T) - designed for high latitudes")


def example_extended_houses() -> None:
    """Example 6: Using houses_ex for sidereal calculations."""
    print("\n" + "=" * 70)
    print("Example 6: Sidereal House Calculations")
    print("=" * 70)

    from libephemeris.constants import FLG_SIDEREAL, SIDM_LAHIRI

    # Location: Chennai, India
    latitude = 13.0827
    longitude = 80.2707

    jd = eph.julday(2024, 1, 14, 6.0)  # Pongal/Makar Sankranti

    print(f"\nLocation: Chennai, India ({latitude}°N, {longitude}°E)")
    print("Time: 2024-01-14 06:00 UT (Pongal)\n")

    # Set sidereal mode first
    eph.set_sid_mode(SIDM_LAHIRI)

    # Calculate tropical houses
    cusps_trop, ascmc_trop = eph.houses(jd, latitude, longitude, ord("W"))

    # Calculate sidereal houses using houses_ex
    cusps_sid, ascmc_sid = eph.houses_ex(
        jd, latitude, longitude, ord("W"), FLG_SIDEREAL
    )

    # Get ayanamsha for reference
    aya = eph.get_ayanamsa_ut(jd)

    print(f"Ayanamsha (Lahiri): {aya:.4f}°\n")

    print(f"{'Bhava':<10}{'Tropical':>14}{'Sidereal':>14}{'Difference':>14}")
    print("-" * 52)

    for i in range(12):
        diff = cusps_trop[i] - cusps_sid[i]
        # Normalize difference
        if diff > 180:
            diff -= 360
        elif diff < -180:
            diff += 360
        print(
            f"  {i + 1:<6}  {cusps_trop[i]:>11.2f}°  {cusps_sid[i]:>11.2f}°  {diff:>11.2f}°"
        )

    print(f"\n{'Angle':<10}{'Tropical':>14}{'Sidereal':>14}")
    print("-" * 38)
    print(f"  ASC     {ascmc_trop[0]:>11.2f}°  {ascmc_sid[0]:>11.2f}°")
    print(f"  MC      {ascmc_trop[1]:>11.2f}°  {ascmc_sid[1]:>11.2f}°")


def example_gauquelin_sectors() -> None:
    """Example 7: Gauquelin sectors for statistical research."""
    print("\n" + "=" * 70)
    print("Example 7: Gauquelin Sectors")
    print("=" * 70)

    print("\nGauquelin sectors divide the diurnal circle into 36 sectors,")
    print("used in astrological research on planetary prominence.\n")

    # Sample data
    latitude = 48.8566  # Paris
    longitude = 2.3522

    jd = eph.julday(2024, 1, 15, 8.0)

    print(f"Location: Paris ({latitude}°N, {longitude}°E)")
    print("Time: 2024-01-15 08:00 UT\n")

    # Get Gauquelin sectors using the G house system
    # Note: houses returns 12 house cusps even for Gauquelin
    # For full 36 sectors, use gauquelin_sector() function
    cusps, ascmc = eph.houses(jd, latitude, longitude, ord("G"))

    print("Gauquelin house cusps (12 primary sectors):")
    print("-" * 45)

    # Show all 12 cusps
    for i in range(12):
        print(f"  Sector {i + 1:2}: {cusps[i]:>10.4f}°")

    # Calculate sector position for a planet
    sun_pos, _ = eph.calc_ut(jd, SUN, FLG_SWIEPH)

    print(f"\nSun longitude: {sun_pos[0]:.4f}°")
    print("\nNote: For the full 36 Gauquelin sectors, use the")
    print("      gauquelin_sector() function.")


def main() -> None:
    """Run all house systems examples."""
    print("\n" + "#" * 70)
    print("# libephemeris House Systems Examples")
    print("#" * 70)

    example_house_systems_overview()
    example_compare_house_systems()
    example_angles_and_vertices()
    example_house_position()
    example_polar_latitudes()
    example_extended_houses()
    example_gauquelin_sectors()

    print("\n" + "#" * 70)
    print("# House systems examples completed!")
    print("#" * 70 + "\n")


if __name__ == "__main__":
    main()
