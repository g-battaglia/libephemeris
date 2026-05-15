#!/usr/bin/env python3
"""
Fixed Stars Example for libephemeris.

This script demonstrates how to work with fixed stars:
1. Looking up stars by name
2. Calculating star positions
3. Getting star magnitudes
4. Finding stars near a given longitude
5. Working with proper motion

Requirements:
    pip install libephemeris

Usage:
    python examples/fixed_stars.py
"""

from __future__ import annotations

import libephemeris as eph
from libephemeris.constants import (
    SUN,
    FLG_SWIEPH,
    FLG_SPEED,
)


# Notable fixed stars with traditional astrological significance
NOTABLE_STARS = [
    ("Aldebaran", "Alpha Tauri - The Eye of the Bull"),
    ("Regulus", "Alpha Leonis - The Heart of the Lion"),
    ("Spica", "Alpha Virginis - The Wheat Sheaf"),
    ("Antares", "Alpha Scorpii - The Heart of the Scorpion"),
    ("Sirius", "Alpha Canis Majoris - The Dog Star"),
    ("Vega", "Alpha Lyrae - The Falling Eagle"),
    ("Altair", "Alpha Aquilae - The Flying Eagle"),
    ("Deneb", "Alpha Cygni - The Tail of the Swan"),
    ("Fomalhaut", "Alpha Piscis Austrini - The Lonely One"),
    ("Arcturus", "Alpha Bootis - The Bear Guard"),
    ("Betelgeuse", "Alpha Orionis - Shoulder of the Giant"),
    ("Rigel", "Beta Orionis - Foot of Orion"),
    ("Polaris", "Alpha Ursae Minoris - The North Star"),
    ("Capella", "Alpha Aurigae - The Little Goat"),
    ("Procyon", "Alpha Canis Minoris - Before the Dog"),
]


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
    sign_num = int(longitude / 30) % 12
    pos_in_sign = longitude % 30
    degrees = int(pos_in_sign)
    minutes = int((pos_in_sign - degrees) * 60)
    return f"{degrees:02d}° {minutes:02d}' {signs[sign_num]}"


def example_star_positions() -> None:
    """Example 1: Calculate positions of notable fixed stars."""
    print("=" * 70)
    print("Example 1: Fixed Star Positions")
    print("=" * 70)

    jd = eph.julday(2024, 1, 1, 0.0)

    print(f"\nFixed star positions at JD {jd} (January 1, 2024):\n")
    print(f"{'Star':<14}{'Longitude':>18}{'Latitude':>12}{'Mag':>6}")
    print("-" * 52)

    for star_name, description in NOTABLE_STARS[:10]:
        try:
            # Get star position
            pos, _, retflag = eph.fixstar_ut(star_name, jd, FLG_SWIEPH)

            # Get magnitude (returns tuple of (magnitude, star_name))
            mag, _ = eph.fixstar_mag(star_name)

            formatted = format_zodiac(pos[0])
            print(f"{star_name:<14}{formatted:>18}{pos[1]:>11.4f}°{mag:>6.2f}")
        except Exception as e:
            print(f"{star_name:<14} - Error: {e}")


def example_star_lookup() -> None:
    """Example 2: Look up stars by partial name."""
    print("\n" + "=" * 70)
    print("Example 2: Star Lookup by Name")
    print("=" * 70)

    jd = eph.julday(2024, 1, 1, 0.0)

    print("\nLooking up stars with partial names (fixstar2_ut):\n")

    # Partial names to search
    search_terms = ["Reg", "Spic", "Sir", "Alp"]

    for term in search_terms:
        try:
            pos, full_name, retflag = eph.fixstar2_ut(term, jd, FLG_SWIEPH)

            print(f"  Search '{term}' -> {full_name}")
            print(f"    Position: {format_zodiac(pos[0])}")
        except Exception as e:
            print(f"  Search '{term}' -> Error: {e}")


def example_star_magnitude() -> None:
    """Example 3: Get star magnitudes and brightness."""
    print("\n" + "=" * 70)
    print("Example 3: Star Magnitudes")
    print("=" * 70)

    print("\nVisual magnitudes of bright stars:\n")
    print("(Lower magnitude = brighter; Sirius is the brightest star in the sky)\n")

    # Stars sorted roughly by brightness
    bright_stars = [
        "Sirius",
        "Canopus",
        "Arcturus",
        "Vega",
        "Capella",
        "Rigel",
        "Procyon",
        "Betelgeuse",
        "Altair",
        "Aldebaran",
    ]

    star_mags = []
    for star_name in bright_stars:
        try:
            mag, _resolved_name = eph.fixstar_mag(star_name)
            star_mags.append((star_name, mag))
        except Exception:
            pass

    # Sort by magnitude (brightest first)
    star_mags.sort(key=lambda x: x[1])

    print(f"{'Star':<14}{'Magnitude':>10}{'Brightness':>15}")
    print("-" * 42)

    for star_name, mag in star_mags:
        # Relative brightness description
        if mag < 0:
            brightness = "Very Bright"
        elif mag < 1:
            brightness = "Bright"
        elif mag < 2:
            brightness = "Moderate"
        else:
            brightness = "Dim"

        print(f"{star_name:<14}{mag:>10.2f}    {brightness:<15}")


def example_proper_motion() -> None:
    """Example 4: Demonstrate proper motion over time."""
    print("\n" + "=" * 70)
    print("Example 4: Proper Motion Over Time")
    print("=" * 70)

    print("\nFixed stars aren't truly fixed - they have proper motion.")
    print("Some stars move noticeably over centuries.\n")

    # Stars with notable proper motion
    high_motion_stars = ["Arcturus", "Sirius", "Procyon", "61Cyg"]

    years = [1900, 2000, 2100, 2200]

    for star_name in high_motion_stars:
        print(f"\n{star_name}:")
        print("-" * 40)

        positions = []
        for year in years:
            jd = eph.julday(year, 1, 1, 0.0)
            try:
                pos, _, retflag = eph.fixstar_ut(star_name, jd, FLG_SWIEPH)
                positions.append((year, pos[0], pos[1]))
                print(f"  {year}: {pos[0]:>10.6f}° lon, {pos[1]:>10.6f}° lat")
            except Exception:
                pass

        # Calculate motion per century
        if len(positions) >= 2:
            first = positions[0]
            last = positions[-1]
            years_diff = last[0] - first[0]
            lon_motion = (
                (last[1] - first[1]) / (years_diff / 100) * 3600
            )  # arcsec/century
            lat_motion = (
                (last[2] - first[2]) / (years_diff / 100) * 3600
            )  # arcsec/century
            print(
                f'  Motion: {lon_motion:+.2f}"/century (lon), {lat_motion:+.2f}"/century (lat)'
            )


def example_conjunction_with_planet() -> None:
    """Example 5: Find stars near a planet's position."""
    print("\n" + "=" * 70)
    print("Example 5: Stars Near Planetary Position")
    print("=" * 70)

    jd = eph.julday(2024, 8, 15, 0.0)

    # Get Sun's position
    sun_pos, _ = eph.calc_ut(jd, SUN, FLG_SWIEPH | FLG_SPEED)
    sun_lon = sun_pos[0]

    print(f"\nSun position on Aug 15, 2024: {format_zodiac(sun_lon)}")
    print("\nStars within 5° of the Sun:\n")

    stars_near = []

    for star_name, description in NOTABLE_STARS:
        try:
            pos, _, retflag = eph.fixstar_ut(star_name, jd, FLG_SWIEPH)
            # Calculate angular separation
            diff = abs(pos[0] - sun_lon)
            if diff > 180:
                diff = 360 - diff

            if diff <= 5:
                stars_near.append((star_name, pos[0], diff))
        except Exception:
            pass

    if stars_near:
        stars_near.sort(key=lambda x: x[2])
        print(f"{'Star':<14}{'Position':>20}{'Distance':>12}")
        print("-" * 48)
        for name, lon, dist in stars_near:
            print(f"{name:<14}{format_zodiac(lon):>20}{dist:>11.2f}°")
    else:
        print("  No notable stars within 5° of the Sun")

    # Also check Regulus which is near Leo
    print("\nNote: Stars near the ecliptic are most likely to be occulted by planets.")


def example_royal_stars() -> None:
    """Example 6: The Royal Stars of Persia."""
    print("\n" + "=" * 70)
    print("Example 6: The Royal Stars of Persia")
    print("=" * 70)

    print("\nThe four Royal Stars mark the four corners of the heavens")
    print(
        "in ancient Persian astrology (circa 3000 BCE they marked the equinoxes and solstices):\n"
    )

    royal_stars = [
        ("Aldebaran", "Watcher of the East", "Archangel Michael"),
        ("Regulus", "Watcher of the North", "Archangel Raphael"),
        ("Antares", "Watcher of the West", "Archangel Uriel"),
        ("Fomalhaut", "Watcher of the South", "Archangel Gabriel"),
    ]

    jd = eph.julday(2024, 1, 1, 0.0)

    print(f"{'Star':<12}{'Guardian':<22}{'Archangel':<18}{'Position':>18}")
    print("-" * 72)

    for star_name, title, archangel in royal_stars:
        try:
            pos, _, retflag = eph.fixstar_ut(star_name, jd, FLG_SWIEPH)
            formatted = format_zodiac(pos[0])
            print(f"{star_name:<12}{title:<22}{archangel:<18}{formatted:>18}")
        except Exception as e:
            print(f"{star_name:<12} Error: {e}")


def example_behenian_stars() -> None:
    """Example 7: The 15 Behenian Fixed Stars."""
    print("\n" + "=" * 70)
    print("Example 7: The Behenian Fixed Stars")
    print("=" * 70)

    print("\nThe 15 Behenian stars are used in medieval astrological magic.")
    print("Each has associated planets, gemstones, and plants.\n")

    # Behenian stars with their ruling planets
    behenian_stars = [
        ("Algol", "Saturn/Jupiter"),
        ("Pleiades", "Moon/Mars"),
        ("Aldebaran", "Mars"),
        ("Capella", "Jupiter/Saturn"),
        ("Sirius", "Jupiter"),
        ("Procyon", "Mercury/Mars"),
        ("Regulus", "Jupiter/Mars"),
        ("Polaris", "Saturn/Venus"),
        ("Algorab", "Saturn/Mars"),
        ("Spica", "Venus/Mercury"),
        ("Arcturus", "Jupiter/Mars"),
        ("Alphecca", "Venus/Mercury"),
        ("Antares", "Mars/Jupiter"),
        ("Vega", "Mercury/Venus"),
        ("Deneb", "Venus/Mercury"),
    ]

    jd = eph.julday(2024, 1, 1, 0.0)

    print(f"{'Star':<12}{'Planets':<16}{'Position':>22}")
    print("-" * 52)

    for star_name, planets in behenian_stars:
        try:
            pos, _, retflag = eph.fixstar_ut(star_name, jd, FLG_SWIEPH)
            formatted = format_zodiac(pos[0])
            print(f"{star_name:<12}{planets:<16}{formatted:>22}")
        except Exception:
            print(f"{star_name:<12}{planets:<16}{'(not found)':>22}")


def main() -> None:
    """Run all fixed stars examples."""
    print("\n" + "#" * 70)
    print("# libephemeris Fixed Stars Examples")
    print("#" * 70)

    example_star_positions()
    example_star_lookup()
    example_star_magnitude()
    example_proper_motion()
    example_conjunction_with_planet()
    example_royal_stars()
    example_behenian_stars()

    print("\n" + "#" * 70)
    print("# Fixed stars examples completed!")
    print("#" * 70 + "\n")


if __name__ == "__main__":
    main()
