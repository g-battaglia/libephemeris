#!/usr/bin/env python3
"""
Sidereal Calculations Example for libephemeris.

This script demonstrates how to work with sidereal (Vedic) astrology:
1. Setting different ayanamsha modes
2. Calculating sidereal planetary positions
3. Comparing tropical vs sidereal positions
4. Working with popular Indian ayanamshas

Requirements:
    pip install libephemeris

Usage:
    python examples/sidereal_calculations.py
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
    # Calculation flags
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_SIDEREAL,
    # Sidereal modes (ayanamshas)
    SIDM_FAGAN_BRADLEY,
    SIDM_LAHIRI,
    SIDM_RAMAN,
    SIDM_USHASHASHI,
    SIDM_KRISHNAMURTI,
    SIDM_DJWHAL_KHUL,
    SIDM_YUKTESHWAR,
    SIDM_JN_BHASIN,
    SIDM_BABYL_KUGLER1,
    SIDM_TRUE_CITRA,
    SIDM_TRUE_REVATI,
    SIDM_TRUE_PUSHYA,
    SIDM_GALCENT_0SAG,
)


# Zodiac signs - same for both tropical and sidereal
SIGNS = [
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

# Vedic names for signs (nakshatras use different divisions)
RASHIS = [
    "Mesha",
    "Vrishabha",
    "Mithuna",
    "Karka",
    "Simha",
    "Kanya",
    "Tula",
    "Vrischika",
    "Dhanu",
    "Makara",
    "Kumbha",
    "Meena",
]


def format_zodiac(longitude: float, use_vedic_names: bool = False) -> str:
    """Convert ecliptic longitude to zodiac notation."""
    signs = RASHIS if use_vedic_names else SIGNS
    sign_num = int(longitude / 30) % 12
    pos_in_sign = longitude % 30
    degrees = int(pos_in_sign)
    minutes = int((pos_in_sign - degrees) * 60)
    return f"{degrees:02d}° {minutes:02d}' {signs[sign_num]}"


def example_ayanamsha_overview() -> None:
    """Example 1: Overview of available ayanamshas."""
    print("=" * 60)
    print("Example 1: Ayanamsha Overview")
    print("=" * 60)

    print("\nThe ayanamsha is the angular difference between the tropical")
    print("and sidereal zodiacs, caused by the precession of equinoxes.")
    print("\nPopular ayanamsha systems:\n")

    jd = eph.julday(2024, 1, 1, 0.0)

    # List of common ayanamshas with descriptions
    ayanamshas = [
        (SIDM_LAHIRI, "Lahiri", "Official Indian government standard"),
        (SIDM_RAMAN, "Raman", "B.V. Raman's ayanamsha"),
        (SIDM_KRISHNAMURTI, "Krishnamurti", "KP astrology system"),
        (SIDM_FAGAN_BRADLEY, "Fagan-Bradley", "Western sidereal standard"),
        (SIDM_TRUE_CITRA, "True Citra", "Spica at 0° Libra"),
        (SIDM_YUKTESHWAR, "Yukteshwar", "Sri Yukteshwar's calculation"),
        (SIDM_JN_BHASIN, "JN Bhasin", "J.N. Bhasin's ayanamsha"),
        (SIDM_USHASHASHI, "Ushashashi", "Based on Surya Siddhanta"),
        (SIDM_DJWHAL_KHUL, "Djwhal Khul", "Theosophical tradition"),
        (SIDM_TRUE_REVATI, "True Revati", "Zeta Piscium at 29°50' Pisces"),
        (SIDM_TRUE_PUSHYA, "True Pushya", "Delta Cancri at 16° Cancer"),
        (SIDM_BABYL_KUGLER1, "Babylonian (Kugler 1)", "Ancient Babylonian"),
        (SIDM_GALCENT_0SAG, "Galactic Center", "Galactic center at 0° Sagittarius"),
    ]

    print(f"{'Ayanamsha':<25} {'Value':>12}  Description")
    print("-" * 65)

    for mode, name, desc in ayanamshas:
        eph.set_sid_mode(mode)
        aya = eph.get_ayanamsa_ut(jd)
        print(f"{name:<25} {aya:>11.6f}°  {desc}")


def example_tropical_vs_sidereal() -> None:
    """Example 2: Compare tropical and sidereal positions."""
    print("\n" + "=" * 60)
    print("Example 2: Tropical vs Sidereal Positions")
    print("=" * 60)

    jd = eph.julday(2024, 3, 21, 12.0)  # Spring equinox 2024

    print(f"\nComparison at JD {jd} (Spring Equinox 2024):\n")

    # Set Lahiri ayanamsha
    eph.set_sid_mode(SIDM_LAHIRI)
    ayanamsha = eph.get_ayanamsa_ut(jd)
    print(f"Lahiri Ayanamsha: {ayanamsha:.6f}°\n")

    planets = [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
        (TRUE_NODE, "Rahu"),
    ]

    print(f"{'Planet':<12} {'Tropical':>22} {'Sidereal (Lahiri)':>22}")
    print("-" * 58)

    for planet_id, name in planets:
        # Tropical position
        pos_trop, _ = eph.calc_ut(jd, planet_id, FLG_SWIEPH | FLG_SPEED)

        # Sidereal position
        pos_sid, _ = eph.calc_ut(
            jd, planet_id, FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL
        )

        trop_str = format_zodiac(pos_trop[0])
        sid_str = format_zodiac(pos_sid[0], use_vedic_names=True)

        print(f"{name:<12} {trop_str:>22} {sid_str:>22}")

    print("\nNote: Ketu (South Node) is 180° opposite from Rahu")


def example_vedic_chart() -> None:
    """Example 3: Calculate a Vedic (Sidereal) natal chart."""
    print("\n" + "=" * 60)
    print("Example 3: Vedic Natal Chart Calculation")
    print("=" * 60)

    # Example birth data: Mumbai, India
    year, month, day, hour = 1990, 5, 15, 10.5
    latitude = 19.0760  # Mumbai latitude
    longitude = 72.8777  # Mumbai longitude

    jd = eph.julday(year, month, day, hour)

    print(
        f"\nBirth data: {year}-{month:02d}-{day:02d} at {int(hour):02d}:{int((hour % 1) * 60):02d}"
    )
    print(f"Location: Mumbai ({latitude}°N, {longitude}°E)")

    # Set Lahiri ayanamsha (most common in India)
    eph.set_sid_mode(SIDM_LAHIRI)
    ayanamsha = eph.get_ayanamsa_ut(jd)
    print(f"Ayanamsha (Lahiri): {ayanamsha:.6f}°\n")

    # Calculate sidereal planets
    print("Grahas (Planets):")
    print("-" * 40)

    planets = [
        (SUN, "Surya (Sun)"),
        (MOON, "Chandra (Moon)"),
        (MERCURY, "Budha (Mercury)"),
        (VENUS, "Shukra (Venus)"),
        (MARS, "Mangal (Mars)"),
        (JUPITER, "Guru (Jupiter)"),
        (SATURN, "Shani (Saturn)"),
        (TRUE_NODE, "Rahu"),
    ]

    flags = FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL

    for planet_id, name in planets:
        pos, _ = eph.calc_ut(jd, planet_id, flags)
        formatted = format_zodiac(pos[0], use_vedic_names=True)
        retro = " (Vakri)" if pos[3] < 0 else ""
        print(f"  {name:20} {formatted}{retro}")

    # Calculate Ketu (opposite of Rahu)
    rahu_pos, _ = eph.calc_ut(jd, TRUE_NODE, flags)
    ketu_lon = (rahu_pos[0] + 180) % 360
    print(f"  {'Ketu':20} {format_zodiac(ketu_lon, use_vedic_names=True)}")

    # Calculate houses (Bhava) - using Whole Sign houses (common in Vedic)
    print("\nBhava Sphutas (House Cusps) - Whole Sign:")
    print("-" * 40)

    # Get ascendant in sidereal
    cusps, ascmc = eph.houses_ex(jd, latitude, longitude, ord("W"), FLG_SIDEREAL)

    for i in range(12):
        print(f"  Bhava {i + 1:2}: {format_zodiac(cusps[i], use_vedic_names=True)}")

    print(f"\n  Lagna (Ascendant): {format_zodiac(ascmc[0], use_vedic_names=True)}")


def example_ayanamsha_change_over_time() -> None:
    """Example 4: Show ayanamsha change over centuries."""
    print("\n" + "=" * 60)
    print("Example 4: Ayanamsha Change Over Time")
    print("=" * 60)

    print("\nThe ayanamsha increases by approximately 50.3 arcseconds per year")
    print("due to the precession of the Earth's axis.\n")

    eph.set_sid_mode(SIDM_LAHIRI)

    years = [1900, 1950, 2000, 2024, 2050, 2100]

    print(f"{'Year':>6}  {'Ayanamsha':>12}  {'Change from 2000':>18}")
    print("-" * 40)

    jd_2000 = eph.julday(2000, 1, 1, 0.0)
    aya_2000 = eph.get_ayanamsa_ut(jd_2000)

    for year in years:
        jd = eph.julday(year, 1, 1, 0.0)
        aya = eph.get_ayanamsa_ut(jd)
        change = aya - aya_2000
        sign = "+" if change >= 0 else ""
        print(f"{year:>6}  {aya:>11.6f}°  {sign}{change:>17.6f}°")

    # Rate of precession
    jd1 = eph.julday(2024, 1, 1, 0.0)
    jd2 = eph.julday(2025, 1, 1, 0.0)
    aya1 = eph.get_ayanamsa_ut(jd1)
    aya2 = eph.get_ayanamsa_ut(jd2)
    annual_precession = aya2 - aya1

    print(f"\nAnnual precession rate: {annual_precession * 3600:.2f} arcseconds")
    print(f"                        ({annual_precession:.6f}°)")


def example_custom_ayanamsha() -> None:
    """Example 5: Set a custom ayanamsha."""
    print("\n" + "=" * 60)
    print("Example 5: Custom Ayanamsha")
    print("=" * 60)

    from libephemeris.constants import SIDM_USER

    jd = eph.julday(2024, 1, 1, 0.0)

    print("\nYou can define a custom ayanamsha by specifying:")
    print("- Reference epoch (Julian Day)")
    print("- Ayanamsha value at that epoch")
    print("- (Optionally) precession rate\n")

    # Set a custom ayanamsha: 23.5° at J2000
    j2000 = 2451545.0  # J2000.0 epoch
    custom_aya = 23.5  # Custom value at J2000

    # SIDM_USER (255) enables user-defined ayanamsha
    eph.set_sid_mode(SIDM_USER, j2000, custom_aya)

    aya = eph.get_ayanamsa_ut(jd)
    print("Custom ayanamsha (23.5° at J2000):")
    print(f"  Value at 2024-01-01: {aya:.6f}°")

    # Compare with standard Lahiri
    eph.set_sid_mode(SIDM_LAHIRI)
    lahiri = eph.get_ayanamsa_ut(jd)
    print("\nLahiri for comparison:")
    print(f"  Value at 2024-01-01: {lahiri:.6f}°")


def example_nakshatra_calculation() -> None:
    """Example 6: Calculate Nakshatra position."""
    print("\n" + "=" * 60)
    print("Example 6: Nakshatra Calculation")
    print("=" * 60)

    # Nakshatras (27 lunar mansions of 13°20' each)
    NAKSHATRAS = [
        "Ashwini",
        "Bharani",
        "Krittika",
        "Rohini",
        "Mrigashira",
        "Ardra",
        "Punarvasu",
        "Pushya",
        "Ashlesha",
        "Magha",
        "Purva Phalguni",
        "Uttara Phalguni",
        "Hasta",
        "Chitra",
        "Swati",
        "Vishakha",
        "Anuradha",
        "Jyeshtha",
        "Mula",
        "Purva Ashadha",
        "Uttara Ashadha",
        "Shravana",
        "Dhanishta",
        "Shatabhisha",
        "Purva Bhadrapada",
        "Uttara Bhadrapada",
        "Revati",
    ]

    jd = eph.julday(2024, 1, 1, 12.0)

    # Set Lahiri ayanamsha
    eph.set_sid_mode(SIDM_LAHIRI)

    print(f"\nNakshatra positions at JD {jd}:\n")

    flags = FLG_SWIEPH | FLG_SIDEREAL

    planets = [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MARS, "Mars"),
        (MERCURY, "Mercury"),
        (JUPITER, "Jupiter"),
        (VENUS, "Venus"),
        (SATURN, "Saturn"),
    ]

    print(f"{'Planet':<10} {'Nakshatra':<20} {'Pada':>6} {'Position':>12}")
    print("-" * 52)

    for planet_id, name in planets:
        pos, _ = eph.calc_ut(jd, planet_id, flags)
        longitude = pos[0]

        # Calculate nakshatra (each is 13°20' = 13.3333°)
        nakshatra_size = 360 / 27  # 13.3333...
        nakshatra_num = int(longitude / nakshatra_size)
        nakshatra_name = NAKSHATRAS[nakshatra_num]

        # Calculate pada (each nakshatra has 4 padas of 3°20' each)
        pos_in_nakshatra = longitude - (nakshatra_num * nakshatra_size)
        pada = int(pos_in_nakshatra / (nakshatra_size / 4)) + 1

        print(f"{name:<10} {nakshatra_name:<20} {pada:>6} {longitude:>11.4f}°")


def main() -> None:
    """Run all sidereal calculation examples."""
    print("\n" + "#" * 60)
    print("# libephemeris Sidereal Calculations Examples")
    print("#" * 60)

    example_ayanamsha_overview()
    example_tropical_vs_sidereal()
    example_vedic_chart()
    example_ayanamsha_change_over_time()
    example_custom_ayanamsha()
    example_nakshatra_calculation()

    print("\n" + "#" * 60)
    print("# Sidereal examples completed!")
    print("#" * 60 + "\n")


if __name__ == "__main__":
    main()
