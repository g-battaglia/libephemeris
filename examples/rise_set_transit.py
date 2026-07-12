#!/usr/bin/env python3
"""
Rise, Set, and Transit Example for libephemeris.

This script demonstrates how to calculate:
1. Sunrise and sunset times
2. Moonrise and moonset
3. Planet rise/set/transit times
4. Twilight times (civil, nautical, astronomical)
5. Meridian transits

Requirements:
    pip install libephemeris

Usage:
    python examples/rise_set_transit.py
"""

from __future__ import annotations

import libephemeris as eph
from libephemeris.constants import (
    # Planet IDs
    SUN,
    MOON,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    # Rise/Set/Transit event types
    CALC_RISE,
    CALC_SET,
    CALC_MTRANSIT,
    BIT_DISC_CENTER,
    BIT_NO_REFRACTION,
    BIT_CIVIL_TWILIGHT,
    BIT_NAUTIC_TWILIGHT,
    BIT_ASTRO_TWILIGHT,
)


def jd_to_time_string(jd: float) -> str:
    """Convert Julian Day to time string (HH:MM)."""
    year, month, day, hour = eph.revjul(jd)
    hours = int(hour)
    minutes = int((hour - hours) * 60)
    return f"{hours:02d}:{minutes:02d}"


def jd_to_datetime_string(jd: float) -> str:
    """Convert Julian Day to datetime string."""
    year, month, day, hour = eph.revjul(jd)
    hours = int(hour)
    minutes = int((hour - hours) * 60)
    return f"{year}-{month:02d}-{day:02d} {hours:02d}:{minutes:02d}"


def example_sunrise_sunset() -> None:
    """Example 1: Calculate sunrise and sunset times."""
    print("=" * 60)
    print("Example 1: Sunrise and Sunset")
    print("=" * 60)

    # Location: San Francisco, CA
    latitude = 37.7749
    longitude = -122.4194
    altitude = 16.0  # meters above sea level

    # Date: Summer solstice 2024
    jd_start = eph.julday(2024, 6, 21, 0.0)

    print(f"\nLocation: San Francisco ({latitude}°N, {longitude}°W)")
    print("Date: 2024-06-21 (Summer Solstice)\n")

    # Calculate sunrise
    flag_rise, tret = eph.rise_trans(
        jd_start, SUN, CALC_RISE, [longitude, latitude, altitude]
    )
    jd_rise = tret[0]

    # Calculate sunset
    flag_set, tret = eph.rise_trans(
        jd_start, SUN, CALC_SET, [longitude, latitude, altitude]
    )
    jd_set = tret[0]

    # Calculate upper meridian transit (solar noon)
    flag_transit, tret = eph.rise_trans(
        jd_start, SUN, CALC_MTRANSIT, [longitude, latitude, altitude]
    )
    jd_transit = tret[0]

    print(f"Sunrise:    {jd_to_time_string(jd_rise)} UT")
    print(f"Solar Noon: {jd_to_time_string(jd_transit)} UT")
    print(f"Sunset:     {jd_to_time_string(jd_set)} UT")

    # Calculate day length
    if jd_set > jd_rise:
        day_length = (jd_set - jd_rise) * 24
        hours = int(day_length)
        minutes = int((day_length - hours) * 60)
        print(f"\nDay length: {hours}h {minutes}m")


def example_twilight_times() -> None:
    """Example 2: Calculate twilight times."""
    print("\n" + "=" * 60)
    print("Example 2: Twilight Times")
    print("=" * 60)

    # Location: London, UK
    latitude = 51.5074
    longitude = -0.1278

    jd_start = eph.julday(2024, 3, 20, 0.0)  # Spring equinox

    print(f"\nLocation: London ({latitude}°N, {longitude}°W)")
    print("Date: 2024-03-20 (Spring Equinox)\n")

    print("Morning Twilight:")
    print("-" * 40)

    # Astronomical dawn (Sun at -18°)
    _, tret = eph.rise_trans(
        jd_start,
        SUN,
        CALC_RISE | BIT_ASTRO_TWILIGHT,
        [longitude, latitude, 0.0],
    )
    jd_astro_dawn = tret[0]
    print(f"  Astronomical dawn (-18°): {jd_to_time_string(jd_astro_dawn)} UT")

    # Nautical dawn (Sun at -12°)
    _, tret = eph.rise_trans(
        jd_start,
        SUN,
        CALC_RISE | BIT_NAUTIC_TWILIGHT,
        [longitude, latitude, 0.0],
    )
    jd_naut_dawn = tret[0]
    print(f"  Nautical dawn (-12°):     {jd_to_time_string(jd_naut_dawn)} UT")

    # Civil dawn (Sun at -6°)
    _, tret = eph.rise_trans(
        jd_start,
        SUN,
        CALC_RISE | BIT_CIVIL_TWILIGHT,
        [longitude, latitude, 0.0],
    )
    jd_civil_dawn = tret[0]
    print(f"  Civil dawn (-6°):         {jd_to_time_string(jd_civil_dawn)} UT")

    # Sunrise
    _, tret = eph.rise_trans(jd_start, SUN, CALC_RISE, [longitude, latitude, 0.0])
    jd_sunrise = tret[0]
    print(f"  Sunrise (0°):             {jd_to_time_string(jd_sunrise)} UT")

    print("\nEvening Twilight:")
    print("-" * 40)

    # Sunset
    _, tret = eph.rise_trans(jd_start, SUN, CALC_SET, [longitude, latitude, 0.0])
    jd_sunset = tret[0]
    print(f"  Sunset (0°):              {jd_to_time_string(jd_sunset)} UT")

    # Civil dusk
    _, tret = eph.rise_trans(
        jd_start,
        SUN,
        CALC_SET | BIT_CIVIL_TWILIGHT,
        [longitude, latitude, 0.0],
    )
    jd_civil_dusk = tret[0]
    print(f"  Civil dusk (-6°):         {jd_to_time_string(jd_civil_dusk)} UT")

    # Nautical dusk
    _, tret = eph.rise_trans(
        jd_start,
        SUN,
        CALC_SET | BIT_NAUTIC_TWILIGHT,
        [longitude, latitude, 0.0],
    )
    jd_naut_dusk = tret[0]
    print(f"  Nautical dusk (-12°):     {jd_to_time_string(jd_naut_dusk)} UT")

    # Astronomical dusk
    _, tret = eph.rise_trans(
        jd_start,
        SUN,
        CALC_SET | BIT_ASTRO_TWILIGHT,
        [longitude, latitude, 0.0],
    )
    jd_astro_dusk = tret[0]
    print(f"  Astronomical dusk (-18°): {jd_to_time_string(jd_astro_dusk)} UT")


def example_moonrise_moonset() -> None:
    """Example 3: Calculate moonrise and moonset."""
    print("\n" + "=" * 60)
    print("Example 3: Moonrise and Moonset")
    print("=" * 60)

    # Location: Tokyo, Japan
    latitude = 35.6762
    longitude = 139.6503

    jd_start = eph.julday(2024, 4, 8, 0.0)  # Total solar eclipse date

    print(f"\nLocation: Tokyo ({latitude}°N, {longitude}°E)")
    print("Date: 2024-04-08\n")

    # Calculate moonrise
    flag, tret = eph.rise_trans(jd_start, MOON, CALC_RISE, [longitude, latitude, 0.0])
    jd_moonrise = tret[0]

    # Calculate moonset
    _, tret = eph.rise_trans(jd_start, MOON, CALC_SET, [longitude, latitude, 0.0])
    jd_moonset = tret[0]

    # Calculate moon transit
    _, tret = eph.rise_trans(jd_start, MOON, CALC_MTRANSIT, [longitude, latitude, 0.0])
    jd_transit = tret[0]

    if flag != -2:  # -2 means circumpolar (never rises/sets)
        print(f"Moonrise:  {jd_to_time_string(jd_moonrise)} UT")
        print(f"Transit:   {jd_to_time_string(jd_transit)} UT")
        print(f"Moonset:   {jd_to_time_string(jd_moonset)} UT")
    else:
        print("Moon is circumpolar (doesn't rise or set)")


def example_planet_visibility() -> None:
    """Example 4: Calculate planet rise/set times for visibility planning."""
    print("\n" + "=" * 60)
    print("Example 4: Planet Visibility Times")
    print("=" * 60)

    # Location: Sydney, Australia
    latitude = -33.8688
    longitude = 151.2093

    jd_start = eph.julday(2024, 7, 15, 0.0)

    print(f"\nLocation: Sydney ({abs(latitude)}°S, {longitude}°E)")
    print("Date: 2024-07-15\n")

    planets = [
        (VENUS, "Venus"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
    ]

    print(f"{'Planet':<10}{'Rise':>10}{'Transit':>10}{'Set':>10}")
    print("-" * 42)

    for planet_id, name in planets:
        try:
            flag_rise, tret = eph.rise_trans(
                jd_start, planet_id, CALC_RISE, [longitude, latitude, 0.0]
            )
            jd_rise = tret[0]

            _, tret = eph.rise_trans(
                jd_start, planet_id, CALC_MTRANSIT, [longitude, latitude, 0.0]
            )
            jd_transit = tret[0]

            _, tret = eph.rise_trans(
                jd_start, planet_id, CALC_SET, [longitude, latitude, 0.0]
            )
            jd_set = tret[0]

            if flag_rise != -2:
                rise_str = jd_to_time_string(jd_rise)
                transit_str = jd_to_time_string(jd_transit)
                set_str = jd_to_time_string(jd_set)
                print(f"{name:<10}{rise_str:>10}{transit_str:>10}{set_str:>10}")
            else:
                print(f"{name:<10}{'(circumpolar)':^30}")
        except Exception as e:
            print(f"{name:<10} Error: {e}")


def example_week_of_sunrise() -> None:
    """Example 5: Sunrise times for a week."""
    print("\n" + "=" * 60)
    print("Example 5: Weekly Sunrise/Sunset Table")
    print("=" * 60)

    # Location: New York City
    latitude = 40.7128
    longitude = -74.0060

    # Start date
    start_year, start_month, start_day = 2024, 6, 17

    print(f"\nLocation: New York City ({latitude}°N, {longitude}°W)")
    print(f"Week of: {start_year}-{start_month:02d}-{start_day:02d}\n")

    print(f"{'Date':<12}{'Sunrise':>10}{'Sunset':>10}{'Day Length':>14}")
    print("-" * 48)

    for day_offset in range(7):
        jd = eph.julday(start_year, start_month, start_day + day_offset, 0.0)
        year, month, day, _ = eph.revjul(jd)

        _, tret = eph.rise_trans(jd, SUN, CALC_RISE, [longitude, latitude, 0.0])
        jd_rise = tret[0]
        _, tret = eph.rise_trans(jd, SUN, CALC_SET, [longitude, latitude, 0.0])
        jd_set = tret[0]

        rise_str = jd_to_time_string(jd_rise)
        set_str = jd_to_time_string(jd_set)

        # Day length
        day_length = (jd_set - jd_rise) * 24
        hours = int(day_length)
        minutes = int((day_length - hours) * 60)
        length_str = f"{hours}h {minutes}m"

        date_str = f"{year}-{month:02d}-{day:02d}"
        print(f"{date_str:<12}{rise_str:>10}{set_str:>10}{length_str:>14}")


def example_no_refraction() -> None:
    """Example 6: Compare with and without atmospheric refraction."""
    print("\n" + "=" * 60)
    print("Example 6: Effect of Atmospheric Refraction")
    print("=" * 60)

    # Location: Denver (high altitude)
    latitude = 39.7392
    longitude = -104.9903
    altitude = 1609.0  # Mile high city

    jd_start = eph.julday(2024, 1, 1, 0.0)

    print(f"\nLocation: Denver ({latitude}°N, {longitude}°W, {altitude}m altitude)")
    print("Date: 2024-01-01\n")

    print("Atmospheric refraction bends light, making the Sun appear")
    print("slightly higher than it actually is, affecting rise/set times.\n")

    # With refraction (default)
    _, tret = eph.rise_trans(jd_start, SUN, CALC_RISE, [longitude, latitude, altitude])
    jd_rise_refr = tret[0]
    _, tret = eph.rise_trans(jd_start, SUN, CALC_SET, [longitude, latitude, altitude])
    jd_set_refr = tret[0]

    # Without refraction
    _, tret = eph.rise_trans(
        jd_start,
        SUN,
        CALC_RISE | BIT_NO_REFRACTION,
        [longitude, latitude, altitude],
    )
    jd_rise_no_refr = tret[0]
    _, tret = eph.rise_trans(
        jd_start,
        SUN,
        CALC_SET | BIT_NO_REFRACTION,
        [longitude, latitude, altitude],
    )
    jd_set_no_refr = tret[0]

    print(f"{'Condition':<24}{'Sunrise':>12}{'Sunset':>12}")
    print("-" * 50)
    print(
        f"{'With refraction':<24}{jd_to_time_string(jd_rise_refr):>12}{jd_to_time_string(jd_set_refr):>12}"
    )
    print(
        f"{'Without refraction':<24}{jd_to_time_string(jd_rise_no_refr):>12}{jd_to_time_string(jd_set_no_refr):>12}"
    )

    # Difference in minutes
    rise_diff = (jd_rise_no_refr - jd_rise_refr) * 24 * 60
    set_diff = (jd_set_refr - jd_set_no_refr) * 24 * 60
    print(f"\nRefraction makes sunrise ~{abs(rise_diff):.1f} min earlier")
    print(f"and sunset ~{abs(set_diff):.1f} min later")


def example_disc_center_vs_limb() -> None:
    """Example 7: Disc center vs upper limb for Sun/Moon."""
    print("\n" + "=" * 60)
    print("Example 7: Disc Center vs Upper Limb")
    print("=" * 60)

    latitude = 45.0
    longitude = 0.0

    jd_start = eph.julday(2024, 3, 20, 0.0)

    print("\nThe Sun and Moon are not point sources. By default,")
    print("rise/set is calculated for the upper limb (edge) of the disc.\n")

    # Upper limb (default)
    _, tret = eph.rise_trans(jd_start, SUN, CALC_RISE, [longitude, latitude, 0.0])
    jd_rise_limb = tret[0]

    # Disc center
    _, tret = eph.rise_trans(
        jd_start, SUN, CALC_RISE | BIT_DISC_CENTER, [longitude, latitude, 0.0]
    )
    jd_rise_center = tret[0]

    print(f"{'Method':<20}{'Sunrise':>15}")
    print("-" * 38)
    print(f"{'Upper limb (default)':<20}{jd_to_time_string(jd_rise_limb):>15}")
    print(f"{'Disc center':<20}{jd_to_time_string(jd_rise_center):>15}")

    diff = (jd_rise_center - jd_rise_limb) * 24 * 60
    print(f"\nDifference: ~{abs(diff):.1f} minutes")
    print("(The Sun's angular diameter is about 32 arcminutes)")


def main() -> None:
    """Run all rise/set/transit examples."""
    print("\n" + "#" * 60)
    print("# libephemeris Rise/Set/Transit Examples")
    print("#" * 60)

    example_sunrise_sunset()
    example_twilight_times()
    example_moonrise_moonset()
    example_planet_visibility()
    example_week_of_sunrise()
    example_no_refraction()
    example_disc_center_vs_limb()

    print("\n" + "#" * 60)
    print("# Rise/Set/Transit examples completed!")
    print("#" * 60 + "\n")


if __name__ == "__main__":
    main()
