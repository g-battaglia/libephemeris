#!/usr/bin/env python3
"""
Thread-Safe Context Example for libephemeris.

This script demonstrates how to use EphemerisContext for:
1. Thread-safe calculations
2. Independent sidereal modes per context
3. Independent topocentric settings per context
4. Parallel calculations in different threads

The EphemerisContext class provides isolated state for ephemeris calculations,
allowing concurrent use in multi-threaded applications without global state
interference.

Requirements:
    pip install libephemeris

Usage:
    python examples/thread_safe_context.py
"""

from __future__ import annotations

import concurrent.futures
import threading
from typing import Any

import libephemeris as eph
from libephemeris import EphemerisContext
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_SIDEREAL,
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
)


def example_basic_context() -> None:
    """Example 1: Basic EphemerisContext usage."""
    print("=" * 60)
    print("Example 1: Basic EphemerisContext Usage")
    print("=" * 60)

    print("\nEphemerisContext provides isolated state for calculations.\n")

    # Create a context
    ctx = EphemerisContext()

    jd = eph.julday(2024, 1, 1, 12.0)

    # Calculate using context methods
    pos, flags = ctx.calc_ut(jd, SUN, FLG_SWIEPH | FLG_SPEED)

    print(f"Sun position using context: {pos[0]:.6f}°")

    # Compare with global function
    pos_global, _ = eph.calc_ut(jd, SUN, FLG_SWIEPH | FLG_SPEED)

    print(f"Sun position using global:  {pos_global[0]:.6f}°")
    print(f"Difference: {abs(pos[0] - pos_global[0]):.10f}°")


def example_independent_sidereal_modes() -> None:
    """Example 2: Different sidereal modes in different contexts."""
    print("\n" + "=" * 60)
    print("Example 2: Independent Sidereal Modes")
    print("=" * 60)

    print("\nEach context can have its own sidereal mode setting.\n")

    jd = eph.julday(2024, 1, 1, 12.0)

    # Context 1: Lahiri ayanamsha
    ctx_lahiri = EphemerisContext()
    ctx_lahiri.set_sid_mode(SIDM_LAHIRI)

    # Context 2: Fagan-Bradley ayanamsha
    ctx_fagan = EphemerisContext()
    ctx_fagan.set_sid_mode(SIDM_FAGAN_BRADLEY)

    # Note: get_ayanamsa_ut is not available on context, but sidereal calculations work
    print("Context 1: Lahiri mode set")
    print("Context 2: Fagan-Bradley mode set")

    # Calculate sidereal Sun position - each context uses its own mode
    pos_lahiri, _ = ctx_lahiri.calc_ut(jd, SUN, FLG_SWIEPH | FLG_SIDEREAL)
    pos_fagan, _ = ctx_fagan.calc_ut(jd, SUN, FLG_SWIEPH | FLG_SIDEREAL)

    print(f"\nSun sidereal position (Lahiri mode):       {pos_lahiri[0]:.6f}°")
    print(f"Sun sidereal position (Fagan-Bradley mode): {pos_fagan[0]:.6f}°")
    print(f"Difference: {abs(pos_lahiri[0] - pos_fagan[0]):.6f}°")


def example_topocentric_settings() -> None:
    """Example 3: Different topocentric settings per context."""
    print("\n" + "=" * 60)
    print("Example 3: Independent Topocentric Settings")
    print("=" * 60)

    print("\nEach context can have its own observer location.\n")

    from libephemeris.constants import FLG_TOPOCTR

    jd = eph.julday(2024, 1, 1, 12.0)

    # Context for New York
    ctx_ny = EphemerisContext()
    ctx_ny.set_topo(lon=-74.0060, lat=40.7128, alt=10)

    # Context for Sydney
    ctx_sydney = EphemerisContext()
    ctx_sydney.set_topo(lon=151.2093, lat=-33.8688, alt=58)

    # Calculate topocentric Moon position (most affected by parallax)
    pos_ny, _ = ctx_ny.calc_ut(jd, MOON, FLG_SWIEPH | FLG_TOPOCTR)
    pos_sydney, _ = ctx_sydney.calc_ut(jd, MOON, FLG_SWIEPH | FLG_TOPOCTR)

    # Geocentric for comparison
    pos_geo, _ = eph.calc_ut(jd, MOON, FLG_SWIEPH)

    print("Moon position (topocentric):")
    print(f"  From New York:   {pos_ny[0]:.6f}° lon, {pos_ny[1]:.6f}° lat")
    print(f"  From Sydney:     {pos_sydney[0]:.6f}° lon, {pos_sydney[1]:.6f}° lat")
    print(f"  Geocentric:      {pos_geo[0]:.6f}° lon, {pos_geo[1]:.6f}° lat")

    diff = abs(pos_ny[0] - pos_sydney[0])
    print(f"\nParallax difference: {diff:.4f}° ({diff * 60:.2f} arcmin)")


def worker_function(
    context: EphemerisContext, jd: float, planet: int, thread_name: str
) -> dict[str, Any]:
    """Worker function for threaded calculations."""
    # Each thread uses its own context
    pos, _ = context.calc_ut(jd, planet, FLG_SWIEPH | FLG_SPEED)

    return {
        "thread": thread_name,
        "planet": planet,
        "longitude": pos[0],
        "latitude": pos[1],
        "thread_id": threading.current_thread().name,
    }


def example_multithreaded_calculations() -> None:
    """Example 4: Multi-threaded calculations with separate contexts."""
    print("\n" + "=" * 60)
    print("Example 4: Multi-threaded Calculations")
    print("=" * 60)

    print("\nUsing separate contexts for concurrent calculations.\n")

    jd = eph.julday(2024, 1, 1, 12.0)

    # Create tasks with different planets and contexts
    tasks = [
        (EphemerisContext(), jd, SUN, "Sun"),
        (EphemerisContext(), jd, MOON, "Moon"),
        (EphemerisContext(), jd, MARS, "Mars"),
    ]

    results = []

    # Use ThreadPoolExecutor for parallel execution
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        futures = [
            executor.submit(worker_function, ctx, jd, planet, name)
            for ctx, jd, planet, name in tasks
        ]

        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())

    print(f"{'Planet':<10}{'Longitude':>12}{'Latitude':>12}{'Thread':>20}")
    print("-" * 56)

    for result in sorted(results, key=lambda x: x["planet"]):
        print(
            f"{result['thread']:<10}"
            f"{result['longitude']:>12.6f}"
            f"{result['latitude']:>12.6f}"
            f"{result['thread_id']:>20}"
        )


def example_context_isolation() -> None:
    """Example 5: Demonstrating context isolation."""
    print("\n" + "=" * 60)
    print("Example 5: Context Isolation")
    print("=" * 60)

    print("\nChanges in one context don't affect others.\n")

    jd = eph.julday(2024, 1, 1, 12.0)

    # Create two contexts
    ctx1 = EphemerisContext()
    ctx2 = EphemerisContext()

    # Set sidereal mode only in ctx1
    ctx1.set_sid_mode(SIDM_LAHIRI)

    # Calculate in both contexts with sidereal flag
    pos1, _ = ctx1.calc_ut(jd, SUN, FLG_SWIEPH | FLG_SIDEREAL)
    pos2, _ = ctx2.calc_ut(jd, SUN, FLG_SWIEPH | FLG_SIDEREAL)

    print("Context 1 (Lahiri set):")
    print(f"  Sun sidereal position: {pos1[0]:.6f}°")

    print("\nContext 2 (default mode):")
    print(f"  Sun sidereal position: {pos2[0]:.6f}°")

    # Note: ctx2 uses default sidereal mode (Fagan-Bradley)
    diff = abs(pos1[0] - pos2[0])
    print(f"\nDifference: {diff:.6f}° (confirms isolation)")


def example_house_calculations() -> None:
    """Example 6: House calculations with context."""
    print("\n" + "=" * 60)
    print("Example 6: House Calculations with Context")
    print("=" * 60)

    print("\nContexts can be used for house calculations too.\n")

    jd = eph.julday(2024, 1, 1, 12.0)
    latitude = 51.5074  # London
    longitude = -0.1278

    # Create context
    ctx = EphemerisContext()

    # Calculate houses using context
    cusps, ascmc = ctx.houses(jd, latitude, longitude, ord("P"))

    print(f"Location: London ({latitude}°N, {longitude}°W)")
    print("\nPlacidus house cusps (using context):")

    for i in range(6):
        print(
            f"  House {i + 1}: {cusps[i]:>10.4f}°    House {i + 7}: {cusps[i + 6]:>10.4f}°"
        )

    print(f"\nAscendant: {ascmc[0]:.4f}°")
    print(f"Midheaven: {ascmc[1]:.4f}°")


def example_sidereal_houses() -> None:
    """Example 7: Sidereal house calculations."""
    print("\n" + "=" * 60)
    print("Example 7: Sidereal Houses with Context")
    print("=" * 60)

    jd = eph.julday(2024, 1, 1, 12.0)
    latitude = 13.0827  # Chennai
    longitude = 80.2707

    print(f"\nLocation: Chennai ({latitude}°N, {longitude}°E)\n")

    # Context for sidereal calculations
    ctx = EphemerisContext()
    ctx.set_sid_mode(SIDM_LAHIRI)

    # Get ayanamsha using global function
    eph.set_sid_mode(SIDM_LAHIRI)
    aya = eph.get_ayanamsa_ut(jd)
    print(f"Lahiri Ayanamsha: {aya:.6f}°")

    # Calculate tropical houses using context
    cusps_trop, ascmc_trop = ctx.houses(jd, latitude, longitude, ord("W"))

    # Calculate sidereal houses using global function with sidereal flag
    cusps_sid, ascmc_sid = eph.houses_ex(
        jd, latitude, longitude, ord("W"), FLG_SIDEREAL
    )

    print("\nWhole Sign Houses (1st house cusp = Lagna):")
    print(f"  Tropical ASC:  {ascmc_trop[0]:>10.4f}°")
    print(f"  Sidereal ASC:  {ascmc_sid[0]:>10.4f}°")
    print(f"  Difference:    {ascmc_trop[0] - ascmc_sid[0]:>10.4f}°")


def example_context_reuse() -> None:
    """Example 8: Reusing a context for multiple calculations."""
    print("\n" + "=" * 60)
    print("Example 8: Reusing Context for Multiple Calculations")
    print("=" * 60)

    print("\nA single context can be reused for multiple calculations.\n")

    ctx = EphemerisContext()

    # Set up context once
    ctx.set_sid_mode(SIDM_LAHIRI)
    ctx.set_topo(lon=77.2090, lat=28.6139, alt=216)  # New Delhi

    print("Calculating planets for a full year (monthly samples):")
    print("-" * 50)

    planets = [(SUN, "Sun"), (MOON, "Moon"), (MARS, "Mars")]

    for month in [1, 4, 7, 10]:  # Quarterly samples
        jd = eph.julday(2024, month, 1, 12.0)
        print(f"\n2024-{month:02d}-01:")

        for planet_id, name in planets:
            pos, _ = ctx.calc_ut(jd, planet_id, FLG_SWIEPH | FLG_SIDEREAL)
            print(f"  {name:8} {pos[0]:>10.4f}° (sidereal)")


def main() -> None:
    """Run all thread-safe context examples."""
    print("\n" + "#" * 60)
    print("# libephemeris Thread-Safe Context Examples")
    print("#" * 60)

    example_basic_context()
    example_independent_sidereal_modes()
    example_topocentric_settings()
    example_multithreaded_calculations()
    example_context_isolation()
    example_house_calculations()
    example_sidereal_houses()
    example_context_reuse()

    print("\n" + "#" * 60)
    print("# Thread-safe context examples completed!")
    print("#" * 60 + "\n")


if __name__ == "__main__":
    main()
