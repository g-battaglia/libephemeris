#!/usr/bin/env python3
"""Verification: LEB vs Skyfield, Uranians, Asteroids, Orbital Elements, Crossings/Stations

Sections covered:
  1.3  LEB vs Skyfield consistency (14 bodies x 100 dates x 3 checks = 4,200)
  1.4  LEB flag modes (14 bodies x 50 dates x 4 modes x 3 checks = 8,400)
  13   Uranians deep (geocentric 9x50x3 + helio 9x30x3 + sidereal 8x3x10x3 = 2,880)
  14   Asteroids vs pyswisseph (5x100x3 + AST_OFFSET mapping = 1,503)
  16   Orbital elements (5x10x6 = 300)
  18   Crossings and stations (equinoxes + solstices + moon crossings + stations ~ 110)

Target: ~17,000+ checks, <30s

Tolerances:
  - Uranians vs pyswisseph: 60" (known difference: different orbital element sources
    for hypothetical Hamburg School bodies; libephemeris uses Skyfield + analytical
    elements while pyswisseph uses Swiss Ephemeris internal tables)
  - Asteroids vs pyswisseph: 6" (Pallas/Juno can differ up to ~5.5")
  - LEB vs Skyfield: 0.005" (same computation path, only Chebyshev approx error)
"""

from __future__ import annotations

import math
import os
import sys
import time
import traceback
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import libephemeris as lib
import swisseph as swe_ref

# Point pyswisseph to its data files
swe_ref.set_ephe_path(_REF_EPHE_PATH)

# ---------------------------------------------------------------------------
# Counters and helpers
# ---------------------------------------------------------------------------
passed = 0
failed = 0
errors: list[str] = []
section_stats: dict[str, dict[str, int]] = {}
# Track max differences per (section, body, coord) for diagnostics
max_diffs: dict[str, float] = defaultdict(float)

_current_section = ""


def set_section(name: str) -> None:
    global _current_section
    _current_section = name
    if name not in section_stats:
        section_stats[name] = {"passed": 0, "failed": 0}


def check(condition: bool, description: str = "") -> None:
    global passed, failed
    if condition:
        passed += 1
        section_stats[_current_section]["passed"] += 1
    else:
        failed += 1
        section_stats[_current_section]["failed"] += 1
        if len(errors) < 100:
            errors.append(f"[{_current_section}] {description}")


def track_max(key: str, value: float) -> None:
    """Track maximum diff for a given key (used in diagnostics)."""
    if value > max_diffs[key]:
        max_diffs[key] = value


def angular_diff(a: float, b: float) -> float:
    """Smallest angular difference in degrees."""
    d = (a - b) % 360.0
    if d > 180.0:
        d = 360.0 - d
    return d


def arcsec(deg: float) -> float:
    """Convert degrees to arcseconds."""
    return abs(deg) * 3600.0


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SEED = 42
JD_MIN = 2415020.5  # ~1900
JD_MAX = 2488069.5  # ~2100
JD_J2000 = 2451545.0  # J2000.0

FLG_SPEED = 256
FLG_SWIEPH = 2
FLG_HELCTR = 8
FLG_J2000 = 32
FLG_NOABERR = 1024
FLG_EQUATORIAL = 2048
FLG_SIDEREAL = 65536

# Bodies for Sec 1.3 / 1.4: Sun(0)-Pluto(9), MeanNode(10), TrueNode(11),
# MeanApog(12), Chiron(15)
LEB_BODIES = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15]
LEB_BODY_NAMES = {
    0: "Sun",
    1: "Moon",
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "Pluto",
    10: "MeanNode",
    11: "TrueNode",
    12: "MeanApog",
    15: "Chiron",
}

# Uranian bodies (Sec 13)
URANIAN_BODIES = [40, 41, 42, 43, 44, 45, 46, 47, 48]
URANIAN_NAMES = {
    40: "Cupido",
    41: "Hades",
    42: "Zeus",
    43: "Kronos",
    44: "Apollon",
    45: "Admetos",
    46: "Vulkanus",
    47: "Poseidon",
    48: "Transpluto",
}

# Asteroids (Sec 14)
ASTEROID_BODIES = [15, 17, 18, 19, 20]
ASTEROID_NAMES = {
    15: "Chiron",
    17: "Ceres",
    18: "Pallas",
    19: "Juno",
    20: "Vesta",
}

# Orbital elements bodies (Sec 16)
ORBITAL_BODIES = [4, 5, 6, 7, 8]  # Mars, Jupiter, Saturn, Uranus, Neptune
ORBITAL_NAMES = {
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
}

# Known semi-major axes (AU) for sanity checks
KNOWN_SEMI_MAJOR = {
    4: 1.524,  # Mars
    5: 5.203,  # Jupiter
    6: 9.537,  # Saturn
    7: 19.19,  # Uranus
    8: 30.07,  # Neptune
}

# Sidereal ayanamsha modes for Sec 13
AYANAMSHA_MODES = [0, 1, 3]  # Fagan-Bradley, Lahiri, Raman

rng = np.random.default_rng(SEED)

# ---------------------------------------------------------------------------
# Date generation helpers
# ---------------------------------------------------------------------------


def random_jds(n: int, jd_min: float = JD_MIN, jd_max: float = JD_MAX) -> list[float]:
    return [float(x) for x in rng.uniform(jd_min, jd_max, n)]


# ===================================================================
# SECTION 1.3: LEB vs Skyfield consistency
# ===================================================================
def run_section_1_3() -> None:
    set_section("Sec1.3 LEB vs Skyfield")
    print(f"\n{'=' * 70}")
    print("SECTION 1.3: LEB vs Skyfield consistency")
    print("  14 bodies x 100 dates x 3 checks = 4,200 checks")
    print(f"{'=' * 70}")

    # Check if LEB is active
    reader = lib.get_leb_reader()
    if reader is None:
        print("  SKIP: No LEB reader configured, skipping LEB vs Skyfield tests")
        return

    N_DATES = 100
    LON_TOL = 0.005  # arcsec
    LAT_TOL = 0.005  # arcsec
    DIST_TOL = 1e-7  # AU

    jds = random_jds(N_DATES)
    flags = FLG_SPEED | FLG_SWIEPH

    for body_id in LEB_BODIES:
        bname = LEB_BODY_NAMES[body_id]

        for jd in jds:
            # Skyfield result
            try:
                lib.set_calc_mode("skyfield")
                sky_res, _ = lib.calc_ut(jd, body_id, flags)
            except Exception as e:
                check(False, f"{bname} jd={jd:.2f} Skyfield crash: {e}")
                continue
            finally:
                lib.set_calc_mode("auto")

            # LEB result (auto mode with LEB configured)
            try:
                leb_res, _ = lib.calc_ut(jd, body_id, flags)
            except Exception as e:
                check(False, f"{bname} jd={jd:.2f} LEB crash: {e}")
                continue

            lon_diff = arcsec(angular_diff(leb_res[0], sky_res[0]))
            lat_diff = arcsec(abs(leb_res[1] - sky_res[1]))
            dist_diff = abs(leb_res[2] - sky_res[2])

            track_max(f"1.3|{bname}|lon", lon_diff)
            track_max(f"1.3|{bname}|lat", lat_diff)
            track_max(f"1.3|{bname}|dist", dist_diff)

            check(
                lon_diff < LON_TOL,
                f'{bname} lon diff={lon_diff:.6f}" > {LON_TOL}" at jd={jd:.2f}',
            )
            check(
                lat_diff < LAT_TOL,
                f'{bname} lat diff={lat_diff:.6f}" > {LAT_TOL}" at jd={jd:.2f}',
            )
            check(
                dist_diff < DIST_TOL,
                f"{bname} dist diff={dist_diff:.2e} > {DIST_TOL} at jd={jd:.2f}",
            )

    print(f"  Done: {section_stats[_current_section]}")


# ===================================================================
# SECTION 1.4: LEB flag modes
# ===================================================================
def run_section_1_4() -> None:
    set_section("Sec1.4 LEB flag modes")
    print(f"\n{'=' * 70}")
    print("SECTION 1.4: LEB flag modes")
    print("  14 bodies x 50 dates x 4 modes x 3 checks = 8,400 checks")
    print(f"{'=' * 70}")

    reader = lib.get_leb_reader()
    if reader is None:
        print("  SKIP: No LEB reader configured")
        return

    N_DATES = 50
    TOL_ARCSEC = 0.01
    DIST_TOL = 1e-7

    jds = random_jds(N_DATES)

    # Four flag modes: SIDEREAL(Lahiri), EQUATORIAL, J2000, NOABERR
    flag_modes = [
        ("SIDEREAL(Lahiri)", FLG_SPEED | FLG_SWIEPH | FLG_SIDEREAL, True),
        ("EQUATORIAL", FLG_SPEED | FLG_SWIEPH | FLG_EQUATORIAL, False),
        ("J2000", FLG_SPEED | FLG_SWIEPH | FLG_J2000, False),
        ("NOABERR", FLG_SPEED | FLG_SWIEPH | FLG_NOABERR, False),
    ]

    for mode_name, flags, needs_sid in flag_modes:
        if needs_sid:
            lib.set_sid_mode(1)  # Lahiri

        for body_id in LEB_BODIES:
            bname = LEB_BODY_NAMES[body_id]

            for jd in jds:
                # Skyfield
                try:
                    lib.set_calc_mode("skyfield")
                    sky_res, _ = lib.calc_ut(jd, body_id, flags)
                except Exception as e:
                    check(False, f"{mode_name} {bname} Skyfield crash: {e}")
                    continue
                finally:
                    lib.set_calc_mode("auto")

                # LEB (auto)
                try:
                    leb_res, _ = lib.calc_ut(jd, body_id, flags)
                except Exception as e:
                    check(False, f"{mode_name} {bname} LEB crash: {e}")
                    continue

                lon_diff = arcsec(angular_diff(leb_res[0], sky_res[0]))
                lat_diff = arcsec(abs(leb_res[1] - sky_res[1]))
                dist_diff = abs(leb_res[2] - sky_res[2])

                track_max(f"1.4|{mode_name}|{bname}|lon", lon_diff)
                track_max(f"1.4|{mode_name}|{bname}|lat", lat_diff)

                check(
                    lon_diff < TOL_ARCSEC,
                    f'{mode_name} {bname} lon={lon_diff:.6f}" at jd={jd:.2f}',
                )
                check(
                    lat_diff < TOL_ARCSEC,
                    f'{mode_name} {bname} lat={lat_diff:.6f}" at jd={jd:.2f}',
                )
                check(
                    dist_diff < DIST_TOL,
                    f"{mode_name} {bname} dist={dist_diff:.2e} at jd={jd:.2f}",
                )

        # Reset sidereal mode
        if needs_sid:
            lib.set_sid_mode(0)

    print(f"  Done: {section_stats[_current_section]}")


# ===================================================================
# SECTION 13: Uranians deep
# ===================================================================
def run_section_13() -> None:
    set_section("Sec13 Uranians deep")
    print(f"\n{'=' * 70}")
    print("SECTION 13: Uranians deep (vs pyswisseph)")
    print("  Geocentric: 9x50x3 = 1,350")
    print("  Heliocentric: 9x30x3 = 810")
    print("  Sidereal: 8x3x10x3 = 720")
    print("  Total = 2,880 checks")
    print(f"{'=' * 70}")

    # Uranians are hypothetical bodies (Hamburg School). libephemeris and
    # pyswisseph use slightly different orbital element sources and
    # geocentric conversion pipelines, so differences up to ~40" in
    # longitude are expected. We use 60" as the pass threshold.
    GEO_LON_TOL = 60.0  # arcsec
    GEO_LAT_TOL = 60.0  # arcsec
    HELIO_LON_TOL = 60.0  # arcsec
    HELIO_LAT_TOL = 60.0  # arcsec
    SID_LON_TOL = 60.0  # arcsec
    SID_LAT_TOL = 60.0  # arcsec
    DIST_TOL = 0.01  # AU

    # --- 13a: Geocentric ---
    N_GEO = 50
    jds_geo = random_jds(N_GEO)
    flags_geo = FLG_SPEED | FLG_SWIEPH

    for body_id in URANIAN_BODIES:
        bname = URANIAN_NAMES[body_id]
        for jd in jds_geo:
            try:
                lib_res, _ = lib.calc_ut(jd, body_id, flags_geo)
            except Exception as e:
                check(False, f"Geo {bname} lib crash jd={jd:.2f}: {e}")
                continue
            try:
                ref_res, _ = swe_ref.calc_ut(jd, body_id, flags_geo)
            except Exception as e:
                check(False, f"Geo {bname} ref crash jd={jd:.2f}: {e}")
                continue

            lon_diff = arcsec(angular_diff(lib_res[0], ref_res[0]))
            lat_diff = arcsec(abs(lib_res[1] - ref_res[1]))
            dist_diff = abs(lib_res[2] - ref_res[2])

            track_max(f"13geo|{bname}|lon", lon_diff)
            track_max(f"13geo|{bname}|lat", lat_diff)
            track_max(f"13geo|{bname}|dist", dist_diff)

            check(
                lon_diff < GEO_LON_TOL,
                f'Geo {bname} lon={lon_diff:.4f}" > {GEO_LON_TOL}" jd={jd:.2f}',
            )
            check(
                lat_diff < GEO_LAT_TOL,
                f'Geo {bname} lat={lat_diff:.4f}" > {GEO_LAT_TOL}" jd={jd:.2f}',
            )
            check(
                dist_diff < DIST_TOL,
                f"Geo {bname} dist={dist_diff:.6f} > {DIST_TOL} jd={jd:.2f}",
            )

    # --- 13b: Heliocentric ---
    N_HELIO = 30
    jds_helio = random_jds(N_HELIO)
    flags_helio = FLG_SPEED | FLG_SWIEPH | FLG_HELCTR

    for body_id in URANIAN_BODIES:
        bname = URANIAN_NAMES[body_id]
        for jd in jds_helio:
            try:
                lib_res, _ = lib.calc_ut(jd, body_id, flags_helio)
            except Exception as e:
                check(False, f"Helio {bname} lib crash jd={jd:.2f}: {e}")
                continue
            try:
                ref_res, _ = swe_ref.calc_ut(jd, body_id, flags_helio)
            except Exception as e:
                check(False, f"Helio {bname} ref crash jd={jd:.2f}: {e}")
                continue

            lon_diff = arcsec(angular_diff(lib_res[0], ref_res[0]))
            lat_diff = arcsec(abs(lib_res[1] - ref_res[1]))
            dist_diff = abs(lib_res[2] - ref_res[2])

            track_max(f"13helio|{bname}|lon", lon_diff)
            track_max(f"13helio|{bname}|lat", lat_diff)
            track_max(f"13helio|{bname}|dist", dist_diff)

            check(
                lon_diff < HELIO_LON_TOL,
                f'Helio {bname} lon={lon_diff:.4f}" > {HELIO_LON_TOL}" jd={jd:.2f}',
            )
            check(
                lat_diff < HELIO_LAT_TOL,
                f'Helio {bname} lat={lat_diff:.4f}" > {HELIO_LAT_TOL}" jd={jd:.2f}',
            )
            check(
                dist_diff < DIST_TOL,
                f"Helio {bname} dist={dist_diff:.6f} > {DIST_TOL} jd={jd:.2f}",
            )

    # --- 13c: Sidereal ---
    # Transpluto (48) excluded from sidereal test: same 8 Hamburg uranians
    SIDEREAL_URANIANS = [40, 41, 42, 43, 44, 45, 46, 47]
    N_SID = 10
    jds_sid = random_jds(N_SID)

    for ayan_mode in AYANAMSHA_MODES:
        ayan_name = {0: "FaganBradley", 1: "Lahiri", 3: "Raman"}[ayan_mode]

        lib.set_sid_mode(ayan_mode)
        swe_ref.set_sid_mode(ayan_mode)

        flags_sid = FLG_SPEED | FLG_SWIEPH | FLG_SIDEREAL

        for body_id in SIDEREAL_URANIANS:
            bname = URANIAN_NAMES[body_id]
            for jd in jds_sid:
                try:
                    lib_res, _ = lib.calc_ut(jd, body_id, flags_sid)
                except Exception as e:
                    check(False, f"Sid({ayan_name}) {bname} lib crash: {e}")
                    continue
                try:
                    ref_res, _ = swe_ref.calc_ut(jd, body_id, flags_sid)
                except Exception as e:
                    check(False, f"Sid({ayan_name}) {bname} ref crash: {e}")
                    continue

                lon_diff = arcsec(angular_diff(lib_res[0], ref_res[0]))
                lat_diff = arcsec(abs(lib_res[1] - ref_res[1]))
                dist_diff = abs(lib_res[2] - ref_res[2])

                track_max(f"13sid|{ayan_name}|{bname}|lon", lon_diff)
                track_max(f"13sid|{ayan_name}|{bname}|lat", lat_diff)

                check(
                    lon_diff < SID_LON_TOL,
                    f'Sid({ayan_name}) {bname} lon={lon_diff:.4f}" jd={jd:.2f}',
                )
                check(
                    lat_diff < SID_LAT_TOL,
                    f'Sid({ayan_name}) {bname} lat={lat_diff:.4f}" jd={jd:.2f}',
                )
                check(
                    dist_diff < DIST_TOL,
                    f"Sid({ayan_name}) {bname} dist={dist_diff:.6f} jd={jd:.2f}",
                )

    # Reset sidereal mode
    lib.set_sid_mode(0)
    swe_ref.set_sid_mode(0)

    # Print max diffs for diagnostics
    print('  Max diffs (Uranians geocentric, lon"):')
    for body_id in URANIAN_BODIES:
        bname = URANIAN_NAMES[body_id]
        k = f"13geo|{bname}|lon"
        print(f'    {bname}: {max_diffs.get(k, 0):.2f}"')

    print(f"  Done: {section_stats[_current_section]}")


# ===================================================================
# SECTION 14: Asteroids vs pyswisseph
# ===================================================================
def run_section_14() -> None:
    set_section("Sec14 Asteroids")
    print(f"\n{'=' * 70}")
    print("SECTION 14: Asteroids vs pyswisseph")
    print("  5 bodies x 100 dates x 3 checks = 1,500 + AST_OFFSET mapping")
    print(f"{'=' * 70}")

    N_DATES = 100
    # Pallas (18) and Juno (19) can differ by up to ~5.5" from pyswisseph
    # due to different ephemeris sources (JPL DE440 vs Swiss Ephemeris internal).
    LON_TOL = 6.0  # arcsec
    LAT_TOL = 6.0  # arcsec
    DIST_TOL = 1e-4  # AU

    jds = random_jds(N_DATES)
    flags = FLG_SPEED | FLG_SWIEPH

    for body_id in ASTEROID_BODIES:
        bname = ASTEROID_NAMES[body_id]
        for jd in jds:
            try:
                lib_res, _ = lib.calc_ut(jd, body_id, flags)
            except Exception as e:
                check(False, f"{bname} lib crash jd={jd:.2f}: {e}")
                continue
            try:
                ref_res, _ = swe_ref.calc_ut(jd, body_id, flags)
            except Exception as e:
                check(False, f"{bname} ref crash jd={jd:.2f}: {e}")
                continue

            lon_diff = arcsec(angular_diff(lib_res[0], ref_res[0]))
            lat_diff = arcsec(abs(lib_res[1] - ref_res[1]))
            dist_diff = abs(lib_res[2] - ref_res[2])

            track_max(f"14|{bname}|lon", lon_diff)
            track_max(f"14|{bname}|lat", lat_diff)
            track_max(f"14|{bname}|dist", dist_diff)

            check(
                lon_diff < LON_TOL,
                f'{bname} lon={lon_diff:.4f}" > {LON_TOL}" jd={jd:.2f}',
            )
            check(
                lat_diff < LAT_TOL,
                f'{bname} lat={lat_diff:.4f}" > {LAT_TOL}" jd={jd:.2f}',
            )
            check(
                dist_diff < DIST_TOL,
                f"{bname} dist={dist_diff:.6f} > {DIST_TOL} jd={jd:.2f}",
            )

    # --- AST_OFFSET mapping checks ---
    # AST_OFFSET + 1 = 10001 should correspond to Ceres (asteroid #1)
    # In pyswisseph: body id 10001 = asteroid #1 = Ceres
    # Compare with CERES (17) for a few dates
    ast_jds = random_jds(3)
    for jd in ast_jds:
        try:
            # pyswisseph: body 10001 = asteroid catalog #1 = Ceres
            ref_ceres_by_offset, _ = swe_ref.calc_ut(jd, 10001, flags)
            ref_ceres_by_id, _ = swe_ref.calc_ut(jd, 17, flags)

            lon_diff = arcsec(angular_diff(ref_ceres_by_offset[0], ref_ceres_by_id[0]))
            check(
                lon_diff < 0.001,
                f'AST_OFFSET Ceres 10001 vs 17 lon={lon_diff:.6f}" jd={jd:.2f}',
            )
        except Exception as e:
            check(False, f"AST_OFFSET mapping check crash: {e}")

    # Print max diffs
    print('  Max diffs (asteroids, lon"):')
    for body_id in ASTEROID_BODIES:
        bname = ASTEROID_NAMES[body_id]
        k = f"14|{bname}|lon"
        print(f'    {bname}: {max_diffs.get(k, 0):.4f}"')

    print(f"  Done: {section_stats[_current_section]}")


# ===================================================================
# SECTION 16: Orbital elements
# ===================================================================
def run_section_16() -> None:
    set_section("Sec16 Orbital elements")
    print(f"\n{'=' * 70}")
    print("SECTION 16: Orbital elements")
    print("  5 bodies x 10 dates x 6 checks = 300")
    print(f"{'=' * 70}")

    N_DATES = 10
    jds = random_jds(N_DATES)

    for body_id in ORBITAL_BODIES:
        bname = ORBITAL_NAMES[body_id]
        known_a = KNOWN_SEMI_MAJOR[body_id]

        for jd in jds:
            # --- Orbital elements ---
            try:
                elems = lib.get_orbital_elements_ut(jd, body_id, FLG_SWIEPH)
            except Exception as e:
                for _ in range(6):
                    check(False, f"{bname} orbital_elements crash jd={jd:.2f}: {e}")
                continue

            # Check we get 50 elements
            check(
                len(elems) == 50,
                f"{bname} expected 50 elements, got {len(elems)} at jd={jd:.2f}",
            )

            semi_major = elems[0]  # a
            eccentricity = elems[1]  # e
            inclination = elems[2]  # i

            # Semi-major axis within 5% of known value
            if known_a > 0:
                rel_err = abs(semi_major - known_a) / known_a
                check(
                    rel_err < 0.05,
                    f"{bname} a={semi_major:.4f} vs known={known_a:.4f} "
                    f"(err={rel_err:.4f}) jd={jd:.2f}",
                )
            else:
                check(True)  # No known value to compare, pass

            # Eccentricity in [0, 1)
            check(
                0 <= eccentricity < 1,
                f"{bname} e={eccentricity:.6f} out of [0,1) jd={jd:.2f}",
            )

            # Inclination in [0, 180)
            check(
                0 <= inclination < 180,
                f"{bname} i={inclination:.4f} out of [0,180) jd={jd:.2f}",
            )

            # --- orbit_max_min_true_distance ---
            try:
                dmax, dmin, dcur = lib.orbit_max_min_true_distance(
                    jd + lib.deltat(jd), body_id, FLG_SWIEPH
                )
            except Exception as e:
                check(False, f"{bname} orbit_max_min crash jd={jd:.2f}: {e}")
                check(False, f"{bname} orbit_max_min crash (max>min) jd={jd:.2f}")
                continue

            # perihelion < aphelion => dmin < dmax
            check(
                dmax > dmin,
                f"{bname} dmax={dmax:.6f} <= dmin={dmin:.6f} jd={jd:.2f}",
            )

            # Current distance should be between min and max (with tolerance)
            tol_factor = 1.05  # 5% tolerance for current distance bound
            check(
                dmin / tol_factor <= dcur <= dmax * tol_factor,
                f"{bname} dcur={dcur:.6f} outside [{dmin:.6f}, {dmax:.6f}] jd={jd:.2f}",
            )

    print(f"  Done: {section_stats[_current_section]}")


# ===================================================================
# SECTION 18: Crossings and stations
# ===================================================================
def run_section_18() -> None:
    set_section("Sec18 Crossings/Stations")
    print(f"\n{'=' * 70}")
    print("SECTION 18: Crossings and stations")
    print("  Equinoxes + Solstices + Moon crossings + Stations")
    print(f"{'=' * 70}")

    SECONDS_PER_DAY = 86400.0
    TIME_TOL_SEC = 60.0  # 60 seconds tolerance
    TIME_TOL_JD = TIME_TOL_SEC / SECONDS_PER_DAY

    # --- 18a: Equinoxes (Vernal equinox, Sun crossing 0 deg) ---
    # 2000-2025, one per year = 26 checks
    for year in range(2000, 2026):
        jd_start = lib.julday(year, 1, 1, 0.0)
        try:
            lib_jd = lib.solcross_ut(0.0, jd_start, FLG_SWIEPH)
        except Exception as e:
            check(False, f"Vernal equinox {year} lib crash: {e}")
            continue
        try:
            ref_jd = swe_ref.solcross_ut(0.0, jd_start, FLG_SWIEPH)
        except Exception as e:
            check(False, f"Vernal equinox {year} ref crash: {e}")
            continue

        diff_sec = abs(lib_jd - ref_jd) * SECONDS_PER_DAY
        track_max("18|equinox|sec", diff_sec)
        check(
            diff_sec < TIME_TOL_SEC,
            f"Vernal equinox {year}: diff={diff_sec:.2f}s > {TIME_TOL_SEC}s",
        )

    # --- 18b: Solstices (90 and 270 degrees) ---
    for year in range(2000, 2026):
        jd_start = lib.julday(year, 1, 1, 0.0)

        # Summer solstice (90 deg)
        try:
            lib_jd = lib.solcross_ut(90.0, jd_start, FLG_SWIEPH)
            ref_jd = swe_ref.solcross_ut(90.0, jd_start, FLG_SWIEPH)
            diff_sec = abs(lib_jd - ref_jd) * SECONDS_PER_DAY
            track_max("18|solstice_summer|sec", diff_sec)
            check(
                diff_sec < TIME_TOL_SEC,
                f"Summer solstice {year}: diff={diff_sec:.2f}s",
            )
        except Exception as e:
            check(False, f"Summer solstice {year} crash: {e}")

        # Winter solstice (270 deg)
        try:
            lib_jd = lib.solcross_ut(270.0, jd_start, FLG_SWIEPH)
            ref_jd = swe_ref.solcross_ut(270.0, jd_start, FLG_SWIEPH)
            diff_sec = abs(lib_jd - ref_jd) * SECONDS_PER_DAY
            track_max("18|solstice_winter|sec", diff_sec)
            check(
                diff_sec < TIME_TOL_SEC,
                f"Winter solstice {year}: diff={diff_sec:.2f}s",
            )
        except Exception as e:
            check(False, f"Winter solstice {year} crash: {e}")

    # --- 18c: Moon crossings ---
    # 12 target longitudes (every 30 deg), starting from J2000
    jd_moon_start = JD_J2000
    for target_deg in range(0, 360, 30):
        try:
            lib_jd = lib.mooncross_ut(float(target_deg), jd_moon_start, FLG_SWIEPH)
        except Exception as e:
            check(False, f"Moon cross {target_deg}deg lib crash: {e}")
            continue
        try:
            ref_jd = swe_ref.mooncross_ut(
                float(target_deg), jd_moon_start, FLG_SWIEPH
            )
        except Exception as e:
            check(False, f"Moon cross {target_deg}deg ref crash: {e}")
            continue

        diff_sec = abs(lib_jd - ref_jd) * SECONDS_PER_DAY
        track_max("18|mooncross|sec", diff_sec)
        check(
            diff_sec < TIME_TOL_SEC,
            f"Moon cross {target_deg}deg: diff={diff_sec:.2f}s",
        )

    # --- 18d: Stations ---
    # Mercury: find 5 stations starting from J2000
    # Mars: find 5 stations starting from J2000
    for planet_id, planet_name, n_stations in [(2, "Mercury", 5), (4, "Mars", 5)]:
        jd_search = JD_J2000
        for i in range(n_stations):
            try:
                jd_station, stype = lib.find_station_ut(planet_id, jd_search)
            except Exception as e:
                check(False, f"{planet_name} station #{i} crash: {e}")
                jd_search += 120  # skip ahead
                continue

            # Verify speed near zero at the station
            try:
                flags_spd = FLG_SPEED | FLG_SWIEPH
                res, _ = lib.calc_ut(jd_station, planet_id, flags_spd)
                speed_lon = res[3]  # speed in longitude (deg/day)

                # At a station, speed should be very close to 0
                check(
                    abs(speed_lon) < 0.05,
                    f"{planet_name} station #{i} ({stype}) speed={speed_lon:.6f} deg/day "
                    f"at jd={jd_station:.4f}",
                )
            except Exception as e:
                check(False, f"{planet_name} station #{i} speed check crash: {e}")

            # Station type should be "SR" or "SD"
            check(
                stype in ("SR", "SD"),
                f"{planet_name} station #{i} type={stype} not SR/SD",
            )

            # Move search past this station
            jd_search = jd_station + 30  # skip 30 days ahead

    print(f"  Done: {section_stats[_current_section]}")


# ===================================================================
# Main
# ===================================================================
def main() -> None:
    global passed, failed

    print("=" * 70)
    print("VERIFICATION: LEB, Uranians, Asteroids, Orbital Elements, Crossings")
    print("=" * 70)
    t0 = time.time()

    # Ensure we start in auto mode (LEB if available)
    lib.set_calc_mode("auto")

    try:
        run_section_1_3()
    except Exception as e:
        print(f"  SECTION 1.3 CRASHED: {e}")
        traceback.print_exc()

    try:
        run_section_1_4()
    except Exception as e:
        print(f"  SECTION 1.4 CRASHED: {e}")
        traceback.print_exc()

    try:
        run_section_13()
    except Exception as e:
        print(f"  SECTION 13 CRASHED: {e}")
        traceback.print_exc()

    try:
        run_section_14()
    except Exception as e:
        print(f"  SECTION 14 CRASHED: {e}")
        traceback.print_exc()

    try:
        run_section_16()
    except Exception as e:
        print(f"  SECTION 16 CRASHED: {e}")
        traceback.print_exc()

    try:
        run_section_18()
    except Exception as e:
        print(f"  SECTION 18 CRASHED: {e}")
        traceback.print_exc()

    elapsed = time.time() - t0

    # ── Summary ──
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}")
    for sec_name, stats in section_stats.items():
        total = stats["passed"] + stats["failed"]
        status = "PASS" if stats["failed"] == 0 else "FAIL"
        print(f"  [{status}] {sec_name}: {stats['passed']}/{total} passed")

    total = passed + failed
    print(f"\n  TOTAL: {passed}/{total} passed, {failed} failed")
    print(f"  Time:  {elapsed:.1f}s")

    if errors:
        print("\n  First failures (max 100):")
        for e in errors[:100]:
            print(f"    - {e}")

    if failed > 0:
        print(f"\n  RESULT: FAILED ({failed} failures)")
        sys.exit(1)
    else:
        print("\n  RESULT: ALL PASSED")
        sys.exit(0)


if __name__ == "__main__":
    main()
