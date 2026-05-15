#!/usr/bin/env python3
"""
Round 14: Deep Asteroid & Minor Body Position Verification
==========================================================

Comprehensive comparison of asteroid/minor body positions between
libephemeris (Skyfield + JPL SPK) and pyswisseph across multiple
body categories, epochs, and flag combinations.

Parts:
  P1: Big 4 + Chiron + Pholus (dedicated IDs) — 10 epochs × 6 bodies × 4 flag sets
  P2: Main belt asteroids — 5 epochs × 6 bodies
  P3: Additional centaurs (Nessus, Asbolus, Chariklo) — 5 epochs
  P4: TNOs (Eris, Sedna, Haumea, Makemake, etc.) — 5 epochs × 9 bodies
  P5: NEAs (Eros, Apophis, etc.) — 5 epochs × 9 bodies
  P6: Astrological asteroids — 5 epochs × 4 bodies
  P7: Multi-flag stress test — Big 4 × 5 epochs × 8 flag combos
  P8: Speed verification — all categories × 3 epochs
"""

from __future__ import annotations

import math
import os
import sys
import time
import traceback

# Force Skyfield mode
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import libephemeris as ephem
from libephemeris import spk_auto
from libephemeris.constants import (
    AST_OFFSET,
    CERES,
    CHIRON,
    JUNO,
    PALLAS,
    PHOLUS,
    VESTA,
    FLG_EQUATORIAL,
    FLG_HELCTR,
    FLG_J2000,
    FLG_NOABERR,
    FLG_NOGDEFL,
    FLG_NONUT,
    FLG_SPEED,
    FLG_TRUEPOS,
    FLG_XYZ,
)

# Set pyswisseph ephemeris path
_EPHE_PATH = os.path.join(os.path.dirname(__file__), "..", "swisseph", "ephe")
swe.set_ephe_path(_EPHE_PATH)

# Enable auto-SPK for all bodies
spk_auto.enable_common_bodies()


# ============================================================================
# Helpers
# ============================================================================


def se_hsys(ch: str) -> bytes:
    return ch.encode("ascii")


def le_hsys(ch: str) -> int:
    return ord(ch)


def angular_diff(a: float, b: float) -> float:
    """Angular difference handling 360-wrap."""
    d = abs(a - b)
    if d > 180:
        d = 360 - d
    return d


def arcsec(deg: float) -> float:
    """Convert degrees to arcseconds."""
    return deg * 3600.0


def fmt_arcsec(deg: float) -> str:
    """Format a degree value as arcseconds string."""
    return f'{arcsec(deg):.3f}"'


# ============================================================================
# Body definitions
# ============================================================================

# Dedicated IDs (15-20)
BIG6 = [
    (CHIRON, "Chiron"),
    (PHOLUS, "Pholus"),
    (CERES, "Ceres"),
    (PALLAS, "Pallas"),
    (JUNO, "Juno"),
    (VESTA, "Vesta"),
]

# Main belt (via AST_OFFSET)
MAIN_BELT = [
    (AST_OFFSET + 10, "Hygiea"),
    (AST_OFFSET + 704, "Interamnia"),
    (AST_OFFSET + 511, "Davida"),
    (AST_OFFSET + 52, "Europa(ast)"),
    (AST_OFFSET + 87, "Sylvia"),
    (AST_OFFSET + 16, "Psyche"),
]

# Additional centaurs
CENTAURS = [
    (AST_OFFSET + 7066, "Nessus"),
    (AST_OFFSET + 8405, "Asbolus"),
    (AST_OFFSET + 10199, "Chariklo"),
]

# Trans-Neptunian Objects
TNOS = [
    (AST_OFFSET + 136199, "Eris"),
    (AST_OFFSET + 90377, "Sedna"),
    (AST_OFFSET + 136108, "Haumea"),
    (AST_OFFSET + 136472, "Makemake"),
    (AST_OFFSET + 90482, "Orcus"),
    (AST_OFFSET + 50000, "Quaoar"),
    (AST_OFFSET + 20000, "Varuna"),
    (AST_OFFSET + 28978, "Ixion"),
    (AST_OFFSET + 225088, "Gonggong"),
]

# Near-Earth Asteroids
NEAS = [
    (AST_OFFSET + 433, "Eros"),
    (AST_OFFSET + 1221, "Amor"),
    (AST_OFFSET + 1566, "Icarus"),
    (AST_OFFSET + 1685, "Toro"),
    (AST_OFFSET + 99942, "Apophis"),
    (AST_OFFSET + 4179, "Toutatis"),
    (AST_OFFSET + 25143, "Itokawa"),
    (AST_OFFSET + 101955, "Bennu"),
    (AST_OFFSET + 162173, "Ryugu"),
]

# Astrological asteroids
ASTRO_AST = [
    (AST_OFFSET + 80, "Sappho"),
    (AST_OFFSET + 55, "Pandora"),
    (AST_OFFSET + 1181, "Lilith(ast)"),
    (AST_OFFSET + 944, "Hidalgo"),
]

# Test epochs
EPOCHS_10 = [
    (1900, 1, 1, 12.0, "1900"),
    (1950, 6, 15, 0.0, "1950"),
    (1970, 1, 1, 0.0, "1970"),
    (1985, 7, 4, 12.0, "1985"),
    (2000, 1, 1, 12.0, "J2000"),
    (2010, 3, 20, 18.0, "2010"),
    (2020, 6, 21, 0.0, "2020"),
    (2024, 11, 5, 12.0, "2024"),
    (2025, 3, 14, 0.0, "Today"),
    (2030, 12, 31, 23.99, "2030"),
]

EPOCHS_5 = EPOCHS_10[::2]  # 1900, 1970, J2000, 2020, Today

EPOCHS_3 = [EPOCHS_10[0], EPOCHS_10[4], EPOCHS_10[8]]  # 1900, J2000, Today


def jd_for(epoch):
    return swe.julday(epoch[0], epoch[1], epoch[2], epoch[3])


# ============================================================================
# Test runner
# ============================================================================


class TestResults:
    def __init__(self, part_name: str):
        self.part_name = part_name
        self.passed = 0
        self.failed = 0
        self.skipped = 0
        self.failures = []
        self.max_lon_diff = 0.0
        self.max_lat_diff = 0.0
        self.max_dist_diff = 0.0
        self.max_speed_diff = 0.0
        self.max_lon_body = ""
        self.max_lat_body = ""

    def record_pass(self, lon_diff, lat_diff, dist_diff, body_name="", speed_diff=0.0):
        self.passed += 1
        if lon_diff > self.max_lon_diff:
            self.max_lon_diff = lon_diff
            self.max_lon_body = body_name
        if lat_diff > self.max_lat_diff:
            self.max_lat_diff = lat_diff
            self.max_lat_body = body_name
        self.max_dist_diff = max(self.max_dist_diff, dist_diff)
        self.max_speed_diff = max(self.max_speed_diff, speed_diff)

    def record_fail(self, msg):
        self.failed += 1
        self.failures.append(msg)

    def record_skip(self, msg):
        self.skipped += 1

    def summary(self):
        total = self.passed + self.failed
        print(f"\n{'=' * 70}")
        print(
            f"  {self.part_name}: {self.passed}/{total} PASSED ({self.skipped} skipped)"
        )
        if self.max_lon_diff > 0:
            print(
                f"  Max lon diff: {fmt_arcsec(self.max_lon_diff)} ({self.max_lon_body})"
            )
        if self.max_lat_diff > 0:
            print(
                f"  Max lat diff: {fmt_arcsec(self.max_lat_diff)} ({self.max_lat_body})"
            )
        if self.max_dist_diff > 0:
            print(f"  Max dist diff: {self.max_dist_diff:.8f} AU")
        if self.max_speed_diff > 0:
            print(f"  Max speed diff: {self.max_speed_diff:.6f} deg/day")
        if self.failures:
            print(f"  FAILURES:")
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


def compare_body(
    jd,
    body_id,
    body_name,
    flags,
    tol_lon,
    tol_lat,
    tol_dist,
    results: TestResults,
    epoch_label="",
    check_speed=False,
    tol_speed=0.01,
):
    """Compare a single body between pyswisseph and libephemeris."""
    label = f"{body_name} @ {epoch_label} flags={flags}"

    try:
        res_se, retflag_se = swe.calc_ut(jd, body_id, flags)
    except Exception as e:
        results.record_skip(f"{label}: pyswisseph error: {e}")
        return

    try:
        res_le, retflag_le = ephem.calc_ut(jd, body_id, flags)
    except Exception as e:
        results.record_fail(f"{label}: libephemeris error: {e}")
        return

    lon_diff = angular_diff(res_se[0], res_le[0])
    lat_diff = abs(res_se[1] - res_le[1])
    dist_diff = abs(res_se[2] - res_le[2])

    speed_diff = 0.0
    speed_ok = True
    if check_speed and (flags & FLG_SPEED):
        speed_diff = abs(res_se[3] - res_le[3])
        # Handle speed wrap
        if speed_diff > 180:
            speed_diff = 360 - speed_diff
        if speed_diff >= tol_speed:
            speed_ok = False

    if lon_diff >= tol_lon:
        results.record_fail(
            f"{label}: LON {res_se[0]:.6f} vs {res_le[0]:.6f} "
            f"diff={fmt_arcsec(lon_diff)} (tol={fmt_arcsec(tol_lon)})"
        )
    elif lat_diff >= tol_lat:
        results.record_fail(
            f"{label}: LAT {res_se[1]:.6f} vs {res_le[1]:.6f} "
            f"diff={fmt_arcsec(lat_diff)} (tol={fmt_arcsec(tol_lat)})"
        )
    elif dist_diff >= tol_dist:
        results.record_fail(
            f"{label}: DIST {res_se[2]:.8f} vs {res_le[2]:.8f} "
            f"diff={dist_diff:.8f} (tol={tol_dist:.8f})"
        )
    elif not speed_ok:
        results.record_fail(
            f"{label}: SPEED {res_se[3]:.6f} vs {res_le[3]:.6f} "
            f"diff={speed_diff:.6f} (tol={tol_speed})"
        )
    else:
        results.record_pass(lon_diff, lat_diff, dist_diff, body_name, speed_diff)


# ============================================================================
# Part 1: Big 6 (dedicated IDs) — 10 epochs × 4 flag sets
# ============================================================================


def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Big 6 (Chiron, Pholus, Ceres, Pallas, Juno, Vesta)")
    print("  10 epochs × 6 bodies × 4 flag sets = 240 tests")
    print('  Tolerance: 1.0" lon/lat, 0.0001 AU dist')
    print("=" * 70)

    results = TestResults("P1: Big 6 Dedicated IDs")
    tol_lon = 1.0 / 3600  # 1 arcsecond
    tol_lat = 1.0 / 3600
    tol_dist = 0.0001

    flag_sets = [
        (FLG_SPEED, "SPEED"),
        (FLG_SPEED | FLG_J2000, "SPEED+J2000"),
        (FLG_SPEED | FLG_EQUATORIAL, "SPEED+EQUAT"),
        (FLG_SPEED | FLG_NONUT, "SPEED+NONUT"),
    ]

    for epoch in EPOCHS_10:
        jd = jd_for(epoch)
        for body_id, body_name in BIG6:
            for flags, flag_label in flag_sets:
                compare_body(
                    jd,
                    body_id,
                    body_name,
                    flags,
                    tol_lon,
                    tol_lat,
                    tol_dist,
                    results,
                    epoch_label=f"{epoch[4]} {flag_label}",
                    check_speed=True,
                    tol_speed=0.005,
                )

    return results.summary(), results


# ============================================================================
# Part 2: Main Belt Asteroids — 5 epochs
# ============================================================================


def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Main Belt Asteroids (Hygiea, Interamnia, Davida, etc.)")
    print("  5 epochs × 6 bodies = 30 tests")
    print('  Tolerance: 2.0" lon/lat')
    print("=" * 70)

    results = TestResults("P2: Main Belt")
    tol_lon = 2.0 / 3600  # 2 arcseconds
    tol_lat = 2.0 / 3600
    tol_dist = 0.001

    for epoch in EPOCHS_5:
        jd = jd_for(epoch)
        for body_id, body_name in MAIN_BELT:
            compare_body(
                jd,
                body_id,
                body_name,
                FLG_SPEED,
                tol_lon,
                tol_lat,
                tol_dist,
                results,
                epoch_label=epoch[4],
                check_speed=True,
                tol_speed=0.01,
            )

    return results.summary(), results


# ============================================================================
# Part 3: Additional Centaurs — 5 epochs
# ============================================================================


def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: Additional Centaurs (Nessus, Asbolus, Chariklo)")
    print("  5 epochs × 3 bodies = 15 tests")
    print('  Tolerance: 2.0" lon/lat')
    print("=" * 70)

    results = TestResults("P3: Centaurs")
    tol_lon = 2.0 / 3600
    tol_lat = 2.0 / 3600
    tol_dist = 0.001

    for epoch in EPOCHS_5:
        jd = jd_for(epoch)
        for body_id, body_name in CENTAURS:
            compare_body(
                jd,
                body_id,
                body_name,
                FLG_SPEED,
                tol_lon,
                tol_lat,
                tol_dist,
                results,
                epoch_label=epoch[4],
                check_speed=True,
                tol_speed=0.01,
            )

    return results.summary(), results


# ============================================================================
# Part 4: TNOs — 5 epochs
# ============================================================================


def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Trans-Neptunian Objects")
    print("  5 epochs × 9 bodies = 45 tests")
    print('  Tolerance: 5.0" lon/lat (distant, slower-moving)')
    print("=" * 70)

    results = TestResults("P4: TNOs")
    tol_lon = 5.0 / 3600  # 5 arcseconds
    tol_lat = 5.0 / 3600
    tol_dist = 0.01

    for epoch in EPOCHS_5:
        jd = jd_for(epoch)
        for body_id, body_name in TNOS:
            compare_body(
                jd,
                body_id,
                body_name,
                FLG_SPEED,
                tol_lon,
                tol_lat,
                tol_dist,
                results,
                epoch_label=epoch[4],
                check_speed=True,
                tol_speed=0.005,
            )

    return results.summary(), results


# ============================================================================
# Part 5: Near-Earth Asteroids — 5 epochs
# ============================================================================


def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Near-Earth Asteroids")
    print("  5 epochs × 9 bodies = 45 tests")
    print('  Tolerance: 5.0" lon/lat (fast-moving, eccentric orbits)')
    print("=" * 70)

    results = TestResults("P5: NEAs")
    tol_lon = 5.0 / 3600
    tol_lat = 5.0 / 3600
    tol_dist = 0.001

    for epoch in EPOCHS_5:
        jd = jd_for(epoch)
        for body_id, body_name in NEAS:
            compare_body(
                jd,
                body_id,
                body_name,
                FLG_SPEED,
                tol_lon,
                tol_lat,
                tol_dist,
                results,
                epoch_label=epoch[4],
                check_speed=True,
                tol_speed=0.05,
            )

    return results.summary(), results


# ============================================================================
# Part 6: Astrological Asteroids — 5 epochs
# ============================================================================


def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Astrological Asteroids (Sappho, Pandora, Lilith, Hidalgo)")
    print("  5 epochs × 4 bodies = 20 tests")
    print('  Tolerance: 2.0" lon/lat')
    print("=" * 70)

    results = TestResults("P6: Astrological Asteroids")
    tol_lon = 2.0 / 3600
    tol_lat = 2.0 / 3600
    tol_dist = 0.001

    for epoch in EPOCHS_5:
        jd = jd_for(epoch)
        for body_id, body_name in ASTRO_AST:
            compare_body(
                jd,
                body_id,
                body_name,
                FLG_SPEED,
                tol_lon,
                tol_lat,
                tol_dist,
                results,
                epoch_label=epoch[4],
                check_speed=True,
                tol_speed=0.01,
            )

    return results.summary(), results


# ============================================================================
# Part 7: Multi-flag stress test — Big 4 × 5 epochs × 8 flag combos
# ============================================================================


def run_part7():
    print("\n" + "=" * 70)
    print("PART 7: Multi-Flag Stress Test")
    print("  Big 4 × 5 epochs × 8 flag combos = 160 tests")
    print('  Tolerance: 1.5" lon/lat')
    print("=" * 70)

    results = TestResults("P7: Multi-Flag Stress")
    tol_lon = 1.5 / 3600
    tol_lat = 1.5 / 3600
    tol_dist = 0.001

    bodies = [
        (CERES, "Ceres"),
        (PALLAS, "Pallas"),
        (CHIRON, "Chiron"),
        (VESTA, "Vesta"),
    ]

    flag_combos = [
        (0, "default"),
        (FLG_SPEED, "SPEED"),
        (FLG_SPEED | FLG_J2000, "SPEED+J2000"),
        (FLG_SPEED | FLG_EQUATORIAL, "SPEED+EQUAT"),
        (FLG_SPEED | FLG_NONUT, "SPEED+NONUT"),
        (FLG_SPEED | FLG_TRUEPOS, "SPEED+TRUEPOS"),
        (FLG_SPEED | FLG_NOABERR, "SPEED+NOABERR"),
        (FLG_SPEED | FLG_NOGDEFL, "SPEED+NOGDEFL"),
    ]

    for epoch in EPOCHS_5:
        jd = jd_for(epoch)
        for body_id, body_name in bodies:
            for flags, flag_label in flag_combos:
                compare_body(
                    jd,
                    body_id,
                    body_name,
                    flags,
                    tol_lon,
                    tol_lat,
                    tol_dist,
                    results,
                    epoch_label=f"{epoch[4]} {flag_label}",
                )

    return results.summary(), results


# ============================================================================
# Part 8: Speed verification — all categories × 3 epochs
# ============================================================================


def run_part8():
    print("\n" + "=" * 70)
    print("PART 8: Speed Component Verification")
    print("  All bodies × 3 epochs, checking lon/lat/dist speeds")
    print("  Tolerance: 0.01 deg/day lon speed, 0.005 deg/day lat speed")
    print("=" * 70)

    results = TestResults("P8: Speed Verification")

    all_bodies = BIG6 + MAIN_BELT + CENTAURS + TNOS[:5] + NEAS[:5] + ASTRO_AST
    tol_lon_speed = 0.01  # deg/day
    tol_lat_speed = 0.005  # deg/day
    tol_dist_speed = 0.001  # AU/day

    for epoch in EPOCHS_3:
        jd = jd_for(epoch)
        for body_id, body_name in all_bodies:
            label = f"{body_name} @ {epoch[4]}"
            try:
                res_se, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
            except Exception as e:
                results.record_skip(f"{label}: pyswisseph error: {e}")
                continue
            try:
                res_le, _ = ephem.calc_ut(jd, body_id, FLG_SPEED)
            except Exception as e:
                results.record_fail(f"{label}: libephemeris error: {e}")
                continue

            # Speed components are indices 3, 4, 5
            lon_speed_diff = abs(res_se[3] - res_le[3])
            if lon_speed_diff > 180:
                lon_speed_diff = 360 - lon_speed_diff
            lat_speed_diff = abs(res_se[4] - res_le[4])
            dist_speed_diff = abs(res_se[5] - res_le[5])

            fail_parts = []
            if lon_speed_diff >= tol_lon_speed:
                fail_parts.append(
                    f"lon_speed: {res_se[3]:.6f} vs {res_le[3]:.6f} "
                    f"diff={lon_speed_diff:.6f}"
                )
            if lat_speed_diff >= tol_lat_speed:
                fail_parts.append(
                    f"lat_speed: {res_se[4]:.6f} vs {res_le[4]:.6f} "
                    f"diff={lat_speed_diff:.6f}"
                )
            if dist_speed_diff >= tol_dist_speed:
                fail_parts.append(
                    f"dist_speed: {res_se[5]:.8f} vs {res_le[5]:.8f} "
                    f"diff={dist_speed_diff:.8f}"
                )

            if fail_parts:
                results.record_fail(f"{label}: {'; '.join(fail_parts)}")
            else:
                results.record_pass(0, 0, 0, body_name, speed_diff=lon_speed_diff)

    return results.summary(), results


# ============================================================================
# Part 9: Heliocentric positions — Big 6 × 5 epochs
# ============================================================================


def run_part9():
    print("\n" + "=" * 70)
    print("PART 9: Heliocentric Positions")
    print("  Big 6 × 5 epochs with FLG_HELCTR")
    print('  Tolerance: 1.0" lon/lat')
    print("=" * 70)

    results = TestResults("P9: Heliocentric")
    tol_lon = 1.0 / 3600
    tol_lat = 1.0 / 3600
    tol_dist = 0.0001

    flags = FLG_SPEED | FLG_HELCTR

    for epoch in EPOCHS_5:
        jd = jd_for(epoch)
        for body_id, body_name in BIG6:
            compare_body(
                jd,
                body_id,
                body_name,
                flags,
                tol_lon,
                tol_lat,
                tol_dist,
                results,
                epoch_label=epoch[4],
                check_speed=True,
                tol_speed=0.005,
            )

    return results.summary(), results


# ============================================================================
# Part 10: J2000 equatorial (XYZ) — Big 4 × 3 epochs
# ============================================================================


def run_part10():
    print("\n" + "=" * 70)
    print("PART 10: J2000 Equatorial + XYZ")
    print("  Big 4 × 3 epochs with FLG_EQUATORIAL|J2000 and FLG_XYZ")
    print('  Tolerance: 1.0" (equat), 0.0001 AU (XYZ)')
    print("=" * 70)

    results = TestResults("P10: Equatorial/XYZ")

    bodies = [
        (CERES, "Ceres"),
        (CHIRON, "Chiron"),
        (PALLAS, "Pallas"),
        (VESTA, "Vesta"),
    ]

    # Equatorial J2000
    flags_eq = FLG_SPEED | FLG_EQUATORIAL | FLG_J2000
    tol_lon = 1.0 / 3600
    tol_lat = 1.0 / 3600

    for epoch in EPOCHS_3:
        jd = jd_for(epoch)
        for body_id, body_name in bodies:
            compare_body(
                jd,
                body_id,
                body_name,
                flags_eq,
                tol_lon,
                tol_lat,
                0.0001,
                results,
                epoch_label=f"{epoch[4]} EQ_J2000",
            )

    # XYZ
    flags_xyz = FLG_SPEED | FLG_XYZ
    for epoch in EPOCHS_3:
        jd = jd_for(epoch)
        for body_id, body_name in bodies:
            label = f"{body_name} @ {epoch[4]} XYZ"
            try:
                res_se, _ = swe.calc_ut(jd, body_id, flags_xyz)
            except Exception as e:
                results.record_skip(f"{label}: pyswisseph: {e}")
                continue
            try:
                res_le, _ = ephem.calc_ut(jd, body_id, flags_xyz)
            except Exception as e:
                results.record_fail(f"{label}: libephemeris: {e}")
                continue

            # XYZ: indices 0,1,2 are x,y,z in AU
            x_diff = abs(res_se[0] - res_le[0])
            y_diff = abs(res_se[1] - res_le[1])
            z_diff = abs(res_se[2] - res_le[2])
            max_diff = max(x_diff, y_diff, z_diff)

            tol_xyz = 0.0001  # AU
            if max_diff >= tol_xyz:
                results.record_fail(
                    f"{label}: XYZ max diff={max_diff:.8f} AU "
                    f"(x={x_diff:.8f}, y={y_diff:.8f}, z={z_diff:.8f})"
                )
            else:
                results.record_pass(max_diff, 0, 0, body_name)

    return results.summary(), results


# ============================================================================
# Main
# ============================================================================


def main():
    print("=" * 70)
    print("ROUND 14: Deep Asteroid & Minor Body Position Verification")
    print("=" * 70)

    start = time.time()
    all_ok = True
    all_results = []

    parts = [
        ("P1", run_part1),
        ("P2", run_part2),
        ("P3", run_part3),
        ("P4", run_part4),
        ("P5", run_part5),
        ("P6", run_part6),
        ("P7", run_part7),
        ("P8", run_part8),
        ("P9", run_part9),
        ("P10", run_part10),
    ]

    for part_name, part_fn in parts:
        try:
            ok, results = part_fn()
            all_results.append((part_name, results))
            if not ok:
                all_ok = False
        except Exception as e:
            print(f"\n  {part_name} CRASHED: {e}")
            traceback.print_exc()
            all_ok = False

    elapsed = time.time() - start

    # Final summary
    print("\n" + "=" * 70)
    print("ROUND 14 FINAL SUMMARY")
    print("=" * 70)

    total_passed = 0
    total_failed = 0
    total_skipped = 0

    for part_name, r in all_results:
        status = "PASS" if r.failed == 0 else "FAIL"
        total = r.passed + r.failed
        print(
            f"  {part_name} {r.part_name}: {r.passed}/{total} "
            f"({r.skipped} skip) [{status}]"
            f"  max_lon={fmt_arcsec(r.max_lon_diff)}"
        )
        total_passed += r.passed
        total_failed += r.failed
        total_skipped += r.skipped

    total = total_passed + total_failed
    print(
        f"\n  TOTAL: {total_passed}/{total} PASSED, "
        f"{total_failed} FAILED, {total_skipped} SKIPPED"
    )
    print(f"  Time: {elapsed:.1f}s")

    if all_ok:
        print("\n  >>> ROUND 14: ALL PASSED <<<")
    else:
        print(f"\n  >>> ROUND 14: {total_failed} FAILURES <<<")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
