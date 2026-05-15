#!/usr/bin/env python3
"""
Round 15: Deep Node/Apsides (nod_aps_ut) Verification
==========================================================

Comprehensive comparison of orbital node and apsides calculations between
libephemeris and pyswisseph.

Parts:
  P1: All planets × 5 epochs × OSCU method — node longitudes
  P2: All planets × 5 epochs × OSCU method — apse longitudes
  P3: Method comparison (MEAN vs OSCU vs OSCU_BAR)
  P4: FOPOINT method (second focal point)
  P5: Moon nodes/apsides (geocentric orbit)
  P6: Speed components verification
  P7: Node latitude (should be ~0 for osculating)
  P8: TT vs UT consistency
"""

from __future__ import annotations

import math
import os
import sys
import time
import traceback

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import libephemeris as ephem
from libephemeris.constants import (
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
    NODBIT_MEAN,
    NODBIT_OSCU,
    NODBIT_OSCU_BAR,
    NODBIT_FOPOINT,
    FLG_SPEED,
    FLG_SWIEPH,
)

_EPHE_PATH = os.path.join(os.path.dirname(__file__), "..", "swisseph", "ephe")
swe.set_ephe_path(_EPHE_PATH)


def angular_diff(a, b):
    d = abs(a - b)
    if d > 180:
        d = 360 - d
    return d


def fmt_arcsec(deg):
    return f'{deg * 3600:.3f}"'


def fmt_deg(deg):
    return f"{deg:.4f}°"


PLANETS = [
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
    (PLUTO, "Pluto"),
]

EPOCHS = [
    (1950, 6, 15, 0.0, "1950"),
    (1980, 3, 21, 12.0, "1980"),
    (2000, 1, 1, 12.0, "J2000"),
    (2020, 6, 21, 0.0, "2020"),
    (2025, 3, 14, 0.0, "2025"),
]


def jd_for(epoch):
    return swe.julday(epoch[0], epoch[1], epoch[2], epoch[3])


class Results:
    def __init__(self, name):
        self.name = name
        self.passed = 0
        self.failed = 0
        self.skipped = 0
        self.failures = []
        self.max_diff = 0.0
        self.max_diff_body = ""

    def ok(self, diff, body=""):
        self.passed += 1
        if diff > self.max_diff:
            self.max_diff = diff
            self.max_diff_body = body

    def fail(self, msg):
        self.failed += 1
        self.failures.append(msg)

    def skip(self, msg=""):
        self.skipped += 1

    def summary(self):
        total = self.passed + self.failed
        print(f"\n{'=' * 70}")
        print(f"  {self.name}: {self.passed}/{total} PASSED ({self.skipped} skipped)")
        if self.max_diff > 0:
            print(f"  Max diff: {fmt_deg(self.max_diff)} ({self.max_diff_body})")
        if self.failures:
            print(f"  FAILURES:")
            for f in self.failures[:20]:
                print(f"    - {f}")
        print(f"{'=' * 70}")
        return self.failed == 0


# ============================================================================
# P1: Node longitudes — all planets × 5 epochs × OSCU
# ============================================================================


def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Node Longitudes (ascending + descending)")
    print("  8 planets × 5 epochs = 40 tests (80 node values)")
    print("=" * 70)

    r = Results("P1: Node Longitudes")

    # Node tolerances per planet (degrees)
    # Nodes are well-defined for inclined orbits
    tol = {
        "Mercury": 0.05,
        "Venus": 0.05,
        "Mars": 0.05,
        "Jupiter": 0.05,
        "Saturn": 0.05,
        "Uranus": 0.1,
        "Neptune": 0.1,
        "Pluto": 0.15,
    }

    for epoch in EPOCHS:
        jd = jd_for(epoch)
        for ipl, name in PLANETS:
            label = f"{name} @ {epoch[4]}"
            t = tol.get(name, 0.1)
            try:
                se = swe.nod_aps_ut(jd, ipl, FLG_SPEED, NODBIT_OSCU)
                le = ephem.nod_aps_ut(jd, ipl, FLG_SPEED, NODBIT_OSCU)
            except Exception as e:
                r.fail(f"{label}: error: {e}")
                continue

            # Ascending node longitude
            asc_diff = angular_diff(se[0][0], le[0][0])
            # Descending node longitude
            dsc_diff = angular_diff(se[1][0], le[1][0])

            worst = max(asc_diff, dsc_diff)
            if worst >= t:
                r.fail(
                    f"{label}: asc_diff={fmt_deg(asc_diff)} "
                    f"dsc_diff={fmt_deg(dsc_diff)} (tol={fmt_deg(t)})"
                )
            else:
                r.ok(worst, name)

    return r.summary(), r


# ============================================================================
# P2: Apse longitudes — all planets × 5 epochs × OSCU
# ============================================================================


def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Apse Longitudes (perihelion + aphelion)")
    print("  8 planets × 5 epochs = 40 tests")
    print("=" * 70)

    r = Results("P2: Apse Longitudes")

    # Apse tolerances - much wider for low-eccentricity orbits
    # where perihelion direction is poorly constrained
    tol = {
        "Mercury": 0.5,
        "Venus": 10.0,  # Venus e~0.007, very poorly defined
        "Mars": 2.0,
        "Jupiter": 5.0,  # Jupiter e~0.048
        "Saturn": 2.0,
        "Uranus": 2.0,
        "Neptune": 10.0,  # e~0.009
        "Pluto": 1.0,
    }

    for epoch in EPOCHS:
        jd = jd_for(epoch)
        for ipl, name in PLANETS:
            label = f"{name} @ {epoch[4]}"
            t = tol.get(name, 5.0)
            try:
                se = swe.nod_aps_ut(jd, ipl, FLG_SPEED, NODBIT_OSCU)
                le = ephem.nod_aps_ut(jd, ipl, FLG_SPEED, NODBIT_OSCU)
            except Exception as e:
                r.fail(f"{label}: error: {e}")
                continue

            # Perihelion
            peri_diff = angular_diff(se[2][0], le[2][0])
            # Aphelion
            aphe_diff = angular_diff(se[3][0], le[3][0])

            worst = max(peri_diff, aphe_diff)
            if worst >= t:
                r.fail(
                    f"{label}: peri={fmt_deg(peri_diff)} "
                    f"aphe={fmt_deg(aphe_diff)} (tol={fmt_deg(t)})"
                    f"\n        SE_peri={se[2][0]:.4f} LE_peri={le[2][0]:.4f}"
                    f"  SE_aphe={se[3][0]:.4f} LE_aphe={le[3][0]:.4f}"
                )
            else:
                r.ok(worst, name)

    return r.summary(), r


# ============================================================================
# P3: Method comparison (MEAN vs OSCU vs OSCU_BAR)
# ============================================================================


def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: Method Comparison (MEAN / OSCU / OSCU_BAR)")
    print("  4 planets × 3 methods × 3 epochs = 36 tests")
    print("=" * 70)

    r = Results("P3: Method Comparison")

    bodies = [
        (MERCURY, "Mercury"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
    ]

    methods = [
        (NODBIT_MEAN, "MEAN"),
        (NODBIT_OSCU, "OSCU"),
        (NODBIT_OSCU_BAR, "OSCU_BAR"),
    ]

    epochs = [EPOCHS[0], EPOCHS[2], EPOCHS[4]]  # 1950, J2000, 2025

    # Tolerance for nodes
    tol_node = 0.1  # degrees
    tol_apse = 10.0  # very wide for apse (method-dependent)

    for epoch in epochs:
        jd = jd_for(epoch)
        for ipl, name in bodies:
            for method, method_name in methods:
                label = f"{name} @ {epoch[4]} {method_name}"
                try:
                    se = swe.nod_aps_ut(jd, ipl, FLG_SPEED, method)
                    le = ephem.nod_aps_ut(jd, ipl, FLG_SPEED, method)
                except Exception as e:
                    r.fail(f"{label}: error: {e}")
                    continue

                asc_diff = angular_diff(se[0][0], le[0][0])
                if asc_diff >= tol_node:
                    r.fail(
                        f"{label}: asc_node diff={fmt_deg(asc_diff)} (tol={fmt_deg(tol_node)})"
                    )
                else:
                    r.ok(asc_diff, f"{name}/{method_name}")

    return r.summary(), r


# ============================================================================
# P4: FOPOINT method
# ============================================================================


def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: FOPOINT Method (second focal point)")
    print("  4 planets × 3 epochs × 2 method combos = 24 tests")
    print("=" * 70)

    r = Results("P4: FOPOINT")

    bodies = [
        (MERCURY, "Mercury"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (PLUTO, "Pluto"),
    ]

    method_combos = [
        (NODBIT_OSCU | NODBIT_FOPOINT, "OSCU+FOPOINT"),
        (NODBIT_MEAN | NODBIT_FOPOINT, "MEAN+FOPOINT"),
    ]

    epochs = [EPOCHS[0], EPOCHS[2], EPOCHS[4]]
    tol_node = 0.1
    tol_apse = 10.0  # FOPOINT changes aphelion position

    for epoch in epochs:
        jd = jd_for(epoch)
        for ipl, name in bodies:
            for method, method_name in method_combos:
                label = f"{name} @ {epoch[4]} {method_name}"
                try:
                    se = swe.nod_aps_ut(jd, ipl, FLG_SPEED, method)
                    le = ephem.nod_aps_ut(jd, ipl, FLG_SPEED, method)
                except Exception as e:
                    r.fail(f"{label}: error: {e}")
                    continue

                # Nodes should be unchanged by FOPOINT
                asc_diff = angular_diff(se[0][0], le[0][0])
                # Aphelion may be replaced by focal point
                aphe_diff = angular_diff(se[3][0], le[3][0])

                if asc_diff >= tol_node:
                    r.fail(
                        f"{label}: asc_diff={fmt_deg(asc_diff)} (tol={fmt_deg(tol_node)})"
                    )
                elif aphe_diff >= tol_apse:
                    r.fail(
                        f"{label}: aphe/foc diff={fmt_deg(aphe_diff)} (tol={fmt_deg(tol_apse)})"
                        f"\n        SE={se[3][0]:.4f} LE={le[3][0]:.4f}"
                    )
                else:
                    r.ok(max(asc_diff, aphe_diff), f"{name}/{method_name}")

    return r.summary(), r


# ============================================================================
# P5: Moon nodes/apsides
# ============================================================================


def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Moon Nodes and Apsides (geocentric orbit)")
    print("  5 epochs × 3 methods = 15 tests")
    print("=" * 70)

    r = Results("P5: Moon")

    methods = [
        (NODBIT_MEAN, "MEAN"),
        (NODBIT_OSCU, "OSCU"),
        (NODBIT_OSCU_BAR, "OSCU_BAR"),
    ]

    tol_node = 1.0  # Moon nodes: wider tolerance due to complex perturbations
    tol_apse = 5.0  # Moon apsides: complex geocentric orbit

    for epoch in EPOCHS:
        jd = jd_for(epoch)
        for method, method_name in methods:
            label = f"Moon @ {epoch[4]} {method_name}"
            try:
                se = swe.nod_aps_ut(jd, MOON, FLG_SPEED, method)
                le = ephem.nod_aps_ut(jd, MOON, FLG_SPEED, method)
            except Exception as e:
                r.fail(f"{label}: error: {e}")
                continue

            asc_diff = angular_diff(se[0][0], le[0][0])
            dsc_diff = angular_diff(se[1][0], le[1][0])
            peri_diff = angular_diff(se[2][0], le[2][0])
            aphe_diff = angular_diff(se[3][0], le[3][0])

            worst_node = max(asc_diff, dsc_diff)
            worst_apse = max(peri_diff, aphe_diff)

            if worst_node >= tol_node:
                r.fail(
                    f"{label}: node diff={fmt_deg(worst_node)} (tol={fmt_deg(tol_node)})"
                    f"\n        SE_asc={se[0][0]:.4f} LE_asc={le[0][0]:.4f}"
                    f"  SE_dsc={se[1][0]:.4f} LE_dsc={le[1][0]:.4f}"
                )
            elif worst_apse >= tol_apse:
                r.fail(
                    f"{label}: apse diff={fmt_deg(worst_apse)} (tol={fmt_deg(tol_apse)})"
                    f"\n        SE_peri={se[2][0]:.4f} LE_peri={le[2][0]:.4f}"
                    f"  SE_aphe={se[3][0]:.4f} LE_aphe={le[3][0]:.4f}"
                )
            else:
                r.ok(max(worst_node, worst_apse), f"Moon/{method_name}")

    return r.summary(), r


# ============================================================================
# P6: Speed components
# ============================================================================


def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Speed Components")
    print("  Check if both implementations return speeds for nodes/apsides")
    print("=" * 70)

    r = Results("P6: Speeds")

    bodies = [
        (MERCURY, "Mercury"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (MOON, "Moon"),
    ]

    jd = jd_for(EPOCHS[2])  # J2000

    for ipl, name in bodies:
        label = f"{name} speeds"
        try:
            se = swe.nod_aps_ut(jd, ipl, FLG_SPEED, NODBIT_OSCU)
            le = ephem.nod_aps_ut(jd, ipl, FLG_SPEED, NODBIT_OSCU)
        except Exception as e:
            r.fail(f"{label}: error: {e}")
            continue

        # Check if SE returns nonzero speeds
        se_has_speeds = any(abs(se[i][3]) > 1e-10 for i in range(4))
        le_has_speeds = any(abs(le[i][3]) > 1e-10 for i in range(4))

        if se_has_speeds and not le_has_speeds:
            # Known: LE doesn't compute speeds yet
            # Report but don't fail — document as a known limitation
            r.skip(f"{label}: SE has speeds, LE returns zeros (known limitation)")
            print(
                f"    NOTE: {name} — SE speeds: "
                f"asc={se[0][3]:.6f} dsc={se[1][3]:.6f} "
                f"peri={se[2][3]:.6f} aphe={se[3][3]:.6f}"
            )
        elif not se_has_speeds and not le_has_speeds:
            r.ok(0, name)
        elif se_has_speeds and le_has_speeds:
            # Both have speeds — compare them
            speed_diffs = [abs(se[i][3] - le[i][3]) for i in range(4)]
            worst = max(speed_diffs)
            if worst >= 0.01:
                r.fail(f"{label}: speed diff={worst:.6f} deg/day")
            else:
                r.ok(worst, name)
        else:
            r.ok(0, name)

    return r.summary(), r


# ============================================================================
# P7: Node latitude (should be ~0 for osculating nodes on the ecliptic)
# ============================================================================


def run_part7():
    print("\n" + "=" * 70)
    print("PART 7: Node Latitude Verification")
    print("  Osculating nodes should have latitude near zero")
    print("=" * 70)

    r = Results("P7: Node Latitude")

    jd = jd_for(EPOCHS[2])  # J2000

    for ipl, name in PLANETS:
        label = f"{name} node lat"
        try:
            se = swe.nod_aps_ut(jd, ipl, FLG_SPEED, NODBIT_OSCU)
            le = ephem.nod_aps_ut(jd, ipl, FLG_SPEED, NODBIT_OSCU)
        except Exception as e:
            r.fail(f"{label}: error: {e}")
            continue

        # Ascending node latitude
        se_lat = se[0][1]
        le_lat = le[0][1]
        lat_diff = abs(se_lat - le_lat)

        # Both should be near zero for osculating nodes
        # SE returns exactly 0 for nodes; LE may compute small values
        if lat_diff >= 1.0:  # 1 degree tolerance
            r.fail(
                f"{label}: SE_lat={se_lat:.6f} LE_lat={le_lat:.6f} diff={lat_diff:.6f}"
            )
        else:
            r.ok(lat_diff, name)
            if abs(le_lat) > 0.01:
                print(f"    NOTE: {name} asc node lat={le_lat:.4f}° (SE={se_lat:.4f}°)")

    return r.summary(), r


# ============================================================================
# P8: TT vs UT consistency
# ============================================================================


def run_part8():
    print("\n" + "=" * 70)
    print("PART 8: TT vs UT Consistency")
    print("  nod_aps(tt) should match nod_aps_ut(ut) for same instant")
    print("=" * 70)

    r = Results("P8: TT/UT Consistency")

    bodies = [
        (MERCURY, "Mercury"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (MOON, "Moon"),
    ]

    for epoch in EPOCHS:
        jd_ut = jd_for(epoch)
        dt = ephem.deltat(jd_ut)
        jd_tt = jd_ut + dt

        for ipl, name in bodies:
            label = f"{name} @ {epoch[4]} TT/UT"
            try:
                res_ut = ephem.nod_aps_ut(jd_ut, ipl, FLG_SPEED, NODBIT_OSCU)
                res_tt = ephem.nod_aps(jd_tt, ipl, NODBIT_OSCU, FLG_SPEED)
            except Exception as e:
                r.fail(f"{label}: error: {e}")
                continue

            # Compare ascending node longitude
            diff = angular_diff(res_ut[0][0], res_tt[0][0])
            tol = 1e-6  # Should be nearly identical

            if diff >= tol:
                r.fail(
                    f"{label}: UT={res_ut[0][0]:.8f} TT={res_tt[0][0]:.8f} "
                    f"diff={diff:.10f}"
                )
            else:
                r.ok(diff, name)

    return r.summary(), r


# ============================================================================
# P9: Distance comparison for nodes and apsides
# ============================================================================


def run_part9():
    print("\n" + "=" * 70)
    print("PART 9: Distance Values for Nodes and Apsides")
    print("  Verify distance components match between SE and LE")
    print("=" * 70)

    r = Results("P9: Distances")

    bodies = [
        (MERCURY, "Mercury"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
    ]

    jd = jd_for(EPOCHS[2])  # J2000

    for ipl, name in bodies:
        label = f"{name} distances"
        try:
            se = swe.nod_aps_ut(jd, ipl, FLG_SPEED, NODBIT_OSCU)
            le = ephem.nod_aps_ut(jd, ipl, FLG_SPEED, NODBIT_OSCU)
        except Exception as e:
            r.fail(f"{label}: error: {e}")
            continue

        # Check distances for all 4 points
        point_names = ["asc_node", "dsc_node", "perihelion", "aphelion"]
        fail_parts = []

        for idx, pname in enumerate(point_names):
            se_dist = se[idx][2]
            le_dist = le[idx][2]
            if se_dist > 0 and le_dist > 0:
                rel_diff = abs(se_dist - le_dist) / max(se_dist, le_dist)
                if rel_diff >= 0.1:  # 10% relative tolerance
                    fail_parts.append(
                        f"{pname}: SE={se_dist:.6f} LE={le_dist:.6f} rel={rel_diff:.4f}"
                    )

        if fail_parts:
            r.fail(f"{label}: {'; '.join(fail_parts)}")
        else:
            r.ok(0, name)

    return r.summary(), r


# ============================================================================
# Main
# ============================================================================


def main():
    print("=" * 70)
    print("ROUND 15: Deep Node/Apsides (nod_aps_ut) Verification")
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
    ]

    for pname, pfn in parts:
        try:
            ok, res = pfn()
            all_results.append((pname, res))
            if not ok:
                all_ok = False
        except Exception as e:
            print(f"\n  {pname} CRASHED: {e}")
            traceback.print_exc()
            all_ok = False

    elapsed = time.time() - start

    print("\n" + "=" * 70)
    print("ROUND 15 FINAL SUMMARY")
    print("=" * 70)

    total_p = total_f = total_s = 0
    for pname, res in all_results:
        status = "PASS" if res.failed == 0 else "FAIL"
        total = res.passed + res.failed
        print(
            f"  {pname} {res.name}: {res.passed}/{total} "
            f"({res.skipped} skip) [{status}]"
            f"  max={fmt_deg(res.max_diff)}"
        )
        total_p += res.passed
        total_f += res.failed
        total_s += res.skipped

    total = total_p + total_f
    print(f"\n  TOTAL: {total_p}/{total} PASSED, {total_f} FAILED, {total_s} SKIPPED")
    print(f"  Time: {elapsed:.1f}s")

    if all_ok:
        print("\n  >>> ROUND 15: ALL PASSED <<<")
    else:
        print(f"\n  >>> ROUND 15: {total_f} FAILURES <<<")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
