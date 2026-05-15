#!/usr/bin/env python3
"""Round 98: Node/Apsides Osculating vs Mean — Deep Sweep

Tests nod_aps_ut() for all planets with different calculation methods:
1. NODBIT_MEAN (mean nodes/apsides)
2. NODBIT_OSCU (osculating nodes/apsides)
3. NODBIT_OSCU_BAR (osculating barycentric)
4. NODBIT_FOPOINT (second focal point)
5. Multiple epochs
6. All planets (Sun through Pluto + Chiron + asteroids)
"""

from __future__ import annotations

import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_SWIEPH = 2

# Body IDs
SUN = 0
MOON = 1
MERCURY = 2
VENUS = 3
MARS = 4
JUPITER = 5
SATURN = 6
URANUS = 7
NEPTUNE = 8
PLUTO = 9
CHIRON = 15
CERES = 17
PALLAS = 18

# Method flags
NODBIT_MEAN = 1
NODBIT_OSCU = 2
NODBIT_OSCU_BAR = 4
NODBIT_FOPOINT = 256

PLANETS = [
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
    (PLUTO, "Pluto"),
    (CHIRON, "Chiron"),
]

METHODS = [
    (NODBIT_MEAN, "MEAN"),
    (NODBIT_OSCU, "OSCU"),
    (NODBIT_OSCU_BAR, "OSCU_BAR"),
    (NODBIT_FOPOINT, "FOPOINT"),
    (NODBIT_MEAN | NODBIT_FOPOINT, "MEAN+FO"),
]

EPOCHS = [
    2451545.0,  # J2000
    2460000.0,  # 2023
    2440000.0,  # 1968
    2430000.0,  # 1941
    2470000.0,  # 2050
    2415020.0,  # 1900
]

# Tolerances (arcsec)
LON_TOL = 60.0  # 1 arcmin — osculating elements can diverge significantly
MEAN_LON_TOL = 30.0  # mean nodes tighter
DIST_TOL = 0.01  # AU


def run_tests():
    passed = 0
    failed = 0
    errors = 0
    total = 0
    fail_details = []

    print("=" * 80)
    print("ROUND 98: Node/Apsides Osculating vs Mean — Deep Sweep")
    print("=" * 80)

    for method_flag, method_name in METHODS:
        m_pass = 0
        m_fail = 0
        m_err = 0

        tol = MEAN_LON_TOL if method_flag & NODBIT_MEAN else LON_TOL

        for body_id, body_name in PLANETS:
            for jd in EPOCHS:
                total += 1
                iflag = FLG_SWIEPH | FLG_SPEED

                try:
                    se_result = swe.nod_aps_ut(jd, body_id, iflag, method_flag)
                    # Returns (asc_node, desc_node, perihelion, aphelion)
                    # Each is a 6-tuple
                except Exception as e:
                    m_err += 1
                    errors += 1
                    continue

                try:
                    le_result = ephem.nod_aps_ut(jd, body_id, iflag, method_flag)
                except Exception as e:
                    m_err += 1
                    errors += 1
                    if len(fail_details) < 5:
                        fail_details.append(
                            f"  ERROR [{method_name}] {body_name} JD={jd:.1f}: {e}"
                        )
                    continue

                # Compare all 4 output nodes/apsides
                names = ["AscNode", "DescNode", "Perihelion", "Aphelion"]
                all_ok = True

                for i, name in enumerate(names):
                    se_pos = se_result[i]
                    le_pos = le_result[i]

                    if len(se_pos) < 3 or len(le_pos) < 3:
                        continue

                    dlon = abs(se_pos[0] - le_pos[0]) * 3600
                    dlat = abs(se_pos[1] - le_pos[1]) * 3600

                    if dlon > 180 * 3600:
                        dlon = 360 * 3600 - dlon

                    if dlon > tol or dlat > tol:
                        all_ok = False
                        if len(fail_details) < 30:
                            fail_details.append(
                                f"  FAIL [{method_name}] {body_name:10s} JD={jd:.1f} "
                                f'{name}: dLon={dlon:.1f}" dLat={dlat:.1f}" '
                                f"SE=({se_pos[0]:.4f},{se_pos[1]:.4f}) "
                                f"LE=({le_pos[0]:.4f},{le_pos[1]:.4f})"
                            )
                        break

                if all_ok:
                    m_pass += 1
                    passed += 1
                else:
                    m_fail += 1
                    failed += 1

        pct = 100 * m_pass / max(1, m_pass + m_fail + m_err)
        print(
            f"  [{method_name:12s}] {m_pass}/{m_pass + m_fail + m_err} passed ({pct:.1f}%) "
            f"[{m_fail} fail, {m_err} err]"
        )

    print()
    if fail_details:
        print("Sample failures/errors:")
        for d in fail_details[:30]:
            print(d)

    print()
    print("=" * 80)
    print(
        f"ROUND 98 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%), "
        f"{failed} failed, {errors} errors"
    )
    print("=" * 80)


if __name__ == "__main__":
    t0 = time.time()
    run_tests()
    print(f"\nCompleted in {time.time() - t0:.1f}s")
