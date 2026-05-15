#!/usr/bin/env python3
"""Round 101: Mean/True Node & Lilith Speed Precision

Tests speed (velocity) values for lunar nodes and Lilith across epochs.
Verifies both position and speed agreement for:
1. TRUE_NODE (11) — true node position + speed
2. MEAN_NODE (10) — mean node position + speed
3. MEAN_APOG (12) — mean Lilith/apogee position + speed
4. OSCU_APOG (13) — osculating Lilith position + speed
"""

from __future__ import annotations
import sys, os, time
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_SWIEPH = 2

BODIES = [
    (10, "MeanNode"),
    (11, "TrueNode"),
    (12, "MeanLilith"),
    (13, "OscuLilith"),
]

# 100 epochs over 100 years
EPOCHS = [2415020.0 + i * 365.25 for i in range(100)]

# Tolerances
LON_TOL = 25.0  # arcsec — MeanNode has known offsets
LAT_TOL = 25.0  # arcsec — MeanLilith latitude ~19" known
SPEED_TOL = 5.0  # arcsec/day — mean node speed is constant vs SE perturbations


def run_tests():
    passed = 0
    failed = 0
    total = 0
    fail_details = []
    flags = FLG_SWIEPH | FLG_SPEED

    print("=" * 80)
    print("ROUND 101: Mean/True Node & Lilith Speed Precision")
    print("=" * 80)

    for body_id, body_name in BODIES:
        b_pass = 0
        b_fail = 0
        worst_lon = 0
        worst_spd = 0

        for jd in EPOCHS:
            total += 1
            try:
                se = swe.calc_ut(jd, body_id, flags)[0]
                le = ephem.calc_ut(jd, body_id, flags)[0]
            except Exception as e:
                failed += 1
                b_fail += 1
                continue

            dlon = abs(se[0] - le[0]) * 3600
            dlat = abs(se[1] - le[1]) * 3600
            dspd = abs(se[3] - le[3]) * 3600

            if dlon > 180 * 3600:
                dlon = 360 * 3600 - dlon

            worst_lon = max(worst_lon, dlon)
            worst_spd = max(worst_spd, dspd)

            if dlon <= LON_TOL and dlat <= LAT_TOL and dspd <= SPEED_TOL:
                passed += 1
                b_pass += 1
            else:
                failed += 1
                b_fail += 1
                if len(fail_details) < 10:
                    fail_details.append(
                        f"  FAIL {body_name:12s} JD={jd:.1f}: "
                        f'dLon={dlon:.2f}" dLat={dlat:.2f}" dSpd={dspd:.2f}"/d '
                        f'SE_spd={se[3] * 3600:.2f}"/d LE_spd={le[3] * 3600:.2f}"/d'
                    )

        pct = 100 * b_pass / max(1, b_pass + b_fail)
        print(
            f"  [{body_name:12s}] {b_pass}/{b_pass + b_fail} passed ({pct:.1f}%) "
            f'worst_lon={worst_lon:.2f}" worst_spd={worst_spd:.2f}"/d'
        )

    print()
    if fail_details:
        print("Sample failures:")
        for d in fail_details:
            print(d)

    print()
    print("=" * 80)
    pct = 100 * passed / max(1, total)
    print(f"ROUND 101 FINAL: {passed}/{total} passed ({pct:.1f}%), {failed} failed")
    print("=" * 80)


if __name__ == "__main__":
    t0 = time.time()
    run_tests()
    print(f"\nCompleted in {time.time() - t0:.1f}s")
