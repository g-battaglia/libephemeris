#!/usr/bin/env python3
"""Round 107: Moon Apsidal Precession (8.85yr cycle)

Verifies Mean Lilith (lunar apogee) position across multiple 8.85-year
apsidal precession cycles. Tests that the precession rate and position
track correctly over 50+ years.
"""

from __future__ import annotations
import sys, os, time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_SWIEPH = 2
MEAN_APOG = 12  # Mean Lilith
OSCU_APOG = 13  # Osculating Lilith

# Sample every 10 days over 50 years (covers ~5.6 apsidal cycles)
EPOCHS = [2430000.0 + i * 10.0 for i in range(1826)]

LON_TOL = 0.5  # arcsec
SPEED_TOL = 1.0  # arcsec/day


def run_tests():
    passed = 0
    failed = 0
    total = 0
    worst_lon = 0
    worst_spd = 0
    flags = FLG_SWIEPH | FLG_SPEED

    print("=" * 80)
    print("ROUND 107: Moon Apsidal Precession (8.85yr cycle)")
    print("=" * 80)

    # PART 1: Mean Lilith positions
    print("\n--- PART 1: Mean Lilith (MEAN_APOG) over 50 years ---")
    p1_pass = 0
    p1_fail = 0

    for jd in EPOCHS:
        total += 1
        try:
            se = swe.calc_ut(jd, MEAN_APOG, flags)[0]
            le = ephem.calc_ut(jd, MEAN_APOG, flags)[0]
        except Exception:
            failed += 1
            p1_fail += 1
            continue

        dlon = abs(se[0] - le[0]) * 3600
        if dlon > 180 * 3600:
            dlon = 360 * 3600 - dlon
        dspd = abs(se[3] - le[3]) * 3600

        worst_lon = max(worst_lon, dlon)
        worst_spd = max(worst_spd, dspd)

        if dlon <= LON_TOL and dspd <= SPEED_TOL:
            passed += 1
            p1_pass += 1
        else:
            failed += 1
            p1_fail += 1

    print(
        f"  Part 1: {p1_pass}/{p1_pass + p1_fail} passed "
        f'(worst lon={worst_lon:.3f}" spd={worst_spd:.3f}"/d)'
    )

    # PART 2: Osculating Lilith positions
    print("\n--- PART 2: Osculating Lilith (OSCU_APOG) over 50 years ---")
    p2_pass = 0
    p2_fail = 0
    worst_lon2 = 0

    for jd in EPOCHS[::5]:  # Every 50 days (faster)
        total += 1
        try:
            se = swe.calc_ut(jd, OSCU_APOG, flags)[0]
            le = ephem.calc_ut(jd, OSCU_APOG, flags)[0]
        except Exception:
            failed += 1
            p2_fail += 1
            continue

        dlon = abs(se[0] - le[0]) * 3600
        if dlon > 180 * 3600:
            dlon = 360 * 3600 - dlon

        worst_lon2 = max(worst_lon2, dlon)

        if dlon <= 2.0:  # Wider tolerance for osculating
            passed += 1
            p2_pass += 1
        else:
            failed += 1
            p2_fail += 1

    print(
        f'  Part 2: {p2_pass}/{p2_pass + p2_fail} passed (worst lon={worst_lon2:.3f}")'
    )

    # PART 3: Precession rate verification
    print("\n--- PART 3: Apsidal Precession Rate ---")
    p3_pass = 0
    p3_fail = 0

    # Mean Lilith advances ~40.69°/year = ~0.1115°/day
    expected_rate = 40.69 / 365.25 * 3600  # arcsec/day

    for jd in EPOCHS[::50]:
        total += 1
        le = ephem.calc_ut(jd, MEAN_APOG, flags)[0]
        le_speed = le[3] * 3600  # arcsec/day

        rate_diff = abs(le_speed - expected_rate)

        if rate_diff < 5.0:  # 5"/day tolerance on ~400"/day rate
            passed += 1
            p3_pass += 1
        else:
            failed += 1
            p3_fail += 1

    print(f"  Part 3: {p3_pass}/{p3_pass + p3_fail} passed")

    pct = 100 * passed / max(1, total)
    print()
    print("=" * 80)
    print(f"ROUND 107 FINAL: {passed}/{total} passed ({pct:.1f}%), {failed} failed")
    print("=" * 80)


if __name__ == "__main__":
    t0 = time.time()
    run_tests()
    print(f"\nCompleted in {time.time() - t0:.1f}s")
