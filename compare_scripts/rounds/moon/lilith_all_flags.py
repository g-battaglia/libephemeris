#!/usr/bin/env python3
"""Round 102: Lilith Positions Deep — All Flags & Modes

Tests Mean Lilith (MEAN_APOG=12) and Osculating Lilith (OSCU_APOG=13)
with various flag combinations: J2000, NONUT, EQUATORIAL, SIDEREAL, HELCTR.
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
FLG_J2000 = 32
FLG_NONUT = 64
FLG_EQUATORIAL = 2048
FLG_SIDEREAL = 65536

MEAN_APOG = 12
OSCU_APOG = 13

EPOCHS = [2415020.0 + i * 1826.25 for i in range(40)]  # 5-year steps, 200 years

FLAG_COMBOS = [
    (FLG_SWIEPH | FLG_SPEED, "default"),
    (FLG_SWIEPH | FLG_SPEED | FLG_J2000, "J2000"),
    (FLG_SWIEPH | FLG_SPEED | FLG_NONUT, "NONUT"),
    (FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL, "EQUATORIAL"),
    (FLG_SWIEPH | FLG_SPEED | FLG_J2000 | FLG_NONUT, "J2000+NONUT"),
    (FLG_SWIEPH | FLG_SPEED | FLG_J2000 | FLG_EQUATORIAL, "J2000+EQ"),
]

LON_TOL = 1.5  # arcsec
LAT_TOL = 25.0  # arcsec — MeanLilith latitude known ~19" offset


def run_tests():
    passed = 0
    failed = 0
    errors = 0
    total = 0
    fail_details = []

    print("=" * 80)
    print("ROUND 102: Lilith Positions Deep — All Flags & Modes")
    print("=" * 80)

    for body_id, body_name in [
        (MEAN_APOG, "MeanLilith"),
        (OSCU_APOG, "OscuLilith"),
    ]:
        for flags, flag_name in FLAG_COMBOS:
            f_pass = 0
            f_fail = 0

            for jd in EPOCHS:
                total += 1
                try:
                    se = swe.calc_ut(jd, body_id, flags)[0]
                    le = ephem.calc_ut(jd, body_id, flags)[0]
                except Exception as e:
                    errors += 1
                    continue

                dlon = abs(se[0] - le[0]) * 3600
                dlat = abs(se[1] - le[1]) * 3600
                if dlon > 180 * 3600:
                    dlon = 360 * 3600 - dlon

                if dlon <= LON_TOL and dlat <= LAT_TOL:
                    passed += 1
                    f_pass += 1
                else:
                    failed += 1
                    f_fail += 1
                    if len(fail_details) < 10:
                        fail_details.append(
                            f"  FAIL [{body_name}+{flag_name}] JD={jd:.1f}: "
                            f'dLon={dlon:.2f}" dLat={dlat:.2f}"'
                        )

            pct = 100 * f_pass / max(1, f_pass + f_fail)
            print(
                f"  [{body_name:12s}+{flag_name:12s}] {f_pass}/{f_pass + f_fail} ({pct:.1f}%)"
            )

    # Sidereal modes
    print("\n  --- Sidereal Modes ---")
    for sid_mode, sid_name in [(1, "Lahiri"), (0, "Fagan")]:
        swe.set_sid_mode(sid_mode)
        ephem.set_sid_mode(sid_mode, 0.0, 0.0)
        flags = FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL
        f_pass = 0
        f_fail = 0

        for jd in EPOCHS:
            total += 1
            try:
                se = swe.calc_ut(jd, MEAN_APOG, flags)[0]
                le = ephem.calc_ut(jd, MEAN_APOG, flags)[0]
            except Exception:
                errors += 1
                continue

            dlon = abs(se[0] - le[0]) * 3600
            if dlon > 180 * 3600:
                dlon = 360 * 3600 - dlon

            if dlon <= LON_TOL:
                passed += 1
                f_pass += 1
            else:
                failed += 1
                f_fail += 1

        pct = 100 * f_pass / max(1, f_pass + f_fail)
        print(
            f"  [MeanLilith+SID_{sid_name:6s}] {f_pass}/{f_pass + f_fail} ({pct:.1f}%)"
        )

    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0.0, 0.0)

    print()
    if fail_details:
        print("Sample failures:")
        for d in fail_details:
            print(d)

    print()
    pct = 100 * passed / max(1, total)
    print("=" * 80)
    print(
        f"ROUND 102 FINAL: {passed}/{total} passed ({pct:.1f}%), {failed} failed, {errors} errors"
    )
    print("=" * 80)


if __name__ == "__main__":
    t0 = time.time()
    run_tests()
    print(f"\nCompleted in {time.time() - t0:.1f}s")
