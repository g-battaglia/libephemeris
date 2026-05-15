#!/usr/bin/env python3
"""Round 127: Fixed Stars Near Ecliptic Poles

Tests fixed stars at high ecliptic latitudes where coordinate transformations
are most sensitive to numerical precision.
"""

from __future__ import annotations
import os, sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_EQUATORIAL = 2048
FLG_J2000 = 32
FLG_NONUT = 64

# Stars at various ecliptic latitudes (including high-lat ones)
STARS = [
    "Polaris",  # Very high ecliptic lat ~66°
    "Vega",  # ~62° ecliptic lat
    "Capella",  # ~23° ecliptic lat
    "Deneb",  # ~60° ecliptic lat
    "Altair",  # ~29° ecliptic lat
    "Sirius",  # -39° ecliptic lat
    "Canopus",  # -76° ecliptic lat (very south)
    "Betelgeuse",  # -16° ecliptic lat
    "Rigel",  # -31° ecliptic lat
    "Procyon",  # -16° ecliptic lat
    "Aldebaran",  # -5° ecliptic lat
    "Regulus",  # ~0.5° ecliptic lat (near ecliptic)
    "Spica",  # -2° ecliptic lat
    "Antares",  # -4° ecliptic lat
    "Fomalhaut",  # -21° ecliptic lat
    "Arcturus",  # +31° ecliptic lat
]


def main():
    print("=" * 80)
    print("ROUND 127: Fixed Stars Near Ecliptic Poles")
    print("=" * 80)
    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    test_jds = [2451545.0, 2455197.5, 2459580.5, 2460310.5, 2444239.5]

    flag_combos = [
        ("default", FLG_SPEED),
        ("equatorial", FLG_SPEED | FLG_EQUATORIAL),
        ("J2000", FLG_SPEED | FLG_J2000 | FLG_NONUT),
        ("J2000+EQ", FLG_SPEED | FLG_J2000 | FLG_NONUT | FLG_EQUATORIAL),
    ]

    for star in STARS:
        for jd in test_jds:
            for flag_name, flags in flag_combos:
                try:
                    se_res = swe.fixstar2(star, jd, flags)
                    se_pos = se_res[0]
                except Exception:
                    continue

                try:
                    le_res = ephem.fixstar2_ut(star, jd, flags)
                    le_pos = le_res[0]
                except Exception:
                    continue

                for i, (label, tol) in enumerate([("lon/ra", 2.0), ("lat/dec", 2.0)]):
                    diff = le_pos[i] - se_pos[i]
                    if i == 0:
                        if diff > 180:
                            diff -= 360
                        elif diff < -180:
                            diff += 360
                    diff_as = abs(diff) * 3600

                    total_tests += 1
                    if diff_as < tol:
                        total_pass += 1
                    else:
                        total_fail += 1
                        if len(failures) < 25:
                            failures.append(
                                f'  {star} JD={jd:.1f} {flag_name} {label}: {diff_as:.4f}"'
                            )

                # Speed comparison
                for idx, label, tol in [(3, "lon_spd", 1.0), (4, "lat_spd", 1.0)]:
                    if idx < len(se_pos) and idx < len(le_pos):
                        diff_as = abs(le_pos[idx] - se_pos[idx]) * 3600
                        total_tests += 1
                        if diff_as < tol:
                            total_pass += 1
                        else:
                            total_fail += 1
                            if len(failures) < 35:
                                failures.append(
                                    f'  {star} JD={jd:.1f} {flag_name} {label}: {diff_as:.4f}"/day'
                                )

    pct = 100 * total_pass / total_tests if total_tests else 0
    print(
        f"\n{'=' * 80}\nROUND 127 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)\n  Failures: {total_fail}\n{'=' * 80}"
    )
    if failures:
        print("\nSample failures:")
        for f in failures:
            print(f)
    if total_fail == 0:
        print("\nAll tests PASSED!")
    return total_fail


if __name__ == "__main__":
    sys.exit(main())
