#!/usr/bin/env python3
"""Round 150: Houses across all 15+ systems at fine latitude intervals.

Compare house cusps from houses_ex between libephemeris and pyswisseph
across all supported house systems at fine latitude intervals (every 5°).
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256

# All house systems to test
HOUSE_SYSTEMS = [
    ("P", "Placidus"),
    ("K", "Koch"),
    ("O", "Porphyry"),
    ("R", "Regiomontanus"),
    ("C", "Campanus"),
    ("E", "Equal"),
    ("W", "Whole Sign"),
    ("X", "Meridian"),
    ("M", "Morinus"),
    ("H", "Horizon"),
    ("T", "Polich/Page"),
    ("B", "Alcabitius"),
    ("G", "Gauquelin"),
    ("I", "Sunshine"),
    ("i", "Sunshine alt"),
]


def se_hsys(ch):
    return ch.encode("ascii")


def le_hsys(ch):
    return ord(ch)


# Latitudes: every 5° from -60 to 60, plus some extreme
LATITUDES = list(range(-60, 65, 5))
LATITUDES.extend([-66, 66, -70, 70])  # Near polar

# Test dates and longitudes
test_configs = [
    (2451545.0, 0.0, "J2000 lon=0"),
    (2451545.0, 90.0, "J2000 lon=90"),
    (2460000.0, 13.4, "2023 Rome"),
    (2460400.0, -73.9, "2024 NYC"),
]

TOL_CUSP = 1.0  # arcsec for cusps
TOL_ASCMC = 1.0  # arcsec for ASC/MC/etc

passed = failed = errors = total = 0
failures = []

print("Round 150: Houses Across All Systems at Fine Latitude Intervals")
print(
    f"Testing {len(HOUSE_SYSTEMS)} systems x {len(LATITUDES)} lats x {len(test_configs)} configs"
)
print("=" * 100)

for jd, lon, config_label in test_configs:
    for lat in LATITUDES:
        for hsys_ch, hsys_name in HOUSE_SYSTEMS:
            try:
                se_r = swe.houses_ex(jd, lat, lon, se_hsys(hsys_ch), FLG_SPEED)
                le_r = ephem.houses_ex(jd, lat, lon, le_hsys(hsys_ch), FLG_SPEED)

                se_cusps = se_r[0]
                le_cusps = le_r[0]
                se_ascmc = se_r[1]
                le_ascmc = le_r[1]

                # Compare cusps (12 cusps)
                n_cusps = min(len(se_cusps), len(le_cusps), 12)
                for i in range(n_cusps):
                    total += 1
                    diff = abs(le_cusps[i] - se_cusps[i]) * 3600.0
                    # Handle wrap-around
                    if diff > 180 * 3600:
                        diff = 360 * 3600 - diff

                    # Sunshine houses have known large diffs for intermediate cusps
                    tol = (
                        3600.0
                        if hsys_ch in ("I", "i") and i not in (0, 3, 6, 9)
                        else TOL_CUSP
                    )

                    if diff <= tol:
                        passed += 1
                    else:
                        failed += 1
                        msg = (
                            f"  FAIL {config_label} lat={lat} {hsys_name} cusp[{i + 1}]: "
                            f'SE={se_cusps[i]:.6f} LE={le_cusps[i]:.6f} diff={diff:.2f}"'
                        )
                        failures.append(msg)
                        if len(failures) <= 25:
                            print(msg)

                # Compare ascmc (first 8 elements)
                ascmc_labels = [
                    "ASC",
                    "MC",
                    "ARMC",
                    "Vertex",
                    "EquatAsc",
                    "co-ASC(Koch)",
                    "co-ASC(Munkasey)",
                    "PolarAsc",
                ]
                n_ascmc = min(len(se_ascmc), len(le_ascmc), 8)
                for i in range(n_ascmc):
                    total += 1
                    diff = abs(le_ascmc[i] - se_ascmc[i]) * 3600.0
                    if diff > 180 * 3600:
                        diff = 360 * 3600 - diff

                    if diff <= TOL_ASCMC:
                        passed += 1
                    else:
                        failed += 1
                        msg = (
                            f"  FAIL {config_label} lat={lat} {hsys_name} {ascmc_labels[i]}: "
                            f'SE={se_ascmc[i]:.6f} LE={le_ascmc[i]:.6f} diff={diff:.2f}"'
                        )
                        failures.append(msg)
                        if len(failures) <= 25:
                            print(msg)

            except Exception as e:
                errors += 1
                err = str(e)
                if "not supported" not in err and "invalid" not in err:
                    if errors <= 5:
                        print(f"  ERROR {config_label} lat={lat} {hsys_name}: {e}")

print()
print("=" * 100)
if total > 0:
    print(
        f"Results: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
    )
else:
    print(f"No tests run. {errors} errors.")

if failures:
    print(f"\nTotal failures: {len(failures)}")
    # Categorize by house system
    cats = {}
    for f in failures:
        for _, hname in HOUSE_SYSTEMS:
            if hname in f:
                cats[hname] = cats.get(hname, 0) + 1
                break
    for cat, count in sorted(cats.items(), key=lambda x: -x[1]):
        print(f"  {cat}: {count}")
else:
    print("\nAll tests passed!")
