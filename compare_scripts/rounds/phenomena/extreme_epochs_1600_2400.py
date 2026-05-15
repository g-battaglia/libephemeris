#!/usr/bin/env python3
"""Round 148: Pheno at extreme dates (far past/future).

Compare pheno_ut results for all major planets at dates spanning
1600 CE to 2400 CE, testing phase angle, phase, elongation, diameter, magnitude.
"""

from __future__ import annotations
import sys, os
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
BODIES = {
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
}
PHENO_LABELS = ["phase_angle", "phase", "elongation", "diameter", "magnitude"]

# Tolerances
TOL = {
    "phase_angle": 60.0,  # arcsec
    "phase": 0.01,  # unitless fraction
    "elongation": 60.0,  # arcsec
    "diameter": 0.5,  # arcsec
    "magnitude": 0.5,  # mag
}

test_dates = []
for year in [
    1600,
    1700,
    1800,
    1850,
    1900,
    1950,
    1980,
    2000,
    2010,
    2020,
    2024,
    2025,
    2030,
    2050,
    2100,
    2150,
    2200,
    2300,
    2400,
]:
    for month in [1, 7]:
        jd = swe.julday(year, month, 1, 12.0)
        test_dates.append((f"{year}-{month:02d}", jd))

passed = failed = errors = total = 0
failures = []

print(f"Round 148: Pheno at Extreme Dates")
print(
    f"Testing {len(BODIES)} bodies x {len(test_dates)} dates = {len(BODIES) * len(test_dates)} combos"
)
print("=" * 90)

for label, jd in test_dates:
    for body, bname in BODIES.items():
        try:
            se_r = swe.pheno_ut(jd, body, FLG_SPEED)
            le_r = ephem.pheno_ut(jd, body, FLG_SPEED)
            le_data = le_r[0]  # LE returns (tuple, retflag)

            for i, pname in enumerate(PHENO_LABELS):
                total += 1
                se_val = se_r[i]
                le_val = le_data[i]

                if pname in ("phase_angle", "elongation"):
                    diff = abs(le_val - se_val) * 3600.0
                    tol = TOL[pname]
                    unit = '"'
                elif pname == "diameter":
                    diff = abs(le_val - se_val) * 3600.0
                    tol = TOL[pname]
                    unit = '"'
                else:
                    diff = abs(le_val - se_val)
                    tol = TOL[pname]
                    unit = ""

                if diff <= tol:
                    passed += 1
                else:
                    failed += 1
                    msg = f"  FAIL {label} {bname} {pname}: SE={se_val:.8f} LE={le_val:.8f} diff={diff:.4f}{unit}"
                    failures.append(msg)
                    if len(failures) <= 30:
                        print(msg)
        except Exception as e:
            errors += 1
            if "not found" not in str(e):
                print(f"  ERROR {label} {bname}: {e}")

print()
print("=" * 90)
print(
    f"Results: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
)
if failures:
    print(f"\nTotal failures: {len(failures)}")
    # Categorize
    cats = {}
    for f in failures:
        for bname in BODIES.values():
            if bname in f:
                cats[bname] = cats.get(bname, 0) + 1
    for cat, count in sorted(cats.items(), key=lambda x: -x[1]):
        print(f"  {cat}: {count}")
else:
    print("\nAll tests passed!")
