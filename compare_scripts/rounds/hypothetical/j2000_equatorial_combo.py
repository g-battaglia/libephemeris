#!/usr/bin/env python3
"""Round 204: Hypothetical bodies J2000+EQUATORIAL combo.

Tests Uranian/hypothetical bodies with combined J2000 and EQUATORIAL flags.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
total = 0
failures = []

URANIANS = [
    ("Cupido", ephem.CUPIDO, swe.CUPIDO),
    ("Hades", ephem.HADES, swe.HADES),
    ("Zeus", ephem.ZEUS, swe.ZEUS),
    ("Kronos", ephem.KRONOS, swe.KRONOS),
    ("Apollon", ephem.APOLLON, swe.APOLLON),
    ("Admetos", ephem.ADMETOS, swe.ADMETOS),
    ("Vulkanus", ephem.VULKANUS, swe.VULKANUS),
    ("Poseidon", ephem.POSEIDON, swe.POSEIDON),
]

FLAG_COMBOS = [
    ("default+SPEED", ephem.FLG_SWIEPH | ephem.FLG_SPEED),
    ("J2000+SPEED", ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_J2000),
    (
        "EQUATORIAL+SPEED",
        ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_EQUATORIAL,
    ),
    (
        "J2000+EQ+SPEED",
        ephem.FLG_SWIEPH
        | ephem.FLG_SPEED
        | ephem.FLG_J2000
        | ephem.FLG_EQUATORIAL,
    ),
    ("NONUT+SPEED", ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_NONUT),
    ("HELCTR+SPEED", ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_HELCTR),
]

TEST_JDS = [2451545.0, 2455197.5, 2458849.5, 2460310.5]


def compare(label, le_body, se_body, jd, flags, tol=60.0):
    global passed, failed, total

    try:
        le_r = ephem.calc_ut(jd, le_body, flags)
        se_r = swe.calc_ut(jd, se_body, flags)
    except Exception:
        return

    # Longitude/RA
    total += 1
    val_diff = abs(le_r[0][0] - se_r[0][0])
    if val_diff > 180:
        val_diff = 360 - val_diff
    val_as = val_diff * 3600
    if val_as <= tol:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} [0]: diff={val_as:.2f}" (tol {tol}")')

    # Latitude/Dec
    total += 1
    lat_as = abs(le_r[0][1] - se_r[0][1]) * 3600
    if lat_as <= tol:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} [1]: diff={lat_as:.2f}"')


if __name__ == "__main__":
    print("=" * 70)
    print("Round 204: Hypothetical Bodies J2000+EQUATORIAL Combo")
    print("=" * 70)

    for name, le_b, se_b in URANIANS:
        print(f"\n--- {name} ---")
        for fname, flags in FLAG_COMBOS:
            for jd in TEST_JDS:
                label = f"{name} {fname} JD={jd:.1f}"
                compare(label, le_b, se_b, jd, flags)

    print(f"\n{'=' * 70}")
    print(
        f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
    )
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
        if len(failures) > 30:
            print(f"  ... and {len(failures) - 30} more")
