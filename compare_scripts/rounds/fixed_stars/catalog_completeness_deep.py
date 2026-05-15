#!/usr/bin/env python3
"""Round 201: Fixed star catalog completeness.

Tests all stars in LE's catalog against SE, checking which stars are
found by both and verifying position agreement.
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
skipped = 0
failures = []

FLAGS = ephem.FLG_SWIEPH | ephem.FLG_SPEED
JD = 2451545.0  # J2000

# Comprehensive star list
STARS = [
    "Aldebaran",
    "Regulus",
    "Spica",
    "Antares",
    "Sirius",
    "Fomalhaut",
    "Vega",
    "Capella",
    "Rigel",
    "Betelgeuse",
    "Procyon",
    "Pollux",
    "Castor",
    "Deneb",
    "Altair",
    "Arcturus",
    "Canopus",
    "Achernar",
    "Hamal",
    "Polaris",
    "Algol",
    "Alcyone",
    "Alphecca",
    "Zubenelgenubi",
    "Zubeneschamali",
    "Unukalhai",
    "Ras Alhague",
    "Scheat",
    "Markab",
    "Alpheratz",
    "Mirach",
    "Almach",
    "Menkar",
    "Bellatrix",
    "Mintaka",
    "Alnilam",
    "Alnitak",
    "Saiph",
    "Adhara",
    "Wezen",
    "Aludra",
    "Mirzam",
    "Naos",
    "Alphard",
    "Denebola",
    "Vindemiatrix",
    "Algorab",
    "Acrux",
    "Mimosa",
    "Agena",
    "Toliman",
    "Izar",
    "Kochab",
    "Cor Caroli",
    "Algenib",
    "Mira",
    "Acubens",
    "Zaniah",
    "Khambalia",
    "Syrma",
    "Foramen",
    "Vertex",
    "Aculeus",
    "Acumen",
    "Sinistra",
    "Spiculum",
    "Nunki",
    "Ascella",
    "Vega",
    "Altair",
    "Albireo",
    "Tarazed",
    "Giedi",
    "Dabih",
    "Nashira",
    "Deneb Algedi",
    "Sadalsuud",
    "Sadalmelik",
    "Fomalhaut",
    "Scheat",
    "Markab",
]

# Remove duplicates
STARS = list(dict.fromkeys(STARS))


def test_star_catalog():
    global passed, failed, total, skipped

    print("=" * 70)
    print("Round 201: Fixed Star Catalog Completeness")
    print("=" * 70)

    found_both = 0
    le_only = 0
    se_only = 0

    for star in STARS:
        le_ok = False
        se_ok = False

        try:
            le_r = ephem.fixstar2_ut(star, JD, FLAGS)
            le_lon = le_r[0][0]
            le_lat = le_r[0][1]
            le_ok = True
        except Exception:
            pass

        try:
            se_r = swe.fixstar2(star, JD, swe.FLG_SWIEPH | swe.FLG_SPEED)
            se_lon = se_r[0][0]
            se_lat = se_r[0][1]
            se_ok = True
        except Exception:
            pass

        if le_ok and se_ok:
            found_both += 1
            total += 1
            lon_diff = abs(le_lon - se_lon)
            if lon_diff > 180:
                lon_diff = 360 - lon_diff
            lon_as = lon_diff * 3600

            if lon_as <= 1.0:
                passed += 1
            else:
                failed += 1
                failures.append(f'  {star} LON: diff={lon_as:.4f}"')

            total += 1
            lat_as = abs(le_lat - se_lat) * 3600
            if lat_as <= 1.0:
                passed += 1
            else:
                failed += 1
                failures.append(f'  {star} LAT: diff={lat_as:.4f}"')
        elif le_ok and not se_ok:
            le_only += 1
        elif not le_ok and se_ok:
            se_only += 1

    print(f"\n  Stars found by both: {found_both}")
    print(f"  Stars in LE only: {le_only}")
    print(f"  Stars in SE only: {se_only}")

    # Test with various flag combinations for the first 20 stars
    print("\n--- Flag Combinations ---")
    flag_combos = [
        (
            "J2000",
            FLAGS | ephem.FLG_J2000,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_J2000,
        ),
        (
            "NONUT",
            FLAGS | ephem.FLG_NONUT,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_NONUT,
        ),
        (
            "EQUATORIAL",
            FLAGS | ephem.FLG_EQUATORIAL,
            swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL,
        ),
    ]

    for star in STARS[:20]:
        for fname, le_f, se_f in flag_combos:
            try:
                le_r = ephem.fixstar2_ut(star, JD, le_f)
                se_r = swe.fixstar2(star, JD, se_f)
            except Exception:
                continue

            total += 1
            lon_diff = abs(le_r[0][0] - se_r[0][0])
            if lon_diff > 180:
                lon_diff = 360 - lon_diff
            lon_as = lon_diff * 3600
            if lon_as <= 1.0:
                passed += 1
            else:
                failed += 1
                failures.append(f'  {star} {fname}: diff={lon_as:.4f}"')


if __name__ == "__main__":
    test_star_catalog()

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
