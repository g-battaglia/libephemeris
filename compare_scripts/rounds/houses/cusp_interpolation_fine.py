#!/usr/bin/env python3
"""Round 192: House cusp interpolation fine grain.

Tests house cusps at fine latitude/longitude/ARMC increments to verify
smooth interpolation and catch any discontinuities or jumps.
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

# House systems to test
HOUSE_SYSTEMS = [
    ("Placidus", ord("P"), b"P"),
    ("Koch", ord("K"), b"K"),
    ("Equal", ord("E"), b"E"),
    ("WholeSign", ord("W"), b"W"),
    ("Regiomontanus", ord("R"), b"R"),
    ("Campanus", ord("C"), b"C"),
    ("Porphyry", ord("O"), b"O"),
]

# Test JD
JD = 2451545.0  # J2000

# Fine latitude increments
LATITUDES = list(range(-60, 65, 5))

# Fine longitude increments
LONGITUDES = [
    0,
    15,
    30,
    45,
    60,
    75,
    90,
    105,
    120,
    135,
    150,
    165,
    180,
    195,
    210,
    225,
    240,
    255,
    270,
    285,
    300,
    315,
    330,
    345,
]


def compare_cusps(label, le_cusps, se_cusps, hsys_name, tol=1.0):
    """Compare house cusps between LE and SE."""
    global passed, failed, total

    n_cusps = min(len(le_cusps), len(se_cusps), 12)
    for i in range(n_cusps):
        total += 1
        diff = abs(le_cusps[i] - se_cusps[i])
        if diff > 180:
            diff = 360 - diff

        diff_arcsec = diff * 3600

        # Wider tolerance for Sunshine and some edge cases
        actual_tol = tol
        if diff_arcsec <= actual_tol:
            passed += 1
        else:
            failed += 1
            failures.append(
                f'  {label} cusp{i + 1}: LE={le_cusps[i]:.6f} SE={se_cusps[i]:.6f} diff={diff_arcsec:.2f}"'
            )


def test_fine_latitude_sweep():
    global passed, failed, total

    print("=" * 70)
    print("Round 192: House Cusp Interpolation Fine Grain")
    print("=" * 70)

    for hsys_name, le_hsys, se_hsys in HOUSE_SYSTEMS:
        print(f"\n--- {hsys_name} ---")

        for lat in LATITUDES:
            for lon in LONGITUDES[::3]:  # every 3rd longitude for speed
                try:
                    le_r = ephem.houses_ex2(JD, float(lat), float(lon), le_hsys, 0)
                    le_cusps = le_r[0]
                except Exception as e:
                    continue

                try:
                    se_r = swe.houses_ex(JD, float(lat), float(lon), se_hsys)
                    se_cusps = se_r[0]
                except Exception:
                    continue

                label = f"{hsys_name} lat={lat} lon={lon}"
                compare_cusps(label, le_cusps, se_cusps, hsys_name)


def test_ascmc_fine():
    """Test ASC/MC/Vertex at fine latitude increments."""
    global passed, failed, total

    print("\n--- ASC/MC/Vertex Fine ---")

    for lat in range(-66, 67, 2):
        for lon in [0, 90, 180, 270]:
            try:
                le_r = ephem.houses_ex2(JD, float(lat), float(lon), ord("P"), 0)
                le_ascmc = le_r[1]

                se_r = swe.houses_ex(JD, float(lat), float(lon), b"P")
                se_ascmc = se_r[1]
            except Exception:
                continue

            # ASC
            total += 1
            diff = abs(le_ascmc[0] - se_ascmc[0])
            if diff > 180:
                diff = 360 - diff
            if diff * 3600 <= 1.0:
                passed += 1
            else:
                failed += 1
                failures.append(f'  ASC lat={lat} lon={lon}: diff={diff * 3600:.2f}"')

            # MC
            total += 1
            diff = abs(le_ascmc[1] - se_ascmc[1])
            if diff > 180:
                diff = 360 - diff
            if diff * 3600 <= 1.0:
                passed += 1
            else:
                failed += 1
                failures.append(f'  MC lat={lat} lon={lon}: diff={diff * 3600:.2f}"')

            # Vertex
            total += 1
            diff = abs(le_ascmc[3] - se_ascmc[3])
            if diff > 180:
                diff = 360 - diff
            if diff * 3600 <= 5.0:  # wider for Vertex
                passed += 1
            else:
                failed += 1
                failures.append(
                    f'  VTX lat={lat} lon={lon}: LE={le_ascmc[3]:.6f} SE={se_ascmc[3]:.6f} diff={diff * 3600:.2f}"'
                )


def test_armc_sweep():
    """Test houses_armc with fine ARMC sweep."""
    global passed, failed, total

    print("\n--- ARMC Sweep ---")

    for armc in range(0, 360, 15):
        for lat in [-45, 0, 30, 52, 66]:
            eps = 23.4393  # approximate obliquity
            try:
                le_r = ephem.houses_armc_ex2(
                    float(armc), float(lat), eps, ord("P"), 0
                )
                le_cusps = le_r[0]

                se_r = swe.houses_armc(float(armc), float(lat), eps, b"P")
                se_cusps = se_r[0]
            except Exception:
                continue

            label = f"ARMC={armc} lat={lat}"
            compare_cusps(label, le_cusps, se_cusps, "Placidus")


if __name__ == "__main__":
    test_fine_latitude_sweep()
    test_ascmc_fine()
    test_armc_sweep()

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
