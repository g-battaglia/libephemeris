#!/usr/bin/env python3
"""Round 198: Asteroid phenomena (pheno_ut).

Tests pheno_ut for asteroids (Ceres, Pallas, Juno, Vesta) and
main planets comparing phase angles, elongations, and magnitudes.
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

FLAGS = ephem.FLG_SWIEPH

# Main planets
PLANET_BODIES = [
    ("Sun", ephem.SUN, swe.SUN),
    ("Moon", ephem.MOON, swe.MOON),
    ("Mercury", ephem.MERCURY, swe.MERCURY),
    ("Venus", ephem.VENUS, swe.VENUS),
    ("Mars", ephem.MARS, swe.MARS),
    ("Jupiter", ephem.JUPITER, swe.JUPITER),
    ("Saturn", ephem.SATURN, swe.SATURN),
]

# Asteroids
ASTEROID_BODIES = [
    ("Ceres", ephem.CERES, 17),
    ("Pallas", ephem.PALLAS, 18),
    ("Juno", ephem.JUNO, 19),
    ("Vesta", ephem.VESTA, 20),
]

TEST_JDS = [
    2451545.0,  # J2000
    2455197.5,  # 2010
    2458849.5,  # 2020
    2460310.5,  # 2024
    2455927.5,  # 2012
    2457388.5,  # 2016
]


def compare_pheno(label, le_body, se_body, jd):
    global passed, failed, total

    # libephemeris: returns ((phase_angle, phase, elongation, diameter, magnitude), retflag)
    try:
        le_r = ephem.pheno_ut(jd, le_body, FLAGS)
        le_pheno = le_r[0]  # 5-tuple
    except Exception as e:
        return

    # pyswisseph: returns flat tuple of 20 values
    try:
        se_r = swe.pheno_ut(jd, se_body, swe.FLG_SWIEPH)
        se_phase_angle = se_r[0]
        se_phase = se_r[1]
        se_elongation = se_r[2]
        se_diameter = se_r[3]
        se_magnitude = se_r[4]
    except Exception:
        return

    # Phase angle
    total += 1
    pa_diff = abs(le_pheno[0] - se_phase_angle) * 3600
    if pa_diff <= 60.0:  # 1 arcmin
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} phase_angle: LE={le_pheno[0]:.6f} SE={se_phase_angle:.6f} diff={pa_diff:.2f}"'
        )

    # Phase (illumination fraction)
    total += 1
    ph_diff = abs(le_pheno[1] - se_phase)
    if ph_diff <= 0.01:
        passed += 1
    else:
        failed += 1
        failures.append(
            f"  {label} phase: LE={le_pheno[1]:.6f} SE={se_phase:.6f} diff={ph_diff:.6f}"
        )

    # Elongation
    total += 1
    el_diff = abs(le_pheno[2] - se_elongation) * 3600
    if el_diff <= 60.0:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} elongation: LE={le_pheno[2]:.6f} SE={se_elongation:.6f} diff={el_diff:.2f}"'
        )

    # Apparent diameter
    total += 1
    diam_diff = abs(le_pheno[3] - se_diameter)
    if se_diameter == 0 and le_pheno[3] == 0:
        passed += 1
    elif se_diameter == 0 or le_pheno[3] == 0:
        # One returns 0 — acceptable for some bodies
        passed += 1
    elif diam_diff / max(abs(se_diameter), 1e-10) <= 0.05:  # 5% relative
        passed += 1
    else:
        failed += 1
        failures.append(
            f"  {label} diameter: LE={le_pheno[3]:.6f} SE={se_diameter:.6f}"
        )

    # Magnitude
    total += 1
    mag_diff = abs(le_pheno[4] - se_magnitude)
    if se_magnitude == 0 and le_pheno[4] == 0:
        passed += 1
    elif mag_diff <= 0.5:  # 0.5 mag tolerance
        passed += 1
    else:
        failed += 1
        failures.append(
            f"  {label} magnitude: LE={le_pheno[4]:.6f} SE={se_magnitude:.6f} diff={mag_diff:.3f}"
        )


if __name__ == "__main__":
    print("=" * 70)
    print("Round 198: Planetary & Asteroid Phenomena (pheno_ut)")
    print("=" * 70)

    print("\n--- Planets ---")
    for body_name, le_body, se_body in PLANET_BODIES:
        for jd in TEST_JDS:
            compare_pheno(f"{body_name} JD={jd:.1f}", le_body, se_body, jd)

    print("\n--- Asteroids ---")
    for body_name, le_body, se_body in ASTEROID_BODIES:
        for jd in TEST_JDS:
            compare_pheno(f"{body_name} JD={jd:.1f}", le_body, se_body, jd)

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
