#!/usr/bin/env python3
"""
Diagnostic: Compare Moon magnitude formula across phase angles.

Finds epochs where Moon is at various phase angles and compares
SE vs LE magnitude output. Focus on high phase angles (>120°).
"""

from __future__ import annotations

import math
import os
import sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import libephemeris as ephem
from libephemeris.constants import *

_EPHE_PATH = os.path.join(os.path.dirname(__file__), "..", "swisseph", "ephe")
swe.set_ephe_path(_EPHE_PATH)

# Scan a year of daily Moon positions to find various phase angles
print("Phase Angle | SE Mag  | LE Mag  | Diff    | SE Phase Angle")
print("-" * 70)

# Start from J2000 and scan daily for 2 years to get good phase angle coverage
jd_start = 2451545.0  # J2000
results = []

for i in range(730):  # 2 years of daily samples
    jd = jd_start + i * 0.5  # every 12 hours for finer sampling
    try:
        se_attr = swe.pheno_ut(jd, MOON, 0)
        le_attr, _ = ephem.pheno_ut(jd, MOON, 0)
    except Exception:
        continue

    se_pa = se_attr[0]  # phase angle
    se_mag = se_attr[4]
    le_pa = le_attr[0]
    le_mag = le_attr[4]
    mag_diff = le_mag - se_mag

    results.append((se_pa, se_mag, le_mag, mag_diff, le_pa))

# Sort by phase angle
results.sort(key=lambda x: x[0])

# Print sampled results at various phase angle ranges
bins = [
    (0, 10),
    (10, 20),
    (20, 40),
    (40, 60),
    (60, 80),
    (80, 100),
    (100, 120),
    (120, 130),
    (130, 140),
    (140, 150),
    (150, 155),
    (155, 160),
    (160, 165),
    (165, 170),
    (170, 175),
    (175, 180),
]

for lo, hi in bins:
    subset = [r for r in results if lo <= r[0] < hi]
    if subset:
        # Pick the one with max abs diff
        worst = max(subset, key=lambda x: abs(x[3]))
        print(
            f"  {worst[0]:6.1f}°   | {worst[1]:7.2f} | {worst[2]:7.2f} | {worst[3]:+7.2f} | {worst[4]:6.1f}°"
        )

print()
print("=" * 70)
print("All samples with phase angle > 140°:")
print("=" * 70)
print(f"{'SE PA':>8} | {'SE Mag':>8} | {'LE Mag':>8} | {'Diff':>8} | {'LE PA':>8}")
print("-" * 55)

high_pa = [r for r in results if r[0] > 140]
for r in high_pa:
    flag = " <<<" if abs(r[3]) > 0.5 else ""
    print(
        f"  {r[0]:6.1f}° | {r[1]:8.2f} | {r[2]:8.2f} | {r[3]:+8.2f} | {r[4]:6.1f}°{flag}"
    )

print()
print("=" * 70)
print("SE Moon magnitude formula reverse-engineering:")
print("For the same phase angles, what does SE actually compute?")
print("=" * 70)

# Try to reverse-engineer SE's formula
# SE formula is likely: V = -12.73 + 0.026*|alpha| + 4e-9*|alpha|^4
# but may cap or use a different formula at high phase angles
# Let's check our formula vs SE at high phase angles
for r in high_pa:
    alpha = r[0]
    # Our formula
    our_darkening = 0.026 * alpha + 4.0e-9 * alpha**4
    our_mag = -12.73 + our_darkening  # ignoring distance correction

    # What darkening does SE use? (removing distance correction)
    # We can't easily extract this without knowing the distance, but let's
    # compare the difference pattern
    print(
        f"  PA={alpha:6.1f}°: SE_mag={r[1]:8.2f}  our_darkening={our_darkening:8.4f}  alpha^4_term={4e-9 * alpha**4:8.4f}"
    )
