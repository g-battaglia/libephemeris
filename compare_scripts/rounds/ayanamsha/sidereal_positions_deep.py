#!/usr/bin/env python3
"""Round 70: Sidereal Mode All Planets Deep Sweep"""

from __future__ import annotations
import sys, os
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)
passed = failed = errors = 0
FLAGS = 256 | 65536  # FLG_SPEED | FLG_SIDEREAL

NAMES = {
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
DATES = [2451545.0 + i * 365.25 for i in range(25)]  # 25 years

# Sidereal modes to test
SID_MODES = [0, 1, 3, 5, 7, 22, 27, 30]
# 0=Fagan/Bradley, 1=Lahiri, 3=Raman, 5=Krishnamurti, 7=Yukteshwar
# 22=Galactic Center 0Sag, 27=True Pushya, 30=Galactic Eq IAU1958

print("=" * 70)
print("ROUND 70: Sidereal Mode All Planets Deep Sweep")
print("=" * 70)

for sid_mode in SID_MODES:
    swe.set_sid_mode(sid_mode)
    ephem.set_sid_mode(sid_mode)
    mode_pass = mode_fail = 0

    for jd in DATES:
        for body in range(10):
            try:
                se = swe.calc_ut(jd, body, FLAGS)
                le = ephem.calc_ut(jd, body, FLAGS)
                # Longitude
                diff_lon = abs(se[0][0] - le[0][0])
                if diff_lon > 180:
                    diff_lon = 360 - diff_lon
                diff_arcsec = diff_lon * 3600
                if diff_arcsec < 20:  # Known ~14" sidereal offset
                    passed += 1
                    mode_pass += 1
                else:
                    failed += 1
                    mode_fail += 1
                    if mode_fail <= 3:
                        print(
                            f"  FAIL mode={sid_mode} {NAMES[body]} jd={jd:.0f}: "
                            f'SE={se[0][0]:.6f} LE={le[0][0]:.6f} diff={diff_arcsec:.1f}"'
                        )
            except Exception as e:
                errors += 1
    print(f"  Mode {sid_mode}: {mode_pass}/{mode_pass + mode_fail} passed")

# Reset sidereal mode
swe.set_sid_mode(0)
ephem.set_sid_mode(0)

print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 70 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
