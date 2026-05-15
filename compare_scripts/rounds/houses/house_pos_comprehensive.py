#!/usr/bin/env python3
"""Round 64: House Position (house_pos) Comprehensive

Compare house_pos() for various house systems, latitudes, and planet positions.
house_pos returns the house position (1.0-12.999) for a given ARMC, geo lat, obliquity, and planet lon/lat.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = failed = errors = 0


def se_hsys(ch):
    return ch.encode("ascii") if isinstance(ch, str) else ch


def le_hsys(ch):
    return ord(ch) if isinstance(ch, str) else ch


HOUSE_SYSTEMS = ["P", "K", "O", "R", "C", "E", "W", "B", "M", "A"]
# P=Placidus, K=Koch, O=Porphyry, R=Regiomontanus, C=Campanus
# E=Equal, W=Whole Sign, B=Alcabitius, M=Morinus, A=Equal from Asc

LATITUDES = [0, 10, 23.44, 30, 40, 45, 50, 55, 60, 65]
ARMCS = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
OBLIQUITY = 23.4393  # approximate J2000

print("=" * 70)
print("ROUND 64: House Position (house_pos) Comprehensive")
print("=" * 70)

# ============================================================
# P1: house_pos for Placidus across ARMC/lat grid
# ============================================================
print("\n=== P1: Placidus house_pos ===")
for lat in LATITUDES:
    for armc in ARMCS:
        for planet_lon in [0, 45, 90, 135, 180, 225, 270, 315]:
            planet_lat = 0.0
            try:
                se_hp = swe.house_pos(
                    armc, lat, OBLIQUITY, (planet_lon, planet_lat), se_hsys("P")
                )
                le_hp = ephem.house_pos(
                    armc, lat, OBLIQUITY, le_hsys("P"), planet_lon, planet_lat
                )
                diff = abs(se_hp - le_hp)
                if diff > 6:
                    diff = 12 - diff  # wrap around house 12->1
                if diff < 0.01:
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL P1 lat={lat} armc={armc} lon={planet_lon}: SE={se_hp:.4f} LE={le_hp:.4f} diff={diff:.4f}"
                    )
            except Exception as e:
                errors += 1
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: house_pos for Koch
# ============================================================
print("\n=== P2: Koch house_pos ===")
for lat in LATITUDES:
    for armc in ARMCS:
        for planet_lon in [0, 60, 120, 180, 240, 300]:
            try:
                se_hp = swe.house_pos(
                    armc, lat, OBLIQUITY, (planet_lon, 0.0), se_hsys("K")
                )
                le_hp = ephem.house_pos(
                    armc, lat, OBLIQUITY, le_hsys("K"), planet_lon, 0.0
                )
                diff = abs(se_hp - le_hp)
                if diff > 6:
                    diff = 12 - diff
                if diff < 0.01:
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL P2 lat={lat} armc={armc} lon={planet_lon}: SE={se_hp:.4f} LE={le_hp:.4f} diff={diff:.4f}"
                    )
            except Exception as e:
                errors += 1
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: house_pos for all systems at specific conditions
# ============================================================
print("\n=== P3: All house systems at lat=45, armc=90 ===")
lat = 45.0
armc = 90.0
for hsys in HOUSE_SYSTEMS:
    for planet_lon in range(0, 360, 15):
        try:
            se_hp = swe.house_pos(
                armc, lat, OBLIQUITY, (planet_lon, 0.0), se_hsys(hsys)
            )
            le_hp = ephem.house_pos(
                armc, lat, OBLIQUITY, le_hsys(hsys), planet_lon, 0.0
            )
            diff = abs(se_hp - le_hp)
            if diff > 6:
                diff = 12 - diff
            tol = 0.05 if hsys in ("B",) else 0.01  # Alcabitius wider tolerance
            if diff < tol:
                passed += 1
            else:
                failed += 1
                print(
                    f"  FAIL P3 {hsys} lon={planet_lon}: SE={se_hp:.4f} LE={le_hp:.4f} diff={diff:.4f}"
                )
        except Exception as e:
            errors += 1
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: house_pos with non-zero planet latitude
# ============================================================
print("\n=== P4: house_pos with planet latitude ===")
for lat in [30, 45, 55]:
    for armc in [0, 90, 180, 270]:
        for planet_lon in [0, 90, 180, 270]:
            for planet_lat in [-5, -2, 0, 2, 5]:
                for hsys in ["P", "K", "R", "C"]:
                    try:
                        se_hp = swe.house_pos(
                            armc,
                            lat,
                            OBLIQUITY,
                            (planet_lon, planet_lat),
                            se_hsys(hsys),
                        )
                        le_hp = ephem.house_pos(
                            armc, lat, OBLIQUITY, le_hsys(hsys), planet_lon, planet_lat
                        )
                        diff = abs(se_hp - le_hp)
                        if diff > 6:
                            diff = 12 - diff
                        if diff < 0.02:
                            passed += 1
                        else:
                            failed += 1
                            print(
                                f"  FAIL P4 {hsys} lat={lat} armc={armc} lon={planet_lon} plat={planet_lat}: SE={se_hp:.4f} LE={le_hp:.4f}"
                            )
                    except Exception as e:
                        errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: house_pos with real planet data
# ============================================================
print("\n=== P5: house_pos with real chart data ===")
FLAGS = 256
jd_test = 2451545.0  # J2000

for lat in [0, 30, 45, 55, 65]:
    for lon in [0, 30, 90]:
        # Get houses first
        try:
            le_houses = ephem.houses_ex(jd_test, lat, lon, le_hsys("P"), 0)
            le_cusps = le_houses[0]
            le_ascmc = le_houses[1]
            armc = le_ascmc[2]
            obl = ephem.calc_ut(jd_test, -1, 0)[0][1]  # mean obliquity

            # Get planet positions and compute house_pos
            for body in [0, 1, 2, 3, 4, 5]:
                planet = ephem.calc_ut(jd_test, body, FLAGS)[0]
                plon, plat = planet[0], planet[1]

                for hsys in ["P", "K", "O", "R", "E"]:
                    try:
                        se_hp = swe.house_pos(
                            armc, lat, obl, (plon, plat), se_hsys(hsys)
                        )
                        le_hp = ephem.house_pos(
                            armc, lat, obl, le_hsys(hsys), plon, plat
                        )
                        diff = abs(se_hp - le_hp)
                        if diff > 6:
                            diff = 12 - diff
                        if diff < 0.02:
                            passed += 1
                        else:
                            failed += 1
                            print(
                                f"  FAIL P5 {hsys} lat={lat} body={body}: SE={se_hp:.4f} LE={le_hp:.4f}"
                            )
                    except:
                        errors += 1
        except Exception as e:
            errors += 1
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P6: house_pos range validation (should be 1.0 to 12.999)
# ============================================================
print("\n=== P6: house_pos range validation ===")
for lat in [0, 30, 60]:
    for armc in range(0, 360, 30):
        for plon in range(0, 360, 10):
            try:
                hp = ephem.house_pos(armc, lat, OBLIQUITY, le_hsys("P"), plon, 0.0)
                if 1.0 <= hp < 13.0:
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL P6 lat={lat} armc={armc} lon={plon}: hp={hp:.4f} (out of range)"
                    )
            except:
                errors += 1
print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P7: Porphyry house_pos (should be simple interpolation)
# ============================================================
print("\n=== P7: Porphyry house_pos ===")
for lat in [0, 20, 40, 60]:
    for armc in range(0, 360, 45):
        for plon in range(0, 360, 30):
            try:
                se_hp = swe.house_pos(armc, lat, OBLIQUITY, (plon, 0.0), se_hsys("O"))
                le_hp = ephem.house_pos(
                    armc, lat, OBLIQUITY, le_hsys("O"), plon, 0.0
                )
                diff = abs(se_hp - le_hp)
                if diff > 6:
                    diff = 12 - diff
                if diff < 0.005:
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL P7 lat={lat} armc={armc} lon={plon}: SE={se_hp:.4f} LE={le_hp:.4f}"
                    )
            except:
                errors += 1
print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P8: Equal house_pos
# ============================================================
print("\n=== P8: Equal house_pos ===")
for lat in [0, 30, 45, 60]:
    for armc in range(0, 360, 30):
        for plon in range(0, 360, 15):
            try:
                se_hp = swe.house_pos(armc, lat, OBLIQUITY, (plon, 0.0), se_hsys("E"))
                le_hp = ephem.house_pos(
                    armc, lat, OBLIQUITY, le_hsys("E"), plon, 0.0
                )
                diff = abs(se_hp - le_hp)
                if diff > 6:
                    diff = 12 - diff
                if diff < 0.005:
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL P8 lat={lat} armc={armc} lon={plon}: SE={se_hp:.4f} LE={le_hp:.4f}"
                    )
            except:
                errors += 1
print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")

print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 64 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
if failed > 0:
    print(f"\n--- FAILURES ({failed}) ---")
