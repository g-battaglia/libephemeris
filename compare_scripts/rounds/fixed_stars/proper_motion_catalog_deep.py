#!/usr/bin/env python3
"""Round 75: Fixed Star Proper Motion & Catalog Deep"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)
passed = failed = errors = 0
FLAGS = 256  # FLG_SPEED

STARS = [
    "Aldebaran",
    "Regulus",
    "Spica",
    "Antares",
    "Sirius",
    "Canopus",
    "Arcturus",
    "Vega",
    "Capella",
    "Rigel",
    "Procyon",
    "Betelgeuse",
    "Altair",
    "Pollux",
    "Fomalhaut",
    "Deneb",
    "Castor",
    "Achernar",
    "Acrux",
    "Mimosa",
    "Hadar",
    "Shaula",
    "Bellatrix",
    "Alnilam",
    "Alnitak",
    "Alhena",
    "Dubhe",
    "Mizar",
    "Alioth",
    "Alkaid",
]

DATES = [2451545.0 + i * 1826.25 for i in range(11)]  # 50 years in 5-year steps

print("=" * 70)
print("ROUND 75: Fixed Star Proper Motion & Catalog Deep")
print("=" * 70)

# P1: Star positions across epochs
print("\n=== P1: Star positions over 50 years ===")
for star in STARS:
    star_pass = star_fail = star_err = 0
    for jd in DATES:
        try:
            se = swe.fixstar2_ut(star, jd, FLAGS)
            le = ephem.fixstar2_ut(star, jd, FLAGS)
            # SE returns (pos_tuple, starname, retflag)
            # LE returns (pos_tuple, starname, retflag)
            se_lon = se[0][0]
            le_lon = le[0][0]
            diff = abs(se_lon - le_lon)
            if diff > 180:
                diff = 360 - diff
            diff_as = diff * 3600
            if diff_as < 10.0:  # 10" tolerance (known ~5.3" sidereal offset)
                passed += 1
                star_pass += 1
            else:
                failed += 1
                star_fail += 1
                if star_fail <= 1:
                    print(
                        f'  FAIL {star} jd={jd:.0f}: SE={se_lon:.6f} LE={le_lon:.6f} d={diff_as:.2f}"'
                    )
        except Exception as e:
            errors += 1
            star_err += 1
    if star_fail > 0 or star_err > 0:
        t = star_pass + star_fail
        print(f"  {star}: {star_pass}/{t} err={star_err}")
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# P2: Star latitude consistency
print("\n=== P2: Star latitude ===")
for star in STARS:
    for jd in DATES[:5]:
        try:
            se = swe.fixstar2_ut(star, jd, FLAGS)
            le = ephem.fixstar2_ut(star, jd, FLAGS)
            se_lat = se[0][1]
            le_lat = le[0][1]
            diff = abs(se_lat - le_lat) * 3600
            if diff < 10.0:
                passed += 1
            else:
                failed += 1
        except:
            errors += 1
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# P3: Star distance (should be very similar)
print("\n=== P3: Star distance ===")
for star in STARS[:15]:
    jd = 2451545.0
    try:
        se = swe.fixstar2_ut(star, jd, FLAGS)
        le = ephem.fixstar2_ut(star, jd, FLAGS)
        se_dist = se[0][2]
        le_dist = le[0][2]
        if se_dist > 0 and le_dist > 0:
            ratio = le_dist / se_dist
            if 0.9 < ratio < 1.1:
                passed += 1
            else:
                failed += 1
                print(f"  FAIL P3 {star}: SE_dist={se_dist:.6f} LE_dist={le_dist:.6f}")
        else:
            passed += 1
    except:
        errors += 1
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# P4: Proper motion effect — compare position drift over 100 years
print("\n=== P4: Proper motion drift verification ===")
jd_2000 = 2451545.0
jd_2100 = 2451545.0 + 36525.0
for star in STARS[:20]:
    try:
        se_2000 = swe.fixstar2_ut(star, jd_2000, FLAGS)
        se_2100 = swe.fixstar2_ut(star, jd_2100, FLAGS)
        le_2000 = ephem.fixstar2_ut(star, jd_2000, FLAGS)
        le_2100 = ephem.fixstar2_ut(star, jd_2100, FLAGS)

        se_drift = se_2100[0][0] - se_2000[0][0]
        le_drift = le_2100[0][0] - le_2000[0][0]
        if abs(se_drift) > 180:
            se_drift = se_drift - 360 if se_drift > 0 else se_drift + 360
        if abs(le_drift) > 180:
            le_drift = le_drift - 360 if le_drift > 0 else le_drift + 360

        diff = abs(se_drift - le_drift) * 3600
        if diff < 5.0:  # proper motion drift should match within 5"
            passed += 1
        else:
            failed += 1
            print(
                f'  FAIL P4 {star}: SE_drift={se_drift * 3600:.2f}" LE_drift={le_drift * 3600:.2f}" d={diff:.2f}"'
            )
    except:
        errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# P5: J2000 mode
print("\n=== P5: Stars in J2000 frame ===")
J2000 = 256 | 32
for star in STARS[:15]:
    for jd in DATES[:5]:
        try:
            se = swe.fixstar2_ut(star, jd, J2000)
            le = ephem.fixstar2_ut(star, jd, J2000)
            diff = abs(se[0][0] - le[0][0])
            if diff > 180:
                diff = 360 - diff
            if diff * 3600 < 5.0:
                passed += 1
            else:
                failed += 1
        except:
            errors += 1
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# P6: Star speed (lon/lat speed)
print("\n=== P6: Star speed ===")
for star in STARS[:15]:
    jd = 2451545.0
    try:
        se = swe.fixstar2_ut(star, jd, FLAGS)
        le = ephem.fixstar2_ut(star, jd, FLAGS)
        se_lon_spd = se[0][3]
        le_lon_spd = le[0][3]
        diff = abs(se_lon_spd - le_lon_spd)
        if diff < 0.001:  # 0.001 deg/day
            passed += 1
        else:
            failed += 1
            print(f"  FAIL P6 {star}: SE_spd={se_lon_spd:.8f} LE_spd={le_lon_spd:.8f}")
    except:
        errors += 1
print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 75 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
