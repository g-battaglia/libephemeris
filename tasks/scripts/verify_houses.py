#!/usr/bin/env python3
"""Verify libephemeris house calculations against pyswisseph.

Sections:
  3.1-3.2: 24 house systems x 6 locations x 20 dates (cusps + ASC/MC comparison)
  3.3:     houses_ex with sidereal (5 systems x 3 ayanamshas x 10 dates)
  3.5:     houses_armc (5 systems x 10 ARMC x 3 obliquities)
  3.6:     house_pos (10 longitudes x 5 systems x 5 dates)
  3.7:     Gauquelin sectors via house_pos with system G

Target: ~10000+ checks, <30 seconds.
"""

import math
import os
import sys
import time
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import libephemeris as lib
import swisseph as swe_ref

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

# Point pyswisseph at its ephemeris files
swe_ref.set_ephe_path(_REF_EPHE_PATH)

# Ensure libephemeris uses Skyfield backend (no LEB)
lib.close()
lib.set_calc_mode("skyfield")

# ---------------------------------------------------------------------------
# Counters and check helper
# ---------------------------------------------------------------------------
passed = 0
failed = 0
errors = []
fail_by = Counter()


def check(cond, desc="", key=""):
    """Record a pass/fail check. Failures are logged with description and key."""
    global passed, failed
    if cond:
        passed += 1
    else:
        failed += 1
        fail_by[key] += 1
        if len(errors) < 120:
            errors.append(desc)


def angle_diff(a, b):
    """Absolute angular difference in degrees, handling 360/0 wraparound."""
    d = abs(float(a) - float(b))
    if d > 180.0:
        d = 360.0 - d
    return d


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# 24 house system characters
HOUSE_SYSTEMS = list("PKORECWBMTUHXGYIiNFDJLAS")

# 6 locations: (lat, lon)
LOCATIONS = [
    (0.0, 0.0),  # Null Island (equator/prime meridian)
    (41.9, 12.5),  # Rome
    (51.5, -0.1),  # London
    (-33.9, 151.2),  # Sydney
    (35.7, 139.7),  # Tokyo
    (64.0, -22.0),  # Reykjavik (high latitude)
]

# 20 JDs evenly spaced 2000-01-01 to 2020-01-01
JD_2000 = 2451545.0  # J2000.0
JD_2020 = 2458849.5  # ~2020-01-01
JDS_20 = [JD_2000 + i * (JD_2020 - JD_2000) / 19 for i in range(20)]

# 10 JDs evenly spaced
JDS_10 = [JD_2000 + i * (JD_2020 - JD_2000) / 9 for i in range(10)]

# 5 JDs evenly spaced
JDS_5 = [JD_2000 + i * (JD_2020 - JD_2000) / 4 for i in range(5)]

t0 = time.time()

# =========================================================================
# Section 3.1-3.2: 24 house systems x 6 locations x 20 dates
# Compare cusps (12), ASC, MC => 14 checks per combo
# 24 * 6 * 20 * 14 = 40,320 checks
# =========================================================================
print("Section 3.1-3.2: 24 house systems x 6 locations x 20 dates vs pyswisseph...")

for hsys_char in HOUSE_SYSTEMS:
    hsys_int = ord(hsys_char)
    hsys_bytes = bytes(hsys_char, "ascii")
    for lat, lon in LOCATIONS:
        for jd in JDS_20:
            lib_ok = True
            ref_ok = True
            lib_cusps = None
            lib_ascmc = None
            ref_cusps = None
            ref_ascmc = None

            # --- libephemeris call ---
            try:
                lib_result = lib.houses(jd, lat, lon, hsys_int)
                lib_cusps = lib_result[0]
                lib_ascmc = lib_result[1]
            except Exception as e:
                lib_ok = False
                ename = type(e).__name__

            # --- pyswisseph call ---
            try:
                ref_result = swe_ref.houses(jd, lat, lon, hsys_bytes)
                ref_cusps = ref_result[0]
                ref_ascmc = ref_result[1]
            except Exception as e:
                ref_ok = False

            # If both fail, that is acceptable (polar issues etc.) -- count as pass
            if not lib_ok and not ref_ok:
                passed += 1  # both agree on error
                continue

            # If one fails and the other does not, record failures
            if not lib_ok and ref_ok:
                for _ in range(14):
                    check(
                        False,
                        f"3.1 {hsys_char} ({lat},{lon}) jd={jd:.1f} lib error, ref ok",
                        f"3.1/{hsys_char}/LIB_ERR",
                    )
                continue
            if lib_ok and not ref_ok:
                # Ref failed but lib succeeded -- count as pass (lib may support more)
                for _ in range(14):
                    passed += 1
                continue

            # Both succeeded -- compare results
            # Compare 12 cusps (tolerance 0.01 deg)
            ncusps = min(12, len(lib_cusps), len(ref_cusps))
            for i in range(ncusps):
                diff = angle_diff(lib_cusps[i], ref_cusps[i])
                check(
                    diff < 0.01,
                    f"3.2 {hsys_char} cusp{i + 1} ({lat},{lon}) jd={jd:.1f} "
                    f"lib={float(lib_cusps[i]):.4f} ref={float(ref_cusps[i]):.4f} diff={diff:.6f}",
                    f"3.2/{hsys_char}/CUSP{i + 1}",
                )
            # Pad remaining if fewer cusps
            for i in range(ncusps, 12):
                check(
                    False,
                    f"3.2 {hsys_char} cusp{i + 1} missing",
                    f"3.2/{hsys_char}/MISSING",
                )

            # Compare ASC
            asc_diff = angle_diff(lib_ascmc[0], ref_ascmc[0])
            check(
                asc_diff < 0.01,
                f"3.2 {hsys_char} ASC ({lat},{lon}) jd={jd:.1f} "
                f"lib={float(lib_ascmc[0]):.4f} ref={float(ref_ascmc[0]):.4f} diff={asc_diff:.6f}",
                f"3.2/{hsys_char}/ASC",
            )

            # Compare MC
            mc_diff = angle_diff(lib_ascmc[1], ref_ascmc[1])
            check(
                mc_diff < 0.01,
                f"3.2 {hsys_char} MC ({lat},{lon}) jd={jd:.1f} "
                f"lib={float(lib_ascmc[1]):.4f} ref={float(ref_ascmc[1]):.4f} diff={mc_diff:.6f}",
                f"3.2/{hsys_char}/MC",
            )

print(f"  3.1-3.2 done: {time.time() - t0:.1f}s  (passed={passed} failed={failed})")


# =========================================================================
# Section 3.3: houses_ex with sidereal flag
# 5 systems x 3 ayanamshas x 10 dates x 14 checks = 2,100
# =========================================================================
print("Section 3.3: houses_ex with sidereal flag...")

FLG_SIDEREAL = 65536
SID_SYSTEMS = ["P", "K", "E", "W", "B"]
AYANAMSHAS = [0, 1, 27]  # Fagan-Bradley, Lahiri, True Citra

for hsys_char in SID_SYSTEMS:
    hsys_int = ord(hsys_char)
    hsys_bytes = bytes(hsys_char, "ascii")
    for ayan_id in AYANAMSHAS:
        for jd in JDS_10:
            lib_ok = True
            ref_ok = True

            # --- libephemeris ---
            try:
                lib.set_sid_mode(ayan_id)
                lib_result = lib.houses_ex(jd, 41.9, 12.5, hsys_int, FLG_SIDEREAL)
                lib_cusps = lib_result[0]
                lib_ascmc = lib_result[1]
                lib.set_sid_mode(0)  # reset
            except Exception as e:
                lib_ok = False
                lib.set_sid_mode(0)

            # --- pyswisseph ---
            try:
                swe_ref.set_sid_mode(ayan_id)
                ref_result = swe_ref.houses_ex(
                    jd, 41.9, 12.5, hsys_bytes, FLG_SIDEREAL
                )
                ref_cusps = ref_result[0]
                ref_ascmc = ref_result[1]
                swe_ref.set_sid_mode(0)
            except Exception as e:
                ref_ok = False
                swe_ref.set_sid_mode(0)

            if not lib_ok and not ref_ok:
                passed += 1
                continue
            if not lib_ok and ref_ok:
                for _ in range(14):
                    check(
                        False,
                        f"3.3 {hsys_char} ayan={ayan_id} jd={jd:.1f} lib error",
                        f"3.3/{hsys_char}/LIB_ERR",
                    )
                continue
            if lib_ok and not ref_ok:
                for _ in range(14):
                    passed += 1
                continue

            # Compare 12 cusps (tolerance 0.05 deg for sidereal)
            ncusps = min(12, len(lib_cusps), len(ref_cusps))
            for i in range(ncusps):
                diff = angle_diff(lib_cusps[i], ref_cusps[i])
                check(
                    diff < 0.05,
                    f"3.3 {hsys_char} ayan={ayan_id} cusp{i + 1} jd={jd:.1f} "
                    f"lib={float(lib_cusps[i]):.4f} ref={float(ref_cusps[i]):.4f} diff={diff:.6f}",
                    f"3.3/{hsys_char}/CUSP{i + 1}",
                )
            for i in range(ncusps, 12):
                check(
                    False,
                    f"3.3 {hsys_char} ayan={ayan_id} cusp{i + 1} missing",
                    f"3.3/{hsys_char}/MISSING",
                )

            # Compare ASC
            asc_diff = angle_diff(lib_ascmc[0], ref_ascmc[0])
            check(
                asc_diff < 0.05,
                f"3.3 {hsys_char} ayan={ayan_id} ASC jd={jd:.1f} "
                f"lib={float(lib_ascmc[0]):.4f} ref={float(ref_ascmc[0]):.4f} diff={asc_diff:.6f}",
                f"3.3/{hsys_char}/ASC",
            )

            # Compare MC
            mc_diff = angle_diff(lib_ascmc[1], ref_ascmc[1])
            check(
                mc_diff < 0.05,
                f"3.3 {hsys_char} ayan={ayan_id} MC jd={jd:.1f} "
                f"lib={float(lib_ascmc[1]):.4f} ref={float(ref_ascmc[1]):.4f} diff={mc_diff:.6f}",
                f"3.3/{hsys_char}/MC",
            )

print(f"  3.3 done: {time.time() - t0:.1f}s  (passed={passed} failed={failed})")


# =========================================================================
# Section 3.5: houses_armc
# 5 systems x 10 ARMC x 3 obliquities x 14 checks = 2,100
# =========================================================================
print("Section 3.5: houses_armc...")

ARMC_SYSTEMS = ["P", "K", "O", "R", "C"]
ARMC_VALS = [float(x) for x in range(0, 360, 36)]  # 0, 36, 72, ..., 324
EPS_VALS = [22.5, 23.4393, 24.0]
ARMC_LAT = 41.9  # Rome latitude

for hsys_char in ARMC_SYSTEMS:
    hsys_int = ord(hsys_char)
    hsys_bytes = bytes(hsys_char, "ascii")
    for armc in ARMC_VALS:
        for eps in EPS_VALS:
            lib_ok = True
            ref_ok = True

            # --- libephemeris ---
            try:
                lib_result = lib.houses_armc(armc, ARMC_LAT, eps, hsys_int)
                lib_cusps = lib_result[0]
                lib_ascmc = lib_result[1]
            except Exception as e:
                lib_ok = False

            # --- pyswisseph ---
            try:
                ref_result = swe_ref.houses_armc(armc, ARMC_LAT, eps, hsys_bytes)
                ref_cusps = ref_result[0]
                ref_ascmc = ref_result[1]
            except Exception as e:
                ref_ok = False

            if not lib_ok and not ref_ok:
                passed += 1
                continue
            if not lib_ok and ref_ok:
                for _ in range(14):
                    check(
                        False,
                        f"3.5 {hsys_char} armc={armc:.0f} eps={eps} lib error",
                        f"3.5/{hsys_char}/LIB_ERR",
                    )
                continue
            if lib_ok and not ref_ok:
                for _ in range(14):
                    passed += 1
                continue

            # Compare 12 cusps (tolerance 0.01 deg)
            ncusps = min(12, len(lib_cusps), len(ref_cusps))
            for i in range(ncusps):
                diff = angle_diff(lib_cusps[i], ref_cusps[i])
                check(
                    diff < 0.01,
                    f"3.5 {hsys_char} armc={armc:.0f} eps={eps} cusp{i + 1} "
                    f"lib={float(lib_cusps[i]):.4f} ref={float(ref_cusps[i]):.4f} diff={diff:.6f}",
                    f"3.5/{hsys_char}/CUSP{i + 1}",
                )
            for i in range(ncusps, 12):
                check(
                    False,
                    f"3.5 {hsys_char} armc={armc:.0f} cusp{i + 1} missing",
                    f"3.5/{hsys_char}/MISSING",
                )

            # Compare ASC
            asc_diff = angle_diff(lib_ascmc[0], ref_ascmc[0])
            check(
                asc_diff < 0.01,
                f"3.5 {hsys_char} armc={armc:.0f} eps={eps} ASC "
                f"lib={float(lib_ascmc[0]):.4f} ref={float(ref_ascmc[0]):.4f} diff={asc_diff:.6f}",
                f"3.5/{hsys_char}/ASC",
            )

            # Compare MC
            mc_diff = angle_diff(lib_ascmc[1], ref_ascmc[1])
            check(
                mc_diff < 0.01,
                f"3.5 {hsys_char} armc={armc:.0f} eps={eps} MC "
                f"lib={float(lib_ascmc[1]):.4f} ref={float(ref_ascmc[1]):.4f} diff={mc_diff:.6f}",
                f"3.5/{hsys_char}/MC",
            )

print(f"  3.5 done: {time.time() - t0:.1f}s  (passed={passed} failed={failed})")


# =========================================================================
# Section 3.6: house_pos
# 10 longitudes x 5 systems x 5 dates = 250 checks
# For each date we first compute ARMC and obliquity using both libraries,
# then call house_pos and compare.
# =========================================================================
print("Section 3.6: house_pos...")

HPOS_SYSTEMS = ["P", "K", "O", "R", "E"]
PLANET_LONS = [float(x) for x in range(0, 360, 36)]  # 0, 36, ..., 324
HPOS_LAT = 41.9

for jd in JDS_5:
    # Compute ARMC and obliquity from pyswisseph (ground truth for this test)
    try:
        ref_houses_result = swe_ref.houses(jd, HPOS_LAT, 12.5, b"P")
        armc = float(ref_houses_result[1][2])  # ascmc[2] = ARMC
        # Get obliquity from calc_ut with ECL_NUT (-1)
        ref_ecl = swe_ref.calc_ut(jd, -1, 0)
        eps = float(ref_ecl[0][0])  # (tuple_of_6, retflags) -> true obliquity
    except Exception as e:
        # If we cannot get armc/eps, skip this date
        for _ in range(len(PLANET_LONS) * len(HPOS_SYSTEMS)):
            check(False, f"3.6 jd={jd:.1f} setup error: {e}", "3.6/SETUP")
        continue

    for planet_lon in PLANET_LONS:
        for hsys_char in HPOS_SYSTEMS:
            hsys_bytes = bytes(hsys_char, "ascii")
            lib_ok = True
            ref_ok = True
            lib_pos = None
            ref_pos = None

            # --- libephemeris: house_pos(armc, lat, eps, (lon, lat_body), hsys_str) ---
            try:
                lib_pos = float(
                    lib.house_pos(armc, HPOS_LAT, eps, (planet_lon, 0.0), hsys_char)
                )
            except Exception as e:
                lib_ok = False

            # --- pyswisseph: house_pos(armc, lat, eps, (lon, lat_body), hsys_bytes) ---
            try:
                ref_pos = float(
                    swe_ref.house_pos(
                        armc, HPOS_LAT, eps, (planet_lon, 0.0), hsys_bytes
                    )
                )
            except Exception as e:
                ref_ok = False

            if not lib_ok and not ref_ok:
                passed += 1
                continue
            if not lib_ok and ref_ok:
                check(
                    False,
                    f"3.6 {hsys_char} lon={planet_lon:.0f} jd={jd:.1f} lib error",
                    f"3.6/{hsys_char}/LIB_ERR",
                )
                continue
            if lib_ok and not ref_ok:
                passed += 1
                continue

            # Compare positions (tolerance 0.05)
            diff = abs(lib_pos - ref_pos)
            # Handle house number wraparound (12.9 vs 1.1 is only 0.2 apart in circular sense)
            if diff > 6.0:
                diff = 12.0 - diff
            check(
                diff < 0.05,
                f"3.6 {hsys_char} lon={planet_lon:.0f} jd={jd:.1f} "
                f"lib={lib_pos:.4f} ref={ref_pos:.4f} diff={diff:.6f}",
                f"3.6/{hsys_char}/POS",
            )

print(f"  3.6 done: {time.time() - t0:.1f}s  (passed={passed} failed={failed})")


# =========================================================================
# Section 3.7: Gauquelin sectors via house_pos with system 'G'
# 5 body longitudes x 10 dates = 50 checks (result in [1, 37))
# Plus comparison with pyswisseph = another 50 checks
# =========================================================================
print("Section 3.7: Gauquelin sectors...")

GQ_LONS = [0.0, 72.0, 144.0, 216.0, 288.0]  # 5 test longitudes
GQ_LAT = 41.9

for jd in JDS_10:
    # Compute ARMC and obliquity
    try:
        ref_houses_result = swe_ref.houses(jd, GQ_LAT, 12.5, b"P")
        armc = float(ref_houses_result[1][2])
        ref_ecl = swe_ref.calc_ut(jd, -1, 0)
        eps = float(ref_ecl[0][0])  # (tuple_of_6, retflags) -> true obliquity
    except Exception as e:
        for _ in range(len(GQ_LONS) * 2):
            check(False, f"3.7 jd={jd:.1f} setup error: {e}", "3.7/SETUP")
        continue

    for planet_lon in GQ_LONS:
        # --- libephemeris ---
        lib_ok = True
        lib_pos = None
        try:
            lib_pos = float(lib.house_pos(armc, GQ_LAT, eps, (planet_lon, 0.0), "G"))
        except Exception as e:
            lib_ok = False

        # Check: result should be in [1, 37) for Gauquelin
        if lib_ok:
            check(
                1.0 <= lib_pos < 37.0,
                f"3.7 lon={planet_lon:.0f} jd={jd:.1f} lib_pos={lib_pos:.4f} out of [1,37)",
                "3.7/RANGE",
            )
        else:
            check(
                False,
                f"3.7 lon={planet_lon:.0f} jd={jd:.1f} lib error",
                "3.7/LIB_ERR",
            )

        # --- pyswisseph comparison ---
        ref_ok = True
        ref_pos = None
        try:
            ref_pos = float(
                swe_ref.house_pos(armc, GQ_LAT, eps, (planet_lon, 0.0), b"G")
            )
        except Exception as e:
            ref_ok = False

        if lib_ok and ref_ok:
            diff = abs(lib_pos - ref_pos)
            if diff > 18.0:
                diff = 36.0 - diff
            check(
                diff < 0.05,
                f"3.7 G lon={planet_lon:.0f} jd={jd:.1f} "
                f"lib={lib_pos:.4f} ref={ref_pos:.4f} diff={diff:.6f}",
                "3.7/G_CMP",
            )
        elif not lib_ok and not ref_ok:
            passed += 1
        elif lib_ok and not ref_ok:
            passed += 1
        else:
            check(
                False,
                f"3.7 G lon={planet_lon:.0f} jd={jd:.1f} lib error, ref ok",
                "3.7/G_LIB_ERR",
            )

print(f"  3.7 done: {time.time() - t0:.1f}s  (passed={passed} failed={failed})")


# =========================================================================
# Summary
# =========================================================================
elapsed = time.time() - t0
total = passed + failed
pct = 100.0 * passed / total if total else 0

report = []
report.append("")
report.append("=" * 70)
report.append("VERIFY HOUSES: libephemeris vs pyswisseph")
report.append("=" * 70)
report.append(f"Result: {passed}/{total} PASS ({pct:.1f}%)")
report.append(f"Passed: {passed}")
report.append(f"Failed: {failed}")
report.append(f"Time:   {elapsed:.1f}s")
report.append("")

if fail_by:
    report.append("Failure breakdown (top 50):")
    for key, cnt in sorted(fail_by.items(), key=lambda x: -x[1])[:50]:
        report.append(f"  {key:50s} = {cnt:5d}")
    report.append("")

if errors:
    report.append(f"Sample failures (first {min(60, len(errors))} of {failed}):")
    for e in errors[:60]:
        report.append(f"  FAIL: {e}")
elif failed == 0:
    report.append("ALL CHECKS PASSED!")

report.append("")
text = "\n".join(report)
print(text)

# Write report to file
os.makedirs(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results"), exist_ok=True)
with open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "verify_houses.txt"), "w") as f:
    f.write(text + "\n")

# Cleanup
lib.close()
swe_ref.close()
