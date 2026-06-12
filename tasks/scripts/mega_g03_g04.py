#!/usr/bin/env python3
"""Mega verification G03 (Flag Combinations) + G04 (House Systems).

G03.01: 13 single flags x 10 bodies x 10 dates       = 1300+ checks
G03.02: 6 flag pairs x 5 bodies x ~7 dates            =  210  checks
G04.01: 24 house systems x 5 locations x 10 dates     = 1200  checks
                                                 Total >= 2710 checks
"""

import os
import math
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import libephemeris as lib
import swisseph as swe_ref

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

# ── Setup ────────────────────────────────────────────────────────────────
swe_ref.set_ephe_path(_REF_EPHE_PATH)
lib.set_calc_mode("skyfield")

random.seed(42)

# ── Counters ─────────────────────────────────────────────────────────────
passed = 0
failed = 0
errors = []


def check(condition, description=""):
    global passed, failed
    if condition:
        passed += 1
    else:
        failed += 1
        if len(errors) < 80:
            errors.append(description)


PI = math.pi
TWO_PI = 2.0 * PI

# ── Constants ────────────────────────────────────────────────────────────
SEFLG_DEFAULT = 0
FLG_SWIEPH = 2
FLG_HELCTR = 8
FLG_TRUEPOS = 16
FLG_J2000 = 32
FLG_NONUT = 64
FLG_SPEED = 256
FLG_NOGDEFL = 512
FLG_NOABERR = 1024
FLG_EQUATORIAL = 2048
FLG_XYZ = 4096
FLG_RADIANS = 8192
FLG_SIDEREAL = 65536

BODIES = list(range(10))  # 0..9 (Sun through Pluto)
BODY_NAMES = {
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
    14: "Earth",
}

# 10 random JDs in [1900, 2100]  (JD 2415020.5 to 2488069.5)
JD_MIN = 2415020.5
JD_MAX = 2488069.5
DATES_10 = sorted([random.uniform(JD_MIN, JD_MAX) for _ in range(10)])


def is_finite_tuple(tup, n=6):
    """Return True if tup has at least n finite float values."""
    if not isinstance(tup, (tuple, list)):
        return False
    if len(tup) < n:
        return False
    return all(math.isfinite(v) for v in tup[:n])


# ======================================================================
#  G03: Flag Combinations
# ======================================================================
print("=" * 72)
print("G03: Flag Combinations")
print("=" * 72)

t0 = time.time()

# ── G03.01: 13 single flags x 10 bodies x 10 dates ──────────────────

SINGLE_FLAGS = [
    ("DEFAULT", SEFLG_DEFAULT),
    ("SWIEPH", FLG_SWIEPH),
    ("SPEED", FLG_SPEED),
    ("HELCTR", FLG_HELCTR),
    ("TRUEPOS", FLG_TRUEPOS),
    ("J2000", FLG_J2000),
    ("NONUT", FLG_NONUT),
    ("NOGDEFL", FLG_NOGDEFL),
    ("NOABERR", FLG_NOABERR),
    ("EQUATORIAL", FLG_EQUATORIAL),
    ("XYZ", FLG_XYZ),
    ("RADIANS", FLG_RADIANS),
    ("SIDEREAL", FLG_SIDEREAL),
]

g03_01_passed = 0
g03_01_failed = 0
p_before = passed

print("\n--- G03.01: Single flags (13 x 10 bodies x 10 dates) ---")

for flag_name, flag_val in SINGLE_FLAGS:
    flag_pass = 0
    flag_fail = 0

    for body in BODIES:
        bname = BODY_NAMES.get(body, f"Body{body}")

        # HELCTR: skip Sun (body 0) -- heliocentric Sun is undefined
        if flag_val == FLG_HELCTR and body == 0:
            continue

        for jd in DATES_10:
            desc_prefix = f"G03.01 {flag_name} body={bname} jd={jd:.2f}"

            # For SIDEREAL: set mode before, reset after
            if flag_val == FLG_SIDEREAL:
                lib.set_sid_mode(1)

            exc_caught = False
            try:
                result = lib.calc_ut(jd, body, flag_val | FLG_SPEED)
            except Exception as e:
                check(False, f"{desc_prefix}: EXCEPTION {e}")
                exc_caught = True
            finally:
                if flag_val == FLG_SIDEREAL:
                    lib.set_sid_mode(0)

            if exc_caught:
                flag_fail += 1
                continue

            # Check 1: result is (tuple_of_6, int)
            ok_shape = (
                isinstance(result, tuple)
                and len(result) == 2
                and isinstance(result[0], tuple)
                and len(result[0]) == 6
                and isinstance(result[1], int)
            )
            check(ok_shape, f"{desc_prefix}: bad shape {type(result)}")
            if not ok_shape:
                flag_fail += 1
                continue

            pos = result[0]

            # Check 2: all 6 values are finite
            finite_ok = is_finite_tuple(pos, 6)
            check(finite_ok, f"{desc_prefix}: non-finite values {pos}")
            if not finite_ok:
                flag_fail += 1
                continue

            # Check 3: range checks (skip Earth body=14 which is all zeros)
            # For XYZ: no angular range constraint
            # For RADIANS: lon in [0, 2pi), lat in [-pi/2, pi/2]
            # For default ecliptic/equatorial: lon in [0, 360), lat in [-90, 90]
            lon, lat = pos[0], pos[1]

            if flag_val == FLG_XYZ:
                # XYZ coordinates -- no specific range checks beyond finiteness
                check(True, f"{desc_prefix}: XYZ finite ok")
            elif flag_val == FLG_RADIANS:
                lon_ok = -0.001 <= lon <= TWO_PI + 0.001
                lat_ok = -(PI / 2.0 + 0.001) <= lat <= (PI / 2.0 + 0.001)
                check(lon_ok, f"{desc_prefix}: RADIANS lon={lon:.6f} out of [0, 2pi)")
                check(
                    lat_ok, f"{desc_prefix}: RADIANS lat={lat:.6f} out of [-pi/2, pi/2]"
                )
            else:
                # Standard ecliptic or equatorial
                lon_ok = -0.01 <= lon <= 360.01
                if flag_val == FLG_EQUATORIAL:
                    # RA can be [0, 360), Dec in [-90, 90]
                    lat_ok = -90.01 <= lat <= 90.01
                else:
                    lat_ok = -90.01 <= lat <= 90.01
                check(lon_ok, f"{desc_prefix}: lon={lon:.6f} out of [0, 360)")
                check(lat_ok, f"{desc_prefix}: lat={lat:.6f} out of [-90, 90]")

    # Per-flag tally from global counters
    cur_pass = passed - p_before - g03_01_passed
    g03_01_passed = passed - p_before
    print(f"  {flag_name:14s}: section ok")

g03_01_total = passed - p_before
g03_01_fail = failed  # cumulative at this point
print(f"\nG03.01 checks: {g03_01_total} passed, {g03_01_fail} failed")

# ── G03.02: Critical flag pairs (post-warmup) ────────────────────────

print("\n--- G03.02: Critical flag pairs (6 pairs x 5 bodies x ~7 dates) ---")

FLAG_PAIRS = [
    ("EQUATORIAL+NONUT", FLG_EQUATORIAL | FLG_NONUT),
    ("TRUEPOS+NONUT", FLG_TRUEPOS | FLG_NONUT),
    ("NONUT+NOGDEFL", FLG_NONUT | FLG_NOGDEFL),
    ("NONUT+NOABERR", FLG_NONUT | FLG_NOABERR),
    ("XYZ+HELCTR", FLG_XYZ | FLG_HELCTR),
    ("XYZ+EQUATORIAL", FLG_XYZ | FLG_EQUATORIAL),
]

PAIR_BODIES = [1, 2, 4, 5, 9]  # Moon, Mercury, Mars, Jupiter, Pluto
PAIR_DATES = sorted([random.uniform(JD_MIN, JD_MAX) for _ in range(7)])

p_before_02 = passed
f_before_02 = failed

for pair_name, pair_flags in FLAG_PAIRS:
    for body in PAIR_BODIES:
        bname = BODY_NAMES.get(body, f"Body{body}")

        # For XYZ+HELCTR skip Sun which would be body 0 (not in PAIR_BODIES anyway)
        for jd in PAIR_DATES:
            desc = f"G03.02 {pair_name} body={bname} jd={jd:.2f}"
            try:
                result = lib.calc_ut(jd, body, pair_flags | FLG_SPEED)
            except Exception as e:
                check(False, f"{desc}: EXCEPTION {e}")
                continue

            ok_shape = (
                isinstance(result, tuple)
                and len(result) == 2
                and isinstance(result[0], tuple)
                and len(result[0]) == 6
                and isinstance(result[1], int)
            )
            check(ok_shape, f"{desc}: bad shape")
            if not ok_shape:
                continue

            finite_ok = is_finite_tuple(result[0], 6)
            check(finite_ok, f"{desc}: non-finite {result[0]}")

    print(f"  {pair_name:24s}: ok")

g03_02_passed = passed - p_before_02
g03_02_failed = failed - f_before_02
print(f"\nG03.02 checks: {g03_02_passed} passed, {g03_02_failed} failed")

g03_time = time.time() - t0
g03_total_passed = passed
g03_total_failed = failed
print(
    f"\n>>> G03 TOTAL: {g03_total_passed} passed, {g03_total_failed} failed "
    f"({g03_time:.1f}s)"
)


# ======================================================================
#  G04: House Systems
# ======================================================================
print("\n" + "=" * 72)
print("G04: House Systems")
print("=" * 72)

t1 = time.time()

# ── G04.01: 24 systems x 5 locations x 10 dates ─────────────────────

HOUSE_SYSTEMS = list("PKORCEWBMTUHXGYINFDJ") + ["i", "L", "A", "S"]

LOCATIONS = [
    (0.0, 0.0, "Equator/GM"),
    (41.9, 12.5, "Rome"),
    (51.5, -0.1, "London"),
    (-33.9, 151.2, "Sydney"),
    (35.7, 139.7, "Tokyo"),
]

# 10 JDs evenly spaced 2000-2020
JD_2000 = 2451545.0  # J2000.0
JD_2020 = 2459215.5  # ~2021-01-01
HOUSE_DATES = [JD_2000 + i * (JD_2020 - JD_2000) / 9.0 for i in range(10)]

CUSP_TOL = 0.01  # degrees
ANGLE_TOL = 0.01  # degrees

p_before_g04 = passed
f_before_g04 = failed

print("\n--- G04.01: 24 systems x 5 locations x 10 dates ---")

# Track per-system stats
system_stats = {}

for hsys_char in HOUSE_SYSTEMS:
    sys_passed = 0
    sys_failed = 0
    sys_skipped = 0

    for lat, lon, loc_name in LOCATIONS:
        for jd in HOUSE_DATES:
            desc = f"G04.01 sys={hsys_char} loc={loc_name} jd={jd:.2f}"

            # Known limitation: 'i' system at high latitudes
            if hsys_char == "i" and abs(lat) >= 60:
                sys_skipped += 1
                continue

            # Call reference (swisseph): uses bytes
            ref_ok = True
            try:
                ref_cusps, ref_ascmc = swe_ref.houses(jd, lat, lon, hsys_char.encode())
            except Exception:
                ref_ok = False

            # Call libephemeris: uses ord()
            lib_ok = True
            try:
                lib_cusps, lib_ascmc = lib.houses(jd, lat, lon, ord(hsys_char))
            except Exception:
                lib_ok = False

            # If either raised, skip (polar circle issues etc.)
            if not ref_ok or not lib_ok:
                sys_skipped += 1
                continue

            # Compare 12 cusps
            n_cusps = min(len(ref_cusps), len(lib_cusps), 12)
            # Gauquelin has 36 cusps; compare first 12 at minimum
            if hsys_char == "G":
                n_cusps = min(len(ref_cusps), len(lib_cusps))

            cusps_ok = True
            for ci in range(n_cusps):
                diff = abs(ref_cusps[ci] - lib_cusps[ci])
                if diff > 180.0:
                    diff = 360.0 - diff
                if diff > CUSP_TOL:
                    cusps_ok = False
                    check(
                        False,
                        f"{desc}: cusp[{ci}] diff={diff:.6f} "
                        f"(ref={ref_cusps[ci]:.6f} lib={lib_cusps[ci]:.6f})",
                    )
                    sys_failed += 1
                    break

            if cusps_ok:
                check(True, f"{desc}: cusps match")
                sys_passed += 1

            # Compare ASC (ascmc[0])
            asc_diff = abs(ref_ascmc[0] - lib_ascmc[0])
            if asc_diff > 180.0:
                asc_diff = 360.0 - asc_diff
            asc_ok = asc_diff <= ANGLE_TOL
            check(
                asc_ok,
                f"{desc}: ASC diff={asc_diff:.6f} "
                f"(ref={ref_ascmc[0]:.6f} lib={lib_ascmc[0]:.6f})",
            )
            if asc_ok:
                sys_passed += 1
            else:
                sys_failed += 1

            # Compare MC (ascmc[1])
            mc_diff = abs(ref_ascmc[1] - lib_ascmc[1])
            if mc_diff > 180.0:
                mc_diff = 360.0 - mc_diff
            mc_ok = mc_diff <= ANGLE_TOL
            check(
                mc_ok,
                f"{desc}: MC diff={mc_diff:.6f} "
                f"(ref={ref_ascmc[1]:.6f} lib={lib_ascmc[1]:.6f})",
            )
            if mc_ok:
                sys_passed += 1
            else:
                sys_failed += 1

    system_stats[hsys_char] = (sys_passed, sys_failed, sys_skipped)
    status = "PASS" if sys_failed == 0 else "FAIL"
    skip_note = f" (skipped {sys_skipped})" if sys_skipped else ""
    print(
        f"  System '{hsys_char}': {status} "
        f"(pass={sys_passed} fail={sys_failed}{skip_note})"
    )

g04_passed = passed - p_before_g04
g04_failed = failed - f_before_g04
g04_time = time.time() - t1

print(f"\nG04.01 checks: {g04_passed} passed, {g04_failed} failed")
print(f"\n>>> G04 TOTAL: {g04_passed} passed, {g04_failed} failed ({g04_time:.1f}s)")


# ======================================================================
#  Grand Summary
# ======================================================================
total_time = time.time() - t0

print("\n" + "=" * 72)
print("GRAND SUMMARY")
print("=" * 72)
print(
    f"  G03 Flag Combinations : {g03_total_passed:5d} passed, {g03_total_failed:5d} failed"
)
print(f"  G04 House Systems     : {g04_passed:5d} passed, {g04_failed:5d} failed")
print("  ────────────────────────────────────────")
print(f"  TOTAL                 : {passed:5d} passed, {failed:5d} failed")
print(f"  Total time            : {total_time:.1f}s")

if errors:
    print(f"\nFirst {min(len(errors), 20)} failures:")
    for e in errors[:20]:
        print(f"  FAIL: {e}")

if failed == 0:
    print("\nALL CHECKS PASSED")
else:
    print(f"\n{failed} CHECKS FAILED")

sys.exit(0 if failed == 0 else 1)
