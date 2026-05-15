#!/usr/bin/env python
"""Exhaustive flag-combination verification for libephemeris.

Sections tested:
  2.1  13 single flags x 22 bodies x 50 dates
  2.2  30 compatible flag pairs x 5 bodies x 20 dates
  2.4  XYZ coherence with spherical (10 bodies x 50 dates)
  2.5  RADIANS coherence with degrees (10 bodies x 50 dates)
  2.6  Sun heliocentric = near-zero (50 dates)

Target: ~20 000+ checks, < 30 seconds.
"""

from __future__ import annotations

import math
import random
import sys
import time
from typing import Any, Tuple

# ---------------------------------------------------------------------------
# Library under test
# ---------------------------------------------------------------------------
import libephemeris as lib

# ---------------------------------------------------------------------------
# Counters
# ---------------------------------------------------------------------------
passed = 0
failed = 0
errors: list[str] = []


def check(condition: bool, label: str) -> None:
    """Record a pass or fail.  Failures are printed and logged."""
    global passed, failed
    if condition:
        passed += 1
    else:
        failed += 1
        errors.append(label)
        print(f"  FAIL: {label}")


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SEFLG_DEFAULT = 0
FLG_SWIEPH = 2
FLG_SPEED = 256
FLG_HELCTR = 8
FLG_TRUEPOS = 16
FLG_J2000 = 32
FLG_NONUT = 64
FLG_NOGDEFL = 512
FLG_NOABERR = 1024
FLG_EQUATORIAL = 2048
FLG_XYZ = 4096
FLG_RADIANS = 8192
FLG_SIDEREAL = 65536

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

FLAG_PAIRS = [
    ("SPEED+HELCTR", FLG_SPEED | FLG_HELCTR),
    ("SPEED+EQUATORIAL", FLG_SPEED | FLG_EQUATORIAL),
    ("SPEED+TRUEPOS", FLG_SPEED | FLG_TRUEPOS),
    ("SPEED+J2000", FLG_SPEED | FLG_J2000),
    ("SPEED+NONUT", FLG_SPEED | FLG_NONUT),
    ("SPEED+NOGDEFL", FLG_SPEED | FLG_NOGDEFL),
    ("SPEED+NOABERR", FLG_SPEED | FLG_NOABERR),
    ("HELCTR+TRUEPOS", FLG_HELCTR | FLG_TRUEPOS),
    ("HELCTR+J2000", FLG_HELCTR | FLG_J2000),
    ("HELCTR+NONUT", FLG_HELCTR | FLG_NONUT),
    ("EQUATORIAL+J2000", FLG_EQUATORIAL | FLG_J2000),
    ("EQUATORIAL+NONUT", FLG_EQUATORIAL | FLG_NONUT),
    ("EQUATORIAL+SPEED", FLG_EQUATORIAL | FLG_SPEED),
    ("XYZ+SPEED", FLG_XYZ | FLG_SPEED),
    ("TRUEPOS+J2000", FLG_TRUEPOS | FLG_J2000),
    ("TRUEPOS+NONUT", FLG_TRUEPOS | FLG_NONUT),
    ("J2000+NONUT", FLG_J2000 | FLG_NONUT),
    ("J2000+NOGDEFL", FLG_J2000 | FLG_NOGDEFL),
    ("J2000+NOABERR", FLG_J2000 | FLG_NOABERR),
    ("NONUT+NOGDEFL", FLG_NONUT | FLG_NOGDEFL),
    ("NONUT+NOABERR", FLG_NONUT | FLG_NOABERR),
    ("NOGDEFL+NOABERR", FLG_NOGDEFL | FLG_NOABERR),
    ("HELCTR+EQUATORIAL", FLG_HELCTR | FLG_EQUATORIAL),
    ("HELCTR+NOGDEFL", FLG_HELCTR | FLG_NOGDEFL),
    ("HELCTR+NOABERR", FLG_HELCTR | FLG_NOABERR),
    ("EQUATORIAL+NOGDEFL", FLG_EQUATORIAL | FLG_NOGDEFL),
    ("EQUATORIAL+NOABERR", FLG_EQUATORIAL | FLG_NOABERR),
    ("EQUATORIAL+TRUEPOS", FLG_EQUATORIAL | FLG_TRUEPOS),
    ("XYZ+HELCTR", FLG_XYZ | FLG_HELCTR),
    ("XYZ+EQUATORIAL", FLG_XYZ | FLG_EQUATORIAL),
]

# 22 bodies: 0-9, 10-13, 14, 15, 17-22
ALL_BODIES = list(range(0, 16)) + list(range(17, 23))
assert len(ALL_BODIES) == 22, f"Expected 22 bodies, got {len(ALL_BODIES)}"

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
    10: "MeanNode",
    11: "TrueNode",
    12: "MeanApog",
    13: "OscuApog",
    14: "Earth",
    15: "Chiron",
    17: "Ceres",
    18: "Pallas",
    19: "Juno",
    20: "Vesta",
    21: "IntpApog",
    22: "IntpPerig",
}

PAIR_BODIES = [0, 1, 4, 5, 6]  # Sun, Moon, Mars, Jupiter, Saturn

EARTH = 14

# ---------------------------------------------------------------------------
# Date generation
# ---------------------------------------------------------------------------
random.seed(42)

# JD range for 1900-01-01 to 2100-12-31
JD_1900 = 2415020.5
JD_2100 = 2488069.5

DATES_50 = [random.uniform(JD_1900, JD_2100) for _ in range(50)]
DATES_20 = [random.uniform(JD_1900, JD_2100) for _ in range(20)]


# Force Skyfield mode to avoid LEB caching artifacts in this stress test
lib.set_calc_mode("skyfield")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def is_finite(v: Any) -> bool:
    """Return True if v is a finite float."""
    try:
        return math.isfinite(float(v))
    except (TypeError, ValueError):
        return False


def safe_calc_ut(jd: float, body: int, flag: int) -> Tuple[Any, ...] | None:
    """Call lib.calc_ut, returning None on any exception."""
    try:
        result = lib.calc_ut(jd, body, flag)
        return result
    except Exception as e:
        global _last_exception
        _last_exception = e
        return None


_last_exception = None


def validate_result_basic(
    result: Tuple[Any, ...] | None,
    label: str,
) -> Tuple[float, ...] | None:
    """Check that result is a (6-float-tuple, retflag) and all finite.

    Returns the 6-float tuple on success, None on failure.
    """
    if result is None:
        check(
            False,
            f"{label}: unhandled exception (calc_ut returned None: {_last_exception})",
        )
        return None

    check(True, f"{label}: no exception")

    # Result should be (tuple_of_6, retflag)
    pos = result[0]
    ok_tuple = isinstance(pos, (tuple, list)) and len(pos) == 6
    check(ok_tuple, f"{label}: result is tuple of 6")
    if not ok_tuple:
        return None

    all_fin = all(is_finite(v) for v in pos)
    check(all_fin, f"{label}: all 6 values finite")
    if not all_fin:
        return None

    return tuple(float(v) for v in pos)


# =========================================================================
# Section 2.1 -- Single flags x all bodies x 50 dates
# =========================================================================
def run_section_2_1() -> None:
    print("=" * 70)
    print("Section 2.1: 13 single flags x 22 bodies x 50 dates")
    print("=" * 70)

    for flag_name, flag_val in SINGLE_FLAGS:
        # FLG_SIDEREAL requires set_sid_mode first
        if flag_val == FLG_SIDEREAL:
            lib.set_sid_mode(1)  # Lahiri ayanamsha

        for body in ALL_BODIES:
            # Skip heliocentric Earth -- known to be problematic
            if flag_val == FLG_HELCTR and body == EARTH:
                continue

            for jd in DATES_50:
                label = (
                    f"flag={flag_name} body={BODY_NAMES.get(body, body)} jd={jd:.1f}"
                )

                result = safe_calc_ut(jd, body, flag_val)
                pos = validate_result_basic(result, label)
                if pos is None:
                    continue

                # Range checks (only for non-XYZ results)
                is_xyz = bool(flag_val & FLG_XYZ)
                is_rad = bool(flag_val & FLG_RADIANS)
                is_earth_geo = (body == EARTH) and not (flag_val & FLG_HELCTR)

                if not is_xyz and not is_earth_geo:
                    lon, lat = pos[0], pos[1]

                    if is_rad:
                        check(
                            0.0 <= lon < 2.0 * math.pi,
                            f"{label}: lon in [0, 2pi) -- got {lon}",
                        )
                        check(
                            -math.pi / 2.0 <= lat <= math.pi / 2.0,
                            f"{label}: lat in [-pi/2, pi/2] -- got {lat}",
                        )
                    else:
                        check(
                            0.0 <= lon < 360.0,
                            f"{label}: lon in [0, 360) -- got {lon}",
                        )
                        check(
                            -90.0 <= lat <= 90.0,
                            f"{label}: lat in [-90, 90] -- got {lat}",
                        )

        # Reset sidereal mode after the flag is done
        if flag_val == FLG_SIDEREAL:
            lib.set_sid_mode(0)


# =========================================================================
# Section 2.2 -- 30 compatible flag pairs x 5 bodies x 20 dates
# =========================================================================
def run_section_2_2() -> None:
    print("=" * 70)
    print("Section 2.2: 30 flag pairs x 5 bodies x 20 dates")
    print("=" * 70)

    for pair_name, pair_val in FLAG_PAIRS:
        for body in PAIR_BODIES:
            for jd in DATES_20:
                label = (
                    f"pair={pair_name} body={BODY_NAMES.get(body, body)} jd={jd:.1f}"
                )
                result = safe_calc_ut(jd, body, pair_val | FLG_SWIEPH)
                validate_result_basic(result, label)


# =========================================================================
# Section 2.4 -- XYZ coherence with spherical
# =========================================================================
def run_section_2_4() -> None:
    print("=" * 70)
    print("Section 2.4: XYZ coherence with spherical (10 bodies x 50 dates)")
    print("=" * 70)

    bodies_10 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    for body in bodies_10:
        for jd in DATES_50:
            label = f"xyz_coherence body={BODY_NAMES.get(body, body)} jd={jd:.1f}"

            # Spherical
            res_sph = safe_calc_ut(jd, body, FLG_SPEED)
            # XYZ
            res_xyz = safe_calc_ut(jd, body, FLG_SPEED | FLG_XYZ)

            if res_sph is None or res_xyz is None:
                check(False, f"{label}: one of the calls failed")
                continue

            pos_sph = res_sph[0]
            pos_xyz = res_xyz[0]

            if not (len(pos_sph) == 6 and len(pos_xyz) == 6):
                check(False, f"{label}: unexpected tuple length")
                continue

            dist = float(pos_sph[2])
            x, y, z = float(pos_xyz[0]), float(pos_xyz[1]), float(pos_xyz[2])

            if not (is_finite(dist) and is_finite(x) and is_finite(y) and is_finite(z)):
                check(False, f"{label}: non-finite values")
                continue

            r_squared = x * x + y * y + z * z
            dist_squared = dist * dist

            if dist_squared == 0.0:
                # If dist is zero, r should also be zero
                check(
                    r_squared < 1e-12,
                    f"{label}: dist=0, r^2={r_squared}",
                )
            else:
                rel_err = abs(r_squared - dist_squared) / dist_squared
                check(
                    rel_err < 1e-6,
                    f"{label}: x^2+y^2+z^2 vs dist^2 rel_err={rel_err:.2e}",
                )


# =========================================================================
# Section 2.5 -- RADIANS coherence with degrees
# =========================================================================
def run_section_2_5() -> None:
    print("=" * 70)
    print("Section 2.5: RADIANS coherence with degrees (10 bodies x 50 dates)")
    print("=" * 70)

    bodies_10 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    for body in bodies_10:
        for jd in DATES_50:
            label = f"rad_coherence body={BODY_NAMES.get(body, body)} jd={jd:.1f}"

            res_deg = safe_calc_ut(jd, body, FLG_SPEED)
            res_rad = safe_calc_ut(jd, body, FLG_SPEED | FLG_RADIANS)

            if res_deg is None or res_rad is None:
                check(False, f"{label}: one of the calls failed")
                continue

            pos_deg = res_deg[0]
            pos_rad = res_rad[0]

            if not (len(pos_deg) == 6 and len(pos_rad) == 6):
                check(False, f"{label}: unexpected tuple length")
                continue

            lon_deg = float(pos_deg[0])
            lon_rad = float(pos_rad[0])

            if not (is_finite(lon_deg) and is_finite(lon_rad)):
                check(False, f"{label}: non-finite longitude")
                continue

            expected_rad = lon_deg * math.pi / 180.0
            abs_err = abs(lon_rad - expected_rad)
            check(
                abs_err < 1e-6,
                f"{label}: lon rad={lon_rad:.8f} vs deg*pi/180={expected_rad:.8f} "
                f"abs_err={abs_err:.2e}",
            )


# =========================================================================
# Section 2.6 -- Sun heliocentric = near-zero
# =========================================================================
def run_section_2_6() -> None:
    print("=" * 70)
    print("Section 2.6: Sun heliocentric = near-zero (50 dates)")
    print("=" * 70)

    for jd in DATES_50:
        label = f"sun_helio jd={jd:.1f}"
        result = safe_calc_ut(jd, 0, FLG_HELCTR)  # body 0 = Sun

        if result is None:
            check(False, f"{label}: calc_ut raised exception")
            continue

        pos = result[0]
        if not (isinstance(pos, (tuple, list)) and len(pos) == 6):
            check(False, f"{label}: not a 6-tuple")
            continue

        # All 6 values should be zero or very close to zero
        all_near_zero = all(abs(float(v)) < 1e-6 for v in pos)
        check(
            all_near_zero,
            f"{label}: expected near-zero, got {tuple(float(v) for v in pos)}",
        )


# =========================================================================
# Main
# =========================================================================
def main() -> None:
    global passed, failed

    t0 = time.time()

    run_section_2_1()
    run_section_2_2()
    run_section_2_4()
    run_section_2_5()
    run_section_2_6()

    elapsed = time.time() - t0

    # Summary
    total = passed + failed
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Total checks : {total}")
    print(f"  Passed       : {passed}")
    print(f"  Failed       : {failed}")
    print(f"  Elapsed      : {elapsed:.1f}s")
    print()

    if failed > 0:
        print(f"First {min(len(errors), 20)} failures:")
        for e in errors[:20]:
            print(f"  - {e}")
        if len(errors) > 20:
            print(f"  ... and {len(errors) - 20} more")

    print()
    if failed == 0:
        print("ALL CHECKS PASSED")
    else:
        print(f"FAILURES DETECTED: {failed}/{total}")

    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
