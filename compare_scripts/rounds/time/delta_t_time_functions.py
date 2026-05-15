#!/usr/bin/env python3
"""
Round 11: Deep Delta-T and Time Functions Audit
=================================================
Compares libephemeris time-related functions against pyswisseph:
  P1: deltat — Delta-T at many epochs (historical, current, future)
  P2: deltat_ex — Extended Delta-T with ephemeris flag
  P3: julday / revjul — Julian Day conversions
  P4: utc_to_jd / jdut1_to_utc — UTC conversions
  P5: date_conversion — Date validation and JD
  P6: day_of_week — Day of week
  P7: Delta-T sweep 1600-2100 — systematic drift detection
  P8: TT-UT consistency (calc_ut vs calc with manual Delta-T)
"""

from __future__ import annotations

import os
import sys
import time
import traceback

import swisseph as swe

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import libephemeris as ephem

EPHE_PATH = os.path.join(os.path.dirname(__file__), "..", "swisseph", "ephe")
if os.path.exists(EPHE_PATH):
    swe.set_ephe_path(EPHE_PATH)

# ============================================================================
# COUNTERS
# ============================================================================

total = 0
passed = 0
failed = 0
skipped = 0
failures = []


def record(name, ok, detail=""):
    global total, passed, failed, failures
    total += 1
    if ok:
        passed += 1
        print(f"  [PASS] {name}: {detail}")
    else:
        failed += 1
        failures.append((name, detail))
        print(f"  [FAIL] {name}: {detail}")


# ============================================================================
# PART 1: deltat at many epochs
# ============================================================================


def test_part1_deltat():
    print("\n" + "=" * 70)
    print("PART 1: deltat — Delta-T at Multiple Epochs")
    print("=" * 70)

    # Test epochs spanning historical to future
    epochs = [
        ("1600", 2305447.5),
        ("1700", 2341972.5),
        ("1800", 2378496.5),
        ("1900", 2415020.5),
        ("1950", 2433282.5),
        ("1970", 2440587.5),
        ("1980", 2444239.5),
        ("1990", 2447892.5),
        ("2000-Jan", 2451544.5),
        ("J2000.0", 2451545.0),
        ("2005", 2453371.5),
        ("2010", 2455197.5),
        ("2015", 2457023.5),
        ("2020", 2458849.5),
        ("2024-Jan", 2460310.5),
        ("2030", 2462502.5),
        ("2050", 2469807.5),
        ("2100", 2488069.5),
    ]

    for epoch_name, jd in epochs:
        test_name = f"P1/deltat/{epoch_name}"
        try:
            dt_se = swe.deltat(jd)
            dt_le = ephem.deltat(jd)

            diff = abs(dt_se - dt_le)
            diff_sec = diff * 86400.0  # Convert days to seconds

            # Tolerance: historical epochs may differ more due to model differences
            # Modern (1900+): within 0.5s, Historical: within 25s, Future: within 10s
            # The 1600 epoch differs by ~20s due to different Delta-T polynomial
            # models (Stephenson/Morrison/Hohenkerk 2016 vs Espenak/Meeus 2006).
            if jd < 2415020.5:  # Before 1900
                tol_sec = 25.0
            elif jd > 2462502.5:  # After 2030
                tol_sec = 10.0
            else:
                tol_sec = 0.5

            ok = diff_sec < tol_sec
            record(
                test_name,
                ok,
                f"SE={dt_se * 86400:.3f}s LE={dt_le * 86400:.3f}s diff={diff_sec:.3f}s",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 2: deltat_ex — Extended Delta-T
# ============================================================================


def test_part2_deltat_ex():
    print("\n" + "=" * 70)
    print("PART 2: deltat_ex — Extended Delta-T with Ephemeris Flag")
    print("=" * 70)

    epochs = [
        ("J2000.0", 2451545.0),
        ("2024-Jan", 2460310.5),
        ("1950", 2433282.5),
    ]

    for epoch_name, jd in epochs:
        test_name = f"P2/deltat_ex/{epoch_name}"
        try:
            # pyswisseph deltat_ex returns a float (despite docs suggesting tuple)
            dt_se = swe.deltat_ex(jd, swe.FLG_SWIEPH)
            if isinstance(dt_se, tuple):
                dt_se = dt_se[0]

            # libephemeris also returns a float
            dt_le = ephem.deltat_ex(jd, ephem.FLG_SWIEPH)
            if isinstance(dt_le, tuple):
                dt_le = dt_le[0]

            diff_sec = abs(dt_se - dt_le) * 86400.0

            ok = diff_sec < 0.5
            record(
                test_name,
                ok,
                f"SE={dt_se * 86400:.3f}s LE={dt_le * 86400:.3f}s diff={diff_sec:.3f}s",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 3: julday / revjul — Julian Day conversions
# ============================================================================


def test_part3_julday():
    print("\n" + "=" * 70)
    print("PART 3: julday / revjul — Julian Day Conversions")
    print("=" * 70)

    # Test dates (year, month, day, hour, calendar)
    test_dates = [
        (2000, 1, 1, 12.0, 1, "J2000.0"),
        (2024, 1, 15, 0.0, 1, "2024-Jan-15"),
        (1900, 1, 1, 0.0, 1, "1900-Jan-1"),
        (1582, 10, 15, 0.0, 1, "Greg_start"),
        (1582, 10, 4, 0.0, 0, "Julian_end"),
        (-4712, 1, 1, 12.0, 0, "JD_epoch"),
        (2100, 12, 31, 0.0, 1, "2100-Dec-31"),
        (1066, 10, 14, 12.0, 0, "Hastings"),
    ]

    for year, month, day, hour, cal, name in test_dates:
        test_name = f"P3/julday/{name}"
        try:
            jd_se = swe.julday(year, month, day, hour, cal)
            jd_le = ephem.julday(year, month, day, hour, cal)

            diff = abs(jd_se - jd_le)
            ok = diff < 1e-10
            record(
                test_name,
                ok,
                f"SE={jd_se:.6f} LE={jd_le:.6f} diff={diff:.2e}",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")

    # Reverse conversion tests
    jds = [2451545.0, 2460310.5, 2415020.5, 2299160.5]
    for jd in jds:
        test_name = f"P3/revjul/JD={jd}"
        try:
            result_se = swe.revjul(jd, 1)  # Gregorian
            result_le = ephem.revjul(jd, 1)

            # Both return (year, month, day, hour)
            ok = (
                result_se[0] == result_le[0]
                and result_se[1] == result_le[1]
                and result_se[2] == result_le[2]
                and abs(result_se[3] - result_le[3]) < 1e-8
            )

            record(
                test_name,
                ok,
                f"SE={result_se} LE={result_le}",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 4: utc_to_jd / jdut1_to_utc
# ============================================================================


def test_part4_utc():
    print("\n" + "=" * 70)
    print("PART 4: utc_to_jd / jdut1_to_utc")
    print("=" * 70)

    # Test UTC dates: (year, month, day, hour, min, sec)
    utc_dates = [
        (2024, 1, 15, 12, 0, 0.0, "2024-Jan-15T12:00"),
        (2000, 1, 1, 12, 0, 0.0, "J2000.0_utc"),
        (1990, 7, 15, 6, 30, 0.0, "1990-Jul-15"),
    ]

    for year, month, day, hour, minute, sec, name in utc_dates:
        test_name = f"P4/utc_to_jd/{name}"
        try:
            # pyswisseph: utc_to_jd(year, month, day, hour, min, sec, gregflag)
            # returns (jd_et, jd_ut1)
            result_se = swe.utc_to_jd(year, month, day, hour, minute, sec, 1)
            jd_et_se, jd_ut1_se = result_se[0], result_se[1]

            result_le = ephem.utc_to_jd(year, month, day, hour, minute, sec, 1)
            jd_et_le, jd_ut1_le = result_le[0], result_le[1]

            diff_et_sec = abs(jd_et_se - jd_et_le) * 86400.0
            diff_ut1_sec = abs(jd_ut1_se - jd_ut1_le) * 86400.0

            ok_et = diff_et_sec < 0.5
            ok_ut1 = diff_ut1_sec < 0.5

            record(
                f"{test_name}/ET",
                ok_et,
                f"SE={jd_et_se:.6f} LE={jd_et_le:.6f} diff={diff_et_sec:.3f}s",
            )
            record(
                f"{test_name}/UT1",
                ok_ut1,
                f"SE={jd_ut1_se:.6f} LE={jd_ut1_le:.6f} diff={diff_ut1_sec:.3f}s",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")

    # Test invalid leap second — both SE and LE should reject it
    test_name = "P4/utc_to_jd/leap_sec_reject"
    se_err = le_err = None
    try:
        swe.utc_to_jd(2024, 6, 30, 23, 59, 60.0, 1)
    except Exception as e:
        se_err = str(e)
    try:
        ephem.utc_to_jd(2024, 6, 30, 23, 59, 60.0, 1)
    except Exception as e:
        le_err = str(e)
    # Both should raise an error (no leap second at that date)
    ok = se_err is not None and le_err is not None
    record(
        test_name,
        ok,
        f"Both reject invalid leap second (SE: {se_err is not None}, LE: {le_err is not None})",
    )

    # Reverse: jdut1_to_utc
    jds = [2460310.5, 2451545.0]
    for jd in jds:
        test_name = f"P4/jdut1_to_utc/JD={jd}"
        try:
            result_se = swe.jdut1_to_utc(jd, 1)  # Gregorian
            result_le = ephem.jdut1_to_utc(jd, 1)

            # Both return (year, month, day, hour, min, sec)
            # Tolerance 0.1s for seconds — different Delta-T models produce
            # slightly different UT1-UTC offsets at the ~0.08s level
            ok = (
                result_se[0] == result_le[0]
                and result_se[1] == result_le[1]
                and result_se[2] == result_le[2]
                and result_se[3] == result_le[3]
                and result_se[4] == result_le[4]
                and abs(result_se[5] - result_le[5]) < 0.1
            )

            record(
                test_name,
                ok,
                f"SE=({result_se[0]},{result_se[1]},{result_se[2]},{result_se[3]}:{result_se[4]}:{result_se[5]:.2f}) "
                f"LE=({result_le[0]},{result_le[1]},{result_le[2]},{result_le[3]}:{result_le[4]}:{result_le[5]:.2f})",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 5: date_conversion
# ============================================================================


def test_part5_date_conversion():
    print("\n" + "=" * 70)
    print("PART 5: date_conversion")
    print("=" * 70)

    # Test valid dates
    valid_dates = [
        (2024, 1, 15, 0.0, b"g", "2024-Jan-15"),
        (2000, 1, 1, 12.0, b"g", "J2000.0"),
        (1582, 10, 15, 0.0, b"g", "Greg_start"),
        (1582, 10, 4, 0.0, b"j", "Julian_end"),
    ]

    for year, month, day, hour, cal, name in valid_dates:
        test_name = f"P5/date_conv/{name}"
        try:
            result_se = swe.date_conversion(year, month, day, hour, cal)
            result_le = ephem.date_conversion(year, month, day, hour, cal)

            # Both return (valid_bool, jd, (year, month, day, hour))
            valid_se, jd_se = result_se[0], result_se[1]
            valid_le, jd_le = result_le[0], result_le[1]

            diff = abs(jd_se - jd_le)
            ok = valid_se == valid_le and diff < 1e-10
            record(
                test_name,
                ok,
                f"valid={valid_se}/{valid_le} JD SE={jd_se:.6f} LE={jd_le:.6f} diff={diff:.2e}",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")

    # Test invalid date — both should return valid=False
    test_name = "P5/date_conv/invalid"
    try:
        result_se = swe.date_conversion(2024, 13, 1, 0.0, b"g")
        result_le = ephem.date_conversion(2024, 13, 1, 0.0, b"g")
        ok = result_se[0] == result_le[0] == False  # noqa: E712
        record(
            test_name,
            ok,
            f"SE valid={result_se[0]} LE valid={result_le[0]} (both reject invalid)",
        )
    except Exception as e:
        record(test_name, True, f"Error handling OK: {e}")


# ============================================================================
# PART 6: day_of_week
# ============================================================================


def test_part6_day_of_week():
    print("\n" + "=" * 70)
    print("PART 6: day_of_week")
    print("=" * 70)

    # Known days
    test_days = [
        (2460310.5, "2024-Jan-15", "Monday"),  # Known Monday
        (2451545.0, "J2000.0", "Saturday"),  # J2000 = Saturday
        (2440587.5, "Unix-epoch", "Thursday"),  # 1970-01-01 = Thursday
        (2415020.5, "1900-Jan-1", "Monday"),  # Known Monday
    ]

    for jd, name, expected_day in test_days:
        test_name = f"P6/dow/{name}"
        try:
            dow_se = swe.day_of_week(jd)
            dow_le = ephem.day_of_week(jd)

            ok = dow_se == dow_le
            days = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
            record(
                test_name,
                ok,
                f"SE={dow_se} LE={dow_le} expected={expected_day}",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 7: Delta-T sweep 1600-2100 — systematic drift detection
# ============================================================================


def test_part7_deltat_sweep():
    print("\n" + "=" * 70)
    print("PART 7: Delta-T Sweep 1600-2100 — Drift Detection")
    print("=" * 70)

    # Dense sweep 1900-2030 (where we expect close agreement)
    max_diff_sec = 0.0
    max_diff_year = 0
    diffs = []

    for year in range(1900, 2031, 5):
        jd = swe.julday(year, 1, 1, 0.0)
        dt_se = swe.deltat(jd)
        dt_le = ephem.deltat(jd)
        diff_sec = abs(dt_se - dt_le) * 86400.0
        diffs.append((year, diff_sec))
        if diff_sec > max_diff_sec:
            max_diff_sec = diff_sec
            max_diff_year = year

    ok = max_diff_sec < 1.0
    record(
        "P7/sweep/1900-2030",
        ok,
        f"max_diff={max_diff_sec:.3f}s at year={max_diff_year}",
    )

    # Historical sweep 1600-1900
    max_diff_hist = 0.0
    max_diff_hist_year = 0
    for year in range(1600, 1901, 25):
        jd = swe.julday(year, 1, 1, 0.0)
        dt_se = swe.deltat(jd)
        dt_le = ephem.deltat(jd)
        diff_sec = abs(dt_se - dt_le) * 86400.0
        if diff_sec > max_diff_hist:
            max_diff_hist = diff_sec
            max_diff_hist_year = year

    ok = max_diff_hist < 30.0  # Historical can diverge more
    record(
        "P7/sweep/1600-1900",
        ok,
        f"max_diff={max_diff_hist:.3f}s at year={max_diff_hist_year}",
    )

    # Future sweep 2030-2100
    max_diff_future = 0.0
    max_diff_future_year = 0
    for year in range(2030, 2101, 5):
        jd = swe.julday(year, 1, 1, 0.0)
        dt_se = swe.deltat(jd)
        dt_le = ephem.deltat(jd)
        diff_sec = abs(dt_se - dt_le) * 86400.0
        if diff_sec > max_diff_future:
            max_diff_future = diff_sec
            max_diff_future_year = year

    ok = max_diff_future < 30.0
    record(
        "P7/sweep/2030-2100",
        ok,
        f"max_diff={max_diff_future:.3f}s at year={max_diff_future_year}",
    )


# ============================================================================
# PART 8: TT-UT consistency
# ============================================================================


def test_part8_tt_ut_consistency():
    print("\n" + "=" * 70)
    print("PART 8: TT-UT Consistency (calc_ut vs calc with Delta-T)")
    print("=" * 70)

    jd_ut = 2460310.5
    flags = swe.FLG_SPEED

    for body_id, body_name in [
        (swe.SUN, "Sun"),
        (swe.MOON, "Moon"),
        (swe.MARS, "Mars"),
    ]:
        test_name = f"P8/tt_ut/{body_name}"
        try:
            # Get position via calc_ut (UT input)
            pos_ut_le, _ = ephem.calc_ut(jd_ut, body_id, flags)

            # Get Delta-T, then use calc (TT input)
            dt_le = ephem.deltat(jd_ut)
            jd_tt = jd_ut + dt_le
            pos_tt_le, _ = ephem.calc(jd_tt, body_id, flags)

            # Both should give identical results
            lon_diff = abs(pos_ut_le[0] - pos_tt_le[0])
            if lon_diff > 180:
                lon_diff = 360 - lon_diff

            ok = lon_diff < 0.00001  # Should be essentially identical
            record(
                test_name,
                ok,
                f"calc_ut={pos_ut_le[0]:.6f} calc(TT)={pos_tt_le[0]:.6f} "
                f'diff={lon_diff:.8f}° ({lon_diff * 3600:.4f}")',
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# MAIN
# ============================================================================


def main():
    global total, passed, failed, skipped, failures

    print("=" * 70)
    print("ROUND 11: Deep Delta-T and Time Functions Audit")
    print("=" * 70)
    t0 = time.time()

    test_part1_deltat()
    test_part2_deltat_ex()
    test_part3_julday()
    test_part4_utc()
    test_part5_date_conversion()
    test_part6_day_of_week()
    test_part7_deltat_sweep()
    test_part8_tt_ut_consistency()

    elapsed = time.time() - t0

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total:   {total}")
    print(f"Passed:  {passed}")
    print(f"Failed:  {failed}")
    print(f"Skipped: {skipped}")
    print(f"Time:    {elapsed:.1f}s")

    if failures:
        print(f"\n--- {len(failures)} FAILURES ---")
        for name, detail in failures:
            print(f"  {name}: {detail}")

    print(f"\nPass rate: {passed}/{total} = {100 * passed / max(total, 1):.1f}%")

    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
