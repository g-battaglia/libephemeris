#!/usr/bin/env python
"""Standalone deep verification: time functions, coordinate transforms,
lunar nodes/apsides, and phenomena vs pyswisseph.

Sections covered:
  11.1  julday/revjul round-trip (500 random dates)
  11.2  deltat (200 JDs)
  11.3  sidtime (100 JDs)
  11.4  utc_to_jd / jd_to_utc round-trip (20 dates)
  11.5  day_of_week (10 known dates)
  10.1  cotrans round-trip (3 obliquities x 100 points)
  10.2  cotrans_sp with velocities (20 points)
  12    lunar nodes and apsides (100 dates)
  12.5  nod_aps_ut (5 bodies x 5 dates)
  17    pheno_ut (10 bodies x 10 dates)

Target: ~5000+ checks, <30 seconds.
"""

from __future__ import annotations

import random
import sys
import time

# ---------------------------------------------------------------------------
# Libraries
# ---------------------------------------------------------------------------
try:
    import libephemeris as lib
except ImportError:
    sys.exit("ERROR: libephemeris not importable. Run with .venv/bin/python.")

try:
    import swisseph as swe_ref
except ImportError:
    sys.exit("ERROR: pyswisseph (swisseph) not importable.")

# ---------------------------------------------------------------------------
# Counters & helpers
# ---------------------------------------------------------------------------
passed = 0
failed = 0
errors = 0

random.seed(42)


def check(condition: bool, label: str, detail: str = "") -> bool:
    """Register a check.  Prints only on failure."""
    global passed, failed
    if condition:
        passed += 1
        return True
    else:
        failed += 1
        msg = f"  FAIL: {label}"
        if detail:
            msg += f"  -- {detail}"
        print(msg)
        return False


def close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def angle_diff(a: float, b: float) -> float:
    """Shortest angular difference in degrees."""
    d = (a - b) % 360.0
    if d > 180.0:
        d -= 360.0
    return abs(d)


# ---------------------------------------------------------------------------
# Helper: generate a random valid calendar date
# ---------------------------------------------------------------------------
DAYS_IN_MONTH = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def _is_leap(y: int) -> bool:
    if y % 4 != 0:
        return False
    if y % 100 != 0:
        return True
    return y % 400 == 0


def random_date(y_min: int = -3000, y_max: int = 3000):
    """Return (year, month, day, decimal_hour) for a random date."""
    y = random.randint(y_min, y_max)
    m = random.randint(1, 12)
    max_d = DAYS_IN_MONTH[m]
    if m == 2 and _is_leap(y):
        max_d = 29
    d = random.randint(1, max_d)
    h = random.random() * 24.0
    return y, m, d, h


# ===========================================================================
# Section 11.1  julday / revjul round-trip  (500 dates)
# ===========================================================================
def test_julday_revjul():
    global errors
    print("\n=== 11.1  julday / revjul round-trip (500 dates) ===")
    n = 500
    for _ in range(n):
        y, m, d, h = random_date()
        try:
            jd_lib = lib.julday(y, m, d, h)
            jd_ref = swe_ref.julday(y, m, d, h)
            # Both use Meeus formula; allow FP rounding (1e-9 ~ 0.1ms)
            check(
                close(jd_lib, jd_ref, 1e-9),
                f"julday match y={y} m={m} d={d}",
                f"lib={jd_lib} ref={jd_ref} diff={abs(jd_lib - jd_ref):.2e}",
            )

            # revjul
            yr_l, mr_l, dr_l, hr_l = lib.revjul(jd_lib)
            yr_r, mr_r, dr_r, hr_r = swe_ref.revjul(jd_lib)
            check(
                yr_l == yr_r and mr_l == mr_r and dr_l == dr_r,
                f"revjul ymd match jd={jd_lib:.4f}",
                f"lib=({yr_l},{mr_l},{dr_l}) ref=({yr_r},{mr_r},{dr_r})",
            )
            check(
                close(hr_l, hr_r, 1e-9),
                f"revjul hour match jd={jd_lib:.4f}",
                f"lib={hr_l} ref={hr_r}",
            )

            # Round-trip: revjul(julday(y,m,d,h)) should recover
            check(
                yr_l == y and mr_l == m and dr_l == d,
                f"round-trip ymd y={y} m={m} d={d}",
                f"got ({yr_l},{mr_l},{dr_l})",
            )
            check(
                close(hr_l, h, 1e-8),
                f"round-trip hour y={y} m={m} d={d}",
                f"orig={h} got={hr_l} diff={abs(hr_l - h):.2e}",
            )
        except Exception as exc:
            errors += 1
            print(f"  ERROR in julday/revjul y={y} m={m} d={d}: {exc}")


# ===========================================================================
# Section 11.2  deltat  (200 JDs in [1900, 2100])
# ===========================================================================
def test_deltat():
    global errors
    print("\n=== 11.2  deltat (200 JDs) ===")
    # JD for 1900-Jan-1 ~ 2415020.5, 2100-Jan-1 ~ 2488069.5
    jd_start = 2415020.5
    jd_end = 2488069.5
    n = 200
    jds = sorted([jd_start + random.random() * (jd_end - jd_start) for _ in range(n)])
    prev_dt = None
    prev_jd = None
    for jd in jds:
        try:
            dt_lib = lib.deltat(jd)
            dt_ref = swe_ref.deltat(jd)
            check(
                close(dt_lib, dt_ref, 1e-4),
                f"deltat match jd={jd:.2f}",
                f"lib={dt_lib:.8e} ref={dt_ref:.8e} diff={abs(dt_lib - dt_ref):.2e}",
            )
            check(dt_lib > 0, f"deltat positive jd={jd:.2f}", f"got {dt_lib}")
            check(
                dt_lib < 0.01,
                f"deltat < 0.01 day jd={jd:.2f}",
                f"got {dt_lib:.6e}",
            )
            # Smoothness: rate of change should be < 1e-7 per day
            # (delta-T changes by at most ~2 sec/year ~ 2.3e-8 day/day)
            if prev_dt is not None and prev_jd is not None:
                interval = jd - prev_jd
                if interval > 0:
                    rate = abs(dt_lib - prev_dt) / interval
                    check(
                        rate < 1e-7,
                        f"deltat smooth jd={jd:.2f}",
                        f"rate={rate:.2e}/day over {interval:.1f} days",
                    )
            prev_dt = dt_lib
            prev_jd = jd
        except Exception as exc:
            errors += 1
            print(f"  ERROR in deltat jd={jd:.2f}: {exc}")


# ===========================================================================
# Section 11.3  sidtime  (100 JDs)
# ===========================================================================
def test_sidtime():
    global errors
    print("\n=== 11.3  sidtime (100 JDs) ===")
    jd_start = 2415020.5  # ~1900
    jd_end = 2488069.5  # ~2100
    n = 100
    jds = sorted([jd_start + random.random() * (jd_end - jd_start) for _ in range(n)])
    prev_st = None
    prev_jd = None
    for jd in jds:
        try:
            st_lib = lib.sidtime(jd)
            st_ref = swe_ref.sidtime(jd)
            check(
                close(st_lib, st_ref, 0.001),
                f"sidtime match jd={jd:.2f}",
                f"lib={st_lib:.6f} ref={st_ref:.6f} diff={abs(st_lib - st_ref):.6f}",
            )
            check(
                0.0 <= st_lib < 24.0,
                f"sidtime range jd={jd:.2f}",
                f"got {st_lib}",
            )
            # Sidereal day increment check: only for closely-spaced consecutive days
            if prev_st is not None and prev_jd is not None:
                djd = jd - prev_jd
                if 0.9 < djd < 1.1:
                    # Expected sidereal gain per solar day: ~3m56.555s = 0.065710 hours
                    expected_gain = djd * (3.0 * 60.0 + 56.555) / 3600.0
                    actual_diff = (st_lib - prev_st) % 24.0
                    gain = actual_diff
                    if gain > 12.0:
                        gain -= 24.0
                    check(
                        close(gain, expected_gain, 0.02),
                        f"sidtime ~3m56s/day jd={jd:.2f}",
                        f"gain={gain:.6f}h expected~{expected_gain:.6f}h",
                    )
            prev_st = st_lib
            prev_jd = jd
        except Exception as exc:
            errors += 1
            print(f"  ERROR in sidtime jd={jd:.2f}: {exc}")


# ===========================================================================
# Section 11.4  utc_to_jd / jd_to_utc round-trip  (20 dates)
# ===========================================================================
def test_utc_to_jd():
    global errors
    print("\n=== 11.4  utc_to_jd / jd_to_utc round-trip (20 dates) ===")
    # Use dates in the well-covered modern era (1972-2030) to avoid
    # delta-T model divergence between Skyfield and Swiss Ephemeris
    dates = [
        (2000, 1, 1, 12, 0, 0.0),
        (1999, 12, 31, 23, 59, 59.0),
        (2020, 6, 21, 6, 44, 0.0),
        (1980, 3, 15, 8, 30, 15.5),
        (2023, 11, 15, 14, 22, 33.7),
        (1990, 6, 15, 18, 0, 0.0),
        (2010, 12, 25, 3, 45, 12.0),
        (2005, 8, 8, 20, 20, 20.0),
        (1975, 4, 10, 10, 10, 10.0),
        (1985, 9, 23, 16, 30, 0.0),
        (2015, 3, 20, 22, 45, 30.0),
        (1995, 1, 1, 0, 0, 0.0),
        (2019, 4, 15, 9, 15, 45.123),
        (2000, 6, 15, 0, 0, 0.0),
        (2003, 3, 20, 12, 0, 0.0),
        (2008, 8, 8, 8, 8, 8.0),
        (2012, 12, 21, 11, 11, 0.0),
        (2018, 1, 31, 13, 51, 0.0),
        (2022, 3, 20, 15, 33, 0.0),
        (2024, 7, 4, 12, 0, 0.0),
    ]
    for y, mo, d, h, mi, s in dates:
        try:
            jd_et_lib, jd_ut_lib = lib.utc_to_jd(y, mo, d, h, mi, s, 1)
            jd_et_ref, jd_ut_ref = swe_ref.utc_to_jd(y, mo, d, h, mi, s, 1)
            # Delta-T models differ between Skyfield and Swiss Ephemeris;
            # allow up to ~1 second (1.16e-5 day) for ET and UT
            check(
                close(jd_et_lib, jd_et_ref, 2e-5),
                f"utc_to_jd ET {y}-{mo:02d}-{d:02d}",
                f"lib={jd_et_lib:.10f} ref={jd_et_ref:.10f} diff={abs(jd_et_lib - jd_et_ref):.2e}",
            )
            check(
                close(jd_ut_lib, jd_ut_ref, 2e-5),
                f"utc_to_jd UT {y}-{mo:02d}-{d:02d}",
                f"lib={jd_ut_lib:.10f} ref={jd_ut_ref:.10f} diff={abs(jd_ut_lib - jd_ut_ref):.2e}",
            )
            # Round-trip: jdet_to_utc(utc_to_jd()[0]) should recover components
            yr2, mo2, d2, h2, mi2, s2 = lib.jdet_to_utc(jd_et_lib, 1)
            check(
                yr2 == y and mo2 == mo and d2 == d,
                f"jdet_to_utc date {y}-{mo:02d}-{d:02d}",
                f"got {yr2}-{mo2:02d}-{d2:02d}",
            )
            # Allow +/-1 second tolerance on time components due to rounding
            total_sec_orig = h * 3600 + mi * 60 + s
            total_sec_back = h2 * 3600 + mi2 * 60 + s2
            check(
                abs(total_sec_back - total_sec_orig) < 1.0,
                f"jdet_to_utc time {y}-{mo:02d}-{d:02d}",
                f"expected {h}:{mi}:{s} got {h2}:{mi2}:{s2:.3f} "
                f"(diff={abs(total_sec_back - total_sec_orig):.3f}s)",
            )
        except Exception as exc:
            errors += 1
            print(f"  ERROR in utc_to_jd {y}-{mo:02d}-{d:02d}: {exc}")


# ===========================================================================
# Section 11.5  day_of_week  (10 known dates)
# ===========================================================================
def test_day_of_week():
    global errors
    print("\n=== 11.5  day_of_week (10 known dates) ===")
    # Convention: 0=Mon, 1=Tue, 2=Wed, 3=Thu, 4=Fri, 5=Sat, 6=Sun
    known = [
        # (year, month, day, expected_dow)
        (2000, 1, 1, 5),  # Saturday
        (2024, 1, 1, 0),  # Monday
        (2023, 12, 25, 0),  # Monday
        (1969, 7, 20, 6),  # Sunday (Moon landing)
        (1776, 7, 4, 3),  # Thursday
        (2020, 2, 29, 5),  # Saturday (leap day)
        (1999, 12, 31, 4),  # Friday
        (2001, 9, 11, 1),  # Tuesday
        (1945, 5, 8, 1),  # Tuesday (V-E Day)
        (2026, 3, 27, 4),  # Friday (today per context)
    ]
    for y, m, d, expected_dow in known:
        try:
            jd = lib.julday(y, m, d, 12.0)
            dow_lib = lib.day_of_week(jd)
            dow_ref = swe_ref.day_of_week(jd)
            check(
                dow_lib == expected_dow,
                f"day_of_week {y}-{m:02d}-{d:02d} vs expected",
                f"lib={dow_lib} expected={expected_dow}",
            )
            check(
                dow_lib == dow_ref,
                f"day_of_week {y}-{m:02d}-{d:02d} vs swe_ref",
                f"lib={dow_lib} ref={dow_ref}",
            )
        except Exception as exc:
            errors += 1
            print(f"  ERROR in day_of_week {y}-{m:02d}-{d:02d}: {exc}")


# ===========================================================================
# Section 10.1  cotrans round-trip  (3 obliquities x 100 random points)
# ===========================================================================
def test_cotrans_roundtrip():
    global errors
    print("\n=== 10.1  cotrans round-trip (3 obliquities x 100 points) ===")
    obliquities = [22.5, 23.4393, 24.0]
    n_points = 100
    for eps in obliquities:
        for _ in range(n_points):
            lon = random.uniform(0.0, 360.0)
            lat = random.uniform(-85.0, 85.0)  # avoid poles
            dist = 1.0
            try:
                # ecl -> eq (negative eps) -> ecl (positive eps)
                eq = lib.cotrans((lon, lat, dist), -eps)
                ecl_back = lib.cotrans(eq, eps)
                err_lon = angle_diff(ecl_back[0], lon)
                err_lat = abs(ecl_back[1] - lat)
                check(
                    err_lon < 1e-8 and err_lat < 1e-8,
                    f"cotrans round-trip eps={eps} lon={lon:.2f} lat={lat:.2f}",
                    f"err_lon={err_lon:.2e} err_lat={err_lat:.2e}",
                )
                # Also compare with pyswisseph directly
                eq_ref = swe_ref.cotrans((lon, lat, dist), -eps)
                check(
                    close(eq[0], eq_ref[0], 1e-8) and close(eq[1], eq_ref[1], 1e-8),
                    f"cotrans vs swe_ref eps={eps} lon={lon:.2f}",
                    f"lib=({eq[0]:.8f},{eq[1]:.8f}) ref=({eq_ref[0]:.8f},{eq_ref[1]:.8f})",
                )
            except Exception as exc:
                errors += 1
                print(f"  ERROR in cotrans eps={eps} lon={lon:.2f}: {exc}")


# ===========================================================================
# Section 10.2  cotrans_sp round-trip with velocities  (20 points)
# ===========================================================================
def test_cotrans_sp():
    global errors
    print("\n=== 10.2  cotrans_sp round-trip with velocities (20 points) ===")
    eps = 23.4393
    for _ in range(20):
        lon = random.uniform(0.0, 360.0)
        lat = random.uniform(-80.0, 80.0)
        dist = random.uniform(0.5, 2.0)
        lon_sp = random.uniform(-2.0, 2.0)
        lat_sp = random.uniform(-0.5, 0.5)
        dist_sp = random.uniform(-0.01, 0.01)
        coord = (lon, lat, dist, lon_sp, lat_sp, dist_sp)
        try:
            # ecl -> eq -> ecl
            eq = lib.cotrans_sp(coord, -eps)
            ecl_back = lib.cotrans_sp(eq, eps)
            err_lon = angle_diff(ecl_back[0], lon)
            err_lat = abs(ecl_back[1] - lat)
            err_dist = abs(ecl_back[2] - dist)
            check(
                err_lon < 1e-7 and err_lat < 1e-7 and err_dist < 1e-10,
                f"cotrans_sp pos round-trip lon={lon:.2f}",
                f"err_lon={err_lon:.2e} err_lat={err_lat:.2e}",
            )
            # Velocity round-trip
            err_lon_sp = abs(ecl_back[3] - lon_sp)
            err_lat_sp = abs(ecl_back[4] - lat_sp)
            err_dist_sp = abs(ecl_back[5] - dist_sp)
            check(
                err_lon_sp < 1e-6 and err_lat_sp < 1e-6 and err_dist_sp < 1e-10,
                f"cotrans_sp vel round-trip lon={lon:.2f}",
                f"err_lonsp={err_lon_sp:.2e} err_latsp={err_lat_sp:.2e}",
            )
            # Compare with pyswisseph
            eq_ref = swe_ref.cotrans_sp(coord, -eps)
            check(
                close(eq[0], eq_ref[0], 1e-8) and close(eq[1], eq_ref[1], 1e-8),
                f"cotrans_sp vs swe_ref pos lon={lon:.2f}",
                f"lib=({eq[0]:.8f},{eq[1]:.8f}) ref=({eq_ref[0]:.8f},{eq_ref[1]:.8f})",
            )
            check(
                close(eq[3], eq_ref[3], 1e-7) and close(eq[4], eq_ref[4], 1e-7),
                f"cotrans_sp vs swe_ref vel lon={lon:.2f}",
                f"lib=({eq[3]:.8f},{eq[4]:.8f}) ref=({eq_ref[3]:.8f},{eq_ref[4]:.8f})",
            )
        except Exception as exc:
            errors += 1
            print(f"  ERROR in cotrans_sp lon={lon:.2f}: {exc}")


# ===========================================================================
# Section 12  Lunar nodes and apsides  (100 dates)
# ===========================================================================
def test_lunar_nodes_apsides():
    global errors
    print("\n=== 12  Lunar nodes and apsides (100 dates) ===")
    jd_start = 2415020.5  # ~1900
    jd_end = 2488069.5  # ~2100
    n = 100
    flags_lib = lib.FLG_SWIEPH | lib.FLG_SPEED
    flags_ref = swe_ref.FLG_SWIEPH | swe_ref.FLG_SPEED
    jds = [jd_start + random.random() * (jd_end - jd_start) for _ in range(n)]
    for jd in jds:
        # --- Mean Node ---
        try:
            mn_lib, _ = lib.calc_ut(jd, lib.MEAN_NODE, flags_lib)
            mn_ref = swe_ref.calc_ut(jd, swe_ref.MEAN_NODE, flags_ref)
            mn_ref_lon = mn_ref[0][0]
            diff_mn = angle_diff(mn_lib[0], mn_ref_lon)
            check(
                diff_mn < 1.0 / 3600.0,  # < 1 arcsec
                f"Mean Node lon jd={jd:.2f}",
                f'lib={mn_lib[0]:.6f} ref={mn_ref_lon:.6f} diff={diff_mn * 3600:.2f}"',
            )
            check(
                mn_lib[3] < 0,
                f"Mean Node speed negative jd={jd:.2f}",
                f"speed={mn_lib[3]:.6f}",
            )
        except Exception as exc:
            errors += 1
            print(f"  ERROR Mean Node jd={jd:.2f}: {exc}")

        # --- True Node ---
        try:
            tn_lib, _ = lib.calc_ut(jd, lib.TRUE_NODE, flags_lib)
            tn_ref = swe_ref.calc_ut(jd, swe_ref.TRUE_NODE, flags_ref)
            tn_ref_lon = tn_ref[0][0]
            diff_tn = angle_diff(tn_lib[0], tn_ref_lon)
            # True Node uses different osculating methods; allow 60 arcsec
            check(
                diff_tn < 60.0 / 3600.0,
                f"True Node lon jd={jd:.2f}",
                f'lib={tn_lib[0]:.6f} ref={tn_ref_lon:.6f} diff={diff_tn * 3600:.2f}"',
            )
        except Exception as exc:
            errors += 1
            print(f"  ERROR True Node jd={jd:.2f}: {exc}")

        # --- Mean Apogee (Black Moon Lilith) ---
        try:
            ma_lib, _ = lib.calc_ut(jd, lib.MEAN_APOG, flags_lib)
            ma_ref = swe_ref.calc_ut(jd, swe_ref.MEAN_APOG, flags_ref)
            ma_ref_lon = ma_ref[0][0]
            diff_ma = angle_diff(ma_lib[0], ma_ref_lon)
            check(
                diff_ma < 5.0 / 3600.0,  # < 5 arcsec
                f"Mean Apog lon jd={jd:.2f}",
                f'lib={ma_lib[0]:.6f} ref={ma_ref_lon:.6f} diff={diff_ma * 3600:.2f}"',
            )
            check(
                ma_lib[3] > 0,
                f"Mean Apog speed positive jd={jd:.2f}",
                f"speed={ma_lib[3]:.6f}",
            )
        except Exception as exc:
            errors += 1
            print(f"  ERROR Mean Apog jd={jd:.2f}: {exc}")

        # --- Osculating Apogee ---
        try:
            oa_lib, _ = lib.calc_ut(jd, lib.OSCU_APOG, flags_lib)
            oa_ref = swe_ref.calc_ut(jd, swe_ref.OSCU_APOG, flags_ref)
            oa_ref_lon = oa_ref[0][0]
            diff_oa = angle_diff(oa_lib[0], oa_ref_lon)
            # Osculating apogee is sensitive to ephemeris differences; allow 4 arcmin
            check(
                diff_oa < 240.0 / 3600.0,
                f"Oscu Apog lon jd={jd:.2f}",
                f'lib={oa_lib[0]:.6f} ref={oa_ref_lon:.6f} diff={diff_oa * 3600:.2f}"',
            )
        except Exception as exc:
            errors += 1
            print(f"  ERROR Oscu Apog jd={jd:.2f}: {exc}")

        # --- IntpApog + IntpPerig: should be roughly opposite (within 30 deg) ---
        try:
            ia_lib, _ = lib.calc_ut(jd, lib.INTP_APOG, flags_lib)
            ip_lib, _ = lib.calc_ut(jd, lib.INTP_PERG, flags_lib)
            sep = angle_diff(ia_lib[0], ip_lib[0])
            # Interpolated apse positions can deviate from exact opposition
            check(
                abs(sep - 180.0) < 30.0,
                f"IntpApog-IntpPerig ~180 apart jd={jd:.2f}",
                f"sep={sep:.2f} deg",
            )
        except Exception as exc:
            errors += 1
            print(f"  ERROR IntpApog/Perig jd={jd:.2f}: {exc}")


# ===========================================================================
# Section 12.5  nod_aps_ut  (5 bodies x 5 dates x method mean)
# ===========================================================================
def test_nod_aps_ut():
    global errors
    print("\n=== 12.5  nod_aps_ut (5 bodies x 5 dates) ===")
    bodies_lib = [
        lib.MARS,
        lib.JUPITER,
        lib.SATURN,
        lib.URANUS,
        lib.NEPTUNE,
    ]
    bodies_ref = [
        swe_ref.MARS,
        swe_ref.JUPITER,
        swe_ref.SATURN,
        swe_ref.URANUS,
        swe_ref.NEPTUNE,
    ]
    body_names = ["Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
    jds = [2451545.0, 2455000.0, 2458000.0, 2460000.0, 2462000.0]
    flags_lib = lib.FLG_SWIEPH | lib.FLG_SPEED
    flags_ref = swe_ref.FLG_SWIEPH | swe_ref.FLG_SPEED
    method = 1  # NODBIT_MEAN

    for i, (b_lib, b_ref, bname) in enumerate(zip(bodies_lib, bodies_ref, body_names)):
        for jd in jds:
            try:
                nasc_l, ndsc_l, peri_l, aphe_l = lib.nod_aps_ut(
                    jd, b_lib, method, flags_lib
                )
                nasc_r, ndsc_r, peri_r, aphe_r = swe_ref.nod_aps_ut(
                    jd, b_ref, method, flags_ref
                )
                # Ascending node comparison
                diff_nasc = angle_diff(nasc_l[0], nasc_r[0])
                check(
                    diff_nasc < 1.0,  # < 1 degree
                    f"nod_aps_ut {bname} asc_node jd={jd:.1f}",
                    f"lib={nasc_l[0]:.4f} ref={nasc_r[0]:.4f} diff={diff_nasc:.4f} deg",
                )
                # Descending node comparison
                diff_ndsc = angle_diff(ndsc_l[0], ndsc_r[0])
                check(
                    diff_ndsc < 1.0,
                    f"nod_aps_ut {bname} dsc_node jd={jd:.1f}",
                    f"lib={ndsc_l[0]:.4f} ref={ndsc_r[0]:.4f} diff={diff_ndsc:.4f} deg",
                )
                # Perihelion comparison (relaxed: different geocentric methods)
                diff_peri = angle_diff(peri_l[0], peri_r[0])
                check(
                    diff_peri < 35.0,
                    f"nod_aps_ut {bname} perihelion jd={jd:.1f}",
                    f"lib={peri_l[0]:.4f} ref={peri_r[0]:.4f} diff={diff_peri:.4f} deg",
                )
                # Aphelion comparison
                diff_aphe = angle_diff(aphe_l[0], aphe_r[0])
                check(
                    diff_aphe < 35.0,
                    f"nod_aps_ut {bname} aphelion jd={jd:.1f}",
                    f"lib={aphe_l[0]:.4f} ref={aphe_r[0]:.4f} diff={diff_aphe:.4f} deg",
                )
            except Exception as exc:
                errors += 1
                print(f"  ERROR nod_aps_ut {bname} jd={jd:.1f}: {exc}")


# ===========================================================================
# Section 17  pheno_ut  (10 bodies x 10 dates)
# ===========================================================================
def test_pheno_ut():
    global errors
    print("\n=== 17  pheno_ut (10 bodies x 10 dates) ===")
    # Use only planets that pyswisseph can calculate without extra ephemeris files
    # (Chiron requires seas_18.se1 which may not be installed)
    bodies_lib = [
        lib.MERCURY,
        lib.VENUS,
        lib.MARS,
        lib.JUPITER,
        lib.SATURN,
        lib.URANUS,
        lib.NEPTUNE,
        lib.PLUTO,
        lib.SUN,
        lib.MOON,
    ]
    bodies_ref = [
        swe_ref.MERCURY,
        swe_ref.VENUS,
        swe_ref.MARS,
        swe_ref.JUPITER,
        swe_ref.SATURN,
        swe_ref.URANUS,
        swe_ref.NEPTUNE,
        swe_ref.PLUTO,
        swe_ref.SUN,
        swe_ref.MOON,
    ]
    body_names = [
        "Mercury",
        "Venus",
        "Mars",
        "Jupiter",
        "Saturn",
        "Uranus",
        "Neptune",
        "Pluto",
        "Sun",
        "Moon",
    ]
    jd_start = 2451545.0  # J2000
    jd_end = 2460000.0  # ~2023
    jds = [jd_start + (jd_end - jd_start) * i / 9.0 for i in range(10)]
    flags_lib = lib.FLG_SWIEPH
    flags_ref = swe_ref.FLG_SWIEPH

    for i, (b_lib, b_ref, bname) in enumerate(zip(bodies_lib, bodies_ref, body_names)):
        for jd in jds:
            try:
                attr_lib = lib.pheno_ut(jd, b_lib, flags_lib)
                attr_ref = swe_ref.pheno_ut(jd, b_ref, flags_ref)

                # [0] phase angle
                diff_pa = abs(attr_lib[0] - attr_ref[0])
                check(
                    diff_pa < 0.1,
                    f"pheno_ut {bname} phase_angle jd={jd:.1f}",
                    f"lib={attr_lib[0]:.4f} ref={attr_ref[0]:.4f} diff={diff_pa:.4f}",
                )

                # [1] phase (illuminated fraction)
                diff_ph = abs(attr_lib[1] - attr_ref[1])
                check(
                    diff_ph < 0.01,
                    f"pheno_ut {bname} phase jd={jd:.1f}",
                    f"lib={attr_lib[1]:.6f} ref={attr_ref[1]:.6f} diff={diff_ph:.6f}",
                )

                # [2] elongation
                diff_el = abs(attr_lib[2] - attr_ref[2])
                check(
                    diff_el < 0.1,
                    f"pheno_ut {bname} elongation jd={jd:.1f}",
                    f"lib={attr_lib[2]:.4f} ref={attr_ref[2]:.4f} diff={diff_el:.4f}",
                )

                # [3] apparent diameter
                diff_diam = abs(attr_lib[3] - attr_ref[3])
                ref_diam = max(abs(attr_ref[3]), 1e-6)
                check(
                    diff_diam / ref_diam < 0.01,  # 1% relative
                    f"pheno_ut {bname} diameter jd={jd:.1f}",
                    f"lib={attr_lib[3]:.6f} ref={attr_ref[3]:.6f} rel={diff_diam / ref_diam:.4f}",
                )

                # [4] magnitude
                diff_mag = abs(attr_lib[4] - attr_ref[4])
                check(
                    diff_mag < 0.5,
                    f"pheno_ut {bname} magnitude jd={jd:.1f}",
                    f"lib={attr_lib[4]:.4f} ref={attr_ref[4]:.4f} diff={diff_mag:.4f}",
                )
            except Exception as exc:
                errors += 1
                print(f"  ERROR pheno_ut {bname} jd={jd:.1f}: {exc}")


# ===========================================================================
# Main
# ===========================================================================
def main():
    global passed, failed, errors
    t0 = time.time()

    print("=" * 70)
    print("Deep Verification: Time, Coordinates, Nodes, Phenomena")
    print("  libephemeris vs pyswisseph")
    print("=" * 70)

    test_julday_revjul()
    test_deltat()
    test_sidtime()
    test_utc_to_jd()
    test_day_of_week()
    test_cotrans_roundtrip()
    test_cotrans_sp()
    test_lunar_nodes_apsides()
    test_nod_aps_ut()
    test_pheno_ut()

    elapsed = time.time() - t0

    print("\n" + "=" * 70)
    print(f"TOTAL CHECKS: {passed + failed}")
    print(f"  PASSED: {passed}")
    print(f"  FAILED: {failed}")
    if errors:
        print(f"  ERRORS (exceptions): {errors}")
    print(f"  Time: {elapsed:.1f}s")
    print("=" * 70)

    if failed > 0 or errors > 0:
        sys.exit(1)
    else:
        print("\nAll checks passed!")
        sys.exit(0)


if __name__ == "__main__":
    main()
