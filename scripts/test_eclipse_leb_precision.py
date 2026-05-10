#!/usr/bin/env python3
"""Compare eclipse calculation precision between LEB and Skyfield paths.

Runs each eclipse function twice — once with LEB active, once forcing Skyfield —
and measures the difference. This validates that the LEB porting doesn't
introduce precision regressions.

Usage:
    uv run python scripts/test_eclipse_leb_precision.py
"""

from __future__ import annotations

import math
import sys
import time

import libephemeris
from libephemeris.constants import (
    SE_CALC_RISE,
    SE_CALC_SET,
    SE_MOON,
    SE_SUN,
    SEFLG_SPEED,
)
from libephemeris.time_utils import swe_julday


def _run_with_mode(mode: str, func, *args, **kwargs):
    """Run a function in a specific calc mode with clean session state."""
    from libephemeris import state

    old_mode = libephemeris.get_calc_mode()
    old_planets = state._PLANETS
    state._PLANETS = None
    libephemeris.set_calc_mode(mode)
    try:
        return func(*args, **kwargs)
    finally:
        libephemeris.set_calc_mode(old_mode)
        state._PLANETS = old_planets


def compare_besselian_elements():
    """Compare Besselian elements LEB vs Skyfield for known eclipses."""
    from libephemeris.eclipse import (
        _calc_gamma,
        _calc_penumbra_limit,
        _calc_umbra_limit,
        calc_besselian_d,
        calc_besselian_l1,
        calc_besselian_l2,
        calc_besselian_mu,
        calc_besselian_x,
        calc_besselian_y,
    )

    eclipses = [
        ("2024-04-08 Total", swe_julday(2024, 4, 8, 18.17)),
        ("2023-10-14 Annular", swe_julday(2023, 10, 14, 17.97)),
        ("2025-03-29 Partial", swe_julday(2025, 3, 29, 10.5)),
        ("2024-10-02 Annular", swe_julday(2024, 10, 2, 18.5)),
        ("2026-02-17 Annular", swe_julday(2026, 2, 17, 12.0)),
    ]

    funcs = [
        ("gamma", _calc_gamma),
        ("penumbra", _calc_penumbra_limit),
        ("umbra", _calc_umbra_limit),
        ("x", calc_besselian_x),
        ("y", calc_besselian_y),
        ("d", calc_besselian_d),
        ("l1", calc_besselian_l1),
        ("l2", calc_besselian_l2),
        ("mu", calc_besselian_mu),
    ]

    print("=" * 70)
    print("BESSELIAN ELEMENTS: LEB vs Skyfield")
    print("=" * 70)

    max_err = 0.0
    all_pass = True

    for ecl_name, jd in eclipses:
        print(f"\n  {ecl_name} (JD {jd:.2f}):")
        for fname, func in funcs:
            try:
                val_leb = _run_with_mode("leb", func, jd)
                val_sf = _run_with_mode("skyfield", func, jd)
                diff = abs(val_leb - val_sf)
                max_err = max(max_err, diff)
                ok = "OK" if diff < 1e-6 else "FAIL"
                if diff >= 1e-6:
                    all_pass = False
                print(f"    {fname:10s}: LEB={val_leb:+12.6f}  SF={val_sf:+12.6f}  diff={diff:.2e}  [{ok}]")
            except Exception as e:
                print(f"    {fname:10s}: ERROR {e}")
                all_pass = False

    print(f"\n  Max error: {max_err:.2e}")
    print("  Target: < 1e-6")
    print(f"  Result: {'PASS' if all_pass else 'FAIL'}")
    return all_pass


def compare_eclipse_timing():
    """Compare eclipse search timing LEB vs Skyfield."""
    from libephemeris.eclipse import lun_eclipse_when, sol_eclipse_when_glob

    print("\n" + "=" * 70)
    print("ECLIPSE TIMING: LEB vs Skyfield")
    print("=" * 70)

    search_dates = [
        ("Solar from 2024-01", swe_julday(2024, 1, 1, 0.0)),
        ("Solar from 2025-01", swe_julday(2025, 1, 1, 0.0)),
        ("Solar from 2026-01", swe_julday(2026, 1, 1, 0.0)),
    ]

    all_pass = True

    for name, jd in search_dates:
        print(f"\n  {name}:")

        t0 = time.perf_counter()
        type_leb, tret_leb = _run_with_mode("leb", sol_eclipse_when_glob, jd)
        dt_leb = time.perf_counter() - t0

        t0 = time.perf_counter()
        type_sf, tret_sf = _run_with_mode("skyfield", sol_eclipse_when_glob, jd)
        dt_sf = time.perf_counter() - t0

        if type_leb == type_sf:
            jd_diff_sec = abs(tret_leb[0] - tret_sf[0]) * 86400.0
            ok = "OK" if jd_diff_sec < 1.0 else "FAIL"
            if jd_diff_sec >= 1.0:
                all_pass = False
            print(f"    sol_glob: type={type_leb}  max_diff={jd_diff_sec:.3f}s  [{ok}]")
            print(f"              LEB: {dt_leb:.3f}s  Skyfield: {dt_sf:.3f}s  speedup: {dt_sf/dt_leb:.1f}x")
        else:
            print(f"    sol_glob: TYPE MISMATCH LEB={type_leb} SF={type_sf}")
            all_pass = False

    print("\n  Lunar timing:")
    for name, jd in search_dates[:2]:
        type_leb, tret_leb = _run_with_mode("leb", lun_eclipse_when, jd)
        type_sf, tret_sf = _run_with_mode("skyfield", lun_eclipse_when, jd)
        if type_leb == type_sf and tret_leb[0] > 0 and tret_sf[0] > 0:
            jd_diff_sec = abs(tret_leb[0] - tret_sf[0]) * 86400.0
            ok = "OK" if jd_diff_sec < 1.0 else "FAIL"
            if jd_diff_sec >= 1.0:
                all_pass = False
            print(f"    lun_when ({name}): max_diff={jd_diff_sec:.3f}s  [{ok}]")
        else:
            print(f"    lun_when ({name}): MISMATCH type LEB={type_leb} SF={type_sf}")
            all_pass = False

    print(f"\n  Result: {'PASS' if all_pass else 'FAIL'}")
    return all_pass


def compare_rise_trans():
    """Compare rise/trans LEB vs Skyfield."""
    from libephemeris.eclipse import rise_trans

    print("\n" + "=" * 70)
    print("RISE/TRANS: LEB vs Skyfield")
    print("=" * 70)

    locations = [
        ("Rome", [12.5, 41.9, 0]),
        ("New York", [-74.0, 40.7, 10]),
        ("Tokyo", [139.7, 35.7, 40]),
        ("Cape Town", [18.4, -33.9, 0]),
        ("Reykjavik", [-21.9, 64.1, 0]),
    ]

    jd = swe_julday(2024, 6, 15, 0.0)
    all_pass = True

    for loc_name, geopos in locations:
        try:
            _, tret_leb = _run_with_mode("leb", rise_trans, jd, SE_SUN, SE_CALC_RISE, geopos)
            _, tret_sf = _run_with_mode("skyfield", rise_trans, jd, SE_SUN, SE_CALC_RISE, geopos)

            if tret_leb[0] > 0 and tret_sf[0] > 0:
                diff_sec = abs(tret_leb[0] - tret_sf[0]) * 86400.0
                ok = "OK" if diff_sec < 2.0 else "FAIL"
                if diff_sec >= 2.0:
                    all_pass = False
                print(f"  {loc_name:15s} sunrise: diff={diff_sec:.3f}s  [{ok}]")
            else:
                print(f"  {loc_name:15s} sunrise: no result")
                all_pass = False
        except Exception as e:
            print(f"  {loc_name:15s} sunrise: ERROR {e}")
            all_pass = False

    print(f"\n  Result: {'PASS' if all_pass else 'FAIL'}")
    return all_pass


def main():
    print("Eclipse LEB vs Skyfield Precision Test")
    print(f"Mode: {libephemeris.get_calc_mode()}")
    print()

    results = []
    results.append(("Besselian elements", compare_besselian_elements()))
    results.append(("Eclipse timing", compare_eclipse_timing()))
    results.append(("Rise/trans", compare_rise_trans()))

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    all_ok = True
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"  {name:25s}: {status}")
        if not passed:
            all_ok = False

    print(f"\n  Overall: {'ALL PASS' if all_ok else 'SOME FAILURES'}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
