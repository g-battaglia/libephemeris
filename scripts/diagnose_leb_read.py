#!/usr/bin/env python3
"""Diagnose LEB file read vs calc at specific JD.

Reads the actual .leb Chebyshev data for Saturn at the worst-case JD
and compares with calc.  This tells us whether the stored Chebyshev
data matches what calc produces, or whether the error is introduced
somewhere else.
"""

from __future__ import annotations

import sys

sys.path.insert(0, ".")

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED
from libephemeris.leb_reader import LEBReader


def check_leb_vs_swe_calc(leb_path: str, body_id: int, body_name: str, jd_ut: float):
    """Read LEB and compare with calc at the same JD."""
    print(f"\n{'=' * 72}")
    print(f"  {body_name} (body {body_id}) at JD_UT = {jd_ut:.6f}")
    print(f"{'=' * 72}")

    reader = LEBReader(leb_path)

    # 1. LEB reader's delta-T conversion
    delta_t = reader.delta_t(jd_ut)
    jd_tt = jd_ut + delta_t
    print(f"  LEB delta_t:  {delta_t:.12f} days = {delta_t * 86400:.6f} seconds")
    print(f"  LEB jd_tt:    {jd_tt:.12f}")

    # 2. Skyfield's delta-T
    from libephemeris.state import get_timescale

    ts = get_timescale()
    t = ts.ut1_jd(jd_ut)
    skyfield_tt = t.tt
    skyfield_delta_t = skyfield_tt - jd_ut
    print(
        f"  Skyfield δT:  {skyfield_delta_t:.12f} days = {skyfield_delta_t * 86400:.6f} seconds"
    )
    print(f"  Skyfield TT:  {skyfield_tt:.12f}")
    print(f"  TT diff:      {(jd_tt - skyfield_tt) * 86400:.6f} seconds")

    # 3. Read LEB at the TT JD
    (lon_leb, lat_leb, dist_leb), (dlon, dlat, ddist) = reader.eval_body(body_id, jd_tt)
    print(f"\n  LEB read at jd_tt={jd_tt:.12f}:")
    print(f"    lon:  {lon_leb:.10f}°")
    print(f"    lat:  {lat_leb:.10f}°")
    print(f"    dist: {dist_leb:.12f} AU")

    # 4. Also read LEB at Skyfield's TT
    (lon_leb2, lat_leb2, dist_leb2), _ = reader.eval_body(body_id, skyfield_tt)
    print(f"\n  LEB read at skyfield_tt={skyfield_tt:.12f}:")
    print(f"    lon:  {lon_leb2:.10f}°")
    print(f"    lat:  {lat_leb2:.10f}°")
    print(f"    dist: {dist_leb2:.12f} AU")

    # 5. calc at both TT values (Skyfield mode, no LEB)
    ephem.set_calc_mode("skyfield")
    ref_ut, _ = ephem.calc_ut(jd_ut, body_id, FLG_SPEED)
    ref_tt, _ = ephem.calc(jd_tt, body_id, FLG_SPEED)
    ref_tt2, _ = ephem.calc(skyfield_tt, body_id, FLG_SPEED)

    print(f"\n  calc_ut(jd_ut={jd_ut:.6f}):")
    print(f"    lon:  {ref_ut[0]:.10f}°")
    print(f"    lat:  {ref_ut[1]:.10f}°")
    print(f"    dist: {ref_ut[2]:.12f} AU")

    print(f"\n  calc(jd_tt={jd_tt:.12f}):")
    print(f"    lon:  {ref_tt[0]:.10f}°")
    print(f"    lat:  {ref_tt[1]:.10f}°")

    print(f"\n  calc(skyfield_tt={skyfield_tt:.12f}):")
    print(f"    lon:  {ref_tt2[0]:.10f}°")
    print(f"    lat:  {ref_tt2[1]:.10f}°")

    # 6. Errors
    def ang_diff(a, b):
        d = abs(a - b)
        if d > 180:
            d = 360 - d
        return d

    # LEB(leb_tt) vs calc_ut (what the compare test does)
    err1_lon = ang_diff(lon_leb, ref_ut[0]) * 3600
    err1_lat = abs(lat_leb - ref_ut[1]) * 3600

    # LEB(leb_tt) vs calc(leb_tt)
    err2_lon = ang_diff(lon_leb, ref_tt[0]) * 3600
    err2_lat = abs(lat_leb - ref_tt[1]) * 3600

    # LEB(skyfield_tt) vs calc_ut
    err3_lon = ang_diff(lon_leb2, ref_ut[0]) * 3600
    err3_lat = abs(lat_leb2 - ref_ut[1]) * 3600

    # LEB(skyfield_tt) vs calc(skyfield_tt)
    err4_lon = ang_diff(lon_leb2, ref_tt2[0]) * 3600
    err4_lat = abs(lat_leb2 - ref_tt2[1]) * 3600

    print("\n  ERRORS (arcsec):")
    print(
        f'    LEB(leb_tt) vs calc_ut:         lon={err1_lon:.6f}"  lat={err1_lat:.6f}"'
    )
    print(
        f'    LEB(leb_tt) vs calc(leb_tt):     lon={err2_lon:.6f}"  lat={err2_lat:.6f}"'
    )
    print(
        f'    LEB(sky_tt) vs calc_ut:          lon={err3_lon:.6f}"  lat={err3_lat:.6f}"'
    )
    print(
        f'    LEB(sky_tt) vs calc(sky_tt):     lon={err4_lon:.6f}"  lat={err4_lat:.6f}"'
    )

    # 7. What the compare test actually does
    print("\n  COMPARE TEST simulation:")
    print("    The test calls calc_ut(jd_ut) in both skyfield and LEB mode")

    # LEB mode: calc_ut → fast_calc_ut → jd_tt = jd_ut + reader.delta_t(jd_ut)
    # → _pipeline_geo_ecliptic → reader.eval_body(ipl, jd_tt)
    # Skyfield mode: calc_ut → t = ts.ut1_jd(jd_ut) → _calc_body(t, ...)
    # The compare test compares these two results

    # So the LEB mode uses reader.delta_t(jd_ut) to get TT
    # And Skyfield uses ts.ut1_jd(jd_ut).tt to get TT internally
    # If these differ, we get a time-shift error

    # Saturn speed in lon is about 0.033 deg/day
    speed_lon = ref_ut[3]  # deg/day
    time_shift_sec = (jd_tt - skyfield_tt) * 86400
    expected_lon_error = abs(speed_lon) * abs(time_shift_sec) / 86400 * 3600
    print(f"    Saturn speed:      {speed_lon:.6f} deg/day")
    print(f"    TT time shift:     {time_shift_sec:.6f} seconds")
    print(f'    Expected lon err:  {expected_lon_error:.6f}" from time shift')

    # 8. Check the actual LEB data stored - what was the generator's pipeline output?
    # The Chebyshev fitting error should be <0.001", so LEB(jd_tt) should be
    # close to what the generator pipeline produced at jd_tt.
    # And from diagnose_pipeline.py, the generator pipeline matches calc(jd_tt) to ~0.0006"
    # So LEB(jd_tt) vs calc(jd_tt) should be < 0.001" + 0.001" ≈ 0.002"
    # But we see err2_lon above...

    reader.close()


def main():
    leb_path = "data/leb/ephemeris_base.leb"

    test_cases = [
        (6, "Saturn", 2501964.8),  # Worst case from compare test
        (6, "Saturn", 2451545.0),  # J2000
        (9, "Pluto", 2501964.8),  # Another failing body
        (8, "Neptune", 2501964.8),  # Another failing body
        (0, "Sun", 2451545.0),  # Control — should be perfect
    ]

    for body_id, body_name, jd_ut in test_cases:
        check_leb_vs_swe_calc(leb_path, body_id, body_name, jd_ut)


if __name__ == "__main__":
    main()
