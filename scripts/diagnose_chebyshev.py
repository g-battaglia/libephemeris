#!/usr/bin/env python3
"""Diagnose Chebyshev fitting error for Saturn segment around JD 2501964.8.

Directly evaluates the generator pipeline at Chebyshev nodes for one segment,
fits, and compares — pinpointing whether the fitting itself is the problem
or the pipeline output.
"""

from __future__ import annotations
import sys
import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval

sys.path.insert(0, ".")

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED
from libephemeris.state import get_planets, get_timescale

# Saturn params: interval=2d, degree=15
BODY_ID = 6
INTERVAL = 2.0
DEGREE = 15
JD_START = 2396758.5  # base tier start
TARGET_JD = 2501964.8


def chebyshev_nodes(n):
    return np.cos(np.pi * (np.arange(n) + 0.5) / n)


def eval_pipeline_at_jds(jds_array):
    """Evaluate the generator's _apply_geo_ecliptic_pipeline at given JDs."""
    from scripts.generate_leb import (
        _eval_body_icrs_vectorized,
        _apply_geo_ecliptic_pipeline,
        _get_spk_jd_range,
    )

    planets = get_planets()
    ts = get_timescale()
    C_AU_DAY = 173.1446326846693

    all_jds = np.array(jds_array)
    earth_pos = _eval_body_icrs_vectorized("earth", all_jds, planets, ts)
    spk_min, spk_max = _get_spk_jd_range(planets)
    clamped = np.clip(all_jds, spk_min + 1.0, spk_max - 1.0)
    t_c = ts.tt_jd(clamped)
    earth_vel = np.asarray(planets["earth"].at(t_c).velocity.au_per_d).T

    target_pos = _eval_body_icrs_vectorized("saturn", all_jds, planets, ts)
    geo = target_pos - earth_pos

    for _ in range(3):
        dist = np.sqrt(np.sum(geo**2, axis=1))
        lt = dist / C_AU_DAY
        retarded = all_jds - lt
        ret_pos = _eval_body_icrs_vectorized("saturn", retarded, planets, ts)
        geo = ret_pos - earth_pos

    result = _apply_geo_ecliptic_pipeline(
        geo, earth_vel, all_jds, ts, earth_pos, lt, planets
    )
    return result  # (N, 3): lon, lat, dist


def eval_swe_calc_at_jds(jds_array):
    """Evaluate calc at given JDs."""
    ephem.set_calc_mode("skyfield")
    results = np.zeros((len(jds_array), 3))
    for i, jd in enumerate(jds_array):
        ref, _ = ephem.calc(float(jd), BODY_ID, FLG_SPEED)
        results[i] = [ref[0], ref[1], ref[2]]
    return results


def ang_diff(a, b):
    d = abs(a - b)
    if d > 180:
        d = 360 - d
    return d


def main():
    # Find the segment containing TARGET_JD
    seg_idx = int((TARGET_JD - JD_START) / INTERVAL)
    seg_start = JD_START + seg_idx * INTERVAL
    seg_end = seg_start + INTERVAL
    mid = 0.5 * (seg_start + seg_end)
    half = 0.5 * (seg_end - seg_start)

    print(f"Saturn segment {seg_idx}: [{seg_start:.6f}, {seg_end:.6f}]")
    print(f"  contains JD {TARGET_JD}")
    print(f"  interval={INTERVAL}d  degree={DEGREE}")

    nodes_01 = chebyshev_nodes(DEGREE + 1)
    node_jds = half * nodes_01 + mid

    print(f"\n--- Evaluating at {DEGREE + 1} Chebyshev nodes ---")

    # Method A: Generator pipeline (vectorized)
    pipeline_vals = eval_pipeline_at_jds(node_jds)

    # Method B: calc (scalar, reference)
    swe_vals = eval_swe_calc_at_jds(node_jds)

    print("\nNode-by-node comparison (pipeline vs calc):")
    for i in range(DEGREE + 1):
        lon_err = ang_diff(pipeline_vals[i, 0], swe_vals[i, 0]) * 3600
        abs(pipeline_vals[i, 1] - swe_vals[i, 1]) * 3600
        print(
            f"  node {i:2d}: JD={node_jds[i]:.6f}  "
            f"pipe_lon={pipeline_vals[i, 0]:.10f}  swe_lon={swe_vals[i, 0]:.10f}  "
            f'diff={lon_err:.6f}"'
        )

    # Fit Chebyshev to pipeline values
    pipe_fit = pipeline_vals.copy()
    pipe_fit[:, 0] = np.degrees(np.unwrap(np.radians(pipe_fit[:, 0])))
    coeffs_pipe = np.zeros((3, DEGREE + 1))
    for c in range(3):
        coeffs_pipe[c] = chebfit(nodes_01, pipe_fit[:, c], DEGREE)

    # Fit Chebyshev to calc values
    swe_fit = swe_vals.copy()
    swe_fit[:, 0] = np.degrees(np.unwrap(np.radians(swe_fit[:, 0])))
    coeffs_swe = np.zeros((3, DEGREE + 1))
    for c in range(3):
        coeffs_swe[c] = chebfit(nodes_01, swe_fit[:, c], DEGREE)

    # Evaluate both fits at test points
    n_test = 50
    print(f"\n--- Verification at {n_test} intermediate points ---")
    print(
        f"{'JD':>16s}  {'pipe_fit_err':>14s}  {'swe_fit_err':>14s}  {'pipe_vs_swe':>14s}"
    )

    max_pipe_err = 0.0
    max_swe_err = 0.0
    max_cross_err = 0.0

    for k in range(n_test):
        frac = (k + 0.5) / n_test
        jd_test = seg_start + frac * (seg_end - seg_start)
        tau = (jd_test - mid) / half

        # Reference: calc at this point
        ref, _ = ephem.calc(float(jd_test), BODY_ID, FLG_SPEED)

        # Pipeline reference at this point
        pipe_ref = eval_pipeline_at_jds([jd_test])[0]

        # Evaluate fitted Chebyshev
        pipe_fitted_lon = float(chebval(tau, coeffs_pipe[0])) % 360.0
        swe_fitted_lon = float(chebval(tau, coeffs_swe[0])) % 360.0

        err_pipe = ang_diff(pipe_fitted_lon, pipe_ref[0]) * 3600
        err_swe = ang_diff(swe_fitted_lon, ref[0]) * 3600
        err_cross = ang_diff(pipe_fitted_lon, ref[0]) * 3600

        if err_pipe > max_pipe_err:
            max_pipe_err = err_pipe
        if err_swe > max_swe_err:
            max_swe_err = err_swe
        if err_cross > max_cross_err:
            max_cross_err = err_cross

        if k % 10 == 0:
            print(
                f'  {jd_test:.6f}  {err_pipe:12.6f}"  {err_swe:12.6f}"  {err_cross:12.6f}"'
            )

    print("\nMax errors:")
    print(f'  Pipeline fit vs pipeline ref:  {max_pipe_err:.6f}"')
    print(f'  calc fit vs calc ref:  {max_swe_err:.6f}"')
    print(f'  Pipeline fit vs calc ref:  {max_cross_err:.6f}"')

    # Also test: what does the LEB file give?
    from libephemeris.leb_reader import LEBReader

    try:
        reader = LEBReader("data/leb/ephemeris_base_planets.leb")
        dt = reader.delta_t(TARGET_JD)
        jd_tt = TARGET_JD + dt
        (leb_lon, leb_lat, leb_dist), _ = reader.eval_body(BODY_ID, jd_tt)
        ref, _ = ephem.calc(float(jd_tt), BODY_ID, FLG_SPEED)
        print(f"\nLEB file read at JD_TT={jd_tt:.12f}:")
        print(
            f"  LEB lon: {leb_lon:.10f}°  calc lon: {ref[0]:.10f}°  "
            f'error: {ang_diff(leb_lon, ref[0]) * 3600:.6f}"'
        )
        reader.close()
    except Exception as e:
        print(f"\nCouldn't read LEB file: {e}")


if __name__ == "__main__":
    main()
