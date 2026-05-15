#!/usr/bin/env python3
"""Diagnose pipeline discrepancy between generator and calc().

Compares the generator's _apply_geo_ecliptic_pipeline() against
Skyfield's observer.at(t).observe(target).apparent() step by step
at specific TT Julian Days for Saturn, Neptune, and Pluto.

This identifies which step introduces the ~1" discrepancy.
"""

from __future__ import annotations

import math
import sys

import numpy as np

# Ensure libephemeris is importable
sys.path.insert(0, ".")

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED
from libephemeris.state import get_planets, get_timescale


def diagnose_body(body_id: int, body_name: str, jd_tt: float) -> None:
    """Compare generator pipeline vs Skyfield apparent() for one body at one JD."""
    print(f"\n{'=' * 72}")
    print(f"  {body_name} (body {body_id}) at JD_TT = {jd_tt:.6f}")
    print(f"{'=' * 72}")

    planets = get_planets()
    ts = get_timescale()
    t = ts.tt_jd(jd_tt)

    # ========================================================================
    # PATH A: Skyfield's observer.at(t).observe(target).apparent()
    # ========================================================================
    from libephemeris.planets import get_planet_target, _PLANET_MAP

    target_name = _PLANET_MAP[body_id]
    target = get_planet_target(planets, target_name)
    observer = planets["earth"]

    # Step 1: Earth barycentric position
    earth_at_t = observer.at(t)
    earth_pos_sky = earth_at_t.position.au  # (3,)
    earth_vel_sky = earth_at_t.velocity.au_per_d  # (3,)

    # Step 2: Light-time corrected astrometric position
    astrometric = earth_at_t.observe(target)
    astro_pos_sky = (
        astrometric.position.au
    )  # (3,) geocentric ICRS, light-time corrected
    light_time_sky = astrometric.light_time  # scalar, days

    # Step 3: Apparent position (deflection + aberration)
    apparent = astrometric.apparent()
    apparent_pos_sky = apparent.position.au  # (3,) after deflection + aberration

    # Step 4: Ecliptic of date
    from skyfield.framelib import ecliptic_frame

    lat_sky, lon_sky, dist_sky = apparent.frame_latlon(ecliptic_frame)
    lon_sky_deg = lon_sky.degrees
    lat_sky_deg = lat_sky.degrees
    dist_sky_au = dist_sky.au

    print("\n  PATH A: Skyfield apparent()")
    print(
        f"    Earth bary pos:  ({earth_pos_sky[0]:.12f}, {earth_pos_sky[1]:.12f}, {earth_pos_sky[2]:.12f})"
    )
    print(
        f"    Earth bary vel:  ({earth_vel_sky[0]:.12f}, {earth_vel_sky[1]:.12f}, {earth_vel_sky[2]:.12f})"
    )
    print(
        f"    Astrometric pos: ({astro_pos_sky[0]:.12f}, {astro_pos_sky[1]:.12f}, {astro_pos_sky[2]:.12f})"
    )
    print(f"    Light-time:      {light_time_sky:.12f} days")
    print(
        f"    Apparent pos:    ({apparent_pos_sky[0]:.12f}, {apparent_pos_sky[1]:.12f}, {apparent_pos_sky[2]:.12f})"
    )
    print(f"    Ecliptic lon:    {lon_sky_deg:.10f}°")
    print(f"    Ecliptic lat:    {lat_sky_deg:.10f}°")
    print(f"    Distance:        {dist_sky_au:.12f} AU")

    # ========================================================================
    # PATH B: Generator pipeline (_apply_geo_ecliptic_pipeline)
    # ========================================================================
    from scripts.generate_leb import (
        _eval_body_icrs_vectorized,
        _apply_geo_ecliptic_pipeline,
        _get_spk_jd_range,
    )

    all_jds = np.array([jd_tt])

    # Step 1: Earth position + velocity
    earth_pos_gen = _eval_body_icrs_vectorized("earth", all_jds, planets, ts)  # (1, 3)

    spk_min, spk_max = _get_spk_jd_range(planets)
    clamped_jds = np.clip(all_jds, spk_min + 1.0, spk_max - 1.0)
    t_clamped = ts.tt_jd(clamped_jds)
    earth_vel_gen = np.asarray(
        planets["earth"].at(t_clamped).velocity.au_per_d
    ).T  # (1, 3)

    # Step 2: Target position
    target_pos_gen = _eval_body_icrs_vectorized(
        target_name, all_jds, planets, ts
    )  # (1, 3)

    # Step 3: Geometric geocentric
    geo_gen = target_pos_gen - earth_pos_gen  # (1, 3)

    # Step 4: Light-time correction (3 iterations)
    C_AU_DAY = 173.1446326846693
    for lt_iter in range(3):
        dist_gen = np.sqrt(np.sum(geo_gen**2, axis=1))
        lt_gen = dist_gen / C_AU_DAY
        retarded_jds = all_jds - lt_gen
        retarded_pos = _eval_body_icrs_vectorized(
            target_name, retarded_jds, planets, ts
        )
        geo_gen = retarded_pos - earth_pos_gen

    print("\n  PATH B: Generator pipeline")
    print(
        f"    Earth bary pos:  ({earth_pos_gen[0, 0]:.12f}, {earth_pos_gen[0, 1]:.12f}, {earth_pos_gen[0, 2]:.12f})"
    )
    print(
        f"    Earth bary vel:  ({earth_vel_gen[0, 0]:.12f}, {earth_vel_gen[0, 1]:.12f}, {earth_vel_gen[0, 2]:.12f})"
    )
    print(
        f"    LT-corrected geo:({geo_gen[0, 0]:.12f}, {geo_gen[0, 1]:.12f}, {geo_gen[0, 2]:.12f})"
    )
    print(f"    Light-time:      {lt_gen[0]:.12f} days")

    # Compare astrometric positions
    diff_astro = np.sqrt(np.sum((geo_gen[0] - astro_pos_sky) ** 2))
    print("\n  COMPARISON: astrometric positions (before deflection/aberration)")
    print(
        f"    Skyfield astro:  ({astro_pos_sky[0]:.12f}, {astro_pos_sky[1]:.12f}, {astro_pos_sky[2]:.12f})"
    )
    print(
        f"    Generator astro: ({geo_gen[0, 0]:.12f}, {geo_gen[0, 1]:.12f}, {geo_gen[0, 2]:.12f})"
    )
    print(
        f"    |diff|:          {diff_astro:.2e} AU = {diff_astro * 206265 / dist_sky_au:.6f} arcsec"
    )

    # Step 5-9: Full pipeline (deflection + aberration + PNM + ecliptic)
    result_gen = _apply_geo_ecliptic_pipeline(
        geo_gen,
        earth_vel_gen,
        all_jds,
        ts,
        earth_pos_gen,
        lt_gen,
        planets,
    )
    lon_gen = result_gen[0, 0]
    lat_gen = result_gen[0, 1]
    dist_gen_final = result_gen[0, 2]

    print("\n  PATH B final:")
    print(f"    Ecliptic lon:    {lon_gen:.10f}°")
    print(f"    Ecliptic lat:    {lat_gen:.10f}°")
    print(f"    Distance:        {dist_gen_final:.12f} AU")

    # ========================================================================
    # PATH C: calc() (what compare tests use as reference)
    # ========================================================================
    ephem.set_calc_mode("skyfield")
    ref, _ = ephem.calc(jd_tt, body_id, FLG_SPEED)
    lon_ref = ref[0]
    lat_ref = ref[1]
    dist_ref = ref[2]

    print(f"\n  PATH C: calc(jd_tt, {body_id}, FLG_SPEED)")
    print(f"    Ecliptic lon:    {lon_ref:.10f}°")
    print(f"    Ecliptic lat:    {lat_ref:.10f}°")
    print(f"    Distance:        {dist_ref:.12f} AU")

    # ========================================================================
    # SUMMARY
    # ========================================================================
    def ang_diff(a, b):
        d = abs(a - b)
        if d > 180:
            d = 360 - d
        return d

    err_sky_vs_ref_lon = ang_diff(lon_sky_deg, lon_ref) * 3600
    err_sky_vs_ref_lat = abs(lat_sky_deg - lat_ref) * 3600
    err_gen_vs_ref_lon = ang_diff(lon_gen, lon_ref) * 3600
    err_gen_vs_ref_lat = abs(lat_gen - lat_ref) * 3600
    err_gen_vs_sky_lon = ang_diff(lon_gen, lon_sky_deg) * 3600
    err_gen_vs_sky_lat = abs(lat_gen - lat_sky_deg) * 3600

    print("\n  ERROR SUMMARY (arcsec):")
    print(
        f'    Skyfield.apparent vs calc:  lon={err_sky_vs_ref_lon:.6f}"  lat={err_sky_vs_ref_lat:.6f}"'
    )
    print(
        f'    Generator vs calc:          lon={err_gen_vs_ref_lon:.6f}"  lat={err_gen_vs_ref_lat:.6f}"'
    )
    print(
        f'    Generator vs Skyfield.apparent:  lon={err_gen_vs_sky_lon:.6f}"  lat={err_gen_vs_sky_lat:.6f}"'
    )

    # ========================================================================
    # DECOMPOSE: Manually apply deflection + aberration to generator's
    # astrometric position using Skyfield's code, and compare
    # ========================================================================
    from scripts.generate_leb import _apply_gravitational_deflection

    # Generator's deflection
    geo_defl_gen = _apply_gravitational_deflection(
        geo_gen,
        earth_pos_gen,
        all_jds,
        lt_gen,
        planets,
        ts,
    )

    # Skyfield's deflection: apparent_pos has deflection+aberration, but
    # we can get deflection-only by examining internals
    # Instead, let's just compare our deflected position
    defl_diff = np.sqrt(np.sum((geo_defl_gen[0] - geo_gen[0]) ** 2))
    print(
        f"\n  DEFLECTION magnitude (generator): {defl_diff:.2e} AU = {defl_diff * 206265 / dist_sky_au:.6f} arcsec"
    )

    # Now apply generator's aberration manually
    p1mag = lt_gen * C_AU_DAY
    vemag = np.sqrt(np.sum(earth_vel_gen**2, axis=1))
    beta = vemag / C_AU_DAY
    dot = np.sum(geo_defl_gen * earth_vel_gen, axis=1)
    safe_denom = np.where(p1mag * vemag > 0, p1mag * vemag, 1.0)
    cosd = dot / safe_denom
    gammai = np.sqrt(1.0 - beta * beta)
    p = beta * cosd
    q = (1.0 + p / (1.0 + gammai)) * lt_gen
    r = 1.0 + p
    safe_r = np.where(np.abs(r) > 1e-30, r, 1.0)
    geo_aber_gen = (
        gammai[:, np.newaxis] * geo_defl_gen + q[:, np.newaxis] * earth_vel_gen
    ) / safe_r[:, np.newaxis]

    # Compare deflected+aberrated position (ICRS) between generator and Skyfield
    diff_aber = np.sqrt(np.sum((geo_aber_gen[0] - apparent_pos_sky) ** 2))
    print(
        f"  DEFLECTED+ABERRATED diff:         {diff_aber:.2e} AU = {diff_aber * 206265 / dist_sky_au:.6f} arcsec"
    )

    # Now check if the difference is in the ICRS apparent position or in the frame rotation
    # Apply Skyfield's frame rotation to the generator's ICRS apparent position
    from skyfield.functions import mxv

    ecl_mat = ecliptic_frame.rotation_at(t)
    gen_aber_ecl = mxv(ecl_mat, geo_aber_gen[0])
    # Convert to spherical
    from skyfield.functions import to_spherical

    d_gen2, lat_gen2, lon_gen2 = to_spherical(gen_aber_ecl)
    lon_gen2_deg = math.degrees(lon_gen2)
    lat_gen2_deg = math.degrees(lat_gen2)
    if lon_gen2_deg < 0:
        lon_gen2_deg += 360.0

    err_gen2_vs_sky_lon = ang_diff(lon_gen2_deg, lon_sky_deg) * 3600
    err_gen2_vs_ref_lon = ang_diff(lon_gen2_deg, lon_ref) * 3600
    err_gen2_vs_ref_lat = abs(lat_gen2_deg - lat_ref) * 3600

    print("\n  GEN ICRS→Skyfield ecliptic frame:")
    print(f"    lon={lon_gen2_deg:.10f}°  lat={lat_gen2_deg:.10f}°")
    print(f'    vs Skyfield.apparent: lon={err_gen2_vs_sky_lon:.6f}"')
    print(
        f'    vs calc:          lon={err_gen2_vs_ref_lon:.6f}"  lat={err_gen2_vs_ref_lat:.6f}"'
    )

    # Also check: generator's own PNM rotation vs Skyfield's ecliptic_frame
    M_arr = np.array(t.M)  # (3, 3, N) or (3, 3) for scalar
    if M_arr.ndim == 2:
        pn_mat = M_arr
    else:
        pn_mat = M_arr  # scalar time → (3, 3)

    mean_obl = float(t._mean_obliquity_radians)
    _dpsi, deps = t._nutation_angles_radians
    eps_true = mean_obl + float(deps)

    cos_eps = math.cos(eps_true)
    sin_eps = math.sin(eps_true)

    # Apply PNM
    geo_eq_gen = pn_mat @ geo_aber_gen[0]
    # Equatorial → ecliptic
    ecl_x = geo_eq_gen[0]
    ecl_y = geo_eq_gen[1] * cos_eps + geo_eq_gen[2] * sin_eps
    ecl_z = -geo_eq_gen[1] * sin_eps + geo_eq_gen[2] * cos_eps

    lon_gen3 = math.degrees(math.atan2(ecl_y, ecl_x)) % 360.0
    lat_gen3 = math.degrees(
        math.asin(ecl_z / math.sqrt(ecl_x**2 + ecl_y**2 + ecl_z**2))
    )

    err_gen3_vs_sky_lon = ang_diff(lon_gen3, lon_sky_deg) * 3600
    err_gen3_vs_ref_lon = ang_diff(lon_gen3, lon_ref) * 3600

    print("\n  GEN ICRS→manual PNM+obliquity:")
    print(f"    lon={lon_gen3:.10f}°  lat={lat_gen3:.10f}°")
    print(f'    vs Skyfield.apparent: lon={err_gen3_vs_sky_lon:.6f}"')
    print(f'    vs calc:          lon={err_gen3_vs_ref_lon:.6f}"')
    print(
        f'    vs Skyfield frame:    lon={ang_diff(lon_gen3, lon_gen2_deg) * 3600:.6f}"'
    )


def main():
    # Test JDs from where errors were observed
    # Saturn worst: JD 2501964.8 (~year 2138)
    # Also test a nearby year for comparison
    test_cases = [
        (6, "Saturn", 2501964.8),
        (6, "Saturn", 2451545.0),  # J2000.0 — should be good
        (6, "Saturn", 2488069.0),  # ~year 2100
        (9, "Pluto", 2501964.8),
        (8, "Neptune", 2501964.8),
        (0, "Sun", 2501964.8),  # Sun should be well-behaved
    ]

    for body_id, body_name, jd_tt in test_cases:
        diagnose_body(body_id, body_name, jd_tt)

    print("\nDone.")


if __name__ == "__main__":
    main()
