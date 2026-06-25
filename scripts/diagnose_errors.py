#!/usr/bin/env python3
"""Decompose LEB errors into: Chebyshev fitting, COB offset, pipeline differences.

For each failing outer planet, measures:
1. Raw Chebyshev error: LEB eval vs Skyfield geometric ICRS (same COB source)
2. COB mismatch: Skyfield geometric with SPK centers vs analytical COB
3. Pipeline error: LEB pipeline using perfect positions vs Skyfield apparent
4. Total end-to-end: LEB pipeline vs calc
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, ".")

os.environ.pop("LIBEPHEMERIS_LEB", None)

import libephemeris as ephem
from libephemeris.constants import (
    EARTH,
    FLG_SPEED,
)
from libephemeris.state import get_planets, get_timescale, get_planet_center_segment
from libephemeris.planets import _PLANET_CENTER_NAIF_IDS
from libephemeris.moon_theories import get_cob_offset
from libephemeris.leb_reader import LEBReader
from libephemeris.fast_calc import fast_calc_tt


BODY_NAMES = {
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "Pluto",
    0: "Sun",
    1: "Moon",
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    14: "Earth",
}

_PLANET_MAP = {
    0: "sun",
    1: "moon",
    2: "mercury",
    3: "venus",
    4: "mars",
    5: "jupiter",
    6: "saturn",
    7: "uranus",
    8: "neptune",
    9: "pluto",
    14: "earth",
}

_BARY_MAP = {
    "jupiter": "jupiter barycenter",
    "saturn": "saturn barycenter",
    "uranus": "uranus barycenter",
    "neptune": "neptune barycenter",
    "pluto": "pluto barycenter",
}


def angular_sep_arcsec(lon1, lat1, lon2, lat2):
    """Great-circle angular separation in arcseconds."""
    lon1r, lat1r = math.radians(lon1), math.radians(lat1)
    lon2r, lat2r = math.radians(lon2), math.radians(lat2)
    a = math.cos(lat2r) * math.sin(lon2r - lon1r)
    b = math.cos(lat1r) * math.sin(lat2r) - math.sin(lat1r) * math.cos(
        lat2r
    ) * math.cos(lon2r - lon1r)
    c = math.sin(lat1r) * math.sin(lat2r) + math.cos(lat1r) * math.cos(
        lat2r
    ) * math.cos(lon2r - lon1r)
    sep = math.atan2(math.sqrt(a * a + b * b), c)
    return abs(sep) * 3600.0 * 180.0 / math.pi


def vec_dist_au(a, b):
    """Euclidean distance between two 3-vectors in AU."""
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def au_to_arcsec_at_dist(au_error, dist_au):
    """Convert positional AU error to approximate angular arcsec at given distance."""
    if dist_au == 0:
        return 0.0
    return (au_error / dist_au) * 206264.806  # radians to arcsec


def get_skyfield_icrs_bary(target_name, jd_tt, planets, ts):
    """Get ICRS barycentric position from Skyfield, handling outer planets."""
    # Try direct
    try:
        target = planets[target_name]
        t = ts.tt_jd(jd_tt)
        pos = target.at(t).position.au
        return (float(pos[0]), float(pos[1]), float(pos[2]))
    except KeyError:
        pass

    # Outer planet: barycenter + center offset
    bary_name = _BARY_MAP.get(target_name)
    if bary_name is None:
        raise ValueError(f"No fallback for {target_name}")

    barycenter = planets[bary_name]
    t = ts.tt_jd(jd_tt)
    bary_pos = barycenter.at(t).position.au

    # Try SPK center
    if target_name in _PLANET_CENTER_NAIF_IDS:
        naif_id = _PLANET_CENTER_NAIF_IDS[target_name]
        seg = get_planet_center_segment(naif_id)
        if seg is not None:
            try:
                offset_pos = seg.at(t).position.au
                return (
                    float(bary_pos[0] + offset_pos[0]),
                    float(bary_pos[1] + offset_pos[1]),
                    float(bary_pos[2] + offset_pos[2]),
                )
            except Exception:
                pass

    # Analytical COB
    offset = get_cob_offset(bary_name, t)
    return (
        float(bary_pos[0] + offset[0]),
        float(bary_pos[1] + offset[1]),
        float(bary_pos[2] + offset[2]),
    )


def get_skyfield_bary_with_spk(target_name, jd_tt, planets, ts):
    """Get position using SPK center offset specifically."""
    bary_name = _BARY_MAP.get(target_name)
    if bary_name is None:
        return None

    barycenter = planets[bary_name]
    t = ts.tt_jd(jd_tt)
    bary_pos = barycenter.at(t).position.au

    if target_name in _PLANET_CENTER_NAIF_IDS:
        naif_id = _PLANET_CENTER_NAIF_IDS[target_name]
        seg = get_planet_center_segment(naif_id)
        if seg is not None:
            try:
                offset_pos = seg.at(t).position.au
                return (
                    float(bary_pos[0] + offset_pos[0]),
                    float(bary_pos[1] + offset_pos[1]),
                    float(bary_pos[2] + offset_pos[2]),
                )
            except Exception:
                return None
    return None


def get_skyfield_bary_with_cob(target_name, jd_tt, planets, ts):
    """Get position using analytical COB correction specifically."""
    bary_name = _BARY_MAP.get(target_name)
    if bary_name is None:
        return None

    barycenter = planets[bary_name]
    t = ts.tt_jd(jd_tt)
    bary_pos = barycenter.at(t).position.au
    offset = get_cob_offset(bary_name, t)
    return (
        float(bary_pos[0] + offset[0]),
        float(bary_pos[1] + offset[1]),
        float(bary_pos[2] + offset[2]),
    )


def diagnose_body(ipl, leb_path, jd_samples, planets, ts):
    """Diagnose error sources for a body at multiple JDs."""
    target_name = _PLANET_MAP[ipl]
    reader = LEBReader(leb_path)

    results = {
        "chebyshev_err_au": [],
        "chebyshev_err_arcsec": [],
        "cob_mismatch_au": [],
        "cob_mismatch_arcsec": [],
        "total_err_arcsec": [],
    }

    for jd_tt in jd_samples:
        jd_tt = float(jd_tt)

        # 1. LEB raw Chebyshev position
        try:
            leb_pos, _ = reader.eval_body(ipl, jd_tt)
        except Exception:
            continue

        # 2. Skyfield geometric ICRS bary (same COB source as generator would use)
        try:
            sky_pos = get_skyfield_icrs_bary(target_name, jd_tt, planets, ts)
        except Exception:
            continue

        cheby_err_au = vec_dist_au(leb_pos, sky_pos)

        # Approximate distance from Earth for angular conversion
        try:
            earth_pos, _ = reader.eval_body(EARTH, jd_tt)
            geo_dist = vec_dist_au(leb_pos, earth_pos)
        except Exception:
            geo_dist = 5.0  # fallback

        cheby_err_arcsec = au_to_arcsec_at_dist(cheby_err_au, geo_dist)
        results["chebyshev_err_au"].append(cheby_err_au)
        results["chebyshev_err_arcsec"].append(cheby_err_arcsec)

        # 3. COB mismatch: SPK center vs analytical COB
        spk_pos = get_skyfield_bary_with_spk(target_name, jd_tt, planets, ts)
        cob_pos = get_skyfield_bary_with_cob(target_name, jd_tt, planets, ts)
        if spk_pos and cob_pos:
            cob_err_au = vec_dist_au(spk_pos, cob_pos)
            cob_err_arcsec = au_to_arcsec_at_dist(cob_err_au, geo_dist)
            results["cob_mismatch_au"].append(cob_err_au)
            results["cob_mismatch_arcsec"].append(cob_err_arcsec)

        # 4. Total end-to-end: LEB fast_calc vs calc
        try:
            leb_result, _ = fast_calc_tt(reader, jd_tt, ipl, FLG_SPEED)
            ref_result, _ = ephem.calc(jd_tt, ipl, FLG_SPEED)
            total_err = angular_sep_arcsec(
                ref_result[0],
                ref_result[1],
                leb_result[0],
                leb_result[1],
            )
            results["total_err_arcsec"].append(total_err)
        except Exception:
            pass

    return results


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--leb", required=True)
    parser.add_argument("--samples", type=int, default=500)
    args = parser.parse_args()

    reader = LEBReader(args.leb)
    jd_start = reader._header.jd_start
    jd_end = reader._header.jd_end

    planets = get_planets()
    ts = get_timescale()

    # Check planet_centers availability
    print("=== Planet Centers Availability ===")
    for name, naif_id in _PLANET_CENTER_NAIF_IDS.items():
        seg = get_planet_center_segment(naif_id)
        status = "AVAILABLE" if seg else "NOT AVAILABLE"
        print(f"  {name:10s} (NAIF {naif_id}): {status}")
    print()

    bodies = [5, 6, 7, 8, 9]  # Outer planets that fail
    jd_samples = np.linspace(jd_start + 1, jd_end - 1, args.samples)

    print(f"Date range: JD {jd_start:.1f} - {jd_end:.1f}")
    print(f"Samples: {args.samples}")
    print()

    hdr = (
        f"{'Body':>10s}  {'ChebyMax':>12s}  {'ChebyMean':>12s}  "
        f"{'COBMax':>12s}  {'COBMean':>12s}  "
        f"{'TotalMax':>12s}  {'TotalMean':>12s}"
    )
    print(hdr)
    print("-" * 90)

    for ipl in bodies:
        name = BODY_NAMES[ipl]
        print(f"  Diagnosing {name}...", end=" ", flush=True)

        results = diagnose_body(ipl, args.leb, jd_samples, planets, ts)

        def stat(arr):
            if not arr:
                return 0.0, 0.0
            return float(np.max(arr)), float(np.mean(arr))

        cheby_max, cheby_mean = stat(results["chebyshev_err_arcsec"])
        cob_max, cob_mean = stat(results["cob_mismatch_arcsec"])
        total_max, total_mean = stat(results["total_err_arcsec"])

        print(
            f"\r{name:>10s}  {cheby_max:>12.6f}  {cheby_mean:>12.6f}  "
            f"{cob_max:>12.6f}  {cob_mean:>12.6f}  "
            f"{total_max:>12.6f}  {total_mean:>12.6f}"
        )

        # Detailed breakdown
        if results["chebyshev_err_au"]:
            print(
                f'           Chebyshev fit:  max={cheby_max:.6f}"  mean={cheby_mean:.6f}"  '
                f"max_au={max(results['chebyshev_err_au']):.2e}"
            )
        if results["cob_mismatch_au"]:
            print(
                f'           COB mismatch:   max={cob_max:.6f}"  mean={cob_mean:.6f}"  '
                f"max_au={max(results['cob_mismatch_au']):.2e}"
            )
        if results["total_err_arcsec"]:
            print(
                f'           Total e2e:      max={total_max:.6f}"  mean={total_mean:.6f}"'
            )
        print()


if __name__ == "__main__":
    main()
