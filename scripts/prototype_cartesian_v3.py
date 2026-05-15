#!/usr/bin/env python3
"""Prototype: test Chebyshev fitting of geocentric ecliptic CARTESIAN vs SPHERICAL.

Hypothesis: ecliptic Cartesian (x, y, z) is smooth even during retrogrades,
while spherical (lon, lat, dist) has cusps from atan2. If true, we can:
1. Store Chebyshev fits of ecliptic Cartesian (x, y, z) — smooth, fits well
2. At runtime: read (x, y, z) → convert to (lon, lat, dist) via atan2/asin — exact

This script tests both representations at Saturn's worst-case retrograde segment.
"""

from __future__ import annotations

import math
import sys

import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval

sys.path.insert(0, ".")

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED
from libephemeris.state import get_planets, get_timescale


def chebyshev_nodes(n):
    return np.cos(np.pi * (np.arange(n) + 0.5) / n)


def eval_apparent_ecliptic_cartesian(body_id, jds):
    """Evaluate Skyfield apparent() → ecliptic-of-date Cartesian at given TT JDs.

    Returns (N, 3) array of (ecl_x, ecl_y, ecl_z) in AU.
    """
    from libephemeris.planets import get_planet_target, _PLANET_MAP
    from skyfield.framelib import ecliptic_frame
    from skyfield.functions import mxv

    planets = get_planets()
    ts = get_timescale()
    target_name = _PLANET_MAP[body_id]
    target = get_planet_target(planets, target_name)
    observer = planets["earth"]

    results = np.zeros((len(jds), 3))
    for i, jd in enumerate(jds):
        t = ts.tt_jd(float(jd))
        apparent = observer.at(t).observe(target).apparent()
        # Get ecliptic-of-date Cartesian
        ecl_mat = ecliptic_frame.rotation_at(t)
        ecl_xyz = mxv(ecl_mat, apparent.position.au)
        results[i] = ecl_xyz
    return results


def eval_apparent_ecliptic_spherical(body_id, jds):
    """Evaluate calc at given TT JDs → (lon, lat, dist)."""
    ephem.set_calc_mode("skyfield")
    results = np.zeros((len(jds), 3))
    for i, jd in enumerate(jds):
        ref, _ = ephem.calc(float(jd), body_id, FLG_SPEED)
        results[i] = [ref[0], ref[1], ref[2]]
    return results


def ang_diff(a, b):
    d = abs(a - b)
    if d > 180:
        d = 360 - d
    return d


def cartesian_to_spherical(x, y, z):
    """Convert ecliptic Cartesian to spherical (lon_deg, lat_deg, dist_au)."""
    dist = math.sqrt(x * x + y * y + z * z)
    lon = math.degrees(math.atan2(y, x)) % 360.0
    lat = math.degrees(math.asin(max(-1.0, min(1.0, z / dist)))) if dist > 0 else 0.0
    return lon, lat, dist


def test_segment(body_id, body_name, seg_start, seg_end, degree, n_test=100):
    """Test Chebyshev fitting for one segment, comparing Cartesian vs Spherical."""
    mid = 0.5 * (seg_start + seg_end)
    half = 0.5 * (seg_end - seg_start)
    nodes_01 = chebyshev_nodes(degree + 1)
    node_jds = half * nodes_01 + mid

    # ── Evaluate at Chebyshev nodes ──
    cart_nodes = eval_apparent_ecliptic_cartesian(body_id, node_jds)
    sph_nodes = eval_apparent_ecliptic_spherical(body_id, node_jds)

    # ── Fit Cartesian ──
    coeffs_cart = np.zeros((3, degree + 1))
    for c in range(3):
        coeffs_cart[c] = chebfit(nodes_01, cart_nodes[:, c], degree)

    # ── Fit Spherical (with longitude unwrapping) ──
    sph_fit = sph_nodes.copy()
    sph_fit[:, 0] = np.degrees(np.unwrap(np.radians(sph_fit[:, 0])))
    coeffs_sph = np.zeros((3, degree + 1))
    for c in range(3):
        coeffs_sph[c] = chebfit(nodes_01, sph_fit[:, c], degree)

    # ── Verify at test points ──
    max_cart_lon_err = 0.0
    max_cart_lat_err = 0.0
    max_sph_lon_err = 0.0
    max_sph_lat_err = 0.0

    test_jds = [
        seg_start + (k + 0.5) / n_test * (seg_end - seg_start) for k in range(n_test)
    ]

    # Reference values at test points
    ref_sph = eval_apparent_ecliptic_spherical(body_id, test_jds)
    ref_cart = eval_apparent_ecliptic_cartesian(body_id, test_jds)

    for k in range(n_test):
        tau = (test_jds[k] - mid) / half

        # Cartesian fit → spherical
        cx = float(chebval(tau, coeffs_cart[0]))
        cy = float(chebval(tau, coeffs_cart[1]))
        cz = float(chebval(tau, coeffs_cart[2]))
        c_lon, c_lat, c_dist = cartesian_to_spherical(cx, cy, cz)

        # Spherical fit
        s_lon = float(chebval(tau, coeffs_sph[0])) % 360.0
        s_lat = float(chebval(tau, coeffs_sph[1]))

        # Errors vs reference
        cart_lon_err = ang_diff(c_lon, ref_sph[k, 0]) * 3600
        cart_lat_err = abs(c_lat - ref_sph[k, 1]) * 3600
        sph_lon_err = ang_diff(s_lon, ref_sph[k, 0]) * 3600
        sph_lat_err = abs(s_lat - ref_sph[k, 1]) * 3600

        max_cart_lon_err = max(max_cart_lon_err, cart_lon_err)
        max_cart_lat_err = max(max_cart_lat_err, cart_lat_err)
        max_sph_lon_err = max(max_sph_lon_err, sph_lon_err)
        max_sph_lat_err = max(max_sph_lat_err, sph_lat_err)

    return max_cart_lon_err, max_cart_lat_err, max_sph_lon_err, max_sph_lat_err


def main():
    JD_START = 2396758.5  # base tier start

    # Test cases: body, name, interval, degree, worst-case JD
    test_cases = [
        (6, "Saturn", 2.0, 15, 2501964.8),
        (6, "Saturn", 2.0, 15, 2451545.0),
        (5, "Jupiter", 0.5, 21, 2451545.0),
        (7, "Uranus", 1.0, 23, 2451545.0),
        (2, "Mercury", 1.0, 17, 2451545.0),
        (8, "Neptune", 4.0, 17, 2451545.0),
        (9, "Pluto", 8.0, 13, 2451545.0),
    ]

    h_cl = 'Cart lon"'
    h_cla = 'Cart lat"'
    h_sl = 'Sph lon"'
    h_sla = 'Sph lat"'
    print(
        f"{'Body':>10s} {'JD':>12s} {'Params':>8s}  "
        f"{h_cl:>10s} {h_cla:>10s}  "
        f"{h_sl:>10s} {h_sla:>10s}  "
        f"{'Winner':>8s}"
    )
    print("-" * 95)

    for body_id, name, interval, degree, target_jd in test_cases:
        seg_idx = int((target_jd - JD_START) / interval)
        seg_start = JD_START + seg_idx * interval
        seg_end = seg_start + interval

        cart_lon, cart_lat, sph_lon, sph_lat = test_segment(
            body_id, name, seg_start, seg_end, degree, n_test=100
        )

        cart_max = max(cart_lon, cart_lat)
        sph_max = max(sph_lon, sph_lat)
        winner = "CART" if cart_max < sph_max else "SPH"

        print(
            f"{name:>10s} {target_jd:>12.1f} {interval:.1f}d/{degree:2d}  "
            f'{cart_lon:10.4f}" {cart_lat:10.4f}"  '
            f'{sph_lon:10.4f}" {sph_lat:10.4f}"  '
            f"{winner:>8s}"
        )

    # Now do a broader scan for Saturn at worst-case: scan multiple segments
    # around retrogrades to find the true maximum
    print(f"\n{'=' * 95}")
    print("Saturn 2d/15: scanning 200 segments around worst-case regions...")
    print(f"{'=' * 95}")

    # Saturn retrogrades every ~378 days. Scan broadly.
    import random

    random.seed(42)

    interval = 2.0
    degree = 15
    body_id = 6

    # Scan ALL segments in the range to find worst Cartesian error
    n_total_segs = int(math.ceil((2506331.5 - JD_START) / interval))

    # Sample 500 random + known worst
    sample_indices = random.sample(range(n_total_segs), min(500, n_total_segs))
    # Add the known worst segment
    worst_idx = int((2501964.8 - JD_START) / interval)
    if worst_idx not in sample_indices:
        sample_indices.append(worst_idx)

    max_cart_lon = 0.0
    max_cart_lat = 0.0
    max_sph_lon = 0.0
    max_sph_lat = 0.0
    worst_cart_jd = 0.0
    worst_sph_jd = 0.0

    for count, si in enumerate(sample_indices):
        seg_start = JD_START + si * interval
        seg_end = seg_start + interval

        cl, cla, sl, sla = test_segment(
            body_id, "Saturn", seg_start, seg_end, degree, n_test=50
        )

        if max(cl, cla) > max(max_cart_lon, max_cart_lat):
            max_cart_lon = cl
            max_cart_lat = cla
            worst_cart_jd = seg_start
        if max(sl, sla) > max(max_sph_lon, max_sph_lat):
            max_sph_lon = sl
            max_sph_lat = sla
            worst_sph_jd = seg_start

        if (count + 1) % 50 == 0:
            print(
                f"  {count + 1}/{len(sample_indices)} segments scanned. "
                f'Cart worst: {max(max_cart_lon, max_cart_lat):.4f}" at JD {worst_cart_jd:.1f}  '
                f'Sph worst: {max(max_sph_lon, max_sph_lat):.4f}" at JD {worst_sph_jd:.1f}'
            )

    print(f"\n  FINAL Saturn 2d/15 ({len(sample_indices)} segments):")
    print(
        f'    Cartesian: max lon={max_cart_lon:.4f}"  max lat={max_cart_lat:.4f}"  worst JD={worst_cart_jd:.1f}'
    )
    print(
        f'    Spherical: max lon={max_sph_lon:.4f}"  max lat={max_sph_lat:.4f}"  worst JD={worst_sph_jd:.1f}'
    )
    print(
        f"    Improvement: {max(max_sph_lon, max_sph_lat) / max(max_cart_lon, max_cart_lat):.1f}x"
    )


if __name__ == "__main__":
    main()
