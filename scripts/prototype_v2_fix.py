#!/usr/bin/env python3
"""Prototype: measure V2 ICRS_BARY pipeline error with and without deflection.

Simulates the runtime _pipeline_icrs from fast_calc.py with Skyfield-sourced
positions (eliminating Chebyshev error), then measures the ecliptic output
against calc(). Tests with:
  A) Current pipeline: no deflection, first-order aberration
  B) +Deflection: full PPN deflection by Sun/Jupiter/Saturn
  C) +Full SR aberration: Skyfield's relativistic formula
  D) B+C combined

If D gives ~0" error, the solution is to upgrade _pipeline_icrs.
"""

from __future__ import annotations

import math
import sys

import numpy as np

sys.path.insert(0, ".")

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED
from libephemeris.state import get_planets, get_timescale


C_AU_DAY = 173.1446326846693


def _first_order_aberration(geo, earth_vel):
    """Current LEB first-order Bradley aberration."""
    dist = math.sqrt(geo[0] ** 2 + geo[1] ** 2 + geo[2] ** 2)
    if dist == 0:
        return geo
    ux, uy, uz = geo[0] / dist, geo[1] / dist, geo[2] / dist
    vx = earth_vel[0] / C_AU_DAY
    vy = earth_vel[1] / C_AU_DAY
    vz = earth_vel[2] / C_AU_DAY
    dot = ux * vx + uy * vy + uz * vz
    ax = ux + vx - ux * dot
    ay = uy + vy - uy * dot
    az = uz + vz - uz * dot
    a_dist = math.sqrt(ax * ax + ay * ay + az * az)
    return (ax / a_dist * dist, ay / a_dist * dist, az / a_dist * dist)


def _full_sr_aberration(geo, earth_vel, light_time):
    """Skyfield's full special-relativistic aberration."""
    p1mag = light_time * C_AU_DAY
    vemag = math.sqrt(earth_vel[0] ** 2 + earth_vel[1] ** 2 + earth_vel[2] ** 2)
    beta = vemag / C_AU_DAY
    dot = geo[0] * earth_vel[0] + geo[1] * earth_vel[1] + geo[2] * earth_vel[2]

    if p1mag * vemag > 0:
        cosd = dot / (p1mag * vemag)
    else:
        cosd = 0.0

    gammai = math.sqrt(1.0 - beta * beta)
    p = beta * cosd
    q = (1.0 + p / (1.0 + gammai)) * light_time
    r = 1.0 + p
    if abs(r) < 1e-30:
        r = 1.0

    return (
        (gammai * geo[0] + q * earth_vel[0]) / r,
        (gammai * geo[1] + q * earth_vel[1]) / r,
        (gammai * geo[2] + q * earth_vel[2]) / r,
    )


def _gravitational_deflection(geo, earth_bary, jd_tt, light_time, planets, ts):
    """Apply PPN gravitational deflection by Sun, Jupiter, Saturn."""
    GS = 1.32712440017987e20
    C = 299792458.0
    AU_M = 149597870700

    deflectors = [
        ("sun", 1.0),
        ("jupiter barycenter", 1047.3486),
        ("saturn barycenter", 3497.898),
    ]

    result = list(geo)
    pmag = math.sqrt(result[0] ** 2 + result[1] ** 2 + result[2] ** 2)
    if pmag == 0:
        return tuple(result)

    t_obs = ts.tt_jd(jd_tt)

    for defl_name, rmass in deflectors:
        deflector = planets[defl_name]

        # Deflector bary position at observation time
        db = deflector.at(t_obs).position.au
        db = (float(db[0]), float(db[1]), float(db[2]))

        # Deflector relative to observer
        gpv = (db[0] - earth_bary[0], db[1] - earth_bary[1], db[2] - earth_bary[2])

        # Unit vector to target
        phat = (result[0] / pmag, result[1] / pmag, result[2] / pmag)

        # Light-time difference: when did light pass closest to deflector
        dlt = (phat[0] * gpv[0] + phat[1] * gpv[1] + phat[2] * gpv[2]) / C_AU_DAY

        # Clamp
        tclose_offset = max(0.0, min(dlt, light_time))
        tclose_jd = jd_tt - tclose_offset

        # Deflector at closest approach time
        t_close = ts.tt_jd(tclose_jd)
        db_close = deflector.at(t_close).position.au
        db_close = (float(db_close[0]), float(db_close[1]), float(db_close[2]))

        # pe = observer - deflector
        pe = (
            earth_bary[0] - db_close[0],
            earth_bary[1] - db_close[1],
            earth_bary[2] - db_close[2],
        )

        # pq = target - deflector (from observer frame)
        pq = (result[0] + pe[0], result[1] + pe[1], result[2] + pe[2])

        qmag = math.sqrt(pq[0] ** 2 + pq[1] ** 2 + pq[2] ** 2)
        emag = math.sqrt(pe[0] ** 2 + pe[1] ** 2 + pe[2] ** 2)

        if qmag == 0 or emag == 0:
            continue

        phat2 = (result[0] / pmag, result[1] / pmag, result[2] / pmag)
        qhat = (pq[0] / qmag, pq[1] / qmag, pq[2] / qmag)
        ehat = (pe[0] / emag, pe[1] / emag, pe[2] / emag)

        pdotq = phat2[0] * qhat[0] + phat2[1] * qhat[1] + phat2[2] * qhat[2]
        qdote = qhat[0] * ehat[0] + qhat[1] * ehat[1] + qhat[2] * ehat[2]
        edotp = ehat[0] * phat2[0] + ehat[1] * phat2[1] + ehat[2] * phat2[2]

        if abs(edotp) > 0.99999999999:
            continue

        fac1 = 2.0 * GS / (C * C * emag * AU_M * rmass)
        fac2 = 1.0 + qdote
        if abs(fac2) < 1e-30:
            continue

        coeff = fac1 / fac2 * pmag
        result[0] += coeff * (pdotq * ehat[0] - edotp * qhat[0])
        result[1] += coeff * (pdotq * ehat[1] - edotp * qhat[1])
        result[2] += coeff * (pdotq * ehat[2] - edotp * qhat[2])

    return tuple(result)


def _icrs_to_ecliptic(geo, jd_tt, ts):
    """ICRS → ecliptic-of-date via PNM + true obliquity (matches Skyfield)."""
    t = ts.tt_jd(jd_tt)
    M = t.M
    # M is (3,3) for scalar time
    pn_mat = np.array(M)

    # PNM rotation
    eq = pn_mat @ np.array(geo)

    # True obliquity
    mean_obl = float(t._mean_obliquity_radians)
    _dpsi, deps = t._nutation_angles_radians
    eps_true = mean_obl + float(deps)
    cos_eps = math.cos(eps_true)
    sin_eps = math.sin(eps_true)

    ecl_x = eq[0]
    ecl_y = eq[1] * cos_eps + eq[2] * sin_eps
    ecl_z = -eq[1] * sin_eps + eq[2] * cos_eps

    dist = math.sqrt(ecl_x**2 + ecl_y**2 + ecl_z**2)
    lon = math.degrees(math.atan2(ecl_y, ecl_x)) % 360.0
    lat = (
        math.degrees(math.asin(max(-1.0, min(1.0, ecl_z / dist)))) if dist > 0 else 0.0
    )

    return lon, lat, dist


def simulate_pipeline(body_id, jd_tt, *, deflection=False, sr_aberration=False):
    """Simulate _pipeline_icrs with optional improvements."""
    planets = get_planets()
    ts = get_timescale()

    from libephemeris.planets import _PLANET_MAP, get_planet_target

    target_name = _PLANET_MAP[body_id]
    target = get_planet_target(planets, target_name)
    observer = planets["earth"]

    t = ts.tt_jd(jd_tt)

    # Earth position + velocity
    earth_at = observer.at(t)
    earth_pos = tuple(float(x) for x in earth_at.position.au)
    earth_vel = tuple(float(x) for x in earth_at.velocity.au_per_d)

    # Target position
    target_pos = tuple(float(x) for x in target.at(t).position.au)

    # Geocentric
    geo = (
        target_pos[0] - earth_pos[0],
        target_pos[1] - earth_pos[1],
        target_pos[2] - earth_pos[2],
    )

    # Light-time (3 iterations)
    lt = 0.0
    for _ in range(3):
        dist = math.sqrt(geo[0] ** 2 + geo[1] ** 2 + geo[2] ** 2)
        lt = dist / C_AU_DAY
        t_ret = ts.tt_jd(jd_tt - lt)
        ret_pos = tuple(float(x) for x in target.at(t_ret).position.au)
        geo = (
            ret_pos[0] - earth_pos[0],
            ret_pos[1] - earth_pos[1],
            ret_pos[2] - earth_pos[2],
        )

    # Deflection
    if deflection:
        geo = _gravitational_deflection(geo, earth_pos, jd_tt, lt, planets, ts)

    # Aberration
    if sr_aberration:
        geo = _full_sr_aberration(geo, earth_vel, lt)
    else:
        geo = _first_order_aberration(geo, earth_vel)

    # ICRS → ecliptic of date
    lon, lat, dist = _icrs_to_ecliptic(geo, jd_tt, ts)
    return lon, lat, dist


def ang_diff(a, b):
    d = abs(a - b)
    if d > 180:
        d = 360 - d
    return d


def main():
    # Bodies to test
    bodies = [
        (0, "Sun"),
        (1, "Moon"),
        (2, "Mercury"),
        (3, "Venus"),
        (4, "Mars"),
        (5, "Jupiter"),
        (6, "Saturn"),
        (7, "Uranus"),
        (8, "Neptune"),
        (9, "Pluto"),
    ]

    # Test JDs — spread across the base tier, including worst-case dates
    test_jds = [
        2396800.0,  # 1850
        2415020.0,  # 1900
        2433282.0,  # 1950
        2451545.0,  # 2000 (J2000)
        2460000.0,  # 2023
        2469807.0,  # 2050
        2488069.0,  # 2100
        2501964.8,  # 2138 (Saturn worst-case)
        2506300.0,  # 2150 (near end)
    ]

    ephem.set_calc_mode("skyfield")

    print("Pipeline comparison: measuring error from missing deflection/aberration")
    print("=" * 110)

    # For each body, test all 4 pipeline variants at all dates
    for body_id, body_name in bodies:
        max_err_A = 0.0  # No deflection, 1st-order aberration
        max_err_B = 0.0  # +Deflection
        max_err_C = 0.0  # +Full SR aberration
        max_err_D = 0.0  # +Both
        worst_jd_A = 0.0

        for jd_tt in test_jds:
            # Reference: calc
            ref, _ = ephem.calc(jd_tt, body_id, FLG_SPEED)

            # A: current pipeline (no defl, 1st-order aber)
            lon_a, lat_a, _ = simulate_pipeline(
                body_id, jd_tt, deflection=False, sr_aberration=False
            )
            err_a = max(ang_diff(lon_a, ref[0]) * 3600, abs(lat_a - ref[1]) * 3600)
            if err_a > max_err_A:
                max_err_A = err_a
                worst_jd_A = jd_tt

            # B: +deflection
            lon_b, lat_b, _ = simulate_pipeline(
                body_id, jd_tt, deflection=True, sr_aberration=False
            )
            err_b = max(ang_diff(lon_b, ref[0]) * 3600, abs(lat_b - ref[1]) * 3600)
            max_err_B = max(max_err_B, err_b)

            # C: +SR aberration (no deflection)
            lon_c, lat_c, _ = simulate_pipeline(
                body_id, jd_tt, deflection=False, sr_aberration=True
            )
            err_c = max(ang_diff(lon_c, ref[0]) * 3600, abs(lat_c - ref[1]) * 3600)
            max_err_C = max(max_err_C, err_c)

            # D: both
            lon_d, lat_d, _ = simulate_pipeline(
                body_id, jd_tt, deflection=True, sr_aberration=True
            )
            err_d = max(ang_diff(lon_d, ref[0]) * 3600, abs(lat_d - ref[1]) * 3600)
            max_err_D = max(max_err_D, err_d)

        print(
            f"  {body_name:>10s}: "
            f'A={max_err_A:8.4f}"  '
            f'B(+defl)={max_err_B:8.4f}"  '
            f'C(+SR)={max_err_C:8.4f}"  '
            f'D(both)={max_err_D:8.4f}"  '
            f"worst_A@JD={worst_jd_A:.1f}"
        )

    print()
    print("Legend:")
    print("  A = current pipeline (no deflection, 1st-order aberration)")
    print("  B = A + gravitational deflection")
    print("  C = A + full SR aberration")
    print("  D = B + C (full fix)")
    print()
    print('If D ≈ 0.000" → upgrading _pipeline_icrs eliminates the conversion error')


if __name__ == "__main__":
    main()
