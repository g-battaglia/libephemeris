"""Measure actual max errors for medium tier LEB V3 across all tolerance categories.

Usage: python scripts/measure_medium_errors.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import libephemeris as ephem
from libephemeris.time_utils import julday

import argparse

_parser = argparse.ArgumentParser()
_parser.add_argument("--tier", default="medium", choices=["base", "medium", "extended"])
_args = _parser.parse_args()
TIER = _args.tier

LEB_PATH = os.path.join(
    os.path.dirname(__file__), "..", "data", "leb", f"ephemeris_{TIER}.leb"
)
LEB_PATH = os.path.abspath(LEB_PATH)

# Constants
FLG_SPEED = 256
FLG_EQUATORIAL = 2048
FLG_J2000 = 32
FLG_HELCTR = 8
FLG_BARYCTR = 16384
FLG_TRUEPOS = 16
FLG_NOABERR = 1024
FLG_SIDEREAL = 64 * 1024

MAIN_PLANETS = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14]
ECLIPTIC_BODIES = [10, 11, 12, 13, 21, 22]
ASTEROIDS = [15, 17, 18, 19, 20]
# Only IDs with a reviewed, independently documented runtime model belong in
# a numerical-error campaign.  Recognised IDs 48, 49, 54, 55, 57, and 58 are
# deliberately names-only compatibility constants and raise UnknownBodyError;
# measuring them would turn the provenance boundary into a misleading test
# failure.  ID 56 is the separately sourced uniform Selena convention rather
# than a row in fictitious_orbits.csv, but it is still a supported local model.
HYPOTHETICAL = [*range(40, 48), 50, 51, 52, 53, 56]

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
    12: "MeanApogee",
    13: "OscuApogee",
    14: "Earth",
    15: "Chiron",
    17: "Ceres",
    18: "Pallas",
    19: "Juno",
    20: "Vesta",
    21: "InterpApogee",
    22: "InterpPerigee",
    40: "Cupido",
    41: "Hades",
    42: "Zeus",
    43: "Kronos",
    44: "Apollon",
    45: "Admetos",
    46: "Vulkanus",
    47: "Poseidon",
    48: "Transpluto",
}


def angular_diff(a: float, b: float) -> float:
    d = abs(a - b)
    if d > 180.0:
        d = 360.0 - d
    return d


def generate_dates(
    n: int, jd_start: float, jd_end: float, margin: float = 30.0
) -> list[float]:
    s = jd_start + margin
    e = jd_end - margin
    step = (e - s) / (n - 1)
    return [s + i * step for i in range(n)]


def skyfield_calc(body: int, jd: float, flags: int = FLG_SPEED):
    ephem.state._LEB_FILE = None
    ephem.state._LEB_READER = None
    ephem.set_precision_tier(TIER)
    ephem.set_calc_mode("skyfield")
    try:
        result, _ = ephem.calc(jd, body, flags)
        return result
    finally:
        ephem.set_calc_mode(None)


def leb_calc(body: int, jd: float, flags: int = FLG_SPEED):
    ephem.state._LEB_FILE = LEB_PATH
    ephem.state._LEB_READER = None
    ephem.set_calc_mode("auto")
    try:
        result, _ = ephem.calc(jd, body, flags)
        return result
    finally:
        ephem.state._LEB_FILE = None
        ephem.state._LEB_READER = None
        ephem.set_calc_mode(None)


# Tier-specific date ranges
_TIER_RANGES = {
    "base": (1859, 2140),  # de440s: 1849-2150
    "medium": (1560, 2640),  # de440: 1550-2650
    "extended": (-4990, 4990),  # de441: -13200 to +17191 (test subset)
}
_yr_start, _yr_end = _TIER_RANGES[TIER]
JD_START = julday(_yr_start, 1, 1, 0.0)
JD_END = julday(_yr_end, 1, 1, 0.0)
# Asteroid SPK coverage
AST_JD_START = julday(1901, 1, 1, 0.0)
AST_JD_END = julday(2099, 1, 1, 0.0)

N_SAMPLES = 100  # Fast but representative


def measure_position(bodies, dates, label):
    """Measure lon/lat/dist errors for bodies."""
    max_lon_err = 0.0
    max_lat_err = 0.0
    max_dist_err = 0.0
    max_lon_body = ""
    max_lat_body = ""
    max_dist_body = ""

    for bid in bodies:
        body_dates = dates
        if bid in {15, 17, 18, 19, 20}:
            body_dates = [d for d in dates if AST_JD_START <= d <= AST_JD_END]

        for jd in body_dates:
            try:
                sf = skyfield_calc(bid, jd)
                lb = leb_calc(bid, jd)
            except Exception:
                continue

            lon_err = angular_diff(sf[0], lb[0]) * 3600  # arcsec
            lat_err = abs(sf[1] - lb[1]) * 3600
            dist_err = abs(sf[2] - lb[2])

            if lon_err > max_lon_err:
                max_lon_err = lon_err
                max_lon_body = f"{BODY_NAMES[bid]} JD={jd:.1f}"
            if lat_err > max_lat_err:
                max_lat_err = lat_err
                max_lat_body = f"{BODY_NAMES[bid]} JD={jd:.1f}"
            if dist_err > max_dist_err:
                max_dist_err = dist_err
                max_dist_body = f"{BODY_NAMES[bid]} JD={jd:.1f}"

    print(f"\n--- {label} ---")
    print(f'  Max lon error:  {max_lon_err:.4f}" ({max_lon_body})')
    print(f'  Max lat error:  {max_lat_err:.4f}" ({max_lat_body})')
    print(f"  Max dist error: {max_dist_err:.2e} AU ({max_dist_body})")
    return max_lon_err, max_lat_err, max_dist_err


def measure_speed(bodies, dates, label):
    """Measure speed errors (lon/lat/dist velocity)."""
    max_slon = 0.0
    max_slat = 0.0
    max_sdist = 0.0
    max_slon_body = ""
    max_slat_body = ""
    max_sdist_body = ""

    for bid in bodies:
        body_dates = dates
        if bid in {15, 17, 18, 19, 20}:
            body_dates = [d for d in dates if AST_JD_START <= d <= AST_JD_END]

        for jd in body_dates:
            try:
                sf = skyfield_calc(bid, jd)
                lb = leb_calc(bid, jd)
            except Exception:
                continue

            slon_err = abs(sf[3] - lb[3])  # deg/day
            slat_err = abs(sf[4] - lb[4])
            sdist_err = abs(sf[5] - lb[5])  # AU/day

            if slon_err > max_slon:
                max_slon = slon_err
                max_slon_body = f"{BODY_NAMES[bid]} JD={jd:.1f}"
            if slat_err > max_slat:
                max_slat = slat_err
                max_slat_body = f"{BODY_NAMES[bid]} JD={jd:.1f}"
            if sdist_err > max_sdist:
                max_sdist = sdist_err
                max_sdist_body = f"{BODY_NAMES[bid]} JD={jd:.1f}"

    print(f"\n--- {label} speed ---")
    print(
        f'  Max speed_lon:  {max_slon:.6f} deg/day = {max_slon * 3600:.2f}"/day ({max_slon_body})'
    )
    print(
        f'  Max speed_lat:  {max_slat:.6f} deg/day = {max_slat * 3600:.2f}"/day ({max_slat_body})'
    )
    print(f"  Max speed_dist: {max_sdist:.2e} AU/day ({max_sdist_body})")
    return max_slon, max_slat, max_sdist


def measure_flags(dates, label, extra_flags):
    """Measure position errors with non-default flags."""
    flags = FLG_SPEED | extra_flags
    test_bodies = [0, 1, 4, 5]  # Sun, Moon, Mars, Jupiter
    max_err = 0.0
    max_body = ""

    for bid in test_bodies:
        for jd in dates:
            try:
                sf = skyfield_calc(bid, jd, flags)
                lb = leb_calc(bid, jd, flags)
            except Exception:
                continue
            lon_err = angular_diff(sf[0], lb[0]) * 3600
            lat_err = abs(sf[1] - lb[1]) * 3600
            err = max(lon_err, lat_err)
            if err > max_err:
                max_err = err
                max_body = f"{BODY_NAMES[bid]} JD={jd:.1f}"

    print(f'  {label}: {max_err:.4f}" ({max_body})')
    return max_err


def measure_ecliptic_per_body(dates):
    """Measure per-body ecliptic errors."""
    print("\n--- Ecliptic bodies (per-body) ---")
    for bid in ECLIPTIC_BODIES:
        max_lon = 0.0
        max_slon = 0.0
        for jd in dates:
            try:
                sf = skyfield_calc(bid, jd)
                lb = leb_calc(bid, jd)
            except Exception:
                continue
            lon_err = angular_diff(sf[0], lb[0]) * 3600
            slon_err = abs(sf[3] - lb[3])
            max_lon = max(max_lon, lon_err)
            max_slon = max(max_slon, slon_err)
        print(
            f'  {BODY_NAMES[bid]:15s}: lon={max_lon:.4f}", speed_lon={max_slon:.6f} deg/day'
        )


def measure_sidereal(dates):
    """Measure sidereal position errors (mode 0 = Fagan-Bradley)."""
    print("\n--- Sidereal (mode 0, Fagan-Bradley) ---")
    test_bodies = [0, 1, 4, 5]
    flags = FLG_SPEED | FLG_SIDEREAL
    max_err = 0.0
    max_body = ""

    ephem.set_sid_mode(0)
    for bid in test_bodies:
        for jd in dates:
            try:
                sf = skyfield_calc(bid, jd, flags)
                lb = leb_calc(bid, jd, flags)
            except Exception:
                continue
            lon_err = angular_diff(sf[0], lb[0]) * 3600
            if lon_err > max_err:
                max_err = lon_err
                max_body = f"{BODY_NAMES[bid]} JD={jd:.1f}"

    print(f'  Max sidereal error: {max_err:.4f}" ({max_body})')
    return max_err


def main():
    print(f"LEB file: {LEB_PATH}")
    print(f"Samples: {N_SAMPLES}")

    # Enable auto SPK download for asteroid Skyfield reference
    ephem.set_auto_spk_download(True)

    dates = generate_dates(N_SAMPLES, JD_START, JD_END)

    # 1. Main planets position
    p_lon, p_lat, p_dist = measure_position(
        MAIN_PLANETS, dates, "Main planets (position)"
    )

    # 2. Main planets speed
    sp_lon, sp_lat, sp_dist = measure_speed(MAIN_PLANETS, dates, "Main planets")

    # 3. Asteroids position
    a_lon, a_lat, a_dist = measure_position(ASTEROIDS, dates, "Asteroids (position)")

    # 4. Asteroids speed
    sa_lon, sa_lat, sa_dist = measure_speed(ASTEROIDS, dates, "Asteroids")

    # 5. Ecliptic bodies
    measure_ecliptic_per_body(dates)

    # 6. Hypothetical
    h_lon, h_lat, h_dist = measure_position(
        HYPOTHETICAL, dates, "Hypothetical (position)"
    )

    # 7. Flag tests
    print("\n--- Flag tests ---")
    eq_err = measure_flags(dates, "EQUATORIAL", FLG_EQUATORIAL)
    j2k_err = measure_flags(dates, "J2000", FLG_J2000)

    # 8. Sidereal
    sid_err = measure_sidereal(dates)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY — suggested TIER_DEFAULTS (2x safety margin)")
    print("=" * 60)
    print(
        f'  POSITION_ARCSEC:          {max(p_lon, p_lat):.4f}" measured → {max(p_lon, p_lat) * 2:.2f}" suggested'
    )
    print(
        f'  ASTEROID_ARCSEC:          {max(a_lon, a_lat):.4f}" measured → {max(a_lon, a_lat) * 2:.4f}" suggested'
    )
    print(
        f'  EQUATORIAL_ARCSEC:        {eq_err:.4f}" measured → {eq_err * 2:.2f}" suggested'
    )
    print(
        f'  J2000_ARCSEC:             {j2k_err:.4f}" measured → {j2k_err * 2:.2f}" suggested'
    )
    print(
        f'  SIDEREAL_ARCSEC:          {sid_err:.4f}" measured → {sid_err * 2:.2f}" suggested'
    )
    print(
        f'  HYPOTHETICAL_ARCSEC:      {max(h_lon, h_lat):.6f}" measured → {max(h_lon, h_lat) * 2:.6f}" suggested'
    )
    print(
        f"  DISTANCE_AU:              {p_dist:.2e} measured → {p_dist * 2:.2e} suggested"
    )
    print(
        f"  SPEED_LON_DEG_DAY:        {sp_lon:.6f} measured → {sp_lon * 2:.6f} suggested"
    )
    print(
        f"  SPEED_LAT_DEG_DAY:        {sp_lat:.6f} measured → {sp_lat * 2:.6f} suggested"
    )
    print(
        f"  SPEED_DIST_AU_DAY:        {sp_dist:.2e} measured → {sp_dist * 2:.2e} suggested"
    )
    print(
        f"  ASTEROID_SPEED_LON:       {sa_lon:.6f} measured → {sa_lon * 2:.6f} suggested"
    )
    print(
        f"  ASTEROID_SPEED_LAT:       {sa_lat:.6f} measured → {sa_lat * 2:.6f} suggested"
    )
    print(
        f"  ASTEROID_SPEED_DIST:      {sa_dist:.2e} measured → {sa_dist * 2:.2e} suggested"
    )


if __name__ == "__main__":
    main()
