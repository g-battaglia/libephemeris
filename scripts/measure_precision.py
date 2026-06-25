#!/usr/bin/env python3
"""Measure end-to-end LEB precision: fast_calc vs calc (Skyfield reference).

Dense sampling across the full date range per body to find worst-case errors.
Reports per-body statistics: mean, P99, max error in arcseconds.
"""

from __future__ import annotations

import math
import os
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, ".")

# Ensure LEB is NOT loaded for reference calc calls.
# We pop the env var AND force skyfield mode so auto-discovery cannot
# silently route reference calls through LEB.
os.environ.pop("LIBEPHEMERIS_LEB", None)

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED
from libephemeris.leb_format import BODY_PARAMS

# Force skyfield mode globally so calc never auto-discovers LEB files
ephem.set_calc_mode("skyfield")

# Suppress MeeusPolynomialWarning spam from ecliptic bodies on extended tier
from libephemeris.lunar import MeeusPolynomialWarning

warnings.filterwarnings("ignore", category=MeeusPolynomialWarning)


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
    12: "MeanApog",
    13: "OscuApog",
    14: "Earth",
    15: "Chiron",
    17: "Ceres",
    18: "Pallas",
    19: "Juno",
    20: "Vesta",
    21: "IntpApog",
    22: "IntpPerg",
    40: "Cupido",
    41: "Hades",
    42: "Zeus",
    43: "Kronos",
    44: "Apollon",
    45: "Admetos",
    46: "Vulkanus",
    47: "Poseidon",
    48: "Isis",
}


def angular_sep_arcsec(lon1, lat1, lon2, lat2):
    """Great-circle angular separation in arcseconds."""
    lon1r, lat1r = math.radians(lon1), math.radians(lat1)
    lon2r, lat2r = math.radians(lon2), math.radians(lat2)
    dlon = lon2r - lon1r
    lat2r - lat1r
    # Vincenty formula for better numerical stability
    a = math.cos(lat2r) * math.sin(dlon)
    b = math.cos(lat1r) * math.sin(lat2r) - math.sin(lat1r) * math.cos(
        lat2r
    ) * math.cos(dlon)
    c = math.sin(lat1r) * math.sin(lat2r) + math.cos(lat1r) * math.cos(
        lat2r
    ) * math.cos(dlon)
    sep = math.atan2(math.sqrt(a * a + b * b), c)
    return abs(sep) * 3600.0 * 180.0 / math.pi


def lon_diff_arcsec(lon1, lon2):
    """Signed longitude difference in arcseconds, handling wrapping."""
    d = lon2 - lon1
    if d > 180:
        d -= 360
    elif d < -180:
        d += 360
    return d * 3600.0


def measure_body(ipl, reader, jd_start, jd_end, n_samples=2000):
    """Measure end-to-end error for a single body.

    Compares fast_calc (LEB) against calc (Skyfield reference).
    The global calc mode must be set to "skyfield" so calc never
    uses LEB auto-discovery.

    Args:
        ipl: Body ID.
        reader: Pre-opened LEBReader instance.
        jd_start: Start of JD range.
        jd_end: End of JD range.
        n_samples: Number of evenly-spaced sample points.

    Returns:
        Dict with per-body error statistics, or None if no valid samples.
    """
    from libephemeris.fast_calc import fast_calc_tt

    # Generate dense sample points
    jds = np.linspace(jd_start, jd_end, n_samples)

    errors_arcsec = []
    errors_lon = []
    errors_lat = []
    errors_dist = []
    worst_jd = 0.0
    worst_err = 0.0

    iflag = FLG_SPEED

    for jd in jds:
        jd_float = float(jd)

        # Reference: calc via Skyfield (calc_mode="skyfield" set globally)
        try:
            ref_result, _ = ephem.calc(jd_float, ipl, iflag)
        except Exception:
            continue

        ref_lon, ref_lat, ref_dist = ref_result[0], ref_result[1], ref_result[2]

        # LEB: fast_calc_tt directly (bypasses global mode)
        try:
            leb_result, _ = fast_calc_tt(reader, jd_float, ipl, iflag)
        except Exception:
            continue

        leb_lon, leb_lat, leb_dist = leb_result[0], leb_result[1], leb_result[2]

        # Angular separation
        sep = angular_sep_arcsec(ref_lon, ref_lat, leb_lon, leb_lat)
        errors_arcsec.append(sep)

        # Component errors
        dlon = lon_diff_arcsec(ref_lon, leb_lon)
        dlat = (leb_lat - ref_lat) * 3600.0
        ddist = abs(leb_dist - ref_dist)

        errors_lon.append(abs(dlon))
        errors_lat.append(abs(dlat))
        errors_dist.append(ddist)

        if sep > worst_err:
            worst_err = sep
            worst_jd = jd_float

    if not errors_arcsec:
        return None

    errors_arcsec = np.array(errors_arcsec)
    errors_lon = np.array(errors_lon)
    errors_lat = np.array(errors_lat)
    errors_dist = np.array(errors_dist)

    return {
        "body": ipl,
        "name": BODY_NAMES.get(ipl, f"body_{ipl}"),
        "n_samples": len(errors_arcsec),
        "mean_arcsec": float(np.mean(errors_arcsec)),
        "p99_arcsec": float(np.percentile(errors_arcsec, 99)),
        "max_arcsec": float(np.max(errors_arcsec)),
        "mean_lon_arcsec": float(np.mean(errors_lon)),
        "max_lon_arcsec": float(np.max(errors_lon)),
        "mean_lat_arcsec": float(np.mean(errors_lat)),
        "max_lat_arcsec": float(np.max(errors_lat)),
        "max_dist_au": float(np.max(errors_dist)),
        "worst_jd": worst_jd,
    }


def _detect_tier(leb_path: str) -> str:
    """Auto-detect precision tier from LEB filename.

    Matches patterns like 'ephemeris_base.leb', 'ephemeris_extended.leb',
    'ephemeris_medium.leb'. Falls back to 'medium' if no match.
    """
    basename = os.path.basename(leb_path).lower()
    for tier in ("base", "medium", "extended"):
        if tier in basename:
            return tier
    return "medium"


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Measure LEB precision: fast_calc vs calc (Skyfield reference)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "The --tier flag ensures the Skyfield reference uses the matching\n"
            "SPK file (de440s/de440/de441). If omitted, the tier is auto-detected\n"
            "from the LEB filename (e.g. ephemeris_extended.leb -> extended)."
        ),
    )
    parser.add_argument("--leb", required=True, help="Path to .leb file")
    parser.add_argument("--samples", type=int, default=2000, help="Samples per body")
    parser.add_argument("--bodies", nargs="*", type=int, help="Specific body IDs")
    parser.add_argument(
        "--group",
        choices=["planets", "asteroids", "ecliptic", "helio", "all"],
        default="all",
        help="Body group to test",
    )
    parser.add_argument(
        "--tier",
        choices=["base", "medium", "extended"],
        default=None,
        help="Precision tier for Skyfield reference (auto-detected from filename if omitted)",
    )
    args = parser.parse_args()

    if not os.path.exists(args.leb):
        print(f"ERROR: LEB file not found: {args.leb}")
        sys.exit(1)

    # Determine and set the precision tier so Skyfield loads the matching SPK
    tier = args.tier or _detect_tier(args.leb)
    ephem.set_precision_tier(tier)
    tier_info = ephem.TIERS[tier]
    print(f"Precision tier: {tier} ({tier_info.ephemeris_file})")

    # Read the LEB file to determine date range
    from libephemeris.leb_reader import LEBReader

    reader = LEBReader(args.leb)
    jd_start = reader._header.jd_start
    jd_end = reader._header.jd_end

    print(f"LEB file: {args.leb}")
    print(f"Date range: JD {jd_start:.1f} - {jd_end:.1f}")
    print(f"Samples per body: {args.samples}")
    print()

    # Determine which bodies to test
    if args.bodies:
        body_ids = args.bodies
    else:
        body_ids = sorted(BODY_PARAMS.keys())
        if args.group == "planets":
            body_ids = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14]
        elif args.group == "asteroids":
            body_ids = [15, 17, 18, 19, 20]
        elif args.group == "ecliptic":
            body_ids = [10, 11, 12, 13, 21, 22]
        elif args.group == "helio":
            body_ids = [40, 41, 42, 43, 44, 45, 46, 47, 48]

    # Filter to bodies in file
    available = set(reader._bodies.keys())
    body_ids = [b for b in body_ids if b in available]

    print(f"Testing {len(body_ids)} bodies: {[BODY_NAMES.get(b, b) for b in body_ids]}")
    print("=" * 100)

    # Coordinate type names for display
    COORD_NAMES = {
        0: "ICRS_BARY",
        1: "ECLIPTIC",
        2: "HELIO_ECL",
        4: "ICRS_SYS",
    }

    # Measure each body
    results = []
    for ipl in body_ids:
        name = BODY_NAMES.get(ipl, f"body_{ipl}")
        params = BODY_PARAMS.get(ipl)
        if params is None:
            print(f"  {name:12s} (id={ipl:2d}) ... SKIPPED (no BODY_PARAMS)")
            continue
        coord_name = COORD_NAMES.get(params[2], "?")

        t0 = time.time()
        print(
            f"  {name:12s} (id={ipl:2d}, {coord_name:10s}, {params[0]}d/{params[1]}) ...",
            end=" ",
            flush=True,
        )

        result = measure_body(ipl, reader, jd_start, jd_end, args.samples)
        elapsed = time.time() - t0

        if result is None:
            print("SKIPPED (no valid samples)")
            continue

        status = "OK" if result["max_arcsec"] < 0.001 else "FAIL"
        print(
            f'max={result["max_arcsec"]:.6f}"  p99={result["p99_arcsec"]:.6f}"  mean={result["mean_arcsec"]:.6f}"  [{status}]  ({elapsed:.1f}s)'
        )
        results.append(result)

    # Summary table
    print()
    print("=" * 100)
    hdr = (
        f"{'Body':>12s}  {'Coord':>10s}  {'Params':>7s}  "
        f"{'Mean':>10s}  {'P99':>10s}  {'Max':>10s}  "
        f"{'MaxLon':>10s}  {'MaxLat':>10s}  "
        f"{'MaxDist':>12s}  {'Status':>6s}"
    )
    print(hdr)
    print("-" * 100)

    all_pass = True
    for r in results:
        ipl = r["body"]
        params = BODY_PARAMS.get(ipl, (0, 0, 0, 0))
        coord_name = COORD_NAMES.get(params[2], "?")
        status = "PASS" if r["max_arcsec"] < 0.001 else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(
            f"{r['name']:>12s}  {coord_name:>10s}  {params[0]}d/{params[1]:>2d}  "
            f"{r['mean_arcsec']:>10.6f}  {r['p99_arcsec']:>10.6f}  {r['max_arcsec']:>10.6f}  "
            f"{r['max_lon_arcsec']:>10.6f}  {r['max_lat_arcsec']:>10.6f}  "
            f"{r['max_dist_au']:>12.2e}  {status:>6s}"
        )

    print("-" * 100)

    # Bodies exceeding target
    failing = [r for r in results if r["max_arcsec"] >= 0.001]
    if failing:
        print(f'\n*** {len(failing)} bodies EXCEED 0.001" target ***')
        for r in failing:
            params = BODY_PARAMS.get(r["body"], (0, 0, 0, 0))
            print(
                f'  {r["name"]:12s}: max={r["max_arcsec"]:.6f}"  '
                f"(worst JD={r['worst_jd']:.6f}, params={params[0]}d/{params[1]})"
            )
    else:
        print(f'\n*** ALL {len(results)} bodies PASS <0.001" target ***')

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
