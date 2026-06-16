#!/usr/bin/env python3
"""JPL Horizons cross-validation for libephemeris.

Queries NASA JPL Horizons for geocentric apparent ecliptic positions
and compares against libephemeris calculations.  This provides independent
ground-truth validation using the authoritative NASA source.

Validation Plan v2, Section 1.

Usage:
    .venv/bin/python3 scripts/horizons_cross_validate.py

Methodology:
    - Horizons ObsEclLon / ObsEclLat = apparent ecliptic of date
    - libephemeris  calc_ut(jd, body)  = apparent ecliptic of date (default)
    - Both referenced from geocentric (Horizons location='500')
    - Both interpret input JD as UT (Universal Time)
    - Internally both convert UT→TDB for computation, but use different
      Delta T models.  This causes growing divergence for dates far from
      the present where Delta T is uncertain/extrapolated.

Tolerance tiers:
    - "core era" (1972-2040): Delta T precisely known from observations
      → tight tolerance (0.5" planets, 1.0" Moon)
    - "historical era" (1900-1972): Delta T well-modeled
      → moderate tolerance (1.0" planets, 2.0" Moon)
    - "extrapolated era" (2040-2100): Delta T extrapolated, models diverge
      → relaxed tolerance (2.0" planets, 20.0" Moon)

Deliverable for validation-plan-v2.md §1.
"""

from __future__ import annotations

import json
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import libephemeris as swe  # noqa: E402
from astroquery.jplhorizons import Horizons  # noqa: E402


# ============================================================================
# Body mappings
# ============================================================================

PLANETARY_BODIES: dict[str, tuple[str, int, str]] = {
    # name: (horizons_id, swe_body_id, id_type)
    "Sun": ("10", swe.SUN, "majorbody"),
    "Moon": ("301", swe.MOON, "majorbody"),
    "Mercury": ("199", swe.MERCURY, "majorbody"),
    "Venus": ("299", swe.VENUS, "majorbody"),
    "Mars": ("499", swe.MARS, "majorbody"),
    "Jupiter": ("599", swe.JUPITER, "majorbody"),
    "Saturn": ("699", swe.SATURN, "majorbody"),
    "Uranus": ("799", swe.URANUS, "majorbody"),
    "Neptune": ("899", swe.NEPTUNE, "majorbody"),
    "Pluto": ("999", swe.PLUTO, "majorbody"),
}

MINOR_BODIES: dict[str, tuple[str, int, str]] = {
    "Chiron": ("2060", swe.CHIRON, "smallbody"),
    "Ceres": ("Ceres", swe.CERES, "smallbody"),
    "Pallas": ("Pallas", swe.PALLAS, "smallbody"),
    "Juno": ("Juno", swe.JUNO, "smallbody"),
    "Vesta": ("Vesta", swe.VESTA, "smallbody"),
}

COB_BODIES: dict[str, tuple[str, str, int]] = {
    # name: (body_center_id, barycenter_id, swe_body_id)
    "Jupiter": ("599", "5", swe.JUPITER),
    "Saturn": ("699", "6", swe.SATURN),
    "Neptune": ("899", "8", swe.NEPTUNE),
    "Pluto": ("999", "9", swe.PLUTO),
}


# ============================================================================
# Era-adaptive tolerances
# ============================================================================

# JD boundaries for tolerance eras
JD_1900 = 2415020.5
JD_1972 = 2441317.5  # UTC leap seconds era begins
JD_2040 = 2466154.5  # ~2040-01-01
JD_2100 = 2488069.5


def _get_tolerance(jd: float, body_name: str) -> float:
    """Return era-adaptive tolerance in arcseconds.

    Delta T is precisely known from observations for 1972-present.
    Before 1972 it relies on historical observations (still good).
    After ~2030 it's extrapolated, and different models diverge.
    The Moon is most sensitive due to its ~0.55"/sec apparent speed.
    """
    is_moon = body_name == "Moon"

    if JD_1972 <= jd <= JD_2040:
        # Core era: DeltaT from observations, < 1s uncertainty
        return 1.0 if is_moon else 0.5
    elif JD_1900 <= jd < JD_1972:
        # Historical era: DeltaT well-modeled, ~few seconds uncertainty
        return 2.0 if is_moon else 1.0
    else:
        # Extrapolated era (>2040 or <1900): DeltaT diverges between models
        # Moon moves ~0.55"/s, so 30s DeltaT difference → ~16" Moon error
        return 20.0 if is_moon else 2.0


# ============================================================================
# Date generation
# ============================================================================


def generate_planetary_dates(n: int = 50) -> list[float]:
    """Generate *n* deterministic dates spanning 1900-2100.

    Includes J2000.0, equinoxes/solstices, and evenly spaced dates.
    """
    key_dates = [
        2415020.5,  # 1900-01-01
        2420000.0,  # 1913-07-07
        2425000.0,  # 1927-01-06
        2430000.0,  # 1940-07-08
        2435000.0,  # 1954-01-08
        2440000.0,  # 1967-07-10
        2441317.5,  # 1972-01-01 (UTC leap seconds start)
        2445000.0,  # 1981-07-12
        2447892.5,  # 1990-01-01
        2450000.0,  # 1995-10-10
        2451179.5,  # 1999-01-01
        2451545.0,  # 2000-01-01.5 (J2000.0)
        2451623.8,  # ~2000-03-20 (vernal equinox)
        2451716.3,  # ~2000-06-21 (summer solstice)
        2451808.8,  # ~2000-09-22 (autumn equinox)
        2451900.1,  # ~2000-12-21 (winter solstice)
        2453005.5,  # 2004-01-01
        2455197.5,  # 2010-01-01
        2457388.5,  # 2016-01-01
        2458849.5,  # 2020-01-01
        2460310.5,  # 2024-01-01
        2460386.8,  # ~2024-03-20 (vernal equinox)
        2460479.7,  # ~2024-06-20 (summer solstice)
        2460572.5,  # ~2024-09-22 (autumn equinox)
        2460664.2,  # ~2024-12-21 (winter solstice)
        2461041.5,  # 2026-01-01
        2465000.0,  # 2036-11-10
        2470000.0,  # 2050-05-12
        2475000.0,  # 2063-11-11
        2480000.0,  # 2077-05-13
        2488069.5,  # 2100-01-01
    ]

    # Fill remaining with reproducible pseudo-random dates
    rng = np.random.default_rng(seed=42)
    while len(key_dates) < n:
        jd = float(rng.uniform(2415020.5, 2488069.5))
        key_dates.append(round(jd, 1))

    return sorted(key_dates[:n])


def generate_minor_body_dates(n: int = 20) -> list[float]:
    """Generate *n* dates for minor bodies (1950-2050, within SPK range)."""
    rng = np.random.default_rng(seed=123)
    jd_start, jd_end = 2433282.5, 2469807.5  # 1950–2050
    key = [2451545.0, 2460311.0]  # J2000.0, 2024-01-01
    while len(key) < n:
        key.append(round(float(rng.uniform(jd_start, jd_end)), 1))
    return sorted(key[:n])


def generate_cob_dates(n: int = 20) -> list[float]:
    """Generate *n* dates for COB tests (1950-2050)."""
    rng = np.random.default_rng(seed=456)
    jd_start, jd_end = 2433282.5, 2469807.5
    key = [2451545.0, 2460311.0]
    while len(key) < n:
        key.append(round(float(rng.uniform(jd_start, jd_end)), 1))
    return sorted(key[:n])


# ============================================================================
# Horizons query helper
# ============================================================================


def query_horizons(
    body_id: str,
    epochs: list[float],
    id_type: str = "majorbody",
) -> list[dict]:
    """Query Horizons for ecliptic positions, returns list of dicts."""
    obj = Horizons(
        id=body_id,
        location="500",  # geocentric
        epochs=epochs,
        id_type=id_type,
    )
    eph = obj.ephemerides(quantities="1,20,31")

    results = []
    for row in eph:
        ecl_lon = row["ObsEclLon"]
        ecl_lat = row["ObsEclLat"]
        # Skip masked / missing values
        if hasattr(ecl_lon, "mask") or str(ecl_lon) == "--":
            continue
        results.append(
            {
                "jd": float(row["datetime_jd"]),
                "ecl_lon": float(ecl_lon),
                "ecl_lat": float(ecl_lat),
                "ra": float(row["RA"]),
                "dec": float(row["DEC"]),
                "dist": (float(row["delta"]) if str(row["delta"]) != "--" else None),
            }
        )
    return results


# ============================================================================
# Comparison logic
# ============================================================================


def _wrap_delta(d: float) -> float:
    """Wrap longitude difference into [-180, +180] degrees."""
    if d > 180:
        d -= 360
    elif d < -180:
        d += 360
    return d


def compare_body(
    name: str,
    horizons_id: str,
    swe_body: int,
    dates: list[float],
    id_type: str = "majorbody",
    fixed_tolerance: float | None = None,
) -> dict | None:
    """Compare libephemeris vs Horizons for a single body.

    Uses era-adaptive tolerance unless fixed_tolerance is given.
    """
    print(f"  {name:18s}", end="", flush=True)

    try:
        horizons_data = query_horizons(horizons_id, dates, id_type=id_type)
    except Exception as e:
        print(f" HORIZONS ERROR: {e}")
        return None

    if not horizons_data:
        print(" NO DATA from Horizons")
        return None

    details: list[dict] = []
    lon_errs: list[float] = []
    lat_errs: list[float] = []

    for h in horizons_data:
        jd = h["jd"]
        tol = fixed_tolerance if fixed_tolerance else _get_tolerance(jd, name)
        try:
            result, _flags = swe.calc_ut(jd, swe_body)
            lib_lon, lib_lat = result[0], result[1]

            dlon = _wrap_delta(lib_lon - h["ecl_lon"]) * 3600  # arcsec
            dlat = (lib_lat - h["ecl_lat"]) * 3600

            lon_errs.append(abs(dlon))
            lat_errs.append(abs(dlat))
            passed = abs(dlon) < tol and abs(dlat) < tol

            details.append(
                {
                    "jd": jd,
                    "horizons_lon": h["ecl_lon"],
                    "horizons_lat": h["ecl_lat"],
                    "lib_lon": round(lib_lon, 8),
                    "lib_lat": round(lib_lat, 8),
                    "dlon_arcsec": round(dlon, 6),
                    "dlat_arcsec": round(dlat, 6),
                    "tolerance_arcsec": tol,
                    "passed": passed,
                }
            )
        except Exception as e:
            details.append({"jd": jd, "error": str(e), "passed": False})

    n_pass = sum(d.get("passed", False) for d in details)
    n_total = len(details)
    max_lon = max(lon_errs) if lon_errs else 0.0
    max_lat = max(lat_errs) if lat_errs else 0.0
    max_err = max(max_lon, max_lat)
    mean_lon = float(np.mean(lon_errs)) if lon_errs else 0.0
    mean_lat = float(np.mean(lat_errs)) if lat_errs else 0.0

    # Core-era stats (1972-2040) — the meaningful precision measurement
    core_lon = [
        abs(d["dlon_arcsec"])
        for d in details
        if "dlon_arcsec" in d and JD_1972 <= d["jd"] <= JD_2040
    ]
    core_lat = [
        abs(d["dlat_arcsec"])
        for d in details
        if "dlat_arcsec" in d and JD_1972 <= d["jd"] <= JD_2040
    ]
    core_max = max(max(core_lon, default=0), max(core_lat, default=0))
    core_mean_lon = float(np.mean(core_lon)) if core_lon else 0.0

    status = "PASS" if n_pass == n_total else "FAIL"
    print(
        f"  {status}  {n_pass:3d}/{n_total:3d}  "
        f'max={max_err:8.4f}"  core_max={core_max:.4f}"  '
        f'core_mean={core_mean_lon:.4f}"'
    )

    return {
        "name": name,
        "horizons_id": horizons_id,
        "swe_body": swe_body,
        "dates_tested": n_total,
        "dates_passed": n_pass,
        "max_lon_err_arcsec": round(max_lon, 6),
        "max_lat_err_arcsec": round(max_lat, 6),
        "mean_lon_err_arcsec": round(mean_lon, 6),
        "mean_lat_err_arcsec": round(mean_lat, 6),
        "core_era_max_err_arcsec": round(core_max, 6),
        "core_era_mean_lon_arcsec": round(core_mean_lon, 6),
        "core_era_n_dates": len(core_lon),
        "pass": n_pass == n_total,
        "details": details,
    }


# ============================================================================
# Main
# ============================================================================


def main() -> int:
    print("=" * 78)
    print("  JPL Horizons Cross-Validation for libephemeris")
    print("=" * 78)
    print()
    print("  Tolerance model: era-adaptive (accounts for Delta T divergence)")
    print('    Core (1972-2040):    0.5" planets, 1.0" Moon')
    print('    Historical (<1972):  1.0" planets, 2.0" Moon')
    print('    Extrapolated (>2040): 2.0" planets, 20.0" Moon')
    print()

    report: dict = {
        "metadata": {
            "generated": datetime.now(timezone.utc).isoformat(),
            "libephemeris_version": getattr(swe, "__version__", "unknown"),
            "horizons_location": "500 (geocentric)",
            "coordinate_system": "apparent ecliptic of date (ObsEclLon/ObsEclLat)",
            "time_scale": "UT (both Horizons observer tables and calc_ut use UT)",
            "methodology": (
                "Both Horizons and calc_ut interpret input JD as UT. "
                "Internally each converts UT→TDB using its own Delta T model. "
                "For dates far from the present (especially >2040), Delta T "
                "models diverge, causing systematic errors proportional to "
                "body speed × DeltaT difference. Era-adaptive tolerances "
                "account for this expected divergence."
            ),
            "tolerance_model": {
                "core_era": '1972-2040: 0.5" planets, 1.0" Moon',
                "historical_era": '<1972: 1.0" planets, 2.0" Moon',
                "extrapolated_era": '>2040: 2.0" planets, 20.0" Moon',
            },
        },
        "sections": {},
    }

    # ---------------------------------------------------------------
    # §1.1  Planetary Positions  (10 bodies × 50 dates)
    # ---------------------------------------------------------------
    print("§1.1  Planetary Positions  (10 bodies × 50 dates)")
    print("-" * 78)
    dates_1_1 = generate_planetary_dates(50)
    print(
        f"  Date range: JD {dates_1_1[0]:.1f} – {dates_1_1[-1]:.1f}"
        f"  ({len(dates_1_1)} dates)"
    )
    print()

    sec_1_1: dict = {}
    for bname, (hid, sid, idt) in PLANETARY_BODIES.items():
        res = compare_body(bname, hid, sid, dates_1_1, id_type=idt)
        if res:
            sec_1_1[bname] = res
        time.sleep(0.3)

    report["sections"]["1.1_planetary_positions"] = sec_1_1

    # ---------------------------------------------------------------
    # §1.2  Lunar Node & Lilith — methodology note
    # ---------------------------------------------------------------
    print()
    print("§1.2  Lunar Node & Lilith")
    print("-" * 78)
    print(
        "  True/Mean Node and Lilith are specialised derived quantities\n"
        "  not directly available as Horizons body positions.\n"
        "  Validated via the reference ephemeris hyper-validation (4400 rounds, 0 FAIL)."
    )
    report["sections"]["1.2_lunar_node_lilith"] = {
        "note": (
            "True/Mean Node and Lilith validated via the reference ephemeris hyper-validation "
            "(4400 rounds, 0 FAIL). These derived astrological quantities are not "
            "directly available as Horizons body positions."
        ),
        "validation_method": "reference_ephemeris_comparison",
    }

    # ---------------------------------------------------------------
    # §1.3  Minor Bodies with SPK  (5 bodies × 20 dates)
    # ---------------------------------------------------------------
    print()
    print("§1.3  Minor Bodies with SPK  (5 bodies × 20 dates)")
    print("-" * 78)
    dates_1_3 = generate_minor_body_dates(20)
    print(
        f"  Date range: JD {dates_1_3[0]:.1f} – {dates_1_3[-1]:.1f}"
        f"  ({len(dates_1_3)} dates)"
    )
    print()

    sec_1_3: dict = {}
    for bname, (hid, sid, idt) in MINOR_BODIES.items():
        res = compare_body(bname, hid, sid, dates_1_3, id_type=idt, fixed_tolerance=1.0)
        if res:
            sec_1_3[bname] = res
        time.sleep(0.3)

    report["sections"]["1.3_minor_bodies"] = sec_1_3

    # ---------------------------------------------------------------
    # §1.4  Outer Planet COB Corrections  (4 bodies × 20 dates)
    # ---------------------------------------------------------------
    print()
    print("§1.4  Outer Planet COB Corrections  (4 bodies × 20 dates)")
    print("-" * 78)
    dates_1_4 = generate_cob_dates(20)
    print(
        f"  Date range: JD {dates_1_4[0]:.1f} – {dates_1_4[-1]:.1f}"
        f"  ({len(dates_1_4)} dates)"
    )
    print()

    sec_1_4: dict = {}
    for bname, (body_id, _bary_id, sid) in COB_BODIES.items():
        res = compare_body(
            f"{bname} (body ctr)",
            body_id,
            sid,
            dates_1_4,
            id_type="majorbody",
            fixed_tolerance=0.5,
        )
        if res:
            sec_1_4[f"{bname}_body_center"] = res
        time.sleep(0.3)

    report["sections"]["1.4_outer_planet_cob"] = sec_1_4

    # ---------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------
    print()
    print("=" * 78)
    print("  SUMMARY")
    print("=" * 78)

    total_pass = 0
    total_tests = 0
    failures: list[str] = []

    for sec_name, sec_data in report["sections"].items():
        if not isinstance(sec_data, dict) or "note" in sec_data:
            continue
        for body_name, body_data in sec_data.items():
            total_pass += body_data["dates_passed"]
            total_tests += body_data["dates_tested"]
            if not body_data["pass"]:
                me = max(
                    body_data["max_lon_err_arcsec"],
                    body_data["max_lat_err_arcsec"],
                )
                failures.append(
                    f"  FAIL  {sec_name}/{body_name}: "
                    f'max_err={me:.4f}" core_max='
                    f'{body_data["core_era_max_err_arcsec"]:.4f}"'
                )

    if failures:
        for f in failures:
            print(f)
    print(f"\n  Total: {total_pass}/{total_tests} date-body comparisons passed")

    all_passed = len(failures) == 0
    if all_passed:
        print("  STATUS: ALL PASS")
    else:
        print(f"  STATUS: {len(failures)} BODY/SECTION FAILURES")

    report["summary"] = {
        "total_comparisons": total_tests,
        "total_passed": total_pass,
        "all_passed": all_passed,
    }

    # Save JSON report
    out = (
        Path(__file__).resolve().parent.parent
        / "data"
        / "horizons_cross_validation.json"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        json.dump(report, fh, indent=2, default=str)
    print(
        f"\n  Report saved to:"
        f" {out.relative_to(Path(__file__).resolve().parent.parent)}"
    )

    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
