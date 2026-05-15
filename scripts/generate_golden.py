"""Generate golden reference files for regression testing.

Creates a JSON file with 100 representative calculations spanning the full
API surface. This file is used by test_golden_regression.py to detect
any unexpected changes in calculation results.

Usage:
    .venv/bin/python3 scripts/generate_golden.py

Output:
    tests/golden/golden_reference.json
"""

from __future__ import annotations

import json
import math
import sys
import time

sys.path.insert(0, ".")
import libephemeris as swe  # noqa: E402

# ── Configuration ─────────────────────────────────────────────────────────────

OUTPUT_FILE = "tests/golden/golden_reference.json"

# Julian dates spanning key eras
JDS = [
    2415020.5,  # 1900-01-01
    2431545.0,  # 1945-01-15
    2440587.5,  # 1970-01-01
    2451545.0,  # 2000-01-01 12:00 (J2000.0)
    2455197.5,  # 2010-01-01
    2459580.5,  # 2022-01-01
    2460676.5,  # 2025-01-15
    2469807.5,  # 2050-01-01
]

BODIES = [
    swe.SUN,  # 0
    swe.MOON,  # 1
    swe.MERCURY,  # 2
    swe.VENUS,  # 3
    swe.MARS,  # 4
    swe.JUPITER,  # 5
    swe.SATURN,  # 6
    swe.URANUS,  # 7
    swe.NEPTUNE,  # 8
    swe.PLUTO,  # 9
    swe.MEAN_NODE,  # 10
    swe.TRUE_NODE,  # 11
]

HOUSE_SYSTEMS = [ord("P"), ord("K"), ord("E"), ord("W")]

LOCATIONS = [
    (12.5, 41.9, "Rome"),
    (139.7, 35.7, "Tokyo"),
    (-74.0, 40.7, "NewYork"),
    (0.0, 0.0, "Equator"),
]


def safe_float(v: float) -> float:
    """Convert to native float, handle non-finite."""
    v = float(v)
    if not math.isfinite(v):
        return 0.0
    return v


def generate_calc_ut_entries() -> list[dict]:
    """Generate calc_ut golden entries: bodies × dates × flag combos."""
    entries = []
    flag_combos = [
        (0, "default"),
        (swe.FLG_SPEED, "speed"),
        (swe.FLG_EQUATORIAL, "equatorial"),
        (swe.FLG_HELCTR, "heliocentric"),
    ]

    # 12 bodies × 8 dates × 1 flag = 96 entries (default flags only for all)
    for body in BODIES:
        for jd in JDS:
            pos, retflag = swe.calc_ut(jd, body, swe.FLG_SPEED)
            entries.append(
                {
                    "type": "calc_ut",
                    "jd": jd,
                    "body": body,
                    "flags": swe.FLG_SPEED,
                    "result": [safe_float(v) for v in pos],
                    "retflag": int(retflag),
                }
            )

    # Additional flag combos for Sun and Moon only (to keep count manageable)
    for body in [swe.SUN, swe.MOON]:
        jd = 2451545.0  # J2000
        for flags, desc in flag_combos:
            if flags == swe.FLG_SPEED:
                continue  # Already covered above
            try:
                pos, retflag = swe.calc_ut(jd, body, flags)
                entries.append(
                    {
                        "type": "calc_ut",
                        "jd": jd,
                        "body": body,
                        "flags": flags,
                        "flags_desc": desc,
                        "result": [safe_float(v) for v in pos],
                        "retflag": int(retflag),
                    }
                )
            except Exception:
                pass

    return entries


def generate_houses_entries() -> list[dict]:
    """Generate houses golden entries."""
    entries = []
    jd = 2451545.0

    for lon, lat, loc_name in LOCATIONS:
        for hsys in HOUSE_SYSTEMS:
            cusps, angles = swe.houses(jd, lat, lon, hsys)
            entries.append(
                {
                    "type": "houses",
                    "jd": jd,
                    "lat": lat,
                    "lon": lon,
                    "location": loc_name,
                    "hsys": chr(hsys),
                    "cusps": [safe_float(v) for v in cusps],
                    "angles": [safe_float(v) for v in angles],
                }
            )

    return entries


def generate_sidereal_entries() -> list[dict]:
    """Generate sidereal position golden entries."""
    entries = []
    jd = 2451545.0
    modes = [
        (swe.SIDM_LAHIRI, "Lahiri"),
        (swe.SIDM_FAGAN_BRADLEY, "FaganBradley"),
    ]

    for mode, mode_name in modes:
        swe.set_sid_mode(mode)
        for body in [swe.SUN, swe.MOON, swe.MARS]:
            pos, retflag = swe.calc_ut(jd, body, swe.FLG_SIDEREAL | swe.FLG_SPEED)
            entries.append(
                {
                    "type": "sidereal",
                    "jd": jd,
                    "body": body,
                    "sid_mode": mode,
                    "sid_mode_name": mode_name,
                    "result": [safe_float(v) for v in pos],
                }
            )

    # Reset to default
    swe.set_sid_mode(swe.SIDM_LAHIRI)
    return entries


def generate_time_entries() -> list[dict]:
    """Generate time conversion golden entries."""
    entries = []

    # julday / revjul roundtrips
    dates = [
        (2000, 1, 1, 12.0),
        (1900, 6, 15, 6.5),
        (2050, 12, 31, 23.99),
    ]
    for y, m, d, h in dates:
        jd = swe.julday(y, m, d, h)
        yr, mr, dr, hr = swe.revjul(jd)
        entries.append(
            {
                "type": "julday",
                "input": [y, m, d, h],
                "jd": safe_float(jd),
                "revjul": [int(yr), int(mr), int(dr), safe_float(hr)],
            }
        )

    # sidtime
    for jd in [2451545.0, 2460676.5]:
        st = swe.sidtime(jd)
        entries.append(
            {
                "type": "sidtime",
                "jd": jd,
                "result": safe_float(st),
            }
        )

    # deltat
    for jd in [2451545.0, 2460676.5]:
        dt = swe.deltat(jd)
        entries.append(
            {
                "type": "deltat",
                "jd": jd,
                "result": safe_float(dt),
            }
        )

    return entries


def generate_eclipse_entries() -> list[dict]:
    """Generate eclipse golden entries."""
    entries = []

    # Solar eclipse
    jd = swe.julday(2024, 4, 1, 0.0)
    ecl_type, times = swe.sol_eclipse_when_glob(jd, ecltype=swe.ECL_TOTAL)
    entries.append(
        {
            "type": "solar_eclipse",
            "search_jd": jd,
            "ecl_type": int(ecl_type),
            "times": [safe_float(t) for t in times],
        }
    )

    # Lunar eclipse
    jd = swe.julday(2025, 3, 1, 0.0)
    ecl_type, times = swe.lun_eclipse_when(jd, ecltype=swe.ECL_TOTAL)
    entries.append(
        {
            "type": "lunar_eclipse",
            "search_jd": jd,
            "ecl_type": int(ecl_type),
            "times": [safe_float(t) for t in times],
        }
    )

    return entries


def main() -> None:
    """Generate the golden reference file."""
    print("Generating golden reference file...")
    start = time.monotonic()

    all_entries: list[dict] = []
    all_entries.extend(generate_calc_ut_entries())
    all_entries.extend(generate_houses_entries())
    all_entries.extend(generate_sidereal_entries())
    all_entries.extend(generate_time_entries())
    all_entries.extend(generate_eclipse_entries())

    elapsed = time.monotonic() - start

    golden = {
        "version": 1,
        "generator": "scripts/generate_golden.py",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "entry_count": len(all_entries),
        "entries": all_entries,
    }

    with open(OUTPUT_FILE, "w") as f:
        json.dump(golden, f, indent=2)

    print(f"Generated {len(all_entries)} golden entries in {elapsed:.2f}s")
    print(f"Saved to {OUTPUT_FILE}")

    # Summary by type
    from collections import Counter

    counts = Counter(e["type"] for e in all_entries)
    for t, c in sorted(counts.items()):
        print(f"  {t}: {c}")


if __name__ == "__main__":
    main()
