# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Generate golden reference files for regression testing.

Creates a JSON file with 100 representative calculations spanning the full
API surface. This file is used by test_golden_regression.py to detect
any unexpected changes in calculation results.

Usage:
    .venv/bin/python3 scripts/generate_golden.py

Output:
    tests/golden/golden_reference.json

Provenance:
    Project-authored regression-fixture generator. It snapshots LibEphemeris's
    own reviewed, hash-pinned runtime behavior; it does not call, persist, or fit
    output from another ephemeris. Golden values detect unintended changes but
    are not scientific sources and cannot justify a production coefficient or
    algorithm.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from pathlib import Path

sys.path.insert(0, ".")
import libephemeris as ephem  # noqa: E402

# ── Configuration ─────────────────────────────────────────────────────────────

OUTPUT_FILE = "tests/golden/golden_reference.json"
REVIEWED_CORE = (
    Path(__file__).resolve().parents[1]
    / "libephemeris"
    / "data"
    / "leb2"
    / "base_core.leb2"
)


def configure_reviewed_core() -> str:
    """Pin generation exclusively to the reviewed clean-room LEB core.

    Planet-center SPKs are optional tier assets whose presence depends on the
    user's data directory.  Golden output must not change when one happens to
    be installed, so this process-local generator explicitly disables their
    lookup.  Gas-giant positions therefore come only from the pinned LEB core.
    """
    from libephemeris import state
    from libephemeris.download import DATA_FILES

    expected = DATA_FILES["base_core.leb2"]["sha256"]
    digest = hashlib.sha256(REVIEWED_CORE.read_bytes()).hexdigest()
    if digest != expected:
        raise RuntimeError(
            "Refusing to generate golden data: bundled base_core.leb2 failed "
            f"SHA-256 verification (expected {expected}, got {digest})"
        )
    ephem.set_leb_file(str(REVIEWED_CORE))
    ephem.set_calc_mode("leb")

    def no_optional_planet_center(naif_id: int, jd: float | None = None) -> None:
        del naif_id, jd
        return None

    state.get_planet_center_segment = no_optional_planet_center
    return digest


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
    ephem.SUN,  # 0
    ephem.MOON,  # 1
    ephem.MERCURY,  # 2
    ephem.VENUS,  # 3
    ephem.MARS,  # 4
    ephem.JUPITER,  # 5
    ephem.SATURN,  # 6
    ephem.URANUS,  # 7
    ephem.NEPTUNE,  # 8
    ephem.PLUTO,  # 9
    ephem.MEAN_NODE,  # 10
    ephem.TRUE_NODE,  # 11
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
        (ephem.FLG_SPEED, "speed"),
        (ephem.FLG_EQUATORIAL, "equatorial"),
        (ephem.FLG_HELCTR, "heliocentric"),
    ]

    # 12 bodies × 8 dates × 1 flag = 96 entries (default flags only for all)
    for body in BODIES:
        for jd in JDS:
            pos, retflag = ephem.calc_ut(jd, body, ephem.FLG_SPEED)
            entries.append(
                {
                    "type": "calc_ut",
                    "jd": jd,
                    "body": body,
                    "flags": ephem.FLG_SPEED,
                    "result": [safe_float(v) for v in pos],
                    "retflag": int(retflag),
                }
            )

    # Additional flag combos for Sun and Moon only (to keep count manageable)
    for body in [ephem.SUN, ephem.MOON]:
        jd = 2451545.0  # J2000
        for flags, desc in flag_combos:
            if flags == ephem.FLG_SPEED:
                continue  # Already covered above
            try:
                pos, retflag = ephem.calc_ut(jd, body, flags)
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
            cusps, angles = ephem.houses(jd, lat, lon, hsys)
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
        (ephem.SIDM_J2000, "J2000"),
        (ephem.SIDM_TRUE_CITRA, "TrueCitra"),
    ]

    for mode, mode_name in modes:
        ephem.set_sid_mode(mode)
        for body in [ephem.SUN, ephem.MOON, ephem.MARS]:
            pos, retflag = ephem.calc_ut(jd, body, ephem.FLG_SIDEREAL | ephem.FLG_SPEED)
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
    ephem.set_sid_mode(ephem.SIDM_J2000)
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
        jd = ephem.julday(y, m, d, h)
        yr, mr, dr, hr = ephem.revjul(jd)
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
        st = ephem.sidtime(jd)
        entries.append(
            {
                "type": "sidtime",
                "jd": jd,
                "result": safe_float(st),
            }
        )

    # deltat
    for jd in [2451545.0, 2460676.5]:
        dt = ephem.deltat(jd)
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
    jd = ephem.julday(2024, 4, 1, 0.0)
    ecl_type, times = ephem.sol_eclipse_when_glob(jd, ecltype=ephem.ECL_TOTAL)
    entries.append(
        {
            "type": "solar_eclipse",
            "search_jd": jd,
            "ecl_type": int(ecl_type),
            "times": [safe_float(t) for t in times],
        }
    )

    # Lunar eclipse
    jd = ephem.julday(2025, 3, 1, 0.0)
    ecl_type, times = ephem.lun_eclipse_when(jd, ecltype=ephem.ECL_TOTAL)
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
    source_sha256 = configure_reviewed_core()

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
        "provenance": (
            "LibEphemeris independent NASA JPL/IAU/IERS runtime output; "
            "no external-reference calls"
        ),
        "source_artifact": "libephemeris/data/leb2/base_core.leb2",
        "source_sha256": source_sha256,
        "optional_planet_centers": "disabled; base_core-only regression",
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
