#!/usr/bin/env python3
"""
Shared diagnostic logic for precision tier scripts.

This module provides the core functionality used by the per-tier diagnostic
scripts (tier_diagnostic_base.py, tier_diagnostic_medium.py,
tier_diagnostic_extended.py).

It calculates positions for all supported celestial bodies and prints a
detailed table showing ecliptic coordinates, equatorial coordinates,
velocities, and whether each body uses DE4xx, SPK, Keplerian fallback,
or analytical computation.
"""

from __future__ import annotations

import argparse
import os
import sys

# Ensure libephemeris can be imported from the project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import libephemeris as eph  # noqa: E402
from libephemeris import state  # noqa: E402
from libephemeris.constants import (  # noqa: E402
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    MEAN_NODE,
    TRUE_NODE,
    CHIRON,
    PHOLUS,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    ERIS,
    SEDNA,
    HAUMEA,
    MAKEMAKE,
    IXION,
    ORCUS,
    QUAOAR,
    VARUNA,
    NESSUS,
    ASBOLUS,
    CHARIKLO,
    GONGGONG,
    APOPHIS,
    HYGIEA,
    EROS,
    FLG_SPEED,
    FLG_EQUATORIAL,
)
from libephemeris.spk_auto import discover_local_spks  # noqa: E402
from libephemeris.state import TIERS, set_precision_tier  # noqa: E402


# =============================================================================
# BODY DEFINITIONS
# =============================================================================

# (ipl, display_name, category)
# category: "major", "lunar", "minor"
BODIES: list[tuple[int, str, str]] = [
    # Major planets
    (SUN, "Sun", "major"),
    (MOON, "Moon", "major"),
    (MERCURY, "Mercury", "major"),
    (VENUS, "Venus", "major"),
    (MARS, "Mars", "major"),
    (JUPITER, "Jupiter", "major"),
    (SATURN, "Saturn", "major"),
    (URANUS, "Uranus", "major"),
    (NEPTUNE, "Neptune", "major"),
    (PLUTO, "Pluto", "major"),
    # Lunar points
    (MEAN_NODE, "Mean Node", "lunar"),
    (TRUE_NODE, "True Node", "lunar"),
    # Minor bodies (from SPK_BODY_NAME_MAP)
    (CHIRON, "Chiron", "minor"),
    (PHOLUS, "Pholus", "minor"),
    (CERES, "Ceres", "minor"),
    (PALLAS, "Pallas", "minor"),
    (JUNO, "Juno", "minor"),
    (VESTA, "Vesta", "minor"),
    (ERIS, "Eris", "minor"),
    (SEDNA, "Sedna", "minor"),
    (HAUMEA, "Haumea", "minor"),
    (MAKEMAKE, "Makemake", "minor"),
    (IXION, "Ixion", "minor"),
    (ORCUS, "Orcus", "minor"),
    (QUAOAR, "Quaoar", "minor"),
    (VARUNA, "Varuna", "minor"),
    (NESSUS, "Nessus", "minor"),
    (ASBOLUS, "Asbolus", "minor"),
    (CHARIKLO, "Chariklo", "minor"),
    (GONGGONG, "Gonggong", "minor"),
    (APOPHIS, "Apophis", "minor"),
    (HYGIEA, "Hygiea", "minor"),
    (EROS, "Eros", "minor"),
]

MAJOR_IPLS = {
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
}
LUNAR_IPLS = {MEAN_NODE, TRUE_NODE}

# Default dates per tier (label, year, month, day, hour)
# Extended to cover full SPK date ranges with regular milestones
DEFAULT_DATES: dict[str, list[tuple[str, int, int, int, float]]] = {
    "base": [
        ("1850-01-01", 1850, 1, 1, 0.0),
        ("1875-01-01", 1875, 1, 1, 0.0),
        ("1900-01-01", 1900, 1, 1, 0.0),
        ("1925-01-01", 1925, 1, 1, 0.0),
        ("1950-01-01", 1950, 1, 1, 0.0),
        ("1975-01-01", 1975, 1, 1, 0.0),
        ("J2000", 2000, 1, 1, 12.0),
        ("2025-01-01", 2025, 1, 1, 0.0),
        ("2050-01-01", 2050, 1, 1, 0.0),
        ("2075-01-01", 2075, 1, 1, 0.0),
        ("2100-01-01", 2100, 1, 1, 0.0),
        ("2125-01-01", 2125, 1, 1, 0.0),
        ("2150-01-01", 2150, 1, 1, 0.0),
    ],
    "medium": [
        ("1550-01-01", 1550, 1, 1, 0.0),
        ("1600-01-01", 1600, 1, 1, 0.0),
        ("1650-01-01", 1650, 1, 1, 0.0),
        ("1700-01-01", 1700, 1, 1, 0.0),
        ("1750-01-01", 1750, 1, 1, 0.0),
        ("1800-01-01", 1800, 1, 1, 0.0),
        ("1850-01-01", 1850, 1, 1, 0.0),
        ("1900-01-01", 1900, 1, 1, 0.0),
        ("J2000", 2000, 1, 1, 12.0),
        ("2050-01-01", 2050, 1, 1, 0.0),
        ("2100-01-01", 2100, 1, 1, 0.0),
        ("2150-01-01", 2150, 1, 1, 0.0),
        ("2200-01-01", 2200, 1, 1, 0.0),
        ("2250-01-01", 2250, 1, 1, 0.0),
        ("2300-01-01", 2300, 1, 1, 0.0),
        ("2350-01-01", 2350, 1, 1, 0.0),
        ("2400-01-01", 2400, 1, 1, 0.0),
        ("2450-01-01", 2450, 1, 1, 0.0),
        ("2500-01-01", 2500, 1, 1, 0.0),
        ("2550-01-01", 2550, 1, 1, 0.0),
        ("2600-01-01", 2600, 1, 1, 0.0),
        ("2650-01-01", 2650, 1, 1, 0.0),
    ],
    "extended": [
        ("year -12800", -12800, 1, 1, 0.0),
        ("year -12000", -12000, 1, 1, 0.0),
        ("year -11200", -11200, 1, 1, 0.0),
        ("year -10400", -10400, 1, 1, 0.0),
        ("year -9600", -9600, 1, 1, 0.0),
        ("year -8800", -8800, 1, 1, 0.0),
        ("year -8000", -8000, 1, 1, 0.0),
        ("year -7200", -7200, 1, 1, 0.0),
        ("year -6400", -6400, 1, 1, 0.0),
        ("year -5600", -5600, 1, 1, 0.0),
        ("year -4800", -4800, 1, 1, 0.0),
        ("year -4000", -4000, 1, 1, 0.0),
        ("year -3200", -3200, 1, 1, 0.0),
        ("year -2400", -2400, 1, 1, 0.0),
        ("year -1600", -1600, 1, 1, 0.0),
        ("year -800", -800, 1, 1, 0.0),
        ("year 0", 0, 1, 1, 0.0),
        ("year 800", 800, 1, 1, 0.0),
        ("J2000", 2000, 1, 1, 12.0),
        ("year 2400", 2400, 1, 1, 0.0),
        ("year 3200", 3200, 1, 1, 0.0),
        ("year 4000", 4000, 1, 1, 0.0),
        ("year 4800", 4800, 1, 1, 0.0),
        ("year 5600", 5600, 1, 1, 0.0),
        ("year 6400", 6400, 1, 1, 0.0),
        ("year 7200", 7200, 1, 1, 0.0),
        ("year 8000", 8000, 1, 1, 0.0),
        ("year 8800", 8800, 1, 1, 0.0),
        ("year 9600", 9600, 1, 1, 0.0),
        ("year 10400", 10400, 1, 1, 0.0),
        ("year 11200", 11200, 1, 1, 0.0),
        ("year 12000", 12000, 1, 1, 0.0),
        ("year 12800", 12800, 1, 1, 0.0),
        ("year 13600", 13600, 1, 1, 0.0),
        ("year 14400", 14400, 1, 1, 0.0),
        ("year 15200", 15200, 1, 1, 0.0),
        ("year 16000", 16000, 1, 1, 0.0),
        ("year 16800", 16800, 1, 1, 0.0),
    ],
}


# =============================================================================
# FORMATTING HELPERS
# =============================================================================


def _format_deg(val: float) -> str:
    """Format degrees with 4 decimal places."""
    return f"{val:11.4f}\u00b0"


def _format_speed(val: float) -> str:
    """Format velocity in degrees/day."""
    return f"{val:+10.6f}"


def _format_dist(val: float) -> str:
    """Format distance in AU."""
    return f"{val:12.6f}"


def _format_dist_speed(val: float) -> str:
    """Format distance speed in AU/day."""
    return f"{val:+12.8f}"


def _format_hms(deg: float) -> str:
    """Convert degrees to hours:minutes:seconds string."""
    h_total = deg / 15.0
    h = int(h_total)
    m_total = (h_total - h) * 60.0
    m = int(m_total)
    s = (m_total - m) * 60.0
    return f"{h:2d}h{m:02d}m{s:05.2f}s"


def _format_dms(deg: float) -> str:
    """Convert degrees to signed degrees:arcmin:arcsec string."""
    sign = "+" if deg >= 0 else "-"
    deg = abs(deg)
    d = int(deg)
    m_total = (deg - d) * 60.0
    m = int(m_total)
    s = (m_total - m) * 60.0
    return f"{sign}{d:02d}\u00b0{m:02d}'{s:05.2f}\""


# =============================================================================
# SOURCE DETECTION
# =============================================================================


def _parse_spk_date(date_str: str) -> float:
    """Parse SPK date string 'YYYY-MM-DD' to Julian Day."""
    parts = date_str.split("-")
    year = int(parts[0])
    month = int(parts[1])
    day = int(parts[2])
    return eph.julday(year, month, day, 0.0)


def _get_source(ipl: int, tier_name: str, jd: float) -> str:
    """Determine the data source used for a body calculation.

    Returns one of: "DE440s", "DE440", "DE441", "Analytical", "SPK", "ASSIST", "Keplerian".
    """
    if ipl in MAJOR_IPLS:
        eph_file = TIERS[tier_name].ephemeris_file
        return eph_file.replace(".bsp", "").upper()
    if ipl in LUNAR_IPLS:
        return "Analytical"
    if ipl in state._SPK_BODY_MAP:
        from libephemeris.spk import get_spk_coverage

        spk_file, _ = state._SPK_BODY_MAP[ipl]
        coverage = get_spk_coverage(spk_file)
        if coverage is not None:
            spk_start, spk_end = coverage
            if spk_start <= jd <= spk_end:
                return "SPK"
    # Check if ASSIST would be used as fallback
    try:
        from libephemeris.rebound_integration import check_assist_available

        if check_assist_available():
            return "ASSIST"
    except ImportError:
        pass
    return "Keplerian"


# =============================================================================
# CALCULATION
# =============================================================================


def _calc_body(jd: float, ipl: int) -> dict:
    """Calculate all available data for a body at a given Julian Day.

    Returns dict with keys:
        ecl: (lon, lat, dist, lon_spd, lat_spd, dist_spd) or None
        equ: (ra, dec, dist, ra_spd, dec_spd, dist_spd) or None
        error: str or None
    """
    result: dict = {"ecl": None, "equ": None, "error": None}

    # Ecliptic coordinates + speed
    try:
        pos_ecl, _ = eph.calc_ut(jd, ipl, FLG_SPEED)
        result["ecl"] = pos_ecl
    except Exception as e:
        result["error"] = str(e)[:60]
        return result

    # Equatorial coordinates + speed
    try:
        pos_equ, _ = eph.calc_ut(jd, ipl, FLG_SPEED | FLG_EQUATORIAL)
        result["equ"] = pos_equ
    except Exception:
        # Equatorial may fail for some bodies; ecliptic is still valid
        pass

    return result


# =============================================================================
# TABLE OUTPUT
# =============================================================================


def _print_separator(char: str = "\u2550", width: int = 160) -> None:
    print(char * width)


def _print_header(tier_name: str) -> None:
    """Print the tier header box."""
    tier = TIERS[tier_name]
    start, end = tier.spk_date_range
    print()
    _print_separator("\u2550")
    print(
        f"  TIER: {tier_name} \u2014 {tier.ephemeris_file}"
        f" \u2014 SPK range: {start} to {end}"
    )
    print(f"  {tier.description}")
    _print_separator("\u2550")


def _print_column_headers() -> None:
    """Print the column header row."""
    print(
        f"  {'Body':<14s}"
        f" {'Lon':>11s}"
        f" {'Lat':>11s}"
        f" {'Dist(AU)':>12s}"
        f" {'LonSpd':>10s}"
        f" {'LatSpd':>10s}"
        f" {'DistSpd':>13s}"
        f"  {'RA':>12s}"
        f" {'Dec':>12s}"
        f" {'RASpd':>10s}"
        f" {'DecSpd':>10s}"
        f"  {'Source':<10s}"
    )
    _print_separator("\u2500")


def _print_group_label(label: str) -> None:
    """Print a group separator label."""
    print(f"  \u2500\u2500 {label} " + "\u2500" * (160 - len(label) - 6))


def _print_body_row(name: str, data: dict, source: str) -> None:
    """Print one body's data row."""
    if data["error"]:
        print(f"  {name:<14s}  ** {data['error']:<140s}  {source}")
        return

    ecl = data["ecl"]
    equ = data["equ"]

    lon = _format_deg(ecl[0])
    lat = _format_deg(ecl[1])
    dist = _format_dist(ecl[2])
    lon_spd = _format_speed(ecl[3])
    lat_spd = _format_speed(ecl[4])
    dist_spd = _format_dist_speed(ecl[5])

    if equ:
        ra = _format_hms(equ[0])
        dec = _format_dms(equ[1])
        ra_spd = _format_speed(equ[3])
        dec_spd = _format_speed(equ[4])
    else:
        ra = "     N/A    "
        dec = "     N/A    "
        ra_spd = "       N/A"
        dec_spd = "       N/A"

    print(
        f"  {name:<14s}"
        f" {lon}"
        f" {lat}"
        f" {dist}"
        f" {lon_spd}"
        f" {lat_spd}"
        f" {dist_spd}"
        f"  {ra}"
        f" {dec}"
        f" {ra_spd}"
        f" {dec_spd}"
        f"  {source}"
    )


def _print_summary(source_counts: dict[str, int], error_count: int) -> None:
    """Print the per-date summary line."""
    parts = []
    for src, count in sorted(source_counts.items()):
        parts.append(f"{count} {src}")
    if error_count > 0:
        parts.append(f"{error_count} Errors")
    print(f"\n  Summary: {' | '.join(parts)}")


# =============================================================================
# DATE PARSING
# =============================================================================


def _parse_date_str(date_str: str) -> tuple[str, int, int, int, float]:
    """Parse a YYYY-MM-DD string into (label, year, month, day, hour)."""
    parts = date_str.split("-")
    if len(parts) != 3:
        raise ValueError(f"Invalid date format: {date_str!r} (expected YYYY-MM-DD)")
    year = int(parts[0])
    month = int(parts[1])
    day = int(parts[2])
    return (date_str, year, month, day, 0.0)


def _parse_cli_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Tier diagnostic: calculate all bodies and show data sources.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python scripts/tier_diagnostic_base.py\n"
            "  python scripts/tier_diagnostic_base.py 1985-06-15\n"
            "  python scripts/tier_diagnostic_base.py 2000-01-01 1950-03-15\n"
            "  python scripts/tier_diagnostic_base.py --jd 2451545.0\n"
        ),
    )
    parser.add_argument(
        "dates",
        nargs="*",
        help="Dates in YYYY-MM-DD format (default: tier-specific set)",
    )
    parser.add_argument(
        "--jd",
        type=float,
        nargs="+",
        metavar="JD",
        help="Dates as Julian Day numbers",
    )
    parser.add_argument(
        "--compact",
        action="store_true",
        help="Omit equatorial coordinates (narrower table)",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download all SPK files before computing (requires network)",
    )
    return parser.parse_args()


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================


def run_diagnostic(tier_name: str) -> None:
    """Run the full tier diagnostic.

    Args:
        tier_name: One of "base", "medium", "extended".
    """
    args = _parse_cli_args()

    # Build date list
    dates: list[tuple[str, float]] = []  # (label, jd)

    if args.dates or args.jd:
        # User-provided dates
        if args.dates:
            for d in args.dates:
                info = _parse_date_str(d)
                jd = eph.julday(info[1], info[2], info[3], info[4])
                dates.append((info[0], jd))
        if args.jd:
            for jd_val in args.jd:
                dates.append((f"JD {jd_val}", jd_val))
    else:
        # Default dates for this tier
        for label, y, m, d, h in DEFAULT_DATES[tier_name]:
            jd = eph.julday(y, m, d, h)
            dates.append((label, jd))

    # Setup
    eph.close()
    set_precision_tier(tier_name)
    eph.set_strict_precision(False)

    if args.download:
        eph.set_auto_spk_download(True)
        print("  Downloading SPK files for all minor bodies...")
        eph.ensure_all_ephemerides(show_progress=True)
        print()
    else:
        eph.set_auto_spk_download(False)

    # Discover local SPK files
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    discovered = discover_local_spks(repo_root)

    # Print header
    _print_header(tier_name)

    if discovered:
        spk_bodies = [name for name, st in discovered.items() if st == "registered"]
        if spk_bodies:
            print(f"\n  Local SPK files discovered: {', '.join(spk_bodies)}")

    # Process each date
    for label, jd in dates:
        print(f"\n  Date: {label} (JD {jd:.1f})")
        _print_separator("\u2500")
        if not args.compact:
            _print_column_headers()

        source_counts: dict[str, int] = {}
        error_count = 0
        current_category = ""

        for ipl, name, category in BODIES:
            # Print group separator
            if category != current_category:
                current_category = category
                group_labels = {
                    "major": "Major planets",
                    "lunar": "Lunar points",
                    "minor": "Minor bodies",
                }
                _print_group_label(group_labels[category])

            # Calculate
            data = _calc_body(jd, ipl)
            source = _get_source(ipl, tier_name, jd)

            if data["error"]:
                error_count += 1
            else:
                source_counts[source] = source_counts.get(source, 0) + 1

            # Print row
            if args.compact:
                _print_body_row_compact(name, data, source)
            else:
                _print_body_row(name, data, source)

        _print_summary(source_counts, error_count)

    # Cleanup
    print()
    _print_separator("\u2550")
    eph.close()


def _print_body_row_compact(name: str, data: dict, source: str) -> None:
    """Print a compact body row (ecliptic only, no equatorial)."""
    if data["error"]:
        print(f"  {name:<14s}  ** {data['error']:<80s}  {source}")
        return

    ecl = data["ecl"]
    lon = _format_deg(ecl[0])
    lat = _format_deg(ecl[1])
    dist = _format_dist(ecl[2])
    lon_spd = _format_speed(ecl[3])
    lat_spd = _format_speed(ecl[4])
    dist_spd = _format_dist_speed(ecl[5])

    print(f"  {name:<14s} {lon} {lat} {dist} {lon_spd} {lat_spd} {dist_spd}  {source}")
