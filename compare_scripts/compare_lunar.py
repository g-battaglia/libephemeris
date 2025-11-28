"""
Comparison script for Lunar calculations: Nodes and Lilith.
Validates LibEphemeris against SwissEphemeris for lunar points.
"""

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *
import sys

sys.path.insert(0, "/Users/giacomo/dev/libephemeris/compare_scripts")
from comparison_utils import *


def compare_lunar_nodes(jd, name, date_str):
    """Compare lunar nodes between implementations."""
    print(f"\n{'=' * 80}")
    print(f"LUNAR NODES - {name} ({date_str})")
    print(f"{'=' * 80}")

    results = {}

    # Mean North Node
    print(f"\n{'Mean North Node':-^80}")
    pos_swe, _ = swe.calc_ut(jd, swe.MEAN_NODE, 0)
    pos_py, _ = ephem.swe_calc_ut(jd, SE_MEAN_NODE, 0)

    diff_lon = angular_diff(pos_swe[0], pos_py[0])
    passed = diff_lon < 0.01

    print(f"SwissEph:     Lon={format_coord(pos_swe[0], 6)}°")
    print(f"LibEphemeris: Lon={format_coord(pos_py[0], 6)}°")
    print(f"Difference:   {format_diff(diff_lon, 8)}° {format_status(passed)}")

    results["mean_node"] = (passed, diff_lon)

    # True North Node
    print(f"\n{'True North Node':-^80}")
    pos_swe, _ = swe.calc_ut(jd, swe.TRUE_NODE, 0)
    pos_py, _ = ephem.swe_calc_ut(jd, SE_TRUE_NODE, 0)

    diff_lon = angular_diff(pos_swe[0], pos_py[0])
    passed = diff_lon < 1.0  # Relaxed tolerance for true node

    print(f"SwissEph:     Lon={format_coord(pos_swe[0], 6)}°")
    print(f"LibEphemeris: Lon={format_coord(pos_py[0], 6)}°")
    print(f"Difference:   {format_diff(diff_lon, 8)}° {format_status(passed)}")

    results["true_node"] = (passed, diff_lon)

    return results


def compare_lilith(jd, name, date_str):
    """Compare Lilith between implementations."""
    print(f"\n{'=' * 80}")
    print(f"LILITH - {name} ({date_str})")
    print(f"{'=' * 80}")

    results = {}

    # Mean Lilith
    print(f"\n{'Mean Lilith (Black Moon)':-^80}")
    pos_swe, _ = swe.calc_ut(jd, swe.MEAN_APOG, 0)
    pos_py, _ = ephem.swe_calc_ut(jd, SE_MEAN_APOG, 0)

    diff_lon = angular_diff(pos_swe[0], pos_py[0])
    passed = diff_lon < 0.1

    print(f"SwissEph:     Lon={format_coord(pos_swe[0], 6)}°")
    print(f"LibEphemeris: Lon={format_coord(pos_py[0], 6)}°")
    print(f"Difference:   {format_diff(diff_lon, 8)}° {format_status(passed)}")

    results["mean_lilith"] = (passed, diff_lon)

    # True Lilith (Osculating Apogee)
    print(f"\n{'True Lilith (Osculating Apogee)':-^80}")
    pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, 0)
    pos_py, _ = ephem.swe_calc_ut(jd, SE_OSCU_APOG, 0)

    diff_lon = angular_diff(pos_swe[0], pos_py[0])
    passed = diff_lon < 5.0  # Very relaxed for osculating apogee

    print(f"SwissEph:     Lon={format_coord(pos_swe[0], 6)}°")
    print(f"LibEphemeris: Lon={format_coord(pos_py[0], 6)}°")
    print(f"Difference:   {format_diff(diff_lon, 8)}° {format_status(passed)}")
    print(
        "\nNote: Osculating apogee can have larger differences due to calculation method."
    )

    results["true_lilith"] = (passed, diff_lon)

    return results


def main():
    print_header("LUNAR CALCULATIONS COMPARISON: LibEphemeris vs SwissEphemeris")

    # Test subjects
    subjects = [
        ("J2000.0", 2000, 1, 1, 12.0),
        ("2024-01-01", 2024, 1, 1, 0.0),
        ("1990-06-15", 1990, 6, 15, 12.0),
        ("2010-12-31", 2010, 12, 31, 23.5),
    ]

    stats = TestStatistics()

    for name, year, month, day, hour in subjects:
        jd = swe.julday(year, month, day, hour)
        date_str = f"{year}-{month:02d}-{day:02d}"

        # Compare nodes
        node_results = compare_lunar_nodes(jd, name, date_str)
        for key, (passed, diff) in node_results.items():
            stats.add_result(passed, diff)

        # Compare Lilith
        lilith_results = compare_lilith(jd, name, date_str)
        for key, (passed, diff) in lilith_results.items():
            stats.add_result(passed, diff)

    # Summary
    stats.print_summary("LUNAR CALCULATIONS SUMMARY")

    if stats.pass_rate() >= 95:
        print("\n✓ Lunar calculations show excellent compatibility!")
    elif stats.pass_rate() >= 80:
        print("\n~ Lunar calculations show good compatibility with some differences.")
    else:
        print("\n✗ Lunar calculations show significant differences.")

    return 0 if stats.pass_rate() >= 80 else 1


if __name__ == "__main__":
    sys.exit(main())
