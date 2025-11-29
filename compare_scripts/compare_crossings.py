"""
Comparison script for Crossing Functions.
Validates solcross_ut and mooncross_ut against SwissEphemeris.
"""

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *
import sys

sys.path.insert(0, "/Users/giacomo/dev/libephemeris/compare_scripts")
from comparison_utils import *


def compare_sun_crossing(target_lon, jd_start, test_name):
    """Compare Sun crossing calculation."""
    print(f"\n{test_name:-^80}")

    # LibEphemeris - now returns float directly, raises exception on error
    try:
        jd_cross_py = ephem.swe_solcross_ut(target_lon, jd_start, 0)
    except Exception as e:
        print(f"ERROR (LibEphemeris): {e}")
        return False, 0.0

    # SwissEphemeris
    jd_cross_swe = swe.solcross_ut(target_lon, jd_start, 0)

    # Calculate time difference in seconds
    diff_seconds = abs(jd_cross_py - jd_cross_swe) * 86400

    # Verify actual Sun position at crossing
    pos_py, _ = ephem.swe_calc_ut(jd_cross_py, SE_SUN, 0)
    pos_swe, _ = swe.calc_ut(jd_cross_swe, swe.SUN, 0)

    lon_diff_py = min(abs(pos_py[0] - target_lon), 360 - abs(pos_py[0] - target_lon))
    lon_diff_swe = min(abs(pos_swe[0] - target_lon), 360 - abs(pos_swe[0] - target_lon))

    # Get date/time for display
    y_py, m_py, d_py, h_py = ephem.swe_revjul(jd_cross_py)
    y_swe, m_swe, d_swe, h_swe = swe.revjul(jd_cross_swe)

    passed = diff_seconds < 60  # Within 1 minute

    print(f"Target Longitude: {target_lon:.2f}°")
    print("\nLibEphemeris:")
    print(f"  Crossing Time: {int(y_py)}-{int(m_py):02d}-{int(d_py):02d} {h_py:.4f}h")
    print(f"  Sun Position:  {pos_py[0]:.8f}° (error: {lon_diff_py:.8f}°)")
    print("\nSwissEphemeris:")
    print(
        f"  Crossing Time: {int(y_swe)}-{int(m_swe):02d}-{int(d_swe):02d} {h_swe:.4f}h"
    )
    print(f"  Sun Position:  {pos_swe[0]:.8f}° (error: {lon_diff_swe:.8f}°)")
    print(f"\nTime Difference: {diff_seconds:.2f} seconds {format_status(passed)}")

    return passed, diff_seconds


def compare_moon_crossing(target_lon, jd_start, test_name):
    """Compare Moon crossing calculation."""
    print(f"\n{test_name:-^80}")

    # LibEphemeris - now returns float directly, raises exception on error
    try:
        jd_cross_py = ephem.swe_mooncross_ut(target_lon, jd_start, 0)
    except Exception as e:
        print(f"ERROR (LibEphemeris): {e}")
        return False, 0.0

    # SwissEphemeris
    jd_cross_swe = swe.mooncross_ut(target_lon, jd_start, 0)

    # Calculate time difference in seconds
    diff_seconds = abs(jd_cross_py - jd_cross_swe) * 86400

    # Verify Moon position
    pos_py, _ = ephem.swe_calc_ut(jd_cross_py, SE_MOON, 0)
    pos_swe, _ = swe.calc_ut(jd_cross_swe, swe.MOON, 0)

    lon_diff_py = min(abs(pos_py[0] - target_lon), 360 - abs(pos_py[0] - target_lon))
    lon_diff_swe = min(abs(pos_swe[0] - target_lon), 360 - abs(pos_swe[0] - target_lon))

    y_py, m_py, d_py, h_py = ephem.swe_revjul(jd_cross_py)
    y_swe, m_swe, d_swe, h_swe = swe.revjul(jd_cross_swe)

    passed = diff_seconds < 120  # Within 2 minutes for Moon

    print(f"Target Longitude: {target_lon:.2f}°")
    print("\nLibEphemeris:")
    print(f"  Crossing Time: {int(y_py)}-{int(m_py):02d}-{int(d_py):02d} {h_py:.4f}h")
    print(f"  Moon Position: {pos_py[0]:.8f}° (error: {lon_diff_py:.8f}°)")
    print("\nSwissEphemeris:")
    print(
        f"  Crossing Time: {int(y_swe)}-{int(m_swe):02d}-{int(d_swe):02d} {h_swe:.4f}h"
    )
    print(f"  Moon Position: {pos_swe[0]:.8f}° (error: {lon_diff_swe:.8f}°)")
    print(f"\nTime Difference: {diff_seconds:.2f} seconds {format_status(passed)}")

    return passed, diff_seconds


def main():
    print_header("CROSSING FUNCTIONS COMPARISON: LibEphemeris vs SwissEphemeris")

    jd_start = swe.julday(2024, 1, 1, 0.0)

    stats = TestStatistics()

    # Test Sun crossings
    print(f"\n{'SUN CROSSINGS':=^80}")

    sun_targets = [0, 30, 60, 90, 120, 150, 180]
    for target in sun_targets:
        test_name = f"Sun crossing {target}° ({target // 30} signs)"
        passed, diff = compare_sun_crossing(target, jd_start, test_name)
        stats.add_result(passed, diff)
        jd_start += 35  # Move forward for next crossing

    # Test Moon crossings
    print(f"\n\n{'MOON CROSSINGS':=^80}")

    jd_start = swe.julday(2024, 11, 1, 0.0)
    moon_targets = [0, 90, 180, 270]
    for target in moon_targets:
        test_name = f"Moon crossing {target}°"
        passed, diff = compare_moon_crossing(target, jd_start, test_name)
        stats.add_result(passed, diff)
        jd_start += 7  # Move forward

    # Summary
    stats.print_summary("CROSSING FUNCTIONS SUMMARY")

    print("\n" + "=" * 80)
    print("PRECISION ASSESSMENT")
    print("=" * 80)
    print("Sun crossings:  Target < 60 seconds difference")
    print("Moon crossings: Target < 120 seconds difference")
    print("=" * 80)

    if stats.pass_rate() >= 90:
        print("\n✓ Crossing functions show excellent agreement!")
        return 0
    else:
        print("\n~ Some crossings show larger differences.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
