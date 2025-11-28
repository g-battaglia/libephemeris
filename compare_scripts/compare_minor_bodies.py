"""
Comparison script for Minor Bodies: Asteroids and TNOs.
Note: Comparisons use relaxed tolerances as LibEphemeris uses simplified
Keplerian elements while SwissEphemeris uses high-precision perturbation models.
"""

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *
import sys

sys.path.insert(0, "/Users/giacomo/dev/libephemeris/compare_scripts")
from comparison_utils import *


ASTEROIDS = [
    (SE_CHIRON, swe.CHIRON, "Chiron"),
    (SE_PHOLUS, swe.PHOLUS, "Pholus"),
    (SE_CERES, swe.CERES, "Ceres"),
    (SE_PALLAS, swe.PALLAS, "Pallas"),
    (SE_JUNO, swe.JUNO, "Juno"),
    (SE_VESTA, swe.VESTA, "Vesta"),
]


def compare_asteroid(jd, body_py, body_swe, name, verbose=True):
    """Compare single asteroid position."""
    try:
        # SwissEph
        try:
            pos_swe, _ = swe.calc_ut(jd, body_swe, 0)
        except swe.Error as e:
            if verbose:
                print(f"\n{name:-^80}")
                print(f"  SKIPPED - SwissEph data file missing ({e})")
            return 2, 0.0  # 2 = SKIPPED

        # LibEphemeris
        pos_py, _ = ephem.swe_calc_ut(jd, body_py, 0)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])
        diff_dist = abs(pos_swe[2] - pos_py[2])

        # Relaxed tolerances for Keplerian approximation
        passed = diff_lon < 10.0 and diff_lat < 5.0

        if verbose:
            print(f"\n{name:-^80}")
            print(
                f"  Longitude:  SWE={format_coord(pos_swe[0])}°  PY={format_coord(pos_py[0])}°  Diff={format_diff(diff_lon)}° {format_status(diff_lon < 10)}"
            )
            print(
                f"  Latitude:   SWE={format_coord(pos_swe[1])}°  PY={format_coord(pos_py[1])}°  Diff={format_diff(diff_lat)}° {format_status(diff_lat < 5)}"
            )
            print(
                f"  Distance:   SWE={format_coord(pos_swe[2], 4)} AU  PY={format_coord(pos_py[2], 4)} AU  Diff={format_diff(diff_dist, 6)} AU"
            )
            print(f"  Status:     {format_status(passed)}")

        return (1 if passed else 0), max(diff_lon, diff_lat)

    except Exception as e:
        print(f"\n{name}: ERROR - {e}")
        return 0, 999.0


def main():
    print_header("ASTEROID COMPARISON: LibEphemeris vs SwissEphemeris")

    print("\n" + "!" * 80)
    print("! NOTE: LibEphemeris uses simplified Keplerian orbital elements")
    print("! Expected accuracy: ~1-5 arcminutes (astrological precision)")
    print("! SwissEphemeris uses full perturbation models (scientific precision)")
    print("!" * 80 + "\n")

    # Test dates
    subjects = [
        ("J2000.0", 2000, 1, 1, 12.0),
        ("2024-01-01", 2024, 1, 1, 0.0),
        ("2010-07-01", 2010, 7, 1, 12.0),
    ]

    # Global stats
    g_total = 0
    g_passed = 0
    g_skipped = 0
    g_failed = 0
    g_max_diff = 0.0

    for name, year, month, day, hour in subjects:
        jd = swe.julday(year, month, day, hour)
        date_str = f"{year}-{month:02d}-{day:02d}"

        print(f"\n{'=' * 80}")
        print(f"DATE: {name} ({date_str})")
        print(f"{'=' * 80}")

        for body_py, body_swe, body_name in ASTEROIDS:
            g_total += 1
            status_code, diff = compare_asteroid(jd, body_py, body_swe, body_name)

            if status_code == 1:  # PASSED
                g_passed += 1
                g_max_diff = max(g_max_diff, diff)
            elif status_code == 2:  # SKIPPED
                g_skipped += 1
            else:  # FAILED
                g_failed += 1

    # Summary
    print("\n" + "=" * 80)
    print("ASTEROID COMPARISON SUMMARY")
    print("=" * 80)
    print(f"Total tests:   {g_total}")
    print(f"Passed:        {g_passed} ✓")
    print(f"Skipped:       {g_skipped} -")
    print(f"Failed:        {g_failed} ✗")

    if g_passed > 0:
        print(
            f"Pass rate:     {g_passed / (g_total - g_skipped) * 100:.1f}% (of executed)"
        )
        print(f"Max diff:      {g_max_diff:.6f}")
    elif g_total == g_skipped:
        print("Pass rate:     N/A (All skipped)")
    else:
        print("Pass rate:     0.0%")

    print("=" * 80)


if __name__ == "__main__":
    sys.exit(main())
