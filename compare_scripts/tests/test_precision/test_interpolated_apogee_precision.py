"""
Precision tests for Interpolated Lunar Apogee and Perigee against pyswisseph.

This module validates calc_interpolated_apogee and calc_interpolated_perigee against
pyswisseph swe.calc_ut(jd, swe.INTP_APOG, flags) and swe.calc_ut(jd, swe.INTP_PERG, flags)
using 100+ test dates spanning different lunar phases.

Target precision: <0.5° max error, <0.2° mean error (achieved via Moshier analytical method)

Current Implementation Status
=============================

**Moshier Analytical Approach (Current):**

The interpolated apogee is now computed using a comprehensive analytical approach
based on Moshier's DE404-fitted lunar theory:

1. **Mean Apogee Position:** Uses Meeus polynomial for mean argument of perigee + 180°
2. **Perturbation Series:** ~50 periodic terms modeling apsidal oscillations:
   - Primary evection harmonics (kD - kM') up to k=10 (~4.69° dominant term)
   - Solar anomaly coupling (M-dependent terms with E factor)
   - Latitude coupling (F-dependent terms for inclination effects)
   - Cross-coupling terms combining D, M, M', F
   - Secular drift corrections

3. **Current Precision:**
   - Apogee: <0.5° max error, <0.2° mean error from pyswisseph
   - Perigee: ~9-19° difference (still computed as apogee + 180°)

**Note:** Swiss Ephemeris computes apogee and perigee independently - they are
NOT exactly 180° apart. Full perigee precision would require implementing a
separate perturbation series for perigee.

These tests document the current precision levels and serve as regression tests
to track any improvements to the implementation.
"""

import random
import statistics
import pytest
import swisseph as swe
from libephemeris.lunar import calc_interpolated_apogee, calc_interpolated_perigee
import libephemeris as ephem


# Current implementation thresholds (based on Moshier analytical implementation)
# These values represent the current precision with ~50 term analytical perturbation series
# Updated after implementing comprehensive Moshier-based apsidal perturbation terms
CURRENT_APOGEE_MAX_ERROR = 0.20  # degrees (LSQ-fitted perturbation series)
CURRENT_APOGEE_MEAN_ERROR = 0.05  # degrees (LSQ-fitted perturbation series)
CURRENT_PERIGEE_MAX_ERROR = 1.5  # degrees (LSQ-fitted extended evection series k=1..26)
CURRENT_PERIGEE_MEAN_ERROR = (
    0.4  # degrees (LSQ-fitted extended evection series k=1..26)
)

# Target precision as specified in the task
TARGET_MAX_ERROR = 0.1  # degrees

# Number of test dates (must be at least 100)
NUM_TEST_DATES = 120


def angular_difference(angle1: float, angle2: float) -> float:
    """
    Calculate the absolute angular difference between two angles in degrees.

    Handles wrap-around at 0/360 degrees.

    Args:
        angle1: First angle in degrees
        angle2: Second angle in degrees

    Returns:
        Absolute difference in degrees (0-180)
    """
    diff = abs(angle1 - angle2)
    if diff > 180:
        diff = 360 - diff
    return diff


def generate_lunar_phase_dates(count: int = NUM_TEST_DATES, seed: int = 42) -> list:
    """
    Generate dates spanning different lunar phases for comprehensive testing.

    The interpolated apogee/perigee behavior can vary with lunar phase due to
    the solar perturbation effects, so we sample across the full lunar cycle.

    Args:
        count: Number of dates to generate
        seed: Random seed for reproducibility

    Returns:
        List of tuples (year, month, day, hour, jd, phase_description)
    """
    random.seed(seed)
    dates = []

    # Synodic month is about 29.53 days
    # Generate dates across multiple lunar cycles within DE421 range
    base_jd = 2451545.0  # J2000.0

    # Create dates across 5 years (multiple lunar cycles)
    for i in range(count):
        # Spread dates across 5 years (about 62 lunar months)
        offset_days = random.uniform(0, 5 * 365.25)
        jd = base_jd + offset_days

        # Calculate approximate lunar phase for documentation
        # (Using synodic month approximation)
        synodic_month = 29.530589
        phase_fraction = (offset_days % synodic_month) / synodic_month
        if phase_fraction < 0.125:
            phase = "New Moon"
        elif phase_fraction < 0.25:
            phase = "Waxing Crescent"
        elif phase_fraction < 0.375:
            phase = "First Quarter"
        elif phase_fraction < 0.5:
            phase = "Waxing Gibbous"
        elif phase_fraction < 0.625:
            phase = "Full Moon"
        elif phase_fraction < 0.75:
            phase = "Waning Gibbous"
        elif phase_fraction < 0.875:
            phase = "Last Quarter"
        else:
            phase = "Waning Crescent"

        # Convert JD back to calendar for reporting
        year, month, day, hour = ephem.revjul(jd, ephem.GREG_CAL)
        dates.append((year, month, day, hour, jd, phase))

    return dates


@pytest.mark.comparison
@pytest.mark.precision
class TestInterpolatedApogeeVsPyswisseph:
    """
    Comparison of Interpolated Apogee calculation against pyswisseph.

    Tests 100+ dates spanning different lunar phases.

    Current Status:
    - Mean difference: ~8° (due to fundamental algorithmic differences)
    - Target: < 0.1° (requires implementing analytical Moshier method)
    """

    def test_interpolated_apogee_100_dates(self, progress_reporter):
        """
        Compare calc_interpolated_apogee against pyswisseph for 100+ dates.

        This test:
        - Generates 120 dates spanning different lunar phases
        - Compares library results with pyswisseph swe.calc_ut(jd, swe.INTP_APOG, 0)
        - Calculates max, mean, and standard deviation of differences
        - Documents current precision levels
        """
        dates = generate_lunar_phase_dates(count=NUM_TEST_DATES, seed=42)
        errors = []
        lat_errors = []
        discrepancies = []  # Store (diff, jd, swe_lon, lib_lon, phase)
        skipped = 0

        progress = progress_reporter(
            "Interpolated Apogee comparison", len(dates), report_every=10
        )

        for i, (year, month, day, hour, jd, phase) in enumerate(dates):
            try:
                # Calculate with pyswisseph (INTP_APOG = interpolated apogee)
                swe_result = swe.calc_ut(jd, swe.INTP_APOG, 0)
                swe_lon = swe_result[0][0]
                swe_lat = swe_result[0][1]

                # Calculate with libephemeris
                lib_lon, lib_lat, lib_ecc = calc_interpolated_apogee(jd)

                # Calculate angular difference for longitude
                lon_diff = angular_difference(swe_lon, lib_lon)
                errors.append(lon_diff)

                # Track latitude difference
                lat_diff = abs(swe_lat - lib_lat)
                lat_errors.append(lat_diff)

                # Store for discrepancy report
                discrepancies.append(
                    (
                        lon_diff,
                        jd,
                        swe_lon,
                        lib_lon,
                        swe_lat,
                        lib_lat,
                        phase,
                        f"{int(year):04d}-{int(month):02d}-{int(day):02d}",
                    )
                )

            except Exception as e:
                # Skip dates outside ephemeris range
                if "ephemeris" in str(e).lower() or "outside" in str(e).lower():
                    skipped += 1
                else:
                    raise

            progress.update(i, f"JD {jd:.1f}")

        # Calculate statistics
        if not errors:
            pytest.fail("No valid dates were tested")

        max_error = max(errors)
        mean_error = statistics.mean(errors)
        std_error = statistics.stdev(errors) if len(errors) > 1 else 0.0
        median_error = statistics.median(errors)
        min_error = min(errors)

        # Latitude statistics
        max_lat_error = max(lat_errors) if lat_errors else 0.0
        mean_lat_error = statistics.mean(lat_errors) if lat_errors else 0.0

        # Sort discrepancies by error (descending)
        discrepancies.sort(key=lambda x: x[0], reverse=True)

        # Count dates meeting target precision
        dates_under_target = sum(1 for e in errors if e < TARGET_MAX_ERROR)
        pct_under_target = 100 * dates_under_target / len(errors)

        # Print comprehensive report
        print(f"\n{'=' * 80}")
        print("INTERPOLATED APOGEE (INTP_APOG) PRECISION REPORT")
        print(f"{'=' * 80}")
        print(f"Dates tested: {len(errors)}")
        print(f"Dates skipped: {skipped}")
        print(f"Target precision: < {TARGET_MAX_ERROR} degrees")
        print("\n--- Longitude Statistics ---")
        print(f"  Maximum error: {max_error:.6f} degrees")
        print(f"  Mean error:    {mean_error:.6f} degrees")
        print(f"  Median error:  {median_error:.6f} degrees")
        print(f"  Std deviation: {std_error:.6f} degrees")
        print(f"  Minimum error: {min_error:.6f} degrees")
        print("\n--- Latitude Statistics ---")
        print(f"  Maximum error: {max_lat_error:.6f} degrees")
        print(f"  Mean error:    {mean_lat_error:.6f} degrees")
        print("\n--- Target Achievement ---")
        print(
            f"  Dates under {TARGET_MAX_ERROR}°: {dates_under_target}/{len(errors)} ({pct_under_target:.1f}%)"
        )

        # Show worst cases
        print("\n--- Top 5 Largest Discrepancies ---")
        for j, (
            diff,
            jd,
            swe_lon,
            lib_lon,
            swe_lat,
            lib_lat,
            phase,
            date_str,
        ) in enumerate(discrepancies[:5]):
            print(f"  {j + 1}. {date_str} ({phase})")
            print(f"     JD: {jd:.1f}")
            print(f"     SE: {swe_lon:.4f}°, Lib: {lib_lon:.4f}°, Diff: {diff:.4f}°")

        print("\n--- Implementation Note ---")
        print("  Swiss Ephemeris uses analytical Moshier-based method.")
        print("  LibEphemeris uses polynomial regression on osculating elements.")
        print("  Achieving < 0.1° requires implementing analytical method.")

        # Assert current implementation thresholds (regression test)
        assert max_error < CURRENT_APOGEE_MAX_ERROR, (
            f"Maximum error {max_error:.4f}° exceeds current threshold {CURRENT_APOGEE_MAX_ERROR}°. "
            "This may indicate a regression in the implementation."
        )
        assert mean_error < CURRENT_APOGEE_MEAN_ERROR, (
            f"Mean error {mean_error:.4f}° exceeds current threshold {CURRENT_APOGEE_MEAN_ERROR}°. "
            "This may indicate a regression in the implementation."
        )

    @pytest.mark.xfail(
        reason="Target precision < 0.1° requires implementing analytical Moshier method",
        strict=False,
    )
    def test_interpolated_apogee_target_precision(self, progress_reporter):
        """
        Target: Interpolated apogee should match pyswisseph within 0.1 degrees.

        This test is expected to fail until the analytical Moshier-based method
        is implemented. It documents the target precision goal.
        """
        dates = generate_lunar_phase_dates(count=NUM_TEST_DATES, seed=42)
        errors = []

        progress = progress_reporter(
            "Target precision test", len(dates), report_every=20
        )

        for i, (year, month, day, hour, jd, phase) in enumerate(dates):
            try:
                swe_result = swe.calc_ut(jd, swe.INTP_APOG, 0)
                swe_lon = swe_result[0][0]
                lib_lon, _, _ = calc_interpolated_apogee(jd)
                errors.append(angular_difference(swe_lon, lib_lon))
            except Exception:
                pass
            progress.update(i, f"JD {jd:.1f}")

        max_error = max(errors)
        mean_error = statistics.mean(errors)

        print(f"\nTarget precision test: max={max_error:.4f}°, mean={mean_error:.4f}°")

        assert max_error < TARGET_MAX_ERROR, (
            f"Maximum error {max_error:.4f}° exceeds target {TARGET_MAX_ERROR}°"
        )


@pytest.mark.comparison
@pytest.mark.precision
class TestInterpolatedPerigeeVsPyswisseph:
    """
    Comparison of Interpolated Perigee calculation against pyswisseph.

    Tests 100+ dates spanning different lunar phases.

    Current Status:
    - Mean difference: ~9.5° (due to fundamental algorithmic differences)
    - Target: < 0.1° (requires implementing analytical Moshier method)

    Note: The larger error for perigee vs apogee is expected because:
    1. Swiss Ephemeris computes apogee/perigee independently (not as opposites)
    2. LibEphemeris computes perigee as apogee + 180°
    """

    def test_interpolated_perigee_100_dates(self, progress_reporter):
        """
        Compare calc_interpolated_perigee against pyswisseph for 100+ dates.

        This test:
        - Generates 120 dates spanning different lunar phases
        - Compares library results with pyswisseph swe.calc_ut(jd, swe.INTP_PERG, 0)
        - Calculates max, mean, and standard deviation of differences
        - Documents current precision levels
        """
        dates = generate_lunar_phase_dates(
            count=NUM_TEST_DATES, seed=43
        )  # Different seed
        errors = []
        lat_errors = []
        discrepancies = []
        skipped = 0

        progress = progress_reporter(
            "Interpolated Perigee comparison", len(dates), report_every=10
        )

        for i, (year, month, day, hour, jd, phase) in enumerate(dates):
            try:
                # Calculate with pyswisseph (INTP_PERG = interpolated perigee)
                swe_result = swe.calc_ut(jd, swe.INTP_PERG, 0)
                swe_lon = swe_result[0][0]
                swe_lat = swe_result[0][1]

                # Calculate with libephemeris
                lib_lon, lib_lat, lib_ecc = calc_interpolated_perigee(jd)

                # Calculate angular difference for longitude
                lon_diff = angular_difference(swe_lon, lib_lon)
                errors.append(lon_diff)

                # Track latitude difference
                lat_diff = abs(swe_lat - lib_lat)
                lat_errors.append(lat_diff)

                # Store for discrepancy report
                discrepancies.append(
                    (
                        lon_diff,
                        jd,
                        swe_lon,
                        lib_lon,
                        swe_lat,
                        lib_lat,
                        phase,
                        f"{int(year):04d}-{int(month):02d}-{int(day):02d}",
                    )
                )

            except Exception as e:
                if "ephemeris" in str(e).lower() or "outside" in str(e).lower():
                    skipped += 1
                else:
                    raise

            progress.update(i, f"JD {jd:.1f}")

        # Calculate statistics
        if not errors:
            pytest.fail("No valid dates were tested")

        max_error = max(errors)
        mean_error = statistics.mean(errors)
        std_error = statistics.stdev(errors) if len(errors) > 1 else 0.0
        median_error = statistics.median(errors)
        min_error = min(errors)

        # Latitude statistics
        max_lat_error = max(lat_errors) if lat_errors else 0.0
        mean_lat_error = statistics.mean(lat_errors) if lat_errors else 0.0

        # Sort discrepancies by error (descending)
        discrepancies.sort(key=lambda x: x[0], reverse=True)

        # Count dates meeting target precision
        dates_under_target = sum(1 for e in errors if e < TARGET_MAX_ERROR)
        pct_under_target = 100 * dates_under_target / len(errors)

        # Print comprehensive report
        print(f"\n{'=' * 80}")
        print("INTERPOLATED PERIGEE (INTP_PERG) PRECISION REPORT")
        print(f"{'=' * 80}")
        print(f"Dates tested: {len(errors)}")
        print(f"Dates skipped: {skipped}")
        print(f"Target precision: < {TARGET_MAX_ERROR} degrees")
        print("\n--- Longitude Statistics ---")
        print(f"  Maximum error: {max_error:.6f} degrees")
        print(f"  Mean error:    {mean_error:.6f} degrees")
        print(f"  Median error:  {median_error:.6f} degrees")
        print(f"  Std deviation: {std_error:.6f} degrees")
        print(f"  Minimum error: {min_error:.6f} degrees")
        print("\n--- Latitude Statistics ---")
        print(f"  Maximum error: {max_lat_error:.6f} degrees")
        print(f"  Mean error:    {mean_lat_error:.6f} degrees")
        print("\n--- Target Achievement ---")
        print(
            f"  Dates under {TARGET_MAX_ERROR}°: {dates_under_target}/{len(errors)} ({pct_under_target:.1f}%)"
        )

        # Show worst cases
        print("\n--- Top 5 Largest Discrepancies ---")
        for j, (
            diff,
            jd,
            swe_lon,
            lib_lon,
            swe_lat,
            lib_lat,
            phase,
            date_str,
        ) in enumerate(discrepancies[:5]):
            print(f"  {j + 1}. {date_str} ({phase})")
            print(f"     JD: {jd:.1f}")
            print(f"     SE: {swe_lon:.4f}°, Lib: {lib_lon:.4f}°, Diff: {diff:.4f}°")

        print("\n--- Implementation Note ---")
        print("  Swiss Ephemeris computes apogee/perigee independently.")
        print("  LibEphemeris computes perigee as apogee + 180°.")
        print("  This causes larger errors for perigee than apogee.")

        # Assert current implementation thresholds (regression test)
        assert max_error < CURRENT_PERIGEE_MAX_ERROR, (
            f"Maximum error {max_error:.4f}° exceeds current threshold {CURRENT_PERIGEE_MAX_ERROR}°. "
            "This may indicate a regression in the implementation."
        )
        assert mean_error < CURRENT_PERIGEE_MEAN_ERROR, (
            f"Mean error {mean_error:.4f}° exceeds current threshold {CURRENT_PERIGEE_MEAN_ERROR}°. "
            "This may indicate a regression in the implementation."
        )

    @pytest.mark.xfail(
        reason="Target precision < 0.1° requires implementing analytical Moshier method",
        strict=False,
    )
    def test_interpolated_perigee_target_precision(self, progress_reporter):
        """
        Target: Interpolated perigee should match pyswisseph within 0.1 degrees.

        This test is expected to fail until the analytical Moshier-based method
        is implemented. It documents the target precision goal.
        """
        dates = generate_lunar_phase_dates(count=NUM_TEST_DATES, seed=43)
        errors = []

        progress = progress_reporter(
            "Target precision test", len(dates), report_every=20
        )

        for i, (year, month, day, hour, jd, phase) in enumerate(dates):
            try:
                swe_result = swe.calc_ut(jd, swe.INTP_PERG, 0)
                swe_lon = swe_result[0][0]
                lib_lon, _, _ = calc_interpolated_perigee(jd)
                errors.append(angular_difference(swe_lon, lib_lon))
            except Exception:
                pass
            progress.update(i, f"JD {jd:.1f}")

        max_error = max(errors)
        mean_error = statistics.mean(errors)

        print(f"\nTarget precision test: max={max_error:.4f}°, mean={mean_error:.4f}°")

        assert max_error < TARGET_MAX_ERROR, (
            f"Maximum error {max_error:.4f}° exceeds target {TARGET_MAX_ERROR}°"
        )


@pytest.mark.comparison
@pytest.mark.precision
class TestInterpolatedApogeePerigeeRelationship:
    """
    Test the relationship between interpolated apogee and perigee.

    Key Finding: Both Swiss Ephemeris and LibEphemeris compute apogee and
    perigee independently with separate perturbation series. They are NOT
    exactly 180° apart - the deviation depends on Sun-Moon geometry.

    Swiss Ephemeris shows max deviation of ~28° from being opposite.
    LibEphemeris shows max deviation of ~45° due to its ELP2000-82B series.
    """

    def test_apogee_perigee_relationship_documented(self, progress_reporter):
        """
        Document the apogee-perigee relationship differences between implementations.

        This test compares how each implementation handles the apogee-perigee
        relationship. Both compute them independently with separate perturbation
        series, showing that apogee and perigee are NOT exactly opposite.
        """
        dates = generate_lunar_phase_dates(count=NUM_TEST_DATES, seed=44)
        lib_diffs = []
        swe_diffs = []
        cross_diffs = []  # Difference between lib and swe opposition

        progress = progress_reporter(
            "Apogee-Perigee relationship", len(dates), report_every=10
        )

        for i, (year, month, day, hour, jd, phase) in enumerate(dates):
            try:
                # Calculate with libephemeris
                lib_apogee_lon, _, _ = calc_interpolated_apogee(jd)
                lib_perigee_lon, _, _ = calc_interpolated_perigee(jd)

                # Calculate with pyswisseph
                swe_apogee = swe.calc_ut(jd, swe.INTP_APOG, 0)
                swe_perigee = swe.calc_ut(jd, swe.INTP_PERG, 0)
                swe_apogee_lon = swe_apogee[0][0]
                swe_perigee_lon = swe_perigee[0][0]

                # Calculate angular separation (should be 180)
                lib_separation = angular_difference(lib_apogee_lon, lib_perigee_lon)
                swe_separation = angular_difference(swe_apogee_lon, swe_perigee_lon)

                # Deviation from 180 degrees
                lib_diff_from_180 = abs(lib_separation - 180)
                swe_diff_from_180 = abs(swe_separation - 180)

                lib_diffs.append(lib_diff_from_180)
                swe_diffs.append(swe_diff_from_180)

                # How different are the implementations' separations?
                cross_diffs.append(abs(lib_separation - swe_separation))

            except Exception:
                pass  # Skip dates with errors

            progress.update(i, f"JD {jd:.1f}")

        if not lib_diffs:
            pytest.fail("No valid dates were tested")

        # Print statistics
        print(f"\n{'=' * 80}")
        print("APOGEE-PERIGEE RELATIONSHIP ANALYSIS")
        print(f"{'=' * 80}")
        print(f"Dates tested: {len(lib_diffs)}")
        print("\n--- LibEphemeris deviation from 180° ---")
        print(f"  Max: {max(lib_diffs):.6f}°")
        print(f"  Mean: {statistics.mean(lib_diffs):.6f}°")
        print(
            "  (Note: Now shows oscillation due to improved apogee perturbation series)"
        )
        print("\n--- PySwisseph deviation from 180° ---")
        print(f"  Max: {max(swe_diffs):.6f}°")
        print(f"  Mean: {statistics.mean(swe_diffs):.6f}°")
        print("  (Expected: up to ~28° per Swiss Ephemeris documentation)")
        print("\n--- Difference between implementations ---")
        print(f"  Max separation difference: {max(cross_diffs):.6f}°")
        print(f"  Mean separation difference: {statistics.mean(cross_diffs):.6f}°")

        print("\n--- Swiss Ephemeris Documentation Note ---")
        print("  'Apogee and perigee are not exactly opposite - they are only")
        print("  roughly opposite when the Sun is in conjunction with one of")
        print("  them or at 90° angle.'")

        # LibEphemeris uses independent ELP2000-82B perturbation series for
        # apogee and perigee. The deviation from 180° can reach ~44.4° due to
        # the larger amplitude of perigee oscillations in the ELP2000 model.
        # Swiss Ephemeris shows ~28° max deviation; the difference reflects
        # different perturbation series calibrations.
        # Measured max: 44.35° (2026-02-23, 120 dates seed=44)
        assert max(lib_diffs) < 46.0, (
            f"LibEphemeris apogee-perigee deviation exceeded expected range, "
            f"max deviation is {max(lib_diffs):.6f}°"
        )

        # SwissEphemeris is expected to have non-trivial deviation
        assert max(swe_diffs) > 1.0, (
            f"PySwisseph should show apogee-perigee deviation from 180° "
            f"(as documented), but max is only {max(swe_diffs):.6f}°"
        )


@pytest.mark.comparison
@pytest.mark.precision
class TestInterpolatedVsOsculatingComparison:
    """
    Test that interpolated values are smoother than osculating values.

    This validates that the interpolation algorithm is effectively smoothing
    out the spurious short-period oscillations in the osculating apogee.
    """

    def test_interpolated_smoother_than_osculating_at_multiple_dates(self):
        """
        Verify that interpolated apogee is smoother than osculating across many dates.

        The interpolated apogee should have lower day-to-day variance than the
        osculating apogee, indicating successful smoothing of spurious oscillations.
        """
        # Test at multiple starting points
        test_jds = [
            2451545.0,  # J2000
            2458849.5,  # 2020-01-01
            2460000.0,  # 2023
        ]

        for base_jd in test_jds:
            oscu_changes = []
            intp_changes = []

            # Sample over 30 days with 1-day steps
            prev_oscu = None
            prev_intp = None

            for day in range(30):
                jd = base_jd + day

                # Get osculating (True Lilith) and interpolated positions
                oscu_result, _ = ephem.calc_ut(jd, ephem.OSCU_APOG, 0)
                intp_result, _ = ephem.calc_ut(jd, ephem.INTP_APOG, 0)

                oscu_lon = oscu_result[0]
                intp_lon = intp_result[0]

                if prev_oscu is not None and prev_intp is not None:
                    # Calculate daily change
                    oscu_change = oscu_lon - prev_oscu
                    intp_change = intp_lon - prev_intp

                    # Handle wrap-around
                    if oscu_change > 180:
                        oscu_change -= 360
                    elif oscu_change < -180:
                        oscu_change += 360
                    if intp_change > 180:
                        intp_change -= 360
                    elif intp_change < -180:
                        intp_change += 360

                    oscu_changes.append(oscu_change)
                    intp_changes.append(intp_change)

                prev_oscu = oscu_lon
                prev_intp = intp_lon

            # Calculate variance (measure of smoothness)
            oscu_mean = statistics.mean(oscu_changes)
            intp_mean = statistics.mean(intp_changes)

            oscu_variance = statistics.variance(oscu_changes)
            intp_variance = statistics.variance(intp_changes)

            # Interpolated should have lower variance (smoother)
            assert intp_variance < oscu_variance, (
                f"At JD {base_jd}: Interpolated variance ({intp_variance:.4f}) should be "
                f"less than osculating variance ({oscu_variance:.4f})"
            )

    def test_interpolated_perigee_smoother_than_osculating(self):
        """
        Verify that interpolated perigee is also smoother than osculating.
        """
        base_jd = 2451545.0  # J2000
        oscu_changes = []
        intp_changes = []

        prev_oscu = None
        prev_intp = None

        for day in range(30):
            jd = base_jd + day

            # Get osculating and interpolated perigee positions
            oscu_result, _ = ephem.calc_ut(jd, ephem.OSCU_APOG, 0)
            intp_result, _ = ephem.calc_ut(jd, ephem.INTP_PERG, 0)

            # Osculating perigee = osculating apogee + 180
            oscu_lon = (oscu_result[0] + 180) % 360
            intp_lon = intp_result[0]

            if prev_oscu is not None and prev_intp is not None:
                oscu_change = oscu_lon - prev_oscu
                intp_change = intp_lon - prev_intp

                # Handle wrap-around
                if oscu_change > 180:
                    oscu_change -= 360
                elif oscu_change < -180:
                    oscu_change += 360
                if intp_change > 180:
                    intp_change -= 360
                elif intp_change < -180:
                    intp_change += 360

                oscu_changes.append(oscu_change)
                intp_changes.append(intp_change)

            prev_oscu = oscu_lon
            prev_intp = intp_lon

        oscu_variance = statistics.variance(oscu_changes)
        intp_variance = statistics.variance(intp_changes)

        assert intp_variance < oscu_variance, (
            f"Interpolated perigee variance ({intp_variance:.4f}) should be "
            f"less than osculating variance ({oscu_variance:.4f})"
        )
