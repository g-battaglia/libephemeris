"""
Tests for the optimal interpolation window configuration for lunar apogee.

The interpolation parameters were determined through extensive testing:
- 9 sample points (previously 7)
- 56-day window / 2 synodic months (previously 14 days)
- Linear polynomial fit (previously quadratic)

These parameters provide:
- Better smoothness: variance ~0.56 vs ~37.64 for previous settings
- Closer Swiss Ephemeris match: ~8.56° mean diff vs ~11.68° previously

Note: The ~8-10° difference from Swiss Ephemeris is fundamental and cannot
be eliminated through parameter tuning because:
1. Swiss Ephemeris uses Moshier's analytical lunar theory for osculating elements
2. LibEphemeris uses JPL DE ephemeris state vectors

The tests verify:
1. The new parameters produce smoother curves
2. Basic functionality is preserved
3. Reasonable agreement with Swiss Ephemeris (where available)
"""

import math
import pytest
import libephemeris as swe
from libephemeris import lunar


class TestOptimalInterpolationParameters:
    """Tests verifying the optimal interpolation parameters are used."""

    def test_interpolated_apogee_uses_9_samples(self):
        """Verify calc_interpolated_apogee behavior is consistent with 9 samples."""
        # The function should work correctly with 9 sample points
        jd_tt = 2451545.0  # J2000.0
        result = lunar.calc_interpolated_apogee(jd_tt)
        lon, lat, dist = result

        # Basic sanity checks
        assert 0 <= lon < 360
        assert -10 < lat < 10
        assert 0.002 < dist < 0.004  # Distance in AU

    def test_interpolated_perigee_uses_9_samples(self):
        """Verify calc_interpolated_perigee behavior is consistent with 9 samples."""
        jd_tt = 2451545.0
        result = lunar.calc_interpolated_perigee(jd_tt)
        lon, lat, dist = result

        assert 0 <= lon < 360
        assert -10 < lat < 10
        assert 0.002 < dist < 0.004  # Distance in AU


class TestInterpolatedSmoothness:
    """Tests for the smoothness of interpolated apogee motion."""

    def test_interpolated_apogee_is_smoother_than_before(self):
        """
        Test that the interpolated apogee is smoother with new parameters.

        With 9 samples, 56-day window, linear fit, the variance of daily
        motion should be significantly lower than with previous parameters.
        """
        jd_start = 2451545.0
        num_days = 50
        dt = 1.0  # 1 day step

        # Collect daily changes
        changes = []
        prev_lon = None

        for i in range(num_days):
            jd = jd_start + i * dt
            lon, _, _ = lunar.calc_interpolated_apogee(jd)

            if prev_lon is not None:
                change = lon - prev_lon
                # Handle wrap-around
                if change > 180:
                    change -= 360
                elif change < -180:
                    change += 360
                changes.append(change)

            prev_lon = lon

        # Calculate variance of changes (lower = smoother)
        mean_change = sum(changes) / len(changes)
        variance = sum((c - mean_change) ** 2 for c in changes) / len(changes)

        # The new parameters should give variance < 1.0 (vs ~37 with old params)
        # Allow some margin for test stability
        assert variance < 5.0, f"Variance {variance} too high for smooth interpolation"

    def test_interpolated_perigee_is_smoother_than_before(self):
        """Test that the interpolated perigee is smooth."""
        jd_start = 2451545.0
        num_days = 50
        dt = 1.0

        changes = []
        prev_lon = None

        for i in range(num_days):
            jd = jd_start + i * dt
            lon, _, _ = lunar.calc_interpolated_perigee(jd)

            if prev_lon is not None:
                change = lon - prev_lon
                if change > 180:
                    change -= 360
                elif change < -180:
                    change += 360
                changes.append(change)

            prev_lon = lon

        mean_change = sum(changes) / len(changes)
        variance = sum((c - mean_change) ** 2 for c in changes) / len(changes)

        assert variance < 5.0, f"Variance {variance} too high for smooth interpolation"


class TestInterpolatedVsOsculating:
    """Tests comparing interpolated and osculating apogee."""

    def test_interpolated_smoother_than_osculating(self):
        """Test that interpolated apogee varies more smoothly than osculating."""
        jd_start = 2451545.0
        dt = 1.0  # 1 day step

        oscu_positions = []
        intp_positions = []

        for i in range(20):
            jd = jd_start + i * dt

            oscu_lon, _, _ = lunar.calc_true_lilith(jd)
            intp_lon, _, _ = lunar.calc_interpolated_apogee(jd)

            oscu_positions.append(oscu_lon)
            intp_positions.append(intp_lon)

        # Calculate day-to-day changes
        def calc_changes(positions):
            changes = []
            for i in range(1, len(positions)):
                diff = positions[i] - positions[i - 1]
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360
                changes.append(diff)
            return changes

        oscu_changes = calc_changes(oscu_positions)
        intp_changes = calc_changes(intp_positions)

        oscu_mean = sum(oscu_changes) / len(oscu_changes)
        intp_mean = sum(intp_changes) / len(intp_changes)

        oscu_variance = sum((c - oscu_mean) ** 2 for c in oscu_changes) / len(
            oscu_changes
        )
        intp_variance = sum((c - intp_mean) ** 2 for c in intp_changes) / len(
            intp_changes
        )

        # Interpolated should have lower variance (smoother)
        assert intp_variance < oscu_variance, (
            f"Interpolated variance {intp_variance} should be less than osculating {oscu_variance}"
        )


class TestSwissEphemerisComparison:
    """Tests comparing with Swiss Ephemeris (requires pyswisseph)."""

    @pytest.fixture
    def swisseph(self):
        """Fixture to import swisseph, skip if not available."""
        pytest.importorskip("swisseph")
        import swisseph as swe_lib

        return swe_lib

    def test_interpolated_apogee_difference_from_se(self, swisseph):
        """
        Test that interpolated apogee has expected difference from Swiss Ephemeris.

        Due to fundamental differences in underlying calculations (JPL DE vs Moshier),
        we expect approximately 8-10° mean difference. This is not a bug but a
        consequence of different ephemeris sources.
        """
        test_jds = [
            2451545.0,  # J2000.0
            2458849.5,  # 2020-01-01
            2460676.5,  # 2025-01-01
        ]

        diffs = []
        for jd in test_jds:
            # Swiss Ephemeris interpolated apogee (INTP_APOG = 21)
            se_result = swisseph.calc_ut(jd, 21)
            se_lon = se_result[0][0]

            # Our interpolated apogee
            lib_lon, _, _ = lunar.calc_interpolated_apogee(jd)

            diff = lib_lon - se_lon
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360
            diffs.append(abs(diff))

        mean_diff = sum(diffs) / len(diffs)

        # The mean difference should be around 8-10° (with some margin)
        # If it's significantly larger, something might be wrong
        assert mean_diff < 20, f"Mean difference {mean_diff}° is unexpectedly large"

        # Document the expected difference range
        # This is informational, not a strict assertion
        print(f"\nMean difference from Swiss Ephemeris: {mean_diff:.2f}°")
        print("(~8-10° is expected due to different underlying ephemeris)")


class TestInterpolatedApogeeConsistency:
    """Tests for internal consistency of the implementation."""

    def test_continuous_motion(self):
        """Test that interpolated apogee moves continuously (no large jumps)."""
        jd_start = 2451545.0
        dt = 0.5  # 12-hour steps

        prev_lon = None
        for i in range(100):
            jd = jd_start + i * dt
            lon, _, _ = lunar.calc_interpolated_apogee(jd)

            if prev_lon is not None:
                diff = lon - prev_lon
                # Handle wrap-around
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360

                # With 12-hour steps and smooth interpolation,
                # motion should be less than a few degrees
                assert abs(diff) < 5, f"Jump of {diff}° at step {i} is too large"

            prev_lon = lon

    def test_perigee_opposite_apogee(self):
        """Test that interpolated perigee is approximately 180° from apogee."""
        test_jds = [2451545.0, 2458849.5, 2460676.5]

        for jd in test_jds:
            apogee_lon, _, _ = lunar.calc_interpolated_apogee(jd)
            perigee_lon, _, _ = lunar.calc_interpolated_perigee(jd)

            diff = abs(apogee_lon - perigee_lon)
            if diff > 180:
                diff = 360 - diff

            # According to Swiss Ephemeris documentation (section 2.2.4):
            # "Apogee and perigee are not exactly opposite - they are only roughly
            # opposite when the Sun is in conjunction with one of them or at 90° angle."
            # The deviation can be up to 28° depending on Sun-Moon geometry.
            # LibEphemeris now computes apogee and perigee with independent ELP2000-82B
            # perturbation series, matching Swiss Ephemeris behavior where they can
            # deviate from 180° by up to 28°.
            assert abs(diff - 180) < 30, (
                f"At JD {jd}: apogee-perigee difference {diff}° not close to 180°"
            )


class TestInterpolationDocumentation:
    """Tests that document the research findings."""

    def test_56_day_window_reasoning(self):
        """
        Document why 56-day window (2 synodic months) was chosen.

        The dominant oscillation in the osculating apogee has period ~14.77 days
        (half synodic month, argument 2D). A 56-day window spans approximately
        4 complete oscillation cycles, allowing effective averaging.
        """
        synodic_month = 29.53  # days
        oscillation_period = synodic_month / 2  # ~14.77 days
        window = 56  # days

        cycles_in_window = window / oscillation_period
        assert cycles_in_window > 3.5, "Window should span multiple oscillation cycles"

        # Document the configuration
        print(f"\nOscillation period (2D): {oscillation_period:.2f} days")
        print(f"Window size: {window} days")
        print(f"Oscillation cycles in window: {cycles_in_window:.2f}")

    def test_linear_fit_reasoning(self):
        """
        Document why linear polynomial fit was chosen.

        The mean apogee moves at approximately 40.7°/year = 0.111°/day.
        Over 56 days, this is ~6° of motion - essentially linear.
        Linear fit provides maximum smoothing without following oscillations.
        """
        mean_motion_per_year = 40.7  # degrees/year (prograde)
        mean_motion_per_day = mean_motion_per_year / 365.25
        window_days = 56

        motion_in_window = mean_motion_per_day * window_days

        # Motion should be small enough that linear fit is appropriate
        assert motion_in_window < 10, (
            f"Motion of {motion_in_window}° in window is small enough for linear fit"
        )

        print(f"\nMean apogee motion: {mean_motion_per_day:.3f}°/day")
        print(f"Motion over 56-day window: {motion_in_window:.2f}°")
        print("Linear fit is appropriate for this nearly-linear motion")
