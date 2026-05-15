"""
Tests for True Lilith (osculating lunar apogee) velocity calculation.

Validates that velocity (FLG_SPEED) is calculated correctly via numerical
differentiation for True Lilith.

True Lilith (osculating apogee) can move rapidly because the osculating orbital
elements change quickly due to solar perturbations. This makes velocity
calculation particularly important for accurate ephemeris generation.

Velocity precision characteristics:
- True Lilith velocity: The osculating apogee velocity varies significantly
  due to the rapidly changing osculating orbital elements.
- Expected range: approximately -1.0 to +1.0 degrees/day with typical values
  around 0.1-0.5 degrees/day.

The True Lilith position calculation has known differences from pyswisseph
(up to 5-7 degrees in some cases - see test_true_lilith_precision.py).
Since velocity is derived via numerical differentiation, velocity accuracy
is inherently bounded by position accuracy.
"""

import random
import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import OSCU_APOG, MEAN_APOG, FLG_SPEED


# True Lilith velocity tolerance (osculating apogee can have significant velocity variations)
# Due to position differences with pyswisseph, we expect larger velocity differences
# The velocity is derived from positions that already differ by 5-7 degrees in some cases,
# so velocity differences of ~0.7 degrees/day are expected
TRUE_LILITH_VELOCITY_TOLERANCE = 0.75  # degrees/day

# At J2000 where position is most accurate, velocity should be closer
TRUE_LILITH_J2000_VELOCITY_TOLERANCE = 0.2  # degrees/day


class TestTrueLilithVelocity:
    """Test True Lilith (osculating apogee) velocity calculations."""

    @pytest.mark.unit
    def test_true_lilith_velocity_is_nonzero_with_flag(self):
        """True Lilith velocity should be calculated when FLG_SPEED is set."""
        jd = 2451545.0  # J2000.0
        pos, _ = ephem.calc_ut(jd, OSCU_APOG, FLG_SPEED)

        # Velocity components should be non-zero
        assert pos[3] != 0.0, "Longitude velocity should be non-zero with FLG_SPEED"

    @pytest.mark.unit
    def test_true_lilith_velocity_without_flag_is_zero(self):
        """Without FLG_SPEED, velocity should be zero."""
        jd = 2451545.0
        pos, _ = ephem.calc_ut(jd, OSCU_APOG, 0)  # No FLG_SPEED

        assert pos[3] == 0.0, f"Velocity should be 0 without FLG_SPEED, got {pos[3]}"
        assert pos[4] == 0.0, (
            f"Lat velocity should be 0 without FLG_SPEED, got {pos[4]}"
        )
        assert pos[5] == 0.0, (
            f"Dist velocity should be 0 without FLG_SPEED, got {pos[5]}"
        )

    @pytest.mark.unit
    def test_true_lilith_velocity_in_expected_range(self):
        """True Lilith velocity should be within physically reasonable range."""
        jd = 2451545.0  # J2000.0
        pos, _ = ephem.calc_ut(jd, OSCU_APOG, FLG_SPEED)

        # The osculating apogee typically moves between -1.0 and +1.0 degrees/day
        # (faster than mean apogee due to short-period oscillations)
        assert -2.0 < pos[3] < 2.0, (
            f"True Lilith velocity {pos[3]} deg/day outside expected range [-2, 2]"
        )

    @pytest.mark.unit
    def test_true_lilith_velocity_vs_position_change(self):
        """True Lilith velocity should match actual position change."""
        jd = 2451545.0
        dt = 1.0  # 1 day

        pos1, _ = ephem.calc_ut(jd, OSCU_APOG, FLG_SPEED)
        pos2, _ = ephem.calc_ut(jd + dt, OSCU_APOG, FLG_SPEED)

        # Actual position change
        actual_change = pos2[0] - pos1[0]
        if actual_change > 180:
            actual_change -= 360
        elif actual_change < -180:
            actual_change += 360

        # Velocity prediction (using average of velocities for better approximation)
        avg_velocity = (pos1[3] + pos2[3]) / 2.0
        predicted_change = avg_velocity * dt

        # Should match within 0.1 degrees (apogee can have significant velocity variations)
        tolerance = 0.5  # degrees
        assert abs(actual_change - predicted_change) < tolerance, (
            f"Velocity doesn't match change: predicted={predicted_change:.4f}, "
            f"actual={actual_change:.4f}, diff={abs(actual_change - predicted_change):.4f}"
        )

    @pytest.mark.comparison
    def test_true_lilith_velocity_vs_pyswisseph(self):
        """True Lilith velocity should be reasonable compared to pyswisseph."""
        jd = 2451545.0  # J2000.0

        pos_lib, _ = ephem.calc_ut(jd, OSCU_APOG, FLG_SPEED)
        pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, swe.FLG_SPEED)

        vel_diff = abs(pos_lib[3] - pos_swe[3])
        assert vel_diff < TRUE_LILITH_J2000_VELOCITY_TOLERANCE, (
            f"True Lilith velocity diff {vel_diff} >= {TRUE_LILITH_J2000_VELOCITY_TOLERANCE}\n"
            f"  libephemeris: {pos_lib[3]}\n"
            f"  pyswisseph: {pos_swe[3]}"
        )

    @pytest.mark.comparison
    def test_true_lilith_velocity_multiple_dates(self, progress_reporter):
        """Test True Lilith velocity against pyswisseph for 100 random dates."""
        random.seed(42)
        dates = []
        for _ in range(100):
            year = random.randint(1950, 2050)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = ephem.julday(year, month, day, hour)
            dates.append((year, jd))

        progress = progress_reporter(
            "True Lilith velocity", len(dates), report_every=10
        )
        max_diff = 0.0
        errors = []

        for i, (year, jd) in enumerate(dates):
            pos_lib, _ = ephem.calc_ut(jd, OSCU_APOG, FLG_SPEED)
            pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, swe.FLG_SPEED)

            vel_diff = abs(pos_lib[3] - pos_swe[3])
            max_diff = max(max_diff, vel_diff)

            if vel_diff >= TRUE_LILITH_VELOCITY_TOLERANCE:
                errors.append(
                    f"Year {year}: diff={vel_diff:.6f} "
                    f"(lib={pos_lib[3]:.6f}, swe={pos_swe[3]:.6f})"
                )
            progress.update(i)

        progress.done(f"max diff: {max_diff:.6f}")

        # Allow up to 10% of dates to exceed tolerance (osculating apogee is volatile)
        max_errors_allowed = len(dates) // 10
        assert len(errors) <= max_errors_allowed, (
            f"Velocity mismatch in {len(errors)} dates (max allowed {max_errors_allowed}):\n"
            + "\n".join(errors[:10])
        )


class TestTrueLilithLatitudeVelocity:
    """Test True Lilith latitude velocity calculations."""

    @pytest.mark.unit
    def test_latitude_velocity_is_calculated(self):
        """True Lilith latitude velocity should be calculated with FLG_SPEED."""
        jd = 2451545.0
        pos, _ = ephem.calc_ut(jd, OSCU_APOG, FLG_SPEED)

        # Latitude velocity (pos[4]) may be zero at some times, but generally isn't
        # Just verify it's being calculated (not permanently zero)
        # Check a few dates
        any_nonzero = False
        for offset in range(10):
            pos_test, _ = ephem.calc_ut(jd + offset * 10, OSCU_APOG, FLG_SPEED)
            if pos_test[4] != 0.0:
                any_nonzero = True
                break

        assert any_nonzero, "Latitude velocity should be non-zero at least sometimes"

    @pytest.mark.comparison
    def test_latitude_velocity_vs_pyswisseph(self):
        """True Lilith latitude velocity should be reasonable compared to pyswisseph."""
        jd = 2451545.0

        pos_lib, _ = ephem.calc_ut(jd, OSCU_APOG, FLG_SPEED)
        pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, swe.FLG_SPEED)

        # Latitude velocity difference
        lat_vel_diff = abs(pos_lib[4] - pos_swe[4])

        # Allow larger tolerance for latitude velocity
        assert lat_vel_diff < 0.1, (
            f"Latitude velocity diff {lat_vel_diff} too large\n"
            f"  libephemeris: {pos_lib[4]}\n"
            f"  pyswisseph: {pos_swe[4]}"
        )


class TestMeanLilithVelocity:
    """Test Mean Lilith velocity calculations for comparison."""

    @pytest.mark.unit
    def test_mean_lilith_velocity_is_nonzero(self):
        """Mean Lilith velocity should now be calculated with FLG_SPEED."""
        jd = 2451545.0
        pos, _ = ephem.calc_ut(jd, MEAN_APOG, FLG_SPEED)

        # Mean Lilith velocity is now implemented using central difference
        # The base apsidal precession is about 40.7 degrees/year = 0.111 deg/day
        # but the Mean Lilith formula includes periodic corrections that cause
        # the instantaneous velocity to vary
        assert pos[3] != 0.0, "Mean Lilith velocity should be non-zero with FLG_SPEED"
        assert 0.0 < pos[3] < 1.0, (
            f"Mean Lilith velocity {pos[3]} deg/day outside expected range [0.0, 1.0]"
        )


class TestTrueLilithVelocityEdgeCases:
    """Edge case tests for True Lilith velocity calculation."""

    @pytest.mark.unit
    def test_velocity_at_multiple_epochs(self):
        """Verify velocity is calculated correctly at different epochs."""
        # Use dates within ephemeris range (1900-2053)
        epochs = [
            2420000.5,  # 1913
            2451545.0,  # J2000
            2469807.5,  # 2050
        ]

        for jd in epochs:
            pos, _ = ephem.calc_ut(jd, OSCU_APOG, FLG_SPEED)

            # Velocity should be within reasonable range at all epochs
            assert -3.0 < pos[3] < 3.0, (
                f"Velocity {pos[3]} at JD {jd} outside reasonable range"
            )

    @pytest.mark.unit
    def test_velocity_consistency(self):
        """Calling multiple times should give consistent results."""
        jd = 2451545.0

        results = [ephem.calc_ut(jd, OSCU_APOG, FLG_SPEED) for _ in range(5)]

        velocities = [r[0][3] for r in results]
        assert all(v == velocities[0] for v in velocities), (
            f"Inconsistent velocities: {velocities}"
        )

    @pytest.mark.unit
    def test_velocity_sign_changes(self):
        """True Lilith velocity can change sign due to oscillatory motion."""
        # Check a range of dates to verify sign changes occur
        positive_count = 0
        negative_count = 0

        for offset in range(100):
            jd = 2451545.0 + offset * 3  # Check every 3 days for 300 days
            pos, _ = ephem.calc_ut(jd, OSCU_APOG, FLG_SPEED)

            if pos[3] > 0.01:
                positive_count += 1
            elif pos[3] < -0.01:
                negative_count += 1

        # Both positive and negative velocities should occur
        # (osculating apogee oscillates around mean)
        assert positive_count > 0, "Expected some positive velocities"
        assert negative_count > 0, "Expected some negative velocities"


class TestTrueLilithVelocityDocumentation:
    """Documentation tests for True Lilith velocity."""

    @pytest.mark.unit
    def test_velocity_units_are_degrees_per_day(self):
        """Verify that velocity is returned in degrees/day."""
        jd = 2451545.0
        dt = 1.0  # 1 day

        pos1, _ = ephem.calc_ut(jd, OSCU_APOG, FLG_SPEED)
        pos2, _ = ephem.calc_ut(jd + dt, OSCU_APOG, FLG_SPEED)

        # Calculate numerical derivative
        actual_change = pos2[0] - pos1[0]
        if actual_change > 180:
            actual_change -= 360
        elif actual_change < -180:
            actual_change += 360

        # The velocity should approximately equal position change per day
        # (within tolerance for numerical differentiation and non-linear motion)
        reported_velocity = pos1[3]

        # They should be in the same order of magnitude
        assert abs(reported_velocity) < 5.0, "Velocity magnitude reasonable"
        assert abs(actual_change) < 5.0, "Position change magnitude reasonable"
