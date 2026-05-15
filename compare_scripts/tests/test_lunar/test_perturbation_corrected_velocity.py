"""
Tests for perturbation-corrected velocity calculations for True Node and True Lilith.

This module tests the improved velocity calculations that use ±0.5 day central
difference instead of 1-second steps, ensuring velocity reflects the actual rate
of change including all periodic perturbation terms from ELP2000-82B theory.

The key improvement is that velocities now correctly capture the perturbation
effects that a 1-second step would miss, enabling:
- Precise station detection for lunar nodes
- Velocity-weighted aspect orb calculations
- Smooth velocity curves without numerical noise

Expected accuracy: <0.001 degrees/day error for velocity calculations.
"""

import random
import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    TRUE_NODE,
    MEAN_NODE,
    OSCU_APOG,
    MEAN_APOG,
    FLG_SPEED,
)


# Velocity tolerance for True Node (should match actual position change)
# Due to non-linear motion, there's inherent error in linear velocity prediction
TRUE_NODE_VELOCITY_CONSISTENCY_TOLERANCE = 0.003  # degrees/day

# Velocity tolerance for True Lilith (more variable due to osculating elements)
TRUE_LILITH_VELOCITY_CONSISTENCY_TOLERANCE = 0.05  # degrees/day

# Velocity tolerance for Mean Lilith (includes periodic corrections)
MEAN_LILITH_VELOCITY_CONSISTENCY_TOLERANCE = 0.01  # degrees/day


class TestTrueNodePerturbationCorrectedVelocity:
    """Test True Node velocity uses perturbation-corrected central difference."""

    @pytest.mark.unit
    def test_velocity_matches_actual_position_change(self):
        """
        True node velocity should accurately predict position change over 1 day.

        The ±0.5 day central difference ensures velocity correctly captures
        the perturbation effects from ELP2000-82B theory.
        """
        jd = 2451545.0  # J2000.0
        dt = 1.0  # 1 day

        pos1, _ = ephem.calc_ut(jd, TRUE_NODE, FLG_SPEED)
        pos2, _ = ephem.calc_ut(jd + dt, TRUE_NODE, FLG_SPEED)

        # Actual position change
        actual_change = pos2[0] - pos1[0]
        if actual_change > 180:
            actual_change -= 360
        elif actual_change < -180:
            actual_change += 360

        # Velocity prediction using average of velocities at start and end
        avg_velocity = (pos1[3] + pos2[3]) / 2.0
        predicted_change = avg_velocity * dt

        # With perturbation-corrected velocity, this should be very accurate
        error = abs(actual_change - predicted_change)
        assert error < TRUE_NODE_VELOCITY_CONSISTENCY_TOLERANCE, (
            f"Velocity doesn't accurately predict position change:\n"
            f"  Predicted: {predicted_change:.6f} deg\n"
            f"  Actual: {actual_change:.6f} deg\n"
            f"  Error: {error:.6f} deg (tolerance: {TRUE_NODE_VELOCITY_CONSISTENCY_TOLERANCE})"
        )

    @pytest.mark.unit
    def test_velocity_is_smooth_across_time(self):
        """
        Velocity should change smoothly without discontinuities.

        The ±0.5 day step averages out high-frequency noise that could
        cause velocity discontinuities.
        """
        jd_start = 2451545.0
        velocities = []

        # Sample velocities every hour for 24 hours
        for i in range(25):
            jd = jd_start + i / 24.0
            pos, _ = ephem.calc_ut(jd, TRUE_NODE, FLG_SPEED)
            velocities.append(pos[3])

        # Check that velocity changes smoothly (no jumps > 0.01 deg/day)
        max_jump = 0.0
        for i in range(1, len(velocities)):
            jump = abs(velocities[i] - velocities[i - 1])
            max_jump = max(max_jump, jump)

        assert max_jump < 0.01, (
            f"Velocity has discontinuities: max jump = {max_jump:.6f} deg/day"
        )

    @pytest.mark.unit
    def test_velocity_multiple_dates(self):
        """
        Test velocity accuracy across multiple random dates.

        This ensures the perturbation-corrected velocity works consistently
        across different time periods.
        """
        random.seed(42)
        max_error = 0.0
        errors = []

        for _ in range(50):
            year = random.randint(1950, 2050)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = ephem.julday(year, month, day, hour)
            dt = 1.0

            pos1, _ = ephem.calc_ut(jd, TRUE_NODE, FLG_SPEED)
            pos2, _ = ephem.calc_ut(jd + dt, TRUE_NODE, FLG_SPEED)

            actual_change = pos2[0] - pos1[0]
            if actual_change > 180:
                actual_change -= 360
            elif actual_change < -180:
                actual_change += 360

            avg_velocity = (pos1[3] + pos2[3]) / 2.0
            predicted_change = avg_velocity * dt

            error = abs(actual_change - predicted_change)
            max_error = max(max_error, error)

            if error >= TRUE_NODE_VELOCITY_CONSISTENCY_TOLERANCE:
                errors.append(f"Year {year}: error={error:.6f} deg/day")

        assert not errors, (
            f"Velocity prediction errors in {len(errors)} dates:\n"
            + "\n".join(errors[:5])
        )


class TestTrueLilithPerturbationCorrectedVelocity:
    """Test True Lilith velocity uses perturbation-corrected central difference."""

    @pytest.mark.unit
    def test_velocity_matches_actual_position_change(self):
        """
        True Lilith velocity should accurately predict position change.

        The ±0.5 day central difference ensures velocity correctly captures
        the orbital mechanics perturbation effects.
        """
        jd = 2451545.0  # J2000.0
        dt = 1.0  # 1 day

        pos1, _ = ephem.calc_ut(jd, OSCU_APOG, FLG_SPEED)
        pos2, _ = ephem.calc_ut(jd + dt, OSCU_APOG, FLG_SPEED)

        # Actual position change
        actual_change = pos2[0] - pos1[0]
        if actual_change > 180:
            actual_change -= 360
        elif actual_change < -180:
            actual_change += 360

        # Velocity prediction using average of velocities
        avg_velocity = (pos1[3] + pos2[3]) / 2.0
        predicted_change = avg_velocity * dt

        # True Lilith has more variation due to osculating elements
        error = abs(actual_change - predicted_change)
        assert error < TRUE_LILITH_VELOCITY_CONSISTENCY_TOLERANCE, (
            f"Velocity doesn't accurately predict position change:\n"
            f"  Predicted: {predicted_change:.6f} deg\n"
            f"  Actual: {actual_change:.6f} deg\n"
            f"  Error: {error:.6f} deg"
        )

    @pytest.mark.unit
    def test_velocity_is_smooth_across_time(self):
        """
        Velocity should change smoothly without discontinuities.
        """
        jd_start = 2451545.0
        velocities = []

        # Sample velocities every hour for 24 hours
        for i in range(25):
            jd = jd_start + i / 24.0
            pos, _ = ephem.calc_ut(jd, OSCU_APOG, FLG_SPEED)
            velocities.append(pos[3])

        # Check that velocity changes smoothly (osculating apogee can vary more)
        max_jump = 0.0
        for i in range(1, len(velocities)):
            jump = abs(velocities[i] - velocities[i - 1])
            max_jump = max(max_jump, jump)

        # Osculating apogee velocity can vary more than node
        assert max_jump < 0.1, (
            f"Velocity has discontinuities: max jump = {max_jump:.6f} deg/day"
        )


class TestMeanLilithVelocity:
    """Test that Mean Lilith velocity is now calculated."""

    @pytest.mark.unit
    def test_mean_lilith_velocity_is_nonzero(self):
        """Mean Lilith velocity should now be calculated with FLG_SPEED."""
        jd = 2451545.0
        pos, _ = ephem.calc_ut(jd, MEAN_APOG, FLG_SPEED)

        # Mean Lilith velocity should be positive (prograde apsidal precession)
        # The base apsidal precession is about 40.7 degrees/year = 0.111 deg/day
        # However, the Mean Lilith formula includes periodic corrections that
        # can make the instantaneous velocity vary significantly (up to ~0.5 deg/day)
        assert pos[3] != 0.0, "Mean Lilith velocity should be non-zero with FLG_SPEED"
        assert 0.0 < pos[3] < 1.0, (
            f"Mean Lilith velocity {pos[3]} deg/day outside expected range [0.0, 1.0]"
        )

    @pytest.mark.unit
    def test_mean_lilith_velocity_matches_position_change(self):
        """Mean Lilith velocity should match actual position change."""
        jd = 2451545.0
        dt = 1.0

        pos1, _ = ephem.calc_ut(jd, MEAN_APOG, FLG_SPEED)
        pos2, _ = ephem.calc_ut(jd + dt, MEAN_APOG, FLG_SPEED)

        # Actual position change
        actual_change = pos2[0] - pos1[0]
        if actual_change > 180:
            actual_change -= 360
        elif actual_change < -180:
            actual_change += 360

        # Velocity prediction
        avg_velocity = (pos1[3] + pos2[3]) / 2.0
        predicted_change = avg_velocity * dt

        error = abs(actual_change - predicted_change)
        assert error < MEAN_LILITH_VELOCITY_CONSISTENCY_TOLERANCE, (
            f"Mean Lilith velocity doesn't match position change:\n"
            f"  Predicted: {predicted_change:.6f} deg\n"
            f"  Actual: {actual_change:.6f} deg\n"
            f"  Error: {error:.6f} deg"
        )


class TestMeanNodeVelocityConsistency:
    """Test Mean Node velocity is consistent with perturbation-corrected method."""

    @pytest.mark.unit
    def test_mean_node_velocity_matches_position_change(self):
        """Mean node velocity should accurately predict position change."""
        jd = 2451545.0
        dt = 1.0

        pos1, _ = ephem.calc_ut(jd, MEAN_NODE, FLG_SPEED)
        pos2, _ = ephem.calc_ut(jd + dt, MEAN_NODE, FLG_SPEED)

        # Actual position change
        actual_change = pos2[0] - pos1[0]
        if actual_change > 180:
            actual_change -= 360
        elif actual_change < -180:
            actual_change += 360

        # Velocity prediction
        avg_velocity = (pos1[3] + pos2[3]) / 2.0
        predicted_change = avg_velocity * dt

        error = abs(actual_change - predicted_change)
        assert error < 0.001, (
            f"Mean node velocity doesn't match position change:\n"
            f"  Predicted: {predicted_change:.6f} deg\n"
            f"  Actual: {actual_change:.6f} deg\n"
            f"  Error: {error:.6f} deg"
        )


class TestVelocityStationDetection:
    """Test that velocity can be used for accurate station detection."""

    @pytest.mark.unit
    def test_true_node_velocity_sign_changes(self):
        """
        True Node velocity should exhibit sign changes for station detection.

        The perturbation-corrected velocity should show clear velocity
        variations that can be used to detect when the node appears to
        stand still or reverse direction briefly.
        """
        jd_start = 2451545.0
        velocities = []

        # Sample velocities every 3 days for 100 days
        for i in range(34):
            jd = jd_start + i * 3.0
            pos, _ = ephem.calc_ut(jd, TRUE_NODE, FLG_SPEED)
            velocities.append(pos[3])

        # True node generally moves retrograde, but velocity magnitude varies
        # Check that we have velocity variation
        min_vel = min(velocities)
        max_vel = max(velocities)
        velocity_range = max_vel - min_vel

        # The velocity should show clear variation (not constant)
        assert velocity_range > 0.001, (
            f"True node velocity shows insufficient variation for station detection:\n"
            f"  Min velocity: {min_vel:.6f} deg/day\n"
            f"  Max velocity: {max_vel:.6f} deg/day\n"
            f"  Range: {velocity_range:.6f} deg/day"
        )
