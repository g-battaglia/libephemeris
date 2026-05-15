"""
Tests for True Lunar Node velocity calculation.

Validates that velocity (FLG_SPEED) is calculated correctly via numerical
differentiation.

Velocity precision characteristics:
- Mean Node: < 0.001 degrees/day (excellent - proves numerical diff is correct)
- True Node at J2000: < 0.001 degrees/day (excellent)
- True Node 1950-2050: < 0.01 degrees/day (limited by position model differences)

The True Node position calculation has known differences from pyswisseph
(up to 0.15° for dates far from J2000 - see test_true_node_precision.py).
Since velocity is derived via numerical differentiation, velocity accuracy
is inherently bounded by position accuracy.
"""

import random
import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import TRUE_NODE, MEAN_NODE, FLG_SPEED


# Mean Node velocity tolerance (should be very precise)
MEAN_NODE_VELOCITY_TOLERANCE = 0.001  # degrees/day

# True Node velocity tolerance at J2000 (where position is most accurate)
TRUE_NODE_J2000_VELOCITY_TOLERANCE = 0.001  # degrees/day

# True Node velocity tolerance for general dates (limited by position accuracy)
# Position can differ by up to 0.15°, so velocity inherits this imprecision
TRUE_NODE_VELOCITY_TOLERANCE = 0.01  # degrees/day


class TestTrueNodeVelocity:
    """Test True Node velocity calculations."""

    @pytest.mark.unit
    def test_true_node_velocity_is_negative(self):
        """True Node moves retrograde (west) at about -0.05 degrees/day."""
        jd = 2451545.0  # J2000.0
        pos, _ = ephem.calc_ut(jd, TRUE_NODE, FLG_SPEED)

        # True node moves retrograde, so velocity should be negative
        # Typical velocity is around -0.05 degrees/day (completes cycle in ~18.6 years)
        assert pos[3] < 0, f"True node should move retrograde, got velocity {pos[3]}"
        assert pos[3] > -0.1, f"True node velocity {pos[3]} too fast (expected ~-0.05)"

    @pytest.mark.unit
    def test_true_node_velocity_vs_position_change(self):
        """True node velocity should match actual position change."""
        jd = 2451545.0
        dt = 1.0  # 1 day

        pos1, _ = ephem.calc_ut(jd, TRUE_NODE, FLG_SPEED)
        pos2, _ = ephem.calc_ut(jd + dt, TRUE_NODE, FLG_SPEED)

        # Actual position change
        actual_change = pos2[0] - pos1[0]
        if actual_change > 180:
            actual_change -= 360
        elif actual_change < -180:
            actual_change += 360

        # Velocity prediction
        predicted_change = pos1[3] * dt

        # Should match within 0.01 degrees (node has slow velocity variations)
        assert abs(actual_change - predicted_change) < 0.01, (
            f"Velocity {pos1[3]} doesn't match change {actual_change}"
        )

    @pytest.mark.comparison
    def test_true_node_velocity_vs_pyswisseph(self):
        """True node velocity should match pyswisseph within 0.001 degrees/day at J2000."""
        jd = 2451545.0  # J2000.0

        pos_lib, _ = ephem.calc_ut(jd, TRUE_NODE, FLG_SPEED)
        pos_swe, _ = swe.calc_ut(jd, swe.TRUE_NODE, swe.FLG_SPEED)

        vel_diff = abs(pos_lib[3] - pos_swe[3])
        assert vel_diff < TRUE_NODE_J2000_VELOCITY_TOLERANCE, (
            f"True node velocity diff {vel_diff} >= {TRUE_NODE_J2000_VELOCITY_TOLERANCE}\n"
            f"  libephemeris: {pos_lib[3]}\n"
            f"  pyswisseph: {pos_swe[3]}"
        )

    @pytest.mark.comparison
    def test_true_node_velocity_multiple_dates(self, progress_reporter):
        """Test true node velocity against pyswisseph for 100 random dates."""
        random.seed(42)
        dates = []
        for _ in range(100):
            year = random.randint(1950, 2050)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = ephem.julday(year, month, day, hour)
            dates.append((year, jd))

        progress = progress_reporter("True node velocity", len(dates), report_every=10)
        max_diff = 0.0
        errors = []

        for i, (year, jd) in enumerate(dates):
            pos_lib, _ = ephem.calc_ut(jd, TRUE_NODE, FLG_SPEED)
            pos_swe, _ = swe.calc_ut(jd, swe.TRUE_NODE, swe.FLG_SPEED)

            vel_diff = abs(pos_lib[3] - pos_swe[3])
            max_diff = max(max_diff, vel_diff)

            if vel_diff >= TRUE_NODE_VELOCITY_TOLERANCE:
                errors.append(
                    f"Year {year}: diff={vel_diff:.6f} "
                    f"(lib={pos_lib[3]:.6f}, swe={pos_swe[3]:.6f})"
                )
            progress.update(i)

        progress.done(f"max diff: {max_diff:.6f}")

        assert not errors, f"Velocity mismatch in {len(errors)} dates:\n" + "\n".join(
            errors[:10]
        )

    @pytest.mark.unit
    def test_true_node_velocity_without_flag_is_zero(self):
        """Without FLG_SPEED, velocity should be zero."""
        jd = 2451545.0
        pos, _ = ephem.calc_ut(jd, TRUE_NODE, 0)  # No FLG_SPEED

        assert pos[3] == 0.0, f"Velocity should be 0 without FLG_SPEED, got {pos[3]}"
        assert pos[4] == 0.0, (
            f"Lat velocity should be 0 without FLG_SPEED, got {pos[4]}"
        )
        assert pos[5] == 0.0, (
            f"Dist velocity should be 0 without FLG_SPEED, got {pos[5]}"
        )


class TestMeanNodeVelocity:
    """Test Mean Node velocity calculations."""

    @pytest.mark.unit
    def test_mean_node_velocity_is_negative(self):
        """Mean Node moves retrograde at about -0.053 degrees/day."""
        jd = 2451545.0
        pos, _ = ephem.calc_ut(jd, MEAN_NODE, FLG_SPEED)

        # Mean node moves retrograde with nearly constant velocity
        assert pos[3] < 0, f"Mean node should move retrograde, got velocity {pos[3]}"
        assert -0.06 < pos[3] < -0.04, (
            f"Mean node velocity {pos[3]} out of expected range"
        )

    @pytest.mark.comparison
    def test_mean_node_velocity_vs_pyswisseph(self):
        """Mean node velocity should match pyswisseph within tolerance."""
        jd = 2451545.0

        pos_lib, _ = ephem.calc_ut(jd, MEAN_NODE, FLG_SPEED)
        pos_swe, _ = swe.calc_ut(jd, swe.MEAN_NODE, swe.FLG_SPEED)

        vel_diff = abs(pos_lib[3] - pos_swe[3])
        assert vel_diff < MEAN_NODE_VELOCITY_TOLERANCE, (
            f"Mean node velocity diff {vel_diff} >= {MEAN_NODE_VELOCITY_TOLERANCE}\n"
            f"  libephemeris: {pos_lib[3]}\n"
            f"  pyswisseph: {pos_swe[3]}"
        )

    @pytest.mark.comparison
    def test_mean_node_velocity_multiple_dates(self, progress_reporter):
        """Test mean node velocity against pyswisseph for 100 random dates."""
        random.seed(42)
        dates = []
        for _ in range(100):
            year = random.randint(1950, 2050)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = ephem.julday(year, month, day, hour)
            dates.append((year, jd))

        progress = progress_reporter("Mean node velocity", len(dates), report_every=10)
        max_diff = 0.0
        errors = []

        for i, (year, jd) in enumerate(dates):
            pos_lib, _ = ephem.calc_ut(jd, MEAN_NODE, FLG_SPEED)
            pos_swe, _ = swe.calc_ut(jd, swe.MEAN_NODE, swe.FLG_SPEED)

            vel_diff = abs(pos_lib[3] - pos_swe[3])
            max_diff = max(max_diff, vel_diff)

            if vel_diff >= MEAN_NODE_VELOCITY_TOLERANCE:
                errors.append(
                    f"Year {year}: diff={vel_diff:.6f} "
                    f"(lib={pos_lib[3]:.6f}, swe={pos_swe[3]:.6f})"
                )
            progress.update(i)

        progress.done(f"max diff: {max_diff:.6f}")

        assert not errors, f"Velocity mismatch in {len(errors)} dates:\n" + "\n".join(
            errors[:10]
        )


class TestSouthNodeVelocity:
    """Test that South Node velocity is correctly derived from North Node."""

    @pytest.mark.unit
    def test_south_node_velocity_equals_north_node(self):
        """South node longitude velocity should equal north node velocity."""
        jd = 2451545.0

        # Using negative body ID for south node
        pos_north, _ = ephem.calc_ut(jd, TRUE_NODE, FLG_SPEED)
        pos_south, _ = ephem.calc_ut(jd, -TRUE_NODE, FLG_SPEED)

        # Longitude velocity should be identical
        assert pos_south[3] == pos_north[3], (
            f"South/North velocity mismatch: {pos_south[3]} != {pos_north[3]}"
        )

        # Latitude velocity should be inverted
        assert pos_south[4] == -pos_north[4], (
            f"South lat velocity should be -North: {pos_south[4]} != {-pos_north[4]}"
        )
