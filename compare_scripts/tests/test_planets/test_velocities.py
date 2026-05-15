"""
Comprehensive tests for planetary velocity calculations.

Tests cover:
- Expected velocity ranges for each planet
- Retrograde motion detection
- Velocity comparison with pyswisseph
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


class TestSunVelocity:
    """Test Sun velocity calculations."""

    @pytest.mark.unit
    def test_sun_velocity_approximately_1_degree_per_day(self):
        """Sun moves approximately 1 degree per day."""
        jd = 2451545.0
        pos, _ = ephem.calc_ut(jd, SUN, FLG_SPEED)

        # Sun velocity should be ~0.9-1.1 degrees/day
        assert 0.9 < pos[3] < 1.1, f"Sun velocity {pos[3]} should be ~1°/day"

    @pytest.mark.unit
    def test_sun_velocity_always_positive(self):
        """Sun never goes retrograde."""
        for jd in [2451545.0, 2455000.0, 2460000.0]:
            pos, _ = ephem.calc_ut(jd, SUN, FLG_SPEED)
            assert pos[3] > 0, "Sun velocity should always be positive"

    @pytest.mark.unit
    def test_sun_velocity_varies_seasonally(self):
        """Sun is faster at perihelion (January) than aphelion (July)."""
        # Perihelion around Jan 3
        jd_jan = ephem.julday(2020, 1, 3, 12.0)
        # Aphelion around Jul 4
        jd_jul = ephem.julday(2020, 7, 4, 12.0)

        pos_jan, _ = ephem.calc_ut(jd_jan, SUN, FLG_SPEED)
        pos_jul, _ = ephem.calc_ut(jd_jul, SUN, FLG_SPEED)

        assert pos_jan[3] > pos_jul[3], "Sun faster at perihelion"


class TestMoonVelocity:
    """Test Moon velocity calculations."""

    @pytest.mark.unit
    def test_moon_velocity_approximately_13_degrees_per_day(self):
        """Moon moves approximately 12-15 degrees per day."""
        jd = 2451545.0
        pos, _ = ephem.calc_ut(jd, MOON, FLG_SPEED)

        assert 11 < pos[3] < 16, f"Moon velocity {pos[3]} should be ~12-15°/day"

    @pytest.mark.unit
    def test_moon_velocity_always_positive(self):
        """Moon never goes retrograde."""
        for jd in [2451545.0, 2455000.0, 2460000.0]:
            pos, _ = ephem.calc_ut(jd, MOON, FLG_SPEED)
            assert pos[3] > 0, "Moon velocity should always be positive"


class TestPlanetVelocityRanges:
    """Test expected velocity ranges for each planet."""

    @pytest.mark.unit
    def test_mercury_velocity_range(self, progress_reporter):
        """Mercury velocity range (includes retrograde)."""
        velocities = []
        days = list(range(2451545, 2451545 + 365))
        progress = progress_reporter("Mercury velocity", len(days), report_every=25)

        for i, jd in enumerate(days):
            pos, _ = ephem.calc_ut(float(jd), MERCURY, FLG_SPEED)
            velocities.append(pos[3])
            progress.update(i)

        min_vel, max_vel = min(velocities), max(velocities)
        # Mercury can go from about -1.5 to +2.2 degrees/day
        assert -2.0 < min_vel < 0, f"Mercury min velocity {min_vel}"
        assert 1.0 < max_vel < 2.5, f"Mercury max velocity {max_vel}"
        progress.done(f"vel range: {min_vel:.2f} to {max_vel:.2f}")

    @pytest.mark.unit
    def test_venus_velocity_range(self, progress_reporter):
        """Venus velocity range (includes retrograde)."""
        velocities = []
        days = list(range(2451545, 2451545 + 500))
        progress = progress_reporter("Venus velocity", len(days), report_every=25)

        for i, jd in enumerate(days):
            pos, _ = ephem.calc_ut(float(jd), VENUS, FLG_SPEED)
            velocities.append(pos[3])
            progress.update(i)

        min_vel, max_vel = min(velocities), max(velocities)
        # Venus can go from about -0.8 to +1.3 degrees/day
        assert -1.5 < min_vel < 0, f"Venus min velocity {min_vel}"
        assert 0.5 < max_vel < 1.5, f"Venus max velocity {max_vel}"
        progress.done(f"vel range: {min_vel:.2f} to {max_vel:.2f}")

    @pytest.mark.unit
    def test_mars_velocity_range(self, progress_reporter):
        """Mars velocity range."""
        velocities = []
        days = list(range(2451545, 2451545 + 700))
        progress = progress_reporter("Mars velocity", len(days), report_every=25)

        for i, jd in enumerate(days):
            pos, _ = ephem.calc_ut(float(jd), MARS, FLG_SPEED)
            velocities.append(pos[3])
            progress.update(i)

        min_vel, max_vel = min(velocities), max(velocities)
        # Mars can go from about -0.6 to +0.8 degrees/day
        assert -1.0 < min_vel, f"Mars min velocity {min_vel}"
        assert max_vel < 1.0, f"Mars max velocity {max_vel}"
        progress.done(f"vel range: {min_vel:.2f} to {max_vel:.2f}")

    @pytest.mark.unit
    def test_outer_planets_slow(self):
        """Outer planets should move slowly."""
        jd = 2451545.0

        for planet_id, name, max_speed in [
            (JUPITER, "Jupiter", 0.25),
            (SATURN, "Saturn", 0.15),
            (URANUS, "Uranus", 0.08),
            (NEPTUNE, "Neptune", 0.05),
            (PLUTO, "Pluto", 0.05),
        ]:
            pos, _ = ephem.calc_ut(jd, planet_id, FLG_SPEED)
            assert abs(pos[3]) < max_speed, (
                f"{name} velocity {pos[3]} exceeds expected {max_speed}"
            )


class TestVelocityVsPositionDifference:
    """Test that velocity matches position change."""

    @pytest.mark.unit
    def test_sun_velocity_matches_position_change(self):
        """Sun velocity should match actual position change."""
        jd = 2451545.0
        dt = 1.0  # 1 day

        pos1, _ = ephem.calc_ut(jd, SUN, FLG_SPEED)
        pos2, _ = ephem.calc_ut(jd + dt, SUN, FLG_SPEED)

        # Actual position change
        actual_change = pos2[0] - pos1[0]
        if actual_change > 180:
            actual_change -= 360
        elif actual_change < -180:
            actual_change += 360

        # Velocity prediction
        predicted_change = pos1[3] * dt

        # Should match within ~0.01 degrees (10% of velocity)
        assert abs(actual_change - predicted_change) < 0.1, (
            f"Velocity {pos1[3]} doesn't match change {actual_change}"
        )


class TestVelocityVsPyswisseph:
    """Compare velocities with pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_velocity_matches_swisseph(self, planet_id, planet_name):
        """Velocity should match pyswisseph."""
        jd = 2451545.0
        tolerance = 0.01  # degrees/day

        pos_lib, _ = ephem.calc_ut(jd, planet_id, FLG_SPEED)
        pos_swe, _ = swe.calc_ut(jd, planet_id, FLG_SPEED)

        # Compare longitude velocity
        vel_diff = abs(pos_lib[3] - pos_swe[3])
        assert vel_diff < tolerance, (
            f"{planet_name} velocity diff {vel_diff} >= {tolerance}"
        )
