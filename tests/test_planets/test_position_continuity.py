"""Tests for position continuity — no discontinuities over time."""

from __future__ import annotations

import math
import pytest
import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    CHIRON,
    FLG_SWIEPH,
    FLG_SPEED,
)

JD_J2000 = 2451545.0


def angle_diff(a: float, b: float) -> float:
    """Minimum angular difference between two angles in degrees."""
    d = (a - b) % 360
    if d > 180:
        d -= 360
    return abs(d)


@pytest.mark.unit
class TestSunContinuity:
    """Sun position should be continuous over time."""

    def test_sun_longitude_smooth_1day(self):
        """Sun longitude changes smoothly day-to-day (~1°/day)."""
        flags = FLG_SWIEPH | FLG_SPEED
        prev, _ = swe.calc_ut(JD_J2000, SUN, flags)
        for day in range(1, 365):
            jd = JD_J2000 + day
            curr, _ = swe.calc_ut(jd, SUN, flags)
            diff = angle_diff(curr[0], prev[0])
            assert diff < 1.2, f"Sun jump {diff}° at day {day}"
            prev = curr

    def test_sun_latitude_near_zero(self):
        """Sun latitude should always be near zero (ecliptic)."""
        flags = FLG_SWIEPH | FLG_SPEED
        for day in range(0, 365, 10):
            jd = JD_J2000 + day
            result, _ = swe.calc_ut(jd, SUN, flags)
            assert abs(result[1]) < 0.01, f"Sun lat {result[1]}° at day {day}"


@pytest.mark.unit
class TestMoonContinuity:
    """Moon position should be continuous over time."""

    def test_moon_longitude_smooth_1hour(self):
        """Moon longitude changes smoothly hour-to-hour (~0.5°/hr)."""
        flags = FLG_SWIEPH | FLG_SPEED
        prev, _ = swe.calc_ut(JD_J2000, MOON, flags)
        for hour in range(1, 48):
            jd = JD_J2000 + hour / 24.0
            curr, _ = swe.calc_ut(jd, MOON, flags)
            diff = angle_diff(curr[0], prev[0])
            assert diff < 0.7, f"Moon jump {diff}° at hour {hour}"
            prev = curr

    def test_moon_latitude_range(self):
        """Moon latitude should stay within ±5.3° (orbital inclination)."""
        flags = FLG_SWIEPH | FLG_SPEED
        for day in range(0, 30):
            jd = JD_J2000 + day
            result, _ = swe.calc_ut(jd, MOON, flags)
            assert abs(result[1]) < 5.4, f"Moon lat {result[1]}° at day {day}"

    def test_moon_distance_range(self):
        """Moon distance should stay within 0.0022-0.0028 AU."""
        flags = FLG_SWIEPH | FLG_SPEED
        for day in range(0, 30):
            jd = JD_J2000 + day
            result, _ = swe.calc_ut(jd, MOON, flags)
            dist = result[2]
            assert 0.0022 < dist < 0.0028, f"Moon dist {dist} AU at day {day}"


@pytest.mark.unit
class TestPlanetContinuity:
    """Planet positions should be continuous (even through retrograde)."""

    @pytest.mark.parametrize(
        "body,name,max_daily_change",
        [
            (MERCURY, "Mercury", 2.5),
            (VENUS, "Venus", 1.5),
            (MARS, "Mars", 1.0),
            (JUPITER, "Jupiter", 0.3),
            (SATURN, "Saturn", 0.15),
        ],
    )
    def test_planet_longitude_smooth(self, body, name, max_daily_change):
        """Planet longitude changes smoothly day-to-day."""
        flags = FLG_SWIEPH | FLG_SPEED
        prev, _ = swe.calc_ut(JD_J2000, body, flags)
        for day in range(1, 90):
            jd = JD_J2000 + day
            curr, _ = swe.calc_ut(jd, body, flags)
            diff = angle_diff(curr[0], prev[0])
            assert diff < max_daily_change, f"{name} jump {diff}° at day {day}"
            prev = curr

    @pytest.mark.parametrize(
        "body,name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_planet_latitude_bounded(self, body, name):
        """Planet latitude should be within expected bounds."""
        flags = FLG_SWIEPH | FLG_SPEED
        max_lat = 0
        for day in range(0, 365, 5):
            jd = JD_J2000 + day
            result, _ = swe.calc_ut(jd, body, flags)
            max_lat = max(max_lat, abs(result[1]))
        # All planets should have latitude < 10° from ecliptic
        assert max_lat < 10, f"{name} max latitude {max_lat}°"


@pytest.mark.unit
class TestNodeContinuity:
    """Node positions should be continuous."""

    def test_mean_node_smooth(self):
        """Mean node longitude changes smoothly day-to-day."""
        flags = FLG_SWIEPH | FLG_SPEED
        prev, _ = swe.calc_ut(JD_J2000, MEAN_NODE, flags)
        for day in range(1, 365):
            jd = JD_J2000 + day
            curr, _ = swe.calc_ut(jd, MEAN_NODE, flags)
            diff = angle_diff(curr[0], prev[0])
            assert diff < 0.1, f"Mean node jump {diff}° at day {day}"
            prev = curr

    def test_true_node_smooth(self):
        """True node longitude changes smoothly day-to-day."""
        flags = FLG_SWIEPH | FLG_SPEED
        prev, _ = swe.calc_ut(JD_J2000, TRUE_NODE, flags)
        for day in range(1, 365):
            jd = JD_J2000 + day
            curr, _ = swe.calc_ut(jd, TRUE_NODE, flags)
            diff = angle_diff(curr[0], prev[0])
            # True node has short-period oscillations; allow larger change
            assert diff < 0.5, f"True node jump {diff}° at day {day}"
            prev = curr


@pytest.mark.unit
class TestSpeedMatchesNumericalDerivative:
    """Verify speed values match numerical derivatives."""

    @pytest.mark.parametrize(
        "body,name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_speed_matches_numerical(self, body, name):
        """Speed from FLG_SPEED should match (pos(t+dt) - pos(t-dt)) / (2*dt)."""
        flags = FLG_SWIEPH | FLG_SPEED
        dt = 0.01  # 0.01 days = ~14 minutes
        jd = JD_J2000
        result, _ = swe.calc_ut(jd, body, flags)
        speed_reported = result[3]  # longitude speed

        r_plus, _ = swe.calc_ut(jd + dt, body, FLG_SWIEPH)
        r_minus, _ = swe.calc_ut(jd - dt, body, FLG_SWIEPH)

        # Handle angle wrapping
        diff = r_plus[0] - r_minus[0]
        if diff > 180:
            diff -= 360
        elif diff < -180:
            diff += 360

        speed_numerical = diff / (2 * dt)
        # Tolerance depends on body — Moon has higher-order terms
        tol = 0.01 if body != MOON else 0.05
        assert speed_reported == pytest.approx(speed_numerical, abs=tol), (
            f"{name} speed {speed_reported} vs numerical {speed_numerical}"
        )

    @pytest.mark.parametrize(
        "body,name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
        ],
    )
    def test_lat_speed_matches_numerical(self, body, name):
        """Latitude speed should match numerical derivative."""
        flags = FLG_SWIEPH | FLG_SPEED
        dt = 0.01
        jd = JD_J2000
        result, _ = swe.calc_ut(jd, body, flags)
        speed_lat = result[4]

        r_plus, _ = swe.calc_ut(jd + dt, body, FLG_SWIEPH)
        r_minus, _ = swe.calc_ut(jd - dt, body, FLG_SWIEPH)

        speed_numerical = (r_plus[1] - r_minus[1]) / (2 * dt)
        tol = 0.01
        assert speed_lat == pytest.approx(speed_numerical, abs=tol), (
            f"{name} lat speed {speed_lat} vs numerical {speed_numerical}"
        )
