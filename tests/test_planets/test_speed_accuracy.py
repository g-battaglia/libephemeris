"""
Speed (velocity) accuracy tests.

Verifies that planetary speeds computed by libephemeris are consistent
with numerical differentiation, and that speed values are physically
plausible for all bodies across many dates.
"""

from __future__ import annotations

import random

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
    URANUS,
    NEPTUNE,
    PLUTO,
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    OSCU_APOG,
    CHIRON,
    FLG_SPEED,
)


def _random_jds(n: int, seed: int = 42) -> list[float]:
    """Generate n random JDs in 1900-2100."""
    rng = random.Random(seed)
    return [rng.uniform(2415020.0, 2488070.0) for _ in range(n)]


# Expected speed ranges (degrees/day) for each body
SPEED_RANGES = {
    SUN: (0.9, 1.1),
    MOON: (10.0, 16.0),
    MERCURY: (-1.5, 2.3),  # Can retrograde
    VENUS: (-0.7, 1.3),  # Can retrograde
    MARS: (-0.5, 0.8),  # Can retrograde
    JUPITER: (-0.15, 0.25),
    SATURN: (-0.12, 0.15),
    URANUS: (-0.05, 0.07),
    NEPTUNE: (-0.03, 0.04),
    PLUTO: (-0.04, 0.04),
    MEAN_NODE: (-0.06, -0.04),  # Always retrograde
    TRUE_NODE: (-0.3, 0.2),  # Oscillates
    MEAN_APOG: (0.10, 0.12),
    OSCU_APOG: (-2.0, 2.5),  # Highly variable
    CHIRON: (-0.1, 0.15),
}

CORE_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
    (PLUTO, "Pluto"),
]

ALL_BODIES = CORE_BODIES + [
    (MEAN_NODE, "MeanNode"),
    (TRUE_NODE, "TrueNode"),
    (MEAN_APOG, "MeanApog"),
    (OSCU_APOG, "OscuApog"),
    (CHIRON, "Chiron"),
]


class TestSpeedPlausibility:
    """Test that speed values are within physically plausible ranges."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", ALL_BODIES)
    def test_speed_in_range_at_j2000(self, body_id: int, body_name: str):
        """Speed at J2000 should be within expected range."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        speed = result[3]
        lo, hi = SPEED_RANGES[body_id]
        assert lo <= speed <= hi, (
            f"{body_name}: speed {speed}°/day outside [{lo}, {hi}]"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", CORE_BODIES)
    def test_speed_range_50_dates(self, body_id: int, body_name: str):
        """Speed stays in range for 50 random dates."""
        jds = _random_jds(50, seed=body_id * 11)
        lo, hi = SPEED_RANGES[body_id]
        for jd in jds:
            result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
            speed = result[3]
            assert lo <= speed <= hi, (
                f"{body_name}: speed {speed}°/day outside range at JD {jd:.1f}"
            )


class TestSpeedNumericalDerivative:
    """Test speed against numerical differentiation."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", CORE_BODIES)
    def test_speed_matches_numerical_derivative(self, body_id: int, body_name: str):
        """Analytical speed should match (pos(t+dt) - pos(t-dt)) / (2*dt)."""
        jd = 2451545.0
        dt = 0.01  # 0.01 day = ~14.4 minutes

        result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        analytical_speed = result[3]

        r_plus, _ = swe.calc_ut(jd + dt, body_id, 0)
        r_minus, _ = swe.calc_ut(jd - dt, body_id, 0)

        # Handle wrapping at 0/360
        num_diff = r_plus[0] - r_minus[0]
        if num_diff > 180:
            num_diff -= 360
        elif num_diff < -180:
            num_diff += 360
        numerical_speed = num_diff / (2 * dt)

        # Tolerance: 0.001°/day for most bodies
        tol = 0.005 if body_id == MOON else 0.001
        assert abs(analytical_speed - numerical_speed) < tol, (
            f"{body_name}: analytical {analytical_speed:.6f} vs "
            f"numerical {numerical_speed:.6f}°/day"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", CORE_BODIES)
    def test_speed_numerical_derivative_20_dates(self, body_id: int, body_name: str):
        """Speed matches numerical derivative at 20 dates per planet."""
        jds = _random_jds(20, seed=body_id * 23)
        dt = 0.01
        tol = 0.01 if body_id == MOON else 0.002

        for jd in jds:
            result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
            analytical = result[3]

            rp, _ = swe.calc_ut(jd + dt, body_id, 0)
            rm, _ = swe.calc_ut(jd - dt, body_id, 0)

            num_diff = rp[0] - rm[0]
            if num_diff > 180:
                num_diff -= 360
            elif num_diff < -180:
                num_diff += 360
            numerical = num_diff / (2 * dt)

            assert abs(analytical - numerical) < tol, (
                f"{body_name} @ JD {jd:.1f}: analytical {analytical:.6f} vs "
                f"numerical {numerical:.6f}°/day"
            )


class TestLatitudeSpeed:
    """Test latitude speed is consistent with numerical derivative."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_lat_speed_numerical_derivative(self, body_id: int, body_name: str):
        """Latitude speed matches numerical derivative."""
        jd = 2451545.0
        dt = 0.01

        result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        analytical = result[4]  # lat speed

        rp, _ = swe.calc_ut(jd + dt, body_id, 0)
        rm, _ = swe.calc_ut(jd - dt, body_id, 0)
        numerical = (rp[1] - rm[1]) / (2 * dt)

        tol = 0.01
        assert abs(analytical - numerical) < tol, (
            f"{body_name}: lat speed analytical {analytical:.6f} vs "
            f"numerical {numerical:.6f}°/day"
        )


class TestDistanceSpeed:
    """Test distance speed is consistent with numerical derivative."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_dist_speed_numerical_derivative(self, body_id: int, body_name: str):
        """Distance speed matches numerical derivative."""
        jd = 2451545.0
        dt = 0.01

        result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        analytical = result[5]  # dist speed

        rp, _ = swe.calc_ut(jd + dt, body_id, 0)
        rm, _ = swe.calc_ut(jd - dt, body_id, 0)
        numerical = (rp[2] - rm[2]) / (2 * dt)

        tol = 1e-5
        assert abs(analytical - numerical) < tol, (
            f"{body_name}: dist speed analytical {analytical:.8f} vs "
            f"numerical {numerical:.8f} AU/day"
        )


class TestSpeedHighVolume:
    """High-volume speed tests — 100 dates per body."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "jd",
        _random_jds(100, seed=5555),
    )
    def test_sun_speed_100_dates(self, jd: float):
        """Sun speed is valid at 100 random dates."""
        result, _ = swe.calc_ut(jd, SUN, FLG_SPEED)
        assert 0.9 < result[3] < 1.1, f"Sun speed {result[3]} at JD {jd}"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "jd",
        _random_jds(100, seed=6666),
    )
    def test_moon_speed_100_dates(self, jd: float):
        """Moon speed is valid at 100 random dates."""
        result, _ = swe.calc_ut(jd, MOON, FLG_SPEED)
        assert 10.0 < result[3] < 16.0, f"Moon speed {result[3]} at JD {jd}"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
        ],
    )
    def test_inner_planet_speed_50_dates(self, body_id: int, body_name: str):
        """Inner planet speed is valid at 50 dates."""
        jds = _random_jds(50, seed=body_id * 31)
        lo, hi = SPEED_RANGES[body_id]
        for jd in jds:
            result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
            assert lo <= result[3] <= hi, f"{body_name} speed {result[3]} at JD {jd}"
