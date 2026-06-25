"""
Tests for speed consistency: numerical derivative vs returned speed.

Verifies that the speed values returned by calc_ut are consistent
with numerical differentiation of positions.
"""

from __future__ import annotations


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
    FLG_SWIEPH,
    FLG_SPEED,
)


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    swe.close()


JD_J2000 = 2451545.0

# Step size for numerical derivative (0.01 days = ~14.4 minutes)
DT = 0.01


def _numerical_speed(body, jd, flags, component=0):
    """Compute numerical derivative via central differences."""
    pos_m, _ = swe.calc_ut(jd - DT, body, flags)
    pos_p, _ = swe.calc_ut(jd + DT, body, flags)
    val_m = pos_m[component]
    val_p = pos_p[component]

    # Handle wrap-around for longitude (component 0)
    if component == 0:
        diff = val_p - val_m
        if diff > 180.0:
            diff -= 360.0
        elif diff < -180.0:
            diff += 360.0
        return diff / (2 * DT)

    return (val_p - val_m) / (2 * DT)


class TestSpeedConsistencyLongitude:
    """Test longitude speed vs numerical derivative."""

    PLANETS = [
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
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("body", PLANETS)
    def test_lon_speed_at_j2000(self, body):
        """Longitude speed matches numerical derivative at J2000."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        returned_speed = pos[3]  # speed in longitude

        num_speed = _numerical_speed(body, JD_J2000, flags, component=0)

        # Tolerance: 0.01 deg/day for most, looser for Moon
        tol = 0.05 if body == MOON else 0.01
        assert returned_speed == pytest.approx(num_speed, abs=tol), (
            f"Body {body}: returned {returned_speed}, numerical {num_speed}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("body", PLANETS)
    def test_lon_speed_at_2020(self, body):
        """Longitude speed matches numerical derivative at 2020."""
        flags = FLG_SWIEPH | FLG_SPEED
        jd = 2458849.5
        pos, _ = swe.calc_ut(jd, body, flags)
        returned_speed = pos[3]

        num_speed = _numerical_speed(body, jd, flags, component=0)

        tol = 0.05 if body == MOON else 0.01
        assert returned_speed == pytest.approx(num_speed, abs=tol)

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body",
        [SUN, MOON, MARS, JUPITER],
    )
    def test_lon_speed_at_multiple_dates(self, body):
        """Speed consistency across several dates."""
        flags = FLG_SWIEPH | FLG_SPEED
        dates = [
            2451545.0,  # J2000
            2455197.5,  # 2010
            2458849.5,  # 2020
            2460000.0,  # 2023
        ]
        for jd in dates:
            pos, _ = swe.calc_ut(jd, body, flags)
            returned = pos[3]
            numerical = _numerical_speed(body, jd, flags, component=0)
            tol = 0.05 if body == MOON else 0.01
            assert returned == pytest.approx(numerical, abs=tol), (
                f"Body {body} at JD {jd}: returned {returned}, numerical {numerical}"
            )


class TestSpeedConsistencyLatitude:
    """Test latitude speed vs numerical derivative."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body",
        [MOON, MERCURY, VENUS, MARS, JUPITER],
    )
    def test_lat_speed(self, body):
        """Latitude speed matches numerical derivative."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        returned_speed = pos[4]  # speed in latitude

        num_speed = _numerical_speed(body, JD_J2000, flags, component=1)

        # Latitude speed is generally smaller, so use relative-ish tolerance
        tol = 0.01
        assert returned_speed == pytest.approx(num_speed, abs=tol), (
            f"Body {body}: lat speed returned {returned_speed}, numerical {num_speed}"
        )


class TestSpeedConsistencyDistance:
    """Test distance speed vs numerical derivative."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body",
        [SUN, MOON, MARS, JUPITER],
    )
    def test_dist_speed(self, body):
        """Distance speed matches numerical derivative."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        returned_speed = pos[5]  # speed in distance

        num_speed = _numerical_speed(body, JD_J2000, flags, component=2)

        # Distance speed in AU/day — tolerance 1e-4 AU/day
        tol = 1e-4
        assert returned_speed == pytest.approx(num_speed, abs=tol), (
            f"Body {body}: dist speed returned {returned_speed}, numerical {num_speed}"
        )


class TestSpeedSign:
    """Test that speed sign makes physical sense."""

    @pytest.mark.unit
    def test_sun_speed_positive(self):
        """Sun longitude speed should always be positive (direct motion)."""
        flags = FLG_SWIEPH | FLG_SPEED
        # Check across several dates
        for jd in [JD_J2000, 2458849.5, 2460000.0]:
            pos, _ = swe.calc_ut(jd, SUN, flags)
            assert pos[3] > 0, f"Sun speed negative at JD {jd}: {pos[3]}"

    @pytest.mark.unit
    def test_sun_speed_magnitude(self):
        """Sun longitude speed should be near 1 deg/day."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, SUN, flags)
        assert 0.9 < pos[3] < 1.1, f"Sun speed {pos[3]} not near 1 deg/day"

    @pytest.mark.unit
    def test_moon_speed_magnitude(self):
        """Moon longitude speed should be near 13 deg/day."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, MOON, flags)
        assert 11.0 < abs(pos[3]) < 16.0, f"Moon speed {pos[3]} not near 13"

    @pytest.mark.unit
    def test_mean_node_speed_negative(self):
        """Mean node speed should be negative (retrograde motion)."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, MEAN_NODE, flags)
        assert pos[3] < 0, f"Mean node speed should be negative: {pos[3]}"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body,min_speed,max_speed",
        [
            (MERCURY, -1.5, 2.5),
            (VENUS, -0.8, 1.5),
            (MARS, -0.5, 0.9),
            (JUPITER, -0.2, 0.3),
            (SATURN, -0.1, 0.2),
        ],
    )
    def test_planet_speed_range(self, body, min_speed, max_speed):
        """Planet speed at J2000 is within expected range."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        assert min_speed < pos[3] < max_speed, (
            f"Body {body} speed {pos[3]} out of range [{min_speed}, {max_speed}]"
        )
