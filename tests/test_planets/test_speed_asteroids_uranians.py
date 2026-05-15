"""
Tests for speed consistency of asteroids and Uranian bodies.

Verifies that speed values returned by calc_ut match numerical
differentiation for Chiron, Ceres, Pallas, Juno, Vesta, and Uranians.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    CHIRON,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    FLG_SWIEPH,
    FLG_SPEED,
)


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    swe.close()


JD_J2000 = 2451545.0
DT = 0.01  # 0.01 days for numerical derivative

# Uranian body IDs
CUPIDO = 40
HADES = 41
ZEUS = 42
KRONOS = 43
APOLLON = 44
ADMETOS = 45
VULKANUS = 46
POSEIDON = 47


def _numerical_speed(body, jd, flags, component=0):
    """Compute numerical derivative via central differences."""
    pos_m, _ = swe.calc_ut(jd - DT, body, flags)
    pos_p, _ = swe.calc_ut(jd + DT, body, flags)
    val_m = pos_m[component]
    val_p = pos_p[component]

    if component == 0:
        diff = val_p - val_m
        if diff > 180.0:
            diff -= 360.0
        elif diff < -180.0:
            diff += 360.0
        return diff / (2 * DT)

    return (val_p - val_m) / (2 * DT)


class TestAsteroidSpeedConsistency:
    """Test speed consistency for main asteroids."""

    ASTEROIDS = [CHIRON, CERES, PALLAS, JUNO, VESTA]

    @pytest.mark.unit
    @pytest.mark.parametrize("body", ASTEROIDS)
    def test_lon_speed(self, body):
        """Asteroid longitude speed matches numerical derivative."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        returned = pos[3]
        numerical = _numerical_speed(body, JD_J2000, flags, 0)
        assert returned == pytest.approx(numerical, abs=0.02), (
            f"Body {body}: returned {returned}, numerical {numerical}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("body", ASTEROIDS)
    def test_lat_speed(self, body):
        """Asteroid latitude speed matches numerical derivative."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        returned = pos[4]
        numerical = _numerical_speed(body, JD_J2000, flags, 1)
        assert returned == pytest.approx(numerical, abs=0.02), (
            f"Body {body}: lat speed returned {returned}, numerical {numerical}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("body", ASTEROIDS)
    def test_dist_speed(self, body):
        """Asteroid distance speed matches numerical derivative."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        returned = pos[5]
        numerical = _numerical_speed(body, JD_J2000, flags, 2)
        assert returned == pytest.approx(numerical, abs=1e-3), (
            f"Body {body}: dist speed returned {returned}, numerical {numerical}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("body", ASTEROIDS)
    def test_speed_magnitude_reasonable(self, body):
        """Asteroid speed magnitude is within physical limits."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        # Asteroids move slower than inner planets
        assert abs(pos[3]) < 1.0, f"Body {body}: lon speed {pos[3]} too fast"


class TestUranianSpeedConsistency:
    """Test speed consistency for Uranian hypothetical bodies."""

    URANIANS = [CUPIDO, HADES, ZEUS, KRONOS, APOLLON, ADMETOS, VULKANUS, POSEIDON]

    @pytest.mark.unit
    @pytest.mark.parametrize("body", URANIANS)
    def test_lon_speed(self, body):
        """Uranian longitude speed matches numerical derivative."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        returned = pos[3]
        numerical = _numerical_speed(body, JD_J2000, flags, 0)
        assert returned == pytest.approx(numerical, abs=0.01), (
            f"Body {body}: returned {returned}, numerical {numerical}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("body", URANIANS)
    def test_speed_magnitude(self, body):
        """Uranian body speed is within expected range."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        # Uranians are very slow-moving (trans-Neptunian orbit periods)
        assert abs(pos[3]) < 0.5, f"Body {body}: lon speed {pos[3]} too fast"

    @pytest.mark.unit
    @pytest.mark.parametrize("body", URANIANS)
    def test_speed_at_2020(self, body):
        """Uranian speed also consistent at a different date."""
        flags = FLG_SWIEPH | FLG_SPEED
        jd = 2458849.5  # 2020
        pos, _ = swe.calc_ut(jd, body, flags)
        returned = pos[3]
        numerical = _numerical_speed(body, jd, flags, 0)
        assert returned == pytest.approx(numerical, abs=0.01), (
            f"Body {body} at 2020: returned {returned}, numerical {numerical}"
        )
