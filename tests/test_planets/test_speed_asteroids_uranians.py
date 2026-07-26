"""Independent speed checks for asteroids and sourced fictitious bodies."""

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
from libephemeris.exceptions import UnknownBodyError


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    swe.close()


JD_J2000 = 2451545.0
DT = 0.01  # 0.01 days for numerical derivative

REVIEWED_FICTITIOUS_IDS = (*range(40, 49), *range(50, 59))
UNVERIFIED_FICTITIOUS_IDS = (49,)


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


class TestFictitiousBodies:
    """Only independently sourced built-ins enter the speed pipeline."""

    @pytest.mark.unit
    def test_harrington_speed_request_returns_finite_state(self) -> None:
        swe.set_calc_mode("skyfield")
        position, _ = swe.calc_ut(JD_J2000, swe.HARRINGTON, FLG_SWIEPH | FLG_SPEED)
        assert len(position) == 6
        assert all(math.isfinite(value) for value in position)

    @pytest.mark.unit
    @pytest.mark.parametrize("body", REVIEWED_FICTITIOUS_IDS)
    def test_reviewed_body_speed_request_returns_finite_state(self, body: int) -> None:
        swe.set_calc_mode("skyfield")
        position, _ = swe.calc_ut(JD_J2000, body, FLG_SWIEPH | FLG_SPEED)
        assert len(position) == 6
        assert all(math.isfinite(value) for value in position)

    @pytest.mark.unit
    @pytest.mark.parametrize("body", UNVERIFIED_FICTITIOUS_IDS)
    def test_unverified_body_speed_request_fails_closed(self, body: int) -> None:
        swe.set_calc_mode("skyfield")
        with pytest.raises(UnknownBodyError) as raised:
            swe.calc_ut(JD_J2000, body, FLG_SWIEPH | FLG_SPEED)
        assert raised.value.body_id == body
