"""FLG_SIDEREAL invariants for asteroid propagation fallback paths.

The ASSIST and Keplerian paths used outside an asteroid's SPK interval must
apply the same frame offset as the direct SPK path. These tests use the native
True Citra definition and check the coordinate identity directly.
"""

from __future__ import annotations

import pytest

import libephemeris as L
from libephemeris.constants import (
    CERES,
    VESTA,
    CHIRON,
    FLG_SWIEPH,
    FLG_SIDEREAL,
    SIDM_TRUE_CITRA,
)


def _signed_angle(angle: float) -> float:
    return (angle + 180.0) % 360.0 - 180.0


@pytest.fixture(autouse=True)
def _true_citra_skyfield():
    L.set_calc_mode("skyfield")
    L.set_sid_mode(SIDM_TRUE_CITRA)
    yield
    L.set_calc_mode("auto")
    L.set_sid_mode(L.SIDM_J2000)


# 2000 = inside SPK window (normal path); 2140 / 1860 = outside -> fallback path.
@pytest.mark.parametrize("year", [2000, 2140, 1860])
@pytest.mark.parametrize("body", [CERES, VESTA, CHIRON])
def test_asteroid_sidereal_applies_ayanamsa(body, year):
    jd = 2451545.0 + (year - 2000) * 365.25
    trop = L.calc(jd, body, FLG_SWIEPH)[0][0]
    sid = L.calc(jd, body, FLG_SWIEPH | FLG_SIDEREAL)[0][0]
    ayan = L.get_ayanamsa(jd)
    # sidereal = tropical - ayanamsa (the defining contract, all bodies)
    residual = _signed_angle((trop - sid) - ayan)
    assert abs(residual) < 0.01, (
        f"body {body} @ {year}: frame residual={residual:.4f} deg"
    )
    # and it actually shifted (guards the silent-tropical regression)
    assert abs(_signed_angle(trop - sid)) > 1.0
