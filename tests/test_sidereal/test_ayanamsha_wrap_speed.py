"""Regression: sidereal speed remains finite across an ayanamsha wrap."""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import (
    FLG_SIDEREAL,
    FLG_SPEED,
    SIDM_GALCENT_COCHRANE,
    VENUS,
)


def _signed_angle(angle: float) -> float:
    return (angle + 180.0) % 360.0 - 180.0


def _true_ayanamsha_signed(jd_ut: float) -> float:
    return _signed_angle(ephem.get_ayanamsa_ex_ut(jd_ut, 0)[1])


def _find_catalog_wrap() -> float:
    """Locate the Cochrane catalog-mode wrap without a stored date fixture."""
    left = ephem.julday(2150, 1, 1, 0.0)
    y_left = _true_ayanamsha_signed(left)
    right = left + 1.0
    limit = ephem.julday(2300, 1, 1, 0.0)
    while right <= limit:
        y_right = _true_ayanamsha_signed(right)
        if y_left == 0.0 or y_left * y_right <= 0.0:
            break
        left, y_left = right, y_right
        right += 1.0
    else:
        raise AssertionError("no Cochrane ayanamsha wrap found in 2150--2300")

    for _ in range(50):
        mid = (left + right) / 2.0
        y_mid = _true_ayanamsha_signed(mid)
        if y_left * y_mid <= 0.0:
            right = mid
        else:
            left, y_left = mid, y_mid
    return (left + right) / 2.0


@pytest.mark.unit
def test_no_speed_spike_at_runtime_located_wrap() -> None:
    """The one-second derivative samples use the shortest angular arc."""
    ephem.set_sid_mode(SIDM_GALCENT_COCHRANE)
    wrap_jd = _find_catalog_wrap()
    flags = FLG_SPEED | FLG_SIDEREAL
    step = 0.5 / 86400.0
    try:
        ephem.set_calc_mode("skyfield")
        for k in range(-6, 7):
            jd = wrap_jd + k * step
            pos, _ = ephem.calc_ut(jd, VENUS, flags)
            assert abs(pos[3]) < 10.0, (
                f"sidereal speed spike {pos[3]:+.1f} deg/day at JD {jd:.8f}"
            )
    finally:
        ephem.set_calc_mode("auto")
