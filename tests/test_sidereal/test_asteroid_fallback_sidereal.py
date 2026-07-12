"""Regression: FLG_SIDEREAL must be applied to asteroids on the fallback paths.

A differential sweep (LEB vs Skyfield) found the asteroid ASSIST / Keplerian
fallback paths (used outside the ~1920-2080 SPK window) applied FLG_EQUATORIAL /
FLG_J2000 / FLG_XYZ but silently dropped FLG_SIDEREAL — returning a TROPICAL
longitude for a sidereal request (~26 deg wrong). Inside the SPK window the
sidereal offset was applied correctly, so the bug only bit at extreme dates.

These check the contract that holds for every other body: sidereal longitude =
(tropical - ayanamsa) mod 360, across both the in-range and the fallback path.
"""

from __future__ import annotations

import pytest

import libephemeris as L
from libephemeris.constants import CERES, VESTA, CHIRON, FLG_SWIEPH, FLG_SIDEREAL


def _wrap(a, b):
    d = abs(a - b) % 360.0
    return min(d, 360.0 - d)


@pytest.fixture(autouse=True)
def _lahiri_skyfield():
    L.set_calc_mode("skyfield")
    L.set_sid_mode(1, 0.0, 0.0)  # Lahiri
    yield
    L.set_calc_mode("auto")
    L.set_sid_mode(0, 0.0, 0.0)


# 2000 = inside SPK window (normal path); 2140 / 1860 = outside -> fallback path.
@pytest.mark.parametrize("year", [2000, 2140, 1860])
@pytest.mark.parametrize("body", [CERES, VESTA, CHIRON])
def test_asteroid_sidereal_applies_ayanamsa(body, year):
    jd = 2451545.0 + (year - 2000) * 365.25
    trop = L.calc(jd, body, FLG_SWIEPH)[0][0]
    sid = L.calc(jd, body, FLG_SWIEPH | FLG_SIDEREAL)[0][0]
    ayan = L.get_ayanamsa_ut(jd)
    # sidereal = tropical - ayanamsa (the defining contract, all bodies)
    assert abs(_wrap(trop, sid) - ayan) < 0.01, (
        f"body {body} @ {year}: trop-sid={_wrap(trop, sid):.4f} vs ayan={ayan:.4f}"
    )
    # and it actually shifted (guards the silent-tropical regression)
    assert _wrap(trop, sid) > 1.0
