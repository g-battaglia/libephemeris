"""Vertex one-sided limits for tiny latitudes (round-W parity fix).

The reference gives tiny NONZERO latitudes their own sign-side limit at the
equator singularity; only exact lat==0.0 gets the special one-sided value
(positive side for every system, negative side for 'H'). Frozen from the
black-box oracle: houses_armc(0, -1e-11, eps, 'P') Vertex = 0.0 while
lat=0.0 gives 180.0; 'H' at +1e-11 gives 180.0 while lat=0.0 gives 360/0-side.
"""

from __future__ import annotations

import libephemeris as le

EPS = 23.4393


def _vertex(lat: float, hsys: str) -> float:
    return le.houses_armc(0.0, lat, EPS, ord(hsys))[1][3]


def test_tiny_negative_latitude_gets_negative_side_limit():
    assert (
        abs(_vertex(-1e-11, "P") - 0.0) < 1e-6
        or abs(_vertex(-1e-11, "P") - 360.0) < 1e-6
    )


def test_exact_zero_latitude_keeps_positive_side_limit():
    assert abs(_vertex(0.0, "P") - 180.0) < 1e-6


def test_h_system_tiny_positive_latitude_gets_positive_side_limit():
    assert abs(_vertex(1e-11, "H") - 180.0) < 1e-6


def test_h_system_exact_zero_keeps_negative_side_limit():
    v = _vertex(0.0, "H")
    assert abs(v - 0.0) < 1e-6 or abs(v - 360.0) < 1e-6
