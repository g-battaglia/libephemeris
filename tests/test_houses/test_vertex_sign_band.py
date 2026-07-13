"""Analytic one-sided Vertex limits around the equator singularity."""

from __future__ import annotations

import pytest

import libephemeris as le

EPS = 23.4393


def _vertex(lat: float, hsys: str) -> float:
    return le.houses_armc(0.0, lat, EPS, ord(hsys))[1][3]


def _angular_difference(a: float, b: float) -> float:
    return abs((a - b + 180.0) % 360.0 - 180.0)


def test_placidus_one_sided_limits_are_antipodal():
    assert _angular_difference(_vertex(-1e-11, "P"), _vertex(1e-11, "P")) > 179.0


def test_exact_zero_selects_placidus_positive_limit():
    assert _angular_difference(_vertex(0.0, "P"), _vertex(1e-11, "P")) < 1e-6


def test_exact_zero_selects_horizontal_positive_limit():
    assert _angular_difference(_vertex(0.0, "H"), _vertex(1e-11, "H")) < 1e-6


@pytest.mark.unit
class TestHousePosIntHsysObjcoordForm:
    """The objcoord-first house_pos form honors an int house code (character
    code, same convention as the 6-arg form) instead of silently falling
    back to Placidus."""

    def test_int_matches_bytes_and_sixarg(self):
        from libephemeris.houses import house_pos

        i = house_pos(0, 23.4, 23.4393, (40.0, 0.0), ord("S"))
        b = house_pos(0, 23.4, 23.4393, (40.0, 0.0), b"S")
        six = house_pos(0, 23.4, 23.4393, ord("S"), 40.0, 0.0)
        assert i == b == six
        # And it is NOT the Placidus value.
        p = house_pos(0, 23.4, 23.4393, (40.0, 0.0), b"P")
        assert i != p

    def test_omitted_hsys_defaults_placidus(self):
        from libephemeris.houses import house_pos

        d = house_pos(0, 23.4, 23.4393, (40.0, 0.0))
        p = house_pos(0, 23.4, 23.4393, (40.0, 0.0), b"P")
        assert d == p
