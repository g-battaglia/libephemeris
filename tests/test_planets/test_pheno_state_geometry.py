"""Independent state-geometry tests for planetary phenomena."""

from __future__ import annotations

import math

import pytest

import libephemeris as le
from libephemeris.exceptions import InputValidationError
from libephemeris.planets import (
    _moon_horizontal_parallax_deg,
    _phenomena_geometry,
)


def _angle(left: tuple[float, ...], right: tuple[float, ...]) -> float:
    cross = (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )
    return math.degrees(
        math.atan2(
            math.sqrt(sum(value * value for value in cross)),
            sum(a * b for a, b in zip(left, right)),
        )
    )


@pytest.mark.parametrize(
    ("body", "sun"),
    [
        ((1.0, 0.0, 0.0), (2.0, 1e-15, 0.0)),
        ((1.0, 0.0, 0.0), (0.0, 1e-15, 0.0)),
        ((0.2, -0.3, 0.4), (1.1, 0.7, -0.2)),
    ],
)
def test_geometry_matches_independent_cross_dot_angles(body, sun) -> None:
    """Phase and elongation use stable vector angles, including endpoints."""
    phase_angle, fraction, elongation = _phenomena_geometry(body, sun)
    body_to_observer = tuple(-value for value in body)
    body_to_sun = tuple(a - b for a, b in zip(sun, body))
    assert phase_angle == pytest.approx(_angle(body_to_observer, body_to_sun))
    assert elongation == pytest.approx(_angle(body, sun))
    assert fraction == pytest.approx((1.0 + math.cos(math.radians(phase_angle))) / 2.0)


def test_apparent_diameter_uses_exact_spherical_geometry() -> None:
    """A finite spherical body uses twice its exact angular radius."""
    from libephemeris.planets import _calc_apparent_diameter

    radius_km = 1737.4
    distance_au = 0.00257
    expected = 2.0 * math.degrees(math.asin(radius_km / (distance_au * 149597870.7)))
    assert _calc_apparent_diameter(radius_km, distance_au) == expected


@pytest.mark.parametrize("distance_au", [0.0024, 0.00257, 0.00272])
def test_horizontal_parallax_matches_earth_radius_relation(distance_au: float) -> None:
    """Geocentric lunar slot 5 follows asin(R_E / Delta)."""
    expected = math.degrees(math.asin(6378.1366 / (distance_au * 149597870.7)))
    assert _moon_horizontal_parallax_deg(distance_au) == expected


@pytest.mark.parametrize("bad", [0.0, -1.0, math.nan, math.inf, "0.0025"])
def test_horizontal_parallax_rejects_invalid_distance(bad) -> None:
    """The pure helper rejects invalid geometry instead of returning a sentinel."""
    with pytest.raises(InputValidationError):
        _moon_horizontal_parallax_deg(bad)


@pytest.mark.parametrize(
    "flags",
    [
        le.FLG_NOABERR,
        le.FLG_NOGDEFL,
        le.FLG_ASTROMETRIC,
        le.FLG_TRUEPOS | le.FLG_NOABERR,
    ],
)
def test_sealed_leb_keeps_reduction_flags_on_the_leb_path(
    flags: int, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reduction flags must not open the direct ephemeris in sealed LEB mode."""
    import libephemeris.planets as planets

    def fail_direct_path():
        raise AssertionError("direct ephemeris path used in sealed LEB mode")

    monkeypatch.setattr(planets, "_get_computation_ephemeris", fail_direct_path)
    result = le.pheno_ut(2451545.0, le.MARS, flags)
    assert len(result) == 20
    assert all(type(value) is float for value in result)


def test_public_sentinels_and_ignored_heliocentric_flag() -> None:
    """Degenerate tuples remain exact and HELCTR does not change phenomena."""
    earth = le.pheno_ut(2451545.0, le.EARTH, 0)
    assert earth == (0.0, 0.0, 0.0, 180.0) + (0.0,) * 16
    nutation = le.pheno_ut(2451545.0, le.ECL_NUT, 0)
    assert all(math.isnan(value) for value in nutation[:3])
    assert nutation[3:] == (0.0,) * 17
    plain = le.pheno_ut(2451545.0, le.MARS, 0)
    heliocentric = le.pheno_ut(2451545.0, le.MARS, le.FLG_HELCTR)
    assert heliocentric == plain
    assert len(plain) == 20 and all(type(value) is float for value in plain)


def test_lunar_parallax_truepos_invariance_and_topocentric_displacement() -> None:
    """Slot 5 keeps its own geometry while topocentric parallax is observer-specific."""
    jd = 2451545.0
    geocentric = le.pheno_ut(jd, le.MOON, 0)
    geometric = le.pheno_ut(jd, le.MOON, le.FLG_TRUEPOS)
    assert geometric[5] == geocentric[5]
    assert all(value == 0.0 for value in le.pheno_ut(jd, le.MARS, 0)[5:])

    le.set_topo(12.5, 41.9, 0.0)
    topocentric = le.pheno_ut(jd, le.MOON, le.FLG_TOPOCTR)
    topocentric_geometric = le.pheno_ut(jd, le.MOON, le.FLG_TOPOCTR | le.FLG_TRUEPOS)
    assert topocentric[5] > 0.0
    assert topocentric_geometric[5] == topocentric[5]
    assert topocentric[5] != geocentric[5]
    le.close()
