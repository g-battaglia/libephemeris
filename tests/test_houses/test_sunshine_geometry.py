"""Geometric invariants for the Sunshine/Treindl cusp construction."""

from __future__ import annotations

import math

import pytest

from libephemeris.exceptions import PolarCircleError


H = __import__("libephemeris.houses", fromlist=["houses"])


def _anchors(armc: float, latitude: float, obliquity: float) -> tuple[float, float]:
    return (
        H._rising_longitude(armc, obliquity, latitude),
        H._armc_to_mc(armc, obliquity),
    )


def _circular_difference(left: float, right: float) -> float:
    return abs((left - right + 180.0) % 360.0 - 180.0)


def test_zero_declination_is_the_same_house_circle_pencil() -> None:
    """An equatorial reference parallel reproduces the rational pencil."""
    for armc, latitude, obliquity in (
        (20.0, 30.0, 5.0),
        (350.0, -30.0, 23.4392911111),
        (200.0, 40.0, 60.0),
    ):
        ascendant, midheaven = _anchors(armc, latitude, obliquity)
        sunshine = H._houses_sunshine(
            armc, latitude, obliquity, ascendant, midheaven, 0.0
        )
        rational = H._houses_regiomontanus(
            armc, latitude, obliquity, ascendant, midheaven
        )
        assert all(
            _circular_difference(left, right) < 1e-12
            for left, right in zip(sunshine[1:], rational[1:])
        )


def test_construction_is_periodic_and_continuous_at_zero_declination() -> None:
    """Whole turns and the two sides of the equatorial limit are continuous."""
    armc, latitude, obliquity = 350.0, -30.0, 23.4392911111
    ascendant, midheaven = _anchors(armc, latitude, obliquity)
    base = H._houses_sunshine(armc, latitude, obliquity, ascendant, midheaven, 0.0)
    rotated = H._houses_sunshine(
        armc + 360.0,
        latitude,
        obliquity,
        ascendant + 360.0,
        midheaven + 360.0,
        0.0,
    )
    below = H._houses_sunshine(armc, latitude, obliquity, ascendant, midheaven, -1e-9)
    above = H._houses_sunshine(armc, latitude, obliquity, ascendant, midheaven, 1e-9)
    assert all(
        _circular_difference(left, right) < 1e-12
        for left, right in zip(base[1:], rotated[1:])
    )
    assert (
        max(
            _circular_difference(left, right)
            for left, right in zip(below[1:], above[1:])
        )
        < 1e-8
    )


def test_zero_offset_is_finite_and_normalized() -> None:
    """The collapsed signed separation has its continuous projection."""
    latitude = 41.9
    declination = 12.0
    obliquity = 23.4392911111
    result = H._sunshine_arc_to_ecliptic(
        0.0,
        True,
        200.0,
        latitude,
        declination,
        math.sin(math.radians(latitude)),
        math.cos(math.radians(latitude)),
        math.cos(math.radians(declination)),
        math.tan(math.radians(declination)),
        math.sin(math.radians(obliquity)),
        math.cos(math.radians(obliquity)),
    )
    assert math.isfinite(result)
    assert 0.0 <= result < 360.0


def test_circumpolar_continuation_and_exact_pole_partition() -> None:
    """Non-crossing paths continue, but an exact geographic pole refuses."""
    armc, latitude, obliquity = 20.0, 80.0, 23.4392911111
    ascendant, midheaven = _anchors(armc, latitude, obliquity)
    cusps = H._houses_sunshine(armc, latitude, obliquity, ascendant, midheaven, 23.0)
    assert len(cusps) == 13
    assert all(math.isfinite(cusp) and 0.0 <= cusp < 360.0 for cusp in cusps)

    near_pole = 89.999999
    ascendant, midheaven = _anchors(armc, near_pole, obliquity)
    assert (
        len(H._houses_sunshine(armc, near_pole, obliquity, ascendant, midheaven, 0.0))
        == 13
    )

    with pytest.raises(PolarCircleError):
        H._houses_sunshine(armc, 90.0, obliquity, 100.0, 280.0, 0.0)


def test_below_horizon_orientation_rotates_only_declared_anchors() -> None:
    """The opposite meridian is reported while the Ascendant stays fixed."""
    armc, latitude, obliquity = 20.0, -80.0, 60.0
    assert H._mc_below_horizon(armc, obliquity, latitude)
    ascendant, midheaven = _anchors(armc, latitude, obliquity)
    cusps = H._houses_sunshine(armc, latitude, obliquity, ascendant, midheaven, 23.0)
    assert cusps[1] == ascendant % 360.0
    assert cusps[7] == (ascendant + 180.0) % 360.0
    assert cusps[10] == (midheaven + 180.0) % 360.0
    assert _circular_difference(cusps[4], midheaven) < 1e-12


@pytest.mark.parametrize(
    "arguments",
    [
        (float("nan"), 0.0, 23.0, 100.0, 200.0, 0.0),
        (0.0, 91.0, 23.0, 100.0, 200.0, 0.0),
        (0.0, 0.0, 90.0, 100.0, 200.0, 0.0),
        (0.0, 0.0, 23.0, 100.0, 200.0, 90.0),
    ],
)
def test_invalid_constructor_inputs_are_rejected(arguments: tuple[float, ...]) -> None:
    """Invalid domains are rejected before spherical evaluation."""
    with pytest.raises(ValueError):
        H._houses_sunshine(*arguments)
