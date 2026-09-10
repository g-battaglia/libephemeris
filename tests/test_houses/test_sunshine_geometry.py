"""Geometric invariants for the Sunshine/Treindl cusp construction."""

from __future__ import annotations

import math

import pytest

from libephemeris.exceptions import CalculationError, PolarCircleError


H = __import__("libephemeris.houses", fromlist=["houses"])


def _anchors(armc: float, latitude: float, obliquity: float) -> tuple[float, float]:
    return (
        H._rising_longitude(armc, obliquity, latitude),
        H._armc_to_mc(armc, obliquity),
    )


def _circular_difference(left: float, right: float) -> float:
    return abs((left - right + 180.0) % 360.0 - 180.0)


def _plane_intersection_expected(
    arc_offset: float,
    is_diurnal: bool,
    armc: float,
    latitude: float,
    declination: float,
    obliquity: float,
) -> tuple[float, float]:
    """Intersect the independently constructed house and ecliptic planes."""
    phi = math.radians(latitude)
    theta = math.radians(armc)
    epsilon = math.radians(obliquity)
    point_ra = math.radians((armc if is_diurnal else armc + 180.0) + arc_offset)
    delta = math.radians(declination)
    north = (
        -math.sin(phi) * math.cos(theta),
        -math.sin(phi) * math.sin(theta),
        math.cos(phi),
    )
    point = (
        math.cos(delta) * math.cos(point_ra),
        math.cos(delta) * math.sin(point_ra),
        math.sin(delta),
    )
    ecliptic_normal = (0.0, -math.sin(epsilon), math.cos(epsilon))

    def cross(left: tuple[float, ...], right: tuple[float, ...]) -> tuple[float, ...]:
        return (
            left[1] * right[2] - left[2] * right[1],
            left[2] * right[0] - left[0] * right[2],
            left[0] * right[1] - left[1] * right[0],
        )

    intersection = cross(ecliptic_normal, cross(north, point))
    longitudes = []
    for direction in (intersection, tuple(-value for value in intersection)):
        ecliptic_y = direction[1] * math.cos(epsilon) + direction[2] * math.sin(epsilon)
        longitudes.append(math.degrees(math.atan2(ecliptic_y, direction[0])) % 360.0)
    return tuple(longitudes)


def test_horizontal_degeneracies_are_finite_and_antipodal() -> None:
    """Horizontal houses preserve paired geometry at parity and plane coincidences."""
    for armc, latitude, obliquity in (
        (0.0, 0.0, 23.4392911),
        (0.0, 0.0, 0.0),
        (0.0, 90.0, 0.0),
        (0.0, -90.0, 0.0),
        (0.0, 0.0, 90.0),
        (180.0, 0.0, 90.0),
    ):
        ascendant, midheaven = _anchors(armc, latitude, obliquity)
        cusps = H._houses_horizontal(armc, latitude, obliquity, ascendant, midheaven)
        assert len(cusps) == 13
        assert all(math.isfinite(value) and 0.0 <= value < 360.0 for value in cusps)
        assert all(
            _circular_difference(cusps[index + 6], cusps[index] + 180.0) < 1e-12
            for index in range(1, 7)
        )


def test_horizontal_poles_share_the_northern_limit() -> None:
    """Both exact poles use the same deterministic horizontal-house ring."""
    armc, obliquity = 123.0, 23.4392911
    north = H._houses_horizontal(
        armc, 90.0, obliquity, *_anchors(armc, 90.0, obliquity)
    )
    south = H._houses_horizontal(
        armc, -90.0, obliquity, *_anchors(armc, -90.0, obliquity)
    )
    assert all(
        _circular_difference(left, right) < 1e-12 for left, right in zip(north, south)
    )


def test_horizontal_whole_turns_leave_the_ring_unchanged() -> None:
    """ARMC and obliquity whole turns preserve every horizontal cusp."""
    armc, latitude, obliquity = 211.0, -41.9, 23.4392911
    base = H._houses_horizontal(
        armc, latitude, obliquity, *_anchors(armc, latitude, obliquity)
    )
    rotated = H._houses_horizontal(
        armc + 360.0,
        latitude,
        obliquity + 360.0,
        *_anchors(armc + 360.0, latitude, obliquity + 360.0),
    )
    assert all(
        _circular_difference(left, right) < 1e-12 for left, right in zip(base, rotated)
    )


@pytest.mark.parametrize(
    ("arc_offset", "is_diurnal"),
    [(20.0, True), (-20.0, True), (20.0, False), (-20.0, False)],
)
def test_helper_selects_one_direct_plane_intersection(
    arc_offset: float, is_diurnal: bool
) -> None:
    """The returned cusp is one antipode of the direct plane intersection."""
    armc = 200.0
    latitude = 41.9
    declination = 10.0
    obliquity = 23.4392911
    actual = H._sunshine_arc_to_ecliptic(
        arc_offset,
        is_diurnal,
        armc,
        latitude,
        declination,
        math.sin(math.radians(latitude)),
        math.cos(math.radians(latitude)),
        math.cos(math.radians(declination)),
        math.tan(math.radians(declination)),
        math.sin(math.radians(obliquity)),
        math.cos(math.radians(obliquity)),
    )
    intersections = _plane_intersection_expected(
        arc_offset, is_diurnal, armc, latitude, declination, obliquity
    )
    assert (
        min(
            _circular_difference(actual, intersection) for intersection in intersections
        )
        < 1e-12
    )


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
    with pytest.raises(PolarCircleError):
        H.houses_armc(armc, 90.0, obliquity, ord("I"), 0.0)


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


def test_exact_meridian_horizon_contact_does_not_rotate() -> None:
    """Exact contact belongs to the unrotated branch without a tolerance."""
    armc = 90.0
    obliquity = 23.44
    latitude = -66.56
    ascendant, midheaven = _anchors(armc, latitude, obliquity)
    cusps = H._houses_sunshine(armc, latitude, obliquity, ascendant, midheaven, 0.0)
    assert cusps[10] == midheaven % 360.0
    assert _circular_difference(cusps[4], midheaven + 180.0) < 1e-12


@pytest.mark.parametrize("error", [IndexError, TypeError, CalculationError])
def test_unrelated_sun_provider_errors_propagate(monkeypatch, error) -> None:
    """Only provider coverage failures select the analytic fallback."""

    def boom(*args, **kwargs):
        raise error("provider failure")

    monkeypatch.setattr(H, "calc_ut", boom)
    with pytest.raises(error, match="provider failure"):
        H._sunshine_sun_declination(2451545.0, 0)


def test_ecliptic_pole_intersection_is_rejected(monkeypatch) -> None:
    """An exactly null homogeneous longitude pair has no longitude."""
    inverse_sines = iter((0.0, math.radians(15.0)))
    monkeypatch.setattr(H.math, "asin", lambda value: next(inverse_sines))
    with pytest.raises(CalculationError, match="longitude is not uniquely defined"):
        H._sunshine_arc_to_ecliptic(
            0.0,
            True,
            0.0,
            0.0,
            0.0,
            math.sin(0.0),
            math.cos(0.0),
            math.cos(0.0),
            math.tan(0.0),
            math.sin(math.radians(75.0)),
            math.cos(math.radians(75.0)),
        )


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
