"""Source-backed checks for independently implemented house constructions."""

from __future__ import annotations

import math

import pytest

from libephemeris.houses import (
    houses_armc,
)
from libephemeris.house_constructions import (
    _cross,
    _dot,
    _house_circle_lon_deg,
    _normalized,
    houses_krusinski,
    houses_savard_a,
    houses_sunshine_makransky,
)


def _angular_distance_deg(first: float, second: float) -> float:
    """Return the unsigned shortest separation of two cyclic longitudes."""
    return abs((first - second + 180.0) % 360.0 - 180.0)


@pytest.mark.parametrize(
    ("axis", "anchor", "eps_deg"),
    [
        ((0.17, 0.91, 0.38), (0.84, -0.22, 0.49), 23.4392911),
        ((-0.73, 0.31, 0.61), (0.12, 0.96, -0.25), 23.4392911),
        ((0.64, 0.57, -0.51), (-0.44, 0.19, 0.88), 17.0),
        ((-0.28, -0.89, 0.36), (0.79, 0.46, 0.41), 29.0),
    ],
)
@pytest.mark.parametrize("same_side", [False, True])
def test_house_circle_branch_is_selected_by_geometric_semicircle(
    axis: tuple[float, float, float],
    anchor: tuple[float, float, float],
    eps_deg: float,
    same_side: bool,
) -> None:
    """The selected ecliptic antipode occupies the requested circle half.

    A longitude/quadrant predicate is not invariant at 0/360 degrees or under
    a polar horizon-frame flip. The implementation instead uses the sign of
    the in-plane coordinate perpendicular to the north/south axis. This test
    verifies that defining invariant directly, without any ephemeris oracle.
    """
    axis_unit = _normalized(axis)
    anchor_unit = _normalized(anchor)
    eps_rad = math.radians(eps_deg)
    longitude = math.radians(
        _house_circle_lon_deg(
            anchor_unit,
            axis_unit,
            eps_rad,
            same_side=same_side,
        )
    )
    intersection = (
        math.cos(longitude),
        math.sin(longitude) * math.cos(eps_rad),
        math.sin(longitude) * math.sin(eps_rad),
    )
    normal = _cross(axis_unit, anchor_unit)
    half_coordinate = _cross(normal, axis_unit)
    product = _dot(intersection, half_coordinate) * _dot(anchor_unit, half_coordinate)

    assert _dot(normal, intersection) == pytest.approx(0.0, abs=2e-15)
    assert (product > 0.0) is same_side


@pytest.mark.parametrize(
    ("armc", "latitude", "ascendant", "midheaven"),
    [
        (12.3, 0.0, 101.2, 13.4),
        (219.7, 41.9, 305.1, 217.4),
        (71.0, -55.0, 159.0, 73.0),
        (359.9, 89.0, 270.1, 359.9),
    ],
)
def test_source_defined_opposite_cusps_remain_antipodal(
    armc: float, latitude: float, ascendant: float, midheaven: float
) -> None:
    """Savard-A and Krusinski explicitly complete six antipodal axes."""
    for cusps in (
        houses_savard_a(armc, latitude, 23.4392911, ascendant, midheaven),
        houses_krusinski(armc, latitude, 23.4392911, ascendant, midheaven),
    ):
        for first, opposite in ((1, 7), (2, 8), (3, 9), (4, 10), (5, 11), (6, 12)):
            separation = (cusps[opposite] - cusps[first]) % 360.0
            assert separation == pytest.approx(180.0, abs=2e-12)


def test_makransky_published_worked_example() -> None:
    """Reproduce the worked Solar House Cusp Algorithm example.

    Source: M. J. Makransky, "The Solar House System" (1990), in
    *Primary Directions: A Primer of Calculation*, pp. 147--148. The author
    gives DEC=-18.38, LAT=51.50 and RAMC=12.37, followed by the eight
    intermediate cusps rounded to whole arcminutes. The scan was distributed
    by the author and is preserved at:

    https://web.archive.org/web/20160308075443id_/http://dearbrutus.com/primarydirections_3.pdf
    """
    cusps = houses_sunshine_makransky(
        armc=12.37,
        lat=51.50,
        eps=23.4392911,
        asc=0.0,
        mc=0.0,
        sun_dec=-18.38,
    )
    published = {
        2: 147.300000,
        3: 167.066667,
        5: 233.100000,
        6: 274.883333,
        8: 321.200000,
        9: 340.550000,
        11: 67.066667,
        12: 105.483333,
    }

    for index, expected in published.items():
        assert cusps[index] == pytest.approx(expected, abs=1.0 / 60.0)


@pytest.mark.parametrize(
    ("latitude", "sun_declination"),
    [
        (-80.0, 8.7),
        (-66.4, 23.3),
        (60.0, -23.3),
        (66.4, -23.3),
    ],
)
def test_makransky_high_latitude_rotation_equivariance(
    latitude: float, sun_declination: float
) -> None:
    """The high-latitude acute-MD branch is equivariant under RA rotation."""
    base = houses_sunshine_makransky(
        17.3,
        latitude,
        0.0,
        0.0,
        0.0,
        sun_declination,
    )
    rotated = houses_sunshine_makransky(
        154.7,
        latitude,
        0.0,
        0.0,
        0.0,
        sun_declination,
    )

    shift = 154.7 - 17.3
    for index in (2, 3, 5, 6, 8, 9, 11, 12):
        delta = (rotated[index] - base[index] - shift + 180.0) % 360.0 - 180.0
        assert delta == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("armc", [float(value) for value in range(0, 360, 30)])
@pytest.mark.parametrize("latitude", [-60.0, -33.9, 0.0, 41.9, 60.0])
def test_makransky_is_continuous_at_zero_solar_declination(
    armc: float, latitude: float
) -> None:
    """The equinoctial solar parallel has one continuous geometric limit.

    Makransky's published first step is

    ``AD = asin(tan(DEC) * tan(LAT))``.

    Therefore ``DEC -> 0`` implies ``AD -> 0`` from both sides, equal
    90-degree diurnal/nocturnal semiarcs, and the same eight division points.
    His later pole and ascensional-difference equations are compositions of
    continuous trigonometric functions in this non-circumpolar domain.  No
    source-defined cusp swap occurs at the equinox.

    This test deliberately surrounds zero instead of comparing against an
    ephemeris oracle.  A wrong antipodal/quadrant branch would differ by about
    180 degrees; the 2e-8-degree bound is more than an order of magnitude wider
    than the observed floating-point motion for a 1e-9-degree perturbation,
    while remaining far too small to conceal such a branch error.
    """
    exact = houses_sunshine_makransky(
        armc,
        latitude,
        23.4392911,
        0.0,
        0.0,
        0.0,
    )

    for nearby_declination in (-1e-9, 1e-9):
        nearby = houses_sunshine_makransky(
            armc,
            latitude,
            23.4392911,
            0.0,
            0.0,
            nearby_declination,
        )
        for index in (2, 3, 5, 6, 8, 9, 11, 12):
            assert _angular_distance_deg(nearby[index], exact[index]) < 2e-8


@pytest.mark.parametrize("armc", [float(value) for value in range(0, 360, 30)])
@pytest.mark.parametrize("sun_declination", [-23.4, 0.0, 23.4])
def test_makransky_equator_is_declination_independent_and_cyclic(
    armc: float, sun_declination: float
) -> None:
    """At latitude zero, the source construction is an ordered RA division.

    With ``LAT = 0``, Makransky's ascensional difference, every house-circle
    pole, and every auxiliary ``Q`` are zero for any non-polar solar
    declination.  All twelve cusps must therefore equal the declination-zero
    solution and retain cusp-number order around the ecliptic.  Obliquity makes
    the ecliptic-longitude gaps alternate around 30 degrees, so equality of
    longitude gaps is *not* the correct invariant; positivity and a complete
    360-degree winding are.

    Exercising :func:`houses_armc` also supplies the geometrically calculated
    Ascendant and Midheaven, rather than arbitrary placeholders for the four
    angular cusps.  This is a source-derived limit test and stores no external
    implementation output.

    The armc entry point returns the continuous construction at every armc,
    including the exact colure hits, so the winding invariant applies to the
    public output directly.
    """
    cusps, _ = houses_armc(
        armc,
        0.0,
        23.4392911,
        ord("i"),
        ascmc9=sun_declination,
    )
    equinoctial_cusps, _ = houses_armc(
        armc,
        0.0,
        23.4392911,
        ord("i"),
        ascmc9=0.0,
    )

    for actual, expected in zip(cusps, equinoctial_cusps, strict=True):
        assert _angular_distance_deg(actual, expected) < 2e-12

    unwrapped = [cusps[0]]
    for longitude in cusps[1:]:
        candidate = longitude
        while candidate <= unwrapped[-1]:
            candidate += 360.0
        unwrapped.append(candidate)

    gaps = [
        unwrapped[index + 1] - unwrapped[index] for index in range(len(unwrapped) - 1)
    ]
    gaps.append(unwrapped[0] + 360.0 - unwrapped[-1])
    assert all(gap > 0.0 for gap in gaps)
    assert sum(gaps) == pytest.approx(360.0, abs=2e-12)


@pytest.mark.parametrize("armc", [30.0, 60.0, 120.0, 150.0, 210.0, 240.0, 300.0, 330.0])
@pytest.mark.parametrize("latitude", [-60.0, -33.9, 0.0, 41.9, 60.0])
def test_sunshine_armc_colure_hits_stay_on_the_construction(
    armc: float, latitude: float
) -> None:
    """Exact colure hits return the construction itself, not a special case.

    With zero declination the division points sit at exact multiples of 30
    degrees from the meridians, so these armc values place one point on the
    equinoctial colure and two on the solstitial colure.  Makransky's
    published construction is continuous through those points, so the armc
    entry point must return it unchanged there — and the value must equal
    the two-sided limit of the same construction.  An earlier revision
    replaced these isolated points with antipode/mirror transforms inferred
    from a grid of external outputs; that reconstruction was removed (see
    docs/comparison/intentional-divergences.md).
    """
    public_cusps, _ = houses_armc(armc, latitude, 23.4392911, ord("i"))
    construction = houses_sunshine_makransky(
        armc,
        latitude,
        23.4392911,
        0.0,
        0.0,
        0.0,
    )

    for index in (2, 3, 5, 6, 8, 9, 11, 12):
        assert (
            _angular_distance_deg(public_cusps[index - 1], construction[index]) < 1e-12
        )

    # The exact-hit value is the two-sided limit: stepping a microdegree to
    # either side moves every division point by less than an arcsecond.
    for offset in (-1e-6, 1e-6):
        nearby, _ = houses_armc(armc + offset, latitude, 23.4392911, ord("i"))
        for index in (2, 3, 5, 6, 8, 9, 11, 12):
            assert (
                _angular_distance_deg(public_cusps[index - 1], nearby[index - 1])
                < 1.0 / 3600.0
            )
