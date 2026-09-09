"""Direct tests for the published planetary magnitude models."""

from __future__ import annotations

import math

import pytest

from libephemeris.constants import (
    JUPITER,
    MARS,
    MERCURY,
    NEPTUNE,
    SATURN,
    URANUS,
    VENUS,
)
from libephemeris.exceptions import InputValidationError, UnknownBodyError
from libephemeris.planets import (
    _calc_neptune_magnitude,
    _calc_planet_magnitude,
    _direction_to_icrs,
    _pole_vector,
    _sub_latitudes,
    _sun_magnitude,
    _uranus_photometric_latitude,
)
from libephemeris.time_utils import julday


J2000 = 2451545.0


def _magnitude(body: int, phase: float, **kwargs: float) -> float:
    return _calc_planet_magnitude(body, 1.0, 1.0, phase, tjd=J2000, **kwargs)


@pytest.mark.parametrize(
    ("body", "phase", "expected"),
    [
        (MERCURY, 0.0, -0.613),
        (VENUS, 0.0, -4.384),
        (MARS, 0.0, -1.601),
        (JUPITER, 0.0, -9.395),
        (NEPTUNE, 0.0, -7.0),
    ],
)
def test_published_zero_phase_constants(
    body: int, phase: float, expected: float
) -> None:
    """Each direct law reproduces its printed zero-phase constant."""
    assert _magnitude(body, phase) == expected


def test_distance_law_is_exactly_inverse_square() -> None:
    """Doubling either distance adds the common published distance modulus."""
    shift = 5.0 * math.log10(2.0)
    for body, phase in (
        (MERCURY, 20.0),
        (VENUS, 20.0),
        (MARS, 20.0),
        (JUPITER, 5.0),
        (SATURN, 1.0),
        (URANUS, 20.0),
        (NEPTUNE, 1.0),
    ):
        kwargs = {
            "geo_lon": 0.0,
            "geo_lat": 0.0,
            "helio_lon": 0.0,
            "helio_lat": 0.0,
        }
        one = _calc_planet_magnitude(body, 1.0, 1.0, phase, tjd=J2000, **kwargs)
        two = _calc_planet_magnitude(body, 2.0, 1.0, phase, tjd=J2000, **kwargs)
        assert two - one == pytest.approx(shift, abs=2e-15)


def test_venus_and_mars_use_first_equation_at_break() -> None:
    """Exact branch boundaries belong to the lower-phase published equations."""
    assert _magnitude(VENUS, 163.7) == pytest.approx(
        -4.384
        - 1.044e-3 * 163.7
        + 3.687e-4 * 163.7**2
        - 2.814e-6 * 163.7**3
        + 8.938e-9 * 163.7**4
    )
    assert _magnitude(MARS, 50.0) == pytest.approx(
        -1.601 + 0.02267 * 50.0 - 0.0001302 * 50.0**2
    )


def test_saturn_opposite_faces_zero_ring_opening() -> None:
    """Opposite Saturn ring faces force the effective beta to zero."""
    opposite = _magnitude(
        SATURN,
        0.0,
        geo_lon=0.0,
        geo_lat=0.0,
        helio_lon=180.0,
        helio_lat=0.0,
    )
    expected = -8.914 - 0.378 * 0.0 + 0.026 * 0.0
    assert opposite == expected


def test_saturn_same_faces_use_geometric_mean_and_mirror_invariance() -> None:
    """Same-face Saturn openings use sqrt(abs(beta_E*beta_S))."""
    first = _magnitude(
        SATURN, 0.0, geo_lon=0.0, geo_lat=0.0, helio_lon=0.0, helio_lat=0.0
    )
    mirrored = _magnitude(
        SATURN, 0.0, geo_lon=180.0, geo_lat=0.0, helio_lon=180.0, helio_lat=0.0
    )
    pole = _pole_vector(40.589, 83.537)
    earth, sun = _sub_latitudes(pole, 0.0, 0.0, 0.0, 0.0, J2000)
    beta = math.degrees(math.sqrt(abs(earth * sun)))
    expected = (
        -8.914
        - 1.825 * math.sin(math.radians(beta))
        - 0.378 * math.sin(math.radians(beta))
    )
    assert first == pytest.approx(expected)
    assert mirrored == pytest.approx(first)


def test_saturn_rejects_ring_opening_at_27_degrees_without_tolerance() -> None:
    """The model refuses the first direction whose effective beta reaches 27°."""
    with pytest.raises(InputValidationError):
        _magnitude(SATURN, 0.0, geo_lon=65.0, helio_lon=65.0)
    assert math.isfinite(_magnitude(SATURN, 0.0, geo_lon=64.0, helio_lon=64.0))


def test_saturn_sub_latitudes_use_date_frame_and_turn_invariance() -> None:
    """Saturn pole latitudes are stable under complete longitude turns."""
    a = _direction_to_icrs(17.0, 12.0, J2000)
    b = _direction_to_icrs(377.0, 12.0, J2000)
    assert a == pytest.approx(b, abs=2e-16)


def test_uranus_pole_frame_and_brightening() -> None:
    """Uranus responds to independently derived pole latitudes."""
    base = _uranus_photometric_latitude(0.0, 0.0, 180.0, 0.0, J2000)
    poleward = _uranus_photometric_latitude(257.311, -15.175, 257.311, -15.175, J2000)
    assert poleward > base
    assert _uranus_photometric_latitude(360.0, 0.0, 540.0, 0.0, J2000) == pytest.approx(
        base
    )


@pytest.mark.parametrize(
    ("year", "expected"),
    [(1979, -6.89), (1980, -6.89), (2000, -6.998), (2001, -7.00)],
)
def test_neptune_calendar_year_boundaries(year: int, expected: float) -> None:
    """Neptune uses the published secular law at Gregorian boundaries."""
    jd = julday(year, 1, 1, 0.0)
    assert _calc_planet_magnitude(NEPTUNE, 1.0, 1.0, 0.0, tjd=jd) == pytest.approx(
        expected
    )


def test_neptune_public_helper_accepts_large_geocentric_phase_independently() -> None:
    """The public-backend helper does not require the rounded phase domain."""
    assert _calc_neptune_magnitude(1.0, 1.0, julday(1979, 1, 1)) == -6.89
    assert _calc_neptune_magnitude(2.0, 1.0, julday(1979, 1, 1)) == pytest.approx(
        -6.89 + 5.0 * math.log10(2.0)
    )


def test_neptune_leap_day_has_calendar_fraction() -> None:
    """Neptune's year fraction distinguishes leap-day from a 365-day year."""
    feb = _calc_planet_magnitude(NEPTUNE, 1.0, 1.0, 0.0, tjd=julday(1999, 2, 28))
    mar = _calc_planet_magnitude(NEPTUNE, 1.0, 1.0, 0.0, tjd=julday(1999, 3, 1))
    assert mar < feb
    assert feb - mar == pytest.approx(0.0054 / 365.0, rel=2e-8)
    leap = _calc_planet_magnitude(NEPTUNE, 1.0, 1.0, 0.0, tjd=julday(2000, 2, 29))
    after = _calc_planet_magnitude(NEPTUNE, 1.0, 1.0, 0.0, tjd=julday(2000, 3, 1))
    assert after == pytest.approx(leap, abs=0.0054 / 366.0)


def test_sun_published_value_and_distance_law() -> None:
    """The shared Sun helper uses Johnson V and the inverse-square law."""
    assert _sun_magnitude(1.0) == -26.76
    assert _sun_magnitude(2.0) - _sun_magnitude(1.0) == pytest.approx(
        5.0 * math.log10(2.0)
    )


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf"), "1"])
def test_sun_rejects_invalid_distance(bad: object) -> None:
    """Sun magnitudes reject invalid values before logarithms."""
    with pytest.raises(InputValidationError):
        _sun_magnitude(bad)


@pytest.mark.parametrize("body", [JUPITER, SATURN, URANUS, NEPTUNE])
def test_direct_model_rejects_large_phase(body: int) -> None:
    """Geocentric model domains are explicit and closed."""
    phase = {JUPITER: 12.0001, SATURN: 6.5, URANUS: 154.0001, NEPTUNE: 1.9001}[body]
    with pytest.raises(InputValidationError):
        _calc_planet_magnitude(body, 1.0, 1.0, phase, tjd=J2000)


def test_direct_model_rejects_invalid_inputs_and_unknown_body() -> None:
    """Invalid types, non-finite values, distances, and body IDs are typed errors."""
    with pytest.raises(UnknownBodyError):
        _calc_planet_magnitude(999, 1.0, 1.0, 0.0)
    for args in (
        ("1", 1.0, 1.0, 0.0),
        (MERCURY, 0.0, 1.0, 0.0),
        (MERCURY, 1.0, -1.0, 0.0),
        (MERCURY, 1.0, 1.0, float("nan")),
    ):
        with pytest.raises(InputValidationError):
            _calc_planet_magnitude(*args)
    with pytest.raises(InputValidationError):
        _calc_planet_magnitude(SATURN, 1.0, 1.0, 0.0, geo_lat=91.0, tjd=J2000)
