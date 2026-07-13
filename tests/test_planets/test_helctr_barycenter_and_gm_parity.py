"""Independent invariants for heliocentric and barycentric calculations.

These regressions intentionally contain no positions or orbital elements
captured from another ephemeris.  Expected values are mathematical identities
or are derived at runtime from the library's independently sourced JPL state.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as lib

_J2000 = 2451545.0


def _ipl(name: str) -> int:
    return getattr(lib, name)


def _signed_angle_delta(after: float, before: float) -> float:
    return (after - before + 180.0) % 360.0 - 180.0


def _norm(vector) -> float:
    return math.sqrt(sum(float(component) ** 2 for component in vector))


@pytest.mark.parametrize("name", ["MERCURY", "JUPITER", "SATURN", "PLUTO"])
def test_helctr_position_is_finite_spherical_state(name: str) -> None:
    pos, _ret = lib.calc(_J2000, _ipl(name), lib.FLG_SWIEPH | lib.FLG_HELCTR)
    assert all(math.isfinite(value) for value in pos)
    assert 0.0 <= pos[0] < 360.0
    assert -90.0 <= pos[1] <= 90.0
    assert pos[2] > 0.0


@pytest.mark.parametrize("name", ["MERCURY", "JUPITER", "PLUTO"])
def test_helctr_speed_is_derivative_of_position(name: str) -> None:
    flags = lib.FLG_SWIEPH | lib.FLG_HELCTR | lib.FLG_SPEED
    pos, _ret = lib.calc(_J2000, _ipl(name), flags)
    h = 1e-3
    p_lo, _ = lib.calc(_J2000 - h, _ipl(name), flags)
    p_hi, _ = lib.calc(_J2000 + h, _ipl(name), flags)
    fd = _signed_angle_delta(p_hi[0], p_lo[0]) / (2.0 * h)
    assert pos[3] == pytest.approx(fd, abs=0.02 / 3600.0)


@pytest.mark.parametrize("name", ["JUPITER", "SATURN", "URANUS", "NEPTUNE", "PLUTO"])
def test_baryctr_position_is_finite_spherical_state(name: str) -> None:
    pos, _ret = lib.calc(_J2000, _ipl(name), lib.FLG_SWIEPH | lib.FLG_BARYCTR)
    assert all(math.isfinite(value) for value in pos)
    assert 0.0 <= pos[0] < 360.0
    assert -90.0 <= pos[1] <= 90.0
    assert pos[2] > 0.0


def test_baryctr_light_time_is_small_and_nonzero() -> None:
    true_pos, _ = lib.calc(
        _J2000, lib.PLUTO, lib.FLG_SWIEPH | lib.FLG_BARYCTR | lib.FLG_TRUEPOS
    )
    apparent, _ = lib.calc(_J2000, lib.PLUTO, lib.FLG_SWIEPH | lib.FLG_BARYCTR)
    delta = abs(_signed_angle_delta(true_pos[0], apparent[0]))
    assert 0.0 < delta < 1.0


def test_jupiter_true_positions_use_a_covered_planet_center() -> None:
    """HELCTR and BARYCTR retain a covered JPL physical center state."""
    from libephemeris.planets import get_planet_target
    from libephemeris.state import (
        get_planet_center_segment,
        get_planets,
        get_timescale,
    )

    previous_mode = lib.get_calc_mode()
    lib.set_calc_mode("skyfield")
    try:
        planets = get_planets()
        segment = get_planet_center_segment(599)
        if segment is None:
            pytest.skip("Jupiter planet-center SPK is not installed")
        spk_segment = segment.spk_segment
        jd_tt = 0.5 * (spk_segment.start_jd + spk_segment.end_jd)
        epoch = get_timescale().tt_jd(jd_tt)
        center = get_planet_target(planets, "jupiter").at(epoch).position.au
        barycenter = planets["jupiter barycenter"].at(epoch).position.au
        sun = planets["sun"].at(epoch).position.au

        helio_expected = _norm(center - sun)
        helio_barycenter = _norm(barycenter - sun)
        bary_expected = _norm(center)
        bary_barycenter = _norm(barycenter)

        helio, _ = lib.calc(
            jd_tt,
            lib.JUPITER,
            lib.FLG_SWIEPH | lib.FLG_HELCTR | lib.FLG_TRUEPOS,
        )
        bary, _ = lib.calc(
            jd_tt,
            lib.JUPITER,
            lib.FLG_SWIEPH | lib.FLG_BARYCTR | lib.FLG_TRUEPOS,
        )
        assert helio[2] == pytest.approx(helio_expected, abs=2e-12)
        assert bary[2] == pytest.approx(bary_expected, abs=2e-12)
        assert abs(helio[2] - helio_barycenter) > 1e-10
        assert abs(bary[2] - bary_barycenter) > 1e-10
    finally:
        lib.set_calc_mode(previous_mode)


def test_mercury_helctr_light_time_uses_sun_to_target_distance() -> None:
    """The retarded epoch solves |Mercury(t-lt)-Sun(t)| = c*lt."""
    from libephemeris.fast_calc import C_LIGHT_AU_DAY
    from libephemeris.planets import get_planet_target
    from libephemeris.state import get_planets, get_timescale

    previous_mode = lib.get_calc_mode()
    lib.set_calc_mode("skyfield")
    try:
        planets = get_planets()
        timescale = get_timescale()
        epoch = timescale.tt_jd(_J2000)
        target = get_planet_target(planets, "mercury")
        observer = planets["sun"].at(epoch).position.au
        relative = target.at(epoch).position.au - observer
        for _ in range(3):
            light_time = _norm(relative) / C_LIGHT_AU_DAY
            emission = timescale.tdb_jd(epoch.whole, epoch.tdb_fraction - light_time)
            relative = target.at(emission).position.au - observer

        apparent, _ = lib.calc(_J2000, lib.MERCURY, lib.FLG_SWIEPH | lib.FLG_HELCTR)
        assert apparent[2] == pytest.approx(_norm(relative), abs=2e-12)
    finally:
        lib.set_calc_mode(previous_mode)


@pytest.mark.parametrize("name", ["JUPITER", "SATURN", "PLUTO"])
def test_baryctr_speed_is_derivative_of_position(name: str) -> None:
    flags = lib.FLG_SWIEPH | lib.FLG_BARYCTR | lib.FLG_SPEED
    pos, _ret = lib.calc(_J2000, _ipl(name), flags)
    h = 1e-3
    p_lo, _ = lib.calc(_J2000 - h, _ipl(name), flags)
    p_hi, _ = lib.calc(_J2000 + h, _ipl(name), flags)
    fd = _signed_angle_delta(p_hi[0], p_lo[0]) / (2.0 * h)
    assert pos[3] == pytest.approx(fd, abs=0.02 / 3600.0)


def test_sun_heliocentric_stays_origin() -> None:
    pos, _ret = lib.calc(_J2000, lib.SUN, lib.FLG_SWIEPH | lib.FLG_HELCTR)
    assert pos[2] == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("name", ["JUPITER", "SATURN", "NEPTUNE", "EARTH"])
def test_osculating_distance_identities(name: str) -> None:
    elements = lib.get_orbital_elements(
        _J2000, _ipl(name), lib.FLG_SWIEPH | lib.FLG_HELCTR
    )
    a, eccentricity, perihelion, aphelion = (
        elements[0],
        elements[1],
        elements[15],
        elements[16],
    )
    assert 0.0 <= eccentricity < 1.0
    assert perihelion == pytest.approx(a * (1.0 - eccentricity), rel=1e-12)
    assert aphelion == pytest.approx(a * (1.0 + eccentricity), rel=1e-12)
    assert perihelion < a < aphelion
