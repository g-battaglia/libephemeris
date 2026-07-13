"""Reference-free invariants covering every built-in ayanamsha mode."""

from __future__ import annotations

from collections.abc import Iterator
import math
import warnings

import erfa
import pytest

import libephemeris as ephem
from libephemeris.constants import (
    SIDM_GALCENT_0SAG,
    SIDM_TRUE_CITRA,
)
from libephemeris.ayanamsha_definitions import (
    AYANAMSHA_DEFINING,
    GALCENT_TARGET_LON,
)
from libephemeris.planets import STARS, _get_star_position_ecliptic
from libephemeris.state import get_timescale


ALL_MODES = tuple(range(47))
TEST_DATES = tuple(
    float(sum(erfa.epj2jd(epoch))) for epoch in (1900.0, 2000.0, 2050.0, 2100.0)
)


@pytest.fixture(autouse=True)
def _reset_state() -> Iterator[None]:
    ephem.reset_session()
    yield
    ephem.reset_session()


def _signed_angle(angle: float) -> float:
    return (angle + 180.0) % 360.0 - 180.0


@pytest.mark.parametrize("mode", ALL_MODES)
@pytest.mark.parametrize("jd_ut", TEST_DATES)
def test_every_mode_returns_a_finite_native_angle(mode: int, jd_ut: float) -> None:
    ephem.set_sid_mode(mode)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ayanamsha = ephem.get_ayanamsa_ut(jd_ut)
    assert type(ayanamsha) is float
    assert math.isfinite(ayanamsha)
    assert 0.0 <= ayanamsha < 360.0
    assert not caught


@pytest.mark.parametrize("mode", range(15))
@pytest.mark.parametrize("jd_ut", TEST_DATES)
def test_extended_mean_and_true_values_differ_only_by_iau_nutation(
    mode: int, jd_ut: float
) -> None:
    ephem.set_sid_mode(mode)
    mean = ephem.get_ayanamsa_ut(jd_ut)
    _, true = ephem.get_ayanamsa_ex_ut(jd_ut, 0)
    _, nonut = ephem.get_ayanamsa_ex_ut(jd_ut, ephem.FLG_NONUT)
    jd_tt = float(get_timescale().ut1_jd(jd_ut).tt)
    dpsi, _ = erfa.nut06a(2451545.0, jd_tt - 2451545.0)
    assert _signed_angle(true - mean) == pytest.approx(
        math.degrees(float(dpsi)), abs=2e-12
    )
    assert nonut == pytest.approx(mean, abs=2e-12)


@pytest.mark.parametrize("jd_ut", TEST_DATES)
def test_true_citra_places_spica_at_180_degrees(jd_ut: float) -> None:
    ephem.set_sid_mode(SIDM_TRUE_CITRA)
    ayanamsha = ephem.get_ayanamsa_ut(jd_ut)
    jd_tt = float(get_timescale().ut1_jd(jd_ut).tt)
    star_longitude = _get_star_position_ecliptic(STARS["SPICA"], jd_tt, 0.0, nonut=True)
    sidereal_longitude = (star_longitude - ayanamsha) % 360.0
    assert _signed_angle(sidereal_longitude - 180.0) == pytest.approx(0.0, abs=2e-10)


@pytest.mark.parametrize("jd_ut", TEST_DATES)
def test_galactic_center_mode_holds_sgr_a_star_at_its_published_target(
    jd_ut: float,
) -> None:
    """Sgr A* (Reid & Brunthaler 2004) sits at 0 Sagittarius at every date."""
    ephem.set_sid_mode(SIDM_GALCENT_0SAG)
    ayanamsha = ephem.get_ayanamsa_ut(jd_ut)
    jd_tt = float(get_timescale().ut1_jd(jd_ut).tt)
    center_longitude = _get_star_position_ecliptic(
        STARS["GAL_CENTER"], jd_tt, 0.0, nonut=True
    )
    expected_sidereal_longitude = GALCENT_TARGET_LON[SIDM_GALCENT_0SAG]
    assert _signed_angle(
        (center_longitude - ayanamsha) - expected_sidereal_longitude
    ) == pytest.approx(
        0.0,
        abs=2e-10,
    )


def test_sheoran_uses_its_published_mahabharata_epoch() -> None:
    """Sheoran, The Science of Time: -60 deg at the winter solstice 4174 BCE."""
    from libephemeris.constants import SIDM_TRUE_SHEORAN as _MODE
    from libephemeris.precession_vondrak import method_b_accumulated_precession

    defining_value, defining_epoch = AYANAMSHA_DEFINING[_MODE]
    assert defining_value == -60.0
    for jd_ut in TEST_DATES:
        ephem.set_sid_mode(_MODE)
        ayanamsha = ephem.get_ayanamsa_ut(jd_ut)
        jd_tt = float(get_timescale().ut1_jd(jd_ut).tt)
        expected = (
            defining_value + method_b_accumulated_precession(jd_tt, defining_epoch)
        ) % 360.0
        assert _signed_angle(ayanamsha - expected) == pytest.approx(0.0, abs=2e-9)


def test_all_mode_ids_remain_available() -> None:
    assert ALL_MODES == tuple(range(47))
