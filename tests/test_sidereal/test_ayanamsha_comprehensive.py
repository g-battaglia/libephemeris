"""Reference-free coverage for every predefined ayanamsha mode."""

from __future__ import annotations

from collections.abc import Iterator
import math
import warnings

import erfa
import pytest

import libephemeris as ephem


NATIVE_MODES = frozenset(range(47))
TEST_DATES = (2415020.0, 2451545.0, 2488070.0)


@pytest.fixture(autouse=True)
def _reset_state() -> Iterator[None]:
    ephem.reset_session()
    yield
    ephem.reset_session()


def _signed_angle(angle: float) -> float:
    return (angle + 180.0) % 360.0 - 180.0


@pytest.mark.parametrize("mode", range(47))
def test_every_predefined_mode_returns_a_native_finite_float(mode: int) -> None:
    ephem.set_sid_mode(mode)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        value = ephem.get_ayanamsa_ut(2451545.0)
    assert type(value) is float
    assert math.isfinite(value)
    assert 0.0 <= value < 360.0
    assert not caught


@pytest.mark.parametrize("mode", sorted(NATIVE_MODES))
@pytest.mark.parametrize("jd_ut", TEST_DATES)
def test_predefined_modes_do_not_warn_or_collapse_to_j2000(
    mode: int, jd_ut: float
) -> None:
    ephem.set_sid_mode(mode)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        value = ephem.get_ayanamsa_ut(jd_ut)
    assert math.isfinite(value)
    assert not caught
    if mode not in {ephem.SIDM_J2000}:
        ephem.set_sid_mode(ephem.SIDM_J2000)
        j2000 = ephem.get_ayanamsa_ut(jd_ut)
        assert abs(_signed_angle(value - j2000)) > 1e-4


@pytest.mark.parametrize(
    ("mode", "epoch"),
    [
        (ephem.SIDM_J2000, float(sum(erfa.epj2jd(2000.0)))),
        (ephem.SIDM_J1900, float(sum(erfa.epj2jd(1900.0)))),
        (ephem.SIDM_B1950, float(sum(erfa.epb2jd(1950.0)))),
    ],
)
def test_standard_epoch_modes_are_zero_at_erfa_epoch(mode: int, epoch: float) -> None:
    ephem.set_sid_mode(mode)
    value = _signed_angle(ephem.get_ayanamsa(epoch))
    assert value == pytest.approx(0.0, abs=2e-7)


@pytest.mark.parametrize("mode", sorted(NATIVE_MODES))
def test_native_modes_do_not_emit_a_warning(mode: int) -> None:
    ephem.set_sid_mode(mode)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        value = ephem.get_ayanamsa_ut(2451545.0)
    assert math.isfinite(value)
    assert not caught


@pytest.mark.parametrize("mode", [0, 17, 18, 34])
def test_sidereal_longitude_uses_the_selected_mean_ayanamsha(mode: int) -> None:
    ephem.set_sid_mode(mode)
    jd_ut = 2451545.0
    tropical, _ = ephem.calc_ut(jd_ut, ephem.SUN, ephem.FLG_NONUT)
    sidereal, _ = ephem.calc_ut(jd_ut, ephem.SUN, ephem.FLG_SIDEREAL | ephem.FLG_NONUT)
    ayanamsha = ephem.get_ayanamsa_ut(jd_ut)
    difference = _signed_angle((tropical[0] - sidereal[0]) - ayanamsha)
    assert difference == pytest.approx(0.0, abs=5e-10)


@pytest.mark.parametrize("mode", range(47))
def test_predefined_modes_have_a_nonempty_name(mode: int) -> None:
    """Every predefined mode 0-46 carries a non-empty descriptive name."""
    assert ephem.get_ayanamsa_name(mode) != ""


@pytest.mark.parametrize("sidmode", [47, 48, 100, 254, ephem.SIDM_USER])
def test_modes_without_a_predefined_name_return_empty_string(sidmode: int) -> None:
    """Measured reference behavior: an id without a predefined name -- the
    unassigned block above mode 46 and the user-defined mode SIDM_USER (255) --
    returns the empty string, not a placeholder such as 'Unknown'."""
    assert ephem.get_ayanamsa_name(sidmode) == ""
