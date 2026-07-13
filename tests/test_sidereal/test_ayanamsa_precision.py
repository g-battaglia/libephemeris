"""Precision invariants for restored predefined ayanamsha definitions."""

from __future__ import annotations

import math
import warnings

import erfa
import pytest

import libephemeris as ephem
from libephemeris.ayanamsha_definitions import AYANAMSHA_DEFINING
from libephemeris.constants import (
    SIDM_FAGAN_BRADLEY,
    SIDM_J2000,
    SIDM_KRISHNAMURTI,
    SIDM_LAHIRI,
    SIDM_RAMAN,
)
from libephemeris.precession_vondrak import method_b_accumulated_precession


J2000 = 2451545.0
TEST_DATES_TT = (2415020.0, J2000, 2488070.0)
FORMULA_SAMPLE_MODES = (
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
    SIDM_RAMAN,
    SIDM_KRISHNAMURTI,
)


def _signed_angle(angle: float) -> float:
    return (angle + 180.0) % 360.0 - 180.0


@pytest.mark.parametrize("sid_mode", FORMULA_SAMPLE_MODES)
@pytest.mark.parametrize("jd_tt", TEST_DATES_TT)
def test_formula_modes_use_their_shared_method_b_definition(
    sid_mode: int, jd_tt: float
) -> None:
    ephem.set_sid_mode(sid_mode)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        actual = ephem.get_ayanamsa(jd_tt)
    defining_value, defining_epoch = AYANAMSHA_DEFINING[sid_mode]
    expected = defining_value + method_b_accumulated_precession(jd_tt, defining_epoch)
    assert _signed_angle(actual - expected) == pytest.approx(0.0, abs=2e-9)
    assert not caught


def test_fagan_bradley_recovers_the_published_b1950_svp_anchor() -> None:
    """Fagan and Firebrace publish SVP 335d57m28.64s at B1950."""
    ephem.set_sid_mode(SIDM_FAGAN_BRADLEY)
    b1950_tt = float(sum(erfa.epb2jd(1950.0)))
    actual = ephem.get_ayanamsa(b1950_tt)
    published = 24.0 + 2.0 / 60.0 + 31.36 / 3600.0
    assert _signed_angle(actual - published) == pytest.approx(0.0, abs=1.5e-4)


def test_raman_matches_his_published_fixed_rate_rule_near_1900() -> None:
    """Raman's Appendix A gives (year - 397) times 50 1/3 arcsec."""
    ephem.set_sid_mode(SIDM_RAMAN)
    j1900_tt = float(sum(erfa.epj2jd(1900.0)))
    actual = ephem.get_ayanamsa(j1900_tt)
    published = (1900.0 - 397.0) * (151.0 / 3.0) / 3600.0
    assert _signed_angle(actual - published) == pytest.approx(0.0, abs=5e-5)


@pytest.mark.parametrize("jd_tt", TEST_DATES_TT)
def test_j2000_mode_uses_vondrak_method_b(jd_tt: float) -> None:
    ephem.set_sid_mode(SIDM_J2000)
    actual = _signed_angle(ephem.get_ayanamsa(jd_tt))
    expected = method_b_accumulated_precession(jd_tt, J2000)
    assert actual == pytest.approx(expected, abs=2e-9)


def test_precession_rate_increases_across_j2000() -> None:
    past = method_b_accumulated_precession(2415020.0, J2000)
    present = method_b_accumulated_precession(J2000, J2000)
    future = method_b_accumulated_precession(2488070.0, J2000)
    assert math.isfinite(past)
    assert future - present > present - past
