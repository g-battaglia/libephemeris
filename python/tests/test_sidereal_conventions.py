"""Declared epoch pairs and live anchors define different sidereal models."""
from data_paths import CORE_LEB
import math
from pathlib import Path

import libephemeris as eph
from libephemeris.ayanamsha_definitions import AYANAMSHA_DEFINING
import pytest


@pytest.fixture(autouse=True)
def configured():
    eph.set_calc_mode('leb')
    eph.set_precision_tier('base')
    eph.set_leb_file(str(CORE_LEB))
    eph.set_delta_t_userdef(64. / 86400.)
    eph.set_iers_delta_t_enabled(False)
    eph.set_network_policy('sealed')
    yield
    eph.set_sid_mode(eph.SIDM_FAGAN_BRADLEY)


@pytest.mark.parametrize('mode', sorted(AYANAMSHA_DEFINING))
def test_named_pair_equals_the_independent_user_realization(mode):
    offset, epoch = AYANAMSHA_DEFINING[mode]
    eph.set_sid_mode(mode)
    assert abs(math.remainder(eph.get_ayanamsa(epoch)-offset, 360.))*3600 <= .001
    for tt in (2396754.5, 2451545., 2460310.5):
        eph.set_sid_mode(mode)
        named = eph.get_ayanamsa(tt)
        eph.set_sid_mode(eph.SIDM_USER, epoch, offset)
        assert named == eph.get_ayanamsa(tt)


@pytest.mark.parametrize('tt', (2447892.5, 2451545., 2460310.5))
def test_aldebaran_is_a_live_direction_constraint(tt):
    eph.set_sid_mode(eph.SIDM_ALDEBARAN_15TAU)
    position, _, _ = eph.fixstar('Aldebaran', tt, eph.FLG_SWIEPH | eph.FLG_SIDEREAL)
    assert abs(math.remainder(position[0]-45., 360.))*3600 <= .001
