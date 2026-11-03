"""Operation-specific public exception classes for issues #31–#33."""
from data_paths import CORE_LEB, DATA_ROOT, LEB_HASHES
import hashlib
import math
from pathlib import Path

import libephemeris as eph
from libephemeris import planetary_moons
from libephemeris.moon_theories import tass17
import pytest

ROOT = Path(__file__).resolve().parents[2]
BASE = CORE_LEB
MEDIUM = DATA_ROOT / 'leb/medium_core.leb2'


@pytest.fixture(autouse=True)
def sealed_stored_configuration():
    eph.set_network_policy('sealed')
    eph.set_auto_spk_download(False)
    eph.set_calc_mode('leb')
    eph.set_precision_tier('base')
    eph.set_leb_file(str(BASE))
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)
    planetary_moons.close_moon_kernels()
    yield
    planetary_moons.close_moon_kernels()
    eph.set_calc_mode('leb')
    eph.set_precision_tier('base')
    eph.set_leb_file(str(BASE))


@pytest.mark.parametrize('function', ['saturn_moon_position', 'saturn_moon_position_velocity'])
@pytest.mark.parametrize('body', [-1, 8])
def test_tass_invalid_indices_keep_the_builtin_value_error(function, body):
    with pytest.raises(ValueError) as caught:
        getattr(tass17, function)(2460311.0, body=body)
    assert type(caught.value) is ValueError
    assert not issubclass(eph.InvalidBodyError, ValueError)


@pytest.mark.parametrize('function,size', [('saturn_moon_position', 3), ('saturn_moon_position_velocity', 6)])
@pytest.mark.parametrize('body', [0, 7])
def test_tass_valid_endpoints_remain_real_model_results(function, size, body):
    actual = getattr(tass17, function)(2460311.0, body=body)
    assert isinstance(actual, tuple) and len(actual) == size
    assert all(math.isfinite(value) for value in actual)
    assert any(value != 0.0 for value in actual[:3])


@pytest.mark.parametrize('body', [9001, 9002, 9003, 9004])
@pytest.mark.parametrize('function', ['calc', 'calc_ut'])
@pytest.mark.parametrize('independent_context', [False, True])
def test_known_galilean_forced_leb_refusal_is_exact_base_error(body, function, independent_context):
    target = eph.EphemerisContext() if independent_context else eph
    if independent_context:
        target.set_leb_file(str(BASE))
    with pytest.raises(eph.Error) as caught:
        getattr(target, function)(2451545.0, body, eph.FLG_SPEED)
    assert type(caught.value) is eph.Error


@pytest.mark.parametrize('function', ['calc', 'calc_ut'])
def test_unknown_body_remains_distinct_from_unavailable_known_moon(function):
    with pytest.raises(eph.UnknownBodyError) as caught:
        getattr(eph, function)(2451545.0, -999999, 2)
    assert type(caught.value) is eph.UnknownBodyError


@pytest.mark.parametrize('function,jd,kwargs', [
    ('planet_occult_when_glob', 2458849.5, {}),
    ('planet_occult_when_loc', 2460310.5, {'lat': 45.0, 'lon': 10.0}),
])
def test_planetary_occultation_exhaustion_is_exact_base_error(function, jd, kwargs):
    if not MEDIUM.is_file():
        pytest.skip('optional medium ephemeris is unavailable')
    assert hashlib.sha256(MEDIUM.read_bytes()).hexdigest() == LEB_HASHES['medium']['core']
    eph.set_precision_tier('medium')
    eph.set_leb_file(str(MEDIUM))
    with pytest.raises(eph.Error) as caught:
        getattr(eph, function)(jd, eph.VENUS, eph.MARS, **kwargs)
    assert type(caught.value) is eph.Error


def test_event_argument_failure_does_not_become_search_exhaustion():
    with pytest.raises(ValueError) as caught:
        eph.planet_occult_when_glob(2458849.5, eph.VENUS, eph.VENUS)
    assert type(caught.value) is ValueError
