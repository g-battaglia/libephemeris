"""Public calc dispatch uses native catalog, observer and context cache models."""
from data_paths import CORE_LEB
from pathlib import Path
import math
import libephemeris as eph
import pytest

BASE = CORE_LEB


@pytest.fixture
def configured():
    eph.set_network_policy('sealed')
    eph.set_auto_spk_download(False)
    eph.set_calc_mode('leb')
    eph.set_precision_tier('base')
    eph.set_leb_file(str(BASE))
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)
    eph.set_sid_mode(eph.SIDM_LAHIRI)
    eph.set_topo(12.4964, 41.9028, 0.)
    eph.clear_angles_cache()
    return eph


@pytest.mark.parametrize('function', ['calc', 'calc_ut'])
@pytest.mark.parametrize('body,name', [(1000001, ',alLeo'), (1000002, ',alVir'), (1000118, ',29Psc')])
def test_numeric_stars_share_the_stellar_engine(configured, function, body, name):
    flags_list = [0, eph.FLG_SPEED, eph.FLG_SPEED | eph.FLG_EQUATORIAL,
                  eph.FLG_SPEED | eph.FLG_XYZ | eph.FLG_J2000,
                  eph.FLG_SPEED | eph.FLG_SIDEREAL | eph.FLG_RADIANS,
                  eph.FLG_SPEED | eph.FLG_TOPOCTR]
    star = eph.fixstar_ut if function == 'calc_ut' else eph.fixstar
    for flags in flags_list:
        actual, returned = getattr(eph, function)(2451545., body, flags)
        expected, _, _ = star(name, 2451545., flags)
        assert actual == expected
        expected_flags = flags | eph.FLG_SWIEPH if function == 'calc_ut' else flags
        if flags & (eph.FLG_J2000 | eph.FLG_SIDEREAL):
            expected_flags |= eph.FLG_NONUT
        assert returned == expected_flags


@pytest.mark.parametrize('function', ['calc', 'calc_ut'])
def test_angle_identifiers_share_the_observer_clock(configured, function):
    jd = 2451545.
    # UT entries share the explicit angle helper. TT clock conversion is
    # covered independently by the matched pinned-Source comparison grid.
    angles = eph.calc_angles(jd, 41.9028, 12.4964) if function == 'calc_ut' else None
    results = []
    for identifier, name in [(1, 'Ascendant'), (2, 'MC'), (3, 'Descendant'),
                             (4, 'IC'), (5, 'Vertex'), (6, 'AntiVertex')]:
        value, returned = getattr(eph, function)(jd, eph.ANGLE_OFFSET + identifier, 0)
        if angles is not None:
            assert abs(value[0] - angles[name]) <= .001 / 3600.
        assert value[1:] == (0.,) * 5
        assert returned == (eph.FLG_SWIEPH if function == 'calc_ut' else 0)
        results.append(value[0])
    assert abs((results[2] - results[0]) % 360. - 180.) < 1e-12
    assert abs((results[3] - results[1]) % 360. - 180.) < 1e-12
    independent = eph.EphemerisContext()
    with pytest.raises(eph.ConfigurationError) as caught:
        getattr(independent, function)(jd, eph.ASCENDANT, 0)
    assert caught.value.missing_config == 'observer_location'


@pytest.mark.parametrize('function', ['calc', 'calc_ut'])
@pytest.mark.parametrize('sun', [90., 270.])
def test_cached_part_identifiers_use_the_owning_context(configured, function, sun):
    cache = dict(Asc=15., Sun=sun, Moon=240., Mercury=130., Venus=110.)
    eph.set_angles_cache(cache)
    expected = eph.calc_all_arabic_parts(cache)
    ids = [eph.PARS_FORTUNAE, eph.PARS_SPIRITUS, eph.PARS_AMORIS, eph.PARS_FIDEI]
    for identifier, name in zip(ids, ['Pars_Fortunae', 'Pars_Spiritus', 'Pars_Amoris', 'Pars_Fidei']):
        value, returned = getattr(eph, function)(2451545., identifier, 0)
        assert value == (expected[name], 0., 0., 0., 0., 0.)
        assert returned == (eph.FLG_SWIEPH if function == 'calc_ut' else 0)
        context = eph.EphemerisContext()
        with pytest.raises(eph.ConfigurationError) as caught:
            getattr(context, function)(2451545., identifier, 0)
        assert caught.value.missing_config == 'angles_cache'
        context.set_angles_cache(dict(cache, Asc=27.))
        own = getattr(context, function)(2451545., identifier, 0)[0]
        assert abs((own[0] - value[0]) % 360. - 12.) < 1e-12
        assert eph.get_angles_cache() == cache


@pytest.mark.parametrize('function', ['calc', 'calc_ut'])
def test_sidereal_output_units_and_frames_are_finite(configured, function):
    base = eph.FLG_SIDEREAL | eph.FLG_SPEED
    for mode in [eph.SIDM_LAHIRI, eph.SIDM_J2000, eph.SIDM_LAHIRI | eph.SIDBIT_ECL_T0,
                 eph.SIDM_LAHIRI | eph.SIDBIT_SSY_PLANE]:
        eph.set_sid_mode(mode)
        for body in [eph.SUN, eph.MOON, eph.MARS, eph.EARTH]:
            for frame in [0, eph.FLG_EQUATORIAL, eph.FLG_J2000 | eph.FLG_EQUATORIAL]:
                flags = base | frame
                degrees = getattr(eph, function)(2451545., body, flags)[0]
                radians = getattr(eph, function)(2451545., body, flags | eph.FLG_RADIANS)[0]
                xyz = getattr(eph, function)(2451545., body, flags | eph.FLG_XYZ)[0]
                assert all(math.isfinite(x) for x in (*degrees, *radians, *xyz))
                for index in [0, 1, 3, 4]:
                    assert math.isclose(radians[index], math.radians(degrees[index]), rel_tol=1e-14, abs_tol=1e-14)
                assert radians[2] == degrees[2] and radians[5] == degrees[5]
                assert math.isclose(math.hypot(*xyz[:3]), degrees[2], rel_tol=1e-14, abs_tol=1e-15)
