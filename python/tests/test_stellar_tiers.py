"""Historical stellar observers and live anchors use the actual date-aware tier."""
from pathlib import Path
import hashlib
import os
import subprocess
import sys

import pytest
from data_paths import DATA_ROOT, LEB_HASHES

ROOT = Path(__file__).resolve().parents[2]
DATA = DATA_ROOT
CORE_HASHES = {tier: groups['core'] for tier, groups in LEB_HASHES.items()}


@pytest.fixture(scope='module')
def mixed_stellar_data(tmp_path_factory):
    root = tmp_path_factory.mktemp('stellar-tiers')
    (root / 'leb').mkdir()
    for tier, expected in CORE_HASHES.items():
        source = DATA / f'leb/{tier}_core.leb2'
        if not source.is_file():
            pytest.skip('the read-only three-tier astronomy dataset is unavailable')
        with source.open('rb') as stream:
            assert hashlib.file_digest(stream, 'sha256').hexdigest() == expected
        (root / 'leb' / source.name).symlink_to(source)
    return root


@pytest.mark.parametrize('mode', ['leb', 'auto'])
def test_historical_stellar_observer_matches_the_selected_explicit_image(mixed_stellar_data, mode):
    environment = {key: value for key, value in os.environ.items() if not key.startswith('LIBEPHEMERIS_')}
    environment.update(
        LIBEPHEMERIS_DATA_DIR=str(mixed_stellar_data),
        LIBEPHEMERIS_CONFIG=str(mixed_stellar_data / 'absent.toml'),
        LIBEPHEMERIS_ENV_FILE=str(mixed_stellar_data / 'absent.env'),
        LIBEPHEMERIS_MODE=mode, LIBEPHEMERIS_PRECISION='extended',
        LIBEPHEMERIS_NETWORK_POLICY='sealed', LIBEPHEMERIS_AUTO_SPK='0',
        LIBEPHEMERIS_IERS_DELTA_T='0', LIBEPHEMERIS_IERS_AUTO_DOWNLOAD='0',
    )
    process = subprocess.run([sys.executable, '-c', r'''
from pathlib import Path
import math
import sys
import libephemeris as e
root, external, mode = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
e.set_calc_mode(mode)
e.set_precision_tier('extended')
e.set_ephe_path(str(external))
e.set_ephemeris_file('de441.bsp')
e.set_auto_spk_download(False)
e.set_iers_delta_t_enabled(False)
e.set_delta_t_userdef(64.0 / 86400.0)
for year, tier in [(-500, 'extended'), (1000, 'extended'), (1700, 'medium'), (2000, 'base')]:
    jd = e.julday(year, 6, 15, 12.0)
    e.set_leb_file(None)
    e.set_sid_mode(0)
    coverage = e.get_body_coverage(e.EARTH, jd)
    assert coverage is not None and coverage.contains(jd)
    assert Path(coverage.data_file).name == f'{tier}_core.leb2', coverage.data_file
    stars = [e.fixstar('Regulus', jd, 258), e.fixstar_ut('Regulus', jd, 258),
             *e.batch_fixstars_ut(['Regulus', 'Aldebaran'], jd, 258)]
    for value, _, _ in stars:
        assert all(math.isfinite(component) for component in value)
    e.set_sid_mode(e.SIDM_ALDEBARAN_15TAU)
    anchors = [e.calc(jd, e.SUN, 258 | e.FLG_SIDEREAL),
               e.calc_ut(jd, e.SUN, 258 | e.FLG_SIDEREAL)]
    e.set_leb_file(str(root / 'leb' / f'{tier}_core.leb2'))
    e.set_sid_mode(0)
    explicit = [e.fixstar('Regulus', jd, 258), e.fixstar_ut('Regulus', jd, 258),
                *e.batch_fixstars_ut(['Regulus', 'Aldebaran'], jd, 258)]
    # This is a provider-selection invariant using identical actual images,
    # not a fitted expected coordinate or a replacement for v3 comparisons.
    assert stars == explicit, (year, mode, 'stellar observer selected a different image')
    e.set_sid_mode(e.SIDM_ALDEBARAN_15TAU)
    explicit_anchors = [e.calc(jd, e.SUN, 258 | e.FLG_SIDEREAL),
                        e.calc_ut(jd, e.SUN, 258 | e.FLG_SIDEREAL)]
    for actual, reference in zip(anchors, explicit_anchors):
        assert actual[1] == reference[1]
        for index, (got, wanted) in enumerate(zip(actual[0], reference[0])):
            delta = abs((got - wanted + 180.0) % 360.0 - 180.0) if index == 0 else abs(got - wanted)
            limit = 0.001 / 3600.0 if index < 2 else (
                abs(wanted) * 1e-8 if index == 2 else
                1e-6 if index < 5 and abs(wanted) <= 1e-6 else abs(wanted) * 1e-6)
            assert delta <= limit, (year, mode, index, delta, limit)
''', str(mixed_stellar_data), str(DATA), mode], env=environment, capture_output=True, text=True, timeout=120)
    assert process.returncode == 0, process.stdout + process.stderr
