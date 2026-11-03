"""Installed Python calls must actually use the selected native LEB backend."""
from data_paths import CORE_LEB
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = CORE_LEB


def calculate(tmp_path, reader, limits=''):
    config = tmp_path / 'reader.toml'
    config.write_text("[libephemeris]\nleb_reader='indexed'\n" + limits)
    environment = {key: value for key, value in os.environ.items() if not key.startswith('LIBEPHEMERIS_')}
    environment.update(LIBEPHEMERIS_CONFIG=str(config), LIBEPHEMERIS_ENV_FILE=str(tmp_path / 'absent.env'),
                       LIBEPHEMERIS_NETWORK_POLICY='sealed', LIBEPHEMERIS_IERS_DELTA_T='0',
                       LIBEPHEMERIS_IERS_AUTO_DOWNLOAD='0', LIBEPHEMERIS_AUTO_SPK='0',
                       LIBEPHEMERIS_DATA_DIR=str(tmp_path), LIBEPHEMERIS_MODE='leb')
    if reader is not None:
        environment['LIBEPHEMERIS_LEB_READER'] = reader
    script = '''
import json,struct,sys
import libephemeris as e
try:
    e.set_leb_file(sys.argv[1])
    results=[]
    for body in [0,1,2,4,9]:
        position, flags=e.calc_ut(2451545.0,body,e.FLG_SPEED)
        results.append(([struct.pack('!d',value).hex() for value in position],flags))
except Exception as error:
    print(json.dumps({'error':type(error).__name__,'message':str(error)}))
else:
    print(json.dumps({'result':results}))
'''
    completed = subprocess.run([sys.executable, '-I', '-c', script, str(FIXTURE)], env=environment,
                               cwd=tmp_path, text=True, capture_output=True, check=True, timeout=60)
    return json.loads(completed.stdout)


def test_indexed_and_owned_public_results_have_identical_binary64_words(tmp_path):
    owned = calculate(tmp_path, 'owned')
    indexed = calculate(tmp_path, 'indexed')
    assert 'result' in owned
    assert indexed == owned


def test_configured_limits_prove_the_indexed_backend_is_really_selected(tmp_path):
    limits = '[libephemeris.leb_limits]\nframe_bytes=0\n'
    result = calculate(tmp_path, 'indexed', limits)
    assert result['error'] == 'ConfigurationError'
    assert 'LEB decoded frame' in result['message']
    # Environment wins over TOML. Owned ignores indexed-only allocation limits.
    assert 'result' in calculate(tmp_path, 'owned', limits)
    assert calculate(tmp_path, None, limits)['error'] == 'ConfigurationError'


def test_invalid_reader_is_not_silently_replaced_by_owned(tmp_path):
    result = calculate(tmp_path, 'typo')
    assert result['error'] == 'ConfigurationError'
    assert 'leb_reader' in result['message']
