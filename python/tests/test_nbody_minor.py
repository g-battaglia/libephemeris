"""Installed Runtime minor-body routing and cache/resource ownership."""

import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys

import pytest
from data_paths import DATA_ROOT
# Fixed identity of the optional local Chiron kernel, not an import from
# the external qualification suite.
CHIRON_SPK_SHA256 = 'd8d93784dbaf078cfdddac4129a88824fa7057fdbb600895da9c8e02dd150b23'


DATA = DATA_ROOT
HASHES = {
    DATA / "de440.bsp": "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5",
    DATA / "assist/linux_p1550p2650.440": "29915576d0a6555766b99485ac3056ee415e86df4fce282611c31afb329ad062",
    DATA / "assist/sb441-n16.bsp": "919d612ce3c72a78fc7158f9120156542d0f21e6b8b052e4c1339c759747fd90",
    DATA / "spk/2060_184911_215002.bsp": CHIRON_SPK_SHA256,
}


@pytest.fixture(scope="module")
def verified_files():
    for path, expected in HASHES.items():
        if not path.is_file():
            pytest.skip("optional local JPL/ASSIST files are unavailable")
        with path.open("rb") as stream:
            assert hashlib.file_digest(stream, "sha256").hexdigest() == expected
    return DATA


SCRIPT = r'''
from concurrent.futures import ThreadPoolExecutor
import json
import math
import os
from pathlib import Path
import threading
import libephemeris as e

root = Path(os.environ["V4_TEST_DATA"])
present = os.environ["V4_TEST_PRESENT"] == "1"
e.set_ephe_path(str(root))
e.set_ephemeris_file("de440.bsp")
e.set_calc_mode("skyfield")
e.set_auto_spk_download(False)
e.set_strict_precision(False)
e.set_iers_delta_t_enabled(False)
e.set_delta_t_userdef(64.0 / 86400.0)
flags = e.FLG_SPEED | e.FLG_SWIEPH | e.FLG_HELCTR | e.FLG_TRUEPOS | e.FLG_J2000 | e.FLG_NONUT
jd = 2451545.0

token = e.start_tracing()
expected = e.calc(jd, e.CHIRON, flags)
source = e.get_trace_results()[e.CHIRON]
assert source == ("ASSIST" if present else "Keplerian"), source
assert expected[0][2] > 0.0 and all(math.isfinite(value) for value in expected[0])
assert token.reset()

# Release cached owners while workers race to acquire or use their own snapshots.
# Every calculation must finish with the same physical result either way.
barrier = threading.Barrier(9)
def calculate(_):
    barrier.wait()
    return e.calc(jd, e.CHIRON, flags)
with ThreadPoolExecutor(max_workers=8) as pool:
    futures = [pool.submit(calculate, index) for index in range(8)]
    barrier.wait()
    e.close()
    observed = [future.result() for future in futures]
assert observed == [expected] * 8
assert e.calc(jd, e.CHIRON, flags) == expected

# close() resets the selected planetary path; configure it again before a
# registered heliocentric kernel needs the Sun origin from the main kernel.
e.set_ephe_path(str(root))
e.set_ephemeris_file("de440.bsp")
# Actual registered coefficients take precedence over either numerical fallback.
e.register_spk_body(e.CHIRON, str(root / "spk/2060_184911_215002.bsp"), 20002060)
token = e.start_tracing()
registered = e.calc(jd, e.CHIRON, flags)
registered_source = e.get_trace_results()[e.CHIRON]
assert registered_source == "SPK", registered_source
assert all(math.isfinite(value) for value in registered[0])
assert token.reset()
e.unregister_spk_body(e.CHIRON)
assert e.calc(jd, e.CHIRON, flags) == expected

e.set_strict_precision(True)
if present:
    assert e.calc(jd, e.CHIRON, flags) == expected
else:
    try:
        e.calc(jd, e.CHIRON, flags)
    except e.SPKRequiredError:
        pass
    else:
        raise AssertionError("strict mode without a precision source must refuse")
print(json.dumps({"source": source, "registered_source": registered_source, "workers": len(observed)}))
'''


@pytest.mark.parametrize("present", [True, False], ids=["assist-present", "assist-absent"])
def test_installed_minor_precision_chain_in_isolated_environment(verified_files, tmp_path, present):
    home = tmp_path / "home"
    home.mkdir()
    absent = tmp_path / "absent-data"
    absent.mkdir()
    environment = {key: value for key, value in os.environ.items() if not key.startswith("LIBEPHEMERIS_")}
    environment.update(
        HOME=str(home), ASSIST_DIR=str(verified_files / "assist" if present else absent),
        LIBEPHEMERIS_NETWORK_POLICY="sealed", LIBEPHEMERIS_AUTO_SPK="0",
        LIBEPHEMERIS_IERS_AUTO_DOWNLOAD="0", LIBEPHEMERIS_IERS_DELTA_T="0",
        LIBEPHEMERIS_CONFIG=str(tmp_path / "absent.toml"),
        LIBEPHEMERIS_ENV_FILE=str(tmp_path / "absent.env"),
        LIBEPHEMERIS_DATA_DIR=str(tmp_path), V4_TEST_DATA=str(verified_files),
        V4_TEST_PRESENT="1" if present else "0",
    )
    process = subprocess.Popen(
        [sys.executable, "-I", "-c", SCRIPT], cwd=tmp_path, env=environment,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, start_new_session=True,
    )
    try:
        stdout, stderr = process.communicate(timeout=180)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()
        raise
    assert process.returncode == 0, stderr
    result = json.loads(stdout.strip())
    assert result == {"source": "ASSIST" if present else "Keplerian", "registered_source": "SPK", "workers": 8}
