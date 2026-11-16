"""Exercise the runtime-owned coefficient cache in isolated Python processes."""

import json
import os
from pathlib import Path
import subprocess
import sys

from data_paths import PUBLISHED_CORE_LEB


def test_process_cache_precedence_resize_and_retained_limit(tmp_path):
    config = tmp_path / "config.toml"
    config.write_text("[libephemeris]\nleb_cache_bytes = 16777216\n", encoding="utf-8")
    env = {name: value for name, value in os.environ.items() if not name.startswith("LIBEPHEMERIS_")}
    env.update(
        LIBEPHEMERIS_CONFIG=str(config),
        LIBEPHEMERIS_ENV_FILE=str(tmp_path / "absent.env"),
        LIBEPHEMERIS_DATA_DIR=str(tmp_path),
    )
    script = """
import json, struct, sys
import libephemeris as eph
from libephemeris import state

assert state.set_leb_cache_bytes is eph.set_leb_cache_bytes
eph.set_calc_mode("leb")
eph.set_network_policy("sealed")
eph.set_leb_file(sys.argv[1])
epochs = [2451545.0, 2451545.5, 2460000.0]
def calculate():
    return [(struct.pack(">6d", *values).hex(), flags)
            for values, flags in (eph.calc_ut(jd, 4, eph.FLG_SPEED) for jd in epochs)]
original = calculate()
retained, limit, entries = eph.get_leb_cache_usage()
assert 0 < retained <= limit and entries > 0
eph.set_leb_cache_bytes(0)
assert eph.get_leb_cache_usage() == (0, 0, 0)
assert calculate() == original
assert eph.get_leb_cache_usage() == (0, 0, 0)
eph.set_leb_cache_bytes(64 * 1024 * 1024)
assert calculate() == original
retained, new_limit, entries = eph.get_leb_cache_usage()
assert 0 < retained <= new_limit and new_limit == 64 * 1024 * 1024 and entries > 0
print(json.dumps((limit, original)))
"""
    cmd = [sys.executable, "-I", "-c", script, str(Path(PUBLISHED_CORE_LEB).resolve())]
    toml_limit, toml_values = json.loads(subprocess.check_output(cmd, cwd=tmp_path, env=env))
    env["LIBEPHEMERIS_LEB_CACHE_MB"] = "32"
    env_limit, env_values = json.loads(subprocess.check_output(cmd, cwd=tmp_path, env=env))
    assert (toml_limit, env_limit) == (16 * 1024**2, 32 * 1024**2)
    assert toml_values == env_values
