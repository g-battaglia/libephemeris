"""Discovery generations and immutable verified images through the installed API."""
from data_paths import CORE_LEB, PUBLISHED_CORE_LEB

from pathlib import Path
import os
import subprocess
import sys

import pytest

ROOT = Path(__file__).parents[2]
FIXTURE = CORE_LEB


def run_isolated(code, data_dir, *args):
    env = {key: value for key, value in os.environ.items() if not key.startswith("LIBEPHEMERIS_")}
    env.update(
        LIBEPHEMERIS_DATA_DIR=str(data_dir),
        LIBEPHEMERIS_CONFIG=str(data_dir / "missing.toml"),
        LIBEPHEMERIS_ENV_FILE=str(data_dir / "missing.env"),
        LIBEPHEMERIS_MODE="leb",
        LIBEPHEMERIS_PRECISION="base",
        LIBEPHEMERIS_NETWORK_POLICY="sealed",
    )
    result = subprocess.run(
        [sys.executable, "-c", code, *(str(arg) for arg in args)],
        env=env, text=True, capture_output=True, timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize('fixture', [CORE_LEB, PUBLISHED_CORE_LEB], ids=['current', 'published'])
def test_warm_discovery_uses_retained_images_and_reload_opens_new_owner(tmp_path, fixture):
    run_isolated(
        """
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import os
import shutil
import sys
import libephemeris as e
from libephemeris.download import DATA_FILES

root, fixture = map(Path, sys.argv[1:])
reviewed = hashlib.sha256(fixture.read_bytes()).hexdigest() == DATA_FILES['base_core.leb2']['sha256']
directory = root / 'leb'
directory.mkdir()
path = directory / 'base_core.leb2'
shutil.copyfile(fixture, path)
e.set_leb_file(None)
old_reader = e.get_leb_reader()
old_state = old_reader.eval_body(e.SUN, 2451545.0)
assert e.get_body_coverage(e.SUN, 2451545.0).reviewed is reviewed
expected = e.calc(2451545.0, e.SUN, e.FLG_SPEED)

hidden = root / 'retired'
directory.rename(hidden)
def warm(_):
    assert e.get_leb_reader().path == str(path)
    assert e.get_body_coverage(e.SUN, 2451545.0).reviewed is reviewed
    return e.calc(2451545.0, e.SUN, e.FLG_SPEED)
with ThreadPoolExecutor(max_workers=8) as workers:
    assert list(workers.map(warm, range(24))) == [expected] * 24
hidden.rename(directory)

replacement = root / 'replacement.leb2'
replacement.write_bytes(fixture.read_bytes() + b'\\0')
os.replace(replacement, path)
assert e.get_body_coverage(e.SUN, 2451545.0).reviewed is reviewed
e.reload_configuration()
assert not e.get_body_coverage(e.SUN, 2451545.0).reviewed
assert old_reader.eval_body(e.SUN, 2451545.0) == old_state
assert e.calc(2451545.0, e.SUN, e.FLG_SPEED) == expected
""",
        tmp_path, tmp_path, fixture,
    )


def test_explicit_missing_and_corrupt_files_never_use_discovered_fallback(tmp_path):
    run_isolated(
        """
from pathlib import Path
import shutil
import sys
import libephemeris as e

root, fixture = map(Path, sys.argv[1:])
(root / 'leb').mkdir()
shutil.copyfile(fixture, root / 'leb/base_core.leb2')
missing = root / 'selected.leb2'
e.set_leb_file(str(missing))
# Tropical geometry consumes no planetary data, even with an explicit path.
expected_houses = e.houses(2451545.0, 41.9, 12.5, ord('O'))
expected_delta = e.deltat(2451545.0)
expected_nutation = e.calc_ut(2451545.0, e.ECL_NUT)
for content in (None, b'corrupt'):
    if content is not None:
        missing.write_bytes(content)
    try:
        e.calc_ut(2451545.0, e.SUN, e.FLG_SWIEPH | e.FLG_SPEED)
    except e.DataNotFoundError:
        pass
    else:
        raise AssertionError('explicit file error was hidden')
shutil.copyfile(fixture, missing)
assert e.get_leb_reader().path == str(missing)
assert e.houses(2451545.0, 41.9, 12.5, ord('O')) == expected_houses
assert e.deltat(2451545.0) == expected_delta
assert e.calc_ut(2451545.0, e.ECL_NUT) == expected_nutation
""",
        tmp_path, tmp_path, FIXTURE,
    )


def test_medium_calculation_excludes_coexisting_extended_images(tmp_path):
    from data_paths import DATA_ROOT
    data = DATA_ROOT
    if not all((data / f"leb/{tier}_core.leb2").is_file() for tier in ("base", "medium", "extended")):
        pytest.skip("read-only mixed-tier data is not available")
    run_isolated(
        """
import math
import libephemeris as e

e.set_precision_tier('medium')
e.set_leb_file(None)
position, flags = e.calc_ut(2451545.0, e.SUN, e.FLG_SPEED)
assert all(math.isfinite(value) for value in position)
assert flags & e.FLG_SPEED
files = e.get_leb_inventory()['files']
assert files and any('medium_core.leb2' in row['path'] for row in files)
assert all('extended_' not in row['path'] for row in files)
""",
        data,
    )


def test_auto_unavailable_explicit_leb_is_retried_after_repair(tmp_path):
    run_isolated(
        """
from pathlib import Path
import shutil
import sys
import libephemeris as e

root, fixture = map(Path, sys.argv[1:])
(root / 'leb').mkdir()
shutil.copyfile(fixture, root / 'leb/base_core.leb2')
selected = root / 'selected.leb2'
e.set_leb_file(str(selected))
e.set_calc_mode('auto')
assert e.get_leb_reader() is None
selected.write_bytes(b'corrupt')
assert e.get_leb_reader() is None
shutil.copyfile(fixture, selected)
assert e.get_leb_reader().path == str(selected)
""",
        tmp_path, tmp_path, FIXTURE,
    )
