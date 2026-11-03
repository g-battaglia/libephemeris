"""The HTML book preserves native help and resolves canonical source links."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace
import os
import re
import subprocess
import sys

import pytest

pytest.importorskip("libephemeris")
ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("documentation_config", ROOT / "docs/conf.py")
CONFIG = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CONFIG)


def test_native_help_is_preserved_verbatim_not_reinterpreted_as_rst():
    original = ["Units: |x| in AU.", "", "```python", "print('hello')", "```"]
    lines = original.copy()
    CONFIG.native_help_text(None, "function", "calc", None, None, lines)
    assert lines[:2] == [".. code-block:: text", ""]
    assert [line.removeprefix("    ") for line in lines[2:-1]] == original


def test_source_links_keep_book_links_but_resolve_project_files():
    app = SimpleNamespace(
        srcdir=str(ROOT),
        env=SimpleNamespace(doc2path=lambda name: str(ROOT / (name + ".md"))),
    )
    source = ["[book](python.md) [manual](../beginners-manual/en/README.md) "
              "[root](../README.md#first-calculation)"]
    CONFIG.source_links(app, "docs/README", source)
    assert source == [
        "[book](python.md) [manual](../beginners-manual/en/README.md) "
        f"[root](https://github.com/g-battaglia/libephemeris/blob/v{CONFIG.release}/README.md#first-calculation)"
    ]


def test_missing_source_destinations_fail_instead_of_suppressing_warnings():
    app = SimpleNamespace(
        srcdir=str(ROOT),
        env=SimpleNamespace(doc2path=lambda name: str(ROOT / (name + ".md"))),
    )
    with pytest.raises(ValueError, match="invalid source link"):
        CONFIG.source_links(app, "docs/README", ["[missing](../missing-root-file.md)"])


def test_migration_boundaries_do_not_depend_on_a_reference_checkout():
    """Keep these user-visible exemptions checked even without the v3 fixture."""
    migration = (ROOT / 'MIGRATION.md').read_text()
    assert 'exception graphs' in migration
    assert 'native named values' in migration.lower()


def test_beginner_examples_run_in_order_without_external_data_or_network(tmp_path):
    """The tutorial is executable instruction, not an untested output transcript."""
    tutorial = (ROOT / 'docs/installation.md').read_text()
    blocks = re.findall(r'```python\n(.*?)```', tutorial, flags=re.DOTALL)
    assert len(blocks) == 2
    program = '\n'.join(blocks) + '\nassert ephe.julday(2024, 4, 8, 12.0) == 2460409.0\n'
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith('LIBEPHEMERIS_')}
    environment.update(LIBEPHEMERIS_DATA_DIR=str(tmp_path),
                       LIBEPHEMERIS_CONFIG=str(tmp_path / 'absent.toml'),
                       LIBEPHEMERIS_ENV_FILE=str(tmp_path / 'absent.env'),
                       LIBEPHEMERIS_NETWORK_POLICY='sealed')
    result = subprocess.run([sys.executable, '-c', program], cwd=tmp_path,
                            env=environment, text=True, capture_output=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'Moon longitude:' in result.stdout and 'Ascendant:' in result.stdout


def test_readme_chart_example_runs_offline_and_isolates_sidereal_state(tmp_path):
    """Exercise the actual landing-page code, including its chart conventions."""
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    blocks = re.findall(r"```python\n(.*?)```", readme, flags=re.DOTALL)
    assert len(blocks) == 1
    # Check units/layouts and this date's motion, not merely successful printing.
    checks = '''
assert set(positions) == {"Sun", "Moon", "Mercury"}
assert all(len(values) == 6 for values in positions.values())
assert all(0 <= values[0] < 360 and values[2] > 0 for values in positions.values())
assert positions["Sun"][3] > 0
assert positions["Mercury"][3] < 0
assert len(cusps) == 12 and all(0 <= cusp < 360 for cusp in cusps)
assert all(0 <= angle < 360 for angle in angles[:2])
assert 20 < (positions["Moon"][0] - moon_sidereal[0]) % 360 < 30
# A context's Lahiri selection must not change the module-level tropical result.
moon_again, _ = ephe.calc_ut(jd_ut1, ephe.MOON, ephe.FLG_SPEED)
assert moon_again == positions["Moon"]
'''
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith("LIBEPHEMERIS_")}
    environment.update(LIBEPHEMERIS_DATA_DIR=str(tmp_path),
                       LIBEPHEMERIS_CONFIG=str(tmp_path / "absent.toml"),
                       LIBEPHEMERIS_ENV_FILE=str(tmp_path / "absent.env"),
                       LIBEPHEMERIS_NETWORK_POLICY="sealed")
    result = subprocess.run([sys.executable, "-c", blocks[0] + checks], cwd=tmp_path,
                            env=environment, text=True, capture_output=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "retrograde" in result.stdout
    assert "Ascendant:" in result.stdout and "Lahiri sidereal longitude:" in result.stdout
