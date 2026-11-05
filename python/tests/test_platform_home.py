"""Native home discovery must not require the Unix HOME variable on Windows."""

import os
import subprocess
import sys

import pytest


@pytest.mark.parametrize("explicit_home", [False, True], ids=["native-profile", "home-override"])
def test_import_and_default_data_root_use_the_user_home(tmp_path, explicit_home):
    home = tmp_path / "user profile with spaces é"
    home.mkdir()
    environment = {
        key: value for key, value in os.environ.items()
        if not key.startswith("LIBEPHEMERIS_") and key != "HOME"
    }
    if explicit_home:
        environment["HOME"] = str(home)
    elif os.name == "nt":
        # PowerShell/cmd normally provide USERPROFILE, not HOME. A temporary
        # profile also verifies that Unicode paths survive native conversion.
        environment["USERPROFILE"] = str(home)
    environment.update(
        LIBEPHEMERIS_CONFIG=str(tmp_path / "absent.toml"),
        LIBEPHEMERIS_ENV_FILE=str(tmp_path / "absent.env"),
        LIBEPHEMERIS_NETWORK_POLICY="sealed",
        LIBEPHEMERIS_IERS_DELTA_T="0",
        LIBEPHEMERIS_AUTO_SPK="0",
        PYTHONUTF8="1",
    )
    program = '''
from pathlib import Path
import sys
import libephemeris as eph
expected = Path(sys.argv[1]) if sys.argv[1] else Path.home()
assert Path(eph.get_data_dir()) == expected / ".libephemeris"
assert Path(eph.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
'''
    result = subprocess.run(
        [sys.executable, "-I", "-c", program, str(home) if explicit_home else ""],
        cwd=tmp_path, env=environment, capture_output=True, text=True,
        encoding="utf-8", timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
