"""Configuration file merging is native and never mutates process variables."""

import os
import libephemeris as eph


def test_dotenv_file_merges_only_library_configuration(tmp_path):
    path = tmp_path / "settings.env"
    path.write_text('LIBEPHEMERIS_MODE="skyfield"\nUNRELATED_NATIVE_TEST=value\n')
    original = dict(os.environ)
    eph.set_calc_mode(None)
    try:
        assert eph.load_dotenv(path, override=True)
        assert eph.get_calc_mode() == "skyfield"
        assert dict(os.environ) == original
        assert not eph.load_dotenv(tmp_path / "missing.env")
    finally:
        eph.reload_configuration()
        eph.set_calc_mode("leb")
