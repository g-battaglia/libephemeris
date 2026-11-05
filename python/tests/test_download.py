"""Native installers publish the bundled coefficient image into owned storage."""
from data_paths import CORE_LEB

from pathlib import Path
import hashlib
import libephemeris as eph
import pytest

SHA256 = "b891de02f82b68346046da67389fd3713aa81ab299287c632cb20cc396397964"


def test_bundled_base_install_and_explicit_activation(tmp_path):
    settings = tmp_path / "settings.env"
    settings.write_text(f"LIBEPHEMERIS_DATA_DIR={tmp_path}\n")
    eph.load_dotenv(settings, override=True)
    eph.set_network_policy("sealed")
    eph.set_calc_mode("leb")
    eph.set_leb_file(str(CORE_LEB))
    old = eph.get_leb_reader().path
    try:
        installed = eph.download_leb_for_tier("base", quiet=True, activate=False)
        assert isinstance(installed, Path)
        assert installed == tmp_path / "leb/base_core.leb2"
        with installed.open("rb") as stream:
            assert hashlib.file_digest(stream, "sha256").hexdigest() == SHA256
        assert eph.get_leb_reader().path == old
        assert eph.download_leb_for_tier("base", quiet=True) == installed
        assert eph.get_leb_reader().path == str(installed)
        for function in (eph.download_for_tier, eph.download_leb_for_tier):
            with pytest.raises(ValueError):
                function("unknown", quiet=True)
    finally:
        eph.set_leb_file(old)
        eph.set_network_policy(None)
        eph.reload_configuration()
