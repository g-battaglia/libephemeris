"""SPK metadata, native registration and retained global state."""

from pathlib import Path
import libephemeris as eph
import pytest


def test_absent_spk_has_no_coverage(tmp_path):
    assert eph.get_spk_coverage(tmp_path / "missing.bsp") is None
    eph.unregister_spk_body(9999999)
    assert eph.get_spk_body_info(9999999) is None


def test_native_registration_accepts_documented_numeric_target_strings():
    from data_paths import DATA_ROOT
    kernel = DATA_ROOT / "de440.bsp"
    if not kernel.is_file():
        pytest.skip("optional DE440 kernel unavailable")
    try:
        eph.register_spk_body(9999999, str(kernel), "399")
        assert eph.get_spk_body_info(9999999) == (str(kernel), 399)
        assert eph.list_spk_bodies()[9999999] == (str(kernel), 399)
        coverage = eph.get_spk_coverage(kernel)
        assert isinstance(coverage, tuple) and coverage[0] < 2451545.0 < coverage[1]
    finally:
        eph.unregister_spk_body(9999999)
