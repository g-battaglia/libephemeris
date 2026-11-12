"""Native discovery leaves unrelated files and absent directories alone."""

import libephemeris as eph
from data_paths import DATA_ROOT


def test_empty_and_non_spk_directories(tmp_path):
    assert eph.discover_local_spks(str(tmp_path / "missing")) == {}
    (tmp_path / "ordinary.txt").write_text("not a kernel")
    assert eph.discover_local_spks(str(tmp_path)) == {}
    assert (tmp_path / "ordinary.txt").read_text() == "not a kernel"


def test_auto_config_owns_fields_and_resolves_a_cache(tmp_path):
    import hashlib
    from libephemeris import spk_auto

    config = spk_auto.AutoSpkConfig(9999999, "Earth center", 399, cache_dir=str(tmp_path))
    digest = hashlib.md5(b"Earth center_2000-01-01_2100-01-01").hexdigest()[:8]
    assert config.get_cache_filename() == f"earth_center_{digest}.bsp"
    assert config.get_cache_dir() == str(tmp_path)
    assert config.get_cache_path() == str(tmp_path / config.get_cache_filename())
    assert config.enabled is True and config.spk_path is None
    assert (config.ipl, config.body_id, config.naif_id) == (9999999, "Earth center", 399)
    config.enabled = False
    config.start = "1990-01-01"
    assert not config.enabled and config.start == "1990-01-01"


def test_enabled_auto_config_has_shared_native_lifetime(tmp_path):
    from libephemeris import spk_auto

    spk_auto.disable_all()
    try:
        spk_auto.enable_auto_spk(9999999, "399", 399, cache_dir=str(tmp_path))
        config = spk_auto.get_auto_spk_config(9999999)
        assert config.body_id == "399"
        assert spk_auto.is_auto_spk_enabled(9999999)
        config.enabled = False
        assert not spk_auto.is_auto_spk_enabled(9999999)
        spk_auto.list_auto_spk_bodies()[9999999].enabled = True
        assert config.enabled
        spk_auto.disable_auto_spk(9999999)
        assert spk_auto.get_auto_spk_config(9999999) is None
        assert config.body_id == "399"  # retained record remains owned
    finally:
        spk_auto.disable_all()


def test_cache_management_preserves_unrelated_files(tmp_path):
    from libephemeris import spk_auto
    import pytest

    unrelated = tmp_path / "ordinary.txt"
    unrelated.write_text("not a kernel")
    assert spk_auto.ensure_cache_dir(str(tmp_path)) == str(tmp_path)
    assert spk_auto.get_cache_path("Earth center", str(tmp_path)) == str(tmp_path / "earth_center.bsp")
    assert spk_auto.cache_info(str(tmp_path)) == dict(cache_dir=str(tmp_path), num_files=0, total_size_mb=0.0, files=[])
    assert spk_auto.list_cached_spk(str(tmp_path)) == []
    assert spk_auto.get_cache_size(str(tmp_path)) == 0.0
    assert spk_auto.clear_spk_cache(str(tmp_path)) == 0
    assert spk_auto.prune_old_cache(30, str(tmp_path)) == 0
    with pytest.raises(ValueError):
        spk_auto.prune_old_cache(0, str(tmp_path))
    with pytest.raises(ValueError):
        spk_auto.auto_get_spk("399", 2451545.0, 2451544.0, str(tmp_path))
    assert unrelated.read_text() == "not a kernel"


def test_cached_download_registers_the_actual_image_without_network(tmp_path):
    import hashlib
    from pathlib import Path
    import pytest
    from libephemeris import spk_auto

    kernel = DATA_ROOT / "de440.bsp"
    if not kernel.is_file():
        pytest.skip("optional DE440 kernel unavailable")
    with kernel.open("rb") as stream:
        assert hashlib.file_digest(stream, "sha256").hexdigest() == "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5"
    spk_auto.disable_all()
    try:
        spk_auto.enable_auto_spk(9999999, "399", 399, start="2000-01-01", end="2001-01-01", cache_dir=str(tmp_path))
        config = spk_auto.get_auto_spk_config(9999999)
        cached = Path(config.get_cache_path())
        cached.symlink_to(kernel)
        assert spk_auto.download_now(9999999) == str(cached)
        assert eph.get_spk_body_info(9999999) == (str(cached), 399)
        assert config.spk_path == str(cached)
        assert spk_auto.try_auto_download(9999999) == str(cached)
        assert spk_auto.clear_cache(9999999) == 1
        assert config.spk_path is None and not cached.exists() and kernel.is_file()
    finally:
        eph.unregister_spk_body(9999999)
        spk_auto.disable_all()


def test_covering_kernel_reused_by_auto_get_spk(tmp_path):
    from pathlib import Path
    import pytest
    from libephemeris import spk_auto

    kernel = DATA_ROOT / "de440.bsp"
    if not kernel.is_file():
        pytest.skip("optional DE440 kernel unavailable")
    cached = tmp_path / "399_200001_200101.bsp"
    cached.symlink_to(kernel)
    old_padding = eph.get_spk_date_padding()
    try:
        eph.set_spk_date_padding(0)
        assert spk_auto.is_spk_cached("399", 2451545.0, 2451910.0, str(tmp_path))
        assert spk_auto.auto_get_spk("399", 2451545.0, 2451910.0, str(tmp_path), 9999999, 399) == str(cached)
        assert eph.get_spk_body_info(9999999) == (str(cached), 399)
    finally:
        eph.unregister_spk_body(9999999)
        eph.set_spk_date_padding(old_padding)


def test_common_settings_are_eight_native_records(tmp_path):
    from libephemeris import spk_auto

    spk_auto.disable_all()
    try:
        spk_auto.enable_common_bodies(cache_dir=str(tmp_path))
        configs = spk_auto.list_auto_spk_bodies()
        assert set(configs) == {eph.CHIRON, eph.PHOLUS, eph.CERES, eph.PALLAS, eph.JUNO, eph.VESTA, eph.ERIS, eph.SEDNA}
        assert all(config.enabled and config.naif_id is not None for config in configs.values())
    finally:
        spk_auto.disable_all()
