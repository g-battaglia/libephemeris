"""IERS lookups and cache lifecycle use owned native row generations."""

from pathlib import Path
import shutil
import libephemeris as eph
import pytest
from data_paths import DATA_ROOT


@pytest.fixture
def cache(tmp_path):
    source = DATA_ROOT / "iers_cache"
    names = ("finals2000A.data", "leap_seconds.dat", "deltat.data")
    if not all((source / name).is_file() for name in names):
        pytest.skip("optional IERS source files unavailable")
    for name in names:
        shutil.copyfile(source / name, tmp_path / name)
    eph.set_iers_cache_dir(str(tmp_path))
    eph.set_iers_auto_download(False)
    eph.set_network_policy("sealed")
    try:
        yield tmp_path
    finally:
        eph.clear_iers_cache()
        eph.set_iers_cache_dir(None)
        eph.set_iers_auto_download(None)
        eph.set_network_policy(None)


def test_native_iers_rows_and_exact_time_identity(cache):
    assert eph.load_iers_data()
    assert eph.get_tai_utc(57754.0) == 37.0
    jd = 2457754.5
    observed = eph.get_observed_delta_t(jd)
    expected = observed if observed is not None else eph.get_tai_utc(jd - 2400000.5) + 32.184 - eph.get_ut1_utc(jd - 2400000.5)
    assert eph.get_delta_t_iers(jd) == expected
    assert eph.is_iers_data_available(jd)
    assert eph.get_iers_data_range()[0] < jd < eph.get_iers_data_range()[1]
    assert eph.get_observed_delta_t_data_range()[0] < jd
    assert eph.is_observed_delta_t_available(jd)


def test_fresh_files_are_kept_without_network_and_delete_is_scoped(cache):
    for function, name in ((eph.download_iers_finals,"finals2000A.data"), (eph.download_leap_seconds,"leap_seconds.dat"), (eph.download_delta_t_data,"deltat.data")):
        assert function() == str(cache / name)
    assert eph.get_iers_cache_dir() == str(cache)
    assert eph.get_iers_auto_download() is False
    assert eph.load_iers_data()
    info = eph.get_iers_cache_info()
    assert info["data_points"] > 0 and info["leap_second_entries"] > 0
    (cache / "unrelated.txt").write_text("keep")
    assert eph.delete_iers_cache_files() == 3
    assert (cache / "unrelated.txt").read_text() == "keep"
    assert eph.delete_iers_cache_files() == 0
