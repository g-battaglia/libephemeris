"""Native configuration, resource metadata and deterministic lifecycle checks."""
from data_paths import CORE_LEB

from pathlib import Path
import inspect

import libephemeris as eph
import pytest

FIXTURE = CORE_LEB


@pytest.fixture(autouse=True)
def configured():
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_leb_file(str(FIXTURE))
    eph.set_delta_t_userdef(64.0 / 86400.0)
    eph.set_iers_delta_t_enabled(False)
    yield
    eph.close()
    eph.set_calc_mode("leb")
    eph.set_leb_file(str(FIXTURE))


def test_setters_and_getters():
    eph.set_tid_acc(-25.0)
    assert eph.get_tid_acc() == -25.0
    eph.set_tid_acc(999999.0)
    assert eph.get_tid_acc() == eph.TIDAL_DEFAULT
    eph.set_lapse_rate(0.006)
    assert eph.get_lapse_rate() == 0.006
    eph.set_lapse_rate(None)
    assert eph.get_lapse_rate() == 0.0065
    eph.set_auto_spk_download(False)
    assert not eph.get_auto_spk_download()
    eph.set_strict_precision(False)
    assert not eph.get_strict_precision()
    eph.set_spk_cache_dir("cache-test")
    assert eph.get_spk_cache_dir() == "cache-test"
    eph.set_spk_date_padding(42)
    assert eph.get_spk_date_padding() == 42
    with pytest.raises(ValueError):
        eph.set_spk_date_padding(-1)
    assert eph.get_spk_date_padding() == 42
    assert eph.get_delta_t_userdef() == 64.0 / 86400.0
    eph.set_delta_t_model("espenak_meeus")
    assert eph.get_delta_t_model() == "espenak_meeus"


def test_tiers_and_signatures():
    rows = eph.list_tiers()
    assert [row.name for row in rows] == ["base", "medium", "extended"]
    assert rows[0].ephemeris_file == "de440s.bsp"
    assert rows[2].spk_date_range == ("1600-01-01", "2500-01-01")
    assert eph.get_spk_date_range_for_tier("medium") == ("1850-01-01", "2150-01-01")
    assert eph.get_spk_date_range_for_tier("extended") == ("1600-01-01", "2500-01-01")
    with pytest.raises(AttributeError):
        rows[0].name = "changed"
    assert list(inspect.signature(eph.PrecisionTier).parameters) == ["name", "ephemeris_file", "spk_date_range", "description"]
    assert inspect.signature(eph.set_topo).parameters["alt"].default == 0.0


def test_reader_ownership_and_real_metadata():
    first = eph.get_leb_reader()
    assert first.path == str(FIXTURE)
    assert first.has_body(0)
    values, rates = first.eval_body(0, 2451545.0)
    assert len(values) == len(rates) == 3
    path, begin, end, number = eph.get_current_file_data(0)
    assert path == first.path and begin <= 2451545.0 <= end and number == 0
    expected = eph.calc_ut(2451545.0, 0, 258)
    first.close()
    assert first.closed
    assert eph.calc_ut(2451545.0, 0, 258) == expected
    assert eph.get_current_file_data(99) == ("", 0.0, 0.0, 0)
    eph.release_data_cache()


def test_horizons_selection_does_not_open_leb(tmp_path):
    eph.set_leb_file(str(tmp_path / "missing.leb2"))
    eph.set_calc_mode("horizons")
    assert eph.get_leb_reader() is None
    assert eph.get_calc_mode() == "horizons"


def test_empty_ephe_path_selects_current_directory_and_none_restores_discovery(tmp_path):
    from pathlib import Path

    eph.set_ephe_path(None)
    default = eph.get_library_path()
    eph.set_ephe_path(str(tmp_path))
    assert eph.get_library_path() == str(tmp_path)
    eph.set_ephe_path("")
    assert eph.get_library_path() == str(Path.cwd())
    eph.release_data_cache()
    reader = eph.get_leb_reader()
    expected = eph.calc_ut(2451545.0, 0, 258)
    expected_vector = reader.eval_body(0, 2451545.0)
    generation = eph.get_leb_reader_generation()
    for path in (".", "", "."):
        eph.set_ephe_path(path)
        assert eph.get_library_path() == str(Path.cwd())
        assert eph.get_leb_reader_generation() == generation
        assert not reader.closed
        assert reader.eval_body(0, 2451545.0) == expected_vector
        assert eph.calc_ut(2451545.0, 0, 258) == expected
        eph.release_data_cache()
    eph.set_ephe_path(None)
    assert eph.get_library_path() == default
    eph.release_data_cache()


def test_cooling_never_acquires_missing_readers(tmp_path):
    eph.close()
    eph.set_leb_file(str(tmp_path / "missing.leb2"))
    assert eph.cool_leb_readers() == 0


def test_reader_generation_tracks_invalidation_without_acquiring_data(tmp_path):
    """Metadata caches can detect changed discovery without opening any file."""
    generation = eph.get_leb_reader_generation()
    assert eph.get_leb_reader_generation() == generation
    eph.set_leb_file(str(tmp_path / "missing-generation-fixture.leb2"))
    changed = eph.get_leb_reader_generation()
    assert changed != generation
    assert eph.cool_leb_readers() == 0
    eph.set_leb_file(str(tmp_path / "missing-generation-fixture.leb2"))
    assert eph.get_leb_reader_generation() == changed


def test_cooling_keeps_live_reader_and_results_valid():
    first = eph.get_leb_reader()
    expected = first.eval_body(0, 2451545.0)
    assert eph.cool_leb_readers(clear_coefficients=False) >= 1
    assert first.eval_body(0, 2451545.0) == expected
    assert eph.cool_leb_readers() >= 1
    assert not first.closed
    assert first.eval_body(0, 2451545.0) == expected
    assert eph.cool_leb_readers() >= 1


def test_reload_retains_explicit_overrides(monkeypatch):
    eph.close()
    eph.set_calc_mode(None)
    monkeypatch.setenv("LIBEPHEMERIS_MODE", "horizons")
    eph.reload_configuration()
    assert eph.get_calc_mode() == "horizons"
    eph.set_calc_mode("leb")
    assert eph.get_calc_mode() == "leb"
    eph.set_calc_mode(None)
    assert eph.get_calc_mode() == "horizons"
    monkeypatch.setenv("LIBEPHEMERIS_MODE", "auto")
    eph.reload_configuration()
    assert eph.get_calc_mode() == "auto"


def test_tier_constant_contains_typed_native_records():
    assert list(eph.TIERS) == ["base", "medium", "extended"]
    for name, tier in eph.TIERS.items():
        assert isinstance(tier, eph.PrecisionTier)
        assert tier.name == name
    assert list(eph.TIERS.values()) == eph.list_tiers()
