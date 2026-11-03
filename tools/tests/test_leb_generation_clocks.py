"""Explicit offline clock configuration, native routing and publication guards."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path
import struct
from types import SimpleNamespace

import libephemeris as engine
from libephemeris import iers_data
from libephemeris.time_utils import deltat
import pytest

from tools.leb_generation import format as wire
from tools.leb_generation import nutation
from tools.leb_generation.generate import auxiliary_sections, generate_leb1, generate_leb2
from tools.leb_generation.options import create_source, source_arguments, validate_source_options
from tools.leb_generation.source import IERS_FILES, NativeSource, check_input_hashes
from tools.tests.test_leb_source_generation import LinearSource


def _metadata_loader(path):
    # Clock tests need no astronomical image. Native clock calls remain real;
    # only unrelated main-image metadata acquisition is replaced.
    return SimpleNamespace(path=str(path), segment_descriptors_tdb=lambda: [])


def _clock_in_spawn(selection):
    engine.get_loader = lambda: _metadata_loader
    source = selection.open()
    result = (deltat(2451544.5), source.input_hashes(), engine.get_delta_t_model(),
              engine.get_iers_delta_t_enabled(), engine.get_iers_auto_download(),
              engine.get_network_policy())
    selection.check_inputs()
    return result


@pytest.fixture
def local_source(tmp_path, monkeypatch):
    kernel = tmp_path / "synthetic.bsp"
    kernel.write_bytes(b"synthetic metadata only")
    monkeypatch.setattr(engine, "get_loader", lambda: _metadata_loader)
    yield kernel
    engine.close()
    engine.set_delta_t_userdef(None)
    engine.set_delta_t_model("smh2016")
    engine.set_iers_delta_t_enabled(False)
    iers_data.set_iers_cache_dir(None)
    engine.set_iers_auto_download(False)
    engine.set_network_policy("sealed")


def _clock_directory(path):
    path.mkdir()
    (path / "deltat.data").write_text("2000 1 1 60.0\n2000 2 1 61.0\n")
    return path


def test_defaults_reset_prior_clock_configuration_without_selecting_local_files(local_source, tmp_path):
    default = NativeSource(local_source, "base")
    expected = deltat(2451544.5)
    cache = _clock_directory(tmp_path / "iers")
    NativeSource(local_source, "base", "espenak_meeus", delta_t_userdef_days=0.125,
                 iers_cache_dir=cache)
    assert deltat(2451544.5) == 0.125
    reopened = NativeSource(local_source, "base")
    assert deltat(2451544.5) == expected
    assert engine.get_delta_t_model() == "smh2016"
    assert engine.get_delta_t_userdef() is None
    assert engine.get_iers_delta_t_enabled() is False
    assert engine.get_iers_auto_download() is False
    assert engine.get_network_policy() == "sealed"
    assert reopened.input_hashes() == default.input_hashes() == {local_source: default.sha256}


@pytest.mark.parametrize("days", [0.001, -0.0])
def test_override_has_native_priority_and_auxiliary_stores_exact_days(local_source, tmp_path, monkeypatch, days):
    cache = _clock_directory(tmp_path / "iers")
    NativeSource(local_source, "base", "espenak_meeus", delta_t_userdef_days=days,
                 iers_cache_dir=cache)
    assert iers_data.get_delta_t_iers(2451544.5) == 60.0
    monkeypatch.setattr(nutation, "fit_nutation", lambda *args: b"")
    data = auxiliary_sections(2451544.5, 2451546.5, tmp_path)[wire.SECTION_DELTA_T]
    assert struct.unpack_from(wire.DELTA_T_HEADER_FMT, data) == (2, 0)
    for index, epoch in enumerate((2451544.5, 2451546.5)):
        actual_epoch, value = struct.unpack_from(wire.DELTA_T_ENTRY_FMT, data,
                                               wire.DELTA_T_HEADER_SIZE + index * wire.DELTA_T_ENTRY_SIZE)
        assert actual_epoch == epoch
        assert struct.pack("d", value) == struct.pack("d", days)


def test_local_iers_routes_observed_rows_and_keeps_native_model_outside_coverage(local_source, tmp_path):
    NativeSource(local_source, "base", "espenak_meeus")
    outside = deltat(2300000.5)
    cache = _clock_directory(tmp_path / "iers")
    source = NativeSource(local_source, "base", "espenak_meeus", iers_cache_dir=cache)
    assert deltat(2451544.5) == 60.0 / 86400.0
    assert iers_data.get_delta_t_iers(2300000.5) is None
    assert deltat(2300000.5) == outside
    assert set(source.clock_hashes) == {cache / name for name in IERS_FILES}
    assert source.clock_hashes[cache / "finals2000A.data"] is None
    assert source.clock_hashes[cache / "leap_seconds.dat"] is None
    check_input_hashes(source.input_hashes().items())


def test_opening_new_clock_directory_cannot_retain_missing_previous_finals(local_source, tmp_path):
    first = tmp_path / "first"
    first.mkdir()
    (first / "finals2000A.data").write_text(f"000101 {51544.0:8.1f}{' ' * 42}I{0.1:10.7f}\n")
    NativeSource(local_source, "base", iers_cache_dir=first)
    assert iers_data.get_delta_t_iers(2451544.5) is not None
    second = tmp_path / "second"
    second.mkdir()
    NativeSource(local_source, "base", iers_cache_dir=second)
    assert iers_data.get_delta_t_iers(2451544.5) is None
    actual = deltat(2451544.5)
    NativeSource(local_source, "base")
    assert deltat(2451544.5) == actual


@pytest.mark.parametrize("days", [None, 0.002])
def test_spawned_selection_reopens_the_same_native_clock_and_all_identities(local_source, tmp_path, days):
    cache = _clock_directory(tmp_path / "iers")
    source = NativeSource(local_source, "medium", "espenak_meeus", delta_t_userdef_days=days,
                          iers_cache_dir=cache)
    expected = deltat(2451544.5)
    selection = source.isolated_selection(0)
    with ProcessPoolExecutor(max_workers=1, mp_context=get_context("spawn")) as pool:
        value, hashes, model, enabled, auto, policy = pool.submit(_clock_in_spawn, selection).result()
    assert value == expected
    assert hashes == source.input_hashes() == dict(selection.hashes)
    assert model == "espenak_meeus" and enabled is True and auto is False and policy == "sealed"
    (cache / "deltat.data").write_text("2000 1 1 62.0\n")
    with pytest.raises(ValueError, match="source identity changed"):
        selection.open()


@pytest.mark.parametrize("generator", [generate_leb1, generate_leb2])
@pytest.mark.parametrize("change", ["replace", "remove", "create"])
def test_changed_clock_inputs_refuse_publication_and_preserve_output(local_source, tmp_path, generator, change):
    cache = _clock_directory(tmp_path / "iers")
    source = NativeSource(local_source, "base", iers_cache_dir=cache)
    synthetic = LinearSource()
    source.require_domain = synthetic.require_domain
    source.states = synthetic.states
    source.values = synthetic.values
    destination = tmp_path / "preserved.leb"
    destination.write_bytes(b"existing destination")

    def mutate(_body):
        if change == "remove":
            (cache / "deltat.data").unlink()
        elif change == "create":
            (cache / "leap_seconds.dat").write_text("41317.0 1 1 1972 10\n")
        else:
            (cache / "deltat.data").write_text("2000 1 1 62.0\n")

    with pytest.raises(ValueError, match="generation source changed"):
        generator(source, destination, [0], 2451545.0, 2451549.0, skip_aux=True,
                  overwrite=True, progress=mutate)
    assert destination.read_bytes() == b"existing destination"
    assert not list(tmp_path.glob(".preserved.leb-*"))


def test_clock_cli_selection_is_explicit_and_reaches_native_source(local_source, tmp_path):
    parser = argparse.ArgumentParser()
    source_arguments(parser)
    defaults = parser.parse_args([])
    assert defaults.delta_t_model == "smh2016" and defaults.delta_t_userdef_days is None
    assert defaults.iers_delta_t is False and defaults.iers_cache_dir is None
    for arguments in (["--iers-delta-t"], ["--iers-cache-dir", str(tmp_path)]):
        with pytest.raises(ValueError, match="must be supplied together"):
            validate_source_options(parser.parse_args(arguments))
    cache = _clock_directory(tmp_path / "iers")
    args = parser.parse_args(["--kernel", str(local_source), "--iers-delta-t", "--iers-cache-dir", str(cache),
                              "--delta-t-userdef-days", "0.003", "--delta-t-model", "espenak_meeus"])
    selected = create_source(args, engine)
    assert deltat(2451544.5) == 0.003
    assert selected.isolated_selection(0).delta_t_userdef_days == 0.003


@pytest.mark.parametrize("days", [float("nan"), float("inf"), -float("inf"), engine.DELTAT_AUTOMATIC])
def test_invalid_explicit_override_is_not_silently_an_automatic_model(local_source, days):
    with pytest.raises(ValueError, match="finite explicit value"):
        NativeSource(local_source, "base", delta_t_userdef_days=days)


@pytest.mark.parametrize("clock_options", [["--delta-t-userdef-days", "0"],
                                          ["--iers-delta-t", "--iers-cache-dir", "."]])
def test_merge_and_structural_verification_cannot_ignore_clock_options(tmp_path, capsys, clock_options):
    from tools.leb_generation import cli, cli_leb1, cli_leb2

    for main, arguments in (
        (cli.main, ["base", "merge", "--output-dir", str(tmp_path)]),
        (cli_leb1.main, ["--merge", "unused.leb", "-o", str(tmp_path / "merged.leb")]),
    ):
        with pytest.raises(SystemExit) as failure:
            main([*arguments, *clock_options])
        assert failure.value.code == 1
        assert "does not accept" in capsys.readouterr().err
    path = tmp_path / "input.leb2"
    generate_leb2(LinearSource(), path, [0], 2451545.0, 2451549.0, skip_aux=True)
    with pytest.raises(SystemExit) as failure:
        cli_leb2.main(["verify", str(path), "--structural-only", *clock_options])
    assert failure.value.code == 1
    assert "does not accept astronomical source selections" in capsys.readouterr().err
