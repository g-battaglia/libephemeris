"""Focused checks for mmap-like byte sources and the lazy reader factory."""

from __future__ import annotations

import math
import random
from pathlib import Path
from types import SimpleNamespace

import pytest

import libephemeris as ephe
import libephemeris.state as state
from libephemeris.leb2_reader import LEB2Reader

_DATA = Path(__file__).parents[1] / "libephemeris" / "data" / "leb2" / "base_core.leb2"


class BytesSource:
    def __init__(self, data):
        self.data = data
        self.closed = False
        self.reads = []

    def __len__(self):
        return len(self.data)

    def __getitem__(self, key):
        if self.closed:
            raise SourceFailure("source unavailable")
        self.reads.append(key)
        return self.data[key]

    def close(self):
        self.closed = True


class SourceFailure(Exception):
    pass


def test_original_bytes_match_native_states_and_boundaries():
    source = BytesSource(_DATA.read_bytes())
    rng = random.Random(17)
    with LEB2Reader(str(_DATA)) as native, LEB2Reader("memory", source) as remote:
        for body, entry in native._bodies.items():
            edge = entry.jd_start + entry.interval_days
            dates = [entry.jd_start, entry.jd_end, edge]
            dates += [math.nextafter(edge, -math.inf), math.nextafter(edge, math.inf)]
            dates += [rng.uniform(entry.jd_start, entry.jd_end) for _ in range(150)]
            for jd in dates:
                assert remote.eval_body(body, jd) == native.eval_body(body, jd)
        for jd in (native._nutation.jd_start, 2451545.0, native._nutation.jd_end):
            assert remote.eval_nutation(jd) == native.eval_nutation(jd)
            assert remote.delta_t(jd) == native.delta_t(jd)
        assert remote._stars == native._stars
        remote.warm(2451545.0, 2451546.0)
        remote.cool()
    assert source.closed


def test_metadata_reads_are_batched():
    source = BytesSource(_DATA.read_bytes())
    with LEB2Reader("memory", source) as reader:
        # Header, directory, bodies, nutation, Delta-T header/table, stars,
        # and two reads per body's chunk index, regardless of entry count.
        assert len(source.reads) <= 7 + 2 * len(reader._bodies)


@pytest.mark.parametrize("mode", ["skyfield", "horizons"])
def test_forced_backend_does_not_load_factory(monkeypatch, mode):
    ephe.set_leb_file(None)
    monkeypatch.setenv("LIBEPHEMERIS_LEB_SOURCE", "test_source:open_reader")
    monkeypatch.setattr(state, "get_calc_mode", lambda: mode)

    def unexpected(_):
        pytest.fail("Factory imported for a forced non-LEB backend")

    monkeypatch.setattr(state.importlib, "import_module", unexpected)
    assert state.get_leb_reader() is None


def test_explicit_file_takes_precedence_over_factory(monkeypatch):
    monkeypatch.setenv("LIBEPHEMERIS_LEB_SOURCE", "test_source:open_reader")
    monkeypatch.setattr(state, "get_calc_mode", lambda: "leb")

    def unexpected(_):
        pytest.fail("Factory imported despite an explicit local file")

    monkeypatch.setattr(state.importlib, "import_module", unexpected)
    ephe.set_leb_file(str(_DATA))
    try:
        assert state.get_leb_reader().eval_body(0, 2451545.0)
    finally:
        ephe.set_leb_file(None)


@pytest.mark.parametrize("mode", ["leb", "auto"])
def test_factory_is_lazy_cached_and_failure_propagates(monkeypatch, mode):
    ephe.set_leb_file(None)
    monkeypatch.delenv("LIBEPHEMERIS_LEB", raising=False)
    monkeypatch.setenv("LIBEPHEMERIS_LEB_SOURCE", "test_source:open_reader")
    monkeypatch.setattr(state, "get_calc_mode", lambda: mode)
    calls = []
    with LEB2Reader("memory", BytesSource(_DATA.read_bytes())) as reader:
        monkeypatch.setattr(
            state.importlib,
            "import_module",
            lambda _: SimpleNamespace(open_reader=lambda: calls.append(1) or reader),
        )
        assert not calls
        assert state.get_leb_reader() is reader
        assert state.get_leb_reader() is reader
        assert calls == [1]
        ephe.set_leb_file(None)
    failure = SourceFailure("source unavailable")

    def fail(_):
        raise failure

    monkeypatch.setattr(state.importlib, "import_module", fail)
    with pytest.raises(SourceFailure) as error:
        state.get_leb_reader()
    assert error.value is failure
    with pytest.raises(SourceFailure):
        ephe.get_leb_inventory()


def test_tier_change_replaces_reader_and_vector_adapter(monkeypatch):
    from libephemeris.leb_vector import get_leb_vector_ephemeris

    payload = _DATA.read_bytes()
    monkeypatch.setenv("LIBEPHEMERIS_LEB_SOURCE", "test_source:open_reader")
    monkeypatch.setattr(
        state.importlib,
        "import_module",
        lambda _: SimpleNamespace(
            open_reader=lambda: LEB2Reader("memory", BytesSource(payload))
        ),
    )
    ephe.set_leb_file(None)
    original_tier = ephe.get_precision_tier()
    try:
        ephe.set_precision_tier("base")
        first = get_leb_vector_ephemeris(state.get_leb_reader())
        ephe.set_precision_tier("medium")
        second = get_leb_vector_ephemeris(state.get_leb_reader())
        assert first is not second
        assert first.reader is not second.reader
        assert second.reader is state.get_leb_reader()
        assert second.reader.eval_body(0, 2451545.0)
    finally:
        ephe.set_precision_tier(original_tier)
        ephe.set_leb_file(None)


@pytest.mark.parametrize("operation", ["body", "nutation"])
def test_source_failure_propagates_during_evaluation(operation):
    source = BytesSource(_DATA.read_bytes())
    with LEB2Reader("memory", source) as reader:
        source.close()
        with pytest.raises(SourceFailure):
            if operation == "body":
                reader.eval_body(0, 2451545.0)
            else:
                reader.eval_nutation(2451545.0)
