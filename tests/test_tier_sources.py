"""Focused checks for original-byte readers and the lazy factory hook."""

from __future__ import annotations

import math
import random
from pathlib import Path
from types import SimpleNamespace

import pytest

import libephemeris as ephe
import libephemeris.state as state
from libephemeris import CoefficientSourceError, Error
from libephemeris.leb2_reader import LEB2Reader
from libephemeris.leb2_remote import RemoteLEB2Reader

_DATA = Path(__file__).parents[1] / "libephemeris" / "data" / "leb2" / "base_core.leb2"


class BytesSource:
    name = "base_core.leb2"
    locator = "memory://tests"

    def __init__(self):
        self.payload = _DATA.read_bytes()
        self.closed = False

    def __len__(self):
        return len(self.payload)

    def read(self, offset, size):
        if self.closed:
            raise RuntimeError("private transport details")
        return self.payload[offset : offset + size]

    def close(self):
        self.closed = True


def test_original_bytes_match_all_native_states_and_boundaries():
    source = BytesSource()
    rng = random.Random(17)
    with (
        LEB2Reader(str(_DATA)) as native,
        RemoteLEB2Reader(source, reviewed=True) as remote,
    ):
        for body, entry in native._bodies.items():
            edge = entry.jd_start + entry.interval_days
            dates = [
                entry.jd_start,
                entry.jd_end,
                edge,
                math.nextafter(edge, -math.inf),
                math.nextafter(edge, math.inf),
            ]
            dates += [rng.uniform(entry.jd_start, entry.jd_end) for _ in range(150)]
            for jd in dates:
                assert remote.eval_body(body, jd) == native.eval_body(body, jd)
        for jd in (native._nutation.jd_start, 2451545.0, native._nutation.jd_end):
            assert remote.eval_nutation(jd) == native.eval_nutation(jd)
            assert remote.delta_t(jd) == native.delta_t(jd)
        assert remote._stars == native._stars
    assert source.closed


@pytest.mark.parametrize("mode", ["leb", "auto"])
def test_factory_is_lazy_cached_and_failures_do_not_fall_back(monkeypatch, mode):
    ephe.set_leb_file(None)
    monkeypatch.delenv("LIBEPHEMERIS_LEB", raising=False)
    monkeypatch.setenv("LIBEPHEMERIS_LEB_SOURCE", "test_source:open_reader")
    monkeypatch.setattr(state, "get_calc_mode", lambda: mode)
    calls = []
    reader = RemoteLEB2Reader(BytesSource(), reviewed=True)
    monkeypatch.setattr(
        state.importlib,
        "import_module",
        lambda _: SimpleNamespace(open_reader=lambda: calls.append(1) or reader),
    )
    assert state.get_leb_reader() is reader
    assert state.get_leb_reader() is reader
    assert calls == [1]
    ephe.set_leb_file(None)
    monkeypatch.setattr(
        state.importlib,
        "import_module",
        lambda _: (_ for _ in ()).throw(RuntimeError("secret")),
    )
    with pytest.raises(CoefficientSourceError) as error:
        state.get_leb_reader()
    assert "secret" not in str(error.value)


@pytest.mark.parametrize("operation", ["body", "nutation"])
def test_transport_and_corruption_are_fatal(operation):
    source = BytesSource()
    with RemoteLEB2Reader(source, reviewed=True) as reader:
        source.closed = True
        with pytest.raises(CoefficientSourceError):
            if operation == "body":
                reader.eval_body(0, 2451545.0)
            else:
                reader.eval_nutation(2451545.0)
    source = BytesSource()
    with RemoteLEB2Reader(source, reviewed=True) as reader:
        chunk = reader._chunk_index[0][0]
        source.payload = (
            source.payload[: chunk.blob_offset]
            + b"x" * chunk.compressed_size
            + source.payload[chunk.blob_offset + chunk.compressed_size :]
        )
        with pytest.raises(CoefficientSourceError):
            reader.eval_body(0, chunk.jd_start)


def test_short_read_and_bad_header_are_source_failures():
    source = BytesSource()
    source.payload = b"bad"
    with pytest.raises(CoefficientSourceError):
        RemoteLEB2Reader(source, reviewed=True)
    assert source.closed


def test_public_source_error_is_plain_exception():
    assert CoefficientSourceError.__bases__ == (Exception,)
    assert not isinstance(CoefficientSourceError("x"), Error)
