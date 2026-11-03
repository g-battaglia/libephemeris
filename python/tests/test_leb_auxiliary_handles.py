"""Public auxiliary handles evaluate synthetic stored data, not global models."""
from data_paths import CORE_LEB

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import struct

import libephemeris as e
from libephemeris import state
from libephemeris.exceptions import ConfigurationError, DataNotFoundError, EphemerisRangeError
import pytest


FIXTURE = CORE_LEB
AFFINE = (2.0 ** -20, 2.0 ** -22, -(2.0 ** -19), 2.0 ** -23)
DELTA_T = ((0.0, 2.0 ** -10), (16.0, 3.0 * 2.0 ** -10))


def auxiliary_image(magic, nutation="affine", delta_t=DELTA_T):
    """Both formats store auxiliary sections uncompressed, after an empty index."""
    payloads = [(0, b"")]
    if nutation is not None:
        count = 0 if nutation == "empty" else 1
        interval = 0.0 if nutation == "bad-interval" else 16.0
        words = AFFINE if nutation != "constant" else (AFFINE[0], 0.0, AFFINE[2], 0.0)
        if nutation == "short-coefficients":
            words = words[:3]
        payload = struct.pack("<dddIIII", 0.0, 16.0, interval, 1, 2, count, 0)
        if count:
            payload += struct.pack(f"<{len(words)}d", *words)
        payloads.append((2, payload))
    if delta_t is not None:
        payload = struct.pack("<II", len(delta_t), 0)
        payload += b"".join(struct.pack("<dd", *row) for row in delta_t)
        payloads.append((3, payload))
    version, flags = (1, 0) if magic == "LEB1" else (2, 1)
    header = struct.pack("<4sIIIdddI20x", magic.encode(), version, len(payloads),
                         0, 0.0, 16.0, 0.0, flags)
    offset = len(header) + 24 * len(payloads)
    directory = []
    for section, payload in payloads:
        directory.append(struct.pack("<IIQQ", section, 0, offset, len(payload)))
        offset += len(payload)
    return header + b"".join(directory) + b"".join(payload for _, payload in payloads)


@pytest.fixture(params=("owned", "indexed"))
def backing_mode(request):
    return request.param


@pytest.fixture(params=("LEB1", "LEB2"))
def reader_factory(request, backing_mode, monkeypatch, tmp_path):
    monkeypatch.setenv("LIBEPHEMERIS_DATA_DIR", str(tmp_path))
    monkeypatch.setenv("LIBEPHEMERIS_CONFIG", str(tmp_path / "absent.toml"))
    monkeypatch.setenv("LIBEPHEMERIS_ENV_FILE", str(tmp_path / "absent.env"))
    monkeypatch.setenv("LIBEPHEMERIS_LEB_READER", backing_mode)
    e.close()
    e.reload_configuration()
    e.set_calc_mode("leb")
    e.set_precision_tier("base")
    readers = []

    def open_reader(**options):
        suffix = "leb" if request.param == "LEB1" else "leb2"
        path = tmp_path / f"auxiliary-{len(readers)}.{suffix}"
        path.write_bytes(auxiliary_image(request.param, **options))
        e.set_leb_file(str(path))
        reader = state.get_leb_reader()
        assert reader is not None and reader.path == str(path)
        readers.append(reader)
        return reader

    yield open_reader
    for reader in readers:
        reader.close()
    e.close()
    monkeypatch.undo()
    e.reload_configuration()
    e.set_calc_mode("leb")
    e.set_precision_tier("base")
    e.set_leb_file(str(FIXTURE))


def test_affine_nutation_radians_and_delta_t_days_are_exact(reader_factory):
    reader = reader_factory()
    assert reader.has_nutation() is True
    assert reader.eval_nutation(0.0) == (3.0 * 2.0 ** -22, -17.0 * 2.0 ** -23)
    assert reader.eval_nutation(8.0) == (2.0 ** -20, -(2.0 ** -19))
    assert reader.eval_nutation(16.0) == (5.0 * 2.0 ** -22, -15.0 * 2.0 ** -23)
    # A stored table is not the configured user override or an analytical model.
    previous = e.get_delta_t_userdef()
    try:
        e.set_delta_t_userdef(9.25)
        assert [reader.delta_t(jd) for jd in (-100.0, 0.0, 8.0, 16.0, 100.0)] == [
            2.0 ** -10, 2.0 ** -10, 2.0 ** -9, 3.0 * 2.0 ** -10, 3.0 * 2.0 ** -10]
    finally:
        e.set_delta_t_userdef(previous)


def test_constant_nutation_and_singleton_clock_remain_constant(reader_factory):
    reader = reader_factory(nutation="constant", delta_t=((8.0, 2.0 ** -11),))
    assert reader.has_nutation()
    for jd in (0.0, 3.0, 8.0, 16.0):
        assert reader.eval_nutation(jd) == (AFFINE[0], AFFINE[2])
    for jd in (-100.0, 8.0, 100.0):
        assert reader.delta_t(jd) == 2.0 ** -11


@pytest.mark.parametrize("nutation", (None, "empty", "bad-interval"))
def test_missing_or_unusable_nutation_does_not_fall_back(reader_factory, nutation):
    reader = reader_factory(nutation=nutation)
    assert reader.has_nutation() is False
    with pytest.raises(DataNotFoundError):
        reader.eval_nutation(8.0)
    assert reader.delta_t(8.0) == 2.0 ** -9


@pytest.mark.parametrize("delta_t", (None, ()))
def test_missing_or_empty_delta_t_does_not_fall_back(reader_factory, delta_t):
    reader = reader_factory(delta_t=delta_t)
    assert reader.has_nutation()
    with pytest.raises(DataNotFoundError):
        reader.delta_t(8.0)
    assert reader.eval_nutation(8.0) == (AFFINE[0], AFFINE[2])


def test_presence_does_not_certify_nutation_coefficient_payload(reader_factory):
    # Keep the short payload last so no other section can supply its missing word.
    reader = reader_factory(nutation="short-coefficients", delta_t=None)
    assert reader.has_nutation()
    with pytest.raises(DataNotFoundError):
        reader.eval_nutation(8.0)


def test_nutation_domain_and_nonfinite_requests_keep_reader_errors(reader_factory):
    reader = reader_factory()
    for jd in (-1.0, 17.0):
        with pytest.raises(EphemerisRangeError) as caught:
            reader.eval_nutation(jd)
        assert caught.value.requested_jd == jd
        assert caught.value.body_id == -1
        assert (caught.value.start_jd, caught.value.end_jd) == (0.0, 16.0)
    for jd in (float("nan"), float("inf"), -float("inf")):
        for evaluate in (reader.eval_nutation, reader.delta_t):
            with pytest.raises(ValueError):
                evaluate(jd)


def test_retained_handles_and_concurrent_auxiliary_reads_are_independent(reader_factory):
    first = reader_factory()
    second = state.get_leb_reader()
    assert second is not None
    first.close()
    for call in (first.has_nutation, lambda: first.eval_nutation(8.0), lambda: first.delta_t(8.0)):
        with pytest.raises(ConfigurationError):
            call()
    # A different selection and closing the runtime do not close retained handles.
    replacement = reader_factory(nutation="constant", delta_t=((8.0, 2.0 ** -11),))
    e.close()
    try:
        def read(_):
            return second.has_nutation(), second.eval_nutation(0.0), second.delta_t(8.0)
        with ThreadPoolExecutor(max_workers=4) as pool:
            assert list(pool.map(read, range(24))) == [
                (True, (3.0 * 2.0 ** -22, -17.0 * 2.0 ** -23), 2.0 ** -9)] * 24
        assert replacement.eval_nutation(0.0) == (AFFINE[0], AFFINE[2])
        assert replacement.delta_t(8.0) == 2.0 ** -11
    finally:
        second.close()
