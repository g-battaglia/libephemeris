# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Focused byte-ownership and parity tests for the private LEB2 snapshot."""

from __future__ import annotations

import builtins
import hashlib
import os
import struct
from dataclasses import replace
from pathlib import Path
from typing import Any, cast

import pytest

import libephemeris._leb2_snapshot as snapshot_module
from libephemeris._leb2_snapshot import (
    _LEB2SnapshotReader,
    _SnapshotClosedError,
    _SnapshotCorruptionError,
    _SnapshotIntegrityError,
    _SnapshotManifestRecord,
    _open_production_snapshot,
    _open_test_snapshot,
)
from libephemeris.constants import EARTH, JUPITER, MOON, SATURN, SUN
from libephemeris.leb2_reader import LEB2Reader
from libephemeris.leb_format import (
    CHUNK_INDEX_HEADER_SIZE,
    LEB2_VERSION_V1,
    SECTION_NUTATION,
    write_chunk_entry,
    write_chunk_index_header,
    write_nutation_header,
)
from tests.test_leb.test_cov100_leb2_reader import JD0, _build_v2_file


def _record(data: bytes, role: str = "synthetic_core.leb2") -> _SnapshotManifestRecord:
    """Preselect a test-only manifest from generated fixture bytes."""
    return _SnapshotManifestRecord(role, len(data), hashlib.sha256(data).hexdigest())


def _bits(values: tuple[float, ...]) -> bytes:
    """Compare exact native float words without persisting state vectors."""
    assert all(type(word) is float for word in values)
    return struct.pack(f"!{len(values)}d", *values)


def _body_bits(reader: LEB2Reader, body_id: int, jd: float) -> bytes:
    pos, vel = reader.eval_body(body_id, jd)
    return _bits(pos + vel)


def _split_bits(reader: LEB2Reader, body_id: int, jd: float, offset: float) -> bytes:
    pos, vel = reader._eval_body_split(body_id, jd, offset)
    return _bits(pos + vel)


def _synthetic(tmp_path: Path, **kwargs: Any) -> tuple[Path, _SnapshotManifestRecord]:
    path = tmp_path / "synthetic_core.leb2"
    _build_v2_file(path, **kwargs)
    return path, _record(path.read_bytes())


def test_synthetic_snapshot_shares_exact_ordinary_evaluation(tmp_path: Path) -> None:
    path, record = _synthetic(
        tmp_path, with_nutation=True, with_delta_t=True, with_stars=True
    )
    with LEB2Reader(str(path)) as ordinary, _open_test_snapshot(path, record) as snap:
        assert snap.production_asset_tag is None
        assert snap.manifest_identity == ("test-only", record.role)
        assert snap.verified_digest == record.sha256
        assert snap.verified_byte_count == record.byte_count
        assert snap._snapshot is snap._mm
        assert snap.path == str(path)
        assert snap.jd_range == ordinary.jd_range
        assert snap.body_coverage(SUN) == ordinary.body_coverage(SUN)
        later_epochs = tuple(
            (chunk.jd_start + chunk.jd_end) / 2.0
            for chunk in ordinary._chunk_index[SUN][1:3]
        )
        for jd in (JD0 + 5.0, *later_epochs, JD0 + 5.0):
            assert _body_bits(snap, SUN, jd) == _body_bits(ordinary, SUN, jd)
        assert all((SUN, index) in snap._chunk_cache for index in range(3))
        for residual in (-(2.0**-40), 0.0, 2.0**-40):
            assert _split_bits(snap, SUN, JD0 + 10.0, residual) == _split_bits(
                ordinary, SUN, JD0 + 10.0, residual
            )
        assert _bits(snap.eval_nutation(JD0 + 5.0)) == _bits(
            ordinary.eval_nutation(JD0 + 5.0)
        )
        assert _bits((snap.delta_t(JD0 + 40.0),)) == _bits(
            (ordinary.delta_t(JD0 + 40.0),)
        )
        assert snap.get_star(1) == ordinary.get_star(1)
        snap.warm(JD0, JD0 + 20.0)
        snap.cool()


@pytest.mark.parametrize("tier", ["base", "medium"])
def test_production_snapshot_bitwise_matches_stable_ordinary(
    tier: str,
) -> None:
    path = Path("data/leb2") / f"{tier}_core.leb2"
    if not path.is_file():
        pytest.skip(f"local {tier} asset is unavailable")
    # Neither reader nor this test changes the ordinary backing bytes.
    with (
        LEB2Reader(str(path)) as ordinary,
        _open_production_snapshot(tier, path) as snap,
    ):
        assert snap.production_asset_tag is not None
        assert snap.production_asset_tag[0:2] == (tier, f"{tier}_core.leb2")
        for body_id in (SUN, MOON, EARTH, JUPITER, SATURN):
            assert ordinary.has_body(body_id) and snap.has_body(body_id)
            for jd in (2451545.0, 2460000.0, 2451545.0):
                assert _body_bits(snap, body_id, jd) == _body_bits(
                    ordinary, body_id, jd
                )
        boundary = ordinary._chunk_index[SUN][0].jd_end
        for offset in (-(2.0**-40), 0.0, 2.0**-40):
            assert _split_bits(snap, SUN, boundary, offset) == _split_bits(
                ordinary, SUN, boundary, offset
            )
        assert ordinary.has_nutation() and snap.has_nutation()
        assert _bits(snap.eval_nutation(2451545.0)) == _bits(
            ordinary.eval_nutation(2451545.0)
        )


@pytest.mark.parametrize("failure", ["short", "extra", "changed", "digest", "size"])
def test_copy_count_and_digest_reject_mutations(tmp_path: Path, failure: str) -> None:
    path, record = _synthetic(tmp_path)
    original = path.read_bytes()
    if failure == "short":
        path.write_bytes(original[:-1])
    elif failure == "extra":
        path.write_bytes(original + b"X")
    elif failure == "changed":
        changed = bytearray(original)
        changed[-1] ^= 1
        path.write_bytes(changed)
    elif failure == "digest":
        record = _SnapshotManifestRecord(record.role, record.byte_count, "0" * 64)
    else:
        record = _SnapshotManifestRecord(
            record.role, record.byte_count - 1, record.sha256
        )
    with pytest.raises(_SnapshotIntegrityError):
        _open_test_snapshot(path, record)


@pytest.mark.parametrize("bad_header", ["magic", "version", "compression"])
def test_matching_test_pin_reaches_parser_rejection(
    tmp_path: Path, bad_header: str
) -> None:
    kwargs: dict[str, object] = {}
    if bad_header == "magic":
        kwargs["magic"] = b"BADLEB2!"
    elif bad_header == "version":
        kwargs["version"] = LEB2_VERSION_V1
    else:
        kwargs["flags"] = 0
    path, record = _synthetic(tmp_path, **kwargs)
    with pytest.raises(_SnapshotCorruptionError, match=bad_header):
        _open_test_snapshot(path, record)


def test_matching_pin_reaches_lazy_corrupt_chunk_and_closes(tmp_path: Path) -> None:
    path, _ = _synthetic(tmp_path)
    with LEB2Reader(str(path)) as ordinary:
        chunk = ordinary._chunk_index[SUN][0]
    changed = bytearray(path.read_bytes())
    changed[chunk.blob_offset : chunk.blob_offset + chunk.compressed_size] = bytes(
        chunk.compressed_size
    )
    path.write_bytes(changed)
    snap = _open_test_snapshot(path, _record(bytes(changed)))
    with pytest.raises(_SnapshotCorruptionError):
        snap.eval_body(SUN, JD0 + 5.0)
    assert snap._snapshot is None
    assert snap._chunk_cache == snap._eval_cache == {}
    with pytest.raises(_SnapshotClosedError):
        snap.eval_body(SUN, JD0 + 5.0)


@pytest.mark.parametrize("fault", ["empty", "uncovered"])
def test_plain_value_error_routing_corruption_closes_with_typed_error(
    tmp_path: Path, fault: str
) -> None:
    path, _ = _synthetic(tmp_path, n_seg=10, chunk_segments=10)
    with LEB2Reader(str(path)) as ordinary:
        body = ordinary._bodies[SUN]
        chunk = ordinary._chunk_index[SUN][0]
    changed = bytearray(path.read_bytes())
    if fault == "empty":
        write_chunk_index_header(changed, body.data_offset, 0, 10.0)
    else:
        write_chunk_entry(
            changed,
            body.data_offset + CHUNK_INDEX_HEADER_SIZE,
            replace(chunk, segment_start=999),
        )
    path.write_bytes(changed)
    snap = _open_test_snapshot(path, _record(bytes(changed)))
    with pytest.raises(_SnapshotCorruptionError):
        snap.eval_body(SUN, JD0 + 5.0)
    assert snap._snapshot is None


def test_absent_body_and_out_of_range_remain_caller_errors(tmp_path: Path) -> None:
    path, record = _synthetic(tmp_path)
    with _open_test_snapshot(path, record) as snap:
        with pytest.raises(KeyError):
            snap.eval_body(999, JD0 + 5.0)
        with pytest.raises(ValueError, match="outside range"):
            snap.eval_body(SUN, JD0 - 1.0)
        assert snap.verified_digest == record.sha256


def test_lazy_nutation_corruption_closes_with_typed_error(tmp_path: Path) -> None:
    path, _ = _synthetic(tmp_path, with_nutation=True)
    with LEB2Reader(str(path)) as ordinary:
        offset = ordinary._sections[SECTION_NUTATION].offset
        nutation = ordinary._nutation
        assert nutation is not None
    changed = bytearray(path.read_bytes())
    write_nutation_header(changed, offset, replace(nutation, interval_days=0.0))
    path.write_bytes(changed)
    snap = _open_test_snapshot(path, _record(bytes(changed)))
    with pytest.raises(_SnapshotCorruptionError):
        snap.eval_nutation(JD0 + 5.0)
    assert snap._snapshot is None


class _TrackedFile:
    """Track descriptor closure and optionally fail the first read."""

    def __init__(self, path: str, fail_read: bool = False) -> None:
        self.handle = builtins.open(path, "rb", buffering=0)
        self.fail_read = fail_read

    def __enter__(self) -> _TrackedFile:
        return self

    def __exit__(self, *args: object) -> None:
        self.handle.close()

    def read(self, count: int) -> bytes:
        if self.fail_read:
            raise OSError("injected read failure")
        return self.handle.read(count)

    @property
    def closed(self) -> bool:
        return self.handle.closed


@pytest.mark.parametrize("failure", ["none", "read", "count", "digest", "parse"])
def test_descriptor_is_closed_before_parse_and_on_every_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    path, record = _synthetic(tmp_path)
    if failure == "count":
        record = _SnapshotManifestRecord(
            record.role, record.byte_count + 1, record.sha256
        )
    elif failure == "digest":
        record = _SnapshotManifestRecord(record.role, record.byte_count, "0" * 64)
    elif failure == "parse":
        data = bytearray(path.read_bytes())
        data[0] ^= 1
        path.write_bytes(data)
        record = _record(bytes(data))
    tracked: list[_TrackedFile] = []

    def tracked_open(local_path: str, mode: str) -> _TrackedFile:
        assert mode == "rb"
        source = _TrackedFile(local_path, failure == "read")
        tracked.append(source)
        return source

    monkeypatch.setattr(snapshot_module, "open", tracked_open, raising=False)
    original_parse = _LEB2SnapshotReader._parse
    parsed: list[_LEB2SnapshotReader] = []

    def checked_parse(reader: _LEB2SnapshotReader) -> None:
        assert tracked and all(source.closed for source in tracked)
        parsed.append(reader)
        original_parse(reader)

    monkeypatch.setattr(_LEB2SnapshotReader, "_parse", checked_parse)
    if failure == "none":
        snap = _open_test_snapshot(path, record)
        assert tracked and tracked[0].closed
        snap.close()
    else:
        expected = OSError if failure == "read" else ValueError
        with pytest.raises(expected):
            _open_test_snapshot(path, record)
        assert tracked and tracked[0].closed
        if failure == "parse":
            assert parsed and parsed[0]._snapshot is None


@pytest.mark.parametrize("already_copied", [True, False])
def test_deterministic_midcopy_mutation_checks_actual_copied_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, already_copied: bool
) -> None:
    path, record = _synthetic(tmp_path)
    monkeypatch.setattr(snapshot_module, "_COPY_CHUNK_SIZE", 32)

    class _MutatingFile(_TrackedFile):
        first = True

        def read(self, count: int) -> bytes:
            part = super().read(count)
            if self.first and part:
                self.first = False
                offset = 0 if already_copied else 40
                with builtins.open(path, "r+b") as target:
                    target.seek(offset)
                    previous = target.read(1)
                    target.seek(offset)
                    target.write(bytes((previous[0] ^ 1,)))
            return part

    monkeypatch.setattr(
        snapshot_module,
        "open",
        lambda local_path, mode: _MutatingFile(local_path),
        raising=False,
    )
    if already_copied:
        with _open_test_snapshot(path, record) as snap:
            assert snap.verified_digest == record.sha256
            assert _body_bits(snap, SUN, JD0 + 5.0)
    else:
        with pytest.raises(_SnapshotIntegrityError):
            _open_test_snapshot(path, record)


def test_replacement_and_inplace_change_do_not_change_live_snapshot(
    tmp_path: Path,
) -> None:
    path, record = _synthetic(
        tmp_path, with_nutation=True, with_delta_t=True, with_stars=True
    )
    with LEB2Reader(str(path)) as ordinary:
        later = ordinary._chunk_index[SUN][2]
        later_epoch = (later.jd_start + later.jd_end) / 2.0
        expected_uncached = _body_bits(ordinary, SUN, later_epoch)
    snap = _open_test_snapshot(path, record)
    before = _body_bits(snap, SUN, JD0 + 5.0)
    replacement = tmp_path / "replacement.leb2"
    replacement.write_bytes(b"bad replacement")
    os.replace(replacement, path)
    assert _body_bits(snap, SUN, JD0 + 5.0) == before
    assert _body_bits(snap, SUN, later_epoch) == expected_uncached
    assert len(_bits(snap.eval_nutation(JD0 + 5.0))) == 16
    assert snap.delta_t(JD0 + 40.0) > 0.0
    assert snap.get_star(1).star_id == 1
    with pytest.raises(_SnapshotIntegrityError):
        _open_test_snapshot(path, record)
    snap.close()

    path, record = _synthetic(tmp_path)
    with LEB2Reader(str(path)) as ordinary:
        later = ordinary._chunk_index[SUN][2]
        later_epoch = (later.jd_start + later.jd_end) / 2.0
        expected_uncached = _body_bits(ordinary, SUN, later_epoch)
    snap = _open_test_snapshot(path, record)
    before = _body_bits(snap, SUN, JD0 + 5.0)
    assert (SUN, 2) not in snap._chunk_cache
    changed = bytearray(path.read_bytes())
    changed[later.blob_offset] ^= 1
    path.write_bytes(changed)
    assert _body_bits(snap, SUN, JD0 + 5.0) == before
    assert _body_bits(snap, SUN, later_epoch) == expected_uncached
    assert (SUN, 2) in snap._chunk_cache
    with pytest.raises(_SnapshotIntegrityError):
        _open_test_snapshot(path, record)
    snap.close()


def test_production_allowlist_and_test_only_identity(tmp_path: Path) -> None:
    path, record = _synthetic(tmp_path)
    with pytest.raises(ValueError):
        _open_production_snapshot("extended", path)
    with pytest.raises(_SnapshotIntegrityError):
        _open_production_snapshot("base", path)
    with _open_test_snapshot(
        path,
        _SnapshotManifestRecord("base_core.leb2", record.byte_count, record.sha256),
    ) as snap:
        assert snap.production_asset_tag is None
        assert snap.manifest_identity == ("test-only", "base_core.leb2")
    synthetic_bytes = path.read_bytes()
    forged_role = _SnapshotManifestRecord(
        "base_core.leb2", record.byte_count, record.sha256
    )
    direct = _LEB2SnapshotReader._from_verified_bytes(
        synthetic_bytes, str(path), forged_role
    )
    with pytest.raises(_SnapshotIntegrityError, match="no factory admission"):
        _ = direct.production_asset_tag
    direct.close()
    base_path = Path("libephemeris/data/leb2/base_core.leb2")
    if base_path.is_file():
        disposable = tmp_path / "other-name.leb2"
        disposable.write_bytes(base_path.read_bytes())
        with _open_production_snapshot("base", disposable) as snap:
            assert snap.production_asset_tag is not None
            assert snap.path == str(disposable)
        with pytest.raises(_SnapshotIntegrityError):
            _open_production_snapshot("medium", disposable)
        disposable.write_bytes(b"changed")
        with pytest.raises(_SnapshotIntegrityError):
            _open_production_snapshot("base", disposable)


@pytest.mark.parametrize(
    "name,value",
    [
        ("_production_tier", "base"),
        ("_manifest_role", "base_core.leb2"),
        (
            "_verified_digest",
            "5d708bdbe3e799e0802ba575984e57a3c5e44720dbfa1b4a01cf826640e0cb82",
        ),
        ("_verified_byte_count", 10_232_283),
    ],
)
def test_test_only_metadata_mutation_cannot_gain_production_tag(
    tmp_path: Path, name: str, value: object
) -> None:
    path, record = _synthetic(tmp_path)
    with _open_test_snapshot(path, record) as snap:
        setattr(snap, name, value)
        assert snap.production_asset_tag is None
        assert snap.manifest_identity == ("test-only", record.role)
        assert snap.verified_digest == record.sha256
        assert snap.verified_byte_count == record.byte_count


def test_fixed_allowlist_and_production_admission_ignore_mutable_fields() -> None:
    with pytest.raises(TypeError):
        cast(Any, snapshot_module._PRODUCTION_MANIFEST)["base"] = (
            _SnapshotManifestRecord("forged", 0, "0" * 64)
        )
    path = Path("data/leb2/base_core.leb2")
    if not path.is_file():
        pytest.skip("local base asset is unavailable")
    with _open_production_snapshot("base", path) as snap:
        approved = snap.production_asset_tag
        assert approved is not None
        setattr(snap, "_production_tier", "medium")
        setattr(snap, "_manifest_role", "medium_core.leb2")
        setattr(snap, "_verified_digest", "0" * 64)
        setattr(snap, "_verified_byte_count", 0)
        assert snap.production_asset_tag == approved
        assert snap.manifest_identity == ("base", "base_core.leb2")
        setattr(snap, "_snapshot", b"not the admitted bytes")
        with pytest.raises(_SnapshotIntegrityError):
            _ = snap.production_asset_tag


def test_postclose_guards_caches_and_context_exit(tmp_path: Path) -> None:
    path, record = _synthetic(
        tmp_path, with_nutation=True, with_delta_t=True, with_stars=True
    )
    with _open_test_snapshot(path, record) as snap:
        snap.eval_body(SUN, JD0 + 5.0)
        assert snap._chunk_cache and snap._eval_cache
    assert snap._snapshot is None and snap._mm is None
    assert not snap._chunk_cache and not snap._eval_cache and not snap._bodies
    snap.close()
    snap.cool()
    queries = (
        lambda: snap.eval_body(SUN, JD0 + 5.0),
        lambda: snap._eval_body_split(SUN, JD0 + 5.0, 2.0**-40),
        lambda: snap.eval_nutation(JD0 + 5.0),
        lambda: snap.delta_t(JD0 + 5.0),
        lambda: snap.get_star(1),
        lambda: snap.has_body(SUN),
        lambda: snap.body_coverage(SUN),
        lambda: snap.has_nutation(),
        lambda: snap._find_chunk(SUN, JD0 + 5.0),
        lambda: snap._read_blob(0, 1, "probe"),
        lambda: snap.jd_range,
        lambda: snap.path,
        lambda: snap.verified_digest,
        lambda: snap.verified_byte_count,
        lambda: snap.manifest_identity,
        lambda: snap.production_asset_tag,
        lambda: snap.warm(JD0, JD0 + 10.0),
    )
    for query in queries:
        with pytest.raises(_SnapshotClosedError):
            query()
