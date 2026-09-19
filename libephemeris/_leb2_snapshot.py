# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Private immutable LEB2-v2 snapshot reader for a later source receipt.

The production factory admits only fixed base/medium core byte identities.
The test factory takes an explicit synthetic record and never grants an
approved asset tag. Both paths copy and hash before the shared LEB2 parser
or evaluator sees any byte. No ordinary ``LEB2Reader(path)`` mode changes.

Provenance:
    Project-authored implementation of the reviewed bounded snapshot-reader
    contract in ``leb2-immutable-snapshot-reader-contract.md`` (SHA-256
    811f128845e87c47bb6ce6a44268873a496b988b2b7e8b94f316de891b985c88).
    The pinned base/medium byte counts and SHA-256 values are the contract's
    computational asset identities. This module does not attest a selected
    tiered child, state result, or astronomical accuracy.
"""

from __future__ import annotations

import hashlib
import os
import struct
import threading
from dataclasses import dataclass
from types import MappingProxyType
from typing import Optional, Tuple

import zstandard as zstd

from .exceptions import LEBCorruptionError
from .leb2_reader import LEB2Reader
from .leb_compression import decompress_body
from .leb_format import LEB2_VERSION, StarEntry

_COPY_CHUNK_SIZE = 65536


@dataclass(frozen=True, slots=True)
class _SnapshotManifestRecord:
    """Expected identity supplied independently of the bytes being read."""

    role: str
    byte_count: int
    sha256: str


_PRODUCTION_MANIFEST = MappingProxyType(
    {
        "base": _SnapshotManifestRecord(
            "base_core.leb2",
            10_232_283,
            "5d708bdbe3e799e0802ba575984e57a3c5e44720dbfa1b4a01cf826640e0cb82",
        ),
        "medium": _SnapshotManifestRecord(
            "medium_core.leb2",
            37_276_175,
            "4d88ec9a79add7e3af9e75ac3ceabe5462a4af440447578a5ccba69a3e0a55b6",
        ),
    }
)


@dataclass(frozen=True, slots=True)
class _SnapshotAdmission:
    """Factory-owned identity of an already copied immutable byte object."""

    record: _SnapshotManifestRecord
    owned_bytes: bytes


@dataclass(frozen=True, slots=True)
class _ProductionAdmission(_SnapshotAdmission):
    """Capability issued only by the fixed production factory."""

    tier: str


class _SnapshotIntegrityError(ValueError):
    """A copied byte count or digest differs from the selected record."""


class _SnapshotCorruptionError(LEBCorruptionError):
    """Accepted bytes cannot be parsed or lazily decoded as LEB2-v2."""


class _SnapshotClosedError(ValueError):
    """A query requires an open immutable snapshot reader."""


def _validate_record(record: _SnapshotManifestRecord) -> None:
    """Reject malformed synthetic records before opening any descriptor."""
    if type(record) is not _SnapshotManifestRecord:
        raise TypeError("record must be _SnapshotManifestRecord")
    if type(record.role) is not str or not record.role:
        raise ValueError("record role must be a nonempty string")
    if type(record.byte_count) is not int or record.byte_count < 0:
        raise ValueError("record byte_count must be a nonnegative integer")
    if (
        type(record.sha256) is not str
        or len(record.sha256) != 64
        or any(c not in "0123456789abcdef" for c in record.sha256)
    ):
        raise ValueError("record sha256 must be a lowercase SHA-256 digest")


def _copy_verified(path: str, record: _SnapshotManifestRecord) -> bytes:
    """Copy at most count+1 bytes and verify that exact immutable sequence."""
    _validate_record(record)
    parts: list[bytes] = []
    remaining = record.byte_count + 1
    # The descriptor is closed on every exit path, before any parser call.
    with open(path, "rb") as source:
        while remaining:
            requested = min(_COPY_CHUNK_SIZE, remaining)
            part = source.read(requested)
            if type(part) is not bytes or len(part) > requested:
                raise _SnapshotIntegrityError("source read violated byte-copy contract")
            if not part:
                break
            parts.append(part)
            remaining -= len(part)
    snapshot = b"".join(parts)
    if len(snapshot) != record.byte_count:
        raise _SnapshotIntegrityError(
            f"snapshot byte count {len(snapshot)} != expected {record.byte_count}"
        )
    digest = hashlib.sha256(snapshot).hexdigest()
    if digest != record.sha256:
        raise _SnapshotIntegrityError("snapshot SHA-256 differs from selected record")
    return snapshot


class _LEB2SnapshotReader(LEB2Reader):
    """Byte-owned LEB2-v2 reader sharing the ordinary parser and evaluator."""

    _snapshot: bytes | None
    _private_decoder: zstd.ZstdDecompressor | None
    _admission: _SnapshotAdmission | None

    def __init__(self) -> None:
        raise TypeError("use a verified snapshot factory")

    @classmethod
    def _from_verified_bytes(
        cls,
        snapshot: bytes,
        path: str,
        record: _SnapshotManifestRecord,
    ) -> _LEB2SnapshotReader:
        """Build only after copy/hash acceptance and descriptor closure."""
        if type(snapshot) is not bytes or len(snapshot) != record.byte_count:
            raise _SnapshotIntegrityError("constructor requires exact verified bytes")
        if hashlib.sha256(snapshot).hexdigest() != record.sha256:
            raise _SnapshotIntegrityError(
                "constructor bytes differ from selected record"
            )
        reader = object.__new__(cls)
        reader._path = path
        reader._snapshot = snapshot
        reader._admission = None
        setattr(reader, "_file", None)
        setattr(reader, "_mm", snapshot)
        reader._cache = {}
        reader._chunk_cache = {}
        reader._decomp_lock = threading.Lock()
        reader._chunk_index = {}
        reader._chunked = False
        reader._eval_cache = {}
        reader._private_decoder = zstd.ZstdDecompressor()
        try:
            reader._parse()
            if reader._header.version != LEB2_VERSION or not reader._chunked:
                raise _SnapshotCorruptionError(
                    "immutable snapshot requires LEB2-v2 version"
                )
        except BaseException as exc:
            reader.close()
            if isinstance(exc, _SnapshotCorruptionError):
                raise
            if isinstance(exc, Exception):
                raise _SnapshotCorruptionError(
                    f"invalid LEB2-v2 snapshot {path!r}: {exc}"
                ) from exc
            raise
        return reader

    def _require_open(self) -> None:
        if self._snapshot is None:
            raise _SnapshotClosedError("immutable LEB2 snapshot reader is closed")

    def _checked_admission(self) -> _SnapshotAdmission:
        """Bind identity claims to the unchanged, factory-admitted bytes."""
        self._require_open()
        admission = self._admission
        if admission is None:
            raise _SnapshotIntegrityError("snapshot has no factory admission")
        owned = self._snapshot
        if (
            type(owned) is not bytes
            or owned is not admission.owned_bytes
            or self._mm is not owned
        ):
            raise _SnapshotIntegrityError("snapshot no longer owns admitted bytes")
        if isinstance(admission, _ProductionAdmission) != (
            type(self) is _ProductionLEB2SnapshotReader
        ):
            raise _SnapshotIntegrityError("snapshot admission class differs")
        if (
            len(owned) != admission.record.byte_count
            or hashlib.sha256(owned).hexdigest() != admission.record.sha256
        ):
            raise _SnapshotIntegrityError("owned bytes differ from admission")
        return admission

    @property
    def verified_digest(self) -> str:
        """Return the SHA-256 of the exact owned bytes while open."""
        return self._checked_admission().record.sha256

    @property
    def verified_byte_count(self) -> int:
        """Return the exact owned byte count while open."""
        return self._checked_admission().record.byte_count

    @property
    def manifest_identity(self) -> tuple[str, str]:
        """Return the logical tier and role; test identity is never approved."""
        admission = self._checked_admission()
        if isinstance(admission, _ProductionAdmission):
            return admission.tier, admission.record.role
        return "test-only", admission.record.role

    @property
    def production_asset_tag(self) -> tuple[str, str, int, str] | None:
        """Return an approved allowlist identity only for production admission."""
        admission = self._checked_admission()
        if not isinstance(admission, _ProductionAdmission):
            return None
        fixed = _PRODUCTION_MANIFEST.get(admission.tier)
        if fixed is None or admission.record != fixed:
            raise _SnapshotIntegrityError("production admission differs from allowlist")
        return (
            admission.tier,
            fixed.role,
            fixed.byte_count,
            fixed.sha256,
        )

    def __enter__(self) -> _LEB2SnapshotReader:
        self._require_open()
        return self

    @property
    def path(self) -> str:
        """Return only a diagnostic label, never an asset identity."""
        self._require_open()
        return self._path

    @property
    def jd_range(self) -> Tuple[float, float]:
        self._require_open()
        return super().jd_range

    def warm(self, jd_start: float, jd_end: float) -> None:
        """A live immutable bytes object needs no mmap prefetch advice."""
        self._require_open()

    def cool(self) -> None:
        """Do nothing for resident bytes, including after close."""

    def has_body(self, body_id: int) -> bool:
        self._require_open()
        return super().has_body(body_id)

    def body_coverage(self, body_id: int) -> Optional[Tuple[float, float]]:
        self._require_open()
        return super().body_coverage(body_id)

    def _find_chunk(self, body_id: int, jd: float, offset: float = 0.0) -> int:
        """Guard direct routing-metadata queries after close."""
        self._require_open()
        return super()._find_chunk(body_id, jd, offset)

    def _read_blob(self, offset: int, size: int, what: str) -> bytes:
        """Read only from live owned bytes, including direct internal calls."""
        self._require_open()
        return super()._read_blob(offset, size, what)

    def _decode_coefficients(
        self,
        compressed: bytes,
        uncompressed_size: int,
        segment_count: int,
        degree: int,
        components: int,
    ) -> bytes:
        """Use this snapshot's decoder while sharing the lossless transforms."""
        self._require_open()
        decoder = self._private_decoder
        if decoder is None:
            raise _SnapshotClosedError("snapshot decoder is closed")
        return decompress_body(
            compressed,
            uncompressed_size,
            segment_count,
            degree,
            components,
            decoder=decoder,
        )

    def _eval_body_split(
        self, body_id: int, jd: float, offset: float
    ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        self._require_open()
        try:
            return super()._eval_body_split(body_id, jd, offset)
        except (LEBCorruptionError, struct.error, IndexError) as exc:
            self.close()
            raise _SnapshotCorruptionError(f"corrupt LEB2-v2 body data: {exc}") from exc
        except KeyError as exc:
            if body_id not in self._bodies:
                raise  # An absent requested body is a caller query, not corruption.
            self.close()
            raise _SnapshotCorruptionError(
                f"corrupt LEB2-v2 body routing: {exc}"
            ) from exc
        except ValueError as exc:
            if not str(exc).startswith("Corrupted LEB2"):
                raise  # Preserve valid out-of-range and invalid-offset errors.
            self.close()
            raise _SnapshotCorruptionError(f"corrupt LEB2-v2 body data: {exc}") from exc

    def has_nutation(self) -> bool:
        self._require_open()
        return super().has_nutation()

    def eval_nutation(self, jd_tt: float) -> Tuple[float, float]:
        self._require_open()
        try:
            return super().eval_nutation(jd_tt)
        except (LEBCorruptionError, struct.error, IndexError) as exc:
            self.close()
            raise _SnapshotCorruptionError(
                f"corrupt LEB2-v2 nutation data: {exc}"
            ) from exc
        except ValueError as exc:
            if not str(exc).startswith("Corrupted LEB2"):
                raise
            self.close()
            raise _SnapshotCorruptionError(
                f"corrupt LEB2-v2 nutation data: {exc}"
            ) from exc

    def delta_t(self, jd: float) -> float:
        self._require_open()
        return super().delta_t(jd)

    def get_star(self, star_id: int) -> StarEntry:
        self._require_open()
        return super().get_star(star_id)

    def close(self) -> None:
        """Idempotently drop owned bytes, parsed metadata, and all caches."""
        self._admission = None
        with self._decomp_lock:
            self._cache.clear()
            self._chunk_cache.clear()
            self._eval_cache.clear()
            self._chunk_index.clear()
            for name in ("_sections", "_bodies", "_stars"):
                mapping = getattr(self, name, None)
                if mapping is not None:
                    mapping.clear()
            for name in ("_delta_t_jds", "_delta_t_vals"):
                values = getattr(self, name, None)
                if values is not None:
                    values.clear()
            setattr(self, "_mm", None)
            self._snapshot = None
            self._private_decoder = None
            self._nutation = None
            setattr(self, "_header", None)


class _ProductionLEB2SnapshotReader(_LEB2SnapshotReader):
    """Factory-only production object class distinct from test snapshots."""


def _open_production_snapshot(
    tier: str, path: str | os.PathLike[str]
) -> _LEB2SnapshotReader:
    """Admit only a fixed base/medium core identity at a local path."""
    if type(tier) is not str or tier not in _PRODUCTION_MANIFEST:
        raise ValueError("snapshot tier must be base or medium")
    record = _PRODUCTION_MANIFEST[tier]
    local_path = os.fspath(path)
    snapshot = _copy_verified(local_path, record)
    reader = _ProductionLEB2SnapshotReader._from_verified_bytes(
        snapshot, local_path, record
    )
    reader._admission = _ProductionAdmission(record, snapshot, tier)
    return reader


def _open_test_snapshot(
    path: str | os.PathLike[str], record: _SnapshotManifestRecord
) -> _LEB2SnapshotReader:
    """Exercise the same path with a synthetic record and no production tag."""
    _validate_record(record)
    local_path = os.fspath(path)
    snapshot = _copy_verified(local_path, record)
    reader = _LEB2SnapshotReader._from_verified_bytes(snapshot, local_path, record)
    reader._admission = _SnapshotAdmission(record, snapshot)
    return reader
