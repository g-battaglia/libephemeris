# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Disconnected in-memory codec for one exact contact request.

This implements only the request half of the reviewed non-normative private
wire design candidate. The caller must supply bytes already received into
memory; this module does not enforce a bound while reading a transport.
It does not start a worker, decode a result, authenticate astronomical
states, or authorize a physical contact. Limits are explicit caller inputs;
no production resource policy has been selected. A synchronous parse alone
does not bound native exact arithmetic or the whole invocation.

Provenance:
    Project-authored canonical encoding from the independently reviewed,
    non-normative B-16 wire design candidate in the private audit. Exact
    rational geometry is reconstructed by the existing private domain
    constructor. No external reference implementation is used.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from fractions import Fraction

from ._finite_source_certificates import _Matrix, _Target, _ValidatedGeometry, _Vector

_MAGIC = b"LEFC1"
_REQUEST_KIND = 1
_HEADER_BYTES = 14
_DIGEST_BYTES = 32
_U32_MAX = (1 << 32) - 1
_U64_MAX = (1 << 64) - 1

_TARGET_TO_TAG = {
    _Target.PENUMBRA: 1,
    _Target.UMBRA: 2,
    _Target.ANTUMBRA: 3,
}
_TAG_TO_TARGET = {tag: target for target, tag in _TARGET_TO_TAG.items()}


class _ContactWireError(ValueError):
    """A request has invalid or noncanonical framing or exact-number bytes."""


@dataclass(frozen=True, slots=True)
class _ContactWireLimits:
    """Explicit size ceilings for a disconnected request codec call."""

    max_frame_bytes: int
    max_integer_bits: int

    def __post_init__(self) -> None:
        if type(self.max_frame_bytes) is not int or not (
            _HEADER_BYTES <= self.max_frame_bytes <= _U64_MAX
        ):
            raise ValueError("max_frame_bytes must fit the U64 frame length")
        if type(self.max_integer_bits) is not int or not (
            1 <= self.max_integer_bits <= 8 * _U32_MAX
        ):
            raise ValueError("max_integer_bits must fit the U32 magnitude length")


@dataclass(frozen=True, slots=True)
class _DecodedContactRequest:
    """One canonical request after size-checked parsing and domain validation.

    The digest distinguishes request frames. It is not a source receipt.
    """

    geometry: _ValidatedGeometry
    target: _Target
    request_sha256: bytes


def _require_digest(value: bytes, name: str) -> bytes:
    """Require an exact 32-byte identity supplied outside parsed payloads."""
    if type(value) is not bytes or len(value) != _DIGEST_BYTES:
        raise TypeError(f"{name} must be exactly 32 bytes")
    return value


def _append(buffer: bytearray, data: bytes, max_body_bytes: int) -> None:
    """Append a fixed field only when it fits the caller's body ceiling."""
    if len(data) > max_body_bytes - len(buffer):
        raise _ContactWireError("request exceeds the frame byte limit")
    buffer.extend(data)


def _append_integer(buffer: bytearray, value: int, limits: _ContactWireLimits) -> None:
    """Append one canonical signed magnitude after bit and byte prechecks."""
    if type(value) is not int:
        raise TypeError("wire integer must be an exact int")
    bits = value.bit_length()
    if bits > limits.max_integer_bits:
        raise _ContactWireError("integer exceeds the bit limit")
    size = (bits + 7) // 8
    if size > _U32_MAX or 5 + size > limits.max_frame_bytes - _HEADER_BYTES - len(
        buffer
    ):
        raise _ContactWireError("integer exceeds the frame byte limit")
    buffer.append(1 if value < 0 else 0)
    buffer.extend(size.to_bytes(4, "big"))
    if size:
        buffer.extend(abs(value).to_bytes(size, "big"))


def _append_fraction(
    buffer: bytearray, value: Fraction, limits: _ContactWireLimits
) -> Fraction:
    """Append and detach one canonical rational after bounded integer writes."""
    if type(value) is not Fraction:
        raise TypeError("wire rational must be an exact Fraction")
    numerator, denominator = value.numerator, value.denominator
    _append_integer(buffer, numerator, limits)
    _append_integer(buffer, denominator, limits)
    if denominator <= 0 or math.gcd(abs(numerator), denominator) != 1:
        raise ValueError("geometry rational must be reduced with positive denominator")
    return Fraction(numerator, denominator)


def _append_vector(
    buffer: bytearray,
    vector: _Vector,
    limits: _ContactWireLimits,
) -> _Vector:
    """Append a three-coordinate vector and return detached rational values."""
    if type(vector) is not tuple or len(vector) != 3:
        raise TypeError("geometry vector must be a three-tuple")
    return (
        _append_fraction(buffer, vector[0], limits),
        _append_fraction(buffer, vector[1], limits),
        _append_fraction(buffer, vector[2], limits),
    )


def _append_matrix(
    buffer: bytearray,
    matrix: _Matrix,
    limits: _ContactWireLimits,
) -> _Matrix:
    """Append exactly three vector rows and retain their detached values."""
    if type(matrix) is not tuple or len(matrix) != 3:
        raise TypeError("geometry matrix must have three rows")
    return (
        _append_vector(buffer, matrix[0], limits),
        _append_vector(buffer, matrix[1], limits),
        _append_vector(buffer, matrix[2], limits),
    )


def _target_tag(target: _Target) -> int:
    """Select the fixed private target tag without enum coercion."""
    if type(target) is not _Target:
        raise TypeError("target must be an exact private _Target")
    return _TARGET_TO_TAG[target]


def _body(
    geometry: _ValidatedGeometry,
    target: _Target,
    wire_sha256: bytes,
    generator_sha256: bytes,
    policy_sha256: bytes,
    limits: _ContactWireLimits,
) -> bytes:
    """Encode base fields and revalidate them before releasing the body."""
    if type(geometry) is not _ValidatedGeometry:
        raise TypeError("geometry must be an exact _ValidatedGeometry")
    max_body_bytes = limits.max_frame_bytes - _HEADER_BYTES
    body = bytearray()
    for name, digest in (
        ("wire_sha256", wire_sha256),
        ("generator_sha256", generator_sha256),
        ("policy_sha256", policy_sha256),
    ):
        _append(body, _require_digest(digest, name), max_body_bytes)
    _append(body, bytes((_target_tag(target),)), max_body_bytes)
    B = _append_vector(body, geometry.B, limits)
    p = _append_vector(body, geometry.p, limits)
    r = _append_fraction(body, geometry.r, limits)
    m = _append_fraction(body, geometry.m, limits)
    A = _append_matrix(body, geometry.A, limits)
    R = _append_fraction(body, geometry.R, limits)
    _ValidatedGeometry(B, p, r, m, A, R)
    return bytes(body)


def _encode_contact_request(
    geometry: _ValidatedGeometry,
    target: _Target,
    *,
    wire_sha256: bytes,
    generator_sha256: bytes,
    policy_sha256: bytes,
    limits: _ContactWireLimits,
) -> bytes:
    """Encode one exact geometry; no worker or astronomical source is read."""
    if type(limits) is not _ContactWireLimits:
        raise TypeError("limits must be an exact _ContactWireLimits")
    body = _body(geometry, target, wire_sha256, generator_sha256, policy_sha256, limits)
    return _MAGIC + bytes((_REQUEST_KIND,)) + len(body).to_bytes(8, "big") + body


class _Cursor:
    """Read a size-checked in-memory frame without copying magnitude slices."""

    __slots__ = ("data", "offset", "limits")

    def __init__(self, data: memoryview, limits: _ContactWireLimits) -> None:
        self.data = data
        self.offset = 0
        self.limits = limits

    def take(self, size: int) -> memoryview:
        if size < 0 or size > len(self.data) - self.offset:
            raise _ContactWireError("truncated request field")
        start = self.offset
        self.offset += size
        return self.data[start : self.offset]

    def integer(self) -> int:
        sign = int(self.take(1)[0])
        size = int.from_bytes(self.take(4), "big")
        if sign not in (0, 1):
            raise _ContactWireError("invalid integer sign")
        if size > len(self.data) - self.offset:
            raise _ContactWireError("truncated integer magnitude")
        if size > (self.limits.max_integer_bits + 7) // 8:
            raise _ContactWireError("integer exceeds the bit limit")
        if size == 0:
            if sign:
                raise _ContactWireError("negative zero is not canonical")
            return 0
        first = int(self.data[self.offset])
        if first == 0:
            raise _ContactWireError("integer magnitude has a leading zero")
        bits = 8 * (size - 1) + first.bit_length()
        if bits > self.limits.max_integer_bits:
            raise _ContactWireError("integer exceeds the bit limit")
        magnitude = int.from_bytes(self.take(size), "big")
        return -magnitude if sign else magnitude

    def fraction(self) -> Fraction:
        numerator = self.integer()
        denominator = self.integer()
        if denominator <= 0 or math.gcd(abs(numerator), denominator) != 1:
            raise _ContactWireError("rational is not reduced with positive denominator")
        return Fraction(numerator, denominator)


def _decode_contact_request(
    frame: bytes,
    *,
    expected_wire_sha256: bytes,
    expected_generator_sha256: bytes,
    expected_policy_sha256: bytes,
    limits: _ContactWireLimits,
) -> _DecodedContactRequest:
    """Parse one size-checked request and reconstruct its exact geometry.

    Expected identities must come from outside the frame. Framing errors
    raise ``_ContactWireError``; geometry validation retains its existing
    invalid-input and out-of-domain exception types.
    """
    if type(limits) is not _ContactWireLimits:
        raise TypeError("limits must be an exact _ContactWireLimits")
    expected = (
        _require_digest(expected_wire_sha256, "expected_wire_sha256"),
        _require_digest(expected_generator_sha256, "expected_generator_sha256"),
        _require_digest(expected_policy_sha256, "expected_policy_sha256"),
    )
    if type(frame) is not bytes:
        raise TypeError("frame must be exact bytes")
    if len(frame) > limits.max_frame_bytes:
        raise _ContactWireError("request exceeds the frame byte limit")
    if len(frame) < _HEADER_BYTES or frame[:5] != _MAGIC:
        raise _ContactWireError("invalid or truncated request header")
    if frame[5] != _REQUEST_KIND:
        raise _ContactWireError("unexpected frame kind")
    size = int.from_bytes(frame[6:_HEADER_BYTES], "big")
    if (
        size > limits.max_frame_bytes - _HEADER_BYTES
        or size != len(frame) - _HEADER_BYTES
    ):
        raise _ContactWireError("request frame length mismatch")
    cursor = _Cursor(memoryview(frame)[_HEADER_BYTES:], limits)
    for name, expected_digest in zip(
        ("wire", "generator", "policy"), expected, strict=True
    ):
        if bytes(cursor.take(_DIGEST_BYTES)) != expected_digest:
            raise _ContactWireError(f"{name} identity mismatch")
    target_tag = int(cursor.take(1)[0])
    target = _TAG_TO_TARGET.get(target_tag)
    if target is None:
        raise _ContactWireError("invalid target tag")
    values = tuple(cursor.fraction() for _ in range(18))
    if cursor.offset != len(cursor.data):
        raise _ContactWireError("trailing request bytes")
    B = (values[0], values[1], values[2])
    p = (values[3], values[4], values[5])
    A = (
        (values[8], values[9], values[10]),
        (values[11], values[12], values[13]),
        (values[14], values[15], values[16]),
    )
    geometry = _ValidatedGeometry(B, p, values[6], values[7], A, values[17])
    return _DecodedContactRequest(geometry, target, hashlib.sha256(frame).digest())
