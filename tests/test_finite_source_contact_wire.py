# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Byte-level fixtures for the disconnected exact contact request codec."""

from __future__ import annotations

import hashlib
from fractions import Fraction as F

import pytest

from libephemeris._finite_source_certificates import (
    _DomainOutcome,
    _OutsideDomain,
    _Target,
    _ValidatedGeometry,
)
from libephemeris._finite_source_contact_wire import (
    _ContactWireError,
    _ContactWireLimits,
    _decode_contact_request,
    _encode_contact_request,
)

_WIRE = b"W" * 32
_GENERATOR = b"G" * 32
_POLICY = b"P" * 32
_LIMITS = _ContactWireLimits(max_frame_bytes=1024, max_integer_bits=32)


def _geometry(
    B: tuple[F, F, F],
    p: tuple[F, F, F],
    r: F,
    m: F,
    diagonal: tuple[F, F, F],
    R: F,
) -> _ValidatedGeometry:
    zero = F(0)
    return _ValidatedGeometry(
        B,
        p,
        r,
        m,
        (
            (diagonal[0], zero, zero),
            (zero, diagonal[1], zero),
            (zero, zero, diagonal[2]),
        ),
        R,
    )


_SMOOTH = _geometry(
    (F(-15), F(10), F(0)),
    (F(-10), F(10), F(0)),
    F(2),
    F(1),
    (F(1), F(1), F(1)),
    F(1),
)


def _encode(geometry: _ValidatedGeometry, target: _Target) -> bytes:
    return _encode_contact_request(
        geometry,
        target,
        wire_sha256=_WIRE,
        generator_sha256=_GENERATOR,
        policy_sha256=_POLICY,
        limits=_LIMITS,
    )


def _decode(frame: bytes, limits: _ContactWireLimits = _LIMITS):
    return _decode_contact_request(
        frame,
        expected_wire_sha256=_WIRE,
        expected_generator_sha256=_GENERATOR,
        expected_policy_sha256=_POLICY,
        limits=limits,
    )


# Independently specified bytes for the smooth fixture's rational fields.
# Each literal is sign:U8, magnitude length:U32, magnitude bytes, repeated
# for numerator and denominator. The header length below is a fixed fixture.
_Q_ZERO = bytes.fromhex("0000000000 000000000101")
_Q_ONE = bytes.fromhex("000000000101 000000000101")
_Q_TWO = bytes.fromhex("000000000102 000000000101")
_Q_TEN = bytes.fromhex("00000000010a 000000000101")
_Q_NEG_TEN = bytes.fromhex("01000000010a 000000000101")
_Q_NEG_THIRTEEN = bytes.fromhex("01000000010d 000000000101")
_Q_NEG_FIFTEEN = bytes.fromhex("01000000010f 000000000101")
_Q_TWENTY = bytes.fromhex("000000000114 000000000101")
_SMOOTH_BODY = (
    _WIRE
    + _GENERATOR
    + _POLICY
    + b"\x01"
    + b"".join(
        (
            _Q_NEG_FIFTEEN,
            _Q_TEN,
            _Q_ZERO,
            _Q_NEG_TEN,
            _Q_TEN,
            _Q_ZERO,
            _Q_TWO,
            _Q_ONE,
            _Q_ONE,
            _Q_ZERO,
            _Q_ZERO,
            _Q_ZERO,
            _Q_ONE,
            _Q_ZERO,
            _Q_ZERO,
            _Q_ZERO,
            _Q_ONE,
            _Q_ONE,
        )
    )
)
_SMOOTH_FRAME = b"LEFC1\x01\x00\x00\x00\x00\x00\x00\x01\x31" + _SMOOTH_BODY


def _with_body(body: bytes) -> bytes:
    return b"LEFC1\x01" + len(body).to_bytes(8, "big") + body


def test_independent_smooth_byte_fixture_and_request_identity() -> None:
    assert len(_SMOOTH_BODY) == 305
    assert _encode(_SMOOTH, _Target.PENUMBRA) == _SMOOTH_FRAME
    decoded = _decode(_SMOOTH_FRAME)
    assert decoded.geometry == _SMOOTH
    assert decoded.target is _Target.PENUMBRA
    assert decoded.request_sha256 == hashlib.sha256(_SMOOTH_FRAME).digest()
    assert _encode(decoded.geometry, decoded.target) == _SMOOTH_FRAME


@pytest.mark.parametrize(
    ("geometry", "target"),
    [
        (_SMOOTH, _Target.ANTUMBRA),
        (
            _geometry(
                (F(-55, 12), F(5, 4), F(0)),
                (F(5, 12), F(5, 4), F(0)),
                F(11, 4),
                F(1, 4),
                (F(1), F(1), F(1)),
                F(1),
            ),
            _Target.PENUMBRA,
        ),
        (
            _geometry(
                (F(-57, 20), F(1), F(0)),
                (F(3, 20), F(1), F(0)),
                F(19, 10),
                F(1, 10),
                (F(4), F(5, 4), F(4)),
                F(9, 10),
            ),
            _Target.UMBRA,
        ),
        (
            _geometry(
                (F(-11), F(0), F(0)),
                (F(-6), F(0), F(0)),
                F(2),
                F(1),
                (F(1), F(1), F(1)),
                F(1),
            ),
            _Target.UMBRA,
        ),
    ],
)
def test_exact_request_round_trip_across_contact_branches(
    geometry: _ValidatedGeometry, target: _Target
) -> None:
    decoded = _decode(_encode(geometry, target))
    assert decoded.geometry == geometry
    assert decoded.target is target
    assert _encode(decoded.geometry, decoded.target) == _encode(geometry, target)


@pytest.mark.parametrize(
    "frame",
    [
        b"",
        _SMOOTH_FRAME[:13],
        b"LEFC2" + _SMOOTH_FRAME[5:],
        _SMOOTH_FRAME[:5] + b"\x02" + _SMOOTH_FRAME[6:],
        _SMOOTH_FRAME[:6] + b"\x00" * 7 + b"\x01" + _SMOOTH_BODY,
        _SMOOTH_FRAME[:-1],
        _SMOOTH_FRAME + b"\x00",
        _with_body(_SMOOTH_BODY[:96] + b"\x00" + _SMOOTH_BODY[97:]),
        _with_body(_SMOOTH_BODY + _Q_ONE),
    ],
)
def test_rejects_wrong_version_kind_length_target_or_extra_fields(frame: bytes) -> None:
    with pytest.raises(_ContactWireError):
        _decode(frame)


@pytest.mark.parametrize("digest_offset", [0, 32, 64])
def test_expected_identities_are_external_and_exact(digest_offset: int) -> None:
    changed = bytearray(_SMOOTH_FRAME)
    changed[14 + digest_offset] ^= 1
    with pytest.raises(_ContactWireError, match="identity mismatch"):
        _decode(bytes(changed))


@pytest.mark.parametrize(
    ("body", "message"),
    [
        (_SMOOTH_BODY[:97] + b"\x02" + _SMOOTH_BODY[98:], "sign"),
        (_SMOOTH_BODY[:121] + b"\x01" + _SMOOTH_BODY[122:], "negative zero"),
        (_SMOOTH_BODY[:102] + b"\x00" + _SMOOTH_BODY[103:], "leading zero"),
        (_SMOOTH_BODY[:108] + b"\x03" + _SMOOTH_BODY[109:], "not reduced"),
        (_SMOOTH_BODY[:103] + b"\x01" + _SMOOTH_BODY[104:], "denominator"),
        (
            _SMOOTH_BODY[:103] + b"\x00\x00\x00\x00\x00" + _SMOOTH_BODY[109:],
            "denominator",
        ),
        (
            _SMOOTH_BODY[:97] + b"\x00\x00\x00\x00\x02\x00\x0f" + _SMOOTH_BODY[103:],
            "leading zero",
        ),
    ],
)
def test_rejects_noncanonical_integer_and_fraction_forms(
    body: bytes, message: str
) -> None:
    with pytest.raises(_ContactWireError, match=message):
        _decode(_with_body(body))


def test_rejects_oversize_before_arithmetic_and_accepts_exact_limits() -> None:
    frame = _SMOOTH_FRAME
    assert len(frame) == 319
    assert _decode(frame, _ContactWireLimits(319, 4)).geometry == _SMOOTH
    with pytest.raises(_ContactWireError, match="frame byte limit"):
        _decode(frame, _ContactWireLimits(318, 4))
    with pytest.raises(_ContactWireError, match="bit limit"):
        _decode(frame, _ContactWireLimits(319, 3))
    with pytest.raises(_ContactWireError, match="bit limit"):
        _encode_contact_request(
            _SMOOTH,
            _Target.PENUMBRA,
            wire_sha256=_WIRE,
            generator_sha256=_GENERATOR,
            policy_sha256=_POLICY,
            limits=_ContactWireLimits(319, 3),
        )
    with pytest.raises(_ContactWireError, match="frame byte limit"):
        _encode_contact_request(
            _SMOOTH,
            _Target.PENUMBRA,
            wire_sha256=_WIRE,
            generator_sha256=_GENERATOR,
            policy_sha256=_POLICY,
            limits=_ContactWireLimits(318, 4),
        )


def test_claimed_magnitude_length_is_rejected_before_reading_it() -> None:
    changed = bytearray(_SMOOTH_FRAME)
    changed[112:116] = b"\xff" * 4
    with pytest.raises(_ContactWireError, match="truncated integer magnitude"):
        _decode(bytes(changed))


def test_geometry_input_and_domain_exceptions_remain_distinct() -> None:
    invalid = _with_body(_SMOOTH_BODY[: -len(_Q_ONE)] + _Q_ZERO)
    with pytest.raises(ValueError, match="R must be positive"):
        _decode(invalid)

    for body, outcome in (
        (
            _SMOOTH_BODY[:179] + _Q_TWO + _SMOOTH_BODY[191:],
            _DomainOutcome.UNSUPPORTED_SLOPE,
        ),
        (
            _SMOOTH_BODY[:132] + _Q_NEG_THIRTEEN + _SMOOTH_BODY[144:],
            _DomainOutcome.EXCLUDED_DOMAIN,
        ),
        (
            _SMOOTH_BODY[: -len(_Q_ONE)] + _Q_TWENTY,
            _DomainOutcome.UNCERTIFIED_DOMAIN,
        ),
    ):
        with pytest.raises(_OutsideDomain) as exc:
            _decode(_with_body(body))
        assert exc.value.outcome is outcome


def test_encoder_revalidates_forged_base_geometry_before_returning_bytes() -> None:
    forged = _geometry(
        (F(-15), F(10), F(0)),
        (F(-10), F(10), F(0)),
        F(2),
        F(1),
        (F(1), F(1), F(1)),
        F(1),
    )
    object.__setattr__(forged, "R", F(0))
    with pytest.raises(ValueError, match="R must be positive"):
        _encode(forged, _Target.PENUMBRA)


def test_invalid_limit_and_call_types_rejected() -> None:
    with pytest.raises(ValueError):
        _ContactWireLimits(13, 8)
    with pytest.raises(ValueError):
        _ContactWireLimits(319, 0)
    with pytest.raises(TypeError):
        _decode_contact_request(
            bytearray(_SMOOTH_FRAME),  # type: ignore[arg-type]
            expected_wire_sha256=_WIRE,
            expected_generator_sha256=_GENERATOR,
            expected_policy_sha256=_POLICY,
            limits=_LIMITS,
        )
