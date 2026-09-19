# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Disconnected, exact post-return input model for finite-source proofs.

Only synthetic native Python float words are accepted here. This module
checks their arithmetic consistency with the fixed product revision and
derives exact rational geometry. It does not authenticate any call, body
role, epoch, state, provider, asset, or successful astronomical result.

Provenance:
    Project-authored implementation of the reviewed bounded input-model
    specification ``finite-source-post-return-input-model.md`` (SHA-256
    e9d3937b290ca2bf7c355b762891cb7de3d3d072063b28638793938aaa837c36)
    for product revision 1147e6a907b35d6d52613e56c1986c094a5897c9.
    Every declared binary64 operation is replayed on fixed constant words;
    the resulting kilometre words are converted individually to Fractions.
"""

from __future__ import annotations

import math
import struct
from dataclasses import dataclass
from fractions import Fraction
from typing import TypeAlias

from libephemeris._finite_source_certificates import (
    _DomainOutcome,
    _ValidatedGeometry,
    _prepare_geometry,
)

_NativeVector: TypeAlias = tuple[float, float, float]
_Vector: TypeAlias = tuple[Fraction, Fraction, Fraction]
_Matrix: TypeAlias = tuple[_Vector, _Vector, _Vector]
_PRODUCT_REVISION = "1147e6a907b35d6d52613e56c1986c094a5897c9"
_IAU_AU_KM = Fraction(1495978707, 10)

# Binary64 words at the fixed product revision. The intermediate and final
# pins guard the native operation chains as well as their source literals.
_AU_KM_BITS = 0x41A1D55D5D666666
_SUN_RADIUS_KM_BITS = 0x41253D8000000000
_MOON_RADIUS_KM_BITS = 0x409B28999999999A
_EARTH_AXIS_KM_BITS = 0x40B8EA23D70A3D71
_ONE_BITS = 0x3FF0000000000000
_FLATTENING_DENOMINATOR_BITS = 0x4072A41A4BDBA0A5
_SUN_RADIUS_AU_BITS = 0x3F730E789D2665C2
_MOON_RADIUS_AU_BITS = 0x3EE85DCDC9F0E0F8
_EARTH_AXIS_AU_BITS = 0x3F065A676F29FFA0
_FLATTENING_BITS = 0x3F6B775F5E758B82


class _InvalidInputModel(ValueError):
    """A finite native word or fixed-revision arithmetic relation is invalid."""


@dataclass(frozen=True, slots=True)
class _NativeConstants:
    """Captured native constant and operation-result words, not a receipt."""

    au_km: float
    sun_radius_km: float
    moon_radius_km: float
    earth_axis_km: float
    flattening_unit: float
    flattening_denominator: float
    sun_radius_au: float
    moon_radius_au: float
    earth_axis_au: float
    r_native: float
    m_native: float
    a_native: float
    f_native: float


@dataclass(frozen=True, slots=True)
class _NativeWords:
    """Captured float words labelled by body; labels are not authenticated."""

    tjd_ut: float
    moon_au: _NativeVector
    moon_km: _NativeVector
    sun_au: _NativeVector
    sun_km: _NativeVector
    constants: _NativeConstants


def _bits(word: float) -> int:
    """Return all 64 IEEE-754 bits, including a zero's sign bit."""
    return struct.unpack("!Q", struct.pack("!d", word))[0]


def _word_from_bits(bits: int) -> float:
    """Reconstruct one fixed-revision binary64 constant word."""
    return struct.unpack("!d", struct.pack("!Q", bits))[0]


def _finite_word(value: object, name: str) -> float:
    """Accept only a finite native Python float without numerical coercion."""
    if type(value) is not float:
        raise TypeError(f"{name} must be a native Python float")
    if not math.isfinite(value):
        raise _InvalidInputModel(f"{name} must be finite")
    return value


def _native_vector(value: object, name: str) -> _NativeVector:
    """Check a captured three-component native float vector."""
    if type(value) is not tuple or len(value) != 3:
        raise TypeError(f"{name} must be a three-tuple of native floats")
    return (
        _finite_word(value[0], f"{name}[0]"),
        _finite_word(value[1], f"{name}[1]"),
        _finite_word(value[2], f"{name}[2]"),
    )


def _equal_word(actual: float, expected: float, name: str) -> None:
    """Reject any differing binary64 encoding, including signed zero."""
    if _bits(actual) != _bits(expected):
        raise _InvalidInputModel(f"{name} differs from fixed native chain")


def _check_constants(constants: _NativeConstants) -> None:
    """Replay and pin every fixed-revision native constant operation."""
    if type(constants) is not _NativeConstants:
        raise TypeError("constants must be _NativeConstants")
    names = tuple(_NativeConstants.__dataclass_fields__)
    for name in names:
        _finite_word(getattr(constants, name), f"constants.{name}")

    c = _word_from_bits(_AU_KM_BITS)
    sun_km = _word_from_bits(_SUN_RADIUS_KM_BITS)
    moon_km = _word_from_bits(_MOON_RADIUS_KM_BITS)
    earth_km = _word_from_bits(_EARTH_AXIS_KM_BITS)
    one = _word_from_bits(_ONE_BITS)
    denominator = _word_from_bits(_FLATTENING_DENOMINATOR_BITS)
    fixed_words = (
        (constants.au_km, c, "au_km"),
        (constants.sun_radius_km, sun_km, "sun_radius_km"),
        (constants.moon_radius_km, moon_km, "moon_radius_km"),
        (constants.earth_axis_km, earth_km, "earth_axis_km"),
        (constants.flattening_unit, one, "flattening_unit"),
        (constants.flattening_denominator, denominator, "flattening_denominator"),
    )
    for actual, expected, name in fixed_words:
        _equal_word(actual, expected, name)

    # These are separate Python binary64 divide/multiply operations. Compare
    # the native replay against pinned words before comparing captured words.
    sun_au = sun_km / c
    moon_au = moon_km / c
    earth_au = earth_km / c
    r_native = sun_au * c
    m_native = moon_au * c
    a_native = earth_au * c
    f_native = one / denominator
    pinned_chain = (
        (sun_au, _SUN_RADIUS_AU_BITS, "sun_radius_au"),
        (moon_au, _MOON_RADIUS_AU_BITS, "moon_radius_au"),
        (earth_au, _EARTH_AXIS_AU_BITS, "earth_axis_au"),
        (r_native, _SUN_RADIUS_KM_BITS, "r_native"),
        (m_native, _MOON_RADIUS_KM_BITS, "m_native"),
        (a_native, _EARTH_AXIS_KM_BITS, "a_native"),
        (f_native, _FLATTENING_BITS, "f_native"),
    )
    for computed, pinned_bits, name in pinned_chain:
        if _bits(computed) != pinned_bits:
            raise _InvalidInputModel(f"native {name} chain changed from fixed revision")
        _equal_word(getattr(constants, name), computed, name)

    if Fraction.from_float(c) != _IAU_AU_KM - Fraction(1, 83886080):
        raise _InvalidInputModel("fixed AU constant rational identity failed")
    if Fraction.from_float(r_native) != Fraction(696000):
        raise _InvalidInputModel("fixed Sun radius identity failed")
    if Fraction.from_float(m_native) - Fraction(173815, 100) != Fraction(
        1, 10995116277760
    ):
        raise _InvalidInputModel("fixed Moon radius identity failed")
    if Fraction.from_float(a_native) - Fraction(6378140, 1000) != Fraction(
        9, 27487790694400
    ):
        raise _InvalidInputModel("fixed Earth axis identity failed")
    if Fraction.from_float(f_native) - Fraction(50000, 14912821) != Fraction(
        3751413, 17193312025252584327479296
    ):
        raise _InvalidInputModel("fixed flattening identity failed")


@dataclass(frozen=True, slots=True, init=False)
class _FiniteSourceInput:
    """Validated words and exact derived geometry for one synthetic input."""

    words: _NativeWords
    B: _Vector
    p: _Vector
    r: Fraction
    m: Fraction
    A: _Matrix
    R: Fraction
    u: _Vector

    def __init__(self, words: _NativeWords) -> None:
        if type(words) is not _NativeWords:
            raise TypeError("words must be _NativeWords")
        _finite_word(words.tjd_ut, "tjd_ut")
        moon_au = _native_vector(words.moon_au, "moon_au")
        moon_km = _native_vector(words.moon_km, "moon_km")
        sun_au = _native_vector(words.sun_au, "sun_au")
        sun_km = _native_vector(words.sun_km, "sun_km")
        _check_constants(words.constants)
        c = words.constants.au_km
        for role, au, km in (("moon", moon_au, moon_km), ("sun", sun_au, sun_km)):
            for i in range(3):
                expected = au[i] * c
                _finite_word(expected, f"{role}_km[{i}] native product")
                _equal_word(km[i], expected, f"{role}_km[{i}]")

        p: _Vector = (
            Fraction.from_float(moon_km[0]),
            Fraction.from_float(moon_km[1]),
            Fraction.from_float(moon_km[2]),
        )
        B: _Vector = (
            Fraction.from_float(sun_km[0]),
            Fraction.from_float(sun_km[1]),
            Fraction.from_float(sun_km[2]),
        )
        r = Fraction.from_float(words.constants.r_native)
        m = Fraction.from_float(words.constants.m_native)
        a = Fraction.from_float(words.constants.a_native)
        f = Fraction.from_float(words.constants.f_native)
        if a <= 0 or not 0 < f < 1:
            raise _InvalidInputModel("native Earth axis and flattening are invalid")
        b = a * (1 - f)
        equatorial = 1 / (a * a)
        polar = 1 / (b * b)
        zero = Fraction(0)
        A: _Matrix = (
            (equatorial, zero, zero),
            (zero, equatorial, zero),
            (zero, zero, polar),
        )
        object.__setattr__(self, "words", words)
        object.__setattr__(self, "p", p)
        object.__setattr__(self, "B", B)
        object.__setattr__(self, "r", r)
        object.__setattr__(self, "m", m)
        object.__setattr__(self, "A", A)
        object.__setattr__(self, "R", a)
        object.__setattr__(self, "u", (p[0] - B[0], p[1] - B[1], p[2] - B[2]))

    @property
    def fixed_product_revision(self) -> str:
        """Name the chain revision; this is not a source receipt."""
        return _PRODUCT_REVISION

    @property
    def native_span(self) -> _NativeVector:
        """Diagnose the legacy rounded subtraction without using it as u."""
        return (
            self.words.moon_km[0] - self.words.sun_km[0],
            self.words.moon_km[1] - self.words.sun_km[1],
            self.words.moon_km[2] - self.words.sun_km[2],
        )

    @property
    def iau_differences(self) -> tuple[_Vector, _Vector]:
        """Exact native-km minus exact-IAU differences, Moon then Sun."""
        moon: _Vector = (
            self.p[0] - Fraction.from_float(self.words.moon_au[0]) * _IAU_AU_KM,
            self.p[1] - Fraction.from_float(self.words.moon_au[1]) * _IAU_AU_KM,
            self.p[2] - Fraction.from_float(self.words.moon_au[2]) * _IAU_AU_KM,
        )
        sun: _Vector = (
            self.B[0] - Fraction.from_float(self.words.sun_au[0]) * _IAU_AU_KM,
            self.B[1] - Fraction.from_float(self.words.sun_au[1]) * _IAU_AU_KM,
            self.B[2] - Fraction.from_float(self.words.sun_au[2]) * _IAU_AU_KM,
        )
        return moon, sun

    def theorem_domain(self) -> _ValidatedGeometry | _DomainOutcome:
        """Delegate theorem-domain outcomes to the disconnected checker."""
        return _prepare_geometry(self.B, self.p, self.r, self.m, self.A, self.R)
