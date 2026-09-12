# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Focused tests for the private B-11 perturber tuple certificate."""

from __future__ import annotations

import hashlib
import struct

import erfa
import pytest
from flint import ctx

from libephemeris.constants import JUPITER, NEPTUNE, SATURN, URANUS
from libephemeris.perturber_elements import (
    _canonical_coefficients_bytes,
    _certified_perturber_tuple,
    _frame_rows,
    _identity,
    _matrix_product,
    _transpose,
)
from libephemeris.planetary_mean_elements import (
    _MeanElementCoefficients,
    _mean_element_coefficients,
)

_IDS = (JUPITER, SATURN, URANUS, NEPTUNE)
_EXPECTED = {
    JUPITER: {
        "a": (5.202603191, 0.0000001913),
        "e": (0.04849485, 0.000163244, -0.0000004719, -0.00000000197),
        "i": (1.303270, -0.0054966, 0.00000465, -0.000000004),
        "node": (100.464441, 1.0209550, 0.00040117, 0.000000569),
        "perihelion": (14.331309, 1.6126668, 0.00103127, -0.000004569),
    },
    SATURN: {
        "a": (9.554909596, -0.0000021389),
        "e": (0.05550862, -0.000346818, -0.0000006456, 0.00000000338),
        "i": (2.488878, -0.0037363, -0.00001516, 0.000000089),
        "node": (113.665524, 0.8770979, -0.00012067, -0.000002380),
        "perihelion": (93.056787, 1.9637694, 0.00083757, 0.000004899),
    },
    URANUS: {
        "a": (19.218446062, -0.0000000372, 0.00000000098),
        "e": (0.04629590, -0.000027337, 0.0000000790, 0.00000000025),
        "i": (0.773196, 0.0007744, 0.00003749, -0.000000092),
        "node": (74.005947, 0.5211258, 0.00133982, 0.000018516),
        "perihelion": (173.005159, 1.4863784, 0.00021450, 0.000000433),
    },
    NEPTUNE: {
        "a": (30.110386869, -0.0000001663, 0.00000000069),
        "e": (0.00898809, 0.000006408, -0.0000000008, -0.00000000005),
        "i": (1.769952, -0.0093082, -0.00000708, 0.000000028),
        "node": (131.784057, 1.1022057, 0.00026006, -0.000000636),
        "perihelion": (48.123691, 1.4262677, 0.00037918, -0.000000003),
    },
}


def test_accessor_is_frozen_and_source_identical() -> None:
    """The private accessor exposes the sole registered coefficient table."""
    for body_id in _IDS:
        row = _mean_element_coefficients(body_id)
        assert isinstance(row, _MeanElementCoefficients)
        assert tuple(row.__dataclass_fields__) == (
            "body_id",
            "a",
            "e",
            "i",
            "node",
            "perihelion",
        )
        assert row.body_id == body_id
        for field, expected in _EXPECTED[body_id].items():
            value = getattr(row, field)
            assert value == expected
            assert all(type(item) is float for item in value)
    assert hashlib.sha256(_canonical_coefficients_bytes()).hexdigest() == (
        "f129d9b73f27d8d9bad8a81e54bd4eb25e375235039fd543e1f3d11402798bf5"
    )
    assert len(_canonical_coefficients_bytes()) == 852
    assert _mean_element_coefficients(0) is None
    with pytest.raises(TypeError):
        _mean_element_coefficients(True)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        _mean_element_coefficients("5")  # type: ignore[arg-type]


def test_canonical_encoding_uses_expected_record_order() -> None:
    """The digest payload has the normative body and field order."""
    payload = _canonical_coefficients_bytes()
    assert payload[:8] == struct.pack(">q", JUPITER)
    assert b"node" in payload and b"perihelion" in payload


def test_j2000_tuple_is_full_polynomial_and_native_float() -> None:
    """At J2000 the relative frame is identity and all five values are floats."""
    value = _certified_perturber_tuple(JUPITER, 2451545.0)
    assert value == (5.202603191, 0.04849485, 1.30327, 100.464441, 273.866868)
    assert all(type(item) is float for item in value)


def test_required_extrapolation_dates_are_certified() -> None:
    """The private gate handles the two conditioning witnesses."""
    for body_id in _IDS:
        for centuries in (-100.0, 10.0):
            value = _certified_perturber_tuple(
                body_id, 2451545.0 + 36525.0 * centuries
            )
            assert len(value) == 5
            assert all(type(item) is float for item in value)


def test_j2000_frame_is_literal_identity() -> None:
    """The exact epoch bypasses ERFA frame construction in both directions."""
    identity = _identity()
    assert _matrix_product(identity, _transpose(identity)) == identity
    assert identity[0][1].is_exact() and identity[0][1].is_zero()


def test_frame_rows_use_direct_erfa_poles(monkeypatch: pytest.MonkeyPatch) -> None:
    """The private frame path directly calls both ERFA pole routines."""
    calls: list[tuple[str, float]] = []
    original_ecl = erfa.ltpecl
    original_equ = erfa.ltpequ

    def ecl(epoch: float):
        calls.append(("ecl", epoch))
        return original_ecl(epoch)

    def equ(epoch: float):
        calls.append(("equ", epoch))
        return original_equ(epoch)

    monkeypatch.setattr(erfa, "ltpecl", ecl)
    monkeypatch.setattr(erfa, "ltpequ", equ)
    rows = _frame_rows(2451545.0 + 36525.0)
    assert [name for name, _ in calls] == ["ecl", "equ"]
    assert all(row[0].is_finite() for row in rows)
    assert _matrix_product(rows, _transpose(rows))


def test_precision_context_is_restored_after_tuple() -> None:
    """The two independent runtime passes do not leak Arb precision."""
    previous = ctx.prec
    _certified_perturber_tuple(SATURN, 2451545.0 + 36525.0 / 4.0)
    assert ctx.prec == previous
