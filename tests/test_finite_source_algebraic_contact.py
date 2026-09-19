# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Focused synthetic checks of the disconnected exact contact verifier."""

from __future__ import annotations

from fractions import Fraction as F

import pytest

from libephemeris._finite_source_algebraic_contact import (
    _AlgebraicContactPayload,
    _RootSign,
    _verify_algebraic_contact,
)
from libephemeris._finite_source_certificates import _Target, _ValidatedGeometry

_I = ((F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(1)))
_ZERO = (F(0),)


def _geometry(B: tuple[F, F, F], p: tuple[F, F, F], r: F, m: F) -> _ValidatedGeometry:
    """Construct a synthetic unit-solid exact geometry."""
    return _ValidatedGeometry(B, p, r, m, _I, F(1))


def _payload(
    tau: F,
    beta: tuple[F, F, F],
    gate: F | None,
    *,
    P: tuple[int, ...] = (0, 1),
    lo: F = F(-1),
    hi: F = F(1),
) -> _AlgebraicContactPayload:
    """Encode rational dual coefficients in one selected algebraic root."""
    return _AlgebraicContactPayload(
        P=P,
        lo=lo,
        hi=hi,
        tau=(tau,),
        beta=((beta[0],), (beta[1],), (beta[2],)),
        gate_multiplier=None if gate is None else (gate,),
    )


def test_smooth_penumbral_contact_has_exact_dual() -> None:
    """The nonsingular smooth contact accepts only its tangent support."""
    geometry = _geometry((F(-15), F(10), F(0)), (F(-10), F(10), F(0)), F(2), F(1))
    payload = _payload(F(2, 25), (F(0), F(0), F(8, 25)), F(0))
    assert _verify_algebraic_contact(geometry, _Target.PENUMBRA, payload)


def test_positive_singular_penumbral_contact_is_not_lost() -> None:
    """A determinant-zero KKT branch has a valid exact contact payload."""
    geometry = _geometry(
        (F(-55, 12), F(5, 4), F(0)),
        (F(5, 12), F(5, 4), F(0)),
        F(11, 4),
        F(1, 4),
    )
    payload = _payload(F(2, 25), (F(0), F(0), F(8, 25)), F(0))
    assert _verify_algebraic_contact(geometry, _Target.PENUMBRA, payload)


def test_core_apex_contact_is_per_target() -> None:
    """The umbral apex touches while the opposite nappe reaches strictly."""
    geometry = _geometry((F(-11), F(0), F(0)), (F(-6), F(0), F(0)), F(2), F(1))
    umbral = _payload(F(2, 5), (F(0), F(0), F(0)), F(0))
    antumbral_bad = _payload(F(-2, 5), (F(0), F(0), F(0)), None)
    assert _verify_algebraic_contact(geometry, _Target.UMBRA, umbral)
    assert not _verify_algebraic_contact(geometry, _Target.ANTUMBRA, antumbral_bad)
    assert geometry.strict_reach(_Target.ANTUMBRA, (F(-9, 10), F(0), F(0)))


def test_irrational_root_and_reducible_polynomial_signs() -> None:
    """Exact gcd and interval refinement decide signs at selected sqrt(2)."""
    root = _RootSign((F(-2), F(0), F(1)), F(1), F(2))
    assert root.sign((F(-2), F(0), F(1))) == 0
    assert root.sign((F(-1), F(1))) == 1
    assert root.sign((F(-2), F(1))) == -1

    reducible = _RootSign((F(6), F(-2), F(-3), F(1)), F(1), F(2))
    assert reducible.sign((F(-2), F(0), F(1))) == 0
    assert reducible.sign((F(-3), F(1))) == -1


def test_higher_degree_root_can_encode_rational_dual() -> None:
    """The verifier does not impose a degree-one producer convention."""
    geometry = _geometry((F(-15), F(10), F(0)), (F(-10), F(10), F(0)), F(2), F(1))
    payload = _payload(
        F(2, 25),
        (F(0), F(0), F(8, 25)),
        F(0),
        P=(-2, 0, 1),
        lo=F(1),
        hi=F(2),
    )
    assert _verify_algebraic_contact(geometry, _Target.PENUMBRA, payload)
    reducible = _payload(
        F(2, 25),
        (F(0), F(0), F(8, 25)),
        F(0),
        P=(6, -2, -3, 1),
        lo=F(1),
        hi=F(2),
    )
    assert _verify_algebraic_contact(geometry, _Target.PENUMBRA, reducible)


@pytest.mark.parametrize(
    ("P", "lo", "hi"),
    [
        ((1, -2, 1), F(0), F(2)),  # Repeated root.
        ((-1, 1), F(1), F(2)),  # Root at left endpoint.
        ((-1, 1), F(2), F(3)),  # No root.
        ((-2, 0, 1), F(-2), F(2)),  # Two roots.
        ((0, 2), F(-1), F(1)),  # Nonprimitive content.
        ((0, -1), F(-1), F(1)),  # Wrong leading sign.
    ],
)
def test_invalid_root_certificate_returns_false(
    P: tuple[int, ...], lo: F, hi: F
) -> None:
    """A valid mathematical contact cannot repair an invalid root encoding."""
    geometry = _geometry((F(-15), F(10), F(0)), (F(-10), F(10), F(0)), F(2), F(1))
    payload = _payload(F(2, 25), (F(0), F(0), F(8, 25)), F(0), P=P, lo=lo, hi=hi)
    assert not _verify_algebraic_contact(geometry, _Target.PENUMBRA, payload)


def test_malformed_payload_is_distinct_from_false_certificate() -> None:
    """Exact coefficient syntax and field degree are structural obligations."""
    geometry = _geometry((F(-15), F(10), F(0)), (F(-10), F(10), F(0)), F(2), F(1))
    valid = _payload(F(2, 25), (F(0), F(0), F(8, 25)), F(0))
    with pytest.raises(TypeError):
        _verify_algebraic_contact(
            geometry,
            _Target.PENUMBRA,
            _AlgebraicContactPayload(
                (0, 1),
                F(-1),
                F(1),
                (0,),  # type: ignore[arg-type]
                valid.beta,
                _ZERO,
            ),
        )
    with pytest.raises(ValueError, match="trailing zero"):
        _verify_algebraic_contact(
            geometry,
            _Target.PENUMBRA,
            _AlgebraicContactPayload(
                (0, 1), F(-1), F(1), (F(2, 25), F(0)), valid.beta, _ZERO
            ),
        )
    with pytest.raises(ValueError, match="degree"):
        _verify_algebraic_contact(
            geometry,
            _Target.PENUMBRA,
            _AlgebraicContactPayload(
                (0, 1), F(-1), F(1), (F(2, 25), F(1)), valid.beta, _ZERO
            ),
        )
    with pytest.raises(TypeError, match="gate multiplier"):
        _verify_algebraic_contact(
            geometry,
            _Target.ANTUMBRA,
            _AlgebraicContactPayload((0, 1), F(-1), F(1), valid.tau, valid.beta, _ZERO),
        )


def test_failed_dual_support_and_target_checks_return_false() -> None:
    """Each independent certificate condition can reject an exact payload."""
    geometry = _geometry((F(-15), F(10), F(0)), (F(-10), F(10), F(0)), F(2), F(1))
    assert not _verify_algebraic_contact(
        geometry,
        _Target.PENUMBRA,
        _payload(F(-2, 25), (F(0), F(0), F(8, 25)), F(0)),
    )
    assert not _verify_algebraic_contact(
        geometry,
        _Target.PENUMBRA,
        _payload(F(2, 25), (F(0), F(0), F(8, 25)), F(-1)),
    )
    assert not _verify_algebraic_contact(
        geometry,
        _Target.PENUMBRA,
        _payload(F(2, 25), (F(0), F(0), F(1)), F(0)),
    )
    assert not _verify_algebraic_contact(
        geometry,
        _Target.PENUMBRA,
        _payload(F(2, 25), (F(0), F(0), F(3, 10)), F(0)),
    )

    separated = _geometry(
        (F(-43, 12), F(2), F(0)),
        (F(17, 12), F(2), F(0)),
        F(11, 4),
        F(1, 4),
    )
    tangent_but_outside = _payload(F(2, 15), (F(0), F(0), F(0)), F(0))
    assert not _verify_algebraic_contact(
        separated, _Target.PENUMBRA, tangent_but_outside
    )


def test_penumbral_full_cone_apex_fails_the_physical_gate() -> None:
    """A valid tangent to the full cone is rejected before the Moon."""
    geometry = _geometry(
        (F(-239, 60), F(4, 5), F(0)),
        (F(61, 60), F(4, 5), F(0)),
        F(11, 4),
        F(1, 4),
    )
    payload = _payload(F(2, 25), (F(0), F(0), F(8, 25)), F(0))
    assert not _verify_algebraic_contact(geometry, _Target.PENUMBRA, payload)


def test_wrong_penumbral_nappe_fails_even_with_tangent_support() -> None:
    """A negative cone side is rejected despite zero squared residual."""
    geometry = _geometry(
        (F(-119, 30), F(63, 80), F(0)),
        (F(31, 30), F(63, 80), F(0)),
        F(11, 4),
        F(1, 4),
    )
    tangent_point = (F(3, 5), F(4, 5), F(0))
    a, gate, h_t, _ = geometry._target_values(_Target.PENUMBRA, tangent_point)
    assert a == F(-1, 4) and a * a == h_t and gate == F(-17, 12)
    payload = _payload(F(2, 25), (F(0), F(0), F(8, 25)), F(0))
    assert not _verify_algebraic_contact(geometry, _Target.PENUMBRA, payload)


def test_derived_geometry_fields_are_recomputed() -> None:
    """A forged derived inverse cannot influence exact contact acceptance."""
    geometry = _geometry((F(-15), F(10), F(0)), (F(-10), F(10), F(0)), F(2), F(1))
    object.__setattr__(geometry, "A_inverse", ((_ZERO, _ZERO, _ZERO),) * 3)
    payload = _payload(F(2, 25), (F(0), F(0), F(8, 25)), F(0))
    assert _verify_algebraic_contact(geometry, _Target.PENUMBRA, payload)
