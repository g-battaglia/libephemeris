# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Focused exact fixtures for the disconnected finite-source generator."""

from __future__ import annotations

from fractions import Fraction as F

import pytest

from libephemeris._finite_source_certificates import (
    _DualPayload,
    _Target,
    _ValidatedGeometry,
)
from libephemeris._finite_source_contact_algebra import (
    _NumberField,
    _compare_distinct_roots,
    _real_roots,
)
from libephemeris._finite_source_contact_extension import _quadratic_real_roots
from libephemeris._finite_source_contact_generator import (
    _CandidateDisposition,
    _CandidatePredicate,
    _ContactBranch,
    _ContactGenerationInvariantError,
    _ContactGenerationStatus,
    _SummaryReason,
    _forms,
    _generate_algebraic_contact,
    _quadratic_coefficients,
    _stationarity_polynomials,
)


def _diagonal(a: F, b: F, c: F) -> tuple[tuple[F, F, F], ...]:
    zero = F(0)
    return ((a, zero, zero), (zero, b, zero), (zero, zero, c))


def _smooth() -> _ValidatedGeometry:
    zero = F(0)
    return _ValidatedGeometry(
        (F(-15), F(10), zero),
        (F(-10), F(10), zero),
        F(2),
        F(1),
        _diagonal(F(1), F(1), F(1)),
        F(1),
    )


def _singular_rational() -> _ValidatedGeometry:
    zero = F(0)
    return _ValidatedGeometry(
        (F(-55, 12), F(5, 4), zero),
        (F(5, 12), F(5, 4), zero),
        F(11, 4),
        F(1, 4),
        _diagonal(F(1), F(1), F(1)),
        F(1),
    )


def _singular_quadratic() -> _ValidatedGeometry:
    zero = F(0)
    return _ValidatedGeometry(
        (F(-57, 20), F(1), zero),
        (F(3, 20), F(1), zero),
        F(19, 10),
        F(1, 10),
        _diagonal(F(4), F(5, 4), F(4)),
        F(9, 10),
    )


def _core_apex() -> _ValidatedGeometry:
    zero = F(0)
    return _ValidatedGeometry(
        (F(-11), zero, zero),
        (F(-6), zero, zero),
        F(2),
        F(1),
        _diagonal(F(1), F(1), F(1)),
        F(1),
    )


def test_rational_smooth_contact_has_complete_private_trace() -> None:
    geometry = _smooth()
    result = _generate_algebraic_contact(geometry, _Target.PENUMBRA)
    assert result.status is _ContactGenerationStatus.CONTACT
    assert result.payload is not None
    assert result.payload.P == (-1, 4600)
    assert result.payload.tau == (F(2, 25),)
    assert result.payload.beta[2] == (F(8, 25),)
    assert tuple(item.summary.branch for item in result.trace) == tuple(_ContactBranch)
    assert tuple(item.summary.completed for item in result.trace) == (True,) * 3
    assert result.trace[0].summary.root_count == 1
    assert result.trace[0].summary.point_count == 1
    assert result.trace[0].attempts[0].disposition is _CandidateDisposition.ACCEPTED
    assert result == _generate_algebraic_contact(_smooth(), _Target.PENUMBRA)

    forms = _forms(geometry, _Target.PENUMBRA)
    C, d, e = _quadratic_coefficients(forms)
    _, _, Pq, Pf = _stationarity_polynomials(geometry, C, d, e)
    _, factors = Pq.gcd(Pf).factor()
    assert len(factors) == 2  # Reducible gcd and repeated negative factor.


def test_positive_singular_root_is_transferred_and_both_points_checked() -> None:
    result = _generate_algebraic_contact(_singular_rational(), _Target.PENUMBRA)
    assert result.status is _ContactGenerationStatus.CONTACT
    assert result.payload is not None
    assert result.payload.P == (-1, 225)
    transfer = result.trace[0].attempts[0]
    assert transfer.disposition is _CandidateDisposition.TRANSFER_SINGULAR
    assert transfer.predicate is _CandidatePredicate.D_ZERO
    assert transfer.point is None
    assert transfer.root == result.trace[1].summary.root
    assert [item.disposition for item in result.trace[1].attempts] == [
        _CandidateDisposition.FAILED_SIDE,
        _CandidateDisposition.ACCEPTED,
    ]
    assert result.trace[1].attempts[0].predicate_value == (F(-9),)


def test_quadratic_singular_contact_uses_canonical_common_field() -> None:
    result = _generate_algebraic_contact(_singular_quadratic(), _Target.PENUMBRA)
    assert result.status is _ContactGenerationStatus.CONTACT
    assert result.payload is not None
    assert result.payload.P == (-61, -360, 1620)
    assert (result.payload.lo, result.payload.hi) == (F(1, 4), F(3, 4))
    assert result.payload.tau == (F(-4, 27), F(4, 3))
    assert result.payload.beta[2] == (F(2, 3),)
    first, second = result.trace[1].attempts
    assert first.ordinal == 0 and second.ordinal == 1
    assert first.disposition is _CandidateDisposition.FAILED_SIDE
    assert first.predicate is _CandidatePredicate.SIDE_NONPOSITIVE
    assert second.disposition is _CandidateDisposition.ACCEPTED
    assert result.trace[1].summary.point_count == 2


def test_common_core_apex_preserves_opposite_target_rejection() -> None:
    geometry = _core_apex()
    umbra = _generate_algebraic_contact(geometry, _Target.UMBRA)
    antumbra = _generate_algebraic_contact(geometry, _Target.ANTUMBRA)
    assert umbra.status is _ContactGenerationStatus.CONTACT
    assert umbra.trace[2].attempts[0].disposition is _CandidateDisposition.ACCEPTED
    assert antumbra.status is _ContactGenerationStatus.NO_CONTACT_EXHAUSTED
    attempt = antumbra.trace[2].attempts[0]
    assert attempt.disposition is _CandidateDisposition.FAILED_APEX_DUAL
    assert attempt.predicate is _CandidatePredicate.TAU_NEGATIVE
    assert attempt.predicate_value == (F(-2, 5),)
    assert attempt.payload is not None
    assert all(item.summary.completed for item in antumbra.trace)


def test_smooth_verifier_disagreement_is_an_invariant_fault(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import libephemeris._finite_source_contact_generator as generator

    monkeypatch.setattr(generator, "_verify_algebraic_contact", lambda *_: False)
    with pytest.raises(_ContactGenerationInvariantError, match="smooth dual"):
        _generate_algebraic_contact(_smooth(), _Target.PENUMBRA)


def test_generator_reconstructs_forged_derived_geometry() -> None:
    geometry = _smooth()
    expected = _generate_algebraic_contact(geometry, _Target.PENUMBRA)
    object.__setattr__(geometry, "A_inverse", _diagonal(F(0), F(0), F(0)))
    assert _generate_algebraic_contact(geometry, _Target.PENUMBRA) == expected


def test_quadratic_extension_and_split_roots_have_exact_fields() -> None:
    field = _NumberField(_real_roots((-2, 0, 1))[1])
    split = _quadratic_real_roots(field.element(1), field.element(0), field.element(-1))
    assert [root_field.degree for root_field, _, _ in split] == [2, 2]
    assert [t.coefficients() for _, _, t in split] == [(F(-1),), (F(1),)]

    extension = _quadratic_real_roots(
        field.element(1), field.element(0), field.element(-3)
    )
    assert [root_field.degree for root_field, _, _ in extension] == [4, 4]
    assert [root_field.root.P for root_field, _, _ in extension] == [
        (1, 0, -10, 0, 1),
        (1, 0, -10, 0, 1),
    ]
    assert extension[0][2].sign() < 0 and extension[1][2].sign() > 0


def test_two_positive_conjugate_roots_are_ordered_exactly() -> None:
    roots = _real_roots((1, -3, 1))
    assert len(roots) == 2
    assert roots[0].sign(roots[0].poly) == 0
    assert _compare_distinct_roots(roots[0], roots[1]) < 0
    assert roots[0].canonical_interval() != roots[1].canonical_interval()


def test_consistent_singular_line_with_negative_discriminant_is_exhausted() -> None:
    zero = F(0)
    geometry = _ValidatedGeometry(
        (F(-10, 3), F(10), zero),
        (F(5, 3), F(10), zero),
        F(2),
        F(1),
        _diagonal(F(1), F(1), F(1)),
        F(1),
    )
    result = _generate_algebraic_contact(geometry, _Target.PENUMBRA)
    assert result.status is _ContactGenerationStatus.NO_CONTACT_EXHAUSTED
    singular = result.trace[1]
    assert singular.summary.reason is _SummaryReason.NO_REAL_QUADRATIC_POINT
    assert singular.summary.root_count == 1
    assert singular.summary.point_count == 0
    assert singular.summary.reason_value is not None
    assert singular.attempts == ()


def test_zero_discriminant_singular_line_yields_one_checked_point() -> None:
    zero = F(0)
    geometry = _ValidatedGeometry(
        (F(-55, 12), F(25, 16), zero),
        (F(5, 12), F(25, 16), zero),
        F(11, 4),
        F(1, 4),
        _diagonal(F(1), F(1), F(1)),
        F(1),
    )
    result = _generate_algebraic_contact(geometry, _Target.PENUMBRA)
    assert result.status is _ContactGenerationStatus.NO_CONTACT_EXHAUSTED
    singular = result.trace[1]
    assert singular.summary.point_count == 1
    assert len(singular.attempts) == 1
    assert singular.attempts[0].point == ((zero,), (F(1),), (zero,))
    assert singular.attempts[0].disposition is _CandidateDisposition.FAILED_F
    assert singular.attempts[0].predicate is _CandidatePredicate.F_NONZERO


def test_wrong_physical_gate_and_inconsistent_singular_line_are_recorded() -> None:
    zero = F(0)
    geometry = _ValidatedGeometry(
        (F(-3), F(1), zero),
        (F(2), F(1), zero),
        F(2),
        F(1),
        _diagonal(F(1), F(1), F(1)),
        F(1),
    )
    result = _generate_algebraic_contact(geometry, _Target.PENUMBRA)
    assert result.status is _ContactGenerationStatus.NO_CONTACT_EXHAUSTED
    smooth = result.trace[0].attempts
    assert len(smooth) == 1
    assert smooth[0].point == ((F(3, 5),), (F(4, 5),), (zero,))
    assert smooth[0].disposition is _CandidateDisposition.FAILED_GATE
    assert smooth[0].predicate is _CandidatePredicate.GATE_NEGATIVE
    assert smooth[0].predicate_value == (F(-4),)
    singular = result.trace[1].summary
    assert singular.reason is _SummaryReason.INCONSISTENT_LINEAR_SYSTEM
    assert singular.reason_value is not None


def test_exact_perturbations_separate_strict_reach_and_miss() -> None:
    zero = F(0)
    A = _diagonal(F(1), F(1), F(1))
    reach_geometry = _ValidatedGeometry(
        (F(-15), F(10) - F(1, 1000), zero),
        (F(-10), F(10) - F(1, 1000), zero),
        F(2),
        F(1),
        A,
        F(1),
    )
    miss_geometry = _ValidatedGeometry(
        (F(-15), F(10) + F(1, 1000), zero),
        (F(-10), F(10) + F(1, 1000), zero),
        F(2),
        F(1),
        A,
        F(1),
    )
    assert reach_geometry.strict_reach(
        _Target.PENUMBRA, (F(299, 500), F(4003, 5000), zero)
    )
    assert miss_geometry.strict_miss(
        _Target.PENUMBRA, _DualPayload(F(2, 25), (zero, zero, F(8, 25)), zero)
    )
    for geometry in (reach_geometry, miss_geometry):
        result = _generate_algebraic_contact(geometry, _Target.PENUMBRA)
        assert result.status is _ContactGenerationStatus.NO_CONTACT_EXHAUSTED
        assert result.payload is None


def test_no_contact_result_is_not_a_miss_certificate() -> None:
    zero = F(0)
    geometry = _ValidatedGeometry(
        (F(-15), F(20), zero),
        (F(-10), F(20), zero),
        F(2),
        F(1),
        _diagonal(F(1), F(1), F(1)),
        F(1),
    )
    result = _generate_algebraic_contact(geometry, _Target.PENUMBRA)
    assert result.status is _ContactGenerationStatus.NO_CONTACT_EXHAUSTED
    assert result.payload is None
    assert all(item.summary.completed for item in result.trace)
    assert result.trace[0].summary.root_count == 0
    assert result.trace[2].summary.reason is _SummaryReason.PENUMBRAL_APEX_EXCLUDED


@pytest.mark.parametrize("target", [_Target.UMBRA, _Target.ANTUMBRA])
def test_core_apex_off_receiver_is_recorded_exactly(target: _Target) -> None:
    zero = F(0)
    geometry = _ValidatedGeometry(
        (F(-12), zero, zero),
        (F(-7), zero, zero),
        F(2),
        F(1),
        _diagonal(F(1), F(1), F(1)),
        F(1),
    )
    result = _generate_algebraic_contact(geometry, target)
    assert result.status is _ContactGenerationStatus.NO_CONTACT_EXHAUSTED
    apex = result.trace[2]
    assert apex.summary.reason is _SummaryReason.APEX_OFF_RECEIVER
    assert apex.summary.reason_value == (F(3),)
    assert apex.summary.root_count == 0
    assert apex.summary.point_count == 0
    assert apex.attempts == ()


def test_malformed_base_geometry_is_revalidated_before_generation() -> None:
    geometry = _smooth()
    object.__setattr__(geometry, "r", F(-1))
    with pytest.raises(ValueError, match="nonnegative"):
        _generate_algebraic_contact(geometry, _Target.PENUMBRA)
    with pytest.raises(TypeError, match="target"):
        _generate_algebraic_contact(geometry, "penumbra")  # type: ignore[arg-type]
