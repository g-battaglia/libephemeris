# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Synthetic exact-oracle tests for the disconnected strict proposer."""

from __future__ import annotations

from dataclasses import replace
from decimal import (
    ROUND_UP,
    Decimal as D,
    Inexact,
    InvalidOperation,
    getcontext,
    localcontext,
)
from fractions import Fraction as F
from typing import Any

import pytest

from libephemeris import _finite_source_strict_proposer as proposer
from libephemeris._finite_source_certificates import (
    _DomainOutcome,
    _DualPayload,
    _OutsideDomain,
    _Target,
    _ValidatedGeometry,
)
from libephemeris._finite_source_contact_generator import (
    _ContactGenerationStatus,
    _generate_algebraic_contact,
)

Z = F(0)
ONE = F(1)
IDENTITY = ((ONE, Z, Z), (Z, ONE, Z), (Z, Z, ONE))
POLICY = proposer._StrictProposalPolicy(16, 60)


def _geometry(
    source: tuple[F, F, F],
    occulter: tuple[F, F, F],
    r: F = F(2),
    m: F = ONE,
    metric: tuple[tuple[F, F, F], tuple[F, F, F], tuple[F, F, F]] = IDENTITY,
    radius: F = ONE,
) -> _ValidatedGeometry:
    """Construct only exact synthetic source/occulter configurations."""
    return _ValidatedGeometry(source, occulter, r, m, metric, radius)


def _near(offset: F) -> _ValidatedGeometry:
    """Move the reviewed penumbral tangency transversely by exact offset."""
    return _geometry((F(-15), F(10) + offset, Z), (F(-10), F(10) + offset, Z))


def _axis() -> _ValidatedGeometry:
    """Put the source and occulter on the negative first coordinate axis."""
    return _geometry((F(-15), Z, Z), (F(-10), Z, Z))


def _attempt_pattern(result: proposer._StrictProposalResult) -> None:
    """Require complete atomic checkpoints in canonical kind order."""
    kinds = tuple(proposer._StrictProposalKind)
    assert len(result.attempts) % 4 == 0
    for start in range(0, len(result.attempts), 4):
        group = result.attempts[start : start + 4]
        assert tuple(item.kind for item in group) == kinds
        assert len({item.iteration for item in group}) == 1
        assert all(type(item.accepted) is bool for item in group)
        assert all(
            isinstance(item.payload, tuple)
            if item.kind
            in (
                proposer._StrictProposalKind.PRIMAL_CURRENT,
                proposer._StrictProposalKind.PRIMAL_AVERAGE,
            )
            else type(item.payload) is _DualPayload
            for item in group
        )


def test_near_contact_strict_reach_and_miss_have_exact_oracles() -> None:
    """Independent rational witnesses identify each side of tangency."""
    reach_geometry = _near(-F(1, 1000))
    miss_geometry = _near(F(1, 1000))
    assert reach_geometry.strict_reach(
        _Target.PENUMBRA, (F(299, 500), F(4003, 5000), Z)
    )
    assert miss_geometry.strict_miss(
        _Target.PENUMBRA, _DualPayload(F(2, 25), (Z, Z, F(8, 25)), Z)
    )
    reach = proposer._propose_strict_certificate(
        reach_geometry, _Target.PENUMBRA, POLICY
    )
    miss = proposer._propose_strict_certificate(miss_geometry, _Target.PENUMBRA, POLICY)
    assert reach.status is proposer._StrictProposalStatus.REACH
    assert miss.status is proposer._StrictProposalStatus.MISS
    assert reach.reason is miss.reason is proposer._StrictProposalReason.ACCEPTED
    assert isinstance(reach.payload, tuple)
    assert type(miss.payload) is _DualPayload
    assert reach.payload[0] != 0 and reach.payload[1] != 0
    assert miss.payload.beta[2] != 0
    assert reach_geometry.strict_reach(_Target.PENUMBRA, reach.payload)
    assert miss_geometry.strict_miss(_Target.PENUMBRA, miss.payload)
    _attempt_pattern(reach)
    _attempt_pattern(miss)
    assert reach == proposer._propose_strict_certificate(
        reach_geometry, _Target.PENUMBRA, POLICY
    )
    assert miss == proposer._propose_strict_certificate(
        miss_geometry, _Target.PENUMBRA, POLICY
    )


def test_gate_only_miss_reaches_zero_tau_endpoint() -> None:
    """Source and occulter downstream give an exact gate separator."""
    geometry = _geometry((F(12), Z, Z), (F(17), Z, Z))
    assert geometry.strict_miss(_Target.PENUMBRA, _DualPayload(Z, (Z, Z, Z), F(1, 3)))
    result = proposer._propose_strict_certificate(geometry, _Target.PENUMBRA, POLICY)
    assert result.status is proposer._StrictProposalStatus.MISS
    assert type(result.payload) is _DualPayload
    assert result.payload.gate_multiplier > 0
    assert geometry.strict_miss(_Target.PENUMBRA, result.payload)
    exact = proposer._exact_fields(geometry, _Target.PENUMBRA)
    with localcontext(proposer._context(30)):
        fields = proposer._decimal_fields(geometry, exact, POLICY.iterations)
        endpoint = proposer._dual_payload((D(0), D(0), D(0), D(0)), exact, fields)
    assert endpoint.tau == 0
    assert endpoint.gate_multiplier == ONE / exact.G0
    assert geometry.strict_miss(_Target.PENUMBRA, endpoint)
    _attempt_pattern(result)


def test_two_core_nappes_remain_separate() -> None:
    """The same axis geometry has an umbral miss and antumbral reach."""
    geometry = _axis()
    assert geometry.strict_miss(_Target.UMBRA, _DualPayload(F(1, 25), (Z, Z, Z), Z))
    assert geometry.strict_reach(_Target.ANTUMBRA, (Z, Z, Z))
    umbra = proposer._propose_strict_certificate(geometry, _Target.UMBRA, POLICY)
    antumbra = proposer._propose_strict_certificate(geometry, _Target.ANTUMBRA, POLICY)
    assert umbra.status is proposer._StrictProposalStatus.MISS
    assert antumbra.status is proposer._StrictProposalStatus.REACH
    assert type(umbra.payload) is _DualPayload
    assert isinstance(antumbra.payload, tuple)
    assert geometry.strict_miss(_Target.UMBRA, umbra.payload)
    assert geometry.strict_reach(_Target.ANTUMBRA, antumbra.payload)
    assert umbra.payload.gate_multiplier >= 0
    _attempt_pattern(umbra)
    _attempt_pattern(antumbra)


def test_coupled_metric_reversed_orientation_stays_exact() -> None:
    """A non-diagonal metric needs the correctly oriented transverse beta."""
    metric = ((F(2), ONE, Z), (ONE, F(2), Z), (Z, Z, ONE))
    geometry = _geometry((F(15), F(10), Z), (F(10), F(10), Z), metric=metric)
    assert geometry.strict_miss(
        _Target.PENUMBRA, _DualPayload(F(2, 25), (Z, Z, -F(8, 25)), Z)
    )
    assert not geometry.strict_miss(
        _Target.PENUMBRA, _DualPayload(F(2, 25), (Z, Z, F(8, 25)), Z)
    )
    assert not geometry.strict_miss(
        _Target.PENUMBRA, _DualPayload(F(2, 25), (Z, Z, Z), Z)
    )
    result = proposer._propose_strict_certificate(geometry, _Target.PENUMBRA, POLICY)
    assert result.status is proposer._StrictProposalStatus.MISS
    assert type(result.payload) is _DualPayload
    assert result.payload.beta[2] < 0
    assert geometry.strict_miss(_Target.PENUMBRA, result.payload)


def test_exact_contact_only_unresolved_with_finite_trace() -> None:
    """An exact tangent cannot acquire a strict certificate by rounding."""
    geometry = _near(Z)
    assert (
        _generate_algebraic_contact(geometry, _Target.PENUMBRA).status
        is _ContactGenerationStatus.CONTACT
    )
    tangent = (F(3, 5), F(4, 5), Z)
    assert geometry.contains(_Target.PENUMBRA, tangent)
    assert not geometry.strict_reach(_Target.PENUMBRA, tangent)
    policy = proposer._StrictProposalPolicy(9, 60)
    result = proposer._propose_strict_certificate(geometry, _Target.PENUMBRA, policy)
    assert result.status is proposer._StrictProposalStatus.UNRESOLVED
    assert result.payload is None
    assert result.reason is proposer._StrictProposalReason.ITERATION_LIMIT
    assert result.iterations_completed == 9
    assert [result.attempts[i].iteration for i in range(0, 20, 4)] == [1, 2, 4, 8, 9]
    assert not any(item.accepted for item in result.attempts)
    _attempt_pattern(result)


def _dummy_fields(gated: bool) -> proposer._DecimalFields:
    """Supply only the fields read by the synthetic projection probe."""
    affine = (D(0), D(0), D(0), D(0))
    row = (D(0), D(0), D(0))
    return proposer._DecimalFields(
        affine,
        affine if gated else None,
        (affine, affine, affine),
        (row, row, row),
        D(1),
        D(1),
        D(1),
    )


def test_projection_zero_unit_and_exterior_norm_branches() -> None:
    """Gated and antumbral projections take all three norm branches."""
    gated = _dummy_fields(True)
    antumbral = _dummy_fields(False)
    assert proposer._project(
        (D("0.5"), D(0), D(0), D(0)), (D(0), D(0), D(0), D(0)), gated
    ) == (D("0.5"), D(0), D(0), D(0))
    assert proposer._project(
        (D("0.5"), D(0), D(0), D(0)), (D(1), D(0), D(0), D(0)), gated
    ) == (D(0), D(0), D(0), D(0))
    assert proposer._project(
        (D("0.5"), D(0), D(0), D(0)), (D(0), D(-1), D(0), D(0)), gated
    ) == (D("0.75"), D("0.75"), D(0), D(0))
    assert proposer._project(
        (D("0.5"), D(0), D(0), D(0)), (D(0), D(-2), D(0), D(0)), gated
    ) == (D(1), D(1), D(0), D(0))
    assert proposer._project((D(0), D(0), D(0)), (D(0), D(0), D(0)), antumbral) == (
        D(0),
        D(0),
        D(0),
    )
    assert proposer._project((D(0), D(0), D(0)), (D(-1), D(0), D(0)), antumbral) == (
        D(1),
        D(0),
        D(0),
    )
    assert proposer._project((D(0), D(0), D(0)), (D(-2), D(0), D(0)), antumbral) == (
        D(1),
        D(0),
        D(0),
    )


def test_exact_primal_inward_and_dual_cone_contractions() -> None:
    """Canonical rational repairs are only proposals, even on boundaries."""
    geometry = _axis()
    point = proposer._primal_payload((D(1), D(0), D(0)), geometry, 4)
    assert point == (F(9999, 10000), Z, Z)
    assert geometry.strict_reach(_Target.ANTUMBRA, point)
    outside = proposer._primal_payload((D(2), D(0), D(0)), geometry, 4)
    assert outside == (F(1, 2), Z, Z)
    assert geometry.strict_reach(_Target.ANTUMBRA, outside)
    exact = proposer._exact_fields(geometry, _Target.PENUMBRA)
    with localcontext(proposer._context(30)):
        fields = proposer._decimal_fields(geometry, exact, 2)
        dual = proposer._dual_payload((D(1), D(2), D(0), D(0)), exact, fields)
    norm_squared = sum((component * component for component in dual.beta), Z)
    assert norm_squared <= dual.tau * dual.tau * exact.H
    assert norm_squared > 0
    assert type(dual) is _DualPayload


def test_caller_decimal_context_cannot_change_trace() -> None:
    """Ambient precision, rounding and Inexact trap do not enter proposals."""
    geometry = _near(-F(1, 1000))
    expected = proposer._propose_strict_certificate(geometry, _Target.PENUMBRA, POLICY)
    with localcontext() as ambient:
        ambient.prec = 2
        ambient.rounding = ROUND_UP
        ambient.traps[Inexact] = True
        before = (ambient.prec, ambient.rounding, ambient.traps.copy())
        actual = proposer._propose_strict_certificate(
            geometry, _Target.PENUMBRA, POLICY
        )
        assert (ambient.prec, ambient.rounding, ambient.traps.copy()) == before
    assert actual == expected
    assert getcontext().prec != 2


def test_policy_and_target_order_precede_geometry_fields() -> None:
    """Each documented input failure occurs before a forged geometry read."""
    unreadable = object.__new__(_ValidatedGeometry)
    with pytest.raises(TypeError, match="target"):
        proposer._propose_strict_certificate(unreadable, "penumbra", POLICY)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="policy"):
        proposer._propose_strict_certificate(unreadable, _Target.PENUMBRA, object())  # type: ignore[arg-type]
    for policy, match, error in (
        (proposer._StrictProposalPolicy(True, 40), "iterations", TypeError),
        (proposer._StrictProposalPolicy(1, False), "decimal_digits", TypeError),
        (proposer._StrictProposalPolicy(0, 1), "iterations", ValueError),
        (proposer._StrictProposalPolicy(1, 1), "decimal_digits", ValueError),
    ):
        with pytest.raises(error, match=match):
            proposer._propose_strict_certificate(unreadable, _Target.PENUMBRA, policy)
    with pytest.raises(TypeError, match="geometry"):
        proposer._propose_strict_certificate(object(), _Target.PENUMBRA, POLICY)  # type: ignore[arg-type]


def test_fresh_geometry_rejects_malformed_and_preserves_domain_outcomes() -> None:
    """Forged derived fields cannot bypass exact reconstruction and gate order."""
    geometry = _near(Z)
    expected = proposer._propose_strict_certificate(geometry, _Target.PENUMBRA, POLICY)
    object.__setattr__(geometry, "A_inverse", ((Z, Z, Z),) * 3)
    object.__setattr__(geometry, "u", (Z, Z, Z))
    assert (
        proposer._propose_strict_certificate(geometry, _Target.PENUMBRA, POLICY)
        == expected
    )
    object.__setattr__(geometry, "r", F(-1))
    with pytest.raises(ValueError, match="nonnegative"):
        proposer._propose_strict_certificate(geometry, _Target.PENUMBRA, POLICY)
    object.__setattr__(geometry, "r", F(2))
    object.__setattr__(geometry, "p", (F(-10), 10, Z))
    with pytest.raises(TypeError, match="Fraction"):
        proposer._propose_strict_certificate(geometry, _Target.PENUMBRA, POLICY)
    for field, value, outcome in (
        ("m", Z, _DomainOutcome.UNSUPPORTED_SLOPE),
        ("p", (F(-14), F(10), Z), _DomainOutcome.EXCLUDED_DOMAIN),
        ("B", (Z, Z, Z), _DomainOutcome.UNCERTIFIED_DOMAIN),
    ):
        candidate = _near(Z)
        object.__setattr__(candidate, field, value)
        with pytest.raises(_OutsideDomain) as caught:
            proposer._propose_strict_certificate(candidate, _Target.PENUMBRA, POLICY)
        assert caught.value.outcome is outcome


def test_decimal_faults_preserve_only_committed_checkpoints(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A fourth provisional fault discards even an earlier accepted point."""
    geometry = _axis()
    original = proposer._dual_payload
    calls = 0

    def fail_second_dual(*args: Any, **kwargs: Any) -> _DualPayload:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise InvalidOperation
        return original(*args, **kwargs)

    monkeypatch.setattr(proposer, "_dual_payload", fail_second_dual)
    result = proposer._propose_strict_certificate(geometry, _Target.ANTUMBRA, POLICY)
    assert result.status is proposer._StrictProposalStatus.UNRESOLVED
    assert result.reason is proposer._StrictProposalReason.DECIMAL_FAULT
    assert result.iterations_completed == 0
    assert result.attempts == ()
    assert calls == 2


def test_projection_fault_retains_completed_checkpoint(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The next dual projection belongs outside the prior iteration."""
    monkeypatch.setattr(
        proposer, "_project", lambda *_args: (_ for _ in ()).throw(InvalidOperation())
    )
    result = proposer._propose_strict_certificate(
        _near(Z), _Target.PENUMBRA, proposer._StrictProposalPolicy(3, 50)
    )
    assert result.status is proposer._StrictProposalStatus.UNRESOLVED
    assert result.reason is proposer._StrictProposalReason.DECIMAL_FAULT
    assert result.iterations_completed == 1
    assert len(result.attempts) == 4
    _attempt_pattern(result)


@pytest.mark.parametrize("negative", [False, True])
def test_invalid_decimal_support_square_is_fault(
    negative: bool,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Numerically impossible positive-definite support is not a miss."""
    original = proposer._decimal_fields

    def bad_inverse(*args: Any, **kwargs: Any) -> proposer._DecimalFields:
        fields = original(*args, **kwargs)
        diagonal = D(-1) if negative else D(0)
        inverse = (
            (diagonal, D(0), D(0)),
            (D(0), diagonal, D(0)),
            (D(0), D(0), diagonal),
        )
        return replace(fields, inverse=inverse)

    monkeypatch.setattr(proposer, "_decimal_fields", bad_inverse)
    result = proposer._propose_strict_certificate(_near(Z), _Target.PENUMBRA, POLICY)
    assert result.status is proposer._StrictProposalStatus.UNRESOLVED
    assert result.reason is proposer._StrictProposalReason.DECIMAL_FAULT
    assert result.iterations_completed == 0
    assert result.attempts == ()


def test_setup_decimal_fault_is_typed_unresolved(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A trapped conversion before iteration leaves an empty trace."""
    monkeypatch.setattr(
        proposer,
        "_decimal_fraction",
        lambda *_: (_ for _ in ()).throw(InvalidOperation()),
    )
    result = proposer._propose_strict_certificate(_near(Z), _Target.PENUMBRA, POLICY)
    assert result.status is proposer._StrictProposalStatus.UNRESOLVED
    assert result.reason is proposer._StrictProposalReason.DECIMAL_FAULT
    assert result.iterations_completed == 0 and result.attempts == ()


def test_verifier_rejection_stays_unresolved_after_final_checkpoint(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """No Decimal estimate is promoted when all exact checks reject."""
    monkeypatch.setattr(_ValidatedGeometry, "strict_reach", lambda *_: False)
    monkeypatch.setattr(_ValidatedGeometry, "strict_miss", lambda *_: False)
    result = proposer._propose_strict_certificate(
        _near(-F(1, 1000)),
        _Target.PENUMBRA,
        proposer._StrictProposalPolicy(3, 50),
    )
    assert result.status is proposer._StrictProposalStatus.UNRESOLVED
    assert result.reason is proposer._StrictProposalReason.ITERATION_LIMIT
    assert result.iterations_completed == 3
    assert [item.iteration for item in result.attempts[::4]] == [1, 2, 3]
    assert all(not item.accepted for item in result.attempts)
    _attempt_pattern(result)


def test_verifier_exception_or_opposite_acceptance_is_invariant(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A broken exact checker cannot return a physical result or fault trace."""
    geometry = _axis()

    def fail(*_args: Any) -> bool:
        raise ValueError("injected verifier error")

    monkeypatch.setattr(_ValidatedGeometry, "strict_reach", fail)
    with pytest.raises(proposer._StrictProposalInvariantError, match="verifier"):
        proposer._propose_strict_certificate(geometry, _Target.ANTUMBRA, POLICY)
    monkeypatch.setattr(_ValidatedGeometry, "strict_reach", lambda *_: True)
    monkeypatch.setattr(_ValidatedGeometry, "strict_miss", lambda *_: True)
    with pytest.raises(proposer._StrictProposalInvariantError, match="both strict"):
        proposer._propose_strict_certificate(geometry, _Target.ANTUMBRA, POLICY)
