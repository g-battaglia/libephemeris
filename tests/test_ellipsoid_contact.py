# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Tests for the private certified ellipsoid contact frame."""

from __future__ import annotations

import dataclasses
import math
from fractions import Fraction

import pytest

from flint import arb, ctx

from libephemeris.ellipsoid_contact import (
    _ConeSection,
    _CentralAxisStatus,
    _CertifiedContactProofs,
    _ContactCapabilityBlock,
    _ContactCapabilityStatus,
    _ContactEqualityProof,
    _ContactSourceTag,
    _CoreShadowClass,
    _ContactEqualityVariant,
    _CrossDualWitness,
    _EllipsoidContactFrame,
    _GlobalReachStatus,
    _ProjectedDualWitness,
    _RegularRootCertificate,
    _WitnessGenerationStatus,
    _ZeroAngleLineRelation,
    _build_contact_inputs,
    _certify_contact_inputs,
    _certify_regular_root_box,
    _core_cone_section,
    _evaluate_central_axis_ray,
    _evaluate_cone_apex,
    _evaluate_nappe_boundary,
    _evaluate_radial_subgradient,
    _evaluate_global_reach,
    _evaluate_regular_contact,
    _generate_global_witness,
    _nappe_boundary_jacobian,
    _pow2_upper_sqrt,
    _primal_level,
    _dual_level,
    _reduce_regular_residuals,
    _regular_contact_jacobian,
)
from libephemeris.intervals import IntervalCertificationError


def _frame(**overrides) -> _EllipsoidContactFrame:
    values = {
        "metric_km_minus_2": (
            (0.25, 0.0, 0.0),
            (0.0, 0.25, 0.0),
            (0.0, 0.0, 0.25),
        ),
        "axis_point_km": (-3.5, 2.0, 6.0),
        "axis_span": (0.0, 3.0, 4.0),
        "penumbra_branch_sign": 1,
        "core_branch_sign": -1,
    }
    values.update(overrides)
    return _EllipsoidContactFrame(**values)


def _dual(
    span=(Fraction(0), Fraction(0), Fraction(1)),
    seed=(Fraction(1), Fraction(0), Fraction(0)),
    scale=Fraction(1),
    multiplier=Fraction(0),
) -> _ProjectedDualWitness:
    return _ProjectedDualWitness(span, seed, scale, multiplier)


@pytest.mark.parametrize(
    ("axis_point", "axis_span", "status"),
    [
        ((0.0, 0.0, -3.0), (0.0, 0.0, 1.0), _CentralAxisStatus.CROSSING),
        ((2.0, 0.0, 0.0), (0.0, 0.0, 1.0), _CentralAxisStatus.TANGENT),
        ((3.0, 0.0, 0.0), (0.0, 0.0, 1.0), _CentralAxisStatus.MISS),
        ((0.0, 0.0, 3.0), (0.0, 0.0, 1.0), _CentralAxisStatus.MISS),
    ],
)
def test_central_axis_ray_respects_downstream_orientation(
    axis_point, axis_span, status
) -> None:
    """The oriented ray distinguishes forward roots, tangency, and behind roots."""
    proof = _evaluate_central_axis_ray(
        _frame(axis_point_km=axis_point, axis_span=axis_span)
    )
    assert proof.status is status


def _source_tag(**overrides) -> _ContactSourceTag:
    values = {
        "provider": "calc_integer",
        "request_flags": 6146,
        "frame": "true_equator_of_date",
        "body_kind": "integer",
        "tjd_ut": 2451545.0,
        "retflag": 6146,
    }
    values.update(overrides)
    return _ContactSourceTag(**values)


def test_capability_carrier_builds_ready_positive_slope_inputs() -> None:
    """Valid positive-slope geometry produces the exact private ready variant."""
    result = _build_contact_inputs(
        _frame(axis_point_km=(0.0, 0.0, -3.0), axis_span=(0.0, 0.0, 1.0)),
        4.0,
        0.9,
        -2.0,
        0.95,
        10.0,
        1.0,
        _source_tag(),
    )
    assert not isinstance(result, _ContactCapabilityBlock)
    assert result.core.shadow_class is _CoreShadowClass.UMBRA
    assert result.central_axis.status is _CentralAxisStatus.CROSSING
    assert tuple(attempt.requested_bits for attempt in result.precision_evidence) == (
        160,
        256,
    )
    assert all(
        attempt.central_axis.status is _CentralAxisStatus.CROSSING
        for attempt in result.precision_evidence
    )


def test_private_carrier_connects_to_both_witness_generators() -> None:
    """Ready producer inputs produce owned penumbral and core proofs privately."""
    inputs = _build_contact_inputs(
        _frame(axis_point_km=(0.0, 0.0, -3.0), axis_span=(0.0, 0.0, 1.0)),
        4.0,
        1.0,
        -2.0,
        1.0,
        10.0,
        1.0,
        _source_tag(),
    )
    assert not isinstance(inputs, _ContactCapabilityBlock)
    result = _certify_contact_inputs(inputs, max_level=1)
    assert isinstance(result, _CertifiedContactProofs)
    assert result.penumbra.source_identity[0] == result.core.source_identity[0]
    assert result.penumbra.status in {
        _WitnessGenerationStatus.REACH,
        _WitnessGenerationStatus.MISS,
        _WitnessGenerationStatus.CONTACT,
    }
    assert result.core.status in {
        _WitnessGenerationStatus.REACH,
        _WitnessGenerationStatus.MISS,
        _WitnessGenerationStatus.CONTACT,
    }


def test_private_certification_rejects_malformed_decisive_result(monkeypatch) -> None:
    """A decisive top-level status cannot bypass proof/evidence ownership."""
    inputs = _build_contact_inputs(
        _frame(axis_point_km=(0.0, 0.0, -3.0), axis_span=(0.0, 0.0, 1.0)),
        4.0,
        1.0,
        -2.0,
        1.0,
        10.0,
        1.0,
        _source_tag(),
    )
    assert not isinstance(inputs, _ContactCapabilityBlock)
    real = _generate_global_witness(inputs.frame, inputs.penumbra, max_level=1)
    malformed = dataclasses.replace(
        real,
        proof=None,
        source_identity=("wrong",),
        accepted_candidate=None,
    )
    monkeypatch.setattr(
        "libephemeris.ellipsoid_contact._generate_global_witness",
        lambda *_args, **_kwargs: malformed,
    )
    result = _certify_contact_inputs(inputs, max_level=1)
    assert not isinstance(result, _CertifiedContactProofs)
    assert result.status is _ContactCapabilityStatus.INVALID_CERTIFICATION


@pytest.mark.parametrize(
    "mutation",
    [
        lambda real: dataclasses.replace(
            real,
            accepted_candidate=dataclasses.replace(
                real.accepted_candidate, kind="bogus"
            ),
        ),
        lambda real: dataclasses.replace(
            real,
            accepted_candidate=dataclasses.replace(
                real.accepted_candidate, payload=("wrong",)
            ),
        ),
        lambda real: dataclasses.replace(
            real,
            work=dataclasses.replace(
                real.work,
                primal_candidates=0,
                dual_candidates=0,
                total_candidates=0,
            ),
        ),
        lambda real: dataclasses.replace(
            real,
            precision_evidence=tuple(
                dataclasses.replace(
                    record,
                    proof=dataclasses.replace(
                        record.proof,
                        reason="wrong precision proof",
                    ),
                )
                for record in real.precision_evidence
            ),
        ),
        lambda real: dataclasses.replace(real, generation_level=99),
        lambda real: dataclasses.replace(real, candidate_index=99),
        lambda real: dataclasses.replace(
            real,
            work=dataclasses.replace(real.work, last_indices=(99, 99, 99)),
        ),
        lambda real: dataclasses.replace(
            real,
            accepted_candidate=dataclasses.replace(real.accepted_candidate, level=99),
        ),
        lambda real: dataclasses.replace(
            real,
            accepted_candidate=dataclasses.replace(
                real.accepted_candidate, indices=(99, 99, 99)
            ),
            work=dataclasses.replace(real.work, last_indices=(99, 99, 99)),
        ),
        lambda real: dataclasses.replace(
            real,
            precision_evidence=tuple(
                dataclasses.replace(
                    record,
                    proof=dataclasses.replace(
                        record.proof,
                        upper_bound=arb(-999),
                    ),
                )
                for record in real.precision_evidence
            ),
        ),
        lambda real: dataclasses.replace(
            real,
            proof=dataclasses.replace(real.proof, upper_bound=arb(-999)),
            precision_evidence=tuple(
                dataclasses.replace(
                    record,
                    proof=dataclasses.replace(
                        record.proof,
                        upper_bound=arb(-999),
                    ),
                )
                for record in real.precision_evidence
            ),
            accepted_candidate=dataclasses.replace(
                real.accepted_candidate,
                precision_evidence=tuple(
                    dataclasses.replace(
                        record,
                        proof=dataclasses.replace(
                            record.proof,
                            upper_bound=arb(-999),
                        ),
                    )
                    for record in real.precision_evidence
                ),
            ),
        ),
    ],
)
def test_private_certification_rejects_forged_candidate_ledger(
    monkeypatch, mutation
) -> None:
    """Kind, payload, counters, and precision records are all authenticated."""
    inputs = _build_contact_inputs(
        _frame(axis_point_km=(0.0, 0.0, -3.0), axis_span=(0.0, 0.0, 1.0)),
        4.0,
        1.0,
        -2.0,
        1.0,
        10.0,
        1.0,
        _source_tag(),
    )
    assert not isinstance(inputs, _ContactCapabilityBlock)
    real = _generate_global_witness(inputs.frame, inputs.penumbra, max_level=1)
    assert real.accepted_candidate is not None
    forged = mutation(real)
    monkeypatch.setattr(
        "libephemeris.ellipsoid_contact._generate_global_witness",
        lambda *_args, **_kwargs: forged,
    )
    result = _certify_contact_inputs(inputs, max_level=1)
    assert not isinstance(result, _CertifiedContactProofs)
    assert result.status is _ContactCapabilityStatus.INVALID_CERTIFICATION


@pytest.mark.parametrize(
    "mutation",
    [
        lambda attempt: dataclasses.replace(
            attempt, level=7, indices=(99, 99, 99), payload=("forged",)
        ),
        lambda attempt: dataclasses.replace(
            attempt,
            precision_evidence=tuple(
                dataclasses.replace(
                    record,
                    proof=dataclasses.replace(
                        record.proof, reason="forged earlier proof"
                    ),
                )
                for record in attempt.precision_evidence
            ),
        ),
    ],
)
def test_private_certification_replays_prior_attempt_ledger(
    monkeypatch, mutation
) -> None:
    """Every stored pre-final attempt must match canonical order and replay."""
    inputs = _build_contact_inputs(
        _frame(
            metric_km_minus_2=(
                (1.0 / 9.0, 0.0, 0.0),
                (0.0, 1.0 / 9.0, 0.0),
                (0.0, 0.0, 0.25),
            ),
            axis_point_km=(-5.0, 0.0, 0.0),
            axis_span=(0.0, 1.0, 1.0),
        ),
        2.0,
        1.0,
        2.0,
        1.0,
        10.0,
        1.0,
        _source_tag(),
    )
    assert not isinstance(inputs, _ContactCapabilityBlock)
    real = _generate_global_witness(inputs.frame, inputs.penumbra, max_level=1)
    assert len(real.attempted_candidates) > 1
    ledger = list(real.attempted_candidates)
    ledger[0] = mutation(ledger[0])
    forged = dataclasses.replace(real, attempted_candidates=tuple(ledger))
    monkeypatch.setattr(
        "libephemeris.ellipsoid_contact._generate_global_witness",
        lambda *_args, **_kwargs: forged,
    )
    result = _certify_contact_inputs(inputs, max_level=1)
    assert not isinstance(result, _CertifiedContactProofs)
    assert result.status is _ContactCapabilityStatus.INVALID_CERTIFICATION


@pytest.mark.parametrize("variant_index", [-1, 1, 99])
def test_private_certification_rejects_forged_equality_variant(
    monkeypatch, variant_index
) -> None:
    """The sole source-selected equality variant has canonical index zero."""
    inputs = _build_contact_inputs(
        _frame(axis_point_km=(0.0, 0.0, -3.0), axis_span=(0.0, 0.0, 1.0)),
        0.0,
        1.0,
        0.0,
        1.0,
        10.0,
        1.0,
        _source_tag(),
    )
    assert not isinstance(inputs, _ContactCapabilityBlock)
    real = _generate_global_witness(inputs.frame, inputs.penumbra, max_level=0)
    assert real.status is _WitnessGenerationStatus.CONTACT
    forged = dataclasses.replace(real, equality_variant_index=variant_index)
    monkeypatch.setattr(
        "libephemeris.ellipsoid_contact._generate_global_witness",
        lambda *_args, **_kwargs: forged,
    )
    result = _certify_contact_inputs(inputs, max_level=0)
    assert not isinstance(result, _CertifiedContactProofs)
    assert result.status is _ContactCapabilityStatus.INVALID_CERTIFICATION


def test_private_certification_keeps_internal_invalid_distinct(monkeypatch) -> None:
    """Generator defects are not mislabeled as malformed source data."""
    inputs = _build_contact_inputs(
        _frame(axis_point_km=(0.0, 0.0, -3.0), axis_span=(0.0, 0.0, 1.0)),
        4.0,
        1.0,
        -2.0,
        1.0,
        10.0,
        1.0,
        _source_tag(),
    )
    assert not isinstance(inputs, _ContactCapabilityBlock)
    real = _generate_global_witness(inputs.frame, inputs.penumbra, max_level=1)
    invalid = dataclasses.replace(
        real,
        status=_WitnessGenerationStatus.INVALID,
        reason="internal defect",
    )
    monkeypatch.setattr(
        "libephemeris.ellipsoid_contact._generate_global_witness",
        lambda *_args, **_kwargs: invalid,
    )
    result = _certify_contact_inputs(inputs, max_level=1)
    assert not isinstance(result, _CertifiedContactProofs)
    assert result.status is _ContactCapabilityStatus.INVALID_CERTIFICATION
    assert result.reason == "internal defect"


@pytest.mark.parametrize(
    ("tag", "source_radius", "status"),
    [
        (_source_tag(frame="j2000"), 10.0, _ContactCapabilityStatus.MIXED_OUTPUT_FRAME),
        (_source_tag(), 0.0, _ContactCapabilityStatus.EXCLUDED_CORE_SLOPE),
        (
            _source_tag(provider="fixstar", body_kind="integer"),
            10.0,
            _ContactCapabilityStatus.INVALID_SOURCE,
        ),
    ],
)
def test_capability_carrier_blocks_unsupported_source_domains(
    tag, source_radius, status
) -> None:
    """Mixed frames, excluded slopes, and inconsistent tags fail closed."""
    result = _build_contact_inputs(
        _frame(axis_point_km=(0.0, 0.0, -3.0), axis_span=(0.0, 0.0, 1.0)),
        4.0,
        0.9,
        -2.0,
        0.95,
        source_radius,
        1.0,
        tag,
    )
    assert isinstance(result, _ContactCapabilityBlock)
    assert result.status is status
    assert result.core is None


@pytest.mark.parametrize(
    "tag",
    [
        _source_tag(request_flags=0),
        _source_tag(retflag=0),
        _source_tag(tjd_ut=math.nan),
        _source_tag(request_flags=True),
    ],
)
def test_capability_carrier_authenticates_every_source_tag_field(tag) -> None:
    """Unauthenticated request, retflag, time, or type stays invalid."""
    result = _build_contact_inputs(
        _frame(axis_point_km=(0.0, 0.0, -3.0), axis_span=(0.0, 0.0, 1.0)),
        4.0,
        0.9,
        -2.0,
        0.95,
        10.0,
        1.0,
        tag,
    )
    assert isinstance(result, _ContactCapabilityBlock)
    assert result.status is _ContactCapabilityStatus.INVALID_SOURCE


@pytest.mark.parametrize("ephemeris", [1, 2, 4])
def test_capability_carrier_accepts_all_owned_ephemeris_selectors(ephemeris) -> None:
    """JPL, accepted selector, and analytic-compatible selector authenticate."""
    flags = ephemeris | 2048 | 4096
    result = _build_contact_inputs(
        _frame(axis_point_km=(0.0, 0.0, -3.0), axis_span=(0.0, 0.0, 1.0)),
        4.0,
        0.9,
        -2.0,
        0.95,
        10.0,
        1.0,
        _source_tag(request_flags=flags, retflag=flags),
    )
    assert not isinstance(result, _ContactCapabilityBlock)


def test_capability_carrier_catches_malformed_penumbral_scalar() -> None:
    """Wrong producer scalar types cannot escape before typed classification."""
    result = _build_contact_inputs(
        _frame(axis_point_km=(0.0, 0.0, -3.0), axis_span=(0.0, 0.0, 1.0)),
        "4.0",  # type: ignore[arg-type]
        0.9,
        -2.0,
        0.95,
        10.0,
        1.0,
        _source_tag(),
    )
    assert isinstance(result, _ContactCapabilityBlock)
    assert result.status is _ContactCapabilityStatus.INVALID_SOURCE
    assert result.precision_evidence == ()


def test_capability_carrier_rejects_same_status_root_disagreement(monkeypatch) -> None:
    """Matching classes with contradictory root signs remain unresolved."""
    calls = 0
    original = _evaluate_central_axis_ray

    def disagree(frame):
        nonlocal calls
        calls += 1
        proof = original(frame)
        if calls == 2:
            return dataclasses.replace(
                proof,
                lower_parameter=arb(-9),
                upper_parameter=arb(9),
            )
        return proof

    monkeypatch.setattr(
        "libephemeris.ellipsoid_contact._evaluate_central_axis_ray", disagree
    )
    result = _build_contact_inputs(
        _frame(axis_point_km=(0.0, 0.0, -3.0), axis_span=(0.0, 0.0, 1.0)),
        4.0,
        0.9,
        -2.0,
        0.95,
        10.0,
        1.0,
        _source_tag(),
    )
    assert isinstance(result, _ContactCapabilityBlock)
    assert result.status is _ContactCapabilityStatus.UNRESOLVED_AXIS
    assert len(result.precision_evidence) == 2


def test_capability_carrier_restores_ambient_precision() -> None:
    """Independent starts restore process-local Arb precision."""
    previous = ctx.prec
    try:
        ctx.prec = 93
        result = _build_contact_inputs(
            _frame(axis_point_km=(0.0, 0.0, -3.0), axis_span=(0.0, 0.0, 1.0)),
            4.0,
            0.9,
            -2.0,
            0.95,
            10.0,
            1.0,
            _source_tag(),
        )
        assert ctx.prec == 93
        assert not isinstance(result, _ContactCapabilityBlock)
    finally:
        ctx.prec = previous


def test_power_of_two_bound_is_smallest_exact_enclosure() -> None:
    """Source-rational bounds use one unique power of two."""
    assert _pow2_upper_sqrt(Fraction(9, 4)) == 2
    assert _pow2_upper_sqrt(Fraction(1, 4)) == Fraction(1, 2)
    with pytest.raises(ValueError):
        _pow2_upper_sqrt(Fraction(0))


def test_canonical_primal_levels_own_each_point_once() -> None:
    """Level-one endpoints survive and recursive duplicates start at level two."""
    bounds = (Fraction(1), Fraction(1), Fraction(1))
    levels = [list(_primal_level(bounds, level)) for level in range(4)]
    assert [len(level) for level in levels] == [1, 26, 98, 604]
    seen = set()
    for level in levels:
        points = {point for _indices, point in level}
        assert not seen & points
        seen |= points
    assert (Fraction(-1),) * 3 in seen
    assert (Fraction(1),) * 3 in seen


def test_canonical_dual_levels_own_each_payload_once() -> None:
    """Dual level one keeps endpoints and skips only physical zero."""
    levels = [
        list(
            _dual_level(
                (Fraction(0), Fraction(0), Fraction(1)),
                Fraction(1),
                Fraction(2),
                level,
            )
        )
        for level in range(4)
    ]
    assert [len(level) for level in levels] == [1, 26, 98, 604]
    seen = set()
    for level in levels:
        payloads = {(vector, multiplier) for _indices, vector, multiplier in level}
        assert not seen & payloads
        seen |= payloads
    assert ((Fraction(-1), Fraction(-1), Fraction(0)), Fraction(0)) in seen


def test_bounded_generator_finds_primal_reach_and_dual_miss() -> None:
    """Canonical stream reaches both approved strict proof paths."""
    reach = _generate_global_witness(
        _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
        _ConeSection(0.0, 0.8, 1),
        max_level=1,
    )
    assert reach.status is _WitnessGenerationStatus.REACH
    miss_frame = _frame(
        metric_km_minus_2=(
            (1.0 / 9.0, 0.0, 0.0),
            (0.0, 1.0 / 9.0, 0.0),
            (0.0, 0.0, 0.25),
        ),
        axis_point_km=(-5.0, 0.0, 0.0),
        axis_span=(0.0, 1.0, 1.0),
    )
    miss = _generate_global_witness(
        miss_frame,
        _ConeSection(1.0, 1.0, 1),
        max_level=1,
    )
    assert miss.status is _WitnessGenerationStatus.MISS
    assert miss.multiplier_bound is not None
    assert miss.multiplier_evidence is not None
    assert miss.multiplier_evidence.power_bound == miss.multiplier_bound
    assert miss.accepted_candidate is not None
    assert miss.accepted_candidate.kind == "dual"
    assert miss.attempted_candidates[-1] is miss.accepted_candidate
    assert {record.bits for record in miss.precision_evidence} == {160, 256}
    assert miss.source_identity[-1] == "witness-generator-v1"


def test_bounded_generator_proves_zero_angle_contact() -> None:
    """Source-selected equality precedes the strict candidate stream."""
    result = _generate_global_witness(
        _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
        _ConeSection(0.0, 1.0, 1),
        max_level=0,
    )
    assert result.status is _WitnessGenerationStatus.CONTACT
    assert result.work.total_candidates == 0
    assert {record.bits for record in result.precision_evidence} == {160, 256}


def test_bounded_generator_selects_exact_zero_angle_tangency() -> None:
    """Exact source discriminant selects TANGENT rather than CROSSING."""
    result = _generate_global_witness(
        _frame(axis_point_km=(2.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
        _ConeSection(0.0, 1.0, 1),
        max_level=0,
    )
    assert result.status is _WitnessGenerationStatus.CONTACT
    assert result.equality_variant_index == 0
    assert result.accepted_candidate is not None
    assert result.accepted_candidate.kind == "equality"
    assert result.proof is not None
    assert result.proof.equality_proof is not None
    assert result.proof.equality_proof.line_relation is _ZeroAngleLineRelation.TANGENT


def test_bounded_generator_returns_typed_invalid_source() -> None:
    """Malformed source data never escapes the private result boundary."""
    result = _generate_global_witness(
        _frame(
            metric_km_minus_2=(
                (-1.0, 0.0, 0.0),
                (0.0, -1.0, 0.0),
                (0.0, 0.0, 1.0),
            )
        ),
        _ConeSection(1.0, 1.0, 1),
        max_level=0,
    )
    assert result.status is _WitnessGenerationStatus.INVALID
    assert result.work.total_candidates == 0
    assert "positive definite" in result.reason


def test_bounded_generator_deduplicates_zero_multiplier_payloads() -> None:
    """Exact payload ownership, not index ownership, drives dual counters."""
    frame = _frame(axis_point_km=(5.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    result = _generate_global_witness(
        frame,
        _ConeSection(1.0, 1.0, 1),
        max_level=1,
        max_primal=1_000,
        max_dual=1_000,
        max_total=2_000,
    )
    assert result.work.dual_candidates < 27


def test_generator_continues_after_ordinary_unresolved_candidate(monkeypatch) -> None:
    """Ordinary non-separation is skipped until a later accepted proof."""
    calls = 0
    original = _evaluate_global_reach

    def delayed(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls <= 2:
            return dataclasses.replace(
                original(*args, **kwargs),
                status=_GlobalReachStatus.UNRESOLVED,
                reason="synthetic unresolved",
            )
        return original(*args, **kwargs)

    monkeypatch.setattr(
        "libephemeris.ellipsoid_contact._evaluate_global_reach", delayed
    )
    result = _generate_global_witness(
        _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
        _ConeSection(0.0, 0.8, 1),
        max_level=1,
    )
    assert result.status is _WitnessGenerationStatus.REACH
    assert len(result.attempted_candidates) >= 2


def test_generator_internal_invalid_is_terminal(monkeypatch) -> None:
    """A malformed generated payload is never skipped as ordinary unresolved."""
    original = _evaluate_global_reach

    def invalid(*args, **kwargs):
        return dataclasses.replace(
            original(*args, **kwargs),
            status=_GlobalReachStatus.INVALID,
            reason="synthetic generator defect",
        )

    monkeypatch.setattr(
        "libephemeris.ellipsoid_contact._evaluate_global_reach", invalid
    )
    result = _generate_global_witness(
        _frame(axis_point_km=(5.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
        _ConeSection(1.0, 1.0, 1),
        max_level=1,
    )
    assert result.status is _WitnessGenerationStatus.INVALID
    assert result.reason == "synthetic generator defect"
    assert result.work.total_candidates == 1


def test_bounded_generator_reports_exact_work_limit() -> None:
    """Operational exhaustion remains unresolved with deterministic counters."""
    result = _generate_global_witness(
        _frame(axis_point_km=(5.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
        _ConeSection(1.0, 1.0, 1),
        max_level=0,
        max_primal=1,
        max_dual=1,
        max_total=1,
    )
    assert result.status is _WitnessGenerationStatus.UNRESOLVED
    assert result.work.total_candidates == 1
    assert result.reason == "witness generation exhausted its candidate limit"


def test_global_reach_keeps_unproved_e1_equality_unresolved() -> None:
    """Matching zero-containing bounds are not an exact contact proof."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    proof = _evaluate_global_reach(
        frame,
        _ConeSection(1.0, 0.8, 1),
        (arb(-1.6), arb(0), arb(1.2)),
        _dual(scale=Fraction(4, 5)),
    )
    assert proof.status is _GlobalReachStatus.UNRESOLVED
    assert proof.equality_proof is None


def _zero_angle_equality(
    witness=None,
    relation=_ZeroAngleLineRelation.CROSSING,
) -> _ContactEqualityProof:
    if witness is None:
        witness = _dual(scale=Fraction(0))
    return _ContactEqualityProof(
        _ContactEqualityVariant.ZERO_ANGLE_LINE,
        witness,
        relation,
    )


@pytest.mark.parametrize("precision", [160, 256])
def test_global_reach_proves_exact_zero_angle_line_contact(precision) -> None:
    """The whole line identity proves exact L=U=0 at both start precisions."""
    previous = ctx.prec
    try:
        ctx.prec = precision
        frame = _frame(
            metric_km_minus_2=(
                (1.0, 0.8, 0.0),
                (0.8, 1.0, 0.0),
                (0.0, 0.0, 1.0),
            ),
            axis_point_km=(0.0, 1.5, 0.0),
            axis_span=(1.0, 0.0, 0.0),
        )
        equality = _ContactEqualityProof(
            _ContactEqualityVariant.ZERO_ANGLE_LINE,
            _dual(
                span=(Fraction(1), Fraction(0), Fraction(0)),
                scale=Fraction(0),
            ),
            _ZeroAngleLineRelation.CROSSING,
        )
        proof = _evaluate_global_reach(
            frame,
            _ConeSection(0.0, 1.0, 1),
            None,
            None,
            equality,
        )
    finally:
        ctx.prec = previous
    assert proof.status is _GlobalReachStatus.CONTACT
    assert proof.lower_bound.is_zero()  # type: ignore[union-attr]
    assert proof.upper_bound.is_zero()  # type: ignore[union-attr]
    assert proof.equality_proof is equality


@pytest.mark.parametrize("precision", [160, 256])
def test_global_reach_proves_exact_nondyadic_axis_tangency(precision) -> None:
    """A raw binary64 rational identity proves non-principal line tangency."""
    previous = ctx.prec
    try:
        ctx.prec = precision
        frame = _frame(
            axis_point_km=(-2.0, 1.0, 1.0),
            axis_span=(0.0, 1.0, 1.0),
        )
        equality = _ContactEqualityProof(
            _ContactEqualityVariant.ZERO_ANGLE_LINE,
            _dual(
                span=(Fraction(0), Fraction(1), Fraction(1)),
                scale=Fraction(0),
            ),
            _ZeroAngleLineRelation.TANGENT,
        )
        proof = _evaluate_global_reach(
            frame,
            _ConeSection(0.0, 1.0, 1),
            None,
            None,
            equality,
        )
    finally:
        ctx.prec = previous
    assert proof.status is _GlobalReachStatus.CONTACT
    assert proof.lower_bound.is_zero()  # type: ignore[union-attr]
    assert proof.upper_bound.is_zero()  # type: ignore[union-attr]


def test_global_reach_rejects_false_tangency_identity() -> None:
    """A tangency tag cannot promote an exact crossing to contact by assertion."""
    equality = _zero_angle_equality(relation=_ZeroAngleLineRelation.TANGENT)
    proof = _evaluate_global_reach(
        _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
        _ConeSection(0.0, 1.0, 1),
        None,
        None,
        equality,
    )
    assert proof.status is _GlobalReachStatus.INVALID
    assert proof.reason == "zero-angle tangency identity does not hold"


@pytest.mark.parametrize(
    ("frame", "cone", "equality", "reason"),
    [
        (
            _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
            _ConeSection(0.0, 0.8, 1),
            _zero_angle_equality(),
            "zero-angle contact requires an exact zero cone angle",
        ),
        (
            _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
            _ConeSection(1.0, 1.0, 1),
            _zero_angle_equality(),
            "zero-angle line contact requires zero cone radius",
        ),
        (
            _frame(axis_point_km=(3.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
            _ConeSection(0.0, 1.0, 1),
            _zero_angle_equality(),
            "zero-angle axis does not reach the ellipsoid",
        ),
        (
            _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
            _ConeSection(0.0, 1.0, 1),
            _zero_angle_equality(witness=_dual()),
            "zero-angle contact requires an exact zero dual vector",
        ),
        (
            _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
            _ConeSection(0.0, 1.0, 1),
            _zero_angle_equality(
                witness=_dual(scale=Fraction(0), multiplier=Fraction(1))
            ),
            "zero-angle contact requires a zero dual multiplier",
        ),
    ],
)
def test_global_reach_rejects_false_zero_angle_equalities(
    frame, cone, equality, reason
) -> None:
    """Every independent equality condition is fail-closed."""
    proof = _evaluate_global_reach(frame, cone, None, None, equality)
    assert proof.status is _GlobalReachStatus.INVALID
    assert proof.reason == reason


@pytest.mark.parametrize(
    ("signed_diameter", "shadow_class", "branch_sign"),
    [
        (-2.0, _CoreShadowClass.UMBRA, -1),
        (2.0, _CoreShadowClass.ANTUMBRA, 1),
        (0.0, _CoreShadowClass.APEX, 1),
        (-0.0, _CoreShadowClass.APEX, 1),
    ],
)
def test_core_cone_keeps_signed_class_separate_from_radius(
    signed_diameter, shadow_class, branch_sign
) -> None:
    """E7 uses one physical radius while retaining umbra/antumbra class."""
    section = _core_cone_section(signed_diameter, 1.0)
    assert section.cone.radius_km == abs(signed_diameter) / 2.0
    assert section.cone.branch_sign == branch_sign
    assert section.signed_diameter_km == signed_diameter
    assert section.shadow_class is shadow_class


def test_e7_signed_core_geometry_is_identical_apart_from_class() -> None:
    """E7's -2/+2 inputs cannot become negative physical cone lengths."""
    umbra = _core_cone_section(-2.0, 1.0)
    antumbra = _core_cone_section(2.0, 1.0)
    assert umbra.cone.radius_km == antumbra.cone.radius_km == 1.0
    assert umbra.cone.cosine == antumbra.cone.cosine == 1.0
    assert umbra.shadow_class is _CoreShadowClass.UMBRA
    assert antumbra.shadow_class is _CoreShadowClass.ANTUMBRA


@pytest.mark.parametrize("value", [math.inf, -math.inf, math.nan, 2, None])
def test_core_cone_rejects_nonfinite_or_nonnative_diameter(value) -> None:
    """Signed classification accepts only finite native runtime values."""
    with pytest.raises(ValueError, match="signed core diameter"):
        _core_cone_section(value, 1.0)


def test_global_reach_proves_e4_miss() -> None:
    """The exact E4 support witness gives the global lower bound +1."""
    root_two = math.sqrt(2.0)
    frame = _frame(
        metric_km_minus_2=(
            (1.0 / 9.0, 0.0, 0.0),
            (0.0, 1.0 / 9.0, 0.0),
            (0.0, 0.0, 0.25),
        ),
        axis_point_km=(-5.0, 0.0, 0.0),
        axis_span=(0.0, root_two, root_two),
    )
    span = tuple(Fraction.from_float(value) for value in frame.axis_span)
    proof = _evaluate_global_reach(
        frame,
        _ConeSection(1.0, 1.0, 1),
        None,
        _ProjectedDualWitness(
            span,
            (Fraction(1), Fraction(0), Fraction(0)),
            Fraction(1) / _dot_fraction(span, span),
            Fraction(0),
        ),
    )
    assert proof.status is _GlobalReachStatus.MISS
    assert proof.lower_bound > 0  # type: ignore[operator]


def _dot_fraction(left, right) -> Fraction:
    return sum((a * b for a, b in zip(left, right)), Fraction(0))


def test_global_reach_uses_strict_primal_upper_bound() -> None:
    """A feasible point strictly inside a zero-radius cone proves reach."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    proof = _evaluate_global_reach(
        frame,
        _ConeSection(0.0, 0.8, 1),
        (arb(0), arb(0), arb(1)),
        None,
    )
    assert proof.status is _GlobalReachStatus.REACH
    assert proof.upper_bound < 0  # type: ignore[operator]


def test_global_reach_classifies_internal_certification_failure_unresolved(
    monkeypatch,
) -> None:
    """A proof-engine failure is not mislabeled as invalid source data."""
    monkeypatch.setattr(
        "libephemeris.ellipsoid_contact._dual_vector",
        lambda *_: (_ for _ in ()).throw(IntervalCertificationError("synthetic")),
    )
    proof = _evaluate_global_reach(
        _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
        _ConeSection(1.0, 0.8, 1),
        None,
        _dual(),
    )
    assert proof.status is _GlobalReachStatus.UNRESOLVED
    assert proof.reason == "synthetic"


def test_cross_dual_payload_proves_same_e4_miss() -> None:
    """The exact cross-product witness implements the second approved variant."""
    frame = _frame(
        metric_km_minus_2=(
            (1.0 / 9.0, 0.0, 0.0),
            (0.0, 1.0 / 9.0, 0.0),
            (0.0, 0.0, 0.25),
        ),
        axis_point_km=(-5.0, 0.0, 0.0),
        axis_span=(0.0, 1.0, 1.0),
    )
    span = tuple(Fraction.from_float(value) for value in frame.axis_span)
    witness = _CrossDualWitness(
        span,
        (Fraction(0), Fraction(0), Fraction(-1)),
        Fraction(-1),
        Fraction(0),
    )
    proof = _evaluate_global_reach(frame, _ConeSection(1.0, 1.0, 1), None, witness)
    assert proof.status is _GlobalReachStatus.MISS
    assert proof.lower_bound > 0  # type: ignore[operator]


def test_global_reach_rejects_nonpoint_primal_and_missing_witnesses() -> None:
    """A primal witness is an exact point and at least one witness is required."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    nonpoint = _evaluate_global_reach(
        frame,
        _ConeSection(0.0, 0.8, 1),
        (arb(0, 1), arb(0), arb(1)),
        None,
    )
    assert nonpoint.status is _GlobalReachStatus.INVALID
    missing = _evaluate_global_reach(frame, _ConeSection(0.0, 0.8, 1), None, None)
    assert missing.status is _GlobalReachStatus.INVALID


@pytest.mark.parametrize(
    "metric",
    [
        (
            (1.0, 0.0, 0.0),
            (0.0, -1.0, 0.0),
            (0.0, 0.0, 1.0),
        ),
        (
            (-1.0, 0.0, 0.0),
            (0.0, -1.0, 0.0),
            (0.0, 0.0, 1.0),
        ),
        (
            (1.0, 2.0, 0.0),
            (2.0, 1.0, 0.0),
            (0.0, 0.0, -1.0),
        ),
    ],
)
def test_global_reach_classifies_non_spd_source_invalid(metric) -> None:
    """Every exactly indefinite source metric is invalid, not unresolved."""
    frame = _frame(
        metric_km_minus_2=metric,
        axis_point_km=(0.0, 0.0, 0.0),
        axis_span=(0.0, 0.0, 1.0),
    )
    proof = _evaluate_global_reach(frame, _ConeSection(1.0, 1.0, 1), None, _dual())
    assert proof.status is _GlobalReachStatus.INVALID
    assert proof.reason == "ellipsoid metric must be positive definite"


@pytest.mark.parametrize(
    "witness",
    [
        _ProjectedDualWitness(
            (Fraction(0),),
            (Fraction(1), Fraction(0), Fraction(0)),
            Fraction(1),
            Fraction(0),
        ),
        _ProjectedDualWitness(
            (Fraction(0), Fraction(0), Fraction(1)),
            (Fraction(1),),
            Fraction(1),
            Fraction(0),
        ),
        _ProjectedDualWitness(
            [Fraction(0), Fraction(0), Fraction(1)],  # type: ignore[arg-type]
            (Fraction(1), Fraction(0), Fraction(0)),
            Fraction(1),
            Fraction(0),
        ),
        _CrossDualWitness(
            (Fraction(0), Fraction(0), Fraction(1)),
            None,  # type: ignore[arg-type]
            Fraction(1),
            Fraction(0),
        ),
    ],
)
def test_global_reach_returns_invalid_for_malformed_dual_shapes(witness) -> None:
    """Malformed exact payload shapes stay inside the typed invalid result."""
    proof = _evaluate_global_reach(
        _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
        _ConeSection(1.0, 1.0, 1),
        None,
        witness,
    )
    assert proof.status is _GlobalReachStatus.INVALID
    assert proof.reason == "dual witness vectors must be exact three-component tuples"


def test_global_reach_rejects_wrong_dual_span_or_norm() -> None:
    """The exact dual payload stays bound to the frame and norm radius."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    wrong_span = _evaluate_global_reach(
        frame,
        _ConeSection(1.0, 0.8, 1),
        (arb(0), arb(0), arb(0)),
        _dual(span=(Fraction(0), Fraction(1), Fraction(0))),
    )
    assert wrong_span.status is _GlobalReachStatus.INVALID
    near_span = _evaluate_global_reach(
        frame,
        _ConeSection(1.0, 0.8, 1),
        None,
        _dual(span=(Fraction(0), Fraction(1), Fraction(1, 10**100))),
    )
    assert near_span.status is _GlobalReachStatus.INVALID
    too_large = _evaluate_global_reach(
        frame,
        _ConeSection(1.0, 0.8, 1),
        None,
        _dual(scale=Fraction(2)),
    )
    assert too_large.status is _GlobalReachStatus.UNRESOLVED


def test_frame_is_private_frozen_and_equal_by_value() -> None:
    """The internal frame is owned immutable data, not ambient metadata."""
    frame = _frame()
    assert frame == _frame()
    with pytest.raises(dataclasses.FrozenInstanceError):
        frame.axis_span = (1.0, 0.0, 0.0)  # type: ignore[misc]


def test_axis_line_derives_unit_direction_and_closest_anchor() -> None:
    """Normalization and orthogonality are construction identities in Arb."""
    line = _frame().certified_axis_line()
    norm_squared = sum((value * value for value in line.direction), 0)
    dot = sum(
        (line.direction[index] * line.anchor_km[index] for index in range(3)),
        0,
    )
    assert norm_squared.contains(1)
    assert dot.contains(0)
    assert line.direction[0].is_zero()
    assert line.direction[1].contains(0.6)
    assert line.direction[2].contains(0.8)
    assert all(value.is_finite() for value in line.anchor_km)


def test_axis_line_is_invariant_under_another_point_on_the_line() -> None:
    """Changing the supplied point along the same line preserves its anchor."""
    first = _frame(axis_point_km=(-3.5, 2.0, 6.0)).certified_axis_line()
    second = _frame(axis_point_km=(-3.5, 5.0, 10.0)).certified_axis_line()
    for left, right in zip(first.anchor_km, second.anchor_km):
        assert left.overlaps(right)


def test_metric_matrix_preserves_exact_binary64_inputs() -> None:
    """The metric enters Arb without decimal re-rounding."""
    metric = _frame().metric_ball_matrix()
    assert metric[0, 0].is_exact()
    assert metric[0, 0] == 0.25
    assert metric[2, 2] == 0.25


@pytest.mark.parametrize("sign", [0, 2, -2, True, 1.0])
def test_branch_signs_are_exactly_plus_or_minus_one(sign) -> None:
    """Nappe selection admits no Boolean, float, zero, or wider integer."""
    with pytest.raises(ValueError, match="sign"):
        _frame(core_branch_sign=sign).validate()


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("axis_span", (0.0, 0.0, 0.0), "squared norm"),
        ("axis_span", (math.inf, 0.0, 1.0), "finite"),
        ("axis_point_km", (0.0, math.nan, 0.0), "finite"),
        ("axis_point_km", [0.0, 0.0, 0.0], "tuple"),
    ],
)
def test_axis_vectors_fail_closed(field, value, message) -> None:
    """Malformed or non-finite source vectors never enter Arb geometry."""
    with pytest.raises((ValueError, IntervalCertificationError), match=message):
        _frame(**{field: value}).validate()


def test_regular_contact_matches_exact_e1_tangency() -> None:
    """E1 satisfies ellipsoid, cone, nappe, and stationarity equations."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    cone = _ConeSection(1.0, 0.8, 1)
    evaluation = _evaluate_regular_contact(
        frame,
        cone,
        (arb(-1.6), arb(0), arb(1.2)),
        arb(1),
    )
    assert evaluation.ellipsoid.contains(0)
    assert evaluation.cone.contains(0)
    assert evaluation.nappe_radius > 0
    assert evaluation.radial_distance > 0
    assert all(component.contains(0) for component in evaluation.stationarity)


def test_regular_contact_jacobian_matches_exact_cylinder_case() -> None:
    """A dyadic cylinder gives an independently explicit 4-by-4 Jacobian."""
    frame = _frame(
        metric_km_minus_2=(
            (0.25, 0.0, 0.0),
            (0.0, 0.25, 0.0),
            (0.0, 0.0, 0.25),
        ),
        axis_point_km=(0.0, 0.0, 0.0),
        axis_span=(0.0, 0.0, 1.0),
    )
    jacobian = _regular_contact_jacobian(
        frame,
        _ConeSection(1.0, 1.0, 1),
        (arb(-2), arb(0), arb(0)),
        arb(1),
    )
    expected = (
        (-1.0, 0.0, 0.0, 0.0),
        (0.5, 0.0, 0.0, -1.0),
        (0.0, 1.0, 0.0, 0.0),
        (0.0, 0.0, 0.5, 0.0),
    )
    assert (jacobian.nrows(), jacobian.ncols()) == (4, 4)
    for row in range(4):
        for column in range(4):
            assert jacobian[row, column].contains(expected[row][column])


def test_regular_contact_jacobian_encloses_e1_finite_difference() -> None:
    """E1 analytic columns contain independent centered-difference estimates."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    cone = _ConeSection(1.0, 0.8, 1)
    point = (-1.6, 0.0, 1.2)
    jacobian = _regular_contact_jacobian(
        frame, cone, tuple(arb(value) for value in point), arb(1)
    )
    step = 1e-6
    for column in range(3):
        low = list(point)
        high = list(point)
        low[column] -= step
        high[column] += step
        lower = _evaluate_regular_contact(
            frame, cone, tuple(arb(value) for value in low), arb(1)
        )
        upper = _evaluate_regular_contact(
            frame, cone, tuple(arb(value) for value in high), arb(1)
        )
        low_values = (lower.ellipsoid,) + lower.stationarity
        high_values = (upper.ellipsoid,) + upper.stationarity
        for row, (before, after) in enumerate(zip(low_values, high_values)):
            difference = (float(after.mid()) - float(before.mid())) / (2 * step)
            assert abs(float(jacobian[row, column].mid()) - difference) < 1e-8


def test_regular_root_box_certifies_e1_kkt_root() -> None:
    """Strict Krawczyk inclusion isolates the regular E1 stationary root."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    cone = _ConeSection(1.0, 0.8, 1)
    box = (
        arb(-1.6, "0.01"),
        arb(0, "0.01"),
        arb(1.2, "0.01"),
        arb(1, "0.01"),
    )
    certificate = _certify_regular_root_box(frame, cone, box)
    assert all(
        outer.contains_interior(inner)
        for outer, inner in zip(certificate.root_box, certificate.image)
    )
    assert certificate.cone_residual.contains(0)


@pytest.mark.parametrize(
    ("axis_point", "axis_span", "expected_contact"),
    [
        ((-4.0, 0.0, 0.0), (0.0, 0.0, 1.0), True),
        ((-4.0, 0.0, 0.0), (0.0, 1.0, 1.0), True),
        ((-5.0, 0.0, 0.0), (0.0, 1.0, 1.0), False),
    ],
)
def test_regular_root_box_covers_e2_e3_and_e4(
    axis_point, axis_span, expected_contact
) -> None:
    """Principal/inclined tangencies and the inclined miss stay distinct."""
    frame = _frame(
        metric_km_minus_2=(
            (1.0 / 9.0, 0.0, 0.0),
            (0.0, 1.0 / 9.0, 0.0),
            (0.0, 0.0, 0.25),
        ),
        axis_point_km=axis_point,
        axis_span=axis_span,
    )
    certificate = _certify_regular_root_box(
        frame,
        _ConeSection(1.0, 1.0, 1),
        (arb(-3, "0.001"), arb(0, "0.001"), arb(0, "0.001"), arb(1.5, "0.001")),
    )
    if expected_contact:
        assert certificate.cone_residual.contains(0)
    else:
        assert certificate.cone_residual > 0


def test_regular_residual_reduction_orders_supplied_candidates_only() -> None:
    """The reducer selects a separated minimum without claiming completeness."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    cone = _ConeSection(1.0, 0.8, 1)
    first = _certify_regular_root_box(
        frame,
        cone,
        (arb(-1.6, "0.01"), arb(0, "0.01"), arb(1.2, "0.01"), arb(1, "0.01")),
    )
    second = _RegularRootCertificate(
        frame,
        cone,
        first.root_box,
        first.image,
        arb(1.0, "0.01"),
    )
    with pytest.raises(IntervalCertificationError, match="sign"):
        _reduce_regular_residuals(cone, (first, second))
    separated = _RegularRootCertificate(
        frame,
        cone,
        first.root_box,
        first.image,
        arb(-0.5, "0.01"),
    )
    reduction = _reduce_regular_residuals(cone, (separated, second))
    assert reduction.minimum is separated.cone_residual
    assert reduction.candidate_count == 2
    assert reduction.sign == 1


def test_regular_residual_reduction_rejects_overlapping_candidates() -> None:
    """Overlapping residual ranges cannot establish one ordered candidate."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    cone = _ConeSection(1.0, 0.8, 1)
    first = _certify_regular_root_box(
        frame,
        cone,
        (arb(-1.6, "0.01"), arb(0, "0.01"), arb(1.2, "0.01"), arb(1, "0.01")),
    )
    with pytest.raises(IntervalCertificationError, match="minimum"):
        _reduce_regular_residuals(cone, (first, first))
    with pytest.raises(ValueError, match="non-empty"):
        _reduce_regular_residuals(cone, ())
    other_cone = _ConeSection(0.5, 0.8, 1)
    with pytest.raises(ValueError, match="one frame and cone"):
        _reduce_regular_residuals(other_cone, (first,))


def test_regular_root_box_rejects_noncontracting_or_negative_lambda_box() -> None:
    """A broad image or multiplier crossing zero cannot certify a root."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    cone = _ConeSection(1.0, 0.8, 1)
    with pytest.raises(IntervalCertificationError, match="strictly inside"):
        _certify_regular_root_box(
            frame,
            cone,
            (arb(-1.6, "0.3"), arb(0, "0.3"), arb(1.2, "0.3"), arb(1.5, "0.3")),
        )
    with pytest.raises(ValueError, match="multiplier"):
        _certify_regular_root_box(
            frame,
            cone,
            (arb(-1.6, "0.01"), arb(0, "0.01"), arb(1.2, "0.01"), arb(0, 1)),
        )


def test_regular_root_certificate_does_not_claim_cone_contact() -> None:
    """A KKT root with a changed cone radius retains a nonzero g enclosure."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    cone = _ConeSection(0.5, 0.8, 1)
    certificate = _certify_regular_root_box(
        frame,
        cone,
        (arb(-1.6, "0.01"), arb(0, "0.01"), arb(1.2, "0.01"), arb(1, "0.01")),
    )
    assert not certificate.cone_residual.contains(0)


def test_nappe_boundary_evaluator_accepts_zero_containing_radius() -> None:
    """The active r=0 equation is returned instead of forced into r>0."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    evaluation = _evaluate_nappe_boundary(
        frame,
        _ConeSection(1.0, 0.8, 1),
        (arb(1), arb(0), arb(-4) / 3),
        arb(0),
        arb(0),
    )
    assert evaluation.nappe_radius.contains(0)
    assert evaluation.radial_distance > 0
    assert all(component.is_finite() for component in evaluation.stationarity)


def test_nappe_boundary_jacobian_matches_centered_cylinder() -> None:
    """A zero-angle active plane gives an independently explicit Jacobian."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    jacobian = _nappe_boundary_jacobian(
        frame,
        _ConeSection(0.0, 1.0, 1),
        (arb(-2), arb(0), arb(0)),
        arb(1),
        arb(0),
    )
    expected = (
        (-1.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.0),
        (0.5, 0.0, 0.0, -1.0, 0.0),
        (0.0, 1.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.5, 0.0, 0.0),
    )
    assert (jacobian.nrows(), jacobian.ncols()) == (5, 5)
    for row in range(5):
        for column in range(5):
            assert jacobian[row, column].contains(expected[row][column])


def test_nappe_boundary_rejects_uncertified_multiplier_or_radial_axis() -> None:
    """The active set keeps non-negative mu and smooth rho contracts."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    with pytest.raises(ValueError, match="nappe multiplier"):
        _evaluate_nappe_boundary(
            frame,
            _ConeSection(1.0, 0.8, 1),
            (arb(1), arb(0), arb(-4) / 3),
            arb(0),
            arb(0, 1),
        )
    with pytest.raises(IntervalCertificationError, match="radial"):
        _evaluate_nappe_boundary(
            frame,
            _ConeSection(0.0, 0.8, 1),
            (arb(0), arb(0), arb(0)),
            arb(0),
            arb(0),
        )


def test_cone_apex_solves_joint_radial_and_nappe_equations() -> None:
    """The nonzero-angle apex is derived directly without a radial offset."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    apex = _evaluate_cone_apex(frame, _ConeSection(1.0, 0.8, 1))
    assert apex.axial_parameter.contains(-4.0 / 3.0)
    assert apex.point_km[0].is_zero()
    assert apex.point_km[1].is_zero()
    assert apex.point_km[2].contains(-4.0 / 3.0)
    assert apex.cone.contains(0)
    assert apex.ellipsoid < 0


def test_cone_apex_is_finite_in_offset_tilted_frame() -> None:
    """Construction identities avoid dependent interval subtraction at the axis."""
    apex = _evaluate_cone_apex(_frame(), _ConeSection(1.0, 0.8, -1))
    assert all(value.is_finite() for value in apex.point_km)
    assert apex.ellipsoid.is_finite()
    assert apex.cone.is_exact()
    assert apex.cone.is_zero()


def test_cone_apex_rejects_zero_or_unresolved_angle() -> None:
    """A cylinder has no finite apex and unresolved sine fails closed."""
    with pytest.raises(ValueError, match="no finite apex"):
        _evaluate_cone_apex(_frame(), _ConeSection(1.0, 1.0, 1))
    previous = ctx.prec
    try:
        ctx.prec = 32
        with pytest.raises(IntervalCertificationError, match="sine"):
            _evaluate_cone_apex(
                _frame(), _ConeSection(1.0, math.nextafter(1.0, 0.0), 1)
            )
    finally:
        ctx.prec = previous


def test_radial_subgradient_evaluates_axis_witness_constraints() -> None:
    """A centered point and zero witness expose exact nonsmooth residuals."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    evaluation = _evaluate_radial_subgradient(
        frame,
        _ConeSection(1.0, 1.0, 1),
        (arb(0), arb(0), arb(0)),
        (arb(0), arb(0), arb(0)),
        arb(0),
        arb(0),
    )
    assert evaluation.radial_squared.is_zero()
    assert evaluation.witness_axis_dot.is_zero()
    assert evaluation.witness_norm_squared.is_zero()
    assert evaluation.nappe_radius == 1
    assert evaluation.ellipsoid == -1
    assert all(component.is_zero() for component in evaluation.stationarity)


def test_radial_subgradient_retains_witness_inequality_residual() -> None:
    """An inadmissible witness is reported, not silently normalized or clamped."""
    evaluation = _evaluate_radial_subgradient(
        _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
        _ConeSection(1.0, 1.0, 1),
        (arb(0), arb(0), arb(0)),
        (arb(2), arb(0), arb(0)),
        arb(0),
        arb(0),
    )
    assert evaluation.witness_norm_squared == 4
    assert not evaluation.witness_norm_squared <= 1


def test_radial_subgradient_rejects_bad_witness_or_multiplier() -> None:
    """Nonsmooth KKT inputs retain exact Arb and non-negative contracts."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    with pytest.raises(ValueError, match="witness"):
        _evaluate_radial_subgradient(
            frame,
            _ConeSection(1.0, 1.0, 1),
            (arb(0), arb(0), arb(0)),
            (arb(0), 0.0, arb(0)),
            arb(0),
            arb(0),
        )
    with pytest.raises(ValueError, match="nappe multiplier"):
        _evaluate_radial_subgradient(
            frame,
            _ConeSection(1.0, 1.0, 1),
            (arb(0), arb(0), arb(0)),
            (arb(0), arb(0), arb(0)),
            arb(0),
            arb(0, 1),
        )


def test_regular_contact_requires_smooth_domain() -> None:
    """Axis and inactive-nappe points stay outside the regular KKT evaluator."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    with pytest.raises(IntervalCertificationError, match="radial"):
        _evaluate_regular_contact(
            frame, _ConeSection(1.0, 1.0, 1), (arb(0), arb(0), arb(0)), arb(0)
        )
    with pytest.raises(IntervalCertificationError, match="nappe"):
        _evaluate_regular_contact(
            frame,
            _ConeSection(0.0, 0.8, 1),
            (arb(1), arb(0), arb(-2)),
            arb(0),
        )


@pytest.mark.parametrize(
    "cone",
    [
        _ConeSection(-1.0, 0.8, 1),
        _ConeSection(1.0, 0.0, 1),
        _ConeSection(1.0, 1.1, 1),
        _ConeSection(1.0, 0.8, 0),
    ],
)
def test_cone_section_rejects_invalid_inputs(cone) -> None:
    """Physical cone data is never clamped or assigned a default branch."""
    with pytest.raises(ValueError):
        cone.validate()


@pytest.mark.parametrize("multiplier", [arb(-1), arb(0, 1), -1.0, None])
def test_regular_contact_rejects_uncertified_multiplier(multiplier) -> None:
    """The KKT multiplier enclosure must be wholly non-negative."""
    with pytest.raises(ValueError, match="multiplier"):
        _evaluate_regular_contact(
            _frame(),
            _ConeSection(1.0, 1.0, 1),
            (arb(-2), arb(0), arb(0)),
            multiplier,
        )


@pytest.mark.parametrize("point", [(arb(-2), 0.0, arb(0)), [arb(-2), arb(0), arb(0)]])
def test_regular_contact_rejects_malformed_point(point) -> None:
    """Wrong point containers and scalar types raise the declared ValueError."""
    with pytest.raises(ValueError, match="contact point"):
        _evaluate_regular_contact(_frame(), _ConeSection(1.0, 1.0, 1), point, arb(0))


def test_metric_rejects_asymmetry_and_nonpositive_definiteness() -> None:
    """The frame does not symmetrize or regularize an invalid ellipsoid."""
    with pytest.raises(ValueError, match="symmetric"):
        _frame(
            metric_km_minus_2=(
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
                (0.0, 0.0, 1.0),
            )
        ).validate()
    with pytest.raises(ValueError, match="positive definite"):
        _frame(
            metric_km_minus_2=(
                (1.0, 0.0, 0.0),
                (0.0, -1.0, 0.0),
                (0.0, 0.0, 1.0),
            )
        ).validate()
