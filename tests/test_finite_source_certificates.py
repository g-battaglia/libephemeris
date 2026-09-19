# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Synthetic exact checks for the disconnected finite-source verifier."""

from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
from fractions import Fraction
from typing import Any, cast

import pytest

from libephemeris._finite_source_certificates import (
    _DomainOutcome,
    _DualPayload,
    _OutsideDomain,
    _Target,
    _ValidatedGeometry,
    _prepare_geometry,
)

F = Fraction
Z = F(0)
IDENTITY = ((F(1), Z, Z), (Z, F(1), Z), (Z, Z, F(1)))
ZERO = (Z, Z, Z)


def geometry(B, p, r, m, A=IDENTITY, R=F(1)) -> _ValidatedGeometry:
    result = _prepare_geometry(B, p, r, m, A, R)
    assert isinstance(result, _ValidatedGeometry)
    return result


def test_antumbral_reach_and_separate_umbral_miss() -> None:
    geo = geometry((F(-12), Z, Z), (F(-7), Z, Z), F(2), F(1))
    assert geo.strict_reach(_Target.ANTUMBRA, ZERO)
    assert geo.strict_miss(_Target.UMBRA, _DualPayload(F(1), ZERO, Z))
    assert not geo.contains(_Target.UMBRA, ZERO)
    assert not geo.strict_reach(_Target.UMBRA, ZERO)


def test_reversed_orientation_has_three_separate_misses() -> None:
    geo = geometry((F(12), Z, Z), (F(17), Z, Z), F(2), F(1))
    assert geo.strict_miss(_Target.PENUMBRA, _DualPayload(F(1), ZERO, Z))
    assert geo.strict_miss(_Target.UMBRA, _DualPayload(Z, ZERO, F(1)))
    assert geo.strict_miss(_Target.ANTUMBRA, _DualPayload(F(1), ZERO, Z))
    for target in _Target:
        assert not geo.contains(target, ZERO)


def test_positive_squared_core_residual_does_not_override_gate() -> None:
    # The Moon is outside K but close to it. The squared umbral residual is
    # positive at the receiver centre while its axial gate is negative.
    geo = geometry((F(1, 10), F(-11), Z), (F(21, 10), Z, Z), F(10), F(1))
    assert geo.u == (F(2), F(11), Z)
    assert geo.U == F(125)
    v = (F(-21, 10), Z, Z)
    w = sum((geo.u[i] * v[i] for i in range(3)), Z)
    assert w == F(-21, 5)
    assert geo.U - (geo.r - geo.m) ** 2 == F(44)
    assert F(231, 10) ** 2 == F(53361, 100)
    a_u, gate, h_t, _ = geo._target_values(_Target.UMBRA, ZERO)
    assert a_u == F(814, 5)
    assert gate == F(-66, 5)
    assert h_t == F(586971, 25)
    assert a_u * a_u > h_t
    assert gate is not None and gate < 0
    assert geo.strict_reach(_Target.PENUMBRA, ZERO)
    assert not geo.contains(_Target.UMBRA, ZERO)


def test_apex_on_receiver_boundary_is_not_strict_reach() -> None:
    geo = geometry((F(-11), Z, Z), (F(-6), Z, Z), F(2), F(1))
    apex = (F(-1), Z, Z)
    assert geo.contains(_Target.UMBRA, apex)
    assert geo.contains(_Target.ANTUMBRA, apex)
    assert not geo.strict_reach(_Target.UMBRA, apex)
    assert not geo.strict_reach(_Target.ANTUMBRA, apex)
    assert geo.strict_reach(_Target.ANTUMBRA, ZERO)
    assert not geo.strict_reach(_Target.UMBRA, ZERO)


def test_coupled_metric_inverse_changes_dual_decision() -> None:
    metric = ((F(2), F(1), Z), (F(1), F(2), Z), (Z, Z, F(1)))
    geo = geometry(
        (F(3, 16), F(-1), F(3)),
        (F(19, 16), Z, F(3)),
        F(1, 2),
        F(1, 4),
        metric,
    )
    payload = _DualPayload(Z, ZERO, F(1))
    assert geo.strict_miss(_Target.PENUMBRA, payload)
    assert geo.A_inverse[0][1] == F(-1, 3)
    n = (F(1), F(1), Z)
    exact_support_squared = sum(
        (
            n[i] * sum((geo.A_inverse[i][j] * n[j] for j in range(3)), Z)
            for i in range(3)
        ),
        Z,
    )
    assert exact_support_squared == F(2, 3)
    assert sum((component * component for component in n), Z) == F(2)
    assert exact_support_squared < F(1) < F(2)


def test_nonzero_beta_orientation_and_exact_dual_coefficients() -> None:
    metric = ((F(2), F(1), Z), (F(1), F(2), Z), (Z, Z, F(1)))
    geo = geometry(
        (F(3, 16), F(-1), F(3)),
        (F(19, 16), Z, F(3)),
        F(1, 2),
        F(1, 4),
        metric,
    )
    assert geo.U - (geo.r + geo.m) ** 2 == F(23, 16)
    payload = _DualPayload(F(1), (Z, Z, F(-1)), Z)
    assert sum((component * component for component in payload.beta), Z) == F(1)
    a_origin, _, _, z_origin = geo._target_values(_Target.PENUMBRA, ZERO)
    assert a_origin == F(-25, 64)
    assert z_origin == (F(-3), F(3), F(19, 16))
    c = a_origin - z_origin[2]
    assert c == F(-101, 64)
    # The oriented cross term adds +x_0-x_1 to the penumbral a-gradient.
    n = (F(7, 4), F(-1, 4), Z)
    assert F(19, 8) == sum(
        (
            n[i] * sum((geo.A_inverse[i][j] * n[j] for j in range(3)), Z)
            for i in range(3)
        ),
        Z,
    )
    assert c * c == F(10201, 4096) > F(19, 8)
    assert geo.strict_miss(_Target.PENUMBRA, payload)
    # Reversing the cross-product term makes F(0)=51/64 positive.
    reversed_payload = _DualPayload(F(1), (Z, Z, F(1)), Z)
    assert a_origin + z_origin[2] == F(51, 64)
    assert not geo.strict_miss(_Target.PENUMBRA, reversed_payload)
    assert not geo.strict_miss(_Target.PENUMBRA, _DualPayload(F(1), (Z, Z, F(-2)), Z))


@pytest.mark.parametrize(
    "r,m", [(F(1), F(1)), (F(1), F(2)), (Z, F(1)), (F(1), Z), (Z, Z)]
)
def test_unsupported_slopes_are_typed_outcomes(r: Fraction, m: Fraction) -> None:
    assert (
        _prepare_geometry((F(-10), Z, Z), (F(10), Z, Z), r, m, IDENTITY, F(1))
        is _DomainOutcome.UNSUPPORTED_SLOPE
    )


def test_proved_exclusion_precedes_uncertified_enclosure() -> None:
    # Source and occulter balls overlap: U=1 <= (r+m)^2=9.
    result = _prepare_geometry((F(5), Z, Z), (F(6), Z, Z), F(2), F(1), IDENTITY, F(1))
    assert result is _DomainOutcome.EXCLUDED_DOMAIN


def test_conservative_source_distance_failure_is_uncertified() -> None:
    metric = ((F(1), Z, Z), (Z, F(100), Z), (Z, Z, F(100)))
    result = _prepare_geometry((Z, Z, F(7, 2)), (Z, Z, F(10)), F(3), F(1), metric, F(1))
    assert result is _DomainOutcome.UNCERTIFIED_DOMAIN
    # K reaches only z=1/10; the true source-to-K distance is 17/5 > 3.
    assert F(7, 2) - F(1, 10) == F(17, 5)


def test_insufficient_enclosure_and_receiver_overlap_are_uncertified() -> None:
    assert (
        _prepare_geometry((F(-12), Z, Z), (F(-7), Z, Z), F(2), F(1), IDENTITY, F(1, 2))
        is _DomainOutcome.UNCERTIFIED_DOMAIN
    )
    assert (
        _prepare_geometry((F(2), Z, Z), (F(10), Z, Z), F(2), F(1), IDENTITY, F(1))
        is _DomainOutcome.UNCERTIFIED_DOMAIN
    )


def test_wrong_dual_norm_gate_and_equality_separator() -> None:
    geo = geometry((F(12), Z, Z), (F(17), Z, Z), F(2), F(1))
    assert not geo.strict_miss(_Target.PENUMBRA, _DualPayload(Z, (F(1), Z, Z), Z))
    assert not geo.strict_miss(_Target.ANTUMBRA, _DualPayload(F(1), ZERO, F(1)))
    assert not geo.strict_miss(_Target.UMBRA, _DualPayload(F(1), ZERO, F(-1)))

    touching = geometry((F(3, 16), Z, F(3)), (F(19, 16), Z, F(3)), F(1, 2), F(1, 4))
    # F=g_p=-1+x_0 touches the ball at its boundary: equality is not MISS.
    assert not touching.strict_miss(_Target.PENUMBRA, _DualPayload(Z, ZERO, F(1)))


@pytest.mark.parametrize("bad", [F(-1), -1, 1.0])
def test_negative_or_nonfraction_radius_is_malformed(bad: object) -> None:
    if type(bad) is Fraction:
        with pytest.raises(ValueError):
            _prepare_geometry(
                (F(-10), Z, Z), (F(10), Z, Z), cast(Any, bad), F(1), IDENTITY, F(1)
            )
    else:
        with pytest.raises(TypeError):
            _prepare_geometry(
                (F(-10), Z, Z), (F(10), Z, Z), cast(Any, bad), F(1), IDENTITY, F(1)
            )


def test_negative_occulter_radius_is_malformed() -> None:
    with pytest.raises(ValueError):
        _prepare_geometry((F(-10), Z, Z), (F(10), Z, Z), F(2), F(-1), IDENTITY, F(1))


@pytest.mark.parametrize(
    "changes",
    [
        {"B": (0, Z, Z)},
        {"p": (Z, Z)},
        {"R": 1},
        {"R": F(0)},
        {"A": ((F(1), F(1), Z), (Z, F(1), Z), (Z, Z, F(1)))},
        {"A": ((F(1), Z, Z), (Z, F(-1), Z), (Z, Z, F(1)))},
    ],
)
def test_malformed_geometry_inputs(changes: dict[str, object]) -> None:
    args: dict[str, object] = {
        "B": (F(-12), Z, Z),
        "p": (F(-7), Z, Z),
        "r": F(2),
        "m": F(1),
        "A": IDENTITY,
        "R": F(1),
    }
    args.update(changes)
    with pytest.raises((TypeError, ValueError)):
        _prepare_geometry(**cast(Any, args))


def test_payload_and_membership_type_errors_are_distinct_from_false() -> None:
    geo = geometry((F(-12), Z, Z), (F(-7), Z, Z), F(2), F(1))
    with pytest.raises(TypeError):
        geo.contains(_Target.ANTUMBRA, cast(Any, (0, Z, Z)))
    with pytest.raises(ValueError):
        geo.contains(_Target.ANTUMBRA, (F(2), Z, Z))
    assert not geo.strict_reach(_Target.ANTUMBRA, (F(1), Z, Z))
    with pytest.raises(TypeError):
        geo.strict_reach(_Target.ANTUMBRA, cast(Any, (0, Z, Z)))
    with pytest.raises(TypeError):
        geo.strict_miss(_Target.UMBRA, _DualPayload(cast(Any, 1), ZERO, Z))
    with pytest.raises(TypeError):
        geo.strict_miss(_Target.UMBRA, _DualPayload(F(1), cast(Any, (0, Z, Z)), Z))
    with pytest.raises(TypeError):
        geo.strict_miss(_Target.UMBRA, _DualPayload(F(1), ZERO, cast(Any, 0)))
    with pytest.raises(TypeError):
        geo.contains(cast(Any, "umbra"), ZERO)
    with pytest.raises(TypeError):
        geo.strict_reach(cast(Any, "umbra"), ZERO)
    with pytest.raises(TypeError):
        geo.strict_miss(cast(Any, "umbra"), _DualPayload(F(1), ZERO, Z))


def test_direct_construction_derives_axis_and_inverse() -> None:
    geo = _ValidatedGeometry((F(-12), Z, Z), (F(-7), Z, Z), F(2), F(1), IDENTITY, F(1))
    assert geo.u == (F(5), Z, Z)
    assert geo.U == F(25)
    assert geo.A_inverse == IDENTITY
    assert geo.strict_miss(_Target.UMBRA, _DualPayload(F(1), ZERO, Z))
    with pytest.raises(AttributeError):
        getattr(geo, "__dict__")
    with pytest.raises(FrozenInstanceError):
        setattr(geo, "U", F(1))


@pytest.mark.parametrize(
    "B,p,r,m,A,R,outcome",
    [
        (
            (F(-10), Z, Z),
            (F(10), Z, Z),
            F(1),
            Z,
            IDENTITY,
            F(1),
            _DomainOutcome.UNSUPPORTED_SLOPE,
        ),
        (
            (F(5), Z, Z),
            (F(6), Z, Z),
            F(2),
            F(1),
            IDENTITY,
            F(1),
            _DomainOutcome.EXCLUDED_DOMAIN,
        ),
        (
            (F(-12), Z, Z),
            (F(-7), Z, Z),
            F(2),
            F(1),
            IDENTITY,
            F(1, 2),
            _DomainOutcome.UNCERTIFIED_DOMAIN,
        ),
    ],
)
def test_direct_construction_rejects_out_of_domain(
    B: object,
    p: object,
    r: object,
    m: object,
    A: object,
    R: object,
    outcome: _DomainOutcome,
) -> None:
    with pytest.raises(_OutsideDomain) as error:
        _ValidatedGeometry(
            cast(Any, B),
            cast(Any, p),
            cast(Any, r),
            cast(Any, m),
            cast(Any, A),
            cast(Any, R),
        )
    assert error.value.outcome is outcome


def test_direct_construction_cannot_supply_forged_derived_values() -> None:
    args = ((F(-12), Z, Z), (F(-7), Z, Z), F(2), F(1), IDENTITY, F(1))
    forged_inverse = ((Z, Z, Z), (Z, Z, Z), (Z, Z, Z))
    with pytest.raises(TypeError):
        cast(Any, _ValidatedGeometry)(*args, (F(99), Z, Z), F(1), forged_inverse)
    with pytest.raises(TypeError):
        cast(Any, _ValidatedGeometry)(*args, A_inverse=forged_inverse)
    geo = _ValidatedGeometry(*args)
    with pytest.raises(TypeError):
        replace(geo, A_inverse=forged_inverse)
