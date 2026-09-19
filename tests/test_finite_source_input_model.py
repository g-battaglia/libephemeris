# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Synthetic native-word arithmetic tests for the finite-source bridge."""

from __future__ import annotations

import math
import struct
import sys
from dataclasses import replace
from decimal import Decimal
from fractions import Fraction
from typing import Any, cast

import pytest

from libephemeris._finite_source_certificates import _DomainOutcome, _ValidatedGeometry
from libephemeris._finite_source_input_model import (
    _FiniteSourceInput,
    _InvalidInputModel,
    _NativeConstants,
    _NativeWords,
)

F = Fraction


def word_bits(value: float) -> int:
    return struct.unpack("!Q", struct.pack("!d", value))[0]


def constants() -> _NativeConstants:
    c = 149597870.7
    sun_radius_km = 696000.0
    moon_radius_km = 1738.15
    earth_axis_km = 6378.140
    sun_radius_au = sun_radius_km / c
    moon_radius_au = moon_radius_km / c
    earth_axis_au = earth_axis_km / c
    return _NativeConstants(
        au_km=c,
        sun_radius_km=sun_radius_km,
        moon_radius_km=moon_radius_km,
        earth_axis_km=earth_axis_km,
        flattening_unit=1.0,
        flattening_denominator=298.25642,
        sun_radius_au=sun_radius_au,
        moon_radius_au=moon_radius_au,
        earth_axis_au=earth_axis_au,
        r_native=sun_radius_au * c,
        m_native=moon_radius_au * c,
        a_native=earth_axis_au * c,
        f_native=1.0 / 298.25642,
    )


def words(
    moon_au: tuple[float, float, float] = (1.0, -0.0, 0.0),
    sun_au: tuple[float, float, float] = (3.0, 0.0, -0.0),
) -> _NativeWords:
    chain = constants()
    return _NativeWords(
        tjd_ut=2451545.0,
        moon_au=moon_au,
        moon_km=(
            moon_au[0] * chain.au_km,
            moon_au[1] * chain.au_km,
            moon_au[2] * chain.au_km,
        ),
        sun_au=sun_au,
        sun_km=(
            sun_au[0] * chain.au_km,
            sun_au[1] * chain.au_km,
            sun_au[2] * chain.au_km,
        ),
        constants=chain,
    )


def test_exact_centres_metric_and_fixed_native_constant_facts() -> None:
    captured = words()
    result = _FiniteSourceInput(captured)
    assert result.fixed_product_revision == "1147e6a907b35d6d52613e56c1986c094a5897c9"
    assert result.p == tuple(F.from_float(value) for value in captured.moon_km)
    assert result.B == tuple(F.from_float(value) for value in captured.sun_km)
    assert result.u == tuple(result.p[i] - result.B[i] for i in range(3))
    assert result.r == F.from_float(captured.constants.r_native)
    assert result.m == F.from_float(captured.constants.m_native)
    a = F.from_float(captured.constants.a_native)
    f = F.from_float(captured.constants.f_native)
    b = a * (1 - f)
    assert result.R == a
    assert result.A == (
        (1 / a**2, F(0), F(0)),
        (F(0), 1 / a**2, F(0)),
        (F(0), F(0), 1 / b**2),
    )
    assert result.A[2][2] > result.A[0][0] > 0
    assert F.from_float(captured.constants.au_km) == F(1495978707, 10) - F(1, 83886080)
    assert result.r - F(696000) == 0
    assert result.m - F(173815, 100) == F(1, 10995116277760)
    assert result.R - F(6378140, 1000) == F(9, 27487790694400)
    assert f - F(50000, 14912821) == F(3751413, 17193312025252584327479296)
    moon_delta, sun_delta = result.iau_differences
    for au, km, delta in (
        (captured.moon_au, captured.moon_km, moon_delta),
        (captured.sun_au, captured.sun_km, sun_delta),
    ):
        for i in range(3):
            assert delta[i] == F.from_float(km[i]) - F.from_float(au[i]) * F(
                1495978707, 10
            )
    assert isinstance(result.theorem_domain(), _ValidatedGeometry)


@pytest.mark.parametrize("role", ["moon", "sun"])
@pytest.mark.parametrize("index", [0, 1, 2])
def test_each_au_word_is_bitwise_checked_against_km(role: str, index: int) -> None:
    captured = words()
    original = getattr(captured, f"{role}_au")
    changed = list(original)
    changed[index] = math.nextafter(original[index], math.inf)
    assert word_bits(changed[index] * captured.constants.au_km) != word_bits(
        getattr(captured, f"{role}_km")[index]
    )
    changed_vector = (changed[0], changed[1], changed[2])
    with pytest.raises(_InvalidInputModel, match=f"{role}_km\\[{index}\\]"):
        _FiniteSourceInput(
            cast(Any, replace)(captured, **{f"{role}_au": changed_vector})
        )


def test_changed_km_word_is_rejected() -> None:
    captured = words()
    moon_km = list(captured.moon_km)
    moon_km[0] = math.nextafter(moon_km[0], math.inf)
    with pytest.raises(_InvalidInputModel, match="moon_km"):
        _FiniteSourceInput(
            replace(captured, moon_km=(moon_km[0], moon_km[1], moon_km[2]))
        )


@pytest.mark.parametrize(
    "field",
    [
        "au_km",
        "sun_radius_km",
        "moon_radius_km",
        "earth_axis_km",
        "flattening_unit",
        "flattening_denominator",
        "sun_radius_au",
        "moon_radius_au",
        "earth_axis_au",
        "r_native",
        "m_native",
        "a_native",
        "f_native",
    ],
)
def test_each_native_constant_or_chain_word_is_pinned(field: str) -> None:
    captured = words()
    old = getattr(captured.constants, field)
    changed = math.nextafter(old, math.inf)
    mutated = replace(captured.constants, **{field: changed})
    with pytest.raises(_InvalidInputModel, match=field):
        _FiniteSourceInput(replace(captured, constants=mutated))


def test_signed_zero_bits_are_retained_and_checked() -> None:
    captured = words()
    result = _FiniteSourceInput(captured)
    assert result.p[1] == result.B[1] == F(0)
    assert word_bits(result.words.moon_au[1]) == word_bits(-0.0)
    assert word_bits(result.words.moon_km[1]) == word_bits(-0.0)
    assert word_bits(result.words.sun_au[1]) == word_bits(0.0)
    assert word_bits(result.words.sun_km[2]) == word_bits(-0.0)
    mutated = replace(captured, moon_km=(captured.moon_km[0], 0.0, 0.0))
    with pytest.raises(_InvalidInputModel, match="moon_km\\[1\\]"):
        _FiniteSourceInput(mutated)


@pytest.mark.parametrize("bad", [1, F(1), Decimal("1")])
def test_non_native_coordinate_words_are_rejected(bad: object) -> None:
    captured = words()
    mutated = replace(captured, moon_au=(cast(Any, bad), -0.0, 0.0))
    with pytest.raises(TypeError, match="native Python float"):
        _FiniteSourceInput(mutated)


@pytest.mark.parametrize("bad", [math.inf, math.nan])
def test_nonfinite_coordinate_words_are_invalid(bad: float) -> None:
    captured = words()
    mutated = replace(captured, moon_au=(bad, -0.0, 0.0))
    with pytest.raises(_InvalidInputModel, match="must be finite"):
        _FiniteSourceInput(mutated)


def test_nonfinite_epoch_and_overflowing_product_are_invalid() -> None:
    with pytest.raises(_InvalidInputModel):
        _FiniteSourceInput(replace(words(), tjd_ut=math.inf))
    captured = words()
    overflowing_au = replace(captured, moon_au=(sys.float_info.max, -0.0, 0.0))
    with pytest.raises(_InvalidInputModel, match="native product"):
        _FiniteSourceInput(overflowing_au)


def test_exact_u_can_differ_from_rounded_native_span() -> None:
    result = _FiniteSourceInput(words(sun_au=(1e20, 0.0, 0.0)))
    native_span_exact = F.from_float(result.native_span[0])
    assert result.u[0] == result.p[0] - result.B[0]
    assert native_span_exact != result.u[0]
    assert result.native_span[0] == -result.words.sun_km[0]


def test_invalid_arithmetic_is_distinct_from_valid_domain_outcomes() -> None:
    excluded = _FiniteSourceInput(words(sun_au=(1.0, 0.0, 0.0)))
    assert excluded.theorem_domain() is _DomainOutcome.EXCLUDED_DOMAIN
    uncertified = _FiniteSourceInput(words(moon_au=(0.0, 0.0, 0.0)))
    assert uncertified.theorem_domain() is _DomainOutcome.UNCERTIFIED_DOMAIN
    bad = replace(words(), moon_km=(0.0, 0.0, 0.0))
    with pytest.raises(_InvalidInputModel):
        _FiniteSourceInput(bad)


@pytest.mark.parametrize("field", ["u", "span", "A", "R"])
def test_independent_derived_inputs_are_not_in_schema(field: str) -> None:
    captured = words()
    with pytest.raises(TypeError):
        cast(Any, _NativeWords)(**{**vars_for_words(captured), field: 0.0})
    with pytest.raises(TypeError):
        cast(Any, _FiniteSourceInput)(captured, **{field: 0.0})


def vars_for_words(captured: _NativeWords) -> dict[str, object]:
    return {
        "tjd_ut": captured.tjd_ut,
        "moon_au": captured.moon_au,
        "moon_km": captured.moon_km,
        "sun_au": captured.sun_au,
        "sun_km": captured.sun_km,
        "constants": captured.constants,
    }


def test_coherent_word_change_can_pass_without_a_source_receipt() -> None:
    captured = words()
    changed_moon = (2.0, captured.moon_au[1], captured.moon_au[2])
    changed_km = (
        changed_moon[0] * captured.constants.au_km,
        changed_moon[1] * captured.constants.au_km,
        changed_moon[2] * captured.constants.au_km,
    )
    changed = replace(captured, moon_au=changed_moon, moon_km=changed_km)
    assert isinstance(_FiniteSourceInput(changed), _FiniteSourceInput)
    assert changed.moon_km[0] != captured.moon_km[0]
