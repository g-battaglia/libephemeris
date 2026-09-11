"""Baseline-free certification tests for Gauquelin sector cusps."""

from __future__ import annotations

import importlib
import math
from fractions import Fraction

import pytest

import libephemeris as ephem

houses_module = importlib.import_module("libephemeris.houses")
from libephemeris.exceptions import CalculationError, PolarCircleError
from libephemeris.houses import (
    _gauge_alpha_to_longitude,
    _gauge_chart_piece,
    _gauge_derivative,
    _gauge_residual,
    _houses_gauquelin,
    house_pos,
)
from libephemeris.intervals import (
    IntervalCertificationError,
    ball_from_float,
    certified_float,
    interval_precision,
    isolate_unique_root,
)


EPS = 23.4392911
ARMC = 17.0
LAT = 48.0
ASC = 123.0
MC = 234.0


def test_strict_ordinary_outputs_and_cardinal_anchors() -> None:
    cusps = _houses_gauquelin(ARMC, LAT, EPS, ASC, MC)
    assert len(cusps) == 37
    assert all(isinstance(value, float) and 0.0 <= value < 360.0 for value in cusps[1:])
    assert cusps[1] == ASC
    assert cusps[10] == MC
    assert cusps[19] == (ASC + 180.0) % 360.0
    assert cusps[28] == (MC + 180.0) % 360.0


def test_all_non_cardinal_roots_are_in_semantic_order() -> None:
    _, angles = ephem.houses_armc(ARMC, LAT, EPS, ord("G"))
    asc, mc = angles[0], angles[1]
    cusps = _houses_gauquelin(ARMC, LAT, EPS, asc, mc)
    for sector in range(1, 19):
        assert cusps[sector + 18] == (cusps[sector] + 180.0) % 360.0
    # The circular clockwise coordinate is the public placement coordinate.
    for sector, cusp in enumerate(cusps[1:], 1):
        placement = house_pos(ARMC, LAT, EPS, "G", cusp, 0.0)
        assert math.isclose(placement, float(sector), abs_tol=3e-12)


def test_local_charts_agree_at_exact_positive_and_negative_seams() -> None:
    with pytest.raises(IntervalCertificationError):
        # A deliberately wide interval crossing two seams must fail closed.
        _gauge_alpha_to_longitude(
            ball_from_float(45.0).union(ball_from_float(135.0)), EPS
        )
    for seam in (45.0, -45.0, 135.0, -135.0):
        alpha = ball_from_float(seam)
        left_k = math.floor(seam / 90.0) * 90
        left = _gauge_chart_piece(alpha, EPS, left_k)
        right = _gauge_chart_piece(alpha, EPS, left_k + 90)
        assert math.isclose(float(left.mid()), float(right.mid()), abs_tol=1e-12)
        assert math.isclose(
            float(_gauge_alpha_to_longitude(alpha, EPS).mid()),
            float(left.mid()),
            abs_tol=1e-12,
        )


def test_eps_zero_midpoint_uses_declared_directional_candidate() -> None:
    assert (
        houses_module._gauquelin_cusp_for_sector(2, 2**-45, 0.0, 0.0) == 350.0 + 2**-44
    )


def test_eps_zero_is_exact_clockwise_ring() -> None:
    asc = math.nextafter(123.0, math.inf)
    cusps = _houses_gauquelin(ARMC, LAT, 0.0, asc, MC)
    expected = Fraction.from_float(asc)
    for sector in range(1, 37):
        expected_value = float((expected - 10 * (sector - 1)) % 360)
        assert math.isclose(cusps[sector], expected_value, abs_tol=1e-12)


def test_equator_keeps_equal_ra_spacing_but_nonlinear_longitude() -> None:
    cusps = _houses_gauquelin(ARMC, 0.0, EPS, ASC, MC)
    assert len({round(cusps[i], 8) for i in range(2, 10)}) == 8
    assert cusps[2] != pytest.approx(ASC - 10.0)


def test_half_turn_and_whole_turn_invariance() -> None:
    base = _houses_gauquelin(ARMC, LAT, EPS, ASC, MC)
    shifted = _houses_gauquelin(ARMC + 360.0, LAT, EPS, ASC + 360.0, MC + 360.0)
    assert all(math.isclose(a, b, abs_tol=3e-13) for a, b in zip(shifted, base))
    for sector in range(1, 19):
        assert base[sector + 18] == (base[sector] + 180.0) % 360.0


def test_tangent_and_exact_pole_zero_obliquity_are_deterministic() -> None:
    first = _houses_gauquelin(270.0, 30.0, 60.0, ASC, MC)
    second = _houses_gauquelin(270.0, 30.0, 60.0, ASC, MC)
    assert first == second
    assert first[1] == ASC
    assert first[2] == (ASC - 10.0) % 360.0
    pole = _houses_gauquelin(270.0, 90.0, 0.0, ASC, MC)
    assert pole == _houses_gauquelin(270.0, 90.0, 0.0, ASC, MC)


def test_first_representable_frame_beyond_tangent_refuses_publicly() -> None:
    with pytest.raises(PolarCircleError):
        ephem.houses_armc(270.0, 30.0, math.nextafter(60.0, math.inf), ord("G"))


def test_private_beyond_polar_fallback_is_explicit() -> None:
    cusps = _houses_gauquelin(270.0, 30.0, math.nextafter(60.0, math.inf), ASC, MC)
    assert cusps[1] == ASC
    assert cusps[2] == (ASC - 10.0) % 360.0


def test_production_uses_both_precisions_and_agrees(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen: list[int] = []
    original = houses_module.isolate_unique_root

    def wrapped(*args, **kwargs):
        seen.append(kwargs["precision"])
        return original(*args, **kwargs)

    monkeypatch.setattr(houses_module, "isolate_unique_root", wrapped)
    houses_module._gauquelin_cusp_for_sector(4, ARMC, LAT, EPS)
    assert seen == [192, 256]


def test_cross_precision_disagreement_fails_closed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original = houses_module._gauge_alpha_to_longitude
    calls = 0

    def disagree(alpha, eps):
        nonlocal calls
        calls += 1
        result = original(alpha, eps)
        if calls == 2:
            return result + ball_from_float(1.0)
        return result

    monkeypatch.setattr(houses_module, "_gauge_alpha_to_longitude", disagree)
    with pytest.raises(CalculationError):
        houses_module._gauquelin_cusp_for_sector(4, ARMC, LAT, EPS)


def test_complete_interval_residual_and_derivative_certificate() -> None:
    coefficient = math.tan(math.radians(LAT)) * math.tan(math.radians(EPS))
    with interval_precision(256):
        armc = ball_from_float(ARMC)
        coeff = ball_from_float(coefficient)
        for sector, bounds in (
            (2, (-180.0, 0.0)),
            (12, (0.0, 180.0)),
            (22, (0.0, 180.0)),
            (32, (180.0, 360.0)),
        ):
            root = isolate_unique_root(
                lambda h: _gauge_residual(sector, armc, h, coeff),
                lambda h: _gauge_derivative(sector, armc, h, coeff),
                *bounds,
                precision=256,
            )
            assert bounds[0] < float(root.mid()) < bounds[1]
            assert not _gauge_derivative(sector, armc, root, coeff).contains(0)
            residual = _gauge_residual(sector, armc, root.mid(), coeff)
            assert abs(float(residual.mid())) < 1e-30


def test_unresolved_nonlinear_midpoint_fails_closed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def unresolved(_value):
        raise IntervalCertificationError("synthetic unresolved midpoint")

    monkeypatch.setattr(houses_module, "certified_float", unresolved)
    with pytest.raises(CalculationError):
        houses_module._gauquelin_cusp_for_sector(2, 0.0, LAT, EPS)


def test_rounding_cell_overlap_fails_closed() -> None:
    candidate = 1.0
    next_value = math.nextafter(candidate, math.inf)
    midpoint = (Fraction.from_float(candidate) + Fraction.from_float(next_value)) / 2
    ambiguous = ball_from_float(float(midpoint)).union(ball_from_float(next_value))
    with pytest.raises(IntervalCertificationError):
        certified_float(ambiguous)


def test_a_zero_sweep_has_no_valid_root_refusals() -> None:
    for lat in (-0.0, 0.0):
        for eps in (5.0, EPS, 60.0):
            for armc in (
                0.0,
                math.nextafter(0.0, -math.inf),
                math.nextafter(0.0, math.inf),
                90.0,
                180.0,
                270.0,
            ):
                for sector in range(2, 37):
                    if sector in (10, 19, 28):
                        continue
                    result = houses_module._gauquelin_cusp_for_sector(
                        sector, armc, lat, eps
                    )
                    assert 0.0 <= result < 360.0


def test_private_signature_is_current_geometry_only() -> None:
    import inspect
    from libephemeris.houses import _gauquelin_cusp_for_sector

    names = list(inspect.signature(_gauquelin_cusp_for_sector).parameters)
    assert names == ["sector", "armc", "lat", "eps"]
    with pytest.raises(CalculationError):
        _gauquelin_cusp_for_sector(1, ARMC, LAT, EPS)
