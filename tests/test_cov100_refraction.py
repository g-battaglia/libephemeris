# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for :mod:`libephemeris.refraction`.

These tests drive the ICAO standard-atmosphere refraction model through
its edge-case branches: the ICAO temperature/pressure profiles, the
modified (observer-supplied) atmosphere, the refraction integral guard,
the Newton-Raphson inverter, and the horizon-dip helper.

The goal is exercising specific source branches deterministically with
fixed inputs; no network access is required.
"""

from __future__ import annotations

import math

from libephemeris import refraction as refr


# ---------------------------------------------------------------------------
# Gauss-Legendre node helper
# ---------------------------------------------------------------------------


def test_gauss_legendre_nodes_symmetry() -> None:
    """A small node set is symmetric and its weights sum to 2."""
    nodes, weights = refr._gauss_legendre_nodes(5)
    assert len(nodes) == 5
    assert len(weights) == 5
    # Nodes are symmetric about 0.
    for k in range(5):
        assert math.isclose(nodes[k], -nodes[4 - k], abs_tol=1e-12)
    # Weights integrate the constant function to 2 over [-1, 1].
    assert math.isclose(sum(weights), 2.0, abs_tol=1e-12)


# ---------------------------------------------------------------------------
# ICAO temperature profile (lines 184-188)
# ---------------------------------------------------------------------------


def test_icao_temperature_across_layers() -> None:
    """Temperature is read from the correct ICAO layer, incl. below sea level."""
    # Troposphere base.
    assert math.isclose(refr._icao_temperature(0.0), 288.15, abs_tol=1e-9)
    # Isothermal tropopause.
    assert math.isclose(refr._icao_temperature(15_000.0), 216.65, abs_tol=1e-9)
    # Stratosphere 1 (positive lapse).
    assert refr._icao_temperature(25_000.0) > 216.65
    # Below the lowest layer base -> falls through to the surface value.
    assert math.isclose(refr._icao_temperature(-100.0), 288.15, abs_tol=1e-9)


# ---------------------------------------------------------------------------
# ICAO pressure profile (lines 193-208)
# ---------------------------------------------------------------------------


def test_icao_pressure_across_layers() -> None:
    """Pressure follows the hydrostatic integral and decreases with altitude."""
    p0 = refr._icao_pressure(0.0)
    assert math.isclose(p0, refr._P0, abs_tol=1e-9)
    # Isothermal layer (lapse == 0): exponential branch.
    p_iso = refr._icao_pressure(15_000.0)
    assert 0.0 < p_iso < p0
    # Above the modeled top: the loop walks every layer then returns at the end.
    p_top = refr._icao_pressure(90_000.0)
    assert 0.0 < p_top < p_iso


# ---------------------------------------------------------------------------
# Refractive index guard (line 215)
# ---------------------------------------------------------------------------


def test_refractive_index_nonpositive_inputs() -> None:
    """Non-positive temperature or pressure returns vacuum index 1.0."""
    assert refr._refractive_index(1000.0, 0.0) == 1.0
    assert refr._refractive_index(0.0, 288.15) == 1.0
    # A valid pair returns an index slightly above 1.
    assert refr._refractive_index(refr._P0, refr._T0) > 1.0


# ---------------------------------------------------------------------------
# Modified temperature fall-through (line 242)
# ---------------------------------------------------------------------------


def test_modified_temperature_above_atmosphere_top() -> None:
    """A height above the modeled top falls through to the 200 K sentinel."""
    t = refr._modified_temperature(100_000.0, 0.0, 288.15, 0.0065)
    assert t == 200.0


def test_modified_temperature_in_tropo_and_strato() -> None:
    """Tropospheric and stratospheric branches return finite temperatures."""
    # Below 11 km: observer-driven linear profile.
    assert refr._modified_temperature(5_000.0, 0.0, 288.15, 0.0065) < 288.15
    # Between 11 and 20 km: i == 1 branch.
    assert refr._modified_temperature(15_000.0, 0.0, 288.15, 0.0065) > 0.0
    # Above 20 km: offset branch.
    assert refr._modified_temperature(25_000.0, 0.0, 288.15, 0.0065) > 0.0


# ---------------------------------------------------------------------------
# Modified pressure isothermal branches (lines 255, 261, 287)
# ---------------------------------------------------------------------------


def test_modified_pressure_isothermal_tropo() -> None:
    """Zero lapse inside the troposphere uses the exponential branch."""
    p = refr._modified_pressure(5_000.0, 0.0, refr._P0, 288.15, 0.0)
    assert 0.0 < p < refr._P0


def test_modified_pressure_isothermal_above_tropo() -> None:
    """Zero lapse with a height above the troposphere walks all layers."""
    # h above the modeled top hits the isothermal P_tropo branch (line 261),
    # iterates every layer, and returns at the final fall-through (line 287).
    p = refr._modified_pressure(100_000.0, 0.0, refr._P0, 288.15, 0.0)
    assert 0.0 < p < refr._P0


def test_modified_pressure_lapse_above_tropo() -> None:
    """Non-zero lapse, height above troposphere: power-law layer integrals."""
    p = refr._modified_pressure(25_000.0, 0.0, refr._P0, 288.15, 0.0065)
    assert 0.0 < p < refr._P0


# ---------------------------------------------------------------------------
# Refraction integral guard: arg <= 0 -> continue (line 445)
# ---------------------------------------------------------------------------


def test_refraction_integral_ray_below_invariant() -> None:
    """Quadrature points the ray cannot reach are skipped (arg <= 0)."""
    # Build a controlled profile where every node has n*r < C so that
    # arg = (n*r)^2 / C^2 - 1 <= 0 for all of them, forcing the continue.
    prof = object.__new__(refr._AtmosphereProfile)
    prof.n_obs = 1.0003
    prof.r_obs = refr._R_EARTH
    prof.r_top = refr._R_EARTH + refr._ATMO_TOP
    prof.L = prof.r_top - prof.r_obs
    n = refr._N_QUAD
    prof._q_r = [prof.r_obs] * n
    prof._q_n = [1.0] * n
    prof._q_dn_dr = [0.0] * n
    # At z = 90 deg, C = n_obs * r_obs > n*r = r_obs, so arg < 0 everywhere.
    result = refr._refraction_integral(90.0, prof)
    assert result == 0.0


def test_refraction_integral_normal_profile() -> None:
    """A genuine profile yields a small positive refraction near the horizon."""
    prof = refr._get_profile(0.0, refr._P0, 288.15, 0.0065)
    r = refr._refraction_integral(85.0, prof)
    assert r > 0.0


# ---------------------------------------------------------------------------
# _trace_ray guards (lines 480, 482)
# ---------------------------------------------------------------------------


def test_trace_ray_nonpositive_zenith() -> None:
    """A source at or above the zenith has zero refraction."""
    assert refr._trace_ray(0.0) == 0.0
    assert refr._trace_ray(-5.0) == 0.0


def test_trace_ray_nonpositive_pressure() -> None:
    """Zero/negative pressure (vacuum) yields zero refraction."""
    assert refr._trace_ray(45.0, obs_pressure=0.0) == 0.0
    assert refr._trace_ray(45.0, obs_pressure=-1.0) == 0.0


def test_trace_ray_below_horizon_extrapolation() -> None:
    """Below the geometric horizon uses linear extrapolation."""
    r = refr._trace_ray(92.0)
    assert r >= 0.0


# ---------------------------------------------------------------------------
# calc_refraction_true_to_app guard (line 541)
# ---------------------------------------------------------------------------


def test_true_to_app_nonpositive_pressure() -> None:
    """Non-positive pressure short-circuits to zero refraction."""
    assert refr.calc_refraction_true_to_app(45.0, pressure=0.0) == 0.0
    assert refr.calc_refraction_true_to_app(45.0, pressure=-10.0) == 0.0


def test_true_to_app_normal() -> None:
    """A standard altitude gives a positive, sub-degree refraction."""
    r = refr.calc_refraction_true_to_app(10.0)
    assert 0.0 < r < 1.0


# ---------------------------------------------------------------------------
# calc_refraction_app_to_true: guard (577), loop exhaustion (582->600)
# ---------------------------------------------------------------------------


def test_app_to_true_nonpositive_pressure() -> None:
    """Non-positive pressure short-circuits to zero refraction."""
    assert refr.calc_refraction_app_to_true(45.0, pressure=0.0) == 0.0
    assert refr.calc_refraction_app_to_true(45.0, pressure=-3.0) == 0.0


def test_app_to_true_loop_exhausts_without_break() -> None:
    """At the horizon the Newton loop runs all iterations (no early break)."""
    # apparent_alt = 0 does not converge to the 1e-12 residual nor hit the
    # degenerate-derivative break within 8 iterations, so the loop falls
    # through to the final return (arc 582->600).
    r = refr.calc_refraction_app_to_true(0.0)
    assert r >= 0.0


def test_app_to_true_normal_convergence() -> None:
    """A high altitude converges and breaks on the residual test."""
    r = refr.calc_refraction_app_to_true(45.0)
    assert 0.0 < r < 0.1


# ---------------------------------------------------------------------------
# calc_refraction_app_to_true: degenerate derivative break (line 597)
# ---------------------------------------------------------------------------


def test_app_to_true_degenerate_derivative(monkeypatch) -> None:
    """A zero apparent-altitude derivative triggers the early break.

    The inverter evaluates the true-to-apparent map twice per iteration:
    first at ``true_est`` (the residual call), then at ``true_est + dt``
    (the forward-difference call). The stateful fake returns ``-0.002``
    then ``-0.003`` so the forward difference is exactly ``-dt`` and the
    apparent-altitude derivative ``dapp_dtrue = 1 + (r_plus - r) / dt``
    collapses to exactly ``0.0`` (line 597). The base value ``-0.002``
    keeps the residual non-zero so the residual break is skipped first.
    """
    calls = {"n": 0}

    def fake_true_to_app(
        true_alt,
        pressure=refr._P0,
        temperature_C=15.0,
        obs_alt=0.0,
        lapse_rate=0.0065,
    ):
        calls["n"] += 1
        # Odd call -> residual evaluation; even call -> derivative step.
        return -0.002 if calls["n"] % 2 == 1 else -0.003

    monkeypatch.setattr(refr, "calc_refraction_true_to_app", fake_true_to_app)
    out = refr.calc_refraction_app_to_true(10.0)
    assert out >= 0.0


# ---------------------------------------------------------------------------
# calc_dip: ratio >= 1.0 (line 647) and d < 0 clamp (line 660)
# ---------------------------------------------------------------------------


def test_calc_dip_nonpositive_altitude() -> None:
    """At or below sea level the dip is zero."""
    assert refr.calc_dip(0.0) == 0.0
    assert refr.calc_dip(-5.0) == 0.0


def test_calc_dip_ratio_rounds_to_one() -> None:
    """A sub-nanometre altitude rounds the radius ratio to exactly 1.0."""
    # obs_alt > 0 passes the sea-level guard, but R/(R+h) rounds to 1.0,
    # so the geometric dip is degenerate and the function returns 0.0.
    assert refr.calc_dip(1e-10) == 0.0


def test_calc_dip_negative_discriminant_clamped() -> None:
    """Extreme conditions drive the refraction discriminant below zero."""
    # Huge lapse + huge pressure + very cold temperature make
    # d = 1 - 1.848*krefr*P/T^2 < 0, which is clamped to 0 (line 660),
    # yielding a zero-magnitude dip.
    d = refr.calc_dip(1000.0, lapse_rate=10.0, atpress=10_000.0, attemp=-200.0)
    assert d == 0.0


def test_calc_dip_standard_conditions() -> None:
    """A 1 km observer under standard conditions sees a negative dip."""
    d = refr.calc_dip(1000.0)
    assert d < 0.0
