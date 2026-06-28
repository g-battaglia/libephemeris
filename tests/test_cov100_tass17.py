# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for ``libephemeris.moon_theories.tass17``.

These tests target the few uncovered branches in the TASS 1.7 theory port:
the empty-longitude-series skip in ``_calc_lon``, the Hyperion branch in
``_calc_tass17_elem``, the zero-eccentricity / zero-inclination /
immediate-Kepler-convergence edge paths in ``_elliptic_to_rectangular``,
the invalid-body guards, and the velocity / Titan wrapper entry points.
"""

from __future__ import annotations

import math

import pytest

import libephemeris.moon_theories.tass17 as tass17


# Deterministic Julian Date (J2000.0, TT) reused across position tests.
J2000_TT = 2451545.0


def test_calc_lon_skips_empty_series(monkeypatch):
    """Cover the 123->119 arc: a body whose longitude series is empty.

    The real ``_BODY_SERIES`` always has a non-empty ``series[1]``, so the
    ``len(series_1) > 0`` guard never falls through. Patch one body to have
    an empty longitude series to exercise the loop-continuation arc.
    """
    original = tass17._BODY_SERIES
    patched = [list(body) for body in original]
    patched[0][1] = []  # Empty longitude series for body 0.
    monkeypatch.setattr(tass17, "_BODY_SERIES", patched)

    lon = tass17._calc_lon(0.0)

    assert len(lon) == 7
    # Body 0 had no longitude terms applied, so it stays at its init value.
    assert lon[0] == 0.0


def test_calc_tass17_elem_hyperion_branch():
    """Cover line 169: Hyperion (body 7) uses ``start_idx = 0``.

    For all other bodies the first longitude multiterm is already consumed
    by ``_calc_lon`` and ``start_idx`` is 1; Hyperion takes the else branch.
    """
    lon = tass17._calc_lon(0.0)
    elem = tass17._calc_tass17_elem(0.0, lon, 7)

    assert len(elem) == 6
    assert all(math.isfinite(v) for v in elem)


def test_elliptic_to_rectangular_zero_eccentricity_and_inclination():
    """Cover lines 238, 246->252, and 285-286 via crafted elements.

    With k = h = 0 the eccentricity is below 1e-12 so ``varpi`` is forced to
    0.0 (line 238); the Kepler iteration converges immediately and breaks
    (246->252); with q = p = 0 the inclination half-angle is below 1e-12 so
    ``cos_Omega`` / ``sin_Omega`` take their default values (lines 285-286).
    """
    # elem = [n, lambda, k, h, q, p]
    elem = [0.5, 1.0, 0.0, 0.0, 0.0, 0.0]

    result = tass17._elliptic_to_rectangular(2.0, elem, 0.0)

    assert len(result) == 6
    assert all(math.isfinite(v) for v in result)


def test_elliptic_to_rectangular_normal_kepler_convergence():
    """Exercise the Kepler break with a small but non-zero eccentricity."""
    elem = [0.5, 1.0, 0.02, 0.01, 0.01, 0.005]

    result = tass17._elliptic_to_rectangular(2.0, elem, 0.0)

    assert len(result) == 6
    assert all(math.isfinite(v) for v in result)


def test_elliptic_to_rectangular_kepler_non_convergence():
    """Cover the 246->252 arc: the Kepler loop runs all 15 iterations.

    A high eccentricity (e = 0.9) makes the fixed-point Kepler iteration
    fail to reach the 1e-14 tolerance within 15 steps, so the loop exits by
    exhaustion (the for-header -> post-loop arc) rather than via ``break``.
    """
    # k = 0.9, h = 0 -> e = 0.9 (still < 1, so beta = sqrt(1 - e^2) is valid).
    elem = [0.5, 1.0, 0.9, 0.0, 0.01, 0.005]

    result = tass17._elliptic_to_rectangular(2.0, elem, 0.0)

    assert len(result) == 6
    assert all(math.isfinite(v) for v in result)


@pytest.mark.parametrize("body", [-1, 8, 100])
def test_saturn_moon_position_invalid_body(body):
    """Cover line 354: invalid body index raises ValueError."""
    with pytest.raises(ValueError, match="Invalid body index"):
        tass17.saturn_moon_position(J2000_TT, body=body)


@pytest.mark.parametrize("body", [-1, 8, 100])
def test_saturn_moon_position_velocity_invalid_body(body):
    """Cover line 393 (within 393-415): invalid body raises ValueError."""
    with pytest.raises(ValueError, match="Invalid body index"):
        tass17.saturn_moon_position_velocity(J2000_TT, body=body)


@pytest.mark.parametrize("body", list(range(8)))
def test_saturn_moon_position_velocity_all_bodies(body):
    """Cover lines 393-415: the full position+velocity pipeline.

    Iterating all eight bodies also exercises the Hyperion branch (body 7)
    inside the elements computation reached from this entry point.
    """
    result = tass17.saturn_moon_position_velocity(J2000_TT, body=body)

    assert len(result) == 6
    assert all(math.isfinite(v) for v in result)
    # Saturn-centric satellite positions are well under 0.03 AU.
    x, y, z = result[0], result[1], result[2]
    assert math.sqrt(x * x + y * y + z * z) < 0.05


def test_titan_position_wrapper():
    """Cover line 425: the backwards-compat Titan wrapper."""
    wrapper = tass17._titan_position(J2000_TT)
    direct = tass17.saturn_moon_position(J2000_TT, body=5)

    assert wrapper == direct
