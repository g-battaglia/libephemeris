# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for :mod:`libephemeris.astrometry`.

These tests exercise both the optional-pyerfa acceleration path and the
pure-numpy analytical fallback path for precession, nutation, obliquity and
aberration helpers. The pyerfa import-time guard (``except ImportError``) is
covered by exec-ing the module source into a throwaway namespace with ``erfa``
mocked absent, without reloading the live module.
"""

from __future__ import annotations

import importlib.util
import math
import sys

import numpy as np
import pytest

import libephemeris.astrometry as astrometry


# ---------------------------------------------------------------------------
# Pure-logic helpers (no erfa branch)
# ---------------------------------------------------------------------------


def test_normalize_angle_negative_wraps_to_positive():
    """Negative angle gets +360 added (line 72)."""
    assert astrometry._normalize_angle(-90.0) == pytest.approx(270.0)


def test_normalize_angle_already_in_range():
    """Non-negative angle is returned unchanged after modulo."""
    assert astrometry._normalize_angle(123.0) == pytest.approx(123.0)


def test_cartesian_to_spherical_zero_vector_returns_origin():
    """A near-zero vector short-circuits to (0, 0, 0) (lines 78-80)."""
    assert astrometry._cartesian_to_spherical(0.0, 0.0, 0.0) == (0.0, 0.0, 0.0)


def test_cartesian_to_spherical_unit_axes():
    """A non-degenerate vector exercises the full conversion (lines 81-84)."""
    lon, lat, dist = astrometry._cartesian_to_spherical(1.0, 0.0, 0.0)
    assert lon == pytest.approx(0.0)
    assert lat == pytest.approx(0.0)
    assert dist == pytest.approx(1.0)

    # +Y axis -> lon 90 deg; +Z axis pushes latitude positive after normalize.
    lon_y, lat_y, _ = astrometry._cartesian_to_spherical(0.0, 1.0, 0.0)
    assert lon_y == pytest.approx(90.0)
    assert lat_y == pytest.approx(0.0)

    # Negative-X gives a negative atan2 longitude that must be normalized (>0).
    lon_neg, _, _ = astrometry._cartesian_to_spherical(-1.0, -1.0, 0.0)
    assert lon_neg == pytest.approx(225.0)


# ---------------------------------------------------------------------------
# Import-time pyerfa guard: exec source with erfa mocked absent (lines 99-100)
# ---------------------------------------------------------------------------


def test_import_without_erfa_sets_has_erfa_false(monkeypatch):
    """Exec the module source with ``erfa`` absent to hit the ImportError pass.

    A throwaway module object is used so the live ``libephemeris.astrometry``
    references held by other tests are never disturbed.
    """
    # Force ``import erfa`` to raise ImportError during exec.
    monkeypatch.setitem(sys.modules, "erfa", None)

    spec = importlib.util.spec_from_file_location(
        "_astrometry_no_erfa", astrometry.__file__
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    assert module._HAS_ERFA is False
    assert module._erfa is None


# ---------------------------------------------------------------------------
# numpy fallback path: force _HAS_ERFA off so analytical branches run
# ---------------------------------------------------------------------------


@pytest.fixture
def no_erfa(monkeypatch):
    """Disable the pyerfa acceleration so analytical fallbacks execute."""
    monkeypatch.setattr(astrometry, "_HAS_ERFA", False)
    monkeypatch.setattr(astrometry, "_erfa", None)


JD = 2451545.0 + 365250.0  # ~J2000 + 1000 years TT, well-defined epoch.


def test_mean_obliquity_numpy_fallback(no_erfa):
    """Analytical obliquity polynomial path (lines 113-117)."""
    eps = astrometry._mean_obliquity(JD)
    # Mean obliquity stays near 23.4 degrees within a few millennia.
    assert 23.0 < eps < 24.0


def test_fundamental_arguments_numpy(no_erfa):
    """Delaunay arguments returned in [0, 2*pi) (lines 127-152)."""
    t = astrometry._jd_to_julian_centuries(JD)
    args = astrometry._fundamental_arguments(t)
    assert len(args) == 5
    for a in args:
        assert 0.0 <= a < astrometry.TWO_PI


def test_nutation_angles_numpy_series(no_erfa):
    """IAU 2000B 77-term numpy series (lines 238-253, 268)."""
    dpsi, deps = astrometry.nutation_angles(JD)
    # Nutation magnitude is tens of arcseconds -> small fraction of a degree.
    assert abs(dpsi) < 0.01
    assert abs(deps) < 0.01

    # Direct call to the numpy implementation (lines 238-253).
    dpsi2, deps2 = astrometry._nutation_angles_numpy(JD)
    assert dpsi2 == pytest.approx(dpsi)
    assert deps2 == pytest.approx(deps)


def test_precession_matrix_numpy_fallback(no_erfa):
    """Analytical precession rotation matrix (lines 287-314)."""
    mat = astrometry._precession_matrix_j2000_to_date(JD)
    assert mat.shape == (3, 3)
    # Rotation matrix must be orthonormal.
    np.testing.assert_allclose(mat @ mat.T, np.eye(3), atol=1e-9)


def test_nutation_matrix_numpy_fallback(no_erfa):
    """Analytical nutation rotation matrix (lines 324-349)."""
    mat = astrometry._nutation_matrix(JD)
    assert mat.shape == (3, 3)
    np.testing.assert_allclose(mat @ mat.T, np.eye(3), atol=1e-9)


def test_precession_nutation_matrix_numpy_fallback(no_erfa):
    """Combined matrix via N @ P numpy path (lines 356-358)."""
    combined = astrometry._precession_nutation_matrix(JD)
    P = astrometry._precession_matrix_j2000_to_date(JD)
    N = astrometry._nutation_matrix(JD)
    np.testing.assert_allclose(combined, N @ P, atol=1e-12)


def test_true_obliquity_numpy(no_erfa):
    """True obliquity = mean + nutation in obliquity (lines 363-365)."""
    eps_mean = astrometry._mean_obliquity(JD)
    _, d_eps = astrometry.nutation_angles(JD)
    assert astrometry._true_obliquity(JD) == pytest.approx(eps_mean + d_eps)


# ---------------------------------------------------------------------------
# erfa acceleration path (lines 109-111, 265-267, 284-285, 319-322, 354-355)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not astrometry._HAS_ERFA, reason="pyerfa not installed")
def test_erfa_paths_are_exercised():
    """When pyerfa is present, the accelerated branches run by default."""
    eps = astrometry._mean_obliquity(JD)
    assert 23.0 < eps < 24.0

    dpsi, deps = astrometry.nutation_angles(JD)
    assert abs(dpsi) < 0.01
    assert abs(deps) < 0.01

    pmat = astrometry._precession_matrix_j2000_to_date(JD)
    assert pmat.shape == (3, 3)

    nmat = astrometry._nutation_matrix(JD)
    assert nmat.shape == (3, 3)

    pnmat = astrometry._precession_nutation_matrix(JD)
    assert pnmat.shape == (3, 3)


# ---------------------------------------------------------------------------
# Aberration zero-magnitude guards (lines 477, 511, 522)
# ---------------------------------------------------------------------------


def test_annual_aberration_cartesian_zero_direction():
    """A near-zero target direction short-circuits to no shift (line 477)."""
    result = astrometry._annual_aberration_cartesian(
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)
    )
    assert result == (0.0, 0.0, 0.0)


def test_annual_aberration_cartesian_nonzero():
    """A normal direction returns a small transverse velocity shift."""
    dx, dy, dz = astrometry._annual_aberration_cartesian(
        (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)
    )
    # Velocity perpendicular to line of sight survives; radial part removed.
    assert dx == pytest.approx(0.0, abs=1e-12)
    assert dy != 0.0


def test_apply_aberration_zero_position_returns_input():
    """Zero-distance position returns the original tuple (lines 509-511)."""
    pos = (0.0, 0.0, 0.0)
    assert astrometry.apply_aberration_to_position(pos, (1.0, 0.0, 0.0)) is pos


def test_apply_aberration_zero_after_correction_returns_input(monkeypatch):
    """Aberration cancelling the unit vector hits the a_mag guard (line 522).

    Patch the Cartesian helper to return exactly the negated unit direction so
    that the aberrated vector collapses to zero magnitude.
    """

    def fake_aberration(target_direction, earth_velocity):
        return tuple(-c for c in target_direction)

    monkeypatch.setattr(astrometry, "_annual_aberration_cartesian", fake_aberration)

    pos = (1.0, 0.0, 0.0)
    assert astrometry.apply_aberration_to_position(pos, (0.0, 0.0, 0.0)) is pos


def test_apply_aberration_normal_position():
    """A real position is aberrated and keeps its distance (lines 513-527)."""
    pos = (1.0, 0.0, 0.0)
    out = astrometry.apply_aberration_to_position(pos, (0.0, 5.0, 0.0))
    # Magnitude (distance) is preserved by construction.
    assert math.sqrt(sum(c * c for c in out)) == pytest.approx(1.0)
    # Position is shifted toward the velocity direction (+Y).
    assert out[1] > 0.0
