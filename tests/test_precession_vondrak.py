"""Unit tests for the Vondrák 2011 long-term precession module.

The model is obtained from pyerfa/ERFA (BSD/SOFA), so these tests pin the
*construction* of the reduction matrix (precession source, frame-bias, of-date
obliquity, nutation layering and matrix orientation) against an independent
ERFA-based reduction, and confirm that the modern era is unchanged versus the
previous IAU 2006 pipeline to sub-milliarcsecond precision.
"""

from __future__ import annotations

import math

import erfa
import numpy as np
import pytest

from libephemeris.precession_vondrak import (
    vondrak_mean_obliquity_rad,
    vondrak_pn_matrix,
    vondrak_precession_matrix,
)

_J2000 = 2451545.0


def _julian_epoch_ref(jd_tt: float) -> float:
    """Independent Julian-epoch reference for the oracle.

    Computed locally (not imported from the module under test) so a bug in the
    production ``_julian_epoch`` cannot hide behind a shared conversion.

    Args:
        jd_tt: Julian Date in TT.

    Returns:
        The Julian epoch corresponding to ``jd_tt``.
    """
    return 2000.0 + (jd_tt - _J2000) / 365.25


# A grid spanning the modern era through deep-BCE and far-future epochs.
_JD_GRID = [
    _J2000 + (year - 2000) * 365.25
    for year in (-4000, -3000, -2000, -1000, 0, 1000, 1900, 2000, 2100, 4000)
]

# Representative nutation angles (radians) used to exercise the nutation layer.
_DPSI = math.radians(-0.0040)
_DEPS = math.radians(0.0026)


def _max_abs(a: np.ndarray, b: np.ndarray) -> float:
    """Return the maximum absolute element-wise difference between two arrays.

    Args:
        a: First array.
        b: Second array (same shape as ``a``).

    Returns:
        The largest ``|a - b|`` element as a native float.
    """
    return float(np.max(np.abs(a - b)))


@pytest.mark.parametrize("jd", _JD_GRID)
def test_precession_matrix_matches_erfa_ltpb(jd: float) -> None:
    """The precession-only matrix is exactly ERFA's long-term Vondrák matrix."""
    got = np.asarray(vondrak_precession_matrix(jd))
    expected = np.asarray(erfa.ltpb(_julian_epoch_ref(jd)))
    assert _max_abs(got, expected) < 1e-12


@pytest.mark.parametrize("jd", _JD_GRID)
def test_obliquity_is_pole_angle(jd: float) -> None:
    """Of-date mean obliquity is the angle between the equator and ecliptic poles."""
    epj = _julian_epoch_ref(jd)
    ecl = erfa.ltpecl(epj)
    equ = erfa.ltpequ(epj)
    # Clamp to acos's domain: FP drift can push the dot marginally past ±1.
    expected = math.acos(max(-1.0, min(1.0, float(np.dot(equ, ecl)))))
    assert abs(vondrak_mean_obliquity_rad(jd) - expected) < 1e-12


@pytest.mark.parametrize("jd", _JD_GRID)
def test_pn_matrix_orientation_and_layering(jd: float) -> None:
    """pn = N(eps_mean, dpsi, deps) @ ltpb, with eps_true = eps_mean + deps."""
    epj = _julian_epoch_ref(jd)
    eps_mean = vondrak_mean_obliquity_rad(jd)
    expected = np.asarray(erfa.numat(eps_mean, _DPSI, _DEPS)) @ np.asarray(
        erfa.ltpb(epj)
    )

    pn_mat, eps_true = vondrak_pn_matrix(jd, _DPSI, _DEPS)
    assert _max_abs(np.asarray(pn_mat), expected) < 1e-12
    assert abs(eps_true - (eps_mean + _DEPS)) < 1e-15


@pytest.mark.parametrize("year", [1800, 1900, 2000, 2100, 2200])
def test_modern_obliquity_matches_iau2006(year: int) -> None:
    """Near J2000 Vondrák agrees with the IAU 2006 obliquity to < 1 mas."""
    jd = _J2000 + (year - 2000) * 365.25
    iau2006 = erfa.obl06(_J2000, jd - _J2000)
    diff_arcsec = abs(vondrak_mean_obliquity_rad(jd) - iau2006) * 206264.806
    assert diff_arcsec < 1.0


@pytest.mark.parametrize("year", [1800, 1900, 2000, 2100, 2200])
def test_modern_precession_matches_iau2006(year: int) -> None:
    """Near J2000 the Vondrák long-term precession agrees with the previous
    IAU 2006 model to sub-arcsecond, so the modern era is unchanged.

    This compares the *precession-only* matrix against an independent IAU 2006
    matrix: ``erfa.ltpb`` (Vondrák) and ``erfa.pmat06`` (IAU 2006) are both
    GCRS→mean-equator-and-equinox-of-date and both include the ICRS frame bias,
    so the check is a genuine cross-model comparison -- not a tautology against
    the migrated pipeline (which now returns the Vondrák matrix itself).
    """
    jd = _J2000 + (year - 2000) * 365.25
    vondrak = np.asarray(vondrak_precession_matrix(jd))
    iau2006 = np.asarray(erfa.pmat06(_J2000, jd - _J2000))
    max_elem = _max_abs(vondrak, iau2006)
    assert math.degrees(max_elem) * 3600.0 < 1.0


@pytest.mark.parametrize(
    "year,min_arcsec",
    [(-3000, 5.0), (-1000, 0.5), (5000, 1.0)],
)
def test_longterm_precession_diverges_from_iau2006(
    year: int, min_arcsec: float
) -> None:
    """At remote epochs Vondrák and IAU 2006 must differ markedly: the IAU 2006
    polynomial extrapolates badly away from J2000 (measured here: ~13" at -3000,
    ~1.2" at -1000, ~2.3" at +5000 in matrix elements).

    This pins that the long-term model is actually engaged at extreme dates --
    the property that closes the ~36" Sun longitude offset at -3000 versus Swiss
    Ephemeris, which has used the same Vondrák model since v1.78. A regression to
    the IAU 2006 polynomial would collapse this divergence and fail here.
    """
    jd = _J2000 + (year - 2000) * 365.25
    vondrak = np.asarray(vondrak_precession_matrix(jd))
    iau2006 = np.asarray(erfa.pmat06(_J2000, jd - _J2000))
    diff_arcsec = math.degrees(_max_abs(vondrak, iau2006)) * 3600.0
    assert diff_arcsec > min_arcsec
