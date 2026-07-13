# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Long-term precession (Vondrák 2011) for apparent-place reduction.

This module is the single source of the precession matrix and of-date mean
obliquity used to reduce raw ICRS positions to the equator/ecliptic of date.
It replaces the IAU 2006 precession polynomials previously used in
``fast_calc.py`` and ``astrometry.py``.

Why Vondrák 2011
----------------
The IAU 2006 precession (Capitaine et al. 2003) is a polynomial only valid for
a few centuries around J2000.0; it diverges rapidly at remote epochs (e.g. ~36"
for the Sun's longitude at year -3000). Long-term ephemeris work instead uses
the model of

    Vondrák, J., Capitaine, N. & Wallace, P. (2011),
    "New precession expressions, valid for long time intervals",
    Astronomy & Astrophysics 534, A22.  DOI: 10.1051/0004-6361/201117274

which is fitted to a numerical integration and stays accurate over ±200,000
years while agreeing with IAU 2006 to sub-milliarcsecond precision near J2000.
Adopting it makes ancient/future positions scientifically correct.

How it is implemented
---------------------
The Vondrák 2011 precession matrix is obtained **from pyerfa/ERFA**, not
transcribed by hand and not derived from any third-party application source:

    * ``erfa.ltpb``   -> long-term precession matrix, ICRS frame bias included
                        (ICRS -> mean equator and equinox of date).
    * ``erfa.ltp``    -> same without the ICRS frame bias.
    * ``erfa.numat``  -> the standard nutation rotation matrix.

ERFA is the IAU SOFA-derived library distributed under a permissive BSD-3-Clause
licence (pyerfa); it is already a hard dependency of libephemeris (we use
``erfa.pmat06``/``nut06a``/``obl06`` elsewhere). Using ERFA's reference Vondrák
routines gives bit-for-bit agreement with the published model.

The nutation angles (``dpsi``, ``deps``) are NOT computed here: callers pass in
their own (LEB Chebyshev nutation on the fast path, IAU 2006/2000A on the
reference path) so that the well-tuned modern nutation is preserved unchanged.
Only the precession and the of-date *mean* obliquity become Vondrák-based.

The of-date mean obliquity is delegated to
:mod:`libephemeris.sidereal_longterm`, which evaluates the Vondrák 2011 obliquity
series directly. This is a single shared realization used by both the position
pipeline and the house cusps, and it replaces the IAU 2006 obliquity polynomial
(which diverges at remote epochs) — the term that closes the Sun's ~36"
longitude gap at year -3000.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Tuple

import erfa
import numpy as np

# J2000.0 epoch as a Julian Date (TT).
_J2000 = 2451545.0

# A 3x3 rotation matrix in the nested-tuple form used throughout fast_calc.py.
Matrix3 = Tuple[
    Tuple[float, float, float],
    Tuple[float, float, float],
    Tuple[float, float, float],
]


def _julian_epoch(jd_tt: float) -> float:
    """Convert a Julian Date (TT) to a Julian epoch (e.g. 2000.0 at J2000.0).

    ERFA's long-term precession routines are parameterised by the Julian epoch.

    Args:
        jd_tt: Julian Date in TT.

    Returns:
        The Julian epoch (Julian years) corresponding to ``jd_tt``.
    """
    return 2000.0 + (jd_tt - _J2000) / 365.25


@lru_cache(maxsize=1024)
def _precession_matrix(jd_tt: float, frame_bias: bool) -> Matrix3:
    """Vondrák precession matrix for ``jd_tt`` (memoised per epoch + bias flag).

    With ``frame_bias`` true (the default everywhere except the FLG_ICRS paths),
    the matrix maps ICRS -> mean equator/equinox of date and includes the ICRS
    frame bias (``erfa.ltpb``). With it false, the bias is omitted (``erfa.ltp``,
    J2000 mean equator -> date), matching the library's FLG_ICRS convention of
    treating ICRS as the J2000 mean frame.

    Args:
        jd_tt: Julian Date in TT.
        frame_bias: Include the ICRS frame bias (``erfa.ltpb``) when true; omit
            it (``erfa.ltp``) when false.

    Returns:
        The 3x3 precession rotation matrix as nested tuples of floats.
    """
    epj = _julian_epoch(jd_tt)
    p = erfa.ltpb(epj) if frame_bias else erfa.ltp(epj)
    return (
        (float(p[0][0]), float(p[0][1]), float(p[0][2])),
        (float(p[1][0]), float(p[1][1]), float(p[1][2])),
        (float(p[2][0]), float(p[2][1]), float(p[2][2])),
    )


def _mean_obliquity_rad(jd_tt: float) -> float:
    """Of-date mean obliquity (radians), Vondrák 2011 (long-term valid).

    Delegates to :func:`libephemeris.sidereal_longterm.mean_obliquity_rad`, which
    evaluates the Vondrák 2011 obliquity series directly. This is the single
    of-date mean-obliquity realization shared by the position pipeline and the
    house cusps, so one chart's bodies and angles sit in a consistent frame.

    Args:
        jd_tt: Julian Date in TT.

    Returns:
        The of-date mean obliquity in radians.
    """
    from . import sidereal_longterm

    return sidereal_longterm.mean_obliquity_rad(jd_tt)


def vondrak_mean_obliquity_rad(jd_tt: float) -> float:
    """Of-date mean obliquity of the ecliptic (radians), Vondrák 2011.

    Long-term-valid replacement for the IAU 2006 obliquity polynomial.

    Args:
        jd_tt: Julian Date in TT.

    Returns:
        The of-date mean obliquity in radians.
    """
    return _mean_obliquity_rad(jd_tt)


def vondrak_mean_obliquity_deg(jd_tt: float) -> float:
    """Of-date mean obliquity of the ecliptic (degrees), Vondrák 2011.

    Args:
        jd_tt: Julian Date in TT.

    Returns:
        The of-date mean obliquity in degrees.
    """
    return math.degrees(_mean_obliquity_rad(jd_tt))


def vondrak_precession_matrix(jd_tt: float, frame_bias: bool = True) -> Matrix3:
    """ICRS -> mean equator/equinox of date matrix (no nutation), Vondrák 2011.

    Used for mean-equator output (e.g. the sidereal/equatorial path). With
    ``frame_bias`` (default) the ICRS frame bias is included, matching the
    previous Fukushima-Williams matrix which carried the bias in its ``gamb``
    angle.

    Args:
        jd_tt: Julian Date in TT.
        frame_bias: Include the ICRS frame bias when true (default).

    Returns:
        The 3x3 precession rotation matrix as nested tuples of floats.
    """
    return _precession_matrix(jd_tt, frame_bias)


def vondrak_pn_matrix(
    jd_tt: float, dpsi: float, deps: float, frame_bias: bool = True
) -> Tuple[Matrix3, float]:
    """Full ICRS -> true-equator-of-date reduction matrix and true obliquity.

    Combines the Vondrák 2011 precession (``erfa.ltpb``/``erfa.ltp``) with the
    caller-supplied nutation angles via the standard nutation matrix
    (``erfa.numat``):

        ``pn = N(eps_mean, dpsi, deps) @ P``

    so that a raw ICRS Cartesian vector ``v`` maps to true-equator-of-date by
    ``pn @ v``. The ecliptic of date is then obtained by rotating about the
    equinox by ``eps_true``.

    Args:
        jd_tt: Julian Date in TT.
        dpsi: Nutation in longitude (radians), from the caller's nutation model.
        deps: Nutation in obliquity (radians), from the caller's nutation model.
        frame_bias: Include the ICRS frame bias (default). Pass false for the
            FLG_ICRS paths, which omit the bias.

    Returns:
        Tuple ``(pn_mat, eps_true_rad)`` where ``pn_mat`` is the nested-tuple
        rotation matrix and ``eps_true_rad = eps_mean + deps``.
    """
    eps_mean_rad = _mean_obliquity_rad(jd_tt)

    # Standard nutation matrix (mean -> true equator of date).
    n = erfa.numat(eps_mean_rad, dpsi, deps)
    p = np.asarray(_precession_matrix(jd_tt, frame_bias))
    pn = n @ p

    pn_mat: Matrix3 = (
        (float(pn[0][0]), float(pn[0][1]), float(pn[0][2])),
        (float(pn[1][0]), float(pn[1][1]), float(pn[1][2])),
        (float(pn[2][0]), float(pn[2][1]), float(pn[2][2])),
    )
    return pn_mat, float(eps_mean_rad + deps)


# =============================================================================
# Method-B ayanamsha precession (long-term, Vondrák ecliptic frame)
# =============================================================================
#
# A caller-defined mean ayanamsha is propagated as
#
#     ayanamsha(t) = ayan_t0 - lambda(t; t0)
#
# where lambda(t; t0) is the ecliptic longitude of the vernal equinox of date
# t, measured on the MEAN ECLIPTIC OF t0 (the epoch that defines the sidereal
# zero point). Both the ecliptic of t0 and the equinox of t are built from the
# Vondrák 2011 long-term ecliptic/equator poles, so the construction follows
# that model over its stated span. This is "Method B": it is *non-additive* in
# t0 because the ecliptic-tilt term means accum(t; t0) is not generally equal to
# accum(t; J2000) - accum(t0; J2000).


@lru_cache(maxsize=4096)
def _ltp_ecliptic_frame(jd_tt: float) -> Tuple[Tuple[float, float, float], ...]:
    """Vondrák 2011 mean ecliptic frame of date (rows in the GCRS basis).

    Built from the ERFA long-term ecliptic pole (``ltpecl``) and equator pole
    (``ltpequ``):

    * ``x`` — vernal equinox of date, the (normalised) cross product of the
      ecliptic and equator poles (ascending node of the equator on the
      ecliptic).
    * ``z`` — ecliptic pole of date.
    * ``y = z x x`` completes the right-handed triad.

    Args:
        jd_tt: Julian Date in TT.

    Returns:
        A 4-tuple ``(row_x, row_y, row_z, equinox)`` of 3-tuples. The first
        three rows map a GCRS unit vector into ecliptic-of-date coordinates;
        ``equinox`` is the vernal-equinox direction of date in GCRS (== row_x).
    """
    epj = _julian_epoch(jd_tt)
    pecl = np.asarray(erfa.ltpecl(epj))
    pequ = np.asarray(erfa.ltpequ(epj))
    x = np.cross(pecl, pequ)
    x = x / np.linalg.norm(x)
    z = pecl
    y = np.cross(z, x)
    xt = (float(x[0]), float(x[1]), float(x[2]))
    return (
        xt,
        (float(y[0]), float(y[1]), float(y[2])),
        (float(z[0]), float(z[1]), float(z[2])),
        xt,
    )


def method_b_accumulated_precession(jd_tt: float, t0_tt: float) -> float:
    """Accumulated precession in ecliptic longitude from ``t0`` to ``jd`` (deg).

    The sidereal reference direction is fixed on the **mean ecliptic of
    ``t0``**; the returned value is the amount by which the vernal equinox of
    ``jd`` has regressed along that fixed ecliptic. Concretely the equinox of
    ``jd`` is projected into the ecliptic-of-``t0`` frame, giving its longitude
    ``lambda``; the accumulated precession is ``-lambda``.

    Both frames use the Vondrák 2011 long-term ecliptic/equator poles, so the
    result is valid over the full Vondrák span (±200,000 yr) and is
    **non-additive** in ``t0``.

    Args:
        jd_tt: Epoch of interest, Julian Date in TT.
        t0_tt: Reference epoch defining the fixed ecliptic, Julian Date in TT.

    Returns:
        Accumulated precession in longitude, in degrees. NOT reduced mod 360,
        so callers can difference two epochs without a 0/360 wrap artefact.
    """
    e0 = _ltp_ecliptic_frame(t0_tt)
    equinox_t = _ltp_ecliptic_frame(jd_tt)[3]
    wx = e0[0][0] * equinox_t[0] + e0[0][1] * equinox_t[1] + e0[0][2] * equinox_t[2]
    wy = e0[1][0] * equinox_t[0] + e0[1][1] * equinox_t[1] + e0[1][2] * equinox_t[2]
    lam = math.degrees(math.atan2(wy, wx))
    return -lam
