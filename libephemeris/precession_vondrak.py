"""Long-term precession (Vondrák 2011) for apparent-place reduction.

This module is the single source of the precession matrix and of-date mean
obliquity used to reduce raw ICRS positions to the equator/ecliptic of date.
It replaces the IAU 2006 precession polynomials previously used in
``fast_calc.py`` and ``astrometry.py``.

Why Vondrák 2011
----------------
The IAU 2006 precession (Capitaine et al. 2003) is a polynomial only valid for
a few centuries around J2000.0; it diverges rapidly at remote epochs (e.g. ~36"
for the Sun's longitude at year -3000). Swiss Ephemeris and modern long-term
work instead use the model of

    Vondrák, J., Capitaine, N. & Wallace, P. (2011),
    "New precession expressions, valid for long time intervals",
    Astronomy & Astrophysics 534, A22.  DOI: 10.1051/0004-6361/201117274

which is fitted to a numerical integration and stays accurate over ±200,000
years while agreeing with IAU 2006 to sub-milliarcsecond precision near J2000.
Adopting it makes ancient/future positions scientifically correct and matches
Swiss Ephemeris's precession (Swiss uses the same model since v1.78).

How it is implemented (provenance / clean-room)
-----------------------------------------------
The Vondrák 2011 model is obtained **from pyerfa/ERFA**, not transcribed by hand
and not derived from any GPL source:

    * ``erfa.ltpb``   -> long-term precession matrix, ICRS frame bias included
                        (ICRS -> mean equator and equinox of date).
    * ``erfa.ltpecl`` -> ecliptic pole unit vector (for the of-date obliquity).
    * ``erfa.ltpequ`` -> equator pole unit vector.
    * ``erfa.numat``  -> the standard nutation rotation matrix.

ERFA is the IAU SOFA-derived library distributed under a permissive BSD-3-Clause
licence (pyerfa); it is already a hard dependency of libephemeris (we use
``erfa.pmat06``/``nut06a``/``obl06`` elsewhere). Using ERFA's reference Vondrák
routines gives bit-for-bit agreement with the published model — and therefore
with Swiss Ephemeris's precession — without copying Swiss Ephemeris's GPL C code.

The nutation angles (``dpsi``, ``deps``) are NOT computed here: callers pass in
their own (LEB Chebyshev nutation on the fast path, IAU 2006/2000A on the
reference path) so that the well-tuned modern nutation is preserved unchanged.
Only the precession and the of-date *mean* obliquity become Vondrák-based.

The of-date mean obliquity is taken as the angle between the Vondrák equator
pole and ecliptic pole (``arccos(equ · ecl)``) rather than the IAU 2006
obliquity polynomial, which also diverges at remote epochs. This is the term
that closes the Sun's ~36" longitude gap at year -3000.
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
    """
    epj = _julian_epoch(jd_tt)
    p = erfa.ltpb(epj) if frame_bias else erfa.ltp(epj)
    return (
        (float(p[0][0]), float(p[0][1]), float(p[0][2])),
        (float(p[1][0]), float(p[1][1]), float(p[1][2])),
        (float(p[2][0]), float(p[2][1]), float(p[2][2])),
    )


@lru_cache(maxsize=1024)
def _mean_obliquity_rad(jd_tt: float) -> float:
    """Of-date mean obliquity (radians) = angle between equator and ecliptic poles."""
    epj = _julian_epoch(jd_tt)
    ecl = erfa.ltpecl(epj)
    equ = erfa.ltpequ(epj)
    cos_eps = float(equ[0] * ecl[0] + equ[1] * ecl[1] + equ[2] * ecl[2])
    return math.acos(max(-1.0, min(1.0, cos_eps)))


def vondrak_mean_obliquity_rad(jd_tt: float) -> float:
    """Of-date mean obliquity of the ecliptic (radians), Vondrák 2011.

    Long-term-valid replacement for the IAU 2006 obliquity polynomial.
    """
    return _mean_obliquity_rad(jd_tt)


def vondrak_mean_obliquity_deg(jd_tt: float) -> float:
    """Of-date mean obliquity of the ecliptic (degrees), Vondrák 2011."""
    return math.degrees(_mean_obliquity_rad(jd_tt))


def vondrak_precession_matrix(jd_tt: float, frame_bias: bool = True) -> Matrix3:
    """ICRS -> mean equator/equinox of date matrix (no nutation), Vondrák 2011.

    Used for mean-equator output (e.g. the sidereal/equatorial path). With
    ``frame_bias`` (default) the ICRS frame bias is included, matching the
    previous Fukushima-Williams matrix which carried the bias in its ``gamb``
    angle.
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
    return pn_mat, eps_mean_rad + deps
