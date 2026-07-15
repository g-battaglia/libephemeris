# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Galilean satellite positions from the Lieske E5 analytical theory.

Implements the E5 ephemeris of the four Galilean satellites -- Io, Europa,
Ganymede and Callisto -- in the algebraic form tabulated by Meeus
(Astronomical Algorithms, 2nd ed., ch. 44): the full periodic series for
each satellite's longitude, latitude and radius vector are evaluated and
assembled into Jupiter-centric Cartesian vectors.

Theory sources:
    * Lieske, J.H. (1998). "Galilean satellite ephemerides E5."
      Astronomy & Astrophysics Supplement Series 129, 205-217.
    * Meeus, J. (1998). Astronomical Algorithms, 2nd ed., Willmann-Bell,
      ch. 44 ("Positions of the satellites of Jupiter").

Frame conventions and project adaptations:
    * The input Julian Date is interpreted directly as JDE (TT/TDB); no
      TT->TDB conversion and no light-time correction are applied here
      (the caller handles light-time upstream).
    * Output vectors are Jupiter-centric, in kilometres, referred to the
      J2000 ecliptic.  The published recipe precesses the B1950-referred
      E5 longitudes to the equinox of date; this module instead applies
      the fixed precession arc B1950.0 -> J2000.0, so the output frame is
      permanently J2000.
    * The orientation angles of Jupiter's equator and orbit used by the
      rotation chain are likewise held fixed at their J2000 values rather
      than being evaluated at the date.
    * Inherent accuracy of the theory in this form is ~100-200 km per
      satellite (about 0.05 arcsec at Jupiter's geocentric distance).

The computation works in degrees end to end, converting to radians only at
each trigonometric call site, and never reduces an argument modulo 360;
every output component is therefore a smooth, pure, deterministic function
of the input epoch (the caller differentiates it numerically to obtain
velocities).
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

__all__ = ["galilean_moon_positions"]


def _to_radians(degrees: float) -> float:
    """Convert degrees to radians as ``degrees * math.pi / 180.0``.

    The multiply-by-pi-then-divide form is deliberate (``math.radians``
    rounds differently in the last ulp): the series arguments reach ~1e7
    degrees at the edges of the supported time span, where that last ulp
    is visible in the output at the ~0.1 m level.

    Args:
        degrees: Angle in degrees.

    Returns:
        The angle in radians.
    """
    return degrees * math.pi / 180.0


# Reference epoch of the E5 arguments; t is measured in days from this JD.
_E5_EPOCH_JD = 2443000.5

# Jupiter's equatorial radius, used to scale the final vectors to km.
_JUPITER_RADIUS_KM = 71492.0

# Fixed precession arc from B1950.0 (JD 2433282.423) to J2000.0
# (JD 2451545.0), in degrees.  Added to the true longitudes and to psi so
# that the output frame is permanently J2000 (project adaptation; the
# published recipe precesses to the equinox of date instead).
_T0 = (2451545.0 - 2433282.423) / 36525
_J2000_PRECESSION_DEG = 1.3966626 * _T0 + 0.0003088 * (_T0 * _T0)

# Orientation of Jupiter, frozen at the J2000 values (project adaptation):
# the inclination of Jupiter's equator on its orbital plane (whose
# 0.0006 deg/century rate is referred to epoch 1900.0 = JD 2415020.5), the
# inclination of Jupiter's orbit on the ecliptic, and the longitude of the
# ascending node of Jupiter's orbit.
_EQUATOR_TILT_DEG = 3.120262 + 0.0006 * (2451545.0 - 2415020.5) / 36525
_ORBIT_TILT_DEG = 1.303267
_ORBIT_NODE_DEG = 100.464407

_COS_EQUATOR_TILT = math.cos(_to_radians(_EQUATOR_TILT_DEG))
_SIN_EQUATOR_TILT = math.sin(_to_radians(_EQUATOR_TILT_DEG))
_COS_ORBIT_TILT = math.cos(_to_radians(_ORBIT_TILT_DEG))
_SIN_ORBIT_TILT = math.sin(_to_radians(_ORBIT_TILT_DEG))
_COS_ORBIT_NODE = math.cos(_to_radians(_ORBIT_NODE_DEG))
_SIN_ORBIT_NODE = math.sin(_to_radians(_ORBIT_NODE_DEG))


@dataclass
class _E5Arguments:
    """Angular arguments of the E5 theory at one epoch, in degrees.

    Attribute names mirror the customary literature notation (Lieske 1998;
    Meeus 1998, ch. 44):

    * ``l1``..``l4``: mean longitudes of Io, Europa, Ganymede, Callisto.
    * ``pi1``..``pi4``: longitudes of the proper perijoves.
    * ``w1``..``w4``: longitudes of the proper nodes on Jupiter's
      equatorial plane (omega_1..omega_4).
    * ``psi``: longitude of the node of Jupiter's equator on the ecliptic.
    * ``phi``: phase of free libration (Phi_lambda).
    * ``g`` / ``gp``: mean anomalies of Jupiter and Saturn (G and G').
    * ``pj``: longitude of Jupiter's perihelion (Pi, constant).
    * ``sigma1``..``sigma4``: periodic longitude corrections (Sigma_i).
    * ``L1``..``L4``: true longitudes L_i = l_i + Sigma_i.

    The correction and true-longitude fields start at zero and are filled
    in once the longitude series have been evaluated; only the latitude
    series read them.
    """

    l1: float
    l2: float
    l3: float
    l4: float
    pi1: float
    pi2: float
    pi3: float
    pi4: float
    w1: float
    w2: float
    w3: float
    w4: float
    psi: float
    phi: float
    g: float
    gp: float
    pj: float
    sigma1: float = 0.0
    sigma2: float = 0.0
    sigma3: float = 0.0
    sigma4: float = 0.0
    L1: float = 0.0
    L2: float = 0.0
    L3: float = 0.0
    L4: float = 0.0


_Vector = tuple[float, float, float]

# A series term is (amplitude, argument builder); the builder returns the
# trigonometric argument in degrees.  The builders are small lambdas so
# that every argument keeps the exact multiplier structure of the printed
# tables -- e.g. "2(l1 - l2)" stays ``2.0 * (a.l1 - a.l2)`` rather than
# being distributed -- which matters for floating-point agreement with the
# recorded baseline.  Lambdas also accommodate the three composite
# latitude arguments that embed a non-integer multiple of a longitude
# correction Sigma_i, which no integer-multiplier table could express.
_SeriesTerm = tuple[float, Callable[[_E5Arguments], float]]
_Series = tuple[_SeriesTerm, ...]

# Sigma_1: periodic corrections to Io's mean longitude (degrees, 21 sine
# terms).
_IO_LONGITUDE_TERMS: _Series = (
    (0.47259, lambda a: 2.0 * (a.l1 - a.l2)),
    (-0.03478, lambda a: a.pi3 - a.pi4),
    (0.01081, lambda a: a.l2 - 2.0 * a.l3 + a.pi3),
    (0.00738, lambda a: a.phi),
    (0.00713, lambda a: a.l2 - 2.0 * a.l3 + a.pi2),
    (-0.00674, lambda a: a.pi1 + a.pi3 - 2.0 * a.pj - 2.0 * a.g),
    (0.00666, lambda a: a.l2 - 2.0 * a.l3 + a.pi4),
    (0.00445, lambda a: a.l1 - a.pi3),
    (-0.00354, lambda a: a.l1 - a.l2),
    (-0.00317, lambda a: 2.0 * a.psi - 2.0 * a.pj),
    (0.00265, lambda a: a.l1 - a.pi4),
    (-0.00186, lambda a: a.g),
    (0.00162, lambda a: a.pi2 - a.pi3),
    (0.00158, lambda a: 4.0 * (a.l1 - a.l2)),
    (-0.00155, lambda a: a.l1 - a.l3),
    (-0.00138, lambda a: a.psi + a.w3 - 2.0 * a.pj - 2.0 * a.g),
    (-0.00115, lambda a: 2.0 * (a.l1 - 2.0 * a.l2 + a.w2)),
    (0.00089, lambda a: a.pi2 - a.pi4),
    (0.00085, lambda a: a.l1 + a.pi3 - 2.0 * a.pj - 2.0 * a.g),
    (0.00083, lambda a: a.w2 - a.w3),
    (0.00053, lambda a: a.psi - a.w2),
)

# Sigma_2: periodic corrections to Europa's mean longitude (degrees, 37
# sine terms).
_EUROPA_LONGITUDE_TERMS: _Series = (
    (1.06476, lambda a: 2.0 * (a.l2 - a.l3)),
    (0.04256, lambda a: a.l1 - 2.0 * a.l2 + a.pi3),
    (0.03581, lambda a: a.l2 - a.pi3),
    (0.02395, lambda a: a.l1 - 2.0 * a.l2 + a.pi4),
    (0.01984, lambda a: a.l2 - a.pi4),
    (-0.01778, lambda a: a.phi),
    (0.01654, lambda a: a.l2 - a.pi2),
    (0.01334, lambda a: a.l2 - 2.0 * a.l3 + a.pi2),
    (0.01294, lambda a: a.pi3 - a.pi4),
    (-0.01142, lambda a: a.l2 - a.l3),
    (-0.01057, lambda a: a.g),
    (-0.00775, lambda a: 2.0 * (a.psi - a.pj)),
    (0.00524, lambda a: 2.0 * (a.l1 - a.l2)),
    (-0.00460, lambda a: a.l1 - a.l3),
    (0.00316, lambda a: a.psi - 2.0 * a.g + a.w3 - 2.0 * a.pj),
    (-0.00203, lambda a: a.pi1 + a.pi3 - 2.0 * a.pj - 2.0 * a.g),
    (0.00146, lambda a: a.psi - a.w3),
    (-0.00145, lambda a: 2.0 * a.g),
    (0.00125, lambda a: a.psi - a.w4),
    (-0.00115, lambda a: a.l1 - 2.0 * a.l3 + a.pi3),
    (-0.00094, lambda a: 2.0 * (a.l2 - a.w2)),
    (0.00086, lambda a: 2.0 * (a.l1 - 2.0 * a.l2 + a.w2)),
    (-0.00086, lambda a: 5.0 * a.gp - 2.0 * a.g + 52.225),
    (-0.00078, lambda a: a.l2 - a.l4),
    (-0.00064, lambda a: 3.0 * a.l3 - 7.0 * a.l4 + 4.0 * a.pi4),
    (0.00064, lambda a: a.pi1 - a.pi4),
    (-0.00063, lambda a: a.l1 - 2.0 * a.l3 + a.pi4),
    (0.00058, lambda a: a.w3 - a.w4),
    (0.00056, lambda a: 2.0 * (a.psi - a.pj - a.g)),
    (0.00056, lambda a: 2.0 * (a.l2 - a.l4)),
    (0.00055, lambda a: 2.0 * (a.l1 - a.l3)),
    (0.00052, lambda a: 3.0 * a.l3 - 7.0 * a.l4 + a.pi3 + 3.0 * a.pi4),
    (-0.00043, lambda a: a.l1 - a.pi3),
    (0.00041, lambda a: 5.0 * (a.l2 - a.l3)),
    (0.00041, lambda a: a.pi4 - a.pj),
    (0.00032, lambda a: a.w2 - a.w3),
    (0.00032, lambda a: 2.0 * (a.l3 - a.g - a.pj)),
)

# Sigma_3: periodic corrections to Ganymede's mean longitude (degrees, 43
# sine terms).
_GANYMEDE_LONGITUDE_TERMS: _Series = (
    (0.16490, lambda a: a.l3 - a.pi3),
    (0.09081, lambda a: a.l3 - a.pi4),
    (-0.06907, lambda a: a.l2 - a.l3),
    (0.03784, lambda a: a.pi3 - a.pi4),
    (0.01846, lambda a: 2.0 * (a.l3 - a.l4)),
    (-0.01340, lambda a: a.g),
    (-0.01014, lambda a: 2.0 * (a.psi - a.pj)),
    (0.00704, lambda a: a.l2 - 2.0 * a.l3 + a.pi3),
    (-0.00620, lambda a: a.l2 - 2.0 * a.l3 + a.pi2),
    (-0.00541, lambda a: a.l3 - a.l4),
    (0.00381, lambda a: a.l2 - 2.0 * a.l3 + a.pi4),
    (0.00235, lambda a: a.psi - a.w3),
    (0.00198, lambda a: a.psi - a.w4),
    (0.00176, lambda a: a.phi),
    (0.00130, lambda a: 3.0 * (a.l3 - a.l4)),
    (0.00125, lambda a: a.l1 - a.l3),
    (-0.00119, lambda a: 5.0 * a.gp - 2.0 * a.g + 52.225),
    (0.00109, lambda a: a.l1 - a.l2),
    (-0.00100, lambda a: 3.0 * a.l3 - 7.0 * a.l4 + 4.0 * a.pi4),
    (0.00091, lambda a: a.w3 - a.w4),
    (0.00080, lambda a: 3.0 * a.l3 - 7.0 * a.l4 + a.pi3 + 3.0 * a.pi4),
    (-0.00075, lambda a: 2.0 * a.l2 - 3.0 * a.l3 + a.pi3),
    (0.00072, lambda a: a.pi1 + a.pi3 - 2.0 * a.pj - 2.0 * a.g),
    (0.00069, lambda a: a.pi4 - a.pj),
    (-0.00058, lambda a: 2.0 * a.l3 - 3.0 * a.l4 + a.pi4),
    (-0.00057, lambda a: a.l3 - 2.0 * a.l4 + a.pi4),
    (0.00056, lambda a: a.l3 + a.pi3 - 2.0 * a.pj - 2.0 * a.g),
    (-0.00052, lambda a: a.l2 - 2.0 * a.l3 + a.pi1),
    (-0.00050, lambda a: a.pi2 - a.pi3),
    (0.00048, lambda a: a.l3 - 2.0 * a.l4 + a.pi3),
    (-0.00045, lambda a: 2.0 * a.l2 - 3.0 * a.l3 + a.pi4),
    (-0.00041, lambda a: a.pi2 - a.pi4),
    (-0.00038, lambda a: 2.0 * a.g),
    (-0.00037, lambda a: a.pi3 - a.pi4 + a.w3 - a.w4),
    (-0.00032, lambda a: 3.0 * a.l3 - 7.0 * a.l4 + 2.0 * a.pi3 + 2.0 * a.pi4),
    (0.00030, lambda a: 4.0 * (a.l3 - a.l4)),
    (0.00029, lambda a: a.l3 + a.pi4 - 2.0 * a.pj - 2.0 * a.g),
    (-0.00028, lambda a: a.w3 + a.psi - 2.0 * a.pj - 2.0 * a.g),
    (0.00026, lambda a: a.l3 - a.pj - a.g),
    (0.00024, lambda a: a.l2 - 3.0 * a.l3 + 2.0 * a.l4),
    (0.00021, lambda a: 2.0 * (a.l3 - a.pj - a.g)),
    (-0.00021, lambda a: a.l3 - a.pi2),
    (0.00017, lambda a: 2.0 * (a.l3 - a.pi3)),
)

# Sigma_4: periodic corrections to Callisto's mean longitude (degrees, 50
# sine terms).
_CALLISTO_LONGITUDE_TERMS: _Series = (
    (0.84287, lambda a: a.l4 - a.pi4),
    (0.03431, lambda a: a.pi4 - a.pi3),
    (-0.03305, lambda a: 2.0 * (a.psi - a.pj)),
    (-0.03211, lambda a: a.g),
    (-0.01862, lambda a: a.l4 - a.pi3),
    (0.01186, lambda a: a.psi - a.w4),
    (0.00623, lambda a: a.l4 + a.pi4 - 2.0 * a.g - 2.0 * a.pj),
    (0.00387, lambda a: 2.0 * (a.l4 - a.pi4)),
    (-0.00284, lambda a: 5.0 * a.gp - 2.0 * a.g + 52.225),
    (-0.00234, lambda a: 2.0 * (a.psi - a.pi4)),
    (-0.00223, lambda a: a.l3 - a.l4),
    (-0.00208, lambda a: a.l4 - a.pj),
    (0.00178, lambda a: a.psi + a.w4 - 2.0 * a.pi4),
    (0.00134, lambda a: a.pi4 - a.pj),
    (0.00125, lambda a: 2.0 * (a.l4 - a.g - a.pj)),
    (-0.00117, lambda a: 2.0 * a.g),
    (-0.00112, lambda a: 2.0 * (a.l3 - a.l4)),
    (0.00107, lambda a: 3.0 * a.l3 - 7.0 * a.l4 + 4.0 * a.pi4),
    (0.00102, lambda a: a.l4 - a.g - a.pj),
    (0.00096, lambda a: 2.0 * a.l4 - a.psi - a.w4),
    (0.00087, lambda a: 2.0 * (a.psi - a.w4)),
    (-0.00085, lambda a: 3.0 * a.l3 - 7.0 * a.l4 + a.pi3 + 3.0 * a.pi4),
    (0.00085, lambda a: a.l3 - 2.0 * a.l4 + a.pi4),
    (-0.00081, lambda a: 2.0 * (a.l4 - a.psi)),
    (0.00071, lambda a: a.l4 + a.pi4 - 2.0 * a.pj - 3.0 * a.g),
    (0.00061, lambda a: a.l1 - a.l4),
    (-0.00056, lambda a: a.psi - a.w3),
    (-0.00054, lambda a: a.l3 - 2.0 * a.l4 + a.pi3),
    (0.00051, lambda a: a.l2 - a.l4),
    (0.00042, lambda a: 2.0 * (a.psi - a.g - a.pj)),
    (0.00039, lambda a: 2.0 * (a.pi4 - a.w4)),
    (0.00036, lambda a: a.psi + a.pj - a.pi4 - a.w4),
    (0.00035, lambda a: 2.0 * a.gp - a.g + 188.37),
    (-0.00035, lambda a: a.l4 - a.pi4 + 2.0 * a.pj - 2.0 * a.psi),
    (-0.00032, lambda a: a.l4 + a.pi4 - 2.0 * a.pj - a.g),
    (0.00030, lambda a: 2.0 * a.gp - 2.0 * a.g + 149.15),
    (0.00029, lambda a: 3.0 * a.l3 - 7.0 * a.l4 + 2.0 * a.pi3 + 2.0 * a.pi4),
    (0.00028, lambda a: a.l4 - a.pi4 + 2.0 * a.psi - 2.0 * a.pj),
    (-0.00028, lambda a: 2.0 * (a.l4 - a.w4)),
    (-0.00027, lambda a: a.pi3 - a.pi4 + a.w3 - a.w4),
    (-0.00026, lambda a: 5.0 * a.gp - 3.0 * a.g + 188.37),
    (0.00025, lambda a: a.w4 - a.w3),
    (-0.00025, lambda a: a.l2 - 3.0 * a.l3 + 2.0 * a.l4),
    (-0.00023, lambda a: 3.0 * (a.l3 - a.l4)),
    (0.00021, lambda a: 2.0 * a.l4 - 2.0 * a.pj - 3.0 * a.g),
    (-0.00021, lambda a: 2.0 * a.l3 - 3.0 * a.l4 + a.pi4),
    (0.00019, lambda a: a.l4 - a.pi4 - a.g),
    (-0.00019, lambda a: 2.0 * a.l4 - a.pi3 - a.pi4),
    (-0.00018, lambda a: a.l4 - a.pi4 + a.g),
    (-0.00016, lambda a: a.l4 + a.pi3 - 2.0 * a.pj - 2.0 * a.g),
)

# tan(B_1): Io's jovicentric latitude series (dimensionless, 7 sine
# terms).  The latitude tables mix true longitudes L_i, mean longitudes
# l_i, and -- in three composite arguments -- non-integer multiples of the
# longitude corrections Sigma_i.
_IO_LATITUDE_TERMS: _Series = (
    (0.0006393, lambda a: a.L1 - a.w1),
    (0.0001825, lambda a: a.L1 - a.w2),
    (0.0000329, lambda a: a.L1 - a.w3),
    (-0.0000311, lambda a: a.L1 - a.psi),
    (0.0000093, lambda a: a.L1 - a.w4),
    (0.0000075, lambda a: 3.0 * a.L1 - 4.0 * a.l2 - 1.9927 * a.sigma1 + a.w2),
    (0.0000046, lambda a: a.L1 + a.psi - 2.0 * a.pj - 2.0 * a.g),
)

# tan(B_2): Europa's jovicentric latitude series (dimensionless, 9 sine
# terms).
_EUROPA_LATITUDE_TERMS: _Series = (
    (0.0081004, lambda a: a.L2 - a.w2),
    (0.0004512, lambda a: a.L2 - a.w3),
    (-0.0003284, lambda a: a.L2 - a.psi),
    (0.0001160, lambda a: a.L2 - a.w4),
    (0.0000272, lambda a: a.l1 - 2.0 * a.l3 + 1.0146 * a.sigma2 + a.w2),
    (-0.0000144, lambda a: a.L2 - a.w1),
    (0.0000143, lambda a: a.L2 + a.psi - 2.0 * a.pj - 2.0 * a.g),
    (0.0000035, lambda a: a.L2 - a.psi + a.g),
    (-0.0000028, lambda a: a.l1 - 2.0 * a.l3 + 1.0146 * a.sigma2 + a.w3),
)

# tan(B_3): Ganymede's jovicentric latitude series (dimensionless, 11 sine
# terms).
_GANYMEDE_LATITUDE_TERMS: _Series = (
    (0.0032402, lambda a: a.L3 - a.w3),
    (-0.0016911, lambda a: a.L3 - a.psi),
    (0.0006847, lambda a: a.L3 - a.w4),
    (-0.0002797, lambda a: a.L3 - a.w2),
    (0.0000321, lambda a: a.L3 + a.psi - 2.0 * a.pj - 2.0 * a.g),
    (0.0000051, lambda a: a.L3 - a.psi + a.g),
    (-0.0000045, lambda a: a.L3 - a.psi - a.g),
    (-0.0000045, lambda a: a.L3 + a.psi - 2.0 * a.pj),
    (0.0000037, lambda a: a.L3 + a.psi - 2.0 * a.pj - 3.0 * a.g),
    (0.0000030, lambda a: 2.0 * a.l2 - 3.0 * a.L3 + 4.03 * a.sigma3 + a.w2),
    (-0.0000021, lambda a: 2.0 * a.l2 - 3.0 * a.L3 + 4.03 * a.sigma3 + a.w3),
)

# tan(B_4): Callisto's jovicentric latitude series (dimensionless, 8 sine
# terms).
_CALLISTO_LATITUDE_TERMS: _Series = (
    (-0.0076579, lambda a: a.L4 - a.psi),
    (0.0044134, lambda a: a.L4 - a.w4),
    (-0.0005112, lambda a: a.L4 - a.w3),
    (0.0000773, lambda a: a.L4 + a.psi - 2.0 * a.pj - 2.0 * a.g),
    (0.0000104, lambda a: a.L4 - a.psi + a.g),
    (-0.0000102, lambda a: a.L4 - a.psi - a.g),
    (0.0000088, lambda a: a.L4 + a.psi - 2.0 * a.pj - 3.0 * a.g),
    (-0.0000038, lambda a: a.L4 + a.psi - 2.0 * a.pj - a.g),
)

# Delta R_1: fractional corrections to Io's radius vector (dimensionless,
# 7 cosine terms).  The radius tables use the mean longitudes only.
_IO_RADIUS_TERMS: _Series = (
    (-0.0041339, lambda a: 2.0 * (a.l1 - a.l2)),
    (-0.0000387, lambda a: a.l1 - a.pi3),
    (-0.0000214, lambda a: a.l1 - a.pi4),
    (0.0000170, lambda a: a.l1 - a.l2),
    (-0.0000131, lambda a: 4.0 * (a.l1 - a.l2)),
    (0.0000106, lambda a: a.l1 - a.l3),
    (-0.0000066, lambda a: a.l1 + a.pi3 - 2.0 * a.pj - 2.0 * a.g),
)

# Delta R_2: fractional corrections to Europa's radius vector
# (dimensionless, 11 cosine terms).
_EUROPA_RADIUS_TERMS: _Series = (
    (0.0093848, lambda a: a.l1 - a.l2),
    (-0.0003116, lambda a: a.l2 - a.pi3),
    (-0.0001744, lambda a: a.l2 - a.pi4),
    (-0.0001442, lambda a: a.l2 - a.pi2),
    (0.0000553, lambda a: a.l2 - a.l3),
    (0.0000523, lambda a: a.l1 - a.l3),
    (-0.0000290, lambda a: 2.0 * (a.l1 - a.l2)),
    (0.0000164, lambda a: 2.0 * (a.l2 - a.w2)),
    (0.0000107, lambda a: a.l1 - 2.0 * a.l3 + a.pi3),
    (-0.0000102, lambda a: a.l2 - a.pi1),
    (-0.0000091, lambda a: 2.0 * (a.l1 - a.l3)),
)

# Delta R_3: fractional corrections to Ganymede's radius vector
# (dimensionless, 10 cosine terms).
_GANYMEDE_RADIUS_TERMS: _Series = (
    (-0.0014388, lambda a: a.l3 - a.pi3),
    (-0.0007919, lambda a: a.l3 - a.pi4),
    (0.0006342, lambda a: a.l2 - a.l3),
    (-0.0001761, lambda a: 2.0 * (a.l3 - a.l4)),
    (0.0000294, lambda a: a.l3 - a.l4),
    (-0.0000156, lambda a: 3.0 * (a.l3 - a.l4)),
    (0.0000156, lambda a: a.l1 - a.l3),
    (-0.0000153, lambda a: a.l1 - a.l2),
    (0.0000070, lambda a: 2.0 * a.l2 - 3.0 * a.l3 + a.pi3),
    (-0.0000051, lambda a: a.l3 + a.pi3 - 2.0 * a.pj - 2.0 * a.g),
)

# Delta R_4: fractional corrections to Callisto's radius vector
# (dimensionless, 16 cosine terms).
_CALLISTO_RADIUS_TERMS: _Series = (
    (-0.0073546, lambda a: a.l4 - a.pi4),
    (0.0001621, lambda a: a.l4 - a.pi3),
    (0.0000974, lambda a: a.l3 - a.l4),
    (-0.0000543, lambda a: a.l4 + a.pi4 - 2.0 * a.pj - 2.0 * a.g),
    (-0.0000271, lambda a: 2.0 * (a.l4 - a.pi4)),
    (0.0000182, lambda a: a.l4 - a.pj),
    (0.0000177, lambda a: 2.0 * (a.l3 - a.l4)),
    (-0.0000167, lambda a: 2.0 * a.l4 - a.psi - a.w4),
    (0.0000167, lambda a: a.psi - a.w4),
    (-0.0000155, lambda a: 2.0 * (a.l4 - a.pj - a.g)),
    (0.0000142, lambda a: 2.0 * (a.l4 - a.psi)),
    (0.0000105, lambda a: a.l1 - a.l4),
    (0.0000092, lambda a: a.l2 - a.l4),
    (-0.0000089, lambda a: a.l4 - a.pj - a.g),
    (-0.0000062, lambda a: a.l4 + a.pi4 - 2.0 * a.pj - 3.0 * a.g),
    (0.0000048, lambda a: 2.0 * (a.l4 - a.w4)),
)


def _sum_sine_series(terms: _Series, a: _E5Arguments) -> float:
    """Accumulate ``amplitude * sin(argument)`` over a table, in order.

    Args:
        terms: Series table of (amplitude, argument builder) rows.
        a: Evaluated angular arguments, in degrees.

    Returns:
        The series sum, in the unit of the amplitudes.
    """
    total = 0.0
    for amplitude, argument in terms:
        total += amplitude * math.sin(_to_radians(argument(a)))
    return total


def _sum_cosine_series(terms: _Series, a: _E5Arguments) -> float:
    """Accumulate ``amplitude * cos(argument)`` over a table, in order.

    Args:
        terms: Series table of (amplitude, argument builder) rows.
        a: Evaluated angular arguments, in degrees.

    Returns:
        The series sum, in the unit of the amplitudes.
    """
    total = 0.0
    for amplitude, argument in terms:
        total += amplitude * math.cos(_to_radians(argument(a)))
    return total


def _fundamental_arguments(t: float) -> _E5Arguments:
    """Evaluate the fundamental angular arguments of the E5 theory.

    Args:
        t: Days elapsed from the E5 reference epoch JD 2443000.5.

    Returns:
        The argument set, in degrees, with the longitude-correction and
        true-longitude fields still zeroed.
    """
    # Fundamental arguments of Lieske's E5 theory (Lieske 1998, A&AS 129).
    # Principal inequality in Jupiter's longitude, folded into G below.
    gamma = 0.33033 * math.sin(_to_radians(163.679 + 0.0010512 * t))
    gamma += 0.03439 * math.sin(_to_radians(34.486 - 0.0161731 * t))
    return _E5Arguments(
        l1=106.07719 + 203.488955790 * t,
        l2=175.73161 + 101.374724735 * t,
        l3=120.55883 + 50.317609207 * t,
        l4=84.44459 + 21.571071177 * t,
        pi1=97.0881 + 0.16138586 * t,
        pi2=154.8663 + 0.04726307 * t,
        pi3=188.1840 + 0.00712734 * t,
        pi4=335.2868 + 0.00184000 * t,
        w1=312.3346 - 0.13279386 * t,
        w2=100.4411 - 0.03263064 * t,
        w3=119.1942 - 0.00717703 * t,
        w4=322.6186 - 0.00175934 * t,
        psi=316.5182 - 0.00000208 * t,
        phi=199.6766 + 0.17379190 * t,
        g=30.23756 + 0.0830925701 * t + gamma,
        gp=31.97853 + 0.0334597339 * t,
        pj=13.469942,
    )


def _j2000_ecliptic_km(
    longitude: float, latitude: float, radius: float, node: float
) -> _Vector:
    """Assemble one satellite vector and rotate it to the J2000 ecliptic.

    The vector is first built in the Jupiter-equatorial frame (X axis
    toward the node ``node``, Z axis toward Jupiter's north pole), then
    carried to the J2000 ecliptic through the fixed rotation chain:
    equator -> Jupiter's orbital plane -> ecliptic -> J2000 equinox.

    Args:
        longitude: Precessed true longitude L' of the satellite (degrees).
        latitude: Jovicentric latitude B of the satellite (radians).
        radius: Radius vector, in units of Jupiter's equatorial radius.
        node: Precessed longitude psi' of the node of Jupiter's equator
            on the ecliptic (degrees).

    Returns:
        The Jupiter-centric position (x, y, z) in km, J2000 ecliptic.
    """
    x = radius * math.cos(_to_radians(longitude - node)) * math.cos(latitude)
    y = radius * math.sin(_to_radians(longitude - node)) * math.cos(latitude)
    z = radius * math.sin(latitude)

    # 1. Tilt Jupiter's equator onto its orbital plane (x-rotation by I).
    y1 = y * _COS_EQUATOR_TILT - z * _SIN_EQUATOR_TILT
    z1 = y * _SIN_EQUATOR_TILT + z * _COS_EQUATOR_TILT
    # 2. Swing the node (z-rotation by psi' - Omega_j).
    swing = node - _ORBIT_NODE_DEG
    cos_swing = math.cos(_to_radians(swing))
    sin_swing = math.sin(_to_radians(swing))
    x2 = x * cos_swing - y1 * sin_swing
    y2 = x * sin_swing + y1 * cos_swing
    # 3. Tilt Jupiter's orbital plane onto the ecliptic (x-rotation by i_j).
    y3 = y2 * _COS_ORBIT_TILT - z1 * _SIN_ORBIT_TILT
    z3 = y2 * _SIN_ORBIT_TILT + z1 * _COS_ORBIT_TILT
    # 4. Rotate to the J2000 equinox (z-rotation by Omega_j).
    x4 = x2 * _COS_ORBIT_NODE - y3 * _SIN_ORBIT_NODE
    y4 = x2 * _SIN_ORBIT_NODE + y3 * _COS_ORBIT_NODE

    return (
        x4 * _JUPITER_RADIUS_KM,
        y4 * _JUPITER_RADIUS_KM,
        z3 * _JUPITER_RADIUS_KM,
    )


def galilean_moon_positions(jd: float) -> tuple[_Vector, _Vector, _Vector, _Vector]:
    """Compute Jupiter-centric positions of the four Galilean satellites.

    Evaluates the Lieske E5 theory (as tabulated by Meeus 1998, ch. 44) at
    the requested epoch and returns the satellite position vectors in the
    J2000 ecliptic frame.

    Args:
        jd: Julian Date, interpreted directly as JDE (TT/TDB).  No
            light-time correction is applied.

    Returns:
        Positions of Io, Europa, Ganymede and Callisto, in that order,
        each as an ``(x, y, z)`` tuple in kilometres, Jupiter-centric,
        referred to the J2000 ecliptic.
    """
    t = jd - _E5_EPOCH_JD
    a = _fundamental_arguments(t)

    # Periodic longitude corrections (degrees) and true longitudes.
    a.sigma1 = _sum_sine_series(_IO_LONGITUDE_TERMS, a)
    a.sigma2 = _sum_sine_series(_EUROPA_LONGITUDE_TERMS, a)
    a.sigma3 = _sum_sine_series(_GANYMEDE_LONGITUDE_TERMS, a)
    a.sigma4 = _sum_sine_series(_CALLISTO_LONGITUDE_TERMS, a)
    a.L1 = a.l1 + a.sigma1
    a.L2 = a.l2 + a.sigma2
    a.L3 = a.l3 + a.sigma3
    a.L4 = a.l4 + a.sigma4

    # Jovicentric latitudes (radians).
    b1 = math.atan(_sum_sine_series(_IO_LATITUDE_TERMS, a))
    b2 = math.atan(_sum_sine_series(_EUROPA_LATITUDE_TERMS, a))
    b3 = math.atan(_sum_sine_series(_GANYMEDE_LATITUDE_TERMS, a))
    b4 = math.atan(_sum_sine_series(_CALLISTO_LATITUDE_TERMS, a))

    # Radius vectors (units of Jupiter's equatorial radius).
    r1 = 5.90569 * (1.0 + _sum_cosine_series(_IO_RADIUS_TERMS, a))
    r2 = 9.39657 * (1.0 + _sum_cosine_series(_EUROPA_RADIUS_TERMS, a))
    r3 = 14.98832 * (1.0 + _sum_cosine_series(_GANYMEDE_RADIUS_TERMS, a))
    r4 = 26.36273 * (1.0 + _sum_cosine_series(_CALLISTO_RADIUS_TERMS, a))

    # Fixed B1950 -> J2000 precession of the true longitudes and of psi.
    # The offset cancels in each satellite longitude (L' - psi'), but psi'
    # -- not psi -- drives the node swing of the rotation chain.
    psi_prime = a.psi + _J2000_PRECESSION_DEG

    io = _j2000_ecliptic_km(a.L1 + _J2000_PRECESSION_DEG, b1, r1, psi_prime)
    europa = _j2000_ecliptic_km(a.L2 + _J2000_PRECESSION_DEG, b2, r2, psi_prime)
    ganymede = _j2000_ecliptic_km(a.L3 + _J2000_PRECESSION_DEG, b3, r3, psi_prime)
    callisto = _j2000_ecliptic_km(a.L4 + _J2000_PRECESSION_DEG, b4, r4, psi_prime)
    return (io, europa, ganymede, callisto)
