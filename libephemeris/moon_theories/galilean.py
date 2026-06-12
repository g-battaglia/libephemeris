"""
Galilean satellite ephemeris based on the E5 theory of Lieske.

This is a port of the high-accuracy method from Jean Meeus's
"Astronomical Algorithms" Chapter 44, which implements Lieske's E5 theory.

The theory provides positions of Jupiter's four Galilean moons:
- Io (satellite I)
- Europa (satellite II)
- Ganymede (satellite III)
- Callisto (satellite IV)

Original PyMeeus implementation by Michael Lutz, Sophie Scholz,
Vittorio Serra, Sebastian Veigl (LGPL-3.0).

Reference:
- Lieske, J.H. (1998) "Galilean satellite ephemerides E5", A&AS 129, 205
- Meeus, J. (1998) "Astronomical Algorithms", 2nd edition, Chapter 44

Precision: ~100-200 km (sufficient for ~0.05 arcsec at Jupiter's distance)
"""

import math
from typing import Tuple

# =============================================================================
# CONSTANTS
# =============================================================================

# Jupiter equatorial radius in km
JUPITER_RADIUS_KM: float = 71492.0

# Mean distances of the satellites in Jupiter radii
A1: float = 5.90569  # Io
A2: float = 9.39657  # Europa
A3: float = 14.98832  # Ganymede
A4: float = 26.36273  # Callisto

# Reference epoch for the theory
# t is measured from JDE 2443000.5
EPOCH_JDE: float = 2443000.5

# AU in km for coordinate conversion
AU_KM: float = 149597870.7

# Obliquity of Jupiter's equator to the ecliptic (degrees)
# Used for transforming from Jupiter equatorial to ecliptic frame
JUPITER_EQUATOR_INCLINATION_DEG: float = 3.12


def _radians(deg: float) -> float:
    """Convert degrees to radians."""
    return deg * math.pi / 180.0


def galilean_moon_positions(
    jd: float,
) -> Tuple[
    Tuple[float, float, float],
    Tuple[float, float, float],
    Tuple[float, float, float],
    Tuple[float, float, float],
]:
    """
    Calculate positions of Jupiter's Galilean moons.

    Uses the E5 theory of Lieske as implemented by Meeus.

    Args:
        jd: Julian Date (TT/TDB)

    Returns:
        Tuple of four (x, y, z) tuples in km, for Io, Europa, Ganymede, Callisto.
        Coordinates are in the J2000 ecliptic frame, relative to Jupiter center.
    """
    # Time since JDE 2443000.5 in days
    t = jd - EPOCH_JDE

    # Mean longitudes of the satellites (degrees)
    l1 = 106.07719 + 203.488955790 * t
    l2 = 175.73161 + 101.374724735 * t
    l3 = 120.55883 + 50.317609207 * t
    l4 = 84.44459 + 21.571071177 * t

    # Longitudes of the perijoves (degrees)
    pi1 = 97.0881 + 0.16138586 * t
    pi2 = 154.8663 + 0.04726307 * t
    pi3 = 188.1840 + 0.00712734 * t
    pi4 = 335.2868 + 0.00184000 * t

    # Longitudes of the nodes on the equatorial plane of Jupiter (degrees)
    omega1 = 312.3346 - 0.13279386 * t
    omega2 = 100.4411 - 0.03263064 * t
    omega3 = 119.1942 - 0.00717703 * t
    omega4 = 322.6186 - 0.00175934 * t

    # Principal inequality in the longitude of Jupiter (degrees)
    GAMMA = 0.33033 * math.sin(_radians(163.679 + 0.0010512 * t)) + 0.03439 * math.sin(
        _radians(34.486 - 0.0161731 * t)
    )

    # Libration for inner satellites (degrees)
    PHI_lambda = 199.6766 + 0.17379190 * t

    # Longitude of the node of the equator of Jupiter on the ecliptic (degrees)
    psi = 316.5182 - 0.00000208 * t

    # Mean anomalies of Jupiter and Saturn (degrees)
    G = 30.23756 + 0.0830925701 * t + GAMMA
    G_ = 31.97853 + 0.0334597339 * t

    # Longitude of the perihelion of Jupiter (degrees)
    PI = 13.469942

    # ==========================================================================
    # PERIODIC TERMS IN LONGITUDE
    # ==========================================================================

    # Satellite 1 (Io)
    sum1 = 0.0
    sum1 += 0.47259 * math.sin(_radians(2 * (l1 - l2)))
    sum1 -= 0.03478 * math.sin(_radians(pi3 - pi4))
    sum1 += 0.01081 * math.sin(_radians(l2 - 2 * l3 + pi3))
    sum1 += 0.00738 * math.sin(_radians(PHI_lambda))
    sum1 += 0.00713 * math.sin(_radians(l2 - 2 * l3 + pi2))
    sum1 -= 0.00674 * math.sin(_radians(pi1 + pi3 - 2 * PI - 2 * G))
    sum1 += 0.00666 * math.sin(_radians(l2 - 2 * l3 + pi4))
    sum1 += 0.00445 * math.sin(_radians(l1 - pi3))
    sum1 -= 0.00354 * math.sin(_radians(l1 - l2))
    sum1 -= 0.00317 * math.sin(_radians(2 * psi - 2 * PI))
    sum1 += 0.00265 * math.sin(_radians(l1 - pi4))
    sum1 -= 0.00186 * math.sin(_radians(G))
    sum1 += 0.00162 * math.sin(_radians(pi2 - pi3))
    sum1 += 0.00158 * math.sin(_radians(4 * (l1 - l2)))
    sum1 -= 0.00155 * math.sin(_radians(l1 - l3))
    sum1 -= 0.00138 * math.sin(_radians(psi + omega3 - 2 * PI - 2 * G))
    sum1 -= 0.00115 * math.sin(_radians(2 * (l1 - 2 * l2 + omega2)))
    sum1 += 0.00089 * math.sin(_radians(pi2 - pi4))
    sum1 += 0.00085 * math.sin(_radians(l1 + pi3 - 2 * PI - 2 * G))
    sum1 += 0.00083 * math.sin(_radians(omega2 - omega3))
    sum1 += 0.00053 * math.sin(_radians(psi - omega2))

    # Satellite 2 (Europa)
    sum2 = 0.0
    sum2 += 1.06476 * math.sin(_radians(2 * (l2 - l3)))
    sum2 += 0.04256 * math.sin(_radians(l1 - 2 * l2 + pi3))
    sum2 += 0.03581 * math.sin(_radians(l2 - pi3))
    sum2 += 0.02395 * math.sin(_radians(l1 - 2 * l2 + pi4))
    sum2 += 0.01984 * math.sin(_radians(l2 - pi4))
    sum2 -= 0.01778 * math.sin(_radians(PHI_lambda))
    sum2 += 0.01654 * math.sin(_radians(l2 - pi2))
    sum2 += 0.01334 * math.sin(_radians(l2 - 2 * l3 + pi2))
    sum2 += 0.01294 * math.sin(_radians(pi3 - pi4))
    sum2 -= 0.01142 * math.sin(_radians(l2 - l3))
    sum2 -= 0.01057 * math.sin(_radians(G))
    sum2 -= 0.00775 * math.sin(_radians(2 * (psi - PI)))
    sum2 += 0.00524 * math.sin(_radians(2 * (l1 - l2)))
    sum2 -= 0.00460 * math.sin(_radians(l1 - l3))
    sum2 += 0.00316 * math.sin(_radians(psi - 2 * G + omega3 - 2 * PI))
    sum2 -= 0.00203 * math.sin(_radians(pi1 + pi3 - 2 * PI - 2 * G))
    sum2 += 0.00146 * math.sin(_radians(psi - omega3))
    sum2 -= 0.00145 * math.sin(_radians(2 * G))
    sum2 += 0.00125 * math.sin(_radians(psi - omega4))
    sum2 -= 0.00115 * math.sin(_radians(l1 - 2 * l3 + pi3))
    sum2 -= 0.00094 * math.sin(_radians(2 * (l2 - omega2)))
    sum2 += 0.00086 * math.sin(_radians(2 * (l1 - 2 * l2 + omega2)))
    sum2 -= 0.00086 * math.sin(_radians(5 * G_ - 2 * G + 52.225))
    sum2 -= 0.00078 * math.sin(_radians(l2 - l4))
    sum2 -= 0.00064 * math.sin(_radians(3 * l3 - 7 * l4 + 4 * pi4))
    sum2 += 0.00064 * math.sin(_radians(pi1 - pi4))
    sum2 -= 0.00063 * math.sin(_radians(l1 - 2 * l3 + pi4))
    sum2 += 0.00058 * math.sin(_radians(omega3 - omega4))
    sum2 += 0.00056 * math.sin(_radians(2 * (psi - PI - G)))
    sum2 += 0.00056 * math.sin(_radians(2 * (l2 - l4)))
    sum2 += 0.00055 * math.sin(_radians(2 * (l1 - l3)))
    sum2 += 0.00052 * math.sin(_radians(3 * l3 - 7 * l4 + pi3 + 3 * pi4))
    sum2 -= 0.00043 * math.sin(_radians(l1 - pi3))
    sum2 += 0.00041 * math.sin(_radians(5 * (l2 - l3)))
    sum2 += 0.00041 * math.sin(_radians(pi4 - PI))
    sum2 += 0.00032 * math.sin(_radians(omega2 - omega3))
    sum2 += 0.00032 * math.sin(_radians(2 * (l3 - G - PI)))

    # Satellite 3 (Ganymede)
    sum3 = 0.0
    sum3 += 0.16490 * math.sin(_radians(l3 - pi3))
    sum3 += 0.09081 * math.sin(_radians(l3 - pi4))
    sum3 -= 0.06907 * math.sin(_radians(l2 - l3))
    sum3 += 0.03784 * math.sin(_radians(pi3 - pi4))
    sum3 += 0.01846 * math.sin(_radians(2 * (l3 - l4)))
    sum3 -= 0.01340 * math.sin(_radians(G))
    sum3 -= 0.01014 * math.sin(_radians(2 * (psi - PI)))
    sum3 += 0.00704 * math.sin(_radians(l2 - 2 * l3 + pi3))
    sum3 -= 0.00620 * math.sin(_radians(l2 - 2 * l3 + pi2))
    sum3 -= 0.00541 * math.sin(_radians(l3 - l4))
    sum3 += 0.00381 * math.sin(_radians(l2 - 2 * l3 + pi4))
    sum3 += 0.00235 * math.sin(_radians(psi - omega3))
    sum3 += 0.00198 * math.sin(_radians(psi - omega4))
    sum3 += 0.00176 * math.sin(_radians(PHI_lambda))
    sum3 += 0.00130 * math.sin(_radians(3 * (l3 - l4)))
    sum3 += 0.00125 * math.sin(_radians(l1 - l3))
    sum3 -= 0.00119 * math.sin(_radians(5 * G_ - 2 * G + 52.225))
    sum3 += 0.00109 * math.sin(_radians(l1 - l2))
    sum3 -= 0.00100 * math.sin(_radians(3 * l3 - 7 * l4 + 4 * pi4))
    sum3 += 0.00091 * math.sin(_radians(omega3 - omega4))
    sum3 += 0.00080 * math.sin(_radians(3 * l3 - 7 * l4 + pi3 + 3 * pi4))
    sum3 -= 0.00075 * math.sin(_radians(2 * l2 - 3 * l3 + pi3))
    sum3 += 0.00072 * math.sin(_radians(pi1 + pi3 - 2 * PI - 2 * G))
    sum3 += 0.00069 * math.sin(_radians(pi4 - PI))
    sum3 -= 0.00058 * math.sin(_radians(2 * l3 - 3 * l4 + pi4))
    sum3 -= 0.00057 * math.sin(_radians(l3 - 2 * l4 + pi4))
    sum3 += 0.00056 * math.sin(_radians(l3 + pi3 - 2 * PI - 2 * G))
    sum3 -= 0.00052 * math.sin(_radians(l2 - 2 * l3 + pi1))
    sum3 -= 0.00050 * math.sin(_radians(pi2 - pi3))
    sum3 += 0.00048 * math.sin(_radians(l3 - 2 * l4 + pi3))
    sum3 -= 0.00045 * math.sin(_radians(2 * l2 - 3 * l3 + pi4))
    sum3 -= 0.00041 * math.sin(_radians(pi2 - pi4))
    sum3 -= 0.00038 * math.sin(_radians(2 * G))
    sum3 -= 0.00037 * math.sin(_radians(pi3 - pi4 + omega3 - omega4))
    sum3 -= 0.00032 * math.sin(_radians(3 * l3 - 7 * l4 + 2 * pi3 + 2 * pi4))
    sum3 += 0.00030 * math.sin(_radians(4 * (l3 - l4)))
    sum3 += 0.00029 * math.sin(_radians(l3 + pi4 - 2 * PI - 2 * G))
    sum3 -= 0.00028 * math.sin(_radians(omega3 + psi - 2 * PI - 2 * G))
    sum3 += 0.00026 * math.sin(_radians(l3 - PI - G))
    sum3 += 0.00024 * math.sin(_radians(l2 - 3 * l3 + 2 * l4))
    sum3 += 0.00021 * math.sin(_radians(2 * (l3 - PI - G)))
    sum3 -= 0.00021 * math.sin(_radians(l3 - pi2))
    sum3 += 0.00017 * math.sin(_radians(2 * (l3 - pi3)))

    # Satellite 4 (Callisto)
    sum4 = 0.0
    sum4 += 0.84287 * math.sin(_radians(l4 - pi4))
    sum4 += 0.03431 * math.sin(_radians(pi4 - pi3))
    sum4 -= 0.03305 * math.sin(_radians(2 * (psi - PI)))
    sum4 -= 0.03211 * math.sin(_radians(G))
    sum4 -= 0.01862 * math.sin(_radians(l4 - pi3))
    sum4 += 0.01186 * math.sin(_radians(psi - omega4))
    sum4 += 0.00623 * math.sin(_radians(l4 + pi4 - 2 * G - 2 * PI))
    sum4 += 0.00387 * math.sin(_radians(2 * (l4 - pi4)))
    sum4 -= 0.00284 * math.sin(_radians(5 * G_ - 2 * G + 52.225))
    sum4 -= 0.00234 * math.sin(_radians(2 * (psi - pi4)))
    sum4 -= 0.00223 * math.sin(_radians(l3 - l4))
    sum4 -= 0.00208 * math.sin(_radians(l4 - PI))
    sum4 += 0.00178 * math.sin(_radians(psi + omega4 - 2 * pi4))
    sum4 += 0.00134 * math.sin(_radians(pi4 - PI))
    sum4 += 0.00125 * math.sin(_radians(2 * (l4 - G - PI)))
    sum4 -= 0.00117 * math.sin(_radians(2 * G))
    sum4 -= 0.00112 * math.sin(_radians(2 * (l3 - l4)))
    sum4 += 0.00107 * math.sin(_radians(3 * l3 - 7 * l4 + 4 * pi4))
    sum4 += 0.00102 * math.sin(_radians(l4 - G - PI))
    sum4 += 0.00096 * math.sin(_radians(2 * l4 - psi - omega4))
    sum4 += 0.00087 * math.sin(_radians(2 * (psi - omega4)))
    sum4 -= 0.00085 * math.sin(_radians(3 * l3 - 7 * l4 + pi3 + 3 * pi4))
    sum4 += 0.00085 * math.sin(_radians(l3 - 2 * l4 + pi4))
    sum4 -= 0.00081 * math.sin(_radians(2 * (l4 - psi)))
    sum4 += 0.00071 * math.sin(_radians(l4 + pi4 - 2 * PI - 3 * G))
    sum4 += 0.00061 * math.sin(_radians(l1 - l4))
    sum4 -= 0.00056 * math.sin(_radians(psi - omega3))
    sum4 -= 0.00054 * math.sin(_radians(l3 - 2 * l4 + pi3))
    sum4 += 0.00051 * math.sin(_radians(l2 - l4))
    sum4 += 0.00042 * math.sin(_radians(2 * (psi - G - PI)))
    sum4 += 0.00039 * math.sin(_radians(2 * (pi4 - omega4)))
    sum4 += 0.00036 * math.sin(_radians(psi + PI - pi4 - omega4))
    sum4 += 0.00035 * math.sin(_radians(2 * G_ - G + 188.37))
    sum4 -= 0.00035 * math.sin(_radians(l4 - pi4 + 2 * PI - 2 * psi))
    sum4 -= 0.00032 * math.sin(_radians(l4 + pi4 - 2 * PI - G))
    sum4 += 0.00030 * math.sin(_radians(2 * G_ - 2 * G + 149.15))
    sum4 += 0.00029 * math.sin(_radians(3 * l3 - 7 * l4 + 2 * pi3 + 2 * pi4))
    sum4 += 0.00028 * math.sin(_radians(l4 - pi4 + 2 * psi - 2 * PI))
    sum4 -= 0.00028 * math.sin(_radians(2 * (l4 - omega4)))
    sum4 -= 0.00027 * math.sin(_radians(pi3 - pi4 + omega3 - omega4))
    sum4 -= 0.00026 * math.sin(_radians(5 * G_ - 3 * G + 188.37))
    sum4 += 0.00025 * math.sin(_radians(omega4 - omega3))
    sum4 -= 0.00025 * math.sin(_radians(l2 - 3 * l3 + 2 * l4))
    sum4 -= 0.00023 * math.sin(_radians(3 * (l3 - l4)))
    sum4 += 0.00021 * math.sin(_radians(2 * l4 - 2 * PI - 3 * G))
    sum4 -= 0.00021 * math.sin(_radians(2 * l3 - 3 * l4 + pi4))
    sum4 += 0.00019 * math.sin(_radians(l4 - pi4 - G))
    sum4 -= 0.00019 * math.sin(_radians(2 * l4 - pi3 - pi4))
    sum4 -= 0.00018 * math.sin(_radians(l4 - pi4 + G))
    sum4 -= 0.00016 * math.sin(_radians(l4 + pi3 - 2 * PI - 2 * G))

    # True longitudes of the satellites (degrees)
    L1 = l1 + sum1
    L2 = l2 + sum2
    L3 = l3 + sum3
    L4 = l4 + sum4

    # ==========================================================================
    # PERIODIC TERMS IN LATITUDE
    # ==========================================================================

    tan_B1 = (
        +0.0006393 * math.sin(_radians(L1 - omega1))
        + 0.0001825 * math.sin(_radians(L1 - omega2))
        + 0.0000329 * math.sin(_radians(L1 - omega3))
        - 0.0000311 * math.sin(_radians(L1 - psi))
        + 0.0000093 * math.sin(_radians(L1 - omega4))
        + 0.0000075 * math.sin(_radians(3 * L1 - 4 * l2 - 1.9927 * sum1 + omega2))
        + 0.0000046 * math.sin(_radians(L1 + psi - 2 * PI - 2 * G))
    )

    tan_B2 = (
        +0.0081004 * math.sin(_radians(L2 - omega2))
        + 0.0004512 * math.sin(_radians(L2 - omega3))
        - 0.0003284 * math.sin(_radians(L2 - psi))
        + 0.0001160 * math.sin(_radians(L2 - omega4))
        + 0.0000272 * math.sin(_radians(l1 - 2 * l3 + 1.0146 * sum2 + omega2))
        - 0.0000144 * math.sin(_radians(L2 - omega1))
        + 0.0000143 * math.sin(_radians(L2 + psi - 2 * PI - 2 * G))
        + 0.0000035 * math.sin(_radians(L2 - psi + G))
        - 0.0000028 * math.sin(_radians(l1 - 2 * l3 + 1.0146 * sum2 + omega3))
    )

    tan_B3 = (
        +0.0032402 * math.sin(_radians(L3 - omega3))
        - 0.0016911 * math.sin(_radians(L3 - psi))
        + 0.0006847 * math.sin(_radians(L3 - omega4))
        - 0.0002797 * math.sin(_radians(L3 - omega2))
        + 0.0000321 * math.sin(_radians(L3 + psi - 2 * PI - 2 * G))
        + 0.0000051 * math.sin(_radians(L3 - psi + G))
        - 0.0000045 * math.sin(_radians(L3 - psi - G))
        - 0.0000045 * math.sin(_radians(L3 + psi - 2 * PI))
        + 0.0000037 * math.sin(_radians(L3 + psi - 2 * PI - 3 * G))
        + 0.0000030 * math.sin(_radians(2 * l2 - 3 * L3 + 4.03 * sum3 + omega2))
        - 0.0000021 * math.sin(_radians(2 * l2 - 3 * L3 + 4.03 * sum3 + omega3))
    )

    tan_B4 = (
        -0.0076579 * math.sin(_radians(L4 - psi))
        + 0.0044134 * math.sin(_radians(L4 - omega4))
        - 0.0005112 * math.sin(_radians(L4 - omega3))
        + 0.0000773 * math.sin(_radians(L4 + psi - 2 * PI - 2 * G))
        + 0.0000104 * math.sin(_radians(L4 - psi + G))
        - 0.0000102 * math.sin(_radians(L4 - psi - G))
        + 0.0000088 * math.sin(_radians(L4 + psi - 2 * PI - 3 * G))
        - 0.0000038 * math.sin(_radians(L4 + psi - 2 * PI - G))
    )

    B1 = math.atan(tan_B1)
    B2 = math.atan(tan_B2)
    B3 = math.atan(tan_B3)
    B4 = math.atan(tan_B4)

    # ==========================================================================
    # PERIODIC TERMS IN RADIUS
    # ==========================================================================

    sum_r1 = (
        -0.0041339 * math.cos(_radians(2 * (l1 - l2)))
        - 0.0000387 * math.cos(_radians(l1 - pi3))
        - 0.0000214 * math.cos(_radians(l1 - pi4))
        + 0.0000170 * math.cos(_radians(l1 - l2))
        - 0.0000131 * math.cos(_radians(4 * (l1 - l2)))
        + 0.0000106 * math.cos(_radians(l1 - l3))
        - 0.0000066 * math.cos(_radians(l1 + pi3 - 2 * PI - 2 * G))
    )

    sum_r2 = (
        +0.0093848 * math.cos(_radians(l1 - l2))
        - 0.0003116 * math.cos(_radians(l2 - pi3))
        - 0.0001744 * math.cos(_radians(l2 - pi4))
        - 0.0001442 * math.cos(_radians(l2 - pi2))
        + 0.0000553 * math.cos(_radians(l2 - l3))
        + 0.0000523 * math.cos(_radians(l1 - l3))
        - 0.0000290 * math.cos(_radians(2 * (l1 - l2)))
        + 0.0000164 * math.cos(_radians(2 * (l2 - omega2)))
        + 0.0000107 * math.cos(_radians(l1 - 2 * l3 + pi3))
        - 0.0000102 * math.cos(_radians(l2 - pi1))
        - 0.0000091 * math.cos(_radians(2 * (l1 - l3)))
    )

    sum_r3 = (
        -0.0014388 * math.cos(_radians(l3 - pi3))
        - 0.0007919 * math.cos(_radians(l3 - pi4))
        + 0.0006342 * math.cos(_radians(l2 - l3))
        - 0.0001761 * math.cos(_radians(2 * (l3 - l4)))
        + 0.0000294 * math.cos(_radians(l3 - l4))
        - 0.0000156 * math.cos(_radians(3 * (l3 - l4)))
        + 0.0000156 * math.cos(_radians(l1 - l3))
        - 0.0000153 * math.cos(_radians(l1 - l2))
        + 0.0000070 * math.cos(_radians(2 * l2 - 3 * l3 + pi3))
        - 0.0000051 * math.cos(_radians(l3 + pi3 - 2 * PI - 2 * G))
    )

    sum_r4 = (
        -0.0073546 * math.cos(_radians(l4 - pi4))
        + 0.0001621 * math.cos(_radians(l4 - pi3))
        + 0.0000974 * math.cos(_radians(l3 - l4))
        - 0.0000543 * math.cos(_radians(l4 + pi4 - 2 * PI - 2 * G))
        - 0.0000271 * math.cos(_radians(2 * (l4 - pi4)))
        + 0.0000182 * math.cos(_radians(l4 - PI))
        + 0.0000177 * math.cos(_radians(2 * (l3 - l4)))
        - 0.0000167 * math.cos(_radians(2 * l4 - psi - omega4))
        + 0.0000167 * math.cos(_radians(psi - omega4))
        - 0.0000155 * math.cos(_radians(2 * (l4 - PI - G)))
        + 0.0000142 * math.cos(_radians(2 * (l4 - psi)))
        + 0.0000105 * math.cos(_radians(l1 - l4))
        + 0.0000092 * math.cos(_radians(l2 - l4))
        - 0.0000089 * math.cos(_radians(l4 - PI - G))
        - 0.0000062 * math.cos(_radians(l4 + pi4 - 2 * PI - 3 * G))
        + 0.0000048 * math.cos(_radians(2 * (l4 - omega4)))
    )

    # Radius vectors in Jupiter radii
    R1 = A1 * (1.0 + sum_r1)
    R2 = A2 * (1.0 + sum_r2)
    R3 = A3 * (1.0 + sum_r3)
    R4 = A4 * (1.0 + sum_r4)

    # ==========================================================================
    # PRECESSION CORRECTION (B1950 -> J2000, fixed)
    # ==========================================================================

    # The E5 longitudes are referenced to the mean equinox of B1950.
    # Precess by the constant arc B1950 -> J2000 so the output frame is
    # the J2000 ecliptic (the COB consumer then applies one fixed
    # rotation to equatorial ICRF).  Precessing to the equinox of date
    # here would silently mix frames at the consumer.
    T0 = (2451545.0 - 2433282.423) / 36525.0

    # Precession in longitude from epoch B1950.0 to J2000.0 (degrees)
    P = 1.3966626 * T0 + 0.0003088 * T0 * T0

    # Corrected longitudes
    L1_corr = L1 + P
    L2_corr = L2 + P
    L3_corr = L3 + P
    L4_corr = L4 + P
    psi_corr = psi + P

    # ==========================================================================
    # CARTESIAN COORDINATES IN JUPITER'S EQUATORIAL FRAME
    # ==========================================================================

    # Position in Jupiter's equatorial frame (Jupiter radii)
    # X towards psi, Y perpendicular, Z towards pole

    def to_cartesian(R, L, B, psi_c):
        """Convert orbital elements to Cartesian in Jupiter's frame."""
        X = R * math.cos(_radians(L - psi_c)) * math.cos(B)
        Y = R * math.sin(_radians(L - psi_c)) * math.cos(B)
        Z = R * math.sin(B)
        return X, Y, Z

    X1, Y1, Z1 = to_cartesian(R1, L1_corr, B1, psi_corr)
    X2, Y2, Z2 = to_cartesian(R2, L2_corr, B2, psi_corr)
    X3, Y3, Z3 = to_cartesian(R3, L3_corr, B3, psi_corr)
    X4, Y4, Z4 = to_cartesian(R4, L4_corr, B4, psi_corr)

    # ==========================================================================
    # TRANSFORM TO J2000 ECLIPTIC FRAME
    # ==========================================================================

    # Jupiter's orbit orientation, held at J2000 to match the J2000
    # equinox the longitudes were precessed to above.
    OMEGA_J = _radians(100.464407)  # Ascending node of Jupiter's orbit
    i_J = _radians(1.303267)  # Inclination of Jupiter's orbit

    # Inclination of Jupiter's equator to orbit (at J2000)
    Inc = _radians(3.120262 + 0.0006 * (2451545.0 - 2415020.5) / 36525.0)

    psi_rad = _radians(psi_corr)

    def transform_to_ecliptic(X, Y, Z):
        """Transform from Jupiter equatorial to J2000 ecliptic."""
        # Rotation sequence to ecliptic frame

        # 1. Rotate by inclination of Jupiter's equator
        A1 = X
        B1 = Y * math.cos(Inc) - Z * math.sin(Inc)
        C1 = Y * math.sin(Inc) + Z * math.cos(Inc)

        # 2. Rotate by angle psi - Omega (position angle of node)
        PHI = psi_rad - OMEGA_J
        A2 = A1 * math.cos(PHI) - B1 * math.sin(PHI)
        B2 = A1 * math.sin(PHI) + B1 * math.cos(PHI)
        C2 = C1

        # 3. Rotate by inclination of Jupiter's orbit
        A3 = A2
        B3 = B2 * math.cos(i_J) - C2 * math.sin(i_J)
        C3 = B2 * math.sin(i_J) + C2 * math.cos(i_J)

        # 4. Rotate to vernal equinox
        A4 = A3 * math.cos(OMEGA_J) - B3 * math.sin(OMEGA_J)
        B4 = A3 * math.sin(OMEGA_J) + B3 * math.cos(OMEGA_J)
        C4 = C3

        return A4, B4, C4

    # Transform each satellite
    x1, y1, z1 = transform_to_ecliptic(X1, Y1, Z1)
    x2, y2, z2 = transform_to_ecliptic(X2, Y2, Z2)
    x3, y3, z3 = transform_to_ecliptic(X3, Y3, Z3)
    x4, y4, z4 = transform_to_ecliptic(X4, Y4, Z4)

    # Convert from Jupiter radii to km
    io_km = (x1 * JUPITER_RADIUS_KM, y1 * JUPITER_RADIUS_KM, z1 * JUPITER_RADIUS_KM)
    europa_km = (x2 * JUPITER_RADIUS_KM, y2 * JUPITER_RADIUS_KM, z2 * JUPITER_RADIUS_KM)
    ganymede_km = (
        x3 * JUPITER_RADIUS_KM,
        y3 * JUPITER_RADIUS_KM,
        z3 * JUPITER_RADIUS_KM,
    )
    callisto_km = (
        x4 * JUPITER_RADIUS_KM,
        y4 * JUPITER_RADIUS_KM,
        z4 * JUPITER_RADIUS_KM,
    )

    return io_km, europa_km, ganymede_km, callisto_km
