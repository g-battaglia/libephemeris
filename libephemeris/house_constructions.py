# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Vector constructions for horizon-anchored house systems.

This module computes the intermediate cusps of four house systems whose
defining geometry anchors great circles on the horizon frame:

- Savard-A (``J``): position circles through the north and south points of
  the horizon, anchored where the prime vertical meets the declination
  parallels at one third and two thirds of the geographic latitude.
  Author-published definition: John Savard, "Astrological House Systems",
  http://www.quadibloc.com/other/as01.htm ("Albategnus" entry: the prime
  vertical is divided by parallels of declination at thirds of the
  east-point-to-zenith declination difference, then great circles from the
  north point to the south point divide the ecliptic).
- Krusinski-Pisa (``U``): the great circle through the Ascendant and the
  zenith is divided into twelve 30-degree arcs starting at the Ascendant,
  and each division point is carried to the ecliptic along its hour circle
  (the great circle through the celestial poles). Public definition:
  Bogdan Krusinski's own 1995 description,
  https://www.astrologia.pl/house-system.html. The author explicitly states
  both the 30-degree division and projection along meridian circles.
- Sunshine/Makransky (``i``): the Sun's diurnal arc is divided into six
  equal parts, and so is its nocturnal arc; the division points are carried
  to the ecliptic along house circles through the north and south points of
  the horizon. Public definition: B. Makransky, "The Sunshine House System"
  (author's article, reproduced at
  https://www.astro.nu/2017/03/03/astrology-house-systems/).
- APC (``Y``): the same construction as the Sunshine system applied to the
  parallel of declination of the Ascendant instead of the Sun's parallel
  (Dutch School of Ram; L. Knegt). The first cusp is the Ascendant itself
  and the tenth house circle is the meridian. Independent construction
  description: Ingmar de Boer, "APC Houses", paragraphs 1--3,
  https://www.ingmardeboer.nl/index.php?title=APC_Houses (the page describes
  the Ascendant Parallel Circle, its separate sixfold divisions above and
  below the horizon, and projection through oblique-ascension planes).

Savard-A, Krusinski-Pisa and APC use project-authored Cartesian vector
constructions: build
the relevant anchor in the equatorial frame of date, form the required
great-circle plane, and intersect it with the ecliptic. Sunshine-Makransky
instead follows Makransky's published spherical-trigonometry construction,
including its explicit upper/lower-meridian convention. Both formulations
select branches from geometry rather than from copied numerical tables.

Frames and conventions:
    - Equatorial right-handed frame of date: x toward the equinox, z toward
      the north celestial pole. A direction with right ascension ``a`` and
      declination ``d`` is ``(cos d cos a, cos d sin a, sin d)``.
    - The zenith has right ascension ``ARMC`` and declination ``lat``.
    - Hour angles are measured westward: a point of a declination parallel
      at hour angle ``H`` has right ascension ``ARMC - H``.
    - All public functions take degrees and return longitudes in degrees in
      ``[0, 360)``, except :func:`apc_cusp` whose scalar inputs are radians
      (its caller works in radians).

Provenance:
    The four geometric definitions and their public locators are stated above;
    ``docs/reference/house-systems.md`` preserves the complete source record.
    Cartesian cross/dot products, oriented-arc branch selection, and degeneracy
    handling are project-authored derivations from those definitions. Branches
    are selected by geometric orientation and anchor invariants, not by stored
    outputs. Epsilon guards are numerical choices documented beside their use.
"""

from __future__ import annotations

import math
from typing import List, Sequence, Tuple

from .exceptions import CalculationError, PolarCircleError

Vec3 = Tuple[float, float, float]

# Guard for degenerate cross products (anchor parallel to the circle axis).
_TINY = 1e-15


def _unit_radec(ra_rad: float, dec_rad: float) -> Vec3:
    """Unit vector of a direction given by right ascension/declination."""
    cd = math.cos(dec_rad)
    return (cd * math.cos(ra_rad), cd * math.sin(ra_rad), math.sin(dec_rad))


def _cross(u: Sequence[float], v: Sequence[float]) -> Vec3:
    """Return the right-handed three-dimensional cross product ``u x v``."""
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def _dot(u: Sequence[float], v: Sequence[float]) -> float:
    """Return the Euclidean dot product of two three-component vectors."""
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def _normalized(v: Sequence[float]) -> Vec3:
    """Return ``v`` scaled to unit length; callers exclude degenerate input."""
    n = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    return (v[0] / n, v[1] / n, v[2] / n)


def _ecliptic_lon_deg(v: Sequence[float], eps_rad: float) -> float:
    """Ecliptic longitude of an equatorial direction, degrees in [0, 360)."""
    lon = math.atan2(v[1] * math.cos(eps_rad) + v[2] * math.sin(eps_rad), v[0])
    return math.degrees(lon) % 360.0


def _hour_circle_lon_deg(v: Sequence[float], eps_rad: float) -> float:
    """Longitude where the hour circle through ``v`` meets the ecliptic.

    The hour circle is the great circle through both celestial poles and
    ``v``; it crosses the ecliptic where the ecliptic has the same right
    ascension as ``v``.
    """
    ra = math.atan2(v[1], v[0])
    lon = math.atan2(math.sin(ra), math.cos(ra) * math.cos(eps_rad))
    return math.degrees(lon) % 360.0


def _horizon_frame(armc_rad: float, lat_rad: float) -> Tuple[Vec3, Vec3, Vec3]:
    """Zenith, east point, and north point of the horizon as unit vectors."""
    zenith = _unit_radec(armc_rad, lat_rad)
    east = _unit_radec(armc_rad + math.pi / 2.0, 0.0)
    north = _cross(zenith, east)
    return zenith, east, north


def _house_circle_lon_deg(
    anchor: Sequence[float],
    axis: Sequence[float],
    eps_rad: float,
    *,
    same_side: bool = True,
) -> float:
    """Ecliptic longitude cut by the great circle through ``axis`` and ``anchor``.

    ``axis`` is the unit vector of the horizon's north point; the circle also
    contains the south point ``-axis``. Of the two ecliptic intersections the
    ``same_side`` selects the intersection on the same half of the circle as
    ``anchor``; when false it selects the antipodal half.  The two halves are
    bounded by the axis endpoints.

    Args:
        anchor: Unit vector the circle must contain.
        axis: Unit vector of the circle axis endpoint (north horizon point).
        eps_rad: True obliquity in radians.
        same_side: Whether to return the ecliptic crossing on the anchor's
            half of the circle rather than the antipodal half.

    Returns:
        Ecliptic longitude in degrees in [0, 360).
    """
    normal = _cross(axis, anchor)
    n2 = _dot(normal, normal)
    if n2 < _TINY:
        # Anchor coincides with the axis: every circle through the axis
        # contains it. Fall back to the hour-circle carry of the anchor.
        longitude = _hour_circle_lon_deg(anchor, eps_rad)
        return longitude if same_side else (longitude + 180.0) % 360.0
    ecl_pole = (0.0, -math.sin(eps_rad), math.cos(eps_rad))
    meet = _cross(normal, ecl_pole)
    m2 = _dot(meet, meet)
    if m2 < _TINY:
        # The house circle lies in the ecliptic plane; the anchor itself is
        # on the ecliptic.
        longitude = _ecliptic_lon_deg(anchor, eps_rad)
        return longitude if same_side else (longitude + 180.0) % 360.0
    # ``half`` lies in the house-circle plane and is perpendicular to the
    # north/south axis.  Consequently dot(point, half) is positive on one
    # open semicircle, negative on the other, and zero only at the two axis
    # endpoints.  Comparing the signs for ``meet`` and ``anchor`` therefore
    # chooses the branch without consulting longitude quadrants: their dot
    # product is positive exactly when both points occupy the same half of
    # the oriented circle.  This invariant is stable across 0/360 longitude,
    # both geographic hemispheres, and polar representations where a test
    # based on east/west longitude would select the antipode.
    half = _cross(normal, axis)
    side_product = _dot(meet, half) * _dot(anchor, half)
    if (side_product < 0.0) == same_side:
        meet = (-meet[0], -meet[1], -meet[2])
    return _ecliptic_lon_deg(meet, eps_rad)


def _semiarc_cusp_hour_angles(day_semiarc_rad: float) -> List[Tuple[int, float]]:
    """Hour angles of the eight intermediate division points of a parallel.

    ``day_semiarc_rad`` is half the diurnal arc of the parallel. Cusps 12/11
    and 9/8 trisect the diurnal semi-arcs east and west of the meridian;
    cusps 2/3 and 5/6 trisect the nocturnal semi-arcs on either side of the
    lower meridian.

    Returns:
        (cusp index, hour angle in radians) pairs, hour angle west-positive.
    """
    day = day_semiarc_rad
    night = math.pi - day
    return [
        (11, -day / 3.0),
        (12, -2.0 * day / 3.0),
        (2, -day - night / 3.0),
        (3, -day - 2.0 * night / 3.0),
        (9, day / 3.0),
        (8, 2.0 * day / 3.0),
        (6, day + night / 3.0),
        (5, day + 2.0 * night / 3.0),
    ]


def _prime_vertical_arc_to_declination(
    lat: float, sin_lat: float, thirds: int
) -> float:
    """Arc along the prime vertical, from the east point, to a parallel.

    On the prime vertical, parametrized by the arc ``t`` from the east point
    toward the zenith, a point has ``sin(dec) = sin(t) sin(lat)``. The
    parallel of declination at ``thirds/3`` of the latitude is therefore
    reached at ``sin(t) = sin(thirds * lat / 3) / sin(lat)``. On the equator
    both sines vanish together and the ratio tends to ``thirds / 3``.

    Args:
        lat: Geographic latitude, degrees.
        sin_lat: ``sin(lat)``, supplied by the caller.
        thirds: 1 or 2, the fraction of the latitude span to reach.

    Returns:
        The arc ``t`` in radians, in ``[-pi/2, pi/2]``.
    """
    if abs(sin_lat) < 1e-10:
        sin_arc = thirds / 3.0
    else:
        sin_arc = math.sin(math.radians(thirds * lat / 3.0)) / sin_lat
    return math.asin(min(1.0, max(-1.0, sin_arc)))


def _prime_vertical_point(
    east: Vec3, zenith: Vec3, arc_rad: float, toward_east: bool
) -> Vec3:
    """Unit vector on the prime vertical at ``arc_rad`` above the horizon.

    The arc is measured from the east point when ``toward_east`` is true and
    from the west point otherwise, in both cases toward the zenith.
    """
    along = math.cos(arc_rad)
    if not toward_east:
        along = -along
    up = math.sin(arc_rad)
    return (
        along * east[0] + up * zenith[0],
        along * east[1] + up * zenith[1],
        along * east[2] + up * zenith[2],
    )


#: Savard-A anchors: ``(cusp, thirds of the latitude span, on the east side)``.
#: The parallel at two thirds anchors cusps 11 (east) and 3 (west), the one at
#: one third anchors cusps 12 (east) and 2 (west).
_SAVARD_ANCHORS: Tuple[Tuple[int, int, bool], ...] = (
    (11, 2, True),
    (3, 2, False),
    (12, 1, True),
    (2, 1, False),
)

#: Cusps Savard-A constructs directly and the opposite house each one is
#: mirrored into.
_OPPOSITE_HOUSES: Tuple[Tuple[int, int], ...] = (
    (1, 7),
    (10, 4),
    (11, 5),
    (12, 6),
    (2, 8),
    (3, 9),
)


def houses_savard_a(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """Savard-A house cusps (system letter ``J``).

    The prime vertical rises from the east point (declination zero) to the
    zenith (declination equal to the latitude). The parallels of declination
    at one third and two thirds of the latitude divide that declination span
    into thirds; where they cross the prime vertical they anchor position
    circles through the north and south points of the horizon. The circle
    anchored at two thirds gives cusps 11/5, the one at one third gives
    cusps 12/6; their west-side mirror images give cusps 3/9 and 2/8.
    :func:`_prime_vertical_arc_to_declination` locates each anchor along the
    prime vertical and :func:`_prime_vertical_point` turns it into a vector.

    Args:
        armc: Right ascension of the midheaven, degrees.
        lat: Geographic latitude, degrees.
        eps: True obliquity, degrees.
        asc: Ascendant longitude, degrees (cusp 1).
        mc: Midheaven longitude, degrees (cusp 10).

    Returns:
        Thirteen floats; index 0 unused, 1-12 are the cusp longitudes.
    """
    eps_rad = math.radians(eps)
    lat_rad = math.radians(lat)
    zenith, east, north = _horizon_frame(math.radians(armc), lat_rad)
    asc_rad = math.radians(asc)
    asc_vec = (
        math.cos(asc_rad),
        math.sin(asc_rad) * math.cos(eps_rad),
        math.sin(asc_rad) * math.sin(eps_rad),
    )
    # At polar latitudes the dispatcher can represent the upper meridian with
    # the antipodal horizon frame.  Tie the east/west cusp numbering to the
    # actual rising point, not to that representational choice: an anchor is
    # on the rising side when its side of the frame is the Ascendant's.
    frame_faces_east = _dot(east, asc_vec) >= 0.0

    cusps = [0.0] * 13
    cusps[1] = asc
    cusps[10] = mc
    sin_lat = math.sin(lat_rad)
    for cusp, thirds, toward_east in _SAVARD_ANCHORS:
        arc = _prime_vertical_arc_to_declination(lat, sin_lat, thirds)
        anchor = _prime_vertical_point(east, zenith, arc, toward_east)
        cusps[cusp] = _house_circle_lon_deg(
            anchor,
            north,
            eps_rad,
            same_side=toward_east == frame_faces_east,
        )
    for built, opposite in _OPPOSITE_HOUSES:
        cusps[opposite] = (cusps[built] + 180.0) % 360.0
    return cusps


def houses_krusinski(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """Krusinski-Pisa house cusps (system letter ``U``).

    The great circle through the Ascendant and the zenith (a vertical
    circle, since the Ascendant lies on the horizon) is divided into twelve
    30-degree arcs starting at the Ascendant and counted toward the zenith:
    the zenith arc is 90 degrees, so cusp 10 is the zenith's division point.
    Each division point is carried to the ecliptic along its hour circle;
    the zenith thereby maps to the midheaven.

    Args:
        armc: Right ascension of the midheaven, degrees.
        lat: Geographic latitude, degrees.
        eps: True obliquity, degrees.
        asc: Ascendant longitude, degrees.
        mc: Midheaven longitude, degrees (unused; cusp 10 is constructed).

    Returns:
        Thirteen floats; index 0 unused, 1-12 are the cusp longitudes.
    """
    eps_rad = math.radians(eps)
    asc_rad = math.radians(asc)
    # Ascendant direction: ecliptic latitude zero at longitude ``asc``.
    a_vec = (
        math.cos(asc_rad),
        math.sin(asc_rad) * math.cos(eps_rad),
        math.sin(asc_rad) * math.sin(eps_rad),
    )
    zenith = _unit_radec(math.radians(armc), math.radians(lat))
    # Component of the zenith orthogonal to the Ascendant: second basis
    # vector of the Ascendant-zenith circle. The Ascendant lies on the
    # horizon, so the two are already orthogonal up to rounding.
    z_dot_a = _dot(zenith, a_vec)
    up = _normalized(
        (
            zenith[0] - z_dot_a * a_vec[0],
            zenith[1] - z_dot_a * a_vec[1],
            zenith[2] - z_dot_a * a_vec[2],
        )
    )
    cusps = [0.0] * 13
    for cusp, arc_deg in (
        (1, 0.0),
        (12, 30.0),
        (11, 60.0),
        (10, 90.0),
        (2, -30.0),
        (3, -60.0),
    ):
        arc = math.radians(arc_deg)
        point = (
            math.cos(arc) * a_vec[0] + math.sin(arc) * up[0],
            math.cos(arc) * a_vec[1] + math.sin(arc) * up[1],
            math.cos(arc) * a_vec[2] + math.sin(arc) * up[2],
        )
        cusps[cusp] = _hour_circle_lon_deg(point, eps_rad)
    for src, opp in ((1, 7), (2, 8), (3, 9), (10, 4), (11, 5), (12, 6)):
        cusps[opp] = (cusps[src] + 180.0) % 360.0
    return cusps


def _makransky_ecliptic_cut(
    w_deg: float,
    pole_deg: float,
    eps_rad: float,
    *,
    eastern_family: bool,
) -> float:
    """Project a Makransky house circle onto the ecliptic.

    This is a compact Cartesian/algebraic replacement for steps 6--9 of the
    published *Solar House Cusp Algorithm* (Makransky, 1990, pp. 147--148).
    Makransky gives a quadrant table involving the auxiliary angles ``M``,
    ``Z`` and ``R``.  Eliminating those auxiliaries gives the standard great-
    circle/ecliptic intersection

    ``tan(lambda) = sin(W) / (cos(eps) cos(W) - sin(eps) tan(p))``.

    ``atan2`` retains the quadrant directly, so no transcription of the
    published case table is needed.  The sign of ``p`` selects which of the
    two orientations of the same house circle belongs to the cusp-number
    family.  Makransky's cusps 11, 12, 2 and 3 use the positive orientation;
    cusps 5, 6, 8 and 9 use the negative orientation.  The magnitude, rather
    than the algebraic sign returned by the preceding southern-hemisphere
    transformation, is therefore intentional.

    Args:
        w_deg: Oblique ascension/descension ``W`` under the circle's own pole.
        pole_deg: Algebraic pole height computed from zenith distance.
        eps_rad: Obliquity of the ecliptic, in radians.
        eastern_family: True for cusps 11, 12, 2 and 3; false for the
            complementary cusp family 5, 6, 8 and 9.

    Returns:
        Ecliptic longitude in degrees in ``[0, 360)``.
    """
    effective_pole = abs(pole_deg) if eastern_family else -abs(pole_deg)
    w_rad = math.radians(w_deg)
    pole_rad = math.radians(effective_pole)
    longitude = math.atan2(
        math.sin(w_rad),
        math.cos(eps_rad) * math.cos(w_rad) - math.sin(eps_rad) * math.tan(pole_rad),
    )
    return math.degrees(longitude) % 360.0


def _makransky_intermediate_cusp(
    division_ra_deg: float,
    working_armc_deg: float,
    lat_deg: float,
    sun_dec_deg: float,
    eps_rad: float,
) -> float:
    """Project one solar-parallel division through the horizon axis.

    Args:
        division_ra_deg: Right ascension of the selected division point.
        working_armc_deg: Upper-meridian right ascension of the working frame.
        lat_deg: Signed geographic latitude.
        sun_dec_deg: Signed solar declination.
        eps_rad: Obliquity in radians.

    Returns:
        The oriented ecliptic intersection in degrees, normalized to
        ``[0, 360)``.

    Raises:
        CalculationError: If the house and ecliptic planes coincide without a
            unique oriented projective limit.
    """
    armc_rad = math.radians(working_armc_deg)
    latitude_rad = math.radians(lat_deg)
    division_ra_rad = math.radians(division_ra_deg)
    declination_rad = math.radians(sun_dec_deg)
    _zenith, _east, north = _horizon_frame(armc_rad, latitude_rad)
    ecliptic_normal = (0.0, -math.sin(eps_rad), math.cos(eps_rad))

    cos_declination = math.cos(declination_rad)
    division = _unit_radec(division_ra_rad, declination_rad)
    # t increases with apparent motion, hence RA is division_ra - t.
    first_derivative = (
        cos_declination * math.sin(division_ra_rad),
        -cos_declination * math.cos(division_ra_rad),
        0.0,
    )
    second_derivative = (
        -cos_declination * math.cos(division_ra_rad),
        -cos_declination * math.sin(division_ra_rad),
        0.0,
    )
    third_derivative = tuple(-value for value in first_derivative)
    fourth_derivative = tuple(-value for value in second_derivative)
    division_derivatives = (
        division,
        first_derivative,
        second_derivative,
        third_derivative,
        fourth_derivative,
    )
    intersections = tuple(
        _cross(ecliptic_normal, _cross(north, derivative))
        for derivative in division_derivatives
    )

    selected = next(
        (vector for vector in intersections if _dot(vector, vector) != 0.0),
        None,
    )
    if selected is None:
        raise CalculationError("Makransky house and ecliptic planes coincide")

    # Derivatives of dot(Q(t), S(t)) at t=0. The first nonzero derivative
    # determines the sign immediately after the exact division in increasing
    # apparent motion, including a zero raw intersection.
    selector_derivatives = []
    for order in range(5):
        selector_derivatives.append(
            math.fsum(
                math.comb(order, left_order)
                * _dot(
                    intersections[left_order],
                    division_derivatives[order - left_order],
                )
                for left_order in range(order + 1)
            )
        )
    selector = next((value for value in selector_derivatives if value != 0.0), None)
    if selector is None:
        raise CalculationError("Makransky cusp has no unique oriented limit")
    if selector < 0.0:
        selected = tuple(-value for value in selected)  # type: ignore[assignment]
    return _ecliptic_lon_deg(selected, eps_rad)


def houses_sunshine_makransky(
    armc: float, lat: float, eps: float, asc: float, mc: float, sun_dec: float
) -> List[float]:
    """Construct Makransky's twelve Sunshine house cusps.

    Args:
        armc: Right ascension of the upper meridian, in degrees.
        lat: Signed geographic latitude, in degrees.
        eps: Obliquity, in degrees.
        asc: Supplied Ascendant anchor, in degrees.
        mc: Supplied Midheaven anchor, in degrees.
        sun_dec: Solar declination, in degrees.

    Returns:
        A mutable 13-element list with index zero unused.

    Raises:
        ValueError: If numeric inputs are non-finite or outside their domains.
        PolarCircleError: If the observer is at a geographic pole or the Sun is
            circumpolar.
        CalculationError: If a cusp plane has no unique oriented intersection.
    """
    values = (armc, lat, eps, asc, mc, sun_dec)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("Makransky construction requires finite angles")
    if not -90.0 <= lat <= 90.0:
        raise ValueError("Makransky latitude must be in [-90, 90] degrees")
    if not 0.0 <= eps < 90.0:
        raise ValueError("Makransky obliquity must be in [0, 90) degrees")
    if not -90.0 < sun_dec < 90.0:
        raise ValueError("Makransky solar declination must be in (-90, 90) degrees")
    if abs(lat) == 90.0:
        raise PolarCircleError(
            "Makransky horizon axis is undefined at a geographic pole",
            latitude=lat,
            house_system="i",
        )

    tangent_product = math.tan(math.radians(sun_dec)) * math.tan(math.radians(lat))
    if abs(tangent_product) >= 1.0:
        raise PolarCircleError(
            f"Sunshine houses undefined: the Sun (declination {sun_dec:.4f}) "
            f"is circumpolar at latitude {lat:.4f}",
            latitude=lat,
            house_system="i",
        )
    ascensional_difference = math.degrees(math.asin(tangent_product))
    night_third = (90.0 - ascensional_difference) / 3.0
    day_third = (90.0 + ascensional_difference) / 3.0

    southern = lat < 0.0
    working_armc = armc + (180.0 if southern else 0.0)
    lower_meridian = working_armc + 180.0
    output_rotation = 180.0 if southern else 0.0
    eps_rad = math.radians(eps)

    cusps = [0.0] * 13
    cusps[1] = float(asc % 360.0)
    cusps[4] = float((mc + 180.0) % 360.0)
    cusps[7] = float((asc + 180.0) % 360.0)
    cusps[10] = float(mc % 360.0)

    for step in (1, 2):
        for cusp, meridian, displacement in (
            (4 - step, lower_meridian, -step * night_third),
            (4 + step, lower_meridian, step * night_third),
            (10 - step, working_armc, -step * day_third),
            (10 + step, working_armc, step * day_third),
        ):
            longitude = _makransky_intermediate_cusp(
                meridian + displacement,
                working_armc,
                lat,
                sun_dec,
                eps_rad,
            )
            cusps[cusp] = float((longitude + output_rotation) % 360.0)
    return cusps


def apc_cusp(index: int, lat_rad: float, eps_rad: float, armc_rad: float) -> float:
    """One APC house cusp (system letter ``Y``).

    The APC construction is the Sunshine construction applied to the
    Ascendant's own parallel of declination: the parallel's diurnal and
    nocturnal arcs are divided into six equal parts each and the division
    points are carried to the ecliptic along house circles through the north
    and south points of the horizon. The Ascendant lies on its parallel at
    the rising hour angle, so cusp 1 reproduces the Ascendant and cusp 10's
    house circle is the meridian.

    Source-to-vector translation:
        Ingmar de Boer's independent construction description divides the
        Ascendant Parallel Circle into six parts above and six below the
        horizon, then carries each division point through a plane of oblique
        ascension to the ecliptic. In the equatorial frame used here, the APC
        is simply the small circle ``declination = declination(Ascendant)``;
        its above/below-horizon arc endpoints have hour angles equal to the
        rising and setting hour angles. ``_semiarc_cusp_hour_angles`` performs
        those two trisects, and a plane of oblique ascension is the great
        circle through the selected point and the north/south horizon axis.
        Intersecting that plane with the ecliptic is therefore exactly the
        construction, expressed without importing another implementation's
        trigonometric case table.

    Args:
        index: House cusp number, 1-12.
        lat_rad: Geographic latitude in radians.
        eps_rad: True obliquity in radians.
        armc_rad: Right ascension of the midheaven in radians.

    Returns:
        Cusp longitude in degrees in [0, 360).
    """
    zenith, east, north = _horizon_frame(armc_rad, lat_rad)
    ecl_pole = (0.0, -math.sin(eps_rad), math.cos(eps_rad))
    # The two intersections of the ecliptic and horizon are antipodal.  The
    # cross-product order alone does not identify the rising intersection at
    # polar latitudes, so select the point in the eastern half of the horizon.
    asc_vec = _normalized(_cross(ecl_pole, zenith))
    if _dot(asc_vec, east) < 0.0:
        asc_vec = (-asc_vec[0], -asc_vec[1], -asc_vec[2])
    dec_asc = math.asin(min(1.0, max(-1.0, asc_vec[2])))
    ra_asc = math.atan2(asc_vec[1], asc_vec[0])
    # Hour angle of the Ascendant, wrapped to (-pi, pi]: it is rising, so
    # its negative is the parallel's diurnal semi-arc.
    ha_asc = math.atan2(math.sin(armc_rad - ra_asc), math.cos(armc_rad - ra_asc))
    day_semiarc = -ha_asc if ha_asc < 0.0 else 2.0 * math.pi - ha_asc

    angles = dict(_semiarc_cusp_hour_angles(day_semiarc))
    angles[1] = -day_semiarc
    angles[7] = day_semiarc
    angles[10] = 0.0
    angles[4] = math.pi
    anchor = _unit_radec(armc_rad - angles[index], dec_asc)
    return _house_circle_lon_deg(anchor, north, eps_rad)
