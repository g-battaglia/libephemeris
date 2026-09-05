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

from .exceptions import PolarCircleError

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
    index: int,
    offset_deg: float,
    armc_deg: float,
    lat_deg: float,
    sun_dec_deg: float,
    eps_rad: float,
) -> float:
    """Compute one non-angular Sunshine cusp from Makransky's construction.

    Provenance and derivation
    -------------------------
    The implementation is a direct re-expression of Bob Makransky's own
    published spherical-trigonometry derivation, not of another ephemeris:

    * M. J. Makransky, *Primary Directions: A Primer of Calculation* (1988),
      pp. 36--39: definitions and derivation of meridian distance, zenith
      distance, pole, ascensional difference ``Q`` and ``W``;
    * M. J. Makransky, "The Solar House System" (November 1990; printed in
      *Mercury Hour*, April 1991), pp. 145--148, especially the cusp algorithm
      on pp. 147--148.

    Author-distributed scans preserved by the Internet Archive:

    * ``primarydirections_1.pdf`` (book pp. 1--48), capture 2016-03-08,
      SHA-256 ``7e8e6c078ae63888605a0377eac187efd5b111a470925cd1944cbd180f04973c``;
    * ``primarydirections_3.pdf`` (including pp. 145--148), same capture,
      SHA-256 ``bd0056bf966f307c3752eb2b177136f78378b88382b166df707a0172e94976d8``.

    The stable capture prefix is
    ``https://web.archive.org/web/20160308075443id_/http://dearbrutus.com/``.
    The scans are cited as evidence only and are not distributed by this
    project.

    Makransky defines every division point by an RA offset ``X`` from the
    upper meridian (day arc) or lower meridian (night arc). His general
    derivation defines meridian distance as the acute distance to the nearer
    upper or lower meridian and states its upper/lower character separately
    (pp. 15, 37--39). Ordinarily ``MD = |X|``. When a two-thirds division lies
    beyond the east/west hour circle, ``|X| > 90 degrees``; changing to the
    nearer meridian produces the complementary distance ``180 - |X|`` and:

    1. changes whether the point is referred to the upper or lower meridian;
    2. measures the acute prime-vertical arc from the 90-degree boundary,
       giving ``|A + F - 90 degrees|``; and
    3. reverses the sign used to form ``W = RA +/- Q``.

    These are the standard spherical-coordinate continuation of the printed
    acute-triangle algorithm across the east/west hour circle; Figure IV-4
    fixes the relevant quadrants. At exactly ``MD = 90 degrees`` the separate
    limiting formula printed on pp. 38 and 147 is used.

    Args:
        index: Cusp number, one of 2, 3, 5, 6, 8, 9, 11 or 12.
        offset_deg: Published RA offset ``X`` from the relevant meridian.
        armc_deg: Working right ascension of the upper meridian.  For southern
            latitudes the caller has already substituted RAIC, as instructed
            on p. 148.
        lat_deg: Signed geographic latitude.
        sun_dec_deg: Signed solar declination.
        eps_rad: Obliquity of the ecliptic, in radians.

    Returns:
        Ecliptic cusp longitude in degrees in ``[0, 360)``.
    """
    nocturnal = index in (2, 3, 5, 6)
    eastern_family = index in (2, 3, 11, 12)

    # Page 147 measures night-arc offsets from RAIC and day-arc offsets from
    # RAMC.  Keeping that distinction explicit makes the formula auditable.
    ra_deg = (armc_deg + (180.0 if nocturnal else 0.0) + offset_deg) % 360.0

    raw_md = abs(offset_deg)
    beyond_east_west = raw_md > 90.0
    md_deg = 180.0 - raw_md if beyond_east_west else raw_md

    # A day-arc point starts as an upper-meridian distance (UMD), while a
    # night-arc point starts as a lower-meridian distance (LMD).  Passing the
    # east/west hour circle swaps those roles (Makransky fig. IV-4).
    is_upper_distance = (not nocturnal) != beyond_east_west
    abs_lat_deg = abs(lat_deg)
    lat_rad = math.radians(lat_deg)
    abs_lat_rad = math.radians(abs_lat_deg)
    dec_rad = math.radians(sun_dec_deg)
    md_rad = math.radians(md_deg)

    if abs(md_deg - 90.0) < 1e-12:
        # Makransky's explicit limiting expression for a point on the hour
        # circle through the east and west horizon points (pp. 38 and 147).
        zenith_distance = 90.0 - math.degrees(
            math.atan(math.sin(abs_lat_rad) * math.tan(dec_rad))
        )
    else:
        # A and B are the two right-spherical-triangle angles on pp. 37--38.
        angle_a = math.degrees(math.atan(math.cos(lat_rad) * math.tan(md_rad)))
        angle_b = math.degrees(math.atan(math.tan(abs_lat_rad) * math.cos(md_rad)))

        # Makransky's four printed C rules collapse to one sign decision:
        # subtract |DEC| exactly when latitude/declination have the same sign
        # and MD is upper, or have opposite signs and MD is lower.
        same_declination_hemisphere = (sun_dec_deg >= 0.0) == (lat_deg >= 0.0)
        c_sign = -1.0 if same_declination_hemisphere == is_upper_distance else 1.0
        triangle_c = angle_b + c_sign * abs(sun_dec_deg)
        angle_f = math.degrees(
            math.atan(
                math.sin(abs_lat_rad)
                * math.sin(md_rad)
                * math.tan(math.radians(triangle_c))
            )
        )
        zenith_distance = angle_a + angle_f

        # For |X| > 90 degrees, A+F is measured past the east/west hour
        # circle.  The acute arc used by the pole formula is its distance from
        # that boundary, not the unreduced angle.
        if beyond_east_west:
            zenith_distance = abs(zenith_distance - 90.0)

    # Formula IV-2 and IV-1 (p. 39).  Clamps only absorb round-off at the
    # mathematical [-1, 1] limits; circumpolar input was rejected by caller.
    pole_argument = math.sin(lat_rad) * math.sin(math.radians(zenith_distance))
    pole_deg = math.degrees(math.asin(max(-1.0, min(1.0, pole_argument))))
    q_argument = math.tan(dec_rad) * math.tan(math.radians(pole_deg))
    q_deg = math.degrees(math.asin(max(-1.0, min(1.0, q_argument))))

    # Formula IV-3 uses RA-Q for the eastern cusp family and RA+Q for the
    # western family.  Beyond the east/west hour circle the celestial
    # quadrant is exchanged, so the sign reverses while the cusp-number
    # family (and therefore the final ecliptic branch) does not.
    q_sign = -1.0 if eastern_family else 1.0
    if beyond_east_west:
        q_sign = -q_sign
    w_deg = (ra_deg + q_sign * q_deg) % 360.0
    return _makransky_ecliptic_cut(
        w_deg,
        pole_deg,
        eps_rad,
        eastern_family=eastern_family,
    )


def houses_sunshine_makransky(
    armc: float, lat: float, eps: float, asc: float, mc: float, sun_dec: float
) -> List[float]:
    """Sunshine house cusps after Makransky (system letter ``i``).

    The Sun's diurnal arc is divided into six equal parts and its nocturnal
    arc into six equal parts; the division points, taken on the Sun's
    parallel of declination, are carried to the ecliptic along house circles
    through the north and south points of the horizon. The horizon and the
    meridian are themselves house circles, so cusps 1, 4, 7 and 10 are the
    four angles. Opposite intermediate cusps are generally not opposite in
    longitude because the diurnal and nocturnal semi-arcs differ.

    Args:
        armc: Right ascension of the midheaven, degrees.
        lat: Geographic latitude, degrees.
        eps: True obliquity, degrees.
        asc: Ascendant longitude, degrees (cusp 1).
        mc: Midheaven longitude, degrees (cusp 10).
        sun_dec: Declination of the Sun, degrees.

    Returns:
        Thirteen floats; index 0 unused, 1-12 are the cusp longitudes.

    Raises:
        PolarCircleError: If the Sun is circumpolar at the given latitude,
            leaving no sunrise or sunset to divide the day by.
    """
    # Two regular limiting cases deserve an explicit note because they are
    # easy to turn accidentally into quadrant special cases:
    #
    # * DEC = 0: the solar parallel is the celestial equator.  Makransky's
    #   ascensional difference is exactly zero, hence NSA = DSA = 90 degrees
    #   and both halves are trisected at 30/60 degrees.  In the subsequent
    #   right-spherical triangles, Q tends continuously to zero.  Nothing in
    #   the published construction changes cusp family or ecliptic branch as
    #   DEC crosses zero.
    # * LAT = 0: the zenith is on the celestial equator.  The pole of each
    #   house circle is zero, so the construction is independent of solar
    #   declination and reduces to carrying the twelve 30-degree RA divisions
    #   to the ecliptic.  The unequal longitude gaps are solely the ordinary
    #   RA-to-ecliptic projection caused by obliquity.
    #
    # Consequently exact zero is evaluated by the same equations as nearby
    # values.  Do not branch on ``sun_dec == 0`` or on the sign bit of a
    # floating-point zero: either would introduce a non-geometric 180-degree
    # jump at the equinox.  Source-neutral regression tests below the public
    # API verify both limits and the cyclic ordering of all twelve cusps.
    tan_product = math.tan(math.radians(sun_dec)) * math.tan(math.radians(lat))
    if abs(tan_product) >= 1.0:
        raise PolarCircleError(
            f"Sunshine houses undefined: the Sun (declination {sun_dec:.4f}) "
            f"is circumpolar at latitude {lat:.4f}",
            latitude=lat,
            house_system="i",
        )
    # Makransky p. 147: AD = asin(tan(DEC) tan(LAT)), then
    # NSA = 90 - AD and DSA = 90 + AD.  The circumpolar guard above ensures
    # the asin argument is strictly inside its domain.
    ascensional_difference = math.degrees(math.asin(tan_product))
    nocturnal_semiarc = 90.0 - ascensional_difference
    diurnal_semiarc = 90.0 + ascensional_difference
    night_third = nocturnal_semiarc / 3.0
    day_third = diurnal_semiarc / 3.0

    # Page 148 instructs southern latitudes to use RAIC in place of RAMC and
    # rotate the resulting non-angular cusps by 180 degrees.  This is an exact
    # coordinate transformation; the angular cusps supplied by the dispatcher
    # are already hemisphere-correct and therefore are not rotated here.
    southern = lat < 0.0
    working_armc = (armc + 180.0) % 360.0 if southern else armc
    output_rotation = 180.0 if southern else 0.0
    eps_rad = math.radians(eps)

    cusps = [0.0] * 13
    cusps[1] = asc
    cusps[7] = (asc + 180.0) % 360.0
    cusps[10] = mc
    cusps[4] = (mc + 180.0) % 360.0
    # Exact offset table printed on p. 147.  It contains the eight
    # intermediate cusps only; the four angles were assigned above.
    offsets = (
        (2, -2.0 * night_third),
        (3, -night_third),
        (5, night_third),
        (6, 2.0 * night_third),
        (8, -2.0 * day_third),
        (9, -day_third),
        (11, day_third),
        (12, 2.0 * day_third),
    )
    for cusp, offset in offsets:
        longitude = _makransky_intermediate_cusp(
            cusp,
            offset,
            working_armc,
            lat,
            sun_dec,
            eps_rad,
        )
        cusps[cusp] = (longitude + output_rotation) % 360.0
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
