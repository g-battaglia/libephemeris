# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
The shadow geometry an eclipse or an occultation is made of.

One spherical light source of radius ``R_s`` casts, past one spherical
occulter of radius ``R_o`` whose centre is ``D`` away, two coaxial cones of
revolution: an inner one, tangent internally to both spheres, and an outer
one, tangent externally. A third body of radius ``R_b`` then sits somewhere
inside or beside them. Three settings of the module share that one picture --
a solar eclipse (Sun, Moon, Earth), an occultation of a planet or a star by
the Moon (the occulted object, Moon, Earth) and a lunar eclipse (Sun, Earth,
Moon) -- and everything a caller needs afterwards is a handful of numbers
about the two cones: how far the axis passes from the shadowed body's centre,
how wide each cone is where it matters, how steep each cone is, and how big
the shadowed body is.

:class:`ShadowGeometry` is the typed carrier of exactly those numbers. It is
a carrier and not a stage of the computation: a producer fills it from
quantities it has already finished computing and never reads a field back to
continue its own arithmetic. Every length in it is in kilometres, the unit is
part of the field name, and the two half-angle sines are derived from the
stored cosines instead of being stored beside them, so that a stored pair can
never drift out of step with its cosine.

Sign convention: the inner cone's diameters are **signed**, negative while
the point of evaluation lies short of the cone's apex (the umbra proper,
hence a total eclipse) and positive beyond it (the antumbra, hence an annular
eclipse). That is the classical sign of the Besselian element ``l2`` and it is
load-bearing: five public entry points publish the total/annular distinction
through it, so a producer that returned an unsigned width would destroy it.

References:
    - Chauvenet, W., "A Manual of Spherical and Practical Astronomy", vol. I
      (5th ed. 1891; repr. Dover 1960), the chapter on eclipses -- the
      geometry of the internally and externally tangent cones.
    - Explanatory Supplement to the Astronomical Almanac, 3rd ed. (2013),
      ch. 11 "Eclipses of the Sun and Moon" -- the Besselian elements, the
      fundamental plane and the tangency conditions.
    - Meeus, J., "Astronomical Algorithms", 2nd ed. (1998), ch. 54 -- the
      same elements in working form and the sign of the umbral radius on the
      fundamental plane.

Provenance:
    Project-authored data type. It carries the published two-cone eclipse
    geometry of Chauvenet and of the Explanatory Supplement (ch. 11) under
    field names that state the physical quantity and its unit; the sign
    convention of the inner cone is the classical Besselian one for the
    element l2. The module holds no ephemeris, no fitted number and no
    algorithm of its own: it declares fields, derives the two half-angle
    sines from the cosines stored beside them, and checks the geometric
    invariants a producer must satisfy. The values themselves are computed by
    libephemeris/eclipse.py from JPL-backed Sun, Moon and Earth states.
"""

from __future__ import annotations

import math
from dataclasses import KW_ONLY, dataclass

__all__ = ["ShadowGeometry"]


@dataclass(frozen=True, slots=True)
class ShadowGeometry:
    """The two shadow cones, the axis and the shadowed body, by name.

    Immutable and equal by value, so one record can be passed down a chain of
    helpers and held across a root-finding loop without any chance of a
    consumer editing it. Only the first field may be given positionally; every
    other one is keyword-only, and the record is deliberately not a sequence:
    it has no length, cannot be unpacked and cannot be indexed. A consumer
    names the field it wants.

    Every numeric attribute is a native Python float and the one flag a native
    bool. Every length is in kilometres -- there is no per-producer unit and no
    unit field.

    Attributes:
        axis_offset_km: Perpendicular distance, in kilometres, between the
            shadow axis and the centre of the shadowed body: the Earth's
            centre for a solar eclipse or an occultation, the Moon's centre
            for a lunar eclipse. This is the classical Besselian
            ``sqrt(x**2 + y**2)`` carried in kilometres rather than in
            equatorial radii. Non-negative.
        umbral_plane_diameter_km: Diameter, in kilometres and **signed**, of
            the inner cone's section on the fundamental plane -- the plane
            through the centre of the shadowed body perpendicular to the
            shadow axis. Twice the classical Besselian ``l2``, with its sign:
            negative short of the cone's apex (umbra), positive beyond it
            (antumbra).
        penumbral_plane_diameter_km: Diameter, in kilometres, of the outer
            cone's section on the same plane: twice the classical Besselian
            ``l1``. Always positive, the outer cone having no apex on the
            shadowed body's side.
        cos_umbral_half_angle: Cosine of the half-angle of the inner cone,
            the one tangent internally to the source and to the occulter.
            Dimensionless, in (0, 1].
        cos_penumbral_half_angle: Cosine of the half-angle of the outer cone,
            the one tangent externally to both. Dimensionless, in (0, 1], and
            never larger than the inner cone's cosine.
        shadowed_radius_km: Radius, in kilometres, of the body that receives
            the shadow: the Earth's **equatorial** radius for a solar eclipse
            or an occultation, the Moon's mean radius for a lunar eclipse.
            Paired with a cone cosine it gives the half-width of that cone at
            the body's limb, which is what a contact condition is made of.
            Positive.
        shadow_misses_body: True when no part of the shadow -- the penumbra
            included -- reaches the shadowed body at this instant, so that the
            surface fields below describe a closest-approach extrapolation
            rather than a real contact.
        umbral_surface_diameter_km: Diameter, in kilometres and signed by the
            same rule as ``umbral_plane_diameter_km``, of the inner cone's
            section at the point where the shadow axis meets the shadowed
            body's surface (or, when the axis misses the body, at the surface
            point of closest approach). ``None`` when the producer resolves no
            such point, as the lunar producer does not.
        penumbral_surface_diameter_km: Diameter, in kilometres, of the outer
            cone's section at that same surface point. Positive, and set
            together with ``umbral_surface_diameter_km`` or not at all.

    Note:
        The two surface diameters are plane sections of the cones taken
        perpendicular to the axis at the surface point, not the oblique
        footprint the shadow paints on the ground: the ellipse actually drawn
        on the surface is wider whenever the source stands away from the
        zenith, and the module has a separate public function for the path
        width proper.

        For a lunar eclipse the Earth's shadow is not a geometric cone -- the
        atmosphere enlarges it. The record takes no position on which
        enlargement convention produced the sections it carries: a consumer
        cannot tell one from the record and must not try to undo or reapply
        one.
    """

    axis_offset_km: float
    _: KW_ONLY
    umbral_plane_diameter_km: float
    penumbral_plane_diameter_km: float
    cos_umbral_half_angle: float
    cos_penumbral_half_angle: float
    shadowed_radius_km: float
    shadow_misses_body: bool
    umbral_surface_diameter_km: "float | None" = None
    penumbral_surface_diameter_km: "float | None" = None

    @property
    def sin_umbral_half_angle(self) -> float:
        """Sine of the inner cone's half-angle, derived from its cosine.

        Both half-angles are acute, so the sine is the non-negative root of
        ``1 - cos**2``. It is derived rather than stored precisely so that it
        can never disagree with the cosine beside it.
        """
        cosine = self.cos_umbral_half_angle
        return math.sqrt(max(0.0, 1.0 - cosine * cosine))

    @property
    def sin_penumbral_half_angle(self) -> float:
        """Sine of the outer cone's half-angle, derived from its cosine."""
        cosine = self.cos_penumbral_half_angle
        return math.sqrt(max(0.0, 1.0 - cosine * cosine))

    def validate(self) -> None:
        """Check the geometric invariants a producer must satisfy.

        Deliberately not called from ``__post_init__``: a contact search
        builds thousands of records per contact time and the hot path stays
        free of checks. Tests and any producer that wants a loud failure call
        it explicitly.

        The invariants are: the outer cone's sections are real widths; the
        outer cone opens at least as wide as the inner one, the two being
        equal exactly when the light source is point-like; the axis offset is
        a distance; the shadowed body has a positive radius; and the two
        surface diameters are set together or not at all.

        Raises:
            ValueError: If any of those does not hold.
        """
        if not self.penumbral_plane_diameter_km > 0.0:
            raise ValueError(
                "the penumbral cone has no apex on the shadowed body's side, so "
                "its section on the fundamental plane must be a positive width, "
                f"not {self.penumbral_plane_diameter_km!r}"
            )
        if not 0.0 < self.cos_penumbral_half_angle <= self.cos_umbral_half_angle <= 1.0:
            raise ValueError(
                "the outer cone opens at least as wide as the inner one, so the "
                "cosines must satisfy 0 < penumbral <= umbral <= 1, not "
                f"{self.cos_penumbral_half_angle!r} and "
                f"{self.cos_umbral_half_angle!r}"
            )
        if self.axis_offset_km < 0.0:
            raise ValueError(
                "the offset of the axis from the body's centre is a distance, "
                f"not {self.axis_offset_km!r}"
            )
        if not self.shadowed_radius_km > 0.0:
            raise ValueError(
                "the shadowed body must have a positive radius, not "
                f"{self.shadowed_radius_km!r}"
            )
        has_umbral = self.umbral_surface_diameter_km is not None
        has_penumbral = self.penumbral_surface_diameter_km is not None
        if has_umbral != has_penumbral:
            raise ValueError(
                "the two surface diameters describe one point and are set "
                "together or not at all"
            )
        if has_penumbral and not self.penumbral_surface_diameter_km > 0.0:  # type: ignore[operator]
            raise ValueError(
                "the penumbral section at the surface point must be a positive "
                f"width, not {self.penumbral_surface_diameter_km!r}"
            )
