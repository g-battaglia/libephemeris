# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Mean lunar node and apogee from IERS fundamental arguments.

This module contains no fitted ephemeris samples and reads no model artifact.
The angular model is built from the ERFA implementations of the IERS 2003
Delaunay arguments ``l``, ``F``, and ``Omega``.  The mean apogee is the point
opposite the mean perigee in a conventional inclined lunar orbit.

ERFA is the BSD-licensed open-source implementation of the IAU SOFA standards.
The argument definitions are published in IERS Conventions (2010), Eq. 5.43:
https://iers-conventions.obspm.fr/content/chapter5/icc5.pdf

Provenance:
    ERFA supplies the published Delaunay arguments; this module combines them
    through explicit mean-node/perigee geometry. Simon et al. (1994) supplies
    the Earth-eccentricity polynomial used by the interpolated model, while the
    IAU 2012 astronomical unit and stated conventional mean lunar orbit supply
    display distances. The half-day derivative interval and range delegation
    are project API choices, not fitted coefficients.
"""

from __future__ import annotations

import math

import erfa

from .constants import MEAN_APOG, MEAN_NODE
from .exceptions import EphemerisRangeError

_J2000_JD = 2451545.0
_JULIAN_CENTURY_DAYS = 36525.0
_RATE_HALF_SPAN_DAYS = 0.5

# Mean Earth-orbit eccentricity through T^2 from Meeus, Astronomical Algorithms
# (2nd ed., 1998), eq. 25.4 (after Simon et al. 1994, A&A 282, 663--683). The
# interpolated-apse harmonic basis uses its ratio to the J2000 value for
# solar-anomaly terms.
EARTH_ECCENTRICITY_J2000 = 0.016708634
EARTH_ECCENTRICITY_T = -0.000042037
EARTH_ECCENTRICITY_T2 = -0.0000001267

# Conventional mean lunar orbit. The mean distance 384400 km and mean
# eccentricity 0.0549 are the conventional values (Meeus, Astronomical
# Algorithms, ch. 47; Astronomical Almanac). Distances are then expressed in
# AU using the exact IAU 2012 (Resolution B2) astronomical unit.
_AU_KM = 149_597_870.7
_MEAN_LUNAR_DISTANCE_KM = 384_400.0
_MEAN_LUNAR_ECCENTRICITY = 0.0549
# Mean inclination of the lunar orbit to the ecliptic. Chapront-Touzé &
# Chapront, "Lunar Tables and Programs from 4000 B.C. to A.D. 8000" (ELP
# 2000-82), give the mean inclination as 5.145396°. The rounded 5.145° value
# left a constant ~1.000077 latitude-scale error (dlat up to ~1.43").
INCLINATION_DEG = 5.145396
MEAN_NODE_DISTANCE_AU = _MEAN_LUNAR_DISTANCE_KM / _AU_KM
MEAN_APOGEE_DISTANCE_AU = (
    _MEAN_LUNAR_DISTANCE_KM * (1.0 + _MEAN_LUNAR_ECCENTRICITY) / _AU_KM
)


def _active_ephemeris_range() -> tuple[str | None, float, float]:
    """Return loaded JPL coverage without forcing a kernel in LEB mode.

    The IERS formula is analytic and is not sampled from a LEB file, so a LEB
    file's own interval must not narrow its validity.  When a JPL kernel is
    already loaded, retain the public direct-evaluator range contract.  When
    LEB is the only active source, leave range enforcement to the surrounding
    ephemeris pipeline and avoid loading a second, unused kernel.
    """
    from .state import get_calc_mode, get_current_file_data, get_leb_reader, get_planets

    path, start_jd, end_jd, _ = get_current_file_data(0)
    # Honor the range only when it belongs to a JPL kernel. Since the
    # provenance fix that made get_current_file_data() report the concrete
    # LEB artifact serving body 0 (with its own interval), this branch also
    # sees LEB paths — and a LEB file's interval must not narrow the
    # analytic validity (the contract above): in auto mode the router falls
    # back to Skyfield/JPL for dates outside every LEB interval, so
    # enforcing the LEB window here rejected requests the pipeline could
    # serve. LEB-derived paths fall through to the mode branch below.
    _is_leb_source = bool(path) and path.endswith((".leb", ".leb2"))
    if start_jd and end_jd and not _is_leb_source:
        return path, float(start_jd), float(end_jd)

    mode = get_calc_mode()

    # Horizons mode has no local-kernel range contract.  These mean points are
    # evaluated from the analytic IERS arguments above, so asking for their
    # range must not instantiate Skyfield's Loader (which would otherwise
    # download de440.bsp merely to answer an unrelated range question).
    if mode == "horizons":
        return None, 0.0, 0.0

    if mode in ("auto", "leb"):
        try:
            reader = get_leb_reader()
        except RuntimeError:
            # A LEB generator can evaluate these standards-based points before
            # its output file exists.  In that bootstrap case there is no
            # active LEB coverage to inspect, so use the configured JPL kernel.
            reader = None
        if reader is not None:
            return str(reader.path), 0.0, 0.0

    get_planets()
    path, start_jd, end_jd, _ = get_current_file_data(0)
    return path, float(start_jd), float(end_jd)


def _check_range(jd: float, body_id: int, body_name: str) -> None:
    """Validate a TT Julian day against the active ephemeris coverage."""
    if not math.isfinite(jd):
        raise EphemerisRangeError(
            f"Julian day {jd} is not finite",
            requested_jd=jd,
            body_id=body_id,
            body_name=body_name,
        )
    path, start_jd, end_jd = _active_ephemeris_range()
    if start_jd and end_jd and not start_jd <= jd <= end_jd:
        raise EphemerisRangeError(
            f"Julian day {jd} is outside the active ephemeris range "
            f"{start_jd} to {end_jd}",
            requested_jd=jd,
            start_jd=start_jd,
            end_jd=end_jd,
            body_id=body_id,
            body_name=body_name,
            ephemeris_file=path,
        )


def _centuries(jd: float) -> float:
    """Convert a TT Julian date to Julian centuries since J2000.0."""
    return (jd - _J2000_JD) / _JULIAN_CENTURY_DAYS


def lunar_delaunay_arguments(jd_tt: float) -> tuple[float, float, float, float]:
    """Return the IERS 2003 lunar arguments ``(D, M, M', F)`` in radians.

    ERFA implements the IAU SOFA fundamental-argument routines published in
    IERS Conventions (2010), Eq. 5.43.  The explicit order used here matches
    the trigonometric basis of the independently generated interpolated-apse
    model: mean elongation, solar anomaly, lunar anomaly, and argument of
    latitude.
    """
    jd = float(jd_tt)
    centuries = _centuries(jd)
    return (
        float(erfa.fad03(centuries)),
        float(erfa.falp03(centuries)),
        float(erfa.fal03(centuries)),
        float(erfa.faf03(centuries)),
    )


def _fundamental_arguments(jd: float) -> tuple[float, float, float]:
    """Return ``(l, F, Omega)`` in radians from ERFA/IERS."""
    centuries = _centuries(jd)
    return (
        float(erfa.fal03(centuries)),
        float(erfa.faf03(centuries)),
        float(erfa.faom03(centuries)),
    )


def _angle_delta(after: float, before: float) -> float:
    """Return the signed principal angular difference in radians."""
    return math.atan2(math.sin(after - before), math.cos(after - before))


def _fundamental_argument_rates(jd: float) -> tuple[float, float, float]:
    """Return smooth IERS argument rates in radians per day.

    The centered one-day span is an independent numerical derivative of the
    standards routine.  It is far shorter than the 18.6-year nodal period and
    long enough to avoid cancellation in double precision.
    """
    before = _fundamental_arguments(jd - _RATE_HALF_SPAN_DAYS)
    after = _fundamental_arguments(jd + _RATE_HALF_SPAN_DAYS)
    span = 2.0 * _RATE_HALF_SPAN_DAYS
    return tuple(
        float(_angle_delta(after[index], before[index]) / span) for index in range(3)
    )  # type: ignore[return-value]


def _plane_node_unchecked(jd: float) -> tuple[float, float]:
    """Evaluate mean apogee-plane and node longitudes without range policy."""
    lunar_anomaly, argument_latitude, node = _fundamental_arguments(jd)
    apogee_plane = node + argument_latitude - lunar_anomaly + math.pi
    return (
        float(math.degrees(apogee_plane) % 360.0),
        float(math.degrees(node) % 360.0),
    )


def mean_lunar_apse_plane_node(jd_tt: float) -> tuple[float, float]:
    """Return mean-apogee plane longitude and mean-node longitude.

    The apogee plane longitude follows directly from the Delaunay identities
    ``F = L - Omega`` and ``l = L - varpi``: therefore
    ``varpi_apogee = Omega + F - l + pi``.
    """
    jd = float(jd_tt)
    _check_range(jd, MEAN_APOG, "mean lunar apogee")
    return _plane_node_unchecked(jd)


def _mean_apogee_cartesian_state(
    jd: float,
) -> tuple[float, float, float, float, float, float]:
    """Return the conventional mean-apogee Cartesian state of date."""
    lunar_anomaly, argument_latitude, node = _fundamental_arguments(jd)
    lunar_anomaly_rate, argument_latitude_rate, node_rate = _fundamental_argument_rates(
        jd
    )
    argument = argument_latitude - lunar_anomaly + math.pi
    argument_rate = argument_latitude_rate - lunar_anomaly_rate

    inclination = math.radians(INCLINATION_DEG)
    cos_i = math.cos(inclination)
    sin_i = math.sin(inclination)
    cos_node = math.cos(node)
    sin_node = math.sin(node)
    cos_argument = math.cos(argument)
    sin_argument = math.sin(argument)

    x_unit = cos_node * cos_argument - sin_node * sin_argument * cos_i
    y_unit = sin_node * cos_argument + cos_node * sin_argument * cos_i
    z_unit = sin_argument * sin_i

    dx_dnode = -y_unit
    dy_dnode = x_unit
    dx_dargument = -cos_node * sin_argument - sin_node * cos_argument * cos_i
    dy_dargument = -sin_node * sin_argument + cos_node * cos_argument * cos_i
    dz_dargument = cos_argument * sin_i

    distance = MEAN_APOGEE_DISTANCE_AU
    return (
        float(distance * x_unit),
        float(distance * y_unit),
        float(distance * z_unit),
        float(distance * (dx_dnode * node_rate + dx_dargument * argument_rate)),
        float(distance * (dy_dnode * node_rate + dy_dargument * argument_rate)),
        float(distance * dz_dargument * argument_rate),
    )


def _cartesian_to_polar_state(
    state: tuple[float, float, float, float, float, float],
) -> tuple[float, float, float, float, float, float]:
    """Convert Cartesian position/velocity to spherical state and rates.

    Inputs are ecliptic Cartesian AU and AU/day; outputs are longitude,
    latitude, distance, degree/day angular rates, and AU/day radial rate.
    """
    x, y, z, vx, vy, vz = state
    rho_squared = x * x + y * y
    rho = math.sqrt(rho_squared)
    distance = math.sqrt(rho_squared + z * z)
    radial_speed = (x * vx + y * vy + z * vz) / distance
    longitude_speed = math.degrees((x * vy - y * vx) / rho_squared)
    latitude_speed = math.degrees((vz * distance - z * radial_speed) / (distance * rho))
    return (
        float(math.degrees(math.atan2(y, x)) % 360.0),
        float(math.degrees(math.atan2(z, rho))),
        float(distance),
        float(longitude_speed),
        float(latitude_speed),
        0.0,
    )


def mean_lunar_apogee_position(jd_tt: float) -> tuple[float, float, float]:
    """Return mean lunar apogee longitude, latitude, and distance."""
    jd = float(jd_tt)
    _check_range(jd, MEAN_APOG, "mean lunar apogee")
    return _mean_lunar_apogee_position_unchecked(jd)


def _mean_lunar_apogee_position_unchecked(
    jd_tt: float,
) -> tuple[float, float, float]:
    """Evaluate the analytic IERS mean apogee without a kernel range check."""
    jd = float(jd_tt)
    return _cartesian_to_polar_state(_mean_apogee_cartesian_state(jd))[:3]


def mean_lunar_apogee_state(
    jd_tt: float, *, speed3: bool = False
) -> tuple[float, float, float, float, float, float]:
    """Return the standards-derived mean lunar apogee polar state.

    ``speed3`` is accepted for the shared backend call signature. This model
    returns its smooth standards-defined instantaneous derivative; public
    SPEED3 callers rebuild the derivative from fully transformed positions.
    """
    del speed3
    jd = float(jd_tt)
    _check_range(jd, MEAN_APOG, "mean lunar apogee")
    return _cartesian_to_polar_state(_mean_apogee_cartesian_state(jd))


def mean_lunar_node_position(jd_tt: float) -> tuple[float, float, float]:
    """Return mean lunar node longitude, zero latitude, and mean distance."""
    jd = float(jd_tt)
    _check_range(jd, MEAN_NODE, "mean lunar node")
    return _mean_lunar_node_position_unchecked(jd)


def _mean_lunar_node_position_unchecked(jd_tt: float) -> tuple[float, float, float]:
    """Evaluate the analytic IERS mean node without a kernel range check."""
    jd = float(jd_tt)
    _, _, node = _fundamental_arguments(jd)
    return (
        float(math.degrees(node) % 360.0),
        0.0,
        float(MEAN_NODE_DISTANCE_AU),
    )


def mean_lunar_node_state(
    jd_tt: float, *, speed3: bool = False
) -> tuple[float, float, float, float, float, float]:
    """Return the ERFA/IERS mean lunar node polar state."""
    del speed3
    jd = float(jd_tt)
    _check_range(jd, MEAN_NODE, "mean lunar node")
    _, _, node = _fundamental_arguments(jd)
    _, _, node_rate = _fundamental_argument_rates(jd)
    return (
        float(math.degrees(node) % 360.0),
        0.0,
        float(MEAN_NODE_DISTANCE_AU),
        float(math.degrees(node_rate)),
        0.0,
        0.0,
    )


__all__ = [
    "INCLINATION_DEG",
    "MEAN_APOGEE_DISTANCE_AU",
    "MEAN_NODE_DISTANCE_AU",
    "lunar_delaunay_arguments",
    "mean_lunar_apogee_position",
    "mean_lunar_apogee_state",
    "mean_lunar_apse_plane_node",
    "mean_lunar_node_position",
    "mean_lunar_node_state",
]
