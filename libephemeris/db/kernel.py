# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Portable coefficient evaluation and read planning with explicit inputs.

No classes, decorators, ambient state, I/O, callbacks, caches or driver objects
belong to this kernel. Functions consume explicit records, binary64 sequences
and maps/vectors; PostgreSQL transport and engine adapters are independent.

Provenance:
    Project-native storage adaptation over IEEE-754 LEB coefficients. Clenshaw
    (1955) recurrences are reused unchanged from leb_reader; no fit or new
    astronomical model is introduced. Arithmetic order matches native readers.
"""

from __future__ import annotations

import math

from ..exceptions import DBDataError
from ..constants import (
    EARTH,
    JUPITER,
    SATURN,
    SUN,
    FLG_BARYCTR,
    FLG_HELCTR,
    FLG_J2000,
    FLG_NOGDEFL,
    FLG_NONUT,
    FLG_TRUEPOS,
)
from ..leb_format import (
    COORD_ECLIPTIC,
    COORD_GEO_ECLIPTIC,
    COORD_HELIO_ECL_RETIRED,
)
from ..leb_reader import _clenshaw, _clenshaw_with_derivative
from .contract import NUTATION_ID, SegmentKey, Series, segment_index


def segment_coordinates(series: Series, jd: float) -> tuple[int, float]:
    """Calculate the segment index and normalized polynomial epoch.

    Args:
        series: Coverage and interval metadata.
        jd: Evaluation epoch in Julian days TT.

    Returns:
        Index and normalized time in [-1, 1]. Preserve this arithmetic order
        across implementations; reassociation changes rounded boundaries.

    Raises:
        DBDataError: Metadata is invalid.
        EphemerisRangeError: Epoch is outside inclusive coverage or non-finite.
    """
    index = segment_index(series, jd)
    start = series.jd_start + index * series.interval_days
    midpoint = start + 0.5 * series.interval_days
    tau = 2.0 * (jd - midpoint) / series.interval_days
    return index, max(-1.0, min(1.0, tau))


def evaluate_body(
    series: Series, coefficients: tuple[float, ...], tau: float
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Evaluate three body components using the existing Clenshaw recurrences.

    Args:
        series: Validated body metadata with exactly three components.
        coefficients: Validated component-major coefficients, passed explicitly.
        tau: Normalized epoch from segment_coordinates.

    Returns:
        Position and velocity in native series units. Angular longitude uses
        Euclidean modulo, not signed remainder.
    """
    width = series.degree + 1
    positions = [0.0, 0.0, 0.0]
    velocities = [0.0, 0.0, 0.0]
    for component in range(3):
        offset = component * width
        value, derivative = _clenshaw_with_derivative(
            coefficients[offset : offset + width], tau
        )
        positions[component] = value
        velocities[component] = derivative * (2.0 / series.interval_days)
    if series.coord_type in (
        COORD_ECLIPTIC,
        COORD_HELIO_ECL_RETIRED,
        COORD_GEO_ECLIPTIC,
    ):
        positions[0] %= 360.0
    return (
        (positions[0], positions[1], positions[2]),
        (velocities[0], velocities[1], velocities[2]),
    )


def evaluate_nutation(
    series: Series, coefficients: tuple[float, ...], tau: float
) -> tuple[float, float]:
    """Evaluate the two stored nutation components without a frame cache.

    Args:
        series: Validated nutation metadata with exactly two components.
        coefficients: Validated component-major coefficients.
        tau: Normalized epoch from segment_coordinates.

    Returns:
        Longitude and obliquity nutation in radians.
    """
    width = series.degree + 1
    return (
        _clenshaw(coefficients[:width], tau),
        _clenshaw(coefficients[width : 2 * width], tau),
    )


def interpolate_delta_t(points: list[tuple[float, float]], jd: float) -> float:
    """Interpolate at most two explicit samples with native endpoint clamping.

    Args:
        points: Ordered bracketing (Julian day, Delta-T days) samples. Outside
            the table's coverage the transport supplies one endpoint sample.
        jd: Requested epoch.

    Returns:
        Delta-T in days, without retaining either samples or the result.

    Raises:
        ValueError: Optional samples are absent or the epoch is non-finite.
        DBDataError: Samples are malformed, non-finite or not ordered.
    """
    if not points:
        raise ValueError("No Delta-T data in this DB dataset")
    if not math.isfinite(jd):
        raise ValueError("Non-finite Delta-T epoch")
    if len(points) > 2 or not all(
        math.isfinite(date) and math.isfinite(days) for date, days in points
    ):
        raise DBDataError("Invalid DB Delta-T samples")
    if len(points) == 1:
        return float(points[0][1])
    (start, before), (end, after) = points
    if end <= start:
        raise DBDataError("DB Delta-T samples are not strictly ordered")
    return float(before + (jd - start) / (end - start) * (after - before))


def required_segment_keys(
    series_by_body: dict[int, Series], jd_tt: float, body_id: int, flags: int
) -> list[SegmentKey]:
    """Plan likely target, observer, deflector and frame inputs without I/O.

    This is only a prefetch plan. Evaluation must still perform an exact lookup
    when an iteration asks for another epoch. A skipped dependency here never
    authorizes extrapolation or fallback to another persistent source.

    Args:
        series_by_body: Metadata owned by the current operation.
        jd_tt: Observation epoch already converted to TT.
        body_id: Requested public body identifier.
        flags: Normalized calculation flags.

    Returns:
        Deterministically ordered body/index pairs including neighboring
        segments. Keys are explicit signed-body/unsigned-index pairs.
    """
    body_ids = {body_id}
    if flags & FLG_HELCTR:
        body_ids.add(SUN)
    elif not flags & FLG_BARYCTR:
        body_ids.add(EARTH)
    if not flags & (FLG_NOGDEFL | FLG_TRUEPOS | FLG_BARYCTR | FLG_HELCTR):
        body_ids.update((SUN, JUPITER, SATURN))
    if not flags & (FLG_NONUT | FLG_J2000):
        body_ids.add(NUTATION_ID)

    keys = []
    for dependency in sorted(body_ids):
        series = series_by_body.get(dependency)
        if series is None or not series.jd_start <= jd_tt <= series.jd_end:
            continue
        index = segment_index(series, jd_tt)
        # Check before subtracting/adding: unsigned segment indices must not
        # underflow at zero or overflow at the final segment.
        if index > 0:
            keys.append((dependency, index - 1))
        keys.append((dependency, index))
        if index < series.segment_count - 1:
            keys.append((dependency, index + 1))
    return keys
