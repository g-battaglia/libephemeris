# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Language-neutral records, byte layout and validation functions.

Provenance:
    Project-native storage metadata. Payloads retain the IEEE-754 binary64 LEB
    component-major layout documented in docs/db-backend.md. Segment selection
    matches the native readers; no astronomical approximation is introduced.
"""

from __future__ import annotations

import math
import struct
from dataclasses import dataclass
from typing import Protocol

from ..exceptions import DBDataError, EphemerisRangeError

SCHEMA_VERSION = 1
NUTATION_ID = -1
SegmentKey = tuple[int, int]


@dataclass(frozen=True, slots=True)
class Series:
    """Immutable data-only record with explicit primitive fields.

    No behavior, inheritance or object serialization belongs to this record.
    Functions below receive an explicit Series argument. Body -1 denotes
    nutation.

    Attributes:
        body_id: Public body identifier (i32), or NUTATION_ID.
        coord_type: Source coordinate-channel identifier (u32).
        segment_count: Stored coefficient segment count (u32).
        jd_start: Inclusive coverage start, in Julian days TT (f64).
        jd_end: Inclusive coverage end, in Julian days TT (f64).
        interval_days: Actual polynomial interval, not its nominal setting (f64).
        degree: Highest Chebyshev degree (u32).
        components: Coordinate channel count: three for bodies, two for nutation.
    """

    body_id: int
    coord_type: int
    segment_count: int
    jd_start: float
    jd_end: float
    interval_days: float
    degree: int
    components: int


def series_fields(
    series: Series,
) -> tuple[int, int, int, float, float, float, int, int]:
    """Return fields in the documented storage-contract record order.

    Args:
        series: Source metadata record.

    Returns:
        Body, coordinate type, count, start, end, interval, degree and components.
        The order is explicit, not inferred by reflection.
    """
    return (
        series.body_id,
        series.coord_type,
        series.segment_count,
        series.jd_start,
        series.jd_end,
        series.interval_days,
        series.degree,
        series.components,
    )


def segment_byte_count(series: Series) -> int:
    """Calculate the exact byte length of one component-major segment.

    Args:
        series: Validated metadata record.

    Returns:
        Number of IEEE-754 binary64 coefficients multiplied by eight.
    """
    return (series.degree + 1) * series.components * 8


def validate_series(series: Series) -> None:
    """Reject metadata that cannot safely describe a polynomial series.

    Args:
        series: Untrusted metadata record to validate.

    Raises:
        DBDataError: Metadata is non-finite or violates format limits.
    """
    finite_values = (series.jd_start, series.jd_end, series.interval_days)
    invalid = (
        not all(math.isfinite(value) for value in finite_values)
        or series.jd_end <= series.jd_start
        or series.interval_days <= 0
        or series.segment_count <= 0
        or not 0 <= series.degree <= 256
        or not 0 <= series.coord_type <= 4
        or series.components != (2 if series.body_id == NUTATION_ID else 3)
    )
    if invalid:
        raise DBDataError("Invalid coefficient series metadata")

    span = series.jd_end - series.jd_start
    capacity = series.segment_count * series.interval_days
    last_start = (series.segment_count - 1) * series.interval_days
    if not all(math.isfinite(value) for value in (span, capacity, last_start)):
        raise DBDataError("Invalid coefficient series coverage")
    # The last fit may extend beyond the advertised endpoint, but every stored
    # segment must be reachable and the coverage cannot extend beyond the grid.
    # Allow only binary64 endpoint/product roundoff, not a relative date margin.
    roundoff = math.ulp(series.jd_start) + math.ulp(series.jd_end) + math.ulp(capacity)
    if span <= last_start or span > capacity + roundoff:
        raise DBDataError("Coefficient series coverage disagrees with segment grid")


def decode_segment(series: Series, payload: bytes) -> tuple[float, ...]:
    """Validate and decode one binary64 component-major payload.

    Args:
        series: Validated metadata describing the payload.
        payload: Exact little-endian coefficient bytes.

    Returns:
        Finite native floats in component/degree order. Decode each eight-byte
        word as little-endian binary64; no native-endian casts.

    Raises:
        DBDataError: Payload has an invalid length or non-finite values.
    """
    if len(payload) != segment_byte_count(series):
        raise DBDataError("Invalid coefficient payload length")
    count = series.components * (series.degree + 1)
    coefficients = struct.unpack(f"<{count}d", payload)
    if not all(math.isfinite(value) for value in coefficients):
        raise DBDataError("Non-finite DB polynomial coefficient")
    return coefficients


def segment_index(series: Series, jd: float) -> int:
    """Select a segment using the same arithmetic as the file readers.

    Args:
        series: Metadata describing inclusive coverage and actual intervals.
        jd: Evaluation epoch in Julian days TT.

    Returns:
        Zero-based segment index. The inclusive final endpoint belongs to the
        last segment, not to a nonexistent next segment.

    Raises:
        DBDataError: Series metadata is invalid.
        EphemerisRangeError: Epoch is non-finite or outside coverage.
    """
    validate_series(series)
    if not math.isfinite(jd) or not series.jd_start <= jd <= series.jd_end:
        raise EphemerisRangeError(
            "Date outside DB series coverage",
            requested_jd=jd,
            start_jd=series.jd_start,
            end_jd=series.jd_end,
            body_id=series.body_id,
        )
    index = int((jd - series.jd_start) / series.interval_days)
    return min(index, series.segment_count - 1)


class CoefficientStore(Protocol):
    """Transport interface, separate from the pure numerical kernel."""

    def metadata(self, dataset_id: str) -> dict[int, Series]:
        """Read metadata from a published dataset.

        Args:
            dataset_id: Immutable dataset UUID.

        Returns:
            Validated series indexed by body identifier.
        """
        ...

    def fetch_segments(
        self, dataset_id: str, keys: list[SegmentKey]
    ) -> dict[SegmentKey, bytes]:
        """Read exactly the requested coefficient payloads in one batch.

        Args:
            dataset_id: Immutable dataset UUID.
            keys: Body identifiers paired with zero-based segment indices.

        Returns:
            Payloads indexed by the requested keys.

        Raises:
            DBDataError: A required segment is missing.
        """
        ...

    def delta_t_points(self, dataset_id: str, jd: float) -> list[tuple[float, float]]:
        """Read the nearest Delta-T samples on each side of an epoch.

        Args:
            dataset_id: Immutable dataset UUID.
            jd: Epoch to bracket.

        Returns:
            Zero, one or two ordered (Julian day, Delta-T days) samples.
        """
        ...

    def star(self, dataset_id: str, star_id: int) -> tuple:
        """Read a single star without retaining the catalog.

        Args:
            dataset_id: Immutable dataset UUID.
            star_id: Source catalog identifier.

        Returns:
            Fields in the order defined by the LEB StarEntry record.
        """
        ...
