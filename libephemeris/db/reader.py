# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Operation-owned reader adapting PostgreSQL records to the LEB engine.

Provenance:
    Storage adaptation over project-native LEB coefficients. Clenshaw (1955)
    polynomial recurrences are reused unchanged from leb_reader; no new fit is
    performed. Coordinate meanings and IEEE-754 payloads retain the registered
    JPL/IAU provenance documented in docs/db-backend.md.
"""

from __future__ import annotations

from typing import Any, Callable

from .contract import (
    CoefficientStore,
    NUTATION_ID,
    SegmentKey,
    Series,
    decode_segment,
    segment_byte_count,
    validate_series,
)
from .kernel import (
    evaluate_body,
    evaluate_nutation,
    interpolate_delta_t,
    required_segment_keys,
    segment_coordinates,
)
from ..exceptions import DBDataError, EphemerisRangeError
from ..leb_format import StarEntry
from ..operations import source_guard


class DBReader:
    """Python engine adapter; ownership is scoped, numerical work is in kernel.

    This class satisfies the existing engine reader interface. It is not a
    numerical-domain class hierarchy: operation-owned maps supply explicit
    arguments to the free functions in contract/kernel.
    """

    source = "DB"
    cacheable = False
    _manifest_verified = False
    tier: str | None = None
    artifact_groups: tuple[str, ...] = ()
    _input_budget: Callable[[int], None] | None = None
    MAX_INPUT_SEGMENTS = 256

    @source_guard
    def __init__(self, store: CoefficientStore, dataset_id: str) -> None:
        """Read dataset headers and establish ownership of transient inputs.

        Args:
            store: Language-neutral storage implementation.
            dataset_id: Explicit immutable dataset UUID.
        """
        self._store = store
        self.dataset_id = dataset_id
        self._series = store.metadata(dataset_id)
        for series in self._series.values():
            validate_series(series)
        self._bodies = {
            body_id: series
            for body_id, series in self._series.items()
            if body_id != NUTATION_ID
        }
        if not self._bodies:
            raise DBDataError("DB dataset contains no body coefficients")
        self._frame_cache: dict[tuple[int, float], Any] = {}
        self._segments: dict[SegmentKey, bytes] = {}
        self._closed = False

    @property
    def jd_range(self) -> tuple[float, float]:
        """Return the union of stored body coverage intervals.

        Returns:
            Inclusive start and end in Julian days TT.
        """
        return (
            min(series.jd_start for series in self._bodies.values()),
            max(series.jd_end for series in self._bodies.values()),
        )

    def has_body(self, body_id: int) -> bool:
        """Report whether a body has a stored channel.

        Args:
            body_id: Public body identifier.

        Returns:
            Whether metadata declares that channel.
        """
        return body_id in self._bodies

    def body_coverage(self, body_id: int) -> tuple[float, float] | None:
        """Read declared coverage without confusing support-body ranges.

        Args:
            body_id: Public body identifier.

        Returns:
            Inclusive TT bounds, or None for an absent channel.
        """
        series = self._bodies.get(body_id)
        return None if series is None else (series.jd_start, series.jd_end)

    def has_nutation(self) -> bool:
        """Report availability of a usable stored nutation series.

        Returns:
            Whether the dataset declares nutation coefficients.
        """
        return NUTATION_ID in self._series

    @source_guard
    def prepare(self, jd_tt: float, body_id: int, flags: int) -> None:
        """Batch likely inputs before the existing calculation pipeline runs.

        The fetch includes neighboring segments for nearby speed and retarded
        epochs. This is an optimization only: eval methods perform exact
        indexed reads if an iteration needs any segment outside this set.

        Args:
            jd_tt: Observation epoch already converted to TT.
            body_id: Requested public body identifier.
            flags: Normalized calculation flags.
        """
        keys = required_segment_keys(self._series, jd_tt, body_id, flags)
        # Large custom inventories still work through exact on-demand reads.
        # Keep long event searches bounded instead of accumulating a dataset.
        self._fetch(keys[: self.MAX_INPUT_SEGMENTS])

    def _fetch(self, keys: list[SegmentKey]) -> None:
        """Fetch missing inputs and verify every payload before retaining it.

        Args:
            keys: Body/segment pairs required by this operation.

        Raises:
            DBDataError: Reader is closed or payload is missing/malformed.
        """
        if self._closed:
            raise DBDataError("DB operation reader is closed")
        unique_keys = list(dict.fromkeys(keys))
        missing = [key for key in unique_keys if key not in self._segments]
        if len(self._segments) + len(missing) > self.MAX_INPUT_SEGMENTS:
            # Discard older operation inputs, not persisted coefficients.
            # Any later evaluation simply reads its exact segment again.
            self._segments.clear()
            missing = unique_keys
        if not missing:
            return
        budget = getattr(self, "_input_budget", None)
        if budget is not None:
            budget(len(missing))
            missing = [key for key in unique_keys if key not in self._segments]
        payloads = self._store.fetch_segments(self.dataset_id, missing)
        if set(payloads) != set(missing):
            raise DBDataError("Missing coefficient payload")
        for key, payload in payloads.items():
            if len(payload) != segment_byte_count(self._series[key[0]]):
                raise DBDataError("Invalid coefficient payload length")
        self._segments.update(payloads)

    def _coefficients(
        self, series: Series, jd: float
    ) -> tuple[tuple[float, ...], float]:
        """Resolve exact coefficients and normalized epoch without extrapolation.

        Args:
            series: Validated metadata for the requested channel.
            jd: Evaluation epoch in Julian days TT.

        Returns:
            Component-major coefficients and normalized time in [-1, 1].
        """
        try:
            index, tau = segment_coordinates(series, jd)
        except EphemerisRangeError:
            # Match the existing reader protocol. Public API guards turn
            # ordinary coverage misses into typed errors; optional deflectors
            # may skip uncovered channels just as they do with file inputs.
            raise ValueError(
                f"JD {jd} outside range [{series.jd_start}, {series.jd_end}] "
                f"for body {series.body_id}"
            ) from None
        key = (series.body_id, index)
        self._fetch([key])
        coefficients = decode_segment(series, self._segments[key])
        return coefficients, tau

    @source_guard
    def eval_body(self, body_id: int, jd: float) -> tuple[tuple, tuple]:
        """Evaluate position and analytic velocity using existing recurrences.

        Args:
            body_id: Public body identifier.
            jd: Evaluation epoch in Julian days TT.

        Returns:
            Position and velocity tuples in the series' native units.

        Raises:
            KeyError: Body channel is absent (the existing reader protocol).
            DBDataError: Persisted coefficients are malformed.
        """
        if self._closed:
            raise DBDataError("DB operation reader is closed")
        series = self._bodies.get(body_id)
        if series is None:
            raise KeyError(f"Body {body_id} not in DB dataset")
        coefficients, tau = self._coefficients(series, jd)
        return evaluate_body(series, coefficients, tau)

    @source_guard
    def eval_nutation(self, jd_tt: float) -> tuple[float, float]:
        """Evaluate stored nutation without a global frame/result cache.

        Args:
            jd_tt: Evaluation epoch in Julian days TT.

        Returns:
            Nutation in longitude and obliquity, in radians.

        Raises:
            ValueError: Optional nutation is absent, as in the native readers.
            DBDataError: A declared nutation payload is corrupt or missing.
        """
        if self._closed:
            raise DBDataError("DB operation reader is closed")
        series = self._series.get(NUTATION_ID)
        if series is None:
            raise ValueError("No nutation data in this DB dataset")
        coefficients, tau = self._coefficients(series, jd_tt)
        return evaluate_nutation(series, coefficients, tau)

    @source_guard
    def delta_t(self, jd: float) -> float:
        """Interpolate the nearest samples, matching LEB endpoint clamping.

        Args:
            jd: Epoch to evaluate.

        Returns:
            Delta-T in days.

        Raises:
            ValueError: Optional samples are absent, as in the native readers.
            DBDataError: A database read or persisted sample is invalid.
        """
        if self._closed:
            raise DBDataError("DB operation reader is closed")
        points = self._store.delta_t_points(self.dataset_id, jd)
        return interpolate_delta_t(points, jd)

    @source_guard
    def get_star(self, star_id: int) -> StarEntry:
        """Read one source-catalog entry without retaining it.

        Args:
            star_id: Source catalog identifier.

        Returns:
            A detached native LEB StarEntry.
        """
        if self._closed:
            raise DBDataError("DB operation reader is closed")
        return StarEntry(*self._store.star(self.dataset_id, star_id))

    def close(self) -> None:
        """Release all ephemeris inputs owned by this operation."""
        self._segments.clear()
        self._frame_cache.clear()
        self._series.clear()
        self._bodies.clear()
        self.artifact_groups = ()
        self._manifest_verified = False
        self.tier = None
        self._closed = True
