# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Segment-backed readers for externally supplied LEB coefficients.

Provenance:
    Evaluation uses the project's Clenshaw (1955) recurrence in ``leb_reader``
    and the project-native LEB metadata and endpoint conventions. Source
    identity, validation, and lifecycle are storage extension contracts, not
    astronomical models; providers supply generator-attested coefficients.
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from bisect import bisect_right
from collections.abc import Mapping, Sequence
from dataclasses import replace

from .exceptions import CoefficientSourceError
from .leb_format import (
    COORD_ECLIPTIC,
    COORD_GEO_ECLIPTIC,
    COORD_HELIO_ECL_RETIRED,
    BodyEntry,
    NutationHeader,
    StarEntry,
)
from .leb_groups import LEB2_GROUPS
from .leb_reader import _clenshaw, _clenshaw_with_derivative

SEGMENT_SOURCE_API_VERSION = 1


def _validate_grid(entry: BodyEntry | NutationHeader, components: int) -> None:
    """Reject invalid source metadata before it can become a fallback signal."""
    numbers = (entry.jd_start, entry.jd_end, entry.interval_days)
    if (
        not all(math.isfinite(value) for value in numbers)
        or entry.jd_end <= entry.jd_start
        or entry.interval_days <= 0
        or entry.segment_count <= 0
        or not 0 <= entry.degree <= 256
        or entry.components != components
    ):
        raise CoefficientSourceError("Invalid coefficient series metadata")
    span = entry.jd_end - entry.jd_start
    capacity = entry.segment_count * entry.interval_days
    last_start = (entry.segment_count - 1) * entry.interval_days
    if not all(math.isfinite(v) for v in (span, capacity, last_start)):
        raise CoefficientSourceError("Invalid coefficient segment grid")
    roundoff = math.ulp(entry.jd_start) + math.ulp(entry.jd_end) + math.ulp(capacity)
    if span <= last_start or span > capacity + roundoff:
        raise CoefficientSourceError("Coefficient range disagrees with segment grid")


class SegmentSource(ABC):
    """Evaluate external segments with the native LEB reader arithmetic.

    Args:
        artifact_name: Canonical ``{tier}_{group}.leb2`` artifact name.
        locator: Public, credential-free source identity (not a connection URL).
        jd_range: Original artifact header range, used for auxiliary selection.
        bodies: Native per-body metadata; file offsets are ignored.
        nutation: Optional native nutation metadata.
        delta_t: Parallel Julian Day and Delta-T (days) tables.
        stars: Optional native star records.
        reviewed: Whether the provider verified the installed manifest pins.

    Providers must raise ``CoefficientSourceError`` for transport failures,
    missing segments, or invalid payloads, never a range/fallback exception.
    ``close`` releases instance caches only; shared pools belong to the provider.
    """

    def __init__(
        self,
        *,
        artifact_name: str,
        locator: str,
        jd_range: tuple[float, float],
        bodies: Mapping[int, BodyEntry],
        nutation: NutationHeader | None = None,
        delta_t: tuple[Sequence[float], Sequence[float]] = ((), ()),
        stars: Mapping[int, StarEntry] | None = None,
        reviewed: bool = False,
    ) -> None:
        valid_names = {
            f"{tier}_{group}.leb2"
            for tier in ("base", "medium", "extended")
            for group in LEB2_GROUPS
        }
        if artifact_name not in valid_names:
            raise CoefficientSourceError("Invalid coefficient artifact name")
        if (
            len(jd_range) != 2
            or not all(math.isfinite(v) for v in jd_range)
            or jd_range[0] >= jd_range[1]
        ):
            raise CoefficientSourceError("Invalid coefficient artifact range")
        self.artifact_name = artifact_name
        self._path = f"{locator.rstrip('/')}/{artifact_name}"
        self._jd_range = tuple(float(v) for v in jd_range)
        self._bodies = {key: replace(value) for key, value in bodies.items()}
        for key, entry in self._bodies.items():
            _validate_grid(entry, 3)
            if key != entry.body_id or not 0 <= entry.coord_type <= 4:
                raise CoefficientSourceError("Invalid coefficient body metadata")
        self._nutation = replace(nutation) if nutation is not None else None
        if self._nutation is not None and self._nutation.segment_count:
            _validate_grid(self._nutation, 2)
        self._delta_t_jds = list(delta_t[0])
        self._delta_t_vals = list(delta_t[1])
        if (
            len(self._delta_t_jds) != len(self._delta_t_vals)
            or not all(math.isfinite(v) for v in self._delta_t_jds)
            or not all(math.isfinite(v) for v in self._delta_t_vals)
            or any(a >= b for a, b in zip(self._delta_t_jds, self._delta_t_jds[1:]))
        ):
            raise CoefficientSourceError("Invalid coefficient Delta-T table")
        self._stars = {key: replace(value) for key, value in (stars or {}).items()}
        for key, star in self._stars.items():
            if key != star.star_id or not all(
                math.isfinite(v)
                for v in (
                    star.ra_j2000,
                    star.dec_j2000,
                    star.pm_ra,
                    star.pm_dec,
                    star.parallax,
                    star.rv,
                    star.magnitude,
                )
            ):
                raise CoefficientSourceError("Invalid coefficient star metadata")
        self._manifest_verified = bool(reviewed)
        self._closed = False
        self._eval_cache: dict[
            tuple[int, float],
            tuple[tuple[float, float, float], tuple[float, float, float]],
        ] = {}

    @property
    def path(self) -> str:
        """Return the public pseudo-path, ending in the native artifact name."""
        return self._path

    @property
    def jd_range(self) -> tuple[float, float]:
        """Return the original artifact header range."""
        return self._jd_range  # type: ignore[return-value]

    def has_body(self, body_id: int) -> bool:
        """Return whether the source contains the requested body."""
        return body_id in self._bodies

    def body_coverage(self, body_id: int) -> tuple[float, float] | None:
        """Return the stored body interval, or None when absent."""
        body = self._bodies.get(body_id)
        return (float(body.jd_start), float(body.jd_end)) if body else None

    @abstractmethod
    def fetch_body_segment(self, body_id: int, idx: int) -> Sequence[float]:
        """Fetch one component-major segment without numerical conversion."""

    @abstractmethod
    def fetch_nutation_segment(self, idx: int) -> Sequence[float]:
        """Fetch one nutation segment, or raise CoefficientSourceError."""

    @staticmethod
    def _coefficients(values: Sequence[float], count: int) -> tuple[float, ...]:
        coeffs = tuple(values)
        if len(coeffs) != count or not all(math.isfinite(v) for v in coeffs):
            raise CoefficientSourceError("Invalid coefficient segment payload")
        return coeffs

    def eval_body(
        self, body_id: int, jd: float
    ) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
        """Evaluate position and analytical velocity in native stored units."""
        if self._closed:
            raise ValueError("LEB reader is closed")
        key = (body_id, jd)
        cached = self._eval_cache.get(key)
        if cached is not None:
            return cached
        if body_id not in self._bodies:
            raise KeyError(f"Body {body_id} not in LEB file")
        body = self._bodies[body_id]
        if not math.isfinite(jd) or jd < body.jd_start or jd > body.jd_end:
            raise ValueError(
                f"JD {jd} outside range [{body.jd_start}, {body.jd_end}] "
                f"for body {body_id}"
            )
        seg_idx = int((jd - body.jd_start) / body.interval_days)
        seg_idx = max(0, min(seg_idx, body.segment_count - 1))
        seg_start = body.jd_start + seg_idx * body.interval_days
        seg_mid = seg_start + 0.5 * body.interval_days
        tau = 2.0 * (jd - seg_mid) / body.interval_days
        if tau > 1.0:
            tau = 1.0
        elif tau < -1.0:
            tau = -1.0
        deg1 = body.degree + 1
        coeffs = self._coefficients(
            self.fetch_body_segment(body_id, seg_idx), body.components * deg1
        )
        pos = []
        vel = []
        scale = 2.0 / body.interval_days
        for c in range(body.components):
            val, deriv = _clenshaw_with_derivative(
                coeffs[c * deg1 : (c + 1) * deg1], tau
            )
            pos.append(val)
            vel.append(deriv * scale)
        if body.coord_type in (
            COORD_ECLIPTIC,
            COORD_HELIO_ECL_RETIRED,
            COORD_GEO_ECLIPTIC,
        ):
            pos[0] = pos[0] % 360.0
        result = tuple(pos), tuple(vel)
        if len(self._eval_cache) > 256:
            self._eval_cache.clear()
        self._eval_cache[key] = result  # type: ignore[assignment]
        return result  # type: ignore[return-value]

    def has_nutation(self) -> bool:
        """Return whether the source contains usable nutation segments."""
        return self._nutation is not None and self._nutation.segment_count > 0

    def eval_nutation(self, jd_tt: float) -> tuple[float, float]:
        """Evaluate nutation angles in radians using native LEB arithmetic."""
        if self._closed:
            raise ValueError("LEB reader is closed")
        if not self.has_nutation():
            raise ValueError("No nutation data in this LEB file")
        nut = self._nutation
        assert nut is not None
        if not math.isfinite(jd_tt) or jd_tt < nut.jd_start or jd_tt > nut.jd_end:
            raise ValueError(
                f"JD {jd_tt} outside nutation range [{nut.jd_start}, {nut.jd_end}]"
            )
        seg_idx = int((jd_tt - nut.jd_start) / nut.interval_days)
        seg_idx = max(0, min(seg_idx, nut.segment_count - 1))
        seg_start = nut.jd_start + seg_idx * nut.interval_days
        seg_mid = seg_start + 0.5 * nut.interval_days
        tau = 2.0 * (jd_tt - seg_mid) / nut.interval_days
        if tau > 1.0:
            tau = 1.0
        elif tau < -1.0:
            tau = -1.0
        deg1 = nut.degree + 1
        coeffs = self._coefficients(
            self.fetch_nutation_segment(seg_idx), nut.components * deg1
        )
        return _clenshaw(coeffs[0:deg1], tau), _clenshaw(coeffs[deg1 : 2 * deg1], tau)

    def delta_t(self, jd: float) -> float:
        """Interpolate the source table, returning Delta-T in days."""
        if not self._delta_t_jds:
            raise ValueError("No Delta-T data in this LEB file")
        jds, vals = self._delta_t_jds, self._delta_t_vals
        if jd <= jds[0]:
            return vals[0]
        if jd >= jds[-1]:
            return vals[-1]
        idx = bisect_right(jds, jd) - 1
        idx = max(0, min(idx, len(jds) - 2))
        span = jds[idx + 1] - jds[idx]
        if span == 0.0:
            return vals[idx]
        t = (jd - jds[idx]) / span
        return vals[idx] + t * (vals[idx + 1] - vals[idx])

    def get_star(self, star_id: int) -> StarEntry:
        """Return a native star record, or raise KeyError when absent."""
        if star_id not in self._stars:
            raise KeyError(f"Star {star_id} not in LEB catalog")
        return self._stars[star_id]

    def warm(self, jd_start: float, jd_end: float) -> None:
        """Do not preload externally supplied data implicitly."""

    def cool(self) -> None:
        """Leave shared provider caches under provider ownership."""

    def on_close(self) -> None:
        """Optional instance-resource hook; must not close shared pools."""

    def close(self) -> None:
        """Close this reader while preserving shared provider resources."""
        if not self._closed:
            self._closed = True
            self._eval_cache.clear()
            self.on_close()

    def __enter__(self) -> SegmentSource:
        return self

    def __exit__(self, *args: object) -> None:
        self.close()
