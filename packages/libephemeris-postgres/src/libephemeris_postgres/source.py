# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""PostgreSQL-backed segment source with bounded page caching."""

from __future__ import annotations

import math
import struct
import threading
from collections import OrderedDict
from collections.abc import Mapping, Sequence
from typing import Any

from libephemeris.segment_source import SegmentSource

from .config import RuntimeConfig
from .pool import get_pool


class PostgresSegmentSource(SegmentSource):
    """Serve immutable coefficient pages from a PostgreSQL dataset."""

    def __init__(
        self,
        *,
        dataset_id: str,
        artifact_no: int,
        artifact_name: str,
        locator: str,
        jd_range: tuple[float, float],
        bodies: Mapping[int, Any],
        nutation: Any | None,
        delta_t: tuple[Sequence[float], Sequence[float]],
        stars: Mapping[int, Any],
        reviewed: bool,
        config: RuntimeConfig,
        page_rows: Mapping[int, int] | None = None,
        store: Any | None = None,
    ) -> None:
        super().__init__(
            artifact_name=artifact_name,
            locator=locator,
            jd_range=jd_range,
            bodies=bodies,
            nutation=nutation,
            delta_t=delta_t,
            stars=stars,
            reviewed=reviewed,
        )
        self.dataset_id = dataset_id
        self.artifact_no = artifact_no
        self._config = config
        self._store = store
        self._page_rows = {
            int(body_id): int(size) for body_id, size in (page_rows or {}).items()
        }
        if any(
            not math.isfinite(float(size)) or size <= 0
            for size in self._page_rows.values()
        ):
            raise self._error("page size must be positive")
        for body_id, size in self._page_rows.items():
            if body_id not in self._bodies and not (
                body_id == -1 and self._nutation is not None
            ):
                raise self._error(f"page metadata references unknown series {body_id}")
        if any(body_id not in self._page_rows for body_id in self._bodies):
            raise self._error("missing page metadata")
        if self._nutation is not None and -1 not in self._page_rows:
            raise self._error("missing nutation page metadata")
        self._pages: OrderedDict[tuple[int, int], bytes] = OrderedDict()
        self._page_lock = threading.Lock()

    def _series(self, body_id: int) -> Any:
        if body_id == -1:
            if self._nutation is None:
                raise self._error("missing nutation series")
            return self._nutation
        try:
            return self._bodies[body_id]
        except KeyError:
            raise self._error(f"missing body {body_id}") from None

    @staticmethod
    def _error(detail: str) -> Exception:
        from libephemeris import CoefficientSourceError

        return CoefficientSourceError(f"invalid coefficient data ({detail})")

    def _cached(self, key: tuple[int, int]) -> bytes | None:
        with self._page_lock:
            value = self._pages.get(key)
            if value is not None:
                self._pages.move_to_end(key)
            return value

    def _put(self, key: tuple[int, int], value: bytes) -> None:
        with self._page_lock:
            self._pages[key] = value
            self._pages.move_to_end(key)
            while len(self._pages) > self._config.cache_pages:
                self._pages.popitem(last=False)

    def _fetch_pages(
        self, page_keys: list[tuple[int, int]]
    ) -> dict[tuple[int, int], bytes]:
        missing = [key for key in page_keys if self._cached(key) is None]
        if not missing:
            return {}
        try:
            if self._store is not None:
                rows = self._store.fetch_pages(missing)
            else:
                pool = get_pool(self._config)
                body_ids = [key[0] for key in missing]
                page_nos = [key[1] for key in missing]
                with pool.connection() as conn, conn.cursor() as cur:
                    cur.execute(
                        """SELECT body_id, page_no, coeffs
                           FROM libephemeris.pages
                           WHERE dataset_id = %s AND artifact_no = %s
                             AND (body_id, page_no) IN
                               (SELECT * FROM unnest(%s::int[], %s::int[]))""",
                        (self.dataset_id, self.artifact_no, body_ids, page_nos),
                    )
                    rows = cur.fetchall()
        except Exception:
            raise self._unavailable() from None
        try:
            payloads = {(int(b), int(p)): bytes(data) for b, p, data in rows}
        except Exception:
            raise self._error("page payload is not bytes") from None
        if set(payloads) != set(missing):
            raise self._error("missing coefficient pages")
        for key, payload in payloads.items():
            self._put(key, payload)
        return payloads

    def _page_for_jd(self, body_id: int, jd: float) -> int:
        entry = self._series(body_id)
        index = int((jd - float(entry.jd_start)) / float(entry.interval_days))
        index = max(0, min(index, int(entry.segment_count) - 1))
        return index // self._page_rows.get(body_id, 32)

    @staticmethod
    def _unavailable() -> Exception:
        from libephemeris import CoefficientSourceError

        return CoefficientSourceError("PostgreSQL coefficient source unavailable")

    def _page(
        self, body_id: int, page_no: int, representative_jd: float | None = None
    ) -> bytes:
        key = (body_id, page_no)
        value = self._cached(key)
        if value is not None:
            return value
        # A miss fetches the matching epoch for every series in this artifact.
        if representative_jd is None:
            representative_jd = float(self._series(body_id).jd_start)
        keys = [
            (other, self._page_for_jd(other, representative_jd))
            for other in self._bodies
        ]
        if self.has_nutation():
            keys.append((-1, self._page_for_jd(-1, representative_jd)))
        keys.append(key)
        payloads = self._fetch_pages(list(dict.fromkeys(keys)))
        value = payloads.get(key) or self._cached(key)
        if value is None:
            # Another thread may evict a page between lookup and prefetch.
            value = self._fetch_pages([key]).get(key)
        if value is None:
            raise self._error(f"page {body_id}/{page_no} is missing")
        return value

    def fetch_body_segment(self, body_id: int, idx: int) -> tuple[float, ...]:
        entry = self._series(body_id)
        page_size = self._page_rows.get(body_id, 32)
        page_no, in_page = divmod(idx, page_size)
        representative_jd = float(entry.jd_start) + idx * float(entry.interval_days)
        payload = self._page(body_id, page_no, representative_jd)
        width = int(entry.components) * (int(entry.degree) + 1) * 8
        offset = in_page * width
        expected = min(page_size, entry.segment_count - page_no * page_size) * width
        if len(payload) != expected or offset < 0 or offset + width > len(payload):
            raise self._error(f"page {body_id}/{page_no} has invalid length")
        values = struct.unpack_from(
            f"<{int(entry.components) * (int(entry.degree) + 1)}d", payload, offset
        )
        if not all(math.isfinite(value) for value in values):
            raise self._error(f"page {body_id}/{page_no} contains non-finite values")
        return values

    def fetch_nutation_segment(self, idx: int) -> tuple[float, ...]:
        return self.fetch_body_segment(-1, idx)

    def on_close(self) -> None:
        with self._page_lock:
            self._pages.clear()
