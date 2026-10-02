# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Optional PostgreSQL transport, independent of the astronomical engine.

Provenance:
    Project-authored SQL transport over the project-native schema documented in
    docs/db-backend.md. Connections and queries define no astronomical model
    or coefficient and never alter mathematical values.
"""

from __future__ import annotations

import os
import threading
from contextlib import contextmanager
from typing import Any, Iterator

from .contract import SCHEMA_VERSION, SegmentKey, Series, validate_series
from ..exceptions import DBDataError, DBError, StarNotFoundError


class PostgresStore:
    """Own connections, never coefficient records or computed states."""

    def __init__(self, dsn: str) -> None:
        """Configure a lazy, process-owned connection pool.

        Args:
            dsn: PostgreSQL connection string. Never included in diagnostics.
        """
        self._dsn = dsn
        self._pid = os.getpid()
        self._pool: Any = None
        self._driver_errors: tuple[type[Exception], ...] = ()
        self._lock = threading.Lock()

    def _get_pool(self) -> Any:
        """Create the optional driver pool on first use.

        Returns:
            A bounded Psycopg connection pool.

        Raises:
            DBError: The extra is missing or pool creation fails.
        """
        with self._lock:
            if self._pool is not None:
                return self._pool
            try:
                import psycopg
                from psycopg_pool import ConnectionPool, PoolTimeout
            except ImportError:
                raise DBError("Install libephemeris[postgres] for DB mode") from None

            self._driver_errors = (psycopg.Error, PoolTimeout)
            try:
                # Parse synchronously so a malformed secret-bearing DSN never
                # reaches the pool's background connection-error logger.
                psycopg.conninfo.conninfo_to_dict(self._dsn)
                self._pool = ConnectionPool(
                    self._dsn,
                    min_size=0,
                    max_size=4,
                    open=True,
                    timeout=10,
                    kwargs={
                        "autocommit": True,
                        "connect_timeout": 10,
                        "options": (
                            "-c default_transaction_read_only=on "
                            "-c statement_timeout=10000"
                        ),
                    },
                )
            except Exception:
                # A malformed DSN can be quoted by the driver: redact it.
                raise DBError("Could not configure PostgreSQL connections") from None
            return self._pool

    @contextmanager
    def connection(self) -> Iterator[Any]:
        """Lease a read-only connection and enforce the DB network boundary.

        Yields:
            One connection. No lease is retained after leaving the context.

        Raises:
            DBError: Driver operation fails or a pool is reused after fork.
        """
        from .policy import require_database

        require_database()
        if os.getpid() != self._pid:
            raise DBError("Create the DB backend after forking workers")
        pool = self._get_pool()
        try:
            with pool.connection() as connection:
                yield connection
        except self._driver_errors:
            # Only transport errors are translated. Application exceptions
            # raised inside the context must retain their original identity.
            raise DBError("PostgreSQL connection, query or timeout failure") from None

    def _rows(self, query: str, parameters: tuple) -> list:
        """Execute one parameterized read with binary result transfer.

        Args:
            query: SQL statement owned by this transport module.
            parameters: Bound SQL values, never interpolated into the text.

        Returns:
            Result records detached from the connection lease.
        """
        with self.connection() as connection:
            return connection.execute(query, parameters, binary=True).fetchall()

    def metadata(self, dataset_id: str) -> dict[int, Series]:
        """Read and validate all small series headers in one round-trip.

        Args:
            dataset_id: Immutable dataset UUID.

        Returns:
            Validated metadata, owned by the calling operation.

        Raises:
            DBDataError: Dataset is unavailable or schema/headers are invalid.
        """
        records = self._rows(
            "SELECT v.version, d.published, s.body_id, s.coord_type, "
            "s.segment_count, s.jd_start, s.jd_end, s.interval_days, "
            "s.degree, s.components "
            "FROM libephemeris.datasets d "
            "CROSS JOIN libephemeris.schema_version v "
            "LEFT JOIN libephemeris.series s USING (dataset_id) "
            "WHERE d.dataset_id = %s ORDER BY s.body_id",
            (dataset_id,),
        )
        if not records or not records[0][1]:
            raise DBDataError("DB dataset is absent or unpublished")
        if any(record[0] != SCHEMA_VERSION for record in records):
            raise DBDataError("Unsupported PostgreSQL ephemeris schema")
        if any(record[2] is None for record in records):
            raise DBDataError("DB dataset contains no coefficient series")

        result = {}
        for record in records:
            series = Series(*record[2:])
            validate_series(series)
            result[series.body_id] = series
        return result

    def dataset_info(self, dataset_id: str) -> tuple[str, dict]:
        """Read publication identity without retaining a manifest globally.

        Args:
            dataset_id: Immutable UUID.

        Returns:
            Declared tier and source-artifact manifest.

        Raises:
            DBDataError: Dataset is absent or unpublished.
        """
        records = self._rows(
            "SELECT tier, manifest FROM libephemeris.datasets "
            "WHERE dataset_id = %s AND published",
            (dataset_id,),
        )
        if len(records) != 1:
            raise DBDataError("DB dataset is absent or unpublished")
        return records[0][0], records[0][1]

    def fetch_segments(
        self, dataset_id: str, keys: list[SegmentKey]
    ) -> dict[SegmentKey, bytes]:
        """Fetch only requested segments using a paired-key SQL join.

        Args:
            dataset_id: Immutable dataset UUID.
            keys: Body/segment pairs. Duplicates are removed before transfer.

        Returns:
            Exact little-endian coefficient bytes indexed by key.

        Raises:
            DBDataError: The published dataset lacks a requested segment.
        """
        if not keys:
            return {}
        unique_keys = list(dict.fromkeys(keys))
        body_ids = [key[0] for key in unique_keys]
        segment_indices = [key[1] for key in unique_keys]
        records = self._rows(
            "SELECT s.body_id, s.segment_index, s.payload "
            "FROM libephemeris.segments s JOIN "
            "unnest(%s::integer[], %s::integer[]) AS k(body_id, segment_index) "
            "USING (body_id, segment_index) WHERE s.dataset_id = %s",
            (body_ids, segment_indices, dataset_id),
        )
        result = {
            (body_id, index): bytes(payload) for body_id, index, payload in records
        }
        if set(result) != set(unique_keys):
            raise DBDataError("Missing coefficient segments in published DB dataset")
        return result

    def delta_t_points(self, dataset_id: str, jd: float) -> list[tuple[float, float]]:
        """Read at most two indexed samples, including endpoint clamping.

        Args:
            dataset_id: Immutable dataset UUID.
            jd: Epoch to bracket.

        Returns:
            Ordered (Julian day, Delta-T days) samples.
        """
        return self._rows(
            "(SELECT jd, days FROM libephemeris.delta_t WHERE dataset_id = %s "
            "AND jd <= %s ORDER BY jd DESC LIMIT 1) UNION "
            "(SELECT jd, days FROM libephemeris.delta_t WHERE dataset_id = %s "
            "AND jd >= %s ORDER BY jd LIMIT 1) ORDER BY jd",
            (dataset_id, jd, dataset_id, jd),
        )

    def star(self, dataset_id: str, star_id: int) -> tuple:
        """Look up one star without loading the catalog.

        Args:
            dataset_id: Immutable dataset UUID.
            star_id: Source catalog identifier.

        Returns:
            Fields in LEB StarEntry order.

        Raises:
            StarNotFoundError: Star is absent from the dataset.
        """
        records = self._rows(
            "SELECT star_id, ra, dec, pm_ra, pm_dec, parallax, rv, magnitude "
            "FROM libephemeris.stars WHERE dataset_id = %s AND star_id = %s",
            (dataset_id, star_id),
        )
        if not records:
            raise StarNotFoundError("Star not stored in DB dataset")
        return tuple(records[0])

    def close(self) -> None:
        """Release pool resources without forgetting connection configuration.

        Subsequent reads create a new pool. Call after in-flight operations
        finish; a closing pool may reject operations concurrently acquiring it.
        """
        with self._lock:
            if self._pool is not None:
                self._pool.close()
                self._pool = None
