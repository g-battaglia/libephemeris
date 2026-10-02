# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Test-only storage doubles using native project artifacts, never snapshots."""

from __future__ import annotations

from bisect import bisect_left
from pathlib import Path

import pytest

import libephemeris as ephemeris
from libephemeris.db import backend
from libephemeris.db.contract import NUTATION_ID, Series
from libephemeris.db.importer import Artifact
from libephemeris.leb2_reader import LEB2Reader

DATASET_ID = "01234567-89ab-cdef-0123-456789abcdef"
BUNDLED_FILE = Path(ephemeris.__file__).parent / "data/leb2/base_core.leb2"


class ArtifactStore:
    """Test source of truth; runtime DB reader never sees the file handle."""

    def __init__(self, reader: LEB2Reader) -> None:
        """Configure a source reader and observable transport counters.

        Args:
            reader: Native artifact opened before selecting DB mode.
        """
        self.reader = reader
        self.metadata_calls = 0
        self.segment_calls = 0
        self.transferred_keys = []

    def metadata(self, dataset_id: str) -> dict[int, Series]:
        """Return newly allocated source headers for each DB operation.

        Args:
            dataset_id: Explicit version, ignored by this single-version double.

        Returns:
            Detached metadata records in the portable format.
        """
        self.metadata_calls += 1
        result = {}
        for body in self.reader._bodies.values():
            result[body.body_id] = Series(
                body.body_id,
                body.coord_type,
                body.segment_count,
                body.jd_start,
                body.jd_end,
                body.interval_days,
                body.degree,
                body.components,
            )
        nutation = self.reader._nutation
        result[NUTATION_ID] = Series(
            NUTATION_ID,
            0,
            nutation.segment_count,
            nutation.jd_start,
            nutation.jd_end,
            nutation.interval_days,
            nutation.degree,
            nutation.components,
        )
        return result

    def fetch_segments(self, dataset_id: str, keys: list) -> dict:
        """Emulate exact indexed DB reads while counting transport calls.

        Args:
            dataset_id: Requested version.
            keys: Body/segment pairs.

        Returns:
            Detached native coefficient payloads, not precomputed positions.
        """
        self.segment_calls += 1
        self.transferred_keys.extend(keys)
        result = {}
        for body_id, index in keys:
            if body_id == NUTATION_ID:
                nutation = self.reader._nutation
                width = (nutation.degree + 1) * nutation.components * 8
                offset = self.reader._nutation_data_offset + index * width
                result[(body_id, index)] = bytes(
                    self.reader._mm[offset : offset + width]
                )
            else:
                body = self.reader._bodies[body_id]
                for chunk_index, chunk in enumerate(self.reader._chunk_index[body_id]):
                    if (
                        chunk.segment_start
                        <= index
                        < chunk.segment_start + chunk.segment_count
                    ):
                        decoded = self.reader._decompress_chunk(body_id, chunk_index)
                        width = (body.degree + 1) * body.components * 8
                        offset = (index - chunk.segment_start) * width
                        result[(body_id, index)] = decoded[offset : offset + width]
                        break
        return result

    def delta_t_points(self, dataset_id: str, jd: float) -> list:
        """Return the same two bracketing samples as the DB query.

        Args:
            dataset_id: Requested version.
            jd: Requested epoch.

        Returns:
            Ordered samples with endpoint clamping.
        """
        dates = self.reader._delta_t_jds
        values = self.reader._delta_t_vals
        index = bisect_left(dates, jd)
        indices = {max(0, index - 1), min(index, len(dates) - 1)}
        return [(dates[i], values[i]) for i in sorted(indices)]

    def star(self, dataset_id: str, star_id: int) -> tuple:
        """Read a detached source star record.

        Args:
            dataset_id: Requested version.
            star_id: Source identifier.

        Returns:
            Portable star fields.
        """
        star = self.reader.get_star(star_id)
        return (
            star.star_id,
            star.ra_j2000,
            star.dec_j2000,
            star.pm_ra,
            star.pm_dec,
            star.parallax,
            star.rv,
            star.magnitude,
        )

    def close(self) -> None:
        """Leave the test-owned file reader under fixture lifetime control."""


@pytest.fixture
def artifact_store():
    """Provide a disposable native artifact transport double.

    Yields:
        A source store and its independently owned file reader.
    """
    with LEB2Reader(str(BUNDLED_FILE)) as reader:
        yield ArtifactStore(reader)


@pytest.fixture
def db_runtime(monkeypatch, artifact_store):
    """Select DB mode without making real connections or persisting outputs.

    Args:
        monkeypatch: Pytest restoration helper.
        artifact_store: Disposable transport double.

    Yields:
        Instrumented store used by public calculation entry points.
    """
    previous_mode = ephemeris.get_calc_mode()
    previous_policy = ephemeris.get_configured_network_policy()
    ephemeris.close()
    backend.set_db_config("postgresql://test.invalid/ephemeris", DATASET_ID)
    monkeypatch.setattr(backend, "PostgresStore", lambda dsn: artifact_store)
    ephemeris.set_calc_mode("db")
    ephemeris.set_network_policy("auto")
    try:
        yield artifact_store
    finally:
        ephemeris.close()
        backend.set_db_config(None)
        ephemeris.set_calc_mode(previous_mode)
        ephemeris.set_network_policy(previous_policy)
