# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Portable byte-layout, endpoint and reader parity tests."""

from __future__ import annotations

import math

import pytest

from libephemeris.db.contract import (
    Series,
    segment_byte_count,
    segment_index,
    series_fields,
    validate_series,
)
from libephemeris.db.reader import DBReader
from libephemeris.exceptions import DBDataError, EphemerisRangeError

from .conftest import DATASET_ID


@pytest.mark.parametrize("jd, expected", [(10.0, 0), (12.0, 1), (14.0, 1)])
def test_segment_endpoint_contract(jd, expected):
    """Check inclusive final endpoint and deterministic segment selection.

    Args:
        jd: Synthetic epoch.
        expected: Expected zero-based index.
    """
    series = Series(0, 0, 2, 10.0, 14.0, 2.0, 2, 3)
    assert segment_index(series, jd) == expected
    assert segment_byte_count(series) == 72
    assert series_fields(series) == (0, 0, 2, 10.0, 14.0, 2.0, 2, 3)


@pytest.mark.parametrize("jd", [9.0, 15.0, math.nan, math.inf, -math.inf])
def test_invalid_epoch_is_a_range_error(jd):
    """Reject non-finite epochs before index arithmetic.

    Args:
        jd: Synthetic invalid epoch.
    """
    with pytest.raises(EphemerisRangeError):
        segment_index(Series(0, 0, 2, 10.0, 14.0, 2.0, 2, 3), jd)


@pytest.mark.parametrize("interval", [0.0, -1.0, math.nan, math.inf])
def test_invalid_metadata_is_not_a_fallback_signal(interval):
    """Classify damaged metadata as a DB data failure.

    Args:
        interval: Invalid polynomial interval.
    """
    with pytest.raises(DBDataError):
        validate_series(Series(0, 0, 2, 10.0, 14.0, interval, 2, 3))


def test_reader_matches_native_artifact_at_boundaries(artifact_store):
    """Compare ephemeral file and DB evaluations without golden outputs.

    Args:
        artifact_store: Native artifact used as a test transport source.
    """
    reader = DBReader(artifact_store, DATASET_ID)
    try:
        for body_id, series in reader._bodies.items():
            boundary = series.jd_start + 200 * series.interval_days
            epochs = [
                2451545.0,
                boundary,
                math.nextafter(boundary, -math.inf),
                math.nextafter(boundary, math.inf),
                series.jd_start,
                series.jd_end,
            ]
            for jd in epochs:
                assert reader.eval_body(body_id, jd) == artifact_store.reader.eval_body(
                    body_id, jd
                )
        assert reader.eval_nutation(2451545.0) == artifact_store.reader.eval_nutation(
            2451545.0
        )
        for jd in [1000000.0, 2451545.0, 4000000.0]:
            assert reader.delta_t(jd) == artifact_store.reader.delta_t(jd)
    finally:
        reader.close()
    assert reader._series == {}
    assert reader._segments == {}


def test_malformed_payload_fails_closed(artifact_store, monkeypatch):
    """Reject a short payload instead of using another backend.

    Args:
        artifact_store: Native transport double.
        monkeypatch: Pytest patch helper.
    """
    reader = DBReader(artifact_store, DATASET_ID)

    def short_payload(dataset_id, keys):
        """Return intentionally truncated inputs.

        Args:
            dataset_id: Requested version.
            keys: Requested segments.

        Returns:
            Invalid test bytes.
        """
        return {key: b"bad" for key in keys}

    monkeypatch.setattr(artifact_store, "fetch_segments", short_payload)
    with pytest.raises(DBDataError, match="length"):
        reader.eval_body(0, 2451545.0)
    with pytest.raises(KeyError):
        reader.eval_body(999999, 2451545.0)
    reader.close()
    with pytest.raises(DBDataError, match="closed"):
        reader.prepare(2451545.0, 0, 0)


def test_event_scan_keeps_only_bounded_operation_inputs(artifact_store):
    """Ensure a long operation cannot accumulate the whole dataset.

    Args:
        artifact_store: Instrumented transport double.
    """
    reader = DBReader(artifact_store, DATASET_ID)
    series = reader._bodies[1]
    for index in range(reader.MAX_INPUT_SEGMENTS + 20):
        jd = series.jd_start + (index + 0.5) * series.interval_days
        reader.eval_body(1, jd)
        assert len(reader._segments) <= reader.MAX_INPUT_SEGMENTS
    reader.close()
