# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Giacomo Battaglia
"""Contract and numerical parity tests for externally supplied LEB segments."""

from __future__ import annotations

import math
import os
import random
from collections.abc import Sequence

import pytest

from libephemeris.exceptions import CoefficientSourceError
from libephemeris.leb_reader import open_leb
from libephemeris.segment_source import SegmentSource


BUNDLED_CORE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
    "libephemeris",
    "data",
    "leb2",
    "base_core.leb2",
)


class FileSegmentSource(SegmentSource):
    """Small test adapter exposing a file reader through the source seam."""

    def __init__(self, reader, *, tier: str = "base") -> None:
        self._reader = reader
        super().__init__(
            artifact_name=f"{tier}_core.leb2",
            locator="file://test-fixture",
            jd_range=reader.header_jd_range,
            bodies=reader.bodies,
            nutation=reader.nutation_header,
            delta_t=reader.delta_t_table,
            stars=reader.stars,
            reviewed=True,
        )

    def fetch_body_segment(self, body_id: int, idx: int) -> Sequence[float]:
        return self._reader.segment_coefficients(body_id, idx)

    def fetch_nutation_segment(self, idx: int) -> Sequence[float]:
        return self._reader.nutation_coefficients(idx)

    def on_close(self) -> None:
        self._reader.close()


def _sample_dates(entry) -> list[float]:
    """Return endpoint and representative segment-boundary dates."""
    points = [entry.jd_start, entry.jd_end]
    if entry.segment_count > 1:
        for index in (1, entry.segment_count // 2, entry.segment_count - 1):
            boundary = entry.jd_start + index * entry.interval_days
            points.extend((math.nextafter(boundary, -math.inf), boundary))
    return points


@pytest.mark.skipif(not os.path.exists(BUNDLED_CORE), reason="bundled core is absent")
def test_segment_source_matches_bundled_reader_bit_for_bit() -> None:
    """The source adapter must preserve native positions and derivatives."""
    reader = open_leb(BUNDLED_CORE)
    source = FileSegmentSource(reader)
    try:
        assert source.path.endswith("base_core.leb2")
        assert source.jd_range == reader.header_jd_range
        assert set(source._bodies) == set(reader.bodies)
        # Exercise the bulk-export seam without decompressing every large
        # body: one segment from each body and one page of each auxiliary.
        for body_id, entry in reader.bodies.items():
            index = min(entry.segment_count - 1, entry.segment_count // 2)
            assert source.fetch_body_segment(
                body_id, index
            ) == reader.segment_coefficients(body_id, index)
            page = next(reader.iter_segment_pages(body_id, page_size=2))
            assert page[0:2] == (0, 0)
            assert (
                len(page[2])
                == min(2, entry.segment_count)
                * (entry.degree + 1)
                * entry.components
                * 8
            )

        rng = random.Random(20261002)
        for body_id, entry in reader.bodies.items():
            dates = _sample_dates(entry)
            dates.extend(rng.uniform(entry.jd_start, entry.jd_end) for _ in range(2))
            for jd in dates:
                assert source.eval_body(body_id, jd) == reader.eval_body(body_id, jd)

        nutation = reader.nutation_header
        assert nutation is not None
        for jd in _sample_dates(nutation):
            assert source.eval_nutation(jd) == reader.eval_nutation(jd)

        jds, values = reader.delta_t_table
        assert len(jds) == len(values)
        for jd in (jds[0], jds[-1], math.nextafter(jds[0], -math.inf)):
            assert source.delta_t(jd) == reader.delta_t(jd)
        for star_id, expected in reader.stars.items():
            assert source.get_star(star_id) == expected
    finally:
        source.close()


def test_export_accessors_match_leb1_fixture(test_leb_file) -> None:
    """The export accessors preserve LEB1 coefficients and auxiliary pages."""
    reader = open_leb(test_leb_file)
    try:
        for body_id, entry in reader.bodies.items():
            assert reader.segment_coefficients(body_id, 0)
            page_no, first, payload = next(reader.iter_segment_pages(body_id, 2))
            assert (page_no, first) == (0, 0)
            assert (
                len(payload)
                == min(2, entry.segment_count)
                * (entry.degree + 1)
                * entry.components
                * 8
            )
    finally:
        reader.close()


def test_source_fetch_payload_errors_are_normalized() -> None:
    """Malformed fetch payloads use the source-specific exception only."""
    from libephemeris.leb_format import BodyEntry, COORD_ICRS_BARY

    class BrokenSource(SegmentSource):
        def fetch_body_segment(self, body_id: int, idx: int) -> Sequence[float]:
            raise CoefficientSourceError("transport detail")

        def fetch_nutation_segment(self, idx: int) -> Sequence[float]:
            raise CoefficientSourceError("transport detail")

    source = BrokenSource(
        artifact_name="medium_core.leb2",
        locator="postgres://dataset",
        jd_range=(1.0, 3.0),
        bodies={
            0: BodyEntry(0, COORD_ICRS_BARY, 1, 1.0, 3.0, 2.0, 0, 3, 0),
        },
    )
    # Providers are responsible for translating transport exceptions at the
    # fetch boundary. The reader must not silently turn them into range misses.
    with pytest.raises(CoefficientSourceError, match="transport detail"):
        source.eval_body(0, 2.0)
    source.close()


def test_source_metadata_and_lifecycle_contract() -> None:
    """Pseudo-paths, reviewed state, no-op warm/cool, and close are stable."""
    from libephemeris.leb_format import BodyEntry, COORD_ICRS_BARY

    class EmptySource(SegmentSource):
        def fetch_body_segment(self, body_id: int, idx: int) -> Sequence[float]:
            return (1.0, 2.0, 3.0)

        def fetch_nutation_segment(self, idx: int) -> Sequence[float]:
            raise CoefficientSourceError("no nutation")

    source = EmptySource(
        artifact_name="extended_exotics.leb2",
        locator="postgres://dataset",
        jd_range=(1.0, 3.0),
        bodies={
            0: BodyEntry(0, COORD_ICRS_BARY, 1, 1.0, 3.0, 2.0, 0, 3, 0),
        },
        reviewed=True,
    )
    assert source.path == "postgres://dataset/extended_exotics.leb2"
    assert source._manifest_verified is True
    source.warm(1.0, 2.0)
    source.cool()
    source.close()
    with pytest.raises(ValueError, match="LEB reader is closed"):
        source.eval_body(0, 2.0)


def test_invalid_artifact_name_is_rejected() -> None:
    """Only canonical tier/group artifact names may enter the seam."""
    from libephemeris.leb_format import BodyEntry, COORD_ICRS_BARY

    class Source(SegmentSource):
        def fetch_body_segment(self, body_id: int, idx: int) -> Sequence[float]:
            return ()

        def fetch_nutation_segment(self, idx: int) -> Sequence[float]:
            return ()

    with pytest.raises(CoefficientSourceError):
        Source(
            artifact_name="secret-host.leb2",
            locator="postgres://dataset",
            jd_range=(1.0, 3.0),
            bodies={
                0: BodyEntry(0, COORD_ICRS_BARY, 1, 1.0, 3.0, 2.0, 0, 3, 0),
            },
        )
