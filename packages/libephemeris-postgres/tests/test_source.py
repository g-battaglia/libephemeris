# SPDX-License-Identifier: AGPL-3.0-only
from __future__ import annotations

import struct

import pytest
from libephemeris import CoefficientSourceError
from libephemeris.leb_format import BodyEntry, COORD_ICRS_BARY

from libephemeris_postgres.config import RuntimeConfig
from libephemeris_postgres.source import PostgresSegmentSource


class FakePages:
    def __init__(self, pages: dict[tuple[int, int], bytes]) -> None:
        self.pages = pages
        self.calls: list[list[tuple[int, int]]] = []

    def fetch_pages(self, keys: list[tuple[int, int]]) -> list[tuple[int, int, bytes]]:
        self.calls.append(keys)
        return [
            (body, page, self.pages[(body, page)])
            for body, page in keys
            if (body, page) in self.pages
        ]


def _source(store: FakePages) -> PostgresSegmentSource:
    body = BodyEntry(4, COORD_ICRS_BARY, 2, 0.0, 2.0, 1.0, 0, 3, 0)
    page = struct.pack("<6d", 1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    store.pages.setdefault((4, 0), page)
    return PostgresSegmentSource(
        dataset_id="dataset",
        artifact_no=0,
        artifact_name="medium_core.leb2",
        locator="postgres://dataset",
        jd_range=(0.0, 2.0),
        bodies={4: body},
        nutation=None,
        delta_t=([0.0], [0.0]),
        stars={},
        reviewed=False,
        config=RuntimeConfig("unused", cache_pages=4),
        page_rows={4: 2},
        store=store,
    )


def test_fetches_and_caches_pages() -> None:
    store = FakePages({})
    source = _source(store)
    assert source.fetch_body_segment(4, 0) == (1.0, 2.0, 3.0)
    assert source.fetch_body_segment(4, 1) == (4.0, 5.0, 6.0)
    assert len(store.calls) == 1


def test_missing_segment_raises_source_error() -> None:
    store = FakePages({})
    source = _source(store)
    store.pages.clear()
    with pytest.raises(CoefficientSourceError):
        source.fetch_body_segment(4, 0)


def test_nonfinite_coefficients_are_rejected() -> None:
    store = FakePages(
        {(4, 0): struct.pack("<6d", 1.0, 2.0, float("nan"), 4.0, 5.0, 6.0)}
    )
    source = _source(store)
    with pytest.raises(CoefficientSourceError):
        source.fetch_body_segment(4, 0)
