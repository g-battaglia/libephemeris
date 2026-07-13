"""Regression tests for the independent JPL planet-center artifact builder."""

from __future__ import annotations

import sys
from types import SimpleNamespace

import numpy as np
import pytest

from scripts import generate_planet_centers_spk as generator


_SECONDS_PER_NOMINAL_YEAR = 365.25 * 86400.0


def _date_to_test_et(date: str) -> float:
    """Map the leading civil year to deterministic test ET seconds."""
    return float(date.split()[0]) * _SECONDS_PER_NOMINAL_YEAR


def _required_medium_windows() -> dict[int, list[tuple[float, float]]]:
    return {
        target_id: [(_date_to_test_et(start), _date_to_test_et(end))]
        for target_id, (_, start, end) in generator.TIER_REQUIRED_COVERAGE_DATES[
            "medium"
        ].items()
    }


def test_reviewed_medium_sources_are_hash_pinned() -> None:
    """Every official input to the reviewed medium build has an exact pin."""
    filenames = {
        url.rsplit("/", 1)[-1]
        for source in generator.TIER_SOURCES["medium"].values()
        for url, _, _ in (source if isinstance(source, list) else [source])
    }
    filenames.add("naif0012.tls")

    assert filenames <= generator.REVIEWED_SOURCE_SHA256.keys()


def test_extract_segment_copies_every_matching_daf_descriptor(monkeypatch) -> None:
    """A target split across DAF segments must retain its complete coverage."""
    descriptors = [
        (np.array([0.0, 10.0]), np.array([599, 5, 1, 2, 0, 0])),
        (np.array([10.0, 20.0]), np.array([599, 5, 1, 2, 0, 0])),
        (np.array([0.0, 20.0]), np.array([501, 5, 1, 2, 0, 0])),
        (np.array([0.0, 20.0]), np.array([599, 0, 1, 2, 0, 0])),
    ]
    state = {"index": -1}
    copied = []

    def _next_descriptor() -> bool:
        state["index"] += 1
        return state["index"] < len(descriptors)

    fake_spice = SimpleNamespace(
        spklef=lambda path: 17,
        spkcov=lambda path, target: object(),
        wncard=lambda cover: 1,
        wnfetd=lambda cover, index: (0.0, 20.0),
        et2utc=lambda et, fmt, precision: str(et),
        dafbfs=lambda handle: state.update(index=-1),
        daffna=_next_descriptor,
        dafgs=lambda: descriptors[state["index"]],
        dafus=lambda raw, nd, ni: raw,
        dafgn=lambda: f"segment-{state['index']}",
        spksub=lambda *args: copied.append(args),
        spkuef=lambda handle: None,
    )
    monkeypatch.setitem(sys.modules, "spiceypy", fake_spice)

    coverage = generator.extract_segment("source.bsp", 23, 599, "jupiter", center_id=5)

    assert coverage == (0.0, 20.0)
    assert len(copied) == 2
    assert [args[3:5] for args in copied] == [(0.0, 10.0), (10.0, 20.0)]


def test_extract_segment_clips_and_tolerates_utc_format_limit(monkeypatch) -> None:
    """Tier clipping must work even when CSPICE cannot format source dates."""
    descriptors = [
        (np.array([-100.0, 10.0]), np.array([699, 6, 1, 2, 0, 0])),
        (np.array([10.0, 100.0]), np.array([699, 6, 1, 2, 0, 0])),
    ]
    state = {"index": -1}
    copied = []

    def _next_descriptor() -> bool:
        state["index"] += 1
        return state["index"] < len(descriptors)

    def _limited_et2utc(et, fmt, precision):
        del fmt, precision
        if et < 0.0 or et > 20.0:
            raise ValueError("YEAROUTOFRANGE")
        return str(et)

    fake_spice = SimpleNamespace(
        spklef=lambda path: 17,
        spkcov=lambda path, target: object(),
        wncard=lambda cover: 1,
        wnfetd=lambda cover, index: (-100.0, 100.0),
        et2utc=_limited_et2utc,
        dafbfs=lambda handle: state.update(index=-1),
        daffna=_next_descriptor,
        dafgs=lambda: descriptors[state["index"]],
        dafus=lambda raw, nd, ni: raw,
        dafgn=lambda: f"segment-{state['index']}",
        spksub=lambda *args: copied.append(args),
        spkuef=lambda handle: None,
    )
    monkeypatch.setitem(sys.modules, "spiceypy", fake_spice)

    coverage = generator.extract_segment(
        "source.bsp",
        23,
        699,
        "saturn",
        clip_start_et=5.0,
        clip_end_et=15.0,
        center_id=6,
    )

    assert coverage == (5.0, 15.0)
    assert [args[3:5] for args in copied] == [(5.0, 10.0), (10.0, 15.0)]


def test_medium_coverage_accepts_all_required_endpoints_with_et_tolerance() -> None:
    """Sub-second endpoint differences must not reject a complete artifact."""
    tolerance = generator._COVERAGE_ENDPOINT_TOLERANCE_ET
    windows = _required_medium_windows()
    for target_id, [(start_et, end_et)] in windows.items():
        windows[target_id] = [(start_et + tolerance / 2.0, end_et - tolerance / 2.0)]

    generator._validate_required_coverage("medium", windows, _date_to_test_et)


@pytest.mark.parametrize(
    ("target_id", "endpoint"),
    [
        (599, "end"),
        (699, "start"),
    ],
)
def test_medium_coverage_rejects_short_endpoint(target_id: int, endpoint: str) -> None:
    """A truncation outside the numeric tolerance must fail closed."""
    windows = _required_medium_windows()
    [(start_et, end_et)] = windows[target_id]
    shortfall = generator._COVERAGE_ENDPOINT_TOLERANCE_ET + 1.0
    if endpoint == "start":
        start_et += shortfall
    else:
        end_et -= shortfall
    windows[target_id] = [(start_et, end_et)]

    with pytest.raises(RuntimeError, match=rf"\({target_id}\) requires ET"):
        generator._validate_required_coverage("medium", windows, _date_to_test_et)


def test_medium_coverage_rejects_internal_gap() -> None:
    """Min/max endpoints cannot conceal a gap between copied descriptors."""
    windows = _required_medium_windows()
    [(start_et, end_et)] = windows[799]
    midpoint = (start_et + end_et) / 2.0
    gap = generator._COVERAGE_ENDPOINT_TOLERANCE_ET + 1.0
    windows[799] = [
        (start_et, midpoint),
        (midpoint + gap, end_et),
    ]

    with pytest.raises(RuntimeError, match=r"Uranus \(799\)"):
        generator._validate_required_coverage("medium", windows, _date_to_test_et)
