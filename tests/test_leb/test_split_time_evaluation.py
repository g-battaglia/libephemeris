# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Synthetic tests for two-part LEB epoch evaluation."""

from __future__ import annotations

import math
import struct
from types import SimpleNamespace
from typing import Any

import pytest

from libephemeris.constants import (
    EARTH,
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_NOABERR,
    FLG_NOGDEFL,
    SUN,
)
from libephemeris.fast_calc import C_LIGHT_AU_DAY, _pipeline_icrs
from libephemeris.leb_composite import CompositeLEBReader, TieredLEBReader
from libephemeris.leb_format import BodyEntry, ChunkEntry, CompressedBodyEntry
from libephemeris.leb_reader import LEBReader
from libephemeris.leb2_reader import LEB2Reader

JD0 = 2451545.0


def _segment_bytes(segments: list[tuple[tuple[float, ...], ...]]) -> bytes:
    values = [
        coeff for segment in segments for component in segment for coeff in component
    ]
    return struct.pack(f"<{len(values)}d", *values)


def _leb1_reader(
    segments: list[tuple[tuple[float, ...], ...]],
    interval: float = 1.0,
    *,
    jd_start: float = JD0,
) -> LEBReader:
    degree = len(segments[0][0]) - 1
    reader = object.__new__(LEBReader)
    reader._mm = _segment_bytes(segments)
    reader._eval_cache = {}
    reader._bodies = {
        SUN: BodyEntry(
            body_id=SUN,
            coord_type=0,
            segment_count=len(segments),
            jd_start=jd_start,
            jd_end=jd_start + len(segments) * interval,
            interval_days=interval,
            degree=degree,
            components=3,
            data_offset=0,
        )
    }
    return reader


def _leb2_reader(
    segments: list[tuple[tuple[float, ...], ...]],
    interval: float = 1.0,
    *,
    chunk_size: int | None = None,
    jd_start: float = JD0,
) -> LEB2Reader:
    degree = len(segments[0][0]) - 1
    raw = _segment_bytes(segments)
    reader = object.__new__(LEB2Reader)
    reader._mm = object()
    reader._eval_cache = {}
    reader._bodies = {
        SUN: CompressedBodyEntry(
            body_id=SUN,
            coord_type=0,
            segment_count=len(segments),
            jd_start=jd_start,
            jd_end=jd_start + len(segments) * interval,
            interval_days=interval,
            degree=degree,
            components=3,
            data_offset=0,
            compressed_size=len(raw),
            uncompressed_size=len(raw),
        )
    }
    if chunk_size is None:
        reader._chunked = False
        reader._decompress_body = lambda body_id: raw
        reader._chunk_index = {}
    else:
        reader._chunked = True
        chunks = []
        chunk_data: dict[int, bytes] = {}
        for chunk_idx, first in enumerate(range(0, len(segments), chunk_size)):
            count = min(chunk_size, len(segments) - first)
            data = _segment_bytes(segments[first : first + count])
            chunk_data[chunk_idx] = data
            chunks.append(
                ChunkEntry(
                    jd_start=jd_start + first * interval,
                    jd_end=jd_start + (first + count) * interval,
                    blob_offset=0,
                    compressed_size=len(data),
                    uncompressed_size=len(data),
                    segment_start=first,
                    segment_count=count,
                )
            )
        reader._chunk_index = {SUN: chunks}
        reader._decompress_chunk = lambda body_id, chunk_idx: chunk_data[chunk_idx]
    return reader


def _constant_segment(value: float, degree: int = 1) -> tuple[tuple[float, ...], ...]:
    coeffs = (value,) + (0.0,) * degree
    return (coeffs, coeffs, coeffs)


@pytest.mark.unit
@pytest.mark.parametrize("reader_factory", [_leb1_reader, _leb2_reader])
def test_split_degree_one_preserves_sub_ulp_position_and_velocity(reader_factory):
    offset = math.ulp(JD0) / 4.0
    segment = (((0.0, 1.0)), ((0.0, 2.0)), ((0.0, -3.0)))
    reader = reader_factory([segment], interval=2.0)

    pos, vel = reader._eval_body_split(SUN, JD0 + 1.0, offset)
    tau = offset

    assert pos == pytest.approx((tau, 2.0 * tau, -3.0 * tau), abs=1e-24)
    assert vel == (1.0, 2.0, -3.0)
    assert pos != reader.eval_body(SUN, JD0 + 1.0)[0]


@pytest.mark.unit
@pytest.mark.parametrize("jd_start", [JD0, -4_000_000.0])
@pytest.mark.parametrize("factory_kwargs", [{}, {"chunk_size": 1}])
@pytest.mark.parametrize(
    "offset",
    [
        math.nextafter(0.0, math.inf),
        math.nextafter(0.0, -math.inf),
        math.ulp(JD0) / 4.0,
        -math.ulp(JD0) / 4.0,
    ],
)
def test_nonbinary_midpoint_extends_legacy_tau_continuously(
    jd_start, factory_kwargs, offset
):
    interval = math.nextafter(64.0 / 7.0, math.inf)
    segments = [_constant_segment(0.0), (((0.0, 1.0)),) * 3]
    reader_factory = _leb2_reader if factory_kwargs else _leb1_reader
    reader = reader_factory(
        segments,
        interval,
        jd_start=jd_start,
        **factory_kwargs,
    )
    seg_start = jd_start + interval
    seg_mid = seg_start + 0.5 * interval
    legacy_pos, legacy_vel = reader.eval_body(SUN, seg_mid)
    split_pos, split_vel = reader._eval_body_split(SUN, seg_mid, offset)
    expected_tau = 2.0 * math.fsum((seg_mid, -seg_mid, offset)) / interval

    assert legacy_pos == (0.0, 0.0, 0.0)
    assert split_pos == pytest.approx((expected_tau,) * 3, abs=5e-324)
    assert split_vel == legacy_vel


@pytest.mark.unit
@pytest.mark.parametrize("jd_start", [JD0, -4_000_000.0])
@pytest.mark.parametrize("factory_kwargs", [{}, {"chunk_size": 1}])
def test_nonbinary_midpoint_degree_two_position_and_velocity_share_legacy_tau(
    jd_start, factory_kwargs
):
    interval = math.nextafter(64.0 / 7.0, math.inf)
    coeffs = (1.0, 2.0, 3.0)
    reader_factory = _leb2_reader if factory_kwargs else _leb1_reader
    reader = reader_factory(
        [_constant_segment(0.0, 2), (coeffs, coeffs, coeffs)],
        interval,
        jd_start=jd_start,
        **factory_kwargs,
    )
    seg_start = jd_start + interval
    seg_mid = seg_start + 0.5 * interval
    offset = math.ulp(seg_mid) / 4.0
    tau = 2.0 * math.fsum((seg_mid, -seg_mid, offset)) / interval
    pos, vel = reader._eval_body_split(SUN, seg_mid, offset)

    expected_pos = 1.0 + 2.0 * tau + 3.0 * (2.0 * tau * tau - 1.0)
    expected_vel = (2.0 + 12.0 * tau) * 2.0 / interval
    assert pos == pytest.approx((expected_pos,) * 3, abs=1e-15)
    assert vel == pytest.approx((expected_vel,) * 3, abs=1e-15)


@pytest.mark.unit
@pytest.mark.parametrize("reader_factory", [_leb1_reader, _leb2_reader])
def test_split_degree_two_uses_same_tau_for_position_and_derivative(reader_factory):
    offset = -math.ulp(JD0) / 8.0
    coeffs = (1.0, 2.0, 3.0)
    reader = reader_factory([(coeffs, coeffs, coeffs)], interval=2.0)

    pos, vel = reader._eval_body_split(SUN, JD0 + 1.0, offset)
    tau = offset
    expected_pos = 1.0 + 2.0 * tau + 3.0 * (2.0 * tau * tau - 1.0)
    expected_vel = 2.0 + 12.0 * tau

    assert pos == pytest.approx((expected_pos,) * 3, abs=1e-15)
    assert vel == pytest.approx((expected_vel,) * 3, abs=1e-15)


@pytest.mark.unit
@pytest.mark.parametrize("reader_factory", [_leb1_reader, _leb2_reader])
def test_split_range_and_internal_segment_boundaries(reader_factory):
    offset = math.ulp(JD0) / 4.0
    reader = reader_factory(
        [_constant_segment(10.0), _constant_segment(20.0), _constant_segment(30.0)]
    )

    assert reader._eval_body_split(SUN, JD0 + 1.0, -offset)[0][0] == 10.0
    assert reader._eval_body_split(SUN, JD0 + 1.0, 0.0)[0][0] == 20.0
    assert reader._eval_body_split(SUN, JD0 + 1.0, offset)[0][0] == 20.0
    assert reader._eval_body_split(SUN, JD0, offset)[0][0] == 10.0
    assert reader._eval_body_split(SUN, JD0 + 3.0, -offset)[0][0] == 30.0
    with pytest.raises(ValueError, match="outside range"):
        reader._eval_body_split(SUN, JD0, -offset)
    with pytest.raises(ValueError, match="outside range"):
        reader._eval_body_split(SUN, JD0 + 3.0, offset)


@pytest.mark.unit
@pytest.mark.parametrize("jd_start", [JD0, -4_000_000.0])
@pytest.mark.parametrize("factory_kwargs", [{}, {"chunk_size": 1}])
@pytest.mark.parametrize("boundary_index", [1, 2])
def test_nonbinary_interval_preserves_both_boundary_faces(
    jd_start, factory_kwargs, boundary_index
):
    interval = math.nextafter(64.0 / 7.0, math.inf)
    rounded_boundary = jd_start + boundary_index * interval
    rounding_residual = math.fsum(
        (jd_start, boundary_index * interval, -rounded_boundary)
    )
    reader_factory = _leb2_reader if factory_kwargs else _leb1_reader
    reader = reader_factory(
        [_constant_segment(10.0), _constant_segment(20.0), _constant_segment(30.0)],
        interval,
        jd_start=jd_start,
        **factory_kwargs,
    )

    if rounding_residual > 0.0:
        before = rounding_residual / 2.0
        after = math.nextafter(rounding_residual, math.inf)
    else:
        before = math.nextafter(rounding_residual, -math.inf)
        after = rounding_residual / 2.0

    def exact_delta(offset):
        return math.fsum(
            (
                rounded_boundary,
                -jd_start,
                offset,
                -boundary_index * interval,
            )
        )

    expected_before = (10.0, 20.0)[boundary_index - 1]
    expected_after = (20.0, 30.0)[boundary_index - 1]
    assert exact_delta(before) < 0.0
    assert exact_delta(after) > 0.0
    assert (
        reader._eval_body_split(SUN, rounded_boundary, before)[0][0] == expected_before
    )
    assert reader._eval_body_split(SUN, rounded_boundary, after)[0][0] == expected_after
    # Reversed repeat calls prove the two offsets do not alias in the cache.
    assert reader._eval_body_split(SUN, rounded_boundary, after)[0][0] == expected_after
    assert (
        reader._eval_body_split(SUN, rounded_boundary, before)[0][0] == expected_before
    )


@pytest.mark.unit
def test_split_chunk_boundary_uses_right_chunk_and_true_side():
    offset = math.ulp(JD0) / 4.0
    reader = _leb2_reader(
        [
            _constant_segment(10.0),
            _constant_segment(11.0),
            _constant_segment(20.0),
            _constant_segment(21.0),
        ],
        chunk_size=2,
    )

    assert reader._eval_body_split(SUN, JD0 + 2.0, -offset)[0][0] == 11.0
    assert reader._eval_body_split(SUN, JD0 + 2.0, 0.0)[0][0] == 20.0
    assert reader._eval_body_split(SUN, JD0 + 2.0, offset)[0][0] == 20.0


@pytest.mark.unit
@pytest.mark.parametrize("reader_factory", [_leb1_reader, _leb2_reader])
def test_split_cache_key_includes_offset_and_zero_is_identical(reader_factory):
    offset = math.ulp(JD0) / 4.0
    segment = (((0.0, 1.0)), ((0.0, 1.0)), ((0.0, 1.0)))
    reader = reader_factory([segment], interval=2.0)
    jd = JD0 + 1.0

    public = reader.eval_body(SUN, jd)
    zero_split = reader._eval_body_split(SUN, jd, 0.0)
    positive = reader._eval_body_split(SUN, jd, offset)
    negative = reader._eval_body_split(SUN, jd, -offset)

    assert zero_split == public
    assert struct.pack("<6d", *(zero_split[0] + zero_split[1])) == struct.pack(
        "<6d", *(public[0] + public[1])
    )
    assert positive[0][0] > public[0][0] > negative[0][0]
    assert reader._eval_body_split(SUN, jd, offset) == positive
    assert reader._eval_body_split(SUN, jd, -offset) == negative


@pytest.mark.unit
@pytest.mark.parametrize("reader_factory", [_leb1_reader, _leb2_reader])
@pytest.mark.parametrize("offset", [math.inf, -math.inf, math.nan])
def test_split_rejects_nonfinite_offset(reader_factory, offset):
    reader = reader_factory([_constant_segment(1.0)])
    with pytest.raises(ValueError, match="finite"):
        reader._eval_body_split(SUN, JD0, offset)


class _DispatchReader:
    def __init__(self, start: float, end: float, value: float) -> None:
        self._bodies = {SUN: SimpleNamespace(jd_start=start, jd_end=end)}
        self._start = start
        self._end = end
        self._value = value
        self.calls: list[tuple[int, float, float]] = []

    def body_coverage(self, body_id: int) -> tuple[float, float]:
        return self._start, self._end

    def eval_body(self, body_id: int, jd: float):
        return ((self._value,) * 3, (0.0,) * 3)

    def _eval_body_split(self, body_id: int, jd: float, offset: float):
        self.calls.append((body_id, jd, offset))
        return ((self._value,) * 3, (0.0,) * 3)


@pytest.mark.unit
def test_tiered_nonbinary_coverage_boundary_keeps_offset_and_cache_distinct():
    interval = math.nextafter(64.0 / 7.0, math.inf)
    boundary = JD0 + interval
    offset = math.ulp(boundary) / 4.0
    base_child = _DispatchReader(JD0, boundary, 1.0)
    medium_child = _DispatchReader(boundary, boundary + interval, 2.0)
    tiered = TieredLEBReader(
        {
            "base": CompositeLEBReader([base_child]),
            "medium": CompositeLEBReader([medium_child]),
        }
    )

    assert tiered._eval_body_split(SUN, boundary, -offset)[0][0] == 1.0
    assert tiered._eval_body_split(SUN, boundary, offset)[0][0] == 2.0
    assert base_child.calls[-1] == (SUN, boundary, -offset)
    assert medium_child.calls[-1] == (SUN, boundary, offset)


@pytest.mark.unit
def test_composite_and_tiered_dispatch_preserve_offset_and_boundary_side():
    offset = math.ulp(JD0) / 4.0
    base_child = _DispatchReader(JD0, JD0 + 1.0, 1.0)
    medium_child = _DispatchReader(JD0 + 1.0, JD0 + 2.0, 2.0)
    base = CompositeLEBReader([base_child])
    medium = CompositeLEBReader([medium_child])
    tiered = TieredLEBReader({"base": base, "medium": medium})

    assert base._eval_body_split(SUN, JD0 + 0.5, offset)[0][0] == 1.0
    assert base_child.calls[-1] == (SUN, JD0 + 0.5, offset)
    assert tiered._eval_body_split(SUN, JD0 + 1.0, -offset)[0][0] == 1.0
    assert tiered._eval_body_split(SUN, JD0 + 1.0, 0.0)[0][0] == 1.0
    assert tiered._eval_body_split(SUN, JD0 + 1.0, offset)[0][0] == 2.0
    assert tiered.selected_tier(SUN, JD0 + 1.0) == "base"


class _PipelineReader:
    def __init__(self) -> None:
        self._bodies: dict[int, Any] = {}
        self.split_calls: list[tuple[int, float, float]] = []

    def has_body(self, body_id: int) -> bool:
        return body_id in (SUN, EARTH)

    def eval_body(self, body_id: int, jd: float):
        if body_id == EARTH:
            return ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
        return ((1.0, 0.0, 0.0), (0.0, 0.0, 0.0))

    def _eval_body_split(self, body_id: int, jd: float, offset: float):
        self.split_calls.append((body_id, jd, offset))
        return self.eval_body(body_id, jd)


class _EvalOnlyPipelineReader(_PipelineReader):
    _eval_body_split = None


@pytest.mark.unit
def test_pipeline_rejects_noncanonical_eval_only_reader_before_losing_residual():
    reader = _EvalOnlyPipelineReader()
    flags = FLG_EQUATORIAL | FLG_J2000 | FLG_NOABERR | FLG_NOGDEFL

    with pytest.raises(TypeError, match="split-epoch"):
        _pipeline_icrs(reader, JD0, SUN, flags, want_xyz=True)


@pytest.mark.unit
def test_pipeline_dispatches_retarded_epoch_as_two_parts():
    reader = _PipelineReader()
    flags = FLG_EQUATORIAL | FLG_J2000 | FLG_NOABERR | FLG_NOGDEFL

    result = _pipeline_icrs(reader, JD0, SUN, flags, want_xyz=True)

    assert result == pytest.approx((1.0, 0.0, 0.0), abs=1e-7)
    rounded_jd = JD0 - 1.0 / C_LIGHT_AU_DAY
    residual = math.fsum((JD0, -1.0 / C_LIGHT_AU_DAY, -rounded_jd))
    assert reader.split_calls
    assert all(call == (SUN, rounded_jd, residual) for call in reader.split_calls)
    assert residual != 0.0
