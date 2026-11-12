"""Offline format production must be accepted by the independent Rust reader."""

import struct
import subprocess
import sys
from pathlib import Path

import libephemeris as engine
import numpy as np
import pytest
from numpy.polynomial.chebyshev import chebval

from tools.leb_generation import compression
from tools.leb_generation import format as wire
from tools.leb_generation.convert import convert_leb1
from tools.leb_generation.writer import BodyData, write_leb2


def linear_container():
    """A synthetic three-component linear polynomial, with no astronomical data."""
    start, end = 2451545.0, 2451549.0
    raw = struct.pack("<6d", 1.0, 0.5, -2.0, -0.25, 3.0, 0.125)
    index_at = wire.HEADER_SIZE + 2 * wire.SECTION_DIR_SIZE
    coefficients_at = index_at + wire.BODY_ENTRY_SIZE
    data = bytearray(coefficients_at + len(raw))
    header = wire.FileHeader(wire.MAGIC, wire.VERSION, 2, 1, start, end, start, 0)
    body = wire.BodyEntry(0, wire.COORD_ICRS_BARY, 1, start, end, 4.0, 1, 3, coefficients_at)
    wire.write_header(data, header)
    for index, section in enumerate((
        wire.SectionEntry(wire.SECTION_BODY_INDEX, index_at, wire.BODY_ENTRY_SIZE),
        wire.SectionEntry(wire.SECTION_CHEBYSHEV, coefficients_at, len(raw)),
    )):
        wire.write_section_dir(data, wire.HEADER_SIZE + index * wire.SECTION_DIR_SIZE, section)
    wire.write_body_entry(data, index_at, body)
    data[coefficients_at:] = raw
    return data, header, body


def native_reader(path):
    """Select a private context without retaining a changed global setting."""
    previous = engine.get_calc_mode()
    try:
        engine.set_calc_mode("leb")
        context = engine.EphemerisContext(str(path.parent))
    finally:
        engine.set_calc_mode(previous)
    context.set_leb_file(str(path))
    return context.get_leb_reader()


def test_written_leb1_is_evaluated_by_the_native_reader(tmp_path):
    data, header, body = linear_container()
    path = tmp_path / "linear.leb"
    path.write_bytes(data)
    assert wire.read_header(data) == header
    assert wire.read_body_entry(data, wire.HEADER_SIZE + 2 * wire.SECTION_DIR_SIZE) == body
    reader = native_reader(path)
    assert reader.body_ids() == [0]
    assert reader.jd_range == (header.jd_start, header.jd_end)
    for jd in (header.jd_start, header.jd_start + 0.75, header.jd_start + 2.0, header.jd_end):
        t = 2.0 * (jd - header.jd_start) / body.interval_days - 1.0
        position, velocity = reader.eval_body(0, jd)
        assert position == [1.0 + 0.5 * t, -2.0 - 0.25 * t, 3.0 + 0.125 * t]
        assert velocity == [0.25, -0.125, 0.0625]


def test_native_reader_rejects_an_out_of_bounds_serialized_payload(tmp_path):
    data, _, body = linear_container()
    body.data_offset = len(data) + 8
    wire.write_body_entry(data, wire.HEADER_SIZE + 2 * wire.SECTION_DIR_SIZE, body)
    path = tmp_path / "invalid-offset.leb"
    path.write_bytes(data)
    # Container admission uses DataNotFoundError; chunk decompression has the
    # separate LEBCorruptionError. This fails before any coefficient access.
    with pytest.raises(engine.DataNotFoundError):
        native_reader(path)


@pytest.mark.parametrize("segment_count", [1, 7])
def test_compression_is_lossless_when_all_mantissa_bits_are_preserved(segment_count):
    values = np.random.default_rng(72).normal(size=(segment_count, 3, 6)).astype("<f8")
    raw = values.tobytes()
    blob = compression.compress_body(raw, segment_count, 5, 3, [52] * 6)
    assert compression.decompress_body(blob, len(raw), segment_count, 5, 3) == raw


def test_chunked_compression_preserves_the_final_partial_chunk():
    values = np.arange(7 * 3 * 4, dtype="<f8").reshape(7, 3, 4)
    chunks = compression.compress_body_chunked(values.tobytes(), 7, 3, 3, [52] * 4, 3)
    assert len(chunks) == 3
    decoded = []
    for (blob, byte_count), count in zip(chunks, (3, 3, 1), strict=True):
        decoded.append(compression.decompress_body(blob, byte_count, count, 3, 3))
    assert b"".join(decoded) == values.tobytes()


def test_lossy_compression_respects_the_polynomial_sum_error_bound():
    degree, segments, components = 9, 4, 3
    target = 1e-10
    coefficients = np.random.default_rng(41).normal(size=(segments, components, degree + 1))
    coefficients *= np.geomspace(1.0, 1e-12, degree + 1)
    coefficients = coefficients.astype("<f8")
    bits = compression.compute_mantissa_bits(coefficients, target)
    blob = compression.compress_body(coefficients.tobytes(), segments, degree, components, bits)
    decoded = compression.decompress_body(blob, coefficients.nbytes, segments, degree, components)
    restored = np.frombuffer(decoded, dtype="<f8").reshape(coefficients.shape)
    # |T_k(t)| <= 1 on [-1,1]: independent polynomial evaluation must respect
    # the sum of the fixed per-coefficient error budgets, including endpoints.
    nodes = np.linspace(-1.0, 1.0, 257)
    for segment in range(segments):
        for component in range(components):
            before = chebval(nodes, coefficients[segment, component])
            after = chebval(nodes, restored[segment, component])
            assert np.max(np.abs(after - before)) <= (degree + 1) * target


def test_small_coefficients_are_never_discarded_above_the_requested_error_budget():
    coefficients = np.array([[[1.0, 1e-35]]], dtype="<f8")
    target = 1e-40
    bits = compression.compute_mantissa_bits(coefficients, target)
    restored = compression.apply_truncation(coefficients, bits)
    assert np.max(np.abs(coefficients - restored)) <= target


@pytest.mark.parametrize("target", [0.0, -1.0, float("nan"), float("inf")])
def test_invalid_coefficient_error_budgets_are_rejected(target):
    with pytest.raises(ValueError):
        compression.compute_mantissa_bits(np.ones((1, 3, 2)), target)


@pytest.mark.parametrize("lossless", [False, True])
def test_conversion_produces_a_native_readable_leb2_without_changing_its_source(tmp_path, lossless):
    data, header, _ = linear_container()
    source = tmp_path / "source.leb"
    output = tmp_path / "converted.leb2"
    source.write_bytes(data)
    convert_leb1(source, output, lossless=lossless)
    assert source.read_bytes() == data
    reader = native_reader(output)
    for jd in (header.jd_start, header.jd_start + 0.75, header.jd_end):
        assert reader.eval_body(0, jd) == native_reader(source).eval_body(0, jd)


@pytest.mark.parametrize("chunk_days", [8.0, 9.5])
def test_streamed_leb2_uses_actual_chunk_boundaries_and_preserves_the_last_segment(tmp_path, chunk_days):
    _, header, entry = linear_container()
    entry.segment_count = 7
    entry.jd_end = entry.jd_start + 28.0
    header.jd_end = entry.jd_end
    coefficients = np.zeros((7, 3, 2), dtype="<f8")
    coefficients[:, :, 0] = np.arange(7)[:, None] + np.array([1.0, 2.0, 3.0])
    coefficients[:, :, 1] = np.array([0.5, 0.25, 0.125])
    output = tmp_path / "chunked.leb2"
    write_leb2(output, [BodyData(entry, coefficients=coefficients)], header,
               chunk_interval_days=chunk_days, lossless=True)
    reader = native_reader(output)
    for segment in range(7):
        jd = entry.jd_start + (segment + 0.375) * entry.interval_days
        position, velocity = reader.eval_body(0, jd)
        assert position == list(coefficients[segment, :, 0] - 0.25 * coefficients[segment, :, 1])
        assert velocity == [0.25, 0.125, 0.0625]
    assert reader.eval_body(0, entry.jd_end)[0] == list(coefficients[-1, :, 0] + coefficients[-1, :, 1])


def test_failed_native_verification_preserves_the_existing_destination(tmp_path):
    _, header, entry = linear_container()
    entry.jd_start = float("nan")  # Invalid body coverage must fail native admission.
    output = tmp_path / "existing.leb2"
    output.write_bytes(b"keep this existing file")
    invalid = BodyData(entry, coefficients=np.ones((1, 3, 2), dtype="<f8"))
    with pytest.raises(ValueError, match="native generated-file verification failed"):
        write_leb2(output, [invalid], header, overwrite=True)
    assert output.read_bytes() == b"keep this existing file"
    assert list(tmp_path.iterdir()) == [output]


def test_conversion_refuses_to_overwrite_its_input(tmp_path):
    data, _, _ = linear_container()
    source = tmp_path / "source.leb"
    source.write_bytes(data)
    with pytest.raises(ValueError, match="output must differ"):
        convert_leb1(source, source, overwrite=True)
    assert source.read_bytes() == data


def test_conversion_cli_is_self_contained_and_preserves_existing_output(tmp_path):
    data, _, _ = linear_container()
    source = tmp_path / "source.leb"
    output = tmp_path / "converted.leb2"
    source.write_bytes(data)
    command = [sys.executable, "-B", "-m", "leb_generation.cli_leb2",
               "convert", str(source), "--output", str(output)]
    subprocess.run(command, capture_output=True, text=True, check=True)
    assert native_reader(output).body_ids() == [0]
    original = output.read_bytes()
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 1
    assert output.read_bytes() == original
    assert source.read_bytes() == data


def test_conversion_cli_rejects_a_truncated_source_without_publishing(tmp_path):
    source = tmp_path / "truncated.leb"
    output = tmp_path / "missing.leb2"
    source.write_bytes(b"LEB1")
    command = [sys.executable, "-B", "-m", "leb_generation.cli_leb2",
               "convert", str(source), "--output", str(output)]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    assert result.returncode == 1
    assert "truncated LEB header" in result.stderr
    assert not output.exists()
