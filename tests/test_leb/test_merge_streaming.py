"""Bounded-memory and atomicity contracts of the LEB1 merge workflow."""

from __future__ import annotations

import builtins
from dataclasses import replace
import io
from pathlib import Path

import pytest

from libephemeris import leb_format as fmt
from libephemeris.leb_reader import LEBReader
from scripts import generate_leb as generator
from tests.leb_build_helpers import write_synthetic_leb


def _payload(path: Path, body_id: int) -> bytes:
    with LEBReader(str(path)) as reader:
        body = reader._bodies[body_id]
        size = body.segment_count * fmt.segment_byte_size(body.degree, body.components)
        return bytes(reader._mm[body.data_offset : body.data_offset + size])


def test_merge_preserves_payloads_and_uses_aux_only_source(tmp_path):
    aux, first, second, output = [
        tmp_path / name for name in ("sun", "first", "second", "merged")
    ]
    write_synthetic_leb(aux, (0,))
    write_synthetic_leb(first, (17,), aux=False)
    write_synthetic_leb(second, (18,), aux=False)
    generator.merge_leb_files(
        [str(first), str(second)], str(output), False, aux_source=str(aux)
    )
    with LEBReader(str(output)) as reader, LEBReader(str(aux)) as source:
        assert set(reader._bodies) == {17, 18}
        assert reader.has_nutation()
        assert reader.delta_t(105) == source.delta_t(105)
        assert reader.eval_nutation(105) == source.eval_nutation(105)
        assert reader.get_star(1) == source.get_star(1)
    assert _payload(first, 17) == _payload(output, 17)
    assert _payload(second, 18) == _payload(output, 18)


def test_empty_aux_placeholders_do_not_hide_real_sections(tmp_path):
    first, second, output = [tmp_path / name for name in ("empty", "real", "merged")]
    write_synthetic_leb(first, (17,), aux=False)
    write_synthetic_leb(second, (0,))
    generator.merge_leb_files([str(first), str(second)], str(output), False)
    with LEBReader(str(output)) as reader:
        assert reader.has_nutation()
        assert reader._delta_t_jds == [100, 110]


@pytest.mark.parametrize(
    "failure", ["duplicate", "range", "section", "payload", "truncated"]
)
def test_rejected_inputs_leave_previous_output_intact(tmp_path, failure):
    first, second, output = [tmp_path / name for name in ("first", "second", "merged")]
    write_synthetic_leb(first, (0,))
    write_synthetic_leb(
        second,
        (0,) if failure == "duplicate" else (17,),
        start=101 if failure == "range" else 100,
    )
    data = bytearray(second.read_bytes())
    if failure == "section":
        section = fmt.read_section_dir(data, fmt.HEADER_SIZE)
        fmt.write_section_dir(
            data, fmt.HEADER_SIZE, replace(section, offset=len(data) + 1)
        )
    elif failure == "payload":
        section = fmt.read_section_dir(data, fmt.HEADER_SIZE)
        body = fmt.read_body_entry(data, section.offset)
        fmt.write_body_entry(
            data, section.offset, replace(body, data_offset=len(data) + 1)
        )
    elif failure == "truncated":
        data = data[:-1]
    second.write_bytes(data)
    output.write_bytes(b"previous output")
    with pytest.raises(ValueError):
        generator.merge_leb_files([str(first), str(second)], str(output), False)
    assert output.read_bytes() == b"previous output"
    assert not list(tmp_path.glob(".leb-merge-*"))


def test_copy_failure_is_atomic(tmp_path, monkeypatch):
    source, output = tmp_path / "source", tmp_path / "output"
    write_synthetic_leb(source)
    output.write_bytes(b"old")

    def fail(*args):
        raise OSError("disk full")

    monkeypatch.setattr(generator, "_copy_leb_payload", fail)
    with pytest.raises(OSError, match="disk full"):
        generator.merge_leb_files([str(source)], str(output), False)
    assert output.read_bytes() == b"old"
    assert not list(tmp_path.glob(".leb-merge-*"))


def test_merge_never_reads_an_entire_large_file(tmp_path, monkeypatch):
    source, output = tmp_path / "source", tmp_path / "output"
    write_synthetic_leb(source, segments=100_000)
    reads = []

    class GuardedFile:
        def __init__(self, stream):
            self.stream = stream

        def __enter__(self):
            return self

        def __exit__(self, *args):
            self.stream.close()

        def __getattr__(self, name):
            return getattr(self.stream, name)

        def read(self, size=-1):
            assert 0 <= size <= 1024 * 1024
            reads.append(size)
            return self.stream.read(size)

    def guarded_open(path, mode):
        return GuardedFile(builtins.open(path, mode))

    monkeypatch.setattr(generator, "open", guarded_open, raising=False)
    generator.merge_leb_files([str(source)], str(output), False)
    assert max(reads) == 1024 * 1024
    assert _payload(source, 0) == _payload(output, 0)


def test_copy_detects_short_reads():
    with pytest.raises(ValueError, match="Short read"):
        generator._copy_leb_payload(io.BytesIO(b"abc"), io.BytesIO(), 0, 5)


def test_aux_source_must_have_real_data_and_matching_range(tmp_path):
    source, aux, output = [tmp_path / name for name in ("source", "aux", "output")]
    write_synthetic_leb(source, (17,))
    write_synthetic_leb(aux, (0,), aux=False)
    with pytest.raises(ValueError, match="Auxiliary source"):
        generator.merge_leb_files(
            [str(source)], str(output), False, aux_source=str(aux)
        )
    write_synthetic_leb(aux, start=101)
    with pytest.raises(ValueError, match="JD range mismatch"):
        generator.merge_leb_files(
            [str(source)], str(output), False, aux_source=str(aux)
        )


def test_merge_of_body_only_inputs_remains_readable(tmp_path):
    source, output = tmp_path / "body", tmp_path / "merged"
    write_synthetic_leb(source, (17,), aux=False)
    generator.merge_leb_files([str(source)], str(output), False)
    with LEBReader(str(output)) as reader:
        assert reader.has_body(17)
        assert not reader.has_nutation()
        assert reader.eval_body(17, 105)[0][0] == pytest.approx(17.099)


def test_empty_partial_does_not_break_legacy_single_body_workflow(tmp_path):
    empty, real, output = [tmp_path / name for name in ("empty", "real", "merged")]
    # Legacy --single may produce an empty partial for a tier-excluded body.
    write_synthetic_leb(empty, (), aux=False)
    write_synthetic_leb(real)
    generator.merge_leb_files([str(empty), str(real)], str(output), False)
    with LEBReader(str(output)) as reader:
        assert set(reader._bodies) == {0}
        assert reader.has_nutation()


def test_source_change_during_copy_prevents_promotion(tmp_path, monkeypatch):
    source, output = tmp_path / "source", tmp_path / "output"
    write_synthetic_leb(source)
    output.write_bytes(b"previous")
    original = generator._copy_leb_payload
    changed = False

    def edit_after_copy(stream, target, offset, size):
        nonlocal changed
        original(stream, target, offset, size)
        if not changed:
            changed = True
            with source.open("ab") as writer:
                writer.write(b"edited")

    monkeypatch.setattr(generator, "_copy_leb_payload", edit_after_copy)
    with pytest.raises(ValueError, match="input changed"):
        generator.merge_leb_files([str(source)], str(output), False)
    assert output.read_bytes() == b"previous"
