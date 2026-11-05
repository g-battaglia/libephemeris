"""Selected maintenance overlays preserve real native-readable LEB1 images."""

import os
from dataclasses import replace

import numpy as np
import pytest

from tools.leb_generation import assemble
from tools.leb_generation import format as wire
from tools.leb_generation.read import read_leb1
from tools.tests.test_leb_generation import native_reader


START, END = 2451545.0, 2451553.0


def body(identifier, value, *, start=START, end=END, segments=2, degree=2):
    entry = wire.BodyEntry(identifier, wire.COORD_ICRS_BARY, segments,
                           start, end, (end - start) / segments, degree, 3, 0)
    coefficients = (np.arange(segments * 3 * (degree + 1), dtype="<f8") / 64 + value)
    coefficients = coefficients.reshape(segments, 3, degree + 1)
    coefficients[0, 0, 1] = -0.0
    coefficients[0, 0, 2] = np.nextafter(0.0, 1.0)
    return entry, coefficients


def fixture(path, bodies, auxiliary=None, *, start=START, end=END):
    """Hand-serialize the sources independently of the assembly writer."""
    auxiliary = dict(auxiliary or {})
    count = 2 + len(auxiliary)
    index_at = wire.HEADER_SIZE + count * wire.SECTION_DIR_SIZE
    coefficients_at = index_at + len(bodies) * wire.BODY_ENTRY_SIZE
    data = bytearray(coefficients_at)
    sections = [wire.SectionEntry(wire.SECTION_BODY_INDEX, index_at,
                                  len(bodies) * wire.BODY_ENTRY_SIZE)]
    for ordinal, (entry, coefficients) in enumerate(bodies):
        wire.write_body_entry(data, index_at + ordinal * wire.BODY_ENTRY_SIZE,
                              replace(entry, data_offset=len(data)))
        data.extend(coefficients.tobytes())
    sections.append(wire.SectionEntry(wire.SECTION_CHEBYSHEV, coefficients_at,
                                      len(data) - coefficients_at))
    for identifier, payload in auxiliary.items():
        sections.append(wire.SectionEntry(identifier, len(data), len(payload)))
        data.extend(payload)
    wire.write_header(data, wire.FileHeader(wire.MAGIC, wire.VERSION, count,
                                           len(bodies), start, end, START - 100, 0))
    for ordinal, section in enumerate(sections):
        wire.write_section_dir(data, wire.HEADER_SIZE + ordinal * wire.SECTION_DIR_SIZE, section)
    path.write_bytes(data)
    return path


def records(path):
    header, bodies, auxiliary = read_leb1(path)
    rows = {row.entry.body_id: (replace(row.entry, data_offset=0),
                                b"".join(values.tobytes() for _, values in row.chunks(1)))
            for row in bodies}
    sections = {}
    for identifier, section in auxiliary.items():
        with section.path.open("rb") as source:
            source.seek(section.offset)
            sections[identifier] = source.read(section.size)
    return header, rows, sections


def sources(tmp_path):
    main = fixture(tmp_path / "main.leb", [body(15, 15), body(0, 1), body(2, 2)])
    partial = fixture(tmp_path / "partial.leb", [body(2, 20), body(1, 10)])
    return main, partial


def test_selected_replacement_and_addition_preserve_unselected_words_and_main_auxiliaries(tmp_path):
    opaque = {wire.SECTION_ORBITAL_ELEMENTS: b"main orbital payload\x00\xff",
              1234: b"opaque main-only data"}
    main = fixture(tmp_path / "main.leb", [body(15, 15, start=START + 1, end=END - 1),
                                          body(2, 2), body(0, 1)], opaque)
    left = fixture(tmp_path / "left.leb", [body(2, 20, segments=4)],
                   {wire.SECTION_ORBITAL_ELEMENTS: b"deliberately unused"})
    right = fixture(tmp_path / "right.leb", [body(1, 10)], {4321: b"not copied"})
    original = {path: path.read_bytes() for path in (main, left, right)}
    output = tmp_path / "new-directory" / "updated.leb"

    assemble.replace_leb1(main, [left, right], output, expected_bodies={1, 2})

    header, actual, auxiliary = records(output)
    old_header, old, _ = records(main)
    assert header == replace(old_header, body_count=4)
    assert list(actual) == [0, 1, 2, 15]
    assert actual == {**old, **records(left)[1], **records(right)[1]}
    assert auxiliary == opaque
    assert original == {path: path.read_bytes() for path in original}
    result = native_reader(output)
    for path in (main, left, right):
        reader = native_reader(path)
        _, rows, _ = read_leb1(path)
        for row in rows:
            if path == main and row.entry.body_id == 2:
                continue
            for epoch in (row.entry.jd_start, (row.entry.jd_start + row.entry.jd_end) / 2,
                          row.entry.jd_end):
                before = np.asarray(reader.eval_body(row.entry.body_id, epoch), dtype="<f8")
                after = np.asarray(result.eval_body(row.entry.body_id, epoch), dtype="<f8")
                assert before.tobytes() == after.tobytes()
        reader.close()
    result.close()


def test_inplace_replacement_requires_explicit_overwrite_and_preserves_parts(tmp_path):
    main, partial = sources(tmp_path)
    main.chmod(0o640)
    original, part_bytes = main.read_bytes(), partial.read_bytes()
    with pytest.raises(FileExistsError):
        assemble.replace_leb1(main, [partial], main, expected_bodies={1, 2})
    assert main.read_bytes() == original
    unselected = records(main)[1][15]
    assemble.replace_leb1(main, [partial], main, expected_bodies={1, 2}, overwrite=True)
    assert records(main)[1][15] == unselected
    assert records(main)[1][2] == records(partial)[1][2]
    assert partial.read_bytes() == part_bytes
    assert main.stat().st_mode & 0o777 == 0o640
    assert sorted(path.name for path in tmp_path.iterdir()) == ["main.leb", "partial.leb"]


@pytest.mark.parametrize("kind", ["missing", "unexpected", "empty", "overlap", "duplicate-row", "domain"])
def test_invalid_replacement_inventory_and_domain_preserve_destination(tmp_path, kind):
    main, partial = sources(tmp_path)
    output = tmp_path / "destination.leb"
    output.write_bytes(b"previous destination")
    inputs, expected = [partial], {1, 2}
    if kind == "missing":
        expected.add(3)
    elif kind == "unexpected":
        expected.remove(1)
    elif kind == "empty":
        expected.clear()
    elif kind == "overlap":
        inputs.append(fixture(tmp_path / "overlap.leb", [body(2, 200)]))
    elif kind == "duplicate-row":
        fixture(partial, [body(2, 20), body(2, 200)])
    elif kind == "domain":
        fixture(partial, [body(1, 10), body(2, 20)], end=np.nextafter(END, np.inf))
    original = {path: path.read_bytes() for path in (main, *inputs, output)}
    with pytest.raises(ValueError):
        assemble.replace_leb1(main, inputs, output, expected_bodies=expected, overwrite=True)
    assert original == {path: path.read_bytes() for path in original}
    assert not list(tmp_path.glob(".*-replace-*"))


@pytest.mark.parametrize("kind", ["same-main", "main-hardlink", "partial-symlink", "output-partial", "output-hardlink", "output-symlink"])
def test_input_aliases_and_partial_destinations_are_refused(tmp_path, kind):
    main, partial = sources(tmp_path)
    output = tmp_path / "updated.leb"
    inputs = [partial]
    if kind == "same-main":
        inputs = [main]
    elif kind == "main-hardlink":
        alias = tmp_path / "main-alias.leb"
        os.link(main, alias)
        inputs = [alias]
    elif kind == "partial-symlink":
        alias = tmp_path / "partial-alias.leb"
        alias.symlink_to(partial)
        inputs.append(alias)
    elif kind == "output-partial":
        output = partial
    elif kind == "output-hardlink":
        os.link(partial, output)
    elif kind == "output-symlink":
        output.symlink_to(partial)
    original = {path: path.read_bytes() for path in (main, partial)}
    with pytest.raises(ValueError):
        assemble.replace_leb1(main, inputs, output, expected_bodies={1, 2}, overwrite=True)
    assert original == {path: path.read_bytes() for path in original}


@pytest.mark.parametrize("failure", ["native", "source"])
def test_late_validation_failures_cannot_replace_main_or_existing_output(tmp_path, monkeypatch, failure):
    main, partial = sources(tmp_path)
    before = main.read_bytes()
    original_verify = assemble.verify_native

    def reject_or_change(path, bodies):
        original_verify(path, bodies)
        if failure == "native":
            raise ValueError("injected native verification failure")
        with partial.open("ab") as stream:
            stream.write(b"changed after reading")

    monkeypatch.setattr(assemble, "verify_native", reject_or_change)
    with pytest.raises(ValueError, match="native verification|source changed"):
        assemble.replace_leb1(main, [partial], main, expected_bodies={1, 2}, overwrite=True)
    assert main.read_bytes() == before
    assert not list(tmp_path.glob(".*-replace-*"))


def test_atomic_no_overwrite_refuses_a_destination_created_after_validation(tmp_path, monkeypatch):
    main, partial = sources(tmp_path)
    output = tmp_path / "updated.leb"
    original_check = assemble.check_input_hashes

    def racing_destination(identities):
        original_check(identities)
        output.write_bytes(b"other writer won")

    monkeypatch.setattr(assemble, "check_input_hashes", racing_destination)
    with pytest.raises(FileExistsError):
        assemble.replace_leb1(main, [partial], output, expected_bodies={1, 2})
    assert output.read_bytes() == b"other writer won"
    assert not list(tmp_path.glob(".*-replace-*"))
