"""Native evaluation of newly fitted coefficients, coverage and atomic assembly."""

from dataclasses import replace
from pathlib import Path
from typing import ClassVar

import numpy as np
import pytest

from tools.leb_generation import format as wire
from tools.leb_generation.assemble import merge_leb1, write_leb1
from tools.leb_generation.convert import convert_all_leb1, convert_leb1
from tools.leb_generation.generate import fit_grid, generate_leb1, generate_leb2, generation_header_start
from tools.leb_generation.groups import LEB2_GROUPS
from tools.leb_generation.read import read_body_index, read_leb1
from tools.leb_generation.source import chain_windows
from tools.leb_generation.verify import verify_source_file
from tools.leb_generation.writer import BodyData, write_leb2
from tools.tests.test_leb_generation import linear_container, native_reader


class LinearSource:
    """Synthetic state with analytic velocity; no reference astronomy data."""

    kernel = Path("/nonexistent-synthetic-source.bsp")
    input_paths = frozenset((kernel,))
    piecewise_grids: ClassVar[dict] = {}

    def input_hashes(self):
        return {}  # This synthetic source has no file-backed inputs.

    def require_domain(self, body_id, start, end):
        if not 2451545.0 <= start < end <= 2451549.0:
            raise ValueError("synthetic source domain exceeded")

    def states(self, body_id, epochs):
        t = np.asarray(epochs) - 2451545.0
        return np.column_stack((1 + 0.01 * t, -2 - 0.02 * t, 3 + 0.03 * t,
                                np.full_like(t, 0.01), np.full_like(t, -0.02), np.full_like(t, 0.03)))

    def values(self, body_id, epochs):
        return self.states(body_id, epochs)[:, :3]


def test_discontinuous_source_rate_refusal_preserves_destination_and_explains_failure(tmp_path):
    class DiscontinuousSource(LinearSource):
        def states(self, body_id, epochs):
            result = super().states(body_id, epochs)
            # A small position jump leaves each side's velocity unchanged.
            # It cannot be represented by one smooth position/rate polynomial.
            result[:, 1] += np.where(np.asarray(epochs) >= 2451546.75, 2e-8, 0.0)
            return result

    output = tmp_path / "existing.leb"
    output.write_bytes(b"preserved")
    with pytest.raises(ValueError, match="fitted Cartesian rate exceeds PLAN relative criterion") as failure:
        generate_leb1(DiscontinuousSource(), output, [0], 2451545.0, 2451549.0,
                      skip_aux=True, overwrite=True)
    diagnostic = str(failure.value)
    for field in ("body 0:", "segment=0", "start_tt=2451545", "interval_days=4",
                  "jd_tt=", "error_au_per_tt_day=", "limit_au_per_tt_day=", "fitted=", "source="):
        assert field in diagnostic
    assert output.read_bytes() == b"preserved"
    assert list(tmp_path.iterdir()) == [output]


@pytest.mark.parametrize("arguments", [
    ["--nbody", "CHIRON"],
    ["--nbody", "CHIRON", "--assist-planets", "/unused-planets"],
    ["--assist-planets", "/unused-planets", "--assist-asteroids", "/unused-asteroids"],
])
def test_nbody_cli_refuses_incomplete_or_unused_force_selections(arguments, capsys):
    from tools.leb_generation.cli_leb1 import main

    with pytest.raises(SystemExit) as result:
        main(arguments)
    assert result.value.code == 1
    assert "--nbody" in capsys.readouterr().err


def test_native_reader_evaluates_new_source_fits_at_both_endpoints_and_interior(tmp_path):
    output = tmp_path / "new.leb"
    checks = generate_leb1(LinearSource(), output, [0], 2451545.0, 2451549.0, skip_aux=True)
    assert checks.positions == checks.cartesian_rates == 9
    reader = native_reader(output)
    for epoch in (2451545.0, 2451545.123456789, 2451547.5, 2451549.0):
        actual = reader.eval_body(0, epoch)
        expected = LinearSource().states(0, [epoch])[0]
        np.testing.assert_allclose(actual[0], expected[:3], rtol=1e-13, atol=0)
        np.testing.assert_allclose(actual[1], expected[3:], rtol=1e-11, atol=0)


def test_direct_leb2_generation_and_source_verification_use_final_native_values(tmp_path):
    output = tmp_path / "new.leb2"
    checks = generate_leb2(LinearSource(), output, [0], 2451545.0, 2451549.0, skip_aux=True)
    assert checks.positions == checks.cartesian_rates == 9
    header, entries = read_body_index(output)
    assert header.magic == wire.LEB2_MAGIC and header.version == 2
    assert [entry.body_id for entry in entries] == [0]
    assert verify_source_file(output, LinearSource(), 17) == {"positions": 17, "cartesian_rates": 17}

    class WrongRateSource(LinearSource):
        def states(self, body_id, epochs):
            result = super().states(body_id, epochs)
            result[:, 3:] *= 2
            return result

    with pytest.raises(ValueError, match="Cartesian rate"):
        verify_source_file(output, WrongRateSource(), 17)


def test_final_leb2_source_validation_failure_preserves_existing_output(tmp_path, monkeypatch):
    from tools.leb_generation import verify

    original = verify.verify_source_file
    checked = []

    def reject_compressed(path, source, samples):
        magic = read_body_index(path)[0].magic
        checked.append(magic)
        if magic == wire.LEB2_MAGIC:
            raise ValueError("final compressed source criterion failed")
        return original(path, source, samples)

    monkeypatch.setattr(verify, "verify_source_file", reject_compressed)
    output = tmp_path / "existing.leb2"
    output.write_bytes(b"preserved")
    with pytest.raises(ValueError, match="final compressed source criterion"):
        generate_leb2(LinearSource(), output, [0], 2451545.0, 2451549.0,
                      skip_aux=True, overwrite=True)
    assert checked == [wire.MAGIC, wire.LEB2_MAGIC]
    assert output.read_bytes() == b"preserved"
    assert list(tmp_path.iterdir()) == [output]


def test_leb2_verify_cli_distinguishes_structure_from_source_precision(tmp_path, capsys):
    from tools.leb_generation.cli_leb2 import main

    output = tmp_path / "new.leb2"
    generate_leb2(LinearSource(), output, [0], 2451545.0, 2451549.0, skip_aux=True)
    assert main(["verify", str(output), "--structural-only"]) == 0
    assert "source precision not checked" in capsys.readouterr().out
    with pytest.raises(SystemExit) as result:
        main(["verify", str(output)])
    assert result.value.code == 1
    assert "--kernel or --data-dir" in capsys.readouterr().err
    with pytest.raises(SystemExit) as result:
        main(["verify", str(output), "--structural-only", "--kernel", "/unused"])
    assert result.value.code == 1
    assert "does not accept astronomical source selections" in capsys.readouterr().err


def test_leb2_verify_cli_checks_unsampled_compressed_payload(tmp_path, capsys):
    from tools.leb_generation.cli_leb2 import main

    output = tmp_path / "corrupt.leb2"
    _, header, body = linear_container()
    values = np.ones((3, 3, 2), dtype="<f8")
    body = replace(body, segment_count=3, interval_days=4.0 / 3.0)
    write_leb2(output, [BodyData(body, coefficients=values)], header, chunk_interval_days=1.5)
    # Keep all container/chunk metadata intact and break the Zstandard magic.
    # Admission alone is insufficient: verification must decompress the chunk.
    entry = read_body_index(output)[1][0]
    with output.open("r+b") as stream:
        stream.seek(entry.data_offset + wire.CHUNK_INDEX_HEADER_SIZE + wire.CHUNK_ENTRY_SIZE)
        chunk = wire.read_chunk_entry(stream.read(wire.CHUNK_ENTRY_SIZE), 0)
        stream.seek(chunk.blob_offset)
        stream.write(b"BAD!")
    with pytest.raises(SystemExit) as result:
        main(["verify", str(output), "--structural-only"])
    assert result.value.code == 1
    assert "verification failed" in capsys.readouterr().err


def test_periodic_longitude_fit_remains_continuous_across_zero(tmp_path):
    class AngularSource(LinearSource):
        def values(self, body_id, epochs):
            t = np.asarray(epochs) - 2451545.0
            return np.column_stack(((359 + t) % 360, np.full_like(t, -2.0), np.full_like(t, 0.003)))

    output = tmp_path / "angular.leb"
    checks = generate_leb1(AngularSource(), output, [10], 2451545.0, 2451549.0, skip_aux=True)
    assert checks.positions == 9 and checks.cartesian_rates == 0
    reader = native_reader(output)
    for epoch in (2451545.0, 2451545.5, 2451546.0, 2451546.5, 2451549.0):
        position, velocity = reader.eval_body(10, epoch)
        expected = AngularSource().values(10, [epoch])[0]
        assert abs((position[0] - expected[0] + 180) % 360 - 180) < 1e-10
        assert abs(velocity[0] - 1.0) < 1e-10


def test_piecewise_source_knots_are_boundaries_for_arbitrary_requested_starts(tmp_path):
    class PiecewiseSource(LinearSource):
        piecewise_grids: ClassVar[dict] = {21: (2451544.5, 2451550.5, 2.0, 365.25)}

        def values(self, body_id, epochs):
            t = np.asarray(epochs) - 2451546.5
            return np.column_stack((100.0 + 0.1 * np.abs(t), np.full_like(t, -2.0), np.full_like(t, 0.003)))

    source = PiecewiseSource()
    output = tmp_path / "piecewise.leb"
    start, end = 2451545.123456789, 2451549.0
    generate_leb1(source, output, [21], start, end, skip_aux=True)
    header, rows, _ = read_leb1(output)
    entry = rows[0].entry
    assert header.jd_start == entry.jd_start <= start
    assert entry.jd_end == end
    assert entry.interval_days == 0.25
    assert ((2451546.5 - entry.jd_start) / entry.interval_days).is_integer()
    reader = native_reader(output)
    # Both sides and the exact derivative break: a polynomial is never fitted
    # across this kink. The requested fractional start is still fully covered.
    epochs = [start, 2451546.49, 2451546.5, 2451546.51, end]
    for epoch in epochs:
        actual, _ = reader.eval_body(21, epoch)
        np.testing.assert_allclose(actual, source.values(21, [epoch])[0], rtol=1e-13, atol=0)


def test_group_header_union_preserves_actual_body_domains_and_merged_words(tmp_path):
    class MixedSource(LinearSource):
        piecewise_grids: ClassVar[dict] = {21: (2451544.5, 2451550.5, 2.0, 365.25)}

        def values(self, body_id, epochs):
            if body_id != 21:
                return super().values(body_id, epochs)
            t = np.asarray(epochs) - 2451546.5
            return np.column_stack((100.0 + 0.1 * np.abs(t), np.full_like(t, -2.0), np.full_like(t, 0.003)))

    source = MixedSource()
    start, end = 2451545.123456789, 2451549.0
    common = generation_header_start(source, {0: (start, end), 21: (start, end)}, start)
    assert common < start
    parts = [tmp_path / "planets.leb", tmp_path / "analytical.leb"]
    for body, path in zip((0, 21), parts):
        generate_leb1(source, path, [body], start, end, skip_aux=True, header_start=common)
    output = tmp_path / "merged.leb"
    merge_leb1(parts, output)
    header, entries = read_body_index(output)
    assert header.jd_start == common
    assert [entry.jd_start for entry in entries] == [start, common]
    merged = native_reader(output)
    import libephemeris as engine
    with pytest.raises(engine.Error):
        merged.eval_body(0, common)
    for body, path in zip((0, 21), parts):
        original = native_reader(path)
        for epoch in (start, 2451546.5, end):
            np.testing.assert_array_equal(merged.eval_body(body, epoch), original.eval_body(body, epoch))


@pytest.mark.parametrize("header_start", [float("nan"), 2451545.5])
def test_invalid_shared_header_preserves_an_existing_output(tmp_path, header_start):
    output = tmp_path / "existing.leb"
    output.write_bytes(b"preserved")
    with pytest.raises(ValueError, match="header start"):
        generate_leb1(LinearSource(), output, [0], 2451545.0, 2451549.0,
                      skip_aux=True, overwrite=True, header_start=header_start)
    assert output.read_bytes() == b"preserved"


@pytest.mark.parametrize("span,preferred", [(1.0, 8.0), (97.0, 8.0), (10000.125, 64.0)])
def test_fit_grid_covers_requested_domain_without_an_extra_segment(span, preferred):
    count, interval = fit_grid(2451545.0, 2451545.0 + span, preferred)
    assert interval <= preferred
    assert np.ceil(span / interval) == count
    assert abs(count * interval - span) <= 4 * np.spacing(span)


def test_source_coverage_intersects_centers_and_preserves_actual_gaps():
    rows = [(301, 3, 1, 0.0, 100.0), (3, 0, 1, 20.0, 40.0), (3, 0, 1, 60.0, 80.0)]
    assert chain_windows(rows, 301) == [(20.0, 40.0), (60.0, 80.0)]
    assert chain_windows(rows, 999) == []
    with pytest.raises(ValueError, match="cyclic"):
        chain_windows([(301, 3, 1, 0.0, 100.0), (3, 301, 1, 0.0, 100.0)], 301)
    with pytest.raises(ValueError, match="J2000"):
        chain_windows([(301, 0, 17, 0.0, 100.0)], 301)


def test_failed_source_fit_preserves_the_existing_destination(tmp_path):
    class InvalidSource(LinearSource):
        def values(self, body_id, epochs):
            return super().values(body_id, epochs) * float("nan")

    output = tmp_path / "existing.leb"
    output.write_bytes(b"existing output")
    with pytest.raises(ValueError, match="invalid native source"):
        generate_leb1(InvalidSource(), output, [0], 2451545.0, 2451549.0, skip_aux=True, overwrite=True)
    assert output.read_bytes() == b"existing output"
    assert list(tmp_path.iterdir()) == [output]


@pytest.mark.parametrize("generate", [generate_leb1, generate_leb2])
def test_generation_cannot_overwrite_a_seed_or_force_image(tmp_path, generate):
    selected_input = tmp_path / "selected-force-file"
    selected_input.write_bytes(b"immutable selected source")

    class ProtectedSource(LinearSource):
        input_paths = frozenset((selected_input, LinearSource.kernel))

    with pytest.raises(ValueError, match="every source and force file"):
        generate(ProtectedSource(), selected_input, [0], 2451545.0, 2451549.0,
                 skip_aux=True, overwrite=True)
    assert selected_input.read_bytes() == b"immutable selected source"




def test_source_domain_is_not_clamped_or_extrapolated(tmp_path):
    output = tmp_path / "outside.leb"
    with pytest.raises(ValueError, match="domain exceeded"):
        generate_leb1(LinearSource(), output, [0], 2451544.0, 2451549.0, skip_aux=True)
    assert not output.exists()


def test_merged_native_containers_keep_all_coefficients(tmp_path):
    _, header, entry = linear_container()
    values = np.array([[[1.0, 0.5], [-2.0, -0.25], [3.0, 0.125]]], dtype="<f8")
    left, right, output = (tmp_path / name for name in ("left.leb", "right.leb", "merged.leb"))
    write_leb1(left, [BodyData(entry, coefficients=values)], header)
    write_leb1(right, [BodyData(replace(entry, body_id=1), coefficients=values * 2)], header)
    originals = [path.read_bytes() for path in (left, right)]
    merge_leb1([left, right], output)
    merged = native_reader(output)
    assert merged.body_ids() == [0, 1]
    for body, source in ((0, left), (1, right)):
        reader = native_reader(source)
        for epoch in (header.jd_start, header.jd_start + 1.25, header.jd_end):
            assert merged.eval_body(body, epoch) == reader.eval_body(body, epoch)
    assert [path.read_bytes() for path in (left, right)] == originals
    _, bodies, _ = read_leb1(output)
    assert {body.entry.body_id for body in bodies} == {0, 1}
    with pytest.raises(ValueError, match="duplicate"):
        merge_leb1([left, left], output, overwrite=True)
    assert merged.eval_body(0, header.jd_start) == native_reader(left).eval_body(0, header.jd_start)


def test_merge_rejects_conflicting_auxiliary_data_and_preserves_the_output(tmp_path):
    _, header, entry = linear_container()
    values = np.ones((1, 3, 2), dtype="<f8")
    left, right, output = (tmp_path / name for name in ("left.leb", "right.leb", "output.leb"))
    # The reserved orbital-elements section is opaque to the reader. A merge
    # must still preserve matching bytes and refuse conflicting snapshots.
    write_leb1(left, [BodyData(entry, coefficients=values)], header, {wire.SECTION_ORBITAL_ELEMENTS: b"first"})
    write_leb1(right, [BodyData(replace(entry, body_id=1), coefficients=values)], header,
               {wire.SECTION_ORBITAL_ELEMENTS: b"other"})
    output.write_bytes(b"original")
    with pytest.raises(ValueError, match="conflicting shared"):
        merge_leb1([left, right], output, overwrite=True)
    assert output.read_bytes() == b"original"


@pytest.mark.parametrize("write", [write_leb1, write_leb2])
def test_extra_validator_failure_cannot_publish_or_replace(tmp_path, write):
    _, header, entry = linear_container()
    values = np.ones((1, 3, 2), dtype="<f8")
    output = tmp_path / "existing.leb"
    output.write_bytes(b"preserved")

    def reject(path):
        assert path != output and native_reader(path).has_body(0)
        raise ValueError("source criterion failed")

    with pytest.raises(ValueError, match="criterion failed"):
        write(output, [BodyData(entry, coefficients=values)], header, overwrite=True, validator=reject)
    assert output.read_bytes() == b"preserved"
    assert list(tmp_path.iterdir()) == [output]


def test_split_conversion_requires_complete_inventories_before_publishing(tmp_path):
    _, header, entry = linear_container()
    values = np.ones((1, 3, 2), dtype="<f8")
    source, directory = tmp_path / "core.leb", tmp_path / "split"
    write_leb1(source, [BodyData(entry, coefficients=values)], header)
    with pytest.raises(ValueError, match="inventory"):
        convert_all_leb1(source, directory, "base")
    assert not directory.exists()
    rows = [BodyData(replace(entry, body_id=body), coefficients=values)
            for body in sorted(LEB2_GROUPS["core"])]
    write_leb1(source, rows, header, overwrite=True)
    outputs = convert_all_leb1(source, directory, "base")
    assert outputs == [directory / "base_core.leb2"]
    reader = native_reader(outputs[0])
    assert set(reader.body_ids()) == LEB2_GROUPS["core"]
    previous = outputs[0].read_bytes()
    with pytest.raises(FileExistsError):
        convert_all_leb1(source, directory, "base")
    assert outputs[0].read_bytes() == previous
