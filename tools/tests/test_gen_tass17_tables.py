"""Standalone lossless table serialization; no TASS astronomical evaluation."""
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys

import pytest

from tools.generators import gen_tass17_tables as generator


ROOT = Path(__file__).resolve().parents[2]
RUST = ROOT / "core/src/tables/tass17_tables.rs"


def catalog():
    return generator.load_catalog(generator.DEFAULT_INPUT)[0]


def words(values):
    return b"".join(struct.pack("<d", value) for value in values)


def test_current_artifact_and_every_periodic_word_index_reference_match_snapshot():
    data, sha = generator.load_catalog(generator.DEFAULT_INPUT)
    source = RUST.read_text()
    assert generator.render_tables(data).encode() == RUST.read_bytes()
    assert sha == hashlib.sha256(generator.DEFAULT_INPUT.read_bytes()).hexdigest()
    term_count, series_count = 0, 0
    for leaf in data["term_lists"]:
        block = re.search(rf'pub static {leaf["name"]}: &\[TassTerm\] =\s*&\[(.*?)\];', source, re.S)
        assert block is not None
        actual = [float(value) for value in re.findall(r"[-+]?(?:\d+\.\d*|\d+)(?:e[-+]?\d+)?", block[1])]
        expected = [value for term in leaf["terms"] for value in term]
        assert words(actual) == words(expected)
        term_count += len(leaf["terms"])
    for group in data["series"]:
        block = re.search(rf'pub static {group["name"]}: &\[TassMultiTerm\] = &\[(.*?)\];', source, re.S)
        assert block is not None
        actual = re.findall(r"\(\[([^]]+)\], ([A-Z]+_[0-3]_\d+)\)", block[1])
        assert [[[int(value) for value in indices.split(",")], name] for indices, name in actual] == group["rows"]
        series_count += len(group["rows"])
    assert (term_count, series_count, len(data["bodies"]), len(data["rotation"])) == (1535, 299, 8, 9)
    for body in data["bodies"]:
        block = re.search(rf'name: "{body["name"]}",\s*mu: ([^,]+),\s*aam: ([^,]+),\s*s0: \[(.*?)\],', source, re.S)
        assert block is not None
        actual = [float(block[1]), float(block[2])] + [float(value) for value in block[3].split(",") if value.strip()]
        assert words(actual) == words([body["mu"], body["aam"]] + body["s0"])
    block = re.search(r"pub static TASS17_TO_VSOP87: &\[f64; 9\] = &\[(.*?)\];", source, re.S)
    assert block is not None
    assert words([float(value) for value in block[1].split(",") if value.strip()]) == words(data["rotation"])


def test_isolated_generator_needs_only_its_snapshot_without_git_or_engine(tmp_path):
    directory = tmp_path / "standalone/tools/generators"
    (directory / "data").mkdir(parents=True)
    script = directory / "gen_tass17_tables.py"
    shutil.copyfile(ROOT / "tools/generators/gen_tass17_tables.py", script)
    shutil.copyfile(generator.DEFAULT_INPUT, directory / "data/tass17.json")
    output = tmp_path / "candidate.rs"
    environment = {**os.environ, "PATH": "", "PYTHONPATH": "", "PYTHONDONTWRITEBYTECODE": "1"}
    command = [sys.executable, "-I", "-B", str(script), "--out", str(output)]
    for _ in range(2):
        result = subprocess.run(command, env=environment, cwd=tmp_path, capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        assert output.read_bytes() == RUST.read_bytes()
    check = subprocess.run(command + ["--check"], env=environment, cwd=tmp_path, capture_output=True, text=True)
    assert check.returncode == 0, check.stderr
    assert "299 term lists, 1535 terms, 299 series rows" in check.stdout
    assert sorted(path.name for path in directory.iterdir()) == ["data", "gen_tass17_tables.py"]


def test_explicit_input_copy_and_json_object_reordering_preserve_identical_output(tmp_path):
    data = catalog()
    # JSON object order and whitespace are not part of the numerical identity;
    # array/body/series order remains intentional.
    source = tmp_path / "input.json"
    source.write_text(json.dumps(data, sort_keys=True))
    output = tmp_path / "out.rs"
    assert generator.main(["--input", str(source), "--out", str(output)]) == 0
    assert output.read_bytes() == RUST.read_bytes()
    assert generator.main(["--input", str(source), "--out", str(output), "--check"]) == 0


@pytest.mark.parametrize("failure", ["schema", "schema-bool", "lineage", "blob", "concat", "license", "note",
                                     "leaf-count", "leaf-order", "leaf-duplicate", "leaf-name", "term-count",
                                     "term-shape", "term-bool", "term-integer", "term-nonfinite", "term-change",
                                     "series-order", "indices-shape", "indices-bool", "indices-range", "reference",
                                     "cross-series", "missing-series", "body-order", "body-seed", "body-mu",
                                     "body-wiring", "rotation-shape", "rotation-nonfinite", "rotation-change", "extra"])
def test_malformed_or_changed_realization_refuses_before_writing(failure, tmp_path):
    data = catalog()
    if failure == "schema":
        data["schema_version"] = 2
    elif failure == "schema-bool":
        data["schema_version"] = True
    elif failure == "lineage":
        data["lineage"]["seed_commit"] = "0"*40
    elif failure == "blob":
        data["lineage"]["seed_blobs_sha256"][next(iter(generator.BLOB_SHAS))] = "g"*64
    elif failure == "concat":
        data["lineage"]["seed_concatenated_sha256"] = "0"*64
    elif failure == "license":
        data["lineage"]["license"] = "unknown"
    elif failure == "note":
        data["lineage"]["note"] = ""
    elif failure == "leaf-count":
        data["term_lists"].pop()
    elif failure == "leaf-order":
        data["term_lists"].reverse()
    elif failure == "leaf-duplicate":
        data["term_lists"][1] = deepcopy(data["term_lists"][0])
    elif failure == "leaf-name":
        data["term_lists"][0]["name"] = "not_a_rust_identifier()"
    elif failure == "term-count":
        data["term_lists"][0]["terms"].append([0.0, 0.0, 0.0])
    elif failure == "term-shape":
        data["term_lists"][0]["terms"][0].pop()
    elif failure == "term-bool":
        data["term_lists"][0]["terms"][0][0] = True
    elif failure == "term-integer":
        data["term_lists"][0]["terms"][0][0] = 0
    elif failure == "term-nonfinite":
        data["term_lists"][0]["terms"][0][0] = float("inf")
    elif failure == "term-change":
        data["term_lists"][0]["terms"][0][0] += 1e-15
    elif failure == "series-order":
        data["series"].reverse()
    elif failure == "indices-shape":
        data["series"][0]["rows"][0][0].pop()
    elif failure == "indices-bool":
        data["series"][0]["rows"][0][0][0] = True
    elif failure == "indices-range":
        data["series"][0]["rows"][0][0][0] = 128
    elif failure == "reference":
        data["series"][0]["rows"][0][1] = "MISSING_0_0"
    elif failure == "cross-series":
        data["series"][0]["rows"][0][1] = data["series"][1]["rows"][0][1]
    elif failure == "missing-series":
        data["series"].pop()
    elif failure == "body-order":
        data["bodies"].reverse()
    elif failure == "body-seed":
        data["bodies"][0]["s0"].pop()
    elif failure == "body-mu":
        data["bodies"][0]["mu"] = 0.0
    elif failure == "body-wiring":
        data["bodies"][0]["series"].reverse()
    elif failure == "rotation-shape":
        data["rotation"].pop()
    elif failure == "rotation-nonfinite":
        data["rotation"][0] = float("nan")
    elif failure == "rotation-change":
        data["rotation"][0] = -data["rotation"][0]
    else:
        data["unexpected"] = "not silently ignored"
    source, output = tmp_path / "bad.json", tmp_path / "out.rs"
    source.write_text(json.dumps(data))
    output.write_text("preserve existing result")
    assert generator.main(["--input", str(source), "--out", str(output)]) == 2
    assert output.read_text() == "preserve existing result"
    assert sorted(path.name for path in tmp_path.iterdir()) == ["bad.json", "out.rs"]


@pytest.mark.parametrize("text", ["null", "[]", '"text"', '{"schema_version":1,"schema_version":1}', '{'])
def test_nonobjects_duplicate_fields_and_broken_json_are_clear_errors(text, tmp_path, capsys):
    source = tmp_path / "bad.json"
    source.write_text(text)
    assert generator.main(["--input", str(source), "--out", str(tmp_path / "out.rs")]) == 2
    assert "TASS17 generator:" in capsys.readouterr().err
    assert not (tmp_path / "out.rs").exists()


def test_check_mismatch_is_read_only_and_missing_output_is_an_input_error(tmp_path):
    output = tmp_path / "out.rs"
    output.write_text("old output")
    assert generator.main(["--out", str(output), "--check"]) == 1
    assert output.read_text() == "old output"
    assert generator.main(["--out", str(tmp_path / "missing.rs"), "--check"]) == 2
    assert not (tmp_path / "missing.rs").exists()


def test_input_aliases_and_hardlinks_refuse_without_destroying_snapshot(tmp_path):
    source = tmp_path / "input.json"
    shutil.copyfile(generator.DEFAULT_INPUT, source)
    original = source.read_bytes()
    alias = tmp_path / "alias"
    alias.hardlink_to(source)
    symlink = tmp_path / "symlink"
    symlink.symlink_to(source)
    for output in (source, alias, symlink):
        assert generator.main(["--input", str(source), "--out", str(output)]) == 2
        assert source.read_bytes() == original


def test_atomic_publication_failure_keeps_previous_output_and_cleans_own_temporary(tmp_path, monkeypatch):
    output = tmp_path / "out.rs"
    output.write_text("previous")
    def refuse(*args):
        raise OSError("controlled publication failure")
    monkeypatch.setattr(generator.os, "replace", refuse)
    assert generator.main(["--out", str(output)]) == 2
    assert output.read_text() == "previous"
    assert [path.name for path in tmp_path.iterdir()] == ["out.rs"]


def test_scalar_literal_keeps_signed_zero_and_ieee_extremes_without_rounding():
    for word in (-0.0, 0.0, float.fromhex("0x0.0000000000001p-1022"), float.fromhex("0x1.fffffffffffffp1023")):
        assert words([float(generator.f64_lit(word))]) == words([word])
