"""The performance gate must reject changes and incomplete captures."""

import copy
import gzip
import json
from pathlib import Path
import shutil
import struct
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bench"))
from bitdiff import compare, encode, outcome  # noqa: E402
from common import ROOT, digest, prepare, rust_source_revisions  # noqa: E402


def capture(path, result, *, truncate=False):
    build = dict.fromkeys(("platform", "release", "machine", "libc", "python", "cpu"), "test")
    rows = [
        {"schema": 1, "spec": {"expected_count": 1, "matrix": [["leb", "base"]]}, "build": build},
        ["leb", "base", 0, "calc_ut", [0], result],
        {"count": 1, "coverage": {"leb/base/calc_ut": 1}},
    ]
    with gzip.open(path, "wt", encoding="utf-8") as file:
        for row in rows[:-1] if truncate else rows:
            print(json.dumps(row), file=file)


@pytest.mark.parametrize("change", ("bit", "zero", "nan", "flag", "warning", "error", "layout"))
def test_exact_comparison(tmp_path, change):
    nan, = struct.unpack(">d", bytes.fromhex("7ff8000000000003"))
    original = {"value": encode(((1., -0., nan), 258)), "warnings": []}
    changed = copy.deepcopy(original)
    values = changed["value"]["tuple"][0]["tuple"]
    if change == "bit":
        values[0]["f64"] = "3ff0000000000001"
    elif change == "zero":
        values[1]["f64"] = "0000000000000000"
    elif change == "nan":
        values[2]["f64"] = "7ff8000000000004"
    elif change == "flag":
        changed["value"]["tuple"][1] = 2
    elif change == "warning":
        changed["warnings"] = [["builtins.UserWarning", "changed"]]
    elif change == "error":
        original = outcome(float, ("invalid",))
        changed = outcome(float, (None,))
    else:
        changed["value"] = {"list": changed["value"]["tuple"]}

    left, right = tmp_path / "baseline.gz", tmp_path / "candidate.gz"
    capture(left, original)
    capture(right, changed)
    assert compare(left, left)["status"] == "identical"
    assert compare(left, right)["differences"] == 1


def test_incomplete_capture(tmp_path):
    left, right = tmp_path / "baseline.gz", tmp_path / "candidate.gz"
    capture(left, {"value": None, "warnings": []})
    capture(right, {"value": None, "warnings": []}, truncate=True)
    with pytest.raises((ValueError, StopIteration)):
        compare(left, right)


def test_staging_keeps_each_run_on_an_independent_inode(tmp_path):
    source = tmp_path / "source"
    original = source / "leb/base_core.leb2"
    original.parent.mkdir(parents=True)
    shutil.copyfile(ROOT / "core/data/leb2/base_core.leb2", original)
    before = original.stat()
    inventory = prepare(source, tmp_path / "one", ("leb",), ("base",))
    prepare(source, tmp_path / "two", ("leb",), ("base",))
    first = tmp_path / "one/leb/base_core.leb2"
    second = tmp_path / "two/leb/base_core.leb2"
    assert not first.samefile(original)
    assert not first.samefile(second)
    assert original.stat().st_ctime_ns == before.st_ctime_ns
    assert digest(first) == digest(second) == inventory["leb/base_core.leb2"]["sha256"]


def test_rust_source_paths_across_binary_read_boundary(tmp_path):
    binary = tmp_path / "extension"
    revision = b"58cf8f9c35b355199727627d07ae34a4e2846d67"
    binary.write_bytes(b"\x00" * (1024 * 1024 - 12) + b"/rustc/" + revision + b"/library"
                       + b"\x00\\rustc\\" + revision + b"\\library")
    assert rust_source_revisions(binary) == [revision.decode("ascii")]
