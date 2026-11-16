"""The performance gate must reject changes and incomplete captures."""

import copy
import gzip
import json
from pathlib import Path
import struct
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bench"))
from bitdiff import compare, encode, outcome  # noqa: E402


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
