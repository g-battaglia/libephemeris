"""The performance gate must reject changes and incomplete captures."""

import copy
import ctypes
import gzip
import json
import mmap
import os
from pathlib import Path
import shutil
import struct
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bench"))
import benchmark  # noqa: E402
from benchmark import linux_resident_bytes, verify_cold_os_pages  # noqa: E402
from bitdiff import compare, encode, outcome  # noqa: E402
from common import (  # noqa: E402
    ROOT, PrerequisiteError, digest, prepare, rust_source_revisions, shared_medium_identity,
)
import wheel_check  # noqa: E402


def capture(path, result, *, truncate=False, api="calc_ut", args=(2451545., 4, 258)):
    build = dict.fromkeys(("platform", "release", "machine", "libc", "python", "cpu"), "test")
    rows = [
        {"schema": 1, "spec": {"expected_count": 1, "matrix": [["leb", "base"]]}, "build": build},
        ["leb", "base", 0, api, encode(args), result],
        {"count": 1, "coverage": {f"leb/base/{api}": 1}},
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
    identical = compare(left, left)
    assert identical["status"] == "identical"
    assert identical["difference_counts"] == {}
    different = compare(left, right)
    assert different["differences"] == 1
    assert different["warning_calls"] == {"baseline": 0, "candidate": int(change == "warning")}
    if change == "warning":
        transitions = {}
    elif change == "error":
        transitions = {"builtins.ValueError -> builtins.TypeError": 1}
    else:
        transitions = {"value -> value": 1}
    assert different["difference_counts"] == {
        "leb/base/calc_ut": {
            "requests": 1,
            "result_transitions": transitions,
            "changed_warnings": int(change == "warning"),
        },
    }


def test_difference_counts_distinguish_value_from_range_error(tmp_path):
    left, right = tmp_path / "left.gz", tmp_path / "right.gz"
    error = "libephemeris.EphemerisRangeError"
    capture(left, {"value": None, "warnings": []})
    capture(right, {"error": error, "warnings": [["builtins.UserWarning", "changed"]]})
    forward = compare(left, right)["difference_counts"]["leb/base/calc_ut"]
    backward = compare(right, left)["difference_counts"]["leb/base/calc_ut"]
    assert forward == {"requests": 1, "result_transitions": {f"value -> {error}": 1}, "changed_warnings": 1}
    assert backward == {"requests": 1, "result_transitions": {f"{error} -> value": 1}, "changed_warnings": 1}


def test_incomplete_capture(tmp_path):
    left, right = tmp_path / "baseline.gz", tmp_path / "candidate.gz"
    capture(left, {"value": None, "warnings": []})
    capture(right, {"value": None, "warnings": []}, truncate=True)
    with pytest.raises((ValueError, StopIteration)):
        compare(left, right)


@pytest.mark.skipif(sys.platform != "linux", reason="Linux shared snapshot mapping")
def test_shared_medium_mapping_reports_kernel_inode():
    fd = os.memfd_create("libephemeris-medium", os.MFD_CLOEXEC)
    try:
        os.write(fd, b"test image")
        with mmap.mmap(fd, 0, access=mmap.ACCESS_READ):
            device, inode = shared_medium_identity()
            assert inode == os.fstat(fd).st_ino
            assert device
    finally:
        os.close(fd)
    with pytest.raises(PrerequisiteError, match="did not map"):
        shared_medium_identity()


def test_position_inventory_uses_captured_requests_including_errors(tmp_path):
    path = tmp_path / "capture.gz"
    capture(path, {"error": "libephemeris.EphemerisRangeError", "warnings": []},
            api="context.calc_pctr", args=(2451545., 9599, 14, 16658))
    assert compare(path, path)["position_inventory"] == {
        "leb/base/context.calc_pctr": {
            "body_ids_by_class": {"unresolved_planetary_moon_id": [9599]},
            "center_ids": [14],
            "requested_flags": [16658],
        },
    }


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


def test_wheel_check_reports_missing_input(tmp_path, monkeypatch):
    dist = tmp_path / "dist"
    dist.mkdir()
    (dist / "candidate.whl").touch()
    output = tmp_path / "reports"
    monkeypatch.setattr(wheel_check, "ROOT", tmp_path)
    monkeypatch.setattr(sys, "argv", ["wheel_check.py", "--output", str(output)])
    monkeypatch.setattr(wheel_check.venv, "create", lambda *_args, **_kwargs: None)

    def pip(args, *, cwd, check):
        assert check and args[1:3] == ["-m", "pip"]
        if args[3] == "download":
            (cwd / "published.whl").touch()

    monkeypatch.setattr(wheel_check.subprocess, "run", pip)
    assert wheel_check.main() == 2
    report = json.loads((output / "comparison.json").read_text())
    assert report["status"] == "missing_prerequisite"
    assert "base_core.leb2" in report["error"]


@pytest.mark.parametrize("workload, flags, actual_flags", (
    ("seq_position", 266, 266), ("time_series", 258, 258), ("topo", 258, 33026),
))
def test_position_preflight_probes_only_the_measured_route(tmp_path, monkeypatch, capsys,
                                                         workload, flags, actual_flags):
    calls = []
    eph = SimpleNamespace(calc_ut=lambda *args: calls.append(args),
                          sidtime=None, get_ayanamsa_ut=None,
                          start_tracing=lambda: SimpleNamespace(reset=lambda: None),
                          get_trace_results=lambda: {4: "Skyfield"})
    monkeypatch.setitem(sys.modules, "libephemeris", eph)
    monkeypatch.setattr(benchmark, "configure", lambda *_, **__: None)
    spec = tmp_path / "job.json"
    spec.write_text(json.dumps(dict(mode="skyfield", tier="base", workload=workload,
                                    dates=[2451545.], flags=flags)))
    benchmark.worker(spec, probe=True)
    assert calls == [(2451545., 4, actual_flags)]
    assert json.loads(capsys.readouterr().out)["routes"] == {"4": "Skyfield"}


def test_rust_source_paths_across_binary_read_boundary(tmp_path):
    binary = tmp_path / "extension"
    revision = b"58cf8f9c35b355199727627d07ae34a4e2846d67"
    binary.write_bytes(b"\x00" * (1024 * 1024 - 12) + b"/rustc/" + revision + b"/library"
                       + b"\x00\\rustc\\" + revision + b"\\library")
    assert rust_source_revisions(binary) == [revision.decode("ascii")]


def test_first_use_samples_get_fresh_staged_files(tmp_path, monkeypatch):
    name = "leb/core.leb2"
    data = tmp_path / "data"
    source = data / name
    source.parent.mkdir(parents=True)
    source.write_bytes(b"same pinned input")
    env = {"LIBEPHEMERIS_DATA_DIR": str(data), "LIBEPHEMERIS_LEB_READER": "indexed"}
    staged = []

    def run_worker(_args, *, env, cwd, text):
        copy = Path(env["LIBEPHEMERIS_DATA_DIR"]) / name
        staged.append(copy)
        assert copy.read_bytes() == source.read_bytes()
        assert not copy.samefile(source)
        assert env["LIBEPHEMERIS_LEB_READER"] == "indexed"
        assert cwd == tmp_path and text
        return '{"first_call_ns": 10}'

    monkeypatch.setattr(benchmark.subprocess, "check_output", run_worker)
    for _ in range(2):
        assert benchmark.timed_sample([], env, tmp_path, [name], "first_use") == {"first_call_ns": 10}
    assert len(set(staged)) == 2 and all(not path.exists() for path in staged)
    assert env["LIBEPHEMERIS_DATA_DIR"] == str(data)


def test_cold_start_refuses_resident_input_pages(tmp_path, monkeypatch):
    (tmp_path / "input").write_bytes(b"sample")
    monkeypatch.setattr(os, "posix_fadvise", lambda *_: None, raising=False)
    monkeypatch.setattr(os, "POSIX_FADV_DONTNEED", 4, raising=False)

    def mincore(_address, _length, residency):
        residency[0] = 1
        return 0

    monkeypatch.setattr(ctypes, "CDLL", lambda *_args, **_kwargs: SimpleNamespace(mincore=mincore))
    with pytest.raises(PrerequisiteError, match="resident bytes after eviction"):
        verify_cold_os_pages(tmp_path, ["input"])


def test_linux_residency_counts_only_low_bits_and_partial_final_page(tmp_path):
    length = mmap.PAGESIZE + 1
    path = tmp_path / "input"
    path.write_bytes(b"\x00" * length)

    def mincore(address, observed_length, residency):
        assert address % mmap.PAGESIZE == 0
        assert observed_length == length
        assert len(residency) == 2
        residency[:] = (254, 1)
        return 0

    with path.open("rb") as stream:
        assert linux_resident_bytes(stream.fileno(), mincore) == mmap.PAGESIZE
