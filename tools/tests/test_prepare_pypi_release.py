"""Exercise release gates offline; publication tests use a fake local uv binary."""

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET
import zipfile

import pytest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("release_preparation", ROOT / "tools/release/prepare_pypi_release.py")
PREPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PREPARE)


@pytest.mark.parametrize("problem", [None, "sha", "running", "conclusion", "missing_job", "failed_job"])
def test_ci_gate_requires_the_exact_successful_native_matrix(problem):
    sha = "a" * 40
    record = {"headSha": sha, "status": "completed", "conclusion": "success",
              "jobs": [{"name": name, "conclusion": "success"} for name in PREPARE.CI_PLATFORMS]}
    if problem == "sha":
        record["headSha"] = "b" * 40
    if problem == "running":
        record["status"] = "in_progress"
    if problem == "conclusion":
        record["conclusion"] = "failure"
    if problem == "missing_job":
        record["jobs"].pop()
    if problem == "failed_job":
        record["jobs"][0]["conclusion"] = "failure"
    if problem is None:
        PREPARE.check_ci_run(record, sha)
    else:
        with pytest.raises(PREPARE.ReleaseError):
            PREPARE.check_ci_run(record, sha)


@pytest.mark.parametrize("kind", ["passed", "failure", "error", "skipped", "empty"])
def test_junit_gate_rejects_failed_empty_or_entirely_skipped_suites(tmp_path, kind):
    suite = ET.Element("testsuite")
    if kind != "empty":
        case = ET.SubElement(suite, "testcase", name="regression")
        if kind != "passed":
            ET.SubElement(case, kind, message="diagnostic")
    path = tmp_path / "report.xml"
    ET.ElementTree(suite).write(path, encoding="utf-8")
    if kind == "passed":
        assert PREPARE.junit_summary(path) == {"passed": 1, "failed": 0, "skipped": []}
    else:
        with pytest.raises(PREPARE.ReleaseError):
            PREPARE.junit_summary(path)


def test_junit_retains_optional_skip_names_and_reasons(tmp_path):
    path = tmp_path / "report.xml"
    path.write_text('<testsuite><testcase name="executed"/><testcase classname="optional" name="input">'
                    '<skipped message="public optional dataset absent"/></testcase></testsuite>')
    assert PREPARE.junit_summary(path) == {
        "passed": 1, "failed": 0,
        "skipped": [{"test": "optional::input", "reason": "public optional dataset absent"}],
    }


def git(root, *arguments):
    """Operate only on each test's temporary local Git repository."""
    return subprocess.check_output(
        ["git", "-c", "user.name=Release Test", "-c", "user.email=release@example.invalid",
         "-c", "commit.gpgsign=false", "-c", "tag.gpgsign=false", "-C", str(root), *arguments],
        text=True, stderr=subprocess.DEVNULL,
    ).strip()


@pytest.mark.parametrize("problem", [None, "changed_file", "extra_file", "head_changed", "missing_tag", "wrong_tag"])
def test_generated_publisher_checks_files_head_and_remote_tag_before_uv(tmp_path, monkeypatch, problem):
    root = tmp_path / "source"
    remote = tmp_path / "remote.git"
    root.mkdir()
    remote.mkdir()
    git(root, "init", "-q")
    git(remote, "init", "--bare", "-q")
    (root / "source.txt").write_text("owned dummy source")
    git(root, "add", "source.txt")
    git(root, "commit", "-qm", "initial")
    git(root, "remote", "add", "origin", str(remote))
    tag = "python-v9.0.0b1"
    if problem != "missing_tag":
        git(root, "tag", "-a", tag, "-m", "test release")
        git(root, "push", "-q", "origin", f"refs/tags/{tag}")
    if problem == "wrong_tag":
        (root / "source.txt").write_text("new owned dummy source")
        git(root, "commit", "-qam", "another commit")
    sha = git(root, "rev-parse", "HEAD")
    work = root / "target/release"
    upload = work / "upload"
    upload.mkdir(parents=True)
    for name in ("linux-x86.whl", "linux-arm.whl", "windows.whl", "macos.whl", "source.tar.gz"):
        (upload / name).write_bytes(b"synthetic payload, not a real distribution")
    record = {"sha": sha, "tag": tag, "version": "9.0.0b1",
              "files": {path.name: PREPARE.sha256(path) for path in upload.iterdir()}}
    (work / "release.json").write_text(json.dumps(record))
    monkeypatch.setattr(PREPARE, "ROOT", root)
    monkeypatch.setattr(PREPARE, "VENV", Path(sys.executable).parent)
    PREPARE.write_publisher(work, sha, tag)
    script = work / "publish.sh"
    assert subprocess.run(["bash", "-n", script]).returncode == 0

    # The fake command only records argv. It cannot reach PyPI or read a token.
    binary = tmp_path / "bin"
    binary.mkdir()
    capture = tmp_path / "uv-arguments.txt"
    uv = binary / "uv"
    uv.write_text('#!/bin/sh\nprintf "%s\\n" "$@" > "$PUBLISH_CAPTURE"\n')
    uv.chmod(0o700)
    environment = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ["PATH"],
                       UV_PUBLISH_TOKEN="synthetic-test-token", PUBLISH_CAPTURE=str(capture))
    if problem == "changed_file":
        (upload / "windows.whl").write_bytes(b"unverified replacement")
    if problem == "extra_file":
        (upload / "unverified.whl").write_bytes(b"extra")
    if problem == "head_changed":
        (root / "source.txt").write_text("modified owned source")
        git(root, "commit", "-qam", "after qualification")
    result = subprocess.run(["bash", script], env=environment, text=True, capture_output=True)
    assert "synthetic-test-token" not in result.stdout + result.stderr
    assert "SyntaxError" not in result.stderr
    if problem is None:
        assert result.returncode == 0, result.stderr
        arguments = capture.read_text().splitlines()
        assert arguments[0] == "publish"
        assert "https://upload.pypi.org/legacy/" in arguments
        assert len([argument for argument in arguments if argument.startswith("upload/")]) == 5
    else:
        assert result.returncode != 0
        assert not capture.exists(), "the publisher reached uv before qualification checks passed"


@pytest.mark.parametrize("problem", [
    None, "target", "platform", "source_report", "metadata_version",
    "junit_failure", "junit_all_skipped", "checksum", "report_version",
])
def test_downloaded_ci_wheels_are_checked_beyond_the_build_json(tmp_path, monkeypatch, problem):
    sha, version = "a" * 40, "9.0.0"
    for index, (platform, (target, tag)) in enumerate(PREPARE.TARGETS.items()):
        directory = tmp_path / ("wheels-" + platform)
        reports = directory / "reports"
        reports.mkdir(parents=True)
        (directory / "dist").mkdir()
        if index == 0 and problem == "platform":
            tag = "unqualified_platform"
        wheel = directory / "dist" / f"libephemeris-{version}-cp312-abi3-{tag}.whl"
        metadata_version = "9.0.1" if index == 0 and problem == "metadata_version" else version
        with zipfile.ZipFile(wheel, "w") as archive:
            archive.writestr(f"libephemeris-{version}.dist-info/METADATA",
                             f"Metadata-Version: 2.4\nName: libephemeris\nVersion: {metadata_version}\n")
        record = {"source_sha": sha, "version": version, "target": target,
                  "wheel": wheel.name, "sha256": PREPARE.sha256(wheel),
                  "tests": {python: {"returncode": 0} for python in PREPARE.PYTHONS}}
        if index == 0:
            if problem == "target":
                record["target"] = "wrong-target"
            if problem == "checksum":
                record["sha256"] = "0" * 64
            if problem == "report_version":
                record["version"] = "9.0.1"
        (reports / "build.json").write_text(json.dumps(record))
        source = "incorrect source" if index == 0 and problem == "source_report" else f"Commit: {sha}\nTarget: {target}\n"
        (reports / "source.txt").write_text(source)
        for python in PREPARE.PYTHONS:
            suite = ET.Element("testsuite")
            case = ET.SubElement(suite, "testcase", name="executed")
            if index == 0 and python == "3.14":
                if problem == "junit_failure":
                    ET.SubElement(case, "failure", message="must not qualify")
                if problem == "junit_all_skipped":
                    ET.SubElement(case, "skipped", message="nothing ran")
            ET.ElementTree(suite).write(reports / f"pytest-{python}.xml", encoding="utf-8")
    # Fixtures stand in for downloaded artifacts; never contact GitHub here.
    monkeypatch.setattr(PREPARE, "run", lambda *args, **kwargs: None)
    if problem is None:
        assert len(PREPARE.ci_wheels("owned-test-run", sha, version, tmp_path)) == 3
    else:
        with pytest.raises(PREPARE.ReleaseError):
            PREPARE.ci_wheels("owned-test-run", sha, version, tmp_path)
