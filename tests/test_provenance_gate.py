"""Tests for the clean-room physical-worktree provenance gate."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

from scripts import check_provenance


def _hits(root: Path) -> set[tuple[str, str]]:
    return {
        (path.relative_to(root).as_posix(), label)
        for path, label, _ in check_provenance._physical_foreign_artifacts(root)
    }


def _git(root: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", "-C", os.fspath(root), *args],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.stdout.strip()


def _init_git_repo(root: Path) -> None:
    root.mkdir(parents=True, exist_ok=True)
    _git(root, "init", "-q")
    _git(root, "config", "user.name", "Provenance Test")
    _git(root, "config", "user.email", "provenance@example.invalid")
    (root / "safe.txt").write_text("independent project data\n")
    _git(root, "add", "--", "safe.txt")
    _git(root, "commit", "-qm", "safe root")


def _history_hits(root: Path) -> set[tuple[str, str]]:
    return {
        (artifact.path, artifact.label)
        for artifact in check_provenance._reachable_history_foreign_artifacts(root)
    }


def test_history_gate_accepts_safe_complete_repository(tmp_path: Path) -> None:
    _init_git_repo(tmp_path)

    assert check_provenance._reachable_history_foreign_artifacts(tmp_path) == []


def test_history_gate_fails_closed_outside_git_repository(tmp_path: Path) -> None:
    assert _history_hits(tmp_path) == {("<repository>", "history-not-git")}


def test_history_gate_catches_committed_then_deleted_path(tmp_path: Path) -> None:
    _init_git_repo(tmp_path)
    forbidden = tmp_path / "libephemeris" / "seorbel.txt"
    forbidden.parent.mkdir()
    forbidden.write_text("payload must not be opened by the gate\n")
    _git(tmp_path, "add", "--", "libephemeris/seorbel.txt")
    _git(tmp_path, "commit", "-qm", "add forbidden path")
    tainted_commit = _git(tmp_path, "rev-parse", "HEAD")
    forbidden.unlink()
    _git(tmp_path, "add", "-u")
    _git(tmp_path, "commit", "-qm", "delete forbidden path")

    artifacts = check_provenance._reachable_history_foreign_artifacts(tmp_path)

    assert [(artifact.path, artifact.label) for artifact in artifacts] == [
        ("libephemeris/seorbel.txt", "foreign-data-file")
    ]
    assert artifacts[0].commits == (tainted_commit,)


def test_history_baseline_waiver_is_exact_and_ancestry_bounded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _init_git_repo(tmp_path)
    safe_baseline = _git(tmp_path, "rev-parse", "HEAD")
    forbidden = tmp_path / "libephemeris" / "seorbel.txt"
    forbidden.parent.mkdir()
    forbidden.write_text("payload must not be opened by the gate\n")
    _git(tmp_path, "add", "--", "libephemeris/seorbel.txt")
    _git(tmp_path, "commit", "-qm", "add forbidden path")
    tainted_commit = _git(tmp_path, "rev-parse", "HEAD")
    artifacts = check_provenance._reachable_history_foreign_artifacts(tmp_path)

    monkeypatch.setattr(
        check_provenance,
        "INHERITED_HISTORY_ALLOWLIST",
        frozenset({("libephemeris/seorbel.txt", "foreign-data-file")}),
    )
    monkeypatch.setattr(
        check_provenance, "INHERITED_HISTORY_BASELINE_COMMITS", (safe_baseline,)
    )
    inherited, actionable = check_provenance._partition_history_baseline(
        tmp_path, artifacts
    )
    assert inherited == []
    assert actionable == artifacts

    monkeypatch.setattr(
        check_provenance,
        "INHERITED_HISTORY_BASELINE_COMMITS",
        (tainted_commit,),
    )
    inherited, actionable = check_provenance._partition_history_baseline(
        tmp_path, artifacts
    )
    assert inherited == artifacts
    assert actionable == []


def test_history_gate_checks_tag_and_custom_ref_reachability(tmp_path: Path) -> None:
    _init_git_repo(tmp_path)
    main_branch = _git(tmp_path, "branch", "--show-current")
    _git(tmp_path, "switch", "-q", "--orphan", "detached-audit")
    forbidden = tmp_path / "vendor" / "swisseph.zip"
    forbidden.parent.mkdir()
    forbidden.write_text("opaque test payload\n")
    _git(tmp_path, "add", "--", "vendor/swisseph.zip")
    _git(tmp_path, "commit", "-qm", "unmerged forbidden path")
    tainted_commit = _git(tmp_path, "rev-parse", "HEAD")
    _git(tmp_path, "tag", "archived-audit", tainted_commit)
    _git(tmp_path, "update-ref", "refs/audit/legacy", tainted_commit)
    _git(tmp_path, "switch", "-q", main_branch)
    _git(tmp_path, "branch", "-D", "detached-audit")

    assert _history_hits(tmp_path) == {("vendor/swisseph.zip", "foreign-source-file")}
    _git(tmp_path, "tag", "-d", "archived-audit")
    assert _history_hits(tmp_path) == {("vendor/swisseph.zip", "foreign-source-file")}


def test_history_gate_ignores_orchestrator_checkpoint_refs(tmp_path: Path) -> None:
    _init_git_repo(tmp_path)
    main_branch = _git(tmp_path, "branch", "--show-current")
    _git(tmp_path, "switch", "-q", "--orphan", "checkpoint-audit")
    forbidden = tmp_path / "libephemeris" / "seorbel.txt"
    forbidden.parent.mkdir()
    forbidden.write_text("checkpoint-only payload\n")
    _git(tmp_path, "add", "--", "libephemeris/seorbel.txt")
    _git(tmp_path, "commit", "-qm", "checkpoint snapshot")
    checkpoint = _git(tmp_path, "rev-parse", "HEAD")
    _git(tmp_path, "update-ref", "refs/t3/checkpoints/test/turn/0", checkpoint)
    _git(tmp_path, "switch", "-q", main_branch)
    _git(tmp_path, "branch", "-D", "checkpoint-audit")

    assert check_provenance._reachable_history_foreign_artifacts(tmp_path) == []


def test_history_gate_rejects_shallow_repository(tmp_path: Path) -> None:
    source = tmp_path / "source"
    shallow = tmp_path / "shallow"
    _init_git_repo(source)
    (source / "second.txt").write_text("second commit\n")
    _git(source, "add", "--", "second.txt")
    _git(source, "commit", "-qm", "second")
    subprocess.run(
        ["git", "clone", "-q", "--depth", "1", source.as_uri(), os.fspath(shallow)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert _history_hits(shallow) == {("<repository>", "history-shallow-repository")}


def test_history_gate_disables_replace_objects(tmp_path: Path) -> None:
    _init_git_repo(tmp_path)
    forbidden = tmp_path / "seorbel.txt"
    forbidden.write_text("historical test payload\n")
    _git(tmp_path, "add", "--", "seorbel.txt")
    _git(tmp_path, "commit", "-qm", "tainted tree")
    tainted_commit = _git(tmp_path, "rev-parse", "HEAD")
    forbidden.unlink()
    _git(tmp_path, "add", "-u")
    _git(tmp_path, "commit", "-qm", "safe replacement tree")
    safe_commit = _git(tmp_path, "rev-parse", "HEAD")
    _git(tmp_path, "replace", tainted_commit, safe_commit)

    # Ordinary Git metadata traversal now hides the original tree.
    assert _git(tmp_path, "ls-tree", "-r", "--name-only", tainted_commit) == "safe.txt"
    assert _history_hits(tmp_path) == {("seorbel.txt", "foreign-data-file")}


def test_history_gate_handles_nul_delimited_odd_filename(tmp_path: Path) -> None:
    _init_git_repo(tmp_path)
    odd_rel = "odd\nfolder/SEFSTARS\tcatalog"
    odd = tmp_path / odd_rel
    odd.parent.mkdir()
    odd.write_text("historical test payload\n")
    _git(tmp_path, "add", "--", odd_rel)
    _git(tmp_path, "commit", "-qm", "odd filename")

    assert _history_hits(tmp_path) == {(odd_rel, "foreign-data-file")}
    assert "\n" not in check_provenance._display_history_path(odd_rel)


def test_physical_gate_catches_ignored_reference_directory(tmp_path: Path) -> None:
    reference_dir = tmp_path / "data" / "reference"
    reference_dir.mkdir(parents=True)

    assert _hits(tmp_path) == {("data/reference", "foreign-data-dir")}


def test_physical_gate_catches_mixed_case_reference_directory(tmp_path: Path) -> None:
    reference_dir = tmp_path / "Data" / "REFERENCE"
    reference_dir.mkdir(parents=True)

    assert _hits(tmp_path) == {("Data/REFERENCE", "foreign-data-dir")}


def test_physical_gate_catches_reference_files_by_name(tmp_path: Path) -> None:
    artifact_dir = tmp_path / "ignored" / "nested"
    artifact_dir.mkdir(parents=True)
    for name in (
        "planet.se1",
        "sefstars.json",
        "sweph.c",
        "sweprivate.h",
        "pyswisseph.py",
    ):
        (artifact_dir / name).touch()

    assert _hits(tmp_path) == {
        ("ignored/nested/planet.se1", "foreign-data-file"),
        ("ignored/nested/sefstars.json", "foreign-data-file"),
        ("ignored/nested/sweph.c", "foreign-source-file"),
        ("ignored/nested/sweprivate.h", "foreign-source-file"),
        ("ignored/nested/pyswisseph.py", "foreign-source-file"),
    }


def test_physical_gate_catches_reference_source_directory(tmp_path: Path) -> None:
    source_dir = tmp_path / "vendor" / "Swiss-Ephemeris-source"
    source_dir.mkdir(parents=True)

    assert _hits(tmp_path) == {("vendor/Swiss-Ephemeris-source", "foreign-source-dir")}


def test_physical_gate_catches_foreign_named_directories(tmp_path: Path) -> None:
    for dirname in ("sefstars", "seorbel", "SEDELTAT"):
        directory = tmp_path / dirname
        directory.mkdir()

    assert _hits(tmp_path) == {
        ("SEDELTAT", "foreign-data-dir"),
        ("sefstars", "foreign-data-dir"),
        ("seorbel", "foreign-data-dir"),
    }


def test_physical_gate_catches_opaque_archives_and_compiled_sources(
    tmp_path: Path,
) -> None:
    for name in (
        "pyswisseph-2.10.3.2.tar.gz",
        "swisseph.zip",
        "swisseph.cpython-312.so",
        "libswe.a",
        "swedll64.dll",
    ):
        (tmp_path / name).touch()

    assert len(_hits(tmp_path)) == 5


def test_archive_path_classifier_rejects_unsafe_and_disguised_names() -> None:
    classify = check_provenance._classify_artifact_name

    assert classify("pkg/data/./reference/planet.bin") == "foreign-data-dir"
    assert classify("pkg/data/x/../reference/planet.bin") == "unsafe-path"
    assert classify("/absolute/pkg/module.py") == "unsafe-path"
    assert classify("pkg/vendor/PySwissEph-2.10/", is_directory=True)
    assert classify("pkg/vendor/swisseph.zip") == "foreign-source-file"


def test_physical_gate_checks_reference_names_inside_external_caches(
    tmp_path: Path,
) -> None:
    for dirname in (".git", ".venv", ".pytest_cache", ".ruff_cache"):
        directory = tmp_path / dirname
        directory.mkdir()
        (directory / "planet.se1").touch()
        if dirname != ".git":
            (directory / "independent-extension.so").touch()

    assert _hits(tmp_path) == {
        (".git/planet.se1", "foreign-data-file"),
        (".pytest_cache/planet.se1", "foreign-data-file"),
        (".ruff_cache/planet.se1", "foreign-data-file"),
        (".venv/planet.se1", "foreign-data-file"),
    }


def test_physical_gate_allows_generic_dependency_reference_directory(
    tmp_path: Path,
) -> None:
    directory = tmp_path / ".venv" / "site-packages" / "package" / "data" / "reference"
    directory.mkdir(parents=True)
    (directory / "independent-fixture.txt").write_text("synthetic fixture\n")

    assert _hits(tmp_path) == set()


def test_physical_gate_traverses_validation_git_names(tmp_path: Path) -> None:
    hidden = tmp_path / "validation" / ".git"
    hidden.mkdir(parents=True)
    (hidden / "sefstars.txt").touch()

    assert _hits(tmp_path) == {("validation/.git/sefstars.txt", "foreign-data-file")}


def test_physical_gate_does_not_prune_nested_git_lookalike(tmp_path: Path) -> None:
    hidden = tmp_path / "pkg" / ".git"
    hidden.mkdir(parents=True)
    (hidden / "sefstars.txt").touch()

    assert _hits(tmp_path) == {("pkg/.git/sefstars.txt", "foreign-data-file")}


def test_physical_gate_does_not_prune_nested_lookalike_environment(
    tmp_path: Path,
) -> None:
    hidden = tmp_path / "libephemeris" / "data" / ".venv"
    hidden.mkdir(parents=True)
    (hidden / "swisseph.zip").touch()

    assert _hits(tmp_path) == {
        ("libephemeris/data/.venv/swisseph.zip", "foreign-source-file")
    }


def test_physical_gate_checks_symlink_target_names_without_following(
    tmp_path: Path,
) -> None:
    target = tmp_path / "target.txt"
    target.write_text("project text")
    link = tmp_path / "linked.txt"
    link.symlink_to(target.name)
    foreign = tmp_path / "foreign-link"
    foreign.symlink_to("outside/swisseph.zip")

    assert _hits(tmp_path) == {("foreign-link", "foreign-symlink-target")}


def test_generated_looking_pyc_and_pycache_artifacts_are_rejected(
    tmp_path: Path,
) -> None:
    cache = tmp_path / "pkg" / "__pycache__"
    cache.mkdir(parents=True)
    (cache / "module.cpython-312.pyc").touch()
    (cache / "module.cpython-310-pytest-8.4.2.pyc").touch()
    (cache / "swisseph.zip").touch()

    assert _hits(tmp_path) == {
        ("pkg/__pycache__/module.cpython-310-pytest-8.4.2.pyc", "compiled-binary"),
        ("pkg/__pycache__/module.cpython-312.pyc", "compiled-binary"),
        ("pkg/__pycache__/swisseph.zip", "foreign-source-file"),
    }


def test_non_generated_pyc_name_is_not_ignored(tmp_path: Path) -> None:
    cache = tmp_path / "pkg" / "__pycache__"
    cache.mkdir(parents=True)
    (cache / "payload.pyc").touch()

    assert _hits(tmp_path) == {("pkg/__pycache__/payload.pyc", "compiled-binary")}


def test_nested_pycache_component_does_not_waive_compiled_file(tmp_path: Path) -> None:
    nested = tmp_path / "pkg" / "__pycache__" / "nested"
    nested.mkdir(parents=True)
    (nested / "module.cpython-312.pyc").touch()

    assert _hits(tmp_path) == {
        ("pkg/__pycache__/nested/module.cpython-312.pyc", "compiled-binary")
    }


def test_physical_gate_rejects_fifo_without_opening_it(tmp_path: Path) -> None:
    fifo = tmp_path / "innocent"
    os.mkfifo(fifo)

    assert _hits(tmp_path) == {("innocent", "special-file")}
    assert check_provenance._content_scan_files(tmp_path) == []


def test_content_scan_prunes_exact_git_payloads(tmp_path: Path) -> None:
    for directory in (tmp_path / ".git", tmp_path / "validation" / ".git"):
        directory.mkdir(parents=True)
        (directory / "opaque-object").write_bytes(b"\xff")
    owned = tmp_path / "owned.txt"
    owned.write_text("owned")

    assert check_provenance._content_scan_files(tmp_path) == [owned]


def test_content_scan_covers_text_regardless_of_suffix_or_validation_name(
    tmp_path: Path,
) -> None:
    paths = [
        tmp_path / "audit.sh",
        tmp_path / "uv.lock",
        tmp_path / "Makefile",
        tmp_path / "page.HTML",
        tmp_path / "config.xml",
        tmp_path / "validation" / "audit.py",
    ]
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("plain UTF-8 project text")
    binary = tmp_path / "image.bin"
    binary.write_bytes(b"\x00\x01\x02")
    partial_kernel = tmp_path / "de440.bsp.download4"
    partial_kernel.write_bytes(b"\xff\x00binary kernel fragment")
    nul_text = tmp_path / "nul.txt"
    nul_text.write_bytes(b"\x00plain UTF-8 project text")
    latin1 = tmp_path / "latin1.txt"
    latin1.write_bytes(b"caf\xe9 PyMeeus")
    cached_text = tmp_path / "pkg" / "__pycache__" / "notes.txt"
    cached_text.parent.mkdir(parents=True)
    cached_text.write_text("plain UTF-8 project text")

    scanned = {
        path.relative_to(tmp_path).as_posix()
        for path in check_provenance._content_scan_files(tmp_path)
    }

    assert scanned == {
        "Makefile",
        "audit.sh",
        "config.xml",
        "latin1.txt",
        "nul.txt",
        "page.HTML",
        "pkg/__pycache__/notes.txt",
        "uv.lock",
        "validation/audit.py",
    }


def test_non_utf8_project_text_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _init_git_repo(tmp_path)
    (tmp_path / "latin1.txt").write_bytes(b"caf\xe9 PyMeeus")
    monkeypatch.setattr(check_provenance, "REPO_ROOT", tmp_path)

    assert check_provenance.main() == 1


def test_history_findings_do_not_hide_worktree_content_findings(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _init_git_repo(tmp_path)
    forbidden = tmp_path / "libephemeris" / "seorbel.txt"
    forbidden.parent.mkdir()
    forbidden.write_text("historical payload must not be opened\n")
    _git(tmp_path, "add", "--", "libephemeris/seorbel.txt")
    _git(tmp_path, "commit", "-qm", "add forbidden path")
    forbidden.unlink()
    _git(tmp_path, "add", "-u")
    _git(tmp_path, "commit", "-qm", "delete forbidden path")
    (tmp_path / "current.txt").write_text("JupiterMoons\n")
    monkeypatch.setattr(check_provenance, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(
        check_provenance,
        "INHERITED_HISTORY_BASELINE_COMMITS",
        (_git(tmp_path, "rev-parse", "HEAD"),),
    )
    monkeypatch.setattr(check_provenance, "INHERITED_HISTORY_ALLOWLIST", frozenset())

    assert check_provenance.main() == 1
    output = capsys.readouterr().out
    assert "history:libephemeris/seorbel.txt" in output
    assert "current.txt:1: [pymeeus]" in output
    assert "1 new reachable-history hit(s)" in output
    assert "1 content hit(s)" in output


def test_inline_waivers_are_line_scoped_and_exact_path_scoped() -> None:
    token = "PyMeeus"

    assert check_provenance._line_hit_labels("docs/history.md", token) == ["pymeeus"]
    assert (
        check_provenance._line_hit_labels(
            "CHANGELOG.md", f"{token} <!-- provenance-implementation-ok -->"
        )
        == []
    )
    assert check_provenance._line_hit_labels(
        "scripts/payload.py",
        f"{token} # provenance-implementation-ok",
    ) == ["pymeeus"]
    assert check_provenance._line_hit_labels(
        "libephemeris/runtime.py",
        f"{token}  # provenance-implementation-ok",
    ) == ["pymeeus"]
    assert check_provenance._line_hit_labels(
        "LIBEPHEMERIS./runtime.py",
        f"{token}  # provenance-implementation-ok",
    ) == ["pymeeus"]
    assert check_provenance._line_hit_labels(
        "docs/history.md", "GPL-3.0 # provenance-implementation-ok"
    ) == ["copyleft"]
    assert (
        check_provenance._line_hit_labels(
            "docs/guides/getting-started.md", "GPL-3.0 # provenance-ok"
        )
        == []
    )
    assert check_provenance._line_hit_labels(
        "scripts/payload.py", "GPL-3.0 # provenance-ok"
    ) == ["copyleft"]


def test_main_fails_before_reading_a_named_foreign_artifact(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _init_git_repo(tmp_path)
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "sefstars.txt").touch()
    monkeypatch.setattr(check_provenance, "REPO_ROOT", tmp_path)

    def _unexpected_read(*args, **kwargs):
        pytest.fail("the name-only gate must return before opening any payload")

    monkeypatch.setattr(Path, "read_text", _unexpected_read)

    assert check_provenance.main() == 1


def test_main_scans_blocked_content_in_shell_files(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _init_git_repo(tmp_path)
    (tmp_path / "audit.sh").write_text("PyMeeus\n")
    monkeypatch.setattr(check_provenance, "REPO_ROOT", tmp_path)

    assert check_provenance.main() == 1
