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


# ---------------------------------------------------------------------------
# Output-fitting narrative gate (6th category)
# ---------------------------------------------------------------------------


def _narrative_hits(root: Path) -> list[tuple[str, str]]:
    hits, _ = check_provenance._reference_narrative_hits(root)
    return [
        (path.relative_to(root).as_posix(), label)
        for path, label, _ in ((p, lb, s) for p, _ln, lb, s in hits)
    ]


def _write_pkg_module(root: Path, name: str, source: str) -> Path:
    pkg = root / "libephemeris"
    pkg.mkdir(parents=True, exist_ok=True)
    path = pkg / name
    path.write_text(source)
    return path


class TestNarrativeGateCatchesRealArchetypes:
    """One fixture per real removed-case archetype the gate must catch."""

    def test_sampling_domain_narrative_is_caught(self, tmp_path: Path) -> None:
        # Colure-fixup archetype: conditional branches justified by a grid
        # sweep of the reference, no unusual literal anywhere.
        _write_pkg_module(
            tmp_path,
            "mod_a.py",
            '''\
def _fixup(armc):
    """Behavioral comparison with the reference API (armc grid 0-355,
    lat -89..89) shows isolated-point outputs at those exact hits."""
    if armc > 180.0:
        return armc - 0.5
    return armc
''',
        )
        assert ("libephemeris/mod_a.py", "reference-derived-evidence") in (
            _narrative_hits(tmp_path)
        )

    def test_agreement_to_tolerance_is_caught(self, tmp_path: Path) -> None:
        # Selector-rule archetype: a rule derived by observing which inputs
        # behaved differently, justified by an agreement claim.
        _write_pkg_module(
            tmp_path,
            "mod_b.py",
            '''\
_PINNED = frozenset({"O", "W", "N"})


def _rate(h):
    """The rule below reproduces the reference output to 0 deg/day
    across latitude 0-78, four epochs."""
    return 0.0 if h in _PINNED else 1.0
''',
        )
        assert ("libephemeris/mod_b.py", "reference-derived-evidence") in (
            _narrative_hits(tmp_path)
        )

    def test_reference_located_constant_is_caught(self, tmp_path: Path) -> None:
        # Eclipse-margin archetype: a bare comment locating an observed
        # transition; a ROUND number, so the numeric-literal gate is blind.
        path = _write_pkg_module(
            tmp_path,
            "mod_c.py",
            "# Measured reference behavior: the transition sits exactly at\n"
            "# max +- 1e-4 day (8.64 s).\n"
            "_MARGIN = 1e-4\n",
        )
        hits, _ack = check_provenance._reference_narrative_hits(tmp_path)
        labels = {label for _p, _ln, label, _s in hits}
        assert "reference-authored-constant" in labels
        # The explicit regression: the high-precision-literal gate must NOT
        # see this round number.
        assert check_provenance._significant_digits("1e-4") < (
            check_provenance.HIGH_PRECISION_SIG_DIGITS
        )
        assert path.exists()

    def test_explicit_fitting_verb_is_caught(self, tmp_path: Path) -> None:
        _write_pkg_module(
            tmp_path,
            "mod_d.py",
            "# Threshold chosen so that our transition instant matches the\n"
            "# reference.\n"
            "_T = 0.5\n",
        )
        labels = {lb for _p, lb in _narrative_hits(tmp_path)}
        assert "reference-derived-evidence" in labels

    def test_round_constant_fires_clause_b_not_the_literal_gate(
        self, tmp_path: Path
    ) -> None:
        _write_pkg_module(
            tmp_path,
            "mod_e.py",
            "# The reference implementation switches branches here.\n_GATE = 0.5\n",
        )
        hits, _ack = check_provenance._reference_narrative_hits(tmp_path)
        assert [(lb) for _p, _ln, lb, _s in hits] == ["reference-authored-constant"]
        literal_hits = check_provenance._high_precision_literal_hits(tmp_path)
        assert literal_hits == []


class TestNarrativeGateLeavesContractsAlone:
    """One fixture per legitimate archetype the gate must NOT touch.

    The annotated statement is an INT binding: contract comments and
    divergence notes are Clause A material, and Clause B deliberately
    stays strict for module-level FLOAT constants that name the
    reference (that is the _PHENO_EARTH-style case the dual-key waiver
    exists for).
    """

    @pytest.mark.parametrize(
        "comment",
        [
            # Retflag-echo contract.
            "# Compatibility contract: the retflag echoes FLG_MOSEPH only\n"
            "# when it is the sole source bit.\n",
            # Error-type contract.
            "# The reference raises TypeError for a multi-character\n# selector.\n",
            # Return-shape contract.
            "# An unknown selector yields an empty tuple rather than the\n"
            "# 12-cusp fallback (the reference rule).\n",
            # Documented divergence magnitude — policy-encouraged.
            "# This leaves the reported speed ~0.2-0.8 arcsec/day off the\n"
            "# reference near the poles.\n",
            # Negated attestation.
            "# No zero point is fitted to compatibility output; the value\n"
            "# is never fitted to external output.\n",
            # Fitting against a PUBLISHED source.
            "# Delaunay series fitted at the DE440 perigee passages over\n"
            "# the medium-kernel interval.\n",
            # Open-source baseline comparison.
            "# Tracks the Skyfield reference to <0.005 arcsec across the\n"
            "# full catalogue.\n",
        ],
        ids=[
            "retflag-echo",
            "error-type",
            "return-shape",
            "divergence-magnitude",
            "negated-attestation",
            "published-source-fit",
            "skyfield-baseline",
        ],
    )
    def test_legitimate_comment_is_not_a_hit(
        self, tmp_path: Path, comment: str
    ) -> None:
        _write_pkg_module(tmp_path, "mod_ok.py", comment + "_X = 1\n")
        assert _narrative_hits(tmp_path) == []


class TestNarrativeGateWaiverMechanism:
    """The dual-key contract-attestation waiver."""

    SOURCE = (
        "# Compatibility contract sentinel: the reference returns this\n"
        "# fixed tuple for the degenerate body.  provenance-contract-ok\n"
        "_SENTINEL = 180.0\n"
    )
    RATIONALE = (
        "degenerate-body sentinel tuple; a shape/echo contract in the same "
        "class as a retflag echo, re-derivable from the public API surface "
        "without reading any numeric output"
    )

    def test_inline_token_alone_does_not_waive(self, tmp_path: Path) -> None:
        _write_pkg_module(tmp_path, "mod_w.py", self.SOURCE)
        hits, ack = check_provenance._reference_narrative_hits(tmp_path)
        assert len(hits) == 1 and ack == []

    def test_registry_alone_does_not_waive(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        source_without_token = self.SOURCE.replace("  provenance-contract-ok", "")
        _write_pkg_module(tmp_path, "mod_w.py", source_without_token)
        monkeypatch.setattr(
            check_provenance,
            "CONTRACT_ATTESTATION_ALLOWLIST",
            frozenset({("libephemeris/mod_w.py", "_SENTINEL", self.RATIONALE)}),
        )
        hits, ack = check_provenance._reference_narrative_hits(tmp_path)
        assert len(hits) == 1 and ack == []

    def test_both_keys_waive_and_are_reported(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        _write_pkg_module(tmp_path, "mod_w.py", self.SOURCE)
        monkeypatch.setattr(
            check_provenance,
            "CONTRACT_ATTESTATION_ALLOWLIST",
            frozenset({("libephemeris/mod_w.py", "_SENTINEL", self.RATIONALE)}),
        )
        hits, ack = check_provenance._reference_narrative_hits(tmp_path)
        assert hits == []
        assert len(ack) == 1
        assert ack[0][2] == "acknowledged-contract-attestation"
        assert "_SENTINEL" in ack[0][3]

    def test_waiver_dies_on_symbol_rename(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        renamed = self.SOURCE.replace("_SENTINEL", "_RENAMED")
        _write_pkg_module(tmp_path, "mod_w.py", renamed)
        monkeypatch.setattr(
            check_provenance,
            "CONTRACT_ATTESTATION_ALLOWLIST",
            frozenset({("libephemeris/mod_w.py", "_SENTINEL", self.RATIONALE)}),
        )
        hits, ack = check_provenance._reference_narrative_hits(tmp_path)
        assert len(hits) == 1 and ack == []

    def test_short_or_vague_rationale_does_not_waive(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        _write_pkg_module(tmp_path, "mod_w.py", self.SOURCE)
        for bad in ("too short", "citation needed " * 10):
            monkeypatch.setattr(
                check_provenance,
                "CONTRACT_ATTESTATION_ALLOWLIST",
                frozenset({("libephemeris/mod_w.py", "_SENTINEL", bad)}),
            )
            hits, ack = check_provenance._reference_narrative_hits(tmp_path)
            assert len(hits) == 1 and ack == []

    def test_generated_module_header_exempts_clause_b(self, tmp_path: Path) -> None:
        _write_pkg_module(
            tmp_path,
            "mod_gen.py",
            "# AUTO-GENERATED by scripts/generate_tables.py\n"
            "# The reference implementation switches branches here.\n"
            "_GATE = 0.5\n",
        )
        hits, _ack = check_provenance._reference_narrative_hits(tmp_path)
        assert all(lb != "reference-authored-constant" for _p, _ln, lb, _s in hits)


class TestNarrativeGateReviewRegressions:
    """Regression fixtures from the adversarial review of the gate."""

    def test_coincide_verb_is_an_agreement_verb(self, tmp_path: Path) -> None:
        _write_pkg_module(
            tmp_path,
            "mod_r1.py",
            "# Hand-adjusted until the outputs coincide with the reference\n"
            "# within 0.5 arcsec across the tested dates.\n"
            "_K = 0.7\n",
        )
        labels = {lb for _p, lb in _narrative_hits(tmp_path)}
        assert "reference-derived-evidence" in labels

    def test_long_parenthetical_negation_is_not_a_hit(self, tmp_path: Path) -> None:
        _write_pkg_module(
            tmp_path,
            "mod_r2.py",
            "# This coefficient was not, contrary to what an earlier revision\n"
            "# of this comment seemed to suggest to some reviewers, tuned to\n"
            "# the reference in any way; it comes straight from Meeus (1998).\n"
            "_C = 0.2\n",
        )
        assert _narrative_hits(tmp_path) == []

    def test_single_copied_output_on_module_constant_is_caught(
        self, tmp_path: Path
    ) -> None:
        # Clause A cannot see a bare copied output (no domain, no
        # tolerance); Clause B is the designed catch for the constant form.
        _write_pkg_module(
            tmp_path,
            "mod_r3.py",
            "# The reference reports 1.5 for this input, so we return 1.5.\n_V = 1.5\n",
        )
        labels = {lb for _p, lb in _narrative_hits(tmp_path)}
        assert "reference-authored-constant" in labels
