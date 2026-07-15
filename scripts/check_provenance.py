#!/usr/bin/env python3
"""Third-party-source provenance sweep (zero-hit CI gate).

Fingerprint greps over the shipped package:

  class 1 — SE source-file references (swehouse, swecl.c, sweph.c, ...)
  class 2 — SE implementation identifiers (swed, dgsect, xs1, mdd, ...)
  class 3 — PyMeeus implementation identifiers (JupiterMoons,
            rectangular_positions_jovian_equatorial, ...)

Classes 1-2 guard against Swiss Ephemeris source identifiers; class 3
against PyMeeus identifiers (``moon_theories/galilean.py`` is an independent
implementation derived from Lieske 1998 / Meeus ch. 44); class 4 flags any
copyleft (L/GPL) license declaration that strays into the tree.

The sweep covers UTF-8 project text across the physical worktree. Name-only
gates run first: one covers the physical tree (including gitignored and
untracked paths), and one covers every commit reachable from every Git ref.
They reject unsafe archive-style paths and block reference-distribution
source/data artifacts without opening them. The history gate disables Git
replace objects and fails closed for non-Git and shallow repositories. Exact
orchestrator checkpoint refs under ``refs/t3/`` are infrastructure and are not
project history. Five path classifications already reachable from the target
and two published pre-review release ancestries are accepted only when every
containing commit belongs to those exact pinned baselines; the same path in any
later commit remains a gate failure.
Symlinks are inspected by link and target name but never followed outside the
worktree. A reachable-history finding does not suppress the independent
physical-worktree content scan; a physical-name finding still stops before any
possibly prohibited payload can be opened. The
gate implementation and its
adversarial filename tests are whole-file allowlisted because they must contain
the blocked tokens. Reviewed legal/history lines use narrow inline waivers;
``libephemeris/`` never accepts waivers. All classes must produce zero hits.

Usage:
    python scripts/check_provenance.py

Provenance:
    Project-authored clean-room policy gate. It scans project paths, text, and
    reachable Git metadata for prohibited provenance fingerprints without
    opening blocked reference-distribution payloads. Pattern literals are audit
    controls, not imported algorithms or data. The script performs no
    astronomical calculation and cannot establish scientific accuracy.
"""

from __future__ import annotations

import os
import re
import stat
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
# The gate and its adversarial tests necessarily contain every blocked token.
# No package/runtime file is ever allowlisted.
ALLOWLIST = frozenset(
    {
        "scripts/check_provenance.py",
        "tests/test_packaging_metadata.py",
        "tests/test_provenance_gate.py",
    }
)

# The target branch already contains these names in commits that predate this
# review. Rewriting that shared ancestry would make a normal PR impossible.
# This is deliberately a path+class+ancestry pin, not a generic history waiver:
# an occurrence in any commit outside the exact baseline ancestry still fails.
INHERITED_HISTORY_BASELINE_COMMITS = (
    "fdd0fe60a7717a45f08f4f0e1823fabba23b8808",  # target main
    "8184e1c2072b615bef77fd0a4a530a4750a49bd1",  # published v3.0.0a5
    "82b44b18dc11dd230b1759e019d63237436afd98",  # published v3.0.0a6
)
INHERITED_HISTORY_ALLOWLIST = frozenset(
    {
        ("docs/reference/swisseph-comparison.md", "foreign-source-file"),
        ("libephemeris/seorbel.txt", "foreign-data-file"),
        ("plans/pyswisseph-api-compat-phase2.md", "foreign-source-file"),
        ("plans/pyswisseph-compat-verification.md", "foreign-source-file"),
        ("tests/test_parse_seorbel.py", "foreign-data-file"),
    }
)

# Documents that legitimately NAME licenses (the GPL nbody extra, historical
# license notes) but must NOT carry SE source-file references, internal
# identifiers or PyMeeus tokens: the root legal/notice files (incl. README,
# the shipped long_description) and the historical release notes.
# These are exempt from the copyleft class ONLY; every other class still
# applies.
LICENSE_NAMING_OK = frozenset(
    {
        "CHANGELOG.md",
        "LICENSE",
        "LICENSING.md",
        "NOTICE.md",
        "RELEASE_NOTES.md",
        "README.md",
        "THIRD_PARTY_NOTICES.md",
    }
)
IMPLEMENTATION_WAIVER_PATHS = frozenset({"CHANGELOG.md", "THIRD_PARTY_NOTICES.md"})
LICENSE_WAIVER_PATHS = frozenset(
    {
        "pyproject.toml",
        "docs/development/architecture-overview.md",
        "docs/guides/getting-started.md",
        "docs/guides/optional-modules.md",
    }
)

SOURCE_FILE_RE = re.compile(
    r"\bswehouse\b|\bswecl\.c\b|\bsweph\.c\b|\bswephlib\b|\bswejpl\b"
    r"|\bswehel\b|\bswemmoon\b|\bswemplan\b|\bswedate\.c\b",
    re.IGNORECASE,
)
IDENTIFIER_RE = re.compile(
    r"\bswed\b|\bdgsect\b|\bxs1\b|\bxh1\b|\bfh1\b|\bmdd\b|\bmdn\b|\badp\b"
    r"|\badmc\b|\bsamc\b|\bacmc\b|\bdfac\b|apc_sector|\bxeq0\b|\bxp0\b"
    r"|\bplaus_iflag\b"
)
# PyMeeus implementation identifiers (confirmed from pymeeus 0.5.12
# JupiterMoons.py). These are PyMeeus's own expression — the independent
# Galilean module uses none of them.
PYMEEUS_RE = re.compile(
    r"\bpymeeus\b|\bJupiterMoons\b|rectangular_positions_jovian_equatorial"
    r"|correct_rectangular_positions|apparent_rectangular_coordinates"
    r"|jupiter_system_angles",
    re.IGNORECASE,
)
# Copyleft license declarations OTHER than AGPL (case-insensitive). Matches
# L/GPL token forms (LGPL, GPL, GPLv2/v3, GPL-2.0/3.0, LGPL-3.0) and the
# spelled-out "GNU/Lesser General Public License". AGPL is the project's
# own license and is NOT flagged. The leading word boundary means the GPL
# inside "AGPL" is not at a token start. Word boundaries also avoid false
# hits on substrings such as the "gPl" inside "PickeringPlanet".
COPYLEFT_RE = re.compile(
    r"\bL?GPL(?:[-\s]?v?[23](?:\.0)?)?\b"
    r"|\blesser\s+general\s+public\b"
    r"|\bGNU\s+general\s+public\b",
    re.IGNORECASE,
)
CLASSES = (
    ("source-file-ref", SOURCE_FILE_RE),
    ("identifier", IDENTIFIER_RE),
    ("pymeeus", PYMEEUS_RE),
    ("copyleft", COPYLEFT_RE),
)
# Data artifacts that must never (re-)enter the tree: the reference
# distribution's ephemeris/orbital-element/star files. This is a
# *filename* gate over the whole repository (text mentions of these names
# in API-compat constants and provenance docs are legitimate; the files
# themselves are not), plus a content scan of the shipped data directory.
FOREIGN_DATA_NAME_RE = re.compile(
    r"(seorbel|sefstars|seleapsec|sedeltat|seasnam)|\.se1(\.|$)",
    re.IGNORECASE,
)
FOREIGN_SOURCE_NAME_RE = re.compile(
    r"^swe[a-z0-9_]*\.(?:c|h)(?:\.|$)",
    re.IGNORECASE,
)
FOREIGN_SOURCE_DIR_RE = re.compile(
    r"^(?:pyswisseph|swisseph|swiss[-_]?ephemeris)(?:[-_.].*)?$",
    re.IGNORECASE,
)
OPAQUE_ARCHIVE_RE = re.compile(
    r"(?:\.zip|\.whl|\.7z|\.rar|\.tgz|\.tbz2?|\.txz|\.tar(?:\.(?:gz|bz2|xz|zst))?)$",
    re.IGNORECASE,
)
COMPILED_BINARY_RE = re.compile(
    r"(?:\.so(?:\.\d+)*|\.dylib|\.dll|\.pyd|\.a|\.lib|\.o|\.obj|\.pyc)$",
    re.IGNORECASE,
)
# Exact infrastructure paths that are external dependency/cache stores rather
# than project artifacts. Exact paths avoid allowing a directory such as
# ``libephemeris/data/.venv`` to become a hiding place.
CONTENT_INFRASTRUCTURE_SKIP_PATHS = frozenset(
    {
        (".git",),
        (".hypothesis",),
        (".venv",),
        (".mypy_cache",),
        (".pytest_cache",),
        (".ruff_cache",),
        (".opencode", "node_modules"),
        ("validation", ".git"),
        ("validation", ".venv"),
        ("validation", ".pytest_cache"),
        ("validation", ".ruff_cache"),
    }
)
# VCS directories are traversed by the physical name gate so that even an
# ignored nested repository cannot hide a prohibited artifact name. Their
# internal object/index payloads remain outside the subsequent content scan:
# inspecting or decompressing repository history is neither necessary nor
# appropriate for this clean-room gate.
NAME_INFRASTRUCTURE_SKIP_PATHS = CONTENT_INFRASTRUCTURE_SKIP_PATHS - {
    (".git",),
    ("validation", ".git"),
}

# Large binary assets are still covered by the name-only artifact gate. They
# are not decoded during the implementation-fingerprint pass.
BINARY_SUFFIXES = frozenset(
    {
        ".7z",
        ".a",
        ".bin",
        ".bsp",
        ".bz2",
        ".coverage",
        ".db",
        ".dll",
        ".dylib",
        ".gif",
        ".gz",
        ".ico",
        ".jpeg",
        ".jpg",
        ".leb",
        ".leb2",
        ".lib",
        ".mov",
        ".mp3",
        ".mp4",
        ".npy",
        ".npz",
        ".o",
        ".obj",
        ".otf",
        ".pdf",
        ".pkl",
        ".png",
        ".pyc",
        ".pyd",
        ".rar",
        ".so",
        ".sqlite",
        ".ttf",
        ".whl",
        ".woff",
        ".woff2",
        ".xz",
        ".zip",
        ".zst",
    }
)
BINARY_FILENAMES = frozenset({".coverage", ".ds_store"})
BINARY_DOWNLOAD_RE = re.compile(r"\.bsp\.download\d*$", re.IGNORECASE)


@dataclass(frozen=True)
class HistoryArtifact:
    """One deduplicated reachable-history name-gate result."""

    path: str
    label: str
    commits: tuple[str, ...] = ()


def _normalize_artifact_parts(
    name: str, *, reject_unsafe: bool = True
) -> tuple[tuple[str, ...], str | None]:
    """Normalize an archive/worktree path using cross-platform name semantics."""
    if "\x00" in name:
        return (), "unsafe-path"
    normalized = name.replace("\\", "/")
    if reject_unsafe and (
        normalized.startswith("/") or re.match(r"^[A-Za-z]:/", normalized)
    ):
        return (), "unsafe-path"

    parts: list[str] = []
    for raw_part in normalized.split("/"):
        if not raw_part or raw_part == ".":
            continue
        canonical = raw_part.rstrip(" .").casefold()
        if canonical == "..":
            if reject_unsafe:
                return (), "unsafe-path"
            continue
        if not canonical:
            return (), "unsafe-path"
        parts.append(canonical)
    return tuple(parts), None


def _canonical_policy_rel(name: str) -> str:
    """Return a case-folded, platform-neutral path for content policy checks."""
    parts, error = _normalize_artifact_parts(name)
    return "" if error is not None else "/".join(parts)


def _is_binary_payload_name(name: str) -> bool:
    """Return whether a safe-named payload is a known project binary."""
    basename = name.replace("\\", "/").rsplit("/", maxsplit=1)[-1]
    canonical_basename = basename.rstrip(" .").casefold()
    return (
        Path(canonical_basename).suffix in BINARY_SUFFIXES
        or canonical_basename in BINARY_FILENAMES
        or bool(BINARY_DOWNLOAD_RE.search(canonical_basename))
    )


def _classify_artifact_name(
    name: str,
    *,
    is_directory: bool = False,
    reject_unsafe: bool = True,
    classify_generic_reference_dir: bool = True,
) -> str | None:
    """Classify a prohibited path from names alone, without opening payloads."""
    parts, error = _normalize_artifact_parts(name, reject_unsafe=reject_unsafe)
    if error is not None:
        return error
    if not parts:
        return None

    if classify_generic_reference_dir:
        for index in range(len(parts) - 1):
            if parts[index : index + 2] == ("data", "reference"):
                return "foreign-data-dir"

    for index, part in enumerate(parts):
        component_is_dir = is_directory or index < len(parts) - 1
        if FOREIGN_DATA_NAME_RE.search(part):
            return "foreign-data-dir" if component_is_dir else "foreign-data-file"
        if FOREIGN_SOURCE_DIR_RE.fullmatch(part):
            return "foreign-source-dir" if component_is_dir else "foreign-source-file"

    basename = parts[-1]
    if FOREIGN_SOURCE_NAME_RE.search(basename):
        return "foreign-source-dir" if is_directory else "foreign-source-file"
    if OPAQUE_ARCHIVE_RE.search(basename):
        return "opaque-archive"
    if COMPILED_BINARY_RE.search(basename):
        return "compiled-binary"
    return None


def _run_git_metadata(root: Path, *args: str) -> subprocess.CompletedProcess[bytes]:
    """Run a read-only Git metadata command with replacements disabled."""
    return subprocess.run(
        ["git", "--no-replace-objects", "-C", os.fspath(root), *args],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _git_failure_detail(result: subprocess.CompletedProcess[bytes]) -> str:
    """Return a bounded, display-safe diagnostic for a failed Git command."""
    detail = result.stderr.decode("utf-8", errors="backslashreplace").strip()
    if not detail:
        detail = f"git exited with status {result.returncode}"
    return detail[:500]


def _display_history_path(path: str) -> str:
    """Escape control characters so one history hit occupies one output line."""
    return path.encode("unicode_escape", errors="backslashreplace").decode("ascii")


def _reachable_history_foreign_artifacts(root: Path) -> list[HistoryArtifact]:
    """Classify prohibited path names in all reachable Git commits.

    This function reads commit/tree metadata only. It deliberately never uses
    ``git show``, ``git cat-file``, or any other command that opens historical
    blob payloads.

    Returns:
        Deduplicated findings. Repository/metadata errors are represented as
        findings so callers fail closed.
    """
    try:
        inside = _run_git_metadata(root, "rev-parse", "--is-inside-work-tree")
    except OSError as error:
        return [HistoryArtifact("<repository>", "history-git-error", (str(error),))]
    if inside.returncode != 0 or inside.stdout.strip() != b"true":
        detail = (
            _git_failure_detail(inside)
            if inside.returncode != 0
            else "not a Git worktree"
        )
        return [HistoryArtifact("<repository>", "history-not-git", (detail,))]

    shallow = _run_git_metadata(root, "rev-parse", "--is-shallow-repository")
    if shallow.returncode != 0:
        return [
            HistoryArtifact(
                "<repository>",
                "history-git-error",
                (_git_failure_detail(shallow),),
            )
        ]
    shallow_state = shallow.stdout.strip()
    if shallow_state == b"true":
        return [
            HistoryArtifact(
                "<repository>",
                "history-shallow-repository",
                ("reachable history is incomplete",),
            )
        ]
    if shallow_state != b"false":
        return [
            HistoryArtifact(
                "<repository>",
                "history-metadata-error",
                (f"unexpected shallow status: {shallow_state!r}",),
            )
        ]

    revisions = _run_git_metadata(root, "rev-list", "--exclude=refs/t3/*", "--all")
    if revisions.returncode != 0:
        return [
            HistoryArtifact(
                "<repository>",
                "history-git-error",
                (_git_failure_detail(revisions),),
            )
        ]

    commit_ids: list[str] = []
    for raw_commit in revisions.stdout.splitlines():
        try:
            commit = raw_commit.decode("ascii", errors="strict")
        except UnicodeDecodeError:
            return [
                HistoryArtifact(
                    "<repository>",
                    "history-metadata-error",
                    (f"non-ASCII revision id: {raw_commit!r}",),
                )
            ]
        if not re.fullmatch(r"[0-9a-fA-F]+", commit):
            return [
                HistoryArtifact(
                    "<repository>",
                    "history-metadata-error",
                    (f"invalid revision id: {commit!r}",),
                )
            ]
        commit_ids.append(commit)

    hits: dict[tuple[str, str], set[str]] = {}
    for commit in commit_ids:
        tree = _run_git_metadata(
            root,
            "ls-tree",
            "-rz",
            "--full-tree",
            "--name-only",
            commit,
            "--",
        )
        if tree.returncode != 0:
            return [
                HistoryArtifact(
                    "<repository>",
                    "history-git-error",
                    (f"{commit}: {_git_failure_detail(tree)}",),
                )
            ]
        if tree.stdout and not tree.stdout.endswith(b"\x00"):
            return [
                HistoryArtifact(
                    "<repository>",
                    "history-metadata-error",
                    (f"{commit}: unterminated ls-tree output",),
                )
            ]
        for raw_path in tree.stdout.split(b"\x00"):
            if not raw_path:
                continue
            try:
                path = raw_path.decode("utf-8", errors="strict")
            except UnicodeDecodeError:
                displayed = repr(raw_path)
                hits.setdefault((displayed, "history-path-decode-error"), set()).add(
                    commit
                )
                continue
            label = _classify_artifact_name(path)
            if label is not None:
                hits.setdefault((path, label), set()).add(commit)

    return [
        HistoryArtifact(path, label, tuple(sorted(commits)))
        for (path, label), commits in sorted(hits.items())
    ]


def _partition_history_baseline(
    root: Path, artifacts: list[HistoryArtifact]
) -> tuple[list[HistoryArtifact], list[HistoryArtifact]]:
    """Split exact inherited target-history names from actionable findings."""
    baseline = _run_git_metadata(root, "rev-list", *INHERITED_HISTORY_BASELINE_COMMITS)
    if baseline.returncode != 0:
        failure = HistoryArtifact(
            "<repository>",
            "history-baseline-error",
            (_git_failure_detail(baseline),),
        )
        return [], [*artifacts, failure]

    baseline_commits = {
        commit.decode("ascii", errors="strict")
        for commit in baseline.stdout.splitlines()
    }
    inherited: list[HistoryArtifact] = []
    actionable: list[HistoryArtifact] = []
    for artifact in artifacts:
        identity = (artifact.path, artifact.label)
        if (
            identity in INHERITED_HISTORY_ALLOWLIST
            and artifact.commits
            and set(artifact.commits) <= baseline_commits
        ):
            inherited.append(artifact)
        else:
            actionable.append(artifact)
    return inherited, actionable


def _skip_infrastructure_directory(
    root: Path,
    path: Path,
    *,
    content_scan: bool = False,
) -> bool:
    """Return whether an exact external/cache directory may be skipped."""
    try:
        parts = path.relative_to(root).parts
    except ValueError:
        return False
    skip_paths = (
        CONTENT_INFRASTRUCTURE_SKIP_PATHS
        if content_scan
        else NAME_INFRASTRUCTURE_SKIP_PATHS
    )
    return parts in skip_paths


def _physical_foreign_artifacts(root: Path) -> list[tuple[Path, str, str]]:
    """Find prohibited reference artifacts by path name without reading them.

    Returns:
        Tuples of ``(path, label, displayed_name)`` for each prohibited path.
    """
    hits: list[tuple[Path, str, str]] = []

    def _walk_error(error: OSError) -> None:
        path = Path(error.filename) if error.filename else root
        hits.append((path, "walk-error", str(error)))

    # External/cache trees may legitimately contain generic archives and
    # compiled dependencies, so the main name gate prunes them. Still inspect
    # their *names* for reference-distribution-specific markers: a local
    # reference-API installation belongs outside the project worktree, and an ignored cache
    # must not become a hiding place for its source or data files. Payloads are
    # never opened here.
    for skipped_parts in sorted(NAME_INFRASTRUCTURE_SKIP_PATHS):
        skipped_root = root.joinpath(*skipped_parts)
        if not skipped_root.is_dir():
            continue
        for current, dirnames, filenames in os.walk(
            skipped_root, followlinks=False, onerror=_walk_error
        ):
            current_path = Path(current)
            dirnames[:] = sorted(dirnames)
            for dirname in tuple(dirnames):
                path = current_path / dirname
                rel = path.relative_to(root)
                label = _classify_artifact_name(
                    rel.as_posix(),
                    is_directory=True,
                    classify_generic_reference_dir=False,
                )
                if label is not None and label.startswith("foreign-"):
                    hits.append((path, label, rel.as_posix()))
                    dirnames.remove(dirname)
            for filename in sorted(filenames):
                path = current_path / filename
                rel = path.relative_to(root)
                label = _classify_artifact_name(
                    rel.as_posix(), classify_generic_reference_dir=False
                )
                if label is not None and label.startswith("foreign-"):
                    hits.append((path, label, rel.as_posix()))

    for current, dirnames, filenames in os.walk(
        root, followlinks=False, onerror=_walk_error
    ):
        current_path = Path(current)
        dirnames[:] = sorted(
            dirname
            for dirname in dirnames
            if not _skip_infrastructure_directory(root, current_path / dirname)
        )

        for dirname in tuple(dirnames):
            path = current_path / dirname
            rel = path.relative_to(root)
            label = _classify_artifact_name(rel.as_posix(), is_directory=True)
            if label is not None:
                hits.append((path, label, rel.as_posix()))
                dirnames.remove(dirname)
                continue
            if path.is_symlink():
                try:
                    target = os.readlink(path)
                except OSError as error:
                    hits.append((path, "symlink-error", str(error)))
                else:
                    target_label = _classify_artifact_name(
                        target, is_directory=True, reject_unsafe=False
                    )
                    if target_label is not None:
                        hits.append((path, "foreign-symlink-target", target))
                dirnames.remove(dirname)

        for filename in sorted(filenames):
            path = current_path / filename
            rel = path.relative_to(root)
            try:
                mode = path.lstat().st_mode
            except OSError as error:
                hits.append((path, "stat-error", str(error)))
                continue
            label = _classify_artifact_name(rel.as_posix())
            if label is not None:
                hits.append((path, label, rel.as_posix()))
                continue
            if stat.S_ISLNK(mode):
                try:
                    target = os.readlink(path)
                except OSError as error:
                    hits.append((path, "symlink-error", str(error)))
                else:
                    target_label = _classify_artifact_name(target, reject_unsafe=False)
                    if target_label is not None:
                        hits.append((path, "foreign-symlink-target", target))
                continue
            if not stat.S_ISREG(mode):
                hits.append((path, "special-file", rel.as_posix()))
    return hits


def _content_scan_files(root: Path) -> list[Path]:
    """Return UTF-8 project text after the physical name gate is clean."""
    files: list[Path] = []
    for current, dirnames, filenames in os.walk(root, followlinks=False):
        current_path = Path(current)
        dirnames[:] = sorted(
            dirname
            for dirname in dirnames
            if not _skip_infrastructure_directory(
                root,
                current_path / dirname,
                content_scan=True,
            )
        )
        for filename in sorted(filenames):
            path = current_path / filename
            if path.is_symlink():
                continue
            try:
                if not stat.S_ISREG(path.lstat().st_mode):
                    # The physical name gate rejects special files before the
                    # content scan. Avoid blocking on a FIFO if this helper is
                    # called independently.
                    continue
            except OSError:
                files.append(path)
                continue
            if _is_binary_payload_name(path.name):
                continue
            try:
                with path.open("rb") as stream:
                    prefix = stream.read(8192)
            except OSError:
                # Keep unreadable project files in the result so the main scan
                # reports a fail-closed content-read-error.
                files.append(path)
                continue
            try:
                prefix.decode("utf-8", errors="strict")
            except UnicodeDecodeError:
                # Unknown non-UTF-8 payloads are not silently treated as binary;
                # main() will report a fail-closed content-read-error.
                files.append(path)
                continue
            files.append(path)
    return sorted(files)


def _line_hit_labels(rel: str, line: str) -> list[str]:
    """Return provenance classes hit by one project-owned text line."""
    exact_rel = rel.replace("\\", "/")
    canonical_rel = _canonical_policy_rel(exact_rel)
    path = Path(canonical_rel)
    license_rel = exact_rel
    if exact_rel in {"PKG-INFO", "libephemeris.egg-info/PKG-INFO"}:
        license_rel = "README.md"
    in_package = bool(path.parts and path.parts[0] == "libephemeris")
    implementation_waived = (
        "provenance-implementation-ok" in line
        and exact_rel in IMPLEMENTATION_WAIVER_PATHS
        and not in_package
    )
    license_waived = (
        "provenance-ok" in line and exact_rel in LICENSE_WAIVER_PATHS and not in_package
    )
    license_naming_ok = (
        license_waived
        or license_rel in LICENSE_NAMING_OK
        or exact_rel.startswith("release-notes/")
    )

    labels: list[str] = []
    for label, pattern in CLASSES:
        if implementation_waived and label in (
            "source-file-ref",
            "identifier",
            "pymeeus",
        ):
            continue
        if label == "copyleft" and license_naming_ok:
            continue
        if pattern.search(line):
            labels.append(label)
    return labels


# ---------------------------------------------------------------------------
# High-precision numeric-literal gate (shipped modules).
#
# A float literal with HIGH_PRECISION_SIG_DIGITS or more significant digits
# must carry a source annotation: unexplained full-precision constants are
# the signature of values read from an external implementation's output.
# A literal is annotated when a citation token appears on its own line or in
# the preceding comment/docstring block, or when the module is an
# AUTO-GENERATED data table whose committed generator documents the full
# derivation (those modules are separately hash-pinned).
# ---------------------------------------------------------------------------

HIGH_PRECISION_SIG_DIGITS = 10

CITATION_TOKEN_RE = re.compile(
    r"IAU|IERS|JPL|DE4\d\d|Hipparcos|Gaia|ERFA|SOFA|Meeus|Simon|Vondr|"
    r"Chapront|Folkner|Park et al|Standish|Lieske|van Leeuwen|XHIP|Capitaine|"
    r"Mathews|Reid|Brunthaler|Burgess|Clark|Neely|Strubell|Harrington|"
    r"Hoyt|published|source:|Source:|References?:|convention|CODATA|VSOP|"
    r"ELP\s?2000|Williams|Boggs|Petit|Luzum|WGS84|Astronomical Almanac|"
    r"Explanatory Supplement|Espenak|Morrison|Stephenson|Bretagnon|Francou|"
    r"Laskar|Newhall|Seidelmann|Kaplan|Wallace|Bowell|Fagan|Firebrace|"
    r"Sheoran|Mardyks|Gil Brand|Wilhelm|Mercier|Huber|Kugler|Britton|"
    r"Lahiri|Krishnamurti|Raman|Yukteswar|De Luce|Bailey|Fiorenza|Valens|"
    r"Kali|Besselian|Julian epoch|frame epoch|golden section|Sgr A|"
    r"derived|measured|fitted at the DE440|generate_|"
    r"speed of light|astronomical unit|AU in km|1 AU|Gauss|gravitational "
    r"constant|synodic|nominal|Horizons|van Gent|Lieske 1998|Allen \(19|Samaha",
    re.IGNORECASE,
)

# A data module is exempt when its header points at the committed generator
# that reproduces it (the derivation lives in the generator; the module is
# separately reviewed/pinned).
GENERATED_HEADER_RE = re.compile(r"AUTO-GENERATED|Generated by scripts/", re.IGNORECASE)

_FLOAT_DIGITS_RE = re.compile(r"^[0-9.]+$")


def _significant_digits(literal: str) -> int:
    """Count significant digits of a float literal (mantissa only)."""
    mantissa = literal.lower().split("e")[0]
    digits = mantissa.replace(".", "").replace("_", "").lstrip("0")
    return len(digits.rstrip("0")) if digits else 0


def _statement_spans(text: str) -> list[tuple[int, int]]:
    """(lineno, end_lineno) spans of every statement in the module."""
    import ast

    try:
        tree = ast.parse(text)
    except SyntaxError:
        return []
    spans: list[tuple[int, int]] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.stmt):
            spans.append((node.lineno, node.end_lineno or node.lineno))
    return spans


def _high_precision_literal_hits(root: Path) -> list[tuple[Path, int, str, str]]:
    """Find unexplained high-precision literals outside generated/vendor code."""
    import io
    import tokenize

    hits: list[tuple[Path, int, str, str]] = []
    package = root / "libephemeris"
    for path in sorted(package.rglob("*.py")):
        rel_parts = path.relative_to(package).parts
        if "vendor" in rel_parts or "__pycache__" in rel_parts:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        head = "\n".join(text.splitlines()[:60])
        if GENERATED_HEADER_RE.search(head):
            continue
        lines = text.splitlines()
        spans = _statement_spans(text)
        try:
            tokens = list(tokenize.generate_tokens(io.StringIO(text).readline))
        except tokenize.TokenizeError:
            continue

        def _annotated(lineno: int) -> bool:
            # Local context: the literal's line and the 7 lines above it.
            lo = max(0, lineno - 8)
            if CITATION_TOKEN_RE.search("\n".join(lines[lo : lineno + 1])):
                return True
            # Statement context: the whole enclosing statement plus up to 12
            # comment lines above it — the natural place for a data table's
            # source note.
            for start, end in spans:
                if start <= lineno <= end:
                    top = max(0, start - 13)
                    header = "\n".join(lines[top:end])
                    if CITATION_TOKEN_RE.search(header):
                        return True
            return False

        flagged_lines: set[int] = set()
        for tok in tokens:
            if tok.type != tokenize.NUMBER:
                continue
            literal = tok.string
            if "." not in literal and "e" not in literal.lower():
                continue
            if literal.lower().endswith("j"):
                continue
            if _significant_digits(literal) < HIGH_PRECISION_SIG_DIGITS:
                continue
            lineno = tok.start[0]
            if lineno in flagged_lines:
                continue
            if _annotated(lineno):
                continue
            flagged_lines.add(lineno)
            hits.append(
                (
                    path,
                    lineno,
                    "unannotated-high-precision-literal",
                    lines[lineno - 1].strip(),
                )
            )
    return hits


def main() -> int:
    """Run every repository independence and provenance policy check."""
    name_hits: list[tuple[Path, int, str, str]] = []

    # Physical name gate: no reference-distribution directory, source file or
    # data file anywhere in the worktree, whether tracked, untracked or ignored.
    # This deliberately inspects path names only and never opens these artifacts.
    for path, label, displayed_name in _physical_foreign_artifacts(REPO_ROOT):
        name_hits.append((path, 0, label, displayed_name))

    # Reachable-history name gate: enumerate commit/tree metadata under every
    # ref with replace objects disabled. Historical blob contents are never
    # opened. Results are deduplicated by path and classification.
    all_history_hits = _reachable_history_foreign_artifacts(REPO_ROOT)
    inherited_history_hits, history_hits = _partition_history_baseline(
        REPO_ROOT, all_history_hits
    )

    # A prohibited physical name may designate an artifact whose payload must
    # never be inspected here, so fail closed before the content scan in that
    # case. Historical findings are commit/tree metadata only: report them, but
    # still scan the current physical worktree so CI exposes both independent
    # classes of failure in one run.
    if name_hits:
        for path, lineno, label, line in name_hits:
            print(f"{path.relative_to(REPO_ROOT)}:{lineno}: [{label}] {line[:100]}")
        for artifact in inherited_history_hits:
            displayed_path = _display_history_path(artifact.path)
            commit_list = ",".join(commit[:12] for commit in artifact.commits[:5])
            remainder = len(artifact.commits) - 5
            suffix = f",+{remainder} more" if remainder > 0 else ""
            print(
                f"history:{displayed_path}:0: [acknowledged-inherited-history] "
                "baselines "
                f"{','.join(commit[:12] for commit in INHERITED_HISTORY_BASELINE_COMMITS)}; "
                f"reachable in {commit_list}{suffix}"
            )
        for artifact in history_hits:
            displayed_path = _display_history_path(artifact.path)
            if artifact.label.startswith("history-"):
                detail = "; ".join(artifact.commits)
            else:
                commit_list = ",".join(commit[:12] for commit in artifact.commits[:5])
                remainder = len(artifact.commits) - 5
                suffix = f",+{remainder} more" if remainder > 0 else ""
                detail = f"reachable in {commit_list}{suffix}"
            print(f"history:{displayed_path}:0: [{artifact.label}] {detail[:500]}")
        print(
            "provenance sweep: "
            f"{len(name_hits)} physical-name hit(s), "
            f"{len(inherited_history_hits)} acknowledged inherited-history hit(s), "
            f"{len(history_hits)} new reachable-history hit(s), "
            "content scan skipped for physical-name safety (gate: 0)"
        )
        return 1

    content_hits: list[tuple[Path, int, str, str]] = []
    files = _content_scan_files(REPO_ROOT)
    scanned = 0
    for path in files:
        rel = path.relative_to(REPO_ROOT).as_posix()
        if rel in ALLOWLIST:
            continue
        scanned += 1
        try:
            text = path.read_text(encoding="utf-8", errors="strict")
        except (OSError, UnicodeDecodeError) as error:
            content_hits.append((path, 0, "content-read-error", str(error)))
            continue
        for lineno, line in enumerate(text.splitlines(), start=1):
            for label in _line_hit_labels(rel, line):
                content_hits.append((path, lineno, label, line.strip()))

    for artifact in inherited_history_hits:
        displayed_path = _display_history_path(artifact.path)
        commit_list = ",".join(commit[:12] for commit in artifact.commits[:5])
        remainder = len(artifact.commits) - 5
        suffix = f",+{remainder} more" if remainder > 0 else ""
        print(
            f"history:{displayed_path}:0: [acknowledged-inherited-history] "
            "baselines "
            f"{','.join(commit[:12] for commit in INHERITED_HISTORY_BASELINE_COMMITS)}; "
            f"reachable in {commit_list}{suffix}"
        )
    for artifact in history_hits:
        displayed_path = _display_history_path(artifact.path)
        if artifact.label.startswith("history-"):
            detail = "; ".join(artifact.commits)
        else:
            commit_list = ",".join(commit[:12] for commit in artifact.commits[:5])
            remainder = len(artifact.commits) - 5
            suffix = f",+{remainder} more" if remainder > 0 else ""
            detail = f"reachable in {commit_list}{suffix}"
        print(f"history:{displayed_path}:0: [{artifact.label}] {detail[:500]}")
    literal_hits = _high_precision_literal_hits(REPO_ROOT)
    for path, lineno, label, line in literal_hits:
        print(f"{path.relative_to(REPO_ROOT)}:{lineno}: [{label}] {line[:100]}")

    for path, lineno, label, line in content_hits:
        print(f"{path.relative_to(REPO_ROOT)}:{lineno}: [{label}] {line[:100]}")
    print(
        "provenance sweep: "
        f"{len(name_hits)} physical-name hit(s), "
        f"{len(inherited_history_hits)} acknowledged inherited-history hit(s), "
        f"{len(history_hits)} new reachable-history hit(s), "
        f"{len(content_hits)} content hit(s) over {scanned} project files, "
        f"{len(literal_hits)} unannotated high-precision literal(s) "
        "in the physical worktree (gate: 0)"
    )
    return 1 if history_hits or content_hits or literal_hits else 0


if __name__ == "__main__":
    sys.exit(main())
