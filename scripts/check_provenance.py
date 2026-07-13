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

The sweep covers UTF-8 project text across the physical worktree. A name-only
gate runs first (including gitignored and untracked paths), rejects unsafe
archive-style paths, and blocks reference-distribution source/data artifacts
without opening them. Symlinks are inspected by link and target name but never
followed outside the worktree. Only after that gate is clean does a content
scan run. The gate implementation and its
adversarial filename tests are whole-file allowlisted because they must contain
the blocked tokens. Reviewed legal/history lines use narrow inline waivers;
``libephemeris/`` never accepts waivers. All classes must produce zero hits.

Usage:
    python scripts/check_provenance.py
"""

from __future__ import annotations

import os
import re
import stat
import sys
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

# Documents that legitimately NAME licenses (the retired AGPL/dual era, the
# optional GPL nbody extra) but must NOT carry SE source-file references,
# internal identifiers or PyMeeus tokens: the root legal/notice files (incl.
# README, the shipped long_description) and the historical release notes.
# These are exempt from the copyleft/agpl classes ONLY; every other class
# still applies.
LICENSE_NAMING_OK = frozenset(
    {
        "CHANGELOG.md",
        "HANDOFF.md",
        "LICENSING.md",
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
# Copyleft license declarations (case-insensitive). Matches L/GPL token
# forms (LGPL, GPL, GPLv2/v3, GPL-2.0/3.0, LGPL-3.0) and the spelled-out
# "GNU/Lesser General Public License", but NOT AGPL (which now appears only
# in historical release notes and the dev-only test oracle, never in this
# project's shipped Apache-2.0 headers): the leading word boundary means the
# GPL inside "AGPL" is not at a token start, and "GNU Affero General Public"
# breaks the GNU...general adjacency. Word boundaries also avoid false hits
# on substrings such as the "gPl" inside "PickeringPlanet".
COPYLEFT_RE = re.compile(
    r"\bL?GPL(?:[-\s]?v?[23](?:\.0)?)?\b"
    r"|\blesser\s+general\s+public\b"
    r"|\bGNU\s+general\s+public\b",
    re.IGNORECASE,
)
# AGPL is its own class: the reference implementation's license. It must
# never appear in shipped code or data; docs may name it only in the
# allowlisted legal/notice files.
AGPL_RE = re.compile(r"\bAGPL\b|\bGNU\s+Affero\b", re.IGNORECASE)
CLASSES = (
    ("source-file-ref", SOURCE_FILE_RE),
    ("identifier", IDENTIFIER_RE),
    ("pymeeus", PYMEEUS_RE),
    ("copyleft", COPYLEFT_RE),
    ("agpl", AGPL_RE),
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
) -> str | None:
    """Classify a prohibited path from names alone, without opening payloads."""
    parts, error = _normalize_artifact_parts(name, reject_unsafe=reject_unsafe)
    if error is not None:
        return error
    if not parts:
        return None

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
        if label in ("copyleft", "agpl") and license_naming_ok:
            continue
        if pattern.search(line):
            labels.append(label)
    return labels


def main() -> int:
    hits: list[tuple[Path, int, str, str]] = []

    # Physical name gate: no reference-distribution directory, source file or
    # data file anywhere in the worktree, whether tracked, untracked or ignored.
    # This deliberately inspects path names only and never opens these artifacts.
    for path, label, displayed_name in _physical_foreign_artifacts(REPO_ROOT):
        hits.append((path, 0, label, displayed_name))

    # Fail closed before any content read. A prohibited name may designate an
    # artifact whose payload must never be inspected in this clean room.
    if hits:
        for path, lineno, label, line in hits:
            print(f"{path.relative_to(REPO_ROOT)}:{lineno}: [{label}] {line[:100]}")
        print(f"provenance sweep: {len(hits)} physical-name hit(s) (gate: 0)")
        return 1

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
            hits.append((path, 0, "content-read-error", str(error)))
            continue
        for lineno, line in enumerate(text.splitlines(), start=1):
            for label in _line_hit_labels(rel, line):
                hits.append((path, lineno, label, line.strip()))

    for path, lineno, label, line in hits:
        print(f"{path.relative_to(REPO_ROOT)}:{lineno}: [{label}] {line[:100]}")
    print(
        f"provenance sweep: {len(hits)} hit(s) over {scanned} project files "
        "in the physical worktree (gate: 0)"
    )
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
