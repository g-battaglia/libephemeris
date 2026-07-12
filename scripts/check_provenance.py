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

The sweep covers ``libephemeris/``, ``scripts/``, ``docs/`` and
``release-notes/`` plus the root legal/packaging files, across every
shipped text format. A name-only gate walks the physical worktree (including
gitignored and untracked paths) and blocks reference-distribution source/data
files and ``data/reference`` directories without opening them; a content scan
flags copyleft markers in ``libephemeris/data/``. Two allowlists: a
full one for the gate and provenance/legal docs that must quote retired-
source identifiers to record what was removed, and a license-naming-only
one (README, release notes) exempt from the copyleft/AGPL classes but still
checked for SE identifiers and PyMeeus tokens. ``libephemeris/`` is never
exempt. All classes must produce zero hits; this gate keeps the tree that
way.

Usage:
    python scripts/check_provenance.py
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCAN_DIRS = ("libephemeris", "scripts", "docs", "release-notes")
# Root-level files that legally represent the package (long_description,
# notices) are scanned individually in addition to the directories.
ROOT_FILES = (
    "pyproject.toml",
    "README.md",
    "NOTICE.md",
    "LICENSING.md",
    "THIRD_PARTY_NOTICES.md",
    "CHANGELOG.md",
    "RELEASE_NOTES.md",
)
# Every text format that ships or documents the package: code, docs, data
# tables, packaging and tool configuration. Binary blobs (.bsp/.leb/.leb2)
# are covered by the dedicated data-artifact class below.
SCAN_SUFFIXES = (
    ".py",
    ".md",
    ".csv",
    ".toml",
    ".cfg",
    ".ini",
    ".txt",
    ".json",
    ".yaml",
    ".yml",
    ".rst",
)

# Files that legitimately record the retired SE / PyMeeus / LGPL history at
# the IMPLEMENTATION level (the gate itself and the provenance-evidence
# docs, which must quote SE source-file names and internal identifiers to
# describe what was removed). These are exempt from ALL classes. Never
# includes any libephemeris/ file — shipped code stays strictly zero-hit.
ALLOWLIST = frozenset(
    {
        "scripts/check_provenance.py",
        "docs/methodology/galilean-e5-spec.md",
        # The legal/notice files may name retired-source identifiers
        # (swehouse.c, PyMeeus, LGPL) to describe what is not present, so
        # every class legitimately appears.
        "NOTICE.md",
        "LICENSING.md",
        "THIRD_PARTY_NOTICES.md",
        "CHANGELOG.md",
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
        "RELEASE_NOTES.md",
        "README.md",
    }
)

SOURCE_FILE_RE = re.compile(
    r"swehouse|swecl\.c|sweph\.c|swephlib|swejpl|swehel"
    r"|swemmoon|swemplan|swedate\.c",
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
# VCS metadata, local virtual environments and interpreter/tool caches are not
# worktree artifacts. Everything else is examined by name, including ignored
# output/build trees and nested tool-created worktrees.
PHYSICAL_SCAN_SKIP_DIRS = frozenset(
    {".git", ".venv", "__pycache__", ".mypy_cache", ".pytest_cache", ".ruff_cache"}
)
DATA_DIR = "libephemeris/data"


def _physical_foreign_artifacts(root: Path) -> list[tuple[Path, str, str]]:
    """Find prohibited reference artifacts by path name without reading them.

    Returns:
        Tuples of ``(path, label, displayed_name)`` for each prohibited path.
    """
    hits: list[tuple[Path, str, str]] = []
    for current, dirnames, filenames in os.walk(root, followlinks=False):
        current_path = Path(current)
        dirnames[:] = sorted(d for d in dirnames if d not in PHYSICAL_SCAN_SKIP_DIRS)

        # A data/reference directory is prohibited even when empty or ignored.
        # Record it once and prune it; its contents need not be named or opened.
        for dirname in tuple(dirnames):
            path = current_path / dirname
            rel = path.relative_to(root)
            folded_parts = tuple(part.casefold() for part in rel.parts[-2:])
            if len(rel.parts) >= 2 and folded_parts == ("data", "reference"):
                hits.append((path, "foreign-data-dir", rel.as_posix()))
                dirnames.remove(dirname)
            elif FOREIGN_SOURCE_DIR_RE.fullmatch(dirname):
                hits.append((path, "foreign-source-dir", rel.as_posix()))
                dirnames.remove(dirname)

        for filename in sorted(filenames):
            path = current_path / filename
            rel = path.relative_to(root)
            if FOREIGN_DATA_NAME_RE.search(filename):
                hits.append((path, "foreign-data-file", filename))
            elif FOREIGN_SOURCE_NAME_RE.search(filename):
                hits.append((path, "foreign-source-file", filename))
    return hits


def main() -> int:
    hits: list[tuple[Path, int, str, str]] = []

    # Physical name gate: no reference-distribution directory, source file or
    # data file anywhere in the worktree, whether tracked, untracked or ignored.
    # This deliberately inspects path names only and never opens these artifacts.
    for path, label, displayed_name in _physical_foreign_artifacts(REPO_ROOT):
        hits.append((path, 0, label, displayed_name))

    # Content scan of the shipped data directory (any text file).
    data_dir = REPO_ROOT / DATA_DIR
    if data_dir.exists():
        for p in sorted(data_dir.rglob("*")):
            if not p.is_file() or p.suffix in (".bsp", ".leb", ".leb2"):
                continue
            try:
                text = p.read_text(encoding="utf-8", errors="strict")
            except (UnicodeDecodeError, OSError):
                continue
            for lineno, line in enumerate(text.splitlines(), start=1):
                if AGPL_RE.search(line) or COPYLEFT_RE.search(line):
                    hits.append((p, lineno, "data-license-marker", line.strip()))
    files = sorted(
        [
            p
            for d in SCAN_DIRS
            for p in (REPO_ROOT / d).rglob("*")
            if p.suffix in SCAN_SUFFIXES and "__pycache__" not in p.parts
        ]
        + [REPO_ROOT / f for f in ROOT_FILES if (REPO_ROOT / f).exists()]
    )
    scanned = 0
    for path in files:
        rel = path.relative_to(REPO_ROOT).as_posix()
        if rel in ALLOWLIST:
            continue
        scanned += 1
        for lineno, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(), start=1
        ):
            # An explicit inline waiver for individually-reviewed lines.
            # It applies ONLY to the license-naming classes (copyleft/agpl
            # — install docs must be able to name an optional extra's GPL
            # license) and never inside libephemeris/.
            waived = "provenance-ok" in line and not path.is_relative_to(
                REPO_ROOT / "libephemeris"
            )
            # The license-naming classes (copyleft/agpl) are expected in:
            # the root legal/notice docs (LICENSE_NAMING_OK), the historical
            # release notes (which record the retired AGPL/dual era), and
            # any line carrying an explicit reviewed waiver. Every OTHER
            # class (SE source-file refs, internal identifiers, PyMeeus)
            # still applies to all of these — a stray SE identifier in
            # README or a release note is still a failure.
            license_naming_ok = (
                waived
                or rel in LICENSE_NAMING_OK
                or path.is_relative_to(REPO_ROOT / "release-notes")
            )
            for label, pattern in CLASSES:
                if label in ("copyleft", "agpl") and license_naming_ok:
                    continue
                if pattern.search(line):
                    hits.append((path, lineno, label, line.strip()))

    for path, lineno, label, line in hits:
        print(f"{path.relative_to(REPO_ROOT)}:{lineno}: [{label}] {line[:100]}")
    print(
        f"provenance sweep: {len(hits)} hit(s) over {scanned} files "
        f"in {'/, '.join(SCAN_DIRS)}/ (gate: 0)"
    )
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
