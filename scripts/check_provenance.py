#!/usr/bin/env python3
"""Third-party-source provenance sweep (zero-hit CI gate).

Re-runs the fingerprint greps recorded in
``docs/methodology/provenance-sweep-2026-06.md`` and
``docs/methodology/galilean-clean-room-2026-06.md`` over the shipped
package:

  class 1 — SE source-file references (swehouse, swecl.c, sweph.c, ...)
  class 2 — SE implementation identifiers (swed, dgsect, xs1, mdd, ...)
  class 3 — PyMeeus implementation identifiers (JupiterMoons,
            rectangular_positions_jovian_equatorial, ...)

Classes 1-2 come from the June 2026 Swiss Ephemeris remediation (WS1);
class 3 from the June 2026 Galilean clean-room rewrite, which replaced the
PyMeeus-adapted ``moon_theories/galilean.py`` with an independent
implementation derived from Lieske 1998 / Meeus ch. 44; class 4 flags any
copyleft (L/GPL) license declaration that strays into the tree.

The sweep covers ``libephemeris/``, ``scripts/`` and ``docs/`` (``*.py``
and ``*.md``). A short allowlist exempts the files that legitimately record
the retired Swiss Ephemeris / PyMeeus history (this script and the
methodology docs); ``libephemeris/`` is never exempt. All classes must
produce zero hits; this gate keeps the tree that way.

Usage:
    python scripts/check_provenance.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCAN_DIRS = ("libephemeris", "scripts", "docs")
# Root-level files that legally represent the package (long_description,
# notices) are scanned individually in addition to the directories.
ROOT_FILES = (
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

# Files that legitimately record the retired SE / PyMeeus / LGPL history
# (the gate itself and the provenance evidence docs). Never includes any
# libephemeris/ file — shipped code stays strictly zero-hit.
ALLOWLIST = frozenset(
    {
        "scripts/check_provenance.py",
        "docs/methodology/provenance-sweep-2026-06.md",
        "docs/methodology/galilean-clean-room-2026-06.md",
        "docs/methodology/galilean-e5-spec.md",
        "docs/methodology/independence-remediation-2026-07.md",
        # Root legal/provenance documents: they exist to name licenses and
        # record the retired-source history, so every class legitimately
        # appears in them.
        "NOTICE.md",
        "LICENSING.md",
        "THIRD_PARTY_NOTICES.md",
        "CHANGELOG.md",
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
# allowlisted provenance/remediation records.
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
    r"^(seorbel|sefstars|seleapsec|sedeltat)(\.|$)|\.se1$",
    re.IGNORECASE,
)
DATA_DIR = "libephemeris/data"


def main() -> int:
    hits: list[tuple[Path, int, str, str]] = []

    # Filename gate: no reference-distribution data file among the files
    # tracked by git (the gitignored local oracle data used by the
    # validation tooling is dev-only and never shipped; anything TRACKED
    # under such a name is a hard failure).
    import subprocess

    tracked = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "ls-files", "-z"],
        capture_output=True,
        text=True,
        check=False,
    ).stdout.split("\0")
    for rel in tracked:
        if not rel:
            continue
        name = rel.rsplit("/", 1)[-1]
        if FOREIGN_DATA_NAME_RE.search(name):
            hits.append((REPO_ROOT / rel, 0, "foreign-data-file", name))

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
            # An explicit inline waiver for individually-reviewed lines
            # (e.g. install docs that must name an optional extra's GPL
            # license). Never valid inside libephemeris/.
            if "provenance-ok" in line and not path.is_relative_to(
                REPO_ROOT / "libephemeris"
            ):
                continue
            for label, pattern in CLASSES:
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
