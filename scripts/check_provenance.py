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
SCAN_SUFFIXES = (".py", ".md")

# Files that legitimately record the retired SE / PyMeeus / LGPL history
# (the gate itself and the provenance evidence docs). Never includes any
# libephemeris/ file — shipped code stays strictly zero-hit.
ALLOWLIST = frozenset(
    {
        "scripts/check_provenance.py",
        "docs/methodology/provenance-sweep-2026-06.md",
        "docs/methodology/galilean-clean-room-2026-06.md",
        "docs/methodology/galilean-e5-spec.md",
    }
)

SOURCE_FILE_RE = re.compile(
    r"swehouse|swecl\.c|sweph\.c|swephlib|swejpl|swehel"
    r"|swemmoon|swemplan|swedate\.c",
    re.IGNORECASE,
)
IDENTIFIER_RE = re.compile(
    r"\bswed\b|\bdgsect\b|\bxs1\b|\bxh1\b|\bfh1\b|\bmdd\b|\bmdn\b|\badp\b"
    r"|\badmc\b|\bsamc\b|\bdfac\b|apc_sector|\bxeq0\b|\bxp0\b"
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
CLASSES = (
    ("source-file-ref", SOURCE_FILE_RE),
    ("identifier", IDENTIFIER_RE),
    ("pymeeus", PYMEEUS_RE),
    ("copyleft", COPYLEFT_RE),
)


def main() -> int:
    hits: list[tuple[Path, int, str, str]] = []
    files = sorted(
        p
        for d in SCAN_DIRS
        for p in (REPO_ROOT / d).rglob("*")
        if p.suffix in SCAN_SUFFIXES and "__pycache__" not in p.parts
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
