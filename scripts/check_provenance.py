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
implementation derived from Lieske 1998 / Meeus ch. 44. All classes must
produce zero hits over ``libephemeris/**/*.py``; this gate keeps the tree
that way.

Usage:
    python scripts/check_provenance.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PACKAGE_ROOT = REPO_ROOT / "libephemeris"

SOURCE_FILE_RE = re.compile(
    r"swehouse|swecl\.c|sweph\.c|swemmoon|swemplan|swedate\.c", re.IGNORECASE
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
CLASSES = (
    ("source-file-ref", SOURCE_FILE_RE),
    ("identifier", IDENTIFIER_RE),
    ("pymeeus", PYMEEUS_RE),
)


def main() -> int:
    hits: list[tuple[Path, int, str, str]] = []
    files = sorted(
        p for p in PACKAGE_ROOT.rglob("*.py") if "__pycache__" not in p.parts
    )
    for path in files:
        for lineno, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(), start=1
        ):
            for label, pattern in CLASSES:
                if pattern.search(line):
                    hits.append((path, lineno, label, line.strip()))

    for path, lineno, label, line in hits:
        print(f"{path.relative_to(REPO_ROOT)}:{lineno}: [{label}] {line[:100]}")
    print(
        f"provenance sweep: {len(hits)} hit(s) over {len(files)} files "
        f"in libephemeris/ (gate: 0)"
    )
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
