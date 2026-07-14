#!/usr/bin/env python3
"""Check (or stamp) SPDX license headers on shipped source files.

Every Python file under ``libephemeris/`` must carry an SPDX identifier in
its first lines.  Owned files are AGPL-3.0-only licensed; vendored/adapted files
keep their upstream identifier (see ``EXCEPTIONS`` and THIRD_PARTY_NOTICES.md).

Usage:
    python scripts/check_spdx_headers.py --check   # CI gate (default)
    python scripts/check_spdx_headers.py --fix     # stamp missing headers

Exit status: 0 if every file carries the expected identifier, 1 otherwise
(--fix only fails on identifier mismatches it refuses to rewrite).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PACKAGE_ROOT = REPO_ROOT / "libephemeris"

PROJECT_LICENSE = "AGPL-3.0-only"
COPYRIGHT_LINE = "# Copyright (c) 2025-2026 Giacomo Battaglia"

# Vendored or adapted third-party code keeps its upstream license and must
# NOT be stamped with the project license.  Paths relative to the repo root.
EXCEPTIONS: dict[str, str] = {
    "libephemeris/vendor/spktype21.py": "MIT",
    "libephemeris/moon_theories/tass17.py": "MIT",
    "libephemeris/moon_theories/tass17_data.py": "MIT",
}
THIRD_PARTY_COPYRIGHT_LINES: dict[str, str] = {
    "libephemeris/vendor/spktype21.py": (
        "# Copyright (c) 2018 The Python Packaging Authority"
    ),
    "libephemeris/moon_theories/tass17.py": ("# Copyright (c) 2005 Johannes Gajdosik"),
    "libephemeris/moon_theories/tass17_data.py": (
        "# Copyright (c) 2005 Johannes Gajdosik"
    ),
}
THIRD_PARTY_NOTE = (
    "# Third-party/adapted code — see file docstring and THIRD_PARTY_NOTICES.md"
)

SPDX_RE = re.compile(r"#\s*SPDX-License-Identifier:\s*(?P<expr>.+?)\s*$")
CODING_RE = re.compile(r"^#.*coding[:=]\s*[-\w.]+")
HEAD_LINES = 6  # the identifier must appear within the first lines


def expected_expression(path: Path) -> str:
    rel = path.relative_to(REPO_ROOT).as_posix()
    return EXCEPTIONS.get(rel, PROJECT_LICENSE)


def find_spdx(lines: list[str]) -> str | None:
    for line in lines[:HEAD_LINES]:
        match = SPDX_RE.match(line.strip())
        if match:
            return match.group("expr").strip()
    return None


def header_lines(path: Path) -> list[str]:
    expr = expected_expression(path)
    lines = [f"# SPDX-License-Identifier: {expr}"]
    if expr == PROJECT_LICENSE:
        lines.append(COPYRIGHT_LINE)
    else:
        rel = path.relative_to(REPO_ROOT).as_posix()
        lines.append(THIRD_PARTY_COPYRIGHT_LINES[rel])
        lines.append(THIRD_PARTY_NOTE)
    return lines


def insertion_index(lines: list[str]) -> int:
    """Index at which the header goes: after shebang and PEP 263 coding line.

    A coding declaration is only honored on line 1 or 2, so the header must
    never push it further down.
    """
    idx = 0
    if lines and lines[0].startswith("#!"):
        idx = 1
    for candidate in range(idx, min(2, len(lines))):
        if CODING_RE.match(lines[candidate]):
            idx = candidate + 1
    return idx


def stamp(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines(keepends=True)
    idx = insertion_index([ln.rstrip("\n") for ln in lines])
    block = "".join(line + "\n" for line in header_lines(path))
    new_text = "".join(lines[:idx]) + block + "".join(lines[idx:])
    path.write_text(new_text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true", help="report problems (default)")
    group.add_argument("--fix", action="store_true", help="stamp missing headers")
    args = parser.parse_args()

    missing: list[Path] = []
    mismatched: list[tuple[Path, str, str]] = []
    missing_copyrights: list[tuple[Path, str]] = []
    files = sorted(
        p for p in PACKAGE_ROOT.rglob("*.py") if "__pycache__" not in p.parts
    )
    for path in files:
        head = path.read_text(encoding="utf-8").splitlines()[:HEAD_LINES]
        found = find_spdx(head)
        expected = expected_expression(path)
        rel = path.relative_to(REPO_ROOT).as_posix()
        expected_copyright = THIRD_PARTY_COPYRIGHT_LINES.get(rel, COPYRIGHT_LINE)
        if found is None:
            missing.append(path)
        elif found != expected:
            mismatched.append((path, found, expected))
        if expected_copyright not in head:
            missing_copyrights.append((path, expected_copyright))

    if args.fix:
        stamped = set(missing)
        for path in missing:
            stamp(path)
            print(f"stamped: {path.relative_to(REPO_ROOT)}")
        missing = []
        missing_copyrights = [
            item for item in missing_copyrights if item[0] not in stamped
        ]

    for path in missing:
        print(f"MISSING SPDX header: {path.relative_to(REPO_ROOT)}")
    for path, found, expected in mismatched:
        print(
            f"MISMATCHED SPDX header: {path.relative_to(REPO_ROOT)}\n"
            f"  found:    {found}\n  expected: {expected}"
        )
    for path, expected in missing_copyrights:
        print(
            f"MISSING copyright notice: {path.relative_to(REPO_ROOT)}\n"
            f"  expected: {expected}"
        )

    problems = len(missing) + len(mismatched) + len(missing_copyrights)
    print(f"{len(files)} files scanned, {problems} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
