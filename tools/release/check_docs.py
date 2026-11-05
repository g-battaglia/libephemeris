#!/usr/bin/env python3
"""Check maintained Markdown links locally, without fetching external sources.

Only documentation, not comparison archives, is traversed. Local destinations
must remain inside the source tree; generated HTML and code fences are ignored.
This catches relocation mistakes without crawling prohibited third-party pages.
"""

from __future__ import annotations

import argparse
import html
import os
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[2]
ROOT_DOCUMENTS = (
    "README.md", "RELEASE_NOTES.md", "MIGRATION.md", "CHANGELOG.md", "CONTRIBUTING.md",
    "CLA.md", "SECURITY.md", "THIRD_PARTY_NOTICES.md",
    "core/README.md", "core/data/README.md", "core/data/leb2/README.md",
    "python/README.md", "python/libephemeris/data/README.md",
    "tools/leb_generation/README.md",
)
LINK = re.compile(r"!?\[[^\]]*\]\(([^\s)]+)(?:\s+\"[^\"]*\")?\)")


def prose(text: str) -> str:
    """Remove fenced code while preserving ordinary Markdown heading/link text."""
    lines = []
    fence = None
    for line in text.splitlines():
        marker = re.match(r"^\s*(`{3,}|~{3,})", line)
        if marker:
            token = marker.group(1)
            if fence is None:
                fence = token
            elif token[0] == fence[0] and len(token) >= len(fence):
                fence = None
            continue
        if fence is None:
            lines.append(line)
    return "\n".join(lines)


def anchors(text: str) -> set[str]:
    """Recognize explicit HTML IDs and GitHub-style ATX heading fragments."""
    result = set(re.findall(r"\bid=[\"']([^\"']+)[\"']", text))
    counts: dict[str, int] = {}
    for line in prose(text).splitlines():
        heading = re.match(r"^#{1,6}\s+(.+?)(?:\s+#+)?\s*$", line)
        if not heading:
            continue
        # Link labels contribute their visible text, not their destination.
        visible = re.sub(r"\[([^\]]+)\]\([^)]*\)", r"\1", heading.group(1))
        visible = re.sub(r"<[^>]+>", "", visible)
        slug = re.sub(r"[^\w\- ]", "", html.unescape(visible).lower()).replace(" ", "-")
        count = counts.get(slug, 0)
        counts[slug] = count + 1
        result.add(slug if count == 0 else f"{slug}-{count}")
    return result


def document_errors(path: Path, root: Path) -> list[str]:
    """Validate local paths/Markdown fragments; never open an external URL."""
    errors = []
    text = prose(path.read_text(encoding="utf-8"))
    if "https://github.com/g-battaglia/libephemeris-rust" in text:
        errors.append(f"{path.relative_to(root)}: depends on retired migration repository")
    for match in LINK.finditer(text):
        target = match.group(1)
        parsed = urlsplit(target)
        if parsed.scheme or parsed.netloc:
            continue
        destination = ((path.parent / unquote(parsed.path)) if parsed.path else path).resolve()
        if not destination.is_relative_to(root.resolve()):
            errors.append(f"{path.relative_to(root)}: link escapes source tree: {target}")
        elif destination.name in {"AGENTS.md", "PLAN.md", "PROGRESS.md"}:
            errors.append(f"{path.relative_to(root)}: public documentation links to process file: {target}")
        elif not destination.exists():
            errors.append(f"{path.relative_to(root)}: missing destination: {target}")
        elif destination.suffix == ".md" and parsed.fragment:
            fragment = unquote(parsed.fragment)
            if fragment not in anchors(destination.read_text(encoding="utf-8")):
                errors.append(f"{path.relative_to(root)}: missing anchor: {target}")
    return errors


def documents(root: Path) -> list[Path]:
    """List maintained Markdown, excluding generated and process archives."""
    paths = [root / name for name in ROOT_DOCUMENTS if (root / name).is_file()]
    # Traverse only maintained books, never sibling worktrees or generated HTML.
    for book in (root / "docs", root / "beginners-manual"):
        paths += [path for path in book.rglob("*.md")
                  if not {"_build", "archive"}.intersection(path.relative_to(book).parts)]
    return sorted(paths)


def relocate_links(text: str, document: Path, root: Path, moves: list[tuple[str, str]]) -> str:
    """Render local links after explicit source moves, preserving all prose.

    More-specific mappings win, so Python examples can move independently of
    Rust examples. Destinations must exist inside the project. Code fences,
    external URLs, labels, fragments and non-link content remain unchanged.
    """
    mappings = sorted(((root / old, root / new) for old, new in moves),
                      key=lambda pair: len(pair[0].parts), reverse=True)

    def replace(match):
        target = match.group(1)
        parsed = urlsplit(target)
        if parsed.scheme or parsed.netloc or not parsed.path:
            return match.group(0)
        original = (document.parent / unquote(parsed.path)).resolve()
        for old, new in mappings:
            if original.is_relative_to(old):
                destination = new / original.relative_to(old)
                if not destination.is_relative_to(root) or not destination.exists():
                    raise ValueError(f"missing relocated destination: {destination}")
                relative = Path(os.path.relpath(destination, document.parent)).as_posix()
                rendered = parsed._replace(path=relative).geturl()
                start = match.start(1) - match.start(0)
                return match.group(0)[:start] + rendered + match.group(0)[start + len(target):]
        return match.group(0)

    result = []
    fence = None
    for line in text.splitlines(keepends=True):
        marker = re.match(r"^\s*(`{3,}|~{3,})", line)
        if marker:
            token = marker.group(1)
            if fence is None:
                fence = token
            elif token[0] == fence[0] and len(token) >= len(fence):
                fence = None
            result.append(line)
        else:
            result.append(line if fence else LINK.sub(replace, line))
    return "".join(result)


def check(root: Path) -> list[str]:
    """Return problems without fetching external source pages."""
    return [error for path in documents(root) for error in document_errors(path, root)]


def main() -> int:
    """Print a concise local-link report suitable for manual release checks."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--relocate", nargs=2, action="append", default=[], metavar=("OLD", "NEW"),
                        help="rewrite links after an explicit repository-relative source move")
    arguments = parser.parse_args()
    root = arguments.root.resolve()
    if arguments.relocate:
        for old, new in arguments.relocate:
            if not (root / old).resolve().is_relative_to(root) or not (root / new).resolve().is_relative_to(root):
                parser.error("source moves must stay inside the project")
        for path in documents(root):
            original = path.read_text(encoding="utf-8")
            rendered = relocate_links(original, path, root, arguments.relocate)
            if rendered != original:
                path.write_text(rendered, encoding="utf-8")
    errors = check(root)
    for error in errors:
        print(error)
    if not errors:
        print("Maintained documentation links passed (local paths and Markdown anchors).")
    return bool(errors)


if __name__ == "__main__":
    raise SystemExit(main())
