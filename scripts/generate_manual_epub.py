#!/usr/bin/env python3
"""Generate LibEphemeris manual EPUB files without pandoc.

Uses ebooklib + markdown for Kobo-compatible EPUB output with proper
chapter splitting, NCX/NAV navigation, and XHTML-safe content.

Requirements (auto-installed):
    pip install ebooklib markdown pyyaml

Usage:
    python scripts/generate_manual_epub.py                  # Both languages
    python scripts/generate_manual_epub.py --lang it        # Italian only
    python scripts/generate_manual_epub.py --lang en        # English only

Output goes to docs/build/.

Provenance:
    Project-authored EPUB packaging over reviewed repository Markdown, using
    the documented ebooklib/Markdown dependencies. Chapter splitting, XHTML
    escaping, and navigation are presentation choices only. This script adds no
    astronomical algorithm, coefficient, or source claim.
"""

from __future__ import annotations

import argparse
import html as html_module
import re
import sys
from pathlib import Path

try:
    import markdown
    import yaml
    from ebooklib import epub
except ImportError:
    print("ERROR: Missing dependencies. Install them with:")
    print("  pip install ebooklib markdown pyyaml")
    sys.exit(1)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DOCS_DIR = PROJECT_ROOT / "docs"
BUILD_DIR = DOCS_DIR / "build"

MANUALS = {
    "it": {
        "src": DOCS_DIR / "manual",
        "out": BUILD_DIR / "libephemeris-manual-it.epub",
        "identifier": "libephemeris-manual-it-2026",
    },
    "en": {
        "src": DOCS_DIR / "manual-eng",
        "out": BUILD_DIR / "libephemeris-manual-en.epub",
        "identifier": "libephemeris-manual-en-2026",
    },
}

# --- CSS optimized for Kobo e-readers and technical manuals ---
STYLE_CSS = """\
body {
    font-family: "Georgia", "Times New Roman", serif;
    line-height: 1.6;
    margin: 5%;
    text-align: justify;
}
h1 {
    font-size: 1.8em;
    text-align: center;
    margin-top: 2em;
    margin-bottom: 0.8em;
    page-break-before: always;
}
h2 {
    font-size: 1.3em;
    margin-top: 1.5em;
    margin-bottom: 0.6em;
    border-bottom: 1px solid #999;
    padding-bottom: 0.2em;
}
h3 {
    font-size: 1.1em;
    margin-top: 1.2em;
    margin-bottom: 0.4em;
}
p {
    text-indent: 0;
    margin-top: 0.5em;
    margin-bottom: 0.5em;
}
em { font-style: italic; }
strong { font-weight: bold; }
code {
    font-family: "Courier New", "Courier", monospace;
    font-size: 0.85em;
    background-color: #f4f4f4;
    padding: 0.1em 0.3em;
    border-radius: 3px;
}
pre {
    font-family: "Courier New", "Courier", monospace;
    font-size: 0.8em;
    background-color: #f4f4f4;
    border: 1px solid #ddd;
    border-radius: 4px;
    padding: 0.8em;
    overflow-x: auto;
    white-space: pre-wrap;
    word-wrap: break-word;
    line-height: 1.4;
    margin: 1em 0;
}
pre code {
    background-color: transparent;
    padding: 0;
    font-size: 1em;
}
blockquote {
    margin: 1em 1.5em;
    padding: 0.5em 1em;
    border-left: 3px solid #999;
    font-style: italic;
    background-color: #f9f9f9;
}
table {
    border-collapse: collapse;
    margin: 1em 0;
    width: 100%;
    font-size: 0.9em;
}
th, td {
    border: 1px solid #ccc;
    padding: 0.4em 0.6em;
    text-align: left;
}
th {
    background-color: #f0f0f0;
    font-weight: bold;
}
tr:nth-child(even) {
    background-color: #f9f9f9;
}
ul, ol {
    margin-left: 1.5em;
    margin-top: 0.5em;
    margin-bottom: 0.5em;
}
li {
    margin-bottom: 0.3em;
}
hr {
    border: none;
    border-top: 1px solid #ccc;
    margin: 2em 0;
}
.title-page {
    text-align: center;
    margin-top: 25%;
}
.title-page h1 {
    font-size: 2.2em;
    page-break-before: avoid;
    border-bottom: none;
}
.title-page p {
    text-align: center;
    font-size: 1.1em;
}
"""


def load_metadata(src_dir: Path) -> dict:
    """Load metadata.yaml from a manual source directory."""
    meta_path = src_dir / "metadata.yaml"
    raw = meta_path.read_text(encoding="utf-8").strip().strip("---").strip()
    return yaml.safe_load(raw)


def get_md_files(src_dir: Path) -> list[Path]:
    """Get all markdown files sorted by filename."""
    return sorted(src_dir.glob("*.md"))


def md_to_html(md_text: str) -> str:
    """Convert markdown to XHTML-safe HTML."""
    raw_html = markdown.markdown(
        md_text,
        extensions=["extra", "codehilite", "tables", "fenced_code"],
        extension_configs={"codehilite": {"css_class": "code", "guess_lang": False}},
    )
    # Convert named HTML entities to UTF-8 characters, but preserve
    # the 5 entities that are valid (and required) in XHTML:
    #   &amp; &lt; &gt; &apos; &quot;
    # html.unescape() would convert &lt; to <, breaking code blocks.
    XML_VALID = {"amp", "lt", "gt", "apos", "quot"}

    def _replace_entity(m: re.Match) -> str:
        name = m.group(1)
        if name in XML_VALID:
            return m.group(0)  # keep as-is
        # Convert to character via html.unescape
        return html_module.unescape(m.group(0))

    return re.sub(r"&([a-zA-Z]+);", _replace_entity, raw_html)


def extract_title(md_text: str) -> str | None:
    """Extract the h1 title from markdown content."""
    match = re.match(r"^#\s+(.+)$", md_text, re.MULTILINE)
    return match.group(1).strip() if match else None


def build_epub(lang: str) -> bool:
    """Build EPUB for the given language. Returns True on success."""
    config = MANUALS[lang]
    src_dir = config["src"]

    if not src_dir.exists():
        print(f"  WARNING: Source directory not found: {src_dir}")
        return False

    metadata = load_metadata(src_dir)
    md_files = get_md_files(src_dir)

    if not md_files:
        print(f"  WARNING: No markdown files found in {src_dir}")
        return False

    book = epub.EpubBook()

    # --- Metadata ---
    book.set_identifier(config["identifier"])
    book.set_title(metadata.get("title", "LibEphemeris"))
    book.set_language(metadata.get("lang", "en"))
    book.add_author(metadata.get("author", "Unknown"))

    for field in ("subtitle", "description"):
        value = metadata.get(field, "")
        if value:
            book.add_metadata("DC", "description", value.strip())

    rights = metadata.get("rights", "")
    if rights:
        book.add_metadata("DC", "rights", rights)

    book.add_metadata("DC", "subject", "Astronomy")
    book.add_metadata("DC", "subject", "Python")
    book.add_metadata("DC", "subject", "Programming")

    # --- CSS ---
    style = epub.EpubItem(
        uid="style",
        file_name="style/default.css",
        media_type="text/css",
        content=STYLE_CSS.encode("utf-8"),
    )
    book.add_item(style)

    # --- Title page ---
    title_html = (
        '<?xml version="1.0" encoding="utf-8"?>\n'
        "<!DOCTYPE html>\n"
        '<html xmlns="http://www.w3.org/1999/xhtml" '
        f'xml:lang="{metadata.get("lang", "en")}" '
        f'lang="{metadata.get("lang", "en")}">\n'
        "<head>\n"
        f"  <title>{metadata.get('title', '')}</title>\n"
        '  <link rel="stylesheet" type="text/css" '
        'href="style/default.css"/>\n'
        "</head>\n"
        "<body>\n"
        '<div class="title-page">\n'
        f"  <h1>{metadata.get('title', '')}</h1>\n"
        f"  <p><em>{metadata.get('subtitle', '')}</em></p>\n"
        "  <p>&#160;</p>\n"
        f"  <p><strong>{metadata.get('author', '')}</strong></p>\n"
        f"  <p>{metadata.get('date', '')}</p>\n"
        "  <p>&#160;</p>\n"
        f"  <p><small>License: {metadata.get('rights', '')}</small></p>\n"
        "</div>\n"
        "</body>\n"
        "</html>"
    )

    title_chap = epub.EpubItem(
        uid="titlepage",
        file_name="titlepage.xhtml",
        media_type="application/xhtml+xml",
        content=title_html.encode("utf-8"),
    )
    book.add_item(title_chap)

    # --- Chapters ---
    toc_entries: list[epub.Link] = []
    spine: list = ["nav", title_chap]
    lang_code = metadata.get("lang", "en").split("-")[0]

    for i, md_file in enumerate(md_files):
        md_text = md_file.read_text(encoding="utf-8")
        title = extract_title(md_text) or md_file.stem
        html_body = md_to_html(md_text)
        file_name = f"chap_{i:02d}.xhtml"

        xhtml = (
            '<?xml version="1.0" encoding="utf-8"?>\n'
            "<!DOCTYPE html>\n"
            '<html xmlns="http://www.w3.org/1999/xhtml" '
            f'xml:lang="{lang_code}" lang="{lang_code}">\n'
            "<head>\n"
            f"  <title>{title}</title>\n"
            '  <link rel="stylesheet" type="text/css" '
            'href="style/default.css"/>\n'
            "</head>\n"
            f"<body>\n{html_body}\n</body>\n"
            "</html>"
        )

        chap = epub.EpubItem(
            uid=f"chap_{i:02d}",
            file_name=file_name,
            media_type="application/xhtml+xml",
            content=xhtml.encode("utf-8"),
        )
        book.add_item(chap)
        spine.append(chap)
        toc_entries.append(epub.Link(file_name, title, f"chap_{i:02d}"))

    # --- Navigation ---
    book.toc = toc_entries
    book.add_item(epub.EpubNcx())
    book.add_item(epub.EpubNav())
    book.spine = spine

    # --- Write ---
    config["out"].parent.mkdir(parents=True, exist_ok=True)
    epub.write_epub(str(config["out"]), book, {})

    size_kb = config["out"].stat().st_size / 1024
    rel_path = config["out"].relative_to(PROJECT_ROOT)
    print(f"  -> {rel_path} ({size_kb:.0f} KB, {len(toc_entries)} chapters)")
    return True


def main() -> None:
    """Generate the selected language editions of the project EPUB manual.

    This entry point performs document packaging only.  It selects reviewed
    Markdown inputs and delegates EPUB construction without changing any
    scientific formula, table, or runtime value.
    """
    parser = argparse.ArgumentParser(
        description="Generate LibEphemeris manual EPUB files (no pandoc required)"
    )
    parser.add_argument(
        "--lang",
        choices=["it", "en", "all"],
        default="all",
        help="Language to build (default: all)",
    )
    args = parser.parse_args()

    langs = ["it", "en"] if args.lang == "all" else [args.lang]
    success = True

    for lang in langs:
        lang_name = {"it": "Italian", "en": "English"}[lang]
        print(f"\n{'=' * 50}")
        print(f"  {lang_name} manual (EPUB)")
        print(f"{'=' * 50}")

        if not build_epub(lang):
            success = False

    print()
    if success:
        print(f"Done. Output in {BUILD_DIR.relative_to(PROJECT_ROOT)}/")
    else:
        print("Some builds failed — see errors above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
