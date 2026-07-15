#!/usr/bin/env python3
"""Build LibEphemeris manuals in EPUB and PDF format.

Requirements (install via Homebrew):
    brew install pandoc tectonic

Usage:
    python scripts/build_manual.py                    # Build all (EPUB + PDF, IT + EN)
    python scripts/build_manual.py --lang it           # Italian only
    python scripts/build_manual.py --lang en           # English only
    python scripts/build_manual.py --format epub       # EPUB only
    python scripts/build_manual.py --format pdf        # PDF only
    python scripts/build_manual.py --lang it --format epub  # Italian EPUB only

Output goes to docs/build/.

Provenance:
    Project-authored documentation build orchestration over the repository's
    reviewed Markdown sources. Pandoc/Tectonic perform presentation conversion;
    this script defines no scientific content, numerical model, or data source.
    Generated manuals inherit the citations and licensing boundaries of their
    source documents.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DOCS_DIR = PROJECT_ROOT / "docs"
BUILD_DIR = DOCS_DIR / "build"

# Emoji-to-text mapping for PDF (LaTeX engines struggle with emoji)
EMOJI_MAP = {
    "\U0001f4d0": "[Theory]",  # 📐
    "\U0001f30d": "[Real life]",  # 🌍
    "\U0001f4bb": "[Code]",  # 💻
    "\u2b50": "[Star]",  # ⭐
    "\U0001f52d": "[Telescope]",  # 🔭
    "\u26a0\ufe0f": "[Warning]",  # ⚠️
    "\u26a0": "[Warning]",  # ⚠ (without variation selector)
    "\U0001f4a1": "[Tip]",  # 💡
    "\U0001f4dd": "[Note]",  # 📝
}

MANUALS = {
    "it": {
        "src": DOCS_DIR / "manual",
        "metadata": DOCS_DIR / "manual" / "metadata.yaml",
        "epub": BUILD_DIR / "libephemeris-manual-it.epub",
        "pdf": BUILD_DIR / "libephemeris-manual-it.pdf",
    },
    "en": {
        "src": DOCS_DIR / "manual-eng",
        "metadata": DOCS_DIR / "manual-eng" / "metadata.yaml",
        "epub": BUILD_DIR / "libephemeris-manual-en.epub",
        "pdf": BUILD_DIR / "libephemeris-manual-en.pdf",
    },
}


def check_dependencies() -> None:
    """Verify that pandoc and tectonic are installed."""
    missing = []
    for tool in ("pandoc", "tectonic"):
        if shutil.which(tool) is None:
            missing.append(tool)
    if missing:
        print(f"ERROR: Missing required tools: {', '.join(missing)}")
        print("Install them with: brew install pandoc tectonic")
        sys.exit(1)


def get_md_files(src_dir: Path) -> list[Path]:
    """Get all markdown files sorted by filename."""
    return sorted(src_dir.glob("*.md"))


def strip_emoji_for_pdf(text: str) -> str:
    """Replace emoji characters with text equivalents for LaTeX compatibility."""
    for emoji, replacement in EMOJI_MAP.items():
        text = text.replace(emoji, replacement)
    # Remove any remaining emoji (broad Unicode ranges)
    # Covers Miscellaneous Symbols, Dingbats, Emoticons, Transport, Supplemental, etc.
    text = re.sub(
        r"[\U0001f600-\U0001f64f"  # Emoticons
        r"\U0001f300-\U0001f5ff"  # Misc Symbols & Pictographs
        r"\U0001f680-\U0001f6ff"  # Transport & Map
        r"\U0001f1e0-\U0001f1ff"  # Flags
        r"\U00002702-\U000027b0"  # Dingbats
        r"\U0000fe00-\U0000fe0f"  # Variation Selectors
        r"\U0001fa00-\U0001fa6f"  # Chess Symbols
        r"\U0001fa70-\U0001faff"  # Symbols Extended-A
        r"\U00002600-\U000026ff"  # Misc Symbols
        r"]+",
        "",
        text,
    )
    return text


def prepare_pdf_sources(md_files: list[Path], tmp_dir: Path) -> list[Path]:
    """Copy markdown files to temp dir with emoji stripped for PDF build."""
    pdf_files = []
    for md_file in md_files:
        content = md_file.read_text(encoding="utf-8")
        cleaned = strip_emoji_for_pdf(content)
        out_path = tmp_dir / md_file.name
        out_path.write_text(cleaned, encoding="utf-8")
        pdf_files.append(out_path)
    return pdf_files


def build_epub(lang: str) -> bool:
    """Build EPUB for the given language. Returns True on success."""
    config = MANUALS[lang]
    md_files = get_md_files(config["src"])

    if not md_files:
        print(f"  WARNING: No markdown files found in {config['src']}")
        return False

    cmd = [
        "pandoc",
        "--metadata-file",
        str(config["metadata"]),
        "--toc",
        "--toc-depth=2",
        "--epub-chapter-level=1",
        "--highlight-style=tango",
        "--wrap=none",
        "-o",
        str(config["epub"]),
    ]
    cmd.extend(str(f) for f in md_files)

    print(f"  Building EPUB ({lang})...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: pandoc failed:\n{result.stderr}")
        return False

    size_kb = config["epub"].stat().st_size / 1024
    print(f"  -> {config['epub'].relative_to(PROJECT_ROOT)} ({size_kb:.0f} KB)")
    return True


def build_pdf(lang: str) -> bool:
    """Build PDF for the given language. Returns True on success."""
    config = MANUALS[lang]
    md_files = get_md_files(config["src"])

    if not md_files:
        print(f"  WARNING: No markdown files found in {config['src']}")
        return False

    # Prepare emoji-free copies for LaTeX
    with tempfile.TemporaryDirectory(prefix="libeph-pdf-") as tmp_dir:
        pdf_files = prepare_pdf_sources(md_files, Path(tmp_dir))

        cmd = [
            "pandoc",
            "--metadata-file",
            str(config["metadata"]),
            "--toc",
            "--toc-depth=2",
            "--pdf-engine=tectonic",
            "--highlight-style=tango",
            "-V",
            "geometry:margin=1in",
            "-V",
            "fontsize=11pt",
            "-V",
            "documentclass=report",
            "-V",
            "colorlinks=true",
            "-V",
            "linkcolor=NavyBlue",
            "-V",
            "urlcolor=NavyBlue",
            "-V",
            "toccolor=black",
            "-V",
            "header-includes=\\usepackage{fvextra}\\DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,commandchars=\\\\\\{\\}}",
            "-o",
            str(config["pdf"]),
        ]
        cmd.extend(str(f) for f in pdf_files)

        print(f"  Building PDF ({lang})...")
        print("  (first run may download LaTeX packages — this is normal)")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  ERROR: pandoc/tectonic failed:\n{result.stderr}")
            return False

    size_kb = config["pdf"].stat().st_size / 1024
    print(f"  -> {config['pdf'].relative_to(PROJECT_ROOT)} ({size_kb:.0f} KB)")
    return True


def main() -> None:
    """Build selected language/format variants of the project manual.

    The CLI validates external document-tool dependencies, creates the build
    directory, and dispatches the deterministic EPUB/PDF assembly pipeline for
    Italian, English, or both. A failed requested artifact causes a non-zero
    process exit; source Markdown files are never modified.
    """
    parser = argparse.ArgumentParser(
        description="Build LibEphemeris manuals (EPUB and PDF)"
    )
    parser.add_argument(
        "--lang",
        choices=["it", "en", "all"],
        default="all",
        help="Language to build (default: all)",
    )
    parser.add_argument(
        "--format",
        choices=["epub", "pdf", "all"],
        default="all",
        dest="fmt",
        help="Output format (default: all)",
    )
    args = parser.parse_args()

    check_dependencies()
    BUILD_DIR.mkdir(parents=True, exist_ok=True)

    langs = ["it", "en"] if args.lang == "all" else [args.lang]
    success = True

    for lang in langs:
        lang_name = "Italian" if lang == "it" else "English"
        print(f"\n{'=' * 50}")
        print(f"  {lang_name} manual")
        print(f"{'=' * 50}")

        if args.fmt in ("epub", "all"):
            if not build_epub(lang):
                success = False

        if args.fmt in ("pdf", "all"):
            if not build_pdf(lang):
                success = False

    print()
    if success:
        print(f"All builds completed. Output in {BUILD_DIR.relative_to(PROJECT_ROOT)}/")
    else:
        print("Some builds failed — see errors above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
