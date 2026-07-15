# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Manual/documentation build commands: EPUB, PDF, pandoc and ebooklib workflows.

Replaces 8 poe tasks: manual:build*, docs:manual:generate*.

Provenance:
    Project-authored documentation-build orchestration with no astronomical
    model, scientific data transformation, or runtime coefficient.
"""

from __future__ import annotations

import subprocess
import sys

import click


def _python(args: list[str]) -> None:
    """Run a python script."""
    sys.exit(subprocess.call([sys.executable, *args]))


@click.group(
    "manual",
    short_help="Build user manuals: EPUB and PDF, Italian and English.",
    help="Build the libephemeris user manual in EPUB and/or PDF format.\n\n"
    "Two workflows are available:\n\n"
    "  Pandoc workflow (build*)      Requires: pandoc + tectonic\n"
    "  Ebooklib workflow (generate*) Requires: ebooklib + markdown (no pandoc)\n\n"
    "Both produce manuals in Italian and English. The ebooklib workflow\n"
    "generates Kobo-compatible EPUBs without external tools.\n\n"
    "  leph manual build             # EPUB + PDF, both languages\n"
    "  leph manual generate-epub     # EPUB only, no pandoc needed",
)
def manual_group() -> None:
    """Manual build commands."""


# ===========================================================================
# Pandoc + tectonic workflow (requires: brew install pandoc tectonic)
# ===========================================================================


@manual_group.command(
    short_help="Build all manuals (EPUB + PDF, Italian + English).",
)
def build() -> None:
    """Build all manuals (EPUB + PDF, Italian + English).

    Requires: pandoc and tectonic (brew install pandoc tectonic).
    """
    _python(["scripts/build_manual.py"])


@manual_group.command(
    "build-epub",
    short_help="Build all manuals in EPUB format only (both languages).",
)
def build_epub() -> None:
    """Build all manuals in EPUB format only (Italian + English).

    Requires: pandoc (brew install pandoc).
    """
    _python(["scripts/build_manual.py", "--format", "epub"])


@manual_group.command(
    "build-pdf",
    short_help="Build all manuals in PDF format only (both languages).",
)
def build_pdf() -> None:
    """Build all manuals in PDF format only (Italian + English).

    Requires: pandoc and tectonic (brew install pandoc tectonic).
    """
    _python(["scripts/build_manual.py", "--format", "pdf"])


@manual_group.command(
    "build-it",
    short_help="Build Italian manual (EPUB + PDF).",
)
def build_it() -> None:
    """Build Italian manual (EPUB + PDF).

    Requires: pandoc and tectonic.
    """
    _python(["scripts/build_manual.py", "--lang", "it"])


@manual_group.command(
    "build-en",
    short_help="Build English manual (EPUB + PDF).",
)
def build_en() -> None:
    """Build English manual (EPUB + PDF).

    Requires: pandoc and tectonic.
    """
    _python(["scripts/build_manual.py", "--lang", "en"])


# ===========================================================================
# ebooklib workflow (no pandoc required, Kobo-compatible)
# ===========================================================================


@manual_group.command(
    "generate-epub",
    short_help="Generate all EPUBs (both languages, no pandoc required).",
)
def generate_epub() -> None:
    """Generate all manual EPUBs (Italian + English, no pandoc required).

    Uses ebooklib + markdown. Kobo-compatible output.
    Requires: pip install ebooklib markdown pyyaml
    """
    _python(["scripts/generate_manual_epub.py"])


@manual_group.command(
    "generate-epub-it",
    short_help="Generate Italian manual EPUB (no pandoc required).",
)
def generate_epub_it() -> None:
    """Generate Italian manual EPUB (no pandoc required)."""
    _python(["scripts/generate_manual_epub.py", "--lang", "it"])


@manual_group.command(
    "generate-epub-en",
    short_help="Generate English manual EPUB (no pandoc required).",
)
def generate_epub_en() -> None:
    """Generate English manual EPUB (no pandoc required)."""
    _python(["scripts/generate_manual_epub.py", "--lang", "en"])
