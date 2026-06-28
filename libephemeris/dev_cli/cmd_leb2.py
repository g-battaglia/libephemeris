# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""LEB2 compressed ephemeris: convert LEB1 -> LEB2, verify accuracy.

Replaces 8 poe tasks: leb2:convert:*, leb2:verify:*.

LEB2 uses error-bounded lossy compression (mantissa truncation + coeff-major
reorder + byte shuffle + zstd) to achieve 4-10x compression while maintaining
<0.001" precision vs LEB1.
"""

from __future__ import annotations

import subprocess
import sys

import click


def _leb2(args: list[str]) -> None:
    """Run the LEB2 generator script."""
    sys.exit(subprocess.call([sys.executable, "scripts/generate_leb2.py", *args]))


# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------


@click.group(
    "leb2",
    short_help='Convert LEB1 to LEB2 compressed format and verify precision (<0.001").',
    help="Convert LEB1 files to LEB2 compressed format and verify precision.\n\n"
    "LEB2 uses error-bounded lossy compression (mantissa truncation + coeff-major\n"
    "reorder + byte shuffle + zstd) to achieve 4-10x smaller files while keeping\n"
    '<0.001" precision vs LEB1. Useful for distributing via PyPI.\n\n'
    "  leph leb2 convert base    # Convert all 4 body groups\n"
    "  leph leb2 verify base     # Verify against LEB1 reference",
)
def leb2_group() -> None:
    """LEB2 compressed format management."""


# ===========================================================================
# leph leb2 convert
# ===========================================================================


@click.group(
    "convert",
    short_help="Convert LEB1 files to LEB2 compressed format (per tier or per group).",
    help="Convert LEB1 files to LEB2 compressed format.\n\nEach tier can be converted as a whole (all 4 body groups: core, asteroids,\napogee, uranians) or one group at a time for finer control.",
)
def convert_group() -> None:
    """LEB2 conversion commands."""


@convert_group.command(
    short_help="Convert base tier LEB1 -> LEB2 (all 4 body groups).",
)
def base() -> None:
    """Convert base tier LEB1 -> LEB2 (all 4 groups: core/asteroids/apogee/uranians)."""
    _leb2(
        [
            "convert-all",
            "data/leb/ephemeris_base.leb",
            "-o",
            "data/leb2/",
            "--tier-name",
            "base",
        ]
    )


@convert_group.command(
    short_help="Convert medium tier LEB1 -> LEB2 (all 4 body groups).",
)
def medium() -> None:
    """Convert medium tier LEB1 -> LEB2 (all 4 groups)."""
    _leb2(
        [
            "convert-all",
            "data/leb/ephemeris_medium.leb",
            "-o",
            "data/leb2/",
            "--tier-name",
            "medium",
        ]
    )


@convert_group.command(
    short_help="Convert extended tier LEB1 -> LEB2 (all 4 body groups).",
)
def extended() -> None:
    """Convert extended tier LEB1 -> LEB2 (all 4 groups)."""
    _leb2(
        [
            "convert-all",
            "data/leb/ephemeris_extended.leb",
            "-o",
            "data/leb2/",
            "--tier-name",
            "extended",
        ]
    )


@convert_group.command(
    "base-core",
    short_help="Convert base tier core group only (14 bodies, ~10.6 MB).",
)
def base_core() -> None:
    """Convert base tier core group only (14 bodies, ~10.6 MB for PyPI)."""
    _leb2(
        [
            "convert",
            "data/leb/ephemeris_base.leb",
            "-o",
            "data/leb2/base_core.leb2",
            "--group",
            "core",
        ]
    )


@convert_group.command(
    "base-asteroids",
    short_help="Convert base tier asteroids group (5 bodies).",
)
def base_asteroids() -> None:
    """Convert base tier asteroids group (Chiron, Ceres, Pallas, Juno, Vesta)."""
    _leb2(
        [
            "convert",
            "data/leb/ephemeris_base.leb",
            "-o",
            "data/leb2/base_asteroids.leb2",
            "--group",
            "asteroids",
        ]
    )


@convert_group.command(
    "base-apogee",
    short_help="Convert base tier apogee group (3 bodies).",
)
def base_apogee() -> None:
    """Convert base tier apogee group (OscuApog, IntpApog, IntpPerig)."""
    _leb2(
        [
            "convert",
            "data/leb/ephemeris_base.leb",
            "-o",
            "data/leb2/base_apogee.leb2",
            "--group",
            "apogee",
        ]
    )


@convert_group.command(
    "base-uranians",
    short_help="Convert base tier uranians group (9 bodies).",
)
def base_uranians() -> None:
    """Convert base tier uranians group (Cupido-Transpluto, 9 bodies)."""
    _leb2(
        [
            "convert",
            "data/leb/ephemeris_base.leb",
            "-o",
            "data/leb2/base_uranians.leb2",
            "--group",
            "uranians",
        ]
    )


leb2_group.add_command(convert_group)


# ===========================================================================
# leph leb2 verify
# ===========================================================================


@click.group(
    "verify",
    short_help='Verify LEB2 precision against LEB1 reference (target: <0.003").',
    help='Verify LEB2 files against LEB1 reference.\n\nSamples random dates and compares LEB2 output to LEB1, reporting\nthe maximum error in arcseconds. Target: <0.003".',
)
def verify_group() -> None:
    """LEB2 verification commands."""


@verify_group.command(
    "base",
    short_help='Verify base tier LEB2 vs LEB1 (<0.003" precision).',
)
def verify_base() -> None:
    """Verify base tier LEB2 core against LEB1 reference (<0.003\" precision)."""
    _leb2(
        [
            "verify",
            "data/leb2/base_core.leb2",
            "--reference",
            "data/leb/ephemeris_base.leb",
            "--samples",
            "500",
        ]
    )


@verify_group.command(
    "medium",
    short_help='Verify medium tier LEB2 vs LEB1 (<0.003" precision).',
)
def verify_medium() -> None:
    """Verify medium tier LEB2 core against LEB1 reference (<0.003\" precision)."""
    _leb2(
        [
            "verify",
            "data/leb2/medium_core.leb2",
            "--reference",
            "data/leb/ephemeris_medium.leb",
            "--samples",
            "500",
        ]
    )


@verify_group.command(
    "extended",
    short_help='Verify extended tier LEB2 vs LEB1 (<0.003" precision).',
)
def verify_extended() -> None:
    """Verify extended tier LEB2 core against LEB1 reference (<0.003\" precision)."""
    _leb2(
        [
            "verify",
            "data/leb2/extended_core.leb2",
            "--reference",
            "data/leb/ephemeris_extended.leb",
            "--samples",
            "500",
        ]
    )


leb2_group.add_command(verify_group)
