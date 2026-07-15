# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""LEB2 compressed ephemeris: convert LEB1 -> LEB2, verify stored data.

Replaces 8 poe tasks: leb2:convert:*, leb2:verify:*.

LEB2 uses error-bounded lossy compression (mantissa truncation + coeff-major
reorder + byte shuffle + zstd) to achieve 4-10x compression. Verification
reports per-component errors in each body's native stored coordinate units.
"""

from __future__ import annotations

import subprocess
import sys

import click

from ..leb_groups import LEB2_GROUPS

_ALL_GROUPS = f"all {len(LEB2_GROUPS)} body groups: {', '.join(LEB2_GROUPS)}"


def _leb2(args: list[str]) -> None:
    """Run the LEB2 generator script."""
    sys.exit(subprocess.call([sys.executable, "scripts/generate_leb2.py", *args]))


# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------


@click.group(
    "leb2",
    short_help="Convert LEB1 to LEB2 and verify native stored components.",
    help=f"""Convert LEB1 files to LEB2 compressed format and verify agreement.

LEB2 uses error-bounded lossy compression (mantissa truncation,
coefficient-major reorder, byte shuffle, and zstd) to achieve 4-10x smaller
files.
Verification compares native stored components with LEB1: Cartesian ICRS
components use AU; ecliptic components use degrees, degrees, and AU.

\b
  leph leb2 convert base    # Convert {_ALL_GROUPS}
  leph leb2 verify base     # Verify the base core companion
""",
)
def leb2_group() -> None:
    """LEB2 compressed format management."""


# ===========================================================================
# leph leb2 convert
# ===========================================================================


@click.group(
    "convert",
    short_help="Convert LEB1 files to LEB2 compressed format (per tier or per group).",
    help=(
        "Convert LEB1 files to LEB2 compressed format.\n\n"
        f"Each tier can be converted as a whole ({_ALL_GROUPS}) or one group "
        "at a time for finer control."
    ),
)
def convert_group() -> None:
    """LEB2 conversion commands."""


@convert_group.command(
    short_help=f"Convert base tier LEB1 -> LEB2 ({_ALL_GROUPS}).",
)
def base() -> None:
    """Convert base-tier LEB1 into all registered LEB2 groups."""
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
    short_help=f"Convert medium tier LEB1 -> LEB2 ({_ALL_GROUPS}).",
)
def medium() -> None:
    """Convert medium-tier LEB1 into all registered LEB2 groups."""
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
    short_help=f"Convert extended tier LEB1 -> LEB2 ({_ALL_GROUPS}).",
)
def extended() -> None:
    """Convert extended-tier LEB1 into all registered LEB2 groups."""
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
    short_help="Convert base tier core group only (14 bodies, ~10.7 MB).",
)
def base_core() -> None:
    """Convert base tier core group only (14 bodies, ~10.7 MB for PyPI)."""
    _leb2(
        [
            "convert",
            "data/leb/ephemeris_base.leb",
            "-o",
            "data/leb2/base_core.leb2",
            "--group",
            "core",
            "--tier",
            "base",
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
            "--tier",
            "base",
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
            "--tier",
            "base",
        ]
    )


@convert_group.command(
    "base-uranians",
    short_help="Retired pending independently sourced elements.",
)
def base_uranians() -> None:
    """Reject conversion of the retired hypothetical-body group."""
    raise click.ClickException(
        "Uranian LEB conversion is retired pending independently sourced "
        "elements and clean-room regeneration."
    )


def _make_group_converter(tier: str, group: str) -> click.Command:
    """Create a tier/group conversion command from the canonical registry."""

    @click.command(
        name=f"{tier}-{group}",
        short_help=f"Convert {tier} tier {group} group only.",
        help=f"Convert the {group} group from the {tier}-tier LEB1 file.",
    )
    def command() -> None:
        _leb2(
            [
                "convert",
                f"data/leb/ephemeris_{tier}.leb",
                "-o",
                f"data/leb2/{tier}_{group}.leb2",
                "--group",
                group,
                "--tier",
                tier,
            ]
        )

    return command


_EXPLICIT_GROUP_COMMANDS = {
    "base-core",
    "base-asteroids",
    "base-apogee",
    "base-uranians",
}
for _tier in ("base", "medium", "extended"):
    for _group in LEB2_GROUPS:
        if f"{_tier}-{_group}" not in _EXPLICIT_GROUP_COMMANDS:
            convert_group.add_command(_make_group_converter(_tier, _group))


leb2_group.add_command(convert_group)


# ===========================================================================
# leph leb2 verify
# ===========================================================================


@click.group(
    "verify",
    short_help="Verify a tier's core LEB2 companion against LEB1.",
    help=(
        "Verify one tier's core LEB2 companion against its LEB1 reference.\n\n"
        "The exact 14-body inventory is required. Random dates are sampled for\n"
        "each body and errors are reported in native stored component units."
    ),
)
def verify_group() -> None:
    """LEB2 verification commands."""


@verify_group.command(
    "base",
    short_help="Verify the base core companion against LEB1.",
)
def verify_base() -> None:
    """Verify ``base_core.leb2`` against the base-tier LEB1 reference."""
    _leb2(
        [
            "verify",
            "data/leb2/base_core.leb2",
            "--reference",
            "data/leb/ephemeris_base.leb",
            "--samples",
            "500",
            "--group",
            "core",
            "--tier",
            "base",
        ]
    )


@verify_group.command(
    "medium",
    short_help="Verify the medium core companion against LEB1.",
)
def verify_medium() -> None:
    """Verify ``medium_core.leb2`` against the medium-tier LEB1 reference."""
    _leb2(
        [
            "verify",
            "data/leb2/medium_core.leb2",
            "--reference",
            "data/leb/ephemeris_medium.leb",
            "--samples",
            "500",
            "--group",
            "core",
            "--tier",
            "medium",
        ]
    )


@verify_group.command(
    "extended",
    short_help="Verify the extended core companion against LEB1.",
)
def verify_extended() -> None:
    """Verify ``extended_core.leb2`` against the extended-tier LEB1 reference."""
    _leb2(
        [
            "verify",
            "data/leb2/extended_core.leb2",
            "--reference",
            "data/leb/ephemeris_extended.leb",
            "--samples",
            "500",
            "--group",
            "core",
            "--tier",
            "extended",
        ]
    )


leb2_group.add_command(verify_group)
