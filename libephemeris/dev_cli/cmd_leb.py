# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""LEB1 binary ephemeris generation and verification commands.

Replaces 28 poe tasks: leb:generate:*, leb:verify:*.

Structure:
  leph leb generate <tier> <mode>   -- generate LEB files
  leph leb verify <tier>            -- verify existing LEB files
"""

from __future__ import annotations

import subprocess
import sys

import click

from libephemeris.leb_groups import LEB1_GROUPS


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_TIERS = ["base", "medium", "extended"]

_TIER_INFO = {
    "base": "de440s, 1850-2150",
    "medium": "de440, 1550-2650",
    "extended": "de441, -5000 to 5000",
}

# Partial files merged into the final {tier}.leb — derived from the canonical
# group partition so a new group (e.g. exotics) is never silently omitted.
_MERGE_FILES = {
    tier: [f"data/leb/ephemeris_{tier}_{g}.leb" for g in LEB1_GROUPS]
    for tier in _TIERS
}


def _gen(args: list[str]) -> None:
    """Run the LEB generator script."""
    sys.exit(subprocess.call([sys.executable, "scripts/generate_leb.py", *args]))


# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------


@click.group(
    "leb",
    short_help="Generate and verify LEB1 binary ephemeris files (Chebyshev polynomials).",
    help="Generate and verify LEB1 binary ephemeris files.\n\n"
    "LEB files store precomputed Chebyshev polynomial approximations of\n"
    "planetary positions. At runtime they provide ~14x speedup over computing\n"
    "positions from DE440 kernels via Skyfield.\n\n"
    "Three precision tiers are available:\n\n"
    "  base      de440s, 1850-2150    (~53 MB)\n"
    "  medium    de440,  1550-2650    (~175 MB, default)\n"
    "  extended  de441, -5000 to 5000 (~1.6 GB)\n\n"
    "Recommended workflow (avoids macOS multiprocessing deadlocks):\n\n"
    "  leph leb generate medium groups   # all body groups + merge\n"
    "  leph leb verify medium            # verify the generated file",
)
def leb_group() -> None:
    """LEB1 binary ephemeris management."""


# ===========================================================================
# leph leb generate — generation commands
# ===========================================================================


@click.group(
    "generate",
    short_help="Generate LEB1 files from Skyfield/DE440 (groups, full, single, body).",
    help="Generate LEB1 binary ephemeris files from Skyfield/DE440 reference data.\n\n"
    "Each tier has multiple generation modes:\n\n"
    "  groups      RECOMMENDED: generate planets/asteroids/analytical separately then merge\n"
    "  full        Generate all bodies at once (may deadlock on macOS)\n"
    "  single      One body at a time (lowest memory, slowest)\n"
    "  body <name> Generate specific body(ies) by name or ID\n\n"
    "Example: leph leb generate medium groups",
)
def generate_group() -> None:
    """LEB generation commands."""


@generate_group.command(
    "all",
    short_help="Generate LEB for all three tiers sequentially.",
)
def generate_all() -> None:
    """Generate LEB files for all three tiers (base + medium + extended), sequentially."""
    for tier in _TIERS:
        click.echo(f"\n{'=' * 60}")
        click.echo(f"  Generating {tier} tier ({_TIER_INFO[tier]})")
        click.echo(f"{'=' * 60}\n")
        ret = subprocess.call(
            [sys.executable, "scripts/generate_leb.py", "--tier", tier, "--verify"]
        )
        if ret != 0:
            sys.exit(ret)


# --- Per-tier generation subgroups ---


def _make_tier_group(tier: str) -> click.Group:
    """Create a generation subgroup for a specific tier."""

    @click.group(
        tier,
        help=f"Generate {tier} tier LEB ({_TIER_INFO[tier]}).",
        short_help=f"Generate {tier} tier LEB ({_TIER_INFO[tier]}).",
    )
    def tier_group() -> None:
        pass

    @tier_group.command(
        "full",
        short_help=f"Generate {tier} tier, all bodies at once + verify.",
    )
    def full() -> None:
        """Generate all bodies at once + verify."""
        _gen(["--tier", tier, "--verify"])

    full.__doc__ = f"Generate {tier} tier LEB file ({_TIER_INFO[tier]}), all bodies at once + verify."

    @tier_group.command(
        short_help=f"Generate {tier} tier planets group (Sun-Pluto, Earth).",
    )
    def planets() -> None:
        """Generate planets group (Sun-Pluto, Earth)."""
        _gen(["--tier", tier, "--group", "planets"])

    planets.__doc__ = f"Generate {tier} tier planets group (Sun-Pluto, Earth)."

    @tier_group.command(
        short_help=f"Generate {tier} tier asteroids group (5 bodies).",
    )
    def asteroids() -> None:
        """Generate asteroids group (Chiron, Ceres, Pallas, Juno, Vesta)."""
        _gen(["--tier", tier, "--group", "asteroids"])

    asteroids.__doc__ = f"Generate {tier} tier asteroids group (Chiron, Ceres-Vesta)."

    @tier_group.command(
        short_help=f"Generate {tier} tier analytical group (nodes, Lilith, Uranians).",
    )
    def analytical() -> None:
        """Generate analytical group (nodes, Lilith, Uranians)."""
        _gen(["--tier", tier, "--group", "analytical"])

    analytical.__doc__ = (
        f"Generate {tier} tier analytical group (nodes, Lilith, Uranians)."
    )

    @tier_group.command(
        short_help=f"Merge {tier} tier partial files into one .leb + verify.",
    )
    def merge() -> None:
        """Merge partial group files into a single .leb + verify."""
        _gen(["--tier", tier, "--merge", *_MERGE_FILES[tier], "--verify"])

    merge.__doc__ = f"Merge {tier} tier partial files into ephemeris_{tier}.leb."

    _groups_desc = " + ".join(LEB1_GROUPS) + " + merge"

    @tier_group.command(
        short_help=f"Full {tier} tier group workflow: {_groups_desc}.",
    )
    def groups() -> None:
        """Full group workflow: all body groups + merge.

        Recommended over 'full' to avoid macOS multiprocessing deadlocks.
        """
        steps = [
            (g, ["--tier", tier, "--group", g]) for g in LEB1_GROUPS
        ]
        steps.append(
            ("merge", ["--tier", tier, "--merge", *_MERGE_FILES[tier], "--verify"])
        )
        for step, step_args in steps:
            click.echo(f"\n--- {tier}/{step} ---\n")
            ret = subprocess.call(
                [sys.executable, "scripts/generate_leb.py", *step_args]
            )
            if ret != 0:
                sys.exit(ret)

    groups.__doc__ = f"Generate {tier} tier via groups ({' + '.join(LEB1_GROUPS)}) then merge."

    @tier_group.command(
        short_help=f"Generate {tier} tier one body at a time (lowest memory).",
    )
    def single() -> None:
        """Generate one body at a time (lowest memory usage) + verify."""
        _gen(["--tier", tier, "--single", "--verify"])

    single.__doc__ = f"Generate {tier} tier one body at a time (lowest memory usage)."

    @tier_group.command(
        short_help=f"Generate {tier} tier for specific body(ies) by name/ID.",
    )
    @click.argument("body")
    def body(body: str) -> None:
        """Generate specific body(ies).

        BODY can be a name (moon, sun) or numeric ID (1, 2), or comma-separated list.
        Examples: 'moon', '1,2,3', 'sun,mercury'.
        """
        _gen(["--tier", tier, "--bodies", body])

    body.__doc__ = (
        f"Generate {tier} tier for specific body(ies) (e.g. 'moon', '1,2,3')."
    )

    return tier_group


for _tier in _TIERS:
    generate_group.add_command(_make_tier_group(_tier))

leb_group.add_command(generate_group)


# ===========================================================================
# leph leb verify — verify existing .leb files
# ===========================================================================


@click.group(
    "verify",
    short_help="Verify existing .leb files against Skyfield reference.",
    help="Verify existing LEB files against Skyfield reference without regenerating.\n\nReads each body from the .leb file, computes the same position via Skyfield,\nand reports the maximum error in arcseconds.",
)
def verify_group() -> None:
    """LEB verification commands."""


def _make_verify_cmd(tier: str) -> click.Command:
    """Create a verify command for a tier."""

    @click.command(
        tier,
        short_help=f"Verify {tier} tier .leb file ({_TIER_INFO[tier]}).",
    )
    def verify_cmd() -> None:
        _gen(["--tier", tier, "--verify-only"])

    verify_cmd.__doc__ = f"Verify existing {tier} tier .leb file ({_TIER_INFO[tier]})."
    return verify_cmd


for _tier in _TIERS:
    verify_group.add_command(_make_verify_cmd(_tier))

leb_group.add_command(verify_group)
