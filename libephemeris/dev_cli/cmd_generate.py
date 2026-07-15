# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Data generation commands for planet-center SPKs and Keplerian elements.

Replaces the planet-center and multi-epoch Keplerian poe tasks.

Provenance:
    Project-authored command dispatch. The invoked generator scripts document
    their JPL inputs, transformations, frames, and generated artifacts; this
    wrapper contains no independent scientific algorithm.
"""

from __future__ import annotations

import subprocess
import sys

import click


def _python(args: list[str]) -> None:
    """Run a python script."""
    sys.exit(subprocess.call([sys.executable, *args]))


@click.group(
    "generate",
    short_help="Generate planet-center SPKs and Keplerian elements.",
    help="Generate derived data files from raw SPK kernels.\n\n"
    "These are one-time generation steps that produce data files used by\n"
    "the library at runtime. Most developers never need to run these unless\n"
    "they are updating the underlying data or adding new bodies.\n\n"
    "  planet-centers   COB-corrected SPKs for sub-arcsecond gas giant positions\n"
    "  keplerian-elements Multi-epoch orbital elements for Keplerian fallback",
)
def generate_group() -> None:
    """Data generation commands."""


# ===========================================================================
# Planet centers SPK
# ===========================================================================


@generate_group.command(
    "planet-centers-base",
    short_help="Build COB-corrected SPK for base tier (1850-2150).",
)
def planet_centers_base() -> None:
    """Build center-of-body corrected SPK for base tier (de440s, 1850-2150).

    Creates SPK files that correct Jupiter, Saturn, Uranus, and Neptune
    positions from barycenter to center-of-body using satellite ephemeris data.
    This provides sub-arcsecond accuracy for gas giant positions.
    Requires: spiceypy (pip install spiceypy). Downloads satellite SPKs from JPL.
    """
    _python(["scripts/generate_planet_centers_spk.py", "--tier", "base"])


@generate_group.command(
    "planet-centers-medium",
    short_help="Build COB-corrected SPK for medium tier (1550-2650).",
)
def planet_centers_medium() -> None:
    """Build center-of-body corrected SPK for medium tier (de440, 1550-2650).

    Same as planet-centers-base but for the wider medium tier date range.
    Requires spiceypy.
    """
    _python(["scripts/generate_planet_centers_spk.py", "--tier", "medium"])


@generate_group.command(
    "planet-centers-extended",
    short_help="Build COB-corrected SPK for extended tier (full range).",
)
def planet_centers_extended() -> None:
    """Build center-of-body corrected SPK for extended tier (de441, full range).

    Requires spiceypy.
    """
    _python(["scripts/generate_planet_centers_spk.py", "--tier", "extended"])


@generate_group.command(
    "planet-centers-all",
    short_help="Build COB-corrected SPKs for all three tiers.",
)
def planet_centers_all() -> None:
    """Build center-of-body corrected SPKs for all three tiers at once.

    Runs base, medium, and extended sequentially. Requires spiceypy.
    """
    _python(["scripts/generate_planet_centers_spk.py", "--all"])


@generate_group.command(
    "planet-centers-spk",
    short_help="(Legacy alias) Same as planet-centers-medium.",
)
def planet_centers_spk() -> None:
    """(Legacy alias) Same as planet-centers-medium."""
    _python(["scripts/generate_planet_centers_spk.py", "--tier", "medium"])


# ===========================================================================
# Keplerian elements
# ===========================================================================


@generate_group.command(
    "keplerian-elements",
    short_help="Generate multi-epoch orbital elements from SPK files.",
)
def keplerian_elements() -> None:
    """Generate multi-epoch Keplerian orbital elements from SPK files.

    Attempts every configured minor body on the reviewed 1650–2450 ten-year
    grid and emits entries only where a local, heliocentric J2000 type-21 SPK
    has coverage. The emitted map fragment is reviewed before incorporation in
    ``multi_epoch_data.py``; the currently shipped table covers 36 of 37
    bodies. These osculating nodes limit two-body propagation distance when an
    SPK or optional numerical integrator is unavailable. Precision is
    body/date dependent and is not represented by a universal degree bound.
    """
    _python(["scripts/generate_multi_epoch_elements.py"])


@generate_group.command(
    "keplerian-dry-run",
    short_help="Show available SPK files and epoch grid (no generation).",
)
def keplerian_dry_run() -> None:
    """Show which SPK files are available and the epoch grid, without generating.

    Useful to check which bodies have SPK coverage before running the
    full generation. No files are created or modified.
    """
    _python(["scripts/generate_multi_epoch_elements.py", "--dry-run"])
