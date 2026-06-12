# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Diagnostic commands: tier positions, data download verification.

Replaces 4 poe tasks: diag:base, diag:medium, diag:extended, diag:download.
"""

from __future__ import annotations

import subprocess
import sys

import click


def _python(args: list[str]) -> None:
    """Run a python script."""
    sys.exit(subprocess.call([sys.executable, *args]))


@click.group(
    "diag",
    short_help="Print body positions per tier, verify data integrity, download missing files.",
    help="Diagnostic tools: print computed positions for every supported body,\nverify data file integrity, and download missing data.\n\nUseful after downloading data or changing tiers to confirm everything works.\n\n  leph diag positions-medium   # Print all body positions for medium tier\n  leph diag download-data      # Download all required files for current tier",
)
def diag_group() -> None:
    """Diagnostic tools."""


@diag_group.command(
    "positions-base",
    short_help="Print positions for every body using base tier (1850-2150).",
)
def positions_base() -> None:
    """Print ecliptic/equatorial positions for every body using base tier (de440s, 1850-2150).

    Shows longitude, latitude, distance, speed, and data source for each
    supported celestial body at the current date. Useful for verifying
    data integrity after download or tier change.
    """
    _python(["scripts/tier_diagnostic_base.py"])


@diag_group.command(
    "positions-medium",
    short_help="Print positions for every body using medium tier (1550-2650).",
)
def positions_medium() -> None:
    """Print ecliptic/equatorial positions for every body using medium tier (de440, 1550-2650)."""
    _python(["scripts/tier_diagnostic_medium.py"])


@diag_group.command(
    "positions-extended",
    short_help="Print positions for every body using extended tier (-13200 to +17191).",
)
def positions_extended() -> None:
    """Print ecliptic/equatorial positions for every body using extended tier (de441, -13200 to +17191)."""
    _python(["scripts/tier_diagnostic_extended.py"])


@diag_group.command(
    "download-data",
    short_help="Download all required data files for the current tier.",
)
def download_data() -> None:
    """Download all required data files for the current tier with progress bars.

    Downloads BSP kernels (DE440/441), asteroid SPKs, delta-T tables,
    leap-second files, and planet-center corrections. Safe to re-run:
    already-downloaded files are skipped.
    """
    _python(
        [
            "-c",
            "import libephemeris as eph; eph.ensure_all_ephemerides(show_progress=True)",
        ]
    )
