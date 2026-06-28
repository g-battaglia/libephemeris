# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Release commands: upload LEB files to GitHub Releases.

Replaces 5 poe tasks: release:leb, release:leb:base/medium/extended, release:leb:dry-run.
Requires: gh CLI authenticated (gh auth login).
"""

from __future__ import annotations

import subprocess
import sys

import click


def _release(args: list[str]) -> None:
    """Run the release script."""
    sys.exit(subprocess.call([sys.executable, "scripts/release_leb.py", *args]))


@click.group(
    "release",
    short_help="Upload LEB files to GitHub Releases and update download hashes.",
    help="Upload LEB binary ephemeris files to GitHub Releases.\n\nCreates or updates a GitHub Release with LEB files for distribution.\nAlso updates the download hash table in the source code so that\nthe production CLI can verify file integrity after download.\n\nRequires: gh CLI authenticated (gh auth login).\n\n  leph release leb 1.0.0            # Upload all tiers\n  leph release leb-dry-run 1.0.0    # Preview without uploading",
)
def release_group() -> None:
    """Release commands."""


@release_group.command(
    "leb",
    short_help="Upload all LEB files to GitHub release and update hashes.",
)
@click.argument("version")
def release_leb(version: str) -> None:
    """Upload all LEB files to GitHub release and update download hashes.

    VERSION is the release version string, e.g. '1.0.0' or '1.1.0'.
    Requires: gh CLI authenticated (gh auth login).
    """
    _release(["--version", version, "--update-hashes"])


@release_group.command(
    "leb-base",
    short_help="Upload base tier LEB to GitHub release and update hashes.",
)
@click.argument("version")
def release_leb_base(version: str) -> None:
    """Upload base tier LEB to GitHub release and update hashes.

    VERSION is the release version string, e.g. '1.0.0'.
    """
    _release(["--version", version, "--tier", "base", "--update-hashes"])


@release_group.command(
    "leb-medium",
    short_help="Upload medium tier LEB to GitHub release and update hashes.",
)
@click.argument("version")
def release_leb_medium(version: str) -> None:
    """Upload medium tier LEB to GitHub release and update hashes.

    VERSION is the release version string, e.g. '1.0.0'.
    """
    _release(["--version", version, "--tier", "medium", "--update-hashes"])


@release_group.command(
    "leb-extended",
    short_help="Upload extended tier LEB to GitHub release and update hashes.",
)
@click.argument("version")
def release_leb_extended(version: str) -> None:
    """Upload extended tier LEB to GitHub release and update hashes.

    VERSION is the release version string, e.g. '1.0.0'.
    """
    _release(["--version", version, "--tier", "extended", "--update-hashes"])


@release_group.command(
    "leb-dry-run",
    short_help="Dry run: show what would be uploaded without uploading.",
)
@click.argument("version")
def release_leb_dry_run(version: str) -> None:
    """Dry run: show what LEB files would be uploaded without uploading.

    VERSION is the release version string, e.g. '1.0.0'.
    """
    _release(["--version", version, "--dry-run"])
