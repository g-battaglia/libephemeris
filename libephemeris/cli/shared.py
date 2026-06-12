# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Shared CLI utilities for both production (libephemeris) and dev (leph) CLIs.

Contains tier metadata, download helpers, and common CLI patterns reused
across both entry points.
"""

from __future__ import annotations

import click

# ---------------------------------------------------------------------------
# Tier metadata — single source of truth for both CLIs
# ---------------------------------------------------------------------------

TIER_INFO = {
    "base": {
        "label": "base",
        "ephemeris": "de440s.bsp (~31 MB)",
        "range": "1850-2150 CE",
        "spk_range": "1850-2150",
        "description": "Lightweight tier for modern-era calculations.",
    },
    "medium": {
        "label": "medium",
        "ephemeris": "de440.bsp (~114 MB)",
        "range": "1550-2650 CE",
        "spk_range": "1900-2100",
        "description": "General purpose tier (default). Covers most historical and future dates.",
    },
    "extended": {
        "label": "extended",
        "ephemeris": "de441.bsp (~3.1 GB)",
        "range": "-13200 to +17191 CE",
        "spk_range": "1600-2500 (JPL Horizons limit)",
        "description": "Full extended range for deep historical and far-future research.",
    },
}

TIERS = list(TIER_INFO.keys())


def tier_download_help(tier: str) -> str:
    """Build detailed help text for a download subcommand."""
    info = TIER_INFO[tier]
    return f"""\
Download all data files for the '{tier}' precision tier.

  Ephemeris:  {info["ephemeris"]}
  Date range: {info["range"]}
  SPK range:  {info["spk_range"]}

{info["description"]}

Downloads:
  1. The ephemeris file ({info["ephemeris"].split(" ")[0]})
  2. planet_centers.bsp precision offsets (~25 MB)
  3. SPK kernels for 21 minor bodies (asteroids, centaurs, TNOs)
"""


def leb_download_help(tier: str) -> str:
    """Build detailed help text for a LEB download subcommand."""
    info = TIER_INFO[tier]
    sizes = {"base": "~53 MB", "medium": "~175 MB", "extended": "(size varies)"}
    return f"""\
Download the precomputed LEB binary ephemeris for the '{tier}' tier.

LEB files contain Chebyshev polynomial approximations for all celestial bodies,
providing ~14x speedup over the Skyfield/JPL pipeline.

  Tier:       {tier} ({info["range"]})
  File:       ephemeris_{tier}.leb ({sizes.get(tier, "")})
  Bodies:     31 (Sun, Moon, planets, nodes, apsides, asteroids)

Files are saved to ~/.libephemeris/leb/ by default.
"""


# ---------------------------------------------------------------------------
# Common click options
# ---------------------------------------------------------------------------


def download_options(f: click.Command) -> click.Command:
    """Add common --force / --no-progress / --quiet flags to a download command."""
    f = click.option(
        "--force", "-f", is_flag=True, help="Force download even if files already exist"
    )(f)
    f = click.option("--no-progress", is_flag=True, help="Disable progress output")(f)
    f = click.option(
        "--quiet", "-q", is_flag=True, help="Suppress all output except errors"
    )(f)
    return f
