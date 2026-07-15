# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Data download commands for local libephemeris development.

The developer CLI downloads only local-development prerequisites and source
inputs. Generated artifacts such as LEB files and planet-centers outputs are
produced by dedicated generation commands, not downloaded here.
"""

from __future__ import annotations

import subprocess
import sys

import click


def _run_python(args: list[str]) -> int:
    """Run a Python module or script and return its exit code."""
    return subprocess.call([sys.executable, *args])


def _run_script_or_raise(args: list[str], description: str) -> None:
    """Run a helper script and raise a Click-friendly error on failure."""
    exit_code = _run_python(args)
    if exit_code != 0:
        raise click.ClickException(f"{description} failed with exit code {exit_code}.")


def _download_spk_tier(tier: str, force: bool = False) -> None:
    """Download the DE kernel and asteroid SPKs for a specific tier."""
    import libephemeris as eph

    eph.set_precision_tier(tier)
    eph.ensure_all_ephemerides(force_download=force, show_progress=True)


def _download_spk_extended(force: bool = False) -> None:
    """Download `de441.bsp` plus max-range minor-body SPKs for extended work."""
    from libephemeris.state import get_loader

    click.echo("Ensuring de441.bsp is available...")
    get_loader()("de441.bsp")

    args = ["scripts/download_max_range_spk.py"]
    if force:
        args.append("--force")
    _run_script_or_raise(args, "Extended-tier SPK bootstrap")


def _download_iers(force: bool = False) -> None:
    """Download all IERS files used for Earth-orientation and Delta T work."""
    from libephemeris.iers_data import (
        download_delta_t_data,
        download_iers_finals,
        download_leap_seconds,
    )

    finals = download_iers_finals(force=force)
    leap_seconds = download_leap_seconds(force=force)
    delta_t = download_delta_t_data(force=force)

    click.echo(f"  finals2000A.data  -> {finals}")
    click.echo(f"  leap_seconds.dat -> {leap_seconds}")
    click.echo(f"  deltat.data      -> {delta_t}")


def _download_planet_center_sources(
    force: bool = False,
    tier: str | None = None,
) -> None:
    """Download cached NAIF source kernels for planet-center generation."""
    args = ["scripts/generate_planet_centers_spk.py", "--download-only"]
    if tier is None:
        args.append("--all")
    else:
        args.extend(["--tier", tier])
    if force:
        args.append("--force")
    _run_script_or_raise(args, "Planet-center source bootstrap")


def _download_assist(force: bool = False) -> None:
    """Download ASSIST n-body data used by local development workflows."""
    from libephemeris.rebound_integration import download_assist_data

    download_assist_data(force=force, show_progress=True, quiet=False)


# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------


@click.group(
    "download",
    short_help="Download source and runtime prerequisites for development.",
    help="Download data files needed for libephemeris development.\n\n"
    "\b\n"
    "Categories:\n"
    "  all                    Full downloadable bootstrap for local development\n"
    "  spk-*        DE kernels + asteroid SPKs for a specific tier\n"
    "  iers                   Earth-orientation and Delta T data\n"
    "  planet-centers-sources Cached NAIF source kernels for later generation\n"
    "  assist                 REBOUND/ASSIST n-body files\n\n"
    "Use `leph download all` to prepare a machine for local development, data\n"
    "generation, verification, and packaging workflows.\n\n"
    "Generated artifacts are intentionally excluded here: use `leph leb ...`,\n"
    "`leph leb2 ...`, and `leph generate planet-centers-*` to build them.\n\n"
    "\b\n"
    "Examples:\n"
    "  leph download all\n"
    "  leph download spk-medium\n"
    "  leph download planet-centers-sources\n"
    "  leph download iers",
)
def download_group() -> None:
    """Data download commands."""


@download_group.command(
    "all",
    short_help="Download the full bootstrap dataset for local development.",
)
@click.option(
    "--force",
    is_flag=True,
    help="Re-download files even if they are already present.",
)
def download_all(force: bool) -> None:
    """Download every downloadable prerequisite needed for local development.

    This bootstrap intentionally excludes generated artifacts.

    \b
    Includes:
      1. DE kernels and asteroid SPKs for base, medium, and extended workflows.
      2. IERS data files used by time-scale and Delta T logic.
      3. Cached NAIF source kernels used later by `leph generate planet-centers-*`.
      4. ASSIST data used by n-body tooling.
    """
    phases = [
        ("SPK base", lambda: _download_spk_tier("base", force=force)),
        ("SPK medium", lambda: _download_spk_tier("medium", force=force)),
        ("SPK extended", lambda: _download_spk_extended(force=force)),
        ("IERS", lambda: _download_iers(force=force)),
        (
            "planet-center sources",
            lambda: _download_planet_center_sources(force=force),
        ),
        ("ASSIST", lambda: _download_assist(force=force)),
    ]

    for label, action in phases:
        click.echo(f"\n== {label} ==")
        action()


# ===========================================================================
# SPK downloads
# ===========================================================================


@download_group.command(
    "spk-base",
    short_help="Download asteroid SPK kernels for base tier (1850-2150).",
)
@click.option(
    "--force",
    is_flag=True,
    help="Re-download files even if they are already present.",
)
def spk_base(force: bool) -> None:
    """Download asteroid SPK21 kernel files for base tier (1850-2150).

    Downloads SPK kernels for Chiron, Ceres, Pallas, Juno, Vesta and other
    minor bodies. Also downloads DE440s if not already present.
    Files are saved to ~/.libephemeris/ by default.
    """
    _download_spk_tier("base", force=force)


@download_group.command(
    "spk-medium",
    short_help="Download asteroid SPK kernels for medium tier (1550-2650).",
)
@click.option(
    "--force",
    is_flag=True,
    help="Re-download files even if they are already present.",
)
def spk_medium(force: bool) -> None:
    """Download asteroid SPK21 kernel files for medium tier (1550-2650).

    Downloads SPK kernels for all supported minor bodies plus DE440.
    """
    _download_spk_tier("medium", force=force)


@download_group.command(
    "spk-extended",
    short_help="Download de441 plus max-range asteroid SPKs for extended work.",
)
@click.option(
    "--force",
    is_flag=True,
    help="Re-download files even if they are already present.",
)
def spk_extended(force: bool) -> None:
    """Download de441 plus max-range asteroid SPK files for extended work.

    The DE441 kernel is downloaded to the shared data directory. Minor-body SPKs
    are then fetched from JPL Horizons with the widest practical range currently
    available to the project.
    """
    _download_spk_extended(force=force)


@download_group.command(
    "iers",
    short_help="Download IERS Earth-orientation and Delta T data.",
)
@click.option(
    "--force",
    is_flag=True,
    help="Re-download files even if they are already present.",
)
def iers(force: bool) -> None:
    """Download all IERS files used by local time and Earth-orientation logic."""
    _download_iers(force=force)


@download_group.command(
    "planet-centers-sources",
    short_help="Download NAIF source kernels for later planet-centers generation.",
)
@click.option(
    "--tier",
    type=click.Choice(["base", "medium", "extended"], case_sensitive=False),
    default=None,
    help="Limit the source download to one tier (default: all tiers).",
)
@click.option(
    "--force",
    is_flag=True,
    help="Re-download files even if they are already present.",
)
def planet_centers_sources(tier: str | None, force: bool) -> None:
    """Download cached source kernels for `leph generate planet-centers-*`."""
    _download_planet_center_sources(force=force, tier=tier)


# ===========================================================================
# ASSIST download
# ===========================================================================


@download_group.command(
    short_help="Download REBOUND/ASSIST n-body data (~120 MB).",
)
@click.option(
    "--force",
    is_flag=True,
    help="Re-download files even if they are already present.",
)
def assist(force: bool) -> None:
    """Download REBOUND/ASSIST n-body integration data (~120 MB).

    Downloads planet ephemeris and asteroid perturber SPK files needed
    for high-precision TNO/asteroid orbit propagation via REBOUND/ASSIST.
    Files are saved to ~/.libephemeris/assist/.
    Requires: pip install libephemeris[nbody]
    """
    _download_assist(force=force)
