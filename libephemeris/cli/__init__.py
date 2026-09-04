# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Command-line interface for libephemeris.

This module provides CLI commands for managing libephemeris data files.
Migrated from argparse to click for tab completion and consistency with
the dev CLI (``leph``).

Usage:
    libephemeris download base          Download data for 'base' tier (1850-2150)
    libephemeris download medium        Download data for 'medium' tier (1550-2650)
    libephemeris download extended      Download data for 'extended' tier (-13200 to +17191)
    libephemeris download leb2-base     Install reviewed LEB2 assets for 'base'
    libephemeris download leb2-medium   Install reviewed LEB2 assets for 'medium'
    libephemeris download leb2-extended Install reviewed LEB2 assets for 'extended'
    libephemeris download iers          Refresh IERS data explicitly
    libephemeris download assist        Download ASSIST n-body data files (~714 MB)
    libephemeris status                 Show comprehensive library and data status
    libephemeris --version              Show version
    libephemeris --help                 Show help

Shell completion (add to your shell profile):

    # zsh
    eval "$(_LIBEPHEMERIS_COMPLETE=zsh_source libephemeris)"

    # bash
    eval "$(_LIBEPHEMERIS_COMPLETE=bash_source libephemeris)"

Provenance:
    Project-authored Click command orchestration. It contains no astronomical
    model or coefficient and delegates downloads, status inspection, and
    initialization to their registered modules. Tier labels and file identities
    are distribution metadata, not scientific values.
"""

from __future__ import annotations

import sys

import click

from .. import __version__
from .shared import SPK_BODY_COUNT, TIER_INFO, leb_download_help, tier_download_help


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _handle_download(func, quiet: bool, **kwargs) -> None:  # type: ignore[no-untyped-def]
    """Run a download function with standard error handling."""
    try:
        func(quiet=quiet, **kwargs)
    except KeyboardInterrupt:
        click.echo("\nDownload cancelled.")
        sys.exit(130)
    except ImportError as e:
        # --quiet suppresses progress, never errors
        click.echo(
            f"Error: {e}\nInstall the nbody extra: pip install libephemeris[nbody]",
            err=True,
        )
        sys.exit(1)
    except (OSError, ValueError, RuntimeError) as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Root CLI group
# ---------------------------------------------------------------------------


_HELP = f"""\
\b
  _ _ _              _                      _    
 | (_) |__  ___ _ __| |_  ___ _ __  ___ _ _(_)___
 | | | '_ \\/ -_) '_ \\ ' \\/ -_) '  \\/ -_) '_| (_-<
 |_|_|_.__/\\___| .__/_||_\\___|_|_|_\\___|_| |_/__/
               |_|   [ Powered by NASA JPL DE440/441 ]
                     [ v{__version__} ]

This CLI manages data files, shows library status, and configures
the calculation backend. It is intended for end-users and CI pipelines.

\b
First-time setup:

  libephemeris download medium       Download data for the default tier
  libephemeris status                Verify everything is installed"""

_EPILOG = """\
\b
Examples:
  libephemeris download medium        Download data for the default tier
  libephemeris download base          Lightweight, modern-era data
  libephemeris download extended      Full range (-13200 to +17191 CE)
  libephemeris download leb2-base     Install reviewed modular LEB2 assets
  libephemeris download iers          Refresh mutable IERS observations
  libephemeris download assist        ASSIST n-body data (~714 MB)
  libephemeris status                 Show comprehensive library and data status

For more information, visit: https://github.com/g-battaglia/libephemeris"""


@click.group(
    name="libephemeris",
    help=_HELP,
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 120},
    epilog=_EPILOG,
)
@click.version_option(__version__, prog_name="libephemeris")
def cli() -> None:
    """Root CLI group."""


# ---------------------------------------------------------------------------
# download subgroup
# ---------------------------------------------------------------------------


@click.group(
    "download",
    short_help="Download tier SPKs, reviewed LEB2 assets, and ASSIST data.",
    help="Download data files required by libephemeris.\n\n"
    "Available data:\n\n"
    "  Tier data (base/medium/extended)   DE440/441 kernels + asteroid SPKs\n"
    "  LEB2 files (leb2-base/medium/ext)  Independently reviewed modular assets\n"
    "  IERS data                          Earth-orientation observations\n"
    "  ASSIST data                        N-body integration files (~714 MB)\n\n"
    "Historical leb-* commands install the pinned core in the current format.\n\n"
    "Most users need only: libephemeris download medium",
)
def download_group() -> None:
    """Download subcommands.

    Entering this explicit provisioning boundary authorizes network access even
    when the calculation configuration itself is sealed.
    """
    from ..net import set_network_policy

    set_network_policy("allow")


# --- Common download options ---


def _download_options(f):  # type: ignore[no-untyped-def]
    """Add --force / --no-progress / --quiet to a download command."""
    f = click.option(
        "--force", "-f", is_flag=True, help="Force download even if files already exist"
    )(f)
    f = click.option("--no-progress", is_flag=True, help="Disable progress output")(f)
    f = click.option(
        "--quiet", "-q", is_flag=True, help="Suppress all output except errors"
    )(f)
    return f


# --- Tier data downloads ---


def _make_tier_download(tier: str) -> click.Command:
    """Create a download command for a tier."""

    @click.command(
        tier,
        short_help=f"Download DE kernel + SPKs for '{tier}' tier ({TIER_INFO[tier]['range']}).",
        help=f"Download data for '{tier}' tier ({TIER_INFO[tier]['range']}).\n\n"
        + tier_download_help(tier),
    )
    @_download_options
    def cmd(force: bool, no_progress: bool, quiet: bool) -> None:
        from ..download import download_for_tier

        _handle_download(
            download_for_tier,
            quiet=quiet,
            tier_name=tier,
            force=force,
            show_progress=not no_progress,
        )

    return cmd


for _tier in TIER_INFO:
    download_group.add_command(_make_tier_download(_tier))


# --- LEB downloads ---


def _make_leb_download(tier: str) -> click.Command:
    """Create the compatibility command that installs a tier's LEB core."""

    @click.command(
        f"leb-{tier}",
        short_help=f"Install the pinned core LEB asset for '{tier}'.",
        help=leb_download_help(tier),
    )
    @click.option(
        "--force", "-f", is_flag=True, help="Reinstall even if already present."
    )
    @click.option(
        "--no-progress", is_flag=True, help="Disable download progress display."
    )
    @click.option("--quiet", "-q", is_flag=True, help="Suppress progress output.")
    def cmd(force: bool, no_progress: bool, quiet: bool) -> None:
        from ..download import download_leb_for_tier

        _handle_download(
            download_leb_for_tier,
            quiet=quiet,
            tier_name=tier,
            force=force,
            show_progress=not no_progress,
            activate=False,
        )

    return cmd


for _tier in TIER_INFO:
    download_group.add_command(_make_leb_download(_tier))


# --- LEB2 downloads ---


def _make_leb2_download(tier: str) -> click.Command:
    """Create a LEB2 download command for a tier."""

    @click.command(
        f"leb2-{tier}",
        short_help=f"Install reviewed LEB2 assets for the '{tier}' tier.",
        help=f"Install reviewed LEB2 assets for the '{tier}' tier.\n\n"
        f"LEB2 uses error-bounded lossy compression (mantissa truncation + zstd)\n"
        f"to achieve 4-10x smaller files.\n\n"
        f"Only SHA-256-pinned entries in the current manifest are used.\n"
        f"The reviewed base core is copied from the wheel; pinned medium and\n"
        f"extended cores are fetched from the data release. Other precomputed\n"
        f"companions are skipped pending artifact-specific review.\n\n"
        f"Files are saved to ~/.libephemeris/leb/ by default.",
    )
    @_download_options
    def cmd(force: bool, no_progress: bool, quiet: bool) -> None:
        from ..download import download_leb2_for_tier

        _handle_download(
            download_leb2_for_tier,
            quiet=quiet,
            tier_name=tier,
            force=force,
            show_progress=not no_progress,
            activate=True,
        )

    return cmd


for _tier in TIER_INFO:
    download_group.add_command(_make_leb2_download(_tier))


# --- Mutable IERS observations (explicit provisioning only) ---


@download_group.command(
    "iers",
    short_help="Refresh cached IERS Earth-orientation observations.",
)
@click.option(
    "--force", "-f", is_flag=True, help="Refresh even when the cache is fresh."
)
@click.option("--quiet", "-q", is_flag=True, help="Suppress progress output.")
def download_iers(force: bool, quiet: bool) -> None:
    """Provision IERS files without enabling calculation-time network access."""
    from ..iers_data import (
        download_delta_t_data,
        download_iers_finals,
        download_leap_seconds,
    )

    downloads = (
        ("finals2000A.data", download_iers_finals),
        ("leap_seconds.dat", download_leap_seconds),
        ("deltat.data", download_delta_t_data),
    )
    try:
        for label, downloader in downloads:
            path = downloader(force=force)
            if not quiet:
                click.echo(f"  [OK] {label}: {path}")
    except KeyboardInterrupt:
        click.echo("\nDownload cancelled.")
        raise click.exceptions.Exit(130)
    except (OSError, ValueError, RuntimeError) as exc:
        raise click.ClickException(str(exc)) from exc


# --- ASSIST download ---


@download_group.command(
    short_help="Download ASSIST n-body data files (~714 MB, requires libephemeris[nbody]).",
)
@_download_options
@click.option("--no-planets", is_flag=True, help="Skip planet ephemeris download")
@click.option("--no-asteroids", is_flag=True, help="Skip asteroid perturbers download")
@click.option(
    "--target-dir",
    type=str,
    default=None,
    help="Directory to save files (default: ~/.libephemeris/assist/)",
)
def assist(
    force: bool,
    no_progress: bool,
    quiet: bool,
    no_planets: bool,
    no_asteroids: bool,
    target_dir: str | None,
) -> None:
    """Download ASSIST n-body data files (~714 MB).

    ASSIST provides sub-arcsecond precision for asteroid orbit propagation by
    including gravitational perturbations from the Sun, Moon, 8 planets, and
    16 massive asteroids.

    \b
    Downloads:
      1. Planet ephemeris (linux_p1550p2650.440, ~98 MB)
      2. Asteroid perturbers (sb441-n16.bsp, ~616 MB)

    Files are saved to ~/.libephemeris/assist/ by default.
    Requires: pip install libephemeris[nbody]
    """
    from ..rebound_integration import download_assist_data

    _handle_download(
        download_assist_data,
        quiet=quiet,
        target_dir=target_dir,
        planets=not no_planets,
        asteroids=not no_asteroids,
        force=force,
        show_progress=not no_progress,
    )


# --- download auto: config-aware smart download ---


@download_group.command(
    "auto",
    short_help="Download only the files required by your current configuration.",
)
@_download_options
def download_auto(force: bool, no_progress: bool, quiet: bool) -> None:
    """Smart download based on your libephemeris-config.toml.

    Reads your configuration (precision tier + calculation mode) and
    downloads exactly the files needed to make it work.  Already-present
    files are skipped unless --force is used.

    \b
    What gets downloaded depends on your config:
      auto mode     best-by-date LEB2 tiers + DE kernel + planet_centers + SPK
      leb mode      best-by-date LEB2 tiers only
      skyfield      DE kernel + planet_centers + SPK
      horizons      planet_centers (optional)

    \b
    Examples:
      libephemeris download auto           Download what your config needs
      libephemeris download auto --force   Re-download everything
    """
    from .init_wizard import _d, _g, _w, _y

    # Read config
    try:
        from ..state import get_calc_mode, get_precision_tier
    except Exception:
        click.echo("Error: could not read library configuration.", err=True)
        sys.exit(1)

    tier = get_precision_tier()
    mode = get_calc_mode()

    if not quiet:
        click.echo()
        click.echo(_w(f"  Downloading files for: {tier} tier, {mode} mode"))
        click.echo()

    try:
        # LEB2 for auto / leb
        if mode in ("auto", "leb"):
            if not quiet:
                click.echo(_d("  LEB2 ephemeris files..."))
            from ..download import download_leb2_for_tier

            tier_order = ("base", "medium", "extended")
            eligible_tiers = tier_order[: tier_order.index(tier) + 1]
            for eligible_tier in eligible_tiers:
                download_leb2_for_tier(
                    tier_name=eligible_tier,
                    force=force,
                    show_progress=not no_progress,
                    quiet=quiet,
                    activate=False,
                )

            # Do not pin the last downloaded core as an explicit single file.
            # Clear any prior activation so state auto-discovery can compose
            # the eligible tiers and select the best source for each date.
            from ..state import set_leb_file

            set_leb_file(None)
            if not quiet:
                click.echo()

        # DE kernel + planet_centers + SPK for auto / skyfield
        if mode in ("auto", "skyfield"):
            if not quiet:
                click.echo(_d("  DE kernel + planet centers + SPK kernels..."))
            from ..download import download_for_tier

            download_for_tier(
                tier_name=tier,
                force=force,
                show_progress=not no_progress,
                quiet=quiet,
            )
        elif mode == "horizons":
            if not quiet:
                click.echo(_d("  Planet centers..."))
            from ..download import _download_planet_centers_for_tier

            _download_planet_centers_for_tier(
                tier_name=tier,
                force=force,
                show_progress=not no_progress,
                quiet=quiet,
            )

        if not quiet:
            click.echo()
            click.echo(f"  {_g('All done!')} Run 'libephemeris status' to verify.")
            click.echo()

    except KeyboardInterrupt:
        click.echo("\n  Download cancelled.")
        sys.exit(130)
    except (OSError, ValueError, RuntimeError) as e:
        click.echo(f"\n  {_y(f'Error: {e}')}", err=True)
        sys.exit(1)


# --- download all: download everything for all tiers/modes ---

_DOWNLOAD_ALL_HELP = f"""Download every data file for complete offline readiness.

Downloads LEB2 files, DE kernels, planet centers, SPK kernels, and IERS
Earth orientation data for ALL three tiers (base, medium, extended).

\b
WARNING: This will download approximately 5-6 GB of data.

\b
Files downloaded:
  - LEB2 ephemeris for base, medium, extended  (~1050 MB)
  - DE kernels: de440s, de440, de441            (~3.2 GB)
  - Planet centers for all tiers                (~320 MB)
  - SPK kernels for {SPK_BODY_COUNT} minor bodies             (~varies)
  - IERS data: finals, leap seconds, delta T    (~3 MB)

\b
Examples:
  libephemeris download all            Download everything
  libephemeris download all --force    Re-download everything
"""


@download_group.command(
    "all",
    short_help="Download ALL data files for every tier and mode (~5 GB + IERS).",
    help=_DOWNLOAD_ALL_HELP,
)
@_download_options
@click.option(
    "--yes",
    "-y",
    is_flag=True,
    help="Skip the ~5 GB download confirmation prompt",
)
def download_all(force: bool, no_progress: bool, quiet: bool, yes: bool) -> None:
    """Download every data file for complete offline readiness."""
    from .init_wizard import _d, _g, _w, _y

    tiers = ["base", "medium", "extended"]

    if not yes:
        if not click.confirm(
            "  About to download ~5-6 GB of data files. Continue?",
            default=False,
        ):
            click.echo("  Cancelled.")
            sys.exit(130)

    if not quiet:
        click.echo()
        click.echo(_w("  Downloading ALL data files for every tier"))
        click.echo(_d("  This will download ~5 GB of data."))
        click.echo()

    try:
        for tier in tiers:
            if not quiet:
                click.echo(_w(f"  ── {tier} tier ──"))
                click.echo()

            # LEB2
            if not quiet:
                click.echo(_d(f"  LEB2 {tier}..."))
            from ..download import download_leb2_for_tier

            download_leb2_for_tier(
                tier_name=tier,
                force=force,
                show_progress=not no_progress,
                quiet=quiet,
                activate=False,
            )

            # DE kernel + planet_centers + SPK
            if not quiet:
                click.echo()
                click.echo(_d(f"  DE kernel + planet centers + SPK for {tier}..."))
            from ..download import download_for_tier

            download_for_tier(
                tier_name=tier,
                force=force,
                show_progress=not no_progress,
                quiet=quiet,
            )

            if not quiet:
                click.echo()

        if not quiet:
            click.echo(_w("  ── IERS data ──"))
            click.echo()
            click.echo(_d("  Earth orientation parameters + leap seconds..."))

        from ..iers_data import (
            download_delta_t_data,
            download_iers_finals,
            download_leap_seconds,
        )

        download_iers_finals(force=force)
        download_leap_seconds(force=force)
        download_delta_t_data(force=force)

        if not quiet:
            click.echo(f"  {_g('All tiers downloaded!')} Full offline readiness.")
            click.echo(_d("  Run 'libephemeris status' to verify."))
            click.echo()

    except KeyboardInterrupt:
        click.echo("\n  Download cancelled.")
        sys.exit(130)
    except (ConnectionError, OSError, ValueError, RuntimeError) as e:
        click.echo(f"\n  {_y(f'Error: {e}')}", err=True)
        sys.exit(1)


cli.add_command(download_group)


# ---------------------------------------------------------------------------
# leb info — inspect the body/range contract stored in an artifact
# ---------------------------------------------------------------------------


@cli.group(
    "leb",
    short_help="Inspect LEB/LEB2 data files.",
)
def leb_group() -> None:
    """Inspect binary ephemeris artifacts without activating them."""


@leb_group.command(
    "info",
    short_help="Show bodies and per-body Julian Day coverage in a LEB file.",
)
@click.argument(
    "path",
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=str),
)
@click.option("--json", "as_json", is_flag=True, help="Output machine-readable JSON.")
def leb_info(path: str, as_json: bool) -> None:
    """Show the file identity and exact coverage encoded in PATH."""
    import json

    from ..inventory import inspect_leb_file

    try:
        info = inspect_leb_file(path)
    except (OSError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc

    if as_json:
        click.echo(json.dumps(info, indent=2))
        return

    click.echo(f"{info['name']} ({info['body_count']} bodies)")
    click.echo(f"Path: {info['path']}")
    if info["size_bytes"] is not None:
        click.echo(f"Size: {info['size_bytes']} bytes")
    click.echo("Body ID       JD start         JD end")
    for body in info["bodies"]:
        click.echo(
            f"{body['body_id']:>7}  {body['jd_start']:>15.6f}  {body['jd_end']:>15.6f}"
        )


# ---------------------------------------------------------------------------
# status command — comprehensive library and data overview
# ---------------------------------------------------------------------------


@cli.command(
    short_help="Show comprehensive library status: version, config, all data files.",
)
@click.option(
    "--json",
    "as_json",
    is_flag=True,
    help="Output status as machine-readable JSON.",
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase detail: -v shows paths, -vv shows full file lists.",
)
def status(as_json: bool, verbose: int) -> None:
    """Show comprehensive library and data file status.

    Displays version, calculation mode, precision tier, LEB file, data directory,
    and the status of all data files: DE kernels, planet center corrections,
    LEB1 binary ephemeris, LEB2 compressed files, SPK asteroid cache,
    ASSIST n-body data, and IERS Earth orientation data.

    \b
    Use --json for machine-readable output (e.g. for CI pipelines).
    Use -v or -vv for more detailed human-readable reports.
    """
    from ..download import print_data_status

    print_data_status(as_json=as_json, verbose=verbose)


# ---------------------------------------------------------------------------
# config command — explain data file locations and environment variables
# ---------------------------------------------------------------------------


@cli.command(
    short_help="Show data file locations, environment variables, and configuration guide.",
)
def config() -> None:
    """Show where data files are stored and how to configure libephemeris.

    Prints every configurable path, environment variable, and Python API
    function, organized by subsystem: data directory, precision tier,
    calculation mode, LEB ephemeris, SPK asteroid cache, IERS data, and
    ASSIST n-body integration.

    \b
    Useful when you need to:
      - Find where files are stored on disk
      - Override default paths via environment variables
      - Configure the library in a .env file or CI pipeline
      - Understand which env vars control which behavior
    """
    import os

    from ..download import get_data_dir

    data_dir = get_data_dir()

    _b = lambda t: click.style(t, bold=True)  # noqa: E731
    _d = lambda t: click.style(t, fg="cyan")  # noqa: E731
    _e = lambda t: click.style(t, fg="green")  # noqa: E731

    click.echo(_b("libephemeris configuration guide"))
    click.echo()

    # --- Data directory ---
    click.echo(_b("Data directory"))
    click.echo(f"  Current:  {data_dir}")
    click.echo("  Default:  ~/.libephemeris")
    click.echo(f"  Env var:  {_e('LIBEPHEMERIS_DATA_DIR')}")
    click.echo("  All downloaded files (DE kernels, planet centers, LEB, SPK)")
    click.echo("  are stored under this directory.")
    click.echo()

    # --- Precision tier ---
    try:
        from ..state import get_precision_tier

        current_tier = get_precision_tier()
    except Exception:
        current_tier = "medium"
    click.echo(_b("Precision tier"))
    click.echo(f"  Current:  {current_tier}")
    click.echo(
        f"  Env var:  {_e('LIBEPHEMERIS_PRECISION')}  (base | medium | extended)"
    )
    click.echo("  Python:   set_precision_tier('medium')")
    click.echo()
    click.echo(f"  {_d('base')}      de440s.bsp (~31 MB)   1849-2150 CE")
    click.echo(f"  {_d('medium')}    de440.bsp  (~114 MB)  1549-2650 CE  (default)")
    click.echo(f"  {_d('extended')}  de441.bsp  (~3.1 GB)  -13200 to +17191 CE")
    click.echo()

    # --- Ephemeris file ---
    eph_file = os.environ.get("LIBEPHEMERIS_EPHEMERIS", "")
    click.echo(_b("Ephemeris file (DE kernel)"))
    click.echo(f"  Current:  {eph_file or '(from tier)'}")
    click.echo(f"  Env var:  {_e('LIBEPHEMERIS_EPHEMERIS')}  (e.g. de441.bsp)")
    click.echo("  Python:   set_ephemeris_file('de441.bsp')")
    click.echo("  Usually set automatically by the precision tier.")
    click.echo(f"  Location: {data_dir}/")
    click.echo()

    # --- Calculation mode ---
    try:
        from ..state import get_calc_mode

        calc_mode = get_calc_mode()
    except Exception:
        calc_mode = "auto"
    click.echo(_b("Calculation mode"))
    click.echo(f"  Current:  {calc_mode}")
    click.echo(
        f"  Env var:  {_e('LIBEPHEMERIS_MODE')}  (auto | skyfield | leb | horizons)"
    )
    click.echo("  Python:   set_calc_mode('leb')")
    click.echo()
    click.echo(
        f"  {_d('auto')}      LEB -> Horizons -> Skyfield fallback chain (default)"
    )
    click.echo(f"  {_d('skyfield')}  Compute from DE kernel via Skyfield (real-time)")
    click.echo(f"  {_d('leb')}       Precomputed Chebyshev polynomials (~14x faster)")
    click.echo(f"  {_d('horizons')}  Query NASA JPL Horizons API (requires internet)")
    click.echo()

    # --- LEB file ---
    leb_path = os.environ.get("LIBEPHEMERIS_LEB", "")
    if not leb_path:
        try:
            from .. import state as _state

            leb_path = getattr(_state, "_LEB_FILE", None) or ""
        except Exception:
            pass
    click.echo(_b("LEB binary ephemeris"))
    click.echo(f"  Current:  {leb_path or '(none)'}")
    click.echo(f"  Env var:  {_e('LIBEPHEMERIS_LEB')}  (path to .leb or .leb2 file)")
    click.echo("  Python:   set_leb_file('path/to/ephemeris_medium.leb')")
    click.echo(f"  Location: {data_dir}/leb/")
    click.echo()
    click.echo("  Core distribution:")
    click.echo("    All tiers have an exact SHA-256-pinned core asset.")
    click.echo("    Existing LEB1 files remain loadable by explicit path.")
    click.echo()
    click.echo("  Canonical LEB2 group names (not all are published):")
    click.echo("    {tier}_core.leb2       core 14 bodies")
    click.echo("    {tier}_asteroids.leb2  Chiron, Ceres, Pallas, Juno, Vesta")
    click.echo(
        "    {tier}_exotics.leb2    31 registry bodies; extended has 23 (no NEAs)"
    )
    click.echo("    {tier}_apogee.leb2     OscuApog, IntpApog, IntpPerig")
    click.echo()
    click.echo("  Install a core: libephemeris download leb-medium")
    click.echo("  Install reviewed groups: libephemeris download leb2-medium")
    click.echo()

    # --- SPK cache ---
    try:
        from ..spk_auto import DEFAULT_AUTO_SPK_DIR
        from ..state import get_spk_cache_dir

        spk_dir = get_spk_cache_dir() or DEFAULT_AUTO_SPK_DIR
    except Exception:
        spk_dir = f"{data_dir}/spk"
    auto_spk = os.environ.get("LIBEPHEMERIS_AUTO_SPK", "")
    click.echo(_b("SPK asteroid cache"))
    click.echo(f"  Current:  {spk_dir}")
    click.echo(f"  Env var:  {_e('LIBEPHEMERIS_SPK_DIR')}  (override cache directory)")
    click.echo("  Python:   set_spk_cache_dir('/custom/path')")
    click.echo("  Auto-download from Horizons on demand:")
    click.echo(f"  Env var:  {_e('LIBEPHEMERIS_AUTO_SPK')}  (1/0, default: 1)")
    click.echo(f"  Current:  {auto_spk or '(default: enabled)'}")
    click.echo()

    # --- Planet centers ---
    click.echo(_b("Planet center corrections"))
    click.echo(f"  Location: {data_dir}/planet_centers_{{tier}}.bsp")
    click.echo("  Downloaded automatically with: libephemeris download <tier>")
    click.echo("  Provides sub-arcsecond precision for Jupiter-Pluto.")
    click.echo()

    # --- IERS ---
    click.echo(_b("IERS Earth orientation data"))
    iers_auto = os.environ.get("LIBEPHEMERIS_IERS_AUTO_DOWNLOAD", "")
    iers_dt = os.environ.get("LIBEPHEMERIS_IERS_DELTA_T", "")
    try:
        from ..iers_data import get_iers_cache_info

        iers_dir = str(get_iers_cache_info()["cache_dir"])
    except (ImportError, OSError, KeyError, ValueError):
        # Resolving the cache path touches the filesystem and the optional
        # config layer; anything else is a defect and should surface.
        iers_dir = f"{data_dir}/iers_cache"
    click.echo(f"  Location: {iers_dir}")
    click.echo("  Python:   set_iers_cache_dir('/custom/path')")
    _cur = f"  [current: {iers_auto}]" if iers_auto else ""
    click.echo(f"  Env var:  {_e('LIBEPHEMERIS_IERS_AUTO_DOWNLOAD')}  (1/0){_cur}")
    _cur = f"  [current: {iers_dt}]" if iers_dt else ""
    click.echo(
        f"  Env var:  {_e('LIBEPHEMERIS_IERS_DELTA_T')}  (1/0, use observed Delta T)"
        f"{_cur}"
    )
    click.echo("  Files:    finals2000A.data, leap_seconds.dat, deltat.data")
    click.echo()

    # --- Delta T model ---
    try:
        from ..state import get_delta_t_model

        deltat_model = get_delta_t_model()
    except Exception:
        deltat_model = "smh2016"
    click.echo(_b("Delta T model"))
    click.echo(f"  Current:  {deltat_model}")
    click.echo(
        f"  Env var:  {_e('LIBEPHEMERIS_DELTAT_MODEL')}  (smh2016 | espenak_meeus)"
    )
    click.echo("  Python:   set_delta_t_model('smh2016')")
    click.echo('  TOML:     deltat_model = "smh2016"')
    click.echo()

    # --- LEB mmap preloading ---
    from .._config_toml import get_bool, get_int

    mmap_preload = get_bool("mmap_preload")
    mmap_start = get_int("mmap_preload_start")
    mmap_end = get_int("mmap_preload_end")
    click.echo(_b("LEB mmap preloading"))
    click.echo(f"  Current:  {mmap_preload if mmap_preload is not None else False}")
    click.echo("  TOML:     mmap_preload = false")
    click.echo(f"            mmap_preload_start = {mmap_start or 1800}")
    click.echo(f"            mmap_preload_end = {mmap_end or 2200}")
    click.echo("  Preloads only the configured year range when a reader opens.")
    click.echo()

    # --- ASSIST ---
    click.echo(_b("ASSIST n-body data"))
    click.echo("  Location: ~/.libephemeris/assist/")
    click.echo("  Files:    linux_p1550p2650.440 (~98 MB)")
    click.echo("            sb441-n16.bsp (~616 MB)")
    click.echo("  Download: libephemeris download assist")
    click.echo("  Requires: pip install libephemeris[nbody]")
    click.echo()

    # --- .env file ---
    env_file_var = os.environ.get("LIBEPHEMERIS_ENV_FILE", "")
    click.echo(_b(".env file"))
    click.echo(f"  Env var:  {_e('LIBEPHEMERIS_ENV_FILE')}  (path to .env file)")
    click.echo("  Default:  ./.env then ~/.libephemeris/.env")
    if env_file_var:
        click.echo(f"  Current:  {env_file_var}")
    click.echo()

    # --- TOML config file ---
    from .._config_toml import get_config_path

    toml_path = get_config_path()
    click.echo(_b("TOML config file"))
    click.echo(f"  Env var:  {_e('LIBEPHEMERIS_CONFIG')}  (path to config file)")
    click.echo(
        "  Default:  ./libephemeris-config.toml then ~/.libephemeris/config.toml"
    )
    click.echo(f"  Current:  {toml_path or '(none found)'}")
    if toml_path:
        from .._config_toml import get_all

        cfg = get_all()
        if cfg:
            click.echo(f"  Values:   {len(cfg)} key(s) loaded")
            for k, v in cfg.items():
                click.echo(f"            {_d(k)} = {v!r}")
    click.echo("  Generate: libephemeris init")
    click.echo()

    # --- Logging ---
    log_level = os.environ.get("LIBEPHEMERIS_LOG_LEVEL", "")
    click.echo(_b("Logging"))
    click.echo(
        f"  Env var:  {_e('LIBEPHEMERIS_LOG_LEVEL')}  (DEBUG | INFO | WARNING | ERROR)"
    )
    click.echo(f"  Current:  {log_level or '(default: WARNING)'}")
    click.echo()

    # --- Other ---
    click.echo(_b("Other settings"))
    strict = os.environ.get("LIBEPHEMERIS_STRICT_PRECISION", "")
    click.echo(
        f"  {_e('LIBEPHEMERIS_STRICT_PRECISION')}  (1/0) "
        f"Raise errors instead of falling back  [{strict or 'default: 1'}]"
    )
    click.echo()

    click.echo("For current file status:  libephemeris status")
    click.echo("Full CLI reference:       see CLI.md")


# ---------------------------------------------------------------------------
# init command — interactive wizard to generate libephemeris-config.toml
# ---------------------------------------------------------------------------


@cli.command(
    short_help="Create a libephemeris-config.toml via interactive wizard.",
)
@click.option(
    "--output",
    "-o",
    default="libephemeris-config.toml",
    show_default=True,
    help="Output file path for the generated config.",
)
@click.option(
    "--force",
    "-f",
    is_flag=True,
    help="Overwrite existing config file without asking.",
)
@click.option(
    "--non-interactive",
    is_flag=True,
    help="Generate config with all defaults (no prompts).",
)
def init(output: str, force: bool, non_interactive: bool) -> None:
    """Generate a libephemeris-config.toml configuration file.

    An adaptive wizard that adjusts its questions based on your choices.
    Generates a fully-documented config file with all options (advanced
    settings are included as commented lines ready to uncomment).

    \b
    Quick start:
      libephemeris init                    Interactive wizard
      libephemeris init -o config.toml     Custom output path
      libephemeris init --non-interactive  All defaults, no prompts
    """
    from .init_wizard import run_wizard

    run_wizard(output=output, force=force, non_interactive=non_interactive)


# ---------------------------------------------------------------------------
# Backward compatibility aliases (download:base -> download base)
# ---------------------------------------------------------------------------
# These allow the old colon-separated syntax to still work for users
# who have scripts or muscle memory using the old argparse-based CLI.

_LEGACY_ALIASES = {
    "download:base": ["download", "base"],
    "download:medium": ["download", "medium"],
    "download:extended": ["download", "extended"],
    "download:leb:base": ["download", "leb-base"],
    "download:leb:medium": ["download", "leb-medium"],
    "download:leb:extended": ["download", "leb-extended"],
    "download:assist": ["download", "assist"],
}


def main(argv: list[str] | None = None) -> None:
    """Main entry point for the CLI.

    Supports both old (colon-separated) and new (space-separated) command syntax.
    """
    args = argv if argv is not None else sys.argv[1:]

    # Rewrite legacy colon-separated commands to space-separated equivalents
    if args:
        first = args[0]
        if first in _LEGACY_ALIASES:
            args = [*_LEGACY_ALIASES[first], *args[1:]]

    cli(args=args, standalone_mode=True)


if __name__ == "__main__":
    main()
