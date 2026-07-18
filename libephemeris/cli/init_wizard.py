# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Interactive configuration wizard for libephemeris.

Progressive-disclosure wizard that adapts questions and descriptions
based on user choices. Generates a fully-documented
libephemeris-config.toml file with download integration.

Usage (from CLI):
    libephemeris init                    Interactive wizard
    libephemeris init --non-interactive  All defaults, no prompts
    libephemeris init -o config.toml     Custom output path

Provenance:
    Project-authored user-interface and configuration serialization. It does
    not evaluate an ephemeris or define scientific constants. Date ranges and
    asset sizes shown to users are release metadata maintained with the tier
    manifests and do not enter calculations.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import click

from ..leb_groups import LEB2_GROUPS


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_DEFAULT_DATA_DIR = os.path.join("~", ".libephemeris")

_DE_KERNELS: Dict[str, tuple] = {
    "base": ("de440s.bsp", 31.0),
    "medium": ("de440.bsp", 114.0),
    "extended": ("de441.bsp", 3100.0),
}

_PC_SIZES: Dict[str, float] = {"base": 25.4, "medium": 191.25, "extended": 222.6}

_LEB2_GROUPS = list(LEB2_GROUPS)
# Only reviewed core assets currently have distributable sizes. Optional
# companion groups remain absent until their clean-room regeneration review.
_LEB2_SIZES: Dict[str, Dict[str, float]] = {
    "base": {"core": 10.8},
    "medium": {"core": 38.3},
    "extended": {"core": 334.9},
}

_TIER_RANGES: Dict[str, str] = {
    "base": "1849\u20132150 CE",
    "medium": "1549\u20132650 CE",
    "extended": "\u221213200 to +17191 CE",
}


# ---------------------------------------------------------------------------
# UI styling
# ---------------------------------------------------------------------------


def _resolve_leb_entry(path: str, precision: str) -> str:
    """Resolve a directory answer for leb_file to an actual LEB file.

    set_leb_file() expects a file; users naturally type the leb/
    directory here. Prefer the current reviewed LEB2 core, then the historical
    merged LEB1 for the chosen tier, then any .leb/.leb2 present; fall back to
    the original answer unchanged.
    """
    if not os.path.isdir(path):
        return path
    candidates = [
        os.path.join(path, f"{precision}_core.leb2"),
        os.path.join(path, f"ephemeris_{precision}.leb"),
    ]
    for cand in candidates:
        if os.path.isfile(cand):
            return cand
    try:
        for name in sorted(os.listdir(path)):
            if name.endswith((".leb", ".leb2")):
                return os.path.join(path, name)
    except OSError:
        pass
    return path


def _toml_str(value: str) -> str:
    """Escape a string for a TOML basic (double-quoted) value.

    Windows paths contain backslashes and values may contain quotes —
    emitting them raw produced invalid TOML escape sequences.
    """
    return value.replace("\\", "\\\\").replace('"', '\\"')


def _g(t: str) -> str:
    """Green bold — success, titles."""
    return click.style(t, fg="green", bold=True)


def _w(t: str) -> str:
    """Bold white — section headers."""
    return click.style(t, bold=True)


def _c(t: str) -> str:
    """Cyan — option values, labels."""
    return click.style(t, fg="cyan")


def _bc(t: str) -> str:
    """Bold cyan — prompt arrow."""
    return click.style(t, fg="cyan", bold=True)


def _d(t: str) -> str:
    """Dim — descriptions, hints."""
    return click.style(t, dim=True)


def _y(t: str) -> str:
    """Yellow — warnings, hints."""
    return click.style(t, fg="yellow")


def _sep() -> None:
    """Print a horizontal separator line."""
    click.echo(_d("\u2500" * 50))


def _size(mb: float) -> str:
    """Format a size in MB to a human-readable string."""
    if mb >= 1024:
        return f"{mb / 1024:.1f} GB"
    return f"{mb:.0f} MB"


# ---------------------------------------------------------------------------
# Filesystem TAB completion
# ---------------------------------------------------------------------------


def _setup_path_completion() -> None:
    """Enable filesystem TAB completion for path prompts."""
    try:
        import readline
        import glob as _glob_mod

        def _completer(text: str, state: int) -> Optional[str]:
            expanded = os.path.expanduser(text)
            if os.path.isdir(expanded) and not expanded.endswith(os.sep):
                expanded += os.sep
            matches = _glob_mod.glob(expanded + "*")
            matches = [m + os.sep if os.path.isdir(m) else m for m in matches]
            if text.startswith("~"):
                home = os.path.expanduser("~")
                matches = [m.replace(home, "~", 1) for m in matches]
            try:
                return matches[state]
            except IndexError:
                return None

        readline.set_completer(_completer)
        readline.set_completer_delims(" \t\n")
        doc = getattr(readline, "__doc__", "") or ""
        if "libedit" in doc:
            readline.parse_and_bind("bind ^I rl_complete")
        else:
            readline.parse_and_bind("tab: complete")
    except ImportError:
        pass


def _teardown_path_completion() -> None:
    """Disable the custom completer."""
    try:
        import readline

        readline.set_completer(None)
    except ImportError:
        pass


# ---------------------------------------------------------------------------
# Arrow-key interactive selector
# ---------------------------------------------------------------------------


def _select(
    options: List[Dict[str, str]],
    default: int = 0,
) -> str:
    """Arrow-key selector for multiple-choice prompts.

    Each *option* must have a ``value`` key and optionally ``description``.
    Returns the selected ``value``.  Falls back to ``click.prompt`` when
    stdin is not a TTY (piped input, CI, Windows without termios).
    """
    can_interact = sys.stdin.isatty() and sys.stdout.isatty()
    if can_interact:
        try:
            __import__("termios")
            __import__("tty")
        except ImportError:
            can_interact = False

    # Fallback: show options list + typed input
    if not can_interact:
        values = [o["value"] for o in options]
        for o in options:
            desc = ""
            if o.get("description"):
                desc = "  " + _d(o["description"])
            click.echo(f"    {_c(o['value'])}{desc}")
        return click.prompt(
            f"  {_bc('>')}",
            type=click.Choice(values, case_sensitive=False),
            default=values[default],
            show_choices=False,
        )

    # Interactive arrow-key selector
    import termios
    import tty

    selected = default
    n = len(options)
    fd = sys.stdin.fileno()
    old = termios.tcgetattr(fd)
    max_w = max(len(o["value"]) for o in options)

    def _render() -> None:
        for i, opt in enumerate(options):
            padded = opt["value"].ljust(max_w)
            if i == selected:
                arrow = _bc("> ")
                label = _bc(padded)
            else:
                arrow = "  "
                label = _c(padded)
            desc = ""
            if opt.get("description"):
                desc = "  " + _d(opt["description"])
            click.echo(f"    {arrow}{label}{desc}")

    click.echo(_d("  (arrow keys to move, Enter to confirm)"))
    click.echo()
    _render()

    try:
        tty.setcbreak(fd)
        while True:
            ch = os.read(fd, 1)
            if ch in (b"\r", b"\n"):
                break
            if ch == b"\x03":  # Ctrl+C
                raise KeyboardInterrupt
            if ch == b"\x1b":
                seq = os.read(fd, 2)
                if seq == b"[A":  # Up
                    selected = (selected - 1) % n
                elif seq == b"[B":  # Down
                    selected = (selected + 1) % n
                else:
                    continue
                sys.stdout.write(f"\x1b[{n}A\x1b[0J")
                sys.stdout.flush()
                _render()
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old)

    # Collapse selector to a single result line
    sys.stdout.write(f"\x1b[{n + 2}A\x1b[0J")
    sys.stdout.flush()
    click.echo(f"  {_bc('>')} {_c(options[selected]['value'])}")

    return options[selected]["value"]


# ---------------------------------------------------------------------------
# Required-files logic
# ---------------------------------------------------------------------------


def _get_required_files(tier: str, mode: str) -> List[Dict[str, Any]]:
    """Determine which data files are needed for *tier* + *mode*."""
    files: List[Dict[str, Any]] = []

    # LEB2 (primary engine for auto / leb)
    if mode in ("auto", "leb"):
        for group in _LEB2_GROUPS:
            size = _LEB2_SIZES.get(tier, {}).get(group, 0)
            files.append(
                {
                    "filename": f"{tier}_{group}.leb2",
                    "size_mb": size,
                    "category": "LEB2",
                    # Every tier has a pinned core; optional companions remain
                    # local-generation-only until separately reviewed.
                    "status": "required" if group == "core" else "optional",
                    "subdir": "leb",
                }
            )

    # DE kernel
    de_file, de_size = _DE_KERNELS[tier]
    if mode == "auto":
        files.append(
            {
                "filename": de_file,
                "size_mb": de_size,
                "category": "DE kernel",
                "status": "fallback",
                "subdir": "",
            }
        )
    elif mode == "skyfield":
        files.append(
            {
                "filename": de_file,
                "size_mb": de_size,
                "category": "DE kernel",
                "status": "required",
                "subdir": "",
            }
        )

    # Planet centers
    pc_size = _PC_SIZES.get(tier, 0)
    files.append(
        {
            "filename": f"planet_centers_{tier}.bsp",
            "size_mb": pc_size,
            "category": "Planet centers",
            "status": (
                "fallback"
                if mode == "auto"
                else "optional"
                if mode == "horizons"
                else "required"
            ),
            "subdir": "",
        }
    )

    return files


def _check_files(files: List[Dict[str, Any]], data_dir: str) -> List[Dict[str, Any]]:
    """Annotate each file dict with ``exists`` and ``path``."""
    data_path = Path(os.path.expanduser(data_dir))
    results: List[Dict[str, Any]] = []
    for f in files:
        if f["subdir"]:
            path = data_path / f["subdir"] / f["filename"]
        else:
            path = data_path / f["filename"]
        results.append({**f, "exists": path.exists(), "path": str(path)})
    return results


# ---------------------------------------------------------------------------
# TOML generation
# ---------------------------------------------------------------------------


def _generate_toml(config: Dict[str, Any]) -> str:
    """Create a complete, well-documented libephemeris-config.toml."""
    lines = [
        "# libephemeris configuration",
        "# Generated by: libephemeris init",
        "#",
        "# Resolution order (highest to lowest priority):",
        "#   1. set_*() function calls (programmatic)",
        "#   2. Environment variables / .env file",
        "#   3. This file (libephemeris-config.toml)",
        "#   4. Built-in defaults",
        "#",
        "# Search order for this file:",
        "#   1. LIBEPHEMERIS_CONFIG env var (explicit path)",
        "#   2. ./libephemeris-config.toml (project directory)",
        "#   3. ~/.libephemeris/config.toml (global)",
        "",
        "[libephemeris]",
        "",
        "# --- Core ---",
    ]

    prec = config.get("precision", "medium")
    mode = config.get("mode", "auto")
    lines.append(f'precision = "{_toml_str(prec)}"')
    lines.append(f'mode = "{_toml_str(mode)}"')
    network_policy = config.get("network_policy", "sealed" if mode == "leb" else "auto")
    lines.append(f'network_policy = "{_toml_str(network_policy)}"')
    lines.append("")

    # --- Minor Bodies ---
    lines.append("# --- Minor Bodies ---")
    auto_spk = config.get("auto_spk", True)
    strict = config.get("strict_precision", auto_spk)
    lines.append(f"auto_spk = {'true' if auto_spk else 'false'}")
    lines.append(f"strict_precision = {'true' if strict else 'false'}")
    lines.append("")

    # --- Time Precision ---
    lines.append("# --- Time Precision ---")
    iers_auto = config.get("iers_auto_download", False)
    iers_dt = config.get("iers_delta_t", iers_auto)
    lines.append(f"iers_auto_download = {'true' if iers_auto else 'false'}")
    lines.append(f"iers_delta_t = {'true' if iers_dt else 'false'}")
    deltat_model = config.get("deltat_model", "smh2016")
    lines.append(f'deltat_model = "{_toml_str(deltat_model)}"')
    lines.append("")

    # --- Advanced ---
    lines.append("# --- Advanced (uncomment to customize) ---")

    log_level = config.get("log_level")
    if log_level and log_level != "WARNING":
        lines.append(f'log_level = "{_toml_str(log_level)}"')
    else:
        lines.append('# log_level = "WARNING"')

    data_dir = config.get("data_dir")
    if data_dir:
        lines.append(f'data_dir = "{_toml_str(data_dir)}"')
    else:
        lines.append(f'# data_dir = "{_DEFAULT_DATA_DIR}"')

    leb_file = config.get("leb_file")
    if leb_file:
        lines.append(f'leb_file = "{_toml_str(leb_file)}"')
    else:
        lines.append('# leb_file = ""')

    spk_dir = config.get("spk_dir")
    if spk_dir:
        lines.append(f'spk_dir = "{_toml_str(spk_dir)}"')
    else:
        lines.append('# spk_dir = ""')

    ephemeris = config.get("ephemeris")
    if ephemeris:
        lines.append(f'ephemeris = "{_toml_str(ephemeris)}"')
    else:
        lines.append('# ephemeris = ""')

    mmap_preload = bool(config.get("mmap_preload", False))
    mmap_preload_start = int(config.get("mmap_preload_start", 1800))
    mmap_preload_end = int(config.get("mmap_preload_end", 2200))
    lines.append(f"mmap_preload = {'true' if mmap_preload else 'false'}")
    lines.append(f"mmap_preload_start = {mmap_preload_start}")
    lines.append(f"mmap_preload_end = {mmap_preload_end}")

    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------


def _run_downloads(tier: str, mode: str) -> None:
    """Download files required for *tier* + *mode*."""
    try:
        if mode in ("auto", "leb"):
            click.echo()
            click.echo(_d("  Downloading LEB2 files..."))
            from ..download import download_leb2_for_tier

            download_leb2_for_tier(
                tier_name=tier,
                force=False,
                show_progress=True,
                quiet=False,
                activate=True,
            )
            click.echo()

        if mode in ("auto", "skyfield"):
            click.echo(_d("  Downloading DE kernel + planet centers + SPK..."))
            from ..download import download_for_tier

            download_for_tier(
                tier_name=tier,
                force=False,
                show_progress=True,
                quiet=False,
            )
        elif mode in ("leb", "horizons"):
            click.echo(_d("  Downloading planet centers..."))
            from ..download import _download_planet_centers_for_tier

            _download_planet_centers_for_tier(
                tier_name=tier,
                force=False,
                show_progress=True,
                quiet=False,
            )

        click.echo()
        click.echo(f"  {_g('Downloads complete!')}")
    except KeyboardInterrupt:
        click.echo("\n  Download cancelled.")
    except (OSError, ValueError, RuntimeError) as e:
        click.echo(f"\n  {_y(f'Download error: {e}')}")
        click.echo(_d("  You can retry with: libephemeris download auto"))


def _show_manual_commands(tier: str, mode: str) -> None:
    """Print the manual download commands for a given mode."""
    cmds: List[str] = []
    if mode in ("auto", "leb"):
        cmds.append(f"libephemeris download leb2-{tier}")
    if mode in ("auto", "skyfield"):
        cmds.append(f"libephemeris download {tier}")
    if mode in ("leb", "horizons"):
        cmds.append(f"libephemeris download {tier}")
    if not cmds:
        cmds.append(f"libephemeris download {tier}")
    # deduplicate
    seen = set()
    for cmd in cmds:
        if cmd not in seen:
            seen.add(cmd)
            click.echo(f"    {_w(cmd)}")


# ---------------------------------------------------------------------------
# Post-wizard: file status + download offer
# ---------------------------------------------------------------------------


def _post_wizard(
    config: Dict[str, Any],
    data_dir: str,
    output: str,
    interactive: bool = True,
) -> None:
    """Show file status, offer downloads, print next steps."""
    tier = config["precision"]
    mode = config["mode"]

    files = _get_required_files(tier, mode)
    checked = _check_files(files, data_dir)

    click.echo()
    click.echo(_w("  Files needed for your configuration:"))
    click.echo()

    missing_required_count = 0
    missing_required_size = 0.0
    missing_fallback_count = 0
    missing_fallback_size = 0.0

    def _status_suffix(status: str) -> str:
        if status == "fallback":
            return _d(" (fallback)")
        if status == "optional":
            return _d(" (optional)")
        return ""

    for f in checked:
        s = _size(f["size_mb"])
        suffix = _status_suffix(f["status"])
        if f["exists"]:
            marker = click.style("[OK]", fg="green")
            click.echo(f"    {marker} {f['filename']:<30s} {s:>8s}{suffix}")
        else:
            marker = click.style("[--]", fg="yellow")
            click.echo(f"    {marker} {f['filename']:<30s} {s:>8s}{suffix}")
            if f["status"] == "required":
                missing_required_count += 1
                missing_required_size += f["size_mb"]
            elif f["status"] == "fallback":
                missing_fallback_count += 1
                missing_fallback_size += f["size_mb"]

    click.echo()

    if missing_required_count == 0 and missing_fallback_count == 0:
        click.echo(f"  {_g('All files for this mode are present!')}")
    else:
        if missing_required_count > 0:
            click.echo(
                f"  {_y(f'Missing required: {missing_required_count} file(s) (~{_size(missing_required_size)})')}"
            )
        if missing_fallback_count > 0:
            click.echo(
                f"  {_y(f'Missing fallback: {missing_fallback_count} file(s) (~{_size(missing_fallback_size)})')}"
            )
        if interactive:
            click.echo()
            if click.confirm(f"  {_bc('>')} Download missing files now?", default=True):
                _run_downloads(tier, mode)
            else:
                click.echo()
                click.echo(f"  {_y('Skipped.')} {_d('Run these commands when ready:')}")
                _show_manual_commands(tier, mode)
        else:
            click.echo()
            click.echo(_d("  To download, run:"))
            _show_manual_commands(tier, mode)

    click.echo()
    _sep()
    click.echo()
    click.echo(_w("  Next steps:"))
    click.echo(_d(f"  \u2022 Commit {output} to version control"))
    click.echo(_d("  \u2022 Run 'libephemeris status' to verify your setup"))
    click.echo(
        _d("  \u2022 Run 'libephemeris download auto' to fetch missing files later")
    )
    click.echo()


# ---------------------------------------------------------------------------
# Main wizard entry point
# ---------------------------------------------------------------------------


def run_wizard(
    output: str = "libephemeris-config.toml",
    force: bool = False,
    non_interactive: bool = False,
) -> None:
    """Run the interactive configuration wizard.

    Creates a ``libephemeris-config.toml`` with adaptive questions,
    filesystem TAB-completion for paths, and integrated downloads.
    """
    out_path = Path(output)

    # --- Guard: file already exists ---
    if out_path.exists() and not force:
        if non_interactive:
            click.echo(f"Error: {output} already exists. Use --force to overwrite.")
            sys.exit(1)
        if not click.confirm(f"{output} already exists. Overwrite?"):
            click.echo("Aborted.")
            return

    # --- Non-interactive: all defaults, write immediately ---
    if non_interactive:
        config: Dict[str, Any] = {
            "precision": "medium",
            "mode": "auto",
            "network_policy": "auto",
            "auto_spk": True,
            "strict_precision": True,
            "iers_auto_download": False,
            "iers_delta_t": False,
        }
        _write_and_finish(out_path, config, output, interactive=False)
        return

    # === INTERACTIVE WIZARD ===
    click.echo()
    click.echo(_g("libephemeris configuration wizard"))
    _sep()
    click.echo()

    # ── Quick setup ──────────────────────────────────────────────
    click.echo(_w("Quick setup"))
    click.echo(_d("  Recommended: medium precision, auto mode, ~/.libephemeris"))
    click.echo(_d("  You can always edit the generated file later."))
    use_defaults = click.confirm(f"  {_bc('>')} Use defaults?", default=True)

    if use_defaults:
        config = {
            "precision": "medium",
            "mode": "auto",
            "network_policy": "auto",
            "auto_spk": True,
            "strict_precision": True,
            "iers_auto_download": False,
            "iers_delta_t": False,
        }
        click.echo()
        _write_and_finish(out_path, config, output, data_dir=_DEFAULT_DATA_DIR)
        return

    config = {}

    # ── 1. Precision tier ────────────────────────────────────────
    click.echo()
    _sep()
    click.echo()
    click.echo(_w("Precision tier"))
    tier = _select(
        [
            {
                "value": "base",
                "description": f"~{_size(_DE_KERNELS['base'][1])} {_TIER_RANGES['base']}",
            },
            {
                "value": "medium",
                "description": f"~{_size(_DE_KERNELS['medium'][1])} {_TIER_RANGES['medium']}",
            },
            {
                "value": "extended",
                "description": f"~{_size(_DE_KERNELS['extended'][1])} {_TIER_RANGES['extended']}",
            },
        ],
        default=1,
    )
    config["precision"] = tier

    # ── 2. Calculation mode ──────────────────────────────────────
    click.echo()
    click.echo(_w("Calculation mode"))
    mode = _select(
        [
            {
                "value": "auto",
                "description": "LEB -> Horizons -> Skyfield (recommended)",
            },
            {"value": "leb", "description": "Precomputed polynomials (~14x faster)"},
            {"value": "skyfield", "description": "DE kernel via Skyfield"},
            {"value": "horizons", "description": "NASA JPL API (requires internet)"},
        ],
        default=0,
    )
    config["mode"] = mode
    config["network_policy"] = "sealed" if mode == "leb" else "auto"

    # ── 3. Minor bodies (adaptive) ───────────────────────────────
    click.echo()
    click.echo(_w("Minor bodies"))

    if mode in ("auto", "leb"):
        click.echo(_d("  LEB files cover the body groups present in your installed"))
        click.echo(_d("  tier companions. Enable auto-download for any *other* minor"))
        click.echo(_d("  body you request (e.g. Eros, Sedna)?"))
        auto_spk_default = mode == "auto"
    elif mode == "skyfield":
        click.echo(_d("  The DE kernel only contains the major planets and the Moon."))
        click.echo(_d("  Enable auto-download of SPK kernels from JPL Horizons for"))
        click.echo(_d("  asteroids, centaurs, and TNOs when requested?"))
        auto_spk_default = True
    else:  # horizons
        click.echo(_d("  Horizons queries NASA live for all supported bodies."))
        click.echo(_d("  SPK kernels are only useful as an offline fallback."))
        auto_spk_default = False

    auto_spk = click.confirm(f"  {_bc('>')} auto_spk", default=auto_spk_default)
    config["auto_spk"] = auto_spk

    if auto_spk:
        click.echo(_d("  Require an ephemeris-grade source for mapped minor bodies?"))
        click.echo(_d("  If disabled, a declared Keplerian approximation may be used."))
        strict = click.confirm(f"  {_bc('>')} strict_precision", default=True)
        config["strict_precision"] = strict
    else:
        # Without auto-download, strict precision raises SPKRequiredError for
        # mapped minor bodies unless local SPK kernels are already
        # provisioned, so the safe default is the declared Keplerian
        # approximation; pre-provisioned offline deployments may still
        # opt in.
        click.echo(
            _d("  Without auto-download, strict precision raises SPKRequiredError")
        )
        click.echo(
            _d("  for mapped minor bodies unless local SPK kernels are provisioned.")
        )
        strict = click.confirm(f"  {_bc('>')} strict_precision", default=False)
        config["strict_precision"] = strict

    # ── 4. Time precision (adaptive) ─────────────────────────────
    click.echo()
    click.echo(_w("Time precision"))
    click.echo(_d("  All modes convert UTC to astronomical time (TT) using Delta T."))
    click.echo(_d("  Built-in tables are accurate through ~2020. Download latest"))
    click.echo(_d("  IERS data (~3 MB) for best precision on recent/future dates?"))
    if mode == "horizons":
        click.echo(
            _d("  Recommended in horizons mode since you\u2019re already online.")
        )

    iers_default = mode == "horizons"
    iers = click.confirm(f"  {_bc('>')} download IERS data", default=iers_default)
    config["iers_auto_download"] = iers
    config["iers_delta_t"] = iers  # derived: download implies use

    if not iers:
        click.echo(_d("  (iers_delta_t skipped \u2014 requires IERS data)"))

    # ── 5. Data directory ────────────────────────────────────────
    click.echo()
    click.echo(_w("Data directory"))
    click.echo(f"  {_y('(press Tab to autocomplete paths)')}")
    click.echo(_d("  Base folder for all downloads (DE kernels, LEB, SPK, IERS)."))

    _setup_path_completion()
    try:
        data_dir = click.prompt(
            f"  {_bc('>')}",
            default=_DEFAULT_DATA_DIR,
        ).strip()
    finally:
        _teardown_path_completion()

    if data_dir and data_dir != _DEFAULT_DATA_DIR:
        config["data_dir"] = data_dir

    effective_data_dir = data_dir or _DEFAULT_DATA_DIR

    # ── 6. Sub-paths (optional expansion) ────────────────────────
    click.echo(_d("  Sub-paths are derived automatically from the data directory."))
    customize = click.confirm(f"  {_bc('>')} Customize sub-paths?", default=False)

    if customize:
        expanded = os.path.expanduser(effective_data_dir)

        _setup_path_completion()
        try:
            # Show sub-paths relevant to the chosen mode
            if mode in ("auto", "leb"):
                click.echo(
                    _d("  Directory containing .leb/.leb2 binary ephemeris files.")
                )
                leb = click.prompt(
                    f"  {_bc('>')} leb_file",
                    default=os.path.join(expanded, "leb"),
                ).strip()
                if leb:
                    config["leb_file"] = _resolve_leb_entry(
                        leb, config.get("precision", "medium")
                    )

            click.echo(_d("  Cache directory for downloaded minor body SPK kernels."))
            spk = click.prompt(
                f"  {_bc('>')} spk_dir",
                default=os.path.join(expanded, "spk"),
            ).strip()
            if spk:
                config["spk_dir"] = spk

            if mode in ("auto", "skyfield"):
                click.echo(_d("  DE kernel filename override (leave empty for auto)."))
                eph = click.prompt(
                    f"  {_bc('>')} ephemeris",
                    default="",
                    show_default=False,
                ).strip()
                if eph:
                    config["ephemeris"] = eph

            # Show paths that are less relevant but still available
            if mode == "skyfield":
                click.echo(_d("  Path to .leb/.leb2 file (not needed in skyfield mode"))
                click.echo(
                    _d("  unless you plan to switch to auto or leb mode later).")
                )
                leb = click.prompt(
                    f"  {_bc('>')} leb_file",
                    default="",
                    show_default=False,
                ).strip()
                if leb:
                    config["leb_file"] = _resolve_leb_entry(
                        leb, config.get("precision", "medium")
                    )

            if mode == "horizons":
                click.echo(_d("  DE kernel and LEB are not needed in horizons mode"))
                click.echo(_d("  unless you plan to switch modes later."))
        finally:
            _teardown_path_completion()

    # ── Write and finish ─────────────────────────────────────────
    click.echo()
    _write_and_finish(out_path, config, output, data_dir=effective_data_dir)


# ---------------------------------------------------------------------------
# Write config + post-wizard flow
# ---------------------------------------------------------------------------


def _write_and_finish(
    out_path: Path,
    config: Dict[str, Any],
    output: str,
    data_dir: str = _DEFAULT_DATA_DIR,
    interactive: bool = True,
) -> None:
    """Write the TOML file and run the post-wizard flow."""
    content = _generate_toml(config)
    try:
        out_path.write_text(content, encoding="utf-8")
    except OSError as e:
        click.echo(f"Error writing {output}: {e}", err=True)
        sys.exit(1)

    _sep()
    click.echo()
    click.echo(f"  {_g('Created')} {output}")

    _post_wizard(config, data_dir, output, interactive=interactive)
