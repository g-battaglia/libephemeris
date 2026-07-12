# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Development status report for the `leph` CLI."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import click

from .. import __version__
from ..constants import SPK_AUTO_DOWNLOAD_BLOCKED
from ..leb_groups import LEB1_GROUPS, LEB2_GROUPS

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_LEB1_DIR = _PROJECT_ROOT / "data" / "leb"
_LEB2_DIR = _PROJECT_ROOT / "data" / "leb2"
_LEB2_GROUPS = list(LEB2_GROUPS)
# Per-tier LEB1 partial paths — derived from the canonical group partition.
_LEB1_PARTIALS = {
    tier: {g: _LEB1_DIR / f"ephemeris_{tier}_{g}.leb" for g in LEB1_GROUPS}
    for tier in ("base", "medium", "extended")
}
_PLANET_CENTER_SOURCE_FILES = {
    "base": ["jup204.bsp", "sat319.bsp", "ura083.bsp", "nep050.bsp", "plu017.bsp"],
    "medium": [
        "jup365.bsp",
        "sat441.bsp",
        "ura184_part-3.bsp",
        "nep105.bsp",
        "plu060.bsp",
    ],
    "extended": [
        "jup365.bsp",
        "sat441xl_part-1.bsp",
        "sat441xl_part-2.bsp",
        "ura111xl-799.bsp",
        "nep097xl-899.bsp",
        "plu060.bsp",
    ],
}
_EXTENDED_SPK_RANGE = ("1600-01-01", "2500-01-01")
_EXPECTED_SPK_FALLBACKS = SPK_AUTO_DOWNLOAD_BLOCKED


def _display_body_name(name: str) -> str:
    """Normalize CLI-facing body names from Horizons-style identifiers."""
    return name.rstrip(";")


def _file_entry(path: Path) -> dict[str, Any]:
    """Return a serializable file entry with existence and size metadata."""
    exists = path.exists()
    return {
        "exists": exists,
        "path": str(path),
        "size": path.stat().st_size if exists else None,
        "name": path.name,
    }


def _status_from_counts(present: int, total: int, *, partial: int = 0) -> str:
    """Collapse per-item counts into `ok`, `partial`, or `missing`."""
    if total == 0:
        return "missing"
    if present == total and partial == 0:
        return "ok"
    if present == 0 and partial == 0:
        return "missing"
    return "partial"


def _resolve_planet_centers_path(tier: str, data_dir: Path) -> Path:
    """Resolve the preferred location of a generated planet-centers artifact."""
    filename = f"planet_centers_{tier}.bsp"
    data_path = data_dir / filename
    if data_path.exists():
        return data_path
    workspace_path = _PROJECT_ROOT / filename
    if workspace_path.exists():
        return workspace_path
    return data_path


def _scan_best_spk_files(cache_dir: Path) -> dict[int, dict[str, Any]]:
    """Map each known minor body to the widest local SPK file available."""
    from ..constants import SPK_BODY_NAME_MAP
    from ..spk import _get_body_name, _get_spk_targets, get_spk_coverage

    naif_to_body: dict[int, tuple[int, str]] = {}
    for ipl, (horizons_id, naif_id) in SPK_BODY_NAME_MAP.items():
        body_name = _display_body_name(_get_body_name(ipl) or horizons_id)
        naif_to_body[naif_id] = (ipl, body_name)
        asteroid_number = naif_id - 2000000
        if asteroid_number > 0:
            naif_to_body.setdefault(asteroid_number + 20000000, (ipl, body_name))

    best_by_ipl: dict[int, dict[str, Any]] = {}
    if not cache_dir.exists():
        return best_by_ipl

    for path in sorted(cache_dir.glob("*.bsp")):
        targets = _get_spk_targets(str(path))
        if not targets:
            continue

        coverage = get_spk_coverage(str(path))
        span = (coverage[1] - coverage[0]) if coverage else -1.0

        for target in targets:
            if target not in naif_to_body:
                continue
            ipl, body_name = naif_to_body[target]
            current = best_by_ipl.get(ipl)
            if current is None or span > current["span"]:
                best_by_ipl[ipl] = {
                    "body": body_name,
                    "path": str(path),
                    "coverage": coverage,
                    "span": span,
                }

    return best_by_ipl


def _collect_spk_status(data_dir: Path) -> dict[str, Any]:
    """Collect development status for DE kernels and asteroid SPK caches."""
    from ..constants import SPK_BODY_NAME_MAP
    from ..spk import _get_body_name
    from ..spk_auto import DEFAULT_AUTO_SPK_DIR, _iso_to_jd
    from ..state import get_spk_cache_dir, get_spk_date_range_for_tier

    cache_dir = Path(get_spk_cache_dir() or DEFAULT_AUTO_SPK_DIR)
    best_by_ipl = _scan_best_spk_files(cache_dir)
    de_files = {
        "base": "de440s.bsp",
        "medium": "de440.bsp",
        "extended": "de441.bsp",
    }
    tiers: dict[str, Any] = {}

    for tier in ["base", "medium", "extended"]:
        if tier == "extended":
            start_date, end_date = _EXTENDED_SPK_RANGE
        else:
            start_date, end_date = get_spk_date_range_for_tier(tier)
        start_jd = _iso_to_jd(start_date)
        end_jd = _iso_to_jd(end_date)

        bodies = []
        ready = 0
        partial = 0
        missing = 0
        expected_fallback = 0
        for ipl, (horizons_id, _) in sorted(
            SPK_BODY_NAME_MAP.items(),
            key=lambda item: _display_body_name(
                _get_body_name(item[0]) or item[1][0]
            ).lower(),
        ):
            body_name = _display_body_name(_get_body_name(ipl) or horizons_id)
            best = best_by_ipl.get(ipl)
            if best is None:
                if ipl in _EXPECTED_SPK_FALLBACKS:
                    body_entry = {
                        "name": body_name,
                        "status": "expected_fallback",
                        "path": None,
                        "coverage": None,
                        "reason": _EXPECTED_SPK_FALLBACKS[ipl],
                    }
                    expected_fallback += 1
                else:
                    body_entry = {
                        "name": body_name,
                        "status": "missing",
                        "path": None,
                        "coverage": None,
                    }
                    missing += 1
            else:
                coverage = best["coverage"]
                is_ready = (
                    coverage is not None
                    and coverage[0] <= start_jd
                    and coverage[1] >= end_jd
                )
                body_entry = {
                    "name": body_name,
                    "status": "ready" if is_ready else "partial",
                    "path": best["path"],
                    "coverage": coverage,
                }
                if is_ready:
                    ready += 1
                else:
                    partial += 1
            bodies.append(body_entry)

        de_entry = _file_entry(data_dir / de_files[tier])
        required_total = len(bodies) - expected_fallback
        tier_status = _status_from_counts(ready, required_total, partial=partial)
        if not de_entry["exists"] and tier_status == "ok":
            tier_status = "partial"
        elif not de_entry["exists"] and ready == 0 and partial == 0:
            tier_status = "missing"

        tiers[tier] = {
            "status": tier_status,
            "de_kernel": de_entry,
            "range": {"start": start_date, "end": end_date},
            "cache_dir": str(cache_dir),
            "total_bodies": len(bodies),
            "required_bodies": required_total,
            "ready_bodies": ready,
            "partial_bodies": partial,
            "missing_bodies": missing,
            "expected_fallback_bodies": expected_fallback,
            "bodies": bodies,
        }

    cache_files = []
    if cache_dir.exists():
        for path in sorted(cache_dir.glob("*.bsp")):
            cache_files.append(_file_entry(path))

    return {"cache_dir": str(cache_dir), "tiers": tiers, "cache_files": cache_files}


def _collect_iers_status() -> dict[str, Any]:
    """Collect development status for IERS support files."""
    from ..iers_data import get_iers_cache_info

    cache_info = get_iers_cache_info()
    cache_dir = Path(str(cache_info.get("cache_dir", "")))
    files = {
        "finals2000A.data": {
            "exists": bool(cache_info.get("finals_exists")),
            "path": str(cache_dir / "finals2000A.data"),
            "age_days": cache_info.get("finals_age_days"),
        },
        "leap_seconds.dat": {
            "exists": bool(cache_info.get("leap_seconds_exists")),
            "path": str(cache_dir / "leap_seconds.dat"),
            "age_days": cache_info.get("leap_seconds_age_days"),
        },
        "deltat.data": {
            "exists": bool(cache_info.get("delta_t_exists")),
            "path": str(cache_dir / "deltat.data"),
            "age_days": cache_info.get("delta_t_age_days"),
        },
    }
    present = sum(1 for info in files.values() if info["exists"])
    return {
        "status": _status_from_counts(present, len(files)),
        "cache_dir": str(cache_dir),
        "present": present,
        "total": len(files),
        "files": files,
    }


def _collect_planet_center_source_status() -> dict[str, Any]:
    """Collect status for downloaded NAIF source kernels used by generation."""
    cache_dir = Path.home() / ".libephemeris" / "tmp" / "planet_centers"
    shared_file = cache_dir / "naif0012.tls"
    shared = {
        "status": "ok" if shared_file.exists() else "missing",
        "present": 1 if shared_file.exists() else 0,
        "total": 1,
        "files": {shared_file.name: _file_entry(shared_file)},
    }

    tiers: dict[str, Any] = {}
    for tier, filenames in _PLANET_CENTER_SOURCE_FILES.items():
        files = {name: _file_entry(cache_dir / name) for name in filenames}
        present = sum(1 for info in files.values() if info["exists"])
        tiers[tier] = {
            "status": _status_from_counts(present, len(files)),
            "present": present,
            "total": len(files),
            "files": files,
        }

    overall_present = shared["present"] + sum(t["present"] for t in tiers.values())
    overall_total = shared["total"] + sum(t["total"] for t in tiers.values())
    overall_status = _status_from_counts(
        overall_present,
        overall_total,
        partial=sum(1 for t in tiers.values() if t["status"] == "partial"),
    )
    return {
        "status": overall_status,
        "cache_dir": str(cache_dir),
        "shared": shared,
        "tiers": tiers,
    }


def _collect_assist_status() -> dict[str, Any]:
    """Collect status for local ASSIST data files."""
    assist_dir = Path.home() / ".libephemeris" / "assist"
    files = {
        "linux_p1550p2650.440": _file_entry(assist_dir / "linux_p1550p2650.440"),
        "sb441-n16.bsp": _file_entry(assist_dir / "sb441-n16.bsp"),
    }
    present = sum(1 for info in files.values() if info["exists"])
    return {
        "status": _status_from_counts(present, len(files)),
        "directory": str(assist_dir),
        "present": present,
        "total": len(files),
        "files": files,
    }


def _collect_leb1_status() -> dict[str, Any]:
    """Collect status for generated LEB1 artifacts in the repository."""
    tiers: dict[str, Any] = {}
    for tier in ["base", "medium", "extended"]:
        final_file = _file_entry(_LEB1_DIR / f"ephemeris_{tier}.leb")
        partials = {
            name: _file_entry(path) for name, path in _LEB1_PARTIALS[tier].items()
        }
        partial_present = sum(1 for info in partials.values() if info["exists"])
        status = (
            "ok"
            if final_file["exists"]
            else _status_from_counts(
                partial_present,
                len(partials),
            )
        )
        tiers[tier] = {
            "status": status,
            "final": final_file,
            "partial_present": partial_present,
            "partial_total": len(partials),
            "partials": partials,
        }
    return tiers


def _collect_leb2_status(data_dir: Path) -> dict[str, Any]:
    """Collect status for LEB2 group files in both repo and user data dir.

    Checks two locations for each group file:
    1. ``data/leb2/`` in the project repo (freshly generated artifacts)
    2. ``~/.libephemeris/leb/`` user data dir (downloaded or installed files)

    Reports the first found for each group.
    """
    user_leb_dir = data_dir / "leb"
    tiers: dict[str, Any] = {}
    for tier in ["base", "medium", "extended"]:
        files: dict[str, dict[str, Any]] = {}
        for group in _LEB2_GROUPS:
            fname = f"{tier}_{group}.leb2"
            repo_path = _LEB2_DIR / fname
            user_path = user_leb_dir / fname
            if repo_path.exists():
                files[group] = _file_entry(repo_path)
            elif user_path.exists():
                files[group] = _file_entry(user_path)
            else:
                # Report repo path as the expected location
                files[group] = _file_entry(repo_path)
        present = sum(1 for info in files.values() if info["exists"])
        tiers[tier] = {
            "status": _status_from_counts(present, len(files)),
            "present": present,
            "total": len(files),
            "files": files,
        }
    return tiers


def _collect_planet_centers_artifacts(data_dir: Path) -> dict[str, Any]:
    """Collect status for generated planet-center output files."""
    tiers: dict[str, Any] = {}
    for tier in ["base", "medium", "extended"]:
        entry = _file_entry(_resolve_planet_centers_path(tier, data_dir))
        tiers[tier] = {"status": "ok" if entry["exists"] else "missing", "file": entry}
    return tiers


def _collect_missing_commands(report: dict[str, Any]) -> list[str]:
    """Suggest the next developer commands needed to fill current gaps."""
    commands: list[str] = []

    for tier, info in report["download_prerequisites"]["spk"]["tiers"].items():
        if info["status"] != "ok":
            commands.append(f"leph download spk-{tier}")

    if report["download_prerequisites"]["iers"]["status"] != "ok":
        commands.append("leph download iers")

    if report["download_prerequisites"]["planet_center_sources"]["status"] != "ok":
        commands.append("leph download planet-centers-sources")

    if report["download_prerequisites"]["assist"]["status"] != "ok":
        commands.append("leph download assist")

    for tier, info in report["generated_artifacts"]["leb1"].items():
        if not info["final"]["exists"]:
            commands.append(f"leph leb generate {tier} groups")

    for tier, info in report["generated_artifacts"]["leb2"].items():
        if (
            info["status"] != "ok"
            and report["generated_artifacts"]["leb1"][tier]["final"]["exists"]
        ):
            commands.append(f"leph leb2 convert {tier}")

    for tier, info in report["generated_artifacts"]["planet_centers"].items():
        if info["status"] != "ok":
            commands.append(f"leph generate planet-centers-{tier}")

    deduped: list[str] = []
    for command in commands:
        if command not in deduped:
            deduped.append(command)
    return deduped


def _build_dev_status_report() -> dict[str, Any]:
    """Build the development-oriented status report consumed by `leph status`."""
    from ..download import get_data_dir

    data_dir = get_data_dir()
    report = {
        "version": str(__version__),
        "project_root": str(_PROJECT_ROOT),
        "data_directory": str(data_dir),
        "download_prerequisites": {
            "spk": _collect_spk_status(data_dir),
            "iers": _collect_iers_status(),
            "planet_center_sources": _collect_planet_center_source_status(),
            "assist": _collect_assist_status(),
        },
        "generated_artifacts": {
            "leb1": _collect_leb1_status(),
            "leb2": _collect_leb2_status(data_dir),
            "planet_centers": _collect_planet_centers_artifacts(data_dir),
        },
    }
    report["missing_commands"] = _collect_missing_commands(report)
    return report


def _status_tag(status: str) -> str:
    """Return a colorized status tag for text rendering."""
    if status in {"ok", "ready", "expected_fallback"}:
        return click.style("[OK]", fg="green")
    if status == "partial":
        return click.style("[~~]", fg="yellow")
    return click.style("[--]", fg="yellow")


def _render_spk_status(report: dict[str, Any], verbose: int) -> None:
    """Render DE/SPK development status."""
    from ..download import _format_size

    click.echo(click.style("Download Prerequisites", bold=True))
    spk_report = report["download_prerequisites"]["spk"]
    for tier in ["base", "medium", "extended"]:
        info = spk_report["tiers"][tier]
        ready = info["ready_bodies"]
        total = info["required_bodies"]
        partial = info["partial_bodies"]
        expected_fallback = info["expected_fallback_bodies"]
        de_name = info["de_kernel"]["name"]
        suffix = f", {partial} partial" if partial else ""
        if expected_fallback:
            suffix += f", {expected_fallback} expected fallback"
        click.echo(
            f"  {_status_tag(info['status'])} spk-{tier}: {de_name} + {ready}/{total} downloadable bodies{suffix}"
        )
        if verbose >= 1:
            click.echo(f"    Range: {info['range']['start']} -> {info['range']['end']}")
            click.echo(f"    Kernel: {info['de_kernel']['path']}")
            click.echo(f"    Cache: {info['cache_dir']}")
            missing = [b["name"] for b in info["bodies"] if b["status"] == "missing"]
            partial_names = [
                b["name"] for b in info["bodies"] if b["status"] == "partial"
            ]
            fallback_names = [
                f"{b['name']} ({b['reason']})"
                for b in info["bodies"]
                if b["status"] == "expected_fallback"
            ]
            if partial_names:
                click.echo(f"    Partial bodies: {', '.join(partial_names)}")
            if missing:
                click.echo(f"    Missing bodies: {', '.join(missing)}")
            if fallback_names:
                label = "body" if len(fallback_names) == 1 else "bodies"
                click.echo(
                    f"    Expected fallback {label}: {', '.join(fallback_names)}"
                )
        if verbose >= 2:
            for body in info["bodies"]:
                label = body["name"]
                if body["status"] == "expected_fallback":
                    label = f"{label} (fallback)"
                click.echo(f"    {_status_tag(body['status'])} {label}")
                if body["path"]:
                    click.echo(f"      Path: {body['path']}")
                if body["coverage"]:
                    start_jd, end_jd = body["coverage"]
                    click.echo(f"      Coverage: JD {start_jd:.1f} -> {end_jd:.1f}")
                if body.get("reason"):
                    click.echo(f"      Reason: {body['reason']}")

    iers = report["download_prerequisites"]["iers"]
    click.echo(
        f"  {_status_tag(iers['status'])} iers: {iers['present']}/{iers['total']} files"
    )
    if verbose >= 1:
        click.echo(f"    Cache: {iers['cache_dir']}")
        for name, info in iers["files"].items():
            if info["exists"]:
                age = info.get("age_days")
                age_text = f" ({age:.0f}d old)" if age is not None else ""
                click.echo(f"    {_status_tag('ok')} {name}{age_text}")
            else:
                click.echo(f"    {_status_tag('missing')} {name}")
            if verbose >= 2:
                click.echo(f"      Path: {info['path']}")

    source_report = report["download_prerequisites"]["planet_center_sources"]
    source_summary = ", ".join(
        f"{tier} {info['present']}/{info['total']}"
        for tier, info in source_report["tiers"].items()
    )
    click.echo(
        f"  {_status_tag(source_report['status'])} planet-centers-sources: "
        f"shared {source_report['shared']['present']}/{source_report['shared']['total']}, {source_summary}"
    )
    if verbose >= 1:
        click.echo(f"    Cache: {source_report['cache_dir']}")
        for tier, info in source_report["tiers"].items():
            missing = [
                name for name, entry in info["files"].items() if not entry["exists"]
            ]
            if missing:
                click.echo(f"    {tier} missing: {', '.join(missing)}")
        if not source_report["shared"]["files"]["naif0012.tls"]["exists"]:
            click.echo("    Shared missing: naif0012.tls")
    if verbose >= 2:
        for tier, info in source_report["tiers"].items():
            click.echo(f"    {tier} source files:")
            for name, entry in info["files"].items():
                click.echo(
                    f"      {_status_tag('ok' if entry['exists'] else 'missing')} {name}"
                )
                click.echo(f"        Path: {entry['path']}")
        shared = source_report["shared"]["files"]["naif0012.tls"]
        click.echo(
            f"    shared: {_status_tag('ok' if shared['exists'] else 'missing')} naif0012.tls"
        )
        click.echo(f"      Path: {shared['path']}")

    assist = report["download_prerequisites"]["assist"]
    click.echo(
        f"  {_status_tag(assist['status'])} assist: {assist['present']}/{assist['total']} files"
    )
    if verbose >= 1:
        click.echo(f"    Directory: {assist['directory']}")
        if verbose >= 2:
            for name, entry in assist["files"].items():
                label = "ok" if entry["exists"] else "missing"
                size_text = (
                    f" ({_format_size(entry['size'])})"
                    if entry["size"] is not None
                    else ""
                )
                click.echo(f"    {_status_tag(label)} {name}{size_text}")
                click.echo(f"      Path: {entry['path']}")
    click.echo()


def _render_generated_artifacts(report: dict[str, Any], verbose: int) -> None:
    """Render generated development artifacts."""
    from ..download import _format_size

    click.echo(click.style("Generated Artifacts", bold=True))
    generated = report["generated_artifacts"]

    for tier in ["base", "medium", "extended"]:
        info = generated["leb1"][tier]
        suffix = ""
        if not info["final"]["exists"] and info["partial_present"]:
            suffix = f" (partials {info['partial_present']}/{info['partial_total']})"
        click.echo(
            f"  {_status_tag(info['status'])} LEB1 {tier}: {info['final']['name']}{suffix}"
        )
        if verbose >= 1:
            click.echo(f"    Path: {info['final']['path']}")
        if verbose >= 2:
            for name, entry in info["partials"].items():
                click.echo(
                    f"    {_status_tag('ok' if entry['exists'] else 'missing')} {name}"
                )
                click.echo(f"      Path: {entry['path']}")

    for tier in ["base", "medium", "extended"]:
        info = generated["leb2"][tier]
        click.echo(
            f"  {_status_tag(info['status'])} LEB2 {tier}: {info['present']}/{info['total']} groups"
        )
        if verbose >= 2:
            for group, entry in info["files"].items():
                size_text = (
                    f" ({_format_size(entry['size'])})"
                    if entry["size"] is not None
                    else ""
                )
                click.echo(
                    f"    {_status_tag('ok' if entry['exists'] else 'missing')} {tier}_{group}.leb2{size_text}"
                )
                click.echo(f"      Path: {entry['path']}")

    for tier in ["base", "medium", "extended"]:
        info = generated["planet_centers"][tier]
        click.echo(
            f"  {_status_tag(info['status'])} Planet centers {tier}: {info['file']['name']}"
        )
        if verbose >= 1:
            click.echo(f"    Path: {info['file']['path']}")
    click.echo()


def _render_missing_commands(report: dict[str, Any]) -> None:
    """Render follow-up commands when development gaps are present."""
    commands = report["missing_commands"]
    if not commands:
        click.echo(click.style("Development Gaps", bold=True))
        click.echo(
            f"  {_status_tag('ok')} no missing development prerequisites detected"
        )
        return

    click.echo(click.style("Development Gaps", bold=True))
    for command in commands:
        click.echo(f"  {command}")


def _print_dev_status(report: dict[str, Any], verbose: int) -> None:
    """Render the development status report as formatted text."""
    click.echo(click.style(f"libephemeris dev {report['version']}", bold=True))
    click.echo()
    click.echo(click.style("Workspace", bold=True))
    click.echo(f"  Project root:     {report['project_root']}")
    click.echo(f"  Data directory:   {report['data_directory']}")
    click.echo()
    _render_spk_status(report, verbose)
    _render_generated_artifacts(report, verbose)
    _render_missing_commands(report)


@click.command(
    "status",
    short_help="Show development prerequisites and generated-artifact status.",
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
    help="Increase detail: -v shows paths and gaps, -vv shows full file lists.",
)
def status(as_json: bool, verbose: int) -> None:
    """Show a development-oriented status recap for local libephemeris work."""
    report = _build_dev_status_report()
    if as_json:
        click.echo(json.dumps(report, indent=2))
        return
    _print_dev_status(report, verbose)
