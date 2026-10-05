"""Offline preflight and frozen-input attestation.

Provenance:
    Project orchestration over local registered NASA/JPL and IAU inputs. Coverage
    checks delegate to existing generator helpers; no new orbit model is defined.
"""

from __future__ import annotations

import importlib
from importlib.metadata import version
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

from .plan import BuildConfig, PROJECT_ROOT, tier_groups, tier_range
from .storage import file_stamp, sha256_file


# ---------------------------------------------------------------------------
# Scientific runtime: explicit public settings, local inputs only
# ---------------------------------------------------------------------------


def spk_target_id(path: str, registry_id: int) -> int:
    """Return the NAIF target stored in a minor-body SPK.

    Horizons kernels may number a body 20000000+N instead of the registry ID, so
    the file is the authority. A kernel holds many segments for one target, hence
    the unique set is inspected.
    """
    from libephemeris.spk import _get_spk_targets

    targets = set(_get_spk_targets(path))
    if registry_id in targets:
        return registry_id
    small_bodies = [target for target in targets if target > 1_000_000]
    if len(small_bodies) != 1:
        raise ValueError(f"Cannot determine the single NAIF target in SPK {path}")
    return small_bodies[0]


def configure_runtime(config: BuildConfig, tier: str, spks: dict[str, str]) -> None:
    """Set the worker's source policy before opening any scientific input.

    Each worker is a fresh Python process. The preflight uses the same settings
    to ensure it inspects the kernels that generation will actually open.
    """
    import libephemeris as ephem
    from libephemeris import _config_toml
    from libephemeris.constants import TIDAL_AUTOMATIC
    from scripts.generate_leb import TIER_CONFIGS, _ASTEROID_NAIF

    _config_toml.load_config(os.devnull)
    os.environ["LIBEPHEMERIS_DATA_DIR"] = config.data_dir
    os.environ["LIBEPHEMERIS_SPK_DIR"] = config.spk_dir
    os.environ["ASSIST_DIR"] = config.assist_dir
    ephem.close()
    ephem.set_calc_mode("skyfield")
    ephem.set_network_policy("sealed")
    ephem.set_precision_tier(tier)
    ephem.set_ephe_path(config.data_dir)
    ephem.set_jpl_file(str(Path(config.data_dir) / TIER_CONFIGS[tier][0]))
    ephem.set_spk_cache_dir(config.spk_dir)
    ephem.set_iers_auto_download(False)
    ephem.set_delta_t_userdef(None)
    ephem.set_tid_acc(TIDAL_AUTOMATIC)
    for body, path in spks.items():
        # Horizons kernels may number a target 20000000+N instead of the registry
        # ID; the file itself is the authority, as in the runtime auto-SPK flow.
        ephem.register_spk_body(
            int(body), path, spk_target_id(path, _ASTEROID_NAIF[int(body)])
        )


def check_provenance() -> None:
    """Run the repository integrity gate before inspecting source payloads."""
    result = subprocess.run(
        [sys.executable, "-B", str(PROJECT_ROOT / "scripts/check_provenance.py")],
        cwd=PROJECT_ROOT,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode:
        report = (result.stdout + result.stderr).strip()
        raise ValueError(
            f"Provenance gate failed; resolve it before building.\n{report}"
        )


def environment_versions(config: BuildConfig) -> dict[str, str]:
    """Import required extensions and record the actual installed versions."""
    distributions = {
        "numpy": "numpy",
        "skyfield": "skyfield",
        "skyfield-data": "skyfield_data",
        "pyerfa": "erfa",
        "zstandard": "zstandard",
        "jplephem": "jplephem",
    }
    if "extended" in config.tiers:
        distributions.update({"rebound": "rebound", "assist": "assist"})
    versions = {
        "python": sys.version,
        "executable": str(Path(sys.executable).resolve()),
    }
    for distribution, module in distributions.items():
        importlib.import_module(module)
        versions[distribution] = version(distribution)
    return versions


# ---------------------------------------------------------------------------
# Discover the finite set of permitted local inputs
# ---------------------------------------------------------------------------


def minor_spk_candidates(config: BuildConfig) -> dict[int, list[Path]]:
    """List registered body-name prefixes without discovering unrelated BSPs."""
    from libephemeris.constants import SPK_BODY_NAME_MAP
    from libephemeris.spk import _sanitize_filename
    from scripts.generate_leb import _ASTEROID_NAIF

    required_bodies = {
        body
        for tier in config.tiers
        for bodies in tier_groups(tier).values()
        for body in bodies
    }
    return {
        body: sorted(
            Path(config.spk_dir).glob(
                f"{_sanitize_filename(str(SPK_BODY_NAME_MAP[body][0]))}_*.bsp"
            )
        )
        for body in sorted(required_bodies & _ASTEROID_NAIF.keys())
    }


def input_paths(config: BuildConfig) -> list[Path]:
    """Enumerate code, models and all candidates that the generators can select.

    Including alternate SPKs detects additions/deletions as well as edits: a
    newly cached wider kernel must not silently change a resumed build. Existing
    LEB artifacts and compiled caches are never inputs to this workflow.
    """
    from scripts.generate_leb import TIER_CONFIGS
    from libephemeris.rebound_integration import _ASSIST_DEFAULT_DIR

    package = PROJECT_ROOT / "libephemeris"
    paths = set(package.rglob("*.py"))
    paths.update((PROJECT_ROOT / "scripts/leb_build").glob("*.py"))
    paths.update(
        PROJECT_ROOT / name
        for name in (
            "scripts/regenerate_leb.py",
            "scripts/generate_leb.py",
            "scripts/generate_leb2.py",
            "scripts/check_provenance.py",
            "pyproject.toml",
            "uv.lock",
        )
    )
    for path in (package / "data").rglob("*"):
        if path.is_file() and path.suffix not in (".leb", ".leb2", ".pyc"):
            paths.add(path)
    data = Path(config.data_dir)
    for tier in config.tiers:
        paths.add(data / TIER_CONFIGS[tier][0])
        centers = data / f"planet_centers_{tier}.bsp"
        if centers.exists():
            paths.add(centers)
    legacy_centers = data / "planet_centers.bsp"
    if "base" in config.tiers and legacy_centers.exists():
        paths.add(legacy_centers)
    for candidates in minor_spk_candidates(config).values():
        paths.update(candidates)
    # These are the only external auxiliary tables read by the default runtime.
    for name in ("finals2000A.data", "leap_seconds.dat", "deltat.data"):
        path = data / "iers_cache" / name
        if path.exists():
            paths.add(path)
    if "extended" in config.tiers:
        paths.add(_ASSIST_DEFAULT_DIR / "linux_m13000p17000.441")
        for directory in (
            Path(config.assist_dir),
            _ASSIST_DEFAULT_DIR,
            PROJECT_ROOT / "data",
        ):
            candidate = directory / "sb441-n16.bsp"
            if candidate.exists():
                paths.add(candidate)
    return sorted({path.resolve(strict=True) for path in paths})


def attest_inputs(config: BuildConfig) -> dict[str, Any]:
    """Hash all potential inputs with their process environment."""
    before = input_stamps(config)
    attestation = {
        "environment": environment_versions(config),
        "files": {path: sha256_file(Path(path)) for path in before},
        "stamps": {path: list(stamp) for path, stamp in before.items()},
    }
    if input_stamps(config) != before:
        raise ValueError("Build inputs changed during source attestation")
    return attestation


def input_stamps(config: BuildConfig) -> dict[str, tuple[int, ...]]:
    """Snapshot file identities for cheap checks before/after each job."""
    return {str(path): file_stamp(path) for path in input_paths(config)}


def select_minor_spks(config: BuildConfig) -> dict[str, dict[str, str]]:
    """Select useful local coverage for fitting and meaningful verification.

    SPK-backed bodies retain their real coverage; a useful partial interval is
    allowed. Extended N-body bodies require the wider Horizons verification
    window as well as their seed. Selections are per tier, so disjoint source
    windows cannot silently replace another tier's selected kernel.
    """
    from libephemeris.minor_bodies import HORIZONS_SPK_JD_MIN, HORIZONS_SPK_JD_MAX
    from libephemeris.exotic_bodies import (
        EXOTIC_EXTENDED_IDS,
        EXOTIC_ASSIST_PERTURBER_IDS,
    )
    from scripts.generate_leb import (
        GENERATION_SPK_PADDING_DAYS,
        _get_asteroid_spk_range,
    )

    candidates = minor_spk_candidates(config)
    selected: dict[str, dict[str, str]] = {}
    for tier in config.tiers:
        selected[tier] = {}
        start, end = tier_range(tier)
        required_start = max(start, HORIZONS_SPK_JD_MIN)
        required_end = min(end, HORIZONS_SPK_JD_MAX)
        bodies = (*tier_groups(tier)["asteroids"], *tier_groups(tier)["exotics"])
        nbody_ids = (
            set(EXOTIC_EXTENDED_IDS) - set(EXOTIC_ASSIST_PERTURBER_IDS)
            if tier == "extended"
            else set()
        )
        for body in bodies:
            covering = []
            for path in candidates[body]:
                coverage = _get_asteroid_spk_range(str(path), body)
                if coverage is None:
                    continue
                if body in nbody_ids:
                    usable = (
                        coverage[0] <= required_start + 30.0
                        and coverage[1] >= required_end - 30.0
                    )
                    span = coverage[1] - coverage[0]
                else:
                    usable_start = max(start, coverage[0] + GENERATION_SPK_PADDING_DAYS)
                    usable_end = min(end, coverage[1] - GENERATION_SPK_PADDING_DAYS)
                    span = usable_end - usable_start
                    # Same 20-year usefulness policy as assemble_leb; actual
                    # fitting grids/ranges are resolved there, never here.
                    usable = span >= 20 * 365.25
                if usable:
                    covering.append((span, str(path.resolve())))
            if not covering:
                raise ValueError(
                    f"Missing usable local SPK for {tier} body {body}. Provision with "
                    "scripts/download_max_range_spk.py --output-dir <spk-dir>. "
                    "The build never downloads or refreshes shared caches."
                )
            selected[tier][str(body)] = max(covering)[1]
    return selected


# ---------------------------------------------------------------------------
# Preflight: provenance, source authentication and resource budgets
# ---------------------------------------------------------------------------


def estimated_body_bytes(tier: str, body: int) -> int:
    """Estimate uncompressed coefficients over the full advertised tier range."""
    from libephemeris.leb_format import BODY_PARAMS, segment_byte_size

    start, end = tier_range(tier)
    interval, degree, _, components = BODY_PARAMS[body]
    return math.ceil((end - start) / interval) * segment_byte_size(degree, components)


def estimated_auxiliary_bytes(tier: str) -> int:
    """Estimate the existing nutation/Delta-T grids plus catalog overhead."""
    from libephemeris.leb_format import DELTA_T_ENTRY_SIZE, segment_byte_size
    from scripts.generate_leb import (
        NUTATION_INTERVAL,
        NUTATION_DEGREE,
        NUTATION_COMPONENTS,
        DELTA_T_INTERVAL,
    )

    start, end = tier_range(tier)
    nutation = math.ceil((end - start) / NUTATION_INTERVAL) * segment_byte_size(
        NUTATION_DEGREE, NUTATION_COMPONENTS
    )
    delta_t = (math.ceil((end - start) / DELTA_T_INTERVAL) + 1) * DELTA_T_ENTRY_SIZE
    catalog_and_headers_allowance = 16 * 1024**2
    return nutation + delta_t + catalog_and_headers_allowance


def estimated_build_bytes(config: BuildConfig) -> int:
    """Budget body checkpoints, group/tier exports, LEB2 and staging pessimistically."""
    total = 0
    for tier in config.tiers:
        raw = sum(
            estimated_body_bytes(tier, body)
            for bodies in tier_groups(tier).values()
            for body in bodies
        )
        # Assume no compression benefit and allow one large staging payload.
        total += raw * 5 + estimated_auxiliary_bytes(tier) * 7
    return total + 2 * 1024**3


def require_disk_space(root: Path, required: int) -> None:
    """Check the destination volume without creating any directory."""
    ancestor = root
    while not ancestor.exists():
        ancestor = ancestor.parent
    free = shutil.disk_usage(ancestor).free
    if free < required:
        raise ValueError(
            f"Insufficient disk space: need {required / 1024**3:.1f} GiB, "
            f"available {free / 1024**3:.1f} GiB"
        )


def preflight(config: BuildConfig) -> tuple[dict[str, Any], dict[str, dict[str, str]]]:
    """Fail before expensive work if integrity, imports or local sources are missing."""
    check_provenance()
    environment_versions(config)
    for directory in (config.data_dir, config.spk_dir):
        if not Path(directory).is_dir():
            raise ValueError(f"Missing provisioned source directory: {directory}")
    spks = select_minor_spks(config)
    for tier in config.tiers:
        configure_runtime(config, tier, {})
        from libephemeris.state import get_current_file_data, get_planets
        from scripts.generate_leb import TIER_CONFIGS, _get_spk_jd_range

        planets = get_planets()
        actual_path, _, _, _ = get_current_file_data(0)
        expected = Path(config.data_dir) / TIER_CONFIGS[tier][0]
        if Path(actual_path).resolve() != expected.resolve():
            raise ValueError(f"Unexpected kernel fallback for {tier}: {actual_path}")
        source_start, source_end = _get_spk_jd_range(planets)
        start, end = tier_range(tier)
        if source_start > start or source_end < end:
            raise ValueError(f"Kernel does not cover the canonical {tier} interval")
        if tier == "extended":
            from scripts.generate_leb import _nbody_coverage_for_range

            _nbody_coverage_for_range(start, end)
    return attest_inputs(config), spks
