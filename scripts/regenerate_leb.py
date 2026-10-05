#!/usr/bin/env python3
"""Regenerate every canonical LEB artifact in an isolated, resumable directory.

Usage:
    uv run python -B scripts/regenerate_leb.py --output-dir /new/build --dry-run
    uv run python -B scripts/regenerate_leb.py --output-dir /new/build --doctor
    uv run python -B scripts/regenerate_leb.py --output-dir /new/build
    uv run python -B scripts/regenerate_leb.py --output-dir /new/build --resume
    uv run python -B scripts/regenerate_leb.py --output-dir /new/build --status

The default plan exports 15 LEB1 (.leb) files and 12 chunked LEB2 (.leb2) files.
Checkpoints are per body and per phase, not inside a body's fitting/integration.
Generation requires provisioned local JPL/SPK/ASSIST inputs and a passing
provenance gate. It never installs, publishes or refreshes existing data.

Module map:
    leb_build.plan        Immutable settings, canonical inventories and jobs.
    leb_build.storage     Output boundaries, durable manifests and POSIX locks.
    leb_build.sources     Offline preflight, resource budgets and input hashes.
    leb_build.validation  Structural artifact acceptance.
    leb_build.runner      Execution, interruption, recovery and scientific calls.

Provenance:
    Project-authored build orchestration over the registered generators. No
    astronomical model or compatibility-reference output is defined here.
"""

from __future__ import annotations

import argparse
from collections import Counter
import os
from pathlib import Path
import sys

# Support both direct execution and imports from the repository test suite.
# Keep orchestration imports from creating provenance-blocking bytecode caches.
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.leb_build.plan import BuildConfig, build_jobs
from scripts.leb_build.runner import (
    BuildInterrupted,
    BuildRunner,
    execute_scientific_job,
    new_manifest,
    recover_checkpoints,
    interruption_signals,
)
from scripts.leb_build.sources import (
    estimated_build_bytes,
    preflight,
    require_disk_space,
)
from scripts.leb_build.storage import (
    build_lock,
    invalidate,
    load_manifest,
    managed_path,
    safe_root,
    save_manifest,
    staging_path,
    validate_job_records,
)


# ---------------------------------------------------------------------------
# Command-line interface and immutable resume settings
# ---------------------------------------------------------------------------


def positive_int(value: str) -> int:
    """Parse positive sample counts with a useful argparse diagnostic."""
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def build_parser() -> argparse.ArgumentParser:
    """Expose operator actions while keeping scientific workers internal."""
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n", 1)[0],
        epilog="POSIX only. Ctrl+C stops the worker; --resume retries the incomplete phase.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="New directory outside the checkout and source caches",
    )
    actions = parser.add_mutually_exclusive_group()
    actions.add_argument(
        "--doctor",
        action="store_true",
        help="Check provenance, dependencies, local inputs and space; no output writes",
    )
    actions.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the entire job plan without checking source payloads or writing files",
    )
    actions.add_argument(
        "--status",
        action="store_true",
        help="Read the manifest summary; no mutation or source access",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Use saved settings and revalidate recorded checkpoints",
    )
    parser.add_argument(
        "--rerun",
        metavar="JOB_ID",
        help="With --resume, invalidate this phase and all dependent phases",
    )
    parser.add_argument(
        "--tiers", help="all (default), or comma-separated base,medium,extended"
    )
    parser.add_argument(
        "--verify-samples",
        type=positive_int,
        default=None,
        help="LEB1 samples per body (default: 500)",
    )
    parser.add_argument(
        "--leb2-verify-samples",
        type=positive_int,
        default=None,
        help="LEB2 samples per body (default: 200)",
    )
    parser.add_argument(
        "--data-dir", type=Path, help="Provisioned planetary/auxiliary source directory"
    )
    parser.add_argument("--spk-dir", type=Path, help="Provisioned minor-body SPK cache")
    parser.add_argument("--worker-job", help=argparse.SUPPRESS)
    parser.add_argument("--worker-output", type=Path, help=argparse.SUPPRESS)
    return parser


def resolve_config(args: argparse.Namespace, saved: BuildConfig | None) -> BuildConfig:
    """Inherit saved options and reject changes rather than mix different builds."""
    from libephemeris.state import _resolve_data_dir, get_spk_cache_dir
    from libephemeris.spk_auto import DEFAULT_AUTO_SPK_DIR
    from libephemeris.rebound_integration import _ASSIST_DEFAULT_DIR
    from scripts.generate_leb import TIER_CONFIGS

    options = (
        saved.to_dict()
        if saved
        else {
            "tiers": list(TIER_CONFIGS),
            "verify_samples": 500,
            "leb2_verify_samples": 200,
            "data_dir": str(Path(_resolve_data_dir()).resolve()),
            "spk_dir": str(Path(get_spk_cache_dir() or DEFAULT_AUTO_SPK_DIR).resolve()),
            "assist_dir": str(
                Path(os.environ.get("ASSIST_DIR", _ASSIST_DEFAULT_DIR)).resolve()
            ),
        }
    )
    if args.tiers is not None:
        tiers = list(TIER_CONFIGS) if args.tiers == "all" else args.tiers.split(",")
        if len(set(tiers)) != len(tiers) or any(
            tier not in TIER_CONFIGS for tier in tiers
        ):
            raise ValueError(
                "--tiers expects all or unique comma-separated canonical tiers"
            )
        options["tiers"] = [tier for tier in TIER_CONFIGS if tier in tiers]
    for name in ("verify_samples", "leb2_verify_samples", "data_dir", "spk_dir"):
        value = getattr(args, name)
        if value is not None:
            options[name] = (
                str(value.expanduser().resolve()) if isinstance(value, Path) else value
            )
    config = BuildConfig.from_dict(options)
    if saved is not None and config != saved:
        raise ValueError(
            "Resume options differ from the manifest; use a new output directory"
        )
    return config


# ---------------------------------------------------------------------------
# Read-only actions and authenticated worker dispatch
# ---------------------------------------------------------------------------


def show_status(manifest: dict) -> None:
    """Describe recorded progress without pretending to revalidate file contents."""
    counts = Counter(record["status"] for record in manifest["jobs"].values())
    print(f"Recorded build status: {manifest['status']}")
    print(
        "Phases: "
        + ", ".join(f"{state}={count}" for state, count in sorted(counts.items()))
    )
    for job_id, record in manifest["jobs"].items():
        if record["status"] in ("failed", "interrupted", "running"):
            print(
                f"  {job_id}: {record['status']}; log: {record.get('log', '(not started)')}"
            )


def run_worker(args: argparse.Namespace) -> None:
    """Accept only the active manifest job and the inherited build-lock inode."""
    root = args.output_dir.resolve()
    manifest = load_manifest(root)
    config = BuildConfig.from_dict(manifest["config"])
    root = safe_root(args.output_dir, config)
    jobs = build_jobs(config)
    validate_job_records(manifest, jobs)
    job = next((job for job in jobs if job.id == args.worker_job), None)
    if job is None or manifest["jobs"][job.id]["status"] != "running":
        raise ValueError("Worker requires a recorded active job")
    lock_fd = int(os.environ["LEB_BUILD_LOCK_FD"])
    inherited = os.fstat(lock_fd)
    expected_lock = managed_path(root, ".lock").stat()
    if (inherited.st_dev, inherited.st_ino) != (
        expected_lock.st_dev,
        expected_lock.st_ino,
    ):
        raise ValueError("Worker did not inherit this build's lock")
    record = manifest["jobs"][job.id]
    target = (
        staging_path(root, job, record["attempt"])
        if job.writes_output
        else managed_path(root, job.output)
    )
    if args.worker_output is None or args.worker_output.absolute() != target:
        raise ValueError("Unexpected worker output path")
    with interruption_signals():
        execute_scientific_job(
            root, config, job, target, manifest["spks"].get(job.tier, {})
        )


# ---------------------------------------------------------------------------
# Public entry point: preflight first, writes only under the exclusive lock
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    """Return shell-compatible status; leave interrupted builds resumable."""
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.rerun and (not args.resume or args.doctor or args.dry_run or args.status):
        parser.error("--rerun requires --resume and an execution action")
    try:
        if args.worker_job:
            run_worker(args)
            return 0
        if args.worker_output:
            raise ValueError("--worker-output is internal to an active worker")
        if args.status:
            recorded_status = load_manifest(args.output_dir.expanduser().resolve())
            show_status(recorded_status)
            return 0
        manifest = (
            load_manifest(args.output_dir.expanduser().resolve())
            if args.resume
            else None
        )
        saved = BuildConfig.from_dict(manifest["config"]) if manifest else None
        config = resolve_config(args, saved)
        root = safe_root(args.output_dir, config)
        jobs = build_jobs(config)
        if manifest:
            validate_job_records(manifest, jobs)
        elif root.exists() and any(root.iterdir()):
            raise ValueError(
                "Output directory is not empty; use --resume with its manifest"
            )
        if args.dry_run:
            print(
                f"{len(jobs)} phases; {sum(job.export for job in jobs)} canonical exports"
            )
            for job in jobs:
                print(f"{job.id:<36} {job.output}")
            return 0
        print("Checking provenance, local sources and input hashes...", flush=True)
        attestation, spks = preflight(config)
        if manifest and (attestation != manifest["inputs"] or spks != manifest["spks"]):
            raise ValueError(
                "Code, environment or source fingerprint changed; use a new directory"
            )
        budget = estimated_build_bytes(config)
        if manifest:
            budget -= sum(
                managed_path(root, job.output).stat().st_size
                for job in jobs
                if job.writes_output
                and manifest["jobs"][job.id]["status"] == "done"
                and managed_path(root, job.output).is_file()
            )
        require_disk_space(root, max(0, budget))
        if args.doctor:
            print(
                f"Preflight passed; {sum(job.export for job in jobs)} exports, conservative total budget {estimated_build_bytes(config) / 1024**3:.1f} GiB"
            )
            print(
                "Fitting still retains one body's arrays in RAM; extended jobs may be large."
            )
            return 0
        root.mkdir(parents=True, exist_ok=True)
        with build_lock(root) as lock_fd:
            # Reread under the lock: another invocation may have completed while
            # this invocation was performing its read-only preflight.
            if args.resume:
                current = load_manifest(root)
                if (
                    current["config"] != config.to_dict()
                    or current["inputs"] != attestation
                    or current["spks"] != spks
                ):
                    raise ValueError("Build configuration changed during preflight")
                manifest = current
                validate_job_records(manifest, jobs)
                recover_checkpoints(root, manifest, jobs)
                if args.rerun:
                    invalidate(manifest, jobs, {args.rerun})
            else:
                if set(path.name for path in root.iterdir()) - {".lock"}:
                    raise ValueError(
                        "Output changed during preflight; refusing to overwrite"
                    )
                manifest = new_manifest(config, jobs, attestation, spks)
            save_manifest(root, manifest)
            runner = BuildRunner(root, config, jobs, manifest, lock_fd)
            runner.run()
        print(f"Complete: {sum(job.export for job in jobs)} verified exports in {root}")
        return 0
    except BuildInterrupted as exc:
        print(
            f"{exc}. Resume using the same --output-dir and --resume.", file=sys.stderr
        )
        return 128 + exc.signum
    except KeyboardInterrupt:
        print("Interrupted; use --resume if a manifest was created.", file=sys.stderr)
        return 130
    except (ValueError, OSError, KeyError, ImportError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
