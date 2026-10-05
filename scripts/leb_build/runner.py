"""Sequential execution, recovery and isolated scientific workers.

Provenance:
    Project-authored process orchestration. Scientific calls delegate unchanged
    to the registered LEB generators and verifiers; this module defines no model.
"""

from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
from typing import Any, Callable, Iterator

from .plan import (
    BuildConfig,
    Job,
    PROJECT_ROOT,
    SCHEMA_VERSION,
    scientific_environment,
    tier_range,
)
from .sources import attest_inputs, configure_runtime, input_stamps, require_disk_space
from .storage import (
    atomic_write,
    fsync_directory,
    invalidate,
    managed_path,
    save_manifest,
    sha256_file,
    staging_path,
)
from .validation import inspect_artifact


# ---------------------------------------------------------------------------
# Manifest creation and conservative recovery
# ---------------------------------------------------------------------------


def utc_now() -> str:
    """Return an unambiguous UTC timestamp for attempts and checkpoints."""
    return datetime.now(timezone.utc).isoformat()


def new_manifest(
    config: BuildConfig,
    jobs: list[Job],
    attestation: dict[str, Any],
    spks: dict[str, dict[str, str]],
) -> dict[str, Any]:
    """Create pending state; file presence alone never grants completion."""
    return {
        "schema": SCHEMA_VERSION,
        "created": utc_now(),
        "status": "pending",
        "config": config.to_dict(),
        "inputs": attestation,
        "spks": spks,
        "jobs": {
            job.id: {"definition": job.to_dict(), "status": "pending", "attempt": 0}
            for job in jobs
        },
    }


def recover_checkpoints(
    root: Path, manifest: dict[str, Any], jobs: list[Job]
) -> set[str]:
    """Rehash/inspect recorded outputs and invalidate interrupted descendants.

    Generation and verification are separate checkpoints. A failed verification
    therefore leaves valid generated coefficients reusable on the next resume.
    Unregistered final files and staging files are never adopted.
    """
    invalid = set()
    inspected = {}
    for job in jobs:
        record = manifest["jobs"][job.id]
        if record["status"] != "done":
            invalid.add(job.id)
            continue
        try:
            artifact = record["artifact"]
            if job.output not in inspected:
                inspected[job.output] = inspect_artifact(
                    managed_path(root, job.output), job
                )
            if inspected[job.output] != artifact:
                invalid.add(job.id)
        except (OSError, ValueError, KeyError):
            invalid.add(job.id)
    return invalidate(manifest, jobs, invalid)


# ---------------------------------------------------------------------------
# Signal and subprocess lifetime management
# ---------------------------------------------------------------------------


class BuildInterrupted(Exception):
    """A handled operator signal, retaining conventional shell exit semantics."""

    def __init__(self, signum: int) -> None:
        super().__init__(f"Interrupted by {signal.Signals(signum).name}")
        self.signum = signum


@contextmanager
def interruption_signals() -> Iterator[None]:
    """Handle the first signal and allow shutdown to finish despite repeats."""
    interrupted = False

    def handle(signum: int, _frame: Any) -> None:
        nonlocal interrupted
        if not interrupted:
            interrupted = True
            raise BuildInterrupted(signum)

    previous = {
        sig: signal.signal(sig, handle) for sig in (signal.SIGINT, signal.SIGTERM)
    }
    try:
        yield
    finally:
        for sig, handler in previous.items():
            signal.signal(sig, handler)


def stop_process_group(process: subprocess.Popen, signum: int) -> None:
    """Stop the active worker group, escalating after a bounded grace period."""
    try:
        os.killpg(process.pid, signum)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()


def run_process(
    command: list[str], env: dict[str, str], log: Path, lock_fd: int
) -> int:
    """Run one isolated process, retain the lock and leave an append-only log."""
    with log.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps({"command": command, "started": utc_now()}) + "\n")
        stream.flush()
        process = subprocess.Popen(
            command,
            cwd=PROJECT_ROOT,
            env=env,
            stdout=stream,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            pass_fds=(lock_fd,),
        )
        try:
            while True:
                try:
                    return process.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    print(f"  still running; log: {log}", flush=True)
        except BaseException as exc:
            signum = exc.signum if isinstance(exc, BuildInterrupted) else signal.SIGTERM
            stop_process_group(process, signum)
            raise


# ---------------------------------------------------------------------------
# Build runner: one job, one staging output, one durable checkpoint
# ---------------------------------------------------------------------------


Executor = Callable[[Job, Path, Path], int]


class BuildRunner:
    """Run an authenticated job plan while keeping recovery state explicit.

    An optional executor supports deterministic tests without launching a full
    scientific build. Production always uses a fresh worker with an inherited
    build lock. The manifest remains authoritative until final completion.
    """

    def __init__(
        self,
        root: Path,
        config: BuildConfig,
        jobs: list[Job],
        manifest: dict[str, Any],
        lock_fd: int,
        executor: Executor | None = None,
    ) -> None:
        self.root = root
        self.config = config
        self.jobs = jobs
        self.manifest = manifest
        self.lock_fd = lock_fd
        self.executor = executor or self._execute_worker
        self.active_job: Job | None = None
        self.stamps = input_stamps(config)
        saved_stamps = manifest["inputs"].get("stamps")
        if saved_stamps is not None and saved_stamps != {
            path: list(stamp) for path, stamp in self.stamps.items()
        }:
            raise ValueError("Build inputs changed after preflight")

    def _execute_worker(self, job: Job, target: Path, log: Path) -> int:
        """Launch the internal worker entry point with explicit file arguments."""
        command = [
            sys.executable,
            "-B",
            str(PROJECT_ROOT / "scripts/regenerate_leb.py"),
            "--output-dir",
            str(self.root),
            "--worker-job",
            job.id,
            "--worker-output",
            str(target),
        ]
        self.manifest["jobs"][job.id]["command"] = command
        save_manifest(self.root, self.manifest)
        env = scientific_environment(self.config, job.tier)
        env["LEB_BUILD_LOCK_FD"] = str(self.lock_fd)
        return run_process(command, env, log, self.lock_fd)

    def _check_inputs_unchanged(self) -> None:
        """Detect source edits, replacement, deletion and new candidates promptly."""
        if input_stamps(self.config) != self.stamps:
            raise ValueError("Build inputs changed; use a new output directory")

    def _job_space_budget(self, job: Job) -> int:
        """Include staging without assuming a favorable compression ratio."""
        from .sources import estimated_auxiliary_bytes, estimated_body_bytes

        if not job.writes_output:
            return 0
        if job.kind in ("merge", "convert"):
            size = sum(
                managed_path(self.root, path).stat().st_size for path in job.inputs
            )
            if job.aux_source:
                size += managed_path(self.root, job.aux_source).stat().st_size
            return size * (2 if job.kind == "convert" else 1) + 64 * 1024**2
        body = job.bodies[0]
        size = estimated_body_bytes(job.tier, body)
        if body == 0:
            size += estimated_auxiliary_bytes(job.tier)
        return size + 64 * 1024**2

    def _start_job(self, job: Job) -> tuple[Path, Path]:
        """Record the attempt before launch; discard only its known stale staging."""
        self._check_inputs_unchanged()
        require_disk_space(self.root, self._job_space_budget(job))
        record = self.manifest["jobs"][job.id]
        if job.writes_output and record["attempt"]:
            staging_path(self.root, job, record["attempt"]).unlink(missing_ok=True)
        record.update(
            status="running", attempt=record["attempt"] + 1, started=utc_now()
        )
        for field in (
            "error",
            "returncode",
            "artifact",
            "finished",
            "command",
            "signal",
        ):
            record.pop(field, None)
        target = (
            staging_path(self.root, job, record["attempt"])
            if job.writes_output
            else managed_path(self.root, job.output)
        )
        target.parent.mkdir(parents=True, exist_ok=True)
        log = managed_path(self.root, f"logs/{job.id}-{record['attempt']}.log")
        log.parent.mkdir(parents=True, exist_ok=True)
        record["log"] = str(log.relative_to(self.root))
        save_manifest(self.root, self.manifest)
        return target, log

    def _finish_job(self, job: Job, target: Path, returncode: int) -> None:
        """Promote structurally accepted bytes and checkpoint scientific PASS."""
        record = self.manifest["jobs"][job.id]
        record["returncode"] = returncode
        if returncode:
            raise ValueError(
                f"Job failed (exit {returncode}): {job.id}; see {record['log']}"
            )
        self._check_inputs_unchanged()
        artifact = inspect_artifact(target, job)
        if job.writes_output:
            with target.open("rb") as stream:
                os.fsync(stream.fileno())
            destination = managed_path(self.root, job.output)
            os.replace(target, destination)
            fsync_directory(destination.parent)
        record.update(status="done", artifact=artifact, finished=utc_now())
        save_manifest(self.root, self.manifest)

    def _record_failure(self, exc: Exception) -> None:
        """Retain reusable completed phases and identify the failed attempt."""
        interrupted = isinstance(exc, BuildInterrupted)
        self.manifest["status"] = "interrupted" if interrupted else "failed"
        if self.active_job is not None:
            job = self.active_job
            record = self.manifest["jobs"][job.id]
            record.update(
                status=self.manifest["status"], error=str(exc), finished=utc_now()
            )
            if isinstance(exc, BuildInterrupted):
                record["signal"] = exc.signum
            if job.writes_output:
                staging_path(self.root, job, record["attempt"]).unlink(missing_ok=True)
        save_manifest(self.root, self.manifest)

    def _complete_build(self) -> None:
        """Rehash all final exports and source contents before declaring complete."""
        self._check_inputs_unchanged()
        if attest_inputs(self.config) != self.manifest["inputs"]:
            raise ValueError("Input fingerprint changed before finalization")
        checksums = []
        for job in self.jobs:
            if job.export:
                expected = self.manifest["jobs"][job.id]["artifact"]["sha256"]
                if sha256_file(managed_path(self.root, job.output)) != expected:
                    raise ValueError(
                        f"Export changed before finalization: {job.output}"
                    )
                checksums.append(f"{expected}  {job.output}\n")
        atomic_write(managed_path(self.root, "checksums.sha256"), "".join(checksums))
        self.manifest.update(status="complete", finished=utc_now())
        save_manifest(self.root, self.manifest)

    def run(self) -> None:
        """Execute pending phases in dependency order; stop on the first failure."""
        self.manifest["status"] = "running"
        self.manifest.pop("finished", None)
        # Previous final checksums no longer attest a build being rerun.
        managed_path(self.root, "checksums.sha256").unlink(missing_ok=True)
        save_manifest(self.root, self.manifest)
        try:
            with interruption_signals():
                for index, job in enumerate(self.jobs, 1):
                    record = self.manifest["jobs"][job.id]
                    if record["status"] == "done":
                        continue
                    if any(
                        self.manifest["jobs"][dep]["status"] != "done"
                        for dep in job.dependencies
                    ):
                        raise ValueError(f"Unsatisfied dependencies for {job.id}")
                    self.active_job = job
                    print(f"[{index}/{len(self.jobs)}] {job.id}", flush=True)
                    target, log = self._start_job(job)
                    self._finish_job(job, target, self.executor(job, target, log))
                    self.active_job = None
                self._complete_build()
        except Exception as exc:
            self._record_failure(exc)
            raise


# ---------------------------------------------------------------------------
# Internal scientific worker (never an alternate public build command)
# ---------------------------------------------------------------------------


def execute_scientific_job(
    root: Path, config: BuildConfig, job: Job, target: Path, spks: dict[str, str]
) -> None:
    """Delegate one phase to existing generators with fixed source settings."""
    from scripts import generate_leb, generate_leb2

    needed_spks = (
        {body: path for body, path in spks.items() if int(body) in job.bodies}
        if job.kind in ("generate", "verify1")
        else {}
    )
    configure_runtime(config, job.tier, needed_spks)
    start, end = tier_range(job.tier)
    if job.kind == "generate":
        generate_leb.assemble_leb(
            output=str(target),
            jd_start=start,
            jd_end=end,
            bodies=list(job.bodies),
            workers=1,
            verbose=True,
            skip_aux=job.bodies != (0,),
            tier=job.tier,
        )
    elif job.kind == "merge":
        generate_leb.merge_leb_files(
            [str(managed_path(root, path)) for path in job.inputs],
            str(target),
            aux_source=str(managed_path(root, job.aux_source))
            if job.aux_source
            else None,
        )
    elif job.kind == "convert":
        generate_leb2.convert_leb1_to_leb2(
            str(managed_path(root, job.inputs[0])),
            str(target),
            group=job.group,
            expected_tier=job.tier,
        )
    elif job.kind == "verify1":
        if not generate_leb.verify_leb(str(target), n_samples=config.verify_samples):
            raise ValueError(f"Scientific verification failed: {job.id}")
    elif job.kind == "verify2":
        if not generate_leb2.verify_leb2(
            str(target),
            reference_leb1=str(managed_path(root, job.inputs[0])),
            n_samples=config.leb2_verify_samples,
            expected_group=job.group,
            expected_tier=job.tier,
        ):
            raise ValueError(f"Scientific verification failed: {job.id}")
    else:
        raise ValueError(f"Unknown job kind: {job.kind}")
