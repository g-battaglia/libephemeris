"""Filesystem boundaries, atomic checkpoints and POSIX locking.

Provenance:
    Project-authored file/process infrastructure using Python's standard library.
    Hashes attest build inputs and outputs; they are not astronomical values.
"""

from __future__ import annotations

from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import stat
import tempfile
from typing import Any, Iterator

from .plan import BuildConfig, Job, PROJECT_ROOT, SCHEMA_VERSION


# ---------------------------------------------------------------------------
# Output boundaries
# ---------------------------------------------------------------------------


def safe_root(path: Path, config: BuildConfig) -> Path:
    """Reject destinations overlapping sources, installed data or the checkout.

    Parent aliases are resolved before comparison. The build directory itself
    and everything managed inside it must be ordinary paths, without symlinks.
    """
    path = path.expanduser().absolute()
    if path.is_symlink():
        raise ValueError("The output directory must not be a symlink")
    root = path.resolve()
    protected = (
        PROJECT_ROOT,
        Path.home() / ".libephemeris",
        Path(config.data_dir),
        Path(config.spk_dir),
        Path(config.assist_dir),
    )
    for source in protected:
        source = source.resolve()
        if root.is_relative_to(source) or source.is_relative_to(root):
            raise ValueError(f"Output overlaps a protected directory: {source}")
    if root.exists() and not root.is_dir():
        raise ValueError("Output must be a directory")
    if root.exists():
        for directory, folders, files in os.walk(root, followlinks=False):
            for name in (*folders, *files):
                item = Path(directory) / name
                if item.is_symlink() or not (item.is_dir() or item.is_file()):
                    raise ValueError(f"Unexpected output path type: {item}")
    return root


def managed_path(root: Path, relative: str) -> Path:
    """Resolve one build-owned path without permitting escapes or symlinks."""
    if root.is_symlink():
        raise ValueError("The build root was replaced with a symlink")
    part = Path(relative)
    if part.is_absolute() or ".." in part.parts or not part.parts:
        raise ValueError(f"Invalid managed path: {relative}")
    target = root / part
    for item in (target, *target.parents):
        if item == root:
            break
        if item.is_symlink():
            raise ValueError(f"Symlink in managed path: {item}")
    if not target.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Managed path escapes build directory: {relative}")
    if target.exists() and not target.is_file():
        raise ValueError(f"Managed file has an unexpected type: {target}")
    return target


# ---------------------------------------------------------------------------
# Durable writes and bounded hashes
# ---------------------------------------------------------------------------


def fsync_directory(directory: Path) -> None:
    """Persist directory entries after a rename on the local filesystem."""
    descriptor = os.open(directory, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_write(path: Path, contents: str) -> None:
    """Replace a small text file only after flush/fsync succeeds."""
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            delete=False,
        ) as stream:
            temporary = Path(stream.name)
            stream.write(contents)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        temporary = None
        fsync_directory(path.parent)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def file_stamp(path: Path) -> tuple[int, ...]:
    """Return metadata used for cheap between-job change detection."""
    info = path.stat()
    if not stat.S_ISREG(info.st_mode):
        raise ValueError(f"Expected a regular input file: {path}")
    return info.st_size, info.st_mtime_ns, info.st_ctime_ns, info.st_ino, info.st_dev


def sha256_file(path: Path) -> str:
    """Hash with bounded memory and reject concurrent file replacement."""
    before = file_stamp(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    if before != file_stamp(path):
        raise ValueError(f"File changed while hashing: {path}")
    return digest.hexdigest()


def save_manifest(root: Path, manifest: dict[str, Any]) -> None:
    """Persist the authoritative build state using strict JSON."""
    atomic_write(
        managed_path(root, "manifest.json"),
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
    )


def load_manifest(root: Path) -> dict[str, Any]:
    """Read a versioned checkpoint; job definitions are validated separately."""
    with managed_path(root, "manifest.json").open(encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, dict) or value.get("schema") != SCHEMA_VERSION:
        raise ValueError("Unsupported or invalid manifest schema")
    if any(
        not isinstance(value.get(key), dict)
        for key in ("config", "inputs", "spks", "jobs")
    ):
        raise ValueError("Missing or invalid manifest records")
    if value.get("status") not in (
        "pending",
        "running",
        "failed",
        "interrupted",
        "complete",
    ):
        raise ValueError("Invalid build status in manifest")
    return value


# ---------------------------------------------------------------------------
# Exclusive lifetime lock and dependency invalidation
# ---------------------------------------------------------------------------


@contextmanager
def build_lock(root: Path) -> Iterator[int]:
    """Hold a kernel lock and expose its descriptor for child inheritance.

    Never unlink the lock inode. A surviving child retains it after SIGKILL of
    its parent, so a second runner cannot start writing the same build.
    """
    path = managed_path(root, ".lock")
    descriptor = os.open(path, os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
    try:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError(
                "Another runner or surviving worker holds this build lock"
            ) from exc
        yield descriptor
    finally:
        # Closing, rather than LOCK_UN, preserves locks held by inherited fds.
        os.close(descriptor)


def validate_job_records(manifest: dict[str, Any], jobs: list[Job]) -> None:
    """Authenticate persisted paths/commands against the current canonical plan."""
    records = manifest.get("jobs")
    if not isinstance(records, dict) or set(records) != {job.id for job in jobs}:
        raise ValueError("Manifest job inventory does not match the build plan")
    for job in jobs:
        record = records[job.id]
        if not isinstance(record, dict) or record.get("definition") != job.to_dict():
            raise ValueError(f"Manifest job definition mismatch: {job.id}")
        if record.get("status") not in (
            "pending",
            "running",
            "done",
            "failed",
            "interrupted",
        ):
            raise ValueError(f"Invalid job state: {job.id}")
        attempt = record.get("attempt", 0)
        if type(attempt) is not int or attempt < 0:
            raise ValueError(f"Invalid attempt count: {job.id}")


def invalidate(manifest: dict[str, Any], jobs: list[Job], seeds: set[str]) -> set[str]:
    """Invalidate a phase and its descendants without deleting completed files.

    Old files remain reviewable; only recomputation plus validation can promote
    them again. Topological ordering makes the dependency walk a single pass.
    """
    if seeds - {job.id for job in jobs}:
        raise ValueError(f"Unknown job IDs: {sorted(seeds)}")
    affected = set(seeds)
    for job in jobs:
        if job.id in affected or affected.intersection(job.dependencies):
            affected.add(job.id)
            record = manifest["jobs"][job.id]
            record["status"] = "pending"
            record.pop("artifact", None)
    if affected:
        manifest["status"] = "pending"
    return affected


def staging_path(root: Path, job: Job, attempt: int) -> Path:
    """Return the one staging name belonging to a recorded job attempt."""
    output = Path(job.output)
    return managed_path(root, str(output.with_name(f".{output.name}.part-{attempt}")))
