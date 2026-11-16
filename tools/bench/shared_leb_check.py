"""Check sealed LEB snapshot lifetime in separate installed-wheel processes."""

import argparse
from contextlib import contextmanager
import json
import os
from pathlib import Path
import selectors
import shutil
import struct
import subprocess
import sys
import tempfile

from common import J2000, PrerequisiteError, configure, dump, environment, shared_medium_identity


def observe():
    import libephemeris as eph

    configure(eph, "leb", "medium")
    values, flags = eph.calc_ut(J2000, 4, 258)
    try:
        image = shared_medium_identity()
    except PrerequisiteError:
        image = None
    print(dump({"bits": struct.pack("<6di", *values, flags).hex(), "image": image}), flush=True)
    sys.stdin.readline()


@contextmanager
def worker(env):
    with subprocess.Popen(
        [sys.executable, "-E", "-s", str(Path(__file__).resolve()), "--worker"], env=env,
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True,
    ) as process:
        try:
            with selectors.DefaultSelector() as selector:
                selector.register(process.stdout, selectors.EVENT_READ)
                if not selector.select(timeout=15):
                    raise PrerequisiteError("shared LEB worker did not become ready")
                line = process.stdout.readline()
            if not line:
                raise PrerequisiteError("shared LEB worker exited before the first position")
            yield process, json.loads(line)
        finally:
            if process.poll() is None:
                try:
                    release(process)
                except (subprocess.TimeoutExpired, BrokenPipeError):
                    process.kill()


def release(process):
    process.stdin.close()
    if process.wait(timeout=5) != 0:
        raise PrerequisiteError("shared LEB worker could not shut down")


def check(data):
    import fcntl

    with tempfile.TemporaryDirectory(prefix="ephemeris-shared-verify-") as tmp:
        registry = Path(tmp) / "registry"
        registry.mkdir(mode=0o700)
        ordinary = environment(data)
        shared = dict(ordinary, LIBEPHEMERIS_LEB_SHARED_DIR=str(registry))
        with worker(ordinary) as (_, baseline):
            if baseline["image"] is not None:
                raise PrerequisiteError("ordinary worker mapped a shared LEB image")
        with worker(shared) as (creator, first), worker(shared) as (_, peer):
            if first["image"] is None or first != peer or first["bits"] != baseline["bits"]:
                raise PrerequisiteError("independent workers failed to share the pinned image")
            release(creator)
            with worker(shared) as (_, replacement):
                if replacement != first:
                    raise PrerequisiteError("shared LEB image was lost when its creator exited")
        with worker(shared) as (_, fresh):
            if fresh["image"] is None or fresh["bits"] != baseline["bits"]:
                raise PrerequisiteError("worker failed to recreate the image after all peers exited")
        with (registry / "medium-v1.lock").open("rb") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            with worker(shared) as (_, blocked):
                if blocked != baseline:
                    raise PrerequisiteError("busy shared LEB lock changed the owned fallback")

        fd = os.memfd_create("libephemeris-medium", os.MFD_CLOEXEC | os.MFD_ALLOW_SEALING)
        try:
            with (data / "leb/medium_core.leb2").open("rb") as src, os.fdopen(os.dup(fd), "wb") as dst:
                shutil.copyfileobj(src, dst)
            offset = os.fstat(fd).st_size - 64
            byte, = os.pread(fd, 1, offset)
            os.pwrite(fd, bytes([byte ^ 1]), offset)
            seals = fcntl.F_SEAL_WRITE | fcntl.F_SEAL_SHRINK | fcntl.F_SEAL_GROW | fcntl.F_SEAL_SEAL
            fcntl.fcntl(fd, fcntl.F_ADD_SEALS, seals)
            (registry / f"medium-v1-peer-{os.getpid()}-{fd}").touch(mode=0o600)
            with worker(shared) as (_, corrupt):
                if corrupt != baseline:
                    raise PrerequisiteError("incorrect sealed bytes changed the owned fallback")
        finally:
            os.close(fd)
    return {"status": "complete", "creator_exit": "retained", "last_holder_exit": "recreated",
            "busy_lock": "owned fallback", "incorrect_sealed_image": "owned fallback"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--data", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if sys.platform != "linux":
        parser.error("shared LEB lifecycle qualification requires Linux")
    if args.worker:
        observe()
        return 0
    if not args.data or not args.output:
        parser.error("provide --data and --output")
    try:
        report = check(args.data.resolve())
    except (OSError, ValueError, PrerequisiteError) as error:
        report = {"status": "incomplete", "error": str(error)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(dump(report) + "\n", encoding="utf-8")
    print(dump(report))
    return 0 if report["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
