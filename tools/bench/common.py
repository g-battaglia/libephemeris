"""Shared inputs and environment for the offline performance tools."""

import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import struct
import subprocess
import sys
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = 1
MODES = ("leb", "auto", "skyfield")
TIERS = ("base", "medium")
LEB_TIERS = (*TIERS, "extended")
IERS_FILES = ("deltat.data", "finals2000A.data", "leap_seconds.dat")
SEED = 20261116
J2000 = 2451545.0
SITE = (12.5, 41.9, 50.0)
MANIFEST = ROOT / "core/data/leb2/data-v4-manifest.json"
KERNELS = {"base": "de440s.bsp", "medium": "de440.bsp"}
SATELLITE_KERNEL = "plu017.bsp"
KERNEL_PINS = {
    "de440s.bsp": "c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2",
    "planet_centers_base.bsp": "a9ec744ff412b095129166587ea0814f81c850faebf92586a738cb5dc103c92a",
    "planet_centers_medium.bsp": "d3c34f5efe9223ef588ec59a8c59c1bd6619b0eab5d5e0b35c353d675efe7b4d",
    "de440.bsp": "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5",
    SATELLITE_KERNEL: "419eb90ad2ce723f8f61e2416c2ce9f01307de51d492ba268926a30cdf359bc0",
}
RUST_SOURCE = re.compile(rb"[/\\]rustc[/\\]([0-9a-f]{40})[/\\]")


class PrerequisiteError(Exception):
    """Required data or runner unavailable."""


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def rust_source_revisions(path):
    """Identify remapped Rust source paths embedded in a compiled binary."""
    revisions = set()
    tail = b""
    with Path(path).open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            contents = tail + chunk
            revisions.update(match.decode("ascii") for match in RUST_SOURCE.findall(contents))
            tail = contents[-48:]
    return sorted(revisions)


def dump(value):
    return json.dumps(value, allow_nan=False, separators=(",", ":"), sort_keys=True)


def command(*args):
    try:
        return subprocess.check_output(args, cwd=ROOT, text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def identity(wheel=None, commit=None, build_options=None):
    import libephemeris as eph
    from libephemeris import _core

    if not Path(eph.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()):
        raise PrerequisiteError("libephemeris must come from the selected interpreter's install")
    extension = Path(_core.__file__)
    extension_sha = digest(extension)
    wheel_info = None
    if wheel is not None:
        with ZipFile(wheel) as archive, archive.open(f"libephemeris/{extension.name}") as stream:
            if hashlib.file_digest(stream, "sha256").hexdigest() != extension_sha:
                raise PrerequisiteError("installed extension differs from the supplied wheel")
        wheel_info = {"name": wheel.name, "sha256": digest(wheel)}
    return {
        "version": eph.__version__,
        "extension_sha256": extension_sha,
        "python": sys.version.split()[0],
        "interpreter": sys.executable,
        "platform": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "libc": platform.libc_ver(),
        "cpu": command("sysctl", "-n", "machdep.cpu.brand_string") or platform.processor(),
        "rustc_on_host": command("rustc", "-Vv"),
        "rust_source_revisions": rust_source_revisions(extension),
        "commit": commit,
        "build_options": build_options,
        "wheel": wheel_info,
    }


def environment(data):
    env = {key: value for key, value in os.environ.items()
           if not key.startswith(("LIBEPHEMERIS_", "PYTHON"))}
    env.update(
        LIBEPHEMERIS_DATA_DIR=str(data),
        LIBEPHEMERIS_CONFIG=str(data / "absent.toml"),
        LIBEPHEMERIS_ENV_FILE=str(data / "absent.env"),
        PYTHONUTF8="1",
    )
    return env


def shared_medium_identity():
    """Read the device and inode of the actual sealed mapping on Linux."""
    for line in Path("/proc/self/maps").read_text().splitlines():
        fields = line.split(maxsplit=5)
        if len(fields) == 6 and fields[5] == "/memfd:libephemeris-medium (deleted)":
            return fields[3], int(fields[4])
    raise PrerequisiteError("worker did not map the requested sealed LEB image")


def configure(eph, mode, tier, *, iers=False):
    eph.set_network_policy("sealed")
    eph.set_calc_mode(mode)
    eph.set_precision_tier(tier)
    eph.set_iers_delta_t_enabled(iers)
    eph.set_auto_spk_download(False)
    eph.set_delta_t_model("smh2016")
    eph.set_sid_mode(0)
    eph.set_topo(*SITE)


def leb_index(path):
    with path.open("rb") as stream:
        magic, version, count, bodies, start, stop, _, _ = struct.unpack(
            "<4sIIIdddI20x", stream.read(64)
        )
        if magic != b"LEB2" or version != 2 or count > 256 or bodies > 10000:
            raise PrerequisiteError(f"unsupported LEB2 header: {path}")
        sections = [struct.unpack("<IIQQ", stream.read(24)) for _ in range(count)]
        offset, length = next((o, size) for kind, _, o, size in sections if kind == 0)
        if length != bodies * 68:
            raise PrerequisiteError(f"invalid LEB2 index length: {path}")
        stream.seek(offset)
        entries = [struct.unpack("<iIIdddIIQQQ", stream.read(68)) for _ in range(bodies)]
        result = []
        for body, _, segments, lower, upper, step, _, _, chunk_offset, _, _ in entries:
            stream.seek(chunk_offset)
            chunks, = struct.unpack("<I", stream.read(4))
            if chunks > segments:
                raise PrerequisiteError(f"invalid LEB2 chunk count: {path}")
            stream.seek(chunk_offset + 16)
            starts = [struct.unpack("<ddQQQII", stream.read(48))[0] for _ in range(chunks)]
            rows = [lower + k * step for k in {1, segments // 2, segments - 1} if 0 < k < segments]
            result.append({"body": body, "start": lower, "stop": upper, "step": step,
                           "segments": segments, "boundaries": sorted(set([lower, upper, *rows, *starts]))})
    if not math.isfinite(start) or not math.isfinite(stop) or start >= stop:
        raise PrerequisiteError(f"invalid LEB2 coverage: {path}")
    return {"range": [start, stop], "bodies": result,
            "decoded_coefficient_bytes": sum(entry[-1] for entry in entries)}


def spk_index(path):
    """Record the actual DAF segment intervals in TDB Julian days."""
    with path.open("rb") as stream:
        header = stream.read(1024)
        order = {b"LTL-IEEE": "<", b"BIG-IEEE": ">"}.get(header[88:96])
        if not header.startswith(b"DAF/SPK") or order is None:
            raise PrerequisiteError(f"unsupported SPK header: {path}")
        nd, ni = struct.unpack_from(order + "ii", header, 8)
        if (nd, ni) != (2, 6):
            raise PrerequisiteError(f"unsupported SPK summary size: {path}")
        record, = struct.unpack_from(order + "i", header, 76)
        visited = set()
        segments = []
        while record:
            if record in visited or not 1 <= record <= path.stat().st_size // 1024:
                raise PrerequisiteError(f"invalid SPK summary chain: {path}")
            visited.add(record)
            stream.seek((record - 1) * 1024)
            raw = stream.read(1024)
            next_record, _, count = struct.unpack_from(order + "ddd", raw)
            if not 0 <= count <= 25 or count != int(count):
                raise PrerequisiteError(f"invalid SPK summary count: {path}")
            for i in range(int(count)):
                lower, upper, target, center, _, _, _, _ = struct.unpack_from(
                    order + "ddiiiiii", raw, 24 + i * 40
                )
                segments.append([target, center, J2000 + lower / 86400, J2000 + upper / 86400])
            record = int(next_record)
        if not segments:
            raise PrerequisiteError(f"empty SPK index: {path}")
        return segments


def prepare(source, destination, modes, tiers, *, companions=False, satellites=False):
    """Stage a whitelist so unrelated local files cannot affect routing."""
    manifest = json.loads(MANIFEST.read_text())["files"]
    names = {SATELLITE_KERNEL} if satellites else set()
    for tier in tiers:
        if {"leb", "auto"}.intersection(modes):
            names.add(f"leb/{tier}_core.leb2")
            if companions:
                names.update(f"leb/{tier}_{group}.leb2" for group in ("apogee", "asteroids", "exotics"))
        if {"skyfield", "auto"}.intersection(modes):
            if tier not in KERNELS:
                raise PrerequisiteError(f"no pinned SPK inventory for {tier}; use LEB mode")
            names.update((KERNELS[tier], f"planet_centers_{tier}.bsp"))

    inventory = {}
    for name in sorted(names):
        path = source / name
        if not path.is_file():
            raise PrerequisiteError(f"missing input: {path}")
        target = destination / name
        target.parent.mkdir(exist_ok=True, parents=True)
        # A hard link changes the original inode's ctime, including while
        # another benchmark is using it. The engine tracks that metadata on
        # its SPK image hints. Each run must own independent input inodes.
        shutil.copyfile(path, target)
        sha = digest(target)
        if name.endswith(".leb2"):
            expected = manifest[path.name]["sha256"]
            if sha != expected:
                raise PrerequisiteError(f"data-v4 SHA-256 mismatch: {path}")
            coverage = leb_index(target)
        else:
            if name in KERNEL_PINS and sha != KERNEL_PINS[name]:
                raise PrerequisiteError(f"SPK SHA-256 mismatch: {path}")
            coverage = {"segments_tdb": spk_index(target)}
        inventory[name] = {"path": str(path.resolve()), "sha256": sha, "bytes": target.stat().st_size,
                           "coverage": coverage}
    return inventory


def prepare_iers(source, data):
    destination = data / "iers_cache"
    destination.mkdir(parents=True)
    inventory = {}
    for name in IERS_FILES:
        path = source / name
        if not path.is_file() or not path.stat().st_size:
            raise PrerequisiteError(f"missing IERS input: {path}")
        target = destination / name
        shutil.copyfile(path, target)
        inventory[name] = {"sha256": digest(target), "bytes": target.stat().st_size}
    return inventory
