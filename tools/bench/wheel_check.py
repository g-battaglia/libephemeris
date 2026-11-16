"""Provision isolated wheel environments and run the CI comparison."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import venv

from bitdiff import PROFILES
from common import (
    KERNEL_PINS, KERNELS, MODES, ROOT, SATELLITE_KERNEL, SCHEMA, TIERS,
    PrerequisiteError, command, digest, dump,
)


MEMORY_CASES = tuple(f"{mode}/{tier}" for mode in MODES for tier in TIERS)


def executables(path):
    windows = os.name == "nt"
    binaries = path / ("Scripts" if windows else "bin")
    suffix = ".exe" if windows else ""
    return binaries / f"python{suffix}", binaries / f"libephemeris{suffix}"


def provision_data(options, root, data, cli, memory_mode, memory_tier):
    shutil.copyfile(ROOT / "core/data/leb2/base_core.leb2", data / "leb/base_core.leb2")
    leb_tiers = set()
    spk_tiers = set()
    if options.profile == "satellites":
        spk_tiers.add("medium")
    elif options.profile != "smoke":
        leb_tiers.update(TIERS)
        spk_tiers.update(TIERS)
    if options.linux_memory or options.linux_startup:
        if memory_mode != "skyfield":
            leb_tiers.add(memory_tier)
        if memory_mode != "leb":
            spk_tiers.add(memory_tier)
    for tier in sorted(leb_tiers):
        print(f"Provisioning the pinned LEB {tier} group", flush=True)
        subprocess.run([str(cli), "download", f"leb2-{tier}", "--data-dir", str(data),
                        "--no-progress"], cwd=root, check=True)
    downloads = []
    for tier in sorted(spk_tiers):
        name = KERNELS[tier]
        downloads.extend((
            (name, f"https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/{name}"),
            (f"planet_centers_{tier}.bsp",
             "https://github.com/g-battaglia/libephemeris/releases/download/data-v2rc1/"
             f"planet_centers_{tier}.bsp"),
        ))
    if options.profile == "satellites":
        downloads.append((
            SATELLITE_KERNEL,
            "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/"
            f"{SATELLITE_KERNEL}",
        ))
    for file_name, url in downloads:
        print(f"Downloading and verifying {file_name}", flush=True)
        path = data / file_name
        with urllib.request.urlopen(url, timeout=120) as src, path.open("wb") as dest:
            shutil.copyfileobj(src, dest, length=1024 * 1024)
        if digest(path) != KERNEL_PINS[file_name]:
            raise PrerequisiteError(f"downloaded {file_name} does not match its pin")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=PROFILES, default="smoke")
    parser.add_argument("--candidate-cache-mib", type=int)
    parser.add_argument("--candidate-positional-centers", action="store_true")
    parser.add_argument("--linux-memory", action="store_true")
    parser.add_argument("--linux-startup", action="store_true")
    parser.add_argument("--linux-memory-case", choices=MEMORY_CASES, default="leb/base")
    parser.add_argument("--output", type=Path, default=ROOT / "target/bench/wheel-check")
    options = parser.parse_args()
    if (options.linux_memory or options.linux_startup) and sys.platform != "linux":
        parser.error("--linux-memory and --linux-startup require Linux")
    memory_mode, memory_tier = options.linux_memory_case.split("/")
    candidate_wheel, = (ROOT / "dist").glob("*.whl")
    with tempfile.TemporaryDirectory(prefix="ephemeris-wheel-check-") as tmp:
        root = Path(tmp)
        data = root / "data"
        (data / "leb").mkdir(parents=True)
        baseline_path = root / "baseline"
        candidate_path = root / "candidate"
        for path in (baseline_path, candidate_path):
            venv.create(path, with_pip=True)
        baseline, cli = executables(baseline_path)
        candidate, _ = executables(candidate_path)
        subprocess.run([str(baseline), "-m", "pip", "download", "--only-binary=:all:", "--no-deps",
                        "--dest", str(root), "libephemeris==4.0.0"], cwd=root, check=True)
        baseline_wheel, = root.glob("*.whl")
        for interpreter, wheel in ((baseline, baseline_wheel), (candidate, candidate_wheel)):
            subprocess.run([str(interpreter), "-m", "pip", "install", str(wheel), "numpy==2.4.6"],
                           cwd=root, check=True)
        try:
            provision_data(options, root, data, cli, memory_mode, memory_tier)
        except (PrerequisiteError, OSError, subprocess.CalledProcessError) as exc:
            options.output.mkdir(parents=True, exist_ok=True)
            report = dump({"schema": SCHEMA, "status": "missing_prerequisite", "error": str(exc)})
            (options.output / "comparison.json").write_text(report + "\n", encoding="utf-8")
            print(report, file=sys.stderr)
            return 2
        candidate_commit = command("git", "rev-parse", "HEAD")
        comparison_command = [
            sys.executable, str(ROOT / "tools/bench/bitdiff.py"), "run", "--profile", options.profile,
            "--data", str(data),
            "--baseline", str(baseline), "--baseline-wheel", str(baseline_wheel),
            "--baseline-commit", "python-v4.0.0",
        ]
        candidate_args = [
            *comparison_command, "--output", str(options.output.resolve()), "--candidate", str(candidate),
            "--candidate-wheel", str(candidate_wheel), "--candidate-commit", candidate_commit,
            "--candidate-build-options", "maturin build --release --locked (wheel matrix)",
        ]
        if options.candidate_cache_mib is not None:
            candidate_args.extend(("--candidate-cache-mib", str(options.candidate_cache_mib)))
        if options.candidate_positional_centers:
            candidate_args.append("--candidate-positional-centers")
        comparison = subprocess.run(candidate_args, cwd=root).returncode
        if options.profile == "auto-edge":
            control = subprocess.run([
                *comparison_command, "--output", str(options.output.resolve() / "published-control"),
                "--candidate", str(baseline), "--candidate-wheel", str(baseline_wheel),
                "--candidate-commit", "python-v4.0.0",
            ], cwd=root).returncode
            comparison = max(comparison, control)
        if comparison or not (options.linux_memory or options.linux_startup):
            return comparison
        diagnostics = 0
        for role, interpreter, wheel, commit in (
            ("baseline", baseline, baseline_wheel, "python-v4.0.0"),
            ("candidate", candidate, candidate_wheel, candidate_commit),
        ):
            benchmark = [
                str(interpreter), str(ROOT / "tools/bench/benchmark.py"),
                "--data", str(data), "--wheel", str(wheel), "--commit", commit,
                "--mode", memory_mode, "--tier", memory_tier, "--companions",
            ]
            if role == "candidate" and options.candidate_cache_mib is not None:
                benchmark.extend(("--leb-cache-mib", str(options.candidate_cache_mib)))
            if role == "candidate" and options.candidate_positional_centers:
                benchmark.append("--positional-centers")
            if options.linux_memory:
                subprocess.run([*benchmark, "--workload", "memory",
                                "--output", str(options.output.resolve() / f"{role}-memory.json")],
                               cwd=root, check=True)
            if options.linux_startup:
                for condition in ("first_use", "page_cache_warm", "page_cache_cold"):
                    startup = subprocess.run([
                        *benchmark, "--workload", "cold_start", "--cold-condition", condition,
                        "--output", str(options.output.resolve() / f"{role}-{condition}.json"),
                    ], cwd=root)
                    diagnostics = max(diagnostics, abs(startup.returncode))
        return diagnostics


if __name__ == "__main__":
    raise SystemExit(main())
