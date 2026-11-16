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
from common import KERNEL_PINS, ROOT, command, digest


def executables(path):
    windows = os.name == "nt"
    binaries = path / ("Scripts" if windows else "bin")
    suffix = ".exe" if windows else ""
    return binaries / f"python{suffix}", binaries / f"libephemeris{suffix}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=PROFILES, default="smoke")
    parser.add_argument("--candidate-cache-mib", type=int)
    parser.add_argument("--candidate-positional-centers", action="store_true")
    parser.add_argument("--leb-reader", choices=("owned", "indexed"))
    parser.add_argument("--linux-memory", action="store_true")
    parser.add_argument("--output", type=Path, default=ROOT / "target/bench/wheel-check")
    options = parser.parse_args()
    if options.linux_memory and sys.platform != "linux":
        parser.error("--linux-memory needs Linux /proc memory accounting")
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
        shutil.copyfile(ROOT / "core/data/leb2/base_core.leb2", data / "leb/base_core.leb2")
        if options.profile != "smoke":
            for tier, name in (("base", "de440s.bsp"), ("medium", "de440.bsp")):
                print(f"Provisioning the pinned LEB {tier} group", flush=True)
                subprocess.run([str(cli), "download", f"leb2-{tier}", "--data-dir", str(data),
                                "--no-progress"], cwd=root, check=True)
                for file_name, url in (
                    (name, f"https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/{name}"),
                    (f"planet_centers_{tier}.bsp",
                     "https://github.com/g-battaglia/libephemeris/releases/download/data-v2rc1/"
                     f"planet_centers_{tier}.bsp"),
                ):
                    print(f"Downloading and verifying {file_name}", flush=True)
                    path = data / file_name
                    with urllib.request.urlopen(url, timeout=120) as src, path.open("wb") as dest:
                        shutil.copyfileobj(src, dest, length=1024 * 1024)
                    if digest(path) != KERNEL_PINS[file_name]:
                        raise ValueError(f"downloaded {file_name} does not match its pin")
        candidate_commit = command("git", "rev-parse", "HEAD")
        comparison_command = [
            sys.executable, str(ROOT / "tools/bench/bitdiff.py"), "run", "--profile", options.profile,
            "--data", str(data),
            "--baseline", str(baseline), "--baseline-wheel", str(baseline_wheel),
            "--baseline-commit", "python-v4.0.0",
        ]
        if options.leb_reader:
            comparison_command.extend(("--leb-reader", options.leb_reader))
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
            return max(comparison, control)
        if comparison or not options.linux_memory:
            return comparison
        for role, interpreter, wheel, commit in (
            ("baseline", baseline, baseline_wheel, "python-v4.0.0"),
            ("candidate", candidate, candidate_wheel, candidate_commit),
        ):
            benchmark = [
                str(interpreter), str(ROOT / "tools/bench/benchmark.py"),
                "--data", str(data), "--wheel", str(wheel), "--commit", commit,
                "--workload", "memory", "--mode", "leb", "--tier", "base",
                "--output", str(options.output.resolve() / f"{role}-memory.json"),
            ]
            if options.leb_reader:
                benchmark.extend(("--leb-reader", options.leb_reader))
            if role == "candidate" and options.candidate_cache_mib is not None:
                benchmark.extend(("--leb-cache-mib", str(options.candidate_cache_mib)))
            subprocess.run(benchmark, cwd=root, check=True)
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
