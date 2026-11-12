#!/usr/bin/env python3
"""Copy resolved Cargo dependencies' actual license and attribution files.

Run before packaging, then use --check to detect stale or incomplete copies.
Cargo metadata is read offline with the lockfile unchanged. The output contains
upstream bytes and Cargo's name, version, license expression and source metadata;
it does not select among alternative licenses or certify legal compliance.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
SCHEMA = "cargo-dependency-licenses/v1"
NOTICE_NAME = re.compile(
    r"^(?:licen[cs]e|copying|copyright|notice|authors|contributors|unlicen[cs]e)"
    r"(?:[-_.].*)?$",
    re.IGNORECASE,
)
LICENSE_DIRECTORIES = {"license", "licenses", "licence", "licences"}


def cargo_metadata(metadata_file: Path | None) -> dict:
    """Read resolved packages without compiling, fetching or modifying the lockfile."""
    if metadata_file is not None:
        return json.loads(metadata_file.read_text(encoding="utf-8"))
    result = subprocess.run(
        [
            "cargo", "metadata", "--locked", "--offline", "--all-features",
            "--format-version", "1",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(result.stdout)


def notice_files(package: dict) -> list[tuple[str, bytes]]:
    """Return upstream license/notice bytes with their original package-relative names."""
    root = Path(package["manifest_path"]).resolve().parent
    candidates = set()
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        relative = path.relative_to(root)
        if NOTICE_NAME.fullmatch(path.name) or any(
            part.lower() in LICENSE_DIRECTORIES for part in relative.parts[:-1]
        ):
            candidates.add(relative)
    if package.get("license_file"):
        path = Path(package["license_file"])
        path = path if path.is_absolute() else root / path
        if not path.resolve().is_relative_to(root):
            raise ValueError(f"{package['name']}: license_file escapes the package")
        if not path.is_file():
            raise ValueError(f"{package['name']}: declared license_file is absent")
        candidates.add(path.relative_to(root))
    if not candidates:
        raise ValueError(
            f"{package['name']} {package['version']}: no upstream license/notice "
            "file found; supply the actual upstream text before packaging"
        )
    files = []
    for relative in sorted(candidates):
        source = root / relative
        if not source.resolve().is_relative_to(root):
            raise ValueError(f"{package['name']}: notice symlink escapes the package")
        content = source.read_bytes()
        if not content:
            raise ValueError(f"{package['name']}: empty upstream notice {relative}")
        files.append((relative.as_posix(), content))
    return files


def bundle_contents(metadata: dict) -> tuple[dict[str, bytes], int]:
    """Build deterministic package metadata and copies, excluding workspace-owned code."""
    workspace = set(metadata["workspace_members"])
    resolved = {node["id"] for node in metadata["resolve"]["nodes"]}
    packages = sorted(
        (p for p in metadata["packages"] if p["id"] in resolved - workspace),
        key=lambda p: (p["name"], p["version"], p.get("source") or ""),
    )
    files: dict[str, bytes] = {}
    entries = []
    directories = set()
    for package in packages:
        directory = f"{package['name']}-{package['version']}"
        if directory in directories:
            raise ValueError(f"duplicate name/version from different sources: {directory}")
        directories.add(directory)
        copied = []
        for relative, content in notice_files(package):
            destination = f"{directory}/{relative}"
            files[destination] = content
            copied.append({"source_path": relative, "path": destination})
        entries.append({
            "name": package["name"],
            "version": package["version"],
            "license": package.get("license"),
            "source": package.get("source"),
            "repository": package.get("repository"),
            "files": copied,
        })
    files["index.json"] = (
        json.dumps({"schema": SCHEMA, "packages": entries}, indent=2, sort_keys=True)
        + "\n"
    ).encode("utf-8")
    return files, len(packages)


def current_contents(output: Path) -> dict[str, bytes]:
    """Read an existing generated bundle without following output symlinks."""
    if output.is_symlink():
        raise ValueError("the output directory must not be a symlink")
    if not output.exists():
        return {}
    current = {}
    for path in output.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"output symlink is not supported: {path}")
        if path.is_file():
            current[path.relative_to(output).as_posix()] = path.read_bytes()
    return current


def write_bundle(output: Path, expected: dict[str, bytes], current: dict[str, bytes]) -> None:
    """Update only this tool's directory and remove superseded dependency notices."""
    if current:
        try:
            schema = json.loads(current["index.json"])["schema"]
        except (KeyError, ValueError, TypeError) as error:
            raise ValueError("refusing to replace an unrecognized nonempty directory") from error
        if schema != SCHEMA:
            raise ValueError("refusing to replace a different output schema")
    output.mkdir(parents=True, exist_ok=True)
    for relative, content in sorted(expected.items()):
        if current.get(relative) != content:
            target = output / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(content)
    for relative in current.keys() - expected.keys():
        (output / relative).unlink()
    for directory in sorted(
        (p for p in output.rglob("*") if p.is_dir()),
        key=lambda p: len(p.parts), reverse=True,
    ):
        if not any(directory.iterdir()):
            directory.rmdir()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, help="read saved Cargo metadata instead")
    parser.add_argument("--output", type=Path, default=ROOT / "core/licenses")
    parser.add_argument("--check", action="store_true", help="compare without writing")
    arguments = parser.parse_args()
    try:
        expected, package_count = bundle_contents(cargo_metadata(arguments.metadata))
        current = current_contents(arguments.output)
        differences = sorted(
            name for name in current.keys() | expected.keys()
            if current.get(name) != expected.get(name)
        )
        if arguments.check and differences:
            print("Dependency license bundle differs: " + ", ".join(differences[:10]), file=sys.stderr)
            return 1
        if not arguments.check:
            write_bundle(arguments.output, expected, current)
        action = "Checked" if arguments.check else "Collected"
        print(f"{action} upstream license/notice files for {package_count} Cargo packages.")
        return 0
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        print(f"License collection failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
