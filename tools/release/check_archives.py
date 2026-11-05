#!/usr/bin/env python3
"""Inspect local crate/wheel/sdist delivery boundaries without extracting files.

This does not publish, verify numerical accuracy or establish cross-platform
support. It rejects development debris and incomplete metadata; Cargo, Maturin
and clean installed/source builds remain separate release checks.
"""

from __future__ import annotations

import argparse
from email.parser import BytesParser
from pathlib import Path, PurePosixPath
import tarfile
import tomllib
import zipfile

# crates.io config/publish_limits.rs defines the default as 10 * 1024 * 1024.
# Check actual compressed archives, not rounded display sizes or extracted bytes.
CRATE_LIMIT = 10 * 1024 * 1024
PROCESS_FILES = {"AGENTS.md", "PLAN.md", "PROGRESS.md"}
CACHE_PARTS = {".git", "target", "dist", "build", "_build", "__pycache__", ".pytest_cache"}


def path_errors(names: list[str]) -> list[str]:
    """Check member paths before any extraction, including ignored build debris."""
    errors = []
    for name in names:
        path = PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts or "\\" in name:
            errors.append(f"unsafe member path: {name}")
        elif path.name in PROCESS_FILES:
            errors.append(f"process file in delivery: {name}")
        elif any(part in CACHE_PARTS or part.startswith(".venv")
                 or part.endswith(".egg-info") for part in path.parts):
            errors.append(f"build/cache member in delivery: {name}")
        elif "validation" in path.parts or "libephemeris-validation-python" in path.parts:
            errors.append(f"local validation repository in delivery: {name}")
        elif any(a == "native" and b in {"rebound", "assist"}
                 for a, b in zip(path.parts, path.parts[1:])):
            errors.append(f"retired implementation in delivery: {name}")
    return errors


def metadata_errors(content: bytes) -> list[str]:
    """Require useful rendered project metadata and its exact current grant."""
    metadata = BytesParser().parsebytes(content)
    errors = []
    if metadata.get("License-Expression") != "AGPL-3.0":
        errors.append("missing or changed project License-Expression")
    if metadata.get("License"):
        errors.append("legacy License and License-Expression must not coexist")
    if metadata.get_content_type() != "text/plain":
        # The package payload is a metadata document, not a MIME multipart body.
        errors.append("unexpected multipart metadata")
    # Maturin includes valid charset/variant parameters on the media type.
    description_type = metadata.get("Description-Content-Type", "").split(";", 1)[0].strip().lower()
    if description_type != "text/markdown":
        errors.append("project README must be declared as Markdown")
    if not metadata.get_payload().strip():
        errors.append("missing project README payload")
    return errors


def check(path: Path) -> list[str]:
    """Inspect a built archive, reporting all discovered boundary violations."""
    errors = []
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            errors += path_errors(names)
            required = {"libephemeris/__init__.py", "libephemeris/py.typed"}
            errors += [f"missing wheel member: {name}" for name in sorted(required - set(names))]
            metadata = [name for name in names if name.endswith(".dist-info/METADATA")]
            if len(metadata) != 1:
                errors.append("wheel must have exactly one distribution METADATA")
            else:
                errors += metadata_errors(archive.read(metadata[0]))
            if not any(name.startswith("libephemeris/_core.")
                       and name.endswith((".so", ".pyd")) for name in names):
                errors.append("wheel has no native _core extension")
            if not any(name.endswith("/licenses/core/licenses/index.json") for name in names):
                errors.append("wheel has no upstream dependency license index")
    else:
        with tarfile.open(path, "r:gz") as archive:
            members = archive.getmembers()
            names = [member.name for member in members]
            errors += path_errors(names)
            errors += [f"archive link is not a standalone input: {m.name}"
                       for m in members if m.issym() or m.islnk()]
            roots = {PurePosixPath(name).parts[0] for name in names if PurePosixPath(name).parts}
            if len(roots) != 1:
                errors.append("source archive must have exactly one root directory")
            relative_names = {PurePosixPath(*PurePosixPath(name).parts[1:]).as_posix() for name in names}
            required = {"Cargo.toml", "COPYING", "THIRD_PARTY_NOTICES.md"}
            if path.suffix == ".crate":
                required.add("src/lib.rs")
                manifest_member = next(member for member in members
                                       if member.name.endswith("/Cargo.toml"))
                with archive.extractfile(manifest_member) as stream:
                    manifest = tomllib.loads(stream.read().decode("utf-8"))
                if manifest["package"]["name"] == "libephemeris-core-data":
                    required.add("leb2/base_core.leb2")
            else:
                required |= {"pyproject.toml", "python/Cargo.toml", "core/Cargo.toml",
                             "core/src/lib.rs", "core/data/leb2/base_core.leb2",
                             "tools/release/check_minimal.py", "tools/leb_generation/pyproject.toml",
                             "docs/installation.md", "docs/conf.py", "docs/index.rst",
                             "beginners-manual/README.md", "beginners-manual/COPYING",
                             "beginners-manual/COPYING.LESSER"}
                # Require both complete books, not merely their landing pages.
                for language in ("it", "en"):
                    required.add(f"beginners-manual/{language}/README.md")
                    required.add(f"beginners-manual/{language}/metadata.yaml")
                    chapters = {PurePosixPath(name).name[:2] for name in relative_names
                                if name.startswith(f"beginners-manual/{language}/")
                                and name.endswith(".md")}
                    for number in range(16):
                        if f"{number:02d}" not in chapters:
                            errors.append(f"missing manual chapter: {language}/{number:02d}")
            errors += [f"missing source member: {name}" for name in sorted(required - relative_names)]
            if path.suffix == ".crate":
                if path.stat().st_size > CRATE_LIMIT:
                    errors.append(f"crate is {path.stat().st_size} bytes; default registry limit is {CRATE_LIMIT}")
            else:
                metadata = [member for member in members if member.name.endswith("/PKG-INFO")]
                if len(metadata) != 1:
                    errors.append("sdist must have exactly one project PKG-INFO")
                else:
                    stream = archive.extractfile(metadata[0])
                    if stream is None:
                        errors.append("sdist PKG-INFO is not a regular file")
                    else:
                        with stream:
                            errors += metadata_errors(stream.read())
    return errors


def main() -> int:
    """Check explicit local artifacts only; never search unrelated directories."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archives", type=Path, nargs="+")
    arguments = parser.parse_args()
    failed = False
    for path in arguments.archives:
        try:
            errors = check(path)
        except (OSError, ValueError, tarfile.TarError, zipfile.BadZipFile) as error:
            errors = [str(error)]
        for error in errors:
            print(f"{path}: {error}")
        if not errors:
            print(f"{path}: delivery boundaries passed")
        failed |= bool(errors)
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
