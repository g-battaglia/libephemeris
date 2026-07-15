#!/usr/bin/env python3
"""Build the distribution artifacts and audit their contents.

Licensing/packaging gate: asserts neither the sdist nor the wheel contains
test code, the dev CLI, or Swiss Ephemeris reference data, and that the
wheel carries the expected payload (package, typing marker, bundled data,
vendored modules).

The wheel is built FROM THE SDIST (python -m build's default two-step), so
the audit is hermetic against stale in-tree ``build/`` directories — a
stale ``build/lib`` is exactly how ``dev_cli/`` once leaked into a wheel.

Usage:
    python scripts/check_wheel_contents.py            # build sdist+wheel
    python scripts/check_wheel_contents.py path.whl   # audit existing wheel
"""

from __future__ import annotations

import argparse
import base64
import configparser
import csv
import hashlib
import io
import math
import os
import re
import stat
import subprocess
import sys
import tarfile
import tempfile
import tomllib
import zipfile
from collections import Counter
from email import policy
from email.parser import BytesParser
from pathlib import Path

from packaging.requirements import InvalidRequirement, Requirement

try:
    from .check_provenance import (
        ALLOWLIST,
        LICENSE_NAMING_OK,
        _classify_artifact_name,
        _canonical_policy_rel,
        _is_binary_payload_name,
        _line_hit_labels,
    )
except ImportError:  # Direct execution: ``python scripts/check_wheel_contents.py``
    from check_provenance import (  # type: ignore[no-redef]
        ALLOWLIST,
        LICENSE_NAMING_OK,
        _classify_artifact_name,
        _canonical_policy_rel,
        _is_binary_payload_name,
        _line_hit_labels,
    )

REPO_ROOT = Path(__file__).resolve().parent.parent
_PROJECT_CONFIG = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())
_PROJECT_METADATA = _PROJECT_CONFIG["project"]
_SETUPTOOLS_CONFIG = _PROJECT_CONFIG["tool"]["setuptools"]
PROJECT_NAME = _PROJECT_METADATA["name"]
PROJECT_VERSION = _PROJECT_METADATA["version"]
PROJECT_LICENSE = _PROJECT_METADATA["license"]
EXPECTED_SDIST_ROOT = f"{PROJECT_NAME}-{PROJECT_VERSION}"
EXPECTED_DIST_INFO = f"{PROJECT_NAME}-{PROJECT_VERSION}.dist-info"
EXPECTED_EGG_INFO = f"{PROJECT_NAME}.egg-info"
EXPECTED_SDIST_FILENAME = f"{EXPECTED_SDIST_ROOT}.tar.gz"
EXPECTED_WHEEL_FILENAME = f"{PROJECT_NAME}-{PROJECT_VERSION}-py3-none-any.whl"
CANONICAL_EXPECTED_SDIST_ROOT = _canonical_policy_rel(EXPECTED_SDIST_ROOT)
DEFAULT_ALLOWED_METADATA_PATHS = frozenset(
    {
        EXPECTED_DIST_INFO,
        f"{EXPECTED_SDIST_ROOT}/{EXPECTED_EGG_INFO}",
    }
)

FORBIDDEN = (
    re.compile(r"(^|/)tests?/", re.IGNORECASE),
    re.compile(r"(^|/)compare_scripts/", re.IGNORECASE),
    re.compile(r"libephemeris/dev_cli(/|\.py)", re.IGNORECASE),
    re.compile(r"\.se1$", re.IGNORECASE),
    re.compile(r"sefstars", re.IGNORECASE),
    re.compile(r"seorbel", re.IGNORECASE),
    re.compile(r"seleapsec", re.IGNORECASE),
    re.compile(r"sedeltat", re.IGNORECASE),
)


def _discover_wheel_project_files() -> tuple[str, ...]:
    """Derive the complete expected project payload from setuptools config."""
    result = subprocess.run(
        ["git", "ls-files", "-z", "--", "libephemeris"],
        cwd=REPO_ROOT,
        capture_output=True,
        check=True,
    )
    tracked = {
        Path(os.fsdecode(item)).as_posix()
        for item in result.stdout.split(b"\0")
        if item
    }
    files: set[str] = set()
    for package in _SETUPTOOLS_CONFIG["packages"]:
        package_root = REPO_ROOT / package.replace(".", "/")
        files.update(
            path.relative_to(REPO_ROOT).as_posix()
            for path in package_root.glob("*.py")
            if path.is_file() and path.relative_to(REPO_ROOT).as_posix() in tracked
        )
    for package, patterns in _SETUPTOOLS_CONFIG["package-data"].items():
        package_root = REPO_ROOT / package.replace(".", "/")
        for pattern in patterns:
            files.update(
                path.relative_to(REPO_ROOT).as_posix()
                for path in package_root.glob(pattern)
                if path.is_file() and path.relative_to(REPO_ROOT).as_posix() in tracked
            )
    return tuple(sorted(files))


WHEEL_REQUIRED = _discover_wheel_project_files()
WHEEL_METADATA_REQUIRED = (
    f"{EXPECTED_DIST_INFO}/METADATA",
    f"{EXPECTED_DIST_INFO}/WHEEL",
    f"{EXPECTED_DIST_INFO}/RECORD",
    f"{EXPECTED_DIST_INFO}/entry_points.txt",
    f"{EXPECTED_DIST_INFO}/top_level.txt",
    f"{EXPECTED_DIST_INFO}/licenses/LICENSE",
    f"{EXPECTED_DIST_INFO}/licenses/LICENSING.md",
    f"{EXPECTED_DIST_INFO}/licenses/NOTICE.md",
    f"{EXPECTED_DIST_INFO}/licenses/THIRD_PARTY_NOTICES.md",
)
SDIST_ROOT_REQUIRED = (
    "LICENSE",
    "LICENSING.md",
    "NOTICE.md",
    "THIRD_PARTY_NOTICES.md",
    "PKG-INFO",
    "pyproject.toml",
    "MANIFEST.in",
    "README.md",
)
SDIST_GENERATED_REQUIRED = (
    f"{EXPECTED_EGG_INFO}/PKG-INFO",
    f"{EXPECTED_EGG_INFO}/SOURCES.txt",
    f"{EXPECTED_EGG_INFO}/dependency_links.txt",
    f"{EXPECTED_EGG_INFO}/entry_points.txt",
    f"{EXPECTED_EGG_INFO}/requires.txt",
    f"{EXPECTED_EGG_INFO}/top_level.txt",
    "setup.cfg",
)
SDIST_REQUIRED = (
    *SDIST_ROOT_REQUIRED,
    *WHEEL_REQUIRED,
    *SDIST_GENERATED_REQUIRED,
)
LEGAL_FILENAMES = (
    "LICENSE",
    "LICENSING.md",
    "NOTICE.md",
    "THIRD_PARTY_NOTICES.md",
)
SDIST_BYTE_REQUIRED = (
    *LEGAL_FILENAMES,
    "pyproject.toml",
    "MANIFEST.in",
    "README.md",
    *WHEEL_REQUIRED,
)
EXPECTED_WHEEL_NAMES = frozenset((*WHEEL_REQUIRED, *WHEEL_METADATA_REQUIRED))
EXPECTED_SDIST_FILE_NAMES = frozenset(
    f"{EXPECTED_SDIST_ROOT}/{name}" for name in SDIST_REQUIRED
)
EXPECTED_SDIST_DIRECTORY_NAMES = frozenset(
    {
        EXPECTED_SDIST_ROOT,
        *(
            "/".join(parts[:index])
            for name in EXPECTED_SDIST_FILE_NAMES
            for parts in (name.split("/"),)
            for index in range(1, len(parts))
        ),
    }
)
EXPECTED_ENTRY_POINTS = (
    "[console_scripts]\n"
    + "".join(
        f"{name} = {target}\n"
        for name, target in sorted(_PROJECT_METADATA["scripts"].items())
    )
).encode()


def _expected_requires_dist() -> tuple[str, ...]:
    """Return normalized dependency metadata derived from pyproject.toml."""
    requirements = [
        str(Requirement(item)) for item in _PROJECT_METADATA["dependencies"]
    ]
    for extra, items in _PROJECT_METADATA["optional-dependencies"].items():
        requirements.extend(
            str(Requirement(f'{item}; extra == "{extra}"')) for item in items
        )
    return tuple(requirements)


EXPECTED_REQUIRES_DIST = _expected_requires_dist()
EXPECTED_PROVIDES_EXTRA = tuple(_PROJECT_METADATA["optional-dependencies"])
PAX_MTIME_RE = re.compile(r"(?:0|[1-9][0-9]*)(?:\.[0-9]+)?")
PAX_MTIME_MAX_LENGTH = 32


def build_artifacts(outdir: Path) -> tuple[Path, Path]:
    """Build sdist + wheel (wheel from the sdist) into outdir."""
    result = subprocess.run(
        [sys.executable, "-m", "build", "--outdir", str(outdir)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        sys.stderr.write(result.stdout[-2000:] + "\n" + result.stderr[-2000:] + "\n")
        raise SystemExit("build failed")
    return next(outdir.glob("*.tar.gz")), next(outdir.glob("*.whl"))


def forbidden_in(
    names: list[str],
    *,
    directory_names: frozenset[str] = frozenset(),
    allowed_metadata_paths: frozenset[str] | None = None,
) -> list[str]:
    """Return archive paths forbidden by packaging or clean-room policy."""
    bad: set[str] = set()
    allowed_metadata_paths = (
        DEFAULT_ALLOWED_METADATA_PATHS
        if allowed_metadata_paths is None
        else allowed_metadata_paths
    )
    canonical_counts = Counter(_canonical_policy_rel(name) for name in names)
    for name in names:
        normalized = name.replace("\\", "/")
        is_directory = normalized.endswith("/") or name in directory_names
        canonical = _canonical_policy_rel(normalized)
        exact_parts = tuple(
            part for part in normalized.split("/") if part and part != "."
        )
        packaging_name = canonical if not is_directory else canonical.rstrip("/") + "/"
        packaging_match = any(pattern.search(packaging_name) for pattern in FORBIDDEN)
        metadata_match = any(
            part.casefold().rstrip(" .").endswith((".dist-info", ".egg-info"))
            and "/".join(exact_parts[: index + 1]) not in allowed_metadata_paths
            for index, part in enumerate(exact_parts)
        )
        clean_room_match = _classify_artifact_name(
            normalized,
            is_directory=is_directory,
        )
        if packaging_match or metadata_match or clean_room_match is not None:
            bad.add(name)
        if canonical and canonical_counts[canonical] > 1:
            bad.add(name)

    return sorted(bad)


def _zip_info_is_dir(info: zipfile.ZipInfo) -> bool:
    """Recognize ZIP directories by name or Unix mode metadata."""
    mode = (info.external_attr >> 16) & 0o177777
    return info.is_dir() or (mode != 0 and stat.S_ISDIR(mode))


def _zip_info_is_special(info: zipfile.ZipInfo) -> bool:
    """Return whether Unix metadata identifies a non-file, non-directory entry."""
    mode = (info.external_attr >> 16) & 0o177777
    file_type = stat.S_IFMT(mode)
    return file_type not in (0, stat.S_IFREG, stat.S_IFDIR)


def _archive_policy_rel(name: str, *, sdist_root: str | None = None) -> str:
    """Map only exact authenticated metadata to a legal-policy namespace."""
    rel = _archive_exact_rel(name, sdist_root=sdist_root)
    parts = tuple(part for part in rel.split("/") if part)
    basename = parts[-1] if parts else rel
    if (
        len(parts) >= 3
        and parts[0] == EXPECTED_DIST_INFO
        and parts[1] == "licenses"
        and basename in LICENSE_NAMING_OK
    ):
        return basename
    if rel == "PKG-INFO" or rel == f"{EXPECTED_DIST_INFO}/METADATA":
        return "README.md"
    return rel


def _archive_exact_rel(name: str, *, sdist_root: str | None = None) -> str:
    """Return a case-sensitive safe relative name for exact allowlists."""
    normalized = name.replace("\\", "/")
    parts = tuple(part for part in normalized.split("/") if part and part != ".")
    if len(parts) > 1 and sdist_root == EXPECTED_SDIST_ROOT and parts[0] == sdist_root:
        parts = parts[1:]
    return "/".join(parts)


def _text_payload_hits(
    name: str, payload: bytes, *, sdist_root: str | None = None
) -> list[tuple[int, str]]:
    """Return provenance hits for one safe-named archive payload."""
    rel = _archive_policy_rel(name, sdist_root=sdist_root)
    exact_rel = _archive_exact_rel(name, sdist_root=sdist_root)
    if exact_rel in ALLOWLIST or _is_binary_payload_name(name):
        return []
    try:
        text = payload.decode("utf-8", errors="strict")
    except UnicodeDecodeError:
        return [(0, "content-decode-error")]
    hits: list[tuple[int, str]] = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        hits.extend((lineno, label) for label in _line_hit_labels(rel, line))
    return hits


def _metadata_identity_hits(payload: bytes) -> list[str]:
    """Return identity mismatches in core metadata bytes."""
    try:
        metadata = BytesParser(policy=policy.default).parsebytes(payload)
    except (TypeError, ValueError):
        return ["metadata-parse-error"]

    if metadata.defects:
        return ["metadata-parse-defect"]

    expected = {
        "Metadata-Version": "2.4",
        "Name": PROJECT_NAME,
        "Version": PROJECT_VERSION,
        "License-Expression": PROJECT_LICENSE,
        "Requires-Python": _PROJECT_METADATA["requires-python"],
    }
    hits = [
        f"metadata-{field.casefold()}-mismatch"
        for field, value in expected.items()
        if [str(item) for item in metadata.get_all(field, [])] != [value]
    ]
    raw_requirements = [str(item) for item in metadata.get_all("Requires-Dist", [])]
    try:
        requirements = [str(Requirement(item)) for item in raw_requirements]
    except InvalidRequirement:
        hits.append("metadata-requires-dist-invalid")
    else:
        if Counter(requirements) != Counter(EXPECTED_REQUIRES_DIST):
            hits.append("metadata-requires-dist-mismatch")
    if Counter(str(item) for item in metadata.get_all("Provides-Extra", [])) != Counter(
        EXPECTED_PROVIDES_EXTRA
    ):
        hits.append("metadata-provides-extra-mismatch")
    if Counter(str(item) for item in metadata.get_all("License-File", [])) != Counter(
        _PROJECT_METADATA["license-files"]
    ):
        hits.append("metadata-license-file-mismatch")
    return hits


def _wheel_metadata_hits(payload: bytes) -> list[str]:
    """Validate the exact structural identity declared by WHEEL."""
    try:
        metadata = BytesParser(policy=policy.default).parsebytes(payload)
    except (TypeError, ValueError):
        return ["wheel-parse-error"]
    if metadata.defects:
        return ["wheel-parse-defect"]

    expected = {
        "Wheel-Version": ["1.0"],
        "Root-Is-Purelib": ["true"],
        "Tag": ["py3-none-any"],
    }
    return [
        f"wheel-{field.casefold()}-mismatch"
        for field, values in expected.items()
        if [str(item) for item in metadata.get_all(field, [])] != values
    ]


def _record_hits(payloads: dict[str, bytes], file_names: list[str]) -> list[str]:
    """Validate a wheel RECORD against every physical file member."""
    record_name = f"{EXPECTED_DIST_INFO}/RECORD"
    payload = payloads.get(record_name)
    if payload is None:
        return []
    try:
        text = payload.decode("utf-8", errors="strict")
        rows = list(csv.reader(io.StringIO(text, newline=""), strict=True))
    except (UnicodeDecodeError, csv.Error):
        return ["record-parse-error"]

    hits: list[str] = []
    records: dict[str, tuple[str, str]] = {}
    for row in rows:
        if len(row) != 3:
            hits.append("record-row-shape")
            continue
        name, digest, size = row
        if name in records:
            hits.append("record-duplicate-path")
            continue
        records[name] = (digest, size)

    if set(records) != set(file_names):
        hits.append("record-member-set-mismatch")
    for name in set(records) & set(file_names):
        digest, size = records[name]
        if name == record_name:
            if digest or size:
                hits.append("record-self-fields-not-empty")
            continue
        if name not in payloads:
            hits.append(f"record-unreadable-member:{name}")
            continue
        member_payload = payloads[name]
        expected_digest = "sha256=" + base64.urlsafe_b64encode(
            hashlib.sha256(member_payload).digest()
        ).rstrip(b"=").decode("ascii")
        if digest != expected_digest:
            hits.append(f"record-hash-mismatch:{name}")
        if size != str(len(member_payload)):
            hits.append(f"record-size-mismatch:{name}")
    return hits


def _requires_txt_hits(payload: bytes) -> list[str]:
    """Validate egg-info requirements against pyproject.toml."""
    try:
        lines = payload.decode("utf-8", errors="strict").splitlines()
    except UnicodeDecodeError:
        return ["requires-txt-decode-error"]

    groups: dict[str, list[str]] = {"": []}
    current = ""
    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            current = line[1:-1]
            if not current or current in groups:
                return ["requires-txt-section-error"]
            groups[current] = []
            continue
        try:
            groups[current].append(str(Requirement(line)))
        except InvalidRequirement:
            return ["requires-txt-invalid-requirement"]

    expected = {
        "": [str(Requirement(item)) for item in _PROJECT_METADATA["dependencies"]],
        **{
            extra: [str(Requirement(item)) for item in items]
            for extra, items in _PROJECT_METADATA["optional-dependencies"].items()
        },
    }
    if set(groups) != set(expected):
        return ["requires-txt-section-mismatch"]
    if any(
        Counter(groups[name]) != Counter(values) for name, values in expected.items()
    ):
        return ["requires-txt-requirement-mismatch"]
    return []


def _sources_txt_hits(payload: bytes) -> list[str]:
    """Validate the sdist source manifest as an exact unique file set."""
    try:
        lines = payload.decode("utf-8", errors="strict").splitlines()
    except UnicodeDecodeError:
        return ["sources-txt-decode-error"]
    expected = set(SDIST_REQUIRED) - {"PKG-INFO", "setup.cfg"}
    if len(lines) != len(set(lines)):
        return ["sources-txt-duplicate"]
    return [] if set(lines) == expected else ["sources-txt-set-mismatch"]


def _setup_cfg_hits(payload: bytes) -> list[str]:
    """Accept only setuptools' inert generated egg_info configuration."""
    try:
        text = payload.decode("utf-8", errors="strict")
        parser = configparser.ConfigParser(interpolation=None)
        parser.read_string(text)
    except (UnicodeDecodeError, configparser.Error):
        return ["setup-cfg-parse-error"]
    if parser.defaults() or parser.sections() != ["egg_info"]:
        return ["setup-cfg-section-mismatch"]
    values = dict(parser.items("egg_info", raw=True))
    return (
        []
        if values == {"tag_build": "", "tag_date": "0"}
        else ["setup-cfg-value-mismatch"]
    )


def _member_pax_headers_are_safe(member: tarfile.TarInfo) -> bool:
    """Allow only setuptools' bounded canonical fractional-mtime PAX header."""
    headers = member.pax_headers
    if not headers:
        return True
    if set(headers) != {"mtime"}:
        return False
    value = headers["mtime"]
    if (
        not isinstance(value, str)
        or not value
        or len(value) > PAX_MTIME_MAX_LENGTH
        or PAX_MTIME_RE.fullmatch(value) is None
    ):
        return False
    try:
        parsed = float(value)
        resolved = float(member.mtime)
    except (TypeError, ValueError, OverflowError):
        return False
    return math.isfinite(parsed) and math.isfinite(resolved) and parsed == resolved


def _required_payload_hits(
    payloads: dict[str, bytes],
    *,
    sdist: bool,
    wheel_file_names: list[str] | None = None,
) -> list[tuple[str, int, str]]:
    """Validate metadata identity and authenticated repository payloads."""
    hits: list[tuple[str, int, str]] = []
    metadata_names = (
        (
            f"{EXPECTED_SDIST_ROOT}/PKG-INFO",
            f"{EXPECTED_SDIST_ROOT}/{EXPECTED_EGG_INFO}/PKG-INFO",
        )
        if sdist
        else (f"{EXPECTED_DIST_INFO}/METADATA",)
    )
    for metadata_name in metadata_names:
        metadata_payload = payloads.get(metadata_name)
        if metadata_payload is not None:
            hits.extend(
                (metadata_name, 0, label)
                for label in _metadata_identity_hits(metadata_payload)
            )

    for filename in LEGAL_FILENAMES:
        member_name = (
            f"{EXPECTED_SDIST_ROOT}/{filename}"
            if sdist
            else f"{EXPECTED_DIST_INFO}/licenses/{filename}"
        )
        payload = payloads.get(member_name)
        if payload is not None and payload != (REPO_ROOT / filename).read_bytes():
            hits.append((member_name, 0, "legal-content-mismatch"))

    expected_payloads = SDIST_BYTE_REQUIRED if sdist else WHEEL_REQUIRED
    for filename in expected_payloads:
        member_name = f"{EXPECTED_SDIST_ROOT}/{filename}" if sdist else filename
        payload = payloads.get(member_name)
        if payload is not None and payload != (REPO_ROOT / filename).read_bytes():
            hits.append((member_name, 0, "repository-content-mismatch"))

    if not sdist:
        wheel_name = f"{EXPECTED_DIST_INFO}/WHEEL"
        wheel_payload = payloads.get(wheel_name)
        if wheel_payload is not None:
            hits.extend(
                (wheel_name, 0, label) for label in _wheel_metadata_hits(wheel_payload)
            )

        entry_points_name = f"{EXPECTED_DIST_INFO}/entry_points.txt"
        entry_points_payload = payloads.get(entry_points_name)
        if (
            entry_points_payload is not None
            and entry_points_payload != EXPECTED_ENTRY_POINTS
        ):
            hits.append((entry_points_name, 0, "entry-points-mismatch"))

        top_level_name = f"{EXPECTED_DIST_INFO}/top_level.txt"
        top_level_payload = payloads.get(top_level_name)
        if (
            top_level_payload is not None
            and top_level_payload != f"{PROJECT_NAME}\n".encode()
        ):
            hits.append((top_level_name, 0, "top-level-mismatch"))

        if wheel_file_names is not None:
            record_name = f"{EXPECTED_DIST_INFO}/RECORD"
            hits.extend(
                (record_name, 0, label)
                for label in _record_hits(payloads, wheel_file_names)
            )
    else:
        generated_payloads = {
            f"{EXPECTED_SDIST_ROOT}/{EXPECTED_EGG_INFO}/dependency_links.txt": b"\n",
            f"{EXPECTED_SDIST_ROOT}/{EXPECTED_EGG_INFO}/entry_points.txt": (
                EXPECTED_ENTRY_POINTS
            ),
            f"{EXPECTED_SDIST_ROOT}/{EXPECTED_EGG_INFO}/top_level.txt": (
                f"{PROJECT_NAME}\n".encode()
            ),
        }
        for member_name, expected_payload in generated_payloads.items():
            payload = payloads.get(member_name)
            if payload is not None and payload != expected_payload:
                hits.append((member_name, 0, "generated-metadata-mismatch"))
        generated_validators = {
            f"{EXPECTED_SDIST_ROOT}/{EXPECTED_EGG_INFO}/requires.txt": (
                _requires_txt_hits
            ),
            f"{EXPECTED_SDIST_ROOT}/{EXPECTED_EGG_INFO}/SOURCES.txt": (
                _sources_txt_hits
            ),
            f"{EXPECTED_SDIST_ROOT}/setup.cfg": _setup_cfg_hits,
        }
        for member_name, validator in generated_validators.items():
            payload = payloads.get(member_name)
            if payload is not None:
                hits.extend((member_name, 0, label) for label in validator(payload))
    return hits


def audit(sdist: Path | None, wheel: Path) -> int:
    problems = 0

    if sdist is not None:
        sdist_identity_hits = (
            []
            if sdist.name == EXPECTED_SDIST_FILENAME
            else [f"expected artifact filename {EXPECTED_SDIST_FILENAME}"]
        )
        with tarfile.open(sdist) as tf:
            members = tf.getmembers()
            sdist_names = [member.name for member in members]
            directory_names = frozenset(
                member.name for member in members if member.isdir()
            )
            canonical_roots = {
                canonical.split("/", maxsplit=1)[0]
                for member in members
                if (canonical := _canonical_policy_rel(member.name))
            }
            exact_roots = {
                exact_parts[0]
                for member in members
                if (
                    exact_parts := tuple(
                        part
                        for part in member.name.replace("\\", "/").split("/")
                        if part and part != "."
                    )
                )
            }
            valid_sdist_root = canonical_roots == {
                CANONICAL_EXPECTED_SDIST_ROOT
            } and exact_roots == {EXPECTED_SDIST_ROOT}
            sdist_root = EXPECTED_SDIST_ROOT if valid_sdist_root else None
            specials = sorted(
                member.name
                for member in members
                if not member.isfile() and not member.isdir()
            )
            archive_metadata_hits = []
            if tf.pax_headers:
                archive_metadata_hits.append(("<archive>", 0, "pax-global-header"))
            archive_metadata_hits.extend(
                (member.name, 0, "pax-member-header")
                for member in members
                if not _member_pax_headers_are_safe(member)
            )
            bad = forbidden_in(
                sdist_names,
                directory_names=directory_names,
                allowed_metadata_paths=frozenset(
                    {f"{EXPECTED_SDIST_ROOT}/{EXPECTED_EGG_INFO}"}
                ),
            )
            bad = sorted(
                set(bad)
                | {
                    name
                    for name in sdist_names
                    if name not in EXPECTED_SDIST_FILE_NAMES
                    and name not in EXPECTED_SDIST_DIRECTORY_NAMES
                }
            )
            missing_sdist = [
                required
                for required in SDIST_REQUIRED
                if f"{EXPECTED_SDIST_ROOT}/{required}" not in sdist_names
            ]
            if not valid_sdist_root:
                missing_sdist.append(f"archive root {EXPECTED_SDIST_ROOT}/")
            content_hits: list[tuple[str, int, str]] = []
            payloads: dict[str, bytes] = {}
            if not bad and not specials and not archive_metadata_hits:
                for member in members:
                    if not member.isfile():
                        continue
                    extracted = tf.extractfile(member)
                    if extracted is None:
                        content_hits.append((member.name, 0, "content-read-error"))
                        continue
                    try:
                        payload = extracted.read()
                    except OSError:
                        content_hits.append((member.name, 0, "content-read-error"))
                        continue
                    payloads[member.name] = payload
                    content_hits.extend(
                        (member.name, lineno, label)
                        for lineno, label in _text_payload_hits(
                            member.name, payload, sdist_root=sdist_root
                        )
                    )
                content_hits.extend(_required_payload_hits(payloads, sdist=True))
            content_hits.extend(archive_metadata_hits)
        for label in sdist_identity_hits:
            print(f"FORBIDDEN sdist identity: {label}")
        for name in bad:
            print(f"FORBIDDEN in sdist: {name}")
        for name in specials:
            print(f"FORBIDDEN special member in sdist: {name}")
        for name, lineno, label in content_hits:
            print(f"FORBIDDEN content in sdist: {name}:{lineno} [{label}]")
        for name in missing_sdist:
            print(f"MISSING from sdist: {name}")
        print(
            f"sdist audit: {sdist.name} — {len(sdist_names)} entries, "
            f"{len(bad)} forbidden names, {len(specials)} special members, "
            f"{len(content_hits)} content hits, {len(missing_sdist)} missing"
        )
        problems += (
            len(sdist_identity_hits)
            + len(bad)
            + len(specials)
            + len(content_hits)
            + len(missing_sdist)
        )

    wheel_identity_hits = (
        []
        if wheel.name == EXPECTED_WHEEL_FILENAME
        else [f"expected artifact filename {EXPECTED_WHEEL_FILENAME}"]
    )
    with zipfile.ZipFile(wheel) as zf:
        infos = zf.infolist()
        names = [info.filename for info in infos]
        directory_names = frozenset(
            info.filename for info in infos if _zip_info_is_dir(info)
        )
        specials = sorted(info.filename for info in infos if _zip_info_is_special(info))
        archive_metadata_hits = []
        if zf.comment:
            archive_metadata_hits.append(("<archive>", 0, "zip-global-comment"))
        for info in infos:
            if info.comment:
                archive_metadata_hits.append((info.filename, 0, "zip-member-comment"))
            if info.extra:
                archive_metadata_hits.append((info.filename, 0, "zip-member-extra"))
        bad = forbidden_in(
            names,
            directory_names=directory_names,
            allowed_metadata_paths=frozenset({EXPECTED_DIST_INFO}),
        )
        bad = sorted(
            set(bad) | {name for name in names if name not in EXPECTED_WHEEL_NAMES}
        )
        content_hits = []
        payloads = {}
        if not bad and not specials and not archive_metadata_hits:
            for info in infos:
                if _zip_info_is_dir(info):
                    continue
                try:
                    payload = zf.read(info)
                except (OSError, RuntimeError, zipfile.BadZipFile):
                    content_hits.append((info.filename, 0, "content-read-error"))
                    continue
                payloads[info.filename] = payload
                content_hits.extend(
                    (info.filename, lineno, label)
                    for lineno, label in _text_payload_hits(info.filename, payload)
                )
            file_names = [info.filename for info in infos if not _zip_info_is_dir(info)]
            content_hits.extend(
                _required_payload_hits(
                    payloads,
                    sdist=False,
                    wheel_file_names=file_names,
                )
            )
        content_hits.extend(archive_metadata_hits)
    missing = [
        required
        for required in (*WHEEL_REQUIRED, *WHEEL_METADATA_REQUIRED)
        if required not in names
    ]
    for label in wheel_identity_hits:
        print(f"FORBIDDEN wheel identity: {label}")
    for name in bad:
        print(f"FORBIDDEN in wheel: {name}")
    for name in specials:
        print(f"FORBIDDEN special member in wheel: {name}")
    for name, lineno, label in content_hits:
        print(f"FORBIDDEN content in wheel: {name}:{lineno} [{label}]")
    for name in missing:
        print(f"MISSING from wheel: {name}")
    print(
        f"wheel audit: {wheel.name} — {len(names)} entries, "
        f"{len(bad)} forbidden names, {len(specials)} special members, "
        f"{len(content_hits)} content hits, {len(missing)} missing"
    )
    problems += (
        len(wheel_identity_hits)
        + len(bad)
        + len(specials)
        + len(content_hits)
        + len(missing)
    )
    return 1 if problems else 0


def main(argv: list[str] | None = None) -> int:
    """Build and audit release artifacts, or audit one existing wheel."""
    parser = argparse.ArgumentParser(
        description="Build and audit rc8 distribution artifacts."
    )
    parser.add_argument(
        "wheel",
        nargs="?",
        type=Path,
        help="existing wheel to audit instead of building fresh artifacts",
    )
    args = parser.parse_args(argv)
    if args.wheel is not None:
        return audit(None, args.wheel)
    with tempfile.TemporaryDirectory() as tmp:
        sdist, wheel = build_artifacts(Path(tmp))
        return audit(sdist, wheel)


if __name__ == "__main__":
    sys.exit(main())
