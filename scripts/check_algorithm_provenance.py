#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Verify complete algorithm, data, and tooling provenance coverage.

The authoritative inventory is
``docs/methodology/provenance-registry.toml``. This gate proves structural
claims that ordinary prose review cannot prove reliably:

* every Python module under ``libephemeris/`` is classified exactly once;
* every Python program under ``scripts/`` is classified exactly once;
* every data asset selected by ``tool.setuptools.package-data`` is recorded;
* every cited source ID and local documentation path resolves;
* generated artifacts identify their generators;
* every module has exactly one substantive ``Provenance:`` section;
* generated modules name their authoritative generator; and
* vendored modules identify their upstream and licensing boundary.
* every project-owned top-level function and class has a docstring.

The check intentionally does not decide whether a scientific citation is
correct or whether an implementation faithfully evaluates its equations.
Those require scientific review and targeted tests. It makes omissions,
duplicate ownership, stale paths, vague placeholders, and undocumented new
files fail closed.

Provenance:
    Project-authored repository-integrity tooling. It parses only project TOML
    and Python syntax and contains no astronomical model or coefficient.
"""

from __future__ import annotations

import ast
import re
import sys
import tomllib
from collections import Counter
from pathlib import Path, PurePosixPath
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
REGISTRY_PATH = REPO_ROOT / "docs/methodology/provenance-registry.toml"
PYPROJECT_PATH = REPO_ROOT / "pyproject.toml"
CHECKER_RELATIVE_PATH = "scripts/check_algorithm_provenance.py"

ALLOWED_CATEGORIES = frozenset(
    {
        "development-tool",
        "generated-scientific-data",
        "infrastructure",
        "numerical-method",
        "public-api",
        "scientific-data",
        "scientific-model",
        "scientific-pipeline",
        "vendored",
        "verification-tool",
    }
)
ALLOWED_DOCSTRING_LEVELS = frozenset({"detailed"})
GENERATED_CATEGORIES = frozenset({"generated-scientific-data"})
SOURCE_REQUIRED_CATEGORIES = frozenset(
    {
        "generated-scientific-data",
        "numerical-method",
        "scientific-data",
        "scientific-model",
        "scientific-pipeline",
        "vendored",
    }
)

# Placeholders that can make an entry look documented while supplying no
# independently recoverable source. This applies to registry prose only, not
# to historical narrative elsewhere in the repository.
VAGUE_SOURCE_RE = re.compile(
    r"\b(?:reference documentation|standard formula|various sources|"
    r"source unknown|unknown source|citation needed|tbd|todo)\b",
    re.IGNORECASE,
)


def _load_toml(path: Path) -> dict[str, Any]:
    """Load one UTF-8-compatible TOML document from ``path``."""
    with path.open("rb") as stream:
        return tomllib.load(stream)


def _python_inventory(root: Path) -> set[str]:
    """Return canonical repository-relative Python paths below ``root``."""
    return {
        path.relative_to(REPO_ROOT).as_posix()
        for path in root.rglob("*.py")
        if "__pycache__" not in path.parts
    }


def _shipped_data_inventory() -> set[str]:
    """Return data assets explicitly selected for the wheel.

    ``tool.setuptools.package-data.libephemeris`` is the release boundary. A
    local ignored download is therefore neither silently accepted nor reported
    as a shipped asset.
    """
    pyproject = _load_toml(PYPROJECT_PATH)
    package_data = pyproject["tool"]["setuptools"]["package-data"]["libephemeris"]
    return {
        f"libephemeris/{entry}"
        for entry in package_data
        if PurePosixPath(entry).parts[:1] == ("data",)
    }


def _valid_repo_path(raw: object, *, label: str, problems: list[str]) -> str | None:
    """Validate and normalize a registry path without allowing traversal."""
    if not isinstance(raw, str) or not raw.strip():
        problems.append(f"{label}: path must be a non-empty string")
        return None
    path = PurePosixPath(raw)
    if path.is_absolute() or ".." in path.parts or raw != path.as_posix():
        problems.append(f"{label}: unsafe or non-canonical path {raw!r}")
        return None
    return raw


def _check_local_document(
    raw: object, *, label: str, problems: list[str]
) -> str | None:
    """Validate a registry documentation path and require an existing file."""
    path = _valid_repo_path(raw, label=label, problems=problems)
    if path is not None and not (REPO_ROOT / path).is_file():
        problems.append(f"{label}: documentation file does not exist: {path}")
    return path


def _markdown_heading_ids(text: str) -> set[str]:
    """Return anchors for the ordinary ASCII headings used in project docs.

    This intentionally implements only the repository's simple GitHub-style
    headings. It is sufficient to make a stale ``docs/file.md#section`` source
    locator fail closed without claiming to emulate every Markdown renderer.
    """
    anchors: set[str] = set()
    for line in text.splitlines():
        match = re.match(r"^#{1,6}\s+(.+?)\s*#*\s*$", line)
        if match is None:
            continue
        heading = match.group(1).strip().lower()
        heading = re.sub(r"[^a-z0-9 _-]", "", heading)
        # GitHub-style slugs remove punctuation such as backticks and dots,
        # replace whitespace with hyphens, and preserve underscores in code
        # identifiers (for example ``arabic_parts.py``).
        anchor = re.sub(r" +", "-", heading).strip("-")
        if anchor:
            anchors.add(anchor)
    return anchors


def _check_local_source_locator(raw: str, *, label: str, problems: list[str]) -> None:
    """Validate a local source file and its optional Markdown anchor."""
    local_path, separator, anchor = raw.partition("#")
    full_path = REPO_ROOT / local_path
    if not full_path.is_file():
        problems.append(f"{label}: local source record does not exist: {raw}")
        return
    if not separator:
        return
    if not anchor:
        problems.append(f"{label}: local source locator has an empty anchor: {raw}")
        return
    try:
        heading_ids = _markdown_heading_ids(full_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError) as exc:
        problems.append(f"{label}: cannot read local source locator {raw}: {exc}")
        return
    if anchor not in heading_ids:
        problems.append(f"{label}: Markdown anchor does not resolve: {raw}")


def _module_docstring(path: str, problems: list[str]) -> str:
    """Parse and return a module docstring, recording syntax/read failures."""
    full_path = REPO_ROOT / path
    try:
        tree = ast.parse(full_path.read_text(encoding="utf-8"), filename=path)
    except (OSError, SyntaxError, UnicodeError) as exc:
        problems.append(f"{path}: cannot parse module for provenance docstring: {exc}")
        return ""
    docstring = ast.get_docstring(tree, clean=False) or ""
    if not docstring.strip():
        problems.append(f"{path}: missing module docstring")
    return docstring


def _provenance_section(path: str, docstring: str, problems: list[str]) -> str:
    """Extract and validate the module-level Google-style provenance block.

    Requiring one exact heading prevents an incidental use of the word
    "provenance" elsewhere in a long docstring from satisfying the gate. The
    body must be indented so its boundary remains mechanically reviewable.
    """
    lines = docstring.splitlines()
    markers = [
        index for index, line in enumerate(lines) if line.strip() == "Provenance:"
    ]
    if len(markers) != 1:
        problems.append(
            f"{path}: module docstring must contain exactly one "
            f"'Provenance:' heading (found {len(markers)})"
        )
        return ""

    body: list[str] = []
    started = False
    for line in lines[markers[0] + 1 :]:
        if not line.strip():
            if started:
                body.append("")
            continue
        if line == line.lstrip():
            break
        started = True
        body.append(line.strip())

    section = " ".join(body).strip()
    if len(section) < 120:
        problems.append(
            f"{path}: Provenance section is too short to explain its source "
            f"boundary ({len(section)} characters; minimum 120)"
        )
    if VAGUE_SOURCE_RE.search(section):
        problems.append(f"{path}: vague placeholder in Provenance section")
    return section


def _is_overload(node: ast.FunctionDef | ast.AsyncFunctionDef) -> bool:
    """Return whether a top-level definition is a typing overload stub."""
    for decorator in node.decorator_list:
        if isinstance(decorator, ast.Name) and decorator.id == "overload":
            return True
        if isinstance(decorator, ast.Attribute) and decorator.attr == "overload":
            return True
    return False


def _check_top_level_docstrings(path: str, problems: list[str]) -> None:
    """Require documentation for every project-owned top-level definition.

    Private top-level helpers often contain the actual numerical branch or
    transformation behind a small public wrapper, so an underscore is not an
    exemption. Nested callbacks and class methods remain governed by their
    enclosing documented algorithm; requiring every tiny closure would reward
    boilerplate rather than improve the source boundary.
    """
    full_path = REPO_ROOT / path
    try:
        tree = ast.parse(full_path.read_text(encoding="utf-8"), filename=path)
    except (OSError, SyntaxError, UnicodeError):
        # ``_module_docstring`` reports the same parse failure with context.
        return
    for node in tree.body:
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and _is_overload(
            node
        ):
            continue
        if not (ast.get_docstring(node, clean=False) or "").strip():
            problems.append(
                f"{path}:{node.lineno}: top-level {type(node).__name__} "
                f"{node.name!r} lacks a docstring"
            )


def _check_project_code_placeholders(
    paths: set[str], vendored_paths: set[str], problems: list[str]
) -> None:
    """Reject non-auditable source placeholders in project-owned Python."""
    for path in sorted(paths - vendored_paths - {CHECKER_RELATIVE_PATH}):
        try:
            source_text = (REPO_ROOT / path).read_text(encoding="utf-8")
        except (OSError, UnicodeError) as exc:
            problems.append(f"{path}: cannot scan source-attribution prose: {exc}")
            continue
        for line_number, line in enumerate(source_text.splitlines(), start=1):
            match = VAGUE_SOURCE_RE.search(line)
            if match is not None:
                problems.append(
                    f"{path}:{line_number}: vague source-attribution placeholder "
                    f"{match.group(0)!r}"
                )


def _check_sources(registry: dict[str, Any], problems: list[str]) -> set[str]:
    """Validate source records and return the usable source-ID inventory."""
    raw_sources = registry.get("sources")
    if not isinstance(raw_sources, dict) or not raw_sources:
        problems.append("registry: [sources] must be a non-empty table")
        return set()

    source_ids: set[str] = set()
    for source_id, record in sorted(raw_sources.items()):
        label = f"source {source_id!r}"
        if not isinstance(source_id, str) or not source_id.strip():
            problems.append(f"{label}: ID must be a non-empty string")
            continue
        source_ids.add(source_id)
        if not isinstance(record, dict):
            problems.append(f"{label}: record must be a table")
            continue
        for field in ("kind", "citation", "url"):
            value = record.get(field)
            if not isinstance(value, str) or not value.strip():
                problems.append(f"{label}: missing non-empty {field!r}")
            elif VAGUE_SOURCE_RE.search(value):
                problems.append(f"{label}: vague placeholder in {field!r}: {value!r}")
        url = record.get("url")
        if isinstance(url, str) and not url.startswith(("https://", "docs/")):
            problems.append(
                f"{label}: URL must be HTTPS or a local docs/ path: {url!r}"
            )
        if isinstance(url, str) and url.startswith("docs/"):
            _check_local_source_locator(url, label=label, problems=problems)
    return source_ids


def _check_scope_component_map(
    registry: dict[str, Any], scope_document: object, problems: list[str]
) -> None:
    """Keep the human component tables exactly aligned with the TOML registry."""
    if not isinstance(scope_document, str):
        return
    path = REPO_ROOT / scope_document
    if not path.is_file():
        return
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        problems.append(f"registry scope_document: cannot read component map: {exc}")
        return
    start_marker = "## Complete runtime-module map"
    end_marker = "## Source register"
    if start_marker not in text or end_marker not in text:
        problems.append(
            "registry scope_document: missing complete component-map boundaries"
        )
        return
    section = text.split(start_marker, 1)[1].split(end_marker, 1)[0]
    documented = re.findall(r"^\| `([^`]+)` \|", section, flags=re.MULTILINE)
    raw_components = registry.get("component", [])
    registered = [
        record.get("id")
        for record in raw_components
        if isinstance(record, dict) and isinstance(record.get("id"), str)
    ]
    for component_id, count in sorted(Counter(documented).items()):
        if count > 1:
            problems.append(
                f"registry scope_document: component {component_id!r} appears "
                f"{count} times in the human map"
            )
    documented_set = set(documented)
    registered_set = set(registered)
    for component_id in sorted(registered_set - documented_set):
        problems.append(
            f"registry scope_document: component missing from human map: "
            f"{component_id!r}"
        )
    for component_id in sorted(documented_set - registered_set):
        problems.append(
            f"registry scope_document: stale component in human map: {component_id!r}"
        )


def _check_component_records(
    registry: dict[str, Any], source_ids: set[str], problems: list[str]
) -> tuple[set[str], dict[str, str], set[str]]:
    """Validate component records and return paths, levels, and vendor paths."""
    raw_components = registry.get("component")
    if not isinstance(raw_components, list) or not raw_components:
        problems.append("registry: at least one [[component]] record is required")
        return set(), {}, set()

    ids: list[str] = []
    covered_paths: list[str] = []
    docstring_levels: dict[str, str] = {}
    vendored_paths: set[str] = set()

    for index, record in enumerate(raw_components, start=1):
        label = f"component #{index}"
        if not isinstance(record, dict):
            problems.append(f"{label}: record must be a table")
            continue

        component_id = record.get("id")
        if not isinstance(component_id, str) or not component_id.strip():
            problems.append(f"{label}: missing non-empty id")
            component_id = f"<invalid-{index}>"
        else:
            ids.append(component_id)
            label = f"component {component_id!r}"

        category = record.get("category")
        if category not in ALLOWED_CATEGORIES:
            problems.append(f"{label}: invalid category {category!r}")

        scientific_content = record.get("scientific_content")
        if not isinstance(scientific_content, bool):
            problems.append(f"{label}: scientific_content must be true or false")

        origin = record.get("origin")
        if not isinstance(origin, str) or len(origin.strip()) < 40:
            problems.append(f"{label}: origin must explain the boundary in detail")
        elif VAGUE_SOURCE_RE.search(origin):
            problems.append(f"{label}: vague placeholder in origin: {origin!r}")

        raw_sources = record.get("sources")
        if not isinstance(raw_sources, list) or not all(
            isinstance(source, str) for source in raw_sources
        ):
            problems.append(f"{label}: sources must be an array of source IDs")
            raw_sources = []
        if scientific_content and not raw_sources:
            problems.append(f"{label}: scientific content requires at least one source")
        if category in SOURCE_REQUIRED_CATEGORIES and not raw_sources:
            problems.append(f"{label}: category {category!r} requires source records")
        for source in raw_sources:
            if source not in source_ids:
                problems.append(f"{label}: unknown source ID {source!r}")

        raw_docs = record.get("documentation")
        if not isinstance(raw_docs, list) or not raw_docs:
            problems.append(
                f"{label}: at least one local documentation path is required"
            )
        else:
            for doc_index, document in enumerate(raw_docs, start=1):
                _check_local_document(
                    document,
                    label=f"{label} documentation #{doc_index}",
                    problems=problems,
                )

        level = record.get("docstring")
        if level not in ALLOWED_DOCSTRING_LEVELS:
            problems.append(f"{label}: invalid docstring level {level!r}")
            level = "registry"

        raw_paths = record.get("paths")
        if not isinstance(raw_paths, list) or not raw_paths:
            problems.append(f"{label}: paths must be a non-empty array")
            continue
        for path_index, raw_path in enumerate(raw_paths, start=1):
            path = _valid_repo_path(
                raw_path, label=f"{label} path #{path_index}", problems=problems
            )
            if path is None:
                continue
            if not path.endswith(".py"):
                problems.append(f"{label}: component path is not Python: {path}")
            if not (REPO_ROOT / path).is_file():
                problems.append(f"{label}: component path does not exist: {path}")
            covered_paths.append(path)
            docstring_levels[path] = level
            if category == "vendored":
                vendored_paths.add(path)

        raw_generators = record.get("generator", [])
        if not isinstance(raw_generators, list) or not all(
            isinstance(path, str) for path in raw_generators
        ):
            problems.append(f"{label}: generator must be an array of script paths")
            raw_generators = []
        if category in GENERATED_CATEGORIES and not raw_generators:
            problems.append(f"{label}: generated category requires generator paths")
        for generator in raw_generators:
            path = _valid_repo_path(
                generator, label=f"{label} generator", problems=problems
            )
            if path is not None and (
                not path.startswith("scripts/") or not (REPO_ROOT / path).is_file()
            ):
                problems.append(f"{label}: generator is not an existing script: {path}")

    duplicate_ids = sorted(item for item, count in Counter(ids).items() if count > 1)
    for component_id in duplicate_ids:
        problems.append(f"registry: duplicate component id {component_id!r}")

    duplicate_paths = sorted(
        item for item, count in Counter(covered_paths).items() if count > 1
    )
    for path in duplicate_paths:
        problems.append(f"registry: Python path classified more than once: {path}")

    return set(covered_paths), docstring_levels, vendored_paths


def _check_asset_records(
    registry: dict[str, Any], source_ids: set[str], problems: list[str]
) -> set[str]:
    """Validate packaged scientific-asset lineage and return registered paths."""
    records = registry.get("asset")
    if not isinstance(records, list):
        problems.append("registry: [[asset]] records are required")
        return set()

    paths: list[str] = []
    for index, record in enumerate(records, start=1):
        label = f"asset #{index}"
        if not isinstance(record, dict):
            problems.append(f"{label}: record must be a table")
            continue
        path = _valid_repo_path(record.get("path"), label=label, problems=problems)
        if path is None:
            continue
        label = f"asset {path!r}"
        paths.append(path)
        if not (REPO_ROOT / path).is_file():
            problems.append(f"{label}: file does not exist")
        category = record.get("category")
        if category not in {"scientific-data", "generated-scientific-data"}:
            problems.append(f"{label}: invalid scientific asset category {category!r}")
        origin = record.get("origin")
        if not isinstance(origin, str) or len(origin.strip()) < 40:
            problems.append(f"{label}: origin must explain the asset lineage")
        raw_sources = record.get("sources")
        if not isinstance(raw_sources, list) or not raw_sources:
            problems.append(f"{label}: at least one source ID is required")
            raw_sources = []
        for source in raw_sources:
            if source not in source_ids:
                problems.append(f"{label}: unknown source ID {source!r}")
        raw_docs = record.get("documentation")
        if not isinstance(raw_docs, list) or not raw_docs:
            problems.append(f"{label}: at least one documentation path is required")
        else:
            for doc_index, document in enumerate(raw_docs, start=1):
                _check_local_document(
                    document,
                    label=f"{label} documentation #{doc_index}",
                    problems=problems,
                )
        generators = record.get("generator", [])
        if not isinstance(generators, list) or not all(
            isinstance(generator, str) for generator in generators
        ):
            problems.append(f"{label}: generator must be an array of script paths")
            generators = []
        if category == "generated-scientific-data" and not generators:
            problems.append(f"{label}: generated asset requires generator paths")
        for generator in generators:
            if not isinstance(generator, str) or not (
                generator.startswith("scripts/") and (REPO_ROOT / generator).is_file()
            ):
                problems.append(f"{label}: invalid generator {generator!r}")

    for path, count in Counter(paths).items():
        if count > 1:
            problems.append(f"registry: data asset classified more than once: {path}")
    return set(paths)


def main() -> int:
    """Run the fail-closed registry/source-tree audit and return shell status."""
    problems: list[str] = []
    try:
        registry = _load_toml(REGISTRY_PATH)
    except (OSError, tomllib.TOMLDecodeError) as exc:
        print(f"cannot load {REGISTRY_PATH.relative_to(REPO_ROOT)}: {exc}")
        return 1

    if registry.get("schema_version") != 1:
        problems.append("registry: schema_version must be 1")
    scope_document = registry.get("scope_document")
    _check_local_document(
        scope_document, label="registry scope_document", problems=problems
    )

    source_ids = _check_sources(registry, problems)
    _check_scope_component_map(registry, scope_document, problems)
    registered_python, docstring_levels, vendored_paths = _check_component_records(
        registry, source_ids, problems
    )
    registered_assets = _check_asset_records(registry, source_ids, problems)

    actual_python = _python_inventory(REPO_ROOT / "libephemeris")
    actual_python |= _python_inventory(REPO_ROOT / "scripts")
    for path in sorted(actual_python - registered_python):
        problems.append(f"inventory: unclassified Python file: {path}")
    for path in sorted(registered_python - actual_python):
        problems.append(f"inventory: registered Python file is outside scope: {path}")

    shipped_assets = _shipped_data_inventory()
    for path in sorted(shipped_assets - registered_assets):
        problems.append(f"inventory: unclassified shipped data asset: {path}")
    for path in sorted(registered_assets - shipped_assets):
        problems.append(
            f"inventory: registered asset is not selected for packaging: {path}"
        )

    component_by_path = {
        path: record
        for record in registry.get("component", [])
        if isinstance(record, dict)
        for path in record.get("paths", [])
        if isinstance(path, str)
    }

    for path in sorted(actual_python):
        docstring = _module_docstring(path, problems)
        if path not in vendored_paths:
            _check_top_level_docstrings(path, problems)
        level = docstring_levels.get(path)
        if level != "detailed":
            continue
        section = _provenance_section(path, docstring, problems)
        record = component_by_path.get(path, {})
        category = record.get("category")

        if category == "generated-scientific-data":
            generators = record.get("generator", [])
            valid_generators = (
                generators
                if isinstance(generators, list)
                and all(isinstance(item, str) for item in generators)
                else []
            )
            if not any(
                PurePosixPath(item).name in docstring for item in valid_generators
            ):
                problems.append(
                    f"{path}: generated module docstring does not name an "
                    "authoritative registered generator"
                )

        if category == "vendored":
            section_lower = section.lower()
            if "upstream" not in section_lower:
                problems.append(f"{path}: vendored provenance must identify upstream")
            if not any(
                marker in section_lower for marker in ("licence", "license", "spdx")
            ):
                problems.append(
                    f"{path}: vendored provenance must identify its licence boundary"
                )

    _check_project_code_placeholders(actual_python, vendored_paths, problems)

    if problems:
        for problem in problems:
            print(f"ERROR: {problem}")
        print(
            f"FAIL: {len(problems)} problem(s); "
            f"{len(actual_python)} Python files and {len(shipped_assets)} assets scanned"
        )
        return 1

    print(
        "OK: algorithm provenance registry covers "
        f"{len(actual_python)} Python files, {len(shipped_assets)} shipped data assets, "
        f"and {len(source_ids)} source records"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
