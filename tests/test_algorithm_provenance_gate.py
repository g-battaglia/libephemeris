# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Tests for the complete algorithm/data provenance inventory gate."""

from __future__ import annotations

from pathlib import Path

from scripts import check_algorithm_provenance as gate


def test_repository_algorithm_provenance_gate_is_green() -> None:
    """The checked-in registry and source tree must satisfy the live gate."""
    assert gate.main() == 0


def test_local_source_locator_requires_a_live_markdown_anchor(
    tmp_path: Path, monkeypatch
) -> None:
    """A source link to an existing file must not hide a stale section name."""
    document = tmp_path / "docs" / "method.md"
    document.parent.mkdir(parents=True)
    document.write_text("# Method\n\n## Public derivation\n", encoding="utf-8")
    monkeypatch.setattr(gate, "REPO_ROOT", tmp_path)

    problems: list[str] = []
    gate._check_local_source_locator(
        "docs/method.md#missing-section", label="source 'TEST'", problems=problems
    )

    assert problems == [
        "source 'TEST': Markdown anchor does not resolve: "
        "docs/method.md#missing-section"
    ]


def test_project_code_placeholder_scan_exempts_identified_vendor_code(
    tmp_path: Path, monkeypatch
) -> None:
    """Project prose fails on empty attribution while vendored prose is preserved."""
    project_path = "libephemeris/model.py"
    vendor_path = "libephemeris/vendor/upstream.py"
    for relative in (project_path, vendor_path):
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# Uses the standard formula.\n", encoding="utf-8")
    monkeypatch.setattr(gate, "REPO_ROOT", tmp_path)

    problems: list[str] = []
    gate._check_project_code_placeholders(
        {project_path, vendor_path}, {vendor_path}, problems
    )

    assert problems == [
        "libephemeris/model.py:1: vague source-attribution placeholder "
        "'standard formula'"
    ]


def test_every_project_top_level_symbol_requires_a_docstring(
    tmp_path: Path, monkeypatch
) -> None:
    """Concrete private and public definitions are checked; overloads are not."""
    module = tmp_path / "libephemeris" / "model.py"
    module.parent.mkdir(parents=True)
    module.write_text(
        "from typing import overload\n\n"
        "@overload\n"
        "def documented(value: int) -> int: ...\n\n"
        "def _private() -> None:\n"
        '    """Document the private implementation boundary too."""\n'
        "    pass\n\n"
        "def undocumented() -> None:\n"
        "    pass\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(gate, "REPO_ROOT", tmp_path)

    problems: list[str] = []
    gate._check_top_level_docstrings("libephemeris/model.py", problems)

    assert problems == [
        "libephemeris/model.py:10: top-level FunctionDef 'undocumented' "
        "lacks a docstring"
    ]


def test_private_top_level_helper_cannot_hide_undocumented_algorithm(
    tmp_path: Path, monkeypatch
) -> None:
    """An underscore does not exempt a module-level numerical helper."""
    module = tmp_path / "libephemeris" / "model.py"
    module.parent.mkdir(parents=True)
    module.write_text("def _undocumented() -> None:\n    pass\n", encoding="utf-8")
    monkeypatch.setattr(gate, "REPO_ROOT", tmp_path)

    problems: list[str] = []
    gate._check_top_level_docstrings("libephemeris/model.py", problems)

    assert problems == [
        "libephemeris/model.py:1: top-level FunctionDef '_undocumented' "
        "lacks a docstring"
    ]
