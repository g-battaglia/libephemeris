"""Regression tests for project and wheel metadata configuration."""

from __future__ import annotations

import tomllib
from pathlib import Path

from scripts.check_wheel_contents import forbidden_in


_PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _pyproject() -> dict:
    """Load the authoritative package configuration."""
    return tomllib.loads((_PROJECT_ROOT / "pyproject.toml").read_text())


def _has_requirement(requirements: list[str], package: str) -> bool:
    """Return whether a requirement list contains the named package."""
    normalized = package.lower()
    return any(
        requirement.lower().startswith(normalized) for requirement in requirements
    )


def test_astroquery_is_optional_catalog_tooling() -> None:
    """Core installs stay lean while catalog and dev installs remain complete."""
    config = _pyproject()["project"]
    extras = config["optional-dependencies"]

    assert not _has_requirement(config["dependencies"], "astroquery")
    assert _has_requirement(extras["stars"], "astroquery")
    assert _has_requirement(extras["dev"], "astroquery")


def test_bundled_data_namespaces_are_explicit_packages() -> None:
    """Setuptools must not rely on ambiguous implicit data-package handling."""
    packages = set(_pyproject()["tool"]["setuptools"]["packages"])

    assert "libephemeris.data" in packages
    assert "libephemeris.data.leb2" in packages


def test_archive_gate_rejects_clean_room_artifact_names_case_insensitively() -> None:
    """Built archives must enforce the same filename gate as the worktree."""
    prohibited = [
        "libephemeris/data/seasnam.txt",
        "libephemeris/data/sweph.c",
        "libephemeris/vendor/sweprivate.h",
        "libephemeris/Data/REFERENCE/planet.bin",
        "libephemeris/vendor/PySwissEph-2.10/module.py",
    ]

    assert forbidden_in(prohibited) == sorted(prohibited)


def test_archive_gate_allows_project_owned_runtime_data() -> None:
    """The clean-room gate must not reject LibEphemeris runtime assets."""
    allowed = [
        "libephemeris/data/fictitious_orbits.csv",
        "libephemeris/data/leb2/base_core.leb2",
        "libephemeris/vendor/spktype21.py",
    ]

    assert forbidden_in(allowed) == []
