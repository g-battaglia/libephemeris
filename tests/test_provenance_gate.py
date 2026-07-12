"""Tests for the clean-room physical-worktree provenance gate."""

from pathlib import Path

from scripts import check_provenance


def _hits(root: Path) -> set[tuple[str, str]]:
    return {
        (path.relative_to(root).as_posix(), label)
        for path, label, _ in check_provenance._physical_foreign_artifacts(root)
    }


def test_physical_gate_catches_ignored_reference_directory(tmp_path: Path) -> None:
    reference_dir = tmp_path / "data" / "reference"
    reference_dir.mkdir(parents=True)

    assert _hits(tmp_path) == {("data/reference", "foreign-data-dir")}


def test_physical_gate_catches_mixed_case_reference_directory(tmp_path: Path) -> None:
    reference_dir = tmp_path / "Data" / "REFERENCE"
    reference_dir.mkdir(parents=True)

    assert _hits(tmp_path) == {("Data/REFERENCE", "foreign-data-dir")}


def test_physical_gate_catches_reference_files_by_name(tmp_path: Path) -> None:
    artifact_dir = tmp_path / "ignored" / "nested"
    artifact_dir.mkdir(parents=True)
    for name in ("planet.se1", "sefstars.json", "sweph.c", "sweprivate.h"):
        (artifact_dir / name).touch()

    assert _hits(tmp_path) == {
        ("ignored/nested/planet.se1", "foreign-data-file"),
        ("ignored/nested/sefstars.json", "foreign-data-file"),
        ("ignored/nested/sweph.c", "foreign-source-file"),
        ("ignored/nested/sweprivate.h", "foreign-source-file"),
    }


def test_physical_gate_catches_reference_source_directory(tmp_path: Path) -> None:
    source_dir = tmp_path / "vendor" / "Swiss-Ephemeris-source"
    source_dir.mkdir(parents=True)

    assert _hits(tmp_path) == {("vendor/Swiss-Ephemeris-source", "foreign-source-dir")}


def test_physical_gate_prunes_metadata_environments_and_caches(tmp_path: Path) -> None:
    for dirname in (".git", ".venv", "__pycache__", ".pytest_cache"):
        directory = tmp_path / dirname
        directory.mkdir()
        (directory / "planet.se1").touch()

    assert _hits(tmp_path) == set()
