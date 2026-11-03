"""Product source and release archives must not require the private harness."""
from pathlib import Path
import tomllib

ROOT = Path(__file__).resolve().parents[2]


def test_scientific_harness_is_external():
    assert not (ROOT / "validation").exists()
    assert not list((ROOT / "tools").glob("*.py")), "maintenance scripts belong in purpose-specific subdirectories"
    # Git's ignore policy matters in a checkout, not in an extracted sdist
    # (which intentionally does not need .git or .gitignore).
    if (ROOT / ".git").exists():
        assert "/libephemeris-validation-python/" in (ROOT / ".gitignore").read_text()
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert "libephemeris-validation-python/**" in project["tool"]["maturin"]["exclude"]
    assert not any(
        item["path"].startswith("validation/")
        for item in project["tool"]["maturin"]["include"]
    )


def test_project_license_copies_are_identical():
    expected = (ROOT / "COPYING").read_bytes()
    for directory in ("core", "core/data", "tools/leb_generation"):
        assert (ROOT / directory / "COPYING").read_bytes() == expected
    notices = (ROOT / "THIRD_PARTY_NOTICES.md").read_bytes()
    for directory in ("core", "core/data"):
        assert (ROOT / directory / "THIRD_PARTY_NOTICES.md").read_bytes() == notices


def test_distribution_versions_are_deliberately_independent():
    core = tomllib.loads((ROOT / "core/Cargo.toml").read_text())["package"]
    data = tomllib.loads((ROOT / "core/data/Cargo.toml").read_text())["package"]
    binding = tomllib.loads((ROOT / "python/Cargo.toml").read_text())
    python = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    tools = tomllib.loads((ROOT / "tools/leb_generation/pyproject.toml").read_text())["project"]
    assert core["version"] == data["version"] == "1.0.0-beta.1"
    assert binding["package"]["version"] == "4.0.0-beta.1"
    assert binding["dependencies"]["libephemeris-core"]["version"] == "=1.0.0-beta.1"
    assert "version" not in python and python["dynamic"] == ["version"]
    assert tools["version"] == "4.0.0b1"
    assert all(item["license"] == "AGPL-3.0" for item in (core, data, binding["package"], python, tools))
