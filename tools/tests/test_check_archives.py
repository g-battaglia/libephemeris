"""Release archive checks reject caches/process files without extracting them."""

import importlib.util
import io
from pathlib import Path
import tarfile
import zipfile

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("check_archives", ROOT / "tools/release/check_archives.py")
CHECKER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECKER)
METADATA = (
    b"Metadata-Version: 2.4\nName: libephemeris\nVersion: 4.0.0b1\n"
    b"License-Expression: AGPL-3.0\nDescription-Content-Type: text/markdown\n\n"
    b"# LibEphemeris\nNative engine.\n"
)


def test_source_path_rejections():
    names = [
        "package/AGENTS.md", "package/docs/_build/page.html", "package/.venv/bin/python",
        "package/tools/application.egg-info/PKG-INFO", "../escape", "/absolute",
        "package/native/rebound/file.c", "package/libephemeris-validation-python/file",
        "package/validation/v4/compare.py",
    ]
    assert len(CHECKER.path_errors(names)) == len(names)
    assert CHECKER.path_errors(["package/src/lib.rs", "package/licenses/index.json"]) == []


def test_metadata_requires_readme_and_unambiguous_license():
    assert CHECKER.metadata_errors(METADATA) == []
    with_parameters = METADATA.replace(b"text/markdown", b"text/markdown; charset=UTF-8; variant=GFM")
    assert CHECKER.metadata_errors(with_parameters) == []
    assert CHECKER.metadata_errors(METADATA.replace(b"Description-Content-Type: text/markdown\n", b""))
    assert CHECKER.metadata_errors(METADATA.replace(b"License-Expression:", b"License:"))
    assert CHECKER.metadata_errors(METADATA.split(b"\n\n")[0] + b"\n\n")


def test_wheel_boundaries_and_native_payload(tmp_path):
    wheel = tmp_path / "candidate.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr("libephemeris/__init__.py", "from ._core import *")
        archive.writestr("libephemeris/py.typed", "")
        archive.writestr("libephemeris/_core.abi3.so", b"native fixture")
        archive.writestr("libephemeris-4.0.0b1.dist-info/METADATA", METADATA)
        archive.writestr("libephemeris-4.0.0b1.dist-info/licenses/core/licenses/index.json", "{}")
    assert CHECKER.check(wheel) == []
    with zipfile.ZipFile(wheel, "a") as archive:
        archive.writestr("libephemeris/__pycache__/cached.pyc", b"not source")
    assert any("build/cache" in error for error in CHECKER.check(wheel))


def test_sdist_requires_both_complete_manuals_and_their_license_notices(tmp_path):
    """A landing page cannot stand in for a missing translation or chapter."""
    required = {
        "Cargo.toml", "COPYING", "THIRD_PARTY_NOTICES.md", "pyproject.toml",
        "python/Cargo.toml", "core/Cargo.toml", "core/src/lib.rs",
        "core/data/leb2/base_core.leb2", "tools/release/check_minimal.py",
        "tools/leb_generation/pyproject.toml", "docs/installation.md",
        "docs/conf.py", "docs/index.rst", "beginners-manual/README.md",
        "beginners-manual/COPYING", "beginners-manual/COPYING.LESSER",
    }
    for language in ("it", "en"):
        required |= {f"beginners-manual/{language}/README.md",
                     f"beginners-manual/{language}/metadata.yaml"}
        required |= {f"beginners-manual/{language}/{number:02d}-chapter.md"
                     for number in range(16)}

    def archive_with(names):
        archive_path = tmp_path / "candidate.tar.gz"
        with tarfile.open(archive_path, "w:gz") as archive:
            for name in sorted(names | {"PKG-INFO"}):
                content = METADATA if name == "PKG-INFO" else b"fixture\n"
                member = tarfile.TarInfo(f"candidate/{name}")
                member.size = len(content)
                archive.addfile(member, io.BytesIO(content))
        return archive_path

    assert CHECKER.check(archive_with(required)) == []
    incomplete = required - {"beginners-manual/it/07-chapter.md", "beginners-manual/COPYING.LESSER"}
    errors = CHECKER.check(archive_with(incomplete))
    assert "missing manual chapter: it/07" in errors
    assert "missing source member: beginners-manual/COPYING.LESSER" in errors
