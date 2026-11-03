"""The installed binary distribution must carry its upstream license texts."""

import importlib.metadata
import json
from pathlib import PurePosixPath


def test_installed_dependency_and_project_licenses():
    distribution = importlib.metadata.distribution("libephemeris")
    assert distribution.metadata["License-Expression"] == "AGPL-3.0"
    files = {str(path).replace("\\", "/") for path in distribution.files or ()}
    suffix = "core/licenses/index.json"
    indices = sorted(path for path in files if path.endswith("/licenses/" + suffix))
    assert len(indices) == 1, "the wheel must contain its Cargo dependency license bundle"
    index_path = indices[0]
    license_root = index_path.removesuffix(suffix)
    index = json.loads(distribution.locate_file(index_path).read_text(encoding="utf-8"))
    assert index["schema"] == "cargo-dependency-licenses/v1"
    packages = index["packages"]
    assert {"pyo3", "pyo3-log", "ureq", "zstd-sys"} <= {p["name"] for p in packages}
    identities = [(p["name"], p["version"], p["source"]) for p in packages]
    assert len(identities) == len(set(identities))
    license_headers = set(distribution.metadata.get_all("License-File") or ())
    for package in packages:
        assert package["version"] and package["files"]
        assert "license" in package and "source" in package
        for notice in package["files"]:
            relative = PurePosixPath(notice["path"])
            assert not relative.is_absolute() and ".." not in relative.parts
            bundled = "core/licenses/" + str(relative)
            assert bundled in license_headers
            assert license_root + bundled in files
            assert distribution.locate_file(license_root + bundled).read_bytes()
    assert "COPYING" in license_headers
    assert b"GNU AFFERO GENERAL PUBLIC LICENSE" in distribution.locate_file(
        license_root + "COPYING"
    ).read_bytes()
    assert "THIRD_PARTY_NOTICES.md" in license_headers
    assert distribution.locate_file(license_root + "THIRD_PARTY_NOTICES.md").read_bytes()
    for removed in ("native/rebound/", "native/assist/"):
        assert not any(path.startswith(removed) for path in license_headers)
        assert not any("/" + removed in "/" + path for path in files)
