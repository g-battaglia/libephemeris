"""Regression tests for project and wheel metadata configuration."""

from __future__ import annotations

import base64
import csv
import hashlib
import io
import stat
import tarfile
import tomllib
import zipfile
from pathlib import Path

import pytest

from scripts.check_wheel_contents import (
    EXPECTED_DIST_INFO,
    EXPECTED_EGG_INFO,
    EXPECTED_ENTRY_POINTS,
    EXPECTED_PROVIDES_EXTRA,
    EXPECTED_REQUIRES_DIST,
    EXPECTED_SDIST_FILENAME,
    EXPECTED_SDIST_ROOT,
    EXPECTED_WHEEL_FILENAME,
    PROJECT_LICENSE,
    PROJECT_NAME,
    PROJECT_VERSION,
    SDIST_REQUIRED,
    WHEEL_METADATA_REQUIRED,
    WHEEL_REQUIRED,
    _member_pax_headers_are_safe,
    _text_payload_hits,
    _zip_info_is_dir,
    _zip_info_is_special,
    audit,
    forbidden_in,
)
from scripts.check_spdx_headers import THIRD_PARTY_COPYRIGHT_LINES


_PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _pyproject() -> dict:
    """Load the authoritative package configuration."""
    return tomllib.loads((_PROJECT_ROOT / "pyproject.toml").read_text())


def _valid_metadata_payload() -> bytes:
    """Build the authenticated install-relevant metadata used by fixtures."""
    lines = [
        "Metadata-Version: 2.4",
        f"Name: {PROJECT_NAME}",
        f"Version: {PROJECT_VERSION}",
        f"License-Expression: {PROJECT_LICENSE}",
        f"Requires-Python: {_pyproject()['project']['requires-python']}",
        *(f"Requires-Dist: {item}" for item in EXPECTED_REQUIRES_DIST),
        *(f"Provides-Extra: {item}" for item in EXPECTED_PROVIDES_EXTRA),
        *(f"License-File: {item}" for item in _pyproject()["project"]["license-files"]),
    ]
    return ("\n".join(lines) + "\n").encode()


def _valid_requires_txt_payload() -> bytes:
    config = _pyproject()["project"]
    lines = [*config["dependencies"]]
    for extra, requirements in config["optional-dependencies"].items():
        lines.extend(("", f"[{extra}]", *requirements))
    return ("\n".join(lines) + "\n").encode()


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


def test_vendored_mit_copyright_notices_are_retained() -> None:
    """The MIT grant's required upstream notices must ship with each file."""
    third_party_notices = (_PROJECT_ROOT / "THIRD_PARTY_NOTICES.md").read_text()

    for rel, comment in THIRD_PARTY_COPYRIGHT_LINES.items():
        notice = comment.removeprefix("# ")
        source_head = (_PROJECT_ROOT / rel).read_text().splitlines()[:6]
        assert comment in source_head
        assert notice in third_party_notices


def test_archive_gate_rejects_clean_room_artifact_names_case_insensitively() -> None:
    """Built archives must enforce the same filename gate as the worktree."""
    prohibited = [
        "libephemeris/data/seasnam.txt",
        "libephemeris/data/sweph.c",
        "libephemeris/vendor/sweprivate.h",
        "libephemeris/Data/REFERENCE/planet.bin",
        "libephemeris/vendor/PySwissEph-2.10/module.py",
        "libephemeris/vendor/PySwissEph-2.10/",
        "./PySwissEph-2.10/",
        "libephemeris/data/./reference/planet.bin",
        "libephemeris/data/x/../reference/planet.bin",
        "libephemeris/vendor/swisseph.zip",
        "libephemeris/vendor/swisseph.cpython-312.so",
        "/absolute/package/module.py",
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


def test_tar_directory_type_is_used_when_name_has_no_trailing_slash() -> None:
    """TarInfo directory metadata closes the final-component name bypass."""
    member = tarfile.TarInfo("pkg/vendor/PySwissEph-2.10")
    member.type = tarfile.DIRTYPE

    assert member.isdir()
    assert forbidden_in([member.name], directory_names=frozenset({member.name})) == [
        member.name
    ]

    for dirname in (
        "pkg/test",
        "pkg/tests",
        "pkg/tests.",
        "pkg/tests ",
        "pkg/compare_scripts",
        "libephemeris/dev_cli",
        "libephemeris/dev_cli.",
        "libephemeris/dev_cli ",
    ):
        test_dir = tarfile.TarInfo(dirname)
        test_dir.type = tarfile.DIRTYPE
        assert forbidden_in(
            [test_dir.name], directory_names=frozenset({test_dir.name})
        ) == [test_dir.name]


def test_zip_unix_metadata_recognizes_directories_and_special_members() -> None:
    """ZIP mode metadata is enforced even when names omit a trailing slash."""
    directory = zipfile.ZipInfo("pkg/vendor/PySwissEph-2.10")
    directory.create_system = 3
    directory.external_attr = (stat.S_IFDIR | 0o755) << 16
    symlink = zipfile.ZipInfo("pkg/link")
    symlink.create_system = 3
    symlink.external_attr = (stat.S_IFLNK | 0o777) << 16

    assert _zip_info_is_dir(directory)
    assert _zip_info_is_special(symlink)
    assert forbidden_in(
        [directory.filename], directory_names=frozenset({directory.filename})
    ) == [directory.filename]


def test_archive_text_payloads_receive_provenance_scan() -> None:
    """Safe-named package text cannot hide blocked implementation tokens."""
    assert _text_payload_hits("libephemeris/runtime.py", b"PyMeeus\n") == [
        (1, "pymeeus")
    ]
    assert _text_payload_hits("libephemeris/README.md", b"GPL-3.0\n") == [
        (1, "copyleft")
    ]
    assert _text_payload_hits("README.md", b"GPL-3.0\n") == []
    assert _text_payload_hits("docs/audit.txt", b"\x00PyMeeus\n") == [(1, "pymeeus")]
    assert _text_payload_hits("docs/latin1.txt", b"caf\xe9 PyMeeus\n") == [
        (0, "content-decode-error")
    ]
    assert _text_payload_hits("libephemeris-extra/README.md", b"GPL-3.0\n") == [
        (1, "copyleft")
    ]
    assert (
        _text_payload_hits(
            "libephemeris-3.0.0rc8/README.md",
            b"GPL-3.0\n",
            sdist_root="libephemeris-3.0.0rc8",
        )
        == []
    )
    assert (
        _text_payload_hits("libephemeris-3.0.0rc8.dist-info/METADATA", b"GPL-3.0\n")
        == []
    )
    assert (
        _text_payload_hits(
            "libephemeris-3.0.0rc8.dist-info/licenses/LICENSING.md",
            b"GPL-3.0\n",
        )
        == []
    )
    assert _text_payload_hits(
        "LIBEPHEMERIS./runtime.py",
        b"PyMeeus # provenance-implementation-ok\n",
    ) == [(1, "pymeeus")]
    assert _text_payload_hits("SCRIPTS/CHECK_PROVENANCE.PY", b"PyMeeus\n") == [
        (1, "pymeeus")
    ]
    assert _text_payload_hits("scripts/check_provenance.py.", b"PyMeeus\n") == [
        (1, "pymeeus")
    ]
    assert _text_payload_hits("README.md.", b"GPL-3.0\n") == [(1, "copyleft")]
    assert _text_payload_hits(
        f"nested/{EXPECTED_DIST_INFO}/METADATA", b"GPL-3.0\n"
    ) == [(1, "copyleft")]
    assert _text_payload_hits(
        f"{EXPECTED_DIST_INFO}/licenses/LICENSING.md.", b"GPL-3.0\n"
    ) == [(1, "copyleft")]


def _write_minimal_valid_wheel(
    path: Path,
    *,
    metadata_payload: bytes | None = None,
    license_payload: bytes | None = None,
    real_repository_payloads: bool = False,
    wheel_payload: bytes | None = None,
    record_payload: bytes | None = None,
    entry_points_payload: bytes | None = None,
    omit: frozenset[str] = frozenset(),
) -> None:
    payloads: dict[str, bytes] = {}
    for name in WHEEL_REQUIRED:
        if name in omit:
            continue
        payloads[name] = (
            (_PROJECT_ROOT / name).read_bytes() if real_repository_payloads else b""
        )
    if metadata_payload is None:
        metadata_payload = _valid_metadata_payload()
    payloads[f"{EXPECTED_DIST_INFO}/METADATA"] = metadata_payload
    payloads[f"{EXPECTED_DIST_INFO}/WHEEL"] = (
        b"Wheel-Version: 1.0\n"
        b"Generator: test\n"
        b"Root-Is-Purelib: true\n"
        b"Tag: py3-none-any\n"
        if wheel_payload is None
        else wheel_payload
    )
    payloads[f"{EXPECTED_DIST_INFO}/entry_points.txt"] = (
        EXPECTED_ENTRY_POINTS if entry_points_payload is None else entry_points_payload
    )
    payloads[f"{EXPECTED_DIST_INFO}/top_level.txt"] = f"{PROJECT_NAME}\n".encode()
    for name in WHEEL_METADATA_REQUIRED:
        basename = name.rsplit("/", maxsplit=1)[-1]
        if "/licenses/" not in name:
            continue
        payload = (_PROJECT_ROOT / basename).read_bytes()
        if basename == "LICENSE" and license_payload is not None:
            payload = license_payload
        payloads[name] = payload

    record_name = f"{EXPECTED_DIST_INFO}/RECORD"
    rows = []
    for name, payload in sorted(payloads.items()):
        digest = "sha256=" + base64.urlsafe_b64encode(
            hashlib.sha256(payload).digest()
        ).rstrip(b"=").decode("ascii")
        rows.append((name, digest, str(len(payload))))
    rows.append((record_name, "", ""))
    record_stream = io.StringIO(newline="")
    csv.writer(record_stream, lineterminator="\n").writerows(rows)
    payloads[record_name] = (
        record_stream.getvalue().encode() if record_payload is None else record_payload
    )

    with zipfile.ZipFile(path, "w") as zf:
        for name, payload in payloads.items():
            zf.writestr(name, payload)


def _add_tar_bytes(tf: tarfile.TarFile, name: str, payload: bytes) -> None:
    info = tarfile.TarInfo(name)
    info.size = len(payload)
    tf.addfile(info, io.BytesIO(payload))


def _write_valid_sdist(path: Path) -> None:
    metadata = _valid_metadata_payload()
    with tarfile.open(path, "w:gz", format=tarfile.PAX_FORMAT) as tf:
        for name in SDIST_REQUIRED:
            if name in {"PKG-INFO", f"{EXPECTED_EGG_INFO}/PKG-INFO"}:
                payload = metadata
            elif name == f"{EXPECTED_EGG_INFO}/entry_points.txt":
                payload = EXPECTED_ENTRY_POINTS
            elif name == f"{EXPECTED_EGG_INFO}/top_level.txt":
                payload = f"{PROJECT_NAME}\n".encode()
            elif name == f"{EXPECTED_EGG_INFO}/dependency_links.txt":
                payload = b"\n"
            elif name == f"{EXPECTED_EGG_INFO}/requires.txt":
                payload = _valid_requires_txt_payload()
            elif name == f"{EXPECTED_EGG_INFO}/SOURCES.txt":
                sources = sorted(set(SDIST_REQUIRED) - {"PKG-INFO", "setup.cfg"})
                payload = "\n".join(sources).encode()
            elif name == "setup.cfg":
                payload = b"[egg_info]\ntag_build = \ntag_date = 0\n"
            elif (_PROJECT_ROOT / name).is_file():
                payload = (_PROJECT_ROOT / name).read_bytes()
            else:
                payload = b""
            _add_tar_bytes(tf, f"{EXPECTED_SDIST_ROOT}/{name}", payload)


def test_archive_audit_requires_authenticated_metadata_and_legal_files(
    tmp_path: Path, capsys
) -> None:
    incomplete = tmp_path / EXPECTED_WHEEL_FILENAME
    with zipfile.ZipFile(incomplete, "w") as zf:
        for name in WHEEL_REQUIRED:
            zf.writestr(name, b"")

    assert audit(None, incomplete) == 1
    output = capsys.readouterr().out
    assert f"MISSING from wheel: {EXPECTED_DIST_INFO}/METADATA" in output
    assert f"MISSING from wheel: {EXPECTED_DIST_INFO}/licenses/LICENSE" in output


def test_archive_audit_rejects_false_metadata_and_legal_payloads(
    tmp_path: Path, capsys
) -> None:
    wheel = tmp_path / EXPECTED_WHEEL_FILENAME
    _write_minimal_valid_wheel(
        wheel,
        metadata_payload=b"Metadata-Version: 2.4\nName: foreign\nVersion: 0\n",
        license_payload=b"not a license\n",
    )

    assert audit(None, wheel) == 1
    output = capsys.readouterr().out
    assert "metadata-name-mismatch" in output
    assert "legal-content-mismatch" in output


def test_archive_audit_rejects_duplicate_names_and_egg_info_in_wheel(
    tmp_path: Path, capsys
) -> None:
    wheel = tmp_path / EXPECTED_WHEEL_FILENAME
    _write_minimal_valid_wheel(wheel)
    with pytest.warns(UserWarning, match="Duplicate name"):
        with zipfile.ZipFile(wheel, "a") as zf:
            zf.writestr("libephemeris/__init__.py", b"")
            zf.writestr("libephemeris.egg-info/PKG-INFO", b"")

    assert audit(None, wheel) == 1
    output = capsys.readouterr().out
    assert "FORBIDDEN in wheel: libephemeris/__init__.py" in output
    assert "FORBIDDEN in wheel: libephemeris.egg-info/PKG-INFO" in output


def test_archive_audit_rejects_wrong_sdist_root(tmp_path: Path, capsys) -> None:
    wheel = tmp_path / EXPECTED_WHEEL_FILENAME
    _write_minimal_valid_wheel(wheel)
    sdist = tmp_path / "wrong.tar.gz"
    with tarfile.open(sdist, "w:gz") as tf:
        member = tarfile.TarInfo("wrong-root/empty.txt")
        member.size = 0
        tf.addfile(member)

    assert audit(sdist, wheel) == 1
    assert "MISSING from sdist: archive root" in capsys.readouterr().out


def test_foreign_metadata_namespace_gets_no_legal_waiver() -> None:
    assert forbidden_in(["foreign-1.0.dist-info/METADATA"]) == [
        "foreign-1.0.dist-info/METADATA"
    ]
    assert _text_payload_hits("foreign-1.0.dist-info/METADATA", b"GPL-3.0\n") == [
        (1, "copyleft")
    ]


def test_untrusted_sdist_root_cannot_strip_package_namespace() -> None:
    assert _text_payload_hits(
        "libephemeris/runtime.py",
        b"PyMeeus # provenance-implementation-ok\n",
        sdist_root="libephemeris",
    ) == [(1, "pymeeus")]
    assert _text_payload_hits(
        f"{EXPECTED_SDIST_ROOT}/{EXPECTED_EGG_INFO}/PKG-INFO.",
        b"GPL-3.0\n",
        sdist_root=EXPECTED_SDIST_ROOT,
    ) == [(1, "copyleft")]


def test_metadata_namespaces_are_allowed_only_at_exact_archive_paths() -> None:
    allowed = [
        f"{EXPECTED_DIST_INFO}/METADATA",
        f"{EXPECTED_SDIST_ROOT}/{EXPECTED_EGG_INFO}/PKG-INFO",
    ]
    prohibited = [
        f"nested/{EXPECTED_DIST_INFO}/METADATA",
        f"{EXPECTED_SDIST_ROOT}/nested/{EXPECTED_EGG_INFO}/PKG-INFO",
        f"{EXPECTED_DIST_INFO.upper()}/METADATA",
        f"{EXPECTED_DIST_INFO}./METADATA",
    ]

    assert forbidden_in(allowed) == []
    assert forbidden_in(prohibited) == sorted(prohibited)


def test_archive_audit_accepts_exact_authenticated_rc8_pair(
    tmp_path: Path,
) -> None:
    wheel = tmp_path / EXPECTED_WHEEL_FILENAME
    sdist = tmp_path / EXPECTED_SDIST_FILENAME
    _write_minimal_valid_wheel(wheel, real_repository_payloads=True)
    _write_valid_sdist(sdist)

    assert audit(sdist, wheel) == 0


def test_archive_audit_rejects_wrong_rc8_artifact_filenames(
    tmp_path: Path, capsys
) -> None:
    wheel = tmp_path / "libephemeris-3.0.0rc7-py3-none-any.whl"
    sdist = tmp_path / "libephemeris-3.0.0rc7.tar.gz"
    _write_minimal_valid_wheel(wheel)
    _write_valid_sdist(sdist)

    assert audit(sdist, wheel) == 1
    output = capsys.readouterr().out
    assert f"expected artifact filename {EXPECTED_WHEEL_FILENAME}" in output
    assert f"expected artifact filename {EXPECTED_SDIST_FILENAME}" in output


def test_archive_audit_rejects_conflicting_duplicate_identity_headers(
    tmp_path: Path, capsys
) -> None:
    metadata = (
        "Metadata-Version: 2.4\n"
        f"Name: {PROJECT_NAME}\n"
        "Name: hostile-name\n"
        f"Version: {PROJECT_VERSION}\n"
        "Version: 0\n"
        f"License-Expression: {PROJECT_LICENSE}\n"
        "License-Expression: GPL-3.0\n"
        f"Requires-Python: {_pyproject()['project']['requires-python']}\n"
    ).encode()
    wheel = tmp_path / EXPECTED_WHEEL_FILENAME
    _write_minimal_valid_wheel(wheel, metadata_payload=metadata)

    assert audit(None, wheel) == 1
    output = capsys.readouterr().out
    assert "metadata-name-mismatch" in output
    assert "metadata-version-mismatch" in output
    assert "metadata-license-expression-mismatch" in output


def test_archive_audit_rejects_unexpected_install_dependency(
    tmp_path: Path, capsys
) -> None:
    metadata = _valid_metadata_payload() + b"Requires-Dist: hostile-package\n"
    wheel = tmp_path / EXPECTED_WHEEL_FILENAME
    _write_minimal_valid_wheel(wheel, metadata_payload=metadata)

    assert audit(None, wheel) == 1
    assert "metadata-requires-dist-mismatch" in capsys.readouterr().out


def test_archive_audit_requires_every_declared_package_module(
    tmp_path: Path, capsys
) -> None:
    omitted = "libephemeris/planets.py"
    assert omitted in WHEEL_REQUIRED
    wheel = tmp_path / EXPECTED_WHEEL_FILENAME
    _write_minimal_valid_wheel(wheel, omit=frozenset({omitted}))

    assert audit(None, wheel) == 1
    assert f"MISSING from wheel: {omitted}" in capsys.readouterr().out


def test_archive_audit_rejects_untracked_extra_package_module(
    tmp_path: Path, capsys
) -> None:
    wheel = tmp_path / EXPECTED_WHEEL_FILENAME
    _write_minimal_valid_wheel(wheel)
    with zipfile.ZipFile(wheel, "a") as zf:
        zf.writestr("libephemeris/untracked_extra.py", b"plain = True\n")

    assert audit(None, wheel) == 1
    assert (
        "FORBIDDEN in wheel: libephemeris/untracked_extra.py" in capsys.readouterr().out
    )


@pytest.mark.parametrize(
    ("kwargs", "expected_label"),
    [
        ({}, "repository-content-mismatch"),
        ({"wheel_payload": b""}, "wheel-wheel-version-mismatch"),
        ({"record_payload": b""}, "record-member-set-mismatch"),
        ({"entry_points_payload": b"[console_scripts]\n"}, "entry-points-mismatch"),
    ],
)
def test_archive_audit_rejects_unauthenticated_or_malformed_wheel_payloads(
    tmp_path: Path,
    capsys,
    kwargs: dict[str, bytes],
    expected_label: str,
) -> None:
    wheel = tmp_path / EXPECTED_WHEEL_FILENAME
    _write_minimal_valid_wheel(wheel, **kwargs)

    assert audit(None, wheel) == 1
    assert expected_label in capsys.readouterr().out


def test_archive_audit_rejects_zip_comment_and_extra_channels(
    tmp_path: Path, capsys
) -> None:
    labels: list[str] = []
    for case in ("global", "member-comment", "member-extra"):
        directory = tmp_path / case
        directory.mkdir()
        wheel = directory / EXPECTED_WHEEL_FILENAME
        _write_minimal_valid_wheel(wheel)
        with zipfile.ZipFile(wheel, "a") as zf:
            if case == "global":
                zf.comment = b"hidden payload"
            else:
                info = zipfile.ZipInfo(f"{case}.txt")
                if case == "member-comment":
                    info.comment = b"hidden payload"
                else:
                    info.extra = b"\xfe\xca\x07\x00payload"
                zf.writestr(info, b"plain")
        assert audit(None, wheel) == 1
        labels.append(capsys.readouterr().out)

    assert "zip-global-comment" in labels[0]
    assert "zip-member-comment" in labels[1]
    assert "zip-member-extra" in labels[2]


def test_archive_audit_rejects_tar_pax_metadata_channels(
    tmp_path: Path, capsys
) -> None:
    wheel = tmp_path / EXPECTED_WHEEL_FILENAME
    _write_minimal_valid_wheel(wheel)

    member_sdist = tmp_path / EXPECTED_SDIST_FILENAME
    with tarfile.open(member_sdist, "w:gz", format=tarfile.PAX_FORMAT) as tf:
        info = tarfile.TarInfo(f"{EXPECTED_SDIST_ROOT}/README.md")
        info.pax_headers = {"comment": "hidden payload"}
        payload = b"plain"
        info.size = len(payload)
        tf.addfile(info, io.BytesIO(payload))
    assert audit(member_sdist, wheel) == 1
    assert "pax-member-header" in capsys.readouterr().out

    global_dir = tmp_path / "global"
    global_dir.mkdir()
    global_sdist = global_dir / EXPECTED_SDIST_FILENAME
    with tarfile.open(
        global_sdist,
        "w:gz",
        format=tarfile.PAX_FORMAT,
        pax_headers={"comment": "hidden payload"},
    ) as tf:
        _add_tar_bytes(tf, f"{EXPECTED_SDIST_ROOT}/README.md", b"plain")
    assert audit(global_sdist, wheel) == 1
    assert "pax-global-header" in capsys.readouterr().out


def test_archive_audit_allows_only_canonical_member_mtime_pax(
    tmp_path: Path, capsys
) -> None:
    wheel = tmp_path / EXPECTED_WHEEL_FILENAME
    _write_minimal_valid_wheel(wheel)
    sdist = tmp_path / EXPECTED_SDIST_FILENAME
    mtime = "1767225600.125"
    with tarfile.open(sdist, "w:gz", format=tarfile.PAX_FORMAT) as tf:
        info = tarfile.TarInfo(f"{EXPECTED_SDIST_ROOT}/README.md")
        info.mtime = float(mtime)
        info.pax_headers = {"mtime": mtime}
        payload = b"plain"
        info.size = len(payload)
        tf.addfile(info, io.BytesIO(payload))

    with tarfile.open(sdist) as tf:
        member = tf.getmembers()[0]
        assert member.pax_headers == {"mtime": mtime}
        assert _member_pax_headers_are_safe(member)

    assert audit(sdist, wheel) == 1  # The intentionally minimal sdist is incomplete.
    assert "pax-member-header" not in capsys.readouterr().out


@pytest.mark.parametrize(
    ("headers", "resolved_mtime"),
    [
        ({"comment": "hidden"}, 1.25),
        ({"mtime": "1.25", "comment": "hidden"}, 1.25),
        ({"mtime": "nan"}, float("nan")),
        ({"mtime": "inf"}, float("inf")),
        ({"mtime": "1e3"}, 1000.0),
        ({"mtime": "+1.25"}, 1.25),
        ({"mtime": "-1.25"}, -1.25),
        ({"mtime": "01.25"}, 1.25),
        ({"mtime": ".25"}, 0.25),
        ({"mtime": "1."}, 1.0),
        ({"mtime": "1.25 "}, 1.25),
        ({"mtime": "1.25"}, 2.5),
        ({"mtime": "9" * 33}, float("9" * 33)),
    ],
)
def test_member_pax_mtime_waiver_rejects_malformed_or_incoherent_values(
    headers: dict[str, str], resolved_mtime: float
) -> None:
    member = tarfile.TarInfo("member")
    member.pax_headers = headers
    member.mtime = resolved_mtime

    assert not _member_pax_headers_are_safe(member)


def test_archive_audit_requires_sdist_package_and_authentic_project_files(
    tmp_path: Path, capsys
) -> None:
    wheel = tmp_path / EXPECTED_WHEEL_FILENAME
    _write_minimal_valid_wheel(wheel)
    sdist = tmp_path / EXPECTED_SDIST_FILENAME
    with tarfile.open(sdist, "w:gz") as tf:
        for name in (
            "LICENSE",
            "LICENSING.md",
            "NOTICE.md",
            "THIRD_PARTY_NOTICES.md",
            "PKG-INFO",
            "pyproject.toml",
        ):
            payload = (
                (_PROJECT_ROOT / name).read_bytes()
                if name != "PKG-INFO"
                else b"Metadata-Version: 2.4\n"
            )
            _add_tar_bytes(tf, f"{EXPECTED_SDIST_ROOT}/{name}", payload)

    assert audit(sdist, wheel) == 1
    output = capsys.readouterr().out
    assert "MISSING from sdist: libephemeris/__init__.py" in output


@pytest.mark.parametrize(
    ("name", "payload", "expected_label"),
    [
        (
            "setup.cfg",
            b"[egg_info]\ntag_build =\ntag_date = 0\n[install]\nprefix = /tmp\n",
            "setup-cfg-section-mismatch",
        ),
        (
            f"{EXPECTED_EGG_INFO}/requires.txt",
            b"hostile-package\n",
            "requires-txt-section-mismatch",
        ),
        (
            f"{EXPECTED_EGG_INFO}/SOURCES.txt",
            b"pyproject.toml\n",
            "sources-txt-set-mismatch",
        ),
    ],
)
def test_archive_audit_authenticates_generated_sdist_control_files(
    tmp_path: Path,
    capsys,
    name: str,
    payload: bytes,
    expected_label: str,
) -> None:
    wheel = tmp_path / EXPECTED_WHEEL_FILENAME
    _write_minimal_valid_wheel(wheel)
    sdist = tmp_path / EXPECTED_SDIST_FILENAME
    with tarfile.open(sdist, "w:gz") as tf:
        _add_tar_bytes(tf, f"{EXPECTED_SDIST_ROOT}/{name}", payload)

    assert audit(sdist, wheel) == 1
    assert expected_label in capsys.readouterr().out


def test_archive_audit_rejects_case_alias_as_second_sdist_root(
    tmp_path: Path, capsys
) -> None:
    wheel = tmp_path / EXPECTED_WHEEL_FILENAME
    _write_minimal_valid_wheel(wheel)
    sdist = tmp_path / EXPECTED_SDIST_FILENAME
    with tarfile.open(sdist, "w:gz") as tf:
        _add_tar_bytes(tf, f"{EXPECTED_SDIST_ROOT}/README.md", b"plain")
        _add_tar_bytes(tf, f"{EXPECTED_SDIST_ROOT.upper()}/extra.txt", b"plain")

    assert audit(sdist, wheel) == 1
    assert "MISSING from sdist: archive root" in capsys.readouterr().out
    (SDIST_REQUIRED,)
