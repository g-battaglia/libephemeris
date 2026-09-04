# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for issue #55: a failed download must never lose data.

Three invariants are pinned here.

1. **Nothing is unlinked before a verified replacement exists.** Every
   downloader streams to a temporary file, checks the digest and the
   structure there, and only then publishes with ``os.replace``. A failure at
   any point leaves the artifact already on disk byte-for-byte intact. The
   files involved are 10 MB to 2.6 GB, so a failed repair used to be an
   expensive loss.

2. **The SPK cache has one address.** The writer (``spk.download_spk``) and
   every reader in ``spk_auto`` resolve the directory through the same
   function. When they disagreed, a downloaded kernel was never found again
   and each process re-downloaded it.

3. **Deleting the IERS cache is scoped and previewable.** The deletion
   honours ``set_iers_cache_dir()``, reuses the same path helpers the
   downloaders use, and can report what it would remove without removing it.
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pytest

import libephemeris.download as dl
from libephemeris import iers_data, spk_auto, state

pytestmark = pytest.mark.usefixtures("allow_mocked_network")


class _FakeHeaders:
    def __init__(self, content_length):
        self._cl = content_length

    def get(self, key, default=None):
        if key == "Content-Length":
            return self._cl if self._cl is not None else default
        return default


class _FakeResponse:
    def __init__(self, data: bytes):
        self._data = data
        self._offset = 0
        self.headers = _FakeHeaders(len(data))

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def read(self, size):
        if self._offset >= len(self._data):
            return b""
        chunk = self._data[self._offset : self._offset + size]
        self._offset += len(chunk)
        return chunk


def _serving(data: bytes):
    def _open_url(req, timeout=None, context=None, **kwargs):
        return _FakeResponse(data)

    return _open_url


def _offline(req, timeout=None, context=None, **kwargs):
    raise OSError("network is unreachable")


# ---------------------------------------------------------------------------
# 1. Atomicity: an existing artifact outlives every failure mode
# ---------------------------------------------------------------------------


def test_failed_validation_leaves_existing_file_untouched(tmp_path, monkeypatch):
    """A structurally invalid payload never replaces a good file."""
    dest = tmp_path / "artifact.bsp"
    dest.write_bytes(b"the good copy")
    monkeypatch.setattr(dl, "open_url", _serving(b"truncated payload"))

    with pytest.raises(ValueError, match="failed validation"):
        dl.download_file(
            "http://example.invalid/a.bsp",
            dest,
            description="a.bsp",
            show_progress=False,
            validator=lambda path: False,
        )

    assert dest.read_bytes() == b"the good copy"
    assert list(tmp_path.glob("*.download")) == []


def test_digest_mismatch_leaves_existing_file_untouched(tmp_path, monkeypatch):
    """A payload that fails the pinned digest never replaces a good file."""
    dest = tmp_path / "artifact.bsp"
    dest.write_bytes(b"the good copy")
    monkeypatch.setattr(dl, "open_url", _serving(b"something else"))

    with pytest.raises(ValueError, match="Hash mismatch"):
        dl.download_file(
            "http://example.invalid/a.bsp",
            dest,
            expected_sha256="0" * 64,
            show_progress=False,
        )

    assert dest.read_bytes() == b"the good copy"
    assert list(tmp_path.glob("*.download")) == []


def test_network_failure_leaves_existing_file_untouched(tmp_path, monkeypatch):
    """An unreachable server never costs the user the copy they had."""
    dest = tmp_path / "artifact.bsp"
    dest.write_bytes(b"the good copy")
    monkeypatch.setattr(dl, "open_url", _offline)

    with pytest.raises(OSError, match="network is unreachable"):
        dl.download_file("http://example.invalid/a.bsp", dest, show_progress=False)

    assert dest.read_bytes() == b"the good copy"
    assert list(tmp_path.glob("*.download")) == []


def test_unexpected_validator_error_leaves_no_temp_file(tmp_path, monkeypatch):
    """A validator is arbitrary caller code and may raise anything.

    The cleanup is unconditional for that reason: a narrow except clause
    would let an unforeseen exception orphan a fully written .download file
    in the data directory, which for these artifacts means hundreds of
    megabytes of dead weight.
    """
    dest = tmp_path / "artifact.bsp"
    dest.write_bytes(b"the good copy")
    monkeypatch.setattr(dl, "open_url", _serving(b"new payload"))

    class _Unexpected(BaseException):
        pass

    def _explodes(path: str) -> bool:
        raise _Unexpected("validator defect")

    with pytest.raises(_Unexpected):
        dl.download_file(
            "http://example.invalid/a.bsp",
            dest,
            show_progress=False,
            validator=_explodes,
        )

    assert dest.read_bytes() == b"the good copy"
    assert list(tmp_path.glob("*.download")) == []


def test_validator_runs_before_publication(tmp_path, monkeypatch):
    """The validator sees the temp file, not the destination.

    Validating after ``os.replace`` is what made a corrupt download destroy
    the previous artifact: by then the damage was already committed.
    """
    dest = tmp_path / "artifact.bsp"
    dest.write_bytes(b"the good copy")
    monkeypatch.setattr(dl, "open_url", _serving(b"new payload"))
    seen: list[bytes] = []

    def _validator(path: str) -> bool:
        assert Path(path) != dest
        seen.append(Path(path).read_bytes())
        assert dest.read_bytes() == b"the good copy"
        return True

    dl.download_file(
        "http://example.invalid/a.bsp",
        dest,
        show_progress=False,
        validator=_validator,
    )

    assert seen == [b"new payload"]
    assert dest.read_bytes() == b"new payload"


def test_bundled_install_validation_preserves_destination(tmp_path, monkeypatch):
    """A bundled resource that fails its structural check installs nothing."""
    package = tmp_path / "package"
    source = package / "data" / "core.leb2"
    source.parent.mkdir(parents=True)
    payload = b"structurally broken core"
    source.write_bytes(payload)
    monkeypatch.setattr(dl.resources, "files", lambda package_name: package)

    dest = tmp_path / "base_core.leb2"
    dest.write_bytes(b"the installed core")

    with pytest.raises(ValueError, match="failed validation"):
        dl._install_bundled_file(
            "data/core.leb2",
            dest,
            expected_sha256=hashlib.sha256(payload).hexdigest(),
            validator=lambda path: False,
        )

    assert dest.read_bytes() == b"the installed core"
    assert list(tmp_path.glob("*.bundled")) == []


def test_leb2_group_survives_a_corrupt_replacement(tmp_path, monkeypatch, capsys):
    """A corrupt LEB2 download leaves the installed group in place."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    installed = leb_dir / "medium_core.leb2"
    installed.write_bytes(b"the installed core")

    monkeypatch.setattr(dl, "open_url", _serving(b"truncated"))
    monkeypatch.setattr(dl, "_is_valid_leb", lambda path: False)
    monkeypatch.setitem(
        dl.DATA_FILES,
        "medium_core.leb2",
        {
            "url": "https://example.invalid/medium_core.leb2",
            "sha256": hashlib.sha256(b"truncated").hexdigest(),
            "size_mb": 1,
            "description": "test",
            "dest_subdir": "leb",
        },
    )

    result = dl.download_leb2_for_tier(
        "medium", groups=["core"], force=True, quiet=False, activate=False
    )

    assert result == []
    assert "[FAIL]" in capsys.readouterr().out
    assert installed.read_bytes() == b"the installed core"


def test_assist_file_survives_a_failed_reverification(tmp_path, monkeypatch):
    """An ASSIST file that fails verification is kept until it is replaced."""
    import libephemeris.rebound_integration as ri

    dest = tmp_path / "linux_p1550p2650.440"
    dest.write_bytes(b"the copy already downloaded")

    monkeypatch.setattr(ri, "_verify_assist_file", lambda path, name=None: False)

    def _boom(**kwargs):
        raise OSError("network is unreachable")

    monkeypatch.setattr(ri, "_download_single_file", _boom)

    with pytest.raises(OSError, match="network is unreachable"):
        ri.download_assist_data(
            target_dir=str(tmp_path), planets=True, asteroids=False, quiet=True
        )

    assert dest.read_bytes() == b"the copy already downloaded"


# ---------------------------------------------------------------------------
# 2. The SPK cache resolves to one directory for readers and writers alike
# ---------------------------------------------------------------------------


class TestSpkCacheDirIsShared:
    """The lookup path and the download path must never diverge."""

    def teardown_method(self):
        state.set_spk_cache_dir(None)

    def test_setter_is_honoured_by_every_consumer(self, tmp_path, monkeypatch):
        monkeypatch.delenv("LIBEPHEMERIS_SPK_DIR", raising=False)
        target = tmp_path / "configured"
        state.set_spk_cache_dir(str(target))

        resolved = spk_auto.resolve_cache_dir()
        assert resolved == str(target)
        assert os.path.dirname(spk_auto.get_cache_path("2060")) == resolved
        assert spk_auto.cache_info()["cache_dir"] == resolved
        assert spk_auto.get_cache_info()["cache_dir"] == resolved
        assert spk_auto.ensure_cache_dir() == resolved

        config = spk_auto.AutoSpkConfig(ipl=15, body_id="2060")
        assert config.get_cache_dir() == resolved

    def test_env_var_is_honoured(self, tmp_path, monkeypatch):
        target = tmp_path / "from-env"
        monkeypatch.setenv("LIBEPHEMERIS_SPK_DIR", str(target))
        state.set_spk_cache_dir(None)

        assert spk_auto.resolve_cache_dir() == str(target)

    def test_data_dir_redirect_keeps_the_cache_inside_it(self, tmp_path, monkeypatch):
        """A redirected data directory owns its own spk/ subdirectory.

        This is where ``download_spk`` has always written; the readers now
        agree instead of looking under the home directory.
        """
        monkeypatch.delenv("LIBEPHEMERIS_SPK_DIR", raising=False)
        monkeypatch.setattr(state, "get_spk_cache_dir", lambda: None)
        monkeypatch.setenv("LIBEPHEMERIS_DATA_DIR", str(tmp_path / "data"))

        assert spk_auto.resolve_cache_dir() == str(tmp_path / "data" / "spk")

    def test_default_stays_under_the_home_directory(self, monkeypatch):
        monkeypatch.delenv("LIBEPHEMERIS_SPK_DIR", raising=False)
        monkeypatch.delenv("LIBEPHEMERIS_DATA_DIR", raising=False)
        monkeypatch.setattr(state, "get_spk_cache_dir", lambda: None)
        monkeypatch.setattr(state, "_resolve_data_dir", lambda: state.DEFAULT_DATA_DIR)

        assert spk_auto.resolve_cache_dir() == os.path.abspath(
            spk_auto.DEFAULT_AUTO_SPK_DIR
        )


# ---------------------------------------------------------------------------
# 3. Deleting the IERS cache is scoped, previewable and side-effect free
# ---------------------------------------------------------------------------


class TestDeleteIersCacheFiles:
    def teardown_method(self):
        iers_data.set_iers_cache_dir(None)

    def test_dry_run_reports_without_deleting(self, tmp_path):
        cache = tmp_path / "iers"
        cache.mkdir()
        iers_data.set_iers_cache_dir(str(cache))
        names = ("finals2000A.data", "leap_seconds.dat", "deltat.data")
        for name in names:
            (cache / name).write_text("x")

        assert iers_data.delete_iers_cache_files(dry_run=True) == 3
        assert all((cache / name).exists() for name in names)

        assert iers_data.delete_iers_cache_files() == 3
        assert not any((cache / name).exists() for name in names)

    def test_only_the_configured_directory_is_touched(self, tmp_path):
        """The override is honoured, so the real cache is never at risk."""
        other = tmp_path / "untouched"
        other.mkdir()
        (other / "finals2000A.data").write_text("keep me")

        cache = tmp_path / "iers"
        cache.mkdir()
        (cache / "finals2000A.data").write_text("delete me")
        iers_data.set_iers_cache_dir(str(cache))

        assert iers_data.delete_iers_cache_files() == 1
        assert (other / "finals2000A.data").read_text() == "keep me"

    def test_does_not_create_the_cache_directory(self, tmp_path):
        """A deletion must not have the side effect of creating a directory."""
        cache = tmp_path / "never-created"
        iers_data.set_iers_cache_dir(str(cache))

        assert iers_data.delete_iers_cache_files(dry_run=True) == 0
        assert iers_data.delete_iers_cache_files() == 0
        assert not cache.exists()

    def test_paths_come_from_the_shared_helpers(self, tmp_path):
        """The removed set cannot drift from the written set."""
        cache = tmp_path / "iers"
        cache.mkdir()
        iers_data.set_iers_cache_dir(str(cache))
        for path in (
            iers_data._get_finals_cache_path(),
            iers_data._get_leap_seconds_cache_path(),
            iers_data._get_delta_t_cache_path(),
        ):
            Path(path).write_text("x")

        assert iers_data.delete_iers_cache_files() == 3
