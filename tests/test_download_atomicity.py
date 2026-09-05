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

2. **The SPK writer honours the configured cache directory.** When
   ``set_spk_cache_dir()`` (or LIBEPHEMERIS_SPK_DIR) names a directory,
   ``download_spk`` writes there — the directory every ``spk_auto`` lookup
   consults. Without an override the file keeps landing in ``<data dir>/spk``,
   so kernels downloaded by earlier versions stay visible.
"""

from __future__ import annotations

import base64
import hashlib
import json
import os
from pathlib import Path

import pytest

import libephemeris.download as dl
from libephemeris import spk, spk_auto, state

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

    def read(self, size=None):
        if self._offset >= len(self._data):
            return b""
        end = len(self._data) if size is None else self._offset + size
        chunk = self._data[self._offset : end]
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
# 2. The SPK writer honours the configured cache directory (issue #55)
# ---------------------------------------------------------------------------


def _horizons_serving(spk_bytes: bytes = b"\x00" * 1024):
    """Stand in for the Horizons SPK-generation endpoint."""
    payload = json.dumps(
        {
            "spk_file_id": "test_id",
            "spk": base64.b64encode(spk_bytes).decode("ascii"),
        }
    ).encode("utf-8")

    def _open_url(req, purpose=None, timeout=None, context=None, **kwargs):
        return _FakeResponse(payload)

    return _open_url


class TestSpkCacheDirIsHonoured:
    """The downloader must write where the configured lookups look."""

    def teardown_method(self):
        state.set_spk_cache_dir(None)

    def test_configured_cache_dir_is_where_the_kernel_lands(
        self, tmp_path, monkeypatch
    ):
        """Issue #55: download_spk() honours set_spk_cache_dir()."""
        monkeypatch.delenv("LIBEPHEMERIS_SPK_DIR", raising=False)
        target = tmp_path / "configured"
        state.set_spk_cache_dir(str(target))
        monkeypatch.setattr(spk, "open_url", _horizons_serving())
        monkeypatch.setattr(spk, "_is_valid_bsp", lambda path: True)

        path = spk.download_spk("2060", "2020-01-01", "2025-01-01")

        assert os.path.dirname(path) == str(target)
        # The same directory every spk_auto lookup resolves to, which is the
        # whole point: a kernel written elsewhere is downloaded again forever.
        assert spk_auto.ensure_cache_dir() == str(target)

    def test_env_var_is_honoured(self, tmp_path, monkeypatch):
        target = tmp_path / "from-env"
        state.set_spk_cache_dir(None)
        monkeypatch.setenv("LIBEPHEMERIS_SPK_DIR", str(target))
        monkeypatch.setattr(spk, "open_url", _horizons_serving())
        monkeypatch.setattr(spk, "_is_valid_bsp", lambda path: True)

        path = spk.download_spk("2060", "2020-01-01", "2025-01-01")

        assert os.path.dirname(path) == str(target)

    def test_without_an_override_the_data_directory_still_owns_the_cache(
        self, tmp_path, monkeypatch
    ):
        """No override: the kernel keeps landing exactly where it used to.

        Relocating this path would hide every kernel already downloaded by an
        installation whose data directory has been redirected.
        """
        monkeypatch.delenv("LIBEPHEMERIS_SPK_DIR", raising=False)
        monkeypatch.setattr(state, "get_spk_cache_dir", lambda: None)
        monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path / "data"))
        monkeypatch.setattr(spk, "open_url", _horizons_serving())
        monkeypatch.setattr(spk, "_is_valid_bsp", lambda path: True)

        path = spk.download_spk("2060", "2020-01-01", "2025-01-01")

        assert os.path.dirname(path) == str(tmp_path / "data" / "spk")


# ---------------------------------------------------------------------------
# 3. A sealed network policy is reported as such, whatever the cache directory
# ---------------------------------------------------------------------------


class TestSealedPolicyBeforeCacheDirectory:
    """download_spk() must not touch the file system before the policy gate.

    Honouring set_spk_cache_dir() put the directory creation ahead of the
    network request; with an uncreatable cache directory a sealed session
    then reported NotADirectoryError instead of NetworkSealedError, and the
    auto-download helpers, which swallow OSError, silently answered None.
    """

    def test_sealed_download_raises_the_typed_error_and_creates_nothing(
        self, tmp_path, monkeypatch
    ):
        from libephemeris.exceptions import NetworkSealedError

        bogus = tmp_path / "not-a-dir"
        bogus.write_text("file, not a directory")
        cache_dir = bogus / "spk"
        from libephemeris import net

        monkeypatch.setenv("LIBEPHEMERIS_SPK_DIR", str(cache_dir))
        net.set_network_policy("sealed")
        try:
            with pytest.raises(NetworkSealedError, match="Horizons SPK generation"):
                spk.download_spk("Ceres", "2000-01-01", "2000-02-01")
        finally:
            net.set_network_policy(None)
        assert not cache_dir.exists()
