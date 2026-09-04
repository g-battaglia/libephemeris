# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for ``libephemeris.download``.

This module performs network / file-IO / optional-dependency work. The
project policy for such modules is full mocking: every reachable line is
driven through fakes (fake urlopen responses, fake SPK/LEB validators,
patched optional dependencies) so no test ever touches the real network.

Each test targets a specific group of branches in ``download.py``:
progress bars, the atomic download pipeline, hash verification, the data
status report, the multi-step initializers, and the tier downloaders.
"""

from __future__ import annotations

import hashlib
import sys
import types
from pathlib import Path

import pytest

import libephemeris.download as dl

pytestmark = pytest.mark.usefixtures("allow_mocked_network")


# ---------------------------------------------------------------------------
# Test doubles
# ---------------------------------------------------------------------------


class _FakeHeaders:
    """Minimal stand-in for ``response.headers`` with ``.get``."""

    def __init__(self, content_length):
        self._cl = content_length

    def get(self, key, default=None):
        if key == "Content-Length":
            return self._cl if self._cl is not None else default
        return default


class _FakeResponse:
    """Context-manager response double whose ``read`` yields chunks."""

    def __init__(self, data: bytes, content_length):
        self._data = data
        self._offset = 0
        self.headers = _FakeHeaders(content_length)

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


def _make_urlopen(data: bytes, content_length):
    """Build a fake ``net.open_url`` that ignores args and returns a response."""

    def _fake_urlopen(req, timeout=None, context=None, **kwargs):
        return _FakeResponse(data, content_length)

    return _fake_urlopen


def _fake_download_file(payload: bytes = b"fresh"):
    """Stand in for ``download_file`` while honouring its atomicity contract.

    The real function writes ``payload`` to a temp file, runs ``validator``
    on it, and only publishes over ``dest_path`` when the check passes. A
    stub that writes straight to ``dest_path`` would silently exempt callers
    from the very guarantee these tests exist to pin down.
    """

    def _download(url, dest_path, validator=None, **kwargs):
        dest_path = Path(dest_path)
        temp_path = dest_path.parent / (dest_path.name + ".fake-download")
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        temp_path.write_bytes(payload)
        if validator is not None and not validator(str(temp_path)):
            temp_path.unlink()
            raise ValueError(
                f"Downloaded file {kwargs.get('description', dest_path.name)} "
                "failed validation - corrupt or incomplete"
            )
        temp_path.replace(dest_path)
        return True

    return _download


# ---------------------------------------------------------------------------
# get_data_dir / _format_size
# ---------------------------------------------------------------------------


def test_get_data_dir_returns_path(tmp_path, monkeypatch):
    """get_data_dir wraps state._get_data_dir in a Path (219-221)."""
    from libephemeris import state

    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    result = dl.get_data_dir()
    assert result == Path(tmp_path)


def test_format_size_all_units():
    """_format_size walks each unit and the TB fall-through (226-230)."""
    assert dl._format_size(512) == "512.0 B"
    assert dl._format_size(2048) == "2.0 KB"
    assert dl._format_size(5 * 1024 * 1024) == "5.0 MB"
    assert dl._format_size(3 * 1024**3) == "3.0 GB"
    # Larger than GB exhausts the loop and hits the TB return.
    assert dl._format_size(2 * 1024**4) == "2.0 TB"


# ---------------------------------------------------------------------------
# SimpleProgressBar
# ---------------------------------------------------------------------------


def test_simple_progress_bar_full_cycle():
    """Exercise __init__, update/_render, percent caching, close (243-283)."""
    import io

    buf = io.StringIO()
    bar = dl.SimpleProgressBar(total=100, description="dl", width=10, file=buf)
    assert bar.current == 0
    assert bar.file is buf

    bar.update(50)  # renders 50%
    # Second update to same percent -> early return on cached percent.
    bar.update(0)
    bar.update(50)  # 100%
    bar.close()

    out = buf.getvalue()
    assert "50%" in out
    assert "100%" in out
    assert out.endswith("\n")


def test_simple_progress_bar_zero_total():
    """_render returns immediately when total == 0 (257-258)."""
    import io

    buf = io.StringIO()
    bar = dl.SimpleProgressBar(total=0, description="dl", file=buf)
    bar.update(10)
    assert buf.getvalue() == ""


def test_simple_progress_bar_default_file(monkeypatch):
    """file defaults to sys.stdout when none provided (247)."""
    bar = dl.SimpleProgressBar(total=10)
    assert bar.file is sys.stdout


# ---------------------------------------------------------------------------
# _get_progress_bar
# ---------------------------------------------------------------------------


def test_get_progress_bar_rich_path(monkeypatch):
    """Inject a fake rich.progress so the RichProgressWrapper is used (299-327)."""
    rich_pkg = types.ModuleType("rich")
    rich_progress = types.ModuleType("rich.progress")

    class _FakeProgress:
        def __init__(self, *columns):
            self.columns = columns
            self.started = False
            self.stopped = False
            self.updates = []

        def start(self):
            self.started = True

        def add_task(self, description, total):
            self.description = description
            self.total = total
            return 7

        def update(self, task, advance):
            self.updates.append((task, advance))

        def stop(self):
            self.stopped = True

    rich_progress.Progress = _FakeProgress
    rich_progress.BarColumn = lambda *a, **k: "bar"
    rich_progress.DownloadColumn = lambda *a, **k: "dl"
    rich_progress.TransferSpeedColumn = lambda *a, **k: "speed"
    rich_progress.TimeRemainingColumn = lambda *a, **k: "time"
    rich_pkg.progress = rich_progress

    monkeypatch.setitem(sys.modules, "rich", rich_pkg)
    monkeypatch.setitem(sys.modules, "rich.progress", rich_progress)

    bar = dl._get_progress_bar(1000, "download")
    bar.update(250)
    bar.close()
    assert bar.task == 7
    assert bar.progress.started is True
    assert bar.progress.stopped is True
    assert bar.progress.updates == [(7, 250)]


def test_get_progress_bar_import_error_fallback(monkeypatch):
    """rich absent -> SimpleProgressBar fallback (328-332)."""
    monkeypatch.setitem(sys.modules, "rich.progress", None)
    bar = dl._get_progress_bar(1000, "download")
    assert isinstance(bar, dl.SimpleProgressBar)


# ---------------------------------------------------------------------------
# download_file
# ---------------------------------------------------------------------------


def test_download_file_with_progress(tmp_path, monkeypatch):
    """total_size>0 + show_progress drives the progress branch (358-416)."""
    data = b"hello world payload" * 10
    monkeypatch.setattr(dl, "open_url", _make_urlopen(data, len(data)))

    captured = {}

    class _Bar:
        def __init__(self):
            self.amount = 0
            self.closed = False

        def update(self, n):
            self.amount += n

        def close(self):
            self.closed = True

    def _fake_get_progress_bar(total, desc):
        captured["bar"] = _Bar()
        return captured["bar"]

    monkeypatch.setattr(dl, "_get_progress_bar", _fake_get_progress_bar)

    dest = tmp_path / "out.bin"
    ok = dl.download_file("http://x", dest, description="X", show_progress=True)
    assert ok is True
    assert dest.read_bytes() == data
    assert captured["bar"].amount == len(data)
    assert captured["bar"].closed is True


def test_download_file_no_progress_zero_length(tmp_path, monkeypatch):
    """total_size==0 and show_progress False -> progress is None (382-386)."""
    data = b"payload"
    monkeypatch.setattr(dl, "open_url", _make_urlopen(data, 0))
    dest = tmp_path / "nop.bin"
    ok = dl.download_file("http://x", dest, show_progress=False)
    assert ok is True
    assert dest.read_bytes() == data


def test_download_file_hash_match(tmp_path, monkeypatch):
    """Matching expected_sha256 passes verification and replaces (404-414)."""
    data = b"verify me"
    digest = hashlib.sha256(data).hexdigest()
    monkeypatch.setattr(dl, "open_url", _make_urlopen(data, len(data)))
    dest = tmp_path / "v.bin"
    ok = dl.download_file("http://x", dest, expected_sha256=digest, show_progress=False)
    assert ok is True
    assert dest.exists()


def test_download_file_hash_mismatch_unlinks_temp(tmp_path, monkeypatch):
    """Mismatch raises ValueError, unlinks temp, leaves no dest (407-411, 418-422)."""
    data = b"verify me"
    monkeypatch.setattr(dl, "open_url", _make_urlopen(data, len(data)))
    dest = tmp_path / "bad.bin"
    with pytest.raises(ValueError, match="Hash mismatch"):
        dl.download_file(
            "http://x", dest, expected_sha256="deadbeef", show_progress=False
        )
    assert not dest.exists()
    # No leftover temp files in the destination directory.
    assert list(tmp_path.glob("*.download")) == []


def test_download_file_oserror_cleans_temp(tmp_path, monkeypatch):
    """An OSError while reading triggers the except cleanup arm (418-422)."""

    class _BoomResponse(_FakeResponse):
        def read(self, size):
            raise OSError("network boom")

    def _boom_urlopen(req, timeout=None, context=None, **kwargs):
        return _BoomResponse(b"", 100)

    monkeypatch.setattr(dl, "open_url", _boom_urlopen)
    dest = tmp_path / "boom.bin"
    with pytest.raises(OSError, match="network boom"):
        dl.download_file("http://x", dest, show_progress=False)
    assert list(tmp_path.glob("*.download")) == []


# ---------------------------------------------------------------------------
# download_planet_centers
# ---------------------------------------------------------------------------


def test_download_planet_centers_already_valid(tmp_path, monkeypatch, capsys):
    """Existing + valid file short-circuits and prints (448-458)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    dest = tmp_path / "planet_centers.bsp"
    dest.write_bytes(b"data")
    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: True)
    monkeypatch.setattr(
        dl,
        "_file_sha256",
        lambda path: dl.DATA_FILES["planet_centers.bsp"]["sha256"],
    )
    result = dl.download_planet_centers(force=False, quiet=False)
    assert result == dest
    out = capsys.readouterr().out
    assert "already exists" in out


def test_download_planet_centers_already_valid_quiet(tmp_path, monkeypatch, capsys):
    """Existing valid file with quiet=True skips the print (arc 455->458)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    dest = tmp_path / "planet_centers.bsp"
    dest.write_bytes(b"data")
    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: True)
    monkeypatch.setattr(
        dl,
        "_file_sha256",
        lambda path: dl.DATA_FILES["planet_centers.bsp"]["sha256"],
    )
    result = dl.download_planet_centers(force=False, quiet=True)
    assert result == dest
    assert capsys.readouterr().out == ""


def test_download_planet_centers_corrupt_then_download(tmp_path, monkeypatch, capsys):
    """Corrupt cache is replaced only once a valid download is verified."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    dest = tmp_path / "planet_centers.bsp"
    dest.write_bytes(b"corrupt")
    # Cached file is corrupt (invalid), the freshly downloaded one is valid.
    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: Path(p).read_bytes() == b"fresh")
    monkeypatch.setattr(dl, "download_file", _fake_download_file(b"fresh"))

    result = dl.download_planet_centers(force=False, quiet=False)
    assert result == dest
    assert dest.read_bytes() == b"fresh"
    out = capsys.readouterr().out
    assert "Downloading planet_centers.bsp" in out
    assert "Downloaded to:" in out


def test_download_planet_centers_keeps_cache_when_download_fails(tmp_path, monkeypatch):
    """A failed re-download leaves the existing artifact byte-for-byte intact.

    This is the data-loss regression: the cached file used to be unlinked
    before the replacement had been fetched, so an offline repair attempt
    ended with no file at all.
    """
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    dest = tmp_path / "planet_centers.bsp"
    dest.write_bytes(b"stale but usable")
    # Cached bytes fail the pinned digest, so the re-download path is taken.
    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: True)
    monkeypatch.setattr(dl, "_file_sha256", lambda path: "0" * 64)

    def _offline(url, dest_path, **kwargs):
        raise OSError("network down")

    monkeypatch.setattr(dl, "download_file", _offline)

    with pytest.raises(OSError, match="network down"):
        dl.download_planet_centers(force=False, quiet=True)

    assert dest.read_bytes() == b"stale but usable"


def test_download_planet_centers_download_invalid(tmp_path, monkeypatch):
    """A structurally invalid download raises and never lands on disk."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: False)
    monkeypatch.setattr(dl, "download_file", _fake_download_file(b"x"))

    with pytest.raises(ValueError, match="failed validation"):
        dl.download_planet_centers(force=True, quiet=True)
    # The corrupt download must not be left behind.
    assert not (tmp_path / "planet_centers.bsp").exists()
    assert list(tmp_path.glob("*.fake-download")) == []


def test_download_planet_centers_invalid_download_keeps_previous(tmp_path, monkeypatch):
    """A corrupt payload never displaces the artifact already installed."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    dest = tmp_path / "planet_centers.bsp"
    dest.write_bytes(b"known good")
    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: False)
    monkeypatch.setattr(dl, "download_file", _fake_download_file(b"truncated"))

    with pytest.raises(ValueError, match="failed validation"):
        dl.download_planet_centers(force=True, quiet=True)

    assert dest.read_bytes() == b"known good"


def test_download_planet_centers_http_404(tmp_path, monkeypatch, capsys):
    """HTTPError 404 prints guidance then re-raises (490-504)."""
    import urllib.error

    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: True)

    def _raise_404(url, dest_path, **kwargs):
        raise urllib.error.HTTPError(url, 404, "Not Found", None, None)

    monkeypatch.setattr(dl, "download_file", _raise_404)
    with pytest.raises(urllib.error.HTTPError):
        dl.download_planet_centers(force=True, quiet=True)
    err = capsys.readouterr().err
    assert "Data file not found" in err


def test_download_planet_centers_http_other(tmp_path, monkeypatch):
    """A non-404 HTTPError re-raises without the 404 branch (490-491, 504)."""
    import urllib.error

    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: True)

    def _raise_500(url, dest_path, **kwargs):
        raise urllib.error.HTTPError(url, 500, "Server Error", None, None)

    monkeypatch.setattr(dl, "download_file", _raise_500)
    with pytest.raises(urllib.error.HTTPError):
        dl.download_planet_centers(force=True, quiet=True)


def test_download_planet_centers_oserror(tmp_path, monkeypatch, capsys):
    """OSError from download_file is reported and re-raised (505-507)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: True)

    def _raise_os(url, dest_path, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(dl, "download_file", _raise_os)
    with pytest.raises(OSError):
        dl.download_planet_centers(force=True, quiet=True)
    err = capsys.readouterr().err
    assert "Error downloading file" in err


# ---------------------------------------------------------------------------
# download_all
# ---------------------------------------------------------------------------


def test_download_all(tmp_path, monkeypatch):
    """download_all delegates to download_planet_centers (525-529)."""
    sentinel = tmp_path / "planet_centers.bsp"
    monkeypatch.setattr(dl, "download_planet_centers", lambda **kwargs: sentinel)
    paths = dl.download_all(force=True, show_progress=False, quiet=True)
    assert paths == [sentinel]


# ---------------------------------------------------------------------------
# check_data_status
# ---------------------------------------------------------------------------


def test_check_data_status(tmp_path, monkeypatch):
    """Covers both the .leb/.leb2 path and the bsp path branches (543-564)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    # Create one bsp file and one leb2 file so both 'exists' arms run.
    (tmp_path / "planet_centers_base.bsp").write_bytes(b"a" * 5)
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    (leb_dir / "base_core.leb2").write_bytes(b"b" * 7)

    status = dl.check_data_status()
    assert status["planet_centers_base.bsp"]["exists"] is True
    assert status["planet_centers_base.bsp"]["size"] == 5
    assert status["base_core.leb2"]["exists"] is True
    assert status["base_core.leb2"]["size"] == 7
    # Reviewed companions are now first-class manifest entries.
    assert status["base_asteroids.leb2"]["exists"] is False


# ---------------------------------------------------------------------------
# print_data_status: exception arms
# ---------------------------------------------------------------------------


def test_print_data_status_exception_arms(tmp_path, monkeypatch, capsys):
    """Force the broad ``except`` arms across print_data_status.

    Hits: version import failure (593-594), calc_mode failure (599-600),
    LEB fallback via _state._LEB_FILE (606-612), DE Loader failure (638-640),
    SPK exception (726-733), ASSIST exception (758-759), IERS exception
    (788-789), and JSON _serialize Path branch (798-800).
    """
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)

    from libephemeris import state

    # calc_mode raises -> except pass arm (599-600).
    monkeypatch.setattr(
        state, "get_calc_mode", lambda: (_ for _ in ()).throw(RuntimeError())
    )
    monkeypatch.setattr(state, "get_precision_tier", lambda: "base")

    # LEB env set to an object whose truthiness check raises, hitting the
    # leb_path except arm (611-612).
    class _BoomEnv(dict):
        def get(self, key, default=None):
            if key == "LIBEPHEMERIS_LEB":
                raise RuntimeError("env boom")
            return super().get(key, default)

    monkeypatch.setattr(dl.os, "environ", _BoomEnv())

    # Force the skyfield Loader import/usage to raise (DE kernel 639-640).
    fake_sf_api = types.ModuleType("skyfield.api")

    def _boom_loader(*a, **k):
        raise RuntimeError("loader boom")

    fake_sf_api.Loader = _boom_loader
    monkeypatch.setitem(sys.modules, "skyfield.api", fake_sf_api)

    # Accessing any attribute of these fake modules raises, exercising the
    # SPK (732) and IERS (788-789) broad except arms.
    class _BoomModule(types.ModuleType):
        def __getattr__(self, name):
            raise RuntimeError("module boom")

    monkeypatch.setitem(
        sys.modules, "libephemeris.spk_auto", _BoomModule("libephemeris.spk_auto")
    )
    monkeypatch.setitem(
        sys.modules, "libephemeris.iers_data", _BoomModule("libephemeris.iers_data")
    )

    # Force the ASSIST block to raise (758-759).
    monkeypatch.setattr(
        dl.os.path,
        "expanduser",
        lambda p: (_ for _ in ()).throw(RuntimeError("no home")),
    )

    # Force version import failure (593-594): remove __version__ so
    # "from . import __version__" raises.
    import libephemeris as _pkg

    monkeypatch.delattr(_pkg, "__version__", raising=False)

    dl.print_data_status(as_json=True, verbose=0)
    out = capsys.readouterr().out
    assert "unknown" in out  # version fell back
    assert '"data_directory"' in out


def test_print_data_status_json_serialize_path(tmp_path, monkeypatch, capsys):
    """A Path value in the result triggers the JSON _serialize arm (797-800)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)

    from libephemeris import state

    monkeypatch.setattr(state, "get_calc_mode", lambda: "leb")
    monkeypatch.setattr(state, "get_precision_tier", lambda: "base")
    monkeypatch.delenv("LIBEPHEMERIS_LEB", raising=False)

    from libephemeris import iers_data

    # finals_age_days flows straight into result without str(); a Path here
    # is not JSON-serializable, so the default=_serialize hook converts it.
    monkeypatch.setattr(
        iers_data,
        "get_iers_cache_info",
        lambda: {
            "cache_dir": str(tmp_path / "iers"),
            "finals_exists": True,
            "finals_age_days": Path("/some/age"),
            "leap_seconds_exists": False,
            "leap_seconds_age_days": None,
            "delta_t_exists": False,
            "delta_t_age_days": None,
        },
    )

    dl.print_data_status(as_json=True, verbose=0)
    out = capsys.readouterr().out
    assert "/some/age" in out


def test_print_data_status_json_serialize_non_path(tmp_path, monkeypatch):
    """A non-Path non-serializable value reaches the `return obj` arm (800).

    _serialize returns the object unchanged; json then re-encounters it and
    raises (TypeError, or ValueError for the circular-reference guard).
    """
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)

    from libephemeris import state

    monkeypatch.setattr(state, "get_calc_mode", lambda: "leb")
    monkeypatch.setattr(state, "get_precision_tier", lambda: "base")
    monkeypatch.delenv("LIBEPHEMERIS_LEB", raising=False)

    from libephemeris import iers_data

    # A set is neither a Path nor JSON-serializable -> _serialize hits 800.
    monkeypatch.setattr(
        iers_data,
        "get_iers_cache_info",
        lambda: {
            "cache_dir": str(tmp_path / "iers"),
            "finals_exists": True,
            "finals_age_days": {1, 2, 3},
            "leap_seconds_exists": False,
            "leap_seconds_age_days": None,
            "delta_t_exists": False,
            "delta_t_age_days": None,
        },
    )

    with pytest.raises((TypeError, ValueError)):
        dl.print_data_status(as_json=True, verbose=0)


def test_print_data_status_de_loader_found(tmp_path, monkeypatch, capsys):
    """A skyfield Loader directory containing the DE kernel hits line 638."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)

    from libephemeris import state

    monkeypatch.setattr(state, "get_calc_mode", lambda: "leb")
    monkeypatch.setattr(state, "get_precision_tier", lambda: "base")
    monkeypatch.delenv("LIBEPHEMERIS_LEB", raising=False)

    # Place DE kernels in a separate "skyfield" dir, not in data_dir, so the
    # `path.exists()` check on data_dir fails and the Loader branch runs.
    sf_dir = tmp_path / "sf"
    sf_dir.mkdir()
    for fn in ("de440s.bsp", "de440.bsp", "de441.bsp"):
        (sf_dir / fn).write_bytes(b"k")

    fake_sf_api = types.ModuleType("skyfield.api")

    class _FakeLoader:
        def __init__(self, directory):
            self.directory = str(sf_dir)

    fake_sf_api.Loader = _FakeLoader
    monkeypatch.setitem(sys.modules, "skyfield.api", fake_sf_api)

    dl.print_data_status(as_json=False, verbose=1)
    out = capsys.readouterr().out
    assert "de440s.bsp" in out


def test_print_data_status_spk_dir_missing(tmp_path, monkeypatch, capsys):
    """SPK cache directory absent hits the else arm (726-731, 909)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)

    from libephemeris import state

    monkeypatch.setattr(state, "get_calc_mode", lambda: "leb")
    monkeypatch.setattr(state, "get_precision_tier", lambda: "base")
    monkeypatch.delenv("LIBEPHEMERIS_LEB", raising=False)

    # Point the SPK cache at a non-existent directory -> else arm + 0 files.
    monkeypatch.setattr(
        state, "get_spk_cache_dir", lambda: str(tmp_path / "no_such_spk")
    )

    dl.print_data_status(as_json=False, verbose=0)
    out = capsys.readouterr().out
    assert "Files:      (none)" in out


def test_print_data_status_text_output(tmp_path, monkeypatch, capsys):
    """Drive the formatted-text path with present + missing files.

    Covers the LEB2 partial-present branch (884), SPK files-present branch,
    ASSIST present + missing arms (921/923), IERS present/missing branches
    (943), and verbose>=2 SPK file listing.
    """
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)

    from libephemeris import state

    monkeypatch.setattr(state, "get_calc_mode", lambda: "leb")
    monkeypatch.setattr(state, "get_precision_tier", lambda: "base")
    monkeypatch.setenv("LIBEPHEMERIS_LEB", str(tmp_path / "leb" / "ephemeris_base.leb"))

    # Create a partial LEB2 set for "base" (present>0 but not all) -> line 884.
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    (leb_dir / "base_core.leb2").write_bytes(b"x" * 9)
    (leb_dir / "base_asteroids.leb2").write_bytes(b"x" * 9)
    # Leave apogee missing so present != total.

    # Provide an SPK cache directory with one file -> files>0 branch.
    spk_dir = tmp_path / "spk"
    spk_dir.mkdir()
    (spk_dir / "one.bsp").write_bytes(b"y" * 11)
    monkeypatch.setattr(state, "get_spk_cache_dir", lambda: str(spk_dir))

    # Provide ASSIST files: one present, one missing -> both arms (921/923).
    fake_assist = tmp_path / "assist"
    fake_assist.mkdir()
    (fake_assist / "linux_p1550p2650.440").write_bytes(b"z" * 13)
    # sb441-n16.bsp deliberately missing.

    real_expanduser = dl.os.path.expanduser

    def _fake_expanduser(p):
        if p == "~/.libephemeris/assist":
            return str(fake_assist)
        return real_expanduser(p)

    monkeypatch.setattr(dl.os.path, "expanduser", _fake_expanduser)

    # Provide IERS info: finals present with age (943), leap_seconds missing.
    from libephemeris import iers_data

    monkeypatch.setattr(
        iers_data,
        "get_iers_cache_info",
        lambda: {
            "cache_dir": str(tmp_path / "iers"),
            "finals_exists": True,
            "finals_age_days": 2.0,
            "leap_seconds_exists": False,
            "leap_seconds_age_days": None,
            "delta_t_exists": True,
            "delta_t_age_days": None,
        },
    )

    dl.print_data_status(as_json=False, verbose=2)
    out = capsys.readouterr().out
    assert "libephemeris" in out
    assert "base_core.leb2" in out
    assert "SPK Asteroid Cache" in out
    assert "ASSIST N-body Data" in out
    assert "IERS Earth Orientation Data" in out


def test_print_data_status_assist_and_iers_empty(tmp_path, monkeypatch, capsys):
    """assist_info / iers_info empty -> '(could not check)' arms (926, 946)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)

    from libephemeris import state

    monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
    monkeypatch.setattr(state, "get_precision_tier", lambda: "base")
    monkeypatch.delenv("LIBEPHEMERIS_LEB", raising=False)

    # Make the ASSIST block raise so assist_info stays empty.
    def _raise_expanduser(p):
        if p == "~/.libephemeris/assist":
            raise RuntimeError("no home")
        return p

    monkeypatch.setattr(dl.os.path, "expanduser", _raise_expanduser)

    # Make the IERS block raise so iers_info stays empty.
    from libephemeris import iers_data

    monkeypatch.setattr(
        iers_data,
        "get_iers_cache_info",
        lambda: (_ for _ in ()).throw(RuntimeError("no iers")),
    )

    dl.print_data_status(as_json=False, verbose=0)
    out = capsys.readouterr().out
    assert "(could not check)" in out


# ---------------------------------------------------------------------------
# init_all
# ---------------------------------------------------------------------------


def test_init_all_full_paths(monkeypatch, capsys):
    """Drive init_all through success + every exception arm (996-1212).

    A small year range keeps the loop short. The fake
    download_spk_from_horizons raises a different error on each call to
    cover: success, START-outside-limits, STOP-outside-limits, and a
    generic failure. A pre-existing cached file covers the skip arm.
    """
    import libephemeris.spk_auto as spk_auto
    from libephemeris import state

    # Step 1: DE440 succeeds.
    monkeypatch.setattr(state, "get_planets", lambda: object())
    # Step 2: planet_centers succeeds.
    monkeypatch.setattr(dl, "download_planet_centers", lambda **kw: Path("/tmp/pc.bsp"))

    # No real sleeping for rate-limit (init_all imports time locally).
    import time as _time

    monkeypatch.setattr(_time, "sleep", lambda s: None)

    monkeypatch.setattr(spk_auto, "ensure_cache_dir", lambda c: "/tmp/spkcache")
    monkeypatch.setattr(
        spk_auto, "_generate_spk_cache_filename", lambda hid, a, b: f"{hid}_{a}_{b}.bsp"
    )
    monkeypatch.setattr(spk_auto, "_iso_to_jd", lambda s: 2451545.0)

    # First chunk: pretend file already exists (skip). Then download path.
    existing = {"first": True}

    def _fake_exists(path):
        if existing["first"]:
            existing["first"] = False
            return True
        return False

    monkeypatch.setattr(dl.os.path, "exists", _fake_exists)

    call_state = {"n": 0}

    def _fake_download_spk(body_id, jd_start, jd_end, output_path, ipl):
        call_state["n"] += 1
        n = call_state["n"]
        if n == 1:
            return None  # success
        if n == 2:
            raise ValueError("START time outside set SPK limits")
        if n == 3:
            raise ValueError("STOP time outside set SPK limits")
        raise ValueError("generic horizons failure")

    monkeypatch.setattr(spk_auto, "download_spk_from_horizons", _fake_download_spk)

    result = dl.init_all(
        cache_dir=None,
        force=False,
        show_progress=True,
        quiet=False,
        start_year=2000,
        end_year=2100,
    )
    out = capsys.readouterr().out
    assert result["de440"] is True
    assert result["planet_centers"] is True
    assert result["spk_skipped"] >= 1
    assert result["spk_success"] >= 1
    assert "Step 1/3" in out
    assert "Step 3/3" in out


def test_init_all_step_failures(monkeypatch, capsys):
    """DE440 + planet_centers failures hit the [FAIL] print arms (1026, 1046)."""
    import libephemeris.spk_auto as spk_auto
    from libephemeris import state

    monkeypatch.setattr(
        state, "get_planets", lambda: (_ for _ in ()).throw(RuntimeError("no de"))
    )
    monkeypatch.setattr(
        dl,
        "download_planet_centers",
        lambda **kw: (_ for _ in ()).throw(OSError("no pc")),
    )

    import time as _time

    monkeypatch.setattr(_time, "sleep", lambda s: None)
    monkeypatch.setattr(spk_auto, "ensure_cache_dir", lambda c: "/tmp/spk2")
    monkeypatch.setattr(
        spk_auto, "_generate_spk_cache_filename", lambda hid, a, b: f"{hid}.bsp"
    )
    monkeypatch.setattr(spk_auto, "_iso_to_jd", lambda s: 2451545.0)
    # All chunks already cached so the SPK loop short-circuits (skip arm).
    monkeypatch.setattr(dl.os.path, "exists", lambda p: True)

    result = dl.init_all(
        force=False,
        quiet=False,
        start_year=2000,
        end_year=2020,
    )
    assert result["de440"] is False
    assert result["planet_centers"] is False
    err = capsys.readouterr().err
    assert "[FAIL] DE440.bsp" in err
    assert "[FAIL] planet_centers.bsp" in err


def test_init_all_quiet_full_loop(monkeypatch):
    """quiet=True path through DE/PC fail, blocked body, and the chunk loop.

    Covers the quiet arcs 1025->1031, 1045->1051, 1089->1094 (blocked body
    quiet), and 1198->1084 (skip the per-body summary, go to next body).
    """
    import libephemeris.spk_auto as spk_auto
    from libephemeris import state

    monkeypatch.setattr(
        state, "get_planets", lambda: (_ for _ in ()).throw(RuntimeError("no de"))
    )
    monkeypatch.setattr(
        dl,
        "download_planet_centers",
        lambda **kw: (_ for _ in ()).throw(OSError("no pc")),
    )

    import time as _time

    monkeypatch.setattr(_time, "sleep", lambda s: None)
    monkeypatch.setattr(spk_auto, "ensure_cache_dir", lambda c: "/tmp/spk5")
    monkeypatch.setattr(
        spk_auto, "_generate_spk_cache_filename", lambda hid, a, b: f"{hid}.bsp"
    )
    monkeypatch.setattr(spk_auto, "_iso_to_jd", lambda s: 2451545.0)
    # Nothing cached -> each chunk attempts a download.
    monkeypatch.setattr(dl.os.path, "exists", lambda p: False)
    monkeypatch.setattr(spk_auto, "download_spk_from_horizons", lambda **kw: None)

    result = dl.init_all(
        force=False,
        quiet=True,
        start_year=2000,
        end_year=2020,
    )
    assert result["de440"] is False
    assert result["planet_centers"] is False
    # A blocked body (e.g. Bennu) was recorded.
    assert "spk_blocked" in result


def test_init_all_keyboard_interrupt(monkeypatch):
    """KeyboardInterrupt during SPK download re-raises (1150-1153)."""
    import libephemeris.spk_auto as spk_auto
    from libephemeris import state

    monkeypatch.setattr(state, "get_planets", lambda: object())
    monkeypatch.setattr(dl, "download_planet_centers", lambda **kw: Path("/tmp/x"))

    import time as _time

    monkeypatch.setattr(_time, "sleep", lambda s: None)
    monkeypatch.setattr(spk_auto, "ensure_cache_dir", lambda c: "/tmp/spk3")
    monkeypatch.setattr(
        spk_auto, "_generate_spk_cache_filename", lambda hid, a, b: f"{hid}.bsp"
    )
    monkeypatch.setattr(spk_auto, "_iso_to_jd", lambda s: 2451545.0)
    monkeypatch.setattr(dl.os.path, "exists", lambda p: False)

    def _interrupt(**kw):
        raise KeyboardInterrupt()

    monkeypatch.setattr(spk_auto, "download_spk_from_horizons", _interrupt)

    with pytest.raises(KeyboardInterrupt):
        dl.init_all(quiet=False, start_year=2000, end_year=2020)


def test_init_all_keyboard_interrupt_quiet(monkeypatch):
    """KeyboardInterrupt with quiet=True skips the stdout write arc (1151->1153)."""
    import libephemeris.spk_auto as spk_auto
    from libephemeris import state

    monkeypatch.setattr(state, "get_planets", lambda: object())
    monkeypatch.setattr(dl, "download_planet_centers", lambda **kw: Path("/tmp/x"))

    import time as _time

    monkeypatch.setattr(_time, "sleep", lambda s: None)
    monkeypatch.setattr(spk_auto, "ensure_cache_dir", lambda c: "/tmp/spk4")
    monkeypatch.setattr(
        spk_auto, "_generate_spk_cache_filename", lambda hid, a, b: f"{hid}.bsp"
    )
    monkeypatch.setattr(spk_auto, "_iso_to_jd", lambda s: 2451545.0)
    monkeypatch.setattr(dl.os.path, "exists", lambda p: False)
    monkeypatch.setattr(
        spk_auto,
        "download_spk_from_horizons",
        lambda **kw: (_ for _ in ()).throw(KeyboardInterrupt()),
    )

    with pytest.raises(KeyboardInterrupt):
        dl.init_all(quiet=True, start_year=2000, end_year=2020)


# ---------------------------------------------------------------------------
# download_for_tier
# ---------------------------------------------------------------------------


def test_download_for_tier_invalid():
    """Unknown tier raises ValueError (1246-1248)."""
    with pytest.raises(ValueError, match="Unknown tier"):
        dl.download_for_tier("nonsense")


def test_download_for_tier_success(monkeypatch, capsys):
    """Full happy path through both steps and the summary print (1244-1308)."""
    import libephemeris.spk_auto as spk_auto
    from libephemeris import state

    monkeypatch.setattr(state, "set_precision_tier", lambda t: None)
    monkeypatch.setattr(
        dl, "_download_planet_centers_for_tier", lambda **kw: Path("/tmp/pc")
    )
    monkeypatch.setattr(
        spk_auto,
        "ensure_all_ephemerides",
        lambda **kw: {
            "summary": {
                "cached": 3,
                "downloaded": 2,
                "fallback": 1,
                "errors": 1,
            }
        },
    )
    result = dl.download_for_tier("base", force=True, quiet=False)
    out = capsys.readouterr().out
    assert "summary" in result
    assert "Download complete:" in out
    assert "SPK fallback:" in out
    assert "SPK errors:" in out


def test_download_for_tier_pc_warning(monkeypatch, capsys):
    """planet_centers step failure hits the WARN arm (1276-1280)."""
    import libephemeris.spk_auto as spk_auto
    from libephemeris import state

    monkeypatch.setattr(state, "set_precision_tier", lambda t: None)
    monkeypatch.setattr(
        dl,
        "_download_planet_centers_for_tier",
        lambda **kw: (_ for _ in ()).throw(ValueError("pc fail")),
    )
    monkeypatch.setattr(
        spk_auto, "ensure_all_ephemerides", lambda **kw: {"summary": {}}
    )
    result = dl.download_for_tier("medium", quiet=False)
    err = capsys.readouterr().err
    assert "[WARN]" in err
    assert isinstance(result, dict)


def test_download_for_tier_non_dict_summary(monkeypatch, capsys):
    """quiet=False + non-dict summary skips the dict block (arc 1300->1308)."""
    import libephemeris.spk_auto as spk_auto
    from libephemeris import state

    monkeypatch.setattr(state, "set_precision_tier", lambda t: None)
    monkeypatch.setattr(
        dl, "_download_planet_centers_for_tier", lambda **kw: Path("/tmp/pc")
    )
    monkeypatch.setattr(
        spk_auto, "ensure_all_ephemerides", lambda **kw: {"summary": "not-a-dict"}
    )
    result = dl.download_for_tier("extended", quiet=False)
    assert "summary" in result
    out = capsys.readouterr().out
    assert "Download complete:" in out
    # The dict-only lines must not appear for a non-dict summary.
    assert "SPK cached:" not in out


def test_download_for_tier_pc_warning_quiet(monkeypatch):
    """planet_centers failure with quiet=True skips the WARN print (arc 1277->1283)."""
    import libephemeris.spk_auto as spk_auto
    from libephemeris import state

    monkeypatch.setattr(state, "set_precision_tier", lambda t: None)
    monkeypatch.setattr(
        dl,
        "_download_planet_centers_for_tier",
        lambda **kw: (_ for _ in ()).throw(ValueError("pc fail")),
    )
    monkeypatch.setattr(
        spk_auto, "ensure_all_ephemerides", lambda **kw: {"summary": {}}
    )
    result = dl.download_for_tier("base", quiet=True)
    assert isinstance(result, dict)


def test_download_for_tier_success_quiet(monkeypatch):
    """planet_centers success with quiet=True skips the [OK] print (arc 1273->1283)."""
    import libephemeris.spk_auto as spk_auto
    from libephemeris import state

    monkeypatch.setattr(state, "set_precision_tier", lambda t: None)
    monkeypatch.setattr(
        dl, "_download_planet_centers_for_tier", lambda **kw: Path("/tmp/pc")
    )
    monkeypatch.setattr(
        spk_auto, "ensure_all_ephemerides", lambda **kw: {"summary": {}}
    )
    result = dl.download_for_tier("base", quiet=True)
    assert isinstance(result, dict)


# ---------------------------------------------------------------------------
# _is_valid_leb / get_leb_dir / get_leb_path_for_tier
# ---------------------------------------------------------------------------


def test_is_valid_leb_true(monkeypatch):
    """open_leb succeeding returns True (1322-1327)."""
    import libephemeris.leb_reader as leb_reader

    class _Reader:
        def close(self):
            pass

    monkeypatch.setattr(leb_reader, "open_leb", lambda p: _Reader())
    assert dl._is_valid_leb("/tmp/x.leb") is True


def test_is_valid_leb_false(monkeypatch):
    """open_leb raising returns False (1328-1329)."""
    import libephemeris.leb_reader as leb_reader

    def _raise(p):
        raise ValueError("bad leb")

    monkeypatch.setattr(leb_reader, "open_leb", _raise)
    assert dl._is_valid_leb("/tmp/x.leb") is False


def test_get_leb_dir_and_path(tmp_path, monkeypatch):
    """get_leb_dir creates dir (1341-1343); get_leb_path_for_tier joins (1355)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    leb_dir = dl.get_leb_dir()
    assert leb_dir == tmp_path / "leb"
    assert leb_dir.exists()
    path = dl.get_leb_path_for_tier("base")
    assert path == tmp_path / "leb" / "ephemeris_base.leb"


# ---------------------------------------------------------------------------
# download_leb_for_tier
# ---------------------------------------------------------------------------


def test_download_leb_for_tier_invalid_tier():
    """Unknown tier raises ValueError (1392-1396)."""
    with pytest.raises(ValueError, match="Unknown tier"):
        dl.download_leb_for_tier("nonsense")


@pytest.mark.parametrize("tier", ["base", "medium", "extended"])
def test_download_leb_for_tier_delegates_to_reviewed_leb2_core(
    tier, tmp_path, monkeypatch
):
    """The compatibility entry point installs the current reviewed format."""
    core = tmp_path / f"{tier}_core.leb2"
    calls = []

    def _download(tier_name, **kwargs):
        calls.append((tier_name, kwargs))
        return [core]

    monkeypatch.setattr(dl, "download_leb2_for_tier", _download)
    result = dl.download_leb_for_tier(
        tier, force=True, show_progress=False, quiet=True, activate=False
    )

    assert result == core
    assert calls == [
        (
            tier,
            {
                "groups": ["core"],
                "force": True,
                "show_progress": False,
                "quiet": True,
                "activate": False,
            },
        )
    ]


def test_download_leb_for_tier_requires_published_core(monkeypatch):
    """A missing manifest entry remains fail-closed."""
    monkeypatch.setattr(dl, "download_leb2_for_tier", lambda *args, **kwargs: [])
    with pytest.raises(RuntimeError, match="SHA-256-pinned"):
        dl.download_leb_for_tier("medium")


@pytest.mark.skip(reason="legacy monolithic release downloads are retired")
def test_download_leb_for_tier_no_hash(monkeypatch):
    """A None sha256 raises RuntimeError (1399-1403)."""
    monkeypatch.setitem(
        dl.DATA_FILES,
        "ephemeris_base.leb",
        {**dl.DATA_FILES["ephemeris_base.leb"], "sha256": None},
    )
    with pytest.raises(RuntimeError, match="no download hash"):
        dl.download_leb_for_tier("base")


@pytest.mark.skip(reason="legacy monolithic release downloads are retired")
def test_download_leb_for_tier_existing_valid_activate(tmp_path, monkeypatch, capsys):
    """Existing valid file short-circuits and activates (1408-1418)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    dest = leb_dir / "ephemeris_base.leb"
    dest.write_bytes(b"leb")
    monkeypatch.setattr(dl, "_is_valid_leb", lambda p: True)

    from libephemeris import state

    activated = {}
    monkeypatch.setattr(state, "set_leb_file", lambda p: activated.update(path=p))

    result = dl.download_leb_for_tier("base", force=False, quiet=False, activate=True)
    assert result == dest
    assert activated["path"] == str(dest)
    out = capsys.readouterr().out
    assert "already exists" in out


@pytest.mark.skip(reason="legacy monolithic release downloads are retired")
def test_download_leb_for_tier_existing_valid_quiet_no_activate(tmp_path, monkeypatch):
    """Existing valid file, quiet=True + activate=False (arcs 1410->1414, 1414->1418)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    dest = leb_dir / "ephemeris_base.leb"
    dest.write_bytes(b"leb")
    monkeypatch.setattr(dl, "_is_valid_leb", lambda p: True)

    from libephemeris import state

    called = {"n": 0}
    monkeypatch.setattr(
        state, "set_leb_file", lambda p: called.update(n=called["n"] + 1)
    )

    result = dl.download_leb_for_tier("base", force=False, quiet=True, activate=False)
    assert result == dest
    # activate=False -> set_leb_file must not be called.
    assert called["n"] == 0


@pytest.mark.skip(reason="legacy monolithic release downloads are retired")
def test_download_leb_for_tier_corrupt_then_download(tmp_path, monkeypatch, capsys):
    """Corrupt cache removed, fresh download validates + activates (1419-1467)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    dest = leb_dir / "ephemeris_base.leb"
    dest.write_bytes(b"corrupt")

    valid_states = {"calls": 0}

    def _is_valid(p):
        # First call (cache check) -> False; subsequent (post-download) -> True.
        valid_states["calls"] += 1
        return valid_states["calls"] > 1

    monkeypatch.setattr(dl, "_is_valid_leb", _is_valid)

    def _fake_download_file(url, dest_path, **kw):
        Path(dest_path).write_bytes(b"fresh-leb")
        return True

    monkeypatch.setattr(dl, "download_file", _fake_download_file)

    from libephemeris import state

    activated = {}
    monkeypatch.setattr(state, "set_leb_file", lambda p: activated.update(path=p))

    result = dl.download_leb_for_tier("base", force=False, quiet=False, activate=True)
    assert result == dest
    assert activated["path"] == str(dest)
    out = capsys.readouterr().out
    assert "Downloaded to:" in out


@pytest.mark.skip(reason="legacy monolithic release downloads are retired")
def test_download_leb_for_tier_corrupt_remove_oserror(tmp_path, monkeypatch):
    """os.remove OSError during corrupt-cache cleanup is swallowed (1422-1424)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    dest = leb_dir / "ephemeris_base.leb"
    dest.write_bytes(b"corrupt")

    valid_states = {"calls": 0}

    def _is_valid(p):
        valid_states["calls"] += 1
        return valid_states["calls"] > 1

    monkeypatch.setattr(dl, "_is_valid_leb", _is_valid)
    monkeypatch.setattr(
        dl.os, "remove", lambda p: (_ for _ in ()).throw(OSError("nope"))
    )
    monkeypatch.setattr(
        dl,
        "download_file",
        lambda url, dest_path, **kw: Path(dest_path).write_bytes(b"f"),
    )
    result = dl.download_leb_for_tier("base", force=False, quiet=True, activate=False)
    assert result == dest


@pytest.mark.skip(reason="legacy monolithic release downloads are retired")
def test_download_leb_for_tier_download_invalid(tmp_path, monkeypatch):
    """Downloaded file failing validation raises + removes (1444-1451)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    # No pre-existing file -> goes straight to download.
    monkeypatch.setattr(dl, "_is_valid_leb", lambda p: False)
    monkeypatch.setattr(
        dl,
        "download_file",
        lambda url, dest_path, **kw: Path(dest_path).write_bytes(b"x"),
    )
    with pytest.raises(ValueError, match="failed LEB validation"):
        dl.download_leb_for_tier("base", force=True, quiet=False, activate=False)


@pytest.mark.skip(reason="legacy monolithic release downloads are retired")
def test_download_leb_for_tier_download_invalid_remove_oserror(tmp_path, monkeypatch):
    """Validation-failure os.remove OSError swallowed before raise (1446-1448)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    monkeypatch.setattr(dl, "_is_valid_leb", lambda p: False)
    monkeypatch.setattr(
        dl,
        "download_file",
        lambda url, dest_path, **kw: Path(dest_path).write_bytes(b"x"),
    )
    monkeypatch.setattr(
        dl.os, "remove", lambda p: (_ for _ in ()).throw(OSError("nope"))
    )
    with pytest.raises(ValueError, match="failed LEB validation"):
        dl.download_leb_for_tier("base", force=True, quiet=True, activate=False)


@pytest.mark.skip(reason="legacy monolithic release downloads are retired")
def test_download_leb_for_tier_download_no_activate(tmp_path, monkeypatch, capsys):
    """Successful download with activate=False skips set_leb_file (1453-1467)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    monkeypatch.setattr(dl, "_is_valid_leb", lambda p: True)
    monkeypatch.setattr(
        dl,
        "download_file",
        lambda url, dest_path, **kw: Path(dest_path).write_bytes(b"x"),
    )
    result = dl.download_leb_for_tier("base", force=True, quiet=False, activate=False)
    assert result == leb_dir / "ephemeris_base.leb"
    out = capsys.readouterr().out
    assert "Downloaded to:" in out


# ---------------------------------------------------------------------------
# download_leb2_for_tier
# ---------------------------------------------------------------------------


def test_download_leb2_for_tier_invalid_tier():
    with pytest.raises(ValueError, match="Unknown tier"):
        dl.download_leb2_for_tier("nonsense")


def test_published_wider_core_pins_match_reviewed_release_bytes():
    # Cumulative data-v3 cores (reviewed regeneration).
    assert dl.DATA_FILES["medium_core.leb2"]["sha256"] == (
        "4d88ec9a79add7e3af9e75ac3ceabe5462a4af440447578a5ccba69a3e0a55b6"
    )
    assert dl.DATA_FILES["extended_core.leb2"]["sha256"] == (
        "ecb8dd43a4934e74324d9c38274a72360ef6821d96607679ef5248d47c9d7afc"
    )
    assert dl.DATA_FILES["medium_core.leb2"]["url"].endswith("data-v3/medium_core.leb2")
    assert dl.DATA_FILES["extended_core.leb2"]["url"].endswith(
        "data-v3/extended_core.leb2"
    )


def test_leb_and_planet_center_release_families_stay_separate():
    """New LEB bytes must not silently repoint unchanged legacy BSP assets."""
    for filename, metadata in dl.DATA_FILES.items():
        url = metadata.get("url")
        if filename.endswith(".leb2") and url is not None:
            assert "/data-v3/" in url
        if filename.startswith("planet_centers") and filename.endswith(".bsp"):
            assert "/data-v2rc1/" in url


def test_regenerated_medium_planet_centers_pin_matches_reviewed_bytes():
    assert dl.DATA_FILES["planet_centers_medium.bsp"]["sha256"] == (
        "d3c34f5efe9223ef588ec59a8c59c1bd6619b0eab5d5e0b35c353d675efe7b4d"
    )
    assert dl.DATA_FILES["planet_centers_medium.bsp"]["size_mb"] == 191.25


def test_install_bundled_file_verifies_and_publishes(tmp_path, monkeypatch):
    """A package resource is hash-checked before atomic publication."""
    package = tmp_path / "package"
    source = package / "data" / "core.leb2"
    source.parent.mkdir(parents=True)
    payload = b"reviewed bundled core"
    source.write_bytes(payload)
    monkeypatch.setattr(dl.resources, "files", lambda package_name: package)

    destination = tmp_path / "cache" / "base_core.leb2"
    dl._install_bundled_file(
        "data/core.leb2",
        destination,
        expected_sha256=hashlib.sha256(payload).hexdigest(),
    )

    assert destination.read_bytes() == payload


def test_install_bundled_file_rejects_hash_mismatch(tmp_path, monkeypatch):
    """A corrupt package resource never replaces an existing destination."""
    package = tmp_path / "package"
    source = package / "data" / "core.leb2"
    source.parent.mkdir(parents=True)
    source.write_bytes(b"unexpected")
    monkeypatch.setattr(dl.resources, "files", lambda package_name: package)
    destination = tmp_path / "base_core.leb2"
    destination.write_bytes(b"known-good")

    with pytest.raises(ValueError, match="Bundled data hash mismatch"):
        dl._install_bundled_file(
            "data/core.leb2",
            destination,
            expected_sha256=hashlib.sha256(b"expected").hexdigest(),
        )

    assert destination.read_bytes() == b"known-good"


def test_download_leb2_for_tier_full(tmp_path, monkeypatch, capsys):
    """Every reviewed companion is installed and the core is activated."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()

    # A stale companion with the right name is replaced because its hash fails.
    (leb_dir / "base_asteroids.leb2").write_bytes(b"ast")

    def _is_valid(p):
        return True

    monkeypatch.setattr(dl, "_is_valid_leb", _is_valid)

    monkeypatch.setattr(
        dl,
        "download_file",
        lambda url, dest_path, **kwargs: Path(dest_path).write_bytes(b"reviewed"),
    )

    from libephemeris import state

    activated = {}
    monkeypatch.setattr(state, "set_leb_file", lambda p: activated.update(path=p))

    # Include a bogus group to hit the "not in DATA_FILES" SKIP arm.
    groups = ["core", "asteroids", "exotics", "apogee", "bogus"]
    downloaded = dl.download_leb2_for_tier(
        "base", groups=groups, force=False, quiet=False, activate=True
    )
    out = capsys.readouterr().out
    assert "[SKIP]" in out
    assert "already exists" not in out
    assert "[FAIL]" not in out
    assert "base_bogus.leb2: not in DATA_FILES" in out
    # core was downloaded -> activation should have happened.
    assert activated["path"].endswith("base_core.leb2")
    assert any(str(p).endswith("base_core.leb2") for p in downloaded)
    assert {p.name for p in downloaded} == {
        "base_core.leb2",
        "base_asteroids.leb2",
        "base_exotics.leb2",
        "base_apogee.leb2",
    }


def test_download_leb2_remote_asset_requires_and_passes_sha_pin(tmp_path, monkeypatch):
    """A reviewed wider-tier manifest entry uses the pinned downloader."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    payload = b"synthetic structurally valid reviewed core"
    digest = hashlib.sha256(payload).hexdigest()
    monkeypatch.setitem(
        dl.DATA_FILES,
        "medium_core.leb2",
        {
            "url": "https://example.invalid/medium_core.leb2",
            "sha256": digest,
            "size_mb": 1,
            "description": "test",
            "dest_subdir": "leb",
        },
    )
    calls = []

    def _download(url, dest_path, **kwargs):
        calls.append((url, kwargs))
        Path(dest_path).write_bytes(payload)
        return True

    monkeypatch.setattr(dl, "download_file", _download)
    monkeypatch.setattr(dl, "_is_valid_leb", lambda path: True)

    result = dl.download_leb2_for_tier(
        "medium", groups=["core"], quiet=True, activate=False
    )

    assert result == [tmp_path / "leb" / "medium_core.leb2"]
    assert len(calls) == 1
    url, kwargs = calls[0]
    assert url == "https://example.invalid/medium_core.leb2"
    assert kwargs["description"] == "medium_core.leb2"
    assert kwargs["show_progress"] is True
    assert kwargs["expected_sha256"] == digest
    # The structural check is handed to the downloader so it runs on the
    # temp file, before the atomic replace.
    assert kwargs["validator"] is dl._is_valid_leb


def test_download_leb2_default_groups_exclude_uranians() -> None:
    """The retired uranians group never re-enters the tier download (3.1.0)."""
    assert "uranians" not in dl.LEB2_GROUPS
    assert not any("uranians" in name for name in dl.DATA_FILES)


def test_download_leb2_for_tier_default_groups_quiet(tmp_path, monkeypatch):
    """groups=None defaults to LEB2_GROUPS; quiet + no activation (1498-1549)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    monkeypatch.setattr(dl, "_is_valid_leb", lambda p: True)
    monkeypatch.setattr(
        dl,
        "download_file",
        lambda url, dest_path, **kw: Path(dest_path).write_bytes(b"d"),
    )

    from libephemeris import state

    monkeypatch.setattr(state, "set_leb_file", lambda p: None)

    downloaded = dl.download_leb2_for_tier(
        "base", groups=None, force=True, quiet=True, activate=False
    )
    # Every canonical group has an independently reviewed distribution asset.
    assert sorted(p.name for p in downloaded) == [
        "base_apogee.leb2",
        "base_asteroids.leb2",
        "base_core.leb2",
        "base_exotics.leb2",
    ]


def test_download_leb2_for_tier_quiet_arcs(tmp_path, monkeypatch):
    """quiet mode installs reviewed companions and skips unknown groups.

    - bogus group skip (quiet)
    - pre-existing hash-mismatched file replaced
    - reviewed core and requested companions installed
    """
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()

    # A cached same-name file must not be trusted without the pinned hash.
    (leb_dir / "base_asteroids.leb2").write_bytes(b"ast")
    monkeypatch.setattr(dl, "_is_valid_leb", lambda p: True)

    monkeypatch.setattr(
        dl,
        "download_file",
        lambda url, dest_path, **kwargs: Path(dest_path).write_bytes(b"reviewed"),
    )

    from libephemeris import state

    monkeypatch.setattr(state, "set_leb_file", lambda p: None)

    groups = ["core", "asteroids", "apogee", "bogus"]
    downloaded = dl.download_leb2_for_tier(
        "base", groups=groups, force=False, quiet=True, activate=True
    )
    assert sorted(p.name for p in downloaded) == [
        "base_apogee.leb2",
        "base_asteroids.leb2",
        "base_core.leb2",
    ]


def test_download_leb2_for_tier_invalid_cache_redownload(tmp_path, monkeypatch):
    """Pre-existing INVALID file falls through to download (arc 1515->1521)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    # core exists but is invalid -> _is_valid_leb False -> download path.
    # The freshly downloaded file (b"d") must validate so the post-download
    # corruption guard accepts it; only the pre-existing b"corrupt" is invalid.
    (leb_dir / "base_core.leb2").write_bytes(b"corrupt")
    monkeypatch.setattr(
        dl, "_is_valid_leb", lambda p: Path(p).read_bytes() != b"corrupt"
    )
    monkeypatch.setattr(
        dl,
        "download_file",
        lambda url, dest_path, **kw: Path(dest_path).write_bytes(b"d"),
    )

    from libephemeris import state

    monkeypatch.setattr(state, "set_leb_file", lambda p: None)

    downloaded = dl.download_leb2_for_tier(
        "base", groups=["core"], force=False, quiet=True, activate=True
    )
    assert len(downloaded) == 1


def test_download_leb2_replaces_valid_but_unreviewed_bundled_core(
    tmp_path, monkeypatch, capsys
):
    """A structurally valid old release core cannot bypass the pinned hash."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    core = leb_dir / "base_core.leb2"
    core.write_bytes(b"old structurally valid core")
    monkeypatch.setattr(dl, "_is_valid_leb", lambda path: True)
    installed = []

    def _install(resource_path, dest_path, *, expected_sha256, validator=None):
        installed.append((resource_path, expected_sha256))
        Path(dest_path).write_bytes(b"reviewed core")

    monkeypatch.setattr(dl, "_install_bundled_file", _install)

    result = dl.download_leb2_for_tier(
        "base", groups=["core"], activate=False, quiet=False
    )

    assert result == [core]
    assert installed == [
        (
            "data/leb2/base_core.leb2",
            "5d708bdbe3e799e0802ba575984e57a3c5e44720dbfa1b4a01cf826640e0cb82",
        )
    ]
    assert "replacing unreviewed cached core" in capsys.readouterr().out


def test_download_leb2_for_tier_activate_core_missing(tmp_path, monkeypatch):
    """core appended but core file not on disk skips activation (arc 1543->1549)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    monkeypatch.setattr(dl, "_is_valid_leb", lambda p: True)
    # The bundled installer "succeeds" without writing the file, so
    # core_path.exists() is False at activation time.
    monkeypatch.setattr(dl, "_install_bundled_file", lambda *args, **kwargs: None)

    from libephemeris import state

    called = {"n": 0}
    monkeypatch.setattr(
        state, "set_leb_file", lambda p: called.update(n=called["n"] + 1)
    )

    downloaded = dl.download_leb2_for_tier(
        "base", groups=["core"], force=True, quiet=True, activate=True
    )
    assert len(downloaded) == 1
    # core_path.exists() is False -> set_leb_file not called.
    assert called["n"] == 0


# ---------------------------------------------------------------------------
# _download_planet_centers_for_tier
# ---------------------------------------------------------------------------


def test_download_planet_centers_for_tier_invalid():
    """Unknown tier raises ValueError (1579-1580)."""
    with pytest.raises(ValueError, match="No planet_centers file"):
        dl._download_planet_centers_for_tier("nonsense")


def test_download_planet_centers_for_tier_existing_valid(tmp_path, monkeypatch, capsys):
    """Existing valid file short-circuits (1584-1588)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    dest = tmp_path / "planet_centers_base.bsp"
    dest.write_bytes(b"pc")
    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: True)
    monkeypatch.setattr(
        dl,
        "_file_sha256",
        lambda path: dl.DATA_FILES["planet_centers_base.bsp"]["sha256"],
    )
    result = dl._download_planet_centers_for_tier("base", force=False, quiet=False)
    assert result == dest
    out = capsys.readouterr().out
    assert "already exists" in out


def test_download_planet_centers_for_tier_existing_valid_quiet(tmp_path, monkeypatch):
    """Existing valid file, quiet=True skips the print (arc 1586->1588)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    dest = tmp_path / "planet_centers_base.bsp"
    dest.write_bytes(b"pc")
    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: True)
    monkeypatch.setattr(
        dl,
        "_file_sha256",
        lambda path: dl.DATA_FILES["planet_centers_base.bsp"]["sha256"],
    )
    result = dl._download_planet_centers_for_tier("base", force=False, quiet=True)
    assert result == dest


def test_download_planet_centers_for_tier_corrupt_then_download(tmp_path, monkeypatch):
    """A corrupt cache is replaced by a verified download, never before it."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    dest = tmp_path / "planet_centers_base.bsp"
    dest.write_bytes(b"corrupt")

    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: Path(p).read_bytes() == b"fresh")
    monkeypatch.setattr(dl, "download_file", _fake_download_file(b"fresh"))

    result = dl._download_planet_centers_for_tier("base", force=False, quiet=True)
    assert result == dest
    assert dest.read_bytes() == b"fresh"


def test_download_planet_centers_for_tier_keeps_cache_when_offline(
    tmp_path, monkeypatch
):
    """Issue #55: a digest-mismatched tier file survives a failed download.

    Reproduces the reported sequence — an artifact from an earlier data
    release, a pinned digest that no longer matches, and no network — and
    pins the invariant that the 191 MB file is still there afterwards.
    """
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    dest = tmp_path / "planet_centers_medium.bsp"
    dest.write_bytes(b"previous data release")
    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: True)
    monkeypatch.setattr(dl, "_file_sha256", lambda path: "0" * 64)

    def _offline(url, dest_path, **kwargs):
        raise OSError("network is unreachable")

    monkeypatch.setattr(dl, "download_file", _offline)

    with pytest.raises(OSError, match="network is unreachable"):
        dl._download_planet_centers_for_tier("medium", force=False, quiet=True)

    assert dest.read_bytes() == b"previous data release"


def test_download_planet_centers_for_tier_download_invalid(tmp_path, monkeypatch):
    """A structurally invalid download raises and leaves nothing behind."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: False)
    monkeypatch.setattr(dl, "download_file", _fake_download_file(b"x"))

    with pytest.raises(ValueError, match="failed validation"):
        dl._download_planet_centers_for_tier("base", force=True, quiet=False)
    assert not (tmp_path / "planet_centers_base.bsp").exists()


def test_download_planet_centers_for_tier_invalid_download_keeps_previous(
    tmp_path, monkeypatch
):
    """A corrupt payload never displaces the tier file already installed."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    dest = tmp_path / "planet_centers_base.bsp"
    dest.write_bytes(b"known good")
    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: False)
    monkeypatch.setattr(dl, "download_file", _fake_download_file(b"truncated"))

    with pytest.raises(ValueError, match="failed validation"):
        dl._download_planet_centers_for_tier("base", force=True, quiet=True)

    assert dest.read_bytes() == b"known good"


def test_download_planet_centers_for_tier_success(tmp_path, monkeypatch, capsys):
    """Happy path: download + validate + final print (1605-1625)."""
    monkeypatch.setattr(dl, "get_data_dir", lambda: tmp_path)
    monkeypatch.setattr(dl, "_is_valid_bsp", lambda p: True)
    monkeypatch.setattr(dl, "download_file", _fake_download_file(b"x"))
    result = dl._download_planet_centers_for_tier("medium", force=True, quiet=False)
    assert result == tmp_path / "planet_centers_medium.bsp"
    out = capsys.readouterr().out
    assert "Downloaded to" in out
