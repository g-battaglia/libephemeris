# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for :mod:`libephemeris.state`.

These tests drive the global-state management module through every
remaining uncovered line / branch arc by monkeypatching helper seams
(TOML config, Skyfield loader, file existence, download, etc.) so that no
network access or real ephemeris files are required. The autouse fixtures
in ``tests/conftest.py`` reset global state between tests.
"""

from __future__ import annotations

import os
from unittest.mock import Mock

import pytest

from libephemeris import state


# ---------------------------------------------------------------------------
# Small fakes/helpers
# ---------------------------------------------------------------------------


class _RaisingCloser:
    """Object whose close()/shutdown() raise to drive except branches."""

    def close(self):
        raise OSError("boom-close")

    def shutdown(self):
        raise OSError("boom-shutdown")


class _SilentCloser:
    """Object with no-op close()/shutdown()."""

    def close(self):
        return None

    def shutdown(self):
        return None


class _FakeLoader:
    """Stand-in for skyfield Loader that never touches the network."""

    def __init__(self, kernel=None):
        self._kernel = kernel if kernel is not None else _SilentCloser()
        self.calls = []

    def __call__(self, name):
        self.calls.append(name)
        return self._kernel

    def timescale(self):
        return "fake-timescale"


# ===========================================================================
# _resolve_data_dir  (line 53)
# ===========================================================================


def test_resolve_data_dir_env_var(monkeypatch, tmp_path):
    """Env var set -> os.path.abspath(env_value) (line 53)."""
    monkeypatch.setenv(state._DATA_DIR_ENV_VAR, str(tmp_path))
    assert state._resolve_data_dir() == os.path.abspath(str(tmp_path))


# ===========================================================================
# get_calc_mode TOML fallback  (line 279)
# ===========================================================================


def test_get_calc_mode_toml_fallback(monkeypatch):
    """No override/env -> TOML returns a valid mode (line 279)."""
    monkeypatch.setattr(state, "_CALC_MODE", None)
    monkeypatch.delenv(state._CALC_MODE_ENV_VAR, raising=False)
    monkeypatch.setattr(
        "libephemeris._config_toml.get_str",
        lambda key: "leb" if key == "mode" else None,
    )
    assert state.get_calc_mode() == "leb"


# ===========================================================================
# set_leb_file  (lines 304-305)
# ===========================================================================


def test_set_leb_file_close_raises(monkeypatch):
    """Closing an existing reader that raises -> debug log (304-305)."""
    monkeypatch.setattr(state, "_LEB_READER", _RaisingCloser())
    monkeypatch.setattr(state, "_LEB_FILE", "/old.leb")
    # Should swallow the OSError from close() and log debug.
    state.set_leb_file("/new.leb")
    assert state._LEB_FILE == "/new.leb"
    assert state._LEB_READER is None


# ===========================================================================
# _discover_leb_file  (lines 336, 341, 349-388)
# ===========================================================================


def test_reviewed_base_core_presence_gate(tmp_path):
    """The packaged core passes; a path with nothing at it does not.

    Since 3.0.0rc15 the gate is presence, not content: an arbitrary file under
    the expected name passes too, deliberately. Only absence is rejected.
    """
    bundled = os.path.join(
        os.path.dirname(state.__file__), "data", "leb2", "base_core.leb2"
    )
    arbitrary = tmp_path / "base_core.leb2"
    arbitrary.write_bytes(b"not the reviewed artifact")

    assert state._is_reviewed_base_core(bundled)
    assert state._is_reviewed_base_core(str(arbitrary))
    assert not state._is_reviewed_base_core(str(tmp_path / "absent.leb2"))


def test_discover_leb_file_reviewed_cached_core(monkeypatch, tmp_path):
    """A cached core is auto-discovered only when its pinned hash matches."""
    monkeypatch.setattr(state, "get_precision_tier", lambda: "base")
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    target = leb_dir / "base_core.leb2"
    target.write_bytes(b"x")
    monkeypatch.setattr(
        state, "_is_reviewed_base_core", lambda path: path == str(target)
    )
    assert state._discover_leb_file() == str(target)


def test_discover_leb_file_reviewed_medium_core(monkeypatch, tmp_path):
    """A wider cached core is trusted only through its tier-specific pin."""
    monkeypatch.setattr(state, "get_precision_tier", lambda: "medium")
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    target = leb_dir / "medium_core.leb2"
    target.write_bytes(b"pinned")
    monkeypatch.setattr(
        state,
        "_is_reviewed_tier_core",
        lambda path, tier: path == str(target) and tier == "medium",
    )

    assert state._discover_leb_file() == str(target)


def test_discover_leb_file_medium_legacy_core_precedes_base_fallback(
    monkeypatch, tmp_path
):
    """The standard legacy core filename retains rc7 auto-discovery."""
    monkeypatch.setattr(state, "get_precision_tier", lambda: "medium")
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    target = leb_dir / "medium_core.leb"
    target.write_bytes(b"x")
    assert state._discover_leb_file() == str(target)


@pytest.mark.parametrize("tier", ["base", "medium"])
def test_discover_leb_file_legacy_ephemeris_name(monkeypatch, tmp_path, tier):
    """The historical ephemeris_{tier}.leb path remains discoverable."""
    monkeypatch.setattr(state, "get_precision_tier", lambda: tier)
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    target = leb_dir / f"ephemeris_{tier}.leb"
    target.write_bytes(b"legacy user file")

    assert state._discover_leb_file() == str(target)


def test_discover_leb_file_bundled_base(monkeypatch, tmp_path):
    """The hash-pinned bundled base core is the safe fallback."""
    monkeypatch.setattr(state, "get_precision_tier", lambda: "base")
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))

    pkg_dir = os.path.dirname(state.__file__)
    bundled = os.path.join(pkg_dir, "data", "leb2", "base_core.leb2")
    monkeypatch.setattr(state, "_is_reviewed_base_core", lambda path: path == bundled)
    assert state._discover_leb_file() == bundled


def test_discover_leb_file_unreviewed_cache_falls_back_to_bundle(monkeypatch, tmp_path):
    """A stale cached core cannot shadow the reviewed package artifact."""
    monkeypatch.setattr(state, "get_precision_tier", lambda: "base")
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    leb_dir = tmp_path / "leb"
    leb_dir.mkdir()
    cached = leb_dir / "base_core.leb2"
    cached.write_bytes(b"unreviewed")
    pkg_dir = os.path.dirname(state.__file__)
    bundled = os.path.join(pkg_dir, "data", "leb2", "base_core.leb2")
    monkeypatch.setattr(state, "_is_reviewed_base_core", lambda path: path == bundled)
    logger = Mock()
    monkeypatch.setattr(state, "get_logger", lambda: logger)

    assert state._discover_leb_file() == bundled
    assert "Ignoring unusable cached LEB core" in logger.warning.call_args.args[0]


def test_discover_leb_file_rejects_bad_bundle(monkeypatch, tmp_path):
    """A bundled file with the wrong hash is never opened."""
    monkeypatch.setattr(state, "get_precision_tier", lambda: "base")
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    monkeypatch.setattr(state, "_is_reviewed_base_core", lambda path: False)
    logger = Mock()
    monkeypatch.setattr(state, "get_logger", lambda: logger)

    assert state._discover_leb_file() is None
    assert (
        "Bundled LEB core failed its pinned SHA-256" in logger.error.call_args.args[0]
    )


def test_reviewed_core_does_not_attach_cached_companions(monkeypatch, tmp_path):
    """Opening the pinned core uses one reader, not prefix auto-discovery."""
    path = tmp_path / "base_core.leb2"
    path.write_bytes(b"reviewed")
    sentinel = object()
    opened = []

    import libephemeris.leb_composite as composite
    import libephemeris.leb_reader as leb_reader

    monkeypatch.setattr(state, "_LEB_FILE", str(path))
    monkeypatch.setattr(state, "_LEB_READER", None)
    monkeypatch.setattr(state, "_is_reviewed_base_core", lambda value: True)
    monkeypatch.setattr(
        leb_reader,
        "open_leb",
        lambda value: opened.append(value) or sentinel,
    )
    monkeypatch.setattr(
        composite.CompositeLEBReader,
        "from_file_with_companions",
        lambda value: pytest.fail("reviewed core attached unreviewed companions"),
    )
    monkeypatch.setattr(state, "_release_when_unused", lambda reader: None)
    monkeypatch.setattr(state, "_maybe_warm_reader", lambda reader: None)

    assert state._get_leb_reader_locked("auto") is sentinel
    assert opened == [str(path)]


def test_reviewed_medium_core_does_not_attach_cached_companions(monkeypatch, tmp_path):
    path = tmp_path / "medium_core.leb2"
    path.write_bytes(b"reviewed")
    sentinel = object()
    opened = []

    import libephemeris.leb_composite as composite
    import libephemeris.leb_reader as leb_reader

    monkeypatch.setattr(state, "_LEB_FILE", str(path))
    monkeypatch.setattr(state, "_LEB_READER", None)
    monkeypatch.setattr(state, "_is_reviewed_core", lambda value: True)
    monkeypatch.setattr(
        leb_reader,
        "open_leb",
        lambda value: opened.append(value) or sentinel,
    )
    monkeypatch.setattr(
        composite.CompositeLEBReader,
        "from_file_with_companions",
        lambda value: pytest.fail("reviewed core attached unreviewed companions"),
    )
    monkeypatch.setattr(state, "_release_when_unused", lambda reader: None)
    monkeypatch.setattr(state, "_maybe_warm_reader", lambda reader: None)

    assert state._get_leb_reader_locked("auto") is sentinel
    assert opened == [str(path)]


# ===========================================================================
# _maybe_warm_reader  (lines 424-440)
# ===========================================================================


class _WarmReader:
    def __init__(self, fail=False):
        self.fail = fail
        self.warmed = None

    def warm(self, jd_start, jd_end):
        if self.fail:
            raise OSError("warm failed")
        self.warmed = (jd_start, jd_end)


def test_maybe_warm_reader_disabled(monkeypatch):
    """mmap_preload false -> early return (line 422 branch)."""
    monkeypatch.setattr("libephemeris._config_toml.get_bool", lambda key: False)
    reader = _WarmReader()
    state._maybe_warm_reader(reader)
    assert reader.warmed is None


def test_maybe_warm_reader_success(monkeypatch):
    """mmap_preload true -> warm() called with JD range (lines 424-437)."""
    monkeypatch.setattr("libephemeris._config_toml.get_bool", lambda key: True)
    monkeypatch.setattr(
        "libephemeris._config_toml.get_int",
        lambda key: {"mmap_preload_start": 1900, "mmap_preload_end": 2000}.get(key),
    )
    reader = _WarmReader()
    state._maybe_warm_reader(reader)
    assert reader.warmed is not None
    assert reader.warmed[0] < reader.warmed[1]


def test_maybe_warm_reader_defaults_and_failure(monkeypatch):
    """mmap_preload true, defaults used, warm() raises (lines 424-425, 438-440)."""
    monkeypatch.setattr("libephemeris._config_toml.get_bool", lambda key: True)
    # get_int returns None -> the `or 1800` / `or 2200` defaults apply.
    monkeypatch.setattr("libephemeris._config_toml.get_int", lambda key: None)
    reader = _WarmReader(fail=True)
    # Failure is swallowed (logged as warning).
    state._maybe_warm_reader(reader)
    assert reader.warmed is None


# ===========================================================================
# release_data_cache  (lines 460-485)
# ===========================================================================


def test_release_data_cache_unavailable(monkeypatch):
    """posix_fadvise unavailable -> early return (line 458-459 branch)."""
    monkeypatch.setattr(state.os, "posix_fadvise", None, raising=False)
    # Should not raise.
    state.release_data_cache()


def test_release_data_cache_no_dir(monkeypatch, tmp_path):
    """fadvise present but data dir missing -> early return (461-462)."""
    monkeypatch.setattr(state.os, "posix_fadvise", lambda *a: None, raising=False)
    monkeypatch.setattr(state.os, "POSIX_FADV_DONTNEED", 4, raising=False)
    monkeypatch.setattr(state, "_resolve_data_dir", lambda: str(tmp_path / "nope"))
    state.release_data_cache()


def test_release_data_cache_walks_files(monkeypatch, tmp_path):
    """fadvise present + real files -> walks and advises (lines 460-485)."""
    fadvise_calls = []

    def fake_fadvise(fd, offset, length, advice):
        fadvise_calls.append((offset, length, advice))

    monkeypatch.setattr(state.os, "posix_fadvise", fake_fadvise, raising=False)
    monkeypatch.setattr(state.os, "POSIX_FADV_DONTNEED", 4, raising=False)

    data_dir = tmp_path / "data"
    sub = data_dir / "sub"
    sub.mkdir(parents=True)
    good = data_dir / "good.bin"
    good.write_bytes(b"abc")
    nested = sub / "nested.bin"
    nested.write_bytes(b"defgh")
    # A FIFO is a non-regular file: os.stat succeeds, S_ISREG is False ->
    # the `continue` at line 471 executes.
    fifo = data_dir / "pipe.fifo"
    try:
        os.mkfifo(str(fifo))
    except (OSError, AttributeError):
        pass

    monkeypatch.setattr(state, "_resolve_data_dir", lambda: str(data_dir))
    state.release_data_cache()
    # At least the two regular files were advised.
    assert len(fadvise_calls) >= 2


def test_release_data_cache_close_oserror(monkeypatch, tmp_path):
    """os.close raising in the finally block -> swallowed (lines 484-485)."""
    monkeypatch.setattr(state.os, "posix_fadvise", lambda *a: None, raising=False)
    monkeypatch.setattr(state.os, "POSIX_FADV_DONTNEED", 4, raising=False)

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "f.bin").write_bytes(b"x")
    monkeypatch.setattr(state, "_resolve_data_dir", lambda: str(data_dir))

    real_close = state.os.close

    def boom_close(fd):
        # Close the fd for real (avoid a descriptor leak) then raise so
        # the except OSError branch in the finally clause is exercised.
        real_close(fd)
        raise OSError("close failed")

    monkeypatch.setattr(state.os, "close", boom_close)
    state.release_data_cache()


def test_release_data_cache_stat_oserror(monkeypatch, tmp_path):
    """os.stat raising OSError inside loop -> continue (lines 472-473)."""
    monkeypatch.setattr(state.os, "posix_fadvise", lambda *a: None, raising=False)
    monkeypatch.setattr(state.os, "POSIX_FADV_DONTNEED", 4, raising=False)

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    f = data_dir / "f.bin"
    f.write_bytes(b"x")

    monkeypatch.setattr(state, "_resolve_data_dir", lambda: str(data_dir))

    real_stat = state.os.stat

    def boom_stat(path, *a, **k):
        # Only the regular file inside the loop raises; directory traversal
        # via os.walk keeps using the real os.stat.
        if str(path).endswith("f.bin"):
            raise OSError("stat failed")
        return real_stat(path, *a, **k)

    monkeypatch.setattr(state.os, "stat", boom_stat)
    state.release_data_cache()


def test_release_data_cache_open_oserror(monkeypatch, tmp_path):
    """os.open raising OSError inside loop -> pass (lines 478-479)."""
    monkeypatch.setattr(state.os, "posix_fadvise", lambda *a: None, raising=False)
    monkeypatch.setattr(state.os, "POSIX_FADV_DONTNEED", 4, raising=False)

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    f = data_dir / "f.bin"
    f.write_bytes(b"x")

    monkeypatch.setattr(state, "_resolve_data_dir", lambda: str(data_dir))

    def boom_open(*a, **k):
        raise OSError("open failed")

    monkeypatch.setattr(state.os, "open", boom_open)
    state.release_data_cache()


# ===========================================================================
# get_horizons_client  (lines 601, 612, 615->614, 618-620, 626-640)
# ===========================================================================


def test_get_horizons_client_mode_horizons(monkeypatch):
    """mode horizons -> _get_or_create_horizons_client (line 601)."""
    monkeypatch.setattr(state, "get_calc_mode", lambda: "horizons")
    sentinel = object()
    monkeypatch.setattr(state, "_get_or_create_horizons_client", lambda: sentinel)
    assert state.get_horizons_client() is sentinel


def test_get_horizons_client_auto_with_local(monkeypatch, tmp_path):
    """mode auto + local ephemeris present -> None (lines 612, 615->616)."""
    monkeypatch.setattr(state, "get_calc_mode", lambda: "auto")
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    monkeypatch.setattr(state, "_get_effective_ephemeris_file", lambda: "de440.bsp")
    monkeypatch.setattr(state, "_EPHEMERIS_PATH", str(tmp_path / "ep"))
    # Make the second candidate (the _EPHEMERIS_PATH one) the existing file.
    ep_dir = tmp_path / "ep"
    ep_dir.mkdir()
    (ep_dir / "de440.bsp").write_bytes(b"x")
    assert state.get_horizons_client() is None


def test_get_horizons_client_auto_no_local(monkeypatch, tmp_path):
    """mode auto + no local ephemeris -> create client (615->614, 618)."""
    monkeypatch.setattr(state, "get_calc_mode", lambda: "auto")
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    monkeypatch.setattr(state, "_get_effective_ephemeris_file", lambda: "de440.bsp")
    monkeypatch.setattr(state, "_EPHEMERIS_PATH", str(tmp_path / "ep"))
    (tmp_path / "ep").mkdir()
    # Neither candidate file exists -> loop iterates and falls through.
    sentinel = object()
    monkeypatch.setattr(state, "_get_or_create_horizons_client", lambda: sentinel)
    assert state.get_horizons_client() is sentinel


def test_get_horizons_client_auto_no_ephemeris_path(monkeypatch, tmp_path):
    """mode auto + _EPHEMERIS_PATH falsy -> skip append (arc 611->614)."""
    monkeypatch.setattr(state, "get_calc_mode", lambda: "auto")
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    monkeypatch.setattr(state, "_get_effective_ephemeris_file", lambda: "de440.bsp")
    monkeypatch.setattr(state, "_EPHEMERIS_PATH", None)
    sentinel = object()
    monkeypatch.setattr(state, "_get_or_create_horizons_client", lambda: sentinel)
    # No file exists in tmp_path -> falls through to client creation.
    assert state.get_horizons_client() is sentinel


def test_get_horizons_client_unknown_mode(monkeypatch):
    """A non-standard mode value -> final return None (line 620)."""
    # Force a mode that is none of skyfield/leb/horizons/auto so we fall
    # through every branch to the trailing `return None`.
    monkeypatch.setattr(state, "get_calc_mode", lambda: "weird")
    assert state.get_horizons_client() is None


def test_get_or_create_horizons_client(monkeypatch):
    """Create singleton client + emit warning once (lines 626-640)."""
    monkeypatch.setattr(state, "_HORIZONS_CLIENT", None)
    monkeypatch.setattr(state, "_HORIZONS_WARNED", False)

    created = []

    class _FakeClient:
        def __init__(self):
            created.append(self)

    import libephemeris.horizons_backend as hb

    monkeypatch.setattr(hb, "HorizonsClient", _FakeClient)

    client = state._get_or_create_horizons_client()
    assert isinstance(client, _FakeClient)
    assert state._HORIZONS_WARNED is True
    # Second call returns cached instance, no new creation.
    again = state._get_or_create_horizons_client()
    assert again is client
    assert len(created) == 1


def test_get_or_create_horizons_client_already_warned(monkeypatch):
    """Client None but already warned -> skip warning (arc 632->640)."""
    monkeypatch.setattr(state, "_HORIZONS_CLIENT", None)
    monkeypatch.setattr(state, "_HORIZONS_WARNED", True)

    class _FakeClient:
        pass

    import libephemeris.horizons_backend as hb

    monkeypatch.setattr(hb, "HorizonsClient", _FakeClient)
    client = state._get_or_create_horizons_client()
    assert isinstance(client, _FakeClient)


def test_get_or_create_horizons_client_lock_race(monkeypatch):
    """Inner ``if _HORIZONS_CLIENT is None`` False -> skip create (628->632)."""
    monkeypatch.setattr(state, "_HORIZONS_CLIENT", None)
    monkeypatch.setattr(state, "_HORIZONS_WARNED", True)
    sentinel = object()
    monkeypatch.setattr(state, "_INIT_LOCK", _SettingLock("_HORIZONS_CLIENT", sentinel))
    assert state._get_or_create_horizons_client() is sentinel


# ===========================================================================
# get_loader  (line 657->659)
# ===========================================================================


def test_get_loader_creates(monkeypatch, tmp_path):
    """_LOADER None -> build inside lock (lines 657-658)."""
    monkeypatch.setattr(state, "_LOADER", None)
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    captured = {}

    class _FakeLoaderCls:
        def __init__(self, d):
            captured["dir"] = d

    monkeypatch.setattr("libephemeris.state.Loader", _FakeLoaderCls)
    loader = state.get_loader()
    assert isinstance(loader, _FakeLoaderCls)
    assert captured["dir"] == str(tmp_path)
    # Second call returns cached (covers 656/659 cached path).
    assert state.get_loader() is loader


class _SettingLock:
    """RLock-like context manager that sets a global on __enter__.

    Simulates the double-checked-locking race: another thread populates the
    singleton between the outer ``is None`` check and acquiring the lock, so
    the inner ``is None`` check is False (the ``->return`` arc).
    """

    def __init__(self, attr, value):
        self._attr = attr
        self._value = value

    def __enter__(self):
        setattr(state, self._attr, self._value)
        return self

    def __exit__(self, *exc):
        return False


def test_get_loader_double_checked_race(monkeypatch):
    """Inner ``if _LOADER is None`` False -> skip build (arc 657->659)."""
    monkeypatch.setattr(state, "_LOADER", None)
    sentinel = object()
    monkeypatch.setattr(state, "_INIT_LOCK", _SettingLock("_LOADER", sentinel))
    assert state.get_loader() is sentinel


# ===========================================================================
# _build_enhanced_timescale fallback (lines 708-709) via get_timescale
# ===========================================================================


def test_build_enhanced_timescale_missing_hist(monkeypatch):
    """historic_deltat.npy missing -> standard timescale fallback (708-709)."""
    np = pytest.importorskip("numpy")

    import skyfield.data
    from pathlib import Path

    data_dir = Path(skyfield.data.__file__).parent

    real_exists = Path.exists

    def fake_exists(self):
        if self.name == "historic_deltat.npy":
            return False
        return real_exists(self)

    monkeypatch.setattr(Path, "exists", fake_exists)

    # Provide a fake loader so we never download.
    fake_loader = _FakeLoader()
    monkeypatch.setattr(state, "get_loader", lambda: fake_loader)

    # Ensure the iers.npz load path works (it exists in skyfield).
    assert (data_dir / "iers.npz").exists()
    assert np is not None
    result = state._build_enhanced_timescale()
    assert result == "fake-timescale"


# ===========================================================================
# get_timescale  (lines 743->756, 746-755)
# ===========================================================================


def test_get_timescale_cached(monkeypatch):
    """_TS already set on outer check -> skip build entirely."""
    sentinel = object()
    monkeypatch.setattr(state, "_TS", sentinel)
    assert state.get_timescale() is sentinel


def test_get_timescale_double_checked_race(monkeypatch):
    """Inner ``if _TS is None`` False -> skip build (arc 743->756)."""
    monkeypatch.setattr(state, "_TS", None)
    sentinel = object()
    monkeypatch.setattr(state, "_INIT_LOCK", _SettingLock("_TS", sentinel))
    assert state.get_timescale() is sentinel


def test_get_timescale_build_fails_fallback(monkeypatch):
    """Enhanced build raises -> warning + standard timescale (746-755)."""
    monkeypatch.setattr(state, "_TS", None)

    def boom():
        raise ValueError("merge failed")

    monkeypatch.setattr(state, "_build_enhanced_timescale", boom)
    fake_loader = _FakeLoader()
    monkeypatch.setattr(state, "get_loader", lambda: fake_loader)
    ts = state.get_timescale()
    assert ts == "fake-timescale"


# ===========================================================================
# _get_current_tier TOML  (line 803)
# ===========================================================================


def test_get_current_tier_toml(monkeypatch):
    """TOML returns a valid tier (line 803)."""
    monkeypatch.setattr(state, "_PRECISION_TIER", None)
    monkeypatch.delenv(state._PRECISION_TIER_ENV_VAR, raising=False)
    monkeypatch.setattr(
        "libephemeris._config_toml.get_str",
        lambda key: "extended" if key == "precision" else None,
    )
    assert state._get_current_tier().name == "extended"


# ===========================================================================
# set_precision_tier close raises  (lines 858-859)
# ===========================================================================


def test_set_precision_tier_close_raises(monkeypatch):
    """_PLANETS.close() raises -> swallowed (lines 858-859)."""
    monkeypatch.setattr(state, "_PLANETS", _RaisingCloser())
    state.set_precision_tier("base")
    assert state._PLANETS is None
    assert state._PRECISION_TIER == "base"


# ===========================================================================
# _get_effective_ephemeris_file TOML  (line 932)
# ===========================================================================


def test_get_effective_ephemeris_file_toml(monkeypatch):
    """No env / no explicit -> TOML ephemeris returned (line 932)."""
    monkeypatch.delenv(state._EPHEMERIS_ENV_VAR, raising=False)
    monkeypatch.setattr(state, "_EPHEMERIS_FILE_EXPLICIT", False)
    monkeypatch.setattr(
        "libephemeris._config_toml.get_str",
        lambda key: "de441.bsp" if key == "ephemeris" else None,
    )
    assert state._get_effective_ephemeris_file() == "de441.bsp"


# ===========================================================================
# get_planets  (lines 961->1011, 969->976, 987-1006)
# ===========================================================================


def test_get_planets_cached(monkeypatch):
    """_PLANETS already set on outer check -> return immediately."""
    sentinel = object()
    monkeypatch.setattr(state, "_PLANETS", sentinel)
    monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
    assert state.get_planets() is sentinel


def test_get_planets_double_checked_race(monkeypatch):
    """Inner ``if _PLANETS is None`` False -> skip load (arc 961->1011)."""
    monkeypatch.setattr(state, "_PLANETS", None)
    sentinel = object()
    monkeypatch.setattr(state, "_INIT_LOCK", _SettingLock("_PLANETS", sentinel))
    monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
    assert state.get_planets() is sentinel


def test_get_planets_custom_path_missing_falls_to_data_dir(monkeypatch, tmp_path):
    """_EPHEMERIS_PATH set but file missing -> data dir path (969->976)."""
    monkeypatch.setattr(state, "_PLANETS", None)
    monkeypatch.setattr(state, "_EPHEMERIS_PATH", str(tmp_path / "custom"))
    monkeypatch.setattr(state, "_get_effective_ephemeris_file", lambda: "de440.bsp")

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "de440.bsp").write_bytes(b"x")
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(data_dir))

    fake_loader = _FakeLoader(kernel="loaded-kernel")
    monkeypatch.setattr(state, "get_loader", lambda: fake_loader)
    monkeypatch.setattr(state, "get_calc_mode", lambda: "auto")

    result = state.get_planets()
    assert result == "loaded-kernel"
    assert fake_loader.calls == [str(data_dir / "de440.bsp")]


def test_get_planets_leb_mode_rejects_existing_kernel(monkeypatch, tmp_path):
    """LEB runtime never opens a local DE kernel, even when it exists."""
    monkeypatch.setattr(state, "_PLANETS", None)
    monkeypatch.setattr(state, "_EPHEMERIS_PATH", None)
    monkeypatch.setattr(state, "_get_effective_ephemeris_file", lambda: "de441.bsp")

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    # de441 missing, but a fallback (de440.bsp) is present.
    (data_dir / "de440.bsp").write_bytes(b"x")
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(data_dir))
    monkeypatch.setattr(state, "get_calc_mode", lambda: "leb")

    fake_loader = _FakeLoader(kernel="fb-kernel")
    monkeypatch.setattr(state, "get_loader", lambda: fake_loader)

    with pytest.raises(RuntimeError, match="disabled"):
        state.get_planets()
    assert fake_loader.calls == []


def test_get_planets_leb_mode_never_downloads_kernel(monkeypatch, tmp_path):
    """LEB runtime fails before a loader can download a DE kernel."""
    monkeypatch.setattr(state, "_PLANETS", None)
    monkeypatch.setattr(state, "_EPHEMERIS_PATH", None)
    monkeypatch.setattr(state, "_get_effective_ephemeris_file", lambda: "de441.bsp")

    data_dir = tmp_path / "data"
    data_dir.mkdir()  # empty: no DE files
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(data_dir))
    monkeypatch.setattr(state, "get_calc_mode", lambda: "leb")

    fake_loader = _FakeLoader(kernel="dl-kernel")
    monkeypatch.setattr(state, "get_loader", lambda: fake_loader)

    with pytest.raises(RuntimeError, match="disabled"):
        state.get_planets()
    assert fake_loader.calls == []


# ===========================================================================
# _get_planet_centers_locked  (1051->1077, 1066->1077, 1067->1066, 1072-1073)
# ===========================================================================


def test_get_planet_centers_locked_already_loaded(monkeypatch):
    """Condition false (loaded + same tier) -> return early (1051->1077)."""
    sentinel = object()
    monkeypatch.setattr(state, "_PLANET_CENTERS", sentinel)
    monkeypatch.setattr(state, "_PLANET_CENTERS_TIER", "base")
    assert state._get_planet_centers_locked("base") is sentinel


def test_get_planet_centers_locked_no_candidates(monkeypatch, tmp_path):
    """No candidate files exist -> loop exhausts, returns None (1066->1077)."""
    monkeypatch.setattr(state, "_PLANET_CENTERS", None)
    monkeypatch.setattr(state, "_PLANET_CENTERS_TIER", None)
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    monkeypatch.setattr(state, "get_loader", lambda: _FakeLoader())
    assert state._get_planet_centers_locked("medium") is None


def test_get_planet_centers_locked_base_legacy_loads_when_pinned(monkeypatch, tmp_path):
    """The legacy destination is accepted only for the exact base artifact."""
    monkeypatch.setattr(state, "_PLANET_CENTERS", None)
    monkeypatch.setattr(state, "_PLANET_CENTERS_TIER", None)
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    # Only the legacy candidate "planet_centers.bsp" exists.
    legacy = tmp_path / "planet_centers.bsp"
    legacy.write_bytes(b"x")
    fake_loader = _FakeLoader(kernel="pc-kernel")
    monkeypatch.setattr(state, "get_loader", lambda: fake_loader)
    monkeypatch.setattr(
        state,
        "_matches_pinned_data_file",
        lambda path, name: name == "planet_centers.bsp",
    )
    result = state._get_planet_centers_locked("base")
    assert result == "pc-kernel"
    assert state._PLANET_CENTERS_TIER == "base"


def test_get_planet_centers_locked_loads_the_named_cache(monkeypatch, tmp_path):
    """A cache under the expected filename is loaded without a content check.

    Content verification moved to install time in 3.0.0rc15, so presence under
    a manifest name is the whole requirement here. A file that is not a usable
    SPK still fails — just later and more precisely, when the loader raises and
    the caller falls through to the next candidate.
    """
    monkeypatch.setattr(state, "_PLANET_CENTERS", None)
    monkeypatch.setattr(state, "_PLANET_CENTERS_TIER", None)
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    candidate = tmp_path / "planet_centers_medium.bsp"
    candidate.write_bytes(b"arbitrary bytes")
    fake_loader = _FakeLoader(kernel="pc-kernel")
    monkeypatch.setattr(state, "get_loader", lambda: fake_loader)

    assert state._get_planet_centers_locked("medium") == "pc-kernel"
    assert fake_loader.calls == [str(candidate)]


def test_rejected_planet_center_check_and_warning_are_cached(monkeypatch, tmp_path):
    """An unchanged rejected SPK is examined and warned about only once.

    The rejection path is now reachable only for a candidate that exists but
    is not a usable regular file, so the predicate is forced here; what is
    under test is the negative cache around it, which still guards the retry.
    """
    monkeypatch.setattr(state, "_PLANET_CENTERS", None)
    monkeypatch.setattr(state, "_PLANET_CENTERS_TIER", None)
    monkeypatch.setattr(state, "_PLANET_CENTER_CACHE_TIER", None)
    monkeypatch.setattr(state, "_PLANET_CENTER_REJECTED_FILES", set())
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    monkeypatch.setattr(state, "get_precision_tier", lambda: "medium")
    monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
    monkeypatch.setattr(state, "get_loader", lambda: _FakeLoader())
    logger = Mock()
    monkeypatch.setattr(state, "get_logger", lambda: logger)
    candidate = tmp_path / "planet_centers_medium.bsp"
    candidate.write_bytes(b"stale planet-center bytes")

    matcher = Mock(return_value=False)
    monkeypatch.setattr(state, "_matches_pinned_data_file", matcher)

    assert state.get_planet_centers() is None
    assert state.get_planet_centers() is None

    assert matcher.call_count == 1
    assert logger.warning.call_count == 1
    assert len(state._PLANET_CENTER_REJECTED_FILES) == 1


def test_replaced_planet_center_file_retries_hash_and_load(monkeypatch, tmp_path):
    """A stat change invalidates the rejection and checks replacement bytes."""
    monkeypatch.setattr(state, "_PLANET_CENTERS", None)
    monkeypatch.setattr(state, "_PLANET_CENTERS_TIER", None)
    monkeypatch.setattr(state, "_PLANET_CENTER_CACHE_TIER", None)
    monkeypatch.setattr(state, "_PLANET_CENTER_REJECTED_FILES", set())
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    monkeypatch.setattr(state, "get_precision_tier", lambda: "medium")
    monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
    fake_loader = _FakeLoader(kernel="replacement-kernel")
    monkeypatch.setattr(state, "get_loader", lambda: fake_loader)
    matcher = Mock(side_effect=[False, True])
    monkeypatch.setattr(state, "_matches_pinned_data_file", matcher)
    candidate = tmp_path / "planet_centers_medium.bsp"
    candidate.write_bytes(b"stale")

    assert state.get_planet_centers() is None
    candidate.write_bytes(b"reviewed replacement with a new stat fingerprint")
    assert state.get_planet_centers() == "replacement-kernel"

    assert matcher.call_count == 2
    assert fake_loader.calls == [str(candidate)]
    assert state._PLANET_CENTER_REJECTED_FILES == set()


def test_planet_center_rejection_cache_resets_on_tier_change_and_close():
    """Tier changes and full teardown clear every planet-center miss cache."""
    fingerprint = ("/tmp/stale.bsp", 123, 456, "planet_centers_medium.bsp")
    state._PLANET_CENTER_REJECTED_FILES.add(fingerprint)
    state._PLANET_CENTER_CACHE_TIER = "medium"
    state._PLANET_CENTER_MISSING.add(599)

    state.set_precision_tier("base")

    assert state._PLANET_CENTER_REJECTED_FILES == set()
    assert state._PLANET_CENTER_CACHE_TIER is None
    assert state._PLANET_CENTER_MISSING == set()

    state._PLANET_CENTER_REJECTED_FILES.add(fingerprint)
    state._PLANET_CENTER_CACHE_TIER = "medium"
    state._PLANET_CENTER_MISSING.add(599)
    state.close()

    assert state._PLANET_CENTER_REJECTED_FILES == set()
    assert state._PLANET_CENTER_CACHE_TIER is None
    assert state._PLANET_CENTER_MISSING == set()


def test_get_planet_centers_locked_load_raises(monkeypatch, tmp_path):
    """Candidate exists but load() raises -> warning, returns None (1072-1073)."""
    monkeypatch.setattr(state, "_PLANET_CENTERS", None)
    monkeypatch.setattr(state, "_PLANET_CENTERS_TIER", None)
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(tmp_path))
    candidate = tmp_path / "planet_centers_medium.bsp"
    candidate.write_bytes(b"x")

    class _BoomLoader:
        def __call__(self, path):
            raise ValueError("bad spk")

    monkeypatch.setattr(state, "get_loader", lambda: _BoomLoader())
    monkeypatch.setattr(state, "_matches_pinned_data_file", lambda path, name: True)
    result = state._get_planet_centers_locked("medium")
    assert result is None


# ===========================================================================
# get_planet_center_segment  (1109, 1113, 1119-1127, 1133-1134)
# ===========================================================================


def test_get_planet_center_segment_negative_cache(monkeypatch):
    """naif_id in missing-cache -> None (line 1109)."""
    monkeypatch.setattr(state, "_PLANET_CENTER_MISSING", {599})
    assert state.get_planet_center_segment(599) is None


def test_get_planet_center_segment_no_centers(monkeypatch):
    """get_planet_centers() returns None -> None (line 1113)."""
    monkeypatch.setattr(state, "_PLANET_CENTER_MISSING", set())
    monkeypatch.setattr(state, "get_planet_centers", lambda: None)
    assert state.get_planet_center_segment(599) is None


class _Seg:
    def __init__(self, target, start_jd=None, end_jd=None, set_attrs=True):
        self.target = target
        if set_attrs:
            self.start_jd = start_jd
            self.end_jd = end_jd


class _Centers:
    def __init__(self, segments):
        self.segments = segments


def test_get_planet_center_segment_in_coverage(monkeypatch):
    """jd within coverage -> returns segment (1119-1123)."""
    monkeypatch.setattr(state, "_PLANET_CENTER_MISSING", set())
    seg = _Seg(599, start_jd=2400000.0, end_jd=2500000.0)
    monkeypatch.setattr(state, "get_planet_centers", lambda: _Centers([seg]))
    result = state.get_planet_center_segment(599, jd=2451545.0)
    assert result is seg


def test_get_planet_center_segment_dispatches_across_all_descriptors(monkeypatch):
    """A split target remains usable after the first descriptor ends."""
    monkeypatch.setattr(state, "_PLANET_CENTER_MISSING", set())
    early = _Seg(599, start_jd=2400000.0, end_jd=2450000.0)
    late = _Seg(599, start_jd=2450000.0, end_jd=2500000.0)
    for seg in (early, late):
        seg.center = 5
        seg.ephemeris = None
    monkeypatch.setattr(state, "get_planet_centers", lambda: _Centers([early, late]))

    stack = state.get_planet_center_segment(599)

    assert stack.segments == [early, late]
    assert state.get_planet_center_segment(599, jd=2490000.0) is late


def test_get_planet_center_segment_outside_coverage(monkeypatch):
    """jd outside coverage -> None (1123)."""
    monkeypatch.setattr(state, "_PLANET_CENTER_MISSING", set())
    seg = _Seg(599, start_jd=2400000.0, end_jd=2450000.0)
    monkeypatch.setattr(state, "get_planet_centers", lambda: _Centers([seg]))
    assert state.get_planet_center_segment(599, jd=2500000.0) is None


def test_get_planet_center_segment_coverage_check_error(monkeypatch):
    """Coverage compare raises TypeError -> debug + return segment (1125-1127)."""
    monkeypatch.setattr(state, "_PLANET_CENTER_MISSING", set())

    # start_jd/end_jd are non-numeric so ``start_jd <= jd <= end_jd``
    # raises TypeError, exercising the except branch.
    seg = _Seg(599, start_jd="not-a-number", end_jd="also-not")
    monkeypatch.setattr(state, "get_planet_centers", lambda: _Centers([seg]))
    result = state.get_planet_center_segment(599, jd=2451545.0)
    assert result is seg


def test_get_planet_center_segment_not_found(monkeypatch):
    """naif_id absent -> added to missing-cache, None (1133-1134)."""
    missing = set()
    monkeypatch.setattr(state, "_PLANET_CENTER_MISSING", missing)
    seg = _Seg(699, start_jd=2400000.0, end_jd=2500000.0)
    monkeypatch.setattr(state, "get_planet_centers", lambda: _Centers([seg]))
    assert state.get_planet_center_segment(599) is None
    assert 599 in missing


# ===========================================================================
# get_sid_mode overloads (1209->exit, 1213->exit) and full path
# ===========================================================================


def test_get_sid_mode_full(monkeypatch):
    """full=True -> tuple form."""
    state.set_sid_mode(5, t0=2451545.0, ayan_t0=23.5)
    mode, t0, ayan = state.get_sid_mode(full=True)
    assert mode == 5
    assert t0 == 2451545.0
    assert ayan == 23.5


def test_get_sid_mode_default():
    """No mode set -> default 0."""
    state.reset_session()
    assert state.get_sid_mode() == 0


# ===========================================================================
# set_ephe_path  (1327-1328, 1340-1342)
# ===========================================================================


def test_set_ephe_path_close_raises_and_discover_raises(monkeypatch, tmp_path):
    """_PLANETS.close() raises (1327-1328) + discover raises (1340-1342)."""
    monkeypatch.setattr(state, "_EPHEMERIS_PATH", None)
    monkeypatch.setattr(state, "_PLANETS", _RaisingCloser())

    real_dir = tmp_path / "ephe"
    real_dir.mkdir()

    import libephemeris.spk_auto as spk_auto

    def boom(path):
        raise OSError("discover failed")

    monkeypatch.setattr(spk_auto, "discover_local_spks", boom)
    state.set_ephe_path(str(real_dir))
    assert state._PLANETS is None
    assert state._EPHEMERIS_PATH == str(real_dir)


# ===========================================================================
# set_ephemeris_file  (1375-1376)
# ===========================================================================


def test_set_ephemeris_file_close_raises(monkeypatch):
    """_PLANETS.close() raises -> swallowed (1375-1376)."""
    monkeypatch.setattr(state, "_PLANETS", _RaisingCloser())
    state.set_ephemeris_file("de441.bsp")
    assert state._PLANETS is None
    assert state._EPHEMERIS_FILE == "de441.bsp"
    assert state._EPHEMERIS_FILE_EXPLICIT is True


# ===========================================================================
# _close_inner  (1706-1709, 1717-1718, 1739-1742, 1748-1749,
#                1784-1785, 1792-1793, 1797-1800, 1808-1809, 1822-1823)
# ===========================================================================


def test_close_inner_all_error_paths(monkeypatch):
    """All close()/shutdown() and ImportError fallbacks (many lines)."""
    monkeypatch.setattr(state, "_HORIZONS_CLIENT", _RaisingCloser())
    monkeypatch.setattr(state, "_LEB_READER", _RaisingCloser())
    monkeypatch.setattr(state, "_PLANETS", _RaisingCloser())
    monkeypatch.setattr(state, "_PLANET_CENTERS", _RaisingCloser())
    monkeypatch.setattr(state, "_SPK_KERNELS", {"a": _RaisingCloser()})
    monkeypatch.setattr(state, "_SPK_TYPE21_KERNELS", {"b": _RaisingCloser()})

    # Force the ImportError fallbacks (1784-1785, 1792-1793, 1822-1823).
    # The relevant statements use ``from . import iers_data`` /
    # ``from ._config_toml import reset`` / ``from .rebound_integration
    # import ...``; intercept __import__ on the fromlist names.
    import builtins

    real_import = builtins.__import__
    blocked = {"iers_data", "reset", "reset_assist_data_cache"}

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if fromlist and any(item in blocked for item in fromlist):
            raise ImportError(f"blocked fromlist {fromlist}")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    # Should swallow everything and reset state.
    state.close()
    assert state._PLANETS is None
    assert state._PLANET_CENTERS is None
    assert state._SPK_KERNELS == {}
    assert state._SPK_TYPE21_KERNELS == {}
    assert state._HORIZONS_CLIENT is None
    assert state._LEB_READER is None


# ===========================================================================
# get_current_file_data  (1878, 1885->1908, 1887->1908, 1893->1897,
#                          1897->1892, 1903, 1905, 1910->1918,
#                          1915-1916, 1920-1922)
# ===========================================================================


def test_get_current_file_data_not_applicable(monkeypatch):
    """ifno not in (0, 1) -> empty tuple."""
    assert state.get_current_file_data(2) == ("", 0.0, 0.0, 0)


def test_get_current_file_data_no_planets(monkeypatch):
    """_PLANETS None: report the active LEB artifact, or empty without one."""
    monkeypatch.setattr(state, "_PLANETS", None)
    result = state.get_current_file_data(0)
    try:
        reader = state.get_leb_reader()
    except RuntimeError:
        reader = None
    if reader is not None:
        path, start, end, denum = result
        assert path and start < end and denum in (440, 441)
    else:
        assert result == ("", 0.0, 0.0, 0)


def test_get_current_file_data_filename_attr_no_spk(monkeypatch):
    """No `path`, uses `filename`; no spk attr -> skip segments (1878, 1885->1908)."""

    class _Planets:
        # No `path` attribute, but a `filename`.
        filename = "/some/de440.bsp"

    monkeypatch.setattr(state, "_PLANETS", _Planets())
    monkeypatch.setattr(state, "_EPHEMERIS_FILE", "de440.bsp")
    path, start, end, denum = state.get_current_file_data(0)
    assert path == "/some/de440.bsp"
    assert start == 0.0
    assert end == 0.0
    assert denum == 440


def test_get_current_file_data_empty_segments(monkeypatch):
    """spk.segments empty -> skip (1887->1908)."""

    class _Spk:
        segments = []

    class _Planets:
        path = "/x/de441.bsp"
        spk = _Spk()

    monkeypatch.setattr(state, "_PLANETS", _Planets())
    monkeypatch.setattr(state, "_EPHEMERIS_FILE", "de441.bsp")
    path, start, end, denum = state.get_current_file_data(0)
    assert start == 0.0
    assert end == 0.0
    assert denum == 441


def test_get_current_file_data_segments_with_dates(monkeypatch):
    """Segments with start/end -> min/max computed; inner-if False arcs.

    The widest segment comes first so the later, narrower segment leaves
    both running extremes unchanged, exercising the ``if seg_start <
    start_jd`` False arc (1895->1897) and the ``if seg_end > end_jd``
    False arc that loops to the next iteration (1899->1892).
    """

    class _SegWide:
        start_jd = 2300000.0
        end_jd = 2500000.0

    class _SegNarrow:
        # Fully inside the wide segment: start not smaller, end not larger.
        start_jd = 2400000.0
        end_jd = 2450000.0

    class _Spk:
        segments = [_SegWide(), _SegNarrow()]

    class _Planets:
        path = "/x/de440.bsp"
        spk = _Spk()

    monkeypatch.setattr(state, "_PLANETS", _Planets())
    monkeypatch.setattr(state, "_EPHEMERIS_FILE", "de440.bsp")
    path, start, end, denum = state.get_current_file_data(0)
    assert start == 2300000.0
    assert end == 2500000.0


def test_get_current_file_data_segments_no_valid_dates(monkeypatch):
    """Segments present but no start/end attrs -> fallback 0.0 (1903, 1905)."""

    class _SegEmpty:
        pass

    class _Spk:
        segments = [_SegEmpty()]

    class _Planets:
        path = "/x/de440.bsp"
        spk = _Spk()

    monkeypatch.setattr(state, "_PLANETS", _Planets())
    monkeypatch.setattr(state, "_EPHEMERIS_FILE", "de440.bsp")
    path, start, end, denum = state.get_current_file_data(0)
    assert start == 0.0
    assert end == 0.0


def test_get_current_file_data_non_de_filename(monkeypatch):
    """Filename not de*.bsp -> denum stays 0 (1910->1918)."""

    class _Planets:
        path = "/x/custom.bsp"

    monkeypatch.setattr(state, "_PLANETS", _Planets())
    monkeypatch.setattr(state, "_EPHEMERIS_FILE", "custom.bsp")
    path, start, end, denum = state.get_current_file_data(0)
    assert denum == 0


def test_get_current_file_data_bad_de_number(monkeypatch):
    """de<garbage>.bsp -> int parse fails -> denum 0 (1915-1916)."""

    class _Planets:
        path = "/x/dexx.bsp"

    monkeypatch.setattr(state, "_PLANETS", _Planets())
    monkeypatch.setattr(state, "_EPHEMERIS_FILE", "dexx.bsp")
    path, start, end, denum = state.get_current_file_data(0)
    assert denum == 0


def test_get_current_file_data_outer_exception(monkeypatch):
    """Accessing _PLANETS.path raises -> outer except returns empty (1920-1922)."""

    class _Planets:
        @property
        def path(self):
            raise ValueError("boom path")

    monkeypatch.setattr(state, "_PLANETS", _Planets())
    monkeypatch.setattr(state, "_EPHEMERIS_FILE", "de440.bsp")
    assert state.get_current_file_data(0) == ("", 0.0, 0.0, 0)


# ===========================================================================
# get_auto_spk_download TOML  (line 2010)
# ===========================================================================


def test_get_auto_spk_download_toml(monkeypatch):
    """TOML returns bool -> returned (line 2010)."""
    monkeypatch.setattr(state, "_AUTO_SPK_DOWNLOAD", None)
    monkeypatch.delenv(state._AUTO_SPK_ENV_VAR, raising=False)
    monkeypatch.setattr(
        "libephemeris._config_toml.get_bool",
        lambda key: False if key == "auto_spk" else None,
    )
    assert state.get_auto_spk_download() is False


# ===========================================================================
# get_spk_cache_dir  (lines 2086, 2092)
# ===========================================================================


def test_get_spk_cache_dir_env(monkeypatch):
    """Env var set -> returned (line 2086)."""
    monkeypatch.setattr(state, "_SPK_CACHE_DIR", None)
    monkeypatch.setenv(state._SPK_CACHE_DIR_ENV_VAR, "/env/cache")
    assert state.get_spk_cache_dir() == "/env/cache"


def test_get_spk_cache_dir_toml(monkeypatch):
    """TOML returns dir -> returned (line 2092)."""
    monkeypatch.setattr(state, "_SPK_CACHE_DIR", None)
    monkeypatch.delenv(state._SPK_CACHE_DIR_ENV_VAR, raising=False)
    monkeypatch.setattr(
        "libephemeris._config_toml.get_str",
        lambda key: "/toml/cache" if key == "spk_dir" else None,
    )
    assert state.get_spk_cache_dir() == "/toml/cache"


# ===========================================================================
# get_iers_delta_t_enabled  (lines 2216, 2218, 2225)
# ===========================================================================


def test_get_iers_delta_t_enabled_env_true(monkeypatch):
    """Env "1" -> True (line 2216)."""
    monkeypatch.setattr(state, "_IERS_DELTA_T_ENABLED", None)
    monkeypatch.setenv(state._IERS_DELTA_T_ENV_VAR, "1")
    assert state.get_iers_delta_t_enabled() is True


def test_get_iers_delta_t_enabled_env_false(monkeypatch):
    """Env "0" -> False (line 2218)."""
    monkeypatch.setattr(state, "_IERS_DELTA_T_ENABLED", None)
    monkeypatch.setenv(state._IERS_DELTA_T_ENV_VAR, "0")
    assert state.get_iers_delta_t_enabled() is False


def test_get_iers_delta_t_enabled_toml(monkeypatch):
    """TOML returns bool -> returned (line 2225)."""
    monkeypatch.setattr(state, "_IERS_DELTA_T_ENABLED", None)
    monkeypatch.delenv(state._IERS_DELTA_T_ENV_VAR, raising=False)
    monkeypatch.setattr(
        "libephemeris._config_toml.get_bool",
        lambda key: True if key == "iers_delta_t" else None,
    )
    assert state.get_iers_delta_t_enabled() is True


# ===========================================================================
# get_strict_precision TOML  (line 2315)
# ===========================================================================


def test_get_strict_precision_toml(monkeypatch):
    """TOML returns bool -> returned (line 2315)."""
    monkeypatch.setattr(state, "_STRICT_PRECISION", None)
    monkeypatch.delenv(state._STRICT_PRECISION_ENV_VAR, raising=False)
    monkeypatch.setattr(
        "libephemeris._config_toml.get_bool",
        lambda key: False if key == "strict_precision" else None,
    )
    assert state.get_strict_precision() is False


# ===========================================================================
# _load_spk_kernel  (lines 2350, 2353, 2359-2360)
# ===========================================================================


def test_load_spk_kernel_cached(monkeypatch):
    """Already-cached path -> return cache (line 2350)."""
    sentinel = object()
    monkeypatch.setattr(state, "_SPK_KERNELS", {"/p.bsp": sentinel})
    assert state._load_spk_kernel("/p.bsp") is sentinel


def test_load_spk_kernel_missing(monkeypatch):
    """Missing file -> FileNotFoundError (line 2353)."""
    monkeypatch.setattr(state, "_SPK_KERNELS", {})
    monkeypatch.setattr(state.os.path, "exists", lambda p: False)
    with pytest.raises(FileNotFoundError):
        state._load_spk_kernel("/missing.bsp")


def test_load_spk_kernel_loads_and_caches(monkeypatch, tmp_path):
    """Valid file -> load, cache, return (lines 2359-2360)."""
    monkeypatch.setattr(state, "_SPK_KERNELS", {})
    f = tmp_path / "body.bsp"
    f.write_bytes(b"x")
    fake_loader = _FakeLoader(kernel="spk-kernel")
    monkeypatch.setattr(state, "get_loader", lambda: fake_loader)
    result = state._load_spk_kernel(str(f))
    assert result == "spk-kernel"
    assert state._SPK_KERNELS[str(f)] == "spk-kernel"


def test_load_spk_kernel_load_raises(monkeypatch, tmp_path):
    """load() raises -> wrapped ValueError (covers except branch)."""
    monkeypatch.setattr(state, "_SPK_KERNELS", {})
    f = tmp_path / "body.bsp"
    f.write_bytes(b"x")

    class _BoomLoader:
        def __call__(self, path):
            raise OSError("bad")

    monkeypatch.setattr(state, "get_loader", lambda: _BoomLoader())
    with pytest.raises(ValueError):
        state._load_spk_kernel(str(f))


# ===========================================================================
# _get_spk_target  (lines 2377-2389)
# ===========================================================================


def test_get_spk_target_not_registered(monkeypatch):
    """ipl not in body map -> None (line 2377-2378)."""
    monkeypatch.setattr(state, "_SPK_BODY_MAP", {})
    assert state._get_spk_target(15) is None


def test_get_spk_target_kernel_missing(monkeypatch):
    """Body mapped but kernel not loaded -> None (lines 2380-2384)."""
    monkeypatch.setattr(state, "_SPK_BODY_MAP", {15: ("/k.bsp", 2060)})
    monkeypatch.setattr(state, "_SPK_KERNELS", {})
    assert state._get_spk_target(15) is None


def test_get_spk_target_success(monkeypatch):
    """Body mapped + kernel present -> kernel[naif_id] (lines 2386-2387)."""

    class _Kernel:
        def __getitem__(self, naif_id):
            assert naif_id == 2060
            return "target-obj"

    monkeypatch.setattr(state, "_SPK_BODY_MAP", {15: ("/k.bsp", 2060)})
    monkeypatch.setattr(state, "_SPK_KERNELS", {"/k.bsp": _Kernel()})
    assert state._get_spk_target(15) == "target-obj"


def test_get_spk_target_keyerror(monkeypatch):
    """kernel[naif_id] raises KeyError -> None (lines 2388-2389)."""

    class _Kernel:
        def __getitem__(self, naif_id):
            raise KeyError(naif_id)

    monkeypatch.setattr(state, "_SPK_BODY_MAP", {15: ("/k.bsp", 2060)})
    monkeypatch.setattr(state, "_SPK_KERNELS", {"/k.bsp": _Kernel()})
    assert state._get_spk_target(15) is None


def test_current_file_data_reports_active_leb():
    """In sealed LEB mode the metadata reports the active LEB artifact.

    The path names the serving core file, the range is the stored union
    coverage of the requested channel, and the DE number is the generating
    kernel's (440 for base/medium tiers, 441 for extended). ifno 2-4 stay
    empty like the JPL path.
    """
    import libephemeris as le

    prev = le.get_calc_mode()
    try:
        le.set_calc_mode("leb")
        le.calc_ut(2451545.0, le.SUN, 0)
        path, start, end, denum = le.get_current_file_data(0)
        assert path.endswith(".leb2") or path.endswith(".leb")
        assert start < 2451545.0 < end
        assert denum in (440, 441)
        assert le.get_current_file_data(2) == ("", 0.0, 0.0, 0)
    finally:
        le.set_calc_mode(prev)
