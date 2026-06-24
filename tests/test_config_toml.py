# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Tests for the TOML configuration loader (libephemeris/_config_toml.py)."""

from __future__ import annotations

import pytest

from libephemeris import _config_toml as cfg

_GOOD_TOML = """
[libephemeris]
precision = "extended"
mode = "leb"
auto_spk = true
strict_precision = false
deltat_model = "espenak_meeus"
mmap_preload_start = 2451545
log_level = "DEBUG"
unknown_key = "ignored"
data_dir = 999
"""


@pytest.fixture(autouse=True)
def _reset_config():
    cfg.reset()
    yield
    cfg.reset()


def _write(tmp_path, text, name="libephemeris-config.toml"):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


class TestFindConfigFile:
    def test_env_override_existing(self, tmp_path, monkeypatch):
        p = _write(tmp_path, _GOOD_TOML)
        monkeypatch.setenv(cfg.CONFIG_FILE_VAR, str(p))
        assert cfg._find_config_file() == p

    def test_env_override_missing_returns_none(self, tmp_path, monkeypatch):
        monkeypatch.setenv(cfg.CONFIG_FILE_VAR, str(tmp_path / "nope.toml"))
        assert cfg._find_config_file() is None

    def test_search_path_hit(self, tmp_path, monkeypatch):
        p = _write(tmp_path, _GOOD_TOML)
        monkeypatch.delenv(cfg.CONFIG_FILE_VAR, raising=False)
        monkeypatch.setattr(cfg, "_SEARCH_PATHS", (lambda: p,))
        assert cfg._find_config_file() == p

    def test_search_path_raising_is_skipped(self, monkeypatch):
        def _boom():
            raise OSError("cwd unavailable")

        monkeypatch.delenv(cfg.CONFIG_FILE_VAR, raising=False)
        monkeypatch.setattr(cfg, "_SEARCH_PATHS", (_boom,))
        assert cfg._find_config_file() is None


class TestLoadConfig:
    def test_no_toml_parser(self, monkeypatch):
        monkeypatch.setattr(cfg, "tomllib", None)
        assert cfg.load_config() is False
        assert cfg._CONFIG_LOADED is True

    def test_explicit_missing_path(self, tmp_path):
        assert cfg.load_config(str(tmp_path / "absent.toml")) is False

    def test_no_file_found(self, monkeypatch):
        monkeypatch.delenv(cfg.CONFIG_FILE_VAR, raising=False)
        monkeypatch.setattr(cfg, "_SEARCH_PATHS", ())
        assert cfg.load_config() is False

    def test_valid_file_loads_and_validates_types(self, tmp_path):
        p = _write(tmp_path, _GOOD_TOML)
        assert cfg.load_config(str(p)) is True
        # Valid keys with correct types are kept...
        assert cfg.get_str("precision") == "extended"
        assert cfg.get_str("mode") == "leb"
        assert cfg.get_str("deltat_model") == "espenak_meeus"
        assert cfg.get_bool("auto_spk") is True
        assert cfg.get_bool("strict_precision") is False
        assert cfg.get_int("mmap_preload_start") == 2451545
        # ...unknown keys (unknown_key) and valid keys with the wrong type
        # (data_dir expects str but is an int) are dropped.
        assert cfg.get_str("unknown_key") is None
        assert cfg.get_str("data_dir") is None

    def test_load_via_search_path(self, tmp_path, monkeypatch):
        # No explicit path -> discovery returns a file -> the load proceeds.
        p = _write(tmp_path, _GOOD_TOML)
        monkeypatch.delenv(cfg.CONFIG_FILE_VAR, raising=False)
        monkeypatch.setattr(cfg, "_SEARCH_PATHS", (lambda: p,))
        assert cfg.load_config() is True
        assert cfg.get_str("mode") == "leb"

    def test_malformed_file_warns_and_returns_false(self, tmp_path):
        p = _write(tmp_path, "this is = = not valid toml [[[")
        with pytest.warns(UserWarning, match="malformed libephemeris config"):
            assert cfg.load_config(str(p)) is False

    def test_open_error_returns_false(self, tmp_path, monkeypatch):
        p = _write(tmp_path, _GOOD_TOML)

        def _raise_open(*_a, **_k):
            raise OSError("permission denied")

        monkeypatch.setattr("builtins.open", _raise_open)
        assert cfg.load_config(str(p)) is False

    def test_section_not_a_dict_returns_false(self, tmp_path):
        p = _write(tmp_path, "libephemeris = 42\n")
        assert cfg.load_config(str(p)) is False


class TestAccessors:
    def test_typed_getters_with_no_config(self):
        # _ensure_loaded triggers a (no-op) load; absent keys -> None.
        assert cfg.get_str("precision") is None
        assert cfg.get_bool("auto_spk") is None
        assert cfg.get_int("mmap_preload_start") is None

    def test_get_int_rejects_bool(self, tmp_path):
        p = _write(tmp_path, "[libephemeris]\nauto_spk = true\n")
        cfg.load_config(str(p))
        # bool is a subclass of int but must not be returned as int.
        assert cfg.get_int("auto_spk") is None

    def test_introspection(self, tmp_path):
        p = _write(tmp_path, _GOOD_TOML)
        cfg.load_config(str(p))
        assert cfg.get_config_path() == str(p)
        assert cfg.get_all()["precision"] == "extended"
        assert cfg.is_loaded() is True

    def test_is_loaded_false_when_empty(self, monkeypatch):
        monkeypatch.delenv(cfg.CONFIG_FILE_VAR, raising=False)
        monkeypatch.setattr(cfg, "_SEARCH_PATHS", ())
        assert cfg.is_loaded() is False
