# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for ``libephemeris.spk_auto``.

This module pushes ``spk_auto`` to full branch coverage using mocks for all
network/file-IO/optional-dependency seams. No real network access is performed:
every call into the ``spk`` download helpers, the ``state`` configuration
accessors and the optional ``astroquery`` import is mocked.
"""

from __future__ import annotations

import os
from unittest.mock import MagicMock, patch

import pytest

from libephemeris import spk_auto
from libephemeris.constants import CHIRON, NAIF_CHIRON
from libephemeris.exceptions import Error


# =============================================================================
# _detect_naif_id_from_file (lines 61-70)
# =============================================================================


class TestDetectNaifIdFromFile:
    """Cover every branch of ``_detect_naif_id_from_file``."""

    def test_targets_lookup_raises(self):
        """Wrapped ``_get_spk_targets`` raising returns None (lines 61-62)."""
        with patch(
            "libephemeris.spk._get_spk_targets",
            side_effect=OSError("boom"),
        ):
            assert spk_auto._detect_naif_id_from_file("/x.bsp") is None

    def test_empty_targets_returns_none(self):
        """Empty target list returns None (line 64)."""
        with patch("libephemeris.spk._get_spk_targets", return_value=[]):
            assert spk_auto._detect_naif_id_from_file("/x.bsp") is None

    def test_preferred_in_targets(self):
        """Preferred NAIF present in file is returned (lines 65-66)."""
        with patch(
            "libephemeris.spk._get_spk_targets",
            return_value=[2002060, 10],
        ):
            assert spk_auto._detect_naif_id_from_file("/x.bsp", 2002060) == 2002060

    def test_single_small_body(self):
        """Single small-body target is returned (lines 67-69)."""
        with patch(
            "libephemeris.spk._get_spk_targets",
            return_value=[10, 2002060],
        ):
            assert spk_auto._detect_naif_id_from_file("/x.bsp") == 2002060

    def test_multiple_small_bodies_returns_none(self):
        """Multiple small-body targets are ambiguous -> None (line 70)."""
        with patch(
            "libephemeris.spk._get_spk_targets",
            return_value=[2002060, 2002061],
        ):
            assert spk_auto._detect_naif_id_from_file("/x.bsp") is None


# =============================================================================
# Default-directory fallthroughs (lines 128, 167, 207)
# =============================================================================


class TestDefaultDirectoryFallbacks:
    """Cover ``cache_dir is None`` -> DEFAULT_AUTO_SPK_DIR branches."""

    def test_ensure_cache_dir_default(self, tmp_path, monkeypatch):
        """ensure_cache_dir uses DEFAULT when state returns None (line 128)."""
        default = str(tmp_path / "def")
        monkeypatch.setattr(spk_auto, "DEFAULT_AUTO_SPK_DIR", default)
        with patch("libephemeris.state.get_spk_cache_dir", return_value=None):
            result = spk_auto.ensure_cache_dir()
        assert result == os.path.abspath(default)
        assert os.path.exists(default)

    def test_get_cache_path_default(self, monkeypatch):
        """get_cache_path uses DEFAULT when cache_dir is None (line 167)."""
        monkeypatch.setattr(spk_auto, "DEFAULT_AUTO_SPK_DIR", "/tmp/cov_def")
        result = spk_auto.get_cache_path("2060")
        assert result == os.path.join("/tmp/cov_def", "2060.bsp")

    def test_cache_info_default_dir(self, tmp_path, monkeypatch):
        """cache_info uses DEFAULT when cache_dir is None (line 207)."""
        default = str(tmp_path / "nope")
        monkeypatch.setattr(spk_auto, "DEFAULT_AUTO_SPK_DIR", default)
        info = spk_auto.cache_info()
        assert info["cache_dir"] == default
        assert info["num_files"] == 0


# =============================================================================
# AutoSpkConfig.get_cache_dir / get_cache_path (lines 267-280, 297)
# =============================================================================


class TestAutoSpkConfigPaths:
    """Cover the three ``get_cache_dir`` branches and ``get_cache_path``."""

    def test_explicit_cache_dir(self, tmp_path):
        """Explicit cache_dir is used and created (lines 270-271, 277-278)."""
        target = str(tmp_path / "explicit")
        config = spk_auto.AutoSpkConfig(ipl=CHIRON, body_id="2060", cache_dir=target)
        assert config.get_cache_dir() == target
        assert os.path.exists(target)

    def test_state_cache_dir(self, tmp_path):
        """State-provided cache_dir is used (lines 272-273)."""
        target = str(tmp_path / "state")
        config = spk_auto.AutoSpkConfig(ipl=CHIRON, body_id="2060")
        with patch("libephemeris.state.get_spk_cache_dir", return_value=target):
            assert config.get_cache_dir() == target
        assert os.path.exists(target)

    def test_default_cache_dir_existing(self, tmp_path, monkeypatch):
        """Default cache_dir branch, directory already exists (lines 275, 277)."""
        default = str(tmp_path / "default")
        os.makedirs(default, exist_ok=True)
        monkeypatch.setattr(spk_auto, "DEFAULT_AUTO_SPK_DIR", default)
        config = spk_auto.AutoSpkConfig(ipl=CHIRON, body_id="2060")
        with patch("libephemeris.state.get_spk_cache_dir", return_value=None):
            assert config.get_cache_dir() == default

    def test_get_cache_path_method(self, tmp_path):
        """get_cache_path joins dir + filename (line 297)."""
        target = str(tmp_path / "cp")
        config = spk_auto.AutoSpkConfig(ipl=CHIRON, body_id="2060", cache_dir=target)
        path = config.get_cache_path()
        assert path.startswith(target)
        assert path.endswith(".bsp")


# =============================================================================
# _download_spk_astroquery (lines 342, 361->364, 366-367)
# =============================================================================


class TestDownloadSpkAstroquery:
    """Cover the download wrapper success and failure paths."""

    def test_empty_output_dir_uses_dot(self, tmp_path, monkeypatch):
        """Bare filename -> output_dir defaults to '.' (line 342)."""
        monkeypatch.chdir(tmp_path)
        produced = tmp_path / "out.bsp"

        def fake_download(**kwargs):
            produced.write_bytes(b"x")
            return str(produced)

        with patch("libephemeris.spk.download_spk", side_effect=fake_download):
            result = spk_auto._download_spk_astroquery(
                body_id="2060", start="2020-01-01", end="2030-01-01",
                output_path="out.bsp",
            )
        assert result == "out.bsp"

    def test_no_rename_when_same_path(self, tmp_path):
        """result_path == output_path skips os.replace (branch 361->364)."""
        out = str(tmp_path / "same.bsp")

        def fake_download(**kwargs):
            with open(out, "wb") as f:
                f.write(b"x")
            return out

        with patch("libephemeris.spk.download_spk", side_effect=fake_download):
            result = spk_auto._download_spk_astroquery(
                body_id="2060", start="2020-01-01", end="2030-01-01",
                output_path=out,
            )
        assert result == out

    def test_download_failure_wrapped(self, tmp_path):
        """download_spk raising is wrapped in ValueError (lines 366-367)."""
        out = str(tmp_path / "fail.bsp")
        with patch(
            "libephemeris.spk.download_spk",
            side_effect=KeyError("nope"),
        ):
            with pytest.raises(ValueError, match="Failed to download SPK"):
                spk_auto._download_spk_astroquery(
                    body_id="2060", start="2020-01-01", end="2030-01-01",
                    output_path=out,
                )


# =============================================================================
# _check_astroquery_available False branch (lines 378-379)
# =============================================================================


class TestCheckAstroquery:
    """Cover the ImportError arm of the availability probe."""

    def test_import_error_returns_false(self):
        """ImportError on astroquery import -> False (lines 378-379)."""
        real_import = __import__

        def fake_import(name, *args, **kwargs):
            if name.startswith("astroquery"):
                raise ImportError("no astroquery")
            return real_import(name, *args, **kwargs)

        with patch.dict("sys.modules", {"astroquery.jplhorizons": None}):
            with patch("builtins.__import__", side_effect=fake_import):
                assert spk_auto._check_astroquery_available() is False

    def test_available_returns_true(self):
        """Present module -> True (lines 374-377)."""
        fake_mod = MagicMock(Horizons=MagicMock())
        with patch.dict("sys.modules", {"astroquery.jplhorizons": fake_mod}):
            assert spk_auto._check_astroquery_available() is True


# =============================================================================
# enable_auto_spk cache_dir resolution + disable_auto_spk (433->436, 460->exit)
# =============================================================================


class TestEnableDisable:
    """Cover the cache_dir fallback in enable and the no-op disable branch."""

    def setup_method(self):
        spk_auto.disable_all()

    def teardown_method(self):
        spk_auto.disable_all()

    def test_enable_resolves_cache_dir_from_state(self):
        """cache_dir None pulls from state (branch 433->436)."""
        with patch(
            "libephemeris.state.get_spk_cache_dir", return_value="/some/dir"
        ):
            spk_auto.enable_auto_spk(ipl=CHIRON, body_id="2060")
        cfg = spk_auto.get_auto_spk_config(CHIRON)
        assert cfg is not None
        assert cfg.cache_dir == "/some/dir"

    def test_disable_unconfigured_is_noop(self):
        """Disabling an unregistered ipl falls through (branch 460->exit)."""
        spk_auto.disable_auto_spk(99999)
        assert spk_auto.get_auto_spk_config(99999) is None


# =============================================================================
# download_now / _ensure_spk_downloaded (lines 527, 543-597)
# =============================================================================


class TestEnsureSpkDownloaded:
    """Cover the full download/register pipeline of _ensure_spk_downloaded."""

    def setup_method(self):
        spk_auto.disable_all()
        from libephemeris import state

        state._SPK_BODY_MAP.clear()

    def teardown_method(self):
        spk_auto.disable_all()
        from libephemeris import state

        state._SPK_BODY_MAP.clear()

    def test_download_now_triggers_download(self, tmp_path):
        """download_now reaches _ensure_spk_downloaded (line 527)."""
        config = spk_auto.AutoSpkConfig(
            ipl=CHIRON, body_id="2060", naif_id=NAIF_CHIRON,
            cache_dir=str(tmp_path),
        )
        spk_auto.enable_auto_spk(
            ipl=CHIRON, body_id="2060", naif_id=NAIF_CHIRON,
            cache_dir=str(tmp_path),
        )

        def fake_download(**kwargs):
            with open(kwargs["output_path"], "wb") as f:
                f.write(b"x")

        with patch.object(
            spk_auto, "_download_spk_astroquery", side_effect=fake_download
        ), patch("libephemeris.spk._get_spk_targets", return_value=[NAIF_CHIRON]), \
                patch("libephemeris.spk.register_spk_body") as mock_reg:
            result = spk_auto.download_now(CHIRON)
        assert result == config.get_cache_path() or result.endswith(".bsp")
        mock_reg.assert_called_once()

    def test_cached_valid_file_used(self, tmp_path):
        """Cached + valid file is reused without download (lines 549-551)."""
        config = spk_auto.AutoSpkConfig(
            ipl=CHIRON, body_id="2060", naif_id=NAIF_CHIRON,
            cache_dir=str(tmp_path),
        )
        cache_path = config.get_cache_path()
        with open(cache_path, "wb") as f:
            f.write(b"x")

        with patch.object(spk_auto, "_is_valid_bsp", return_value=True), \
                patch.object(spk_auto, "_download_spk_astroquery") as mock_dl, \
                patch(
                    "libephemeris.spk._get_spk_targets", return_value=[NAIF_CHIRON]
                ), patch("libephemeris.spk.register_spk_body"):
            result = spk_auto._ensure_spk_downloaded(config)
        assert result == cache_path
        assert config.spk_path == cache_path
        mock_dl.assert_not_called()

    def test_corrupted_cache_redownload(self, tmp_path):
        """Corrupted cached file is removed + re-downloaded (lines 552-560)."""
        config = spk_auto.AutoSpkConfig(
            ipl=CHIRON, body_id="2060", naif_id=NAIF_CHIRON,
            cache_dir=str(tmp_path),
        )
        cache_path = config.get_cache_path()
        with open(cache_path, "wb") as f:
            f.write(b"corrupt")

        def fake_download(**kwargs):
            with open(kwargs["output_path"], "wb") as f:
                f.write(b"fresh")

        with patch.object(spk_auto, "_is_valid_bsp", return_value=False), \
                patch.object(
                    spk_auto, "_download_spk_astroquery", side_effect=fake_download
                ) as mock_dl, \
                patch(
                    "libephemeris.spk._get_spk_targets", return_value=[NAIF_CHIRON]
                ), patch("libephemeris.spk.register_spk_body"):
            result = spk_auto._ensure_spk_downloaded(config)
        assert result == cache_path
        mock_dl.assert_called_once()

    def test_corrupted_cache_remove_oserror(self, tmp_path):
        """os.remove failure on corrupted cache is swallowed (lines 558-559)."""
        config = spk_auto.AutoSpkConfig(
            ipl=CHIRON, body_id="2060", naif_id=NAIF_CHIRON,
            cache_dir=str(tmp_path),
        )
        cache_path = config.get_cache_path()
        with open(cache_path, "wb") as f:
            f.write(b"corrupt")

        def fake_download(**kwargs):
            with open(kwargs["output_path"], "wb") as f:
                f.write(b"fresh")

        with patch.object(spk_auto, "_is_valid_bsp", return_value=False), \
                patch("os.remove", side_effect=OSError("locked")), \
                patch.object(
                    spk_auto, "_download_spk_astroquery", side_effect=fake_download
                ), \
                patch(
                    "libephemeris.spk._get_spk_targets", return_value=[NAIF_CHIRON]
                ), patch("libephemeris.spk.register_spk_body"):
            result = spk_auto._ensure_spk_downloaded(config)
        assert result == cache_path

    def test_naif_id_detected_from_file(self, tmp_path):
        """Detected NAIF overrides config (lines 578-581)."""
        config = spk_auto.AutoSpkConfig(
            ipl=CHIRON, body_id="2060", naif_id=None,
            cache_dir=str(tmp_path),
        )

        def fake_download(**kwargs):
            with open(kwargs["output_path"], "wb") as f:
                f.write(b"x")

        with patch.object(
            spk_auto, "_download_spk_astroquery", side_effect=fake_download
        ), patch(
            "libephemeris.spk._get_spk_targets", return_value=[20002060]
        ), patch("libephemeris.spk.register_spk_body") as mock_reg:
            spk_auto._ensure_spk_downloaded(config)
        assert config.naif_id == 20002060
        mock_reg.assert_called_once()

    def test_naif_id_deduced(self, tmp_path):
        """NAIF deduced from body_id when none detected (lines 582-583)."""
        config = spk_auto.AutoSpkConfig(
            ipl=CHIRON, body_id="2060", naif_id=None,
            cache_dir=str(tmp_path),
        )

        def fake_download(**kwargs):
            with open(kwargs["output_path"], "wb") as f:
                f.write(b"x")

        with patch.object(
            spk_auto, "_download_spk_astroquery", side_effect=fake_download
        ), patch.object(
            spk_auto, "_detect_naif_id_from_file", return_value=None
        ), patch("libephemeris.spk.register_spk_body") as mock_reg:
            spk_auto._ensure_spk_downloaded(config)
        assert config.naif_id == NAIF_CHIRON
        mock_reg.assert_called_once()

    def test_naif_id_cannot_deduce_raises(self, tmp_path):
        """Undeducible NAIF raises ValueError (lines 584-588)."""
        config = spk_auto.AutoSpkConfig(
            ipl=CHIRON, body_id="NoNumber", naif_id=None,
            cache_dir=str(tmp_path),
        )

        def fake_download(**kwargs):
            with open(kwargs["output_path"], "wb") as f:
                f.write(b"x")

        with patch.object(
            spk_auto, "_download_spk_astroquery", side_effect=fake_download
        ), patch.object(
            spk_auto, "_detect_naif_id_from_file", return_value=None
        ):
            with pytest.raises(ValueError, match="Cannot deduce NAIF ID"):
                spk_auto._ensure_spk_downloaded(config)

    def test_inlock_recheck_finds_file(self, tmp_path):
        """File appearing between checks skips download (branch 566->573)."""
        config = spk_auto.AutoSpkConfig(
            ipl=CHIRON, body_id="2060", naif_id=NAIF_CHIRON,
            cache_dir=str(tmp_path),
        )
        cache_path = config.get_cache_path()

        real_exists = os.path.exists
        state = {"count": 0}

        def exists_side(p):
            if p == cache_path:
                state["count"] += 1
                # First check (outer + cached-check) misses, then file "appears"
                # for the in-lock re-check at line 566.
                return state["count"] > 2
            return real_exists(p)

        with patch("os.path.exists", side_effect=exists_side), \
                patch.object(spk_auto, "_download_spk_astroquery") as mock_dl, \
                patch(
                    "libephemeris.spk._get_spk_targets", return_value=[NAIF_CHIRON]
                ), patch("libephemeris.spk.register_spk_body"):
            result = spk_auto._ensure_spk_downloaded(config, force=False)
        assert result == cache_path
        mock_dl.assert_not_called()

    def test_register_skipped_for_same_path_rereg_for_different(self, tmp_path):
        """register_spk_body is skipped only when the SAME path is already
        registered; a different already-registered path is re-registered
        (branch 594) so a force re-download to a new path takes effect."""
        from libephemeris import state

        config = spk_auto.AutoSpkConfig(
            ipl=CHIRON, body_id="2060", naif_id=NAIF_CHIRON,
            cache_dir=str(tmp_path),
        )
        resolved = config.get_cache_path()

        def fake_download(**kwargs):
            with open(kwargs["output_path"], "wb") as f:
                f.write(b"x")

        try:
            # (a) Same path already registered -> register is skipped.
            state._SPK_BODY_MAP[CHIRON] = (resolved, NAIF_CHIRON)
            with patch.object(
                spk_auto, "_download_spk_astroquery", side_effect=fake_download
            ), patch(
                "libephemeris.spk._get_spk_targets", return_value=[NAIF_CHIRON]
            ), patch("libephemeris.spk.register_spk_body") as mock_reg:
                spk_auto._ensure_spk_downloaded(config)
            mock_reg.assert_not_called()

            # (b) A different (stale) path registered -> re-register.
            state._SPK_BODY_MAP[CHIRON] = ("stale_old.bsp", NAIF_CHIRON)
            with patch.object(
                spk_auto, "_download_spk_astroquery", side_effect=fake_download
            ), patch(
                "libephemeris.spk._get_spk_targets", return_value=[NAIF_CHIRON]
            ), patch("libephemeris.spk.register_spk_body") as mock_reg:
                spk_auto._ensure_spk_downloaded(config)
            mock_reg.assert_called_once()
        finally:
            state._SPK_BODY_MAP.clear()


# =============================================================================
# try_auto_download (lines 617-622)
# =============================================================================


class TestTryAutoDownload:
    """Cover the success and except arms of try_auto_download."""

    def setup_method(self):
        spk_auto.disable_all()

    def teardown_method(self):
        spk_auto.disable_all()

    def test_success(self):
        """Returns ensure result on success (lines 617-618)."""
        cfg = spk_auto.AutoSpkConfig(ipl=CHIRON, body_id="2060")
        with patch.object(spk_auto, "get_auto_spk_config", return_value=cfg), \
                patch.object(
                    spk_auto, "_ensure_spk_downloaded", return_value="/p.bsp"
                ):
            assert spk_auto.try_auto_download(CHIRON) == "/p.bsp"

    def test_failure_returns_none(self):
        """Exception in ensure is swallowed -> None (lines 619-622)."""
        cfg = spk_auto.AutoSpkConfig(ipl=CHIRON, body_id="2060")
        with patch.object(spk_auto, "get_auto_spk_config", return_value=cfg), \
                patch.object(
                    spk_auto, "_ensure_spk_downloaded", side_effect=Error("x")
                ):
            assert spk_auto.try_auto_download(CHIRON) is None

    def test_disabled_returns_none(self):
        """Disabled config returns None (line 615)."""
        cfg = spk_auto.AutoSpkConfig(ipl=CHIRON, body_id="2060")
        cfg.enabled = False
        with patch.object(spk_auto, "get_auto_spk_config", return_value=cfg):
            assert spk_auto.try_auto_download(CHIRON) is None


# =============================================================================
# clear_cache with files (lines 642-654)
# =============================================================================


class TestClearCache:
    """Cover clear_cache removal of per-body and all-body files."""

    def setup_method(self):
        spk_auto.disable_all()

    def teardown_method(self):
        spk_auto.disable_all()

    def test_clear_single_body(self, tmp_path):
        """Removes the spk_path for a single body (lines 642-644)."""
        f = tmp_path / "chiron.bsp"
        f.write_bytes(b"x")
        cfg = spk_auto.AutoSpkConfig(ipl=CHIRON, body_id="2060")
        cfg.spk_path = str(f)
        spk_auto._AUTO_SPK_REGISTRY[CHIRON] = cfg
        deleted = spk_auto.clear_cache(CHIRON)
        assert deleted == 1
        assert not f.exists()
        assert cfg.spk_path is None

    def test_clear_all_bodies(self, tmp_path):
        """Removes all registered files; OSError swallowed (lines 646-654)."""
        f1 = tmp_path / "a.bsp"
        f1.write_bytes(b"x")
        cfg1 = spk_auto.AutoSpkConfig(ipl=CHIRON, body_id="2060")
        cfg1.spk_path = str(f1)
        spk_auto._AUTO_SPK_REGISTRY[CHIRON] = cfg1

        f2 = tmp_path / "b.bsp"
        f2.write_bytes(b"x")
        cfg2 = spk_auto.AutoSpkConfig(ipl=99998, body_id="999")
        cfg2.spk_path = str(f2)
        spk_auto._AUTO_SPK_REGISTRY[99998] = cfg2

        # A third config with no spk_path exercises the False arc (648->647).
        cfg3 = spk_auto.AutoSpkConfig(ipl=99997, body_id="998")
        cfg3.spk_path = None
        spk_auto._AUTO_SPK_REGISTRY[99997] = cfg3

        # f2 removal raises OSError to hit lines 653-654.
        orig_remove = os.remove

        def selective_remove(p):
            if p == str(f2):
                raise OSError("locked")
            return orig_remove(p)

        with patch("os.remove", side_effect=selective_remove):
            deleted = spk_auto.clear_cache()
        assert deleted == 1
        assert cfg1.spk_path is None


# =============================================================================
# get_cache_info with files (lines 680-683)
# =============================================================================


class TestGetCacheInfo:
    """Cover get_cache_info listing real files."""

    def test_with_files(self, tmp_path, monkeypatch):
        """Lists .bsp files and computes size (lines 680-683)."""
        monkeypatch.setattr(spk_auto, "DEFAULT_AUTO_SPK_DIR", str(tmp_path))
        (tmp_path / "x.bsp").write_bytes(b"x" * 2048)
        (tmp_path / "ignore.txt").write_bytes(b"x")
        info = spk_auto.get_cache_info()
        assert info["num_files"] == 1
        assert "x.bsp" in info["files"]


# =============================================================================
# list_cached_spk filename fallback parsing (783->794, 791-792)
# =============================================================================


class TestListCachedSpk:
    """Cover the filename-parse fallback branches of list_cached_spk."""

    def test_short_filename_no_jd(self, tmp_path):
        """Coverage read fails, filename has < 3 parts (branch 783->794)."""
        (tmp_path / "chiron.bsp").write_bytes(b"x")
        with patch(
            "libephemeris.spk.get_spk_coverage",
            side_effect=ValueError("bad"),
        ):
            result = spk_auto.list_cached_spk(str(tmp_path))
        assert len(result) == 1
        assert result[0]["jd_start"] is None

    def test_filename_jd_parse_error(self, tmp_path):
        """Coverage fails + non-numeric trailing parts (lines 791-792)."""
        (tmp_path / "chiron_aa_bb.bsp").write_bytes(b"x")
        with patch(
            "libephemeris.spk.get_spk_coverage",
            side_effect=ValueError("bad"),
        ):
            result = spk_auto.list_cached_spk(str(tmp_path))
        assert result[0]["jd_start"] is None

    def test_filename_jd_parse_success(self, tmp_path):
        """Coverage fails but filename JDs parse (lines 784-790)."""
        (tmp_path / "2060_2450000_2470000.bsp").write_bytes(b"x")
        with patch(
            "libephemeris.spk.get_spk_coverage",
            side_effect=ValueError("bad"),
        ):
            result = spk_auto.list_cached_spk(str(tmp_path))
        assert result[0]["jd_start"] == 2450000.0


# =============================================================================
# clear_spk_cache (lines 852-854, 859->858)
# =============================================================================


class TestClearSpkCache:
    """Cover the OSError-on-remove and registry-cleanup branches."""

    def setup_method(self):
        spk_auto.disable_all()

    def teardown_method(self):
        spk_auto.disable_all()

    def test_remove_oserror_and_registry_cleanup(self, tmp_path):
        """OSError on remove is swallowed; registry cleared (lines 852-859)."""
        good = tmp_path / "good.bsp"
        good.write_bytes(b"x")
        bad = tmp_path / "bad.bsp"
        bad.write_bytes(b"x")

        # Register a config whose spk_path will no longer exist.
        cfg = spk_auto.AutoSpkConfig(ipl=CHIRON, body_id="2060")
        cfg.spk_path = str(tmp_path / "removed.bsp")  # never created
        spk_auto._AUTO_SPK_REGISTRY[CHIRON] = cfg

        # A second config with spk_path=None exercises the False arc (859->858).
        cfg2 = spk_auto.AutoSpkConfig(ipl=99996, body_id="997")
        cfg2.spk_path = None
        spk_auto._AUTO_SPK_REGISTRY[99996] = cfg2

        orig_remove = os.remove

        def selective_remove(p):
            if p == str(bad):
                raise OSError("locked")
            return orig_remove(p)

        with patch("os.remove", side_effect=selective_remove):
            deleted = spk_auto.clear_spk_cache(str(tmp_path))
        assert deleted == 1
        assert cfg.spk_path is None


# =============================================================================
# get_cache_size OSError (lines 915-916)
# =============================================================================


class TestGetCacheSize:
    """Cover the OSError arm of get_cache_size."""

    def test_getsize_oserror(self, tmp_path):
        """getsize raising is swallowed (lines 915-916)."""
        (tmp_path / "x.bsp").write_bytes(b"x")
        with patch("os.path.getsize", side_effect=OSError("gone")):
            size = spk_auto.get_cache_size(str(tmp_path))
        assert size == 0.0


# =============================================================================
# prune_old_cache (lines 998-1000, 1005->1004)
# =============================================================================


class TestPruneOldCache:
    """Cover prune removal, OSError arm and registry cleanup."""

    def setup_method(self):
        spk_auto.disable_all()

    def teardown_method(self):
        spk_auto.disable_all()

    def test_invalid_max_age(self):
        """Non-positive max_age raises ValueError (lines 963-964)."""
        with pytest.raises(ValueError, match="positive integer"):
            spk_auto.prune_old_cache(0)

    def test_prune_remove_oserror_and_registry(self, tmp_path):
        """Old file whose remove fails is swallowed; registry cleared (998-1006)."""
        f = tmp_path / "old.bsp"
        f.write_bytes(b"x")
        # Make the file very old so it passes the age threshold.
        os.utime(str(f), (0.0, 0.0))

        cfg = spk_auto.AutoSpkConfig(ipl=CHIRON, body_id="2060")
        cfg.spk_path = str(tmp_path / "vanished.bsp")  # never created
        spk_auto._AUTO_SPK_REGISTRY[CHIRON] = cfg

        with patch("os.remove", side_effect=OSError("locked")):
            deleted = spk_auto.prune_old_cache(1, str(tmp_path))
        assert deleted == 0  # remove failed
        assert cfg.spk_path is None  # registry cleared (1005->1006)

    def test_prune_removes_old_file_and_keeps_existing(self, tmp_path):
        """Old file removed; existing config spk_path kept (branch 1005->1004)."""
        f = tmp_path / "old.bsp"
        f.write_bytes(b"x")
        os.utime(str(f), (0.0, 0.0))

        # Config whose spk_path still exists -> False branch (1005->1004).
        keep = tmp_path / "keep.bsp"
        keep.write_bytes(b"x")
        cfg = spk_auto.AutoSpkConfig(ipl=CHIRON, body_id="2060")
        cfg.spk_path = str(keep)
        spk_auto._AUTO_SPK_REGISTRY[CHIRON] = cfg

        deleted = spk_auto.prune_old_cache(1, str(tmp_path))
        assert deleted == 1
        assert not f.exists()
        assert cfg.spk_path == str(keep)


# =============================================================================
# _jd_to_iso_date / _iso_to_jd (lines 1104, 1148-1160)
# =============================================================================


class TestJdConversions:
    """Cover the Julian/Gregorian branch and the inverse helper."""

    def test_jd_to_iso_julian_branch(self):
        """z < 2299161 hits the Julian-calendar branch (line 1104)."""
        # JD 2299160 is before the Gregorian cutover.
        result = spk_auto._jd_to_iso_date(2299160.0)
        assert isinstance(result, str)
        assert len(result.split("-")) == 3

    def test_iso_to_jd_basic(self):
        """_iso_to_jd converts a January date (lines 1148-1160)."""
        assert spk_auto._iso_to_jd("2000-01-01") == 2451544.5

    def test_iso_to_jd_jan_feb_branch(self):
        """month <= 2 triggers year/month adjustment (lines 1152-1154)."""
        assert spk_auto._iso_to_jd("2000-02-01") == 2451575.5

    def test_iso_to_jd_march(self):
        """month > 2 skips the adjustment branch."""
        jd = spk_auto._iso_to_jd("2000-03-01")
        assert jd > 2451575.5


# =============================================================================
# _find_covering_spk filename parsing (1222, 1230->1220, 1247-1248)
# =============================================================================


class TestFindCoveringSpk:
    """Cover non-bsp skip, short-name skip and JD-parse failure."""

    def test_skips_non_bsp(self, tmp_path):
        """Non-.bsp files are skipped (line 1222)."""
        (tmp_path / "2060_1_2.txt").write_bytes(b"x")
        result = spk_auto._find_covering_spk("2060", 1.0, 2.0, str(tmp_path))
        assert result is None

    def test_short_filename_skipped(self, tmp_path):
        """Filename with < 3 parts is skipped (branch 1230->1220)."""
        (tmp_path / "2060_x.bsp").write_bytes(b"x")
        result = spk_auto._find_covering_spk("2060", 1.0, 2.0, str(tmp_path))
        assert result is None

    def test_jd_parse_value_error(self, tmp_path):
        """Non-numeric JD parts continue past the file (lines 1247-1248)."""
        (tmp_path / "2060_aa_bb.bsp").write_bytes(b"x")
        result = spk_auto._find_covering_spk("2060", 1.0, 2.0, str(tmp_path))
        assert result is None


# =============================================================================
# auto_get_spk padding + register-cached branches (1409-1410, 1421, 1437-1447)
# =============================================================================


class TestAutoGetSpk:
    """Cover padding, registration-of-cached and the locked re-checks."""

    def setup_method(self):
        spk_auto.disable_all()
        from libephemeris import state

        state._SPK_BODY_MAP.clear()

    def teardown_method(self):
        spk_auto.disable_all()
        from libephemeris import state

        state._SPK_BODY_MAP.clear()

    def test_date_padding_applied(self, tmp_path):
        """Positive padding widens the requested range (lines 1409-1410)."""
        captured = {}

        def fake_download(**kwargs):
            captured.update(kwargs)
            with open(kwargs["output_path"], "wb") as f:
                f.write(b"x")

        with patch("libephemeris.state.get_spk_date_padding", return_value=10), \
                patch.object(
                    spk_auto, "_download_spk_astroquery", side_effect=fake_download
                ):
            spk_auto.auto_get_spk("2060", 2458849.5, 2462502.5, str(tmp_path))
        # With padding the start date moves earlier than 2020-01-01.
        assert captured["start"] < "2020-01-01"

    def test_register_existing_covering(self, tmp_path):
        """Existing covering file is registered when ipl given (line 1421)."""
        f = tmp_path / "2060_2450000_2470000.bsp"
        f.write_bytes(b"x")
        # _find_covering_spk reads real segment metadata via get_spk_coverage;
        # the dummy file has none, so patch it to report a covering span.
        with patch(
            "libephemeris.spk.get_spk_coverage",
            return_value=(2450000.0, 2470000.0),
        ), patch.object(spk_auto, "_is_valid_bsp", return_value=True), \
                patch.object(
                    spk_auto, "_register_spk_after_download"
                ) as mock_reg:
            result = spk_auto.auto_get_spk(
                "2060", 2458849.5, 2462502.5, str(tmp_path), ipl=CHIRON
            )
        assert result == str(f)
        mock_reg.assert_called_once()

    def test_locked_output_path_valid(self, tmp_path):
        """Output path exists+valid inside lock -> register cached (1435-1439)."""
        # Make _find_covering_spk miss but the exact output_path exist & valid.
        out = tmp_path / "2060_2458849_2462502.bsp"
        out.write_bytes(b"x")
        with patch.object(
            spk_auto, "_find_covering_spk", return_value=None
        ), patch.object(spk_auto, "_is_valid_bsp", return_value=True), \
                patch.object(
                    spk_auto, "_register_spk_after_download"
                ) as mock_reg, \
                patch.object(spk_auto, "_download_spk_astroquery") as mock_dl:
            result = spk_auto.auto_get_spk(
                "2060", 2458849.5, 2462502.5, str(tmp_path), ipl=CHIRON
            )
        assert result == str(out)
        mock_dl.assert_not_called()
        mock_reg.assert_called_once()

    def test_locked_output_path_valid_no_ipl(self, tmp_path):
        """Output exists+valid in lock, ipl None -> no register (branch 1437->1439)."""
        out = tmp_path / "2060_2458849_2462502.bsp"
        out.write_bytes(b"x")
        with patch.object(
            spk_auto, "_find_covering_spk", return_value=None
        ), patch.object(spk_auto, "_is_valid_bsp", return_value=True), \
                patch.object(
                    spk_auto, "_register_spk_after_download"
                ) as mock_reg, \
                patch.object(spk_auto, "_download_spk_astroquery") as mock_dl:
            result = spk_auto.auto_get_spk(
                "2060", 2458849.5, 2462502.5, str(tmp_path)
            )
        assert result == str(out)
        mock_dl.assert_not_called()
        mock_reg.assert_not_called()

    def test_locked_covering_appears_no_ipl(self, tmp_path):
        """Covering appears in lock, ipl None -> no register (branch 1445->1447)."""
        covering = tmp_path / "2060_2450000_2470000.bsp"
        calls = {"n": 0}

        def find_side_effect(body_id, jd_start, jd_end, cache_dir):
            calls["n"] += 1
            if calls["n"] == 1:
                return None
            covering.write_bytes(b"x")
            return str(covering)

        with patch.object(
            spk_auto, "_find_covering_spk", side_effect=find_side_effect
        ), patch.object(
            spk_auto, "_is_valid_bsp", return_value=False
        ), patch.object(
            spk_auto, "_register_spk_after_download"
        ) as mock_reg, patch.object(
            spk_auto, "_download_spk_astroquery"
        ) as mock_dl:
            result = spk_auto.auto_get_spk(
                "2060", 2458849.5, 2462502.5, str(tmp_path)
            )
        assert result == str(covering)
        mock_dl.assert_not_called()
        mock_reg.assert_not_called()

    def test_locked_covering_appears(self, tmp_path):
        """Covering file appears on the in-lock re-check (lines 1442-1447)."""
        out_name = "2060_2458849_2462502.bsp"
        covering = tmp_path / "2060_2450000_2470000.bsp"

        calls = {"n": 0}

        def find_side_effect(body_id, jd_start, jd_end, cache_dir):
            calls["n"] += 1
            if calls["n"] == 1:
                return None  # first (pre-lock) check misses
            covering.write_bytes(b"x")
            return str(covering)

        with patch.object(
            spk_auto, "_find_covering_spk", side_effect=find_side_effect
        ), patch.object(
            spk_auto, "_is_valid_bsp", return_value=False
        ), patch.object(
            spk_auto, "_register_spk_after_download"
        ) as mock_reg, patch.object(
            spk_auto, "_download_spk_astroquery"
        ) as mock_dl:
            result = spk_auto.auto_get_spk(
                "2060", 2458849.5, 2462502.5, str(tmp_path), ipl=CHIRON
            )
        assert result == str(covering)
        assert out_name  # name unused otherwise
        mock_dl.assert_not_called()
        mock_reg.assert_called_once()


# =============================================================================
# _register_spk_after_download detected NAIF (lines 1491, 1504)
# =============================================================================


class TestRegisterSpkAfterDownload:
    """Cover the detected-NAIF override and successful registration."""

    def setup_method(self):
        from libephemeris import state

        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def teardown_method(self):
        from libephemeris import state

        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def test_detected_naif_used_and_registered(self, tmp_path):
        """Detected NAIF overrides None and register is called (1491, 1504)."""
        f = tmp_path / "x.bsp"
        f.write_bytes(b"x")
        with patch.object(
            spk_auto, "_detect_naif_id_from_file", return_value=20002060
        ), patch("libephemeris.spk.register_spk_body") as mock_reg:
            spk_auto._register_spk_after_download(str(f), "2060", CHIRON, None)
        mock_reg.assert_called_once()
        assert mock_reg.call_args[0][2] == 20002060

    def test_deduced_naif_registers(self, tmp_path):
        """naif_id None, detected None, deduced succeeds (branch 1496->1504)."""
        f = tmp_path / "x.bsp"
        f.write_bytes(b"x")
        with patch.object(
            spk_auto, "_detect_naif_id_from_file", return_value=None
        ), patch("libephemeris.spk.register_spk_body") as mock_reg:
            spk_auto._register_spk_after_download(str(f), "2060", CHIRON, None)
        mock_reg.assert_called_once()
        assert mock_reg.call_args[0][2] == NAIF_CHIRON


# =============================================================================
# download_spk_from_horizons (1580, 1591, 1602, 1617->1620, 1622-1632)
# =============================================================================


class TestDownloadSpkFromHorizons:
    """Cover validation, dir creation, location mapping and error arms."""

    def setup_method(self):
        from libephemeris import state

        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def teardown_method(self):
        from libephemeris import state

        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def test_invalid_range_raises(self, tmp_path):
        """jd_end <= jd_start raises ValueError (line 1580)."""
        with pytest.raises(ValueError, match="must be greater than"):
            spk_auto.download_spk_from_horizons(
                "2060", 2462502.5, 2458849.5, str(tmp_path / "x.bsp")
            )

    def test_creates_output_dir_and_renames(self, tmp_path):
        """Missing output dir is created; result renamed (1591, 1617->1620)."""
        out = tmp_path / "sub" / "chiron.bsp"
        produced = tmp_path / "produced.bsp"

        def fake_download(**kwargs):
            produced.write_bytes(b"x")
            return str(produced)

        with patch("libephemeris.spk.download_spk", side_effect=fake_download):
            result = spk_auto.download_spk_from_horizons(
                "2060", 2458849.5, 2462502.5, str(out)
            )
        assert result == str(out)
        assert out.exists()

    def test_no_rename_when_same(self, tmp_path):
        """result_path == output_path skips rename (branch 1617->1620)."""
        out = tmp_path / "same.bsp"

        def fake_download(**kwargs):
            out.write_bytes(b"x")
            return str(out)

        with patch("libephemeris.spk.download_spk", side_effect=fake_download):
            result = spk_auto.download_spk_from_horizons(
                "2060", 2458849.5, 2462502.5, str(out)
            )
        assert result == str(out)

    def test_non_at_location_passthrough(self, tmp_path):
        """A non-'@' location is passed through unchanged (line 1602)."""
        out = tmp_path / "loc.bsp"
        captured = {}

        def fake_download(**kwargs):
            captured.update(kwargs)
            out.write_bytes(b"x")
            return str(out)

        with patch("libephemeris.spk.download_spk", side_effect=fake_download):
            spk_auto.download_spk_from_horizons(
                "2060", 2458849.5, 2462502.5, str(out), location="500@10"
            )
        assert captured["center"] == "500@10"

    def test_value_error_wrapped(self, tmp_path):
        """ValueError from download is re-wrapped (lines 1622-1626)."""
        out = tmp_path / "v.bsp"
        with patch(
            "libephemeris.spk.download_spk", side_effect=ValueError("bad body")
        ):
            with pytest.raises(ValueError, match="Failed to download SPK"):
                spk_auto.download_spk_from_horizons(
                    "2060", 2458849.5, 2462502.5, str(out)
                )

    def test_connection_error_wrapped(self, tmp_path):
        """ConnectionError is re-wrapped (lines 1627-1630)."""
        out = tmp_path / "c.bsp"
        with patch(
            "libephemeris.spk.download_spk",
            side_effect=ConnectionError("offline"),
        ):
            with pytest.raises(ConnectionError, match="Network error"):
                spk_auto.download_spk_from_horizons(
                    "2060", 2458849.5, 2462502.5, str(out)
                )

    def test_generic_error_wrapped(self, tmp_path):
        """OSError/KeyError fall into the generic arm (lines 1631-1632)."""
        out = tmp_path / "g.bsp"
        with patch(
            "libephemeris.spk.download_spk", side_effect=OSError("disk")
        ):
            with pytest.raises(ValueError, match="Failed to download SPK"):
                spk_auto.download_spk_from_horizons(
                    "2060", 2458849.5, 2462502.5, str(out)
                )


# =============================================================================
# discover_local_spks (lines 1674-1781)
# =============================================================================


class TestDiscoverLocalSpks:
    """Cover the local SPK discovery scan and registration logic."""

    def setup_method(self):
        from libephemeris import state

        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def teardown_method(self):
        from libephemeris import state

        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

    def test_path_not_dir(self, tmp_path):
        """Non-directory path returns empty (lines 1680-1682)."""
        result = spk_auto.discover_local_spks(str(tmp_path / "missing"))
        assert result == {}

    def test_no_bsp_files(self, tmp_path):
        """Directory with no .bsp files returns empty (lines 1686-1688)."""
        (tmp_path / "readme.txt").write_bytes(b"x")
        result = spk_auto.discover_local_spks(str(tmp_path))
        assert result == {}

    def test_register_new_body(self, tmp_path):
        """A matching NAIF target is registered (lines 1755-1763)."""
        (tmp_path / "chiron.bsp").write_bytes(b"x")
        with patch(
            "libephemeris.spk._get_spk_targets", return_value=[NAIF_CHIRON]
        ), patch("libephemeris.spk.register_spk_body") as mock_reg:
            result = spk_auto.discover_local_spks(str(tmp_path))
        assert mock_reg.called
        assert "registered" in result.values()

    def test_targets_read_error_skipped(self, tmp_path):
        """File whose targets cannot be read is skipped (lines 1715-1716)."""
        (tmp_path / "chiron.bsp").write_bytes(b"x")
        with patch(
            "libephemeris.spk._get_spk_targets",
            side_effect=ValueError("bad"),
        ):
            result = spk_auto.discover_local_spks(str(tmp_path))
        assert result == {}

    def test_empty_targets_skipped(self, tmp_path):
        """File with no targets is skipped (lines 1718-1719)."""
        (tmp_path / "chiron.bsp").write_bytes(b"x")
        with patch("libephemeris.spk._get_spk_targets", return_value=[]):
            result = spk_auto.discover_local_spks(str(tmp_path))
        assert result == {}

    def test_unknown_target_skipped(self, tmp_path):
        """A NAIF not in the reverse map is skipped (branch 1722)."""
        (tmp_path / "x.bsp").write_bytes(b"x")
        with patch(
            "libephemeris.spk._get_spk_targets", return_value=[999999999]
        ):
            result = spk_auto.discover_local_spks(str(tmp_path))
        assert result == {}

    def test_register_error_recorded(self, tmp_path):
        """register_spk_body failure is recorded as error (lines 1764-1771)."""
        (tmp_path / "chiron.bsp").write_bytes(b"x")
        with patch(
            "libephemeris.spk._get_spk_targets", return_value=[NAIF_CHIRON]
        ), patch(
            "libephemeris.spk.register_spk_body",
            side_effect=ValueError("reg fail"),
        ):
            result = spk_auto.discover_local_spks(str(tmp_path))
        assert any(str(v).startswith("error:") for v in result.values())

    def test_already_registered_wider_reregisters(self, tmp_path):
        """Wider new coverage triggers re-registration (lines 1728-1748)."""
        from libephemeris import state

        (tmp_path / "chiron.bsp").write_bytes(b"x")
        state._SPK_BODY_MAP[CHIRON] = ("existing.bsp", NAIF_CHIRON)

        def coverage_side(path):
            if path == "existing.bsp":
                return (2451545.0, 2451645.0)  # narrow span
            return (2451545.0, 2461545.0)  # wide span

        with patch(
            "libephemeris.spk._get_spk_targets", return_value=[NAIF_CHIRON]
        ), patch(
            "libephemeris.spk.get_spk_coverage", side_effect=coverage_side
        ), patch("libephemeris.spk.register_spk_body") as mock_reg:
            result = spk_auto.discover_local_spks(str(tmp_path))
        mock_reg.assert_called()
        assert result.get(_body_name()) == "registered"

    def test_already_registered_not_wider(self, tmp_path):
        """Equal/narrower coverage stays already_registered (lines 1751-1753)."""
        from libephemeris import state

        (tmp_path / "chiron.bsp").write_bytes(b"x")
        state._SPK_BODY_MAP[CHIRON] = ("existing.bsp", NAIF_CHIRON)

        with patch(
            "libephemeris.spk._get_spk_targets", return_value=[NAIF_CHIRON]
        ), patch(
            "libephemeris.spk.get_spk_coverage",
            return_value=(2451545.0, 2461545.0),
        ), patch("libephemeris.spk.register_spk_body") as mock_reg:
            result = spk_auto.discover_local_spks(str(tmp_path))
        mock_reg.assert_not_called()
        assert result.get(_body_name()) == "already_registered"

    def test_already_registered_coverage_none(self, tmp_path):
        """Falsy coverage skips the wider-span check (branch 1733->1751)."""
        from libephemeris import state

        (tmp_path / "chiron.bsp").write_bytes(b"x")
        state._SPK_BODY_MAP[CHIRON] = ("existing.bsp", NAIF_CHIRON)

        with patch(
            "libephemeris.spk._get_spk_targets", return_value=[NAIF_CHIRON]
        ), patch(
            "libephemeris.spk.get_spk_coverage", return_value=None
        ), patch("libephemeris.spk.register_spk_body") as mock_reg:
            result = spk_auto.discover_local_spks(str(tmp_path))
        mock_reg.assert_not_called()
        assert result.get(_body_name()) == "already_registered"

    def test_already_registered_body_already_recorded(self, tmp_path):
        """Second file for an already-recorded body skips re-record (1751->1753)."""
        from libephemeris import state

        # Two files both resolving to CHIRON via two NAIF targets each.
        (tmp_path / "a.bsp").write_bytes(b"x")
        (tmp_path / "b.bsp").write_bytes(b"x")
        state._SPK_BODY_MAP[CHIRON] = ("existing.bsp", NAIF_CHIRON)

        with patch(
            "libephemeris.spk._get_spk_targets", return_value=[NAIF_CHIRON]
        ), patch(
            "libephemeris.spk.get_spk_coverage",
            return_value=(2451545.0, 2461545.0),
        ), patch("libephemeris.spk.register_spk_body") as mock_reg:
            result = spk_auto.discover_local_spks(str(tmp_path))
        mock_reg.assert_not_called()
        assert result.get(_body_name()) == "already_registered"

    def test_reverse_map_low_naif_and_collision(self, tmp_path):
        """Crafted map hits asteroid_number<=0 and collision arcs (1698,1700)."""
        # Two entries: one with naif_id <= 2_000_000 (asteroid_number <= 0 ->
        # branch 1698->1692 skips the horizons mapping) and two entries whose
        # horizons_naif collide (branch 1700->1692 skips the duplicate).
        crafted = {
            CHIRON: ("2060", NAIF_CHIRON),
            777: ("low", 1500000),  # naif <= 2_000_000
            888: ("dupA", 2000060),  # horizons_naif = 20000060
            999: ("dupB", 2000060),  # same naif -> same horizons_naif (collision)
        }
        (tmp_path / "x.bsp").write_bytes(b"x")
        with patch(
            "libephemeris.constants.SPK_BODY_NAME_MAP", crafted
        ), patch(
            "libephemeris.spk._get_spk_targets", return_value=[999999999]
        ), patch("libephemeris.spk._get_body_name", return_value=None):
            result = spk_auto.discover_local_spks(str(tmp_path))
        assert result == {}

    def test_already_registered_coverage_error(self, tmp_path):
        """Coverage read failure falls through to already_registered (1749-1753)."""
        from libephemeris import state

        (tmp_path / "chiron.bsp").write_bytes(b"x")
        state._SPK_BODY_MAP[CHIRON] = ("existing.bsp", NAIF_CHIRON)

        with patch(
            "libephemeris.spk._get_spk_targets", return_value=[NAIF_CHIRON]
        ), patch(
            "libephemeris.spk.get_spk_coverage",
            side_effect=ValueError("bad"),
        ), patch("libephemeris.spk.register_spk_body") as mock_reg:
            result = spk_auto.discover_local_spks(str(tmp_path))
        mock_reg.assert_not_called()
        assert result.get(_body_name()) == "already_registered"


def _body_name() -> str:
    """Resolve the body name discover_local_spks records for CHIRON."""
    from libephemeris import spk
    from libephemeris.constants import SPK_BODY_NAME_MAP

    horizons_id = SPK_BODY_NAME_MAP[CHIRON][0]
    return spk._get_body_name(CHIRON) or horizons_id


# =============================================================================
# ensure_all_ephemerides (lines 1820-1961)
# =============================================================================


class TestEnsureAllEphemerides:
    """Cover the ensure-all orchestration across cached/download/error arms."""

    def setup_method(self):
        from libephemeris import state

        state._SPK_BODY_MAP.clear()

    def teardown_method(self):
        from libephemeris import state

        state._SPK_BODY_MAP.clear()

    def test_full_run_eph_cached_and_spks(self, tmp_path):
        """Ephemeris file found in search path -> cached, SPKs downloaded."""
        from libephemeris import state

        tier = state._get_current_tier()
        eph_file = tier.ephemeris_file

        # Put the ephemeris file under _EPHEMERIS_PATH so it is "cached".
        eph_dir = tmp_path / "eph"
        eph_dir.mkdir()
        (eph_dir / eph_file).write_bytes(b"x")

        with patch.object(state, "_EPHEMERIS_PATH", str(eph_dir)), \
                patch("libephemeris.spk.download_and_register_spk") as mock_dl, \
                patch("sys.stdout") as mock_stdout:
            mock_stdout.isatty.return_value = False
            result = spk_auto.ensure_all_ephemerides(show_progress=True)

        assert result["tier"] == tier.name
        assert result["summary"]["cached"] >= 1
        assert mock_dl.called
        assert result["summary"]["downloaded"] >= 1

    def test_eph_download_path_and_progress(self, tmp_path):
        """Ephemeris not found -> loader path; progress printed on tty."""
        from libephemeris import state

        loader = MagicMock()

        with patch.object(state, "_EPHEMERIS_PATH", None), \
                patch("os.path.exists", return_value=False), \
                patch("libephemeris.state.get_loader", return_value=loader), \
                patch("libephemeris.spk.download_and_register_spk"), \
                patch("sys.stdout") as mock_stdout:
            mock_stdout.isatty.return_value = True
            result = spk_auto.ensure_all_ephemerides(force_download=True)

        loader.assert_called_once()
        assert result["summary"]["downloaded"] >= 1
        assert mock_stdout.write.called

    def test_eph_download_error(self, tmp_path):
        """Loader failure records an ephemeris error (lines 1880-1882)."""
        from libephemeris import state

        def loader(_):
            raise OSError("download failed")

        with patch.object(state, "_EPHEMERIS_PATH", None), \
                patch("os.path.exists", return_value=False), \
                patch("libephemeris.state.get_loader", return_value=lambda f: loader(f)), \
                patch("libephemeris.spk.download_and_register_spk"), \
                patch("sys.stdout") as mock_stdout:
            mock_stdout.isatty.return_value = False
            result = spk_auto.ensure_all_ephemerides()

        assert "error:" in str(result["ephemeris"]["status"])
        assert result["summary"]["errors"] >= 1

    def test_spk_already_registered_cached(self, tmp_path):
        """Registered bodies count as cached and skip download (lines 1919-1922)."""
        from libephemeris import state
        from libephemeris.constants import SPK_BODY_NAME_MAP

        # Register every downloadable body so they're all "cached".
        for ipl in SPK_BODY_NAME_MAP:
            state._SPK_BODY_MAP[ipl] = ("f.bsp", 1)

        # Ephemeris found to avoid loader path.
        tier = state._get_current_tier()
        eph_dir = tmp_path / "eph"
        eph_dir.mkdir()
        (eph_dir / tier.ephemeris_file).write_bytes(b"x")

        with patch.object(state, "_EPHEMERIS_PATH", str(eph_dir)), \
                patch("libephemeris.spk.download_and_register_spk") as mock_dl, \
                patch("sys.stdout") as mock_stdout:
            mock_stdout.isatty.return_value = False
            result = spk_auto.ensure_all_ephemerides()

        mock_dl.assert_not_called()
        assert result["summary"]["cached"] >= 1

    def test_spk_download_failure_fallback(self, tmp_path):
        """SPK download failure records a fallback (lines 1935-1943)."""
        from libephemeris import state

        tier = state._get_current_tier()
        eph_dir = tmp_path / "eph"
        eph_dir.mkdir()
        (eph_dir / tier.ephemeris_file).write_bytes(b"x")

        with patch.object(state, "_EPHEMERIS_PATH", str(eph_dir)), \
                patch(
                    "libephemeris.spk.download_and_register_spk",
                    side_effect=ValueError("no spk"),
                ), patch("sys.stdout") as mock_stdout:
            mock_stdout.isatty.return_value = False
            result = spk_auto.ensure_all_ephemerides()

        assert result["summary"]["fallback"] >= 1

    def test_blocked_body_fallback(self, tmp_path):
        """A blocked body is reported as fallback (lines 1912-1916)."""
        from libephemeris import state
        from libephemeris.constants import BENNU

        tier = state._get_current_tier()
        eph_dir = tmp_path / "eph"
        eph_dir.mkdir()
        (eph_dir / tier.ephemeris_file).write_bytes(b"x")

        # Patch the map so only the blocked body is processed.
        with patch.object(state, "_EPHEMERIS_PATH", str(eph_dir)), \
                patch(
                    "libephemeris.constants.SPK_BODY_NAME_MAP",
                    {BENNU: ("101955", 2101955)},
                ), patch(
                    "libephemeris.spk.download_and_register_spk"
                ) as mock_dl, patch("sys.stdout") as mock_stdout:
            mock_stdout.isatty.return_value = False
            result = spk_auto.ensure_all_ephemerides()

        mock_dl.assert_not_called()
        assert result["summary"]["fallback"] >= 1
