# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for ``libephemeris.spk``.

This module exercises every reachable line/branch of the SPK kernel support
module using full mocks. It never touches the real network: all HTTP access
(``urllib.request.urlopen``), SPK file IO (``jplephem``/Skyfield kernels), and
optional-dependency imports (the vendored ``spktype21``) are mocked or replaced
with small in-memory doubles.

The real DE440 kernel (already downloaded in the dev environment) is reused as a
type 2/3 "fake" SPK kernel where a concrete Skyfield ``VectorFunction`` is
needed; no network access is required for that path either.
"""

from __future__ import annotations

import base64
import json
import sys
import types
import urllib.error
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from libephemeris import spk, state
from libephemeris.constants import (
    CHIRON,
    FLG_HELCTR,
    FLG_J2000,
    FLG_NOABERR,
    FLG_NONUT,
    FLG_SIDEREAL,
    FLG_SPEED,
    FLG_TRUEPOS,
    NAIF_ASTEROID_OFFSET,
)
from libephemeris.exceptions import EphemerisRangeError


# ---------------------------------------------------------------------------
# Test doubles
# ---------------------------------------------------------------------------


class _FakeResponse:
    """Minimal context-manager double for ``urlopen`` return values."""

    def __init__(self, payload: dict):
        self._raw = json.dumps(payload).encode("utf-8")

    def read(self):
        return self._raw

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _mock_urlopen(payload: dict):
    """Build a ``urlopen`` mock returning ``payload`` as JSON once called."""
    mock = MagicMock()
    mock.return_value = _FakeResponse(payload)
    return mock


def _valid_spk_payload(size: int = 1024) -> dict:
    """Canned Horizons response with valid base64-encoded SPK bytes."""
    return {
        "spk_file_id": "test_id",
        "spk": base64.b64encode(b"\x00" * size).decode("ascii"),
    }


class _FakeSegment:
    def __init__(self, target=None, data_type=None, start_jd=None, end_jd=None):
        if target is not None:
            self.target = target
        if data_type is not None:
            self.data_type = data_type
        if start_jd is not None:
            self.start_jd = start_jd
        if end_jd is not None:
            self.end_jd = end_jd


class _FakeType21Kernel:
    """Double for an spktype21 kernel object."""

    def __init__(self, fail=False):
        self._fail = fail

    def compute_type21(self, center, naif_id, jd):
        if self._fail:
            raise ValueError("boom")
        # Return a plausible heliocentric position (km) and velocity (km/s).
        return ([1.0e8, 2.0e8, 3.0e7], [1.0, 2.0, 0.5])


# ---------------------------------------------------------------------------
# Autouse cleanup
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _clean_spk_state():
    """Save/restore SPK global maps and the module-global spktype21 handle."""
    saved_map = dict(state._SPK_BODY_MAP)
    saved_kernels = dict(state._SPK_KERNELS)
    saved_t21 = dict(state._SPK_TYPE21_KERNELS)
    saved_module = spk._spktype21_module

    state._SPK_BODY_MAP.clear()
    state._SPK_KERNELS.clear()
    state._SPK_TYPE21_KERNELS.clear()

    yield

    state._SPK_BODY_MAP.clear()
    state._SPK_BODY_MAP.update(saved_map)
    state._SPK_KERNELS.clear()
    state._SPK_KERNELS.update(saved_kernels)
    state._SPK_TYPE21_KERNELS.clear()
    state._SPK_TYPE21_KERNELS.update(saved_t21)
    spk._spktype21_module = saved_module


# ---------------------------------------------------------------------------
# _get_spktype21 (lines 73->80, 78-79)
# ---------------------------------------------------------------------------


class TestGetSpktype21:
    def test_import_error_marks_unavailable(self):
        """Absent vendored module -> ImportError arm sets sentinel, returns None."""
        spk._spktype21_module = None
        with patch.dict(sys.modules, {"libephemeris.vendor.spktype21": None}):
            result = spk._get_spktype21()
        assert result is None
        # Sentinel is now False (cached unavailable).
        assert spk._spktype21_module is False

    def test_already_loaded_skips_import(self):
        """Non-None cached module short-circuits the import (73->80)."""
        sentinel = object()
        spk._spktype21_module = sentinel
        assert spk._get_spktype21() is sentinel

    def test_present_module_import(self):
        """Inject a fake module so the success arm is exercised."""
        spk._spktype21_module = None
        fake_cls = object()
        fake_mod = types.ModuleType("libephemeris.vendor.spktype21")
        fake_mod.SPKType21 = fake_cls
        with patch.dict(sys.modules, {"libephemeris.vendor.spktype21": fake_mod}):
            result = spk._get_spktype21()
        assert result is fake_cls


# ---------------------------------------------------------------------------
# _deduce_naif_id (150->153) and _find_naif_id_for_asteroid (205-209)
# ---------------------------------------------------------------------------


class TestNaifHelpers:
    def test_deduce_with_explicit_number(self):
        """Passing asteroid_number skips the extract branch (150->153)."""
        assert spk._deduce_naif_id("ignored", 2060) == 2060 + NAIF_ASTEROID_OFFSET

    def test_find_naif_horizons(self):
        targets = [20002060, 10, 0]
        with patch.object(spk, "_get_spk_targets", return_value=targets):
            assert spk._find_naif_id_for_asteroid("f.bsp", 2060) == 20002060

    def test_find_naif_legacy(self):
        targets = [2002060, 10, 0]
        with patch.object(spk, "_get_spk_targets", return_value=targets):
            assert spk._find_naif_id_for_asteroid("f.bsp", 2060) == 2002060

    def test_find_naif_neither(self):
        with patch.object(spk, "_get_spk_targets", return_value=[10, 0]):
            assert spk._find_naif_id_for_asteroid("f.bsp", 2060) is None


# ---------------------------------------------------------------------------
# _get_horizons_id_for_body None path (line 270)
# ---------------------------------------------------------------------------


class TestHorizonsIdForBody:
    def test_unknown_body_returns_none(self):
        assert spk._get_horizons_id_for_body(-99999) is None


# ---------------------------------------------------------------------------
# download_spk (lines 400-404, 438, 443-445, 455-456, 466-467, 472-475,
# 485-491, 530)
# ---------------------------------------------------------------------------


class TestDownloadSpk:
    def test_cached_corrupt_redownload(self, tmp_path):
        """Existing-but-corrupt cached file -> warn, remove, re-download."""
        existing = tmp_path / "2060_202001_202501.bsp"
        existing.write_bytes(b"\x00" * 16)

        # First _is_valid_bsp call (cached) -> False, second (downloaded) -> True.
        with (
            patch("urllib.request.urlopen", _mock_urlopen(_valid_spk_payload())),
            patch(
                "libephemeris.spk._is_valid_bsp",
                side_effect=[False, True],
            ),
        ):
            path = spk.download_spk(
                "2060", "2020-01-01", "2025-01-01", path=str(tmp_path)
            )
        assert path.endswith(".bsp")

    def test_cached_corrupt_remove_oserror(self, tmp_path):
        """os.remove raising OSError is swallowed (lines 403-404)."""
        existing = tmp_path / "2060_202001_202501.bsp"
        existing.write_bytes(b"\x00" * 16)

        with (
            patch("urllib.request.urlopen", _mock_urlopen(_valid_spk_payload())),
            patch(
                "libephemeris.spk._is_valid_bsp",
                side_effect=[False, True],
            ),
            patch("libephemeris.spk.os.remove", side_effect=OSError("denied")),
        ):
            path = spk.download_spk(
                "2060", "2020-01-01", "2025-01-01", path=str(tmp_path)
            )
        assert path.endswith(".bsp")

    def test_error_key_in_response(self, tmp_path):
        """{'error': ...} response -> ValueError (line 438)."""
        with patch("urllib.request.urlopen", _mock_urlopen({"error": "bad command"})):
            with pytest.raises(ValueError, match="Horizons API error"):
                spk.download_spk("2060", "2020-01-01", "2025-01-01", path=str(tmp_path))

    def test_result_key_in_response(self, tmp_path):
        """Missing spk keys but 'result' present -> ValueError (443-444)."""
        with patch(
            "urllib.request.urlopen",
            _mock_urlopen({"result": "No matches found.\n" * 50}),
        ):
            with pytest.raises(ValueError, match="Horizons error"):
                spk.download_spk("2060", "2020-01-01", "2025-01-01", path=str(tmp_path))

    def test_unexpected_format(self, tmp_path):
        """Missing spk keys and no 'result' -> ValueError (445)."""
        with patch("urllib.request.urlopen", _mock_urlopen({"unexpected": 1})):
            with pytest.raises(ValueError, match="Unexpected Horizons response"):
                spk.download_spk("2060", "2020-01-01", "2025-01-01", path=str(tmp_path))

    def test_base64_decode_failure(self, tmp_path):
        """Invalid base64 SPK payload -> ValueError (455-456)."""
        payload = {"spk_file_id": "x", "spk": "!!!not-base64!!!"}
        with patch("urllib.request.urlopen", _mock_urlopen(payload)):
            with patch(
                "base64.b64decode",
                side_effect=ValueError("bad pad"),
            ):
                with pytest.raises(ValueError, match="Failed to decode SPK"):
                    spk.download_spk(
                        "2060", "2020-01-01", "2025-01-01", path=str(tmp_path)
                    )

    def test_downloaded_data_fails_validation(self, tmp_path):
        """Downloaded bytes fail _is_valid_bsp -> unlink + ValueError (466-475)."""
        with (
            patch("urllib.request.urlopen", _mock_urlopen(_valid_spk_payload())),
            patch("libephemeris.spk._is_valid_bsp", return_value=False),
        ):
            with pytest.raises(ValueError, match="failed validation"):
                spk.download_spk("2060", "2020-01-01", "2025-01-01", path=str(tmp_path))
        # Temp file should have been cleaned up.
        leftovers = list(tmp_path.glob("*.download"))
        assert leftovers == []

    def test_os_replace_failure_unlinks_temp(self, tmp_path):
        """os.replace raising leaves temp_path -> except unlinks it (472-475)."""
        with (
            patch("urllib.request.urlopen", _mock_urlopen(_valid_spk_payload())),
            patch("libephemeris.spk._is_valid_bsp", return_value=True),
            patch(
                "libephemeris.spk.os.replace",
                side_effect=OSError("cannot replace"),
            ),
        ):
            with pytest.raises(OSError, match="cannot replace"):
                spk.download_spk("2060", "2020-01-01", "2025-01-01", path=str(tmp_path))
        # The temp file was unlinked by the except block.
        assert list(tmp_path.glob("*.download")) == []

    def test_http_400_bad_body(self, tmp_path):
        """HTTP 400 -> parse error JSON and raise ValueError (485-491)."""

        class _HTTP400(urllib.error.HTTPError):
            def __init__(self):
                super().__init__(
                    url="http://x", code=400, msg="Bad", hdrs=None, fp=None
                )

            def read(self):
                return json.dumps({"message": "no such body"}).encode("utf-8")

        with patch("urllib.request.urlopen", side_effect=_HTTP400()):
            with pytest.raises(ValueError, match="Invalid body"):
                spk.download_spk("2060", "2020-01-01", "2025-01-01", path=str(tmp_path))

    def test_http_400_unparseable_body(self, tmp_path):
        """HTTP 400 with unparseable JSON falls back to str(e) (488-489)."""

        class _HTTP400Bad(urllib.error.HTTPError):
            def __init__(self):
                super().__init__(
                    url="http://x", code=400, msg="Bad", hdrs=None, fp=None
                )

            def read(self):
                return b"not json"

        with patch("urllib.request.urlopen", side_effect=_HTTP400Bad()):
            with pytest.raises(ValueError, match="Invalid body"):
                spk.download_spk("2060", "2020-01-01", "2025-01-01", path=str(tmp_path))

    def test_http_500_retries_then_fail(self, tmp_path):
        """HTTP 500 retries (time.sleep patched) then raises ConnectionError."""
        err = urllib.error.HTTPError(
            url="http://x", code=500, msg="ISE", hdrs=None, fp=None
        )
        with (
            patch("urllib.request.urlopen", side_effect=err),
            patch("libephemeris.spk.time.sleep"),
        ):
            with pytest.raises(ConnectionError):
                spk.download_spk("2060", "2020-01-01", "2025-01-01", path=str(tmp_path))

    def test_url_error_retries_then_success(self, tmp_path):
        """URLError twice then a valid response succeeds (retry loop)."""
        good = _FakeResponse(_valid_spk_payload())
        side = [
            urllib.error.URLError("conn 1"),
            urllib.error.URLError("conn 2"),
            good,
        ]
        with (
            patch("urllib.request.urlopen", side_effect=side),
            patch("libephemeris.spk._is_valid_bsp", return_value=True),
            patch("libephemeris.spk.time.sleep"),
        ):
            path = spk.download_spk(
                "2060", "2020-01-01", "2025-01-01", path=str(tmp_path)
            )
        assert path.endswith(".bsp")

    def test_url_error_all_fail(self, tmp_path):
        """URLError on every attempt -> ConnectionError (network branch)."""
        with (
            patch(
                "urllib.request.urlopen",
                side_effect=urllib.error.URLError("down"),
            ),
            patch("libephemeris.spk.time.sleep"),
        ):
            with pytest.raises(ConnectionError, match="Network error"):
                spk.download_spk("2060", "2020-01-01", "2025-01-01", path=str(tmp_path))

    def test_default_path_branch(self, tmp_path):
        """path=None uses _get_data_dir()/spk (line 379-382)."""
        with (
            patch(
                "libephemeris.state._get_data_dir",
                return_value=str(tmp_path),
            ),
            patch("urllib.request.urlopen", _mock_urlopen(_valid_spk_payload())),
            patch("libephemeris.spk._is_valid_bsp", return_value=True),
        ):
            path = spk.download_spk("Chiron", "2020-01-01", "2025-01-01")
        assert str(tmp_path) in path


# ---------------------------------------------------------------------------
# register_spk_body (lines 568, 572-586, 609-634, 655)
# ---------------------------------------------------------------------------


class TestRegisterSpkBody:
    def test_string_naif_id_and_relative_resolution(self, tmp_path):
        """str naif_id converted (568); relative path resolved in spk subdir."""
        spk_dir = tmp_path / "spk"
        spk_dir.mkdir()
        rel_name = "chiron.bsp"
        full = spk_dir / rel_name
        full.write_bytes(b"\x00" * 16)

        with (
            patch(
                "libephemeris.state._get_data_dir",
                return_value=str(tmp_path),
            ),
            patch("libephemeris.spk._detect_spk_type", return_value=2),
            patch("libephemeris.state._load_spk_kernel"),
        ):
            # Kernel present, target lookup OK.
            fake_kernel = MagicMock()
            fake_kernel.__getitem__ = MagicMock(return_value=object())
            state._SPK_KERNELS[str(full)] = fake_kernel
            spk.register_spk_body(CHIRON, rel_name, "2002060")

        assert CHIRON in state._SPK_BODY_MAP
        assert state._SPK_BODY_MAP[CHIRON][1] == 2002060

    def test_relative_found_via_cwd_fallthrough(self, tmp_path):
        """Relative name absent in data dirs but present via plain exists check.

        Exercises the for-else with no break, where ``os.path.exists(spk_file)``
        on the bare relative name returns True (branch 582->592).
        """
        rel_name = "cwd_only.bsp"

        def _fake_exists(p):
            # Both data-dir candidates are absent; the bare relative name exists.
            if p == rel_name:
                return True
            return False

        with (
            patch(
                "libephemeris.state._get_data_dir",
                return_value=str(tmp_path),
            ),
            patch("libephemeris.spk.os.path.exists", side_effect=_fake_exists),
            patch("libephemeris.spk._detect_spk_type", return_value=2),
            patch("libephemeris.state._load_spk_kernel"),
        ):
            # No cached kernel -> validation skipped, registration proceeds.
            spk.register_spk_body(CHIRON, rel_name, 2002060)
        assert CHIRON in state._SPK_BODY_MAP

    def test_type2_naif_not_found_no_names(self, tmp_path):
        """Kernel without names() -> skip names list, still raise (632->634)."""
        f = tmp_path / "t2.bsp"
        f.write_bytes(b"\x00" * 16)

        class _KNoNames:
            def __getitem__(self, key):
                raise KeyError(key)

        with (
            patch("libephemeris.spk._detect_spk_type", return_value=2),
            patch("libephemeris.state._load_spk_kernel"),
        ):
            state._SPK_KERNELS[str(f)] = _KNoNames()
            with pytest.raises(ValueError, match="not found in SPK kernel"):
                spk.register_spk_body(CHIRON, str(f), 2002060)

    def test_relative_not_found_raises(self, tmp_path):
        """Relative file not found anywhere -> SPKNotFoundError (582-586)."""
        from libephemeris.exceptions import SPKNotFoundError

        with patch("libephemeris.state._get_data_dir", return_value=str(tmp_path)):
            with pytest.raises(SPKNotFoundError):
                spk.register_spk_body(CHIRON, "missing.bsp", 2002060)

    def test_type21_kernel_none_raises(self, tmp_path):
        """Type 21 with kernel load failing -> ValueError (609-612)."""
        f = tmp_path / "t21.bsp"
        f.write_bytes(b"\x00" * 16)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch("libephemeris.spk._load_type21_kernel", return_value=None),
        ):
            with pytest.raises(ValueError, match="Failed to load type 21"):
                spk.register_spk_body(CHIRON, str(f), 2002060)

    def test_type21_kernel_ok(self, tmp_path):
        """Type 21 with kernel loaded -> registers (skips validation)."""
        f = tmp_path / "t21.bsp"
        f.write_bytes(b"\x00" * 16)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch(
                "libephemeris.spk._load_type21_kernel",
                return_value=_FakeType21Kernel(),
            ),
        ):
            spk.register_spk_body(CHIRON, str(f), 2002060)
        assert CHIRON in state._SPK_BODY_MAP

    def test_type2_naif_found_via_str(self, tmp_path):
        """Type 2/3: int key KeyError, str key OK (628-629 success)."""
        f = tmp_path / "t2.bsp"
        f.write_bytes(b"\x00" * 16)

        class _K:
            def __getitem__(self, key):
                if isinstance(key, int):
                    raise KeyError(key)
                return object()  # string lookup succeeds

        with (
            patch("libephemeris.spk._detect_spk_type", return_value=2),
            patch("libephemeris.state._load_spk_kernel"),
        ):
            state._SPK_KERNELS[str(f)] = _K()
            spk.register_spk_body(CHIRON, str(f), 2002060)
        assert CHIRON in state._SPK_BODY_MAP

    def test_type2_naif_not_found_raises(self, tmp_path):
        """Type 2/3: both int and str lookups KeyError -> ValueError (631-637)."""
        f = tmp_path / "t2.bsp"
        f.write_bytes(b"\x00" * 16)

        class _K:
            def __getitem__(self, key):
                raise KeyError(key)

            def names(self):
                return ["A", "B", "C"]

        with (
            patch("libephemeris.spk._detect_spk_type", return_value=2),
            patch("libephemeris.state._load_spk_kernel"),
        ):
            state._SPK_KERNELS[str(f)] = _K()
            with pytest.raises(ValueError, match="not found in SPK kernel"):
                spk.register_spk_body(CHIRON, str(f), 2002060)

    def test_type2_kernel_none_skips_validation(self, tmp_path):
        """Type 2/3 with kernel None -> registration still proceeds (621 False)."""
        f = tmp_path / "t2.bsp"
        f.write_bytes(b"\x00" * 16)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=2),
            patch("libephemeris.state._load_spk_kernel"),
        ):
            # No kernel in _SPK_KERNELS -> .get() returns None.
            spk.register_spk_body(CHIRON, str(f), 2002060)
        assert CHIRON in state._SPK_BODY_MAP

    def test_unregister_existing(self):
        """unregister removes a registered body (line 655)."""
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 2002060)
        spk.unregister_spk_body(CHIRON)
        assert CHIRON not in state._SPK_BODY_MAP


# ---------------------------------------------------------------------------
# get_spk_coverage (lines 701-708, 711, 726->725, 737-739, 748-764)
# ---------------------------------------------------------------------------


class TestGetSpkCoverage:
    def test_relative_resolution_and_not_found(self, tmp_path):
        """Relative path not found anywhere -> None (701-708, 711)."""
        with patch("libephemeris.state._get_data_dir", return_value=str(tmp_path)):
            assert spk.get_spk_coverage("missing.bsp") is None

    def test_relative_resolved_in_data_dir(self, tmp_path):
        """Relative file found in data dir resolves then computes (706-708)."""
        f = tmp_path / "found.bsp"
        f.write_bytes(b"\x00" * 16)
        with (
            patch(
                "libephemeris.state._get_data_dir",
                return_value=str(tmp_path),
            ),
            patch("libephemeris.spk._detect_spk_type", return_value=2),
        ):
            # Force the type 2/3 kernel-None branch (returns None).
            assert spk.get_spk_coverage("found.bsp") is None

    def test_type21_coverage(self, tmp_path):
        """Type 21: iterate jplephem segments and compute (726->725, 735-736)."""
        f = tmp_path / "t21.bsp"
        f.write_bytes(b"\x00" * 16)

        fake_spk = MagicMock()
        fake_spk.segments = [
            _FakeSegment(start_jd=2451545.0, end_jd=2470000.0),
            _FakeSegment(start_jd=2440000.0, end_jd=2460000.0),
            _FakeSegment(),  # no start_jd/end_jd -> 726->725 skip branch
        ]
        fake_mod = types.ModuleType("jplephem.spk")
        fake_mod.SPK = MagicMock()
        fake_mod.SPK.open = MagicMock(return_value=fake_spk)

        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch.dict(sys.modules, {"jplephem.spk": fake_mod}),
        ):
            cov = spk.get_spk_coverage(str(f))
        assert cov == (2440000.0, 2470000.0)

    def test_type21_no_valid_segments(self, tmp_path):
        """Type 21: segments lack start/end -> start_jd None -> None (735->739)."""
        f = tmp_path / "t21.bsp"
        f.write_bytes(b"\x00" * 16)

        fake_spk = MagicMock()
        # Segments with no start_jd/end_jd attributes at all.
        fake_spk.segments = [_FakeSegment(), _FakeSegment()]
        fake_mod = types.ModuleType("jplephem.spk")
        fake_mod.SPK = MagicMock()
        fake_mod.SPK.open = MagicMock(return_value=fake_spk)

        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch.dict(sys.modules, {"jplephem.spk": fake_mod}),
        ):
            assert spk.get_spk_coverage(str(f)) is None

    def test_type21_coverage_exception(self, tmp_path):
        """Type 21: jplephem raising -> None (737-739)."""
        f = tmp_path / "t21.bsp"
        f.write_bytes(b"\x00" * 16)

        fake_mod = types.ModuleType("jplephem.spk")
        fake_mod.SPK = MagicMock()
        fake_mod.SPK.open = MagicMock(side_effect=OSError("nope"))

        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch.dict(sys.modules, {"jplephem.spk": fake_mod}),
        ):
            assert spk.get_spk_coverage(str(f)) is None

    def test_type2_load_failure_returns_none(self, tmp_path):
        """Type 2/3: kernel not cached, load raises ValueError -> None (744-746)."""
        f = tmp_path / "t2.bsp"
        f.write_bytes(b"\x00" * 16)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=2),
            patch(
                "libephemeris.state._load_spk_kernel",
                side_effect=ValueError("bad"),
            ),
        ):
            assert spk.get_spk_coverage(str(f)) is None

    def test_type2_cached_kernel_none(self, tmp_path):
        """Type 2/3: cached kernel is None -> return None (line 750)."""
        f = tmp_path / "t2.bsp"
        f.write_bytes(b"\x00" * 16)
        state._SPK_KERNELS[str(f)] = None  # cached but None
        with patch("libephemeris.spk._detect_spk_type", return_value=2):
            assert spk.get_spk_coverage(str(f)) is None

    def test_type2_kernel_no_spk_attr(self, tmp_path):
        """Type 2/3: kernel lacks .spk -> fall through to None (753->764)."""
        f = tmp_path / "t2.bsp"
        f.write_bytes(b"\x00" * 16)

        class _NoSpk:
            pass

        state._SPK_KERNELS[str(f)] = _NoSpk()
        with patch("libephemeris.spk._detect_spk_type", return_value=2):
            assert spk.get_spk_coverage(str(f)) is None

    def test_type2_empty_segments(self, tmp_path):
        """Type 2/3: empty segments list -> None (755->764)."""
        f = tmp_path / "t2.bsp"
        f.write_bytes(b"\x00" * 16)
        kernel = MagicMock()
        kernel.spk.segments = []
        state._SPK_KERNELS[str(f)] = kernel
        with patch("libephemeris.spk._detect_spk_type", return_value=2):
            assert spk.get_spk_coverage(str(f)) is None

    def test_type2_coverage_from_segments(self, tmp_path):
        """Type 2/3: compute coverage from kernel.spk.segments (748-760)."""
        f = tmp_path / "t2.bsp"
        f.write_bytes(b"\x00" * 16)

        kernel = MagicMock()
        kernel.spk.segments = [
            _FakeSegment(start_jd=2451545.0, end_jd=2470000.0),
            _FakeSegment(start_jd=2440000.0, end_jd=2460000.0),
        ]
        state._SPK_KERNELS[str(f)] = kernel
        with patch("libephemeris.spk._detect_spk_type", return_value=2):
            cov = spk.get_spk_coverage(str(f))
        assert cov == (2440000.0, 2470000.0)

    def test_type2_segments_exception(self, tmp_path):
        """Type 2/3: segment float conversion raising -> None (761-764)."""
        f = tmp_path / "t2.bsp"
        f.write_bytes(b"\x00" * 16)

        bad_seg = MagicMock()
        type(bad_seg).start_jd = property(
            lambda self: (_ for _ in ()).throw(ValueError("x"))
        )
        kernel = MagicMock()
        kernel.spk.segments = [bad_seg]
        state._SPK_KERNELS[str(f)] = kernel
        with patch("libephemeris.spk._detect_spk_type", return_value=2):
            assert spk.get_spk_coverage(str(f)) is None


# ---------------------------------------------------------------------------
# _detect_spk_type (lines 792->791, 794)
# ---------------------------------------------------------------------------


class TestDetectSpkType:
    def test_detect_type2(self, tmp_path):
        """Segments without type-21 -> returns 2 (792->791 loop, 794)."""
        f = tmp_path / "t.bsp"
        f.write_bytes(b"\x00" * 16)

        fake_spk = MagicMock()
        fake_spk.segments = [
            _FakeSegment(data_type=2),
            _FakeSegment(data_type=3),
        ]
        fake_mod = types.ModuleType("jplephem.spk")
        fake_mod.SPK = MagicMock()
        fake_mod.SPK.open = MagicMock(return_value=fake_spk)
        with patch.dict(sys.modules, {"jplephem.spk": fake_mod}):
            assert spk._detect_spk_type(str(f)) == 2

    def test_detect_type21(self, tmp_path):
        """A type-21 segment -> returns 21."""
        f = tmp_path / "t.bsp"
        f.write_bytes(b"\x00" * 16)

        fake_spk = MagicMock()
        fake_spk.segments = [_FakeSegment(data_type=21)]
        fake_mod = types.ModuleType("jplephem.spk")
        fake_mod.SPK = MagicMock()
        fake_mod.SPK.open = MagicMock(return_value=fake_spk)
        with patch.dict(sys.modules, {"jplephem.spk": fake_mod}):
            assert spk._detect_spk_type(str(f)) == 21


# ---------------------------------------------------------------------------
# _load_type21_kernel (lines 819-822, 828-830)
# ---------------------------------------------------------------------------


class TestLoadType21Kernel:
    def test_cache_hit(self):
        cached = object()
        state._SPK_TYPE21_KERNELS["cached.bsp"] = cached
        assert spk._load_type21_kernel("cached.bsp") is cached

    def test_spktype21_unavailable(self):
        """SPKType21 None -> warn, return None (819-822)."""
        with patch("libephemeris.spk._get_spktype21", return_value=None):
            assert spk._load_type21_kernel("x.bsp") is None

    def test_open_success(self, tmp_path):
        f = tmp_path / "t21.bsp"
        f.write_bytes(b"\x00" * 16)
        fake_kernel = _FakeType21Kernel()
        fake_cls = MagicMock()
        fake_cls.open = MagicMock(return_value=fake_kernel)
        with patch("libephemeris.spk._get_spktype21", return_value=fake_cls):
            result = spk._load_type21_kernel(str(f))
        assert result is fake_kernel
        assert state._SPK_TYPE21_KERNELS[str(f)] is fake_kernel

    def test_open_failure(self, tmp_path):
        """SPKType21.open raising -> warn, return None (828-830)."""
        f = tmp_path / "t21.bsp"
        f.write_bytes(b"\x00" * 16)
        fake_cls = MagicMock()
        fake_cls.open = MagicMock(side_effect=ValueError("corrupt"))
        with patch("libephemeris.spk._get_spktype21", return_value=fake_cls):
            assert spk._load_type21_kernel(str(f)) is None


# ---------------------------------------------------------------------------
# _SpkType21Target (lines 928-938, 944-947)
# ---------------------------------------------------------------------------


class TestSpkType21Target:
    def test_at_array_and_scalar(self):
        """Exercise both scalar and array time paths plus at() (928-947)."""
        planets = state.get_planets()
        sun = planets["sun"]
        kernel = _FakeType21Kernel()
        target = spk._SpkType21Target(kernel, 20002060, sun)

        ts = state.get_timescale()

        # Scalar time path.
        t_scalar = ts.tt_jd(2451545.0)
        icrf = target.at(t_scalar)
        assert icrf is not None
        assert repr(target).startswith("<SpkType21Target")

        # Array time path (928-938).
        t_array = ts.tt_jd([2451545.0, 2451546.0, 2451547.0])
        pos, vel, a, b = target._at(t_array)
        assert pos.shape == (3, 3)
        assert vel.shape == (3, 3)
        assert a is None and b is None


# ---------------------------------------------------------------------------
# get_spk_type21_target (lines 998, 1000->1008, 1002->1008, 1006, 1010)
# ---------------------------------------------------------------------------


class TestGetSpkType21Target:
    def test_not_registered(self):
        assert spk.get_spk_type21_target(CHIRON) is None

    def test_not_type21(self, tmp_path):
        """Registered but not type 21 -> None (line 998)."""
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 20002060)
        with patch("libephemeris.spk._detect_spk_type", return_value=2):
            assert spk.get_spk_type21_target(CHIRON) is None

    def test_jd_outside_coverage(self):
        """jd_tt outside coverage -> None (1002->1008 -> 1006)."""
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 20002060)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch(
                "libephemeris.spk.get_spk_coverage",
                return_value=(2451545.0, 2460000.0),
            ),
        ):
            assert spk.get_spk_type21_target(CHIRON, jd_tt=2400000.0) is None

    def test_jd_within_coverage_kernel_none(self):
        """jd within coverage but kernel load fails -> None (1010)."""
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 20002060)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch(
                "libephemeris.spk.get_spk_coverage",
                return_value=(2451545.0, 2460000.0),
            ),
            patch("libephemeris.spk._load_type21_kernel", return_value=None),
        ):
            assert spk.get_spk_type21_target(CHIRON, jd_tt=2455000.0) is None

    def test_jd_none_success(self):
        """jd_tt None skips coverage; returns a target (1000->1008)."""
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 20002060)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch(
                "libephemeris.spk._load_type21_kernel",
                return_value=_FakeType21Kernel(),
            ),
        ):
            target = spk.get_spk_type21_target(CHIRON)
        assert isinstance(target, spk._SpkType21Target)

    def test_coverage_none_branch(self):
        """jd_tt provided but coverage None -> proceeds (1002->1008)."""
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 20002060)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch("libephemeris.spk.get_spk_coverage", return_value=None),
            patch(
                "libephemeris.spk._load_type21_kernel",
                return_value=_FakeType21Kernel(),
            ),
        ):
            target = spk.get_spk_type21_target(CHIRON, jd_tt=2455000.0)
        assert isinstance(target, spk._SpkType21Target)


# ---------------------------------------------------------------------------
# _calc_type21_position (lines 1097, 1109-1111, 1126-1130, 1135-1137,
# 1151-1152, 1163->1177, 1190, 1195->1201, 1220-1221, 1225, 1233-1234)
# ---------------------------------------------------------------------------


class TestCalcType21Position:
    def _run(self, iflag, kernel=None):
        kernel = kernel or _FakeType21Kernel()
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        return spk._calc_type21_position(kernel, 20002060, t, iflag)

    def test_geocentric_full(self):
        """Geocentric path: light-time, aberration, precession, nutation,
        speed (1097, 1109-1111, 1163->1177, 1190, 1220-1221, 1225)."""
        result = self._run(FLG_SPEED)
        assert result is not None
        lon, lat, dist = result[0], result[1], result[2]
        assert 0 <= lon < 360
        assert dist > 0

    def test_heliocentric(self):
        """Heliocentric path: light-time on Sun->body, no aberration
        (1124, 1158-1160)."""
        result = self._run(FLG_HELCTR | FLG_SPEED)
        assert result is not None

    def test_heliocentric_no_speed(self):
        """Heliocentric, no aberration, no speed -> earth_vel else (1101)."""
        result = self._run(FLG_HELCTR)
        assert result is not None

    def test_speed_at_ecliptic_pole(self):
        """Body at the ecliptic pole -> xy_sq == 0 speed branch (1228-1229)."""

        import math as _math

        # ICRS km coords that rotate to ecliptic (0, 0, Z): purely along the
        # ecliptic +Z axis so x == y == 0 exactly after the rotation, forcing
        # xy_sq == 0.
        _au_km = 149597870.7
        _e = _math.radians(23.4392911)
        _z = 3.0
        _y = -_z * _math.tan(_e)
        _pole_pos = [0.0, _y * _au_km, _z * _au_km]

        class _PoleKernel:
            def compute_type21(self, c, naif, jd):
                return (list(_pole_pos), [0.0, 0.0, 0.0])

        # Heliocentric so pos_final == pos_ecl (stays at the pole).
        result = self._run(FLG_HELCTR | FLG_SPEED, kernel=_PoleKernel())
        assert result is not None
        # x == y == 0 -> longitude speed is zero.
        assert result[3] == 0.0
        assert result[4] == 0.0

    def test_truepos_no_speed(self):
        """FLG_TRUEPOS skips light-time iteration (1135-1137 single compute)
        and no FLG_SPEED -> 1233-1234."""
        result = self._run(FLG_TRUEPOS)
        assert result is not None
        # No speed requested -> zeros.
        assert result[3] == 0.0 and result[4] == 0.0 and result[5] == 0.0

    def test_j2000_no_precession_no_nutation(self):
        """FLG_J2000 disables precession+nutation (1195->1201)."""
        result = self._run(FLG_J2000)
        assert result is not None

    def test_nonut(self):
        """FLG_NONUT disables nutation but keeps precession."""
        result = self._run(FLG_NONUT)
        assert result is not None

    def test_sidereal(self):
        """FLG_SIDEREAL applies ayanamsa correction."""
        result = self._run(FLG_SIDEREAL | FLG_SPEED)
        assert result is not None

    def test_light_time_compute_fails(self):
        """compute_type21 raising during light-time loop -> None (1113-1116)."""
        kernel = _FakeType21Kernel(fail=True)
        assert self._run(FLG_SPEED, kernel=kernel) is None

    def test_no_lighttime_compute_fails(self):
        """compute_type21 raising on the no-light-time branch -> None
        (1136-1138)."""
        kernel = _FakeType21Kernel(fail=True)
        assert self._run(FLG_TRUEPOS, kernel=kernel) is None

    def test_final_compute_fails(self):
        """Final compute raising -> None (1143-1145).

        First three light-time computes succeed, the final one fails.
        """
        calls = {"n": 0}

        class _K:
            def compute_type21(self, c, naif, jd):
                calls["n"] += 1
                # 3 iterations in the light-time loop, then the final compute.
                if calls["n"] >= 4:
                    raise ValueError("final fail")
                return ([1.0e8, 2.0e8, 3.0e7], [1.0, 2.0, 0.5])

        assert self._run(FLG_SPEED, kernel=_K()) is None


# ---------------------------------------------------------------------------
# download_and_register_spk (lines 1299->1328, 1302->1306, 1308, 1315-1318,
# 1321)
# ---------------------------------------------------------------------------


class TestDownloadAndRegister:
    def test_naif_provided_skips_deduction(self, tmp_path):
        """naif_id given -> skip whole deduction block (1299->1328)."""
        with (
            patch(
                "libephemeris.spk.download_spk",
                return_value=str(tmp_path / "f.bsp"),
            ),
            patch("libephemeris.spk.register_spk_body") as reg,
        ):
            out = spk.download_and_register_spk(
                "2060", CHIRON, "2020-01-01", "2025-01-01", naif_id=2002060
            )
        assert out == str(tmp_path / "f.bsp")
        reg.assert_called_once()

    def test_naif_from_asteroid_number(self, tmp_path):
        """Asteroid number extracted, found in file (1302->? success, 1308)."""
        with (
            patch(
                "libephemeris.spk.download_spk",
                return_value=str(tmp_path / "f.bsp"),
            ),
            patch(
                "libephemeris.spk._find_naif_id_for_asteroid",
                return_value=20002060,
            ),
            patch("libephemeris.spk.register_spk_body") as reg,
        ):
            spk.download_and_register_spk("2060", CHIRON, "2020-01-01", "2025-01-01")
        # register called with discovered naif_id
        assert reg.call_args[0][2] == 20002060

    def test_naif_legacy_fallback(self, tmp_path):
        """Number present but not found in file -> legacy deduce (1314-1316)."""
        with (
            patch(
                "libephemeris.spk.download_spk",
                return_value=str(tmp_path / "f.bsp"),
            ),
            patch(
                "libephemeris.spk._find_naif_id_for_asteroid",
                return_value=None,
            ),
            patch("libephemeris.spk.register_spk_body") as reg,
        ):
            spk.download_and_register_spk("2060", CHIRON, "2020-01-01", "2025-01-01")
        assert reg.call_args[0][2] == 2060 + NAIF_ASTEROID_OFFSET

    def test_naif_scan_unique_target(self, tmp_path):
        """Name body, single unique target in file scan (1315-1318, 1321)."""
        with (
            patch(
                "libephemeris.spk.download_spk",
                return_value=str(tmp_path / "f.bsp"),
            ),
            patch(
                "libephemeris.spk._get_spk_targets",
                return_value=[0, 10, 20000001],
            ),
            patch("libephemeris.spk.register_spk_body") as reg,
        ):
            spk.download_and_register_spk("Ceres;", CHIRON, "2020-01-01", "2025-01-01")
        assert reg.call_args[0][2] == 20000001

    def test_naif_cannot_deduce(self, tmp_path):
        """No number, scan ambiguous -> ValueError (1328-1333)."""
        with (
            patch(
                "libephemeris.spk.download_spk",
                return_value=str(tmp_path / "f.bsp"),
            ),
            patch(
                "libephemeris.spk._get_spk_targets",
                return_value=[20000001, 20000002],
            ),
        ):
            with pytest.raises(ValueError, match="Cannot deduce NAIF ID"):
                spk.download_and_register_spk(
                    "SomeName", CHIRON, "2020-01-01", "2025-01-01"
                )


# ---------------------------------------------------------------------------
# calc_spk_body_position (lines 1382->1391, 1385-1387, 1397-1506)
# ---------------------------------------------------------------------------


class _FakeTargetWrapper:
    """Wraps a real Skyfield VectorFunction so the kernel[naif_id] works."""

    def __init__(self, vf):
        self._vf = vf

    def at(self, t):
        return self._vf.at(t)


class _FakeKernel:
    def __init__(self, vf):
        self._vf = vf

    def __getitem__(self, naif_id):
        return self._vf


class _FakeKernelKeyError:
    def __getitem__(self, naif_id):
        raise KeyError(naif_id)


class _WrapTarget:
    """Type 2/3 target whose ecliptic longitude crosses the 0/360 seam.

    The very-large distance makes the observer (Earth) offset negligible, so the
    relative vector's longitude equals the configured value to high precision.
    Sampling at t-dt / t+dt yields a |delta| near 360 degrees, which forces the
    longitude-wrap corrections in calc_spk_body_position's speed block.
    """

    def __init__(self, t0_tt, lon_before, lon_after, dist=1.0e6):
        self._t0 = t0_tt
        self._before = lon_before
        self._after = lon_after
        self._dist = dist

    def at(self, t):
        import math

        from skyfield.framelib import ecliptic_frame
        from skyfield.positionlib import ICRF

        lon_deg = self._before if t.tt < self._t0 else self._after
        lon = math.radians(lon_deg)
        ecl = np.array([self._dist * math.cos(lon), self._dist * math.sin(lon), 0.0])
        rot = ecliptic_frame.rotation_at(t)
        icrs = rot.T.dot(ecl)
        return ICRF(icrs, velocity_au_per_d=np.zeros(3), t=t, center=0)


class TestCalcSpkBodyPosition:
    def test_type21_out_of_coverage_raises(self):
        """Type 21, jd outside coverage -> EphemerisRangeError (1385-1398)."""
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 20002060)
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch(
                "libephemeris.spk.get_spk_coverage",
                return_value=(2460000.0, 2470000.0),
            ),
        ):
            with pytest.raises(EphemerisRangeError):
                spk.calc_spk_body_position(t, CHIRON, 0)

    def test_type21_success(self):
        """Type 21 within coverage, calc returns a result (1399-1403)."""
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 20002060)
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch(
                "libephemeris.spk.get_spk_coverage",
                return_value=(2440000.0, 2470000.0),
            ),
            patch(
                "libephemeris.spk._load_type21_kernel",
                return_value=_FakeType21Kernel(),
            ),
        ):
            result = spk.calc_spk_body_position(t, CHIRON, FLG_SPEED)
        assert result is not None

    def test_type21_calc_returns_none(self):
        """Type 21 kernel loads but calc returns None -> None (1402->1405)."""
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 20002060)
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch("libephemeris.spk.get_spk_coverage", return_value=None),
            patch(
                "libephemeris.spk._load_type21_kernel",
                return_value=_FakeType21Kernel(),
            ),
            patch("libephemeris.spk._calc_type21_position", return_value=None),
        ):
            assert spk.calc_spk_body_position(t, CHIRON, 0) is None

    def test_type21_kernel_none(self):
        """Type 21 kernel load fails -> None (1404-1405)."""
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 20002060)
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=21),
            patch("libephemeris.spk.get_spk_coverage", return_value=None),
            patch("libephemeris.spk._load_type21_kernel", return_value=None),
        ):
            assert spk.calc_spk_body_position(t, CHIRON, 0) is None

    def test_type2_kernel_none(self):
        """Type 2/3 with no cached kernel -> None (1408-1410)."""
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 4)
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        with patch("libephemeris.spk._detect_spk_type", return_value=2):
            assert spk.calc_spk_body_position(t, CHIRON, 0) is None

    def test_type2_out_of_coverage_raises(self):
        """Type 2/3 jd outside coverage -> EphemerisRangeError (1416-1422)."""
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 4)
        state._SPK_KERNELS["f.bsp"] = MagicMock()
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=2),
            patch(
                "libephemeris.spk.get_spk_coverage",
                return_value=(2460000.0, 2470000.0),
            ),
        ):
            with pytest.raises(EphemerisRangeError):
                spk.calc_spk_body_position(t, CHIRON, 0)

    def test_type2_in_coverage(self):
        """Type 2/3 with jd inside a real coverage range (1416->1425)."""
        planets = state.get_planets()
        mars = planets[4]
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 4)
        state._SPK_KERNELS["f.bsp"] = _FakeKernel(mars)
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=2),
            patch(
                "libephemeris.spk.get_spk_coverage",
                return_value=(2440000.0, 2470000.0),
            ),
        ):
            result = spk.calc_spk_body_position(t, CHIRON, 0)
        assert result is not None

    def test_type2_speed_longitude_wrap_positive(self):
        """speed_lon overshoots +180/(2dt) -> subtract correction (1504-1505)."""
        t0 = 2451545.0
        # before-seam longitude ~0.002, after-seam ~359.998 -> next-prev large +.
        target = _WrapTarget(t0, lon_before=0.002, lon_after=359.998)
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 4)
        state._SPK_KERNELS["f.bsp"] = _FakeKernel(target)
        ts = state.get_timescale()
        t = ts.tt_jd(t0)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=2),
            patch("libephemeris.spk.get_spk_coverage", return_value=None),
        ):
            result = spk.calc_spk_body_position(
                t, CHIRON, FLG_SPEED | FLG_TRUEPOS | FLG_NOABERR
            )
        assert result is not None

    def test_type2_speed_longitude_wrap_negative(self):
        """speed_lon undershoots -180/(2dt) -> add correction (1506-1507)."""
        t0 = 2451545.0
        target = _WrapTarget(t0, lon_before=359.998, lon_after=0.002)
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 4)
        state._SPK_KERNELS["f.bsp"] = _FakeKernel(target)
        ts = state.get_timescale()
        t = ts.tt_jd(t0)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=2),
            patch("libephemeris.spk.get_spk_coverage", return_value=None),
        ):
            result = spk.calc_spk_body_position(
                t, CHIRON, FLG_SPEED | FLG_TRUEPOS | FLG_NOABERR
            )
        assert result is not None

    def test_type2_target_keyerror(self):
        """Type 2/3 kernel[naif_id] KeyError -> None (1426-1428)."""
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 999999)
        state._SPK_KERNELS["f.bsp"] = _FakeKernelKeyError()
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=2),
            patch("libephemeris.spk.get_spk_coverage", return_value=None),
        ):
            assert spk.calc_spk_body_position(t, CHIRON, 0) is None

    def test_type2_geocentric_full(self):
        """Type 2/3 full geocentric pipeline using DE440 Mars as the target
        (1430-1516, FLG_SPEED)."""
        planets = state.get_planets()
        mars = planets[4]  # Mars barycenter VectorFunction
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 4)
        state._SPK_KERNELS["f.bsp"] = _FakeKernel(mars)
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=2),
            patch("libephemeris.spk.get_spk_coverage", return_value=None),
        ):
            result = spk.calc_spk_body_position(t, CHIRON, FLG_SPEED)
        assert result is not None
        assert 0 <= result[0] < 360

    def test_type2_heliocentric_truepos_noaberr(self):
        """Type 2/3 heliocentric + truepos + noaberr branches (1437, 1451-1452,
        1464 break)."""
        planets = state.get_planets()
        mars = planets[4]
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 4)
        state._SPK_KERNELS["f.bsp"] = _FakeKernel(mars)
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=2),
            patch("libephemeris.spk.get_spk_coverage", return_value=None),
        ):
            result = spk.calc_spk_body_position(
                t, CHIRON, FLG_HELCTR | FLG_TRUEPOS | FLG_NOABERR
            )
        assert result is not None

    def test_type2_sidereal(self):
        """Type 2/3 sidereal correction branch (1510-1514)."""
        planets = state.get_planets()
        mars = planets[4]
        state._SPK_BODY_MAP[CHIRON] = ("f.bsp", 4)
        state._SPK_KERNELS["f.bsp"] = _FakeKernel(mars)
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        with (
            patch("libephemeris.spk._detect_spk_type", return_value=2),
            patch("libephemeris.spk.get_spk_coverage", return_value=None),
        ):
            result = spk.calc_spk_body_position(t, CHIRON, FLG_SIDEREAL)
        assert result is not None

    def test_not_registered(self):
        """Unregistered body -> None (1376-1377)."""
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        assert spk.calc_spk_body_position(t, CHIRON, 0) is None


# Reference numpy so linters don't flag the import as unused.
_ = np.array([0.0])
