"""Regression tests for two SPK fixes.

Fix 1 (spk.py ``_calc_type21_position``): the FLG_SPEED longitude/latitude
speed for type-21 bodies (asteroids/TNOs) is now the true time-derivative of
the *of-date* longitude/latitude. Previously the speed was computed
analytically in the fixed J2000 ecliptic while the position was precessed to
the equinox of date and had nutation added, so the of-date frame-rotation rate
(precession + nutation of the equinox, ~0.2"/day) was missing from the reported
motion.

Fix 2 (spk_auto.py ``_find_covering_spk``): cache coverage detection now reads
the real time span from each kernel's segment metadata instead of parsing the
filename, so all on-disk naming schemes are recognised:

    - body_<jdstart>_<jdend>.bsp  (_generate_spk_cache_filename)
    - body_<YYYYMM>_<YYYYMM>.bsp  (download_spk)
    - body_<md5hash>.bsp          (AutoSpkConfig.get_cache_filename)

The previous filename-parsing logic recognised only the first form, so
YYYYMM- and hash-named kernels were never seen as covering and triggered
redundant re-downloads.

Tests that need a real type-21 kernel or a shippable .bsp are guarded with
``pytest.skip`` so they no-op when those files are unavailable (e.g. CI).
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Optional

import pytest

from libephemeris import spk, spk_auto, state
from libephemeris.constants import FLG_J2000, FLG_SPEED


# ---------------------------------------------------------------------------
# Helpers for locating usable kernels on disk (skip-guarded, never network).
# ---------------------------------------------------------------------------


def _find_any_valid_bsp() -> Optional[str]:
    """Return the path to any shipped/cached valid .bsp, or None.

    Prefers the repository's de440s.bsp (wide coverage, always type 2/3),
    then any cached kernel. Used as a stand-in kernel for the filename-scheme
    coverage tests -- its data content is irrelevant, only that
    get_spk_coverage() can read a real time span from it.
    """
    candidates = [
        Path(__file__).resolve().parents[1] / "de440s.bsp",
        Path.cwd() / "de440s.bsp",
    ]
    data_dir = None
    try:
        from libephemeris.state import _get_data_dir

        data_dir = _get_data_dir()
    except Exception:
        data_dir = None
    if data_dir:
        candidates.append(Path(data_dir) / "de440s.bsp")

    for c in candidates:
        if c.is_file():
            return str(c)
    return None


def _find_type21_kernel() -> Optional[tuple[str, int, tuple[float, float]]]:
    """Discover a cached SPK type-21 kernel.

    Returns (path, small_body_naif_id, (start_jd, end_jd)) or None if no
    type-21 kernel is available. Honours the ``LIBEPHEMERIS_TEST_TYPE21_BSP``
    environment override.
    """
    override = os.environ.get("LIBEPHEMERIS_TEST_TYPE21_BSP")
    search_files: list[str] = []
    if override and os.path.isfile(override):
        search_files.append(override)

    search_dirs = [spk_auto.DEFAULT_AUTO_SPK_DIR]
    try:
        cache_dir = state.get_spk_cache_dir()
        if cache_dir:
            search_dirs.append(cache_dir)
    except Exception:
        pass

    for d in search_dirs:
        if not d or not os.path.isdir(d):
            continue
        for name in sorted(os.listdir(d)):
            if name.startswith("._") or not name.endswith(".bsp"):
                continue
            search_files.append(os.path.join(d, name))

    for path in search_files:
        try:
            if spk._detect_spk_type(path) != 21:
                continue
            targets = spk._get_spk_targets(path)
            small = [t for t in targets if t > 1_000_000]
            if not small:
                continue
            coverage = spk.get_spk_coverage(path)
            if coverage is None:
                continue
            return (path, small[0], coverage)
        except Exception:
            continue
    return None


# ---------------------------------------------------------------------------
# Fix 2: filename-scheme-independent coverage recognition.
# ---------------------------------------------------------------------------


class TestFindCoveringSpkSchemes:
    """_find_covering_spk must recognise coverage regardless of filename."""

    def _make_cache(self, tmpdir: str, filename: str, src_bsp: str) -> str:
        """Place a stand-in kernel under ``filename`` in ``tmpdir``.

        Uses a symlink when possible (fast; avoids copying a large kernel) and
        falls back to a real copy on platforms without symlink support.
        """
        dest = os.path.join(tmpdir, filename)
        try:
            os.symlink(os.path.abspath(src_bsp), dest)
        except (OSError, NotImplementedError):
            import shutil

            shutil.copy(src_bsp, dest)
        return dest

    def _covered_range(self, src_bsp: str) -> tuple[float, float]:
        cov = spk.get_spk_coverage(src_bsp)
        assert cov is not None, "stand-in kernel must expose coverage"
        start, end = cov
        # Request a range comfortably inside the kernel's real coverage.
        span = end - start
        return (start + 0.25 * span, start + 0.75 * span)

    def test_yyyymm_named_file_recognised(self):
        """A YYYYMM_YYYYMM-named cached kernel is recognised as covering.

        This is the schema written by download_spk(); its numeric tokens
        (e.g. 200001) are NOT Julian Days, so the old filename-parsing path
        wrongly rejected it and re-downloaded.
        """
        src = _find_any_valid_bsp()
        if src is None:
            pytest.skip("no shippable .bsp available to build a cache fixture")

        jd_start, jd_end = self._covered_range(src)
        with tempfile.TemporaryDirectory() as tmp:
            # Tokens are calendar YYYYMM, deliberately unrelated to the real
            # (much larger) Julian-Day coverage of the kernel.
            self._make_cache(tmp, "testbody_200001_201001.bsp", src)
            found = spk_auto._find_covering_spk("testbody", jd_start, jd_end, tmp)
            assert found is not None
            assert os.path.basename(found) == "testbody_200001_201001.bsp"

    def test_md5_hash_named_file_recognised(self):
        """A body_<md5hash>.bsp cached kernel is recognised as covering.

        This is the schema written by AutoSpkConfig; the hash token is not
        numeric, so the old path skipped it entirely.
        """
        src = _find_any_valid_bsp()
        if src is None:
            pytest.skip("no shippable .bsp available to build a cache fixture")

        jd_start, jd_end = self._covered_range(src)
        with tempfile.TemporaryDirectory() as tmp:
            self._make_cache(tmp, "testbody_deadbeef.bsp", src)
            found = spk_auto._find_covering_spk("testbody", jd_start, jd_end, tmp)
            assert found is not None
            assert os.path.basename(found) == "testbody_deadbeef.bsp"

    def test_jd_range_named_file_still_recognised(self):
        """The original body_<jdstart>_<jdend>.bsp schema keeps working."""
        src = _find_any_valid_bsp()
        if src is None:
            pytest.skip("no shippable .bsp available to build a cache fixture")

        cov = spk.get_spk_coverage(src)
        assert cov is not None
        start, end = cov
        fname = f"testbody_{int(start)}_{int(end)}.bsp"
        jd_start, jd_end = self._covered_range(src)
        with tempfile.TemporaryDirectory() as tmp:
            self._make_cache(tmp, fname, src)
            found = spk_auto._find_covering_spk("testbody", jd_start, jd_end, tmp)
            assert found is not None
            assert os.path.basename(found) == fname

    def test_out_of_coverage_returns_none(self):
        """A request outside the kernel's real coverage is not matched."""
        src = _find_any_valid_bsp()
        if src is None:
            pytest.skip("no shippable .bsp available to build a cache fixture")

        cov = spk.get_spk_coverage(src)
        assert cov is not None
        start, end = cov
        with tempfile.TemporaryDirectory() as tmp:
            self._make_cache(tmp, "testbody_200001_201001.bsp", src)
            # Request a range that starts before the kernel's real coverage.
            found = spk_auto._find_covering_spk("testbody", start - 10_000.0, end, tmp)
            assert found is None

    def test_other_body_not_matched(self):
        """Only files whose name prefix matches the body are considered."""
        src = _find_any_valid_bsp()
        if src is None:
            pytest.skip("no shippable .bsp available to build a cache fixture")

        jd_start, jd_end = self._covered_range(src)
        with tempfile.TemporaryDirectory() as tmp:
            self._make_cache(tmp, "otherbody_200001_201001.bsp", src)
            found = spk_auto._find_covering_spk("testbody", jd_start, jd_end, tmp)
            assert found is None


# ---------------------------------------------------------------------------
# Fix 1: type-21 of-date speed self-consistency.
# ---------------------------------------------------------------------------


class TestType21SpeedSelfConsistency:
    """The reported FLG_SPEED motion equals d(of-date position)/dt."""

    _TEST_IPL = 987654  # arbitrary id, only used to register the kernel

    @pytest.fixture()
    def registered_type21(self):
        info = _find_type21_kernel()
        if info is None:
            pytest.skip("no cached SPK type-21 kernel available")
        path, naif_id, coverage = info
        if spk._get_spktype21() is None:
            pytest.skip("spktype21 support not installed")
        spk.register_spk_body(self._TEST_IPL, path, naif_id)
        try:
            yield coverage
        finally:
            state._SPK_BODY_MAP.pop(self._TEST_IPL, None)

    def _lonlat_at(self, jd_tt: float, iflag: int) -> tuple[float, float, float]:
        ts = state.get_timescale()
        res = spk.calc_spk_body_position(ts.tt_jd(jd_tt), self._TEST_IPL, iflag)
        assert res is not None, "type-21 position computation returned None"
        return res[0], res[1], res[2]

    @staticmethod
    def _unwrap_deg(delta: float) -> float:
        if delta > 180.0:
            delta -= 360.0
        elif delta < -180.0:
            delta += 360.0
        return delta

    def test_speed_matches_ofdate_finite_difference(self, registered_type21):
        """Reported dlon/dlat == central finite-diff of of-date lon/lat.

        The internal speed uses a 1 s step; here we cross-check with an
        independent, larger step (~0.02 d), so agreement to <0.01"/day is a
        genuine self-consistency check, not a tautology.
        """
        start_jd, end_jd = registered_type21
        jd_mid = 0.5 * (start_jd + end_jd)
        ts = state.get_timescale()

        res = spk.calc_spk_body_position(ts.tt_jd(jd_mid), self._TEST_IPL, FLG_SPEED)
        assert res is not None
        _, _, _, dlon, dlat, ddist = res

        h = 0.02  # days
        lon_p, lat_p, dist_p = self._lonlat_at(jd_mid - h, 0)
        lon_n, lat_n, dist_n = self._lonlat_at(jd_mid + h, 0)

        fd_lon = self._unwrap_deg(lon_n - lon_p) / (2.0 * h)
        fd_lat = (lat_n - lat_p) / (2.0 * h)
        fd_dist = (dist_n - dist_p) / (2.0 * h)

        # Tolerances in the reported units (deg/day for lon/lat). 0.01"/day.
        tol_deg_per_day = 0.01 / 3600.0
        assert abs(dlon - fd_lon) < tol_deg_per_day, (
            f"dlon {dlon} vs finite-diff {fd_lon} "
            f"({abs(dlon - fd_lon) * 3600:.5f} arcsec/day)"
        )
        assert abs(dlat - fd_lat) < tol_deg_per_day, (
            f"dlat {dlat} vs finite-diff {fd_lat} "
            f"({abs(dlat - fd_lat) * 3600:.5f} arcsec/day)"
        )
        # Distance rate is much smaller in magnitude; use an absolute AU/day tol.
        assert abs(ddist - fd_dist) < 1e-6

    def test_frame_rotation_rate_present(self, registered_type21):
        """Of-date longitude speed differs from the J2000-frame speed.

        The whole bug was that the two were identical (the of-date speed was
        just the J2000 analytic derivative). The difference must now be a
        precession/nutation-scale rate: a few tenths of an arcsec/day. If the
        frame-rotation term were dropped again, this difference collapses to ~0.
        """
        start_jd, end_jd = registered_type21
        jd_mid = 0.5 * (start_jd + end_jd)
        ts = state.get_timescale()

        of_date = spk.calc_spk_body_position(
            ts.tt_jd(jd_mid), self._TEST_IPL, FLG_SPEED
        )
        j2000 = spk.calc_spk_body_position(
            ts.tt_jd(jd_mid), self._TEST_IPL, FLG_SPEED | FLG_J2000
        )
        assert of_date is not None and j2000 is not None

        diff_arcsec_per_day = abs(of_date[3] - j2000[3]) * 3600.0
        # General precession in longitude is ~0.14"/day; nutation modulates it.
        assert 0.02 < diff_arcsec_per_day < 0.5, (
            f"frame-rotation contribution {diff_arcsec_per_day:.4f} arcsec/day "
            "outside the expected precession/nutation band"
        )
