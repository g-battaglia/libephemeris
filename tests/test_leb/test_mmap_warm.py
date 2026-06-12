"""Tests for mmap warm() and _madvise_ranges."""

from __future__ import annotations

import mmap
import os
import tempfile

import pytest

from libephemeris.leb_format import _madvise_ranges, _PAGE_SIZE, _PAGE_MASK


# ---------------------------------------------------------------------------
# _madvise_ranges
# ---------------------------------------------------------------------------


class TestMadviseRanges:
    def _make_mmap(self, size: int = 64 * 1024):
        """Create a temporary mmap for testing."""
        fd, path = tempfile.mkstemp()
        try:
            os.write(fd, b"\x00" * size)
            os.lseek(fd, 0, os.SEEK_SET)
            mm = mmap.mmap(fd, 0, access=mmap.ACCESS_READ)
        finally:
            os.close(fd)
            os.unlink(path)
        return mm

    def test_empty_ranges(self):
        mm = self._make_mmap()
        _madvise_ranges(mm, [])
        mm.close()

    def test_single_range(self):
        mm = self._make_mmap()
        _madvise_ranges(mm, [(100, 200)])
        mm.close()

    def test_adjacent_ranges_merged(self):
        mm = self._make_mmap()
        _madvise_ranges(mm, [(100, 50), (150, 50), (200, 50)])
        mm.close()

    def test_non_overlapping_ranges(self):
        size = 64 * 1024
        mm = self._make_mmap(size)
        _madvise_ranges(mm, [(0, 100), (size - 100, 100)])
        mm.close()

    def test_range_clamped_to_mmap_length(self):
        size = 8 * 1024
        mm = self._make_mmap(size)
        _madvise_ranges(mm, [(size - 10, 1000)])
        mm.close()

    def test_page_alignment(self):
        assert _PAGE_SIZE > 0
        assert _PAGE_SIZE & (_PAGE_SIZE - 1) == 0, "PAGE_SIZE must be power of 2"
        assert _PAGE_MASK == ~(_PAGE_SIZE - 1)


# ---------------------------------------------------------------------------
# LEBReader.warm()
# ---------------------------------------------------------------------------


class TestLEBReaderWarm:
    def test_warm_no_error(self, leb_reader):
        jd_start, jd_end = leb_reader.jd_range
        leb_reader.warm(jd_start, jd_end)

    def test_warm_subset(self, leb_reader):
        jd_start, jd_end = leb_reader.jd_range
        mid = (jd_start + jd_end) / 2
        leb_reader.warm(mid - 100, mid + 100)

    def test_warm_outside_range(self, leb_reader):
        jd_start, _ = leb_reader.jd_range
        leb_reader.warm(jd_start - 100000, jd_start - 99000)

    def test_no_madvise_willneed_in_init(self):
        """Verify that __init__ does not call madvise(MADV_WILLNEED)."""
        import libephemeris.leb_reader as mod
        import inspect

        source = inspect.getsource(mod.LEBReader.__init__)
        assert "MADV_WILLNEED" not in source


# ---------------------------------------------------------------------------
# LEB2Reader.warm()
# ---------------------------------------------------------------------------


class TestLEB2ReaderWarm:
    @pytest.fixture
    def leb2_file(self):
        """Find an available LEB2 file for testing."""
        data_dir = os.path.join(
            os.path.dirname(__file__), "..", "..", "data", "leb2"
        )
        for name in ("base_core.leb2", "medium_core.leb2", "extended_core.leb2"):
            path = os.path.join(data_dir, name)
            if os.path.exists(path):
                return path
        pytest.skip("No LEB2 data files available")

    @pytest.fixture
    def leb2_reader(self, leb2_file):
        from libephemeris.leb2_reader import LEB2Reader

        reader = LEB2Reader(leb2_file)
        yield reader
        reader.close()

    def test_warm_no_error(self, leb2_reader):
        jd_start, jd_end = leb2_reader.jd_range
        mid = (jd_start + jd_end) / 2
        leb2_reader.warm(mid - 365, mid + 365)

    def test_warm_modern_range(self, leb2_reader):
        from libephemeris.state import _year_to_jd

        leb2_reader.warm(_year_to_jd(1800), _year_to_jd(2200))

    def test_warm_outside_range(self, leb2_reader):
        jd_start, _ = leb2_reader.jd_range
        leb2_reader.warm(jd_start - 100000, jd_start - 99000)

    def test_no_madvise_willneed_in_init(self):
        """Verify that __init__ does not call madvise(MADV_WILLNEED)."""
        import libephemeris.leb2_reader as mod
        import inspect

        source = inspect.getsource(mod.LEB2Reader.__init__)
        assert "MADV_WILLNEED" not in source


# ---------------------------------------------------------------------------
# CompositeLEBReader.warm()
# ---------------------------------------------------------------------------


class TestCompositeLEBReaderWarm:
    @pytest.fixture
    def composite_reader(self):
        data_dir = os.path.join(
            os.path.dirname(__file__), "..", "..", "data", "leb2"
        )
        core = os.path.join(data_dir, "base_core.leb2")
        if not os.path.exists(core):
            for name in ("medium_core.leb2", "extended_core.leb2"):
                core = os.path.join(data_dir, name)
                if os.path.exists(core):
                    break
            else:
                pytest.skip("No LEB2 data files available")

        from libephemeris.leb_composite import CompositeLEBReader

        reader = CompositeLEBReader.from_file_with_companions(core)
        yield reader
        reader.close()

    def test_warm_delegates(self, composite_reader):
        from libephemeris.state import _year_to_jd

        composite_reader.warm(_year_to_jd(1900), _year_to_jd(2100))

    def test_warm_covers_all_readers(self, composite_reader):
        jd_start, jd_end = composite_reader.jd_range
        composite_reader.warm(jd_start, jd_end)


# ---------------------------------------------------------------------------
# TOML config (mmap_preload default)
# ---------------------------------------------------------------------------


class TestMmapPreloadConfig:
    def test_default_is_false(self):
        from libephemeris._config_toml import get_bool, reset

        reset()
        assert get_bool("mmap_preload") is None


# ---------------------------------------------------------------------------
# cool() — page cache release
# ---------------------------------------------------------------------------


class TestLEBReaderCool:
    def test_cool_no_error(self, leb_reader):
        leb_reader.cool()

    def test_cool_idempotent(self, leb_reader):
        leb_reader.cool()
        leb_reader.cool()

    def test_cool_after_close(self, test_leb_file):
        from libephemeris.leb_reader import LEBReader

        reader = LEBReader(test_leb_file)
        reader.close()
        reader.cool()

    def test_cool_does_not_clear_eval_cache(self, leb_reader):
        jd_start, jd_end = leb_reader.jd_range
        mid = (jd_start + jd_end) / 2
        leb_reader.eval_body(0, mid)
        assert len(leb_reader._eval_cache) > 0
        leb_reader.cool()
        assert len(leb_reader._eval_cache) > 0

    def test_calc_works_after_cool(self, leb_reader):
        jd_start, jd_end = leb_reader.jd_range
        mid = (jd_start + jd_end) / 2
        pos1, vel1 = leb_reader.eval_body(0, mid)
        leb_reader.cool()
        pos2, vel2 = leb_reader.eval_body(0, mid)
        assert pos1 == pos2
        assert vel1 == vel2


class TestLEB2ReaderCool:
    @pytest.fixture
    def leb2_file(self):
        data_dir = os.path.join(
            os.path.dirname(__file__), "..", "..", "data", "leb2"
        )
        for name in ("base_core.leb2", "medium_core.leb2", "extended_core.leb2"):
            path = os.path.join(data_dir, name)
            if os.path.exists(path):
                return path
        pytest.skip("No LEB2 data files available")

    @pytest.fixture
    def leb2_reader(self, leb2_file):
        from libephemeris.leb2_reader import LEB2Reader

        reader = LEB2Reader(leb2_file)
        yield reader
        reader.close()

    def test_cool_no_error(self, leb2_reader):
        leb2_reader.cool()

    def test_cool_idempotent(self, leb2_reader):
        leb2_reader.cool()
        leb2_reader.cool()

    def test_cool_after_close(self, leb2_file):
        from libephemeris.leb2_reader import LEB2Reader

        reader = LEB2Reader(leb2_file)
        reader.close()
        reader.cool()

    def test_cool_does_not_clear_chunk_cache(self, leb2_reader):
        jd_start, jd_end = leb2_reader.jd_range
        mid = (jd_start + jd_end) / 2
        for body_id in list(leb2_reader._bodies.keys())[:1]:
            leb2_reader.eval_body(body_id, mid)
        has_cache = len(leb2_reader._chunk_cache) > 0 or len(leb2_reader._eval_cache) > 0
        assert has_cache
        leb2_reader.cool()
        has_cache_after = len(leb2_reader._chunk_cache) > 0 or len(leb2_reader._eval_cache) > 0
        assert has_cache_after


class TestCompositeLEBReaderCool:
    @pytest.fixture
    def composite_reader(self):
        data_dir = os.path.join(
            os.path.dirname(__file__), "..", "..", "data", "leb2"
        )
        for name in ("base_core.leb2", "medium_core.leb2", "extended_core.leb2"):
            core = os.path.join(data_dir, name)
            if os.path.exists(core):
                from libephemeris.leb_composite import CompositeLEBReader

                reader = CompositeLEBReader.from_file_with_companions(core)
                yield reader
                reader.close()
                return
        pytest.skip("No LEB2 data files available")

    def test_cool_delegates(self, composite_reader):
        composite_reader.cool()

    def test_cool_idempotent(self, composite_reader):
        composite_reader.cool()
        composite_reader.cool()


# ---------------------------------------------------------------------------
# release_data_cache()
# ---------------------------------------------------------------------------


class TestReleaseDataCache:
    def test_no_error(self, tmp_path, monkeypatch):
        from libephemeris.state import release_data_cache

        (tmp_path / "test.leb2").write_bytes(b"\x00" * 1024)
        monkeypatch.setenv("LIBEPHEMERIS_DATA_DIR", str(tmp_path))
        release_data_cache()

    def test_no_error_missing_dir(self, monkeypatch):
        from libephemeris.state import release_data_cache

        monkeypatch.setenv("LIBEPHEMERIS_DATA_DIR", "/nonexistent/path")
        release_data_cache()

    def test_no_error_without_posix_fadvise(self, monkeypatch):
        from libephemeris.state import release_data_cache

        monkeypatch.delattr(os, "posix_fadvise", raising=False)
        release_data_cache()

    def test_skips_non_regular_files(self, tmp_path, monkeypatch):
        from libephemeris.state import release_data_cache

        (tmp_path / "regular.dat").write_bytes(b"\x00" * 512)
        if not hasattr(os, "mkfifo"):
            pytest.skip("os.mkfifo is not supported on this platform")
        os.mkfifo(str(tmp_path / "fifo.dat"))
        monkeypatch.setenv("LIBEPHEMERIS_DATA_DIR", str(tmp_path))
        release_data_cache()
        os.unlink(str(tmp_path / "fifo.dat"))
