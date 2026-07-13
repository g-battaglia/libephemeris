"""Tests for LEB2Reader and CompositeLEBReader."""

import glob
import os
import struct
import tempfile

import numpy as np
import pytest

from libephemeris.leb_compression import compress_body, compute_mantissa_bits
from libephemeris.leb_format import (
    segment_byte_size,
)
from libephemeris.leb_reader import LEBReader, open_leb

# Skip all tests if base LEB file doesn't exist
LEB1_PATH = os.path.join(
    os.path.dirname(__file__), "..", "..", "data", "leb", "ephemeris_base.leb"
)
LEB2_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "data", "leb2")

requires_leb1 = pytest.mark.skipif(
    not os.path.isfile(LEB1_PATH), reason="Base LEB1 file not available"
)
requires_leb2 = pytest.mark.skipif(
    not os.path.isdir(LEB2_DIR), reason="LEB2 files not available"
)


@pytest.fixture(scope="module")
def single_tier_dir(tmp_path_factory):
    """One tier's independently generated LEB2 group files only."""
    src = sorted(glob.glob(os.path.join(LEB2_DIR, "base_*.leb2")))
    if not src:
        pytest.skip("base_*.leb2 files not found")
    d = tmp_path_factory.mktemp("leb2_single_tier")
    for s in src:
        os.symlink(os.path.abspath(s), os.path.join(d, os.path.basename(s)))
    return str(d)


def _create_mini_leb2(body_ids=(0,)):
    """Create a minimal LEB2 file from the base LEB1 file for testing."""
    reader1 = LEBReader(LEB1_PATH)

    body_entries = []
    for bid in body_ids:
        if not reader1.has_body(bid):
            continue
        entry = reader1._bodies[bid]
        seg_size = segment_byte_size(entry.degree, entry.components)
        raw_size = seg_size * entry.segment_count
        raw = struct.pack_from_mmap(reader1._mm, entry.data_offset, raw_size)
        arr = np.frombuffer(raw, dtype=np.float64).reshape(
            entry.segment_count, entry.components, entry.degree + 1
        )
        bits = compute_mantissa_bits(arr)
        blob = compress_body(
            raw, entry.segment_count, entry.degree, entry.components, bits
        )
        body_entries.append((entry, blob, raw_size))

    reader1.close()
    return body_entries


def struct_pack_from_mmap(mm, offset, size):
    return bytes(mm[offset : offset + size])


@requires_leb1
class TestOpenLeb:
    def test_open_leb1(self):
        reader = open_leb(LEB1_PATH)
        assert isinstance(reader, LEBReader)
        reader.close()

    def test_open_invalid_magic(self):
        with tempfile.NamedTemporaryFile(suffix=".leb", delete=False) as f:
            f.write(b"XXXX" + b"\x00" * 60)
            f.flush()
            path = f.name
        try:
            with pytest.raises(ValueError, match="Unknown LEB format"):
                open_leb(path)
        finally:
            os.unlink(path)


@requires_leb2
class TestLEB2Reader:
    def test_open_and_read(self):
        from libephemeris.leb2_reader import LEB2Reader

        path = os.path.join(LEB2_DIR, "base_core.leb2")
        if not os.path.isfile(path):
            pytest.skip("base_core.leb2 not found")

        reader = LEB2Reader(path)
        assert reader.has_body(0)  # Sun
        assert reader.has_body(1)  # Moon

        pos, vel = reader.eval_body(0, 2451545.0)
        assert len(pos) == 3
        assert len(vel) == 3
        assert abs(pos[0]) < 0.1  # Sun ICRS x near 0

        reader.close()

    def test_context_manager(self):
        from libephemeris.leb2_reader import LEB2Reader

        path = os.path.join(LEB2_DIR, "base_core.leb2")
        if not os.path.isfile(path):
            pytest.skip("base_core.leb2 not found")

        with LEB2Reader(path) as reader:
            pos, _ = reader.eval_body(0, 2451545.0)
            assert len(pos) == 3

    def test_nutation(self):
        from libephemeris.leb2_reader import LEB2Reader

        path = os.path.join(LEB2_DIR, "base_core.leb2")
        if not os.path.isfile(path):
            pytest.skip("base_core.leb2 not found")

        with LEB2Reader(path) as reader:
            dpsi, deps = reader.eval_nutation(2451545.0)
            assert abs(dpsi) < 0.001  # radians
            assert abs(deps) < 0.001

    def test_delta_t(self):
        from libephemeris.leb2_reader import LEB2Reader

        path = os.path.join(LEB2_DIR, "base_core.leb2")
        if not os.path.isfile(path):
            pytest.skip("base_core.leb2 not found")

        with LEB2Reader(path) as reader:
            dt = reader.delta_t(2451545.0)
            assert 0 < dt < 0.01  # days

    def test_body_not_found(self):
        from libephemeris.leb2_reader import LEB2Reader

        path = os.path.join(LEB2_DIR, "base_core.leb2")
        if not os.path.isfile(path):
            pytest.skip("base_core.leb2 not found")

        with LEB2Reader(path) as reader:
            with pytest.raises(KeyError):
                reader.eval_body(999, 2451545.0)

    @requires_leb1
    def test_precision_vs_leb1(self):
        """Verify LEB2 results are within 0.003\" of LEB1 for core bodies."""
        from libephemeris.leb2_reader import LEB2Reader

        path = os.path.join(LEB2_DIR, "base_core.leb2")
        if not os.path.isfile(path):
            pytest.skip("base_core.leb2 not found")

        reader1 = LEBReader(LEB1_PATH)
        reader2 = LEB2Reader(path)

        rng = np.random.default_rng(42)
        jds = rng.uniform(2396760, 2506330, 100)

        for bid in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14]:
            for jd in jds:
                pos1, _ = reader1.eval_body(bid, float(jd))
                pos2, _ = reader2.eval_body(bid, float(jd))
                err = max(abs(a - b) for a, b in zip(pos1, pos2))
                err_arcsec = err * 206265
                assert err_arcsec < 0.01, (
                    f'Body {bid} at JD {jd}: error {err_arcsec:.4f}" > 0.01"'
                )

        reader1.close()
        reader2.close()


requires_leb2_asteroids = pytest.mark.skipif(
    not os.path.isfile(os.path.join(LEB2_DIR, "base_asteroids.leb2")),
    reason="independently generated asteroid companion is not available",
)
requires_leb2_all_companions = pytest.mark.skipif(
    not all(
        os.path.isfile(os.path.join(LEB2_DIR, f"base_{group}.leb2"))
        for group in ("asteroids", "apogee")
    ),
    reason="all independently generated LEB2 companion files are not available",
)


@requires_leb2
class TestCompositeLEBReader:
    def test_from_directory(self, single_tier_dir):
        from libephemeris.leb_composite import CompositeLEBReader

        reader = CompositeLEBReader.from_directory(single_tier_dir)
        assert reader.has_body(0)  # Sun (core always present)
        reader.close()

    def test_from_directory_mixed_tiers_raises(self):
        """A directory containing multiple tiers cannot be merged."""
        from libephemeris.leb_composite import CompositeLEBReader

        tiers = {
            os.path.basename(path).split("_", maxsplit=1)[0]
            for path in glob.glob(os.path.join(LEB2_DIR, "*.leb2"))
        }
        if len(tiers) < 2:
            pytest.skip("local LEB2 directory contains only one reviewed tier")
        with pytest.raises(ValueError, match="mixed tiers"):
            CompositeLEBReader.from_directory(LEB2_DIR)

    @requires_leb2_all_companions
    def test_from_directory_all_groups(self, single_tier_dir):
        from libephemeris.leb_composite import CompositeLEBReader

        reader = CompositeLEBReader.from_directory(single_tier_dir)
        assert reader.has_body(15)  # Chiron (asteroids)
        assert reader.has_body(13)  # OscuApog (apogee)
        reader.close()

    def test_from_file_with_companions(self):
        from libephemeris.leb_composite import CompositeLEBReader

        path = os.path.join(LEB2_DIR, "base_core.leb2")
        if not os.path.isfile(path):
            pytest.skip("base_core.leb2 not found")
        reader = CompositeLEBReader.from_file_with_companions(path)
        assert reader.has_body(0)  # from primary
        reader.close()

    @requires_leb2_asteroids
    def test_from_file_with_companions_all(self):
        from libephemeris.leb_composite import CompositeLEBReader

        path = os.path.join(LEB2_DIR, "base_core.leb2")
        if not os.path.isfile(path):
            pytest.skip("base_core.leb2 not found")
        reader = CompositeLEBReader.from_file_with_companions(path)
        assert reader.has_body(15)  # from companion
        reader.close()

    @requires_leb2_asteroids
    def test_eval_body_cross_files(self, single_tier_dir):
        from libephemeris.leb_composite import CompositeLEBReader

        reader = CompositeLEBReader.from_directory(single_tier_dir)
        jd = 2451545.0

        # Bodies from different files
        pos_sun, _ = reader.eval_body(0, jd)  # core
        pos_chi, _ = reader.eval_body(15, jd)  # asteroids

        assert len(pos_sun) == 3
        assert len(pos_chi) == 3

        reader.close()

    def test_missing_body_raises(self, single_tier_dir):
        from libephemeris.leb_composite import CompositeLEBReader

        reader = CompositeLEBReader.from_directory(single_tier_dir)
        with pytest.raises(KeyError):
            reader.eval_body(999, 2451545.0)
        reader.close()
