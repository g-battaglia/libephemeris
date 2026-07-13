"""
Tests for the LEB file format module (leb_format.py).

Tests struct sizes, round-trip serialization, and body parameter table.
"""

from __future__ import annotations

import struct

import pytest

from libephemeris.leb_format import (
    BODY_ENTRY_FMT,
    BODY_ENTRY_SIZE,
    BODY_PARAMS,
    COORD_ECLIPTIC,
    COORD_HELIO_ECL,
    COORD_ICRS_BARY,
    COORD_ICRS_BARY_SYSTEM,
    DELTA_T_ENTRY_FMT,
    DELTA_T_ENTRY_SIZE,
    DELTA_T_HEADER_FMT,
    DELTA_T_HEADER_SIZE,
    HEADER_FMT,
    HEADER_SIZE,
    MAGIC,
    NUTATION_HEADER_FMT,
    NUTATION_HEADER_SIZE,
    SECTION_DIR_FMT,
    SECTION_DIR_SIZE,
    STAR_ENTRY_FMT,
    STAR_ENTRY_SIZE,
    VERSION,
    BodyEntry,
    FileHeader,
    NutationHeader,
    SectionEntry,
    StarEntry,
    read_body_entry,
    read_header,
    read_nutation_header,
    read_section_dir,
    read_star_entry,
    segment_byte_size,
    write_body_entry,
    write_header,
    write_nutation_header,
    write_section_dir,
    write_star_entry,
)


class TestStructSizes:
    """Verify struct format strings produce expected sizes."""

    @pytest.mark.unit
    def test_header_size(self):
        """Header struct must be exactly 64 bytes."""
        assert struct.calcsize(HEADER_FMT) == HEADER_SIZE

    @pytest.mark.unit
    def test_section_dir_size(self):
        """Section directory entry must be exactly 24 bytes."""
        assert struct.calcsize(SECTION_DIR_FMT) == SECTION_DIR_SIZE

    @pytest.mark.unit
    def test_body_entry_size(self):
        """Body index entry struct size must match BODY_ENTRY_SIZE."""
        assert struct.calcsize(BODY_ENTRY_FMT) == BODY_ENTRY_SIZE

    @pytest.mark.unit
    def test_star_entry_size(self):
        """Star catalog entry must be exactly 64 bytes."""
        assert struct.calcsize(STAR_ENTRY_FMT) == STAR_ENTRY_SIZE

    @pytest.mark.unit
    def test_nutation_header_size(self):
        """Nutation header struct size must match."""
        assert struct.calcsize(NUTATION_HEADER_FMT) == NUTATION_HEADER_SIZE

    @pytest.mark.unit
    def test_delta_t_entry_size(self):
        """Delta-T entry must be exactly 16 bytes (2 x f64)."""
        assert struct.calcsize(DELTA_T_ENTRY_FMT) == DELTA_T_ENTRY_SIZE
        assert DELTA_T_ENTRY_SIZE == 16

    @pytest.mark.unit
    def test_delta_t_header_size(self):
        """Delta-T header must be exactly 8 bytes."""
        assert struct.calcsize(DELTA_T_HEADER_FMT) == DELTA_T_HEADER_SIZE


class TestRoundTrip:
    """Test serialization round-trips for all data structures."""

    @pytest.mark.unit
    def test_header_round_trip(self):
        """FileHeader should survive write -> read cycle."""
        header = FileHeader(
            magic=MAGIC,
            version=VERSION,
            section_count=6,
            body_count=25,
            jd_start=2415020.5,
            jd_end=2488069.5,
            generation_epoch=2460000.5,
            flags=0,
        )
        buf = bytearray(HEADER_SIZE)
        write_header(buf, header)
        recovered = read_header(bytes(buf))

        assert recovered.magic == header.magic
        assert recovered.version == header.version
        assert recovered.section_count == header.section_count
        assert recovered.body_count == header.body_count
        assert recovered.jd_start == header.jd_start
        assert recovered.jd_end == header.jd_end
        assert recovered.generation_epoch == header.generation_epoch
        assert recovered.flags == header.flags

    @pytest.mark.unit
    def test_section_dir_round_trip(self):
        """SectionEntry should survive write -> read cycle."""
        entry = SectionEntry(
            section_id=1,
            offset=12345678,
            size=9876543,
        )
        buf = bytearray(SECTION_DIR_SIZE)
        write_section_dir(buf, 0, entry)
        recovered = read_section_dir(bytes(buf), 0)

        assert recovered.section_id == entry.section_id
        assert recovered.offset == entry.offset
        assert recovered.size == entry.size

    @pytest.mark.unit
    def test_body_entry_round_trip(self):
        """BodyEntry should survive write -> read cycle."""
        entry = BodyEntry(
            body_id=0,
            coord_type=COORD_ICRS_BARY,
            segment_count=2282,
            jd_start=2415020.5,
            jd_end=2488069.5,
            interval_days=32.0,
            degree=13,
            components=3,
            data_offset=1024,
        )
        buf = bytearray(BODY_ENTRY_SIZE)
        write_body_entry(buf, 0, entry)
        recovered = read_body_entry(bytes(buf), 0)

        assert recovered.body_id == entry.body_id
        assert recovered.coord_type == entry.coord_type
        assert recovered.segment_count == entry.segment_count
        assert recovered.jd_start == entry.jd_start
        assert recovered.jd_end == entry.jd_end
        assert recovered.interval_days == entry.interval_days
        assert recovered.degree == entry.degree
        assert recovered.components == entry.components
        assert recovered.data_offset == entry.data_offset

    @pytest.mark.unit
    def test_star_entry_round_trip(self):
        """StarEntry should survive write -> read cycle."""
        entry = StarEntry(
            star_id=42,
            ra_j2000=186.649563,
            dec_j2000=12.354321,
            pm_ra=0.001234,
            pm_dec=-0.000567,
            parallax=0.0089,
            rv=-12.5,
            magnitude=1.35,
        )
        buf = bytearray(STAR_ENTRY_SIZE)
        write_star_entry(buf, 0, entry)
        recovered = read_star_entry(bytes(buf), 0)

        assert recovered.star_id == entry.star_id
        assert abs(recovered.ra_j2000 - entry.ra_j2000) < 1e-10
        assert abs(recovered.dec_j2000 - entry.dec_j2000) < 1e-10
        assert abs(recovered.pm_ra - entry.pm_ra) < 1e-10
        assert abs(recovered.pm_dec - entry.pm_dec) < 1e-10
        assert abs(recovered.parallax - entry.parallax) < 1e-10
        assert abs(recovered.rv - entry.rv) < 1e-10
        assert abs(recovered.magnitude - entry.magnitude) < 1e-10

    @pytest.mark.unit
    def test_nutation_header_round_trip(self):
        """NutationHeader should survive write -> read cycle."""
        header = NutationHeader(
            jd_start=2415020.5,
            jd_end=2488069.5,
            interval_days=32.0,
            degree=16,
            components=2,
            segment_count=2282,
            reserved=0,
        )
        buf = bytearray(NUTATION_HEADER_SIZE)
        write_nutation_header(buf, 0, header)
        recovered = read_nutation_header(bytes(buf), 0)

        assert recovered.jd_start == header.jd_start
        assert recovered.jd_end == header.jd_end
        assert recovered.interval_days == header.interval_days
        assert recovered.degree == header.degree
        assert recovered.components == header.components
        assert recovered.segment_count == header.segment_count

    @pytest.mark.unit
    def test_multiple_body_entries(self):
        """Multiple body entries written sequentially should be readable."""
        entries = [
            BodyEntry(0, COORD_ICRS_BARY, 100, 2451545.0, 2488069.5, 32.0, 13, 3, 0),
            BodyEntry(1, COORD_ICRS_BARY, 200, 2451545.0, 2488069.5, 8.0, 13, 3, 4000),
            BodyEntry(10, COORD_ECLIPTIC, 150, 2451545.0, 2488069.5, 8.0, 13, 3, 8000),
        ]
        buf = bytearray(BODY_ENTRY_SIZE * len(entries))
        for i, entry in enumerate(entries):
            write_body_entry(buf, i * BODY_ENTRY_SIZE, entry)

        for i, expected in enumerate(entries):
            recovered = read_body_entry(bytes(buf), i * BODY_ENTRY_SIZE)
            assert recovered.body_id == expected.body_id
            assert recovered.coord_type == expected.coord_type
            assert recovered.data_offset == expected.data_offset


class TestBodyParams:
    """Test the BODY_PARAMS table."""

    @pytest.mark.unit
    def test_all_standard_planets_present(self):
        """Sun through Pluto and Earth must be in BODY_PARAMS."""
        for body_id in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14]:
            assert body_id in BODY_PARAMS, f"Body {body_id} missing from BODY_PARAMS"

    @pytest.mark.unit
    def test_lunar_nodes_present(self):
        """Lunar nodes and Lilith variants must be in BODY_PARAMS."""
        for body_id in [10, 11, 12, 13, 21, 22]:
            assert body_id in BODY_PARAMS, f"Body {body_id} missing from BODY_PARAMS"

    @pytest.mark.unit
    def test_asteroids_present(self):
        """Major asteroids must be in BODY_PARAMS."""
        for body_id in [15, 17, 18, 19, 20]:
            assert body_id in BODY_PARAMS, f"Body {body_id} missing from BODY_PARAMS"

    @pytest.mark.unit
    def test_reviewed_hypothetical_policy(self):
        """Fictitious-body models must never enter persisted LEB data."""
        assert not set(range(40, 59)).intersection(BODY_PARAMS)

    @pytest.mark.unit
    def test_coord_types_valid(self):
        """All coord types must be one of the defined types."""
        valid_types = {
            COORD_ICRS_BARY,
            COORD_ICRS_BARY_SYSTEM,
            COORD_ECLIPTIC,
            COORD_HELIO_ECL,
        }
        for body_id, params in BODY_PARAMS.items():
            _, _, coord_type, _ = params
            assert coord_type in valid_types, (
                f"Body {body_id} has invalid coord_type {coord_type}"
            )

    @pytest.mark.unit
    def test_components_always_3(self):
        """All bodies should have 3 components."""
        for body_id, params in BODY_PARAMS.items():
            _, _, _, components = params
            assert components == 3, f"Body {body_id} has {components} components"

    @pytest.mark.unit
    def test_degree_range(self):
        """Polynomial degree should be in range [7, 17]."""
        for body_id, params in BODY_PARAMS.items():
            _, degree, _, _ = params
            assert 7 <= degree <= 17, (
                f"Body {body_id} has degree {degree}, expected 7-17"
            )


class TestSegmentByteSize:
    """Test segment_byte_size calculation."""

    @pytest.mark.unit
    def test_degree_13_components_3(self):
        """Degree 13 with 3 components: 3 * 14 * 8 = 336 bytes."""
        assert segment_byte_size(13, 3) == 336

    @pytest.mark.unit
    def test_degree_15_components_3(self):
        """Degree 15 with 3 components: 3 * 16 * 8 = 384 bytes."""
        assert segment_byte_size(15, 3) == 384

    @pytest.mark.unit
    def test_degree_11_components_3(self):
        """Degree 11 with 3 components: 3 * 12 * 8 = 288 bytes."""
        assert segment_byte_size(11, 3) == 288

    @pytest.mark.unit
    def test_degree_9_components_3(self):
        """Degree 9 with 3 components: 3 * 10 * 8 = 240 bytes."""
        assert segment_byte_size(9, 3) == 240

    @pytest.mark.unit
    def test_nutation_segment(self):
        """Nutation: degree 16, 2 components: 2 * 17 * 8 = 272 bytes."""
        assert segment_byte_size(16, 2) == 272


class TestMagicAndVersion:
    """Test format constants."""

    @pytest.mark.unit
    def test_magic_is_4_bytes(self):
        assert len(MAGIC) == 4
        assert MAGIC == b"LEB1"

    @pytest.mark.unit
    def test_version_is_1(self):
        assert VERSION == 1


class _FakeMM:
    """Minimal mmap stand-in for madvise helper tests."""

    def __init__(self, length, raise_on_madvise=False):
        self._len = length
        self.raise_on_madvise = raise_on_madvise
        self.calls = []

    def __len__(self):
        return self._len

    def madvise(self, *args):
        self.calls.append(args)
        if self.raise_on_madvise:
            raise OSError("madvise failed")


class TestCompressedFormatRoundTrips:
    """LEB2 write-side packers (used by the generator, not the reader)."""

    def test_compressed_body_entry_round_trip(self):
        from libephemeris.leb_format import (
            COMPRESSED_BODY_ENTRY_SIZE,
            CompressedBodyEntry,
            read_compressed_body_entry,
            write_compressed_body_entry,
        )

        entry = CompressedBodyEntry(
            body_id=7,
            coord_type=1,
            segment_count=42,
            jd_start=2451545.0,
            jd_end=2469545.0,
            interval_days=32.0,
            degree=13,
            components=3,
            data_offset=4096,
            compressed_size=1234,
            uncompressed_size=5678,
        )
        buf = bytearray(COMPRESSED_BODY_ENTRY_SIZE)
        write_compressed_body_entry(buf, 0, entry)
        assert read_compressed_body_entry(buf, 0) == entry

    def test_chunk_index_header_round_trip(self):
        from libephemeris.leb_format import (
            CHUNK_INDEX_HEADER_SIZE,
            read_chunk_index_header,
            write_chunk_index_header,
        )

        buf = bytearray(CHUNK_INDEX_HEADER_SIZE)
        write_chunk_index_header(buf, 0, chunk_count=5, chunk_interval_days=3650.0)
        count, interval = read_chunk_index_header(buf, 0)
        assert count == 5
        assert interval == 3650.0

    def test_chunk_entry_round_trip(self):
        from libephemeris.leb_format import (
            CHUNK_ENTRY_SIZE,
            ChunkEntry,
            read_chunk_entry,
            write_chunk_entry,
        )

        entry = ChunkEntry(
            jd_start=2451545.0,
            jd_end=2455197.0,
            blob_offset=8192,
            compressed_size=999,
            uncompressed_size=4096,
            segment_start=10,
            segment_count=20,
        )
        buf = bytearray(CHUNK_ENTRY_SIZE)
        write_chunk_entry(buf, 0, entry)
        assert read_chunk_entry(buf, 0) == entry


class TestMadviseHelpers:
    def test_madvise_ranges_skips_nonpositive_length(self):
        from libephemeris.leb_format import _PAGE_SIZE, _madvise_ranges

        mm = _FakeMM(10)
        # offset past the (tiny) mapping -> page-aligned start exceeds len ->
        # computed length <= 0 -> the range is skipped, madvise never called.
        _madvise_ranges(mm, [(2 * _PAGE_SIZE, 8)])
        assert mm.calls == []

    def test_madvise_ranges_swallows_errors(self):
        from libephemeris.leb_format import _madvise_ranges

        mm = _FakeMM(1_000_000, raise_on_madvise=True)
        # Should not raise even though madvise() throws.
        _madvise_ranges(mm, [(0, 100)])
        assert mm.calls  # madvise was attempted

    def test_madvise_dontneed_none_is_noop(self):
        from libephemeris.leb_format import _madvise_dontneed

        _madvise_dontneed(None)  # no error

    def test_madvise_dontneed_success_and_error(self):
        from libephemeris.leb_format import _madvise_dontneed

        ok = _FakeMM(100)
        _madvise_dontneed(ok)
        assert ok.calls

        bad = _FakeMM(100, raise_on_madvise=True)
        _madvise_dontneed(bad)  # error swallowed
