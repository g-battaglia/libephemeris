"""Native handles return the selected file's raw star records in stored units."""
from data_paths import CORE_LEB

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import struct

import libephemeris as e
from libephemeris import state
from libephemeris.exceptions import ConfigurationError, DataNotFoundError
import pytest


FIXTURE = CORE_LEB
FIELDS = ("star_id", "ra_j2000_deg", "dec_j2000_deg", "pm_ra_deg_per_yr",
          "pm_dec_deg_per_yr", "parallax_arcsec", "rv_km_per_s", "magnitude")
# Even an ID shared with the built-in catalog must return these selected words.
# The synthetic magnitude is deliberately unrelated to Spica's catalog value.
ROWS = ((e.SPICA_STAR, 123.25, -45.5, 2.0 ** -20, -(2.0 ** -19), 0.125, -17.5, 12.75),
        (-271828, 271.5, 12.25, -0.0, 2.0 ** -25, 0.0625, 31.25, -1.5))


def star_image(magic, rows=ROWS, trailing=b""):
    """Section 4 is raw `<iddddddd4s` in both formats, without a row header."""
    payloads = [(0, b"")]
    if rows is not None:
        stars = b"".join(struct.pack("<i7d4s", *row, b"\xa1\xb2\xc3\xd4") for row in rows)
        payloads.append((4, stars + trailing))
    version, flags = (1, 0) if magic == "LEB1" else (2, 1)
    header = struct.pack("<4sIIIdddI20x", magic.encode(), version, len(payloads),
                         0, 0.0, 16.0, 0.0, flags)
    offset = len(header) + 24 * len(payloads)
    directory = []
    for section, payload in payloads:
        directory.append(struct.pack("<IIQQ", section, 0, offset, len(payload)))
        offset += len(payload)
    return header + b"".join(directory) + b"".join(payload for _, payload in payloads)


def record_values(record):
    return tuple(getattr(record, field) for field in FIELDS)


def assert_stored_words(record, expected):
    actual = record_values(record)
    assert actual[0] == expected[0]
    assert struct.pack("<7d", *actual[1:]) == struct.pack("<7d", *expected[1:])


@pytest.fixture(params=("owned", "indexed"))
def backing_mode(request):
    return request.param


@pytest.fixture(params=("LEB1", "LEB2"))
def reader_factory(request, backing_mode, monkeypatch, tmp_path):
    monkeypatch.setenv("LIBEPHEMERIS_DATA_DIR", str(tmp_path))
    monkeypatch.setenv("LIBEPHEMERIS_CONFIG", str(tmp_path / "absent.toml"))
    monkeypatch.setenv("LIBEPHEMERIS_ENV_FILE", str(tmp_path / "absent.env"))
    monkeypatch.setenv("LIBEPHEMERIS_LEB_READER", backing_mode)
    e.close()
    e.reload_configuration()
    e.set_calc_mode("leb")
    e.set_precision_tier("base")
    readers = []

    def open_reader(**options):
        suffix = "leb" if request.param == "LEB1" else "leb2"
        path = tmp_path / f"stars-{len(readers)}.{suffix}"
        path.write_bytes(star_image(request.param, **options))
        e.set_leb_file(str(path))
        reader = state.get_leb_reader()
        assert reader is not None and reader.path == str(path)
        readers.append(reader)
        return reader

    yield open_reader
    for reader in readers:
        reader.close()
    e.close()
    monkeypatch.undo()
    e.reload_configuration()
    e.set_calc_mode("leb")
    e.set_precision_tier("base")
    e.set_leb_file(str(FIXTURE))


def test_selected_star_values_keep_stored_units_and_are_immutable(reader_factory):
    reader = reader_factory()
    assert not hasattr(e, "LebStarRecord")
    for row in ROWS:
        record = reader.get_star(row[0])
        assert isinstance(record, state.LebStarRecord)
        assert_stored_words(record, row)
        for field, value in zip(FIELDS, row):
            with pytest.raises(AttributeError):
                setattr(record, field, value)


@pytest.mark.parametrize("rows", (None, (), (ROWS[1],)))
def test_absent_catalog_empty_catalog_or_missing_id_does_not_fall_back(reader_factory, rows):
    reader = reader_factory(rows=rows)
    with pytest.raises(DataNotFoundError):
        reader.get_star(e.SPICA_STAR)


def test_last_complete_duplicate_wins_and_partial_tail_is_not_a_record(reader_factory):
    replacement = (*ROWS[0][:-1], 7.5)
    reader = reader_factory(rows=(ROWS[0], replacement),
                            trailing=struct.pack("<i", ROWS[1][0]) + b"\xff" * 59)
    assert_stored_words(reader.get_star(e.SPICA_STAR), replacement)
    with pytest.raises(DataNotFoundError):
        reader.get_star(ROWS[1][0])


def test_raw_ieee_record_is_not_physically_validated_or_normalized(reader_factory):
    words = (0x8000000000000000, 0x0000000000000001, 0x7FF8000000001234,
             0xFFF0000000000000, 0x7FF0000000000000, 0, 0xBFF8000000000000)
    values = struct.unpack("<7d", struct.pack("<7Q", *words))
    row = (123456789, *values)
    reader = reader_factory(rows=(row,))
    assert_stored_words(reader.get_star(row[0]), row)


def test_records_and_retained_handles_outlive_close_with_concurrent_reads(reader_factory):
    first = reader_factory()
    saved = first.get_star(e.SPICA_STAR)
    second = state.get_leb_reader()
    assert second is not None
    first.close()
    for star_id in (e.SPICA_STAR, 999999):
        with pytest.raises(ConfigurationError):
            first.get_star(star_id)
    replacement = reader_factory(rows=((*ROWS[0][:-1], 7.5),))
    e.close()
    try:
        with ThreadPoolExecutor(max_workers=4) as pool:
            records = list(pool.map(second.get_star, [row[0] for row in ROWS] * 12))
        for record, expected in zip(records, ROWS * 12):
            assert_stored_words(record, expected)
        assert replacement.get_star(e.SPICA_STAR).magnitude == 7.5
        second.close()
        replacement.close()
        assert_stored_words(saved, ROWS[0])
        for record, expected in zip(records, ROWS * 12):
            assert_stored_words(record, expected)
    finally:
        second.close()
