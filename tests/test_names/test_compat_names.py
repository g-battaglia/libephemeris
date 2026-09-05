# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Tests for the display-string tables in ``libephemeris.compat_names``.

Each table is compared entry by entry with the public function that answers
from it; the lookup rules the entry points keep for themselves (selector
folding, the low-byte mask, the unknown-selector default) are exercised over
their whole domain; and every table literal is parsed from the source so a
duplicated key, which a dict literal would silently collapse, is caught.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

import libephemeris as le
from libephemeris import compat_names, constants, contrib
from libephemeris.compat_names import (
    AYANAMSHA_NAMES,
    BODY_DISPLAY_NAMES,
    CONTRIB_HOUSE_SYSTEM_NAMES,
    HOUSE_SYSTEM_NAMES,
)

TABLES: dict[str, dict] = {
    "HOUSE_SYSTEM_NAMES": HOUSE_SYSTEM_NAMES,
    "AYANAMSHA_NAMES": AYANAMSHA_NAMES,
    "BODY_DISPLAY_NAMES": BODY_DISPLAY_NAMES,
    "CONTRIB_HOUSE_SYSTEM_NAMES": CONTRIB_HOUSE_SYSTEM_NAMES,
}

EXPECTED_SIZES = {
    "HOUSE_SYSTEM_NAMES": 26,
    "AYANAMSHA_NAMES": 47,
    "BODY_DISPLAY_NAMES": 42,
    "CONTRIB_HOUSE_SYSTEM_NAMES": 25,
}

#: The six SIDBIT_* projection flags, each alone and all together, plus the
#: bare mode as the control.
SIDBIT_WORDS = (
    0,
    256,
    512,
    1024,
    2048,
    4096,
    8192,
    256 | 512 | 1024 | 2048 | 4096 | 8192,
)


def _fold_house_selector(char: str) -> str:
    """The case rule ``house_name`` applies before the table lookup.

    ASCII lowercase letters select like their uppercase form, except ``'i'``
    (a system of its own); ``'g'`` is folded to Gauquelin by ``house_name``
    itself.
    """
    if char == "i":
        return char
    if "a" <= char <= "z":
        return char.upper()
    return char


# ---------------------------------------------------------------------------
# table shape
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("name", sorted(TABLES))
def test_table_size(name: str) -> None:
    """Each table carries exactly the entries the public function answers."""
    assert len(TABLES[name]) == EXPECTED_SIZES[name]


@pytest.mark.unit
@pytest.mark.parametrize("name", sorted(TABLES))
def test_table_values_are_plain_non_empty_strings(name: str) -> None:
    """Every label is a native ``str`` with content."""
    for key, value in TABLES[name].items():
        assert type(value) is str, (name, key, value)
        assert value, (name, key)


@pytest.mark.unit
@pytest.mark.parametrize("name", ["HOUSE_SYSTEM_NAMES", "CONTRIB_HOUSE_SYSTEM_NAMES"])
def test_house_tables_are_keyed_by_single_characters(name: str) -> None:
    """House tables are keyed by the one-character selector."""
    for key in TABLES[name]:
        assert type(key) is str and len(key) == 1, (name, key)


@pytest.mark.unit
@pytest.mark.parametrize("name", ["AYANAMSHA_NAMES", "BODY_DISPLAY_NAMES"])
def test_identifier_tables_are_keyed_by_ints(name: str) -> None:
    """Identifier tables are keyed by the public integer identifiers."""
    for key in TABLES[name]:
        assert type(key) is int, (name, key)


@pytest.mark.unit
def test_ayanamsha_table_covers_the_predefined_modes_exactly() -> None:
    """The predefined modes are 0-46 and nothing else."""
    assert set(AYANAMSHA_NAMES) == set(range(47))
    assert constants.SIDM_USER not in AYANAMSHA_NAMES


@pytest.mark.unit
def test_all_lists_the_four_tables() -> None:
    """The module exports exactly its four tables."""
    assert sorted(compat_names.__all__) == sorted(TABLES)


# ---------------------------------------------------------------------------
# no duplicate keys in the source literals
# ---------------------------------------------------------------------------


def _source_keys(table_name: str) -> list[object]:
    """Keys of the ``table_name`` dict literal, read from the module source.

    A duplicated key in a dict literal keeps only the last value at import
    time, so the check has to look at the source rather than at the dict.
    Constant keys are resolved through ``libephemeris.constants``.
    """
    source = Path(compat_names.__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    for node in tree.body:
        if not isinstance(node, ast.AnnAssign):
            continue
        if not (isinstance(node.target, ast.Name) and node.target.id == table_name):
            continue
        assert isinstance(node.value, ast.Dict), table_name
        keys: list[object] = []
        for key in node.value.keys:
            if isinstance(key, ast.Constant):
                keys.append(key.value)
            elif isinstance(key, ast.Name):
                keys.append(getattr(constants, key.id))
            else:  # pragma: no cover - the table literal is not that shape
                raise AssertionError(f"{table_name}: unexpected key node {key!r}")
        return keys
    raise AssertionError(f"{table_name}: dict literal not found in the source")


@pytest.mark.unit
@pytest.mark.parametrize("name", sorted(TABLES))
def test_source_literal_has_no_duplicate_keys(name: str) -> None:
    """No key appears twice in the source literal of a table."""
    keys = _source_keys(name)
    duplicates = sorted({repr(k) for k in keys if keys.count(k) > 1})
    assert not duplicates, f"{name}: duplicated keys {duplicates}"
    assert len(keys) == len(TABLES[name])
    assert list(TABLES[name]) == keys


# ---------------------------------------------------------------------------
# house_name
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("selector,label", sorted(HOUSE_SYSTEM_NAMES.items()))
def test_house_name_matches_table_entry(selector: str, label: str) -> None:
    """``house_name`` answers the table entry for int, bytes and str forms."""
    assert le.house_name(ord(selector)) == label
    assert le.house_name(selector.encode("latin-1")) == label
    assert le.house_name(selector) == label


@pytest.mark.unit
def test_house_name_whole_byte_domain_is_the_folded_table() -> None:
    """Over 0-255, ``house_name`` is the table after the case fold.

    ``'g'`` is folded to ``'G'`` by ``house_name`` (the tuple shape of
    ``houses()`` does not fold it); every other selector outside the table
    answers the empty string.
    """
    for code in range(256):
        char = _fold_house_selector(chr(code))
        if char == "g":
            char = "G"
        expected = HOUSE_SYSTEM_NAMES.get(char, "")
        assert le.house_name(code) == expected, code
        assert le.house_name(bytes([code])) == expected, code
        assert le.house_name(chr(code)) == expected, code


@pytest.mark.unit
def test_house_name_gauquelin_lowercase_folds_to_uppercase() -> None:
    """``'g'`` names Gauquelin sectors although it is not a selector alias."""
    assert le.house_name(ord("g")) == HOUSE_SYSTEM_NAMES["G"]
    assert "g" not in HOUSE_SYSTEM_NAMES


@pytest.mark.unit
@pytest.mark.parametrize("selector", [-1, -2, -100, 0x110000, 2**40])
def test_house_name_int_outside_chr_range_is_unknown(selector: int) -> None:
    """An int ``chr()`` rejects folds to the unknown-selector default."""
    assert le.house_name(selector) == ""


@pytest.mark.unit
@pytest.mark.parametrize("selector", [b"", b"AB", b"\xff\xfe", "", "AB"])
def test_house_name_multi_character_selector_is_refused(selector) -> None:
    """A selector that is not one character is refused with a TypeError."""
    with pytest.raises(TypeError, match="single character"):
        le.house_name(selector)


# ---------------------------------------------------------------------------
# get_ayanamsa_name
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("sidmode,label", sorted(AYANAMSHA_NAMES.items()))
def test_get_ayanamsa_name_matches_table_entry(sidmode: int, label: str) -> None:
    """``get_ayanamsa_name`` answers the table entry for each predefined mode."""
    assert le.get_ayanamsa_name(sidmode) == label


@pytest.mark.unit
def test_get_ayanamsa_name_low_byte_domain_is_the_table() -> None:
    """Over ids 0-300 under every SIDBIT word, the answer is the masked table."""
    for sidmode in range(301):
        expected = AYANAMSHA_NAMES.get(sidmode & 0xFF, "")
        for bits in SIDBIT_WORDS:
            assert le.get_ayanamsa_name(sidmode | bits) == expected, (sidmode, bits)


@pytest.mark.unit
def test_get_ayanamsa_name_user_mode_and_unassigned_block_are_empty() -> None:
    """``SIDM_USER`` and the unassigned ids above 46 have no label."""
    assert le.get_ayanamsa_name(constants.SIDM_USER) == ""
    for sidmode in range(47, 255):
        assert le.get_ayanamsa_name(sidmode) == "", sidmode


# ---------------------------------------------------------------------------
# get_planet_name
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("body_id,label", sorted(BODY_DISPLAY_NAMES.items()))
def test_get_planet_name_matches_table_entry(body_id: int, label: str) -> None:
    """``get_planet_name`` answers the table entry for each listed body."""
    assert le.get_planet_name(body_id) == label


@pytest.mark.unit
def test_get_planet_name_below_the_asteroid_block_is_the_table() -> None:
    """Below the moon block, every identifier not in the table is unnamed."""
    for body_id in range(-100, 9000):
        expected = BODY_DISPLAY_NAMES.get(body_id, "")
        assert le.get_planet_name(body_id) == expected, body_id


@pytest.mark.unit
def test_body_table_lists_planets_lunar_points_and_fictitious_bodies() -> None:
    """The table spans the planets, the lunar points and the fictitious ids."""
    assert set(range(0, 23)) <= set(BODY_DISPLAY_NAMES)
    assert set(range(40, 59)) <= set(BODY_DISPLAY_NAMES)


# ---------------------------------------------------------------------------
# contrib.house_system_name
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("code,label", sorted(CONTRIB_HOUSE_SYSTEM_NAMES.items()))
def test_house_system_name_matches_table_entry(code: str, label: str) -> None:
    """``contrib.house_system_name`` answers the table entry for each code."""
    assert contrib.house_system_name(code) == label


@pytest.mark.unit
def test_house_system_name_ascii_domain_is_the_table_without_folding() -> None:
    """Over the ASCII range, the answer is the table as given, else Unknown."""
    for code in range(128):
        char = chr(code)
        expected = CONTRIB_HOUSE_SYSTEM_NAMES.get(char, "Unknown")
        assert contrib.house_system_name(char) == expected, code


@pytest.mark.unit
@pytest.mark.parametrize("code", ["", "PK", "p", "j", "J", "z"])
def test_house_system_name_unlisted_code_is_unknown(code: str) -> None:
    """A str the table does not list answers ``"Unknown"``."""
    assert contrib.house_system_name(code) == "Unknown"


@pytest.mark.unit
@pytest.mark.parametrize("code", [ord("P"), b"P", None, 1.0, ["P"]])
def test_house_system_name_rejects_non_str(code: object) -> None:
    """A non-str argument raises ``TypeError`` naming the offending type."""
    with pytest.raises(TypeError, match="expects a str code, not"):
        contrib.house_system_name(code)  # type: ignore[arg-type]
