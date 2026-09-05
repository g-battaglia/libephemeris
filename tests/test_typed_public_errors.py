# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for issue #60: no untyped exception escapes a public API.

Five surfaces let a standard-library exception out where the typed hierarchy
in ``libephemeris.exceptions`` is the documented contract, or answered one
input differently from an equivalent one:

* ``lun_occult_where()`` raised ``KeyError`` from an internal reader lookup
  for a negative body id. The reference answers a negative id exactly as it
  answers the Sun, so the entry points now do the same; an unknown body
  raises ``IllegalBodyError``, the documented type that is both an
  ``UnknownBodyError`` and a ``ValueError``.
* Every entry point taking a house-system selector raised ``ValueError``
  from ``chr()`` for an int outside the character range while an equally
  unmapped int such as 999 fell through to the default. The reference
  answers an unmapped selector with the default (Placidus cusps; the
  cusp-interpolated position in ``house_pos``; ``""`` in ``house_name``).
* ``parse_orbital_elements("")`` raised ``IsADirectoryError`` because
  ``Path("")`` is ``.``, which exists. Special files such as ``/dev/null``
  must keep parsing.
* ``date_conversion()``, ``get_saros_number()`` and ``get_inex_number()``
  normalised with ``.lower()`` before validating, so a non-string leaked
  ``AttributeError``.
* ``rise_trans()`` named a private class in the traceback.
"""

from __future__ import annotations

import os

import pytest

import libephemeris as eph
from libephemeris.exceptions import IllegalBodyError, UnknownBodyError
from libephemeris.time_utils import date_conversion as date_conversion_impl

_JD = 2451545.0
_GEOPOS = [12.0, 42.0, 0.0]


# ---------------------------------------------------------------------------
# lun_occult_* : a negative id is answered as the Sun; an unknown body raises
# the documented contract type
# ---------------------------------------------------------------------------


_OCCULT_ENTRY_POINTS = [
    ("where", lambda body: eph.lun_occult_where(_JD, body)),
    ("when_glob", lambda body: eph.lun_occult_when_glob(_JD, body)),
    ("when_loc", lambda body: eph.lun_occult_when_loc(_JD, body, _GEOPOS)),
]


@pytest.fixture(scope="module")
def sun_answers():
    """The Sun's answer from each entry point, computed once."""
    return {label: call(eph.SUN) for label, call in _OCCULT_ENTRY_POINTS}


@pytest.mark.unit
@pytest.mark.parametrize(
    ("label", "call"),
    _OCCULT_ENTRY_POINTS,
    ids=[name for name, _ in _OCCULT_ENTRY_POINTS],
)
@pytest.mark.parametrize("body", [-1, -2, -1000])
def test_negative_body_is_answered_as_the_sun(label, call, body, sun_answers):
    """ECL_NUT (-1) is a pseudo-target; the reference answers it as body 0."""
    assert call(body) == sun_answers[label]


@pytest.mark.unit
def test_negative_body_no_longer_raises_a_bare_key_error():
    """The reported symptom: the call answers instead of leaking a KeyError."""
    retflag, geopos, attr = eph.lun_occult_where(_JD, -1)
    assert isinstance(retflag, int)
    assert len(geopos) >= 2
    assert len(attr) >= 1


@pytest.mark.unit
@pytest.mark.parametrize(
    ("label", "call"),
    _OCCULT_ENTRY_POINTS,
    ids=[name for name, _ in _OCCULT_ENTRY_POINTS],
)
def test_an_unknown_body_raises_the_documented_contract(label, call):
    """The docstrings promise both UnknownBodyError and ValueError."""
    with pytest.raises(IllegalBodyError) as excinfo:
        call(999)
    exc = excinfo.value
    assert isinstance(exc, UnknownBodyError)
    assert isinstance(exc, ValueError)
    assert exc.body_id == 999
    assert isinstance(exc.__cause__, UnknownBodyError)


@pytest.mark.unit
def test_a_real_occultation_target_still_works():
    """The fold must not touch a body that can actually be occulted."""
    retflag, geopos, attr = eph.lun_occult_where(_JD, eph.VENUS)
    assert isinstance(retflag, int)
    assert len(geopos) >= 2
    assert len(attr) >= 1


# ---------------------------------------------------------------------------
# house-system selectors : every unmapped selector behaves the same way
# ---------------------------------------------------------------------------


# Entry points where the selector picks the system whose cusps are returned.
_PRIMARY_CUSP_ENTRY_POINTS = [
    ("houses", lambda h: eph.houses(_JD, 51.5, -0.12, h)),
    ("houses_ex", lambda h: eph.houses_ex(_JD, 51.5, -0.12, h, 0)),
    ("houses_ex2", lambda h: eph.houses_ex2(_JD, 51.5, -0.12, h, 0)),
    ("houses_armc", lambda h: eph.houses_armc(100.0, 41.9, 23.4, h)),
    ("houses_armc_ex2", lambda h: eph.houses_armc_ex2(100.0, 41.9, 23.4, h)),
    (
        "houses_with_fallback",
        lambda h: eph.houses_with_fallback(_JD, 51.5, -0.12, h),
    ),
    (
        "houses_armc_with_fallback",
        lambda h: eph.houses_armc_with_fallback(100.0, 41.9, 23.4, h),
    ),
]

# The fallback selector is converted up front too; at this latitude the
# fallback is never used, so the answer stays the primary (Placidus) one.
_FALLBACK_CUSP_ENTRY_POINTS = [
    (
        "houses_with_fallback/fallback_hsys",
        lambda h: eph.houses_with_fallback(_JD, 51.5, -0.12, ord("P"), h),
    ),
    (
        "houses_armc_with_fallback/fallback_hsys",
        lambda h: eph.houses_armc_with_fallback(100.0, 41.9, 23.4, ord("P"), h),
    ),
]

_CUSP_ENTRY_POINTS = _PRIMARY_CUSP_ENTRY_POINTS + _FALLBACK_CUSP_ENTRY_POINTS

_HOUSE_POS_ENTRY_POINTS = [
    ("house_pos/5-arg", lambda h: eph.house_pos(100.0, 41.9, 23.4, (10.0, 0.0), h)),
    ("house_pos/6-arg", lambda h: eph.house_pos(100.0, 41.9, 23.4, h, 10.0, 0.0)),
]

_SELECTOR_ENTRY_POINTS = _CUSP_ENTRY_POINTS + _HOUSE_POS_ENTRY_POINTS

# -1 and 0x110000 sit outside chr()'s range; 999 is inside it but unmapped.
_UNMAPPED_SELECTORS = [-1, 999, 0x110000]


def _ids(pairs):
    return [label for label, _ in pairs]


@pytest.mark.unit
@pytest.mark.parametrize(
    ("label", "call"), _SELECTOR_ENTRY_POINTS, ids=_ids(_SELECTOR_ENTRY_POINTS)
)
@pytest.mark.parametrize("hsys", _UNMAPPED_SELECTORS)
def test_every_entry_point_treats_an_unmapped_int_as_an_unknown_selector(
    label, call, hsys
):
    """An int chr() rejects answers exactly like an unmapped byte selector."""
    assert call(hsys) == call(b"?")


@pytest.mark.unit
@pytest.mark.parametrize(
    ("label", "call"), _CUSP_ENTRY_POINTS, ids=_ids(_CUSP_ENTRY_POINTS)
)
@pytest.mark.parametrize("hsys", _UNMAPPED_SELECTORS)
def test_cusp_entry_points_fall_back_to_placidus(label, call, hsys):
    """The default for an unknown selector is Placidus, as in the reference."""
    assert call(hsys) == call(b"P")


@pytest.mark.unit
@pytest.mark.parametrize(
    ("label", "call"), _SELECTOR_ENTRY_POINTS, ids=_ids(_SELECTOR_ENTRY_POINTS)
)
def test_a_valid_selector_is_accepted_in_every_form(label, call):
    """An int code, a byte, and a str select the same system."""
    assert call(ord("K")) == call(b"K") == call("K")


@pytest.mark.unit
@pytest.mark.parametrize(
    ("label", "call"),
    _PRIMARY_CUSP_ENTRY_POINTS,
    ids=_ids(_PRIMARY_CUSP_ENTRY_POINTS),
)
def test_a_valid_selector_still_selects_its_own_system(label, call):
    """Koch cusps differ from Placidus cusps at this latitude."""
    assert call(ord("K")) != call(b"P")


@pytest.mark.unit
@pytest.mark.parametrize("hsys", [-1, -2, 999, 0, 0x110000])
def test_house_name_returns_empty_for_any_unmapped_selector(hsys):
    """-1 used to raise while 999 returned "" — the same kind of input."""
    assert eph.house_name(hsys) == ""


@pytest.mark.unit
def test_house_name_still_resolves_known_selectors():
    assert eph.house_name(ord("P")) == "Placidus"
    assert eph.house_name(b"P") == "Placidus"
    assert eph.house_name("P") == "Placidus"
    assert eph.house_name(ord("K")) == "Koch"


# ---------------------------------------------------------------------------
# parse_orbital_elements : a directory is not a file, a special file is
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("path", ["", ".", "..", "/"])
def test_non_file_paths_raise_the_documented_error(path):
    """Path("") is ".", which exists, so the guard was simply bypassed."""
    with pytest.raises(FileNotFoundError):
        eph.parse_orbital_elements(path)


@pytest.mark.unit
def test_a_directory_argument_raises_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        eph.parse_orbital_elements(str(tmp_path))


@pytest.mark.unit
def test_missing_file_still_raises_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        eph.parse_orbital_elements(str(tmp_path / "absent.txt"))


@pytest.mark.unit
@pytest.mark.skipif(not os.path.exists("/dev/null"), reason="no /dev/null")
def test_a_special_file_is_still_read():
    """/dev/null exists and is not a directory: it parses to an empty list."""
    assert eph.parse_orbital_elements("/dev/null") == []


@pytest.mark.unit
def test_an_empty_regular_file_parses_to_an_empty_list(tmp_path):
    path = tmp_path / "empty.txt"
    path.write_text("", encoding="utf-8")
    assert eph.parse_orbital_elements(str(path)) == []


@pytest.mark.unit
def test_a_real_file_still_parses(tmp_path):
    path = tmp_path / "orbits.txt"
    path.write_text(
        "# comment\n2451545.0,2451545.0,0.0,10.0,0.1,0.0,0.0,0.0,MyPlanet\n",
        encoding="utf-8",
    )
    elements = eph.parse_orbital_elements(str(path))
    assert [element.name for element in elements] == ["MyPlanet"]


# ---------------------------------------------------------------------------
# String arguments validated before they are normalised
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("calendar", [123, None, 1.5, ["g"], object()])
def test_date_conversion_rejects_a_non_string_calendar(calendar):
    """.lower() ran before the check, so a non-string leaked AttributeError."""
    with pytest.raises(ValueError):
        eph.date_conversion(2000, 1, 1, 0.0, calendar)
    with pytest.raises(ValueError):
        date_conversion_impl(2000, 1, 1, 0.0, calendar)


@pytest.mark.unit
@pytest.mark.parametrize("eclipse_type", [123, None, 1.5, ["solar"], object()])
def test_saros_and_inex_reject_a_non_string_type(eclipse_type):
    with pytest.raises(ValueError):
        eph.get_saros_number(_JD, eclipse_type)
    with pytest.raises(ValueError):
        eph.get_inex_number(_JD, eclipse_type)


@pytest.mark.unit
def test_valid_string_arguments_are_unaffected():
    assert eph.date_conversion(2000, 1, 1, 0.0, "g")[0] is True
    assert eph.date_conversion(2000, 1, 1, 0.0, b"g")[0] is True
    assert date_conversion_impl(1582, 10, 15, 12.0, "j") == (1582, 10, 5, 12.0)
    assert isinstance(eph.get_saros_number(_JD, "solar"), int)
    assert isinstance(eph.get_inex_number(_JD, "lunar"), int)


@pytest.mark.unit
@pytest.mark.parametrize("bad", ["x", "gregorian", ""])
def test_an_unknown_calendar_string_still_raises_value_error(bad):
    with pytest.raises(ValueError):
        eph.date_conversion(2000, 1, 1, 0.0, bad)


# ---------------------------------------------------------------------------
# rise_trans : the traceback names a public class
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_illegal_rise_body_error_is_public():
    """A user reading a traceback must not meet a private class name."""
    with pytest.raises(IllegalBodyError) as excinfo:
        eph.rise_trans(_JD, -5, eph.CALC_RISE, _GEOPOS)
    exc = excinfo.value
    assert not type(exc).__name__.startswith("_")
    assert type(exc).__module__ == "libephemeris.exceptions"
    assert exc.body_id == -5


@pytest.mark.unit
def test_illegal_body_error_keeps_both_contracts():
    """Existing callers catch UnknownBodyError or ValueError; both still work."""
    assert issubclass(IllegalBodyError, UnknownBodyError)
    assert issubclass(IllegalBodyError, ValueError)
    assert eph.IllegalBodyError is IllegalBodyError
    with pytest.raises(UnknownBodyError):
        eph.rise_trans(_JD, -5, eph.CALC_RISE, _GEOPOS)
    with pytest.raises(ValueError):
        eph.rise_trans_true_hor(_JD, -5, eph.CALC_RISE, _GEOPOS, 0.0, 0.0, 0.0)
