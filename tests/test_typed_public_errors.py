# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for issue #60: no untyped exception escapes a public API.

Five entry points let a standard-library exception out where the typed
hierarchy in ``libephemeris.exceptions`` is the documented contract, or where
a neighbouring input already had a defined answer:

* ``lun_occult_where()`` raised ``KeyError`` from an internal reader lookup.
* ``house_name()`` raised ``ValueError`` from ``chr()`` for a negative
  selector while returning ``""`` for every other unmapped one.
* ``parse_orbital_elements("")`` raised ``IsADirectoryError`` because
  ``Path("")`` is ``.``, which exists.
* ``date_conversion()``, ``get_saros_number()`` and ``get_inex_number()``
  normalised with ``.lower()`` before validating, so a non-string leaked
  ``AttributeError``.
* ``rise_trans()`` named a private class in the traceback.
"""

from __future__ import annotations

import pytest

import libephemeris as eph
from libephemeris.exceptions import Error, IllegalBodyError, UnknownBodyError
from libephemeris.time_utils import date_conversion as date_conversion_impl

_JD = 2451545.0
_GEOPOS = [12.0, 42.0, 0.0]


# ---------------------------------------------------------------------------
# lun_occult_* : a pseudo-target is an unknown body, not a KeyError
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize(
    ("label", "call"),
    [
        ("where", lambda body: eph.lun_occult_where(_JD, body)),
        ("when_glob", lambda body: eph.lun_occult_when_glob(_JD, body)),
        ("when_loc", lambda body: eph.lun_occult_when_loc(_JD, body, _GEOPOS)),
    ],
    ids=lambda v: v if isinstance(v, str) else "",
)
@pytest.mark.parametrize("body", [-1, -2, -10])
def test_negative_body_is_an_unknown_body(label, call, body):
    """ECL_NUT (-1) answers obliquity, not a position; it cannot be occulted."""
    with pytest.raises(UnknownBodyError):
        call(body)


@pytest.mark.unit
def test_negative_body_no_longer_raises_a_bare_key_error():
    """The reported symptom, stated as the type it must not be."""
    with pytest.raises(Error) as excinfo:
        eph.lun_occult_where(_JD, -1)
    assert not isinstance(excinfo.value, KeyError)


@pytest.mark.unit
def test_a_real_occultation_target_still_works():
    """The guard must not reject a body that can actually be occulted."""
    retflag, geopos, attr = eph.lun_occult_where(_JD, eph.VENUS)
    assert isinstance(retflag, int)
    assert len(geopos) >= 2
    assert len(attr) >= 1


# ---------------------------------------------------------------------------
# house_name / houses* : every unmapped selector behaves the same way
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("hsys", [-1, -2, 999, 0, 1114112])
def test_house_name_returns_empty_for_any_unmapped_selector(hsys):
    """-1 used to raise while 999 returned "" — the same kind of input."""
    assert eph.house_name(hsys) == ""


@pytest.mark.unit
def test_house_name_still_resolves_known_selectors():
    assert eph.house_name(ord("P")) == "Placidus"
    assert eph.house_name(b"P") == "Placidus"
    assert eph.house_name(ord("K")) == "Koch"


@pytest.mark.unit
@pytest.mark.parametrize("hsys", [-1, 999])
def test_house_calculations_treat_out_of_range_like_any_unknown_selector(hsys):
    """The three cusp entry points carried the identical chr() crash."""
    reference = eph.houses(_JD, 51.5, -0.12, 999)
    assert eph.houses(_JD, 51.5, -0.12, hsys) == reference

    reference_armc = eph.houses_armc(100.0, 41.9, 23.4, 999)
    assert eph.houses_armc(100.0, 41.9, 23.4, hsys) == reference_armc

    reference_ex = eph.houses_ex(_JD, 51.5, -0.12, 999, 0)
    assert eph.houses_ex(_JD, 51.5, -0.12, hsys, 0) == reference_ex


# ---------------------------------------------------------------------------
# parse_orbital_elements : a directory is not a file
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
