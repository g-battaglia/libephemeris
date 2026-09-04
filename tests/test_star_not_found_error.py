# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for issue #59: fixed-star lookups raise StarNotFoundError.

``StarNotFoundError`` was exported, documented and never raised. Every
fixed-star failure came back as the base ``Error``, so the ``except
StarNotFoundError`` clause the documentation shows could not run, and a
caller had to fall back to matching message text.

Every failure mode of both families now raises it, carrying the identifier
that failed and the search form that was attempted.
"""

from __future__ import annotations

import pytest

import libephemeris as eph
from libephemeris.exceptions import DataNotFoundError, Error, StarNotFoundError

_JD = 2451545.0

_UT_ENTRY_POINTS = [
    ("fixstar_ut", lambda name: eph.fixstar_ut(name, _JD)),
    ("fixstar2_ut", lambda name: eph.fixstar2_ut(name, _JD)),
]
_ET_ENTRY_POINTS = [
    ("fixstar", lambda name: eph.fixstar(name, _JD)),
    ("fixstar2", lambda name: eph.fixstar2(name, _JD)),
]
_MAG_ENTRY_POINTS = [
    ("fixstar_mag", eph.fixstar_mag),
    ("fixstar2_mag", eph.fixstar2_mag),
]
_ALL_POSITION_ENTRY_POINTS = _UT_ENTRY_POINTS + _ET_ENTRY_POINTS


@pytest.mark.unit
@pytest.mark.parametrize(
    ("label", "call"),
    _ALL_POSITION_ENTRY_POINTS,
    ids=lambda v: v if isinstance(v, str) else "",
)
def test_unknown_name_raises_the_typed_error(label, call):
    """The reported case, on every position entry point."""
    with pytest.raises(StarNotFoundError):
        call("NoSuchStar12345")


@pytest.mark.unit
@pytest.mark.parametrize(
    ("label", "call"), _MAG_ENTRY_POINTS, ids=lambda v: v if isinstance(v, str) else ""
)
def test_unknown_name_raises_the_typed_error_from_the_magnitude_helpers(label, call):
    """The magnitude helpers share the same resolvers and the same contract."""
    with pytest.raises(StarNotFoundError):
        call("NoSuchStar12345")


@pytest.mark.unit
@pytest.mark.parametrize(
    ("search", "expected_type"),
    [
        ("NoSuchStar12345", "name"),
        (",zzNoSuchNomenclature", "nomenclature"),
        ("99999", "sequential"),
        ("ZzzzNoSuchPrefix%", "wildcard"),
        ("", "empty"),
    ],
)
def test_every_v2_search_form_reports_its_kind(search, expected_type):
    """The v2 family covers all five documented search forms."""
    with pytest.raises(StarNotFoundError) as excinfo:
        eph.fixstar2_ut(search, _JD)
    assert excinfo.value.search_type == expected_type
    assert excinfo.value.star_id == search


@pytest.mark.unit
@pytest.mark.parametrize("search", ["NoSuchStar12345", "", "99999"])
def test_v1_family_reports_the_identifier_that_failed(search):
    """The v1 family carries the same attributes."""
    with pytest.raises(StarNotFoundError) as excinfo:
        eph.fixstar_ut(search, _JD)
    assert excinfo.value.star_id == search
    assert excinfo.value.search_type


@pytest.mark.unit
def test_the_documented_except_clause_now_runs():
    """The exact shape from the class docstring and the issue."""
    caught = None
    try:
        eph.fixstar_ut("NoSuchStar12345", _JD)
    except StarNotFoundError as exc:
        caught = exc
    except Error:  # pragma: no cover - the defect this test exists for
        pytest.fail("fixed-star lookup still raises the base Error")
    assert caught is not None
    assert caught.star_id == "NoSuchStar12345"


@pytest.mark.unit
def test_the_error_stays_inside_the_documented_hierarchy():
    """Existing `except Error` callers keep working unchanged."""
    assert issubclass(StarNotFoundError, DataNotFoundError)
    assert issubclass(StarNotFoundError, Error)
    with pytest.raises(Error):
        eph.fixstar_ut("NoSuchStar12345", _JD)
    with pytest.raises(DataNotFoundError):
        eph.fixstar_ut("NoSuchStar12345", _JD)


@pytest.mark.unit
def test_batch_lookup_raises_the_typed_error():
    """The batch helper resolves names itself and must agree."""
    with pytest.raises(StarNotFoundError):
        eph.batch_fixstars_ut(["Regulus", "NoSuchStar12345"], _JD)


@pytest.mark.unit
def test_batch_lookup_still_skips_when_asked():
    """skip_errors catches the subclass exactly as it caught the base."""
    results = eph.batch_fixstars_ut(
        ["Regulus", "NoSuchStar12345"], _JD, skip_errors=True
    )
    # The unresolved entry keeps its slot as None, as it always has.
    assert len(results) == 2
    assert results[0] is not None
    assert results[1] is None


@pytest.mark.unit
def test_known_stars_are_unaffected():
    """Narrowing the error type must not disturb a successful lookup."""
    position, name, _flags = eph.fixstar_ut("Regulus", _JD)
    assert name
    assert len(position) == 6
    magnitude, mag_name = eph.fixstar_mag("Regulus")
    assert isinstance(magnitude, float)
    assert mag_name
