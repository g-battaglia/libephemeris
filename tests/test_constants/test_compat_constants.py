"""Literal contracts for compatibility filenames, paths, and tidal constants."""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris import constants

pytestmark = pytest.mark.unit

FILENAME_CONSTANTS = {
    "FNAME_DE200": "de200.eph",
    "FNAME_DE403": "de403.eph",
    "FNAME_DE404": "de404.eph",
    "FNAME_DE405": "de405.eph",
    "FNAME_DE406": "de406.eph",
    "FNAME_DE431": "de431.eph",
    "FNAME_DFT": "de431.eph",
    "FNAME_DFT2": "de406.eph",
    "SE_FNAME_DE431": "de431.eph",
}

TIDAL_CONSTANTS = {
    "TIDAL_DE200": -23.8946,
    "TIDAL_DE403": -25.580,
    "TIDAL_DE404": -25.580,
    "TIDAL_DE405": -25.826,
    "TIDAL_DE406": -25.826,
    "TIDAL_DE421": -25.85,
    "TIDAL_DE422": -25.85,
    "TIDAL_DE430": -25.82,
    "TIDAL_DE431": -25.80,
    "TIDAL_DE440": -25.936,
    "TIDAL_DE441": -25.936,
    "TIDAL_WILLIAMS_BOGGS_2016": -25.97,
    "TIDAL_DEFAULT": -25.80,
    "TIDAL_AUTOMATIC": 999999,
    "TIDAL_26": -26.0,
    "TIDAL_JPLEPH": -25.8,
    "TIDAL_MOSEPH": -25.58,
    "TIDAL_STEPHENSON_2016": -25.85,
    "TIDAL_SWIEPH": -25.8,
}


@pytest.mark.parametrize(("name", "expected"), sorted(FILENAME_CONSTANTS.items()))
def test_filename_constant_literal(name: str, expected: str) -> None:
    """Every public compatibility filename retains its exact text."""
    value = getattr(constants, name)
    assert type(value) is str
    assert value == expected
    assert getattr(ephem, name) == value
    assert name in constants.__all__
    assert name in ephem.__all__


def test_ephemeris_search_path_literal() -> None:
    """The compatibility search-path constant retains its exact punctuation."""
    assert constants.EPHE_PATH == ".:/users/ephe2/:/users/ephe/"
    assert ephem.EPHE_PATH == constants.EPHE_PATH
    assert "EPHE_PATH" in constants.__all__
    assert "EPHE_PATH" in ephem.__all__


@pytest.mark.parametrize(("name", "expected"), sorted(TIDAL_CONSTANTS.items()))
def test_tidal_constant_literal(name: str, expected: float | int) -> None:
    """Every named tidal constant retains its exact native value and type."""
    value = getattr(constants, name)
    assert type(value) is type(expected)
    assert value == expected
    assert getattr(ephem, name) == value
    assert name in constants.__all__
    assert name in ephem.__all__


def test_named_alias_relationships_are_exact() -> None:
    """Aliases remain identities rather than independently drifting values."""
    assert constants.FNAME_DFT == constants.FNAME_DE431
    assert constants.FNAME_DFT2 == constants.FNAME_DE406
    assert constants.SE_FNAME_DE431 == constants.FNAME_DE431
    assert constants.TIDAL_DEFAULT == constants.TIDAL_DE431
    assert constants.TIDAL_DE403 == constants.TIDAL_DE404
    assert constants.TIDAL_DE405 == constants.TIDAL_DE406
    assert constants.TIDAL_DE421 == constants.TIDAL_DE422
    assert constants.TIDAL_DE440 == constants.TIDAL_DE441
