# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for issue #57: impossible type masks must stay refused.

Two eclipse geometries cannot occur: a partial eclipse is never central, and
a hybrid (annular-total) eclipse always has a central line. The guards that
refuse them compared the whole mask for equality, so setting one unrelated
bit carried the same impossible request straight past them and into the
search.

The guards now test the bits. The refusal is scoped to masks that name no
realisable alternative type, so a filter asking for totals *and* partials
keeps working -- a central total exists, and the impossible part simply
never matches.
"""

from __future__ import annotations

import pytest

from libephemeris.constants import (
    ECL_ANNULAR,
    ECL_ANNULAR_TOTAL,
    ECL_CENTRAL,
    ECL_NONCENTRAL,
    ECL_PARTIAL,
    ECL_TOTAL,
)
from libephemeris.eclipse import (
    _normalize_occultation_filter,
    _sol_glob_reject_impossible,
)
from libephemeris.exceptions import Error


@pytest.mark.unit
@pytest.mark.parametrize(
    "mask",
    [
        ECL_PARTIAL | ECL_CENTRAL,
        ECL_PARTIAL | ECL_CENTRAL | ECL_NONCENTRAL,
    ],
    ids=["pair", "pair-plus-noncentral"],
)
def test_central_partial_is_refused_whatever_else_is_set(mask):
    """The reported case: one extra bit used to smuggle the pair through."""
    with pytest.raises(Error, match="central partial eclipses do not exist"):
        _sol_glob_reject_impossible(mask)


@pytest.mark.unit
@pytest.mark.parametrize(
    "mask",
    [
        ECL_ANNULAR_TOTAL | ECL_NONCENTRAL,
        ECL_ANNULAR_TOTAL | ECL_NONCENTRAL | ECL_CENTRAL,
    ],
    ids=["pair", "pair-plus-central"],
)
def test_noncentral_hybrid_is_refused_whatever_else_is_set(mask):
    """The second impossible pair has the same shape and the same defect."""
    with pytest.raises(Error, match="non-central hybrid"):
        _sol_glob_reject_impossible(mask)


@pytest.mark.unit
@pytest.mark.parametrize(
    "mask",
    [
        0,
        ECL_PARTIAL,
        ECL_PARTIAL | ECL_NONCENTRAL,
        ECL_TOTAL,
        ECL_TOTAL | ECL_CENTRAL,
        ECL_ANNULAR | ECL_CENTRAL,
        ECL_ANNULAR_TOTAL | ECL_CENTRAL,
        # A realisable alternative type is named, so the mask describes
        # something even though it also contains an impossible pair.
        ECL_TOTAL | ECL_PARTIAL | ECL_CENTRAL,
        ECL_TOTAL | ECL_ANNULAR_TOTAL | ECL_NONCENTRAL,
        ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_CENTRAL | ECL_NONCENTRAL,
    ],
)
def test_realisable_masks_are_still_accepted(mask):
    """Narrowing the guard must not start refusing legitimate filters."""
    _sol_glob_reject_impossible(mask)


@pytest.mark.unit
@pytest.mark.parametrize(
    "ecltype",
    [
        ECL_PARTIAL | ECL_CENTRAL,
        ECL_PARTIAL | ECL_CENTRAL | ECL_NONCENTRAL,
    ],
    ids=["pair", "pair-plus-noncentral"],
)
def test_occultation_filter_refuses_central_partial(ecltype):
    """The occultation normaliser carried the identical defect."""
    with pytest.raises(Error, match="central partial eclipses do not exist"):
        _normalize_occultation_filter(ecltype, "Venus", is_sun=False)


@pytest.mark.unit
def test_occultation_filter_keeps_expanding_an_everything_request():
    """A filter naming several types is not an impossible request.

    ``ECL_TOTAL | ECL_PARTIAL | ECL_CENTRAL | ECL_NONCENTRAL`` contains the
    impossible pair but also asks for totals, which are central; refusing it
    would break the explicit spelling of "everything".
    """
    wanted = _normalize_occultation_filter(
        ECL_TOTAL | ECL_PARTIAL | ECL_CENTRAL | ECL_NONCENTRAL,
        "Venus",
        is_sun=False,
    )
    assert wanted & ECL_TOTAL
    assert wanted & ECL_PARTIAL
    assert wanted & ECL_CENTRAL
    assert wanted & ECL_NONCENTRAL


@pytest.mark.unit
def test_occultation_empty_filter_still_expands_to_everything():
    """The empty filter is expanded after the guard, so it is unaffected."""
    wanted = _normalize_occultation_filter(0, "Venus", is_sun=False)
    assert wanted == (ECL_TOTAL | ECL_PARTIAL | ECL_NONCENTRAL | ECL_CENTRAL)
