"""The eclipse/occultation type-filter tables and the where/how result layout.

Pure functions only, no ephemeris is opened. The tables spell out what the
public searches rely on: which classified eclipses a ``sol_eclipse_when_glob``
mask accepts, how a shorthand occultation filter expands, how a conjunction
seek measures its longitude gaps, and the slot layout of the ``geopos`` and
``attr`` tuples the where/how functions return.
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
    ECL_VISIBLE,
)
from libephemeris.eclipse import (
    _LocalCircumstances,
    _longitude_gap_ahead,
    _longitude_gap_nearest,
    _normalize_occultation_filter,
    _shadow_center_geopos,
    _sol_glob_accepts,
    _sol_glob_reject_impossible,
)
from libephemeris.exceptions import Error

CLASSES = (ECL_TOTAL, ECL_ANNULAR, ECL_PARTIAL, ECL_ANNULAR_TOTAL)
CENTRALITIES = (ECL_CENTRAL, ECL_NONCENTRAL)
BOTH = ECL_CENTRAL | ECL_NONCENTRAL
#: Every type a classified eclipse can carry: one centrality, one class.
CLASSIFIED = tuple(g | c for g in CENTRALITIES for c in CLASSES)
FILTER_BITS = BOTH | ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL


def accepted_by(mask: int) -> set[int]:
    return {rf for rf in CLASSIFIED if _sol_glob_accepts(mask, rf)}


class TestSolarGlobFilter:
    def test_zero_accepts_everything(self):
        assert accepted_by(0) == set(CLASSIFIED)

    @pytest.mark.parametrize("mask", [ECL_CENTRAL, ECL_NONCENTRAL, BOTH])
    def test_centrality_alone_matches_nothing(self, mask):
        assert accepted_by(mask) == set()

    @pytest.mark.parametrize("cls", CLASSES)
    def test_single_class_takes_either_centrality(self, cls):
        assert accepted_by(cls) == {ECL_CENTRAL | cls, ECL_NONCENTRAL | cls}

    @pytest.mark.parametrize(
        "mask",
        [
            ECL_TOTAL | ECL_ANNULAR,
            ECL_TOTAL | ECL_PARTIAL,
            ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL,
        ],
    )
    def test_several_classes_without_centrality_match_nothing(self, mask):
        assert accepted_by(mask) == set()

    def test_class_with_centrality_needs_both(self):
        assert accepted_by(ECL_TOTAL | ECL_ANNULAR | ECL_CENTRAL) == {
            ECL_CENTRAL | ECL_TOTAL,
            ECL_CENTRAL | ECL_ANNULAR,
        }
        assert accepted_by(ECL_PARTIAL | ECL_NONCENTRAL) == {
            ECL_NONCENTRAL | ECL_PARTIAL
        }
        assert accepted_by(ECL_TOTAL | BOTH) == {
            ECL_CENTRAL | ECL_TOTAL,
            ECL_NONCENTRAL | ECL_TOTAL,
        }

    def test_bits_outside_the_filter_are_ignored(self):
        assert accepted_by(ECL_VISIBLE) == set()
        assert accepted_by(ECL_TOTAL | ECL_VISIBLE) == accepted_by(ECL_TOTAL)

    @pytest.mark.parametrize(
        "mask",
        [
            ECL_CENTRAL | ECL_PARTIAL,
            ECL_NONCENTRAL | ECL_ANNULAR_TOTAL,
            ECL_CENTRAL | ECL_PARTIAL | ECL_VISIBLE,
        ],
    )
    def test_impossible_filters_are_refused(self, mask):
        with pytest.raises(Error, match="cannot occur"):
            _sol_glob_reject_impossible(mask)

    def test_every_other_filter_passes(self):
        refused = {ECL_CENTRAL | ECL_PARTIAL, ECL_NONCENTRAL | ECL_ANNULAR_TOTAL}
        for mask in range(FILTER_BITS + 1):
            if mask not in refused:
                _sol_glob_reject_impossible(mask)


class TestOccultationFilter:
    EVERY_COMMON = ECL_TOTAL | ECL_PARTIAL | BOTH
    EVERY_SOLAR = EVERY_COMMON | ECL_ANNULAR | ECL_ANNULAR_TOTAL

    @pytest.mark.parametrize(
        ("ecltype", "is_sun", "expected"),
        [
            (0, True, EVERY_SOLAR),
            (0, False, EVERY_COMMON),
            (ECL_TOTAL, False, ECL_TOTAL | BOTH),
            (ECL_TOTAL | ECL_CENTRAL, False, ECL_TOTAL | BOTH),
            (ECL_PARTIAL, False, ECL_PARTIAL | ECL_NONCENTRAL),
            (ECL_PARTIAL | ECL_NONCENTRAL, True, ECL_PARTIAL | ECL_NONCENTRAL),
            (ECL_ANNULAR, True, ECL_ANNULAR | BOTH),
            (ECL_ANNULAR_TOTAL, True, ECL_ANNULAR_TOTAL | BOTH),
            # solar-only classes are dropped for any other body
            (ECL_ANNULAR | ECL_TOTAL, False, ECL_TOTAL | BOTH),
            (ECL_ANNULAR | ECL_PARTIAL, False, ECL_PARTIAL | ECL_NONCENTRAL),
            # ... and an empty remainder means "every class this body shows"
            (ECL_ANNULAR | ECL_ANNULAR_TOTAL, False, EVERY_COMMON),
            # a centrality with no class left is returned as it stands
            (ECL_ANNULAR_TOTAL | ECL_CENTRAL, False, ECL_CENTRAL),
            (ECL_CENTRAL, False, ECL_CENTRAL),
        ],
    )
    def test_expansion(self, ecltype, is_sun, expected):
        assert _normalize_occultation_filter(ecltype, 2, is_sun) == expected

    @pytest.mark.parametrize("is_sun", [True, False])
    def test_central_partial_is_refused_for_every_body(self, is_sun):
        with pytest.raises(Error, match="central partial event"):
            _normalize_occultation_filter(ECL_PARTIAL | ECL_CENTRAL, 2, is_sun)

    @pytest.mark.parametrize(
        "ecltype",
        [
            ECL_ANNULAR,
            ECL_ANNULAR | ECL_CENTRAL,
            ECL_ANNULAR | ECL_NONCENTRAL,
            ECL_ANNULAR | BOTH,
            ECL_ANNULAR_TOTAL,
        ],
    )
    def test_annular_only_is_refused_for_a_non_solar_body(self, ecltype):
        with pytest.raises(Error, match="annular event of body 2"):
            _normalize_occultation_filter(ecltype, 2, False)
        # the same request is legitimate for the Sun
        solar = _normalize_occultation_filter(ecltype, 0, True)
        assert solar & (ECL_ANNULAR | ECL_ANNULAR_TOTAL)
        assert solar & BOTH


class TestConjunctionGaps:
    def test_ahead_heads_the_way_of_the_search(self):
        assert _longitude_gap_ahead(10.0, 350.0, 1) == 20.0
        assert _longitude_gap_ahead(350.0, 10.0, 1) == 340.0
        assert _longitude_gap_ahead(10.0, 350.0, -1) == -340.0
        assert _longitude_gap_ahead(350.0, 10.0, -1) == -20.0

    def test_ahead_at_conjunction(self):
        assert _longitude_gap_ahead(100.0, 100.0, 1) == 0.0
        assert _longitude_gap_ahead(100.0, 100.0, -1) == -360.0

    def test_nearest_folds_into_half_a_turn(self):
        assert _longitude_gap_nearest(10.0, 350.0) == 20.0
        assert _longitude_gap_nearest(350.0, 10.0) == -20.0
        assert _longitude_gap_nearest(190.0, 0.0) == -170.0
        assert _longitude_gap_nearest(180.0, 0.0) == 180.0


class TestResultLayout:
    def test_shadow_center_geopos(self):
        geopos = _shadow_center_geopos(12.5, 41.9)
        assert geopos == (12.5, 41.9) + (0.0,) * 8
        assert all(type(v) is float for v in geopos)

    def test_field_order_is_slot_order(self):
        assert _LocalCircumstances._fields == (
            "magnitude",
            "diameter_ratio",
            "obscuration",
            "core_shadow_km",
            "azimuth",
            "true_altitude",
            "apparent_altitude",
            "separation_deg",
            "nasa_magnitude",
            "saros_series",
            "saros_member",
        )

    def test_from_how_core_names_the_slots(self):
        core = [float(i) for i in range(20)]
        core[3] = 0.0  # the core leaves this slot to its callers
        circumstances = _LocalCircumstances.from_how_core(core, -123.0)
        assert circumstances.core_shadow_km == -123.0
        assert circumstances.as_attr() == (
            (0.0, 1.0, 2.0, -123.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0) + (0.0,) * 9
        )
        assert len(circumstances.as_attr()) == 20

    def test_without_eclipse_keeps_the_sky(self):
        core = [float(i) for i in range(20)]
        cleared = _LocalCircumstances.from_how_core(core, -123.0).without_eclipse()
        assert cleared.as_attr() == (
            (0.0, 0.0, 0.0, 0.0, 4.0, 5.0, 6.0, 7.0, 0.0, 0.0, 0.0) + (0.0,) * 9
        )
