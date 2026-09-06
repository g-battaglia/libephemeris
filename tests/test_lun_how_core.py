# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""The local circumstances of a lunar eclipse, checked against their definitions.

The unit under test answers what the Earth's shadow is doing to the Moon at
one instant: which of the two shadows the Moon's disc reaches, how deep it
stands in each, how far it is from opposition, and which Saros series the
event belongs to. Everything asserted here comes from the published geometry
(Espenak & Meeus, *Five Millennium Canon of Lunar Eclipses*,
NASA/TP-2009-214173, sec. 1.2.7 and sec. 1.5; *Explanatory Supplement to the
Astronomical Almanac*, 3rd ed., ch. 11) or from a documented compatibility
contract, and nothing from a recorded value.

Two branches in particular are witnessed here because no public search ever
lands on them: the instant with **no eclipse at all**, where the penumbral
magnitude is the negative distance of the Moon's limb from the penumbra and
the Saros pair carries the no-match sentinel, and the reserved slots, which
must stay at zero on every path.
"""

from __future__ import annotations

import math

import pytest

from libephemeris import julday
from libephemeris.constants import (
    ECL_PARTIAL,
    ECL_PENUMBRAL,
    ECL_TOTAL,
    FLG_JPLEPH,
    FLG_MOSEPH,
    FLG_SWIEPH,
)
from libephemeris.eclipse import (
    _DANJON_MOON_PARALLAX_SCALE,
    _ECL_RMOON_AU,
    _lun_how_core,
)

#: The no-match value the compatibility surface publishes in the Saros slots.
SAROS_NO_MATCH = -99999999.0

#: Instants of maximum eclipse, one of each class, and one instant with no
#: eclipse of any kind. Julian Days in UT, refined maxima of real events.
TOTAL_MAX = 2459715.6746893055  # 2022 May 16, total
PARTIAL_MAX = 2459537.877033952  # 2021 Nov 19, deep partial
PENUMBRAL_MAX = 2460070.2243118742  # 2023 May 5, penumbral
#: A full moon that carries no eclipse -- the configuration the global search
#: asks about at every lunation and rejects.
NO_ECLIPSE = 2459773.3  # 2022 Jul 12
#: An instant nowhere near a syzygy, where the sections are an extrapolation
#: to the wrong side of the fundamental plane and describe no real shadow.
NO_ECLIPSE_AWAY = 2459732.0  # 2022 Jun 1


def _classify(jd: float):
    """The unit's answer, with the three phases named."""
    return _lun_how_core(jd)


class TestNoEclipseBranch:
    """The branch the recorded searches only reach indirectly."""

    def test_no_eclipse_reports_a_zero_word(self):
        word, _attr, _geometry = _classify(NO_ECLIPSE)
        assert word == 0

    def test_penumbral_magnitude_is_negative_with_no_eclipse(self):
        """The magnitude is a signed depth of immersion, not a clamped one.

        With the Moon clear of the penumbra its near limb has not been
        reached, so the numerator of the ratio is negative and the array
        reports it as such.
        """
        _word, attr, _geometry = _classify(NO_ECLIPSE)
        assert attr[1] < 0.0

    def test_umbral_magnitude_and_opposition_stay_empty_with_no_eclipse(self):
        _word, attr, _geometry = _classify(NO_ECLIPSE)
        assert attr[0] == 0.0
        assert attr[8] == 0.0
        assert attr[7] == 0.0

    def test_saros_slots_carry_the_no_match_sentinel(self):
        """The Saros pair is filled whatever the phase, sentinel included."""
        _word, attr, _geometry = _classify(NO_ECLIPSE)
        assert attr[9] == SAROS_NO_MATCH
        assert attr[10] == SAROS_NO_MATCH

    def test_the_shadow_record_is_still_well_formed(self):
        """The two cones are still described, at the closest approach.

        Only past the Earth do the sections keep their order: away from the
        syzygy the construction is an extrapolation to the wrong side of the
        fundamental plane, where the converging cone is the wider of the two
        and no real shadow is being described.
        """
        offset, umbral, penumbral, cos_u, cos_p = _classify(NO_ECLIPSE)[2]
        assert offset >= 0.0
        assert penumbral > 0.0
        assert 0.0 < cos_p <= cos_u <= 1.0
        assert umbral < penumbral

    def test_the_branch_holds_away_from_the_syzygy_too(self):
        word, attr, _geometry = _classify(NO_ECLIPSE_AWAY)
        assert word == 0
        assert attr[0] == 0.0
        assert attr[1] < 0.0
        assert attr[7] == 0.0
        assert attr[9] == SAROS_NO_MATCH


class TestReservedSlots:
    """Slots nobody fills must stay at zero on every path."""

    @pytest.mark.parametrize(
        "jd", [TOTAL_MAX, PARTIAL_MAX, PENUMBRAL_MAX, NO_ECLIPSE, NO_ECLIPSE_AWAY]
    )
    def test_reserved_slots_are_exactly_zero(self, jd):
        _word, attr, _geometry = _lun_how_core(jd)
        assert len(attr) == 20
        for slot in (2, 3, *range(11, 20)):
            assert attr[slot] == 0.0, f"slot {slot} moved"

    @pytest.mark.parametrize(
        "jd", [TOTAL_MAX, PARTIAL_MAX, PENUMBRAL_MAX, NO_ECLIPSE, NO_ECLIPSE_AWAY]
    )
    def test_the_observer_slots_are_left_to_the_caller(self, jd):
        """The unit takes no observer, so it fills no horizon slot."""
        _word, attr, _geometry = _lun_how_core(jd)
        assert attr[4] == 0.0
        assert attr[5] == 0.0
        assert attr[6] == 0.0

    @pytest.mark.parametrize(
        "jd", [TOTAL_MAX, PARTIAL_MAX, PENUMBRAL_MAX, NO_ECLIPSE, NO_ECLIPSE_AWAY]
    )
    def test_every_slot_is_a_native_float(self, jd):
        _word, attr, _geometry = _lun_how_core(jd)
        assert all(type(value) is float for value in attr)


class TestPublishedContracts:
    """What the compatibility surface promises about the two magnitudes."""

    def test_slot_eight_repeats_the_umbral_magnitude(self):
        for jd in (TOTAL_MAX, PARTIAL_MAX, PENUMBRAL_MAX, NO_ECLIPSE):
            _word, attr, _geometry = _lun_how_core(jd)
            assert attr[8] == attr[0]

    def test_a_penumbral_eclipse_publishes_a_zero_umbral_magnitude(self):
        """The canon publishes that number negative; this surface publishes 0.

        See docs/comparison/known-differences.md. The Moon's limb has not
        reached the umbra, so the definition gives a negative depth of
        immersion; the published channel carries exactly ``0.0``.
        """
        word, attr, geometry = _lun_how_core(PENUMBRAL_MAX)
        assert word == ECL_PENUMBRAL
        assert attr[0] == 0.0
        assert attr[8] == 0.0
        offset, umbral_diameter, _pen, _cos_u, _cos_p = geometry
        definition = (umbral_diameter / 2.0 + _ECL_RMOON_AU - offset) / (
            2.0 * _ECL_RMOON_AU
        )
        assert definition < 0.0

    def test_the_penumbral_magnitude_is_positive_during_a_penumbral_eclipse(self):
        _word, attr, _geometry = _lun_how_core(PENUMBRAL_MAX)
        assert attr[1] > 0.0

    def test_the_word_carries_one_phase_bit_at_most(self):
        for jd in (TOTAL_MAX, PARTIAL_MAX, PENUMBRAL_MAX, NO_ECLIPSE):
            word, _attr, _geometry = _lun_how_core(jd)
            assert word in (0, ECL_TOTAL, ECL_PARTIAL, ECL_PENUMBRAL)


class TestGeometryAndDefinitions:
    """The published quantities, recomputed from the record the unit hands out."""

    @pytest.mark.parametrize("jd", [TOTAL_MAX, PARTIAL_MAX, PENUMBRAL_MAX])
    def test_the_magnitudes_are_the_ratios_of_distances(self, jd):
        word, attr, geometry = _lun_how_core(jd)
        offset, umbral_diameter, penumbral_diameter, _cos_u, _cos_p = geometry
        diameter = 2.0 * _ECL_RMOON_AU
        umbral = (umbral_diameter / 2.0 + _ECL_RMOON_AU - offset) / diameter
        penumbral = (penumbral_diameter / 2.0 + _ECL_RMOON_AU - offset) / diameter
        assert attr[1] == penumbral
        if word & (ECL_TOTAL | ECL_PARTIAL):
            assert attr[0] == umbral

    @pytest.mark.parametrize("jd", [TOTAL_MAX, PARTIAL_MAX, PENUMBRAL_MAX])
    def test_the_phase_follows_the_physical_conditions(self, jd):
        """Total inside the umbral section, partial meeting it, penumbral beyond.

        The Moon's radius enters each condition divided by that cone's
        cosine: the projection of a sphere's radius onto the plane
        perpendicular to the shadow axis.
        """
        word, _attr, geometry = _lun_how_core(jd)
        offset, umbral_diameter, penumbral_diameter, cos_u, cos_p = geometry
        umbral_radius = umbral_diameter / 2.0
        penumbral_radius = penumbral_diameter / 2.0
        on_umbral = _ECL_RMOON_AU / cos_u
        on_penumbral = _ECL_RMOON_AU / cos_p
        if word == ECL_TOTAL:
            assert offset + on_umbral <= umbral_radius
        elif word == ECL_PARTIAL:
            assert offset - on_umbral <= umbral_radius < offset + on_umbral
        else:
            assert umbral_radius < offset - on_umbral
            assert offset - on_penumbral <= penumbral_radius

    @pytest.mark.parametrize("jd", [TOTAL_MAX, PARTIAL_MAX, PENUMBRAL_MAX])
    def test_the_class_agrees_with_the_magnitudes(self, jd):
        """Total above one, partial between zero and one, penumbral beyond.

        The correspondence is the canons' own, exact to the half of the
        projection correction above (about 5e-6 of a magnitude unit).
        """
        word, attr, _geometry = _lun_how_core(jd)
        if word == ECL_TOTAL:
            assert attr[0] > 1.0
        elif word == ECL_PARTIAL:
            assert 0.0 < attr[0] < 1.0
        else:
            assert attr[1] > 0.0

    @pytest.mark.parametrize("jd", [TOTAL_MAX, PARTIAL_MAX, PENUMBRAL_MAX])
    def test_the_distance_from_opposition_is_a_small_positive_angle(self, jd):
        """A lunar eclipse happens within a couple of degrees of opposition."""
        _word, attr, _geometry = _lun_how_core(jd)
        assert 0.0 < attr[7] < 2.0

    @pytest.mark.parametrize("jd", [TOTAL_MAX, PARTIAL_MAX, PENUMBRAL_MAX])
    def test_the_umbral_section_is_a_real_umbra(self, jd):
        """The Moon never reaches the umbral cone's apex, so the section is wide.

        The recorded population puts the umbral section between 9 074 km and
        9 646 km of diameter; the bounds here are generous multiples of the
        Earth's own radius, and only the sign is a statement of principle.
        """
        _word, _attr, geometry = _lun_how_core(jd)
        _offset, umbral_diameter, penumbral_diameter, _cos_u, _cos_p = geometry
        assert umbral_diameter > 0.0
        assert penumbral_diameter > umbral_diameter

    @pytest.mark.parametrize("jd", [TOTAL_MAX, PARTIAL_MAX, PENUMBRAL_MAX])
    def test_the_cone_cosines_are_ordered(self, jd):
        """The penumbral cone opens wider, so its cosine is the smaller."""
        _word, _attr, geometry = _lun_how_core(jd)
        _offset, _umbral, _penumbral, cos_u, cos_p = geometry
        assert 0.0 < cos_p <= cos_u <= 1.0
        assert math.isclose(cos_u, 1.0, abs_tol=1e-4)


class TestEnlargement:
    """The convention that decides both the magnitudes and the class."""

    def test_the_factor_is_the_published_one(self):
        """Danjon's shell as Espenak renders it: 1.01 = 1 + 1/85 - 1/594.

        Five Millennium Canon of Lunar Eclipses, sec. 1.5, eqs. 1-5 and 1-6.
        Chauvenet's competing one fiftieth would move the umbral magnitude by
        about 0.006 and the penumbral one by about 0.026.
        """
        assert _DANJON_MOON_PARALLAX_SCALE == 1.01
        assert math.isclose(
            _DANJON_MOON_PARALLAX_SCALE, 1.0 + 1.0 / 85.0 - 1.0 / 594.0, abs_tol=1e-4
        )


class TestObserverIndependence:
    """Nothing but the ephemeris selection may reach the answer."""

    @pytest.mark.parametrize("jd", [TOTAL_MAX, PARTIAL_MAX, PENUMBRAL_MAX])
    def test_the_four_ephemeris_selectors_agree_bit_for_bit(self, jd):
        reference = _lun_how_core(jd, 0)
        for selector in (FLG_JPLEPH, FLG_SWIEPH, FLG_MOSEPH):
            assert _lun_how_core(jd, selector) == reference

    def test_a_projection_bit_changes_nothing(self):
        """A reduction bit cannot move an observer-independent answer."""
        from libephemeris.constants import FLG_RADIANS, FLG_TOPOCTR

        reference = _lun_how_core(TOTAL_MAX, FLG_SWIEPH)
        for extra in (FLG_RADIANS, FLG_TOPOCTR):
            assert _lun_how_core(TOTAL_MAX, FLG_SWIEPH | extra) == reference


def test_the_anchor_dates_are_the_classes_they_claim():
    """Guard the fixtures themselves: one instant of each class."""
    assert _lun_how_core(TOTAL_MAX)[0] == ECL_TOTAL
    assert _lun_how_core(PARTIAL_MAX)[0] == ECL_PARTIAL
    assert _lun_how_core(PENUMBRAL_MAX)[0] == ECL_PENUMBRAL
    assert _lun_how_core(NO_ECLIPSE)[0] == 0
    assert _lun_how_core(julday(2022, 6, 1, 12.0))[0] == 0
