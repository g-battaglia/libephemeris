"""Parity between the calc_ut()/calc() fixed-star dispatch and fixstar_ut().

Fixed stars requested through the generic calc_ut()/calc() entry points (via
the FIXSTAR_OFFSET id range) delegate to the same computation as the
dedicated fixstar_ut(). This pins two properties that a past regression
broke:

1. Dispatch is BY ID, not by resolving the traditional name — Flamsteed
   names starting with digits (e.g. "29Psc") otherwise parse as sequential
   catalog numbers and return the wrong star.
2. Every flag (FLG_EQUATORIAL, FLG_J2000, FLG_XYZ, FLG_SPEED, FLG_SIDEREAL)
   gets identical handling on both paths.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris import fixed_stars

_FLAG_SETS = [
    ("tropical-ecliptic", 0),
    ("equatorial", ephem.FLG_EQUATORIAL),
    ("j2000", ephem.FLG_J2000),
    ("xyz", ephem.FLG_XYZ),
    ("speed", ephem.FLG_SPEED),
    ("sidereal", ephem.FLG_SIDEREAL),
    ("equatorial-speed", ephem.FLG_EQUATORIAL | ephem.FLG_SPEED),
    ("j2000-xyz", ephem.FLG_J2000 | ephem.FLG_XYZ),
]

# Sample of stars with UNAMBIGUOUS names for the by-name reference: fixstar_ut
# takes a name string, so a digit-leading Flamsteed name (e.g. "29Psc") would
# resolve as a sequential catalog number there too (a name-resolution
# ambiguity shared with the reference API, independent of the id dispatch).
# The digit-leading ids get their own by-id reference below.
_UNAMBIGUOUS_IDS = [
    ephem.REGULUS,
    ephem.SPICA_STAR,
] + [e.id for e in fixed_stars.STAR_CATALOG if not e.name[:1].isdigit()][:8]


@pytest.mark.unit
class TestCalcUtFixedStarDispatchParity:
    """calc_ut()/calc() must match fixstar_ut() for the same star and flags."""

    @pytest.mark.parametrize("flag_label,flags", _FLAG_SETS)
    def test_calc_ut_matches_fixstar_ut(self, flag_label, flags, standard_jd):
        """Every sampled star matches fixstar_ut() across the flag matrix."""
        if flags & ephem.FLG_SIDEREAL:
            ephem.set_sid_mode(ephem.SIDM_LAHIRI, 0, 0)
        for star_id in _UNAMBIGUOUS_IDS:
            name = fixed_stars.get_canonical_star_name(star_id)
            assert name is not None
            ref_pos, _ref_name, _rf = fixed_stars.fixstar_ut(name, standard_jd, flags)
            got_pos, _got_rf = ephem.calc_ut(standard_jd, star_id, flags)
            for i, (a, b) in enumerate(zip(ref_pos, got_pos)):
                assert a == pytest.approx(b, abs=1e-9), (
                    f"{flag_label} star_id={star_id} ({name}) component {i}: "
                    f"fixstar_ut={a} calc_ut={b}"
                )

    def test_digit_leading_name_returns_correct_star(self, standard_jd):
        """A digit-leading Flamsteed id must not resolve as a catalog number.

        Regression guard: dispatching by the traditional name '29Psc' made
        _resolve_star_ref read it as sequential catalog number 29, returning
        a different star entirely.
        """
        digit_ids = [e.id for e in fixed_stars.STAR_CATALOG if e.name[:1].isdigit()]
        if not digit_ids:
            pytest.skip("catalog has no digit-leading Flamsteed names")
        for star_id in digit_ids[:20]:
            name = fixed_stars.get_canonical_star_name(star_id)
            # jd_tt: calc_ut takes UT, calc_fixed_star_position takes TT
            jd_tt = standard_jd + ephem.deltat(standard_jd)
            direct = fixed_stars.calc_fixed_star_position(star_id, jd_tt)[0]
            via_calc = ephem.calc_ut(standard_jd, star_id, 0)[0][0]
            assert via_calc == pytest.approx(direct, abs=1e-9), (
                f"star_id={star_id} ({name}): calc_ut dispatched to the wrong "
                f"star (got {via_calc}, expected {direct})"
            )
