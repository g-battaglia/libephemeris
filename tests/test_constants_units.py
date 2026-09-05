"""The AU conversion factors are the correctly rounded doubles of their definitions."""

from __future__ import annotations

import math

import libephemeris as ephem
from libephemeris import constants


def test_au_in_km_is_the_exact_iau_2012_value():
    assert constants.AUNIT_TO_KM == 149597870.7
    assert ephem.AUNIT_TO_KM == constants.AUNIT_TO_KM


def test_au_per_light_year_is_derived_from_c_and_the_julian_year():
    expected = 149597870700.0 / (299792458.0 * 365.25 * 86400.0)
    assert constants.AUNIT_TO_LIGHTYEAR == expected
    # The former decimal literal was a truncated copy, 5.9e-14 away.
    assert math.isclose(
        constants.AUNIT_TO_LIGHTYEAR, 1.5812507409819728e-05, rel_tol=1e-13
    )
    assert constants.AUNIT_TO_LIGHTYEAR != 1.5812507409819728e-05


def test_au_per_parsec_is_pi_over_648000():
    assert constants.AUNIT_TO_PARSEC == math.pi / 648000.0
    # 1 pc = 648000/pi au (IAU 2015 B2): the reciprocal must round-trip.
    assert math.isclose(
        1.0 / constants.AUNIT_TO_PARSEC, 648000.0 / math.pi, rel_tol=1e-15
    )
    assert math.isclose(constants.AUNIT_TO_PARSEC, 4.848136811095274e-06, rel_tol=1e-13)
