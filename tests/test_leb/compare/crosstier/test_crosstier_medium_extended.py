"""
Cross-Tier Consistency: Medium vs Extended.

Validates that the medium and extended LEB files produce consistent results
in their overlapping date range (1560-2640).
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED
from libephemeris.exceptions import EphemerisRangeError

from tests.test_leb.compare.conftest import (
    ECLIPTIC_BODIES,
    HYPOTHETICAL_BODIES,
    ICRS_PLANETS,
    lon_error_arcsec,
)

from .conftest import TOLS_CROSS, CrossTierHelper

# Asteroids may have different SPK coverage per tier, test cautiously
CROSSTIER_BODIES = ICRS_PLANETS + ECLIPTIC_BODIES + HYPOTHETICAL_BODIES


class TestMediumExtendedPosition:
    """Position consistency between medium and extended tier."""

    @pytest.mark.leb_compare_crosstier
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", CROSSTIER_BODIES)
    def test_longitude_consistency(
        self,
        cross_medium_extended: CrossTierHelper,
        medium_ext_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Longitude from medium matches extended within cross-tier tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in medium_ext_dates:
            try:
                medium, _ = cross_medium_extended.tier_a(
                    ephem.calc_ut, jd, body_id, FLG_SPEED
                )
                extended, _ = cross_medium_extended.tier_b(
                    ephem.calc_ut, jd, body_id, FLG_SPEED
                )
            except (KeyError, ValueError, EphemerisRangeError):
                continue

            err = lon_error_arcsec(medium[0], extended[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_CROSS.POSITION_ARCSEC, (
            f'{body_name}: medium-extended lon diff = {max_err:.4f}" '
            f"at JD {worst_jd:.1f}"
        )


class TestMediumExtendedSpeed:
    """Speed consistency between medium and extended tier."""

    @pytest.mark.leb_compare_crosstier
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", CROSSTIER_BODIES)
    def test_speed_consistency(
        self,
        cross_medium_extended: CrossTierHelper,
        medium_ext_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Speed from medium matches extended within cross-tier tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in medium_ext_dates:
            try:
                medium, _ = cross_medium_extended.tier_a(
                    ephem.calc_ut, jd, body_id, FLG_SPEED
                )
                extended, _ = cross_medium_extended.tier_b(
                    ephem.calc_ut, jd, body_id, FLG_SPEED
                )
            except (KeyError, ValueError, EphemerisRangeError):
                continue

            err = abs(medium[3] - extended[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_CROSS.SPEED_LON_DEG_DAY, (
            f"{body_name}: medium-extended speed diff = {max_err:.6f} deg/day "
            f"at JD {worst_jd:.1f}"
        )


class TestMediumExtendedDistance:
    """Distance consistency between medium and extended tier."""

    @pytest.mark.leb_compare_crosstier
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_distance_consistency(
        self,
        cross_medium_extended: CrossTierHelper,
        medium_ext_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Distance from medium matches extended within cross-tier tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in medium_ext_dates:
            medium, _ = cross_medium_extended.tier_a(
                ephem.calc_ut, jd, body_id, FLG_SPEED
            )
            extended, _ = cross_medium_extended.tier_b(
                ephem.calc_ut, jd, body_id, FLG_SPEED
            )

            err = abs(medium[2] - extended[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_CROSS.DISTANCE_AU, (
            f"{body_name}: medium-extended dist diff = {max_err:.2e} AU "
            f"at JD {worst_jd:.1f}"
        )
