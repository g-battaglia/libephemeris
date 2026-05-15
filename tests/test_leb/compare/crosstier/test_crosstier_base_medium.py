"""
Cross-Tier Consistency: Base vs Medium.

Validates that the base and medium LEB files produce consistent results
in their overlapping date range (1860-2140).
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED
from libephemeris.exceptions import EphemerisRangeError

from tests.test_leb.compare.conftest import (
    ALL_LEB_BODIES,
    ICRS_PLANETS,
    lon_error_arcsec,
)

from .conftest import TOLS_CROSS, CrossTierHelper


class TestBaseMediumPosition:
    """Position consistency between base and medium tier."""

    @pytest.mark.leb_compare_crosstier
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ALL_LEB_BODIES)
    def test_longitude_consistency(
        self,
        cross_base_medium: CrossTierHelper,
        base_medium_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Longitude from base matches medium within cross-tier tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_medium_dates:
            try:
                base, _ = cross_base_medium.tier_a(
                    ephem.calc_ut, jd, body_id, FLG_SPEED
                )
                medium, _ = cross_base_medium.tier_b(
                    ephem.calc_ut, jd, body_id, FLG_SPEED
                )
            except (KeyError, ValueError, EphemerisRangeError):
                continue

            err = lon_error_arcsec(base[0], medium[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_CROSS.POSITION_ARCSEC, (
            f'{body_name}: base-medium lon diff = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestBaseMediumSpeed:
    """Speed consistency between base and medium tier."""

    @pytest.mark.leb_compare_crosstier
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ALL_LEB_BODIES)
    def test_speed_consistency(
        self,
        cross_base_medium: CrossTierHelper,
        base_medium_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Speed from base matches medium within cross-tier tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_medium_dates:
            try:
                base, _ = cross_base_medium.tier_a(
                    ephem.calc_ut, jd, body_id, FLG_SPEED
                )
                medium, _ = cross_base_medium.tier_b(
                    ephem.calc_ut, jd, body_id, FLG_SPEED
                )
            except (KeyError, ValueError, EphemerisRangeError):
                continue

            err = abs(base[3] - medium[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_CROSS.SPEED_LON_DEG_DAY, (
            f"{body_name}: base-medium speed diff = {max_err:.6f} deg/day "
            f"at JD {worst_jd:.1f}"
        )


class TestBaseMediumDistance:
    """Distance consistency between base and medium tier."""

    @pytest.mark.leb_compare_crosstier
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_distance_consistency(
        self,
        cross_base_medium: CrossTierHelper,
        base_medium_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Distance from base matches medium within cross-tier tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_medium_dates:
            base, _ = cross_base_medium.tier_a(
                ephem.calc_ut, jd, body_id, FLG_SPEED
            )
            medium, _ = cross_base_medium.tier_b(
                ephem.calc_ut, jd, body_id, FLG_SPEED
            )

            err = abs(base[2] - medium[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_CROSS.DISTANCE_AU, (
            f"{body_name}: base-medium dist diff = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )
