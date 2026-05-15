"""
Cross-Tier Consistency: All Three Tiers.

Validates that base, medium, and extended LEB files all produce consistent
results in their common overlap range (1860-2140).

This is the ultimate consistency check: same body, same date, three different
LEB files should all agree within tolerance.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED

from tests.test_leb.compare.conftest import (
    ICRS_PLANETS,
    ECLIPTIC_BODIES,
    HYPOTHETICAL_BODIES,
    CompareHelper,
    leb_file_path,
    lon_error_arcsec,
)

from .conftest import TOLS_CROSS

# Representative bodies from each pipeline
ALL_TIER_BODIES = [
    (0, "Sun"),
    (1, "Moon"),
    (4, "Mars"),
    (5, "Jupiter"),
    (6, "Saturn"),
    (9, "Pluto"),
    (10, "MeanNode"),
    (11, "TrueNode"),
    (12, "MeanApogee"),
    (40, "Cupido"),
    (48, "Transpluto"),
]


class TestAllTierConsistency:
    """Three-way consistency for representative bodies."""

    @pytest.mark.leb_compare_crosstier
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ALL_TIER_BODIES)
    def test_three_tier_longitude(
        self,
        leb_file_base: str,
        leb_file: str,
        leb_file_extended: str,
        all_overlap_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """All three tiers agree on longitude within tolerance."""
        helper_base = CompareHelper(leb_file_base)
        helper_medium = CompareHelper(leb_file)
        helper_extended = CompareHelper(leb_file_extended)
        helper_base.setup()

        try:
            max_bm_err = 0.0
            max_me_err = 0.0
            max_be_err = 0.0

            for jd in all_overlap_dates:
                base, _ = helper_base.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)
                medium, _ = helper_medium.leb(
                    ephem.calc_ut, jd, body_id, FLG_SPEED
                )
                extended, _ = helper_extended.leb(
                    ephem.calc_ut, jd, body_id, FLG_SPEED
                )

                bm = lon_error_arcsec(base[0], medium[0])
                me = lon_error_arcsec(medium[0], extended[0])
                be = lon_error_arcsec(base[0], extended[0])

                max_bm_err = max(max_bm_err, bm)
                max_me_err = max(max_me_err, me)
                max_be_err = max(max_be_err, be)

            assert max_bm_err < TOLS_CROSS.POSITION_ARCSEC, (
                f'{body_name}: base-medium = {max_bm_err:.4f}"'
            )
            assert max_me_err < TOLS_CROSS.POSITION_ARCSEC, (
                f'{body_name}: medium-extended = {max_me_err:.4f}"'
            )
            assert max_be_err < TOLS_CROSS.POSITION_ARCSEC, (
                f'{body_name}: base-extended = {max_be_err:.4f}"'
            )
        finally:
            helper_base.teardown()


class TestAllTierDeltaT:
    """Delta-T consistency across all three tiers."""

    @pytest.mark.leb_compare_crosstier
    @pytest.mark.slow
    def test_deltat_consistency(
        self,
        leb_file_base: str,
        leb_file: str,
        leb_file_extended: str,
        all_overlap_dates: list[float],
    ):
        """Delta-T values agree across all three tiers."""
        helper_base = CompareHelper(leb_file_base)
        helper_medium = CompareHelper(leb_file)
        helper_extended = CompareHelper(leb_file_extended)
        helper_base.setup()

        try:
            max_err = 0.0
            for jd in all_overlap_dates:
                dt_base = helper_base.leb(ephem.deltat, jd)
                dt_medium = helper_medium.leb(ephem.deltat, jd)
                dt_extended = helper_extended.leb(ephem.deltat, jd)

                bm = abs(dt_base - dt_medium) * 86400.0  # seconds
                me = abs(dt_medium - dt_extended) * 86400.0
                max_err = max(max_err, bm, me)

            assert max_err < TOLS_CROSS.DELTAT_SEC, (
                f"Delta-T cross-tier diff = {max_err:.4f}s"
            )
        finally:
            helper_base.teardown()
