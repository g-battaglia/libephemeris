"""
LEB vs Skyfield Comparison: Uranian Planets + Transpluto (Medium Tier).

Validates all 9 Pipeline C (heliocentric analytical) bodies
across the medium tier range (1560-2640).
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED

from tests.test_leb.compare.conftest import (
    HYPOTHETICAL_BODIES,
    CompareHelper,
    lon_error_arcsec,
)

from .conftest import TOLS_MEDIUM


class TestMediumHypotheticalPosition:
    """Hypothetical body position precision."""

    @pytest.mark.leb_compare_medium
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", HYPOTHETICAL_BODIES)
    def test_position(
        self,
        compare: CompareHelper,
        medium_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Hypothetical body longitude and latitude match Skyfield."""
        max_lon_err = 0.0
        max_lat_err = 0.0
        worst_jd = 0.0

        for jd in medium_dates_300:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            lon_err = lon_error_arcsec(ref[0], leb[0])
            lat_err = abs(ref[1] - leb[1]) * 3600.0

            if lon_err > max_lon_err:
                max_lon_err = lon_err
                worst_jd = jd
            if lat_err > max_lat_err:
                max_lat_err = lat_err

        assert max_lon_err < TOLS_MEDIUM.HYPOTHETICAL_ARCSEC, (
            f'{body_name}: max lon error = {max_lon_err:.4f}" at JD {worst_jd:.1f}'
        )
        assert max_lat_err < TOLS_MEDIUM.HYPOTHETICAL_ARCSEC, (
            f'{body_name}: max lat error = {max_lat_err:.4f}"'
        )


class TestMediumHypotheticalSpeed:
    """Hypothetical body velocity precision."""

    @pytest.mark.leb_compare_medium
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", HYPOTHETICAL_BODIES)
    def test_speed(
        self,
        compare: CompareHelper,
        medium_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        """Hypothetical body speed matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in medium_dates_300:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            err = abs(ref[3] - leb[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_MEDIUM.SPEED_LON_DEG_DAY, (
            f"{body_name}: max speed error = {max_err:.6f} deg/day at JD {worst_jd:.1f}"
        )
