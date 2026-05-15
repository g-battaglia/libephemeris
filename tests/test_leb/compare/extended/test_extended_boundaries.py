"""
LEB vs Skyfield Comparison: Range Boundaries (Extended Tier).

Tests precision at the extreme edges of the extended LEB file range.
Validates that Chebyshev polynomial clamping and segment boundaries
do not introduce excessive error at the start/end of coverage.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import SUN, MOON, MARS, JUPITER, FLG_SPEED

from tests.test_leb.compare.conftest import (
    ICRS_PLANETS,
    CompareHelper,
    lon_error_arcsec,
)

from .conftest import TOLS_EXT

BOUNDARY_PLANETS = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
]


class TestBoundaryPosition:
    """Planet position precision at range boundaries."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", BOUNDARY_PLANETS)
    def test_boundary_longitude(
        self,
        compare: CompareHelper,
        ext_boundary_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Planet longitude stays within tolerance at range boundaries."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_boundary_dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.POSITION_ARCSEC, (
            f'{body_name}: max boundary lon error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestBoundarySpeed:
    """Planet speed precision at range boundaries."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", BOUNDARY_PLANETS)
    def test_boundary_speed(
        self,
        compare: CompareHelper,
        ext_boundary_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Planet speed stays within tolerance at range boundaries."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_boundary_dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            err = abs(ref[3] - leb[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.SPEED_LON_DEG_DAY, (
            f"{body_name}: max boundary speed error = {max_err:.6f} deg/day "
            f"at JD {worst_jd:.1f}"
        )


class TestBoundaryAllPlanets:
    """All ICRS planets at boundaries (broader check)."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_boundary_all_components(
        self,
        compare: CompareHelper,
        ext_boundary_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """All components stay within tolerance at range boundaries."""
        max_lon_err = 0.0
        max_dist_err = 0.0

        for jd in ext_boundary_dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            lon_err = lon_error_arcsec(ref[0], leb[0])
            dist_err = abs(ref[2] - leb[2])

            max_lon_err = max(max_lon_err, lon_err)
            max_dist_err = max(max_dist_err, dist_err)

        assert max_lon_err < TOLS_EXT.POSITION_ARCSEC, (
            f'{body_name}: boundary lon error {max_lon_err:.4f}"'
        )
        assert max_dist_err < TOLS_EXT.DISTANCE_AU, (
            f"{body_name}: boundary dist error {max_dist_err:.2e} AU"
        )
