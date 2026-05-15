"""
LEB vs Skyfield Comparison: Future Dates (Extended Tier).

Stress-tests LEB precision for future dates (3000-5000).
Validates that Chebyshev approximation error does not grow excessively
for dates far into the future.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED

from tests.test_leb.compare.conftest import (
    ICRS_PLANETS,
    ECLIPTIC_BODIES,
    ECLIPTIC_TOLERANCES,
    CompareHelper,
    lon_error_arcsec,
)

from .conftest import TOLS_EXT, NUTATION_FLOOR_ARCSEC


class TestFuturePlanetPosition:
    """Planet position precision in future dates (3000-5000)."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_future_longitude(
        self,
        compare: CompareHelper,
        ext_future_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Planet longitude stays within tolerance for future dates."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_future_dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.POSITION_ARCSEC, (
            f'{body_name}: max future lon error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestFuturePlanetSpeed:
    """Planet speed precision in future dates."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_future_speed(
        self,
        compare: CompareHelper,
        ext_future_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Planet speed stays within tolerance for future dates."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_future_dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            err = abs(ref[3] - leb[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.SPEED_LON_DEG_DAY, (
            f"{body_name}: max future speed error = {max_err:.6f} deg/day "
            f"at JD {worst_jd:.1f}"
        )


class TestFutureEcliptic:
    """Ecliptic body precision in future dates."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ECLIPTIC_BODIES)
    def test_future_ecliptic_longitude(
        self,
        compare: CompareHelper,
        ext_future_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Ecliptic body longitude stays within tolerance for future dates."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_future_dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        per_body = ECLIPTIC_TOLERANCES.get(body_id, {}).get(
            "lon", TOLS_EXT.ECLIPTIC_ARCSEC
        )
        tol = max(per_body, NUTATION_FLOOR_ARCSEC)
        assert max_err < tol, (
            f'{body_name}: max future lon error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )
