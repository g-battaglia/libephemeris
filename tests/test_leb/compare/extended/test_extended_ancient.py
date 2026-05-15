"""
LEB vs Skyfield Comparison: Ancient Dates (Extended Tier).

Stress-tests LEB precision for ancient dates (-5000 to -1000).
Validates that Chebyshev approximation error does not grow excessively
for dates far from the modern epoch.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import SUN, MOON, MARS, JUPITER, FLG_SPEED

from tests.test_leb.compare.conftest import (
    ICRS_PLANETS,
    ECLIPTIC_BODIES,
    ECLIPTIC_TOLERANCES,
    CompareHelper,
    lon_error_arcsec,
)

from .conftest import TOLS_EXT, NUTATION_FLOOR_ARCSEC

# Core planets to stress-test in ancient dates
ANCIENT_PLANETS = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
]


class TestAncientPlanetPosition:
    """Planet position precision in ancient dates (-5000 to -1000)."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_ancient_longitude(
        self,
        compare: CompareHelper,
        ext_ancient_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Planet longitude stays within tolerance for ancient dates."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_ancient_dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.POSITION_ARCSEC, (
            f'{body_name}: max ancient lon error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestAncientPlanetSpeed:
    """Planet speed precision in ancient dates."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_ancient_speed(
        self,
        compare: CompareHelper,
        ext_ancient_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Planet speed stays within tolerance for ancient dates."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_ancient_dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            err = abs(ref[3] - leb[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.SPEED_LON_DEG_DAY, (
            f"{body_name}: max ancient speed error = {max_err:.6f} deg/day "
            f"at JD {worst_jd:.1f}"
        )


class TestAncientEcliptic:
    """Ecliptic body precision in ancient dates."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ECLIPTIC_BODIES)
    def test_ancient_ecliptic_longitude(
        self,
        compare: CompareHelper,
        ext_ancient_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Ecliptic body longitude stays within tolerance for ancient dates."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_ancient_dates:
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
            f'{body_name}: max ancient lon error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )
