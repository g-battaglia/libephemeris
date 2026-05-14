"""
LEB vs Skyfield Comparison: Flag Combinations.

Validates all flag combinations supported natively by LEB.
This is the most important test for correctness of the coordinate transform pipeline.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MARS,
    SE_JUPITER,
    SEFLG_SPEED,
    SEFLG_EQUATORIAL,
    SEFLG_J2000,
    SEFLG_HELCTR,
    SEFLG_BARYCTR,
    SEFLG_TRUEPOS,
    SEFLG_NOABERR,
    SEFLG_TOPOCTR,
    SEFLG_XYZ,
    SEFLG_RADIANS,
    SEFLG_NONUT,
)

from .conftest import (
    TOLS,
    CompareHelper,
    lon_error_arcsec,
)

FLAG_BODIES = [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
]

FLAG_COMBINATIONS = [
    (SEFLG_SPEED, "default"),
    (SEFLG_SPEED | SEFLG_EQUATORIAL, "equatorial"),
    (SEFLG_SPEED | SEFLG_J2000, "J2000"),
    (SEFLG_SPEED | SEFLG_EQUATORIAL | SEFLG_J2000, "equatorial_J2000"),
    (SEFLG_SPEED | SEFLG_HELCTR, "heliocentric"),
    (SEFLG_SPEED | SEFLG_BARYCTR, "barycentric"),
    (SEFLG_SPEED | SEFLG_TRUEPOS, "truepos"),
    (SEFLG_SPEED | SEFLG_NOABERR, "noaberr"),
]


class TestFlagCombinations:
    """All 8 flag combinations for 4 planets."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", FLAG_BODIES)
    @pytest.mark.parametrize("flags,flag_name", FLAG_COMBINATIONS)
    def test_flags(
        self,
        compare: CompareHelper,
        test_dates_50: list[float],
        body_id: int,
        body_name: str,
        flags: int,
        flag_name: str,
    ):
        """All 6 components match within tolerance for each flag combo."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in test_dates_50:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, flags)

            lon_err = lon_error_arcsec(ref[0], leb[0])
            lat_err = abs(ref[1] - leb[1]) * 3600.0
            err = max(lon_err, lat_err)

            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS.EQUATORIAL_ARCSEC, (
            f'{body_name} {flag_name}: max error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestFallbackFlags:
    """Flags now handled natively in LEB — results should match Skyfield closely."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("body_id,body_name", [(SE_SUN, "Sun"), (SE_MARS, "Mars")])
    def test_topoctr_identical(
        self,
        compare: CompareHelper,
        test_dates_20: list[float],
        body_id: int,
        body_name: str,
    ):
        """SEFLG_TOPOCTR LEB vs Skyfield within sub-arcsecond tolerance."""
        flags = SEFLG_SPEED | SEFLG_TOPOCTR
        ephem.set_topo(41.9, 12.5, 0.0)

        for jd in test_dates_20:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, flags)

            assert ref[0] == pytest.approx(leb[0], abs=1e-6), (
                f"{body_name} TOPOCTR lon differs at JD {jd}"
            )

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("body_id,body_name", [(SE_SUN, "Sun"), (SE_MARS, "Mars")])
    def test_xyz_identical(
        self,
        compare: CompareHelper,
        test_dates_20: list[float],
        body_id: int,
        body_name: str,
    ):
        """SEFLG_XYZ LEB vs Skyfield — native LEB path, sub-arcsecond precision."""
        flags = SEFLG_SPEED | SEFLG_XYZ

        for jd in test_dates_20:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, flags)

            for i in range(6):
                assert ref[i] == pytest.approx(leb[i], abs=1e-4), (
                    f"{body_name} XYZ component {i} differs at JD {jd}"
                )

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("body_id,body_name", [(SE_SUN, "Sun"), (SE_MARS, "Mars")])
    def test_radians_identical(
        self,
        compare: CompareHelper,
        test_dates_20: list[float],
        body_id: int,
        body_name: str,
    ):
        """SEFLG_RADIANS LEB vs Skyfield — native LEB path."""
        flags = SEFLG_SPEED | SEFLG_RADIANS

        for jd in test_dates_20:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, flags)

            for i in range(6):
                assert ref[i] == pytest.approx(leb[i], abs=1e-4), (
                    f"{body_name} RADIANS component {i} differs at JD {jd}"
                )

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("body_id,body_name", [(SE_SUN, "Sun"), (SE_MARS, "Mars")])
    def test_nonut_identical(
        self,
        compare: CompareHelper,
        test_dates_20: list[float],
        body_id: int,
        body_name: str,
    ):
        """SEFLG_NONUT LEB vs Skyfield within sub-arcsecond tolerance."""
        flags = SEFLG_SPEED | SEFLG_NONUT

        for jd in test_dates_20:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, flags)

            assert ref[0] == pytest.approx(leb[0], abs=1e-7), (
                f"{body_name} NONUT lon differs at JD {jd}"
            )
