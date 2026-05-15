"""
LEB vs Skyfield Comparison: Flag Combinations (Extended Tier).

Validates all flag combinations supported natively by LEB
across the extended tier range (-5000 to 5000).
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    JUPITER,
    FLG_SPEED,
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_HELCTR,
    FLG_BARYCTR,
    FLG_TRUEPOS,
    FLG_NOABERR,
)

from tests.test_leb.compare.conftest import (
    CompareHelper,
    lon_error_arcsec,
)

from .conftest import TOLS_EXT

FLAG_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
]

FLAG_COMBINATIONS = [
    (FLG_SPEED, "default"),
    (FLG_SPEED | FLG_EQUATORIAL, "equatorial"),
    (FLG_SPEED | FLG_J2000, "J2000"),
    (FLG_SPEED | FLG_EQUATORIAL | FLG_J2000, "equatorial_J2000"),
    (FLG_SPEED | FLG_HELCTR, "heliocentric"),
    (FLG_SPEED | FLG_BARYCTR, "barycentric"),
    (FLG_SPEED | FLG_TRUEPOS, "truepos"),
    (FLG_SPEED | FLG_NOABERR, "noaberr"),
]


class TestExtFlagCombinations:
    """All 8 flag combinations for 4 planets."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", FLAG_BODIES)
    @pytest.mark.parametrize("flags,flag_name", FLAG_COMBINATIONS)
    def test_flags(
        self,
        compare: CompareHelper,
        ext_dates_100: list[float],
        body_id: int,
        body_name: str,
        flags: int,
        flag_name: str,
    ):
        """All 6 components match within tolerance for each flag combo."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_100:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            lon_err = lon_error_arcsec(ref[0], leb[0])
            lat_err = abs(ref[1] - leb[1]) * 3600.0
            err = max(lon_err, lat_err)

            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.EQUATORIAL_ARCSEC, (
            f'{body_name} {flag_name}: max error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestExtFlagVelocity:
    """Velocity under different flag combinations (missing from extended tier)."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", FLAG_BODIES)
    @pytest.mark.parametrize(
        "flags,flag_name",
        [
            (FLG_SPEED | FLG_EQUATORIAL, "equatorial"),
            (FLG_SPEED | FLG_J2000, "J2000"),
            (FLG_SPEED | FLG_HELCTR, "heliocentric"),
            (FLG_SPEED | FLG_BARYCTR, "barycentric"),
        ],
    )
    def test_flag_speed(
        self,
        compare: CompareHelper,
        ext_dates_100: list[float],
        body_id: int,
        body_name: str,
        flags: int,
        flag_name: str,
    ):
        """Speed matches within tolerance for each flag combo."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_100:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = abs(ref[3] - leb[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.SPEED_LON_DEG_DAY, (
            f"{body_name} {flag_name}: max speed error = {max_err:.6f} deg/day "
            f"at JD {worst_jd:.1f}"
        )

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", FLAG_BODIES)
    @pytest.mark.parametrize(
        "flags,flag_name",
        [
            (FLG_SPEED | FLG_EQUATORIAL, "equatorial"),
            (FLG_SPEED | FLG_J2000, "J2000"),
            (FLG_SPEED | FLG_EQUATORIAL | FLG_J2000, "equatorial_J2000"),
        ],
    )
    def test_flag_lat_speed(
        self,
        compare: CompareHelper,
        ext_dates_100: list[float],
        body_id: int,
        body_name: str,
        flags: int,
        flag_name: str,
    ):
        """Latitude speed matches within tolerance for each flag combo."""
        max_err = 0.0

        for jd in ext_dates_100:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = abs(ref[4] - leb[4])
            max_err = max(max_err, err)

        assert max_err < TOLS_EXT.SPEED_LAT_DEG_DAY, (
            f"{body_name} {flag_name}: max lat speed error = {max_err:.6f} deg/day"
        )

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", FLAG_BODIES)
    @pytest.mark.parametrize(
        "flags,flag_name",
        [
            (FLG_SPEED, "default"),
            (FLG_SPEED | FLG_HELCTR, "heliocentric"),
            (FLG_SPEED | FLG_BARYCTR, "barycentric"),
        ],
    )
    def test_flag_dist_speed(
        self,
        compare: CompareHelper,
        ext_dates_100: list[float],
        body_id: int,
        body_name: str,
        flags: int,
        flag_name: str,
    ):
        """Distance speed matches within tolerance for each flag combo."""
        max_err = 0.0

        for jd in ext_dates_100:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = abs(ref[5] - leb[5])
            max_err = max(max_err, err)

        assert max_err < TOLS_EXT.SPEED_DIST_AU_DAY, (
            f"{body_name} {flag_name}: max dist speed error = {max_err:.2e} AU/day"
        )
