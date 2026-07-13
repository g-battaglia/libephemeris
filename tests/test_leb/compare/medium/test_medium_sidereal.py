"""
LEB vs Skyfield Comparison: Sidereal Mode (Medium Tier).

Validates sidereal positions across predefined and custom ayanamsha modes
across the medium tier range (1560-2640).
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
    FLG_SIDEREAL,
)

from tests.test_leb.compare.conftest import (
    SIDEREAL_MODE_SWEEP,
    CompareHelper,
    generate_test_dates,
    lon_error_arcsec,
    year_to_jd,
)

from .conftest import TOLS_MEDIUM

SIDEREAL_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
]


@pytest.fixture(scope="module")
def medium_sidereal_dates() -> list[float]:
    """50 dates in 1560-2640 range for sidereal tests."""
    return generate_test_dates(50, year_to_jd(1560), year_to_jd(2640))


class TestMediumSiderealPosition:
    """Sidereal position precision across the public mode sweep."""

    @pytest.mark.leb_compare_medium
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", SIDEREAL_BODIES)
    @pytest.mark.parametrize("sid_mode", SIDEREAL_MODE_SWEEP)
    def test_sidereal_position(
        self,
        compare: CompareHelper,
        medium_sidereal_dates: list[float],
        body_id: int,
        body_name: str,
        sid_mode: int,
    ):
        """Sidereal longitude matches Skyfield within tolerance."""
        flags = FLG_SPEED | FLG_SIDEREAL
        max_err = 0.0
        worst_jd = 0.0

        for jd in medium_sidereal_dates:
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_MEDIUM.SIDEREAL_ARCSEC, (
            f'{body_name} sidereal mode {sid_mode}: max error = {max_err:.4f}" '
            f"at JD {worst_jd:.1f}"
        )


class TestMediumSiderealSpeed:
    """Sidereal speed precision (includes precession correction)."""

    @pytest.mark.leb_compare_medium
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", SIDEREAL_BODIES)
    @pytest.mark.parametrize("sid_mode", [0, 1, 2, 3])
    def test_sidereal_speed(
        self,
        compare: CompareHelper,
        medium_sidereal_dates: list[float],
        body_id: int,
        body_name: str,
        sid_mode: int,
    ):
        """Sidereal speed matches Skyfield within tolerance."""
        flags = FLG_SPEED | FLG_SIDEREAL
        max_err = 0.0

        for jd in medium_sidereal_dates:
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = abs(ref[3] - leb[3])
            max_err = max(max_err, err)

        assert max_err < TOLS_MEDIUM.SPEED_LON_DEG_DAY, (
            f"{body_name} sidereal mode {sid_mode}: max speed error = "
            f"{max_err:.6f} deg/day"
        )
