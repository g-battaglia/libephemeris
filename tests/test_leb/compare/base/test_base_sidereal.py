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

from .conftest import TOLS_BASE

SIDEREAL_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
]


@pytest.fixture(scope="module")
def base_sidereal_dates() -> list[float]:
    return generate_test_dates(50, year_to_jd(1860), year_to_jd(2140))


class TestBaseSiderealPrecision:
    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", SIDEREAL_BODIES)
    @pytest.mark.parametrize("sid_mode", SIDEREAL_MODE_SWEEP)
    def test_sidereal_all(
        self,
        compare: CompareHelper,
        base_sidereal_dates: list[float],
        body_id: int,
        body_name: str,
        sid_mode: int,
    ):
        flags = FLG_SPEED | FLG_SIDEREAL
        max_lon_err = 0.0
        max_speed_err = 0.0
        worst_lon_jd = 0.0

        for jd in base_sidereal_dates:
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            lon_err = lon_error_arcsec(ref[0], leb[0])
            speed_err = abs(ref[3] - leb[3])

            if lon_err > max_lon_err:
                max_lon_err = lon_err
                worst_lon_jd = jd
            max_speed_err = max(max_speed_err, speed_err)

        assert max_lon_err < TOLS_BASE.SIDEREAL_ARCSEC, (
            f'{body_name} sidereal mode {sid_mode}: max error = {max_lon_err:.4f}" '
            f"at JD {worst_lon_jd:.1f}"
        )
        assert max_speed_err < TOLS_BASE.SPEED_LON_DEG_DAY, (
            f"{body_name} sidereal mode {sid_mode}: max speed error = "
            f"{max_speed_err:.6f} deg/day"
        )
