"""
LEB vs Skyfield Comparison: Sidereal Mode.

Validates sidereal positions for all 27 formula-based ayanamsha modes.
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

from .conftest import (
    TOLS,
    FORMULA_SIDEREAL_MODES,
    STAR_BASED_SIDEREAL_MODES,
    CompareHelper,
    lon_error_arcsec,
    year_to_jd,
    generate_test_dates,
)

SIDEREAL_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
]


@pytest.fixture(scope="module")
def sidereal_dates():
    """20 dates in 1900-2100 range for sidereal tests."""
    return generate_test_dates(20, year_to_jd(1900), year_to_jd(2100))


class TestSiderealPosition:
    """Sidereal position precision for formula-based modes."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", SIDEREAL_BODIES)
    @pytest.mark.parametrize("sid_mode", FORMULA_SIDEREAL_MODES)
    def test_sidereal_position(
        self,
        compare: CompareHelper,
        sidereal_dates: list[float],
        body_id: int,
        body_name: str,
        sid_mode: int,
    ):
        """Sidereal longitude matches Skyfield within tolerance."""
        flags = FLG_SPEED | FLG_SIDEREAL
        max_err = 0.0
        worst_jd = 0.0

        for jd in sidereal_dates:
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS.SIDEREAL_ARCSEC, (
            f'{body_name} sidereal mode {sid_mode}: max error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestSiderealSpeed:
    """Sidereal speed precision (includes precession correction)."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", SIDEREAL_BODIES)
    @pytest.mark.parametrize("sid_mode", [0, 1, 2, 3])  # Test a subset for speed
    def test_sidereal_speed(
        self,
        compare: CompareHelper,
        sidereal_dates: list[float],
        body_id: int,
        body_name: str,
        sid_mode: int,
    ):
        """Sidereal speed matches Skyfield within tolerance."""
        flags = FLG_SPEED | FLG_SIDEREAL
        max_err = 0.0

        for jd in sidereal_dates:
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = abs(ref[3] - leb[3])
            max_err = max(max_err, err)

        assert max_err < TOLS.SPEED_LON_DEG_DAY, (
            f"{body_name} sidereal mode {sid_mode}: max speed error = {max_err:.6f} deg/day"
        )


class TestStarBasedFallback:
    """Star-based sidereal modes should produce identical results (both fallback)."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("body_id,body_name", [(SUN, "Sun")])
    @pytest.mark.parametrize("sid_mode", STAR_BASED_SIDEREAL_MODES)
    def test_star_based_identical(
        self,
        compare: CompareHelper,
        sidereal_dates: list[float],
        body_id: int,
        body_name: str,
        sid_mode: int,
    ):
        """Star-based modes produce identical results (both fallback to Skyfield)."""
        flags = FLG_SPEED | FLG_SIDEREAL

        for jd in sidereal_dates[:5]:  # Fewer dates for fallback modes
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            assert ref[0] == pytest.approx(leb[0], rel=1e-10), (
                f"{body_name} star-based mode {sid_mode} differs at JD {jd}"
            )
