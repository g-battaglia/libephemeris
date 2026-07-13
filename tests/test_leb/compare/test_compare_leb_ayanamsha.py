"""
LEB vs Skyfield Comparison: Ayanamsha Values.

Validates get_ayanamsa_ut() LEB dispatch directly.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem

from .conftest import (
    TOLS,
    ADDITIONAL_SIDEREAL_MODES,
    SIDEREAL_MODE_SWEEP,
    CompareHelper,
    year_to_jd,
    generate_test_dates,
)


@pytest.fixture(scope="module")
def ayanamsha_dates():
    """50 dates in 1560-2640 range for ayanamsha tests."""
    return generate_test_dates(50, year_to_jd(1560), year_to_jd(2640))


class TestAyanamshaValues:
    """Direct ayanamsha value comparison across the public mode sweep."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("sid_mode", SIDEREAL_MODE_SWEEP)
    def test_ayanamsha_value(
        self, compare: CompareHelper, ayanamsha_dates: list[float], sid_mode: int
    ):
        """Ayanamsha value matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ayanamsha_dates:
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            ref = compare.skyfield(ephem.get_ayanamsa_ut, jd)
            leb = compare.leb(ephem.get_ayanamsa_ut, jd)

            err = abs(ref - leb) * 3600.0  # Convert to arcsec
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS.AYANAMSHA_ARCSEC, (
            f'Ayanamsha mode {sid_mode}: max error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestAyanamshaConsistency:
    """Verify ayanamsha consistency with sidereal positions."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("sid_mode", [0, 1, 2, 3])
    def test_sidereal_offset_leb_vs_skyfield(
        self, compare: CompareHelper, ayanamsha_dates: list[float], sid_mode: int
    ):
        """LEB and Skyfield agree on the (tropical - sidereal) offset.

        The sidereal offset = (tropical_lon - sidereal_lon) mod 360.
        LEB and Skyfield must produce identical offsets since both use
        the same selected ayanamsha definition. We compare the offsets rather
        than comparing against get_ayanamsa_ex_ut() because the
        nutation component in the ayanamsha is handled differently
        (known limitation: ~17" dpsi*cos(eps) architectural difference).
        """
        from libephemeris.constants import SUN, FLG_SPEED, FLG_SIDEREAL

        max_err = 0.0

        for jd in ayanamsha_dates[:10]:
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            # LEB: tropical - sidereal offset
            trop_leb, _ = compare.leb(ephem.calc_ut, jd, SUN, FLG_SPEED)
            sid_leb, _ = compare.leb(ephem.calc_ut, jd, SUN, FLG_SPEED | FLG_SIDEREAL)
            offset_leb = (trop_leb[0] - sid_leb[0]) % 360.0

            # Skyfield: tropical - sidereal offset
            trop_sky, _ = compare.skyfield(ephem.calc_ut, jd, SUN, FLG_SPEED)
            sid_sky, _ = compare.skyfield(
                ephem.calc_ut, jd, SUN, FLG_SPEED | FLG_SIDEREAL
            )
            offset_sky = (trop_sky[0] - sid_sky[0]) % 360.0

            err = abs(offset_leb - offset_sky) * 3600.0
            max_err = max(max_err, err)

        assert max_err < TOLS.AYANAMSHA_ARCSEC, (
            f'Mode {sid_mode}: LEB vs Skyfield sidereal offset diff = {max_err:.4f}"'
        )


class TestAdditionalAyanamshaModes:
    """Additional modes should produce identical results across backends."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("sid_mode", ADDITIONAL_SIDEREAL_MODES[:4])
    def test_additional_mode_identical(
        self, compare: CompareHelper, ayanamsha_dates: list[float], sid_mode: int
    ):
        """Catalog and compatibility ayanamshas share one implementation."""
        for jd in ayanamsha_dates[:5]:
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            ref = compare.skyfield(ephem.get_ayanamsa_ut, jd)
            leb = compare.leb(ephem.get_ayanamsa_ut, jd)

            assert ref == pytest.approx(leb, rel=1e-10), (
                f"Additional mode {sid_mode} differs at JD {jd}"
            )
