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
    VENUS,
    FLG_SPEED,
    FLG_SIDEREAL,
    SIDM_GALCENT_COCHRANE,
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
                f"{body_name} star-based mode {sid_mode} position differs at JD {jd}"
            )
            # Speed slot [3]: the sidereal longitude speed must subtract the
            # TRUE ayanamsa drift (the anchor star's apparent-longitude drift,
            # annual-aberration-dominated), not the precession-polynomial rate.
            # The LEB path once used the polynomial here and was ~0.36"/day off
            # the Skyfield backend for these modes; this assertion guards it.
            assert ref[3] == pytest.approx(leb[3], abs=1e-6), (
                f"{body_name} star-based mode {sid_mode} speed differs at JD {jd}"
            )

    @pytest.mark.leb_compare
    def test_star_based_speed_no_spike_at_ayanamsha_wrap(self, compare: CompareHelper):
        """No speed spike where the star-anchored ayanamsha wraps at 0/360.

        SIDM_GALCENT_COCHRANE's ayanamsha crosses 0 deg inside the medium
        tier (~year 2225). The speed correction central-differences the
        mod-360 ayanamsha at +/-1 s, so for ~2 s around each crossing the
        two samples straddle the wrap: without shortest-arc unwrapping the
        raw -360 deg delta reported ~1.6e7 deg/day (360 deg / 2 s). The two
        backends wrap at different instants (the LEB fast path wraps the
        Skyfield-delegated anchor value BEFORE adding nutation; the Skyfield
        path wraps the final true ayanamsha), so both crossings are scanned
        on both backends. Crossing JDs bisected to sub-second precision
        against each path's own ayanamsha curve.
        """
        flags = FLG_SPEED | FLG_SIDEREAL
        wrap_jds = (2533810.013278473, 2533851.698196034)  # LEB, Skyfield
        step = 0.5 / 86400.0
        ephem.set_sid_mode(SIDM_GALCENT_COCHRANE, 2451545.0, 0.0)
        for anchor in wrap_jds:
            for k in range(-6, 7):
                jd = anchor + k * step
                ref, _ = compare.skyfield(ephem.calc_ut, jd, VENUS, flags)
                leb, _ = compare.leb(ephem.calc_ut, jd, VENUS, flags)
                assert abs(ref[3]) < 10.0, (
                    f"Skyfield speed spike {ref[3]:+.1f} deg/day at JD {jd:.8f}"
                )
                assert abs(leb[3]) < 10.0, (
                    f"LEB speed spike {leb[3]:+.1f} deg/day at JD {jd:.8f}"
                )
                assert ref[3] == pytest.approx(leb[3], abs=1e-6), (
                    f"backend speed mismatch at JD {jd:.8f}"
                )
