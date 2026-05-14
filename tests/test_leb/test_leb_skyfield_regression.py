"""LEB vs Skyfield regression tests.

Each test computes the same value in LEB mode and Skyfield mode, then
asserts the difference is within a tight tolerance.  This catches
precision regressions introduced by future code changes.
"""

from __future__ import annotations

import math

import pytest

import libephemeris
from libephemeris.constants import (
    SE_MARS,
    SE_MEAN_NODE,
    SE_MOON,
    SE_SUN,
    SE_VENUS,
    SEFLG_EQUATORIAL,
    SEFLG_J2000,
    SEFLG_NOABERR,
    SEFLG_NOGDEFL,
    SEFLG_NONUT,
    SEFLG_RADIANS,
    SEFLG_SIDEREAL,
    SEFLG_SPEED,
    SEFLG_TOPOCTR,
    SEFLG_TRUEPOS,
    SEFLG_XYZ,
)
from libephemeris.planets import swe_calc_ut, swe_pheno_ut
from libephemeris.time_utils import swe_julday

# Test date well inside LEB base range (1849-2150)
JD_TEST = swe_julday(2024, 6, 15, 12.0)
# Eclipse date (2024-04-08 total solar eclipse)
JD_ECLIPSE = swe_julday(2024, 4, 8, 18.17)
# Rome geopos
GEOPOS_ROME = (12.5, 41.9, 0.0)


def _arcsec(deg: float) -> float:
    return deg * 3600.0


def _run_leb(func, *args, **kwargs):
    libephemeris.set_calc_mode("leb")
    try:
        return func(*args, **kwargs)
    finally:
        libephemeris.set_calc_mode("auto")


def _run_skyfield(func, *args, **kwargs):
    libephemeris.set_calc_mode("skyfield")
    try:
        return func(*args, **kwargs)
    finally:
        libephemeris.set_calc_mode("auto")


# =========================================================================
# 1. Fixed Stars
# =========================================================================


class TestFixedStarsRegression:

    def _compare_star(self, star_name, flags, tol_arcsec=0.1):
        from libephemeris.fixed_stars import swe_fixstar_ut

        pos_leb, _, _ = _run_leb(swe_fixstar_ut, star_name, JD_TEST, flags)
        pos_sf, _, _ = _run_skyfield(swe_fixstar_ut, star_name, JD_TEST, flags)
        diff_lon = abs(pos_leb[0] - pos_sf[0]) * 3600
        if diff_lon > 648000:
            diff_lon = 1296000 - diff_lon
        diff_lat = abs(pos_leb[1] - pos_sf[1]) * 3600
        assert diff_lon < tol_arcsec, f"lon diff {diff_lon}\" > {tol_arcsec}\""
        assert diff_lat < tol_arcsec, f"lat diff {diff_lat}\" > {tol_arcsec}\""

    def test_star_default_flags(self):
        self._compare_star("Regulus", SEFLG_SPEED)

    def test_star_j2000(self):
        self._compare_star("Regulus", SEFLG_J2000)

    def test_star_noaberr(self):
        self._compare_star("Regulus", SEFLG_NOABERR)

    def test_star_equatorial(self):
        self._compare_star("Regulus", SEFLG_EQUATORIAL | SEFLG_SPEED)

    def test_star_speed_velocity(self):
        from libephemeris.fixed_stars import swe_fixstar_ut

        pos_leb, _, _ = _run_leb(swe_fixstar_ut, "Regulus", JD_TEST, SEFLG_SPEED)
        pos_sf, _, _ = _run_skyfield(swe_fixstar_ut, "Regulus", JD_TEST, SEFLG_SPEED)
        assert abs(pos_leb[3] - pos_sf[3]) < 0.001

    def test_star_no_parallax(self):
        self._compare_star("Mira", SEFLG_SPEED)

    def test_star_high_proper_motion(self):
        self._compare_star("Arcturus", SEFLG_SPEED, tol_arcsec=0.1)

    def test_star_batch(self):
        from libephemeris.fixed_stars import swe_batch_fixstars_ut

        stars = ["Regulus", "Spica", "Aldebaran", "Antares", "Sirius",
                 "Vega", "Arcturus", "Fomalhaut", "Betelgeuse", "Rigel"]
        results_leb = _run_leb(swe_batch_fixstars_ut, stars, JD_TEST, SEFLG_SPEED)
        results_sf = _run_skyfield(swe_batch_fixstars_ut, stars, JD_TEST, SEFLG_SPEED)
        for i, (r_leb, r_sf) in enumerate(zip(results_leb, results_sf)):
            assert r_leb is not None and r_sf is not None
            diff = abs(r_leb[0][0] - r_sf[0][0]) * 3600
            if diff > 648000:
                diff = 1296000 - diff
            assert diff < 0.1, f"{stars[i]}: lon diff {diff}\" > 0.1\""


# =========================================================================
# 2. SEFLG Flag Combinations
# =========================================================================


class TestFlagCombinations:

    def test_xyz_pipeline_a_distance(self):
        pos_xyz, _ = _run_leb(swe_calc_ut, JD_TEST, SE_SUN, SEFLG_XYZ)
        pos_sph, _ = _run_leb(swe_calc_ut, JD_TEST, SE_SUN, SEFLG_SPEED)
        r = math.sqrt(pos_xyz[0] ** 2 + pos_xyz[1] ** 2 + pos_xyz[2] ** 2)
        assert abs(r - pos_sph[2]) < 1e-8

    def test_xyz_pipeline_b_distance(self):
        pos_xyz, _ = _run_leb(swe_calc_ut, JD_TEST, SE_MEAN_NODE, SEFLG_XYZ)
        pos_sph, _ = _run_leb(swe_calc_ut, JD_TEST, SE_MEAN_NODE, SEFLG_SPEED)
        r = math.sqrt(pos_xyz[0] ** 2 + pos_xyz[1] ** 2 + pos_xyz[2] ** 2)
        assert abs(r - pos_sph[2]) < 1e-6

    def test_xyz_sidereal(self):
        flags = SEFLG_XYZ | SEFLG_SIDEREAL | SEFLG_SPEED
        pos_leb, _ = _run_leb(swe_calc_ut, JD_TEST, SE_SUN, flags)
        pos_sf, _ = _run_skyfield(swe_calc_ut, JD_TEST, SE_SUN, flags)
        for i in range(3):
            assert abs(pos_leb[i] - pos_sf[i]) < 1e-6, f"component {i} diff"

    def test_xyz_radians_noop(self):
        pos_xyz, _ = _run_leb(swe_calc_ut, JD_TEST, SE_SUN, SEFLG_XYZ)
        pos_xyzr, _ = _run_leb(swe_calc_ut, JD_TEST, SE_SUN, SEFLG_XYZ | SEFLG_RADIANS)
        for i in range(3):
            assert pos_xyz[i] == pos_xyzr[i], f"XYZ+RADIANS changed component {i}"

    def test_nonut_vs_skyfield(self):
        pos_leb, _ = _run_leb(swe_calc_ut, JD_TEST, SE_SUN, SEFLG_NONUT | SEFLG_SPEED)
        pos_sf, _ = _run_skyfield(swe_calc_ut, JD_TEST, SE_SUN, SEFLG_NONUT | SEFLG_SPEED)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01

    def test_nonut_equatorial(self):
        flags = SEFLG_NONUT | SEFLG_EQUATORIAL | SEFLG_SPEED
        pos_leb, _ = _run_leb(swe_calc_ut, JD_TEST, SE_SUN, flags)
        pos_sf, _ = _run_skyfield(swe_calc_ut, JD_TEST, SE_SUN, flags)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01

    def test_nonut_sidereal_mean_ayanamsha(self):
        libephemeris.set_sid_mode(1)  # Lahiri
        try:
            flags = SEFLG_NONUT | SEFLG_SIDEREAL | SEFLG_SPEED
            pos_leb, _ = _run_leb(swe_calc_ut, JD_TEST, SE_SUN, flags)
            pos_sf, _ = _run_skyfield(swe_calc_ut, JD_TEST, SE_SUN, flags)
            assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01
        finally:
            libephemeris.set_sid_mode(0)

    def test_topoctr_moon(self):
        libephemeris.set_topo(*GEOPOS_ROME)
        try:
            flags = SEFLG_TOPOCTR | SEFLG_SPEED
            pos_leb, _ = _run_leb(swe_calc_ut, JD_TEST, SE_MOON, flags)
            pos_sf, _ = _run_skyfield(swe_calc_ut, JD_TEST, SE_MOON, flags)
            assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 1.0
            assert _arcsec(abs(pos_leb[1] - pos_sf[1])) < 1.0
        finally:
            from libephemeris import state
            state._TOPO = None

    def test_topoctr_pipeline_b_rejects(self):
        libephemeris.set_topo(*GEOPOS_ROME)
        try:
            pos, _ = swe_calc_ut(JD_TEST, SE_MEAN_NODE, SEFLG_TOPOCTR)
            assert 0 < pos[0] < 360
        finally:
            from libephemeris import state
            state._TOPO = None

    def test_icrs_fallback(self):
        from libephemeris.constants import SEFLG_ICRS
        pos, _ = swe_calc_ut(JD_TEST, SE_SUN, SEFLG_ICRS | SEFLG_SPEED)
        assert 0 < pos[0] < 360


# =========================================================================
# 3. swe_pheno_ut
# =========================================================================


class TestPhenoRegression:

    def test_pheno_mars(self):
        r_leb = _run_leb(swe_pheno_ut, JD_TEST, SE_MARS, 0)
        r_sf = _run_skyfield(swe_pheno_ut, JD_TEST, SE_MARS, 0)
        assert abs(r_leb[0] - r_sf[0]) < 0.1  # phase angle
        assert abs(r_leb[2] - r_sf[2]) < 0.1  # elongation
        assert abs(r_leb[4] - r_sf[4]) < 0.1  # magnitude

    def test_pheno_moon_phase_angle(self):
        r_leb = _run_leb(swe_pheno_ut, JD_TEST, SE_MOON, 0)
        r_sf = _run_skyfield(swe_pheno_ut, JD_TEST, SE_MOON, 0)
        assert abs(r_leb[0] - r_sf[0]) < 0.1

    def test_pheno_truepos(self):
        r_default = _run_leb(swe_pheno_ut, JD_TEST, SE_MARS, 0)
        r_truepos = _run_leb(swe_pheno_ut, JD_TEST, SE_MARS, SEFLG_TRUEPOS)
        assert r_default[0] != r_truepos[0]

    def test_pheno_noaberr_skips_leb(self):
        r_leb = _run_leb(swe_pheno_ut, JD_TEST, SE_MARS, SEFLG_NOABERR)
        r_sf = _run_skyfield(swe_pheno_ut, JD_TEST, SE_MARS, SEFLG_NOABERR)
        assert abs(r_leb[0] - r_sf[0]) < 0.001

    def test_pheno_unsupported_body(self):
        result = _run_leb(swe_pheno_ut, JD_TEST, SE_MEAN_NODE, 0)
        assert result == (0.0,) * 20

    def test_pheno_sun(self):
        result = _run_leb(swe_pheno_ut, JD_TEST, SE_SUN, 0)
        assert result[0] == 0.0
        assert result[2] == 0.0


# =========================================================================
# 4. Heliacal
# =========================================================================


class TestHeliacalRegression:

    @pytest.mark.slow
    def test_heliacal_pheno_venus(self):
        from libephemeris.heliacal import heliacal_pheno_ut

        r_leb, _ = _run_leb(heliacal_pheno_ut, JD_TEST, 41.9, 12.5, body=SE_VENUS)
        r_sf, _ = _run_skyfield(heliacal_pheno_ut, JD_TEST, 41.9, 12.5, body=SE_VENUS)
        # dret[1] = body altitude, dret[2] = geo altitude
        assert abs(r_leb[1] - r_sf[1]) < 0.5
        assert abs(r_leb[2] - r_sf[2]) < 0.5

    def test_heliacal_azimuth_convention(self):
        from libephemeris.heliacal import heliacal_pheno_ut

        r_leb, _ = _run_leb(heliacal_pheno_ut, JD_TEST, 41.9, 12.5, body=SE_VENUS)
        r_sf, _ = _run_skyfield(heliacal_pheno_ut, JD_TEST, 41.9, 12.5, body=SE_VENUS)
        # dret[3] = body azimuth, dret[5] = sun azimuth
        assert abs(r_leb[3] - r_sf[3]) < 5.0, "Body azimuth off by > 5°"
        assert abs(r_leb[5] - r_sf[5]) < 5.0, "Sun azimuth off by > 5°"


# =========================================================================
# 5. Eclipse
# =========================================================================


class TestEclipseRegression:

    def test_gamma_leb_vs_skyfield(self):
        from libephemeris.eclipse import _calc_gamma

        g_leb = _run_leb(_calc_gamma, JD_ECLIPSE)
        g_sf = _run_skyfield(_calc_gamma, JD_ECLIPSE)
        assert abs(g_leb - g_sf) < 1e-6

    def test_besselian_all_elements(self):
        from libephemeris.eclipse import (
            calc_besselian_d, calc_besselian_l1, calc_besselian_l2,
            calc_besselian_mu, calc_besselian_x, calc_besselian_y,
        )
        for func in [calc_besselian_x, calc_besselian_y, calc_besselian_d,
                     calc_besselian_l1, calc_besselian_l2, calc_besselian_mu]:
            v_leb = _run_leb(func, JD_ECLIPSE)
            v_sf = _run_skyfield(func, JD_ECLIPSE)
            assert abs(v_leb - v_sf) < 1e-6, f"{func.__name__}: {abs(v_leb - v_sf)}"

    def test_sol_eclipse_timing(self):
        from libephemeris.eclipse import sol_eclipse_when_glob

        _, tret_leb = _run_leb(sol_eclipse_when_glob, swe_julday(2024, 1, 1, 0.0))
        _, tret_sf = _run_skyfield(sol_eclipse_when_glob, swe_julday(2024, 1, 1, 0.0))
        diff_sec = abs(tret_leb[0] - tret_sf[0]) * 86400
        assert diff_sec < 0.1, f"Eclipse max time diff: {diff_sec:.3f}s"

    def test_eclipse_topo_state_preserved(self):
        from libephemeris.eclipse import swe_sol_eclipse_how

        libephemeris.set_topo(10.0, 45.0, 100.0)
        from libephemeris.state import get_topo
        topo_before = get_topo()

        swe_sol_eclipse_how(JD_ECLIPSE, [0.0, 0.0, 0.0])

        topo_after = get_topo()
        if topo_before is not None and topo_after is not None:
            assert abs(topo_after.longitude.degrees - topo_before.longitude.degrees) < 0.001
            assert abs(topo_after.latitude.degrees - topo_before.latitude.degrees) < 0.001

        from libephemeris import state
        state._TOPO = None

    def test_rise_trans_sunrise(self):
        from libephemeris.constants import SE_CALC_RISE
        from libephemeris.eclipse import rise_trans

        _, tret_leb = _run_leb(rise_trans, JD_TEST, SE_SUN, SE_CALC_RISE, [12.5, 41.9, 0])
        _, tret_sf = _run_skyfield(rise_trans, JD_TEST, SE_SUN, SE_CALC_RISE, [12.5, 41.9, 0])
        diff_sec = abs(tret_leb[0] - tret_sf[0]) * 86400
        assert diff_sec < 1.0, f"Sunrise diff: {diff_sec:.3f}s"


# =========================================================================
# 6. angular_separation edge cases
# =========================================================================


class TestAngularSeparation:

    def test_separation_zero(self):
        from libephemeris.utils import angular_separation
        assert angular_separation(45.0, 30.0, 45.0, 30.0) == 0.0

    def test_separation_antipodal(self):
        from libephemeris.utils import angular_separation
        assert abs(angular_separation(0.0, 90.0, 0.0, -90.0) - 180.0) < 1e-10

    def test_separation_small_angle(self):
        from libephemeris.utils import angular_separation
        tiny = 0.001 / 3600.0  # 0.001 arcsec in degrees
        sep = angular_separation(100.0, 45.0, 100.0 + tiny, 45.0)
        expected = tiny * math.cos(math.radians(45.0))
        assert abs(sep - expected) < 1e-9, f"Small angle: got {sep}, expected {expected}"


# =========================================================================
# 7. Fallback edge cases
# =========================================================================


class TestFallbackEdgeCases:

    def test_star_out_of_range_fallback(self):
        from libephemeris.fixed_stars import swe_fixstar_ut

        jd_far = swe_julday(1800, 1, 1, 0.0)
        pos, _, _ = swe_fixstar_ut("Regulus", jd_far, SEFLG_SPEED)
        assert 0 < pos[0] < 360

    def test_pheno_out_of_range_fallback(self):
        jd_far = swe_julday(1800, 1, 1, 0.0)
        result = swe_pheno_ut(jd_far, SE_MARS, 0)
        assert result[2] > 0  # elongation

    def test_fixstar_batch_out_of_range(self):
        from libephemeris.fixed_stars import swe_batch_fixstars_ut

        jd_far = swe_julday(1800, 1, 1, 0.0)
        results = swe_batch_fixstars_ut(["Regulus", "Spica"], jd_far, 0)
        assert all(r is not None for r in results)

    @pytest.mark.slow
    def test_heliacal_out_of_range_fallback(self):
        from libephemeris.heliacal import heliacal_pheno_ut

        jd_far = swe_julday(1800, 6, 15, 12.0)
        result, _ = heliacal_pheno_ut(jd_far, 41.9, 12.5, body=SE_VENUS)
        assert isinstance(result, tuple)
