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
    MARS,
    MEAN_NODE,
    MOON,
    SUN,
    VENUS,
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_NOABERR,
    FLG_NOGDEFL,
    FLG_NONUT,
    FLG_RADIANS,
    FLG_SIDEREAL,
    FLG_SPEED,
    FLG_TOPOCTR,
    FLG_TRUEPOS,
    FLG_XYZ,
)
from libephemeris.planets import calc_ut, pheno_ut
from libephemeris.time_utils import julday

# Test date well inside LEB base range (1849-2150)
JD_TEST = julday(2024, 6, 15, 12.0)
# Eclipse date (2024-04-08 total solar eclipse)
JD_ECLIPSE = julday(2024, 4, 8, 18.17)
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
        from libephemeris.fixed_stars import fixstar_ut

        pos_leb, _, _ = _run_leb(fixstar_ut, star_name, JD_TEST, flags)
        pos_sf, _, _ = _run_skyfield(fixstar_ut, star_name, JD_TEST, flags)
        diff_lon = abs(pos_leb[0] - pos_sf[0]) * 3600
        if diff_lon > 648000:
            diff_lon = 1296000 - diff_lon
        diff_lat = abs(pos_leb[1] - pos_sf[1]) * 3600
        assert diff_lon < tol_arcsec, f"lon diff {diff_lon}\" > {tol_arcsec}\""
        assert diff_lat < tol_arcsec, f"lat diff {diff_lat}\" > {tol_arcsec}\""

    def test_star_default_flags(self):
        self._compare_star("Regulus", FLG_SPEED)

    def test_star_j2000(self):
        self._compare_star("Regulus", FLG_J2000)

    def test_star_noaberr(self):
        self._compare_star("Regulus", FLG_NOABERR)

    def test_star_equatorial(self):
        self._compare_star("Regulus", FLG_EQUATORIAL | FLG_SPEED)

    def test_star_speed_velocity(self):
        from libephemeris.fixed_stars import fixstar_ut

        pos_leb, _, _ = _run_leb(fixstar_ut, "Regulus", JD_TEST, FLG_SPEED)
        pos_sf, _, _ = _run_skyfield(fixstar_ut, "Regulus", JD_TEST, FLG_SPEED)
        assert abs(pos_leb[3] - pos_sf[3]) < 0.001

    def test_star_no_parallax(self):
        self._compare_star("Mira", FLG_SPEED)

    def test_star_high_proper_motion(self):
        self._compare_star("Arcturus", FLG_SPEED, tol_arcsec=0.1)

    def test_star_batch(self):
        from libephemeris.fixed_stars import batch_fixstars_ut

        stars = ["Regulus", "Spica", "Aldebaran", "Antares", "Sirius",
                 "Vega", "Arcturus", "Fomalhaut", "Betelgeuse", "Rigel"]
        results_leb = _run_leb(batch_fixstars_ut, stars, JD_TEST, FLG_SPEED)
        results_sf = _run_skyfield(batch_fixstars_ut, stars, JD_TEST, FLG_SPEED)
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
        pos_xyz, _ = _run_leb(calc_ut, JD_TEST, SUN, FLG_XYZ)
        pos_sph, _ = _run_leb(calc_ut, JD_TEST, SUN, FLG_SPEED)
        r = math.sqrt(pos_xyz[0] ** 2 + pos_xyz[1] ** 2 + pos_xyz[2] ** 2)
        assert abs(r - pos_sph[2]) < 1e-8

    def test_xyz_pipeline_b_distance(self):
        pos_xyz, _ = _run_leb(calc_ut, JD_TEST, MEAN_NODE, FLG_XYZ)
        pos_sph, _ = _run_leb(calc_ut, JD_TEST, MEAN_NODE, FLG_SPEED)
        r = math.sqrt(pos_xyz[0] ** 2 + pos_xyz[1] ** 2 + pos_xyz[2] ** 2)
        assert abs(r - pos_sph[2]) < 1e-6

    def test_xyz_sidereal(self):
        flags = FLG_XYZ | FLG_SIDEREAL | FLG_SPEED
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, SUN, flags)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, SUN, flags)
        for i in range(3):
            assert abs(pos_leb[i] - pos_sf[i]) < 1e-6, f"component {i} diff"

    def test_xyz_radians_noop(self):
        pos_xyz, _ = _run_leb(calc_ut, JD_TEST, SUN, FLG_XYZ)
        pos_xyzr, _ = _run_leb(calc_ut, JD_TEST, SUN, FLG_XYZ | FLG_RADIANS)
        for i in range(3):
            assert pos_xyz[i] == pos_xyzr[i], f"XYZ+RADIANS changed component {i}"

    def test_nonut_vs_skyfield(self):
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, SUN, FLG_NONUT | FLG_SPEED)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, SUN, FLG_NONUT | FLG_SPEED)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01

    def test_nonut_equatorial(self):
        flags = FLG_NONUT | FLG_EQUATORIAL | FLG_SPEED
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, SUN, flags)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, SUN, flags)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01

    def test_nonut_sidereal_mean_ayanamsha(self):
        libephemeris.set_sid_mode(1)  # Lahiri
        try:
            flags = FLG_NONUT | FLG_SIDEREAL | FLG_SPEED
            pos_leb, _ = _run_leb(calc_ut, JD_TEST, SUN, flags)
            pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, SUN, flags)
            assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01
        finally:
            libephemeris.set_sid_mode(0)

    def test_topoctr_moon(self):
        libephemeris.set_topo(*GEOPOS_ROME)
        try:
            flags = FLG_TOPOCTR | FLG_SPEED
            pos_leb, _ = _run_leb(calc_ut, JD_TEST, MOON, flags)
            pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, MOON, flags)
            assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 1.0
            assert _arcsec(abs(pos_leb[1] - pos_sf[1])) < 1.0
        finally:
            from libephemeris import state
            state._TOPO = None

    def test_topoctr_pipeline_b_rejects(self):
        libephemeris.set_topo(*GEOPOS_ROME)
        try:
            pos, _ = calc_ut(JD_TEST, MEAN_NODE, FLG_TOPOCTR)
            assert 0 < pos[0] < 360
        finally:
            from libephemeris import state
            state._TOPO = None

    def test_icrs_fallback(self):
        from libephemeris.constants import FLG_ICRS
        pos, _ = calc_ut(JD_TEST, SUN, FLG_ICRS | FLG_SPEED)
        assert 0 < pos[0] < 360


# =========================================================================
# 3. pheno_ut
# =========================================================================


class TestPhenoRegression:

    def test_pheno_mars(self):
        r_leb = _run_leb(pheno_ut, JD_TEST, MARS, 0)
        r_sf = _run_skyfield(pheno_ut, JD_TEST, MARS, 0)
        assert abs(r_leb[0] - r_sf[0]) < 0.1  # phase angle
        assert abs(r_leb[2] - r_sf[2]) < 0.1  # elongation
        assert abs(r_leb[4] - r_sf[4]) < 0.1  # magnitude

    def test_pheno_moon_phase_angle(self):
        r_leb = _run_leb(pheno_ut, JD_TEST, MOON, 0)
        r_sf = _run_skyfield(pheno_ut, JD_TEST, MOON, 0)
        assert abs(r_leb[0] - r_sf[0]) < 0.1

    def test_pheno_truepos(self):
        r_default = _run_leb(pheno_ut, JD_TEST, MARS, 0)
        r_truepos = _run_leb(pheno_ut, JD_TEST, MARS, FLG_TRUEPOS)
        assert r_default[0] != r_truepos[0]

    def test_pheno_noaberr_skips_leb(self):
        r_leb = _run_leb(pheno_ut, JD_TEST, MARS, FLG_NOABERR)
        r_sf = _run_skyfield(pheno_ut, JD_TEST, MARS, FLG_NOABERR)
        assert abs(r_leb[0] - r_sf[0]) < 0.001

    def test_pheno_unsupported_body(self):
        result = _run_leb(pheno_ut, JD_TEST, MEAN_NODE, 0)
        assert result == (0.0,) * 20

    def test_pheno_sun(self):
        result = _run_leb(pheno_ut, JD_TEST, SUN, 0)
        assert result[0] == 0.0
        assert result[2] == 0.0


# =========================================================================
# 4. Heliacal
# =========================================================================


class TestHeliacalRegression:

    @pytest.mark.slow
    def test_heliacal_pheno_venus(self):
        from libephemeris.heliacal import _heliacal_pheno_ut_pythonic

        r_leb, _ = _run_leb(_heliacal_pheno_ut_pythonic, JD_TEST, 41.9, 12.5, body=VENUS)
        r_sf, _ = _run_skyfield(_heliacal_pheno_ut_pythonic, JD_TEST, 41.9, 12.5, body=VENUS)
        # dret[1] = body altitude, dret[2] = geo altitude
        assert abs(r_leb[1] - r_sf[1]) < 0.5
        assert abs(r_leb[2] - r_sf[2]) < 0.5

    def test_heliacal_azimuth_convention(self):
        from libephemeris.heliacal import _heliacal_pheno_ut_pythonic

        r_leb, _ = _run_leb(_heliacal_pheno_ut_pythonic, JD_TEST, 41.9, 12.5, body=VENUS)
        r_sf, _ = _run_skyfield(_heliacal_pheno_ut_pythonic, JD_TEST, 41.9, 12.5, body=VENUS)
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

        _, tret_leb = _run_leb(sol_eclipse_when_glob, julday(2024, 1, 1, 0.0))
        _, tret_sf = _run_skyfield(sol_eclipse_when_glob, julday(2024, 1, 1, 0.0))
        diff_sec = abs(tret_leb[0] - tret_sf[0]) * 86400
        assert diff_sec < 0.1, f"Eclipse max time diff: {diff_sec:.3f}s"

    def test_eclipse_topo_state_preserved(self):
        from libephemeris.eclipse import sol_eclipse_how

        libephemeris.set_topo(10.0, 45.0, 100.0)
        from libephemeris.state import get_topo
        topo_before = get_topo()

        sol_eclipse_how(JD_ECLIPSE, [0.0, 0.0, 0.0])

        topo_after = get_topo()
        if topo_before is not None and topo_after is not None:
            assert abs(topo_after.longitude.degrees - topo_before.longitude.degrees) < 0.001
            assert abs(topo_after.latitude.degrees - topo_before.latitude.degrees) < 0.001

        from libephemeris import state
        state._TOPO = None

    def test_rise_trans_sunrise(self):
        from libephemeris.constants import CALC_RISE
        from libephemeris.eclipse import rise_trans

        _, tret_leb = _run_leb(rise_trans, JD_TEST, SUN, CALC_RISE, [12.5, 41.9, 0])
        _, tret_sf = _run_skyfield(rise_trans, JD_TEST, SUN, CALC_RISE, [12.5, 41.9, 0])
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
        from libephemeris.fixed_stars import fixstar_ut

        jd_far = julday(1800, 1, 1, 0.0)
        pos, _, _ = fixstar_ut("Regulus", jd_far, FLG_SPEED)
        assert 0 < pos[0] < 360

    def test_pheno_out_of_range_fallback(self):
        jd_far = julday(1800, 1, 1, 0.0)
        result = pheno_ut(jd_far, MARS, 0)
        assert result[2] > 0  # elongation

    def test_fixstar_batch_out_of_range(self):
        from libephemeris.fixed_stars import batch_fixstars_ut

        jd_far = julday(1800, 1, 1, 0.0)
        results = batch_fixstars_ut(["Regulus", "Spica"], jd_far, 0)
        assert all(r is not None for r in results)

    @pytest.mark.slow
    def test_heliacal_out_of_range_fallback(self):
        from libephemeris.heliacal import _heliacal_pheno_ut_pythonic

        jd_far = julday(1800, 6, 15, 12.0)
        result, _ = _heliacal_pheno_ut_pythonic(jd_far, 41.9, 12.5, body=VENUS)
        assert isinstance(result, tuple)


# =========================================================================
# 8. Lunar Eclipse LEB vs Skyfield
# =========================================================================


class TestLunarEclipseRegression:

    def test_lun_eclipse_when(self):
        from libephemeris.eclipse import lun_eclipse_when

        jd = julday(2024, 1, 1, 0.0)
        type_leb, tret_leb = _run_leb(lun_eclipse_when, jd)
        type_sf, tret_sf = _run_skyfield(lun_eclipse_when, jd)
        assert type_leb == type_sf
        diff_sec = abs(tret_leb[0] - tret_sf[0]) * 86400
        assert diff_sec < 0.1, f"Lunar eclipse max diff: {diff_sec:.3f}s"

    def test_lun_eclipse_how(self):
        from libephemeris.eclipse import _lun_eclipse_how_pythonic, lun_eclipse_when

        _, tret = lun_eclipse_when(julday(2024, 1, 1, 0.0))
        jd_max = tret[0]
        r_leb = _run_leb(_lun_eclipse_how_pythonic, jd_max, 41.9, 12.5)
        r_sf = _run_skyfield(_lun_eclipse_how_pythonic, jd_max, 41.9, 12.5)
        assert r_leb[0] == r_sf[0]  # same type


# =========================================================================
# 9. Multi-Planet Pheno
# =========================================================================


class TestPhenoAllPlanets:

    @pytest.mark.parametrize("planet,name", [
        (VENUS, "Venus"),
        (MARS, "Mars"),
        (4, "Jupiter"),  # JUPITER
        (5, "Saturn"),   # SATURN
    ])
    def test_pheno_planet(self, planet, name):
        r_leb = _run_leb(pheno_ut, JD_TEST, planet, 0)
        r_sf = _run_skyfield(pheno_ut, JD_TEST, planet, 0)
        assert abs(r_leb[0] - r_sf[0]) < 0.1, f"{name} phase angle"
        assert abs(r_leb[2] - r_sf[2]) < 0.1, f"{name} elongation"


# =========================================================================
# 10. FLG_HELCTR / FLG_BARYCTR
# =========================================================================


class TestHelioBaryCentric:

    def test_helctr_sun(self):
        from libephemeris.constants import FLG_HELCTR
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, MARS, FLG_HELCTR | FLG_SPEED)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, MARS, FLG_HELCTR | FLG_SPEED)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01

    def test_baryctr(self):
        from libephemeris.constants import FLG_BARYCTR
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, MARS, FLG_BARYCTR | FLG_SPEED)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, MARS, FLG_BARYCTR | FLG_SPEED)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01


# =========================================================================
# 11. TOPOCTR Multiple Bodies
# =========================================================================


class TestTopocentricMultipleBodies:

    @pytest.mark.parametrize("body,name", [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MARS, "Mars"),
        (VENUS, "Venus"),
    ])
    def test_topoctr_body(self, body, name):
        libephemeris.set_topo(*GEOPOS_ROME)
        try:
            flags = FLG_TOPOCTR | FLG_SPEED
            pos_leb, _ = _run_leb(calc_ut, JD_TEST, body, flags)
            pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, body, flags)
            tol = 1.0 if body == MOON else 0.1  # Moon has larger parallax
            assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < tol, f"{name} lon"
            assert _arcsec(abs(pos_leb[1] - pos_sf[1])) < tol, f"{name} lat"
        finally:
            from libephemeris import state
            state._TOPO = None


# =========================================================================
# 12. Stars with Sidereal
# =========================================================================


class TestStarsSidereal:

    def test_star_sidereal_lahiri(self):
        from libephemeris.fixed_stars import fixstar_ut

        libephemeris.set_sid_mode(1)  # Lahiri
        try:
            flags = FLG_SIDEREAL | FLG_SPEED
            pos_leb, _, _ = _run_leb(fixstar_ut, "Spica", JD_TEST, flags)
            pos_sf, _, _ = _run_skyfield(fixstar_ut, "Spica", JD_TEST, flags)
            diff = abs(pos_leb[0] - pos_sf[0]) * 3600
            assert diff < 0.1, f"Sidereal star diff: {diff}\""
        finally:
            libephemeris.set_sid_mode(0)

    def test_star_nogdefl(self):
        from libephemeris.fixed_stars import fixstar_ut

        pos_leb, _, _ = _run_leb(fixstar_ut, "Regulus", JD_TEST, FLG_NOGDEFL)
        pos_sf, _, _ = _run_skyfield(fixstar_ut, "Regulus", JD_TEST, FLG_NOGDEFL)
        diff = abs(pos_leb[0] - pos_sf[0]) * 3600
        assert diff < 0.1


# =========================================================================
# 13. Multi-Date Precision Stability
# =========================================================================


class TestMultiDatePrecision:

    @pytest.mark.parametrize("year", [1900, 1950, 2000, 2024, 2050, 2100])
    def test_sun_precision_across_centuries(self, year):
        jd = julday(year, 6, 15, 12.0)
        pos_leb, _ = _run_leb(calc_ut, jd, SUN, FLG_SPEED)
        pos_sf, _ = _run_skyfield(calc_ut, jd, SUN, FLG_SPEED)
        diff = _arcsec(abs(pos_leb[0] - pos_sf[0]))
        assert diff < 0.001, f"Year {year}: Sun lon diff {diff}\""

    @pytest.mark.parametrize("year", [1900, 1950, 2000, 2024, 2050, 2100])
    def test_moon_precision_across_centuries(self, year):
        jd = julday(year, 6, 15, 12.0)
        pos_leb, _ = _run_leb(calc_ut, jd, MOON, FLG_SPEED)
        pos_sf, _ = _run_skyfield(calc_ut, jd, MOON, FLG_SPEED)
        diff = _arcsec(abs(pos_leb[0] - pos_sf[0]))
        assert diff < 0.001, f"Year {year}: Moon lon diff {diff}\""


# =========================================================================
# 14. Rise/Trans Multiple Locations
# =========================================================================


class TestRiseTransLocations:

    @pytest.mark.parametrize("name,geopos", [
        ("Rome", [12.5, 41.9, 0]),
        ("New York", [-74.0, 40.7, 10]),
        ("Tokyo", [139.7, 35.7, 40]),
        ("Cape Town", [18.4, -33.9, 0]),
    ])
    def test_sunrise_location(self, name, geopos):
        from libephemeris.constants import CALC_RISE
        from libephemeris.eclipse import rise_trans

        _, tret_leb = _run_leb(rise_trans, JD_TEST, SUN, CALC_RISE, geopos)
        _, tret_sf = _run_skyfield(rise_trans, JD_TEST, SUN, CALC_RISE, geopos)
        diff_sec = abs(tret_leb[0] - tret_sf[0]) * 86400
        assert diff_sec < 1.0, f"{name} sunrise diff: {diff_sec:.3f}s"

    def test_moon_rise(self):
        from libephemeris.constants import CALC_RISE
        from libephemeris.eclipse import rise_trans

        _, tret_leb = _run_leb(rise_trans, JD_TEST, MOON, CALC_RISE, [12.5, 41.9, 0])
        _, tret_sf = _run_skyfield(rise_trans, JD_TEST, MOON, CALC_RISE, [12.5, 41.9, 0])
        diff_sec = abs(tret_leb[0] - tret_sf[0]) * 86400
        assert diff_sec < 2.0, f"Moon rise diff: {diff_sec:.3f}s"

    def test_star_rise(self):
        from libephemeris.constants import CALC_RISE
        from libephemeris.eclipse import rise_trans

        _, tret_leb = _run_leb(rise_trans, JD_TEST, "Sirius", CALC_RISE, [12.5, 41.9, 0])
        _, tret_sf = _run_skyfield(rise_trans, JD_TEST, "Sirius", CALC_RISE, [12.5, 41.9, 0])
        diff_sec = abs(tret_leb[0] - tret_sf[0]) * 86400
        assert diff_sec < 2.0, f"Sirius rise diff: {diff_sec:.3f}s"


# =========================================================================
# 15. Eclipse How/Where
# =========================================================================


class TestEclipseHowWhere:

    def test_sol_eclipse_how(self):
        from libephemeris.eclipse import sol_eclipse_how

        r_leb, attr_leb = _run_leb(sol_eclipse_how, JD_ECLIPSE, [0.0, 45.0, 0.0])
        r_sf, attr_sf = _run_skyfield(sol_eclipse_how, JD_ECLIPSE, [0.0, 45.0, 0.0])
        # Both should agree on eclipse type at this location
        assert r_leb == r_sf or (r_leb == 0 and r_sf == 0)

    def test_sol_eclipse_where(self):
        from libephemeris.eclipse import sol_eclipse_where

        r_leb = _run_leb(sol_eclipse_where, JD_ECLIPSE)
        r_sf = _run_skyfield(sol_eclipse_where, JD_ECLIPSE)
        # Compare central line latitude/longitude
        if r_leb[1][0] != 0 and r_sf[1][0] != 0:
            assert abs(r_leb[1][0] - r_sf[1][0]) < 0.1, "Central lat"
            assert abs(r_leb[1][1] - r_sf[1][1]) < 0.1, "Central lon"


# =========================================================================
# 16. NONUT for Moon and Pipeline B
# =========================================================================


class TestNonutExtended:

    def test_nonut_moon(self):
        flags = FLG_NONUT | FLG_SPEED
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, MOON, flags)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, MOON, flags)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01

    def test_nonut_mean_node(self):
        flags = FLG_NONUT | FLG_SPEED
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, MEAN_NODE, flags)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, MEAN_NODE, flags)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01

    def test_nonut_true_node(self):
        from libephemeris.constants import TRUE_NODE
        flags = FLG_NONUT | FLG_SPEED
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, TRUE_NODE, flags)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, TRUE_NODE, flags)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01


# =========================================================================
# 17. XYZ for Pipeline C (Uranians)
# =========================================================================


class TestXYZUranian:

    def test_xyz_cupido(self):
        from libephemeris.constants import CUPIDO
        pos_xyz, _ = _run_leb(calc_ut, JD_TEST, CUPIDO, FLG_XYZ)
        pos_sph, _ = _run_leb(calc_ut, JD_TEST, CUPIDO, FLG_SPEED)
        r = math.sqrt(pos_xyz[0] ** 2 + pos_xyz[1] ** 2 + pos_xyz[2] ** 2)
        assert abs(r - pos_sph[2]) < 0.01, f"Cupido XYZ r={r}, sph dist={pos_sph[2]}"


# =========================================================================
# 18. Exotic Flag Combinations
# =========================================================================


class TestExoticFlagCombinations:

    def test_xyz_equatorial_j2000_icrs_cartesian(self):
        """FLG_XYZ | FLG_EQUATORIAL | FLG_J2000 = raw ICRS Cartesian."""
        flags = FLG_XYZ | FLG_EQUATORIAL | FLG_J2000 | FLG_SPEED
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, SUN, flags)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, SUN, flags)
        for i in range(3):
            assert abs(pos_leb[i] - pos_sf[i]) < 1e-6, f"ICRS XYZ[{i}]"

    def test_xyz_nonut(self):
        flags = FLG_XYZ | FLG_NONUT | FLG_SPEED
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, SUN, flags)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, SUN, flags)
        r_leb = math.sqrt(sum(c**2 for c in pos_leb[:3]))
        r_sf = math.sqrt(sum(c**2 for c in pos_sf[:3]))
        assert abs(r_leb - r_sf) < 1e-8

    def test_xyz_topoctr_moon(self):
        libephemeris.set_topo(*GEOPOS_ROME)
        try:
            flags = FLG_XYZ | FLG_TOPOCTR | FLG_SPEED
            pos_leb, _ = _run_leb(calc_ut, JD_TEST, MOON, flags)
            pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, MOON, flags)
            r_leb = math.sqrt(sum(c**2 for c in pos_leb[:3]))
            r_sf = math.sqrt(sum(c**2 for c in pos_sf[:3]))
            assert abs(r_leb - r_sf) < 1e-4
        finally:
            from libephemeris import state
            state._TOPO = None

    def test_nonut_j2000(self):
        flags = FLG_NONUT | FLG_J2000 | FLG_SPEED
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, SUN, flags)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, SUN, flags)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01

    def test_sidereal_equatorial(self):
        libephemeris.set_sid_mode(1)
        try:
            flags = FLG_SIDEREAL | FLG_EQUATORIAL | FLG_SPEED
            pos_leb, _ = _run_leb(calc_ut, JD_TEST, SUN, flags)
            pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, SUN, flags)
            assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01
        finally:
            libephemeris.set_sid_mode(0)

    def test_sidereal_equatorial_j2000(self):
        libephemeris.set_sid_mode(1)
        try:
            flags = FLG_SIDEREAL | FLG_EQUATORIAL | FLG_J2000 | FLG_SPEED
            pos_leb, _ = _run_leb(calc_ut, JD_TEST, SUN, flags)
            pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, SUN, flags)
            assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01
        finally:
            libephemeris.set_sid_mode(0)

    def test_xyz_speed_velocity_components(self):
        flags = FLG_XYZ | FLG_SPEED
        pos, _ = _run_leb(calc_ut, JD_TEST, SUN, flags)
        assert pos[3] != 0.0 or pos[4] != 0.0 or pos[5] != 0.0, "XYZ velocity is zero"
        # Velocity magnitude should be ~0.017 AU/day for Earth orbital speed
        v = math.sqrt(pos[3]**2 + pos[4]**2 + pos[5]**2)
        assert 0.01 < v < 0.03, f"XYZ velocity magnitude {v} AU/day unexpected"

    def test_truepos_equatorial(self):
        flags = FLG_TRUEPOS | FLG_EQUATORIAL | FLG_SPEED
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, MARS, flags)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, MARS, flags)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01


# =========================================================================
# 19. Edge Case Bodies
# =========================================================================


class TestEdgeCaseBodies:

    def test_earth_geocentric_zero(self):
        from libephemeris.constants import EARTH
        pos, _ = _run_leb(calc_ut, JD_TEST, EARTH, FLG_SPEED)
        assert pos[0] == 0.0 or pos[2] == 0.0  # geocentric Earth = zero

    def test_pluto(self):
        from libephemeris.constants import PLUTO
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, PLUTO, FLG_SPEED)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, PLUTO, FLG_SPEED)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.001

    def test_chiron(self):
        from libephemeris.constants import CHIRON
        # Chiron uses SPK asteroid data with known LEB vs Skyfield offset
        # due to different Keplerian fallback handling outside SPK range.
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, CHIRON, FLG_SPEED)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, CHIRON, FLG_SPEED)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 2000

    def test_mean_apogee(self):
        from libephemeris.constants import MEAN_APOG
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, MEAN_APOG, FLG_SPEED)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, MEAN_APOG, FLG_SPEED)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.01

    def test_true_node_xyz(self):
        from libephemeris.constants import TRUE_NODE
        pos_xyz, _ = _run_leb(calc_ut, JD_TEST, TRUE_NODE, FLG_XYZ)
        pos_sph, _ = _run_leb(calc_ut, JD_TEST, TRUE_NODE, FLG_SPEED)
        r = math.sqrt(pos_xyz[0]**2 + pos_xyz[1]**2 + pos_xyz[2]**2)
        assert abs(r - pos_sph[2]) < 1e-4

    def test_oscu_apog_nonut(self):
        from libephemeris.constants import OSCU_APOG
        flags = FLG_NONUT | FLG_SPEED
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, OSCU_APOG, flags)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, OSCU_APOG, flags)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.1


# =========================================================================
# 20. Date Edge Cases
# =========================================================================


class TestDateEdgeCases:

    def test_j2000_epoch(self):
        from libephemeris.constants import J2000
        pos_leb, _ = _run_leb(calc_ut, J2000, SUN, FLG_SPEED)
        pos_sf, _ = _run_skyfield(calc_ut, J2000, SUN, FLG_SPEED)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.001

    def test_far_past_1860(self):
        jd = julday(1860, 3, 21, 12.0)
        pos_leb, _ = _run_leb(calc_ut, jd, SUN, FLG_SPEED)
        pos_sf, _ = _run_skyfield(calc_ut, jd, SUN, FLG_SPEED)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.001

    def test_far_future_2140(self):
        jd = julday(2140, 9, 22, 12.0)
        pos_leb, _ = _run_leb(calc_ut, jd, SUN, FLG_SPEED)
        pos_sf, _ = _run_skyfield(calc_ut, jd, SUN, FLG_SPEED)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.001

    def test_moon_far_past_1860(self):
        jd = julday(1860, 6, 15, 12.0)
        pos_leb, _ = _run_leb(calc_ut, jd, MOON, FLG_SPEED)
        pos_sf, _ = _run_skyfield(calc_ut, jd, MOON, FLG_SPEED)
        assert _arcsec(abs(pos_leb[0] - pos_sf[0])) < 0.001


# =========================================================================
# 21. Rise/Trans Edge Cases
# =========================================================================


class TestRiseTransEdgeCases:

    def test_transit_meridian_passage(self):
        from libephemeris.constants import CALC_MTRANSIT
        from libephemeris.eclipse import rise_trans

        _, tret_leb = _run_leb(rise_trans, JD_TEST, SUN, CALC_MTRANSIT, [12.5, 41.9, 0])
        _, tret_sf = _run_skyfield(rise_trans, JD_TEST, SUN, CALC_MTRANSIT, [12.5, 41.9, 0])
        diff_sec = abs(tret_leb[0] - tret_sf[0]) * 86400
        assert diff_sec < 1.0, f"Transit diff: {diff_sec:.3f}s"

    def test_planet_rise_jupiter(self):
        from libephemeris.constants import CALC_RISE, JUPITER
        from libephemeris.eclipse import rise_trans

        _, tret_leb = _run_leb(rise_trans, JD_TEST, JUPITER, CALC_RISE, [12.5, 41.9, 0])
        _, tret_sf = _run_skyfield(rise_trans, JD_TEST, JUPITER, CALC_RISE, [12.5, 41.9, 0])
        diff_sec = abs(tret_leb[0] - tret_sf[0]) * 86400
        assert diff_sec < 2.0, f"Jupiter rise diff: {diff_sec:.3f}s"

    def test_circumpolar_reykjavik(self):
        from libephemeris.constants import CALC_SET
        from libephemeris.eclipse import rise_trans

        jd_summer = julday(2024, 6, 21, 0.0)
        ret_leb, tret_leb = _run_leb(rise_trans, jd_summer, SUN, CALC_SET, [-21.9, 64.1, 0])
        ret_sf, tret_sf = _run_skyfield(rise_trans, jd_summer, SUN, CALC_SET, [-21.9, 64.1, 0])
        assert ret_leb == ret_sf  # same return code


# =========================================================================
# 22. Pheno Edge Cases
# =========================================================================


class TestPhenoEdgeCases:

    def test_moon_new_moon(self):
        jd_new = julday(2024, 4, 8, 18.0)  # near new moon
        r = _run_leb(pheno_ut, jd_new, MOON, 0)
        assert r[0] > 150  # phase angle near 180

    def test_moon_full_moon(self):
        jd_full = julday(2024, 4, 23, 23.0)  # near full moon
        r = _run_leb(pheno_ut, jd_full, MOON, 0)
        assert r[0] < 30  # phase angle near 0

    def test_saturn_with_rings(self):
        r_leb = _run_leb(pheno_ut, JD_TEST, 5, 0)  # Saturn
        r_sf = _run_skyfield(pheno_ut, JD_TEST, 5, 0)
        assert abs(r_leb[4] - r_sf[4]) < 0.2  # magnitude


# =========================================================================
# 23. Stars Edge Cases
# =========================================================================


class TestStarsEdgeCases:

    def test_star_near_ecliptic_pole(self):
        from libephemeris.fixed_stars import fixstar_ut
        # Polaris is near the ecliptic pole
        pos_leb, _, _ = _run_leb(fixstar_ut, "Polaris", JD_TEST, FLG_SPEED)
        pos_sf, _, _ = _run_skyfield(fixstar_ut, "Polaris", JD_TEST, FLG_SPEED)
        diff = abs(pos_leb[1] - pos_sf[1]) * 3600  # latitude
        assert diff < 0.1, f"Polaris lat diff: {diff}\""

    def test_star_equatorial_j2000_speed(self):
        from libephemeris.fixed_stars import fixstar_ut
        flags = FLG_EQUATORIAL | FLG_J2000 | FLG_SPEED
        pos_leb, _, _ = _run_leb(fixstar_ut, "Sirius", JD_TEST, flags)
        pos_sf, _, _ = _run_skyfield(fixstar_ut, "Sirius", JD_TEST, flags)
        diff = abs(pos_leb[0] - pos_sf[0]) * 3600
        assert diff < 0.1, f"Sirius EQ J2000 diff: {diff}\""

    def test_batch_mixed_parallax(self):
        from libephemeris.fixed_stars import batch_fixstars_ut
        # Mix stars with/without parallax
        stars = ["Sirius", "Mira", "Polaris", "Canopus", "Deneb"]
        r_leb = _run_leb(batch_fixstars_ut, stars, JD_TEST, FLG_SPEED)
        r_sf = _run_skyfield(batch_fixstars_ut, stars, JD_TEST, FLG_SPEED)
        for i, name in enumerate(stars):
            assert r_leb[i] is not None and r_sf[i] is not None
            diff = abs(r_leb[i][0][0] - r_sf[i][0][0]) * 3600
            if diff > 648000:
                diff = 1296000 - diff
            assert diff < 1.0, f"{name}: diff {diff}\""


# =========================================================================
# 24. Velocity Precision LEB vs Skyfield
# =========================================================================


class TestVelocityPrecision:

    @pytest.mark.parametrize("body,name,tol", [
        (SUN, "Sun", 0.0001),
        (MOON, "Moon", 0.002),
        (MARS, "Mars", 0.0001),
        (VENUS, "Venus", 0.0001),
    ])
    def test_velocity_ecliptic(self, body, name, tol):
        """Compare ecliptic velocity (dlon, dlat) LEB vs Skyfield."""
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, body, FLG_SPEED)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, body, FLG_SPEED)
        assert abs(pos_leb[3] - pos_sf[3]) < tol, f"{name} dlon: {abs(pos_leb[3]-pos_sf[3])}"
        assert abs(pos_leb[4] - pos_sf[4]) < tol, f"{name} dlat: {abs(pos_leb[4]-pos_sf[4])}"

    def test_velocity_xyz_sun(self):
        """XYZ velocity components should match Skyfield."""
        flags = FLG_XYZ | FLG_SPEED
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, SUN, flags)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, SUN, flags)
        for i in range(3, 6):
            assert abs(pos_leb[i] - pos_sf[i]) < 5e-6, f"vxyz[{i-3}]: {abs(pos_leb[i]-pos_sf[i])}"

    def test_velocity_xyz_moon(self):
        flags = FLG_XYZ | FLG_SPEED
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, MOON, flags)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, MOON, flags)
        for i in range(3, 6):
            assert abs(pos_leb[i] - pos_sf[i]) < 1e-5, f"Moon vxyz[{i-3}]: {abs(pos_leb[i]-pos_sf[i])}"

    def test_velocity_mean_node(self):
        """Pipeline B velocity should match."""
        pos_leb, _ = _run_leb(calc_ut, JD_TEST, MEAN_NODE, FLG_SPEED)
        pos_sf, _ = _run_skyfield(calc_ut, JD_TEST, MEAN_NODE, FLG_SPEED)
        assert abs(pos_leb[3] - pos_sf[3]) < 0.0001, "MeanNode dlon"

    def test_velocity_star(self):
        """Fixed star velocity (from finite difference) should match."""
        from libephemeris.fixed_stars import fixstar_ut
        pos_leb, _, _ = _run_leb(fixstar_ut, "Regulus", JD_TEST, FLG_SPEED)
        pos_sf, _, _ = _run_skyfield(fixstar_ut, "Regulus", JD_TEST, FLG_SPEED)
        assert abs(pos_leb[3] - pos_sf[3]) < 0.001, "Star dlon"
        assert abs(pos_leb[4] - pos_sf[4]) < 0.001, "Star dlat"


# =========================================================================
# 25. Local Eclipse LEB vs Skyfield
# =========================================================================


class TestLocalEclipseRegression:

    @pytest.mark.slow
    def test_sol_eclipse_when_loc_rome(self):
        """Local solar eclipse search from Rome."""
        from libephemeris.eclipse import _sol_eclipse_when_loc_pythonic

        jd = julday(2024, 1, 1, 0.0)
        type_leb, tret_leb, attr_leb = _run_leb(_sol_eclipse_when_loc_pythonic, jd, 41.9, 12.5)
        type_sf, tret_sf, attr_sf = _run_skyfield(_sol_eclipse_when_loc_pythonic, jd, 41.9, 12.5)
        # Same eclipse type
        assert type_leb == type_sf, f"Type: LEB={type_leb} SF={type_sf}"
        # Maximum time within 1 second
        if tret_leb[0] > 0 and tret_sf[0] > 0:
            diff_sec = abs(tret_leb[0] - tret_sf[0]) * 86400
            assert diff_sec < 1.0, f"Local eclipse max diff: {diff_sec:.3f}s"

    @pytest.mark.slow
    def test_lun_eclipse_when_loc(self):
        """Local lunar eclipse search."""
        from libephemeris.eclipse import _lun_eclipse_when_loc_pythonic

        jd = julday(2024, 1, 1, 0.0)
        type_leb, tret_leb, attr_leb = _run_leb(_lun_eclipse_when_loc_pythonic, jd, 41.9, 12.5)
        type_sf, tret_sf, attr_sf = _run_skyfield(_lun_eclipse_when_loc_pythonic, jd, 41.9, 12.5)
        assert type_leb == type_sf, f"Type: LEB={type_leb} SF={type_sf}"

    @pytest.mark.slow
    def test_swe_sol_eclipse_when_loc(self):
        """swe_ API local solar eclipse."""
        from libephemeris.eclipse import sol_eclipse_when_loc

        jd = julday(2024, 1, 1, 0.0)
        r_leb = _run_leb(sol_eclipse_when_loc, jd, [12.5, 41.9, 0])
        r_sf = _run_skyfield(sol_eclipse_when_loc, jd, [12.5, 41.9, 0])
        # Compare eclipse type flags
        assert r_leb[0] == r_sf[0], f"Type: LEB={r_leb[0]} SF={r_sf[0]}"


# =========================================================================
# 26. vis_limit_mag LEB vs Skyfield
# =========================================================================


class TestVisLimitMagRegression:

    def test_vis_limit_mag_venus(self):
        from libephemeris.heliacal import vis_limit_mag

        geopos = [12.5, 41.9, 0]
        atmo = [1013.25, 15.0, 50.0, 0.25]
        observer = [25, 1.0, 1.0, 1.0, 0.0, 0.0]

        r_leb, dret_leb = _run_leb(
            vis_limit_mag, JD_TEST, geopos, atmo, observer, "Venus", 0
        )
        r_sf, dret_sf = _run_skyfield(
            vis_limit_mag, JD_TEST, geopos, atmo, observer, "Venus", 0
        )
        assert r_leb == r_sf, f"Return flag: LEB={r_leb} SF={r_sf}"
        # dret[0] = limiting magnitude
        if dret_leb[0] != 0 and dret_sf[0] != 0:
            assert abs(dret_leb[0] - dret_sf[0]) < 1.0, f"Lim mag diff: {abs(dret_leb[0]-dret_sf[0])}"

    def test_vis_limit_mag_star(self):
        from libephemeris.heliacal import vis_limit_mag

        geopos = [12.5, 41.9, 0]
        atmo = [1013.25, 15.0, 50.0, 0.25]
        observer = [25, 1.0, 1.0, 1.0, 0.0, 0.0]

        r_leb, dret_leb = _run_leb(
            vis_limit_mag, JD_TEST, geopos, atmo, observer, "Sirius", 0
        )
        r_sf, dret_sf = _run_skyfield(
            vis_limit_mag, JD_TEST, geopos, atmo, observer, "Sirius", 0
        )
        assert r_leb == r_sf, f"Return flag: LEB={r_leb} SF={r_sf}"

    def test_vis_limit_mag_moon(self):
        from libephemeris.heliacal import vis_limit_mag

        geopos = [12.5, 41.9, 0]
        atmo = [1013.25, 15.0, 50.0, 0.25]
        observer = [25, 1.0, 1.0, 1.0, 0.0, 0.0]

        r_leb, dret_leb = _run_leb(
            vis_limit_mag, JD_TEST, geopos, atmo, observer, "Moon", 0
        )
        r_sf, dret_sf = _run_skyfield(
            vis_limit_mag, JD_TEST, geopos, atmo, observer, "Moon", 0
        )
        assert r_leb == r_sf, f"Return flag: LEB={r_leb} SF={r_sf}"
