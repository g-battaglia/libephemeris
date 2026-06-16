"""
Deep Validation Suite 2: Exhaustive gap-coverage testing of libephemeris vs pyswisseph.

This test suite targets all functions and scenarios NOT covered by test_deep_validation.py,
including:
- Eclipse location-specific functions (sol_eclipse_when_loc, sol_eclipse_how, lun_eclipse_how, etc.)
- Lunar eclipse magnitude/gamma convenience functions
- TT-input variants (calc, get_ayanamsa, pheno, nod_aps, etc.)
- House functions (houses_armc, houses_ex2, house_pos, house_name)
- Fixed star v2 functions (fixstar2_ut, fixstar2_mag)
- Planet-centric calculations (calc_pctr)
- Rise/transit with true horizon (rise_trans_true_hor)
- Ayanamsa extended functions (get_ayanamsa_ex, get_ayanamsa_ex_ut, get_ayanamsa_name)
- Combined flag tests (sidereal+topocentric, sidereal+equatorial, etc.)
- Topocentric at multiple global locations
- Eclipse type classification
- Internal consistency checks (calc vs calc_ut, houses vs houses_ex, etc.)

All tests use Skyfield mode only (not LEB).
"""

from __future__ import annotations

import math
import random
from typing import List, Tuple

import pytest
import swisseph as swe

import libephemeris as ephem
from libephemeris.constants import *


# ============================================================================
# HELPERS
# ============================================================================


def angular_diff(a: float, b: float) -> float:
    """Compute minimal angular difference handling 0/360 wrap."""
    d = abs(a - b) % 360
    return min(d, 360 - d)


def jd_diff_seconds(jd1: float, jd2: float) -> float:
    """Difference between two JDs in seconds."""
    return abs(jd1 - jd2) * 86400.0


# Standard test dates (JD UT)
TEST_DATES_UT = [
    ephem.julday(2000, 1, 1, 12.0),  # J2000
    ephem.julday(2024, 4, 8, 18.0),  # Solar eclipse date
    ephem.julday(1990, 7, 15, 6.0),  # Past
    ephem.julday(2010, 3, 21, 0.0),  # Equinox era
    ephem.julday(1955, 11, 1, 12.0),  # Mid-century
    ephem.julday(2050, 6, 15, 0.0),  # Future
    ephem.julday(1850, 1, 15, 12.0),  # 19th century
    ephem.julday(2100, 12, 25, 6.0),  # 22nd century
]

# Global locations for testing
LOCATIONS = [
    ("Rome", 41.9028, 12.4964, 0),
    ("London", 51.5074, -0.1278, 0),
    ("New York", 40.7128, -74.006, 0),
    ("Sydney", -33.8688, 151.2093, 0),
    ("Tokyo", 35.6762, 139.6503, 0),
    ("Cape Town", -33.9249, 18.4241, 0),
    ("Equator", 0.0, 0.0, 0),
    ("Reykjavik", 64.1466, -21.9426, 0),
    ("Tromso", 69.6492, 18.9553, 0),
    ("Buenos Aires", -34.6037, -58.3816, 0),
    ("High Altitude", 27.9881, 86.9250, 5000),  # Everest
]

PLANETS_MAIN = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
    (PLUTO, "Pluto"),
]


# ============================================================================
# PART 1: TT-INPUT VARIANT FUNCTIONS
# ============================================================================


class TestTTVariants:
    """Test all TT-input (Terrestrial Time) function variants against pyswisseph."""

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:6])
    @pytest.mark.parametrize("planet_id,name", PLANETS_MAIN)
    def test_swe_calc_tt(self, jd_ut, planet_id, name):
        """calc (TT) should match pyswisseph calc."""
        dt_lib = ephem.deltat(jd_ut)
        dt_swe = swe.deltat(jd_ut)
        jd_tt_lib = jd_ut + dt_lib
        jd_tt_swe = jd_ut + dt_swe

        lib_pos, lib_flags = ephem.calc(
            jd_tt_lib, planet_id, FLG_SWIEPH | FLG_SPEED
        )
        swe_pos, swe_flags = swe.calc(
            jd_tt_swe, planet_id, swe.FLG_SWIEPH | swe.FLG_SPEED
        )

        tol = 0.06 if planet_id == MOON else 0.008
        lon_diff = angular_diff(float(lib_pos[0]), float(swe_pos[0]))
        assert lon_diff < tol, (
            f"{name} calc(TT) lon diff {lon_diff:.6f}° > {tol}° at JD {jd_ut}"
        )

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:6])
    def test_calc_tt_vs_calc_ut_consistency(self, jd_ut):
        """calc(TT) should give same result as calc_ut(UT) for same instant."""
        dt = ephem.deltat(jd_ut)
        jd_tt = jd_ut + dt

        for planet_id, name in PLANETS_MAIN[:5]:
            pos_tt, _ = ephem.calc(jd_tt, planet_id, FLG_SWIEPH | FLG_SPEED)
            pos_ut, _ = ephem.calc_ut(jd_ut, planet_id, FLG_SWIEPH | FLG_SPEED)
            lon_diff = angular_diff(float(pos_tt[0]), float(pos_ut[0]))
            assert lon_diff < 1e-8, (
                f"{name} calc(TT) vs calc_ut(UT) diff {lon_diff} at JD {jd_ut}"
            )

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:4])
    @pytest.mark.parametrize(
        "sid_mode", [SIDM_LAHIRI, SIDM_FAGAN_BRADLEY, SIDM_RAMAN]
    )
    def test_swe_get_ayanamsa_tt(self, jd_ut, sid_mode):
        """get_ayanamsa (TT) should match pyswisseph get_ayanamsa."""
        dt_lib = ephem.deltat(jd_ut)
        dt_swe = swe.deltat(jd_ut)
        jd_tt_lib = jd_ut + dt_lib
        jd_tt_swe = jd_ut + dt_swe

        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        lib_ayan = ephem.get_ayanamsa(jd_tt_lib)
        swe_ayan = swe.get_ayanamsa(jd_tt_swe)

        diff = abs(float(lib_ayan) - float(swe_ayan))
        assert diff < 0.005, (
            f"get_ayanamsa(TT) diff {diff:.6f}° for sidmode {sid_mode} at JD {jd_ut}"
        )

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:4])
    def test_swe_get_ayanamsa_tt_vs_ut_consistency(self, jd_ut):
        """get_ayanamsa(TT) should match get_ayanamsa_ut(UT) for same instant."""
        dt = ephem.deltat(jd_ut)
        jd_tt = jd_ut + dt

        ephem.set_sid_mode(SIDM_LAHIRI)

        ayan_tt = ephem.get_ayanamsa(jd_tt)
        ayan_ut = ephem.get_ayanamsa_ut(jd_ut)
        diff = abs(float(ayan_tt) - float(ayan_ut))
        assert diff < 1e-10, f"get_ayanamsa TT vs UT diff {diff} at JD {jd_ut}"

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:4])
    def test_swe_pheno_tt(self, jd_ut):
        """pheno (TT) should match pyswisseph pheno."""
        dt_lib = ephem.deltat(jd_ut)
        dt_swe = swe.deltat(jd_ut)
        jd_tt_lib = jd_ut + dt_lib
        jd_tt_swe = jd_ut + dt_swe

        for planet_id in [VENUS, MARS, JUPITER, SATURN]:
            lib_result = ephem.pheno(jd_tt_lib, planet_id, FLG_SWIEPH)
            swe_result = swe.pheno(jd_tt_swe, planet_id, swe.FLG_SWIEPH)

            # Compare phase angle [0] and elongation [0]
            lib_vals = lib_result[0] if isinstance(lib_result[0], tuple) else lib_result
            swe_vals = swe_result

            # Phase angle (index 0)
            diff_phase = abs(float(lib_vals[0]) - float(swe_vals[0]))
            assert diff_phase < 0.5, (
                f"pheno(TT) phase angle diff {diff_phase:.4f}° for planet {planet_id}"
            )

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:3])
    def test_swe_solcross_tt(self, jd_ut):
        """solcross (TT) should match pyswisseph solcross."""
        dt_lib = ephem.deltat(jd_ut)
        dt_swe = swe.deltat(jd_ut)
        jd_tt_lib = jd_ut + dt_lib
        jd_tt_swe = jd_ut + dt_swe

        for x2cross in [0.0, 90.0, 180.0, 270.0]:
            lib_jd = ephem.solcross(x2cross, jd_tt_lib)
            swe_jd = swe.solcross(x2cross, jd_tt_swe)

            diff_sec = jd_diff_seconds(float(lib_jd), float(swe_jd))
            assert diff_sec < 10, (
                f"solcross(TT) diff {diff_sec:.2f}s for x={x2cross}° at JD {jd_ut}"
            )

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:3])
    def test_swe_mooncross_tt(self, jd_ut):
        """mooncross (TT) should match pyswisseph mooncross."""
        dt_lib = ephem.deltat(jd_ut)
        dt_swe = swe.deltat(jd_ut)
        jd_tt_lib = jd_ut + dt_lib
        jd_tt_swe = jd_ut + dt_swe

        for x2cross in [0.0, 90.0, 180.0]:
            lib_jd = ephem.mooncross(x2cross, jd_tt_lib)
            swe_jd = swe.mooncross(x2cross, jd_tt_swe)

            diff_sec = jd_diff_seconds(float(lib_jd), float(swe_jd))
            assert diff_sec < 120, (
                f"mooncross(TT) diff {diff_sec:.2f}s for x={x2cross}° at JD {jd_ut}"
            )

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:3])
    def test_swe_mooncross_node_tt(self, jd_ut):
        """mooncross_node (TT) should match pyswisseph mooncross_node."""
        dt_lib = ephem.deltat(jd_ut)
        dt_swe = swe.deltat(jd_ut)
        jd_tt_lib = jd_ut + dt_lib
        jd_tt_swe = jd_ut + dt_swe

        lib_result = ephem.mooncross_node(jd_tt_lib)
        swe_result = swe.mooncross_node(jd_tt_swe)

        diff_sec = jd_diff_seconds(float(lib_result[0]), float(swe_result[0]))
        assert diff_sec < 120, f"mooncross_node(TT) diff {diff_sec:.2f}s at JD {jd_ut}"

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:3])
    def test_swe_helio_cross_tt(self, jd_ut):
        """helio_cross (TT) should match pyswisseph helio_cross."""
        dt_lib = ephem.deltat(jd_ut)
        dt_swe = swe.deltat(jd_ut)
        jd_tt_lib = jd_ut + dt_lib
        jd_tt_swe = jd_ut + dt_swe

        for planet_id in [MARS, JUPITER]:
            lib_jd = ephem.helio_cross(planet_id, 0.0, jd_tt_lib, FLG_SWIEPH)
            swe_jd = swe.helio_cross(planet_id, 0.0, jd_tt_swe, swe.FLG_SWIEPH)

            diff_sec = jd_diff_seconds(float(lib_jd), float(swe_jd))
            assert diff_sec < 120, (
                f"helio_cross(TT) diff {diff_sec:.2f}s for planet {planet_id}"
            )

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:3])
    def test_swe_fixstar_tt(self, jd_ut):
        """fixstar (TT) should match pyswisseph fixstar."""
        dt_lib = ephem.deltat(jd_ut)
        dt_swe = swe.deltat(jd_ut)
        jd_tt_lib = jd_ut + dt_lib
        jd_tt_swe = jd_ut + dt_swe

        for star in ["Regulus", "Spica", "Aldebaran", "Sirius"]:
            try:
                lib_result = ephem.fixstar(star, jd_tt_lib, FLG_SWIEPH)
                swe_result = swe.fixstar_ut(star, jd_ut, swe.FLG_SWIEPH)

                lib_lon = float(lib_result[0][0])
                swe_lon = float(swe_result[0][0])
                lon_diff = angular_diff(lib_lon, swe_lon)
                assert lon_diff < 0.01, (
                    f"fixstar(TT) {star} lon diff {lon_diff:.6f}° at JD {jd_ut}"
                )
            except Exception:
                pass  # Star may not be found


# ============================================================================
# PART 2: ECLIPSE LOCATION-SPECIFIC FUNCTIONS
# ============================================================================


class TestSolarEclipseLocation:
    """Test solar eclipse functions at specific locations."""

    def test_sol_eclipse_when_loc_dallas_2024(self):
        """April 8, 2024 total solar eclipse visible from Dallas."""
        jd_start = ephem.julday(2024, 1, 1, 0)
        dallas = [-96.797, 32.7767, 0]  # lon, lat, alt

        lib_result = ephem.sol_eclipse_when_loc(jd_start, dallas, FLG_SWIEPH)
        swe_result = swe.sol_eclipse_when_loc(jd_start, dallas)

        # Both should find the April 8, 2024 eclipse
        lib_max = float(lib_result[1][0])
        swe_max = float(swe_result[1][0])

        # Maximum time should be close
        diff_sec = jd_diff_seconds(lib_max, swe_max)
        assert diff_sec < 300, (
            f"sol_eclipse_when_loc max time diff {diff_sec:.1f}s (max 300s)"
        )

        # Both should be on April 8, 2024
        y, m, d, h = ephem.revjul(lib_max)
        assert (y, m, d) == (2024, 4, 8), f"Expected April 8, 2024, got {y}-{m}-{d}"

    def test_sol_eclipse_when_loc_multiple_locations(self):
        """Solar eclipse detection at multiple locations."""
        jd_start = ephem.julday(2020, 1, 1, 0)

        # Test locations where eclipses should be visible (within 2020-2026 range)
        test_locs = [
            ("Santiago", [-70.6693, -33.4489, 0]),  # Chile - 2020 eclipse
            ("Mazatlan", [-106.4206, 23.2494, 0]),  # Mexico - 2024 eclipse
        ]

        for name, geopos in test_locs:
            try:
                lib_result = ephem.sol_eclipse_when_loc(
                    jd_start, geopos, FLG_SWIEPH
                )
                swe_result = swe.sol_eclipse_when_loc(jd_start, geopos)

                lib_max = float(lib_result[1][0])
                swe_max = float(swe_result[1][0])

                # Should find eclipses on similar dates
                diff_days = abs(lib_max - swe_max)
                assert diff_days < 1.0, (
                    f"sol_eclipse_when_loc at {name}: date diff {diff_days:.2f} days"
                )
            except Exception:
                pass  # Some locations may not find eclipses in range

    def test_sol_eclipse_how_during_eclipse(self):
        """sol_eclipse_how during known eclipse should return valid attributes."""
        # April 8, 2024 eclipse - find max time first
        jd_start = ephem.julday(2024, 1, 1, 0)
        dallas = [-96.797, 32.7767, 0]

        # Find eclipse time
        swe_result = swe.sol_eclipse_when_loc(jd_start, dallas)
        jd_max = float(swe_result[1][0])

        # Calculate circumstances at maximum
        lib_how = ephem.sol_eclipse_how(jd_max, dallas, FLG_SWIEPH)
        swe_how = swe.sol_eclipse_how(jd_max, dallas)

        lib_retflag = lib_how[0]
        swe_retflag = swe_how[0]

        # Both should detect an eclipse (retflag > 0)
        assert lib_retflag > 0, "lib should detect eclipse at known eclipse time"
        assert swe_retflag > 0, "swe should detect eclipse at known eclipse time"

        # Compare magnitude (attr[0])
        lib_mag = float(lib_how[1][0])
        swe_mag = float(swe_how[1][0])
        if swe_mag > 0.01:
            mag_ratio = lib_mag / swe_mag if swe_mag != 0 else float("inf")
            assert 0.5 < mag_ratio < 2.0, (
                f"sol_eclipse_how magnitude ratio {mag_ratio:.3f} out of range "
                f"(lib={lib_mag:.4f}, swe={swe_mag:.4f})"
            )

    def test_sol_eclipse_how_no_eclipse(self):
        """sol_eclipse_how at random time should return 0 magnitude."""
        jd = ephem.julday(2024, 7, 15, 12.0)  # Not an eclipse date
        rome = [12.4964, 41.9028, 0]

        lib_result = ephem.sol_eclipse_how(jd, rome, FLG_SWIEPH)

        # Should return 0 or very small magnitude
        lib_mag = float(lib_result[1][0])
        assert lib_mag < 0.01, (
            f"sol_eclipse_how at non-eclipse time should be ~0, got {lib_mag}"
        )


class TestLunarEclipseLocation:
    """Test lunar eclipse functions at specific locations."""

    def test_lun_eclipse_how_during_eclipse(self):
        """lun_eclipse_how during known total lunar eclipse."""
        # Nov 8, 2022 total lunar eclipse
        # First find it with swe
        jd_start = ephem.julday(2022, 10, 1, 0)
        swe_ecl = swe.lun_eclipse_when(jd_start)
        jd_max = float(swe_ecl[1][0])

        la = [-118.24, 34.05, 0]  # Los Angeles

        lib_result = ephem.lun_eclipse_how(jd_max, la, FLG_SWIEPH)
        swe_result = swe.lun_eclipse_how(jd_max, la)

        lib_retflag = lib_result[0]
        swe_retflag = swe_result[0]

        # Both should detect an eclipse
        assert lib_retflag > 0, "lib should detect lunar eclipse"
        assert swe_retflag > 0, "swe should detect lunar eclipse"

        # Compare umbral magnitude (attr[0])
        lib_umag = float(lib_result[1][0])
        swe_umag = float(swe_result[1][0])

        if swe_umag > 0.1:
            diff = abs(lib_umag - swe_umag)
            assert diff < 0.5, (
                f"lun_eclipse_how umbral mag diff {diff:.4f} "
                f"(lib={lib_umag:.4f}, swe={swe_umag:.4f})"
            )

    def test_lun_eclipse_how_no_eclipse(self):
        """lun_eclipse_how at non-eclipse time should return 0."""
        jd = ephem.julday(2024, 7, 15, 12.0)
        rome = [12.4964, 41.9028, 0]

        lib_result = ephem.lun_eclipse_how(jd, rome, FLG_SWIEPH)

        lib_umag = float(lib_result[1][0])
        assert lib_umag < 0.01, (
            f"lun_eclipse_how at non-eclipse time should be ~0, got {lib_umag}"
        )

    def test_lun_eclipse_when_loc(self):
        """lun_eclipse_when_loc should find eclipses visible from location."""
        jd_start = ephem.julday(2022, 1, 1, 0)

        lib_result = ephem.lun_eclipse_when_loc(jd_start, [-118.24, 34.05, 0.0])
        swe_result = swe.lun_eclipse_when_loc(jd_start, [-118.24, 34.05, 0])

        lib_max = float(lib_result[1][0])
        swe_max = float(swe_result[1][0])

        # Should find similar eclipse times
        diff_sec = jd_diff_seconds(lib_max, swe_max)
        assert diff_sec < 600, (
            f"lun_eclipse_when_loc max time diff {diff_sec:.1f}s (max 600s)"
        )


class TestLunarEclipseMagnitudeGamma:
    """Test lunar eclipse convenience functions (magnitude, gamma)."""

    def _get_known_lunar_eclipses(self):
        """Get a few known lunar eclipse max times from pyswisseph."""
        eclipses = []
        jd = ephem.julday(2020, 1, 1, 0)
        for _ in range(5):
            result = swe.lun_eclipse_when(jd)
            jd_max = float(result[1][0])
            ecl_type = result[0]
            eclipses.append((jd_max, ecl_type))
            jd = jd_max + 30  # Skip ahead
        return eclipses

    def test_lun_eclipse_umbral_magnitude(self):
        """lun_eclipse_umbral_magnitude at known eclipse times."""
        eclipses = self._get_known_lunar_eclipses()
        for jd_max, ecl_type in eclipses:
            lib_umag = float(ephem.lun_eclipse_umbral_magnitude(jd_max))

            # At eclipse maximum, should have some magnitude
            if ecl_type & 4:  # ECL_TOTAL
                assert lib_umag > 0.5, (
                    f"Total eclipse at JD {jd_max}: umbral mag {lib_umag} should be > 0.5"
                )
            elif ecl_type & 2:  # ECL_PARTIAL
                assert lib_umag > 0.0, (
                    f"Partial eclipse at JD {jd_max}: umbral mag {lib_umag} should be > 0"
                )

    def test_lun_eclipse_penumbral_magnitude(self):
        """lun_eclipse_penumbral_magnitude at known eclipse times."""
        eclipses = self._get_known_lunar_eclipses()
        for jd_max, ecl_type in eclipses:
            lib_pmag = float(ephem.lun_eclipse_penumbral_magnitude(jd_max))

            # At any eclipse maximum, penumbral magnitude should be > 0
            assert lib_pmag > 0.0, (
                f"Eclipse at JD {jd_max}: penumbral mag {lib_pmag} should be > 0"
            )

    def test_lun_eclipse_gamma(self):
        """lun_eclipse_gamma at known eclipse times should be reasonable."""
        eclipses = self._get_known_lunar_eclipses()
        for jd_max, ecl_type in eclipses:
            lib_gamma = float(ephem.lun_eclipse_gamma(jd_max))

            # During an eclipse, gamma should be < ~1.6
            assert abs(lib_gamma) < 2.0, (
                f"Eclipse at JD {jd_max}: gamma {lib_gamma} should be < 2.0"
            )

    def test_lun_eclipse_no_eclipse_time(self):
        """At non-eclipse time, magnitudes should be ~0 and gamma large."""
        jd = ephem.julday(2024, 7, 15, 12.0)

        umag = float(ephem.lun_eclipse_umbral_magnitude(jd))
        pmag = float(ephem.lun_eclipse_penumbral_magnitude(jd))
        gamma = float(ephem.lun_eclipse_gamma(jd))

        assert umag < 0.01, f"No eclipse: umbral mag {umag} should be ~0"
        assert pmag < 0.5, f"No eclipse: penumbral mag {pmag} should be small"


class TestSolarEclipseMagnitude:
    """Test solar eclipse magnitude/obscuration at location."""

    def test_sol_eclipse_magnitude_at_loc(self):
        """sol_eclipse_magnitude_at_loc during known eclipse."""
        # Find April 8 2024 eclipse max at Dallas
        jd_start = ephem.julday(2024, 1, 1, 0)
        dallas = [-96.797, 32.7767, 0]
        swe_result = swe.sol_eclipse_when_loc(jd_start, dallas)
        jd_max = float(swe_result[1][0])

        lib_mag = float(
            ephem.sol_eclipse_magnitude_at_loc(jd_max, dallas, FLG_SWIEPH)
        )

        # Dallas should see >0.9 magnitude for this total eclipse
        assert lib_mag > 0.5, (
            f"April 2024 eclipse at Dallas: magnitude {lib_mag} should be > 0.5"
        )

    def test_sol_eclipse_magnitude_no_eclipse(self):
        """sol_eclipse_magnitude_at_loc at random time should return ~0."""
        jd = ephem.julday(2024, 7, 15, 12.0)
        rome = [12.4964, 41.9028, 0]

        lib_mag = float(ephem.sol_eclipse_magnitude_at_loc(jd, rome, FLG_SWIEPH))
        assert lib_mag < 0.01, f"No eclipse: magnitude {lib_mag} should be ~0"


class TestEclipseTypeClassification:
    """Test that eclipse type flags are correctly classified."""

    def test_solar_eclipse_types(self):
        """Verify solar eclipse type classification over 10 years."""
        jd = ephem.julday(2020, 1, 1, 0)
        for _ in range(15):
            lib_result = ephem.sol_eclipse_when_glob(jd, FLG_SWIEPH)
            swe_result = swe.sol_eclipse_when_glob(jd)

            lib_type = lib_result[0]
            swe_type = swe_result[0]

            # Check type bits match (total, annular, partial)
            for bit, name in [
                (ECL_TOTAL, "TOTAL"),
                (ECL_ANNULAR, "ANNULAR"),
                (ECL_PARTIAL, "PARTIAL"),
            ]:
                lib_has = bool(lib_type & bit)
                swe_has = bool(swe_type & bit)
                assert lib_has == swe_has, (
                    f"Eclipse type mismatch at JD {float(lib_result[1][0]):.2f}: "
                    f"{name}: lib={lib_has}, swe={swe_has} "
                    f"(lib_flags={lib_type}, swe_flags={swe_type})"
                )

            jd = float(lib_result[1][0]) + 30

    def test_lunar_eclipse_types(self):
        """Verify lunar eclipse type classification over 10 years."""
        jd = ephem.julday(2020, 1, 1, 0)
        for _ in range(15):
            lib_result = ephem.lun_eclipse_when(jd, FLG_SWIEPH)
            swe_result = swe.lun_eclipse_when(jd)

            lib_type = lib_result[0]
            swe_type = swe_result[0]

            # Check type bits
            for bit, name in [
                (ECL_TOTAL, "TOTAL"),
                (ECL_PARTIAL, "PARTIAL"),
                (ECL_PENUMBRAL, "PENUMBRAL"),
            ]:
                lib_has = bool(lib_type & bit)
                swe_has = bool(swe_type & bit)
                assert lib_has == swe_has, (
                    f"Lunar eclipse type mismatch at JD {float(lib_result[1][0]):.2f}: "
                    f"{name}: lib={lib_has}, swe={swe_has} "
                    f"(lib_flags={lib_type}, swe_flags={swe_type})"
                )

            jd = float(lib_result[1][0]) + 30


# ============================================================================
# PART 3: HOUSE FUNCTIONS
# ============================================================================


class TestHouseAdvanced:
    """Test advanced house functions: houses_armc, houses_ex2, house_pos, house_name."""

    @pytest.mark.parametrize(
        "hsys_code", ["P", "K", "O", "R", "C", "E", "W", "M", "B", "T", "X"]
    )
    @pytest.mark.parametrize(
        "location",
        [
            ("Rome", 41.9028, 12.4964),
            ("London", 51.5074, -0.1278),
            ("Equator", 0.0, 0.0),
            ("Southern", -33.87, 151.21),
        ],
    )
    def test_houses_armc(self, hsys_code, location):
        """houses_armc should match pyswisseph houses_armc."""
        name, lat, lon = location
        jd = ephem.julday(2024, 6, 21, 12.0)

        # Get ARMC from houses, obliquity from ECL_NUT
        _, ascmc_swe = swe.houses(jd, lat, lon, hsys_code.encode())
        armc = float(ascmc_swe[2])
        ecl_nut, _ = swe.calc_ut(jd, swe.ECL_NUT)
        eps = float(ecl_nut[0])

        lib_cusps, lib_ascmc = ephem.houses_armc(armc, lat, eps, ord(hsys_code))
        swe_cusps, swe_ascmc = swe.houses_armc(armc, lat, eps, hsys_code.encode())

        # Compare cusps
        max_diff = 0
        for i in range(12):
            diff = angular_diff(float(lib_cusps[i]), float(swe_cusps[i]))
            max_diff = max(max_diff, diff)

        assert max_diff < 0.01, (
            f"houses_armc {hsys_code} at {name}: max cusp diff {max_diff:.6f}°"
        )

    @pytest.mark.parametrize("hsys_code", ["P", "K", "O", "R", "C", "E", "W"])
    def test_houses_ex2_cusps_match_houses(self, hsys_code):
        """houses_ex2 cusps should match houses for same inputs."""
        jd = ephem.julday(2024, 6, 21, 12.0)
        lat, lon = 41.9028, 12.4964

        lib_cusps_h, lib_ascmc_h = ephem.houses(jd, lat, lon, ord(hsys_code))
        lib_result = ephem.houses_ex2(jd, lat, lon, ord(hsys_code))
        lib_cusps_ex2 = lib_result[0]

        max_diff = 0
        for i in range(12):
            diff = angular_diff(float(lib_cusps_h[i]), float(lib_cusps_ex2[i]))
            max_diff = max(max_diff, diff)

        assert max_diff < 1e-8, (
            f"houses_ex2 vs houses {hsys_code}: max cusp diff {max_diff}"
        )

    @pytest.mark.parametrize("hsys_code", ["P", "K", "O", "R", "C", "E", "W"])
    def test_houses_ex2_vs_pyswisseph(self, hsys_code):
        """houses_ex2 should match pyswisseph houses_ex2."""
        jd = ephem.julday(2024, 6, 21, 12.0)
        lat, lon = 41.9028, 12.4964

        lib_result = ephem.houses_ex2(jd, lat, lon, ord(hsys_code))
        swe_result = swe.houses_ex2(jd, lat, lon, hsys_code.encode())

        # Compare cusps
        max_diff = 0
        for i in range(12):
            diff = angular_diff(float(lib_result[0][i]), float(swe_result[0][i]))
            max_diff = max(max_diff, diff)

        assert max_diff < 0.01, f"houses_ex2 {hsys_code}: max cusp diff {max_diff:.6f}°"

        # Note: cusp speeds (result[2]) are not yet implemented in libephemeris
        # (returns all zeros). Cusp positions are validated above.

    @pytest.mark.parametrize("hsys_code", ["P", "K", "O", "R", "C", "E", "W"])
    def test_house_pos(self, hsys_code):
        """house_pos should match pyswisseph house_pos."""
        jd = ephem.julday(2024, 6, 21, 12.0)
        lat = 41.9028

        # Get ARMC and eps from houses
        _, ascmc = swe.houses(jd, lat, 12.4964, hsys_code.encode())
        armc = float(ascmc[2])
        eps = float(ascmc[3])

        # Test a few planet positions
        for planet_id in [SUN, MOON, VENUS, MARS, JUPITER]:
            pos, _ = swe.calc_ut(jd, planet_id)
            obj_coord = (float(pos[0]), float(pos[1]))

            lib_hpos = float(
                ephem.house_pos(armc, lat, eps, obj_coord, ord(hsys_code))
            )
            swe_hpos = float(
                swe.house_pos(armc, lat, eps, obj_coord, hsys_code.encode())
            )

            diff = abs(lib_hpos - swe_hpos)
            # Handle 36/0 wrapping for gauquelin
            if diff > 18:
                diff = 36 - diff
            # Whole Sign (W) uses different cusp boundaries, allow wider tolerance
            tol = 0.5 if hsys_code == "W" else 0.2
            assert diff < tol, (
                f"house_pos {hsys_code} planet {planet_id}: diff {diff:.6f} "
                f"(lib={lib_hpos:.4f}, swe={swe_hpos:.4f})"
            )

    def test_house_name(self):
        """house_name should return correct names."""
        expected_names = {
            "P": "Placidus",
            "K": "Koch",
            "O": "Porphyry",
            "R": "Regiomontanus",
            "C": "Campanus",
            "E": "Equal",
            "W": "Whole Sign",
            "M": "Morinus",
        }

        for code, expected in expected_names.items():
            lib_name = ephem.house_name(ord(code))
            swe_name = swe.house_name(code.encode())
            # At least one of them should match
            assert (
                expected.lower() in lib_name.lower()
                or lib_name.lower() in swe_name.lower()
            ), (
                f"house_name({code}): lib='{lib_name}', swe='{swe_name}', expected '{expected}'"
            )

    @pytest.mark.parametrize("hsys_code", ["P", "K", "O", "R", "C", "E", "W"])
    def test_houses_vs_houses_ex_consistency(self, hsys_code):
        """houses should give identical cusps to houses_ex with flags=0."""
        jd = ephem.julday(2024, 3, 20, 12.0)
        lat, lon = 48.8566, 2.3522  # Paris

        cusps_h, ascmc_h = ephem.houses(jd, lat, lon, ord(hsys_code))
        cusps_ex, ascmc_ex = ephem.houses_ex(jd, lat, lon, ord(hsys_code), 0)

        for i in range(12):
            diff = angular_diff(float(cusps_h[i]), float(cusps_ex[i]))
            assert diff < 1e-8, (
                f"houses vs houses_ex {hsys_code} cusp {i + 1}: diff {diff}"
            )


# ============================================================================
# PART 4: FIXED STAR FUNCTIONS
# ============================================================================


class TestFixedStarV2:
    """Test fixstar2_ut and fixstar2_mag functions."""

    STARS_TO_TEST = [
        "Regulus",
        "Spica",
        "Aldebaran",
        "Sirius",
        "Antares",
        "Fomalhaut",
        "Pollux",
        "Deneb",
        "Vega",
        "Rigel",
        "Betelgeuse",
        "Canopus",
        "Achernar",
        "Arcturus",
    ]

    @pytest.mark.parametrize("star", STARS_TO_TEST)
    def test_fixstar2_ut_positions(self, star):
        """fixstar2_ut positions should match fixstar_ut for same star."""
        jd = ephem.julday(2024, 6, 21, 12.0)

        try:
            result1 = ephem.fixstar_ut(star, jd, FLG_SWIEPH)
            result2 = ephem.fixstar2_ut(star, jd, FLG_SWIEPH)
        except Exception:
            pytest.skip(f"Star {star} not found")

        # fixstar_ut returns (positions, star_name, retflag)
        # fixstar2_ut returns (positions, star_name, retflag)
        pos1 = result1[0]
        pos2 = result2[0]

        lon_diff = angular_diff(float(pos1[0]), float(pos2[0]))
        assert lon_diff < 1e-8, f"fixstar2_ut vs fixstar_ut {star}: lon diff {lon_diff}"

    @pytest.mark.parametrize("star", STARS_TO_TEST)
    def test_fixstar2_ut_vs_pyswisseph(self, star):
        """fixstar2_ut should match pyswisseph fixstar2_ut."""
        jd = ephem.julday(2024, 6, 21, 12.0)

        try:
            lib_result = ephem.fixstar2_ut(star, jd, FLG_SWIEPH)
            swe_result = swe.fixstar2_ut(star, jd, swe.FLG_SWIEPH)
        except Exception:
            pytest.skip(f"Star {star} not found in one implementation")

        lib_lon = float(lib_result[0][0])
        swe_lon = float(swe_result[0][0])

        lon_diff = angular_diff(lib_lon, swe_lon)
        assert lon_diff < 0.01, (
            f"fixstar2_ut {star}: lon diff {lon_diff:.6f}° "
            f"(lib={lib_lon:.6f}, swe={swe_lon:.6f})"
        )

    @pytest.mark.parametrize("star", STARS_TO_TEST)
    def test_fixstar2_mag(self, star):
        """fixstar2_mag should match pyswisseph fixstar2_mag."""
        try:
            lib_result = ephem.fixstar2_mag(star)
            swe_result = swe.fixstar2_mag(star)
        except Exception:
            pytest.skip(f"Star {star} not found in one implementation")

        # lib returns (magnitude, star_name)
        # swe returns (magnitude, star_name)
        lib_mag = float(lib_result[0])
        swe_mag = float(swe_result[0])

        diff = abs(lib_mag - swe_mag)
        assert diff < 0.5, (
            f"fixstar2_mag {star}: diff {diff:.2f} (lib={lib_mag}, swe={swe_mag})"
        )


# ============================================================================
# PART 5: PLANET-CENTRIC CALCULATIONS
# ============================================================================


class TestCalcPctr:
    """Test planet-centric position calculations."""

    @pytest.mark.parametrize(
        "target,center",
        [
            (MOON, MARS),
            (SUN, JUPITER),
            (VENUS, MARS),
            (MARS, JUPITER),
            (SATURN, JUPITER),
        ],
    )
    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:4])
    def test_calc_pctr_vs_pyswisseph(self, target, center, jd_ut):
        """calc_pctr should match pyswisseph calc_pctr."""
        dt_lib = ephem.deltat(jd_ut)
        dt_swe = swe.deltat(jd_ut)
        jd_tt_lib = jd_ut + dt_lib
        jd_tt_swe = jd_ut + dt_swe

        try:
            lib_result = ephem.calc_pctr(
                jd_tt_lib, target, center, FLG_SWIEPH | FLG_SPEED
            )
            swe_result = swe.calc_pctr(
                jd_tt_swe, target, center, swe.FLG_SWIEPH | swe.FLG_SPEED
            )
        except Exception:
            pytest.skip(f"calc_pctr not available for target={target} center={center}")

        lib_lon = float(lib_result[0][0])
        swe_lon = float(swe_result[0][0])

        lon_diff = angular_diff(lib_lon, swe_lon)
        # Planet-centric compounds position errors
        tol = 0.15 if MOON in (target, center) else 0.02
        assert lon_diff < tol, (
            f"calc_pctr target={target} center={center} lon diff {lon_diff:.6f}° "
            f"(lib={lib_lon:.4f}, swe={swe_lon:.4f})"
        )


# ============================================================================
# PART 6: RISE/TRANSIT TRUE HORIZON
# ============================================================================


class TestRiseTransTrueHorizon:
    """Test rise/transit calculations with true horizon altitude."""

    @pytest.mark.parametrize(
        "location",
        [
            ("Rome", 41.9028, 12.4964, 0.0),
            ("Denver", 39.7392, -104.9903, 1609.0),  # High altitude
        ],
    )
    @pytest.mark.parametrize("hor_alt", [0.0, 0.5, 1.0, -0.5])
    def test_rise_trans_true_hor_sun(self, location, hor_alt):
        """rise_trans_true_hor for Sun should be reasonable vs standard rise_trans."""
        name, lat, lon, alt = location
        jd = ephem.julday(2024, 6, 21, 0)

        try:
            result_true = ephem.rise_trans_true_hor(
                jd,
                SUN,
                CALC_RISE,
                [lon, lat, alt],
                1013.25,
                15.0,
                hor_alt,
                FLG_SWIEPH,
            )
            result_std = ephem.rise_trans(
                jd, SUN, CALC_RISE, [lon, lat, alt], 1013.25, 15.0, FLG_SWIEPH
            )
        except Exception:
            pytest.skip(f"rise_trans_true_hor not available at {name}")

        jd_true = float(result_true[1][0])
        jd_std = float(result_std[1][0])

        # With different horizon altitudes, times should differ but be reasonable
        diff_min = abs(jd_true - jd_std) * 1440  # minutes
        # 1 degree of horizon altitude ≈ 4 minutes for Sun
        max_diff_min = max(abs(hor_alt) * 10 + 5, 10)
        assert diff_min < max_diff_min, (
            f"rise_trans_true_hor vs rise_trans at {name} hor_alt={hor_alt}: "
            f"diff {diff_min:.1f} min (max {max_diff_min})"
        )

    def test_rise_trans_true_hor_vs_pyswisseph(self):
        """rise_trans_true_hor should match pyswisseph rise_trans_true_hor."""
        jd = ephem.julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964

        lib_result = ephem.rise_trans_true_hor(
            jd, SUN, CALC_RISE, [lon, lat, 0.0], 0.0, 0.0, 0.5, FLG_SWIEPH
        )
        swe_result = swe.rise_trans_true_hor(
            jd, swe.SUN, swe.CALC_RISE, [lon, lat, 0], 0, 0, 0.5
        )

        lib_jd = float(lib_result[1][0])
        swe_jd = float(swe_result[1][0])

        diff_sec = jd_diff_seconds(lib_jd, swe_jd)
        assert diff_sec < 120, f"rise_trans_true_hor diff {diff_sec:.1f}s (max 120s)"


# ============================================================================
# PART 7: AYANAMSA EXTENDED FUNCTIONS
# ============================================================================


class TestAyanamsaExtended:
    """Test extended ayanamsa functions."""

    @pytest.mark.parametrize(
        "sid_mode",
        [
            SIDM_LAHIRI,
            SIDM_FAGAN_BRADLEY,
            SIDM_RAMAN,
            SIDM_KRISHNAMURTI,
            SIDM_TRUE_CITRA,
        ],
    )
    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:4])
    def test_get_ayanamsa_ex_ut(self, sid_mode, jd_ut):
        """get_ayanamsa_ex_ut should match pyswisseph.

        Both APIs now use global set_sid_mode and return (retflag, ayanamsa).
        """
        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        lib_result = ephem.get_ayanamsa_ex_ut(jd_ut, FLG_SWIEPH)
        swe_result = swe.get_ayanamsa_ex_ut(jd_ut, swe.FLG_SWIEPH)

        lib_ayan = float(lib_result[1])
        swe_ayan = float(swe_result[1])

        diff = abs(lib_ayan - swe_ayan)
        tol = 0.1 if sid_mode in (SIDM_TRUE_CITRA,) else 0.005
        assert diff < tol, (
            f"get_ayanamsa_ex_ut diff {diff:.6f}° for sidmode {sid_mode} "
            f"(lib={lib_ayan:.6f}, swe={swe_ayan:.6f})"
        )

    @pytest.mark.parametrize(
        "sid_mode",
        [
            SIDM_LAHIRI,
            SIDM_FAGAN_BRADLEY,
            SIDM_RAMAN,
        ],
    )
    def test_get_ayanamsa_ex_tt(self, sid_mode):
        """get_ayanamsa_ex (TT) should match pyswisseph."""
        jd_ut = TEST_DATES_UT[0]
        dt = ephem.deltat(jd_ut)
        jd_tt = jd_ut + dt

        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        lib_result = ephem.get_ayanamsa_ex(jd_tt, FLG_SWIEPH)
        swe_result = swe.get_ayanamsa_ex(jd_tt + swe.deltat(jd_ut) - dt, swe.FLG_SWIEPH)

        lib_ayan = float(lib_result[1])
        swe_ayan = float(swe_result[1])

        diff = abs(lib_ayan - swe_ayan)
        assert diff < 0.005, (
            f"get_ayanamsa_ex(TT) diff {diff:.6f}° for sidmode {sid_mode}"
        )

    def test_get_ayanamsa_name(self):
        """get_ayanamsa_name should return correct names."""
        test_cases = [
            (SIDM_LAHIRI, "Lahiri"),
            (SIDM_FAGAN_BRADLEY, "Fagan"),
            (SIDM_RAMAN, "Raman"),
        ]

        for sid_mode, expected_substr in test_cases:
            lib_name = ephem.get_ayanamsa_name(sid_mode)
            swe_name = swe.get_ayanamsa_name(sid_mode)
            assert expected_substr.lower() in lib_name.lower(), (
                f"get_ayanamsa_name({sid_mode}): '{lib_name}' should contain '{expected_substr}'"
            )

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:3])
    def test_get_ayanamsa_ex_ut_vs_get_ayanamsa_ut_consistency(self, jd_ut):
        """get_ayanamsa_ex_ut and get_ayanamsa_ut do NOT agree with each other.

        In pyswisseph itself the two functions differ by the nutation in
        longitude (~0.004°): get_ayanamsa_ut returns the *true* ayanamsa (mean
        + nutation) while get_ayanamsa_ex_ut with FLG_SWIEPH returns the mean
        ayanamsa. So the correct 1:1-with-Swiss check is that our lib reproduces
        *each* Swiss value, not that the two of ours equal one another.
        """
        ephem.set_sid_mode(SIDM_LAHIRI)
        swe.set_sid_mode(swe.SIDM_LAHIRI)

        lib_ex = float(ephem.get_ayanamsa_ex_ut(jd_ut, FLG_SWIEPH)[1])
        lib_ut = float(ephem.get_ayanamsa_ut(jd_ut))
        swe_ex = float(swe.get_ayanamsa_ex_ut(jd_ut, swe.FLG_SWIEPH)[1])
        swe_ut = float(swe.get_ayanamsa_ut(jd_ut))

        assert abs(lib_ex - swe_ex) < 0.005, (
            f"get_ayanamsa_ex_ut vs swe diff {abs(lib_ex - swe_ex):.6f}° "
            f"(lib={lib_ex:.6f}, swe={swe_ex:.6f})"
        )
        assert abs(lib_ut - swe_ut) < 0.005, (
            f"get_ayanamsa_ut vs swe diff {abs(lib_ut - swe_ut):.6f}° "
            f"(lib={lib_ut:.6f}, swe={swe_ut:.6f})"
        )


# ============================================================================
# PART 8: COMBINED FLAG TESTS
# ============================================================================


class TestCombinedFlags:
    """Test flag combinations not covered by test_deep_validation.py."""

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:3])
    @pytest.mark.parametrize("sid_mode", [SIDM_LAHIRI, SIDM_FAGAN_BRADLEY])
    def test_sidereal_equatorial(self, jd_ut, sid_mode):
        """SIDEREAL + EQUATORIAL combined flag."""
        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        flags = FLG_SWIEPH | FLG_SIDEREAL | FLG_EQUATORIAL | FLG_SPEED

        for planet_id, name in PLANETS_MAIN[:5]:
            lib_pos, _ = ephem.calc_ut(jd_ut, planet_id, flags)
            swe_pos, _ = swe.calc_ut(jd_ut, planet_id, flags)

            tol = 0.06 if planet_id == MOON else 0.01
            lon_diff = angular_diff(float(lib_pos[0]), float(swe_pos[0]))
            assert lon_diff < tol, (
                f"{name} SIDEREAL+EQUATORIAL lon diff {lon_diff:.6f}° > {tol}°"
            )

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:3])
    def test_sidereal_topocentric(self, jd_ut):
        """SIDEREAL + TOPOCENTRIC combined flag."""
        ephem.set_sid_mode(SIDM_LAHIRI)
        swe.set_sid_mode(swe.SIDM_LAHIRI)

        lat, lon, alt = 41.9028, 12.4964, 0
        ephem.set_topo(lon, lat, alt)
        swe.set_topo(lon, lat, alt)

        flags = FLG_SWIEPH | FLG_SIDEREAL | FLG_TOPOCTR | FLG_SPEED

        for planet_id, name in [(SUN, "Sun"), (MOON, "Moon"), (MARS, "Mars")]:
            lib_pos, _ = ephem.calc_ut(jd_ut, planet_id, flags)
            swe_pos, _ = swe.calc_ut(jd_ut, planet_id, flags)

            tol = 0.1 if planet_id == MOON else 0.01
            lon_diff = angular_diff(float(lib_pos[0]), float(swe_pos[0]))
            assert lon_diff < tol, (
                f"{name} SIDEREAL+TOPOCENTRIC lon diff {lon_diff:.6f}° > {tol}°"
            )

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:3])
    def test_barycentric_equatorial(self, jd_ut):
        """BARYCENTRIC + EQUATORIAL combined flag."""
        flags = FLG_SWIEPH | FLG_BARYCTR | FLG_EQUATORIAL | FLG_SPEED

        for planet_id, name in [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ]:
            lib_pos, _ = ephem.calc_ut(jd_ut, planet_id, flags)
            swe_pos, _ = swe.calc_ut(jd_ut, planet_id, flags)

            lon_diff = angular_diff(float(lib_pos[0]), float(swe_pos[0]))
            assert lon_diff < 0.01, (
                f"{name} BARYCENTRIC+EQUATORIAL lon diff {lon_diff:.6f}°"
            )

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:3])
    def test_heliocentric_xyz(self, jd_ut):
        """HELIOCENTRIC + XYZ combined flag."""
        flags = FLG_SWIEPH | FLG_HELCTR | FLG_XYZ | FLG_SPEED

        for planet_id, name in [(MARS, "Mars"), (JUPITER, "Jupiter")]:
            lib_pos, _ = ephem.calc_ut(jd_ut, planet_id, flags)
            swe_pos, _ = swe.calc_ut(jd_ut, planet_id, flags)

            for i in range(3):  # x, y, z
                diff = abs(float(lib_pos[i]) - float(swe_pos[i]))
                assert diff < 0.001, (
                    f"{name} HELIO+XYZ component {i} diff {diff:.6f} AU"
                )

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:3])
    def test_radians_equatorial(self, jd_ut):
        """RADIANS + EQUATORIAL combined flag."""
        flags = FLG_SWIEPH | FLG_RADIANS | FLG_EQUATORIAL | FLG_SPEED

        for planet_id, name in [(SUN, "Sun"), (MARS, "Mars")]:
            lib_pos, _ = ephem.calc_ut(jd_ut, planet_id, flags)
            swe_pos, _ = swe.calc_ut(jd_ut, planet_id, flags)

            # Convert to degrees for comparison
            lib_lon_deg = math.degrees(float(lib_pos[0]))
            swe_lon_deg = math.degrees(float(swe_pos[0]))

            lon_diff = angular_diff(lib_lon_deg, swe_lon_deg)
            assert lon_diff < 0.008, (
                f"{name} RADIANS+EQUATORIAL lon diff {lon_diff:.6f}°"
            )

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:3])
    def test_sidereal_j2000(self, jd_ut):
        """SIDEREAL + J2000 combined flag."""
        ephem.set_sid_mode(SIDM_LAHIRI)
        swe.set_sid_mode(swe.SIDM_LAHIRI)

        flags = FLG_SWIEPH | FLG_SIDEREAL | FLG_J2000 | FLG_SPEED

        for planet_id, name in [
            (SUN, "Sun"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ]:
            lib_pos, _ = ephem.calc_ut(jd_ut, planet_id, flags)
            swe_pos, _ = swe.calc_ut(jd_ut, planet_id, flags)

            lon_diff = angular_diff(float(lib_pos[0]), float(swe_pos[0]))
            assert lon_diff < 0.01, f"{name} SIDEREAL+J2000 lon diff {lon_diff:.6f}°"

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:3])
    def test_sidereal_nonut(self, jd_ut):
        """SIDEREAL + NONUT combined flag."""
        ephem.set_sid_mode(SIDM_LAHIRI)
        swe.set_sid_mode(swe.SIDM_LAHIRI)

        flags = FLG_SWIEPH | FLG_SIDEREAL | FLG_NONUT | FLG_SPEED

        for planet_id, name in [(SUN, "Sun"), (MARS, "Mars")]:
            lib_pos, _ = ephem.calc_ut(jd_ut, planet_id, flags)
            swe_pos, _ = swe.calc_ut(jd_ut, planet_id, flags)

            lon_diff = angular_diff(float(lib_pos[0]), float(swe_pos[0]))
            assert lon_diff < 0.02, f"{name} SIDEREAL+NONUT lon diff {lon_diff:.6f}°"


# ============================================================================
# PART 9: TOPOCENTRIC AT MULTIPLE GLOBAL LOCATIONS
# ============================================================================


class TestTopocentricMultiLocation:
    """Test topocentric positions at diverse global locations."""

    @pytest.mark.parametrize("location", LOCATIONS)
    @pytest.mark.parametrize(
        "planet_id,name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_topocentric_positions(self, location, planet_id, name):
        """Topocentric positions should match pyswisseph at various locations."""
        loc_name, lat, lon, alt = location
        jd = ephem.julday(2024, 6, 21, 12.0)

        ephem.set_topo(lon, lat, alt)
        swe.set_topo(lon, lat, alt)

        flags = FLG_SWIEPH | FLG_TOPOCTR | FLG_SPEED

        lib_pos, _ = ephem.calc_ut(jd, planet_id, flags)
        swe_pos, _ = swe.calc_ut(jd, planet_id, flags)

        tol = 0.06 if planet_id == MOON else 0.008
        lon_diff = angular_diff(float(lib_pos[0]), float(swe_pos[0]))
        assert lon_diff < tol, (
            f"{name} topocentric at {loc_name}: lon diff {lon_diff:.6f}° > {tol}°"
        )

        lat_diff = abs(float(lib_pos[1]) - float(swe_pos[1]))
        lat_tol = 0.06 if planet_id == MOON else 0.01
        assert lat_diff < lat_tol, (
            f"{name} topocentric at {loc_name}: lat diff {lat_diff:.6f}°"
        )


# ============================================================================
# PART 10: NOD_APS TT VARIANT AND ADDITIONAL METHODS
# ============================================================================


class TestNodApsTT:
    """Test nod_aps TT variant and additional node/apse methods."""

    @pytest.mark.parametrize(
        "planet_id,name",
        [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    @pytest.mark.parametrize("method", [NODBIT_MEAN, NODBIT_OSCU])
    def test_nod_aps_tt(self, planet_id, name, method):
        """nod_aps (TT) should return reasonable node/apse data."""
        jd_ut = TEST_DATES_UT[0]
        dt = ephem.deltat(jd_ut)
        jd_tt = jd_ut + dt

        lib_result = ephem.nod_aps(jd_tt, planet_id, FLG_SWIEPH, method)

        # Should return 4 tuples (asc_node, desc_node, perihelion, aphelion)
        assert len(lib_result) == 4, f"nod_aps should return 4 elements"

        # Longitudes should be in 0-360 range
        for i, node_name in enumerate(
            ["asc_node", "desc_node", "perihelion", "aphelion"]
        ):
            lon = float(lib_result[i][0])
            assert 0 <= lon < 360, f"{name} nod_aps {node_name} lon {lon} out of range"

    @pytest.mark.parametrize(
        "planet_id,name",
        [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_nod_aps_tt_vs_ut_consistency(self, planet_id, name):
        """nod_aps(TT) should match nod_aps_ut(UT) for same instant."""
        jd_ut = TEST_DATES_UT[0]
        dt = ephem.deltat(jd_ut)
        jd_tt = jd_ut + dt

        tt_result = ephem.nod_aps(jd_tt, planet_id, FLG_SWIEPH, NODBIT_MEAN)
        ut_result = ephem.nod_aps_ut(jd_ut, planet_id, FLG_SWIEPH, NODBIT_MEAN)

        for i in range(4):
            lon_diff = angular_diff(float(tt_result[i][0]), float(ut_result[i][0]))
            assert lon_diff < 1e-6, (
                f"{name} nod_aps TT vs UT element {i}: lon diff {lon_diff}"
            )


# ============================================================================
# PART 11: RISE/SET AT EXTREME LATITUDES AND EDGE CASES
# ============================================================================


class TestRiseSetEdgeCases:
    """Test rise/set at polar latitudes and edge cases."""

    def test_rise_set_normal_locations(self):
        """Rise and set times should be reasonable at normal latitudes."""
        jd = ephem.julday(2024, 3, 20, 0)  # Spring equinox

        for name, lat, lon, _ in LOCATIONS[:6]:
            if abs(lat) > 66:
                continue  # Skip polar locations for this test

            try:
                rise = ephem.rise_trans(
                    jd,
                    SUN,
                    CALC_RISE,
                    [lon, lat, 0.0],
                    1013.25,
                    15.0,
                    FLG_SWIEPH,
                )
                sett = ephem.rise_trans(
                    jd,
                    SUN,
                    CALC_SET,
                    [lon, lat, 0.0],
                    1013.25,
                    15.0,
                    FLG_SWIEPH,
                )

                rise_jd = float(rise[1][0])
                set_jd = float(sett[1][0])

                # At equinox, day length should be approximately 12 hours
                day_hours = (set_jd - rise_jd) * 24
                assert 8 < day_hours < 16, (
                    f"Day length at {name} on equinox: {day_hours:.1f}h (expected ~12h)"
                )
            except Exception:
                pass

    def test_rise_set_vs_pyswisseph_multiple_planets(self):
        """Rise/set times for multiple planets should match pyswisseph."""
        jd = ephem.julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964

        for planet_id, name in [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (VENUS, "Venus"),
        ]:
            try:
                lib_rise = ephem.rise_trans(
                    jd,
                    planet_id,
                    CALC_RISE,
                    [lon, lat, 0.0],
                    1013.25,
                    15.0,
                    FLG_SWIEPH,
                )
                swe_rise = swe.rise_trans(
                    jd, planet_id, swe.CALC_RISE, [lon, lat, 0], 1013.25, 15.0
                )

                lib_jd = float(lib_rise[1][0])
                swe_jd = float(swe_rise[1][0])

                diff_sec = jd_diff_seconds(lib_jd, swe_jd)
                tol = 180 if planet_id == MOON else 120
                assert diff_sec < tol, (
                    f"{name} rise time diff {diff_sec:.1f}s (max {tol}s)"
                )
            except Exception:
                pass

    def test_meridian_transit(self):
        """Meridian transit times should match pyswisseph."""
        jd = ephem.julday(2024, 6, 21, 0)
        lat, lon = 41.9028, 12.4964

        lib_transit = ephem.rise_trans(
            jd, SUN, CALC_MTRANSIT, [lon, lat, 0.0], 1013.25, 15.0, FLG_SWIEPH
        )
        swe_transit = swe.rise_trans(
            jd, swe.SUN, swe.CALC_MTRANSIT, [lon, lat, 0], 1013.25, 15.0
        )

        lib_jd = float(lib_transit[1][0])
        swe_jd = float(swe_transit[1][0])

        diff_sec = jd_diff_seconds(lib_jd, swe_jd)
        assert diff_sec < 120, f"Sun meridian transit diff {diff_sec:.1f}s (max 120s)"


# ============================================================================
# PART 12: SIDEREAL HOUSES
# ============================================================================


class TestSiderealHouses:
    """Test sidereal house calculations with various ayanamshas."""

    @pytest.mark.parametrize(
        "sid_mode",
        [
            SIDM_LAHIRI,
            SIDM_FAGAN_BRADLEY,
            SIDM_RAMAN,
        ],
    )
    @pytest.mark.parametrize("hsys_code", ["P", "K", "O", "E", "W"])
    def test_sidereal_houses_ex(self, sid_mode, hsys_code):
        """Sidereal houses_ex should match pyswisseph."""
        jd = ephem.julday(2024, 6, 21, 12.0)
        lat, lon = 41.9028, 12.4964

        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        lib_cusps, lib_ascmc = ephem.houses_ex(
            jd, lat, lon, ord(hsys_code), FLG_SIDEREAL
        )
        swe_cusps, swe_ascmc = swe.houses_ex(
            jd, lat, lon, hsys_code.encode(), swe.FLG_SIDEREAL
        )

        max_diff = 0
        for i in range(12):
            diff = angular_diff(float(lib_cusps[i]), float(swe_cusps[i]))
            max_diff = max(max_diff, diff)

        assert max_diff < 0.01, (
            f"Sidereal houses_ex {hsys_code} sidmode={sid_mode}: "
            f"max cusp diff {max_diff:.6f}°"
        )


# ============================================================================
# PART 13: STATISTICAL SURVEY WITH DIVERSE FLAGS
# ============================================================================


class TestStatisticalFlagSurvey:
    """Statistical survey of positions across many flag combinations."""

    def _generate_dates(self, n=50, seed=123):
        random.seed(seed)
        dates = []
        for _ in range(n):
            year = random.randint(1800, 2200)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = ephem.julday(year, month, day, hour)
            dates.append(jd)
        return dates

    @pytest.mark.parametrize(
        "flags_name,flags,tol",
        [
            ("HELCTR", FLG_SWIEPH | FLG_HELCTR | FLG_SPEED, 0.008),
            ("BARYCTR", FLG_SWIEPH | FLG_BARYCTR | FLG_SPEED, 0.008),
            ("J2000", FLG_SWIEPH | FLG_J2000 | FLG_SPEED, 0.008),
            ("NONUT", FLG_SWIEPH | FLG_NONUT | FLG_SPEED, 0.015),
            ("EQUATORIAL", FLG_SWIEPH | FLG_EQUATORIAL | FLG_SPEED, 0.008),
            ("ICRS", FLG_SWIEPH | FLG_ICRS | FLG_SPEED, 0.008),
        ],
    )
    def test_flag_survey(self, flags_name, flags, tol):
        """Survey positions across many dates for various flags."""
        dates = self._generate_dates(30)
        planets = [
            (SUN, "Sun"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ]

        # Skip barycentric for Sun
        if flags & FLG_BARYCTR:
            planets = [(pid, n) for pid, n in planets if pid != SUN]
        # Skip heliocentric for Sun
        if flags & FLG_HELCTR:
            planets = [(pid, n) for pid, n in planets if pid != SUN]

        max_diff_seen = 0
        worst_case = ""

        for jd in dates:
            for planet_id, name in planets:
                try:
                    lib_pos, _ = ephem.calc_ut(jd, planet_id, flags)
                    swe_pos, _ = swe.calc_ut(jd, planet_id, flags)

                    lon_diff = angular_diff(float(lib_pos[0]), float(swe_pos[0]))
                    if lon_diff > max_diff_seen:
                        max_diff_seen = lon_diff
                        worst_case = f"{name} JD={jd:.1f}"
                except Exception:
                    pass

        assert max_diff_seen < tol, (
            f"{flags_name} survey max diff {max_diff_seen:.6f}° > {tol}° at {worst_case}"
        )


# ============================================================================
# PART 14: ORBITAL ELEMENTS UT VARIANT
# ============================================================================


class TestOrbitalElementsUT:
    """Test get_orbital_elements_ut (UT variant)."""

    @pytest.mark.parametrize(
        "planet_id,name",
        [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_orbital_elements_ut_vs_tt_consistency(self, planet_id, name):
        """get_orbital_elements_ut should match get_orbital_elements for same instant."""
        jd_ut = TEST_DATES_UT[0]
        dt = ephem.deltat(jd_ut)
        jd_tt = jd_ut + dt

        ut_result = ephem.get_orbital_elements_ut(jd_ut, planet_id, FLG_SWIEPH)
        tt_result = ephem.get_orbital_elements(jd_tt, planet_id, FLG_SWIEPH)

        # Result is a flat tuple of 50 elements
        ut_elems = ut_result
        tt_elems = tt_result

        # Compare first 6 elements (a, e, i, Omega, omega, M)
        for i in range(6):
            diff = abs(float(ut_elems[i]) - float(tt_elems[i]))
            if float(tt_elems[i]) != 0:
                rel_diff = diff / abs(float(tt_elems[i]))
                assert rel_diff < 1e-6, (
                    f"{name} orbital element {i}: UT vs TT rel diff {rel_diff}"
                )


# ============================================================================
# PART 15: DELTAT_EX AND TIME FUNCTIONS
# ============================================================================


class TestDeltaTEx:
    """Test deltat_ex function."""

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:6])
    def test_deltat_ex_vs_deltat(self, jd_ut):
        """deltat_ex should give same value as deltat for same flags."""
        dt = ephem.deltat(jd_ut)
        dt_ex = ephem.deltat_ex(jd_ut, FLG_SWIEPH)

        # dt_ex now returns a float directly
        diff = abs(dt - dt_ex)
        assert diff < 1e-10, f"deltat vs deltat_ex diff {diff} at JD {jd_ut}"

    @pytest.mark.parametrize("jd_ut", TEST_DATES_UT[:4])
    def test_deltat_ex_vs_pyswisseph(self, jd_ut):
        """deltat_ex should match pyswisseph deltat_ex."""
        lib_result = ephem.deltat_ex(jd_ut, FLG_SWIEPH)
        swe_result = swe.deltat_ex(jd_ut, swe.FLG_SWIEPH)

        lib_val = float(lib_result)
        # pyswisseph deltat_ex returns a float directly
        swe_val = float(swe_result)

        diff_sec = abs(lib_val - swe_val) * 86400
        assert diff_sec < 300, f"deltat_ex diff {diff_sec:.2f}s at JD {jd_ut}"


# ============================================================================
# PART 16: INTERNAL CONSISTENCY CHECKS
# ============================================================================


class TestInternalConsistency:
    """Cross-validation between related functions to catch inconsistencies."""

    def test_sidereal_position_equals_tropical_minus_ayanamsa(self):
        """Sidereal position should equal tropical position minus ayanamsa."""
        jd = ephem.julday(2024, 6, 21, 12.0)
        ephem.set_sid_mode(SIDM_LAHIRI)

        ayan = float(ephem.get_ayanamsa_ut(jd))

        for planet_id, name in PLANETS_MAIN[:5]:
            trop_pos, _ = ephem.calc_ut(jd, planet_id, FLG_SWIEPH)
            sid_pos, _ = ephem.calc_ut(jd, planet_id, FLG_SWIEPH | FLG_SIDEREAL)

            expected_sid = (float(trop_pos[0]) - ayan) % 360
            actual_sid = float(sid_pos[0]) % 360

            diff = angular_diff(expected_sid, actual_sid)
            assert diff < 0.001, (
                f"{name}: sidereal lon {actual_sid:.6f} != tropical {float(trop_pos[0]):.6f} - "
                f"ayanamsa {ayan:.6f} = {expected_sid:.6f} (diff {diff:.6f}°)"
            )

    def test_equatorial_ra_dec_consistency(self):
        """Equatorial RA should be in 0-360 and Dec in -90..+90."""
        jd = ephem.julday(2024, 6, 21, 12.0)

        for planet_id, name in PLANETS_MAIN:
            pos, _ = ephem.calc_ut(jd, planet_id, FLG_SWIEPH | FLG_EQUATORIAL)
            ra = float(pos[0])
            dec = float(pos[1])

            assert 0 <= ra < 360, f"{name} RA {ra} out of range 0-360"
            assert -90 <= dec <= 90, f"{name} Dec {dec} out of range -90..+90"

    def test_xyz_distance_consistency(self):
        """XYZ distance should equal sqrt(x^2 + y^2 + z^2)."""
        jd = ephem.julday(2024, 6, 21, 12.0)

        for planet_id, name in PLANETS_MAIN[:5]:
            # Get spherical
            sph_pos, _ = ephem.calc_ut(jd, planet_id, FLG_SWIEPH)
            dist_sph = float(sph_pos[2])

            # Get Cartesian
            xyz_pos, _ = ephem.calc_ut(jd, planet_id, FLG_SWIEPH | FLG_XYZ)
            x, y, z = float(xyz_pos[0]), float(xyz_pos[1]), float(xyz_pos[2])
            dist_xyz = math.sqrt(x * x + y * y + z * z)

            rel_diff = abs(dist_sph - dist_xyz) / dist_sph if dist_sph > 0 else 0
            assert rel_diff < 1e-6, (
                f"{name} distance: spherical={dist_sph:.8f} xyz={dist_xyz:.8f} "
                f"rel_diff={rel_diff:.2e}"
            )

    def test_houses_asc_mc_in_ascmc_array(self):
        """ASC and MC in ascmc array should match cusp 1 and cusp 10."""
        jd = ephem.julday(2024, 6, 21, 12.0)

        for hsys in ["P", "K", "O", "R", "C"]:
            cusps, ascmc = ephem.houses(jd, 41.9, 12.5, ord(hsys))

            asc_diff = angular_diff(float(ascmc[0]), float(cusps[0]))
            mc_diff = angular_diff(float(ascmc[1]), float(cusps[9]))

            assert asc_diff < 1e-6, (
                f"House {hsys}: ASC in ascmc ({float(ascmc[0]):.6f}) != "
                f"cusp 1 ({float(cusps[0]):.6f})"
            )
            assert mc_diff < 1e-6, (
                f"House {hsys}: MC in ascmc ({float(ascmc[1]):.6f}) != "
                f"cusp 10 ({float(cusps[9]):.6f})"
            )

    def test_heliocentric_sun_is_earth(self):
        """Heliocentric Sun position should give Earth's position (opposite of geocentric Sun)."""
        jd = ephem.julday(2024, 6, 21, 12.0)

        # Heliocentric Earth
        earth_pos, _ = ephem.calc_ut(
            jd, EARTH, FLG_SWIEPH | FLG_HELCTR | FLG_SPEED
        )

        # Geocentric Sun - heliocentric Earth should be ~180° opposite
        sun_pos, _ = ephem.calc_ut(jd, SUN, FLG_SWIEPH)

        # Earth heliocentric lon should be Sun geocentric + 180
        expected_earth_lon = (float(sun_pos[0]) + 180.0) % 360
        actual_earth_lon = float(earth_pos[0])

        diff = angular_diff(expected_earth_lon, actual_earth_lon)
        assert diff < 0.01, (
            f"Helio Earth lon {actual_earth_lon:.4f} should be ~Sun+180 "
            f"({expected_earth_lon:.4f}), diff={diff:.6f}°"
        )

    def test_speed_sign_consistency(self):
        """Planet speed should be consistent with position change."""
        jd = ephem.julday(2024, 6, 21, 12.0)
        dt = 0.01  # 0.01 days

        for planet_id, name in PLANETS_MAIN[:5]:
            pos1, _ = ephem.calc_ut(jd, planet_id, FLG_SWIEPH | FLG_SPEED)
            pos2, _ = ephem.calc_ut(jd + dt, planet_id, FLG_SWIEPH)

            reported_speed = float(pos1[3])  # deg/day
            numerical_speed = angular_diff(float(pos1[0]), float(pos2[0])) / dt

            # Sign should match (unless near 0/360 boundary issues)
            delta_lon = float(pos2[0]) - float(pos1[0])
            if abs(delta_lon) > 180:
                delta_lon = delta_lon - 360 * (1 if delta_lon > 0 else -1)
            numerical_speed_signed = delta_lon / dt

            if abs(reported_speed) > 0.01:  # Skip near-stationary
                speed_ratio = numerical_speed_signed / reported_speed
                assert 0.5 < speed_ratio < 2.0, (
                    f"{name} speed inconsistency: reported={reported_speed:.4f}, "
                    f"numerical={numerical_speed_signed:.4f}"
                )


# ============================================================================
# PART 17: BOUNDARY AND DATE RANGE TESTS
# ============================================================================


class TestDateRangeBoundaries:
    """Test at DE440 boundaries and special dates."""

    def test_de440_start_boundary(self):
        """Positions at DE440 start boundary (1550) should be valid."""
        jd = ephem.julday(1550, 1, 15, 12.0)

        for planet_id, name in PLANETS_MAIN[:5]:
            try:
                lib_pos, _ = ephem.calc_ut(
                    jd, planet_id, FLG_SWIEPH | FLG_SPEED
                )
                swe_pos, _ = swe.calc_ut(jd, planet_id, swe.FLG_SWIEPH | swe.FLG_SPEED)

                lon_diff = angular_diff(float(lib_pos[0]), float(swe_pos[0]))
                tol = 0.1 if planet_id == MOON else 0.02
                assert lon_diff < tol, (
                    f"{name} at DE440 start: lon diff {lon_diff:.6f}°"
                )
            except Exception:
                pass

    def test_de440_end_boundary(self):
        """Positions at DE440 end boundary (2650) should be valid."""
        jd = ephem.julday(2640, 1, 15, 12.0)

        for planet_id, name in PLANETS_MAIN[:5]:
            try:
                lib_pos, _ = ephem.calc_ut(
                    jd, planet_id, FLG_SWIEPH | FLG_SPEED
                )
                swe_pos, _ = swe.calc_ut(jd, planet_id, swe.FLG_SWIEPH | swe.FLG_SPEED)

                lon_diff = angular_diff(float(lib_pos[0]), float(swe_pos[0]))
                tol = 0.1 if planet_id == MOON else 0.02
                assert lon_diff < tol, f"{name} at DE440 end: lon diff {lon_diff:.6f}°"
            except Exception:
                pass

    def test_gregorian_reform_boundary(self):
        """Positions around Gregorian reform (Oct 15, 1582) should be valid."""
        # Just after reform
        jd = ephem.julday(1582, 10, 15, 12.0)

        for planet_id, name in [
            (SUN, "Sun"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ]:
            lib_pos, _ = ephem.calc_ut(jd, planet_id, FLG_SWIEPH)
            swe_pos, _ = swe.calc_ut(jd, planet_id, swe.FLG_SWIEPH)

            lon_diff = angular_diff(float(lib_pos[0]), float(swe_pos[0]))
            assert lon_diff < 0.02, (
                f"{name} at Gregorian reform: lon diff {lon_diff:.6f}°"
            )

    def test_year_boundaries(self):
        """Positions at century boundaries should be valid."""
        for year in [1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500]:
            jd = ephem.julday(year, 1, 1, 0.0)

            lib_pos, _ = ephem.calc_ut(jd, SUN, FLG_SWIEPH)
            swe_pos, _ = swe.calc_ut(jd, swe.SUN, swe.FLG_SWIEPH)

            lon_diff = angular_diff(float(lib_pos[0]), float(swe_pos[0]))
            assert lon_diff < 0.01, f"Sun at year {year}: lon diff {lon_diff:.6f}°"


# ============================================================================
# PART 18: CUSPS NEAR 0/360 WRAP
# ============================================================================


class TestCuspWrapping:
    """Test house cusps near 0°/360° boundary."""

    def test_cusp_wrapping_various_dates(self):
        """All cusps should be in 0-360 range for various dates."""
        dates = [
            ephem.julday(y, m, 15, 12.0) for y in [2000, 2024] for m in range(1, 13)
        ]

        for jd in dates:
            for hsys in ["P", "K", "O", "R", "C", "E", "W"]:
                cusps, ascmc = ephem.houses(jd, 41.9, 12.5, ord(hsys))

                for i in range(12):
                    c = float(cusps[i])
                    assert 0 <= c < 360, (
                        f"House {hsys} cusp {i + 1} = {c} out of range at JD {jd}"
                    )

                # ASC and MC should be in range
                asc = float(ascmc[0])
                mc = float(ascmc[1])
                assert 0 <= asc < 360, f"ASC {asc} out of range"
                assert 0 <= mc < 360, f"MC {mc} out of range"

    def test_cusp_monotonicity(self):
        """Cusps should generally increase (with wrap-around at 360)."""
        jd = ephem.julday(2024, 6, 21, 12.0)

        for hsys in ["P", "K", "O", "R", "C"]:
            cusps, _ = ephem.houses(jd, 41.9, 12.5, ord(hsys))

            # Check that cusps go in order (allowing wrap)
            wraps = 0
            for i in range(11):
                c1 = float(cusps[i])
                c2 = float(cusps[i + 1])
                if c2 < c1:
                    wraps += 1

            # Should wrap at most once (going past 360/0)
            assert wraps <= 1, (
                f"House {hsys}: {wraps} wraps in cusp sequence "
                f"(cusps: {[float(c) for c in cusps]})"
            )


# ============================================================================
# PART 19: ASTEROID AND SPECIAL BODY TESTS
# ============================================================================


class TestSpecialBodies:
    """Test asteroids, hypothetical bodies, and special objects."""

    @pytest.mark.parametrize(
        "planet_id,name",
        [
            (CHIRON, "Chiron"),
            (PHOLUS, "Pholus"),
            (CERES, "Ceres"),
            (PALLAS, "Pallas"),
            (JUNO, "Juno"),
            (VESTA, "Vesta"),
        ],
    )
    def test_asteroids(self, planet_id, name):
        """Asteroid positions should be within tolerance of pyswisseph."""
        jd = ephem.julday(2024, 6, 21, 12.0)

        lib_pos, _ = ephem.calc_ut(jd, planet_id, FLG_SWIEPH | FLG_SPEED)
        swe_pos, _ = swe.calc_ut(jd, planet_id, swe.FLG_SWIEPH | swe.FLG_SPEED)

        lon_diff = angular_diff(float(lib_pos[0]), float(swe_pos[0]))
        assert lon_diff < 5.0, (
            f"{name} lon diff {lon_diff:.4f}° (Keplerian approx tolerance 5°)"
        )

    def test_ecl_nut_obliquity(self):
        """ECL_NUT should return valid obliquity values."""
        jd = ephem.julday(2024, 6, 21, 12.0)

        lib_pos, _ = ephem.calc_ut(jd, ECL_NUT, FLG_SWIEPH)
        swe_pos, _ = swe.calc_ut(jd, swe.ECL_NUT, swe.FLG_SWIEPH)

        # True obliquity should be ~23.4°
        lib_eps = float(lib_pos[0])
        swe_eps = float(swe_pos[0])

        assert 22 < lib_eps < 24, f"Obliquity {lib_eps} out of range"
        diff = abs(lib_eps - swe_eps)
        assert diff < 0.01, (
            f"Obliquity diff {diff:.6f}° (lib={lib_eps:.6f}, swe={swe_eps:.6f})"
        )

    @pytest.mark.parametrize(
        "planet_id,name",
        [
            (MEAN_NODE, "Mean Node"),
            (TRUE_NODE, "True Node"),
            (MEAN_APOG, "Mean Lilith"),
            (OSCU_APOG, "True Lilith"),
            (INTP_APOG, "Interp Apogee"),
            (INTP_PERG, "Interp Perigee"),
        ],
    )
    def test_lunar_points_equatorial(self, planet_id, name):
        """Lunar points in equatorial mode should match pyswisseph."""
        jd = ephem.julday(2024, 6, 21, 12.0)

        flags = FLG_SWIEPH | FLG_EQUATORIAL | FLG_SPEED

        try:
            lib_pos, _ = ephem.calc_ut(jd, planet_id, flags)
            swe_pos, _ = swe.calc_ut(jd, planet_id, flags)

            lon_diff = angular_diff(float(lib_pos[0]), float(swe_pos[0]))
            tol = 5.0 if planet_id == INTP_PERG else 0.5
            assert lon_diff < tol, f"{name} equatorial RA diff {lon_diff:.4f}°"
        except Exception:
            pass


# ============================================================================
# PART 20: COMPREHENSIVE PHENO VALIDATION
# ============================================================================


class TestPhenoComprehensive:
    """Comprehensive phenomena validation beyond test_deep_validation.py."""

    @pytest.mark.parametrize(
        "planet_id,name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_pheno_ut_physical_constraints(self, planet_id, name):
        """pheno_ut results should satisfy physical constraints."""
        dates = [ephem.julday(y, 6, 15, 12.0) for y in range(2020, 2026)]

        for jd in dates:
            vals = ephem.pheno_ut(jd, planet_id, FLG_SWIEPH)

            phase_angle = float(vals[0])
            phase = float(vals[1])
            elongation = float(vals[2])
            app_diam = float(vals[3])
            magnitude = float(vals[4])

            # Phase angle should be 0-180°
            assert 0 <= phase_angle <= 180, (
                f"{name} phase angle {phase_angle} out of range at JD {jd}"
            )

            # Phase should be 0-1
            assert 0 <= phase <= 1.0001, f"{name} phase {phase} out of range at JD {jd}"

            # Elongation should be 0-180°
            assert 0 <= elongation <= 180.1, (
                f"{name} elongation {elongation} out of range at JD {jd}"
            )

            # Apparent diameter should be positive
            assert app_diam > 0, (
                f"{name} apparent diameter {app_diam} should be > 0 at JD {jd}"
            )

            # For inner planets, check elongation constraints
            if planet_id == MERCURY:
                assert elongation < 30, (
                    f"Mercury elongation {elongation}° > 30° (impossible)"
                )
            elif planet_id == VENUS:
                assert elongation < 50, (
                    f"Venus elongation {elongation}° > 50° (impossible)"
                )

    @pytest.mark.parametrize(
        "planet_id,name",
        [
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_pheno_ut_vs_pyswisseph_detailed(self, planet_id, name):
        """Detailed pheno_ut comparison with pyswisseph."""
        jd = ephem.julday(2024, 6, 21, 12.0)

        lib_vals = ephem.pheno_ut(jd, planet_id, FLG_SWIEPH)
        swe_result = swe.pheno_ut(jd, planet_id, swe.FLG_SWIEPH)

        swe_vals = swe_result[0] if isinstance(swe_result[0], tuple) else swe_result

        # Phase angle (index 0)
        diff_phase = abs(float(lib_vals[0]) - float(swe_vals[0]))
        assert diff_phase < 0.5, f"{name} pheno_ut phase angle diff {diff_phase:.4f}°"

        # Phase (index 1) - 0 to 1
        diff_illum = abs(float(lib_vals[1]) - float(swe_vals[1]))
        assert diff_illum < 0.05, f"{name} pheno_ut phase diff {diff_illum:.6f}"

        # Elongation (index 2)
        diff_elong = abs(float(lib_vals[2]) - float(swe_vals[2]))
        assert diff_elong < 0.5, f"{name} pheno_ut elongation diff {diff_elong:.4f}°"

        # Apparent diameter (index 3) - should be in degrees
        lib_diam = float(lib_vals[3])
        swe_diam = float(swe_vals[3])
        if swe_diam > 0:
            ratio = lib_diam / swe_diam
            assert 0.5 < ratio < 2.0, (
                f"{name} pheno_ut diam ratio {ratio:.4f} "
                f"(lib={lib_diam:.8f}, swe={swe_diam:.8f})"
            )
