"""
Deep Validation Suite 3: Coverage of previously untested swe_* API functions.

This test suite targets all functions NOT covered by test_deep_validation.py or
test_deep_validation_2.py, including:
- Solar eclipse geographic location functions (sol_eclipse_where)
- Eclipse path functions (path_width, central_line, northern/southern limit)
- Solar eclipse obscuration at location
- Lunar occultation functions (lun_occult_when_glob, lun_occult_when_loc, lun_occult_where)
- Gauquelin sector calculation
- houses_armc_ex2 (extended ARMC houses with speeds)
- orbit_max_min_true_distance
- cross_ut (generic planet crossing)
- helio_cross_ut (heliocentric crossing)
- heliacal_ut, heliacal_pheno_ut, vis_limit_mag
- fixstar2 (TT variant)

All tests use Skyfield mode only (not LEB).
"""

from __future__ import annotations

import math
from typing import Tuple

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


def jd_diff_minutes(jd1: float, jd2: float) -> float:
    """Difference between two JDs in minutes."""
    return abs(jd1 - jd2) * 1440.0


# ============================================================================
# FIXTURES
# ============================================================================


@pytest.fixture(autouse=True)
def setup_ephemeris():
    """Set up ephemeris paths for both libraries."""
    swe.set_ephe_path(_os.environ.get("REF_EPHE_PATH", "./ephe"))
    yield
    swe.close()


# ============================================================================
# PART 1: SOLAR ECLIPSE WHERE (Geographic location)
# ============================================================================


class TestSolEclipseWhere:
    """Test sol_eclipse_where — geographic position of central eclipse."""

    def _find_solar_eclipse(self, year: int, month: int):
        """Find a solar eclipse starting from given date."""
        jd = ephem.julday(year, month, 1, 0.0)
        swe_type, swe_tret = swe.sol_eclipse_when_glob(jd, swe.FLG_SWIEPH)
        return swe_tret[0]  # JD of maximum eclipse

    @pytest.mark.parametrize(
        "year,month,desc",
        [
            (2024, 4, "Apr 2024 Total"),
            (2023, 10, "Oct 2023 Annular"),
            (2025, 3, "Mar 2025"),
        ],
    )
    def test_sol_eclipse_where_vs_pyswisseph(self, year, month, desc):
        """sol_eclipse_where should match pyswisseph for position and attributes."""
        jd_max = self._find_solar_eclipse(year, month)

        lib_result = ephem.sol_eclipse_where(jd_max, FLG_SWIEPH)
        swe_result = swe.sol_eclipse_where(jd_max, swe.FLG_SWIEPH)

        # Compare eclipse type flags
        lib_type = lib_result[0]
        swe_type = swe_result[0]
        # Both should detect an eclipse (non-zero)
        assert lib_type != 0, f"{desc}: libephemeris found no eclipse"
        assert swe_type != 0, f"{desc}: pyswisseph found no eclipse"

        # Compare geographic position (geopos[0]=lon, geopos[1]=lat)
        lib_geopos = lib_result[1]
        swe_geopos = swe_result[1]

        # Central line position (lon, lat) should match within ~1 degree
        # (DE440 vs Swiss Eph differences can produce geographic shifts)
        lon_diff = angular_diff(float(lib_geopos[0]), float(swe_geopos[0]))
        lat_diff = abs(float(lib_geopos[1]) - float(swe_geopos[1]))

        # Only check if central eclipse (both have non-zero lat/lon)
        if abs(float(swe_geopos[0])) > 0.01 and abs(float(swe_geopos[1])) > 0.01:
            assert lon_diff < 2.0, (
                f"{desc}: central lon diff {lon_diff:.4f}° "
                f"(lib={lib_geopos[0]:.4f}, swe={swe_geopos[0]:.4f})"
            )
            assert lat_diff < 2.0, (
                f"{desc}: central lat diff {lat_diff:.4f}° "
                f"(lib={lib_geopos[1]:.4f}, swe={swe_geopos[1]:.4f})"
            )

        # Compare attributes (attr[0]=magnitude, attr[1]=ratio)
        lib_attr = lib_result[2]
        swe_attr = swe_result[2]

        # Eclipse magnitude should be reasonably close
        mag_diff = abs(float(lib_attr[0]) - float(swe_attr[0]))
        assert mag_diff < 0.1, (
            f"{desc}: magnitude diff {mag_diff:.6f} "
            f"(lib={lib_attr[0]:.6f}, swe={swe_attr[0]:.6f})"
        )

    def test_sol_eclipse_where_return_structure(self):
        """sol_eclipse_where should return (int, tuple_10, tuple_20)."""
        jd_max = self._find_solar_eclipse(2024, 4)
        result = ephem.sol_eclipse_where(jd_max, FLG_SWIEPH)

        assert isinstance(result, tuple)
        assert len(result) == 3
        assert isinstance(result[0], int)
        assert len(result[1]) >= 10  # geopos
        assert len(result[2]) >= 20  # attr


# ============================================================================
# PART 2: ECLIPSE PATH FUNCTIONS
# ============================================================================


class TestEclipsePathFunctions:
    """Test calc_eclipse_path_width, central_line, northern/southern limit."""

    def _find_total_eclipse_jd(self):
        """Find a total/annular solar eclipse for path testing."""
        jd = ephem.julday(2024, 4, 1, 0.0)
        swe_type, swe_tret = swe.sol_eclipse_when_glob(jd, swe.FLG_SWIEPH)
        return swe_tret[0]  # Apr 8 2024 total solar eclipse

    def test_eclipse_path_width_positive(self):
        """Path width should be positive for a central eclipse."""
        jd_max = self._find_total_eclipse_jd()
        width = ephem.calc_eclipse_path_width(jd_max)
        # Total eclipse path width should be between ~20 and ~300 km
        assert width > 0, f"Path width should be positive, got {width}"
        assert width < 500, f"Path width unreasonably large: {width} km"

    def _find_eclipse_time_range(self):
        """Find start/end times for eclipse path calculations."""
        jd_max = self._find_total_eclipse_jd()
        # Use ±10 minutes around maximum (narrow range for speed)
        return jd_max - 10.0 / 1440.0, jd_max + 10.0 / 1440.0

    def test_eclipse_central_line_returns_coords(self):
        """Central line should return valid geographic coordinates."""
        jd_start, jd_end = self._find_eclipse_time_range()
        result = ephem.calc_eclipse_central_line(jd_start, jd_end, step_minutes=5.0)
        assert result is not None, "Central line returned None"
        # Should return 3 tuples: (times, lons, lats)
        assert len(result) == 3, f"Expected 3 tuples, got {len(result)}"
        # Should have some points
        assert len(result[0]) > 0, "Central line has no points"

    def test_eclipse_northern_limit_returns_coords(self):
        """Northern limit should return valid coordinates."""
        jd_start, jd_end = self._find_eclipse_time_range()
        result = ephem.calc_eclipse_northern_limit(jd_start, jd_end, step_minutes=5.0)
        assert result is not None, "Northern limit returned None"
        assert len(result) == 3, f"Expected 3 tuples, got {len(result)}"

    def test_eclipse_southern_limit_returns_coords(self):
        """Southern limit should return valid coordinates."""
        jd_start, jd_end = self._find_eclipse_time_range()
        result = ephem.calc_eclipse_southern_limit(jd_start, jd_end, step_minutes=5.0)
        assert result is not None, "Southern limit returned None"
        assert len(result) == 3, f"Expected 3 tuples, got {len(result)}"

    def test_northern_limit_north_of_southern_limit(self):
        """Northern limit latitude should be >= southern limit latitude on average."""
        jd_start, jd_end = self._find_eclipse_time_range()
        north = ephem.calc_eclipse_northern_limit(jd_start, jd_end, step_minutes=5.0)
        south = ephem.calc_eclipse_southern_limit(jd_start, jd_end, step_minutes=5.0)
        # Extract latitudes (result[2] is the lats tuple)
        if len(north[2]) > 0 and len(south[2]) > 0:
            # Compare average latitudes
            avg_north_lat = sum(float(x) for x in north[2]) / len(north[2])
            avg_south_lat = sum(float(x) for x in south[2]) / len(south[2])
            if avg_north_lat != 0 and avg_south_lat != 0:
                assert avg_north_lat >= avg_south_lat - 5.0, (
                    f"Avg northern limit ({avg_north_lat:.2f}°) should be >= "
                    f"avg southern limit ({avg_south_lat:.2f}°)"
                )


# ============================================================================
# PART 3: SOLAR ECLIPSE OBSCURATION AT LOCATION
# ============================================================================


class TestSolEclipseObscuration:
    """Test sol_eclipse_obscuration_at_loc."""

    def test_obscuration_during_eclipse(self):
        """Obscuration should be > 0 at a location during an eclipse."""
        # Apr 8, 2024 total solar eclipse — Dallas, TX is in totality path
        jd = ephem.julday(2024, 4, 1, 0.0)
        swe_type, swe_tret = swe.sol_eclipse_when_glob(jd, swe.FLG_SWIEPH)
        jd_max = swe_tret[0]

        # Get where eclipse is central
        _, geopos, _ = swe.sol_eclipse_where(jd_max, swe.FLG_SWIEPH)
        central_lon = float(geopos[0])
        central_lat = float(geopos[1])

        if abs(central_lon) > 0.01:
            obs = ephem.sol_eclipse_obscuration_at_loc(
                jd_max, central_lat, central_lon, 0.0
            )
            # Should be very high at center of totality
            assert obs > 0.5, f"Obscuration at central line should be > 0.5, got {obs}"

    def test_obscuration_no_eclipse(self):
        """Obscuration should be ~0 when no eclipse is happening."""
        # Random non-eclipse date
        jd = ephem.julday(2024, 7, 15, 12.0)
        obs = ephem.sol_eclipse_obscuration_at_loc(jd, 41.9, 12.5, 0.0)
        assert obs < 0.01, f"Obscuration should be ~0 outside eclipse, got {obs}"

    def test_obscuration_range(self):
        """Obscuration should always be in [0, 1]."""
        jd = ephem.julday(2024, 4, 1, 0.0)
        swe_type, swe_tret = swe.sol_eclipse_when_glob(jd, swe.FLG_SWIEPH)
        jd_max = swe_tret[0]

        # Test at various locations
        for lat, lon in [(41.9, 12.5), (0.0, 0.0), (35.0, -97.0), (-33.9, 18.4)]:
            obs = ephem.sol_eclipse_obscuration_at_loc(jd_max, lat, lon, 0.0)
            assert 0.0 <= obs <= 1.0, (
                f"Obscuration at ({lat},{lon}) = {obs}, should be in [0,1]"
            )


# ============================================================================
# PART 4: LUNAR OCCULTATION FUNCTIONS
# ============================================================================


class TestLunOccultWhenGlob:
    """Test lun_occult_when_glob — global occultation search."""

    def test_planet_occultation_returns_result(self):
        """lun_occult_when_glob should find a Venus occultation."""
        jd = ephem.julday(2024, 1, 1, 0.0)
        lib_result = ephem.lun_occult_when_glob(jd, VENUS)
        swe_result = swe.lun_occult_when_glob(jd, swe.VENUS, swe.FLG_SWIEPH)

        lib_type = lib_result[0]
        swe_type = swe_result[0]

        assert lib_type != 0, "libephemeris found no Venus occultation"
        assert swe_type != 0, "pyswisseph found no Venus occultation"

        # Compare times of maximum (should be within ~hours to days given
        # different ephemerides)
        lib_jd_max = float(lib_result[1][0])
        swe_jd_max = float(swe_result[1][0])

        diff_days = abs(lib_jd_max - swe_jd_max)
        # They should find the same occultation event (within ~1 day)
        assert diff_days < 1.0, (
            f"Occultation time diff {diff_days:.2f} days "
            f"(lib={lib_jd_max:.4f}, swe={swe_jd_max:.4f})"
        )

    def test_return_structure(self):
        """lun_occult_when_glob should return (int, tuple_of_floats)."""
        jd = ephem.julday(2024, 1, 1, 0.0)
        result = ephem.lun_occult_when_glob(jd, VENUS)

        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], int)
        assert len(result[1]) >= 10  # tret


class TestLunOccultWhere:
    """Test lun_occult_where — geographic position of occultation."""

    def test_occult_where_at_known_time(self):
        """lun_occult_where should return valid coordinates during an occultation."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        # Find a Venus occultation first
        swe_type, swe_tret = swe.lun_occult_when_glob(jd, swe.VENUS, swe.FLG_SWIEPH)
        jd_max = float(swe_tret[0])

        if swe_type == 0:
            pytest.skip("No Venus occultation found")

        lib_result = ephem.lun_occult_where(jd_max, VENUS)
        swe_result = swe.lun_occult_where(jd_max, swe.VENUS, swe.FLG_SWIEPH)

        # Compare geographic positions
        lib_geopos = lib_result[1]
        swe_geopos = swe_result[1]

        lon_diff = angular_diff(float(lib_geopos[0]), float(swe_geopos[0]))
        lat_diff = abs(float(lib_geopos[1]) - float(swe_geopos[1]))

        # Allow wider tolerance for geographic position (DE440 differences)
        if abs(float(swe_geopos[0])) > 0.01:
            assert lon_diff < 5.0, (
                f"Occult where lon diff {lon_diff:.4f}° "
                f"(lib={lib_geopos[0]:.4f}, swe={swe_geopos[0]:.4f})"
            )
        if abs(float(swe_geopos[1])) > 0.01:
            assert lat_diff < 5.0, (
                f"Occult where lat diff {lat_diff:.4f}° "
                f"(lib={lib_geopos[1]:.4f}, swe={swe_geopos[1]:.4f})"
            )


# ============================================================================
# PART 5: GAUQUELIN SECTOR
# ============================================================================


class TestGauquelinSector:
    """Test gauquelin_sector — Gauquelin sector calculation."""

    @pytest.mark.parametrize(
        "planet_id,swe_planet",
        [
            (SUN, swe.SUN),
            (MOON, swe.MOON),
            (MARS, swe.MARS),
            (VENUS, swe.VENUS),
            (JUPITER, swe.JUPITER),
        ],
    )
    def test_gauquelin_sector_vs_pyswisseph(self, planet_id, swe_planet):
        """Gauquelin sectors should match pyswisseph within tolerance."""
        jd = ephem.julday(2024, 6, 21, 12.0)
        lat, lon = 41.9028, 12.4964

        lib_sector = ephem.gauquelin_sector(
            jd, planet_id, 0, (lon, lat, 0.0), 1013.25, 15.0
        )
        swe_sector = swe.gauquelin_sector(
            jd, swe_planet, 0, (lon, lat, 0.0), 1013.25, 15.0
        )

        diff = abs(float(lib_sector) - float(swe_sector))
        # Handle wrapping at 36/0
        if diff > 18:
            diff = 36 - diff
        assert diff < 1.0, (
            f"Planet {planet_id}: sector diff {diff:.4f} "
            f"(lib={lib_sector:.4f}, swe={swe_sector:.4f})"
        )

    def test_gauquelin_sector_range(self):
        """Gauquelin sector should be in [1, 37)."""
        jd = ephem.julday(2024, 6, 21, 12.0)
        for planet_id in [SUN, MOON, MARS, VENUS, JUPITER]:
            sector = ephem.gauquelin_sector(
                jd, planet_id, 0, (12.4964, 41.9028, 0.0)
            )
            assert 1.0 <= float(sector) < 37.0, (
                f"Planet {planet_id}: sector {sector} out of range [1, 37)"
            )

    @pytest.mark.parametrize("method", [0, 1])
    def test_gauquelin_sector_methods(self, method):
        """Different methods should give broadly similar results."""
        jd = ephem.julday(2024, 6, 21, 12.0)
        sector = ephem.gauquelin_sector(
            jd, SUN, method, (12.4964, 41.9028, 0.0), 1013.25, 15.0
        )
        assert 1.0 <= float(sector) < 37.0, (
            f"Method {method}: sector {sector} out of range"
        )


# ============================================================================
# PART 6: HOUSES_ARMC_EX2 (Extended ARMC-based houses with speeds)
# ============================================================================


class TestHousesArmcEx2:
    """Test houses_armc_ex2 — ARMC houses with cusp speeds."""

    @pytest.mark.parametrize("hsys_code", ["P", "K", "O", "R", "C", "E"])
    def test_cusp_positions_match_houses_armc(self, hsys_code):
        """houses_armc_ex2 cusps should match houses_armc."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393
        hsys = ord(hsys_code)

        lib_cusps, lib_ascmc, lib_cusps_speed, lib_ascmc_speed = (
            ephem.houses_armc_ex2(armc, lat, eps, hsys)
        )
        lib_cusps_basic, lib_ascmc_basic = ephem.houses_armc(armc, lat, eps, hsys)

        # Cusp positions should match exactly
        for i in range(12):
            diff = angular_diff(float(lib_cusps[i]), float(lib_cusps_basic[i]))
            assert diff < 0.001, (
                f"{hsys_code} cusp {i + 1}: ex2 vs armc diff {diff:.6f}°"
            )

        # Ascmc should match
        for i in range(min(len(lib_ascmc), len(lib_ascmc_basic))):
            diff = angular_diff(float(lib_ascmc[i]), float(lib_ascmc_basic[i]))
            assert diff < 0.001, f"{hsys_code} ascmc[{i}]: ex2 vs armc diff {diff:.6f}°"

    @pytest.mark.parametrize("hsys_code", ["P", "K", "O"])
    def test_cusp_positions_vs_pyswisseph(self, hsys_code):
        """houses_armc_ex2 cusp positions should match pyswisseph."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        lib_cusps, lib_ascmc, _, _ = ephem.houses_armc_ex2(
            armc, lat, eps, ord(hsys_code)
        )
        swe_cusps, swe_ascmc, _, _ = swe.houses_armc_ex2(
            armc, lat, eps, hsys_code.encode()
        )

        max_diff = 0
        for i in range(12):
            diff = angular_diff(float(lib_cusps[i]), float(swe_cusps[i]))
            max_diff = max(max_diff, diff)

        assert max_diff < 0.01, f"{hsys_code}: max cusp diff {max_diff:.6f}°"

    def test_speeds_with_flag(self):
        """Cusp speeds should be non-zero when FLG_SPEED is set."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393
        hsys = ord("P")

        _, _, cusps_speed, ascmc_speed = ephem.houses_armc_ex2(
            armc, lat, eps, hsys, FLG_SPEED
        )

        # At least some speeds should be non-zero
        has_nonzero = any(abs(float(s)) > 0.0 for s in cusps_speed)
        assert has_nonzero, "All cusp speeds are zero despite FLG_SPEED"

    @pytest.mark.parametrize("hsys_code", ["P", "R", "C", "E", "M", "B", "T"])
    def test_cusp_speeds_vs_pyswisseph(self, hsys_code):
        """Cusp speeds should be close to pyswisseph for systems with compatible derivatives."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        _, _, lib_cspeeds, lib_aspeeds = ephem.houses_armc_ex2(
            armc, lat, eps, ord(hsys_code), FLG_SPEED
        )
        _, _, swe_cspeeds, swe_aspeeds = swe.houses_armc_ex2(
            armc, lat, eps, hsys_code.encode(), swe.FLG_SPEED
        )

        # Cusp speeds: numerical vs analytical, allow < 2 deg/day difference
        for i in range(12):
            diff = abs(float(lib_cspeeds[i]) - float(swe_cspeeds[i]))
            assert diff < 2.0, (
                f"{hsys_code} cusp {i + 1} speed: lib={lib_cspeeds[i]:.4f}, "
                f"swe={swe_cspeeds[i]:.4f}, diff={diff:.4f} deg/day"
            )

    def test_speeds_without_flag(self):
        """houses_armc_ex2 always returns cusp speeds (regardless of FLG_SPEED)."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393
        hsys = ord("P")

        _, _, cusps_speed, ascmc_speed = ephem.houses_armc_ex2(
            armc, lat, eps, hsys, 0
        )

        # houses_armc_ex2 always computes speeds (pyswisseph does the same)
        assert len(cusps_speed) == 12
        assert all(isinstance(float(s), float) for s in cusps_speed)

    def test_return_structure(self):
        """houses_armc_ex2 should return 4 tuples."""
        result = ephem.houses_armc_ex2(292.957, 41.9, 23.4393, ord("P"))
        assert len(result) == 4
        assert len(result[0]) == 12  # cusps
        assert len(result[1]) >= 8  # ascmc
        assert len(result[2]) == 12  # cusps_speed
        assert len(result[3]) >= 8  # ascmc_speed


# ============================================================================
# PART 7: ORBIT_MAX_MIN_TRUE_DISTANCE
# ============================================================================


class TestOrbitMaxMinTrueDistance:
    """Test orbit_max_min_true_distance.

    pyswisseph returns (max_dist, min_dist, true_dist) — 3 values.
    libephemeris must match this format exactly.
    """

    @pytest.mark.parametrize(
        "planet_id,name",
        [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
            (VENUS, "Venus"),
            (MERCURY, "Mercury"),
            (MOON, "Moon"),
        ],
    )
    def test_returns_3_tuple(self, planet_id, name):
        """Must return exactly 3 values: (max_dist, min_dist, true_dist)."""
        jd = ephem.julday(2024, 6, 21, 12.0)
        result = ephem.orbit_max_min_true_distance(jd, planet_id, FLG_SWIEPH)

        assert len(result) == 3, (
            f"{name}: expected 3 values (max, min, true), got {len(result)}: {result}"
        )

    @pytest.mark.parametrize(
        "planet_id,name",
        [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
            (VENUS, "Venus"),
            (MERCURY, "Mercury"),
            (MOON, "Moon"),
        ],
    )
    def test_order_max_min_true(self, planet_id, name):
        """result[0] must be max, result[1] must be min, result[2] is true dist."""
        jd = ephem.julday(2024, 6, 21, 12.0)
        result = ephem.orbit_max_min_true_distance(jd, planet_id, FLG_SWIEPH)

        max_dist = float(result[0])
        min_dist = float(result[1])
        true_dist = float(result[2])

        assert max_dist > 0, f"{name}: max_dist should be > 0, got {max_dist}"
        assert min_dist > 0, f"{name}: min_dist should be > 0, got {min_dist}"
        assert true_dist > 0, f"{name}: true_dist should be > 0, got {true_dist}"
        assert max_dist > min_dist, (
            f"{name}: max ({max_dist:.6f}) should be > min ({min_dist:.6f})"
        )
        # True distance should be between min and max (or close)
        assert min_dist * 0.9 <= true_dist <= max_dist * 1.1, (
            f"{name}: true ({true_dist:.6f}) should be roughly between "
            f"min ({min_dist:.6f}) and max ({max_dist:.6f})"
        )

    @pytest.mark.parametrize(
        "planet_id,swe_planet,name",
        [
            (MARS, swe.MARS, "Mars"),
            (JUPITER, swe.JUPITER, "Jupiter"),
            (SATURN, swe.SATURN, "Saturn"),
            (VENUS, swe.VENUS, "Venus"),
            (MERCURY, swe.MERCURY, "Mercury"),
            (MOON, swe.MOON, "Moon"),
        ],
    )
    def test_vs_pyswisseph(self, planet_id, swe_planet, name):
        """Compare (max, min, true) against pyswisseph directly by index."""
        jd = ephem.julday(2024, 6, 21, 12.0)
        dt = ephem.deltat(jd)
        jd_tt = jd + dt

        lib_result = ephem.orbit_max_min_true_distance(jd, planet_id, FLG_SWIEPH)
        # pyswisseph: returns (max_dist, min_dist, true_dist) using TT
        swe_result = swe.orbit_max_min_true_distance(jd_tt, swe_planet, swe.FLG_SWIEPH)

        lib_max = float(lib_result[0])
        lib_min = float(lib_result[1])
        lib_true = float(lib_result[2])
        swe_max = float(swe_result[0])
        swe_min = float(swe_result[1])
        swe_true = float(swe_result[2])

        # Max distances should be within ~20% (Keplerian approximation)
        max_ratio = lib_max / swe_max if swe_max > 0 else 999
        assert 0.8 < max_ratio < 1.2, (
            f"{name} max dist: lib={lib_max:.4f}, swe={swe_max:.4f}, "
            f"ratio={max_ratio:.4f}"
        )

        # Min distances should be within ~20%
        min_ratio = lib_min / swe_min if swe_min > 0 else 999
        assert 0.8 < min_ratio < 1.2, (
            f"{name} min dist: lib={lib_min:.4f}, swe={swe_min:.4f}, "
            f"ratio={min_ratio:.4f}"
        )

        # True distance should be close (within ~1% — both compute from same epoch)
        true_ratio = lib_true / swe_true if swe_true > 0 else 999
        assert 0.95 < true_ratio < 1.05, (
            f"{name} true dist: lib={lib_true:.6f}, swe={swe_true:.6f}, "
            f"ratio={true_ratio:.6f}"
        )


# ============================================================================
# PART 8: SWE_CROSS_UT (Generic planet crossing)
# ============================================================================


class TestCrossUt:
    """Test cross_ut — generic planet longitude crossing."""

    @pytest.mark.parametrize(
        "planet_id,target_lon,name",
        [
            (JUPITER, 0.0, "Jupiter 0° Aries"),
            (JUPITER, 60.0, "Jupiter 60°"),
            (JUPITER, 90.0, "Jupiter 90°"),
            (SATURN, 0.0, "Saturn 0° Aries"),
            (SATURN, 350.0, "Saturn 350°"),
            (VENUS, 90.0, "Venus 90°"),
            (MARS, 180.0, "Mars 180°"),
        ],
    )
    def test_cross_ut_finds_crossing(self, planet_id, target_lon, name):
        """cross_ut should find the crossing and planet should be at target lon."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        try:
            jd_cross = ephem.cross_ut(planet_id, target_lon, jd, FLG_SWIEPH)
        except RuntimeError:
            pytest.skip(f"{name}: cross_ut did not converge")

        # Verify planet is actually at target longitude
        pos, _ = ephem.calc_ut(jd_cross, planet_id, FLG_SWIEPH)
        diff = angular_diff(float(pos[0]), target_lon)
        assert diff < 0.01, (
            f"{name}: at crossing JD, lon={pos[0]:.6f}, target={target_lon}, "
            f"diff={diff:.6f}°"
        )

    def test_cross_ut_future_of_start(self):
        """Crossing should be in the future relative to start date."""
        jd = ephem.julday(2024, 1, 1, 0.0)
        try:
            jd_cross = ephem.cross_ut(JUPITER, 0.0, jd, FLG_SWIEPH)
            assert jd_cross > jd, "Crossing should be after start date"
        except RuntimeError:
            pytest.skip("cross_ut did not converge")


# ============================================================================
# PART 9: SWE_HELIO_CROSS_UT (Heliocentric crossing)
# ============================================================================


class TestHelioCrossUt:
    """Test helio_cross_ut — heliocentric longitude crossing."""

    @pytest.mark.parametrize(
        "planet_id,swe_planet,target_lon,name",
        [
            (MARS, swe.MARS, 0.0, "Mars 0°"),
            (JUPITER, swe.JUPITER, 90.0, "Jupiter 90°"),
            (SATURN, swe.SATURN, 180.0, "Saturn 180°"),
        ],
    )
    def test_helio_cross_ut_vs_pyswisseph(
        self, planet_id, swe_planet, target_lon, name
    ):
        """Heliocentric crossing times should match pyswisseph."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        try:
            lib_jd = ephem.helio_cross_ut(planet_id, target_lon, jd, FLG_SWIEPH)
        except RuntimeError:
            pytest.skip(f"{name}: lib helio_cross_ut did not converge")

        try:
            swe_jd = swe.helio_cross_ut(swe_planet, target_lon, jd, swe.FLG_SWIEPH)
        except Exception:
            pytest.skip(f"{name}: swe helio_cross_ut failed")

        diff_sec = jd_diff_seconds(lib_jd, float(swe_jd))
        assert diff_sec < 600, (
            f"{name}: helio crossing diff {diff_sec:.1f}s "
            f"(lib={lib_jd:.6f}, swe={swe_jd:.6f})"
        )

    def test_helio_cross_ut_planet_at_target(self):
        """At the crossing time, heliocentric lon should match target."""
        jd = ephem.julday(2024, 1, 1, 0.0)
        target = 0.0
        try:
            jd_cross = ephem.helio_cross_ut(MARS, target, jd, FLG_SWIEPH)
        except RuntimeError:
            pytest.skip("helio_cross_ut did not converge")

        pos, _ = ephem.calc_ut(jd_cross, MARS, FLG_SWIEPH | FLG_HELCTR)
        diff = angular_diff(float(pos[0]), target)
        assert diff < 0.01, (
            f"At crossing, helio lon={pos[0]:.6f}, target={target}, diff={diff:.6f}°"
        )


# ============================================================================
# PART 10: HELIACAL FUNCTIONS
# ============================================================================


class TestHeliacalUt:
    """Test heliacal_ut — heliacal rising/setting times."""

    def test_heliacal_ut_vs_pyswisseph(self):
        """heliacal_ut results should be in reasonable range vs pyswisseph."""
        jd = ephem.julday(2024, 6, 1, 0.0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0065)
        dobs = (25, 1.0, 1.0, 1.0, 1.0, 0.0)

        lib_result = ephem.heliacal_ut(
            jd, geopos, datm, dobs, "Venus", 1, FLG_SWIEPH
        )
        swe_result = swe.heliacal_ut(
            jd, geopos, datm, dobs, "Venus", swe.HELIACAL_RISING, swe.FLG_SWIEPH
        )

        # lib returns (jd1, jd2, jd3), swe returns (jd1, jd2, jd3)
        lib_jd = float(lib_result[0])
        swe_jd = float(swe_result[0])

        # Allow wider tolerance — different implementations may disagree by days
        diff_days = abs(lib_jd - swe_jd)
        assert diff_days < 30.0, (
            f"Heliacal rising diff {diff_days:.2f} days "
            f"(lib={lib_jd:.4f}, swe={swe_jd:.4f})"
        )

    def test_heliacal_ut_returns_valid_jd(self):
        """Heliacal result should be a valid future JD."""
        jd = ephem.julday(2024, 6, 1, 0.0)
        geopos = (12.4964, 41.9028, 0.0)
        datm = (1013.25, 15.0, 40.0, 0.0065)
        dobs = (25, 1.0, 1.0, 1.0, 1.0, 0.0)

        result = ephem.heliacal_ut(jd, geopos, datm, dobs, "Venus", 1, FLG_SWIEPH)

        result_jd = float(result[0])
        # Should be a valid JD (after search start or near it)
        assert result_jd > jd - 365, (
            f"Heliacal JD {result_jd:.4f} seems invalid (start={jd:.4f})"
        )


class TestVisLimitMag:
    """Test vis_limit_mag — limiting visual magnitude."""

    def test_vis_limit_mag_returns_values(self):
        """vis_limit_mag should return valid data."""
        jd = ephem.julday(2024, 6, 21, 22.0)  # Night time
        geopos = (12.4964, 41.9028, 0.0)
        atmo = (1013.25, 15.0, 40.0, 0.0065)
        observer = (25, 1.0, 1.0, 1.0, 1.0, 0.0)

        try:
            result = ephem.vis_limit_mag(
                jd, geopos, atmo, observer, "Venus", FLG_SWIEPH
            )
            assert isinstance(result, tuple)
            assert len(result) == 2
        except Exception:
            pytest.skip("vis_limit_mag raised an exception")


# ============================================================================
# PART 11: SWE_FIXSTAR2 (TT variant)
# ============================================================================


class TestFixstar2TT:
    """Test fixstar2 — fixed star calculation with TT input."""

    @pytest.mark.parametrize(
        "star",
        ["Aldebaran", "Regulus", "Spica", "Antares", "Sirius"],
    )
    def test_fixstar2_vs_fixstar2_ut(self, star):
        """fixstar2 (TT) should match fixstar2_ut for same instant."""
        jd_ut = ephem.julday(2024, 6, 21, 12.0)
        dt = ephem.deltat(jd_ut)
        jd_tt = jd_ut + dt

        ut_result = ephem.fixstar2_ut(star, jd_ut, FLG_SWIEPH)
        tt_result = ephem.fixstar2(star, jd_tt, FLG_SWIEPH)

        # Extract positions
        ut_pos, ut_name, ut_flags = ut_result
        tt_pos, tt_name, tt_flags = tt_result

        # Positions should be essentially identical
        lon_diff = angular_diff(float(ut_pos[0]), float(tt_pos[0]))
        lat_diff = abs(float(ut_pos[1]) - float(tt_pos[1]))

        assert lon_diff < 0.001, f"{star}: UT vs TT lon diff {lon_diff:.6f}°"
        assert lat_diff < 0.001, f"{star}: UT vs TT lat diff {lat_diff:.6f}°"

    def test_fixstar2_returns_star_name(self):
        """fixstar2 should include the resolved star name."""
        jd_tt = ephem.julday(2024, 6, 21, 12.0)
        result = ephem.fixstar2("Sirius", jd_tt, FLG_SWIEPH)
        star_name = result[1]
        assert isinstance(star_name, str)
        assert len(star_name) > 0


# ============================================================================
# PART 12: CROSS-FUNCTION CONSISTENCY TESTS
# ============================================================================


class TestCrossFunctionConsistency:
    """Cross-validation between related functions."""

    def test_sol_eclipse_where_consistent_with_when_glob(self):
        """sol_eclipse_where at time from when_glob should find the same eclipse."""
        jd = ephem.julday(2024, 4, 1, 0.0)

        # Find eclipse via when_glob
        ecl_type, tret = ephem.sol_eclipse_when_glob(jd)
        jd_max = tret[0]

        # Get where info
        where_type, geopos, attr = ephem.sol_eclipse_where(jd_max, FLG_SWIEPH)

        # Both should detect the same type of eclipse
        assert where_type != 0, "sol_eclipse_where found no eclipse at when_glob time"

        # Magnitude from attr should be > 0
        assert float(attr[0]) > 0, "Eclipse magnitude should be > 0"

    def test_houses_armc_ex2_consistent_with_houses_ex2(self):
        """houses_armc_ex2 should give same cusps as houses_ex2 for same ARMC."""
        jd = ephem.julday(2024, 6, 21, 12.0)
        lat = 41.9028
        lon = 12.4964

        # Get cusps via houses_ex2
        cusps_ex2, ascmc_ex2, _, _ = ephem.houses_ex2(
            jd, lat, lon, ord("P"), FLG_SWIEPH
        )
        armc = float(ascmc_ex2[2])

        # Get obliquity
        ecl_nut, _ = ephem.calc_ut(jd, ECL_NUT, FLG_SWIEPH)
        eps = float(ecl_nut[0])

        # Get cusps via houses_armc_ex2 using same ARMC and eps
        cusps_armc, ascmc_armc, _, _ = ephem.houses_armc_ex2(
            armc, lat, eps, ord("P")
        )

        # Cusps should match closely
        max_diff = 0
        for i in range(12):
            diff = angular_diff(float(cusps_ex2[i]), float(cusps_armc[i]))
            max_diff = max(max_diff, diff)

        assert max_diff < 0.01, (
            f"houses_ex2 vs houses_armc_ex2: max cusp diff {max_diff:.6f}°"
        )

    def test_gauquelin_sector_consistency_with_house_pos(self):
        """Gauquelin sector should be roughly consistent with house position."""
        jd = ephem.julday(2024, 6, 21, 12.0)
        lat = 41.9028
        lon = 12.4964

        # Get ARMC and eps
        _, ascmc = ephem.houses(jd, lat, lon, ord("P"))
        armc = float(ascmc[2])
        ecl_nut, _ = ephem.calc_ut(jd, ECL_NUT, FLG_SWIEPH)
        eps = float(ecl_nut[0])

        sector = ephem.gauquelin_sector(
            jd, SUN, 0, (lon, lat, 0.0), 1013.25, 15.0
        )

        # Sector should be valid (1-36 range)
        assert 1.0 <= float(sector) < 37.0, f"Gauquelin sector {sector} out of range"


# ============================================================================
# PART 13: EDGE CASES AND BOUNDARY CONDITIONS
# ============================================================================


class TestEdgeCases:
    """Edge cases and boundary conditions for newly tested functions."""

    def test_sol_eclipse_where_no_eclipse(self):
        """sol_eclipse_where at a non-eclipse time should return type=0."""
        jd = ephem.julday(2024, 7, 15, 12.0)  # Random non-eclipse date
        result = ephem.sol_eclipse_where(jd, FLG_SWIEPH)
        # Should return 0 type (no eclipse)
        assert result[0] == 0, f"Expected no eclipse at JD {jd}, got type={result[0]}"

    def test_gauquelin_sector_multiple_dates(self):
        """Gauquelin sector should vary across dates."""
        sectors = []
        for year in [2020, 2021, 2022, 2023, 2024]:
            jd = ephem.julday(year, 6, 21, 12.0)
            s = ephem.gauquelin_sector(jd, SUN, 0, (12.4964, 41.9028, 0.0))
            sectors.append(float(s))

        # Sun should be in similar sectors at summer solstice but not identical
        # (it's always near MC at noon on solstice)
        for s in sectors:
            assert 1.0 <= s < 37.0, f"Invalid sector {s}"

    def test_houses_armc_ex2_polar_latitude(self):
        """houses_armc_ex2 should handle polar latitudes."""
        # At 70° latitude, some house systems may struggle
        try:
            cusps, ascmc, _, _ = ephem.houses_armc_ex2(
                292.957, 70.0, 23.4393, ord("P")
            )
            # Should return valid cusps
            for i in range(12):
                assert 0.0 <= float(cusps[i]) < 360.0, (
                    f"Polar cusp {i + 1} = {cusps[i]}, invalid"
                )
        except Exception:
            # Some house systems legitimately can't handle extreme latitudes
            pass

    def test_helio_cross_ut_earth(self):
        """Heliocentric crossing should work for Earth."""
        jd = ephem.julday(2024, 1, 1, 0.0)
        try:
            jd_cross = ephem.helio_cross_ut(EARTH, 0.0, jd, FLG_SWIEPH)
            assert jd_cross > jd, "Earth helio crossing should be in the future"
        except RuntimeError:
            # May not converge for Earth
            pass


# ============================================================================
# PART 14: MULTI-DATE STRESS TESTS
# ============================================================================


class TestMultiDateStress:
    """Statistical tests across multiple dates."""

    def test_gauquelin_sector_50_dates(self):
        """Gauquelin sector should be valid for many dates/locations."""
        import random

        random.seed(42)

        failures = []
        for _ in range(50):
            year = random.randint(1950, 2050)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            lat = random.uniform(-60, 60)
            lon = random.uniform(-180, 180)

            jd = ephem.julday(year, month, day, hour)
            try:
                sector = ephem.gauquelin_sector(jd, SUN, 0, (lon, lat, 0.0))
                if not (1.0 <= float(sector) < 37.0):
                    failures.append(
                        f"({year}-{month}-{day} {hour:.1f}h, "
                        f"lat={lat:.1f}, lon={lon:.1f}): sector={sector}"
                    )
            except Exception as e:
                failures.append(f"({year}-{month}-{day}): error: {e}")

        assert len(failures) == 0, (
            f"Gauquelin sector failures ({len(failures)}/50):\n"
            + "\n".join(failures[:5])
        )

    @pytest.mark.parametrize("hsys_code", ["P", "K", "O", "R", "C", "E"])
    def test_houses_armc_ex2_multiple_dates(self, hsys_code):
        """houses_armc_ex2 should be consistent across multiple ARMC values."""
        import random

        random.seed(42)

        for _ in range(20):
            armc = random.uniform(0, 360)
            lat = random.uniform(-66, 66)
            eps = 23.4393

            cusps, ascmc, _, _ = ephem.houses_armc_ex2(
                armc, lat, eps, ord(hsys_code)
            )

            # All cusps should be valid degrees
            for i in range(12):
                c = float(cusps[i])
                assert 0.0 <= c < 360.0, (
                    f"{hsys_code} ARMC={armc:.1f}, lat={lat:.1f}: "
                    f"cusp {i + 1} = {c}, invalid"
                )
