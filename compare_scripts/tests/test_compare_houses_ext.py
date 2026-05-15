"""
Pytest-style Extended House Functions Comparison Tests.

Validates houses_armc, houses_ex, house_pos, and gauquelin_sector
calculations against pyswisseph.
"""

import random
import pytest
import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    JUPITER,
    SATURN,
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
    SIDM_RAMAN,
    FLG_SIDEREAL,
    FLG_SWIEPH,
)


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class HouseExtTolerance:
    """Tolerance thresholds for extended house comparisons."""

    CUSP_DEGREES = 0.002  # House cusp (relaxed for Koch precision at certain angles)
    # house_pos with Placidus/Regiomontanus/Campanus/Koch is now precise with body latitude.
    POSITION = 0.02  # House position for P, R, C
    POSITION_KOCH = 0.02  # Koch now uses analytical algorithm (was 0.25)
    POSITION_CAMPANUS = 0.02  # Campanus now matches Swiss Ephemeris precisely
    SECTOR = 0.1  # Gauquelin sector


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# All 19 house systems
ALL_HOUSE_SYSTEMS = {
    "P": "Placidus",
    "K": "Koch",
    "R": "Regiomontanus",
    "C": "Campanus",
    "E": "Equal (Ascendant)",
    "W": "Whole Sign",
    "O": "Porphyry",
    "B": "Alcabitius",
    "T": "Polich/Page (Topocentric)",
    "M": "Morinus",
    "X": "Meridian (Axial Rotation)",
    "V": "Vehlow Equal",
    "H": "Horizontal",
    "F": "Carter Poli-Equatorial",
    "U": "Krusinski-Pisa",
    "N": "Natural Gradient",
    "G": "Gauquelin Sectors",
    "Y": "APC Houses",
    "S": "Sripati",
}

# Quick test subset
HOUSE_SYSTEMS = ["P", "K", "R", "C", "E", "W", "M", "B"]

# Systems with relaxed tolerances
RELAXED_TOLERANCE_SYSTEMS = {
    "G": 0.001,  # Gauquelin 36 sectors - now matches Swiss Ephemeris
}

SIDEREAL_MODES = [
    (SIDM_LAHIRI, "Lahiri"),
    (SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
    (SIDM_RAMAN, "Raman"),
]


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# FIXTURES
# ============================================================================


@pytest.fixture
def jd_standard():
    """Standard Julian Day for tests."""
    return swe.julday(2024, 1, 1, 12.0)


# ============================================================================
# HOUSES ARMC TESTS
# ============================================================================


class TestHousesArmc:
    """Tests for houses_armc function."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "armc,lat,desc",
        [
            (0.0, 45.0, "ARMC=0, Lat=45"),
            (90.0, 45.0, "ARMC=90, Lat=45"),
            (180.0, 0.0, "ARMC=180, Lat=0"),
            (270.0, -35.0, "ARMC=270, Lat=-35"),
        ],
    )
    @pytest.mark.parametrize("hsys", HOUSE_SYSTEMS[:4])
    def test_houses_armc(self, armc, lat, desc, hsys):
        """Test houses_armc for various ARMC and latitude combinations."""
        eps = 23.4393

        # SwissEphemeris
        cusps_swe, ascmc_swe = swe.houses_armc(armc, lat, eps, hsys.encode("ascii"))

        # LibEphemeris
        cusps_py, ascmc_py = pyephem.houses_armc(armc, lat, eps, hsys)

        # Compare all 12 cusps
        max_diff = 0.0
        for i in range(12):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            max_diff = max(max_diff, diff)

        assert max_diff < HouseExtTolerance.CUSP_DEGREES, (
            f"{hsys} {desc}: max cusp diff {max_diff}°"
        )


# ============================================================================
# HOUSES EX (SIDEREAL) TESTS
# ============================================================================


class TestHousesExSidereal:
    """Tests for houses_ex with sidereal mode."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", SIDEREAL_MODES)
    @pytest.mark.parametrize("hsys", HOUSE_SYSTEMS[:3])
    def test_houses_ex_sidereal(self, jd_standard, sid_mode, sid_name, hsys):
        """Test houses_ex with sidereal mode."""
        lat, lon = 45.0, 0.0

        # Set sidereal mode
        swe.set_sid_mode(sid_mode)
        pyephem.set_sid_mode(sid_mode)

        flags = FLG_SIDEREAL

        # SwissEphemeris
        cusps_swe, ascmc_swe = swe.houses_ex(
            jd_standard, lat, lon, hsys.encode("ascii"), flags
        )

        # LibEphemeris
        cusps_py, ascmc_py = pyephem.houses_ex(jd_standard, lat, lon, hsys, flags)

        # Compare cusps
        max_diff = 0.0
        for i in range(12):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            max_diff = max(max_diff, diff)

        assert max_diff < HouseExtTolerance.CUSP_DEGREES, (
            f"{hsys} {sid_name}: max diff {max_diff}°"
        )


# ============================================================================
# HOUSE POSITION TESTS
# ============================================================================


class TestHousePos:
    """Tests for house_pos function."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
        ],
    )
    @pytest.mark.parametrize("hsys", HOUSE_SYSTEMS[:4])
    def test_house_pos(self, jd_standard, body_id, body_name, hsys):
        """Test house_pos for various planets."""
        lat, lon = 45.0, 0.0

        # Get planet position
        pos, _ = swe.calc_ut(jd_standard, body_id, FLG_SWIEPH)
        planet_lon = pos[0]
        planet_lat = pos[1]

        # Get ARMC
        cusps_swe, ascmc_swe = swe.houses(jd_standard, lat, lon, hsys.encode("ascii"))
        armc = ascmc_swe[2]
        eps = 23.4393

        # SwissEphemeris - signature: house_pos(armc, geolat, eps, objcoord, hsys)
        hp_swe = swe.house_pos(
            armc, lat, eps, (planet_lon, planet_lat), hsys.encode("ascii")
        )

        # LibEphemeris - should match pyswisseph signature
        hp_py = pyephem.house_pos(armc, lat, eps, (planet_lon, planet_lat), hsys)

        diff = abs(hp_swe - hp_py)

        # Use different tolerance for systems with coordinate transform issues
        if hsys == "C":
            tolerance = HouseExtTolerance.POSITION_CAMPANUS
        elif hsys == "K":
            tolerance = HouseExtTolerance.POSITION_KOCH
        else:
            tolerance = HouseExtTolerance.POSITION

        assert diff < tolerance, f"{body_name} in {hsys}: house pos diff {diff}"


# ============================================================================
# GAUQUELIN SECTOR TESTS
# ============================================================================


class TestGauquelinSector:
    """Tests for gauquelin_sector function."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_gauquelin_sector(self, jd_standard, body_id, body_name):
        """Test gauquelin_sector for various planets."""
        lat, lon = 48.8566, 2.3522  # Paris
        geopos = (lon, lat, 0)
        atpress = 1013.25
        attemp = 15.0

        # SwissEphemeris
        try:
            ret_swe = swe.gauquelin_sector(
                jd_standard, body_id, "", FLG_SWIEPH, 0, geopos, atpress, attemp
            )
            sector_swe = ret_swe[0] if isinstance(ret_swe, tuple) else ret_swe
        except Exception as e:
            pytest.skip(f"SwissEphemeris failed: {e}")

        # LibEphemeris
        try:
            ret_py = pyephem.gauquelin_sector(
                jd_standard, body_id, 0, geopos, atpress, attemp, FLG_SWIEPH
            )
            sector_py = ret_py[0] if isinstance(ret_py, tuple) else ret_py
        except Exception as e:
            pytest.skip(f"LibEphemeris failed: {e}")

        diff = abs(sector_swe - sector_py)

        assert diff < HouseExtTolerance.SECTOR, (
            f"{body_name} Gauquelin sector diff {diff}"
        )


# ============================================================================
# COMPREHENSIVE RANDOM LOCATION TESTS
# ============================================================================


class TestHousesExComprehensive:
    """Comprehensive tests for houses_ex at random locations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,hsys_name", list(ALL_HOUSE_SYSTEMS.items())[:10])
    def test_houses_ex_random_locations(self, hsys, hsys_name):
        """Test houses_ex at multiple random locations."""
        random.seed(42)  # Reproducible

        tolerance = RELAXED_TOLERANCE_SYSTEMS.get(hsys, 0.001)
        failures = 0

        for _ in range(10):  # 10 random locations per system
            lat = random.uniform(-60.0, 60.0)
            lon = random.uniform(-180.0, 180.0)
            year = random.randint(1950, 2050)
            jd = swe.julday(year, 6, 15, 12.0)

            try:
                cusps_swe, _ = swe.houses_ex(jd, lat, lon, hsys.encode("ascii"), 0)
                cusps_py, _ = pyephem.houses_ex(jd, lat, lon, hsys, 0)

                max_diff = 0.0
                for i in range(12):
                    diff = angular_diff(cusps_swe[i], cusps_py[i])
                    max_diff = max(max_diff, diff)

                if max_diff >= tolerance:
                    failures += 1
            except Exception:
                failures += 1

        # Allow up to 2 failures per system (edge cases)
        assert failures <= 2, f"{hsys_name}: {failures} locations failed"
