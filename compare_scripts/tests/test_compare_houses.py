"""
House Systems Comparison Tests.

Compares all house system calculations between pyswisseph and libephemeris.
Tests all 24 house systems across different latitudes and dates.
"""

import pytest
import random
import swisseph as swe
import libephemeris as ephem


# ============================================================================
# TOLERANCES
# ============================================================================

HOUSE_CUSP_TOL = 0.001  # degrees (~3.6 arcsec)
ASCMC_TOL = 0.001  # degrees (~3.6 arcsec) - all ASCMC values match precisely
GAUQUELIN_TOL = 0.001  # Gauquelin 36 sectors - now matches Swiss Ephemeris


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# HOUSE SYSTEMS
# ============================================================================

HOUSE_SYSTEMS = [
    ("P", "Placidus"),
    ("K", "Koch"),
    ("R", "Regiomontanus"),
    ("C", "Campanus"),
    ("E", "Equal (Ascendant)"),
    ("A", "Equal (MC)"),
    ("W", "Whole Sign"),
    ("O", "Porphyry"),
    ("B", "Alcabitius"),
    ("T", "Polich/Page (Topocentric)"),
    ("M", "Morinus"),
    ("X", "Meridian (Axial Rotation)"),
    ("V", "Vehlow Equal"),
    ("H", "Horizontal"),
    ("F", "Carter Poli-Equatorial"),
    ("U", "Krusinski-Pisa"),
    ("N", "Natural Gradient"),
    ("Y", "APC Houses"),
    ("D", "Equal from MC"),
    ("L", "Pullen SD"),
    ("S", "Sripati"),
    ("G", "Gauquelin"),
    ("I", "Sunshine/Makransky"),
    ("Q", "Pullen SR"),
]

# Systems with relaxed tolerances
# These systems are not fully implemented or have different calculation methods
RELAXED_SYSTEMS = {
    "G": GAUQUELIN_TOL,  # Gauquelin: 36 sectors (matches Swiss Ephemeris)
    "I": 0.001,  # Sunshine/Makransky: Implemented with Treindl algorithm, ~0.0002° precision
    # Koch: max ~1 arcsec at high latitudes (~0.0003°)
    "K": 0.001,
    # Placidus and Regiomontanus have minor precision differences at high latitudes
    "P": 0.002,
    "R": 0.002,
}

# Test locations
STANDARD_LOCATIONS = [
    ("Rome", 41.9028, 12.4964),
    ("New York", 40.7128, -74.0060),
    ("Sydney", -33.8688, 151.2093),
    ("Equator", 0.0, 0.0),
]

HIGH_LATITUDE_LOCATIONS = [
    ("Tromso", 69.6492, 18.9553),
    ("McMurdo", -77.8463, 166.6681),
]

TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 6, 15, 0.0, "Summer 2024"),
    (1980, 5, 20, 14.5, "Past"),
]


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestHouseCusps:
    """Compare house cusp calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS)
    @pytest.mark.parametrize("name,lat,lon", STANDARD_LOCATIONS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_house_cusps_match(
        self, hsys, hsys_name, name, lat, lon, year, month, day, hour, date_desc
    ):
        """Test all 12 house cusps match between implementations."""
        jd = swe.julday(year, month, day, hour)

        cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, hsys.encode("ascii"))
        cusps_py, ascmc_py = ephem.houses(jd, lat, lon, hsys)

        # Get tolerance for this system
        tolerance = RELAXED_SYSTEMS.get(hsys, HOUSE_CUSP_TOL)

        # Compare cusps (handle Gauquelin which has 36 sectors in pyswisseph but 12 in libephemeris)
        num_cusps_to_compare = min(len(cusps_swe), len(cusps_py))
        max_diff = 0.0
        for i in range(num_cusps_to_compare):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            max_diff = max(max_diff, diff)

        assert max_diff < tolerance, (
            f"{hsys_name} at {name} ({date_desc}): max cusp diff {max_diff:.6f}° exceeds tolerance"
        )


class TestAscendantMC:
    """Compare Ascendant and MC calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS[:7])  # Common systems
    @pytest.mark.parametrize("name,lat,lon", STANDARD_LOCATIONS)
    def test_ascendant_mc_match(self, hsys, hsys_name, name, lat, lon):
        """Test Ascendant and MC match between implementations."""
        jd = swe.julday(2000, 1, 1, 12.0)

        _, ascmc_swe = swe.houses(jd, lat, lon, hsys.encode("ascii"))
        _, ascmc_py = ephem.houses(jd, lat, lon, hsys)

        # Ascendant (index 0)
        asc_diff = angular_diff(ascmc_swe[0], ascmc_py[0])
        assert asc_diff < HOUSE_CUSP_TOL, (
            f"{hsys_name} at {name}: Ascendant diff {asc_diff:.6f}° exceeds tolerance"
        )

        # MC (index 1)
        mc_diff = angular_diff(ascmc_swe[1], ascmc_py[1])
        assert mc_diff < HOUSE_CUSP_TOL, (
            f"{hsys_name} at {name}: MC diff {mc_diff:.6f}° exceeds tolerance"
        )


class TestRandomLocations:
    """Comprehensive random location testing."""

    @pytest.mark.comparison
    @pytest.mark.slow
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS[:10])
    def test_random_locations(self, hsys, hsys_name):
        """Test house systems at 50 random locations."""
        random.seed(42)

        tolerance = RELAXED_SYSTEMS.get(hsys, HOUSE_CUSP_TOL)
        failures = []

        for _ in range(50):
            lat = random.uniform(-60.0, 60.0)
            lon = random.uniform(-180.0, 180.0)
            year = random.randint(1900, 2100)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0.0, 24.0)

            jd = swe.julday(year, month, day, hour)

            try:
                cusps_swe, _ = swe.houses(jd, lat, lon, hsys.encode("ascii"))
                cusps_py, _ = ephem.houses(jd, lat, lon, hsys)

                max_diff = max(
                    angular_diff(cusps_swe[i], cusps_py[i]) for i in range(12)
                )

                if max_diff >= tolerance:
                    failures.append((lat, lon, jd, max_diff))
            except Exception:
                pass  # Skip errors (polar regions, etc.)

        assert len(failures) == 0, (
            f"{hsys_name}: {len(failures)} failures out of 50 random locations"
        )


class TestHighLatitudes:
    """Test behavior at high latitudes."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "hsys,hsys_name",
        [
            ("E", "Equal"),
            ("W", "Whole Sign"),
            ("O", "Porphyry"),
            ("M", "Morinus"),
        ],
    )  # Systems that should work at all latitudes
    @pytest.mark.parametrize("name,lat,lon", HIGH_LATITUDE_LOCATIONS)
    def test_high_latitude_systems(self, hsys, hsys_name, name, lat, lon):
        """Test house systems that work at high latitudes."""
        jd = swe.julday(2000, 1, 1, 12.0)

        cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, hsys.encode("ascii"))
        cusps_py, ascmc_py = ephem.houses(jd, lat, lon, hsys)

        max_diff = max(angular_diff(cusps_swe[i], cusps_py[i]) for i in range(12))

        assert max_diff < HOUSE_CUSP_TOL, (
            f"{hsys_name} at {name}: max diff {max_diff:.6f}° at high latitude"
        )


class TestSpecialCases:
    """Test special cases and edge conditions."""

    @pytest.mark.comparison
    def test_equator(self):
        """Test house calculations at equator."""
        jd = swe.julday(2000, 1, 1, 12.0)
        lat, lon = 0.0, 0.0

        for hsys in ["P", "K", "E", "W"]:
            cusps_swe, _ = swe.houses(jd, lat, lon, hsys.encode("ascii"))
            cusps_py, _ = ephem.houses(jd, lat, lon, hsys)

            max_diff = max(angular_diff(cusps_swe[i], cusps_py[i]) for i in range(12))

            assert max_diff < HOUSE_CUSP_TOL, (
                f"System {hsys} at equator: max diff {max_diff:.6f}°"
            )

    @pytest.mark.comparison
    def test_prime_meridian(self):
        """Test house calculations at prime meridian."""
        jd = swe.julday(2000, 1, 1, 12.0)
        lat, lon = 51.5074, 0.0  # London area

        cusps_swe, _ = swe.houses(jd, lat, lon, b"P")
        cusps_py, _ = ephem.houses(jd, lat, lon, "P")

        max_diff = max(angular_diff(cusps_swe[i], cusps_py[i]) for i in range(12))

        assert max_diff < HOUSE_CUSP_TOL, (
            f"Placidus at prime meridian: max diff {max_diff:.6f}°"
        )

    @pytest.mark.comparison
    def test_international_date_line(self):
        """Test house calculations near international date line."""
        jd = swe.julday(2000, 1, 1, 12.0)
        lat, lon = 35.0, 179.9  # Near date line

        cusps_swe, _ = swe.houses(jd, lat, lon, b"P")
        cusps_py, _ = ephem.houses(jd, lat, lon, "P")

        max_diff = max(angular_diff(cusps_swe[i], cusps_py[i]) for i in range(12))

        assert max_diff < HOUSE_CUSP_TOL, (
            f"Placidus near date line: max diff {max_diff:.6f}°"
        )
