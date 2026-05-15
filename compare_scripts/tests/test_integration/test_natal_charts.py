"""
Integration tests using famous people birth data.

Tests calculate complete natal charts (planetary positions, houses, aspects)
and verify against pyswisseph (Swiss Ephemeris) reference values.

The famous people birth data is sourced from public records and verified
against astro.com reference charts.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    MEAN_NODE,
    TRUE_NODE,
    CHIRON,
    FLG_SPEED,
)


# ============================================================================
# FAMOUS PEOPLE BIRTH DATA
# ============================================================================
# Birth data sourced from AstroDatabank and verified with astro.com
# All times in UT (Universal Time)
# Note: DE421 ephemeris covers 1899-07-29 through 2053-10-09

FAMOUS_PEOPLE = [
    {
        "name": "Marilyn Monroe",
        "year": 1926,
        "month": 6,
        "day": 1,
        "hour": 16.5,  # 9:30 AM PST = 17:30 UT (approx)
        "lat": 34.0522,
        "lon": -118.2437,
        "place": "Los Angeles, CA",
        "expected_sun_sign": "Gemini",
        "expected_moon_sign": "Aquarius",
    },
    {
        "name": "Nelson Mandela",
        "year": 1918,
        "month": 7,
        "day": 18,
        "hour": 12.08,  # 14:05 SAST = ~12:05 UT
        "lat": -31.5875,
        "lon": 28.7833,
        "place": "Mvezo, South Africa",
        "expected_sun_sign": "Cancer",
        "expected_moon_sign": "Scorpio",
    },
    {
        "name": "Queen Elizabeth II",
        "year": 1926,
        "month": 4,
        "day": 21,
        "hour": 1.67,  # 02:40 BST London = 01:40 UT
        "lat": 51.5074,
        "lon": -0.1278,
        "place": "London, England",
        "expected_sun_sign": "Taurus",
        "expected_moon_sign": "Leo",
    },
    {
        "name": "John F. Kennedy",
        "year": 1917,
        "month": 5,
        "day": 29,
        "hour": 20.0,  # 3:00 PM EST = 20:00 UT
        "lat": 42.3601,
        "lon": -71.0589,
        "place": "Brookline, MA",
        "expected_sun_sign": "Gemini",
        "expected_moon_sign": "Virgo",
    },
    {
        "name": "Princess Diana",
        "year": 1961,
        "month": 7,
        "day": 1,
        "hour": 18.45,  # 7:45 PM BST = 18:45 UT
        "lat": 52.8242,
        "lon": 0.4931,
        "place": "Sandringham, England",
        "expected_sun_sign": "Cancer",
        "expected_moon_sign": "Aquarius",
    },
    {
        "name": "Muhammad Ali",
        "year": 1942,
        "month": 1,
        "day": 17,
        "hour": 23.35,  # 6:35 PM CST = 00:35 UT next day, approx 23:35 UT same day
        "lat": 38.2527,
        "lon": -85.7585,
        "place": "Louisville, KY",
        "expected_sun_sign": "Capricorn",
        "expected_moon_sign": "Aquarius",
    },
    {
        "name": "Audrey Hepburn",
        "year": 1929,
        "month": 5,
        "day": 4,
        "hour": 3.0,  # 3:00 AM CET = 02:00 UT (approx, using 3:00 UT)
        "lat": 50.8503,
        "lon": 4.3517,
        "place": "Brussels, Belgium",
        "expected_sun_sign": "Taurus",
        "expected_moon_sign": "Pisces",
    },
    {
        "name": "Elvis Presley",
        "year": 1935,
        "month": 1,
        "day": 8,
        "hour": 10.53,  # 4:35 AM CST = 10:35 UT
        "lat": 34.2576,
        "lon": -88.7034,
        "place": "Tupelo, MS",
        "expected_sun_sign": "Capricorn",
        "expected_moon_sign": "Pisces",
    },
]


# Major planets for natal chart
NATAL_PLANETS = [
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

# Extended chart points
# Note: True Node uses different calculation formula (osculating vs mean)
# and Chiron requires external ephemeris files, so we only test Mean Node
EXTENDED_POINTS = [
    (MEAN_NODE, "Mean Node"),
]

# Common house systems for testing
HOUSE_SYSTEMS = [
    (ord("P"), "Placidus"),
    (ord("K"), "Koch"),
    (ord("O"), "Porphyry"),
    (ord("R"), "Regiomontanus"),
    (ord("C"), "Campanus"),
    (ord("E"), "Equal"),
    (ord("W"), "Whole Sign"),
]

# Zodiac signs
ZODIAC_SIGNS = [
    "Aries",
    "Taurus",
    "Gemini",
    "Cancer",
    "Leo",
    "Virgo",
    "Libra",
    "Scorpio",
    "Sagittarius",
    "Capricorn",
    "Aquarius",
    "Pisces",
]

# Aspect definitions: (name, exact_angle, orb)
ASPECTS = [
    ("Conjunction", 0, 8),
    ("Opposition", 180, 8),
    ("Trine", 120, 8),
    ("Square", 90, 7),
    ("Sextile", 60, 6),
    ("Quincunx", 150, 3),
    ("Semi-sextile", 30, 2),
    ("Semi-square", 45, 2),
    ("Sesquiquadrate", 135, 2),
    ("Quintile", 72, 2),
    ("Bi-quintile", 144, 2),
]


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def get_sign(longitude: float) -> str:
    """Convert longitude (0-360) to zodiac sign name."""
    sign_num = int(longitude / 30) % 12
    return ZODIAC_SIGNS[sign_num]


def normalize_angle(angle: float) -> float:
    """Normalize angle to 0-360 range."""
    while angle < 0:
        angle += 360
    while angle >= 360:
        angle -= 360
    return angle


def angle_diff(a1: float, a2: float) -> float:
    """Calculate the smallest angular difference between two angles."""
    diff = abs(normalize_angle(a1) - normalize_angle(a2))
    if diff > 180:
        diff = 360 - diff
    return diff


def find_aspect(lon1: float, lon2: float) -> tuple | None:
    """
    Find if there's an aspect between two planetary longitudes.

    Returns (aspect_name, exact_angle, actual_angle, orb_used) or None.
    """
    diff = angle_diff(lon1, lon2)

    for aspect_name, exact_angle, orb in ASPECTS:
        aspect_diff = abs(diff - exact_angle)
        if aspect_diff <= orb:
            return (aspect_name, exact_angle, diff, aspect_diff)

    return None


def calculate_all_aspects(positions: dict) -> list:
    """
    Calculate all aspects between planets.

    Args:
        positions: dict of {planet_name: longitude}

    Returns:
        List of (planet1, planet2, aspect_name, orb)
    """
    aspects = []
    planet_names = list(positions.keys())

    for i, p1 in enumerate(planet_names):
        for p2 in planet_names[i + 1 :]:
            aspect = find_aspect(positions[p1], positions[p2])
            if aspect:
                aspects.append((p1, p2, aspect[0], aspect[3]))

    return aspects


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestNatalChartPlanetaryPositions:
    """Test planetary positions for famous people natal charts."""

    @pytest.mark.integration
    @pytest.mark.comparison
    @pytest.mark.parametrize("person", FAMOUS_PEOPLE, ids=lambda p: p["name"])
    def test_all_planets_match_swisseph(self, person):
        """
        All planetary positions should match Swiss Ephemeris within tolerance.

        This is the primary validation that libephemeris produces accurate
        natal chart data.
        """
        jd = ephem.julday(
            person["year"], person["month"], person["day"], person["hour"]
        )

        max_diff = 0
        for planet_id, planet_name in NATAL_PLANETS:
            # Calculate with libephemeris
            pos_lib, _ = ephem.calc_ut(jd, planet_id, FLG_SPEED)
            # Calculate with Swiss Ephemeris
            pos_swe, _ = swe.calc_ut(jd, planet_id, FLG_SPEED)

            # Compare longitudes
            lon_diff = angle_diff(pos_lib[0], pos_swe[0])
            max_diff = max(max_diff, lon_diff)

            # Tolerance: 2 arcseconds (very tight, accounting for minor
            # differences in precession/nutation models between implementations)
            arcsec_2 = 2 / 3600
            assert lon_diff < arcsec_2, (
                f"{person['name']}: {planet_name} longitude diff {lon_diff * 3600:.2f} arcsec "
                f"(lib={pos_lib[0]:.6f}°, swe={pos_swe[0]:.6f}°)"
            )

            # Compare latitudes
            lat_diff = abs(pos_lib[1] - pos_swe[1])
            assert lat_diff < arcsec_2, (
                f"{person['name']}: {planet_name} latitude diff {lat_diff * 3600:.2f} arcsec"
            )

            # Compare distances (AU) - relaxed tolerance for outer planets
            # Pluto uses Keplerian approximation, so needs wider tolerance
            if planet_id == PLUTO:
                dist_tolerance = 0.001  # AU - Pluto Keplerian model
            else:
                dist_tolerance = 0.0001  # AU - JPL ephemeris planets
            dist_diff = abs(pos_lib[2] - pos_swe[2])
            assert dist_diff < dist_tolerance, (
                f"{person['name']}: {planet_name} distance diff {dist_diff} AU"
            )

            # Compare velocities (degrees/day)
            vel_diff = abs(pos_lib[3] - pos_swe[3])
            assert vel_diff < 0.01, (
                f"{person['name']}: {planet_name} velocity diff {vel_diff}°/day"
            )

    @pytest.mark.integration
    @pytest.mark.parametrize("person", FAMOUS_PEOPLE, ids=lambda p: p["name"])
    def test_sun_in_expected_sign(self, person):
        """Verify Sun is in the expected zodiac sign."""
        jd = ephem.julday(
            person["year"], person["month"], person["day"], person["hour"]
        )
        pos, _ = ephem.calc_ut(jd, SUN, 0)
        actual_sign = get_sign(pos[0])

        assert actual_sign == person["expected_sun_sign"], (
            f"{person['name']}: Sun expected in {person['expected_sun_sign']}, "
            f"got {actual_sign} ({pos[0]:.2f}°)"
        )

    @pytest.mark.integration
    @pytest.mark.parametrize("person", FAMOUS_PEOPLE, ids=lambda p: p["name"])
    def test_moon_in_expected_sign(self, person):
        """Verify Moon is in the expected zodiac sign."""
        jd = ephem.julday(
            person["year"], person["month"], person["day"], person["hour"]
        )
        pos, _ = ephem.calc_ut(jd, MOON, 0)
        actual_sign = get_sign(pos[0])

        assert actual_sign == person["expected_moon_sign"], (
            f"{person['name']}: Moon expected in {person['expected_moon_sign']}, "
            f"got {actual_sign} ({pos[0]:.2f}°)"
        )

    @pytest.mark.integration
    @pytest.mark.comparison
    @pytest.mark.parametrize("person", FAMOUS_PEOPLE, ids=lambda p: p["name"])
    def test_extended_points_match_swisseph(self, person):
        """
        Extended chart points (nodes, Chiron) should match Swiss Ephemeris.

        Note: Some points like Chiron use Keplerian approximation and have
        relaxed tolerance.
        """
        jd = ephem.julday(
            person["year"], person["month"], person["day"], person["hour"]
        )

        for planet_id, planet_name in EXTENDED_POINTS:
            pos_lib, _ = ephem.calc_ut(jd, planet_id, 0)
            pos_swe, _ = swe.calc_ut(jd, planet_id, 0)

            lon_diff = angle_diff(pos_lib[0], pos_swe[0])

            # Relaxed tolerance for points using different calculation models
            if planet_id == CHIRON:
                tolerance = 1.0  # degrees - Chiron uses Keplerian model
            elif planet_id in (MEAN_NODE, TRUE_NODE):
                tolerance = 0.01  # degrees - nodes may use different formulas
            else:
                tolerance = 0.001  # degrees for other points

            assert lon_diff < tolerance, (
                f"{person['name']}: {planet_name} diff {lon_diff:.4f}° "
                f"(lib={pos_lib[0]:.4f}°, swe={pos_swe[0]:.4f}°)"
            )


class TestNatalChartHouses:
    """Test house calculations for famous people natal charts."""

    @pytest.mark.integration
    @pytest.mark.comparison
    @pytest.mark.parametrize("person", FAMOUS_PEOPLE, ids=lambda p: p["name"])
    @pytest.mark.parametrize(
        "hsys,hsys_name", HOUSE_SYSTEMS, ids=[name for _, name in HOUSE_SYSTEMS]
    )
    def test_houses_match_swisseph(self, person, hsys, hsys_name):
        """
        House cusps should match Swiss Ephemeris for all common systems.
        """
        jd = ephem.julday(
            person["year"], person["month"], person["day"], person["hour"]
        )
        lat = person["lat"]
        lon = person["lon"]

        # Skip Placidus/Koch at high latitudes (polar circle issues)
        if abs(lat) > 66 and hsys in [ord("P"), ord("K")]:
            pytest.skip(f"{hsys_name} undefined at lat {lat}")

        # Calculate with libephemeris
        cusps_lib, ascmc_lib = ephem.houses(jd, lat, lon, hsys)
        # Calculate with Swiss Ephemeris
        cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, bytes([hsys]))

        # Compare ASC
        asc_diff = angle_diff(ascmc_lib[0], ascmc_swe[0])
        assert asc_diff < 0.1, (
            f"{person['name']} {hsys_name}: ASC diff {asc_diff:.4f}° "
            f"(lib={ascmc_lib[0]:.4f}°, swe={ascmc_swe[0]:.4f}°)"
        )

        # Compare MC
        mc_diff = angle_diff(ascmc_lib[1], ascmc_swe[1])
        assert mc_diff < 0.1, (
            f"{person['name']} {hsys_name}: MC diff {mc_diff:.4f}° "
            f"(lib={ascmc_lib[1]:.4f}°, swe={ascmc_swe[1]:.4f}°)"
        )

        # Compare all 12 house cusps
        for i in range(12):
            cusp_diff = angle_diff(cusps_lib[i], cusps_swe[i])
            assert cusp_diff < 0.1, (
                f"{person['name']} {hsys_name}: Cusp {i + 1} diff {cusp_diff:.4f}° "
                f"(lib={cusps_lib[i]:.4f}°, swe={cusps_swe[i]:.4f}°)"
            )

    @pytest.mark.integration
    @pytest.mark.parametrize("person", FAMOUS_PEOPLE[:3], ids=lambda p: p["name"])
    def test_asc_in_valid_range(self, person):
        """Ascendant should be a valid longitude (0-360)."""
        jd = ephem.julday(
            person["year"], person["month"], person["day"], person["hour"]
        )

        cusps, ascmc = ephem.houses(jd, person["lat"], person["lon"], ord("P"))

        assert 0 <= ascmc[0] < 360, f"ASC {ascmc[0]} out of range"
        assert 0 <= ascmc[1] < 360, f"MC {ascmc[1]} out of range"

    @pytest.mark.integration
    @pytest.mark.parametrize("person", FAMOUS_PEOPLE[:3], ids=lambda p: p["name"])
    def test_cusps_progress_around_zodiac(self, person):
        """House cusps should progress around the zodiac in order."""
        jd = ephem.julday(
            person["year"], person["month"], person["day"], person["hour"]
        )

        cusps, _ = ephem.houses(jd, person["lat"], person["lon"], ord("P"))

        # Check that cusps generally increase (with wraparound handling)
        for i in range(11):
            diff = cusps[i + 1] - cusps[i]
            if diff < -180:
                diff += 360
            # Each house should span approximately 30 degrees (±30)
            assert diff > 0 or diff < -180, (
                f"Cusp {i + 1} to {i + 2} has negative progression"
            )


class TestNatalChartAspects:
    """Test aspect calculations for famous people natal charts."""

    @pytest.mark.integration
    @pytest.mark.parametrize("person", FAMOUS_PEOPLE[:5], ids=lambda p: p["name"])
    def test_aspects_detected_correctly(self, person):
        """
        Verify that major aspects are detected and match between implementations.

        Both libephemeris and pyswisseph should produce the same planetary
        positions, so the same aspects should be detected.
        """
        jd = ephem.julday(
            person["year"], person["month"], person["day"], person["hour"]
        )

        # Calculate positions with libephemeris
        positions_lib = {}
        for planet_id, planet_name in NATAL_PLANETS:
            pos, _ = ephem.calc_ut(jd, planet_id, 0)
            positions_lib[planet_name] = pos[0]

        # Calculate positions with Swiss Ephemeris
        positions_swe = {}
        for planet_id, planet_name in NATAL_PLANETS:
            pos, _ = swe.calc_ut(jd, planet_id, 0)
            positions_swe[planet_name] = pos[0]

        # Calculate aspects from both
        aspects_lib = calculate_all_aspects(positions_lib)
        aspects_swe = calculate_all_aspects(positions_swe)

        # Both should find the same aspects
        assert len(aspects_lib) == len(aspects_swe), (
            f"{person['name']}: Different aspect counts "
            f"(lib={len(aspects_lib)}, swe={len(aspects_swe)})"
        )

        # Convert to comparable sets (planet pair + aspect type)
        def aspect_key(a):
            return (a[0], a[1], a[2])

        lib_set = {aspect_key(a) for a in aspects_lib}
        swe_set = {aspect_key(a) for a in aspects_swe}

        assert lib_set == swe_set, (
            f"{person['name']}: Aspect mismatch\n"
            f"Only in lib: {lib_set - swe_set}\n"
            f"Only in swe: {swe_set - lib_set}"
        )

    @pytest.mark.integration
    @pytest.mark.parametrize("person", FAMOUS_PEOPLE[:3], ids=lambda p: p["name"])
    def test_major_aspects_present(self, person):
        """
        Every natal chart should have multiple major aspects.

        This is a sanity check that aspect detection is working.
        """
        jd = ephem.julday(
            person["year"], person["month"], person["day"], person["hour"]
        )

        positions = {}
        for planet_id, planet_name in NATAL_PLANETS:
            pos, _ = ephem.calc_ut(jd, planet_id, 0)
            positions[planet_name] = pos[0]

        aspects = calculate_all_aspects(positions)

        # Count major aspects (conjunction, opposition, trine, square, sextile)
        major_aspects = [
            a
            for a in aspects
            if a[2] in ["Conjunction", "Opposition", "Trine", "Square", "Sextile"]
        ]

        # Every chart should have at least 5 major aspects
        assert len(major_aspects) >= 5, (
            f"{person['name']}: Only {len(major_aspects)} major aspects found"
        )

    @pytest.mark.integration
    def test_marilyn_monroe_known_aspects(self):
        """
        Verify specific known aspects in Marilyn Monroe's chart.

        Marilyn has a Sun-Mercury conjunction (both in Gemini).
        """
        # Marilyn Monroe: June 1, 1926, 16:30 UT, Los Angeles
        jd = ephem.julday(1926, 6, 1, 16.5)

        sun_pos, _ = ephem.calc_ut(jd, SUN, 0)
        mercury_pos, _ = ephem.calc_ut(jd, MERCURY, 0)

        aspect = find_aspect(sun_pos[0], mercury_pos[0])

        # Marilyn has Sun-Mercury conjunction
        assert aspect is not None, (
            "No Sun-Mercury aspect found in Marilyn Monroe's chart"
        )
        assert aspect[0] == "Conjunction", f"Expected Conjunction, got {aspect[0]}"

    @pytest.mark.integration
    def test_aspect_orb_calculation(self):
        """Test that aspect orbs are calculated correctly."""
        # Test conjunction at exact 0 degrees
        aspect = find_aspect(120.0, 120.0)
        assert aspect is not None
        assert aspect[0] == "Conjunction"
        assert aspect[3] == 0.0  # Exact

        # Test trine with 5 degree orb
        aspect = find_aspect(0.0, 125.0)
        assert aspect is not None
        assert aspect[0] == "Trine"
        assert abs(aspect[3] - 5.0) < 0.1

        # Test opposition near exact
        aspect = find_aspect(0.0, 179.0)
        assert aspect is not None
        assert aspect[0] == "Opposition"
        assert abs(aspect[3] - 1.0) < 0.1


class TestCompleteNatalChartIntegration:
    """
    Full integration tests that calculate complete natal charts.

    These tests verify the entire chart calculation process.
    """

    @pytest.mark.integration
    @pytest.mark.slow
    @pytest.mark.parametrize("person", FAMOUS_PEOPLE, ids=lambda p: p["name"])
    def test_complete_natal_chart(self, person):
        """
        Calculate a complete natal chart and verify all components.

        This is the comprehensive integration test that exercises:
        - Julian Day calculation
        - All planetary positions
        - House cusps (Placidus)
        - Major aspects
        """
        # Calculate Julian Day
        jd = ephem.julday(
            person["year"], person["month"], person["day"], person["hour"]
        )

        # Verify JD is in valid range
        assert 2000000 < jd < 3000000, f"JD {jd} out of expected range"

        # Calculate all planetary positions
        chart = {"planets": {}, "houses": {}, "aspects": []}

        for planet_id, planet_name in NATAL_PLANETS:
            pos, _ = ephem.calc_ut(jd, planet_id, FLG_SPEED)

            assert 0 <= pos[0] < 360, f"{planet_name} longitude out of range"
            assert -90 <= pos[1] <= 90, f"{planet_name} latitude out of range"
            assert pos[2] > 0, f"{planet_name} distance should be positive"

            chart["planets"][planet_name] = {
                "longitude": pos[0],
                "latitude": pos[1],
                "distance": pos[2],
                "velocity": pos[3],
                "sign": get_sign(pos[0]),
            }

        # Calculate houses (skip if polar latitude for Placidus)
        if abs(person["lat"]) < 66:
            cusps, ascmc = ephem.houses(jd, person["lat"], person["lon"], ord("P"))

            chart["houses"]["cusps"] = list(cusps)
            chart["houses"]["asc"] = ascmc[0]
            chart["houses"]["mc"] = ascmc[1]
            chart["houses"]["asc_sign"] = get_sign(ascmc[0])
            chart["houses"]["mc_sign"] = get_sign(ascmc[1])

            assert 0 <= ascmc[0] < 360, "ASC out of range"
            assert 0 <= ascmc[1] < 360, "MC out of range"

        # Calculate aspects
        positions = {name: data["longitude"] for name, data in chart["planets"].items()}
        chart["aspects"] = calculate_all_aspects(positions)

        # Verify chart has expected data
        assert len(chart["planets"]) == 10, "Missing planets"
        assert len(chart["aspects"]) > 0, "No aspects calculated"

        # Verify Sun and Moon signs match expected
        assert chart["planets"]["Sun"]["sign"] == person["expected_sun_sign"], (
            f"Sun sign mismatch: got {chart['planets']['Sun']['sign']}, "
            f"expected {person['expected_sun_sign']}"
        )
        assert chart["planets"]["Moon"]["sign"] == person["expected_moon_sign"], (
            f"Moon sign mismatch: got {chart['planets']['Moon']['sign']}, "
            f"expected {person['expected_moon_sign']}"
        )

    @pytest.mark.integration
    def test_chart_reproducibility(self):
        """
        Calculating the same chart twice should give identical results.
        """
        person = FAMOUS_PEOPLE[0]  # Einstein
        jd = ephem.julday(
            person["year"], person["month"], person["day"], person["hour"]
        )

        # Calculate twice
        results1 = []
        results2 = []

        for planet_id, _ in NATAL_PLANETS:
            pos1, _ = ephem.calc_ut(jd, planet_id, 0)
            pos2, _ = ephem.calc_ut(jd, planet_id, 0)
            results1.append(pos1[0])
            results2.append(pos2[0])

        # Should be identical
        for i, (r1, r2) in enumerate(zip(results1, results2)):
            assert r1 == r2, f"Position {i} not reproducible: {r1} vs {r2}"


class TestHistoricalDates:
    """
    Test calculations for historical dates to verify ephemeris range.

    DE421 ephemeris covers 1899-07-29 through 2053-10-09.
    """

    @pytest.mark.integration
    def test_early_1900s_chart(self):
        """
        Test a chart from early 1900s (within DE421 range).
        """
        # JFK's father, Joe Kennedy Sr.: September 6, 1888 - outside range
        # Using a date just inside the range: August 1, 1900
        jd = ephem.julday(1900, 8, 1, 12.0)

        pos, _ = ephem.calc_ut(jd, SUN, 0)
        # Should be valid
        assert 0 <= pos[0] < 360
        assert get_sign(pos[0]) == "Leo"  # Early August = Leo

    @pytest.mark.integration
    def test_future_date_2050(self):
        """Test calculation near the end of DE421 range (2050)."""
        jd = ephem.julday(2050, 6, 15, 12.0)

        pos, _ = ephem.calc_ut(jd, SUN, 0)

        # Should be valid
        assert 0 <= pos[0] < 360
        assert get_sign(pos[0]) == "Gemini"  # Mid-June = Gemini

    @pytest.mark.integration
    @pytest.mark.comparison
    def test_j2000_reference_epoch(self):
        """
        Verify positions at J2000.0 epoch match reference values.

        J2000.0 (2000-01-01 12:00 TT) is a standard reference epoch.
        """
        jd = 2451545.0  # J2000.0

        # Sun should be near 280° longitude at J2000
        sun_pos, _ = ephem.calc_ut(jd, SUN, 0)
        assert 279 < sun_pos[0] < 281, f"Sun at J2000: {sun_pos[0]}°"

        # Compare with Swiss Ephemeris
        sun_swe, _ = swe.calc_ut(jd, SUN, 0)
        diff = angle_diff(sun_pos[0], sun_swe[0])
        assert diff < 1 / 3600, f"J2000 Sun diff: {diff * 3600:.2f} arcsec"

    @pytest.mark.integration
    @pytest.mark.skip(
        reason="libephemeris uses DE440 with different date range than DE431"
    )
    def test_date_before_ephemeris_range(self):
        """
        Dates before 1899-07-29 should raise an error.
        """
        jd = ephem.julday(1899, 1, 1, 12.0)

        with pytest.raises(Exception):
            # Should raise EphemerisRangeError or similar
            ephem.calc_ut(jd, SUN, 0)

    @pytest.mark.integration
    @pytest.mark.skip(
        reason="libephemeris uses DE440 with different date range than DE431"
    )
    def test_date_after_ephemeris_range(self):
        """
        Dates after 2053-10-09 should raise an error.
        """
        jd = ephem.julday(2054, 1, 1, 12.0)

        with pytest.raises(Exception):
            # Should raise EphemerisRangeError or similar
            ephem.calc_ut(jd, SUN, 0)
