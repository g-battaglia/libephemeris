"""
Aspect matching precision tests.

Tests verify that libephemeris and pyswisseph produce identical aspect counts
and matching for natal charts. This is a derived test that depends on planetary
position accuracy (fixes #5 Moon speed and #6 sidereal precision).

The key validation is that both implementations find the same aspects with
the same orbs, proving that planetary positions are sufficiently close.
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
    TRUE_NODE,
)


# ============================================================================
# TEST DATA - Specific persons from KR_FIX.md issue #8
# ============================================================================

ASPECT_TEST_PEOPLE = [
    {
        "name": "john_lennon",
        "year": 1940,
        "month": 10,
        "day": 9,
        "hour": 7.0,  # 07:00 UT
        "lat": 53.4084,
        "lon": -2.9916,
        "place": "Liverpool, UK",
    },
    {
        "name": "johnny_depp",
        "year": 1963,
        "month": 6,
        "day": 9,
        "hour": 3.75,  # 03:45 UT (adjusted from PST)
        "lat": 37.7749,
        "lon": -122.4194,
        "place": "San Francisco, CA",
    },
    {
        "name": "paul_mccartney",
        "year": 1942,
        "month": 6,
        "day": 18,
        "hour": 14.5,  # 14:30 UT
        "lat": 53.4084,
        "lon": -2.9916,
        "place": "Liverpool, UK",
    },
]

# Major aspects used in traditional astrology (subset for testing)
MAJOR_ASPECTS = [
    ("Conjunction", 0, 8),
    ("Opposition", 180, 8),
    ("Trine", 120, 8),
    ("Square", 90, 8),
    ("Sextile", 60, 6),
]

# Planets for aspect calculation (excluding Chiron to avoid SPK requirement)
ASPECT_PLANETS = [
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
    (TRUE_NODE, "TrueNode"),
]


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def normalize_angle(angle: float) -> float:
    """Normalize angle to 0-360 range."""
    angle = angle % 360
    if angle < 0:
        angle += 360
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

    for aspect_name, exact_angle, orb in MAJOR_ASPECTS:
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


class TestNatalAspectsCount:
    """Test that aspect counts match between libephemeris and Swiss Ephemeris."""

    @pytest.mark.integration
    @pytest.mark.comparison
    @pytest.mark.parametrize("person", ASPECT_TEST_PEOPLE, ids=lambda p: p["name"])
    def test_natal_aspects_count(self, person):
        """
        Verify aspect counts match between libephemeris and Swiss Ephemeris.

        This test validates that planetary positions are accurate enough that
        aspects at the edge of orb tolerances are detected identically.
        """
        jd = ephem.julday(
            person["year"], person["month"], person["day"], person["hour"]
        )

        # Calculate positions with libephemeris
        positions_lib = {}
        for planet_id, planet_name in ASPECT_PLANETS:
            pos, _ = ephem.calc_ut(jd, planet_id, 0)
            positions_lib[planet_name] = pos[0]

        # Calculate positions with Swiss Ephemeris
        positions_swe = {}
        for planet_id, planet_name in ASPECT_PLANETS:
            pos, _ = swe.calc_ut(jd, planet_id, 0)
            positions_swe[planet_name] = pos[0]

        # Calculate aspects
        aspects_lib = calculate_all_aspects(positions_lib)
        aspects_swe = calculate_all_aspects(positions_swe)

        assert len(aspects_lib) == len(aspects_swe), (
            f"{person['name']}: Aspect count mismatch "
            f"(lib={len(aspects_lib)}, swe={len(aspects_swe)})"
        )


class TestNatalAspectsMatch:
    """Test that exact aspect sets match between implementations."""

    @pytest.mark.integration
    @pytest.mark.comparison
    @pytest.mark.parametrize("person", ASPECT_TEST_PEOPLE, ids=lambda p: p["name"])
    def test_natal_aspects_match(self, person):
        """
        Verify that the exact same aspects are detected by both implementations.

        Each aspect is identified by the pair of planets and the aspect type.
        """
        jd = ephem.julday(
            person["year"], person["month"], person["day"], person["hour"]
        )

        # Calculate positions with libephemeris
        positions_lib = {}
        for planet_id, planet_name in ASPECT_PLANETS:
            pos, _ = ephem.calc_ut(jd, planet_id, 0)
            positions_lib[planet_name] = pos[0]

        # Calculate positions with Swiss Ephemeris
        positions_swe = {}
        for planet_id, planet_name in ASPECT_PLANETS:
            pos, _ = swe.calc_ut(jd, planet_id, 0)
            positions_swe[planet_name] = pos[0]

        # Calculate aspects
        aspects_lib = calculate_all_aspects(positions_lib)
        aspects_swe = calculate_all_aspects(positions_swe)

        # Convert to comparable sets (planet pair + aspect type)
        def aspect_key(a):
            return (a[0], a[1], a[2])

        lib_set = {aspect_key(a) for a in aspects_lib}
        swe_set = {aspect_key(a) for a in aspects_swe}

        assert lib_set == swe_set, (
            f"{person['name']}: Aspect mismatch\n"
            f"Only in libephemeris: {lib_set - swe_set}\n"
            f"Only in Swiss Ephemeris: {swe_set - lib_set}"
        )


class TestNatalAspectByIndex:
    """Test that aspects can be matched by index (same ordering)."""

    @pytest.mark.integration
    @pytest.mark.comparison
    @pytest.mark.parametrize("person", ASPECT_TEST_PEOPLE, ids=lambda p: p["name"])
    def test_aspect_orbs_match(self, person):
        """
        Verify that aspect orbs are very close between implementations.

        The orb difference should be less than 1 arcsecond for each aspect.
        """
        jd = ephem.julday(
            person["year"], person["month"], person["day"], person["hour"]
        )

        # Calculate positions with libephemeris
        positions_lib = {}
        for planet_id, planet_name in ASPECT_PLANETS:
            pos, _ = ephem.calc_ut(jd, planet_id, 0)
            positions_lib[planet_name] = pos[0]

        # Calculate positions with Swiss Ephemeris
        positions_swe = {}
        for planet_id, planet_name in ASPECT_PLANETS:
            pos, _ = swe.calc_ut(jd, planet_id, 0)
            positions_swe[planet_name] = pos[0]

        # Calculate aspects
        aspects_lib = calculate_all_aspects(positions_lib)
        aspects_swe = calculate_all_aspects(positions_swe)

        # Match aspects by planet pair and type
        lib_dict = {(a[0], a[1], a[2]): a[3] for a in aspects_lib}
        swe_dict = {(a[0], a[1], a[2]): a[3] for a in aspects_swe}

        for key in lib_dict:
            if key in swe_dict:
                orb_diff = abs(lib_dict[key] - swe_dict[key])
                # TrueNode has larger position differences (up to 15 arcsec)
                # Other planets should be within 2 arcseconds
                if "TrueNode" in key[0] or "TrueNode" in key[1]:
                    tolerance = 20 / 3600  # 20 arcseconds for node aspects
                else:
                    tolerance = 2 / 3600  # 2 arcseconds for planet aspects
                assert orb_diff < tolerance, (
                    f"{person['name']}: {key} orb diff {orb_diff * 3600:.2f} arcsec"
                )


class TestAspectMovement:
    """Test that aspects behave correctly as time progresses."""

    @pytest.mark.integration
    @pytest.mark.parametrize("person", ASPECT_TEST_PEOPLE, ids=lambda p: p["name"])
    def test_aspect_stability(self, person):
        """
        Verify aspects are stable over small time increments.

        An aspect found at time T should also be found at T +/- 1 hour
        (unless it's at the very edge of the orb).
        """
        jd = ephem.julday(
            person["year"], person["month"], person["day"], person["hour"]
        )

        # Calculate positions at base time
        positions_base = {}
        for planet_id, planet_name in ASPECT_PLANETS:
            pos, _ = ephem.calc_ut(jd, planet_id, 0)
            positions_base[planet_name] = pos[0]

        aspects_base = calculate_all_aspects(positions_base)

        # Check aspects at +1 hour
        jd_plus1 = jd + 1 / 24
        positions_plus1 = {}
        for planet_id, planet_name in ASPECT_PLANETS:
            pos, _ = ephem.calc_ut(jd_plus1, planet_id, 0)
            positions_plus1[planet_name] = pos[0]

        aspects_plus1 = calculate_all_aspects(positions_plus1)

        # Most aspects should remain (only tight orbs might change)
        base_set = {(a[0], a[1], a[2]) for a in aspects_base}
        plus1_set = {(a[0], a[1], a[2]) for a in aspects_plus1}

        # At least 80% of aspects should remain stable
        stable_count = len(base_set & plus1_set)
        stability_ratio = stable_count / len(base_set) if base_set else 1.0

        assert stability_ratio >= 0.8, (
            f"{person['name']}: Aspect stability {stability_ratio:.0%} < 80%"
        )


class TestSynastryAspects:
    """Test synastry aspect calculations between two charts."""

    @pytest.mark.integration
    @pytest.mark.comparison
    def test_synastry_aspects_john_paul(self):
        """
        Test synastry aspects between John Lennon and Paul McCartney.

        Both implementations should find the same inter-chart aspects.
        """
        john = ASPECT_TEST_PEOPLE[0]  # john_lennon
        paul = ASPECT_TEST_PEOPLE[2]  # paul_mccartney

        jd_john = ephem.julday(
            john["year"], john["month"], john["day"], john["hour"]
        )
        jd_paul = ephem.julday(
            paul["year"], paul["month"], paul["day"], paul["hour"]
        )

        # Calculate positions with libephemeris
        positions_john_lib = {}
        positions_paul_lib = {}
        for planet_id, planet_name in ASPECT_PLANETS:
            pos, _ = ephem.calc_ut(jd_john, planet_id, 0)
            positions_john_lib[planet_name] = pos[0]
            pos, _ = ephem.calc_ut(jd_paul, planet_id, 0)
            positions_paul_lib[planet_name] = pos[0]

        # Calculate positions with Swiss Ephemeris
        positions_john_swe = {}
        positions_paul_swe = {}
        for planet_id, planet_name in ASPECT_PLANETS:
            pos, _ = swe.calc_ut(jd_john, planet_id, 0)
            positions_john_swe[planet_name] = pos[0]
            pos, _ = swe.calc_ut(jd_paul, planet_id, 0)
            positions_paul_swe[planet_name] = pos[0]

        # Calculate synastry aspects
        def synastry_aspects(pos1, pos2):
            aspects = []
            for p1_name, p1_lon in pos1.items():
                for p2_name, p2_lon in pos2.items():
                    aspect = find_aspect(p1_lon, p2_lon)
                    if aspect:
                        aspects.append(
                            (f"John_{p1_name}", f"Paul_{p2_name}", aspect[0])
                        )
            return set(aspects)

        lib_synastry = synastry_aspects(positions_john_lib, positions_paul_lib)
        swe_synastry = synastry_aspects(positions_john_swe, positions_paul_swe)

        assert lib_synastry == swe_synastry, (
            f"Synastry mismatch:\n"
            f"Only in lib: {lib_synastry - swe_synastry}\n"
            f"Only in swe: {swe_synastry - lib_synastry}"
        )

    @pytest.mark.integration
    @pytest.mark.comparison
    def test_synastry_aspect_count(self):
        """Test that synastry produces a reasonable number of aspects."""
        john = ASPECT_TEST_PEOPLE[0]
        paul = ASPECT_TEST_PEOPLE[2]

        jd_john = ephem.julday(
            john["year"], john["month"], john["day"], john["hour"]
        )
        jd_paul = ephem.julday(
            paul["year"], paul["month"], paul["day"], paul["hour"]
        )

        # Calculate positions
        positions_john = {}
        positions_paul = {}
        for planet_id, planet_name in ASPECT_PLANETS:
            pos, _ = ephem.calc_ut(jd_john, planet_id, 0)
            positions_john[planet_name] = pos[0]
            pos, _ = ephem.calc_ut(jd_paul, planet_id, 0)
            positions_paul[planet_name] = pos[0]

        # Calculate synastry aspects
        synastry_count = 0
        for p1_name, p1_lon in positions_john.items():
            for p2_name, p2_lon in positions_paul.items():
                if find_aspect(p1_lon, p2_lon):
                    synastry_count += 1

        # Two charts with 11 planets each should have many synastry aspects
        # 11 * 11 = 121 possible pairs, expect at least 20 aspects
        assert synastry_count >= 20, (
            f"Only {synastry_count} synastry aspects found, expected >= 20"
        )


class TestDifdeg2nConsistency:
    """Test that difdeg2n is consistent with pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "p1,p2",
        [
            (180, 0),
            (0, 180),
            (270, 90),
            (90, 270),
            (10, 20),
            (350, 10),
            (10, 350),
            (0.0001, 359.9999),
        ],
    )
    def test_difdeg2n_matches_swisseph(self, p1, p2):
        """Test that difdeg2n matches pyswisseph exactly."""
        lib_result = ephem.difdeg2n(p1, p2)
        swe_result = swe.difdeg2n(p1, p2)

        assert lib_result == pytest.approx(swe_result, abs=1e-10), (
            f"difdeg2n({p1}, {p2}): lib={lib_result}, swe={swe_result}"
        )

    def test_difdeg2n_180_returns_negative(self):
        """Test that difdeg2n returns -180 for 180 degree separation."""
        # This is the pyswisseph behavior
        assert ephem.difdeg2n(180, 0) == pytest.approx(-180.0)
        assert ephem.difdeg2n(0, 180) == pytest.approx(-180.0)
        assert ephem.difdeg2n(270, 90) == pytest.approx(-180.0)
