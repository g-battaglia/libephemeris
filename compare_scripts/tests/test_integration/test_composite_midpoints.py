"""
Composite midpoint precision tests.

Tests verify that libephemeris and pyswisseph produce composite midpoints
with acceptable precision. This is a derived test that depends on planetary
position accuracy (fixes #5 Moon speed, #6 sidereal precision, and #7
perspective precision).

The key validation is that composite midpoints calculated from individual
chart positions match between both implementations, proving that the underlying
planetary positions are sufficiently accurate.

Tolerance notes:
- Most planets: 0.0003° (about 1 arcsec) - achievable with current precision
- TrueNode: 0.002° (about 7 arcsec) - has larger inherent calculation differences
  between implementations due to orbital element model differences (same tolerance
  as in test_aspect_matching.py which uses 20 arcsec for TrueNode aspects)

Reference: KR_FIX.md issue #9
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
# TEST DATA - Pairs of people for composite chart testing
# ============================================================================

COMPOSITE_TEST_PAIRS = [
    {
        "name": "lennon_mccartney",
        "person1": {
            "name": "john_lennon",
            "year": 1940,
            "month": 10,
            "day": 9,
            "hour": 7.0,
            "lat": 53.4084,
            "lon": -2.9916,
            "place": "Liverpool, UK",
        },
        "person2": {
            "name": "paul_mccartney",
            "year": 1942,
            "month": 6,
            "day": 18,
            "hour": 14.5,
            "lat": 53.4084,
            "lon": -2.9916,
            "place": "Liverpool, UK",
        },
    },
    {
        "name": "depp_mccartney",
        "person1": {
            "name": "johnny_depp",
            "year": 1963,
            "month": 6,
            "day": 9,
            "hour": 3.75,
            "lat": 37.7749,
            "lon": -122.4194,
            "place": "San Francisco, CA",
        },
        "person2": {
            "name": "paul_mccartney",
            "year": 1942,
            "month": 6,
            "day": 18,
            "hour": 14.5,
            "lat": 53.4084,
            "lon": -2.9916,
            "place": "Liverpool, UK",
        },
    },
    {
        "name": "lennon_depp",
        "person1": {
            "name": "john_lennon",
            "year": 1940,
            "month": 10,
            "day": 9,
            "hour": 7.0,
            "lat": 53.4084,
            "lon": -2.9916,
            "place": "Liverpool, UK",
        },
        "person2": {
            "name": "johnny_depp",
            "year": 1963,
            "month": 6,
            "day": 9,
            "hour": 3.75,
            "lat": 37.7749,
            "lon": -122.4194,
            "place": "San Francisco, CA",
        },
    },
]

# Planets for composite calculation (excluding TrueNode for main tests)
COMPOSITE_PLANETS = [
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

# TrueNode tested separately with relaxed tolerance
TRUE_NODE = (TRUE_NODE, "TrueNode")

# Planets most affected by precision issues per KR_FIX.md #9
SENSITIVE_PLANETS = ["Moon", "Jupiter", "Uranus", "Sun"]


# ============================================================================
# TOLERANCES
# ============================================================================

# Standard tolerance for planets: 0.0003° (about 1 arcsec)
# This accounts for the composite midpoint inheriting half of each position's error
PLANET_TOLERANCE = 0.0003

# Relaxed tolerance for TrueNode: 0.002° (about 7 arcsec)
# TrueNode calculations have larger inherent differences between implementations
# due to different orbital element models (consistent with test_aspect_matching.py)
TRUE_NODE_TOLERANCE = 0.002


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def deg_midpoint(a: float, b: float) -> float:
    """
    Calculate the midpoint between two angles in degrees.

    Uses the shorter arc (same algorithm as libephemeris.deg_midp).
    """
    a = a % 360.0
    b = b % 360.0

    diff = b - a
    if diff > 180.0:
        diff -= 360.0
    elif diff < -180.0:
        diff += 360.0

    midp = a + diff / 2.0
    return midp % 360.0


def angle_diff(a1: float, a2: float) -> float:
    """Calculate the smallest angular difference between two angles."""
    a1 = a1 % 360.0
    a2 = a2 % 360.0
    diff = abs(a1 - a2)
    if diff > 180:
        diff = 360 - diff
    return diff


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestCompositeExpectedData:
    """
    Test that composite midpoints match between libephemeris and Swiss Ephemeris.

    These tests verify fix #9 from KR_FIX.md: Composite midpoint precision.
    Uses tiered tolerance: 0.0003° for planets, 0.002° for TrueNode.
    """

    @pytest.mark.integration
    @pytest.mark.comparison
    @pytest.mark.parametrize("pair", COMPOSITE_TEST_PAIRS, ids=lambda p: p["name"])
    def test_composite_positions_match(self, pair):
        """
        Verify composite midpoint positions match between implementations.

        Each planet's composite position is the midpoint of the two natal positions.
        Both implementations should calculate the same midpoints within tolerance.
        """
        p1, p2 = pair["person1"], pair["person2"]

        jd1 = ephem.julday(p1["year"], p1["month"], p1["day"], p1["hour"])
        jd2 = ephem.julday(p2["year"], p2["month"], p2["day"], p2["hour"])

        errors = []

        for planet_id, planet_name in COMPOSITE_PLANETS:
            # Calculate positions with libephemeris
            pos1_lib, _ = ephem.calc_ut(jd1, planet_id, 0)
            pos2_lib, _ = ephem.calc_ut(jd2, planet_id, 0)

            # Calculate positions with Swiss Ephemeris
            pos1_swe, _ = swe.calc_ut(jd1, planet_id, 0)
            pos2_swe, _ = swe.calc_ut(jd2, planet_id, 0)

            # Calculate composite midpoints
            composite_lib = deg_midpoint(pos1_lib[0], pos2_lib[0])
            composite_swe = deg_midpoint(pos1_swe[0], pos2_swe[0])

            # Compare
            diff = angle_diff(composite_lib, composite_swe)

            if diff > PLANET_TOLERANCE:
                errors.append(
                    f"{planet_name}: lib={composite_lib:.6f}, swe={composite_swe:.6f}, "
                    f"diff={diff:.6f} ({diff * 3600:.2f} arcsec)"
                )

        assert not errors, (
            f"{pair['name']}: Composite midpoint precision exceeded {PLANET_TOLERANCE}deg:\n"
            + "\n".join(errors)
        )

    @pytest.mark.integration
    @pytest.mark.comparison
    @pytest.mark.parametrize("pair", COMPOSITE_TEST_PAIRS, ids=lambda p: p["name"])
    def test_true_node_composite_precision(self, pair):
        """
        Test TrueNode composite precision with relaxed tolerance.

        TrueNode has larger inherent differences between implementations due to
        different orbital element models (consistent with test_aspect_matching.py).
        """
        p1, p2 = pair["person1"], pair["person2"]

        jd1 = ephem.julday(p1["year"], p1["month"], p1["day"], p1["hour"])
        jd2 = ephem.julday(p2["year"], p2["month"], p2["day"], p2["hour"])

        planet_id, planet_name = TRUE_NODE

        # Calculate positions with libephemeris
        pos1_lib, _ = ephem.calc_ut(jd1, planet_id, 0)
        pos2_lib, _ = ephem.calc_ut(jd2, planet_id, 0)

        # Calculate positions with Swiss Ephemeris
        pos1_swe, _ = swe.calc_ut(jd1, planet_id, 0)
        pos2_swe, _ = swe.calc_ut(jd2, planet_id, 0)

        # Calculate composite midpoints
        composite_lib = deg_midpoint(pos1_lib[0], pos2_lib[0])
        composite_swe = deg_midpoint(pos1_swe[0], pos2_swe[0])

        diff = angle_diff(composite_lib, composite_swe)

        assert diff <= TRUE_NODE_TOLERANCE, (
            f"{pair['name']} {planet_name}: diff={diff:.6f} ({diff * 3600:.2f} arcsec) "
            f"> tolerance {TRUE_NODE_TOLERANCE}deg"
        )

    @pytest.mark.integration
    @pytest.mark.comparison
    @pytest.mark.parametrize("pair", COMPOSITE_TEST_PAIRS, ids=lambda p: p["name"])
    def test_sensitive_planets_precision(self, pair):
        """
        Test precision for planets most affected by precision issues.

        Moon, Jupiter, Uranus, and Sun are mentioned in KR_FIX.md #9 as the
        bodies most affected by composite precision issues.
        """
        p1, p2 = pair["person1"], pair["person2"]

        jd1 = ephem.julday(p1["year"], p1["month"], p1["day"], p1["hour"])
        jd2 = ephem.julday(p2["year"], p2["month"], p2["day"], p2["hour"])

        max_diff = 0.0
        max_planet = None

        for planet_id, planet_name in COMPOSITE_PLANETS:
            if planet_name not in SENSITIVE_PLANETS:
                continue

            # Calculate positions with libephemeris
            pos1_lib, _ = ephem.calc_ut(jd1, planet_id, 0)
            pos2_lib, _ = ephem.calc_ut(jd2, planet_id, 0)

            # Calculate positions with Swiss Ephemeris
            pos1_swe, _ = swe.calc_ut(jd1, planet_id, 0)
            pos2_swe, _ = swe.calc_ut(jd2, planet_id, 0)

            # Calculate composite midpoints
            composite_lib = deg_midpoint(pos1_lib[0], pos2_lib[0])
            composite_swe = deg_midpoint(pos1_swe[0], pos2_swe[0])

            diff = angle_diff(composite_lib, composite_swe)
            if diff > max_diff:
                max_diff = diff
                max_planet = planet_name

        assert max_diff <= PLANET_TOLERANCE, (
            f"{pair['name']}: Sensitive planet {max_planet} exceeded tolerance: "
            f"diff={max_diff:.6f} ({max_diff * 3600:.2f} arcsec)"
        )

    @pytest.mark.integration
    @pytest.mark.comparison
    @pytest.mark.parametrize("pair", COMPOSITE_TEST_PAIRS, ids=lambda p: p["name"])
    def test_composite_via_deg_midp_function(self, pair):
        """
        Test composite calculation using libephemeris.deg_midp directly.

        Verifies that the deg_midp function produces results matching
        Swiss Ephemeris-derived midpoints.
        """
        p1, p2 = pair["person1"], pair["person2"]

        jd1 = ephem.julday(p1["year"], p1["month"], p1["day"], p1["hour"])
        jd2 = ephem.julday(p2["year"], p2["month"], p2["day"], p2["hour"])

        for planet_id, planet_name in COMPOSITE_PLANETS:
            # Calculate positions with libephemeris
            pos1_lib, _ = ephem.calc_ut(jd1, planet_id, 0)
            pos2_lib, _ = ephem.calc_ut(jd2, planet_id, 0)

            # Use libephemeris deg_midp function
            composite_lib = ephem.deg_midp(pos1_lib[0], pos2_lib[0])

            # Calculate with Swiss Ephemeris and manual midpoint
            pos1_swe, _ = swe.calc_ut(jd1, planet_id, 0)
            pos2_swe, _ = swe.calc_ut(jd2, planet_id, 0)
            composite_swe = deg_midpoint(pos1_swe[0], pos2_swe[0])

            diff = angle_diff(composite_lib, composite_swe)

            assert diff <= PLANET_TOLERANCE, (
                f"{pair['name']} {planet_name}: deg_midp diff={diff:.6f} "
                f"({diff * 3600:.2f} arcsec) > {PLANET_TOLERANCE}"
            )


class TestCompositeStatistics:
    """Statistical tests for composite midpoint precision across all pairs."""

    @pytest.mark.integration
    @pytest.mark.comparison
    def test_overall_composite_precision(self):
        """
        Verify overall composite precision statistics.

        All composite midpoints should be within tolerance (0.0003° for planets,
        0.002° for TrueNode), with expected precision around 0.0001° (0.36 arcsec).
        """
        planet_diffs = []
        node_diffs = []

        for pair in COMPOSITE_TEST_PAIRS:
            p1, p2 = pair["person1"], pair["person2"]

            jd1 = ephem.julday(p1["year"], p1["month"], p1["day"], p1["hour"])
            jd2 = ephem.julday(p2["year"], p2["month"], p2["day"], p2["hour"])

            # Test planets
            for planet_id, planet_name in COMPOSITE_PLANETS:
                pos1_lib, _ = ephem.calc_ut(jd1, planet_id, 0)
                pos2_lib, _ = ephem.calc_ut(jd2, planet_id, 0)

                pos1_swe, _ = swe.calc_ut(jd1, planet_id, 0)
                pos2_swe, _ = swe.calc_ut(jd2, planet_id, 0)

                composite_lib = deg_midpoint(pos1_lib[0], pos2_lib[0])
                composite_swe = deg_midpoint(pos1_swe[0], pos2_swe[0])

                diff = angle_diff(composite_lib, composite_swe)
                planet_diffs.append((pair["name"], planet_name, diff))

            # Test TrueNode
            planet_id, planet_name = TRUE_NODE
            pos1_lib, _ = ephem.calc_ut(jd1, planet_id, 0)
            pos2_lib, _ = ephem.calc_ut(jd2, planet_id, 0)
            pos1_swe, _ = swe.calc_ut(jd1, planet_id, 0)
            pos2_swe, _ = swe.calc_ut(jd2, planet_id, 0)
            composite_lib = deg_midpoint(pos1_lib[0], pos2_lib[0])
            composite_swe = deg_midpoint(pos1_swe[0], pos2_swe[0])
            diff = angle_diff(composite_lib, composite_swe)
            node_diffs.append((pair["name"], planet_name, diff))

        # Analyze planet results
        planet_max_diff = max(d[2] for d in planet_diffs)
        planet_avg_diff = sum(d[2] for d in planet_diffs) / len(planet_diffs)
        planet_under_tolerance = sum(
            1 for d in planet_diffs if d[2] <= PLANET_TOLERANCE
        )

        # All planets should be under tolerance
        assert planet_under_tolerance == len(planet_diffs), (
            f"Only {planet_under_tolerance}/{len(planet_diffs)} planets under "
            f"{PLANET_TOLERANCE}deg tolerance.\n"
            f"Max diff: {planet_max_diff:.6f} ({planet_max_diff * 3600:.2f} arcsec)\n"
            f"Avg diff: {planet_avg_diff:.6f} ({planet_avg_diff * 3600:.2f} arcsec)"
        )

        # All TrueNodes should be under their tolerance
        node_under_tolerance = sum(1 for d in node_diffs if d[2] <= TRUE_NODE_TOLERANCE)
        assert node_under_tolerance == len(node_diffs), (
            f"Only {node_under_tolerance}/{len(node_diffs)} TrueNodes under "
            f"{TRUE_NODE_TOLERANCE}deg tolerance."
        )

        # Log statistics for reference
        print("\nComposite precision statistics:")
        print(f"  Planets tested: {len(planet_diffs)}")
        print(
            f"  Planet max diff: {planet_max_diff:.6f} deg ({planet_max_diff * 3600:.2f} arcsec)"
        )
        print(
            f"  Planet avg diff: {planet_avg_diff:.6f} deg ({planet_avg_diff * 3600:.2f} arcsec)"
        )
        print(
            f"  All planets under {PLANET_TOLERANCE}deg: {planet_under_tolerance == len(planet_diffs)}"
        )

    @pytest.mark.integration
    @pytest.mark.comparison
    def test_no_systematic_bias(self):
        """
        Verify no significant systematic bias in composite calculations.

        The differences should be approximately distributed around zero.
        """
        signed_diffs = []

        for pair in COMPOSITE_TEST_PAIRS:
            p1, p2 = pair["person1"], pair["person2"]

            jd1 = ephem.julday(p1["year"], p1["month"], p1["day"], p1["hour"])
            jd2 = ephem.julday(p2["year"], p2["month"], p2["day"], p2["hour"])

            for planet_id, planet_name in COMPOSITE_PLANETS:
                pos1_lib, _ = ephem.calc_ut(jd1, planet_id, 0)
                pos2_lib, _ = ephem.calc_ut(jd2, planet_id, 0)

                pos1_swe, _ = swe.calc_ut(jd1, planet_id, 0)
                pos2_swe, _ = swe.calc_ut(jd2, planet_id, 0)

                composite_lib = deg_midpoint(pos1_lib[0], pos2_lib[0])
                composite_swe = deg_midpoint(pos1_swe[0], pos2_swe[0])

                # Signed difference (lib - swe)
                # Handle wraparound
                diff = composite_lib - composite_swe
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360

                signed_diffs.append(diff)

        # Count positive vs negative differences
        positive = sum(1 for d in signed_diffs if d > 0)
        negative = sum(1 for d in signed_diffs if d < 0)
        total = len(signed_diffs)

        # Mean should be small (within 0.0001 degree)
        mean_diff = sum(signed_diffs) / len(signed_diffs)

        assert abs(mean_diff) < 0.0001, (
            f"Systematic bias detected: mean={mean_diff:.7f} deg\n"
            f"Positive: {positive}/{total}, Negative: {negative}/{total}"
        )
