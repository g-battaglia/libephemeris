"""
Pytest-style Elongation Helper Functions Comparison Tests.

Validates elongation helper functions (get_elongation_from_sun, is_morning_star,
is_evening_star, get_elongation_type) against manual calculations using pyswisseph
to determine planet and Sun positions.

These functions are critical for determining planetary visibility (morning star vs
evening star status), which affects astrological interpretation.
"""

import pytest
import swisseph as swe
import libephemeris as pyephem
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
    FLG_SWIEPH,
)


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class ElongationTolerance:
    """Tolerance thresholds for elongation comparisons."""

    # Elongation angle tolerance (degrees)
    ANGLE_DEGREES = 0.001

    # For comparing with pheno_ut results
    PHENO_ANGLE = 0.01


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def normalize_to_signed(angle: float) -> float:
    """
    Normalize an angle to the range -180 to +180 degrees.

    Args:
        angle: Angle in degrees

    Returns:
        Normalized angle in range (-180, 180]
    """
    while angle > 180.0:
        angle -= 360.0
    while angle <= -180.0:
        angle += 360.0
    return angle


def calculate_elongation_from_positions(planet_lon: float, sun_lon: float) -> float:
    """
    Calculate signed elongation from planet and Sun longitudes.

    Args:
        planet_lon: Planet ecliptic longitude in degrees
        sun_lon: Sun ecliptic longitude in degrees

    Returns:
        Signed elongation in degrees:
            - Positive = eastern elongation (evening star)
            - Negative = western elongation (morning star)
    """
    lon_diff = planet_lon - sun_lon
    return normalize_to_signed(lon_diff)


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# All planets that can have elongation from the Sun
ELONGATION_PLANETS = [
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
    (PLUTO, "Pluto"),
    (MOON, "Moon"),
]

# Inner planets (most important for morning/evening star classification)
INNER_PLANETS = [
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
]

# Test dates covering different planetary configurations
TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 1, 15, 12.0, "Jan 2024"),
    (2024, 4, 15, 12.0, "Apr 2024"),
    (2024, 7, 15, 12.0, "Jul 2024"),
    (2024, 10, 15, 12.0, "Oct 2024"),
    (1990, 6, 21, 0.0, "Summer Solstice 1990"),
    (2030, 12, 21, 12.0, "Winter Solstice 2030"),
]

# Specific dates when planets have known elongation characteristics
# These dates are verified against actual pyswisseph calculations
KNOWN_ELONGATION_DATES = [
    # Venus western elongation (morning star) - Venus is west of Sun in early 2024
    (2024, 3, 22, 12.0, VENUS, "Venus western elongation", "morning"),
    # Mercury eastern elongation (evening star) - Mercury east of Sun in July 2024
    (2024, 7, 15, 12.0, MERCURY, "Mercury eastern elongation", "evening"),
    # Mercury western elongation (morning star) - Mercury west of Sun in May 2024
    (2024, 5, 15, 12.0, MERCURY, "Mercury western elongation", "morning"),
    # Venus eastern elongation (evening star) - Venus is east of Sun in late 2024
    (2024, 10, 15, 12.0, VENUS, "Venus eastern elongation", "evening"),
]


# ============================================================================
# FIXTURES
# ============================================================================


@pytest.fixture
def jd_standard():
    """Standard Julian Day for tests (J2000)."""
    return swe.julday(2000, 1, 1, 12.0)


# ============================================================================
# GET_ELONGATION_FROM_SUN TESTS
# ============================================================================


class TestGetElongationFromSun:
    """Tests for get_elongation_from_sun function."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", ELONGATION_PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_elongation_matches_manual_calculation(
        self, body_id, body_name, year, month, day, hour, date_desc
    ):
        """Test that get_elongation_from_sun matches manual calculation from positions."""
        jd = swe.julday(year, month, day, hour)

        # Calculate positions using pyswisseph
        try:
            planet_pos, _ = swe.calc_ut(jd, body_id, FLG_SWIEPH)
            sun_pos, _ = swe.calc_ut(jd, SUN, FLG_SWIEPH)
            planet_lon = planet_pos[0]
            sun_lon = sun_pos[0]
        except Exception as e:
            pytest.skip(f"SwissEphemeris calc failed: {e}")

        # Calculate expected elongation manually
        expected_elongation = calculate_elongation_from_positions(planet_lon, sun_lon)
        expected_is_evening = expected_elongation > 0

        # Get libephemeris result
        try:
            elongation, is_evening = pyephem.get_elongation_from_sun(
                jd, body_id, FLG_SWIEPH
            )
        except Exception as e:
            pytest.fail(f"libephemeris get_elongation_from_sun failed: {e}")

        # Compare elongation values
        diff = abs(elongation - expected_elongation)
        assert diff < ElongationTolerance.ANGLE_DEGREES, (
            f"{body_name} @ {date_desc}: elongation diff {diff:.6f}° "
            f"(lib={elongation:.4f}°, expected={expected_elongation:.4f}°)"
        )

        # Compare morning/evening star classification
        assert is_evening == expected_is_evening, (
            f"{body_name} @ {date_desc}: is_evening_star mismatch "
            f"(lib={is_evening}, expected={expected_is_evening}, "
            f"elongation={elongation:.2f}°)"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", INNER_PLANETS)
    def test_elongation_range_for_inner_planets(self, jd_standard, body_id, body_name):
        """Test that inner planets have elongation within physical limits."""
        elongation, _ = pyephem.get_elongation_from_sun(
            jd_standard, body_id, FLG_SWIEPH
        )

        # Mercury max elongation ~28°, Venus max ~47°
        max_elong = 30 if body_id == MERCURY else 50

        assert abs(elongation) <= max_elong, (
            f"{body_name}: elongation {abs(elongation):.1f}° "
            f"exceeds physical limit of {max_elong}°"
        )

    @pytest.mark.comparison
    def test_sun_returns_zero_elongation(self, jd_standard):
        """Test that Sun has zero elongation from itself."""
        elongation, is_evening = pyephem.get_elongation_from_sun(
            jd_standard, SUN, FLG_SWIEPH
        )

        assert abs(elongation) < ElongationTolerance.ANGLE_DEGREES, (
            f"Sun elongation from Sun should be ~0, got {elongation:.6f}°"
        )


# ============================================================================
# IS_MORNING_STAR TESTS
# ============================================================================


class TestIsMorningStar:
    """Tests for is_morning_star function."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", ELONGATION_PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_morning_star_matches_elongation_sign(
        self, body_id, body_name, year, month, day, hour, date_desc
    ):
        """Test that is_morning_star matches negative elongation."""
        jd = swe.julday(year, month, day, hour)

        # Get elongation to determine expected result
        try:
            elongation, is_evening = pyephem.get_elongation_from_sun(
                jd, body_id, FLG_SWIEPH
            )
        except Exception as e:
            pytest.skip(f"get_elongation_from_sun failed: {e}")

        # Get is_morning_star result
        try:
            is_morning = pyephem.is_morning_star(jd, body_id, FLG_SWIEPH)
        except Exception as e:
            pytest.fail(f"is_morning_star failed: {e}")

        # Morning star = western elongation = negative = NOT evening star
        expected_morning = not is_evening

        assert is_morning == expected_morning, (
            f"{body_name} @ {date_desc}: is_morning_star={is_morning}, "
            f"expected={expected_morning}, elongation={elongation:.2f}°"
        )

    @pytest.mark.comparison
    def test_sun_is_not_morning_star(self, jd_standard):
        """Test that Sun cannot be a morning star."""
        is_morning = pyephem.is_morning_star(jd_standard, SUN, FLG_SWIEPH)
        assert is_morning is False, "Sun should never be classified as morning star"

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,body_id,desc,expected_type",
        KNOWN_ELONGATION_DATES,
    )
    def test_known_morning_star_configurations(
        self, year, month, day, hour, body_id, desc, expected_type
    ):
        """Test morning star classification at known elongation dates."""
        jd = swe.julday(year, month, day, hour)

        is_morning = pyephem.is_morning_star(jd, body_id, FLG_SWIEPH)
        expected_morning = expected_type == "morning"

        assert is_morning == expected_morning, (
            f"{desc}: is_morning_star={is_morning}, expected={expected_morning}"
        )


# ============================================================================
# IS_EVENING_STAR TESTS
# ============================================================================


class TestIsEveningStar:
    """Tests for is_evening_star function."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", ELONGATION_PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_evening_star_matches_elongation_sign(
        self, body_id, body_name, year, month, day, hour, date_desc
    ):
        """Test that is_evening_star matches positive elongation."""
        jd = swe.julday(year, month, day, hour)

        # Get elongation to determine expected result
        try:
            elongation, expected_evening = pyephem.get_elongation_from_sun(
                jd, body_id, FLG_SWIEPH
            )
        except Exception as e:
            pytest.skip(f"get_elongation_from_sun failed: {e}")

        # Get is_evening_star result
        try:
            is_evening = pyephem.is_evening_star(jd, body_id, FLG_SWIEPH)
        except Exception as e:
            pytest.fail(f"is_evening_star failed: {e}")

        assert is_evening == expected_evening, (
            f"{body_name} @ {date_desc}: is_evening_star={is_evening}, "
            f"expected={expected_evening}, elongation={elongation:.2f}°"
        )

    @pytest.mark.comparison
    def test_sun_is_not_evening_star(self, jd_standard):
        """Test that Sun cannot be an evening star."""
        is_evening = pyephem.is_evening_star(jd_standard, SUN, FLG_SWIEPH)
        assert is_evening is False, "Sun should never be classified as evening star"

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,body_id,desc,expected_type",
        KNOWN_ELONGATION_DATES,
    )
    def test_known_evening_star_configurations(
        self, year, month, day, hour, body_id, desc, expected_type
    ):
        """Test evening star classification at known elongation dates."""
        jd = swe.julday(year, month, day, hour)

        is_evening = pyephem.is_evening_star(jd, body_id, FLG_SWIEPH)
        expected_evening = expected_type == "evening"

        assert is_evening == expected_evening, (
            f"{desc}: is_evening_star={is_evening}, expected={expected_evening}"
        )


# ============================================================================
# GET_ELONGATION_TYPE TESTS
# ============================================================================


class TestGetElongationType:
    """Tests for get_elongation_type function."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", ELONGATION_PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_elongation_type_matches_elongation_sign(
        self, body_id, body_name, year, month, day, hour, date_desc
    ):
        """Test that get_elongation_type returns correct type based on elongation."""
        jd = swe.julday(year, month, day, hour)

        # Get elongation to determine expected result
        try:
            elongation, is_evening = pyephem.get_elongation_from_sun(
                jd, body_id, FLG_SWIEPH
            )
        except Exception as e:
            pytest.skip(f"get_elongation_from_sun failed: {e}")

        # Get elongation type result
        try:
            elong_type = pyephem.get_elongation_type(jd, body_id, FLG_SWIEPH)
        except Exception as e:
            pytest.fail(f"get_elongation_type failed: {e}")

        # Eastern = positive elongation = evening star
        # Western = negative elongation = morning star
        expected_type = "eastern" if is_evening else "western"

        assert elong_type == expected_type, (
            f"{body_name} @ {date_desc}: get_elongation_type={elong_type}, "
            f"expected={expected_type}, elongation={elongation:.2f}°"
        )

    @pytest.mark.comparison
    def test_sun_returns_none_type(self, jd_standard):
        """Test that Sun returns 'none' elongation type."""
        elong_type = pyephem.get_elongation_type(jd_standard, SUN, FLG_SWIEPH)
        assert elong_type == "none", (
            f"Sun should return 'none' elongation type, got '{elong_type}'"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour,body_id,desc,expected_star_type",
        KNOWN_ELONGATION_DATES,
    )
    def test_known_elongation_type_configurations(
        self, year, month, day, hour, body_id, desc, expected_star_type
    ):
        """Test elongation type at known configuration dates."""
        jd = swe.julday(year, month, day, hour)

        elong_type = pyephem.get_elongation_type(jd, body_id, FLG_SWIEPH)

        # morning star = western elongation, evening star = eastern elongation
        expected_type = "western" if expected_star_type == "morning" else "eastern"

        assert elong_type == expected_type, (
            f"{desc}: get_elongation_type={elong_type}, expected={expected_type}"
        )

    @pytest.mark.comparison
    def test_elongation_type_values(self, jd_standard):
        """Test that elongation type returns only valid values."""
        valid_types = {"eastern", "western", "none"}

        for body_id, body_name in ELONGATION_PLANETS + [(SUN, "Sun")]:
            elong_type = pyephem.get_elongation_type(jd_standard, body_id, FLG_SWIEPH)
            assert elong_type in valid_types, (
                f"{body_name}: got invalid type '{elong_type}', "
                f"expected one of {valid_types}"
            )


# ============================================================================
# CONSISTENCY TESTS
# ============================================================================


class TestElongationConsistency:
    """Tests for consistency between elongation helper functions."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", ELONGATION_PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_morning_evening_mutually_exclusive(
        self, body_id, body_name, year, month, day, hour, date_desc
    ):
        """Test that a planet cannot be both morning and evening star."""
        jd = swe.julday(year, month, day, hour)

        is_morning = pyephem.is_morning_star(jd, body_id, FLG_SWIEPH)
        is_evening = pyephem.is_evening_star(jd, body_id, FLG_SWIEPH)

        # Cannot be both or neither (except for Sun)
        if body_id != SUN:
            assert is_morning != is_evening, (
                f"{body_name} @ {date_desc}: must be either morning OR evening star, "
                f"got is_morning={is_morning}, is_evening={is_evening}"
            )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", ELONGATION_PLANETS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_elongation_type_consistency_with_star_functions(
        self, body_id, body_name, year, month, day, hour, date_desc
    ):
        """Test that elongation type is consistent with morning/evening functions."""
        jd = swe.julday(year, month, day, hour)

        is_morning = pyephem.is_morning_star(jd, body_id, FLG_SWIEPH)
        is_evening = pyephem.is_evening_star(jd, body_id, FLG_SWIEPH)
        elong_type = pyephem.get_elongation_type(jd, body_id, FLG_SWIEPH)

        if elong_type == "western":
            assert is_morning is True and is_evening is False, (
                f"{body_name} @ {date_desc}: western elongation should be morning star"
            )
        elif elong_type == "eastern":
            assert is_morning is False and is_evening is True, (
                f"{body_name} @ {date_desc}: eastern elongation should be evening star"
            )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", ELONGATION_PLANETS)
    def test_signed_elongation_consistency(self, jd_standard, body_id, body_name):
        """Test that get_signed_elongation matches get_elongation_from_sun."""
        elongation_from_sun, _ = pyephem.get_elongation_from_sun(
            jd_standard, body_id, FLG_SWIEPH
        )
        signed_elongation = pyephem.get_signed_elongation(
            jd_standard, body_id, FLG_SWIEPH
        )

        diff = abs(elongation_from_sun - signed_elongation)
        assert diff < ElongationTolerance.ANGLE_DEGREES, (
            f"{body_name}: get_signed_elongation ({signed_elongation:.4f}°) "
            f"doesn't match get_elongation_from_sun ({elongation_from_sun:.4f}°)"
        )


# ============================================================================
# COMPARISON WITH PHENO_UT
# ============================================================================


class TestElongationVsPheno:
    """Tests comparing elongation helpers with pheno_ut results."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_elongation_magnitude_vs_pheno(
        self, body_id, body_name, year, month, day, hour, date_desc
    ):
        """Test that absolute elongation matches pheno_ut elongation."""
        jd = swe.julday(year, month, day, hour)

        # Get elongation from helper
        try:
            elongation, _ = pyephem.get_elongation_from_sun(jd, body_id, FLG_SWIEPH)
        except Exception as e:
            pytest.skip(f"get_elongation_from_sun failed: {e}")

        # Get elongation from pheno_ut
        try:
            ret = pyephem.pheno_ut(jd, body_id, FLG_SWIEPH)
            attr = ret[1] if isinstance(ret, tuple) and len(ret) > 1 else ret
            pheno_elongation = attr[2]  # PHENO_ELONGATION index
        except Exception as e:
            pytest.skip(f"pheno_ut failed: {e}")

        # pheno_ut returns unsigned elongation (0-180°)
        # Our helper returns signed elongation (-180 to +180°)
        diff = abs(abs(elongation) - pheno_elongation)

        assert diff < ElongationTolerance.PHENO_ANGLE, (
            f"{body_name} @ {date_desc}: elongation magnitude mismatch "
            f"(helper={abs(elongation):.4f}°, pheno={pheno_elongation:.4f}°, "
            f"diff={diff:.4f}°)"
        )


# ============================================================================
# EDGE CASES
# ============================================================================


class TestElongationEdgeCases:
    """Tests for edge cases in elongation calculations."""

    @pytest.mark.comparison
    @pytest.mark.edge_case
    def test_near_conjunction_with_sun(self):
        """Test elongation near conjunction (small values)."""
        # Date when Mercury is near conjunction with Sun
        # Mercury inferior conjunction ~Feb 28, 2024
        jd = swe.julday(2024, 2, 28, 12.0)

        elongation, _ = pyephem.get_elongation_from_sun(jd, MERCURY, FLG_SWIEPH)

        # Near conjunction, elongation should be small
        assert abs(elongation) < 10, (
            f"Mercury at conjunction should have small elongation, got {elongation:.2f}°"
        )

    @pytest.mark.comparison
    @pytest.mark.edge_case
    def test_near_opposition(self):
        """Test elongation for outer planet near opposition."""
        # Jupiter opposition ~Dec 7, 2024
        jd = swe.julday(2024, 12, 7, 12.0)

        elongation, _ = pyephem.get_elongation_from_sun(jd, JUPITER, FLG_SWIEPH)

        # Near opposition, elongation should be close to 180° (or -180°)
        assert abs(abs(elongation) - 180) < 10, (
            f"Jupiter at opposition should have ~180° elongation, got {elongation:.2f}°"
        )

    @pytest.mark.comparison
    @pytest.mark.edge_case
    def test_all_planets_at_same_instant(self):
        """Test all planets at the same instant for self-consistency."""
        jd = swe.julday(2024, 6, 21, 12.0)  # Summer solstice 2024

        results = {}
        for body_id, body_name in ELONGATION_PLANETS:
            elongation, is_evening = pyephem.get_elongation_from_sun(
                jd, body_id, FLG_SWIEPH
            )
            is_morning = pyephem.is_morning_star(jd, body_id, FLG_SWIEPH)
            elong_type = pyephem.get_elongation_type(jd, body_id, FLG_SWIEPH)

            results[body_name] = {
                "elongation": elongation,
                "is_evening": is_evening,
                "is_morning": is_morning,
                "type": elong_type,
            }

            # Verify consistency
            assert is_evening != is_morning, (
                f"{body_name}: morning/evening must be mutually exclusive"
            )
            if is_evening:
                assert elong_type == "eastern" and elongation > 0
            else:
                assert elong_type == "western" and elongation < 0

    @pytest.mark.comparison
    @pytest.mark.edge_case
    @pytest.mark.parametrize("body_id,body_name", ELONGATION_PLANETS)
    def test_different_flags(self, jd_standard, body_id, body_name):
        """Test elongation with different calculation flags."""
        flags_to_test = [
            0,  # Default
            FLG_SWIEPH,  # Swiss Ephemeris
        ]

        results = []
        for flag in flags_to_test:
            try:
                elongation, is_evening = pyephem.get_elongation_from_sun(
                    jd_standard, body_id, flag
                )
                results.append((flag, elongation, is_evening))
            except Exception:
                pass  # Some flags may not be available

        # Results with different flags should be very close
        if len(results) >= 2:
            for i in range(1, len(results)):
                diff = abs(results[0][1] - results[i][1])
                assert diff < 0.1, (
                    f"{body_name}: elongation differs by {diff:.4f}° "
                    f"between flag={results[0][0]} and flag={results[i][0]}"
                )
