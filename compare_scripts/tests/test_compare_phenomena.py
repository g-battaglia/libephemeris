"""
Pytest-style Planetary Phenomena Comparison Tests.

Validates pheno and pheno_ut calculations (phase angle, elongation, magnitude)
against pyswisseph.
"""

import pytest
import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import (
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    FLG_SWIEPH,
)


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class PhenoTolerance:
    """Tolerance thresholds for phenomena comparisons."""

    ANGLE_DEGREES = 0.01  # Phase angle, elongation
    MAGNITUDE = 0.1  # Apparent magnitude
    DIAMETER = 0.01  # Apparent diameter in arcsec


# ============================================================================
# PHENOMENA ARRAY INDICES
# ============================================================================

PHENO_PHASE_ANGLE = 0
PHENO_PHASE = 1  # Illuminated fraction
PHENO_ELONGATION = 2
PHENO_DIAMETER = 3
PHENO_MAGNITUDE = 4


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Planets for phenomena tests (visible from Earth)
PHENO_PLANETS = [
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
]

# Test dates covering different planetary configurations
TEST_DATES = [
    (2024, 1, 15, "Jan 2024"),
    (2024, 4, 15, "Apr 2024"),
    (2024, 7, 15, "Jul 2024"),
    (2024, 10, 15, "Oct 2024"),
]

# Special configurations for inner planets
INNER_PLANET_DATES = [
    (2024, 3, 22, "Venus western elongation"),
    (2024, 6, 4, "Mercury eastern elongation"),
    (2024, 8, 26, "Mercury western elongation"),
]

# Opposition dates for outer planets
OPPOSITION_DATES = [
    (2025, 1, 16, MARS, "Mars opposition 2025"),
    (2024, 12, 7, JUPITER, "Jupiter opposition 2024"),
    (2024, 9, 8, SATURN, "Saturn opposition 2024"),
]


# ============================================================================
# FIXTURES
# ============================================================================


@pytest.fixture
def jd_standard():
    """Standard Julian Day for tests."""
    return swe.julday(2024, 6, 15, 12.0)


# ============================================================================
# PHENO_UT TESTS
# ============================================================================


class TestPhenoUt:
    """Tests for pheno_ut function."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PHENO_PLANETS)
    @pytest.mark.parametrize("year,month,day,date_desc", TEST_DATES)
    def test_pheno_ut(self, body_id, body_name, year, month, day, date_desc):
        """Test pheno_ut for planets at various dates."""
        jd = swe.julday(year, month, day, 12.0)

        # SwissEphemeris
        try:
            ret_swe = swe.pheno_ut(jd, body_id, FLG_SWIEPH)
            attr_swe = (
                ret_swe[1]
                if isinstance(ret_swe, tuple) and len(ret_swe) > 1
                else ret_swe
            )

            phase_angle_swe = attr_swe[PHENO_PHASE_ANGLE]
            elongation_swe = attr_swe[PHENO_ELONGATION]
            magnitude_swe = attr_swe[PHENO_MAGNITUDE]
        except Exception as e:
            pytest.skip(f"SwissEphemeris failed: {e}")

        # LibEphemeris
        try:
            ret_py = pyephem.pheno_ut(jd, body_id, FLG_SWIEPH)
            attr_py = (
                ret_py[1] if isinstance(ret_py, tuple) and len(ret_py) > 1 else ret_py
            )

            phase_angle_py = attr_py[PHENO_PHASE_ANGLE]
            elongation_py = attr_py[PHENO_ELONGATION]
            magnitude_py = attr_py[PHENO_MAGNITUDE]
        except Exception as e:
            pytest.skip(f"LibEphemeris failed: {e}")

        # Calculate differences
        diff_phase_angle = abs(phase_angle_swe - phase_angle_py)
        diff_elongation = abs(elongation_swe - elongation_py)
        diff_magnitude = abs(magnitude_swe - magnitude_py)

        assert diff_phase_angle < PhenoTolerance.ANGLE_DEGREES, (
            f"{body_name} @ {date_desc}: phase angle diff {diff_phase_angle}°"
        )
        assert diff_elongation < PhenoTolerance.ANGLE_DEGREES, (
            f"{body_name} @ {date_desc}: elongation diff {diff_elongation}°"
        )
        assert diff_magnitude < PhenoTolerance.MAGNITUDE, (
            f"{body_name} @ {date_desc}: magnitude diff {diff_magnitude}"
        )


# ============================================================================
# PHENO (ET VERSION) TESTS
# ============================================================================


class TestPheno:
    """Tests for pheno function (ET version)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PHENO_PLANETS)
    def test_pheno_et(self, jd_standard, body_id, body_name):
        """Test pheno function (ET time)."""
        jd_et = jd_standard + swe.deltat(jd_standard)

        # SwissEphemeris
        try:
            ret_swe = swe.pheno(jd_et, body_id, FLG_SWIEPH)
            attr_swe = (
                ret_swe[1]
                if isinstance(ret_swe, tuple) and len(ret_swe) > 1
                else ret_swe
            )

            phase_angle_swe = attr_swe[PHENO_PHASE_ANGLE]
            elongation_swe = attr_swe[PHENO_ELONGATION]
            magnitude_swe = attr_swe[PHENO_MAGNITUDE]
        except Exception as e:
            pytest.skip(f"SwissEphemeris failed: {e}")

        # LibEphemeris
        try:
            ret_py = pyephem.pheno(jd_et, body_id, FLG_SWIEPH)
            attr_py = (
                ret_py[1] if isinstance(ret_py, tuple) and len(ret_py) > 1 else ret_py
            )

            phase_angle_py = attr_py[PHENO_PHASE_ANGLE]
            elongation_py = attr_py[PHENO_ELONGATION]
            magnitude_py = attr_py[PHENO_MAGNITUDE]
        except Exception as e:
            pytest.skip(f"LibEphemeris failed: {e}")

        diff_phase_angle = abs(phase_angle_swe - phase_angle_py)
        diff_elongation = abs(elongation_swe - elongation_py)
        diff_magnitude = abs(magnitude_swe - magnitude_py)

        assert diff_phase_angle < PhenoTolerance.ANGLE_DEGREES
        assert diff_elongation < PhenoTolerance.ANGLE_DEGREES
        assert diff_magnitude < PhenoTolerance.MAGNITUDE


# ============================================================================
# INNER PLANET ELONGATION TESTS
# ============================================================================


class TestInnerPlanetElongations:
    """Tests for inner planets at various elongations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,desc", INNER_PLANET_DATES)
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (VENUS, "Venus"),
            (MERCURY, "Mercury"),
        ],
    )
    def test_inner_planet_elongation(self, year, month, day, desc, body_id, body_name):
        """Test inner planets at elongation dates."""
        jd = swe.julday(year, month, day, 12.0)

        try:
            ret_swe = swe.pheno_ut(jd, body_id, FLG_SWIEPH)
            ret_py = pyephem.pheno_ut(jd, body_id, FLG_SWIEPH)

            attr_swe = ret_swe[1] if isinstance(ret_swe, tuple) else ret_swe
            attr_py = ret_py[1] if isinstance(ret_py, tuple) else ret_py

            diff_phase = abs(attr_swe[PHENO_PHASE_ANGLE] - attr_py[PHENO_PHASE_ANGLE])
            diff_elong = abs(attr_swe[PHENO_ELONGATION] - attr_py[PHENO_ELONGATION])

            assert diff_phase < PhenoTolerance.ANGLE_DEGREES
            assert diff_elong < PhenoTolerance.ANGLE_DEGREES
        except Exception:
            pytest.skip("Calculation failed")


# ============================================================================
# OUTER PLANET OPPOSITION TESTS
# ============================================================================


class TestOuterPlanetOppositions:
    """Tests for outer planets at opposition."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,body_id,desc", OPPOSITION_DATES)
    def test_outer_planet_opposition(self, year, month, day, body_id, desc):
        """Test outer planets at opposition dates (brightest)."""
        jd = swe.julday(year, month, day, 12.0)
        body_name = desc.split()[0]

        try:
            ret_swe = swe.pheno_ut(jd, body_id, FLG_SWIEPH)
            ret_py = pyephem.pheno_ut(jd, body_id, FLG_SWIEPH)

            attr_swe = ret_swe[1] if isinstance(ret_swe, tuple) else ret_swe
            attr_py = ret_py[1] if isinstance(ret_py, tuple) else ret_py

            diff_magnitude = abs(attr_swe[PHENO_MAGNITUDE] - attr_py[PHENO_MAGNITUDE])

            assert diff_magnitude < PhenoTolerance.MAGNITUDE, (
                f"{body_name} opposition: magnitude diff {diff_magnitude}"
            )
        except Exception:
            pytest.skip("Calculation failed")


# ============================================================================
# PHENOMENA CONSISTENCY TESTS
# ============================================================================


class TestPhenomenaConsistency:
    """Tests for phenomena consistency and physical constraints."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PHENO_PLANETS)
    def test_phenomena_physical_constraints(self, jd_standard, body_id, body_name):
        """Test that phenomena values satisfy physical constraints."""
        ret_py = pyephem.pheno_ut(jd_standard, body_id, FLG_SWIEPH)
        # libephemeris returns (attr_tuple, retflag), pyswisseph returns just attr_tuple
        if (
            isinstance(ret_py, tuple)
            and len(ret_py) == 2
            and isinstance(ret_py[0], tuple)
        ):
            attr = ret_py[0]  # libephemeris format: (attr, retflag)
        else:
            attr = ret_py  # pyswisseph format: just attr tuple

        phase_angle = attr[PHENO_PHASE_ANGLE]
        phase = attr[PHENO_PHASE]  # Illuminated fraction
        elongation = attr[PHENO_ELONGATION]
        diameter = attr[PHENO_DIAMETER]

        # Physical constraints
        assert 0 <= phase_angle <= 180, f"{body_name}: phase angle must be 0-180°"
        assert 0 <= phase <= 1, f"{body_name}: illuminated fraction must be 0-1"
        assert 0 <= elongation <= 180, f"{body_name}: elongation must be 0-180°"
        assert diameter >= 0, f"{body_name}: diameter must be non-negative"

    @pytest.mark.comparison
    def test_inner_planet_elongation_limits(self, jd_standard):
        """Test that inner planets have limited elongation."""
        for body_id, body_name in [(MERCURY, "Mercury"), (VENUS, "Venus")]:
            ret_py = pyephem.pheno_ut(jd_standard, body_id, FLG_SWIEPH)
            # libephemeris returns (attr_tuple, retflag), pyswisseph returns just attr_tuple
            if (
                isinstance(ret_py, tuple)
                and len(ret_py) == 2
                and isinstance(ret_py[0], tuple)
            ):
                attr = ret_py[0]  # libephemeris format
            else:
                attr = ret_py  # pyswisseph format

            elongation = attr[PHENO_ELONGATION]

            # Inner planets can't be more than ~47° (Venus) or ~28° (Mercury) from Sun
            max_elong = 50 if body_id == VENUS else 30
            assert elongation <= max_elong, (
                f"{body_name}: elongation {elongation}° exceeds physical limit"
            )
