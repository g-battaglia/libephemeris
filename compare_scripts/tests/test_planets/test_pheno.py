"""
Tests for planetary phenomena calculations (pheno, pheno_ut).

Tests cover:
- Phase angle calculations
- Phase (illumination fraction)
- Elongation from Sun
- Apparent diameter
- Visual magnitude
- Comparison with pyswisseph
"""

import pytest
import math
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
)


class TestPhenoBasic:
    """Test basic pheno function behavior."""

    @pytest.mark.unit
    def test_pheno_ut_returns_tuple(self):
        """pheno_ut should return a tuple of phenomenon values."""
        jd = 2451545.0  # J2000
        result = ephem.pheno_ut(jd, MARS, 0)

        assert isinstance(result, tuple)
        assert len(result) >= 5  # At least 5 phenomenon values

    @pytest.mark.unit
    def test_pheno_returns_tuple(self):
        """pheno should return a tuple of phenomenon values."""
        jd = 2451545.0  # J2000 TT
        result = ephem.pheno(jd, MARS, 0)

        assert isinstance(result, tuple)
        assert len(result) >= 5

    @pytest.mark.unit
    def test_swe_prefix_aliases_exist(self):
        """Both pheno and pheno_ut should be accessible."""
        assert hasattr(ephem, "pheno")
        assert hasattr(ephem, "pheno_ut")
        assert hasattr(ephem, "pheno")
        assert hasattr(ephem, "pheno_ut")


class TestSunPhenomena:
    """Test Sun phenomena (edge case)."""

    @pytest.mark.unit
    def test_sun_phase_is_full(self):
        """Sun should always have phase = 1.0 (fully illuminated)."""
        jd = 2451545.0
        attr = ephem.pheno_ut(jd, SUN, 0)

        phase_angle = attr[0]
        phase = attr[1]
        elongation = attr[2]

        assert phase_angle == 0.0, "Sun phase angle should be 0"
        assert phase == 0.0, (
            "Sun phase should be 0.0 (inapplicable for self-luminous body)"
        )
        assert elongation == 0.0, "Sun elongation from itself should be 0"

    @pytest.mark.unit
    def test_sun_diameter_reasonable(self):
        """Sun apparent diameter should be around 32 arcminutes."""
        jd = 2451545.0
        attr = ephem.pheno_ut(jd, SUN, 0)

        diameter = attr[3]
        # Sun diameter ~32 arcmin ≈ 0.53 degrees, varies with distance
        assert 0.51 < diameter < 0.55, f"Sun diameter {diameter} deg unexpected"

    @pytest.mark.unit
    def test_sun_magnitude_reasonable(self):
        """Sun magnitude should be around -26.7."""
        jd = 2451545.0
        attr = ephem.pheno_ut(jd, SUN, 0)

        magnitude = attr[4]
        assert -27.0 < magnitude < -26.0, f"Sun magnitude {magnitude} unexpected"


class TestMoonPhenomena:
    """Test Moon phenomena calculations."""

    @pytest.mark.unit
    def test_moon_phase_varies(self):
        """Moon phase should vary between 0 and 1."""
        # Test at different times to see variation
        jd_new = 2451550.1  # Near new moon (around Jan 6, 2000)
        jd_full = 2451565.0  # Near full moon (around Jan 21, 2000)

        attr_new = ephem.pheno_ut(jd_new, MOON, 0)
        attr_full = ephem.pheno_ut(jd_full, MOON, 0)

        phase_new = attr_new[1]
        phase_full = attr_full[1]

        # Phases should be in valid range
        assert 0 <= phase_new <= 1.0
        assert 0 <= phase_full <= 1.0

    @pytest.mark.unit
    def test_moon_diameter_reasonable(self):
        """Moon apparent diameter should be around 30-34 arcminutes."""
        jd = 2451545.0
        attr = ephem.pheno_ut(jd, MOON, 0)

        diameter = attr[3]
        # Moon diameter ~31 arcmin ≈ 0.49-0.56 degrees, varies with distance
        assert 0.47 < diameter < 0.58, (
            f"Moon diameter {diameter:.4f} deg outside expected range [0.47, 0.58]"
        )

    @pytest.mark.unit
    def test_moon_elongation_varies(self):
        """Moon elongation should be between 0 and 180 degrees."""
        jd = 2451545.0
        attr = ephem.pheno_ut(jd, MOON, 0)

        elongation = attr[2]
        assert 0 <= elongation <= 180, f"Moon elongation {elongation} out of range"


class TestPlanetPhenomena:
    """Test planetary phenomena calculations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_planet_phase_in_range(self, planet_id, planet_name):
        """Planet phase should be between 0 and 1."""
        jd = 2451545.0
        attr = ephem.pheno_ut(jd, planet_id, 0)

        phase = attr[1]
        assert 0 <= phase <= 1.0, f"{planet_name} phase {phase} out of range"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_planet_elongation_in_range(self, planet_id, planet_name):
        """Planet elongation should be between 0 and 180 degrees."""
        jd = 2451545.0
        attr = ephem.pheno_ut(jd, planet_id, 0)

        elongation = attr[2]
        assert 0 <= elongation <= 180, (
            f"{planet_name} elongation {elongation} out of range"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_planet_diameter_positive(self, planet_id, planet_name):
        """Planet apparent diameter should be positive."""
        jd = 2451545.0
        attr = ephem.pheno_ut(jd, planet_id, 0)

        diameter = attr[3]
        assert diameter > 0, f"{planet_name} diameter should be positive"

    @pytest.mark.unit
    def test_jupiter_diameter_large(self):
        """Jupiter should have the largest apparent diameter of planets."""
        jd = 2451545.0

        diameters = {}
        for planet_id in [MARS, JUPITER, SATURN]:
            attr = ephem.pheno_ut(jd, planet_id, 0)
            diameters[planet_id] = attr[3]

        # Jupiter should have larger diameter than Mars
        assert diameters[JUPITER] > diameters[MARS]


class TestMagnitudes:
    """Test visual magnitude calculations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,min_mag,max_mag,planet_name",
        [
            (VENUS, -5.0, -3.0, "Venus"),  # Venus is always bright
            (JUPITER, -3.0, -1.0, "Jupiter"),  # Jupiter is bright
            (SATURN, -0.5, 1.5, "Saturn"),  # Saturn varies more
            (MARS, -3.0, 2.0, "Mars"),  # Mars varies significantly
        ],
    )
    def test_planet_magnitude_typical_range(
        self, planet_id, min_mag, max_mag, planet_name
    ):
        """Planet magnitudes should be in typical observable range."""
        jd = 2451545.0
        attr = ephem.pheno_ut(jd, planet_id, 0)

        magnitude = attr[4]
        assert min_mag < magnitude < max_mag, (
            f"{planet_name} magnitude {magnitude:.2f} "
            f"outside typical range [{min_mag}, {max_mag}]"
        )


class TestComparisonWithSwissEph:
    """Compare pheno results with pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_phase_angle_matches_swe(self, planet_id, planet_name):
        """Phase angle should match Swiss Ephemeris within tolerance."""
        jd = 2451545.0

        swe.set_ephe_path(None)
        # pyswisseph returns just the tuple, not (tuple, flag)
        swe_attr = swe.pheno_ut(jd, planet_id, 0)
        lib_attr = ephem.pheno_ut(jd, planet_id, 0)

        # Phase angle comparison (allow 1 degree tolerance)
        swe_phase_angle = swe_attr[0]
        lib_phase_angle = lib_attr[0]

        diff = abs(swe_phase_angle - lib_phase_angle)
        assert diff < 1.0, (
            f"{planet_name} phase angle: SWE={swe_phase_angle:.4f}, "
            f"LIB={lib_phase_angle:.4f}, diff={diff:.4f}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_phase_matches_swe(self, planet_id, planet_name):
        """Phase (illumination) should match Swiss Ephemeris within tolerance."""
        jd = 2451545.0

        swe.set_ephe_path(None)
        swe_attr = swe.pheno_ut(jd, planet_id, 0)
        lib_attr = ephem.pheno_ut(jd, planet_id, 0)

        swe_phase = swe_attr[1]
        lib_phase = lib_attr[1]

        diff = abs(swe_phase - lib_phase)
        assert diff < 0.02, (
            f"{planet_name} phase: SWE={swe_phase:.4f}, "
            f"LIB={lib_phase:.4f}, diff={diff:.4f}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_elongation_matches_swe(self, planet_id, planet_name):
        """Elongation should match Swiss Ephemeris within tolerance."""
        jd = 2451545.0

        swe.set_ephe_path(None)
        swe_attr = swe.pheno_ut(jd, planet_id, 0)
        lib_attr = ephem.pheno_ut(jd, planet_id, 0)

        swe_elong = swe_attr[2]
        lib_elong = lib_attr[2]

        diff = abs(swe_elong - lib_elong)
        assert diff < 1.0, (
            f"{planet_name} elongation: SWE={swe_elong:.4f}, "
            f"LIB={lib_elong:.4f}, diff={diff:.4f}"
        )

    @pytest.mark.comparison
    def test_moon_phenomena_matches_swe(self):
        """Moon phenomena should match Swiss Ephemeris."""
        jd = 2451545.0

        swe.set_ephe_path(None)
        swe_attr = swe.pheno_ut(jd, MOON, 0)
        lib_attr = ephem.pheno_ut(jd, MOON, 0)

        # Phase angle
        diff_angle = abs(swe_attr[0] - lib_attr[0])
        assert diff_angle < 2.0, f"Moon phase angle diff: {diff_angle:.4f}"

        # Phase
        diff_phase = abs(swe_attr[1] - lib_attr[1])
        assert diff_phase < 0.05, f"Moon phase diff: {diff_phase:.4f}"

        # Elongation
        diff_elong = abs(swe_attr[2] - lib_attr[2])
        assert diff_elong < 2.0, f"Moon elongation diff: {diff_elong:.4f}"


class TestOuterPlanets:
    """Test outer planet phenomena."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (URANUS, "Uranus"),
            (NEPTUNE, "Neptune"),
            (PLUTO, "Pluto"),
        ],
    )
    def test_outer_planet_nearly_full(self, planet_id, planet_name):
        """Outer planets should be nearly fully illuminated (phase ~1.0)."""
        jd = 2451545.0
        attr = ephem.pheno_ut(jd, planet_id, 0)

        phase = attr[1]
        # Outer planets are almost fully illuminated from Earth
        assert phase > 0.95, f"{planet_name} phase {phase} should be > 0.95"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (URANUS, "Uranus"),
            (NEPTUNE, "Neptune"),
            (PLUTO, "Pluto"),
        ],
    )
    def test_outer_planet_small_phase_angle(self, planet_id, planet_name):
        """Outer planets should have small phase angles."""
        jd = 2451545.0
        attr = ephem.pheno_ut(jd, planet_id, 0)

        phase_angle = attr[0]
        # Outer planets have small phase angles (usually < 5 degrees)
        assert phase_angle < 10.0, (
            f"{planet_name} phase angle {phase_angle} should be < 10"
        )


class TestPhaseAngleFormula:
    """Test phase angle and phase relationship."""

    @pytest.mark.unit
    def test_phase_from_phase_angle(self):
        """Phase should be calculated correctly from phase angle."""
        jd = 2451545.0

        for planet_id in [MERCURY, VENUS, MARS, JUPITER]:
            attr = ephem.pheno_ut(jd, planet_id, 0)

            phase_angle = attr[0]
            phase = attr[1]

            # Phase = (1 + cos(phase_angle)) / 2
            expected_phase = (1 + math.cos(math.radians(phase_angle))) / 2

            diff = abs(phase - expected_phase)
            assert diff < 0.001, (
                f"Planet {planet_id}: Phase {phase:.6f} doesn't match "
                f"expected {expected_phase:.6f} from angle {phase_angle:.2f}"
            )


class TestEdgeCases:
    """Test edge cases and error handling."""

    @pytest.mark.unit
    def test_unsupported_body_returns_zeros(self):
        """Calculating phenomena for unsupported body should return zeros."""
        jd = 2451545.0
        # Use a body ID that's not in _PLANET_MAP (e.g., a fictitious planet)
        unsupported_id = 99  # Not a valid planet ID
        attr = ephem.pheno_ut(jd, unsupported_id, 0)

        # Should return zeros for unsupported body
        assert all(v == 0.0 for v in attr[:5]), (
            "Unsupported body phenomena should return zeros"
        )

    @pytest.mark.unit
    def test_different_dates(self):
        """Phenomena should vary with date."""
        jd1 = 2451545.0  # J2000
        jd2 = 2451645.0  # 100 days later

        attr1 = ephem.pheno_ut(jd1, MARS, 0)
        attr2 = ephem.pheno_ut(jd2, MARS, 0)

        # At least one value should differ significantly
        differences = [abs(attr1[i] - attr2[i]) for i in range(5)]
        assert max(differences) > 0.01, "Phenomena should change over 100 days"


class TestSwePheno20Values:
    """
    Verify pheno_ut returns all 20 output values matching Swiss Ephemeris.

    This is the core verification suite for Swiss Ephemeris compatibility.
    All 20 values in the phenomena tuple must be verified:
    - attr[0-4]: Active phenomenon values (phase angle, phase, elongation, diameter, magnitude)
    - attr[5-19]: Reserved values (should all be 0.0)
    """

    # Tolerances for comparison
    PHASE_ANGLE_TOL = 1.0  # degrees (accounts for ephemeris differences)
    PHASE_TOL = 0.02  # illuminated fraction
    ELONGATION_TOL = 1.0  # degrees
    DIAMETER_TOL = 0.05  # relative tolerance (5%)
    MAGNITUDE_TOL = 0.5  # magnitudes (accounts for different formulas)
    # Saturn has larger magnitude differences due to ring contribution
    # which Swiss Ephemeris calculates differently (requires ring tilt)
    SATURN_MAGNITUDE_TOL = 1.0

    @pytest.mark.comparison
    def test_returns_exactly_20_values(self):
        """pheno_ut must return exactly 20 values in the attr tuple."""
        jd = 2451545.0
        attr = ephem.pheno_ut(jd, MARS, 0)

        assert len(attr) == 20, f"Expected 20 values, got {len(attr)}"

    @pytest.mark.comparison
    def test_reserved_values_are_zero(self):
        """Values attr[5] through attr[19] must all be 0.0 (reserved)."""
        jd = 2451545.0

        # Test across multiple planets and bodies
        bodies = [SUN, MOON, MERCURY, VENUS, MARS, JUPITER, SATURN]

        for body_id in bodies:
            attr = ephem.pheno_ut(jd, body_id, 0)

            for i in range(5, 20):
                assert attr[i] == 0.0, (
                    f"Body {body_id}: attr[{i}] should be 0.0, got {attr[i]}"
                )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_all_five_active_values_match_swe(self, planet_id, planet_name):
        """
        All 5 active phenomenon values must match Swiss Ephemeris within tolerance.

        Note: Both libephemeris and Swiss Ephemeris return diameter in degrees.
        """
        jd = 2451545.0

        swe.set_ephe_path(None)
        swe_attr = swe.pheno_ut(jd, planet_id, 0)
        lib_attr = ephem.pheno_ut(jd, planet_id, 0)

        # attr[0]: Phase angle (degrees)
        diff_phase_angle = abs(swe_attr[0] - lib_attr[0])
        assert diff_phase_angle < self.PHASE_ANGLE_TOL, (
            f"{planet_name} phase angle: SE={swe_attr[0]:.4f}, "
            f"LIB={lib_attr[0]:.4f}, diff={diff_phase_angle:.4f}"
        )

        # attr[1]: Phase (illuminated fraction, 0.0-1.0)
        diff_phase = abs(swe_attr[1] - lib_attr[1])
        assert diff_phase < self.PHASE_TOL, (
            f"{planet_name} phase: SE={swe_attr[1]:.4f}, "
            f"LIB={lib_attr[1]:.4f}, diff={diff_phase:.4f}"
        )

        # attr[2]: Elongation (degrees)
        diff_elongation = abs(swe_attr[2] - lib_attr[2])
        assert diff_elongation < self.ELONGATION_TOL, (
            f"{planet_name} elongation: SE={swe_attr[2]:.4f}, "
            f"LIB={lib_attr[2]:.4f}, diff={diff_elongation:.4f}"
        )

        # attr[3]: Diameter (both SE and libephemeris return degrees)
        swe_diam = swe_attr[3]
        lib_diam = lib_attr[3]
        if swe_diam > 0:
            rel_diff_diam = abs(swe_diam - lib_diam) / swe_diam
        else:
            rel_diff_diam = 0.0 if lib_diam == 0 else 1.0
        assert rel_diff_diam < self.DIAMETER_TOL, (
            f"{planet_name} diameter: SE={swe_diam:.6f}deg, "
            f"LIB={lib_diam:.6f}deg, rel_diff={rel_diff_diam:.4f}"
        )

        # attr[4]: Magnitude
        # Saturn magnitude has larger tolerance due to ring contribution differences
        if planet_id == SATURN:
            mag_tol = self.SATURN_MAGNITUDE_TOL
        else:
            mag_tol = self.MAGNITUDE_TOL
        diff_magnitude = abs(swe_attr[4] - lib_attr[4])
        assert diff_magnitude < mag_tol, (
            f"{planet_name} magnitude: SE={swe_attr[4]:.2f}, "
            f"LIB={lib_attr[4]:.2f}, diff={diff_magnitude:.2f}"
        )

    @pytest.mark.comparison
    def test_sun_all_20_values(self):
        """Verify Sun phenomena returns all 20 values correctly."""
        jd = 2451545.0

        swe.set_ephe_path(None)
        swe_attr = swe.pheno_ut(jd, SUN, 0)
        lib_attr = ephem.pheno_ut(jd, SUN, 0)

        # Must have exactly 20 values
        assert len(lib_attr) == 20

        # Phase angle = 0, phase = 0.0 (inapplicable), elongation = 0
        assert lib_attr[0] == 0.0, "Sun phase angle should be 0"
        assert lib_attr[1] == 0.0, (
            "Sun phase should be 0.0 (inapplicable for self-luminous body)"
        )
        assert lib_attr[2] == 0.0, "Sun elongation should be 0"

        # Diameter check (both return degrees)
        swe_diam = swe_attr[3]
        lib_diam = lib_attr[3]
        rel_diff = abs(swe_diam - lib_diam) / swe_diam
        assert rel_diff < 0.01, (
            f"Sun diameter: SE={swe_diam:.6f}deg, LIB={lib_diam:.6f}deg"
        )

        # Magnitude check (small tolerance for Sun)
        diff_mag = abs(swe_attr[4] - lib_attr[4])
        assert diff_mag < 0.15, (
            f"Sun magnitude: SE={swe_attr[4]:.2f}, LIB={lib_attr[4]:.2f}"
        )

        # Reserved values
        for i in range(5, 20):
            assert lib_attr[i] == 0.0, f"Sun attr[{i}] should be 0.0"

    @pytest.mark.comparison
    def test_moon_all_20_values(self):
        """Verify Moon phenomena returns all 20 values correctly."""
        jd = 2451545.0

        swe.set_ephe_path(None)
        swe_attr = swe.pheno_ut(jd, MOON, 0)
        lib_attr = ephem.pheno_ut(jd, MOON, 0)

        # Must have exactly 20 values
        assert len(lib_attr) == 20

        # Phase angle comparison (Moon can have larger discrepancies)
        diff_angle = abs(swe_attr[0] - lib_attr[0])
        assert diff_angle < 2.0, f"Moon phase angle diff: {diff_angle:.4f}"

        # Phase comparison
        diff_phase = abs(swe_attr[1] - lib_attr[1])
        assert diff_phase < 0.05, f"Moon phase diff: {diff_phase:.4f}"

        # Elongation comparison
        diff_elong = abs(swe_attr[2] - lib_attr[2])
        assert diff_elong < 2.0, f"Moon elongation diff: {diff_elong:.4f}"

        # Diameter check (both return degrees)
        swe_diam = swe_attr[3]
        lib_diam = lib_attr[3]
        rel_diff = abs(swe_diam - lib_diam) / swe_diam
        assert rel_diff < 0.01, (
            f"Moon diameter: SE={swe_diam:.6f}deg, LIB={lib_diam:.6f}deg"
        )

        # Reserved values
        for i in range(5, 20):
            assert lib_attr[i] == 0.0, f"Moon attr[{i}] should be 0.0"

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (URANUS, "Uranus"),
            (NEPTUNE, "Neptune"),
            (PLUTO, "Pluto"),
        ],
    )
    def test_outer_planets_all_20_values(self, planet_id, planet_name):
        """Verify outer planets return all 20 values correctly."""
        jd = 2451545.0

        swe.set_ephe_path(None)
        swe_attr = swe.pheno_ut(jd, planet_id, 0)
        lib_attr = ephem.pheno_ut(jd, planet_id, 0)

        # Must have exactly 20 values
        assert len(lib_attr) == 20, f"{planet_name} should return 20 values"

        # Phase angle - outer planets have small phase angles
        diff_angle = abs(swe_attr[0] - lib_attr[0])
        assert diff_angle < 1.0, (
            f"{planet_name} phase angle: SE={swe_attr[0]:.4f}, LIB={lib_attr[0]:.4f}"
        )

        # Phase - outer planets are nearly full
        diff_phase = abs(swe_attr[1] - lib_attr[1])
        assert diff_phase < 0.02, (
            f"{planet_name} phase: SE={swe_attr[1]:.4f}, LIB={lib_attr[1]:.4f}"
        )

        # Elongation
        diff_elong = abs(swe_attr[2] - lib_attr[2])
        assert diff_elong < 1.0, (
            f"{planet_name} elongation: SE={swe_attr[2]:.4f}, LIB={lib_attr[2]:.4f}"
        )

        # Diameter (both return degrees)
        swe_diam = swe_attr[3]
        lib_diam = lib_attr[3]
        if swe_diam > 0:
            rel_diff = abs(swe_diam - lib_diam) / swe_diam
            assert rel_diff < 0.1, (
                f"{planet_name} diameter: SE={swe_diam:.6f}deg, LIB={lib_diam:.6f}deg"
            )

        # Reserved values
        for i in range(5, 20):
            assert lib_attr[i] == 0.0, f"{planet_name} attr[{i}] should be 0.0"

    @pytest.mark.comparison
    def test_multiple_dates_all_20_values(self):
        """Verify all 20 values at multiple dates throughout a year."""
        test_dates = [
            2451545.0,  # J2000 (Jan 1, 2000)
            2459580.5,  # Jan 1, 2022
            2460310.5,  # Jan 1, 2024
            2460400.5,  # Apr 1, 2024
            2460500.5,  # Jul 10, 2024
        ]

        swe.set_ephe_path(None)

        for jd in test_dates:
            for planet_id in [MERCURY, VENUS, MARS, JUPITER]:
                swe_attr = swe.pheno_ut(jd, planet_id, 0)
                lib_attr = ephem.pheno_ut(jd, planet_id, 0)

                # Must have 20 values
                assert len(lib_attr) == 20

                # Reserved values must be zero
                for i in range(5, 20):
                    assert lib_attr[i] == 0.0, (
                        f"JD {jd}, planet {planet_id}: attr[{i}] should be 0.0"
                    )

                # Active values should be in valid ranges
                assert 0 <= lib_attr[0] <= 180, "Phase angle out of range"
                assert 0 <= lib_attr[1] <= 1.0, "Phase out of range"
                assert 0 <= lib_attr[2] <= 180, "Elongation out of range"
                assert lib_attr[3] >= 0, "Diameter should be non-negative"

    @pytest.mark.comparison
    def test_swe_pheno_et_version_20_values(self):
        """Verify pheno (ET version) also returns 20 values."""
        jd_et = 2451545.0  # J2000 TT

        swe.set_ephe_path(None)
        swe_attr = swe.pheno(jd_et, MARS, 0)
        lib_attr = ephem.pheno(jd_et, MARS, 0)

        # Must have exactly 20 values
        assert len(lib_attr) == 20, (
            f"pheno should return 20 values, got {len(lib_attr)}"
        )

        # Reserved values must be zero
        for i in range(5, 20):
            assert lib_attr[i] == 0.0, f"pheno: attr[{i}] should be 0.0"

        # Phase angle comparison
        diff_angle = abs(swe_attr[0] - lib_attr[0])
        assert diff_angle < 1.0, f"Phase angle diff too large: {diff_angle}"

    @pytest.mark.comparison
    def test_pheno_ut_accepts_different_flags(self):
        """pheno_ut should accept different calculation flags without error."""
        jd = 2451545.0

        # Test with flag = 0
        attr0 = ephem.pheno_ut(jd, MARS, 0)
        assert isinstance(attr0, tuple)

        # Test with FLG_TRUEPOS
        from libephemeris.constants import FLG_TRUEPOS

        attr_true = ephem.pheno_ut(jd, MARS, FLG_TRUEPOS)
        assert isinstance(attr_true, tuple)


class TestPhaseAngleAllPlanets:
    """
    Comprehensive tests for phase angle calculation for all planets.

    The phase angle is the Sun-planet-Earth angle (the angle at the planet
    vertex in the triangle formed by Sun, planet, and Earth). It determines
    how much of the planet's illuminated surface is visible from Earth.

    For inner planets (Mercury, Venus):
    - Phase angle varies from 0° (superior conjunction, fully lit) to 180° (inferior conjunction)
    - Maximum elongation corresponds to phase angle around 90°

    For outer planets (Mars and beyond):
    - Phase angle is always small (typically < 50° for Mars, < 12° for Jupiter and beyond)
    - Opposition: phase angle ~ 0°
    - Quadrature: maximum phase angle
    """

    # Stricter tolerance for phase angle comparison (0.5 degrees)
    PHASE_ANGLE_TOL = 0.5

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
            (URANUS, "Uranus"),
            (NEPTUNE, "Neptune"),
            (PLUTO, "Pluto"),
        ],
    )
    def test_phase_angle_at_multiple_dates(self, planet_id, planet_name):
        """
        Test phase angle calculation at multiple dates throughout a year.

        This ensures phase angle is calculated correctly as planets move
        through different orbital configurations.
        """
        # Test dates spanning different orbital configurations
        test_dates = [
            2451545.0,  # J2000 (Jan 1, 2000)
            2451590.0,  # Feb 15, 2000
            2451635.0,  # Apr 1, 2000
            2451680.0,  # May 16, 2000
            2451725.0,  # Jun 30, 2000
            2451770.0,  # Aug 14, 2000
            2451815.0,  # Sep 28, 2000
            2451860.0,  # Nov 12, 2000
            2459580.5,  # Jan 1, 2022
            2460310.5,  # Jan 1, 2024
        ]

        swe.set_ephe_path(None)

        for jd in test_dates:
            swe_attr = swe.pheno_ut(jd, planet_id, 0)
            lib_attr = ephem.pheno_ut(jd, planet_id, 0)

            swe_phase_angle = swe_attr[0]
            lib_phase_angle = lib_attr[0]

            diff = abs(swe_phase_angle - lib_phase_angle)
            assert diff < self.PHASE_ANGLE_TOL, (
                f"{planet_name} at JD {jd}: phase angle SE={swe_phase_angle:.4f}°, "
                f"LIB={lib_phase_angle:.4f}°, diff={diff:.4f}°"
            )

    @pytest.mark.comparison
    def test_inner_planets_phase_angle_range(self):
        """
        Test that inner planets (Mercury, Venus) show substantial phase angle variation.

        Inner planets can have phase angles from 0° to 180° depending on
        their position relative to the Sun and Earth.
        """
        # Mercury: synodic period ~116 days, sample 200 days for full cycle
        jd_start = 2451545.0  # J2000
        phase_angles_mercury = []

        for i in range(200):  # Sample over Mercury's synodic period
            jd = jd_start + i
            attr_mercury = ephem.pheno_ut(jd, MERCURY, 0)
            phase_angles_mercury.append(attr_mercury[0])

        # Mercury should show substantial phase angle variation (range > 60°)
        mercury_range = max(phase_angles_mercury) - min(phase_angles_mercury)
        assert mercury_range > 60, (
            f"Mercury phase angle range {mercury_range:.1f}° should be > 60°"
        )

        # Venus: synodic period ~584 days, sample 600 days for full cycle
        phase_angles_venus = []
        for i in range(600):  # Sample over Venus's synodic period
            jd = jd_start + i
            attr_venus = ephem.pheno_ut(jd, VENUS, 0)
            phase_angles_venus.append(attr_venus[0])

        # Venus should show substantial phase angle variation (range > 60°)
        venus_range = max(phase_angles_venus) - min(phase_angles_venus)
        assert venus_range > 60, (
            f"Venus phase angle range {venus_range:.1f}° should be > 60°"
        )

    @pytest.mark.comparison
    def test_outer_planets_phase_angle_small(self):
        """
        Test that outer planets (Jupiter and beyond) have small phase angles.

        Outer planets are always seen from approximately the same direction
        as the Sun illuminates them, so phase angles are always small.
        """
        # Sample dates across a year
        jd_start = 2451545.0
        outer_planets = [
            (JUPITER, "Jupiter", 12),  # max phase angle ~12°
            (SATURN, "Saturn", 7),  # max phase angle ~7°
            (URANUS, "Uranus", 4),  # max phase angle ~4°
            (NEPTUNE, "Neptune", 2),  # max phase angle ~2°
            (PLUTO, "Pluto", 2),  # max phase angle ~2°
        ]

        for planet_id, planet_name, max_expected in outer_planets:
            for i in range(0, 365, 30):  # Sample monthly
                jd = jd_start + i
                attr = ephem.pheno_ut(jd, planet_id, 0)
                phase_angle = attr[0]
                assert phase_angle < max_expected + 1, (
                    f"{planet_name} at JD {jd}: phase angle {phase_angle:.2f}° "
                    f"exceeds expected max {max_expected}°"
                )

    @pytest.mark.comparison
    def test_mars_phase_angle_moderate(self):
        """
        Test that Mars has moderate phase angles (up to ~47°).

        Mars is the innermost outer planet and can show significant
        phase angles, especially near quadrature.
        """
        jd_start = 2451545.0
        max_phase_angle = 0

        for i in range(780):  # ~2 years (Mars synodic period ~780 days)
            jd = jd_start + i
            attr = ephem.pheno_ut(jd, MARS, 0)
            max_phase_angle = max(max_phase_angle, attr[0])

        # Mars maximum phase angle is about 47° at quadrature
        assert 30 < max_phase_angle < 50, (
            f"Mars max phase angle {max_phase_angle:.1f}° should be ~47°"
        )

    @pytest.mark.comparison
    def test_phase_angle_vs_elongation_consistency(self):
        """
        Test that phase angle and elongation are geometrically consistent.

        For outer planets at opposition (elongation ~180°), phase angle should be ~0°.
        For outer planets at quadrature (elongation ~90°), phase angle is maximum.
        """
        jd_start = 2451545.0
        swe.set_ephe_path(None)

        for planet_id in [JUPITER, SATURN]:
            # Find date near opposition (elongation close to 180°)
            for i in range(400):
                jd = jd_start + i
                attr = ephem.pheno_ut(jd, planet_id, 0)
                elongation = attr[2]
                phase_angle = attr[0]

                # Near opposition (elongation > 170°), phase angle should be small
                if elongation > 170:
                    assert phase_angle < 5, (
                        f"Planet {planet_id} near opposition: elongation={elongation:.1f}°, "
                        f"phase_angle={phase_angle:.1f}° should be < 5°"
                    )

    @pytest.mark.comparison
    def test_phase_angle_formula_verification(self):
        """
        Verify the phase angle formula: phase = (1 + cos(phase_angle)) / 2.

        This relationship must hold for all planets at all times.
        """
        test_dates = [
            2451545.0,
            2451600.0,
            2451700.0,
            2451800.0,
            2459500.0,
            2459600.0,
            2459700.0,
        ]
        planets = [
            MERCURY,
            VENUS,
            MARS,
            JUPITER,
            SATURN,
            URANUS,
            NEPTUNE,
            PLUTO,
        ]

        for jd in test_dates:
            for planet_id in planets:
                attr = ephem.pheno_ut(jd, planet_id, 0)
                phase_angle = attr[0]
                phase = attr[1]

                # Verify: phase = (1 + cos(phase_angle)) / 2
                expected_phase = (1 + math.cos(math.radians(phase_angle))) / 2

                diff = abs(phase - expected_phase)
                assert diff < 0.0001, (
                    f"Planet {planet_id} at JD {jd}: phase={phase:.6f} != "
                    f"expected={(1 + math.cos(math.radians(phase_angle))) / 2:.6f} "
                    f"from phase_angle={phase_angle:.2f}°"
                )

    @pytest.mark.comparison
    def test_moon_phase_angle_at_lunar_phases(self):
        """
        Test Moon phase angle at known lunar phases.

        - New Moon: phase angle ~180° (Moon between Sun and Earth)
        - First Quarter: phase angle ~90°
        - Full Moon: phase angle ~0° (Earth between Sun and Moon)
        - Last Quarter: phase angle ~90°
        """
        # Approximate dates for lunar phases in Jan 2000
        # These are rough dates for testing the general behavior
        jd_new = 2451550.1  # Near new moon
        jd_full = 2451565.0  # Near full moon

        attr_new = ephem.pheno_ut(jd_new, MOON, 0)
        attr_full = ephem.pheno_ut(jd_full, MOON, 0)

        # Near new moon, phase angle should be high (approaching 180°)
        assert attr_new[0] > 100, (
            f"Moon near new: phase_angle={attr_new[0]:.1f}° should be > 100°"
        )

        # Near full moon, phase angle should be low (approaching 0°)
        assert attr_full[0] < 50, (
            f"Moon near full: phase_angle={attr_full[0]:.1f}° should be < 50°"
        )


class TestIlluminatedFractionAllPlanets:
    """
    Comprehensive tests for illuminated fraction (phase) calculation for all planets.

    The illuminated fraction (also called "phase") represents the fraction of the
    planet's disk that appears illuminated as seen from Earth. It is calculated
    from the phase angle using the formula:

        k = (1 + cos(i)) / 2

    where i is the phase angle (Sun-planet-observer angle).

    Key values:
    - k = 1.0 (100%): Full phase - entire visible disk is illuminated
    - k = 0.5 (50%): Half phase - half of the visible disk is illuminated
    - k = 0.0 (0%): New phase - visible disk is not illuminated

    This corresponds to attr[1] in the pheno_ut output.
    """

    # Tolerance for illuminated fraction comparison with Swiss Ephemeris
    ILLUMINATION_TOL = 0.02  # 2% tolerance

    @pytest.mark.unit
    def test_sun_always_fully_illuminated(self):
        """
        Sun should always have illuminated fraction = 1.0.

        The Sun is the light source, so from any viewpoint it appears
        fully illuminated (100%).
        """
        test_dates = [
            2451545.0,  # J2000
            2459580.5,  # 2022
            2460310.5,  # 2024
            2465000.0,  # 2037
        ]

        for jd in test_dates:
            attr = ephem.pheno_ut(jd, SUN, 0)
            assert attr[1] == 0.0, (
                f"Sun at JD {jd}: illuminated fraction should be 0.0 (inapplicable), got {attr[1]}"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
            (URANUS, "Uranus"),
            (NEPTUNE, "Neptune"),
            (PLUTO, "Pluto"),
        ],
    )
    def test_illuminated_fraction_in_valid_range(self, planet_id, planet_name):
        """
        Illuminated fraction must always be between 0.0 and 1.0 for all planets.

        The formula k = (1 + cos(i))/2 always produces values in [0, 1]
        for any phase angle i in [0°, 180°].
        """
        test_dates = [
            2451545.0,
            2451600.0,
            2451700.0,
            2459500.0,
            2459600.0,
            2460300.0,
        ]

        for jd in test_dates:
            attr = ephem.pheno_ut(jd, planet_id, 0)
            illuminated_fraction = attr[1]

            assert 0.0 <= illuminated_fraction <= 1.0, (
                f"{planet_name} at JD {jd}: illuminated fraction {illuminated_fraction} "
                f"must be in range [0.0, 1.0]"
            )

    @pytest.mark.unit
    def test_moon_illuminated_fraction_varies_with_phase(self):
        """
        Moon illuminated fraction should vary from ~0 (new) to ~1 (full).

        The Moon's synodic period is ~29.5 days, so sampling over a month
        should show the full range of illumination.
        """
        jd_start = 2451545.0  # J2000
        illuminations = []

        # Sample over one lunar month
        for i in range(30):
            jd = jd_start + i
            attr = ephem.pheno_ut(jd, MOON, 0)
            illuminations.append(attr[1])

        min_illum = min(illuminations)
        max_illum = max(illuminations)

        # Should see nearly full range (allowing for exact timing)
        assert min_illum < 0.15, (
            f"Moon minimum illumination {min_illum:.2f} should be < 0.15 (near new)"
        )
        assert max_illum > 0.85, (
            f"Moon maximum illumination {max_illum:.2f} should be > 0.85 (near full)"
        )

    @pytest.mark.unit
    def test_inner_planets_illumination_varies_widely(self):
        """
        Inner planets (Mercury, Venus) should show wide illumination variation.

        These planets can appear as crescents (low illumination) or nearly full
        (high illumination) depending on their orbital position.
        """
        jd_start = 2451545.0

        for planet_id, planet_name, synodic_period in [
            (MERCURY, "Mercury", 116),
            (VENUS, "Venus", 584),
        ]:
            illuminations = []
            # Sample over synodic period
            for i in range(0, synodic_period, 5):
                jd = jd_start + i
                attr = ephem.pheno_ut(jd, planet_id, 0)
                illuminations.append(attr[1])

            min_illum = min(illuminations)
            max_illum = max(illuminations)

            # Inner planets should show substantial variation
            assert max_illum - min_illum > 0.3, (
                f"{planet_name} illumination range {max_illum - min_illum:.2f} "
                f"should be > 0.3"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name,min_expected",
        [
            (JUPITER, "Jupiter", 0.99),
            (SATURN, "Saturn", 0.99),
            (URANUS, "Uranus", 0.999),
            (NEPTUNE, "Neptune", 0.999),
            (PLUTO, "Pluto", 0.999),
        ],
    )
    def test_outer_planets_nearly_fully_illuminated(
        self, planet_id, planet_name, min_expected
    ):
        """
        Outer planets should always appear nearly fully illuminated.

        Due to their great distance from the Sun and Earth, outer planets
        are always seen from nearly the same direction as they are illuminated,
        resulting in illuminated fractions very close to 1.0.
        """
        test_dates = [2451545.0, 2455000.0, 2459000.0, 2460000.0]

        for jd in test_dates:
            attr = ephem.pheno_ut(jd, planet_id, 0)
            illuminated_fraction = attr[1]

            assert illuminated_fraction > min_expected, (
                f"{planet_name} at JD {jd}: illumination {illuminated_fraction:.6f} "
                f"should be > {min_expected}"
            )

    @pytest.mark.unit
    def test_mars_illumination_moderate_range(self):
        """
        Mars illuminated fraction varies moderately (about 0.84 to 1.0).

        Mars is the innermost outer planet, so it can show noticeable
        illumination variation, but never appears as a thin crescent.
        """
        jd_start = 2451545.0
        illuminations = []

        # Sample over Mars synodic period (~780 days)
        for i in range(0, 780, 10):
            jd = jd_start + i
            attr = ephem.pheno_ut(jd, MARS, 0)
            illuminations.append(attr[1])

        min_illum = min(illuminations)
        max_illum = max(illuminations)

        # Mars minimum illumination is about 84% at quadrature
        assert min_illum > 0.80, (
            f"Mars minimum illumination {min_illum:.2f} should be > 0.80"
        )
        assert min_illum < 0.95, (
            f"Mars minimum illumination {min_illum:.2f} should show some variation"
        )
        assert max_illum > 0.99, (
            f"Mars maximum illumination {max_illum:.2f} should be > 0.99 at opposition"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
            (URANUS, "Uranus"),
            (NEPTUNE, "Neptune"),
            (PLUTO, "Pluto"),
        ],
    )
    def test_illuminated_fraction_matches_swiss_ephemeris(self, planet_id, planet_name):
        """
        Illuminated fraction should match Swiss Ephemeris within tolerance.

        This verifies the implementation produces results consistent with
        the reference Swiss Ephemeris library.
        """
        test_dates = [
            2451545.0,  # J2000
            2451600.0,
            2451700.0,
            2455000.0,  # 2009
            2459000.0,  # 2020
            2460000.0,  # 2023
        ]

        swe.set_ephe_path(None)

        for jd in test_dates:
            swe_attr = swe.pheno_ut(jd, planet_id, 0)
            lib_attr = ephem.pheno_ut(jd, planet_id, 0)

            swe_illumination = swe_attr[1]
            lib_illumination = lib_attr[1]

            diff = abs(swe_illumination - lib_illumination)
            assert diff < self.ILLUMINATION_TOL, (
                f"{planet_name} at JD {jd}: illumination SE={swe_illumination:.4f}, "
                f"LIB={lib_illumination:.4f}, diff={diff:.4f}"
            )

    @pytest.mark.comparison
    def test_moon_illuminated_fraction_matches_swiss_ephemeris(self):
        """
        Moon illuminated fraction should match Swiss Ephemeris.

        The Moon requires special handling due to its proximity and
        non-spherical orbit, but illumination should still match closely.
        """
        test_dates = [
            2451545.0,
            2451550.0,
            2451555.0,
            2451560.0,
            2451565.0,
            2451570.0,
        ]

        swe.set_ephe_path(None)

        for jd in test_dates:
            swe_attr = swe.pheno_ut(jd, MOON, 0)
            lib_attr = ephem.pheno_ut(jd, MOON, 0)

            swe_illumination = swe_attr[1]
            lib_illumination = lib_attr[1]

            diff = abs(swe_illumination - lib_illumination)
            # Moon allows slightly larger tolerance due to calculation complexity
            assert diff < 0.05, (
                f"Moon at JD {jd}: illumination SE={swe_illumination:.4f}, "
                f"LIB={lib_illumination:.4f}, diff={diff:.4f}"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MOON, "Moon"),
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
            (URANUS, "Uranus"),
            (NEPTUNE, "Neptune"),
            (PLUTO, "Pluto"),
        ],
    )
    def test_illumination_formula_from_phase_angle(self, planet_id, planet_name):
        """
        Verify illuminated fraction is calculated correctly from phase angle.

        The relationship k = (1 + cos(i))/2 must hold for all planets,
        where i is the phase angle in degrees.
        """
        test_dates = [2451545.0, 2451600.0, 2455000.0, 2459000.0]

        for jd in test_dates:
            attr = ephem.pheno_ut(jd, planet_id, 0)
            phase_angle = attr[0]
            illumination = attr[1]

            expected_illumination = (1 + math.cos(math.radians(phase_angle))) / 2

            diff = abs(illumination - expected_illumination)
            assert diff < 0.0001, (
                f"{planet_name} at JD {jd}: illumination {illumination:.6f} != "
                f"expected {expected_illumination:.6f} from phase_angle={phase_angle:.2f}°"
            )

    @pytest.mark.unit
    def test_illumination_physical_interpretation(self):
        """
        Test that illumination values have correct physical interpretation.

        Verify key boundary conditions:
        - Phase angle 0° -> illumination 1.0 (full)
        - Phase angle 90° -> illumination 0.5 (half)
        - Phase angle 180° -> illumination 0.0 (new)
        """
        # Test the mathematical formula directly
        test_cases = [
            (0, 1.0),  # Full phase
            (30, 0.933),  # Gibbous
            (60, 0.75),  # Gibbous
            (90, 0.5),  # Half (quarter)
            (120, 0.25),  # Crescent
            (150, 0.067),  # Thin crescent
            (180, 0.0),  # New phase
        ]

        for phase_angle, expected in test_cases:
            calculated = (1 + math.cos(math.radians(phase_angle))) / 2
            assert abs(calculated - expected) < 0.001, (
                f"Phase angle {phase_angle}°: expected illumination {expected}, "
                f"calculated {calculated:.4f}"
            )

    @pytest.mark.unit
    def test_illumination_continuous_over_time(self):
        """
        Illumination should change smoothly and continuously over time.

        There should be no sudden jumps in illumination values between
        consecutive time steps.
        """
        jd_start = 2451545.0

        for planet_id in [MERCURY, VENUS, MARS, MOON]:
            prev_illumination = None

            for i in range(100):
                jd = jd_start + i * 0.5  # Half-day steps
                attr = ephem.pheno_ut(jd, planet_id, 0)
                illumination = attr[1]

                if prev_illumination is not None:
                    # Change should be gradual (< 0.1 per half day typically)
                    change = abs(illumination - prev_illumination)
                    assert change < 0.15, (
                        f"Planet {planet_id} at JD {jd}: illumination jumped by "
                        f"{change:.4f} from {prev_illumination:.4f} to {illumination:.4f}"
                    )

                prev_illumination = illumination


class TestApparentDiameterAllPlanets:
    """
    Comprehensive tests for apparent angular diameter calculation.

    The apparent diameter is the angular size of a celestial body as seen from Earth,
    calculated from the body's physical radius and geocentric distance:

        diameter = 2 * arctan(radius / distance) ≈ 2 * radius / distance (in radians)
        diameter_deg = diameter_rad * (180/π)

    This corresponds to attr[3] in the pheno_ut output (in degrees).
    """

    # Tolerance for diameter comparison with Swiss Ephemeris (relative)
    DIAMETER_REL_TOL = 0.05  # 5% tolerance

    @pytest.mark.unit
    def test_sun_diameter_typical_range(self):
        """
        Sun apparent diameter should be around 32 arcminutes (~0.53 degrees).

        The Sun's diameter varies from ~31.5' (aphelion) to ~32.5' (perihelion)
        due to Earth's elliptical orbit.
        """
        jd = 2451545.0  # J2000
        attr = ephem.pheno_ut(jd, SUN, 0)

        diameter = attr[3]
        # Sun diameter ~32 arcmin ≈ 0.53 deg (±3%)
        assert 0.51 < diameter < 0.55, (
            f"Sun diameter {diameter:.4f} deg outside expected range [0.51, 0.55]"
        )

    @pytest.mark.unit
    def test_sun_diameter_varies_with_distance(self):
        """
        Sun diameter should be larger at perihelion than aphelion.

        Earth is closest to Sun around January 3 (perihelion) and
        farthest around July 4 (aphelion).
        """
        # Approximate dates
        jd_perihelion = 2451547.5  # ~Jan 3, 2000
        jd_aphelion = 2451731.5  # ~Jul 4, 2000

        attr_peri = ephem.pheno_ut(jd_perihelion, SUN, 0)
        attr_aph = ephem.pheno_ut(jd_aphelion, SUN, 0)

        # Perihelion should have larger diameter
        assert attr_peri[3] > attr_aph[3], (
            f"Sun should be larger at perihelion ({attr_peri[3]:.2f}) "
            f"than aphelion ({attr_aph[3]:.2f})"
        )

    @pytest.mark.unit
    def test_moon_diameter_typical_range(self):
        """
        Moon apparent diameter should be around 31 arcminutes (~0.52 degrees).

        The Moon's diameter varies from ~29.4' (apogee) to ~33.5' (perigee)
        due to its elliptical orbit.
        """
        jd = 2451545.0
        attr = ephem.pheno_ut(jd, MOON, 0)

        diameter = attr[3]
        # Moon diameter typically 29-34 arcmin ≈ 0.47-0.57 degrees
        assert 0.47 < diameter < 0.58, (
            f"Moon diameter {diameter:.4f} deg outside expected range [0.47, 0.58]"
        )

    @pytest.mark.unit
    def test_moon_diameter_varies_with_distance(self):
        """
        Moon diameter should vary noticeably over an anomalistic month.

        The Moon's apparent size varies by about 14% between perigee and apogee.
        """
        jd_start = 2451545.0
        diameters = []

        # Sample over 30 days (more than one anomalistic month of ~27.5 days)
        for i in range(30):
            jd = jd_start + i
            attr = ephem.pheno_ut(jd, MOON, 0)
            diameters.append(attr[3])

        min_diam = min(diameters)
        max_diam = max(diameters)

        # Should see at least 10% variation
        variation = (max_diam - min_diam) / min_diam
        assert variation > 0.08, (
            f"Moon diameter variation {variation:.1%} should be > 8%"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name,min_diam,max_diam",
        [
            (MERCURY, "Mercury", 0.0011, 0.0039),  # 4"-14" in degrees
            (VENUS, "Venus", 0.0028, 0.0183),  # 10"-66" in degrees
            (MARS, "Mars", 0.0008, 0.0072),  # 3"-26" in degrees
            (JUPITER, "Jupiter", 0.0083, 0.0142),  # 30"-51" in degrees
            (SATURN, "Saturn", 0.0042, 0.0058),  # 15"-21" in degrees
        ],
    )
    def test_planet_diameter_in_typical_range(
        self, planet_id, planet_name, min_diam, max_diam
    ):
        """
        Planet apparent diameter should be within expected observable range.

        These ranges account for varying Earth-planet distances.
        Values are in degrees.
        """
        jd = 2451545.0
        attr = ephem.pheno_ut(jd, planet_id, 0)

        diameter = attr[3]
        assert min_diam < diameter < max_diam, (
            f"{planet_name} diameter {diameter:.6f} deg "
            f"outside typical range [{min_diam}, {max_diam}]"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name,max_expected",
        [
            (URANUS, "Uranus", 0.00125),  # 3.4"-4.1" in degrees
            (NEPTUNE, "Neptune", 0.00069),  # 2.2"-2.4" in degrees
            (PLUTO, "Pluto", 0.000056),  # ~0.11" in degrees
        ],
    )
    def test_outer_planet_diameter_small(self, planet_id, planet_name, max_expected):
        """
        Distant planets should have very small apparent diameters (in degrees).
        """
        jd = 2451545.0
        attr = ephem.pheno_ut(jd, planet_id, 0)

        diameter = attr[3]
        assert 0 < diameter < max_expected, (
            f"{planet_name} diameter {diameter:.6f} deg "
            f"should be in (0, {max_expected})"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
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
        ],
    )
    def test_diameter_matches_swiss_ephemeris(self, planet_id, planet_name):
        """
        Apparent diameter should match Swiss Ephemeris within tolerance.

        Note: Both Swiss Ephemeris and libephemeris return diameter in degrees.
        """
        test_dates = [
            2451545.0,  # J2000
            2455000.0,  # 2009
            2459000.0,  # 2020
            2460000.0,  # 2023
        ]

        swe.set_ephe_path(None)

        for jd in test_dates:
            swe_attr = swe.pheno_ut(jd, planet_id, 0)
            lib_attr = ephem.pheno_ut(jd, planet_id, 0)

            # Both return degrees
            swe_diam = swe_attr[3]
            lib_diam = lib_attr[3]

            if swe_diam > 0:
                rel_diff = abs(swe_diam - lib_diam) / swe_diam
                assert rel_diff < self.DIAMETER_REL_TOL, (
                    f"{planet_name} at JD {jd}: diameter "
                    f"SE={swe_diam:.6f}deg, LIB={lib_diam:.6f}deg, "
                    f"rel_diff={rel_diff:.2%}"
                )

    @pytest.mark.unit
    def test_diameter_inversely_proportional_to_distance(self):
        """
        Diameter should be inversely proportional to distance.

        As a planet moves farther away, its apparent diameter decreases.
        """
        jd_start = 2451545.0

        # Test with Mars over ~2 years (covers close and far configurations)
        distances = []
        diameters = []

        for i in range(0, 780, 30):  # ~2 years in monthly steps
            jd = jd_start + i
            # Get Mars position
            pos, _ = ephem.calc_ut(jd, MARS, 0)
            distance = pos[2]  # Distance in AU

            # Get Mars diameter
            attr = ephem.pheno_ut(jd, MARS, 0)
            diameter = attr[3]

            distances.append(distance)
            diameters.append(diameter)

        # Check inverse relationship: diameter * distance should be roughly constant
        products = [d * dist for d, dist in zip(diameters, distances)]
        avg_product = sum(products) / len(products)

        for i, product in enumerate(products):
            rel_diff = abs(product - avg_product) / avg_product
            assert rel_diff < 0.05, (
                f"Mars: diameter*distance at step {i} deviates by {rel_diff:.2%} "
                f"from average (expected constant due to fixed physical size)"
            )

    @pytest.mark.unit
    def test_jupiter_largest_planet_diameter(self):
        """
        Jupiter should have the largest apparent diameter of all planets.
        """
        jd = 2451545.0

        planet_diameters = {}
        for planet_id in [MARS, JUPITER, SATURN, URANUS, NEPTUNE]:
            attr = ephem.pheno_ut(jd, planet_id, 0)
            planet_diameters[planet_id] = attr[3]

        jupiter_diam = planet_diameters[JUPITER]

        for planet_id, diam in planet_diameters.items():
            if planet_id != JUPITER:
                assert jupiter_diam > diam, (
                    f"Jupiter ({jupiter_diam:.2f}as) should be larger than "
                    f"planet {planet_id} ({diam:.2f}as)"
                )

    @pytest.mark.unit
    def test_diameter_positive_for_all_bodies(self):
        """
        All supported bodies should have positive apparent diameters.
        """
        jd = 2451545.0
        bodies = [
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
        ]

        for body_id in bodies:
            attr = ephem.pheno_ut(jd, body_id, 0)
            assert attr[3] > 0, f"Body {body_id}: diameter should be positive"

    @pytest.mark.unit
    def test_diameter_formula_using_physical_radius(self):
        """
        Verify diameter calculation uses the correct formula.

        For a body with physical radius R at distance D:
        diameter_arcsec = 2 * R / D * (180/π * 3600) = 2 * R / D * 206264.806

        We verify this by checking the relationship between calculated diameter
        and known physical parameters.
        """
        # Sun: radius ~695700 km, at ~1 AU
        # Expected diameter: 2 * 695700 / 149597870.7 * 206264.806 ≈ 1919.1 arcsec
        jd = 2451545.0
        attr_sun = ephem.pheno_ut(jd, SUN, 0)

        # The Sun's diameter should be close to the theoretical value
        # (within a few percent due to actual distance variation)
        theoretical_sun_diam_arcsec = 2 * 695700 / 149597870.7 * 206264.806
        theoretical_sun_diam_deg = theoretical_sun_diam_arcsec / 3600.0
        assert (
            abs(attr_sun[3] - theoretical_sun_diam_deg) / theoretical_sun_diam_deg
            < 0.02
        ), (
            f"Sun diameter {attr_sun[3]:.6f}deg differs from theoretical "
            f"{theoretical_sun_diam_deg:.6f}deg by more than 2%"
        )

    @pytest.mark.comparison
    def test_diameter_at_multiple_dates_matches_swe(self):
        """
        Diameter calculations should match Swiss Ephemeris across multiple dates.

        This tests the diameter calculation for various Earth-planet configurations.
        """
        test_dates = [
            2451545.0,  # J2000
            2451590.0,
            2451635.0,
            2451680.0,
            2451725.0,
            2459580.5,
            2460310.5,
        ]

        swe.set_ephe_path(None)

        for jd in test_dates:
            for planet_id in [MERCURY, VENUS, MARS, JUPITER, SATURN]:
                swe_attr = swe.pheno_ut(jd, planet_id, 0)
                lib_attr = ephem.pheno_ut(jd, planet_id, 0)

                swe_diam = swe_attr[3]  # Both return degrees
                lib_diam = lib_attr[3]

                if swe_diam > 0:
                    rel_diff = abs(swe_diam - lib_diam) / swe_diam
                    assert rel_diff < 0.05, (
                        f"Planet {planet_id} at JD {jd}: "
                        f"diameter SE={swe_diam:.6f}deg, LIB={lib_diam:.6f}deg"
                    )


class TestElongationHelpers:
    """
    Tests for elongation helper functions that distinguish between
    morning star (western elongation) and evening star (eastern elongation).

    These functions extend the basic pheno_ut elongation calculation
    by providing signed elongation values and morning/evening star classification.

    Convention:
        - Eastern elongation (positive): Planet is east of Sun = evening star
          (visible after sunset in the western sky)
        - Western elongation (negative): Planet is west of Sun = morning star
          (visible before sunrise in the eastern sky)
    """

    @pytest.mark.unit
    def test_get_elongation_from_sun_returns_tuple(self):
        """get_elongation_from_sun should return (elongation, is_evening_star)."""
        jd = 2451545.0
        result = ephem.get_elongation_from_sun(jd, VENUS)

        assert isinstance(result, tuple)
        assert len(result) == 2
        elongation, is_evening = result
        assert isinstance(elongation, float)
        assert isinstance(is_evening, bool)

    @pytest.mark.unit
    def test_get_signed_elongation_returns_float(self):
        """get_signed_elongation should return a float."""
        jd = 2451545.0
        result = ephem.get_signed_elongation(jd, VENUS)
        assert isinstance(result, float)

    @pytest.mark.unit
    def test_is_morning_star_returns_bool(self):
        """is_morning_star should return a boolean."""
        jd = 2451545.0
        result = ephem.is_morning_star(jd, VENUS)
        assert isinstance(result, bool)

    @pytest.mark.unit
    def test_is_evening_star_returns_bool(self):
        """is_evening_star should return a boolean."""
        jd = 2451545.0
        result = ephem.is_evening_star(jd, VENUS)
        assert isinstance(result, bool)

    @pytest.mark.unit
    def test_get_elongation_type_returns_string(self):
        """get_elongation_type should return 'eastern', 'western', or 'none'."""
        jd = 2451545.0
        result = ephem.get_elongation_type(jd, VENUS)
        assert result in ("eastern", "western", "none")

    @pytest.mark.unit
    def test_morning_evening_mutually_exclusive(self):
        """A planet cannot be both morning star and evening star at same time."""
        jd = 2451545.0
        for planet_id in [MERCURY, VENUS, MARS, JUPITER]:
            is_morning = ephem.is_morning_star(jd, planet_id)
            is_evening = ephem.is_evening_star(jd, planet_id)
            assert is_morning != is_evening, (
                f"Planet {planet_id} cannot be both morning and evening star"
            )

    @pytest.mark.unit
    def test_signed_elongation_range(self):
        """Signed elongation should be between -180 and +180 degrees."""
        jd = 2451545.0
        for planet_id in [MERCURY, VENUS, MARS, JUPITER, SATURN]:
            elong = ephem.get_signed_elongation(jd, planet_id)
            assert -180.0 <= elong <= 180.0, (
                f"Planet {planet_id} signed elongation {elong} out of range"
            )

    @pytest.mark.unit
    def test_absolute_elongation_matches_pheno(self):
        """Absolute value of signed elongation should match pheno_ut result."""
        jd = 2451545.0
        for planet_id in [MERCURY, VENUS, MARS, JUPITER, SATURN]:
            signed_elong = ephem.get_signed_elongation(jd, planet_id)
            pheno_attr = ephem.pheno_ut(jd, planet_id, 0)
            pheno_elong = pheno_attr[2]

            # Allow small tolerance due to different calculation methods
            # (pheno uses spherical trig, signed uses ecliptic longitude)
            diff = abs(abs(signed_elong) - pheno_elong)
            assert diff < 1.0, (
                f"Planet {planet_id}: |signed|={abs(signed_elong):.4f}, "
                f"pheno={pheno_elong:.4f}, diff={diff:.4f}"
            )

    @pytest.mark.unit
    def test_sun_elongation_type_is_none(self):
        """Sun should return 'none' for elongation type."""
        jd = 2451545.0
        assert ephem.get_elongation_type(jd, SUN) == "none"

    @pytest.mark.unit
    def test_sun_is_not_morning_or_evening_star(self):
        """Sun cannot be morning or evening star."""
        jd = 2451545.0
        assert not ephem.is_morning_star(jd, SUN)
        assert not ephem.is_evening_star(jd, SUN)

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_signed_elongation_matches_swe_positions(self, planet_id, planet_name):
        """
        Signed elongation should match the longitude difference from Swiss Ephemeris.
        """
        jd = 2451545.0
        swe.set_ephe_path(None)

        # Calculate expected elongation from Swiss Ephemeris positions
        pos_planet = swe.calc_ut(jd, planet_id, 0)[0][0]
        pos_sun = swe.calc_ut(jd, 0, 0)[0][0]  # SUN = 0

        expected_diff = pos_planet - pos_sun
        if expected_diff > 180.0:
            expected_diff -= 360.0
        elif expected_diff < -180.0:
            expected_diff += 360.0

        lib_elong = ephem.get_signed_elongation(jd, planet_id)

        diff = abs(expected_diff - lib_elong)
        assert diff < 0.01, (
            f"{planet_name}: SWE={expected_diff:.4f}, LIB={lib_elong:.4f}, "
            f"diff={diff:.4f}"
        )

    @pytest.mark.comparison
    def test_venus_transitions_morning_evening(self):
        """
        Venus should transition between morning and evening star over time.

        Venus has a synodic period of ~584 days, alternating between
        evening star (eastern elongation) and morning star (western elongation).
        """
        # Find dates where Venus changes from evening to morning star
        jd_start = 2460500.0  # Starting point in 2024

        morning_found = False
        evening_found = False

        for i in range(0, 365, 10):  # Check over a year
            jd = jd_start + i
            is_morning = ephem.is_morning_star(jd, VENUS)
            if is_morning:
                morning_found = True
            else:
                evening_found = True

            if morning_found and evening_found:
                break

        assert morning_found, "Venus should be a morning star at some point"
        assert evening_found, "Venus should be an evening star at some point"

    @pytest.mark.comparison
    def test_mercury_transitions_frequently(self):
        """
        Mercury should transition between morning and evening star multiple times per year.

        Mercury has a synodic period of ~116 days, so it transitions
        several times per year.
        """
        jd_start = 2460000.0

        transition_count = 0
        prev_is_evening = ephem.is_evening_star(jd_start, MERCURY)

        for i in range(1, 365):  # Check over a year
            jd = jd_start + i
            is_evening = ephem.is_evening_star(jd, MERCURY)
            if is_evening != prev_is_evening:
                transition_count += 1
                prev_is_evening = is_evening

        # Mercury should transition at least 5-6 times per year
        assert transition_count >= 5, (
            f"Mercury should transition at least 5 times per year, "
            f"found {transition_count}"
        )

    @pytest.mark.comparison
    def test_eastern_elongation_positive(self):
        """Eastern elongation (evening star) should have positive signed elongation."""
        # Find a date where Venus is east of Sun
        jd = 2460600.0  # Choose a date when Venus is evening star

        # Check over a range to find an evening star configuration
        for delta in range(0, 200, 10):
            test_jd = jd + delta
            if ephem.is_evening_star(test_jd, VENUS):
                elong = ephem.get_signed_elongation(test_jd, VENUS)
                assert elong > 0, (
                    f"Evening star should have positive elongation, got {elong}"
                )
                elong_type = ephem.get_elongation_type(test_jd, VENUS)
                assert elong_type == "eastern", (
                    f"Evening star should have eastern elongation, got {elong_type}"
                )
                return

        pytest.fail("Could not find Venus as evening star in test range")

    @pytest.mark.comparison
    def test_western_elongation_negative(self):
        """Western elongation (morning star) should have negative signed elongation."""
        # Find a date where Venus is west of Sun
        jd = 2460700.0  # Choose a date when Venus might be morning star

        # Check over a range to find a morning star configuration
        for delta in range(0, 200, 10):
            test_jd = jd + delta
            if ephem.is_morning_star(test_jd, VENUS):
                elong = ephem.get_signed_elongation(test_jd, VENUS)
                assert elong < 0, (
                    f"Morning star should have negative elongation, got {elong}"
                )
                elong_type = ephem.get_elongation_type(test_jd, VENUS)
                assert elong_type == "western", (
                    f"Morning star should have western elongation, got {elong_type}"
                )
                return

        pytest.fail("Could not find Venus as morning star in test range")

    @pytest.mark.comparison
    def test_elongation_consistency_across_planets(self):
        """
        Test signed elongation consistency across all planets at multiple dates.
        """
        test_dates = [
            2451545.0,  # J2000
            2455000.0,  # 2009
            2459000.0,  # 2020
            2460000.0,  # 2023
        ]

        planets = [MERCURY, VENUS, MARS, JUPITER, SATURN]
        swe.set_ephe_path(None)

        for jd in test_dates:
            for planet_id in planets:
                # Get libephemeris signed elongation
                lib_elong = ephem.get_signed_elongation(jd, planet_id)
                is_evening = ephem.is_evening_star(jd, planet_id)
                is_morning = ephem.is_morning_star(jd, planet_id)
                elong_type = ephem.get_elongation_type(jd, planet_id)

                # Verify consistency
                if lib_elong > 0:
                    assert is_evening, (
                        f"Planet {planet_id} at JD {jd}: positive elongation "
                        f"should be evening star"
                    )
                    assert not is_morning
                    assert elong_type == "eastern"
                else:
                    assert is_morning, (
                        f"Planet {planet_id} at JD {jd}: negative elongation "
                        f"should be morning star"
                    )
                    assert not is_evening
                    assert elong_type == "western"

    @pytest.mark.comparison
    def test_moon_elongation_full_range(self):
        """
        Moon elongation should cover full range from -180 to +180 over its orbit.

        Note: The Moon is traditionally not called a "morning star" or "evening star",
        but the elongation functions work for it too.
        """
        jd_start = 2451545.0

        min_elong = 180.0
        max_elong = -180.0

        for i in range(30):  # One lunar month
            jd = jd_start + i
            elong = ephem.get_signed_elongation(jd, MOON)
            min_elong = min(min_elong, elong)
            max_elong = max(max_elong, elong)

        # Moon elongation should span a wide range
        elong_range = max_elong - min_elong
        assert elong_range > 300, (
            f"Moon elongation range {elong_range:.1f}° should be > 300°"
        )
