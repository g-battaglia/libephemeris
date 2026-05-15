"""
Utility Functions Comparison Tests.

Compares utility functions between pyswisseph and libephemeris:
- Angle normalization functions
- Degree/centiseconds conversion
- String formatting functions
- Atmospheric refraction (refrac_extended)
"""

import pytest
import swisseph as swe
import libephemeris as ephem


# ============================================================================
# TOLERANCES
# ============================================================================

ANGLE_TOL = 1e-10
REFRACTION_TOL = 0.01  # degrees (refraction can vary slightly between implementations)
DIP_TOL = 0.001  # degrees (dip of horizon tolerance)


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestDegnorm:
    """Compare degnorm function."""

    TEST_VALUES = [
        0.0,
        90.0,
        180.0,
        270.0,
        360.0,
        -90.0,
        -180.0,
        -360.0,
        450.0,
        720.0,
        -720.0,
        359.9999,
        0.0001,
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("angle", TEST_VALUES)
    def test_degnorm(self, angle):
        """Test degnorm normalizes to 0-360."""
        result_swe = swe.degnorm(angle)
        result_py = ephem.degnorm(angle)

        diff = abs(result_swe - result_py)

        assert diff < ANGLE_TOL, f"degnorm({angle}): diff {diff:.15f} exceeds tolerance"

        # Result should be in [0, 360)
        assert 0 <= result_py < 360, (
            f"degnorm({angle}) = {result_py}, expected [0, 360)"
        )


class TestRadnorm:
    """Compare radnorm function."""

    import math

    PI = math.pi

    TEST_VALUES = [
        0.0,
        PI / 2,
        PI,
        3 * PI / 2,
        2 * PI,
        -PI / 2,
        -PI,
        -2 * PI,
        3 * PI,
        4 * PI,
        -4 * PI,
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("angle", TEST_VALUES)
    def test_radnorm(self, angle):
        """Test radnorm normalizes to 0-2pi."""
        result_swe = swe.radnorm(angle)
        result_py = ephem.radnorm(angle)

        diff = abs(result_swe - result_py)

        assert diff < ANGLE_TOL, f"radnorm({angle}): diff {diff:.15f} exceeds tolerance"


class TestCsnorm:
    """Compare csnorm function (centiseconds normalization)."""

    TEST_VALUES = [
        0,
        360 * 3600 * 100,
        -360 * 3600 * 100,
        180 * 3600 * 100,
        -180 * 3600 * 100,
        90 * 3600 * 100,
        450 * 3600 * 100,
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("cs", TEST_VALUES)
    def test_csnorm(self, cs):
        """Test csnorm normalizes centiseconds."""
        result_swe = swe.csnorm(cs)
        result_py = ephem.csnorm(cs)

        diff = abs(result_swe - result_py)

        assert diff < 1, f"csnorm({cs}): diff {diff} exceeds tolerance"


class TestD2l:
    """Compare d2l function (degrees to centiseconds)."""

    TEST_VALUES = [0.0, 1.0, 45.0, 90.0, 180.0, 270.0, 359.9999]

    @pytest.mark.comparison
    @pytest.mark.parametrize("deg", TEST_VALUES)
    def test_d2l(self, deg):
        """Test d2l conversion."""
        result_swe = swe.d2l(deg)
        result_py = ephem.d2l(deg)

        diff = abs(result_swe - result_py)

        assert diff < 1, f"d2l({deg}): diff {diff} exceeds tolerance"


class TestDifcsn:
    """Compare difcsn function (centisecond difference)."""

    TEST_PAIRS = [
        (0, 360 * 3600 * 100),
        (90 * 3600 * 100, 270 * 3600 * 100),
        (180 * 3600 * 100, 0),
        (350 * 3600 * 100, 10 * 3600 * 100),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("cs1,cs2", TEST_PAIRS)
    def test_difcsn(self, cs1, cs2):
        """Test difcsn calculates shortest arc difference."""
        result_swe = swe.difcsn(cs1, cs2)
        result_py = ephem.difcsn(cs1, cs2)

        diff = abs(result_swe - result_py)

        assert diff < 1, f"difcsn({cs1}, {cs2}): diff {diff} exceeds tolerance"


class TestDifdegn:
    """Compare difdegn function (degree difference)."""

    TEST_PAIRS = [
        (0.0, 360.0),
        (90.0, 270.0),
        (180.0, 0.0),
        (350.0, 10.0),
        (1.0, 359.0),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("deg1,deg2", TEST_PAIRS)
    def test_difdegn(self, deg1, deg2):
        """Test difdegn calculates shortest arc difference in degrees."""
        result_swe = swe.difdegn(deg1, deg2)
        result_py = ephem.difdegn(deg1, deg2)

        diff = abs(result_swe - result_py)

        assert diff < ANGLE_TOL, (
            f"difdegn({deg1}, {deg2}): diff {diff:.15f} exceeds tolerance"
        )


class TestSplitDeg:
    """Compare split_deg function."""

    TEST_VALUES = [
        (0.0, 0),
        (45.5, 0),
        (123.456789, 0),
        (359.999999, 0),
        (-45.5, 0),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("deg,roundflag", TEST_VALUES)
    def test_split_deg(self, deg, roundflag):
        """Test split_deg breaks degrees into components."""
        result_swe = swe.split_deg(deg, roundflag)
        result_py = ephem.split_deg(deg, roundflag)

        # Compare degree component
        assert result_swe[0] == result_py[0], (
            f"split_deg({deg}): degrees differ {result_swe[0]} vs {result_py[0]}"
        )
        # Compare minutes
        assert result_swe[1] == result_py[1], (
            f"split_deg({deg}): minutes differ {result_swe[1]} vs {result_py[1]}"
        )
        # Compare seconds (with small tolerance)
        diff_sec = abs(result_swe[2] - result_py[2])
        assert diff_sec < 0.001, f"split_deg({deg}): seconds diff {diff_sec:.6f}"


class TestRefracExtended:
    """Compare refrac_extended function with pyswisseph.

    Tests atmospheric refraction calculations with extended parameters:
    - Various altitude angles (horizon to zenith)
    - Different atmospheric pressures
    - Different temperatures
    - Different observer elevations
    - Various lapse rates
    """

    # Altitude angles to test (degrees)
    ALTITUDE_ANGLES = [
        0.0,  # Horizon - maximum refraction
        5.0,  # Very low altitude
        10.0,  # Low altitude
        30.0,  # Medium altitude
        45.0,  # Mid-sky
        60.0,  # Higher altitude
        90.0,  # Zenith - minimal refraction
    ]

    # Atmospheric pressures (hPa/mbar)
    PRESSURES = [
        900.0,  # Low pressure (high altitude location)
        1013.25,  # Standard sea level pressure
        1050.0,  # High pressure
    ]

    # Temperatures (Celsius)
    TEMPERATURES = [
        -20.0,  # Very cold
        0.0,  # Freezing
        15.0,  # Standard
        35.0,  # Hot
    ]

    # Observer elevations (meters)
    OBSERVER_ELEVATIONS = [
        0.0,  # Sea level
        100.0,  # Low elevation
        500.0,  # Medium elevation
        1000.0,  # High elevation
        5000.0,  # Very high (mountain)
    ]

    # Lapse rates (K/m)
    LAPSE_RATES = [
        0.003,  # Low lapse rate
        0.0065,  # Standard lapse rate
        0.010,  # High lapse rate
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("altitude", ALTITUDE_ANGLES)
    def test_refrac_extended_various_altitudes(self, altitude):
        """Compare refrac_extended at various altitude angles."""
        result_swe = swe.refrac_extended(
            altitude, 0.0, 1013.25, 15.0, 0.0065, swe.TRUE_TO_APP
        )
        result_py = ephem.refrac_extended(
            altitude, 0.0, 1013.25, 15.0, 0.0065, ephem.TRUE_TO_APP
        )

        # Compare returned altitude
        diff_alt = abs(result_swe[0] - result_py[0])
        tolerance = max(REFRACTION_TOL, abs(result_swe[1][2]) * 0.15)
        assert diff_alt < tolerance, (
            f"refrac_extended alt={altitude}: altitude diff {diff_alt:.6f} exceeds tolerance"
        )

        # Compare refraction value
        diff_refrac = abs(result_swe[1][2] - result_py[1][2])
        assert diff_refrac < tolerance, (
            f"refrac_extended alt={altitude}: refraction diff {diff_refrac:.6f}"
        )

        # Compare dip (should be 0 at sea level)
        assert result_py[1][3] == 0.0, "Dip should be 0 at sea level"

    @pytest.mark.comparison
    @pytest.mark.parametrize("pressure", PRESSURES)
    def test_refrac_extended_various_pressures(self, pressure):
        """Compare refrac_extended with different atmospheric pressures."""
        altitude = 10.0  # Use 10 degrees for visible effect

        result_swe = swe.refrac_extended(
            altitude, 0.0, pressure, 15.0, 0.0065, swe.TRUE_TO_APP
        )
        result_py = ephem.refrac_extended(
            altitude, 0.0, pressure, 15.0, 0.0065, ephem.TRUE_TO_APP
        )

        # Compare refraction values
        diff_refrac = abs(result_swe[1][2] - result_py[1][2])
        tolerance = max(REFRACTION_TOL, abs(result_swe[1][2]) * 0.15)
        assert diff_refrac < tolerance, (
            f"refrac_extended pressure={pressure}: refraction diff {diff_refrac:.6f}"
        )

        # Compare returned altitude
        diff_alt = abs(result_swe[0] - result_py[0])
        assert diff_alt < tolerance, (
            f"refrac_extended pressure={pressure}: altitude diff {diff_alt:.6f}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("temperature", TEMPERATURES)
    def test_refrac_extended_various_temperatures(self, temperature):
        """Compare refrac_extended with different temperatures."""
        altitude = 10.0  # Use 10 degrees for visible effect

        result_swe = swe.refrac_extended(
            altitude, 0.0, 1013.25, temperature, 0.0065, swe.TRUE_TO_APP
        )
        result_py = ephem.refrac_extended(
            altitude, 0.0, 1013.25, temperature, 0.0065, ephem.TRUE_TO_APP
        )

        # Compare refraction values
        diff_refrac = abs(result_swe[1][2] - result_py[1][2])
        tolerance = max(REFRACTION_TOL, abs(result_swe[1][2]) * 0.15)
        assert diff_refrac < tolerance, (
            f"refrac_extended temp={temperature}: refraction diff {diff_refrac:.6f}"
        )

        # Compare returned altitude
        diff_alt = abs(result_swe[0] - result_py[0])
        assert diff_alt < tolerance, (
            f"refrac_extended temp={temperature}: altitude diff {diff_alt:.6f}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("geoalt", OBSERVER_ELEVATIONS)
    def test_refrac_extended_various_observer_elevations(self, geoalt):
        """Compare refrac_extended dip at various observer elevations."""
        result_swe = swe.refrac_extended(
            0.0, geoalt, 1013.25, 15.0, 0.0065, swe.TRUE_TO_APP
        )
        result_py = ephem.refrac_extended(
            0.0, geoalt, 1013.25, 15.0, 0.0065, ephem.TRUE_TO_APP
        )

        # Compare dip values (most important for elevated observers)
        diff_dip = abs(result_swe[1][3] - result_py[1][3])
        assert diff_dip < DIP_TOL, (
            f"refrac_extended geoalt={geoalt}: dip diff {diff_dip:.6f} "
            f"(lib={result_py[1][3]:.6f}, swe={result_swe[1][3]:.6f})"
        )

        # Compare refraction values
        diff_refrac = abs(result_swe[1][2] - result_py[1][2])
        tolerance = max(REFRACTION_TOL, abs(result_swe[1][2]) * 0.15)
        assert diff_refrac < tolerance, (
            f"refrac_extended geoalt={geoalt}: refraction diff {diff_refrac:.6f}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("lapse_rate", LAPSE_RATES)
    def test_refrac_extended_various_lapse_rates(self, lapse_rate):
        """Compare refrac_extended with different atmospheric lapse rates."""
        # Use elevated observer where lapse rate affects dip calculation
        result_swe = swe.refrac_extended(
            0.0, 1000.0, 1013.25, 15.0, lapse_rate, swe.TRUE_TO_APP
        )
        result_py = ephem.refrac_extended(
            0.0, 1000.0, 1013.25, 15.0, lapse_rate, ephem.TRUE_TO_APP
        )

        # Compare dip values (affected by lapse rate)
        diff_dip = abs(result_swe[1][3] - result_py[1][3])
        assert diff_dip < DIP_TOL, (
            f"refrac_extended lapse={lapse_rate}: dip diff {diff_dip:.6f} "
            f"(lib={result_py[1][3]:.6f}, swe={result_swe[1][3]:.6f})"
        )

    @pytest.mark.comparison
    def test_refrac_extended_app_to_true(self):
        """Compare refrac_extended in APP_TO_TRUE mode."""
        altitudes = [5.0, 10.0, 30.0, 45.0, 60.0]

        for altitude in altitudes:
            result_swe = swe.refrac_extended(
                altitude, 0.0, 1013.25, 15.0, 0.0065, swe.APP_TO_TRUE
            )
            result_py = ephem.refrac_extended(
                altitude, 0.0, 1013.25, 15.0, 0.0065, ephem.APP_TO_TRUE
            )

            # Compare returned true altitude
            diff_alt = abs(result_swe[0] - result_py[0])
            tolerance = max(REFRACTION_TOL, abs(altitude - result_swe[0]) * 0.15)
            assert diff_alt < tolerance, (
                f"refrac_extended APP_TO_TRUE alt={altitude}: diff {diff_alt:.6f}"
            )

    @pytest.mark.comparison
    def test_refrac_extended_combined_conditions(self):
        """Test refrac_extended with realistic combined atmospheric conditions.

        Simulates real-world scenarios for astronomical observations:
        - Mountain observatory at high altitude with cold, thin air
        - Sea level tropical observation with hot, humid conditions
        """
        # Mountain observatory scenario: high altitude, cold, low pressure
        # Note: At extreme conditions, implementations may differ slightly
        # Use a larger tolerance for these edge cases
        EXTREME_DIP_TOL = 0.1  # degrees - acceptable for extreme conditions

        mountain_swe = swe.refrac_extended(
            10.0, 3000.0, 700.0, -5.0, 0.0065, swe.TRUE_TO_APP
        )
        mountain_py = ephem.refrac_extended(
            10.0, 3000.0, 700.0, -5.0, 0.0065, ephem.TRUE_TO_APP
        )

        diff_dip = abs(mountain_swe[1][3] - mountain_py[1][3])
        assert diff_dip < EXTREME_DIP_TOL, f"Mountain scenario dip diff: {diff_dip:.6f}"

        diff_refrac = abs(mountain_swe[1][2] - mountain_py[1][2])
        tolerance = max(REFRACTION_TOL, abs(mountain_swe[1][2]) * 0.15)
        assert diff_refrac < tolerance, (
            f"Mountain scenario refraction diff: {diff_refrac:.6f}"
        )

        # Tropical sea level scenario: hot, high pressure
        tropical_swe = swe.refrac_extended(
            10.0, 0.0, 1020.0, 30.0, 0.0065, swe.TRUE_TO_APP
        )
        tropical_py = ephem.refrac_extended(
            10.0, 0.0, 1020.0, 30.0, 0.0065, ephem.TRUE_TO_APP
        )

        diff_refrac = abs(tropical_swe[1][2] - tropical_py[1][2])
        tolerance = max(REFRACTION_TOL, abs(tropical_swe[1][2]) * 0.15)
        assert diff_refrac < tolerance, (
            f"Tropical scenario refraction diff: {diff_refrac:.6f}"
        )

    @pytest.mark.comparison
    def test_refrac_extended_zero_pressure(self):
        """Test that zero pressure disables refraction in both implementations."""
        result_swe = swe.refrac_extended(30.0, 0.0, 0.0, 15.0, 0.0065, swe.TRUE_TO_APP)
        result_py = ephem.refrac_extended(
            30.0, 0.0, 0.0, 15.0, 0.0065, ephem.TRUE_TO_APP
        )

        # With zero pressure, libephemeris returns input altitude exactly
        assert result_py[0] == 30.0, "Zero pressure should return input altitude"
        assert result_py[1][2] == 0.0, "Zero pressure should have zero refraction"

        # Note: pyswisseph may still apply a minimal adjustment even with zero pressure
        # Verify both implementations produce similar results
        diff = abs(result_swe[0] - result_py[0])
        assert diff < 0.01, (
            f"Zero pressure altitude diff: {diff:.6f} "
            f"(swe={result_swe[0]:.6f}, lib={result_py[0]:.6f})"
        )
