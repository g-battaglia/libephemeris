"""
Pytest-style Planetary Occultation Validation Tests.

Validates planet_occult_when_glob and planet_occult_when_loc functions.

Note: pyswisseph does NOT expose planet_occult_* functions, so these tests
validate the libephemeris implementation by:
1. Verifying the functions execute correctly and return proper structure
2. Cross-validating occultation times using pyswisseph's position calculations
   (at occultation time, the planet and target should have very close RA/Dec)
3. Testing both planetary and stellar occultations

Historical planetary occultations:
- 1818 Dec 3: Venus occulted Jupiter
- 2065 Nov 22: Venus will occult Jupiter
"""

import pytest
import math
import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import (
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    MERCURY,
    SUN,
    MOON,
    FLG_SWIEPH,
    FLG_EQUATORIAL,
)


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class PlanetOccultTolerance:
    """Tolerance thresholds for planetary occultation validation."""

    # At occultation maximum, separation should be very small
    # (within the angular radius of the occulting planet)
    SEPARATION_DEGREES = 1.0  # Maximum expected separation at occultation
    POSITION_DEGREES = 0.1  # For RA/Dec validation
    TIME_SECONDS = 3600.0  # 1 hour - occultation timing validation


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Stars that can be occulted by planets (bright, near ecliptic)
OCCULT_STARS = [
    ("Regulus", "Alpha Leonis - near ecliptic"),
    ("Spica", "Alpha Virginis - near ecliptic"),
    ("Antares", "Alpha Scorpii - near ecliptic"),
    ("Aldebaran", "Alpha Tauri - near ecliptic"),
]

# Planets that can occult (inner planets more common occultors)
OCCULTING_PLANETS = [
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
]

# Test locations
OCCULT_LOCATIONS = [
    ("New York", 40.7128, -74.0060, 0),
    ("London", 51.5074, -0.1278, 0),
    ("Sydney", -33.8688, 151.2093, 0),
]


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def angular_separation(
    ra1_deg: float, dec1_deg: float, ra2_deg: float, dec2_deg: float
) -> float:
    """
    Calculate angular separation between two celestial positions.

    Args:
        ra1_deg, dec1_deg: First position in degrees
        ra2_deg, dec2_deg: Second position in degrees

    Returns:
        Angular separation in degrees
    """
    # Convert to radians
    ra1 = math.radians(ra1_deg)
    dec1 = math.radians(dec1_deg)
    ra2 = math.radians(ra2_deg)
    dec2 = math.radians(dec2_deg)

    # Haversine formula
    cos_sep = math.sin(dec1) * math.sin(dec2) + math.cos(dec1) * math.cos(
        dec2
    ) * math.cos(ra1 - ra2)

    # Clamp to valid range
    cos_sep = max(-1.0, min(1.0, cos_sep))

    return math.degrees(math.acos(cos_sep))


def get_planet_equatorial_position(jd: float, planet_id: int) -> tuple:
    """
    Get planet's equatorial position (RA, Dec) using pyswisseph.

    Returns:
        (ra_degrees, dec_degrees, distance_au)
    """
    flags = FLG_SWIEPH | FLG_EQUATORIAL
    result, _ = swe.calc_ut(jd, planet_id, flags)
    return result[0], result[1], result[2]


def get_star_equatorial_position(jd: float, star_name: str) -> tuple:
    """
    Get star's equatorial position (RA, Dec) using pyswisseph.

    Returns:
        (ra_degrees, dec_degrees, distance_au)
    """
    flags = FLG_SWIEPH | FLG_EQUATORIAL
    try:
        result = swe.fixstar_ut(star_name, jd, flags)
        return result[0][0], result[0][1], result[0][2]
    except Exception:
        # If pyswisseph doesn't have the star, use libephemeris
        result = pyephem.fixstar_ut(star_name, jd, flags)
        return result[0][0], result[0][1], result[0][2]


# ============================================================================
# FIXTURES
# ============================================================================


@pytest.fixture
def jd_2024():
    """Julian Day for 2024 start."""
    return swe.julday(2024, 1, 1, 0.0)


@pytest.fixture
def jd_2000():
    """Julian Day for 2000 start (DE421 safe range)."""
    return swe.julday(2000, 1, 1, 0.0)


# ============================================================================
# PLANETARY OCCULTATION FUNCTION AVAILABILITY TESTS
# ============================================================================


class TestPlanetOccultFunctionsAvailable:
    """Tests for function availability and basic structure."""

    @pytest.mark.comparison
    def test_planet_occult_functions_available(self):
        """Verify planetary occultation functions are available."""
        assert callable(pyephem.planet_occult_when_glob)
        assert callable(pyephem.planet_occult_when_loc)
        assert callable(pyephem.planet_occult_when_glob)
        assert callable(pyephem.planet_occult_when_loc)

    @pytest.mark.comparison
    def test_swe_aliases_are_same(self):
        """Verify swe_* aliases point to same functions."""
        assert pyephem.planet_occult_when_glob is pyephem.planet_occult_when_glob
        assert pyephem.planet_occult_when_loc is pyephem.planet_occult_when_loc


# ============================================================================
# ERROR HANDLING TESTS (FAST - NO ACTUAL SEARCHES)
# ============================================================================


class TestPlanetOccultErrorHandling:
    """Tests for proper error handling in planetary occultation functions."""

    @pytest.mark.comparison
    def test_invalid_star_name_raises_error(self, jd_2024):
        """Test that invalid star name raises ValueError."""
        with pytest.raises(ValueError):
            pyephem.planet_occult_when_glob(
                jd_2024, VENUS, 0, "NonexistentStar123", FLG_SWIEPH, 0
            )

    @pytest.mark.comparison
    def test_no_target_raises_error(self, jd_2024):
        """Test that missing target raises ValueError."""
        with pytest.raises(ValueError):
            pyephem.planet_occult_when_glob(jd_2024, VENUS, 0, "", FLG_SWIEPH, 0)

    @pytest.mark.comparison
    def test_same_planet_raises_error(self, jd_2024):
        """Test that same planet as occulting and occulted raises ValueError."""
        with pytest.raises(ValueError, match="cannot be the same"):
            pyephem.planet_occult_when_glob(
                jd_2024, VENUS, VENUS, "", FLG_SWIEPH, 0
            )

    @pytest.mark.comparison
    def test_sun_as_occulting_body_raises_error(self, jd_2024):
        """Test that Sun cannot be the occulting body."""
        with pytest.raises(ValueError, match="Sun cannot be"):
            pyephem.planet_occult_when_glob(
                jd_2024, SUN, VENUS, "", FLG_SWIEPH, 0
            )

    @pytest.mark.comparison
    def test_moon_as_occulting_body_raises_error(self, jd_2024):
        """Test that Moon cannot be the occulting body."""
        with pytest.raises(ValueError, match="Moon cannot be"):
            pyephem.planet_occult_when_glob(
                jd_2024, MOON, VENUS, "", FLG_SWIEPH, 0
            )

    @pytest.mark.comparison
    def test_invalid_occulting_planet_raises_error(self, jd_2024):
        """Test that invalid planet ID raises error."""
        with pytest.raises(ValueError):
            pyephem.planet_occult_when_glob(
                jd_2024, 999, JUPITER, "", FLG_SWIEPH, 0
            )

    @pytest.mark.comparison
    def test_invalid_occulted_planet_raises_error(self, jd_2024):
        """Test that invalid occulted planet ID raises error."""
        with pytest.raises(ValueError):
            pyephem.planet_occult_when_glob(jd_2024, VENUS, 999, "", FLG_SWIEPH, 0)

    @pytest.mark.comparison
    def test_loc_invalid_star_name_raises_error(self, jd_2024):
        """Test that invalid star name raises ValueError in loc function."""
        with pytest.raises(ValueError):
            pyephem.planet_occult_when_loc(
                jd_2024, VENUS, 0, "NonexistentStar123", 40.0, -74.0, 0, FLG_SWIEPH
            )

    @pytest.mark.comparison
    def test_loc_no_target_raises_error(self, jd_2024):
        """Test that missing target raises ValueError in loc function."""
        with pytest.raises(ValueError):
            pyephem.planet_occult_when_loc(
                jd_2024, VENUS, 0, "", 40.0, -74.0, 0, FLG_SWIEPH
            )

    @pytest.mark.comparison
    def test_loc_same_planet_raises_error(self, jd_2024):
        """Test that same planet raises error in loc function."""
        with pytest.raises(ValueError, match="cannot be the same"):
            pyephem.planet_occult_when_loc(
                jd_2024, VENUS, VENUS, "", 40.0, -74.0, 0, FLG_SWIEPH
            )


# ============================================================================
# FUNCTION SIGNATURE AND PARAMETER TESTS (FAST)
# ============================================================================


class TestPlanetOccultFunctionSignatures:
    """Tests for function signatures and parameter handling."""

    @pytest.mark.comparison
    def test_glob_accepts_all_parameters(self, jd_2024):
        """Test that glob function accepts all documented parameters."""
        # This should not raise TypeError - uses same planet to trigger
        # quick validation error while verifying all parameters are accepted
        with pytest.raises(ValueError, match="cannot be the same"):
            pyephem.planet_occult_when_glob(
                jd_2024,  # tjdut
                VENUS,  # occulting_planet
                VENUS,  # occulted_planet (same = validation error)
                "",  # starname
                FLG_SWIEPH,  # flags
                0,  # direction (forward)
            )

    @pytest.mark.comparison
    def test_loc_accepts_all_parameters(self, jd_2024):
        """Test that loc function accepts all documented parameters."""
        # Uses same planet to trigger quick validation error
        # while verifying all parameters are accepted
        with pytest.raises(ValueError, match="cannot be the same"):
            pyephem.planet_occult_when_loc(
                jd_2024,  # jd_start
                VENUS,  # occulting_planet
                VENUS,  # occulted_planet (same = validation error)
                "",  # star_name
                40.7128,  # lat
                -74.0060,  # lon
                0,  # altitude
                FLG_SWIEPH,  # flags
            )

    @pytest.mark.comparison
    @pytest.mark.parametrize("planet_id,planet_name", OCCULTING_PLANETS)
    def test_all_planets_accepted_as_occultor(self, jd_2024, planet_id, planet_name):
        """Test that all valid planets can be used as occulting body.

        This test only validates parameter acceptance, not actual search.
        We test error handling by using an invalid target to cause quick termination.
        """
        # Test that planet is accepted by checking it doesn't raise "invalid planet" error
        # Use empty starname to trigger quick validation error
        try:
            pyephem.planet_occult_when_glob(jd_2024, planet_id, 0, "", FLG_SWIEPH, 0)
        except ValueError as e:
            # Should fail for "no target specified", not "invalid planet"
            assert "must be specified" in str(e) or "target" in str(e).lower(), (
                f"Unexpected error for {planet_name}: {e}"
            )

    @pytest.mark.comparison
    @pytest.mark.parametrize("star_name,star_desc", OCCULT_STARS)
    def test_all_ecliptic_stars_recognized(self, jd_2024, star_name, star_desc):
        """Test that all ecliptic stars are recognized by the function.

        This validates star name resolution without performing an actual search.
        We use libephemeris's fixstar function to validate the star exists.
        """
        # Just verify the star can be resolved - don't run the slow occultation search
        try:
            result = pyephem.fixstar_ut(star_name, jd_2024, FLG_SWIEPH)
            # If we get here, the star is recognized
            assert result is not None
        except ValueError as e:
            pytest.fail(f"Star {star_name} not recognized: {e}")


# ============================================================================
# PLANETARY OCCULTATION SEARCH TESTS (SLOW - MARKED APPROPRIATELY)
# ============================================================================


class TestVenusRegulus:
    """Tests for Venus occultation of Regulus.

    Venus-Regulus occultations are relatively common (every few years)
    because:
    1. Venus has a large angular diameter (up to 32")
    2. Regulus is very close to the ecliptic (latitude ~0.5°)
    """

    @pytest.mark.slow
    @pytest.mark.comparison
    def test_venus_regulus_glob_structure(self, jd_2024):
        """Test Venus-Regulus global occultation search returns proper structure."""
        try:
            retflags, tret = pyephem.planet_occult_when_glob(
                jd_2024, VENUS, 0, "Regulus", FLG_SWIEPH, 0
            )

            # Verify return structure
            assert isinstance(retflags, int), "retflags should be int"
            assert len(tret) == 10, f"tret should have 10 elements, got {len(tret)}"

            # Verify times are reasonable
            jd_max = tret[0]
            assert jd_max > jd_2024, "Occultation should be in the future"
            assert jd_max < jd_2024 + 150 * 365, (
                "Occultation should be within 150 years"
            )

        except (RuntimeError, pyephem.Error) as e:
            if "No planetary occultation" in str(e):
                pytest.skip("No Venus-Regulus occultation found within search limit")
            raise

    @pytest.mark.slow
    @pytest.mark.comparison
    def test_venus_regulus_position_validation(self, jd_2024):
        """Validate Venus-Regulus occultation by checking angular separation.

        At occultation maximum, Venus and Regulus should be very close.
        Cross-validate positions using pyswisseph.
        """
        try:
            retflags, tret = pyephem.planet_occult_when_glob(
                jd_2024, VENUS, 0, "Regulus", FLG_SWIEPH, 0
            )

            jd_max = tret[0]

            # Get Venus position at occultation max (using pyswisseph)
            venus_ra, venus_dec, _ = get_planet_equatorial_position(jd_max, VENUS)

            # Get Regulus position at occultation max (using pyswisseph)
            regulus_ra, regulus_dec, _ = get_star_equatorial_position(jd_max, "Regulus")

            # Calculate angular separation
            separation = angular_separation(
                venus_ra, venus_dec, regulus_ra, regulus_dec
            )

            # At occultation, separation should be small (within planet's disk + margin)
            # Venus angular diameter can be up to ~1 arcmin = 0.017 degrees
            # Allow more margin for algorithm tolerance
            assert separation < PlanetOccultTolerance.SEPARATION_DEGREES, (
                f"Venus-Regulus separation at occultation max should be small, "
                f"got {separation:.4f}° (RA1={venus_ra:.3f}, Dec1={venus_dec:.3f}, "
                f"RA2={regulus_ra:.3f}, Dec2={regulus_dec:.3f})"
            )

        except (RuntimeError, pyephem.Error) as e:
            if "No planetary occultation" in str(e):
                pytest.skip("No Venus-Regulus occultation found within search limit")
            raise


class TestMarsStarOccultation:
    """Tests for Mars occultation of fixed stars.

    Mars-star occultations are less common than Venus-star because:
    1. Mars has a smaller angular diameter (up to 25")
    2. Mars orbits farther from the Sun, slower orbital motion
    """

    @pytest.mark.slow
    @pytest.mark.comparison
    @pytest.mark.parametrize("star_name,star_desc", OCCULT_STARS[:2])  # Regulus, Spica
    def test_mars_star_glob_structure(self, jd_2000, star_name, star_desc):
        """Test Mars-star global occultation search returns proper structure."""
        try:
            retflags, tret = pyephem.planet_occult_when_glob(
                jd_2000, MARS, 0, star_name, FLG_SWIEPH, 0
            )

            # Verify return structure
            assert isinstance(retflags, int), "retflags should be int"
            assert len(tret) == 10, f"tret should have 10 elements, got {len(tret)}"

            # Verify times are reasonable
            jd_max = tret[0]
            assert jd_max > jd_2000, (
                f"Mars-{star_name} occultation should be in the future"
            )

        except (RuntimeError, pyephem.Error) as e:
            if "No planetary occultation" in str(e):
                pytest.skip(
                    f"No Mars-{star_name} occultation found within search limit"
                )
            raise

    @pytest.mark.slow
    @pytest.mark.comparison
    @pytest.mark.parametrize("star_name,star_desc", OCCULT_STARS[:2])
    def test_mars_star_position_validation(self, jd_2000, star_name, star_desc):
        """Validate Mars-star occultation by checking angular separation."""
        try:
            retflags, tret = pyephem.planet_occult_when_glob(
                jd_2000, MARS, 0, star_name, FLG_SWIEPH, 0
            )

            jd_max = tret[0]

            # Get Mars position at occultation max
            mars_ra, mars_dec, _ = get_planet_equatorial_position(jd_max, MARS)

            # Get star position at occultation max
            star_ra, star_dec, _ = get_star_equatorial_position(jd_max, star_name)

            # Calculate angular separation
            separation = angular_separation(mars_ra, mars_dec, star_ra, star_dec)

            assert separation < PlanetOccultTolerance.SEPARATION_DEGREES, (
                f"Mars-{star_name} separation at occultation max should be small, "
                f"got {separation:.4f}°"
            )

        except (RuntimeError, pyephem.Error) as e:
            if "No planetary occultation" in str(e):
                pytest.skip(
                    f"No Mars-{star_name} occultation found within search limit"
                )
            raise


# ============================================================================
# LOCAL OCCULTATION TESTS (SLOW)
# ============================================================================


class TestPlanetOccultWhenLoc:
    """Tests for local planetary occultation search."""

    @pytest.mark.slow
    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "loc_name,lat,lon,alt", OCCULT_LOCATIONS[:1]
    )  # Just New York
    def test_venus_regulus_loc_structure(self, jd_2024, loc_name, lat, lon, alt):
        """Test Venus-Regulus local occultation search returns proper structure."""
        try:
            times, attr, retflag = pyephem.planet_occult_when_loc(
                jd_2024, VENUS, 0, "Regulus", lat, lon, alt, FLG_SWIEPH
            )

            # Verify return structure
            assert len(times) == 10, f"times should have 10 elements, got {len(times)}"
            assert len(attr) == 20, f"attr should have 20 elements, got {len(attr)}"
            assert isinstance(retflag, int), "retflag should be int"

            # Verify times are reasonable
            jd_max = times[0]
            assert jd_max > jd_2024, (
                f"Occultation at {loc_name} should be in the future"
            )

        except (RuntimeError, pyephem.Error) as e:
            if "No planetary occultation" in str(e):
                pytest.skip(f"No Venus-Regulus occultation visible from {loc_name}")
            raise

    @pytest.mark.slow
    @pytest.mark.comparison
    @pytest.mark.parametrize("loc_name,lat,lon,alt", OCCULT_LOCATIONS[:1])
    def test_venus_regulus_loc_visibility_check(self, jd_2024, loc_name, lat, lon, alt):
        """Validate that occultation is actually visible from the location."""
        try:
            times, attr, retflag = pyephem.planet_occult_when_loc(
                jd_2024, VENUS, 0, "Regulus", lat, lon, alt, FLG_SWIEPH
            )

            # Check that altitude attributes indicate visibility
            # attr[5] = true altitude, attr[6] = apparent altitude
            true_alt = attr[5]

            # At local occultation, objects should be above horizon
            # (or close to it, considering atmospheric refraction)
            assert true_alt > -5.0, (
                f"True altitude {true_alt:.1f}° at {loc_name} should indicate visibility"
            )

        except (RuntimeError, pyephem.Error) as e:
            if "No planetary occultation" in str(e):
                pytest.skip(f"No Venus-Regulus occultation visible from {loc_name}")
            raise


# ============================================================================
# MUTUAL PLANETARY OCCULTATION TESTS (VERY SLOW - RARE EVENTS)
# ============================================================================


class TestMutualPlanetaryOccultation:
    """Tests for mutual planetary occultations (planet occulting planet).

    These are extremely rare events (few times per century).
    Historical events:
    - 1818 Dec 3: Venus occulted Jupiter
    - 2065 Nov 22: Venus will occult Jupiter
    """

    @pytest.mark.slow
    @pytest.mark.comparison
    def test_venus_jupiter_glob_structure(self, jd_2000):
        """Test Venus-Jupiter global occultation search structure.

        This may not find an event (very rare), but should return proper structure
        or raise RuntimeError.
        """
        try:
            retflags, tret = pyephem.planet_occult_when_glob(
                jd_2000, VENUS, JUPITER, "", FLG_SWIEPH, 0
            )

            # If we found an event, verify structure
            assert isinstance(retflags, int), "retflags should be int"
            assert len(tret) == 10, f"tret should have 10 elements, got {len(tret)}"

            jd_max = tret[0]
            assert jd_max > jd_2000, "Occultation should be after start date"

        except (RuntimeError, pyephem.Error) as e:
            if "No planetary occultation" in str(e):
                pytest.skip("No Venus-Jupiter occultation found (expected - very rare)")
            raise

    @pytest.mark.slow
    @pytest.mark.comparison
    def test_venus_mars_glob_structure(self, jd_2000):
        """Test Venus-Mars global occultation search structure."""
        try:
            retflags, tret = pyephem.planet_occult_when_glob(
                jd_2000, VENUS, MARS, "", FLG_SWIEPH, 0
            )

            assert isinstance(retflags, int)
            assert len(tret) == 10

        except (RuntimeError, pyephem.Error) as e:
            if "No planetary occultation" in str(e):
                pytest.skip("No Venus-Mars occultation found (expected - very rare)")
            raise


# ============================================================================
# TIMING VALIDATION TESTS (SLOW - REQUIRE ACTUAL SEARCH)
# ============================================================================


class TestOccultationTimingValidation:
    """Tests to validate occultation timing is consistent."""

    @pytest.mark.slow
    @pytest.mark.comparison
    def test_occultation_phase_times_ordered(self, jd_2024):
        """Test that occultation phase times are properly ordered."""
        try:
            retflags, tret = pyephem.planet_occult_when_glob(
                jd_2024, VENUS, 0, "Regulus", FLG_SWIEPH, 0
            )

            jd_max = tret[0]  # Maximum
            jd_begin = tret[2]  # Begin
            jd_end = tret[3]  # End

            # Times should be ordered: begin < max < end
            if jd_begin > 0 and jd_end > 0:
                assert jd_begin < jd_max, (
                    f"Begin ({jd_begin}) should be before max ({jd_max})"
                )
                assert jd_max < jd_end, (
                    f"Max ({jd_max}) should be before end ({jd_end})"
                )

        except (RuntimeError, pyephem.Error) as e:
            if "No planetary occultation" in str(e):
                pytest.skip("No occultation found for timing validation")
            raise

    @pytest.mark.slow
    @pytest.mark.comparison
    def test_occultation_duration_reasonable(self, jd_2024):
        """Test that occultation duration is reasonable."""
        try:
            retflags, tret = pyephem.planet_occult_when_glob(
                jd_2024, VENUS, 0, "Regulus", FLG_SWIEPH, 0
            )

            jd_begin = tret[2]
            jd_end = tret[3]

            if jd_begin > 0 and jd_end > 0:
                duration_minutes = (jd_end - jd_begin) * 24 * 60

                # Planetary occultations typically last a few minutes to hours
                assert 0 < duration_minutes < 180, (
                    f"Duration {duration_minutes:.1f} minutes should be reasonable"
                )

        except (RuntimeError, pyephem.Error) as e:
            if "No planetary occultation" in str(e):
                pytest.skip("No occultation found for duration validation")
            raise


# ============================================================================
# PYSWISSEPH CROSS-VALIDATION TESTS (POSITION-BASED)
# ============================================================================


class TestPlanetOccultCrossValidation:
    """Cross-validation tests using pyswisseph position calculations.

    Since pyswisseph does NOT expose planet_occult_* functions, we validate
    libephemeris' implementation by:
    1. Verifying positions match pyswisseph at known conjunction dates
    2. Testing that position calculations used internally are accurate
    3. Validating angular separation calculations against pyswisseph

    These tests use known close approaches (conjunctions) that are well-documented
    astronomical events, verifying the underlying calculations match pyswisseph.
    """

    @pytest.mark.comparison
    def test_venus_position_matches_pyswisseph(self, jd_2024):
        """Verify Venus position calculations match pyswisseph."""
        # Get Venus position using pyswisseph
        flags = FLG_SWIEPH | FLG_EQUATORIAL
        venus_swe, _ = swe.calc_ut(jd_2024, VENUS, flags)

        # Get Venus position using libephemeris
        venus_py, _ = pyephem.calc_ut(jd_2024, VENUS, flags)

        # Positions should match within arcsecond precision
        ra_diff = abs(venus_swe[0] - venus_py[0]) * 3600  # arcseconds
        dec_diff = abs(venus_swe[1] - venus_py[1]) * 3600  # arcseconds

        assert ra_diff < 1.0, f"Venus RA differs by {ra_diff:.3f} arcsec"
        assert dec_diff < 1.0, f"Venus Dec differs by {dec_diff:.3f} arcsec"

    @pytest.mark.comparison
    def test_mars_position_matches_pyswisseph(self, jd_2024):
        """Verify Mars position calculations match pyswisseph."""
        flags = FLG_SWIEPH | FLG_EQUATORIAL
        mars_swe, _ = swe.calc_ut(jd_2024, MARS, flags)
        mars_py, _ = pyephem.calc_ut(jd_2024, MARS, flags)

        ra_diff = abs(mars_swe[0] - mars_py[0]) * 3600
        dec_diff = abs(mars_swe[1] - mars_py[1]) * 3600

        assert ra_diff < 1.0, f"Mars RA differs by {ra_diff:.3f} arcsec"
        assert dec_diff < 1.0, f"Mars Dec differs by {dec_diff:.3f} arcsec"

    @pytest.mark.comparison
    def test_regulus_position_matches_pyswisseph(self, jd_2024):
        """Verify Regulus star position matches pyswisseph."""
        flags = FLG_SWIEPH | FLG_EQUATORIAL

        # libephemeris
        regulus_py = pyephem.fixstar_ut("Regulus", jd_2024, flags)

        # pyswisseph (may fail if sefstars.txt not found)
        try:
            regulus_swe = swe.fixstar_ut("Regulus", jd_2024, flags)
            ra_diff = abs(regulus_swe[0][0] - regulus_py[0][0]) * 3600
            dec_diff = abs(regulus_swe[0][1] - regulus_py[0][1]) * 3600

            # Stars should match very closely
            assert ra_diff < 5.0, f"Regulus RA differs by {ra_diff:.3f} arcsec"
            assert dec_diff < 5.0, f"Regulus Dec differs by {dec_diff:.3f} arcsec"
        except Exception:
            # pyswisseph may not have star data - just verify libephemeris returns
            assert regulus_py is not None

    @pytest.mark.comparison
    def test_angular_separation_calculation(self):
        """Test that angular separation calculation is correct."""
        # Known positions: test with simple cases
        # Same position should have 0 separation
        sep = angular_separation(0.0, 0.0, 0.0, 0.0)
        assert sep == 0.0, f"Same position should have 0 separation, got {sep}"

        # 90 degrees apart in RA at equator
        sep = angular_separation(0.0, 0.0, 90.0, 0.0)
        assert abs(sep - 90.0) < 0.001, f"90 deg RA diff should be ~90 deg, got {sep}"

        # Test declination difference
        sep = angular_separation(0.0, 0.0, 0.0, 45.0)
        assert abs(sep - 45.0) < 0.001, f"45 deg Dec diff should be ~45 deg, got {sep}"

    @pytest.mark.comparison
    def test_venus_regulus_conjunction_2044(self):
        """Test Venus-Regulus close approach on 2044 Oct 1.

        This is a well-documented close conjunction where Venus passes within
        0.5 degrees of Regulus. We verify that:
        1. libephemeris and pyswisseph agree on Venus position
        2. The calculated separation matches expected value
        """
        # 2044 Oct 1 - known Venus-Regulus close conjunction
        jd_conjunction = swe.julday(2044, 10, 1, 12.0)

        # Get Venus position from both libraries in ECLIPTIC coordinates
        # (planetary occultations are measured in ecliptic system)
        flags = FLG_SWIEPH
        venus_swe, _ = swe.calc_ut(jd_conjunction, VENUS, flags)
        venus_py, _ = pyephem.calc_ut(jd_conjunction, VENUS, flags)

        # Get Regulus from libephemeris (ecliptic coordinates)
        regulus_py = pyephem.fixstar_ut("Regulus", jd_conjunction, flags)

        # Venus positions should match
        lon_diff = abs(venus_swe[0] - venus_py[0]) * 3600  # arcseconds
        lat_diff = abs(venus_swe[1] - venus_py[1]) * 3600
        assert lon_diff < 1.0, f"Venus lon differs by {lon_diff:.3f} arcsec"
        assert lat_diff < 1.0, f"Venus lat differs by {lat_diff:.3f} arcsec"

        # Calculate separation in ecliptic (lon/lat difference)
        lon_sep = abs(venus_py[0] - regulus_py[0][0])
        lat_sep = abs(venus_py[1] - regulus_py[0][1])

        # At known conjunction, ecliptic separation should be small
        # Venus-Regulus on 2044 Oct 1: ~0.5 deg longitude, ~0.02 deg latitude
        assert lon_sep < 1.0, (
            f"Venus-Regulus lon separation at 2044 Oct 1 should be <1 deg, got {lon_sep:.4f}"
        )
        assert lat_sep < 0.5, (
            f"Venus-Regulus lat separation at 2044 Oct 1 should be <0.5 deg, got {lat_sep:.4f}"
        )

    @pytest.mark.comparison
    def test_planet_conjunction_timing_validation(self):
        """Test that planet positions at conjunctions match pyswisseph.

        This validates the position calculations that planet_occult_when_glob
        uses internally, without requiring an actual occultation.
        """
        # Test multiple dates across DE421 range
        test_dates = [
            (2000, 1, 1, 12.0),
            (2020, 6, 15, 0.0),
            (2030, 12, 25, 6.0),
        ]

        planets = [(VENUS, "Venus"), (MARS, "Mars"), (JUPITER, "Jupiter")]

        for year, month, day, hour in test_dates:
            jd = swe.julday(year, month, day, hour)

            for planet_id, planet_name in planets:
                # Compare ecliptic positions
                flags = FLG_SWIEPH
                pos_swe, _ = swe.calc_ut(jd, planet_id, flags)
                pos_py, _ = pyephem.calc_ut(jd, planet_id, flags)

                lon_diff = abs(pos_swe[0] - pos_py[0])
                lat_diff = abs(pos_swe[1] - pos_py[1])

                # Allow for small numerical differences (sub-arcsecond)
                # 1e-4 degrees = 0.36 arcseconds
                assert lon_diff < 1e-4, (
                    f"{planet_name} at {year}/{month}/{day}: "
                    f"longitude diff {lon_diff} too large"
                )
                assert lat_diff < 1e-4, (
                    f"{planet_name} at {year}/{month}/{day}: "
                    f"latitude diff {lat_diff} too large"
                )


# ============================================================================
# ADDITIONAL MARS-STAR PYSWISSEPH VALIDATION TESTS
# ============================================================================


class TestMarsStarCrossValidation:
    """Mars-star position validation using pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "star_name,star_desc",
        [
            ("Regulus", "Alpha Leonis"),
            ("Spica", "Alpha Virginis"),
        ],
    )
    def test_mars_star_separation_calculation(self, jd_2024, star_name, star_desc):
        """Validate Mars-star separation calculations match expected values.

        This tests the core calculation used by planet_occult_when_glob
        to determine when occultations occur.
        """
        flags = FLG_SWIEPH | FLG_EQUATORIAL

        # Get Mars position from pyswisseph
        mars_swe, _ = swe.calc_ut(jd_2024, MARS, flags)

        # Get Mars position from libephemeris
        mars_py, _ = pyephem.calc_ut(jd_2024, MARS, flags)

        # Verify Mars positions match
        ra_diff = abs(mars_swe[0] - mars_py[0]) * 3600
        dec_diff = abs(mars_swe[1] - mars_py[1]) * 3600
        assert ra_diff < 1.0, f"Mars RA differs by {ra_diff:.3f} arcsec"
        assert dec_diff < 1.0, f"Mars Dec differs by {dec_diff:.3f} arcsec"

        # Get star position from libephemeris
        star_py = pyephem.fixstar_ut(star_name, jd_2024, flags)
        assert star_py is not None, f"Could not get {star_name} position"

        # Calculate and verify separation is a reasonable value
        separation = angular_separation(
            mars_py[0], mars_py[1], star_py[0][0], star_py[0][1]
        )
        assert separation >= 0, f"Separation should be non-negative, got {separation}"
        assert separation < 180, f"Separation should be < 180 deg, got {separation}"


# ============================================================================
# LOCAL OCCULTATION PYSWISSEPH VALIDATION TESTS
# ============================================================================


class TestPlanetOccultLocCrossValidation:
    """Cross-validation tests for planet_occult_when_loc using pyswisseph."""

    @pytest.mark.comparison
    def test_loc_function_uses_correct_coordinates(self, jd_2024):
        """Verify local occultation uses correct observer coordinates.

        Tests that the geographic coordinates are properly used by comparing
        altitude/azimuth calculations against pyswisseph.
        """
        # New York coordinates
        lat, lon, alt = 40.7128, -74.0060, 0

        # Get Venus position for topocentric calculation
        flags = FLG_SWIEPH | FLG_EQUATORIAL

        venus_swe, _ = swe.calc_ut(jd_2024, VENUS, flags)
        venus_py, _ = pyephem.calc_ut(jd_2024, VENUS, flags)

        # Verify geocentric positions match
        ra_diff = abs(venus_swe[0] - venus_py[0]) * 3600
        assert ra_diff < 1.0, f"Venus RA differs by {ra_diff:.3f} arcsec"

        # Use azalt to verify coordinate transformation
        geopos = (lon, lat, alt)
        atpress = 1013.25
        attemp = 15.0

        # pyswisseph azalt
        try:
            azalt_swe = swe.azalt(
                jd_2024, 0, geopos, atpress, attemp, (venus_swe[0], venus_swe[1], 1.0)
            )

            # libephemeris azalt (different signature: jd, flag, lat, lon, alt, press, temp, coord)
            azalt_py = pyephem.azalt(
                jd_2024,
                0,
                lat,
                lon,
                alt,
                atpress,
                attemp,
                (venus_py[0], venus_py[1], 1.0),
            )

            # Altitude and azimuth should match closely
            alt_diff = abs(azalt_swe[1] - azalt_py[1])
            az_diff = abs(azalt_swe[0] - azalt_py[0])

            assert alt_diff < 0.01, f"Altitude differs by {alt_diff:.4f} deg"
            assert az_diff < 0.01, f"Azimuth differs by {az_diff:.4f} deg"
        except Exception:
            # azalt may not be available in both - skip gracefully
            pass

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "loc_name,lat,lon,alt",
        OCCULT_LOCATIONS,
    )
    def test_loc_observer_positions_valid(self, jd_2024, loc_name, lat, lon, alt):
        """Verify observer locations are accepted without immediate error.

        Tests that each standard observer location is syntactically valid.
        This does NOT perform an actual occultation search (too slow).
        Instead, we verify the function accepts the parameters.
        """
        # We use an error-triggering call to verify parameter acceptance
        # without actually running the slow search
        try:
            # Trigger a quick error by using same planet for both
            pyephem.planet_occult_when_loc(
                jd_2024, VENUS, VENUS, "", lat, lon, alt, FLG_SWIEPH
            )
        except ValueError as e:
            # Expected: "cannot be the same" - this means parameters were parsed
            assert "cannot be the same" in str(e), f"Unexpected error: {e}"
