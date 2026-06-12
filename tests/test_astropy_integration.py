"""
Tests for astropy integration evaluation module.

This module tests the astropy_integration module which evaluates how
astropy.coordinates and astropy.time could supplement libephemeris calculations.

The tests verify:
1. Basic functionality (check_astropy_available, evaluate_integration_potential)
2. Time scale conversions and Delta T comparisons
3. Coordinate transformations (ICRS, Galactic, Alt/Az)
4. Geodetic utilities
5. Barycentric corrections
"""

import pytest

from libephemeris.astropy_integration import (
    check_astropy_available,
    evaluate_integration_potential,
)

# Conditionally import the rest only if astropy is available
pytestmark = pytest.mark.skipif(
    not check_astropy_available(),
    reason="Astropy is not installed. Install with: pip install astropy",
)


# =============================================================================
# BASIC FUNCTIONALITY TESTS
# =============================================================================


def test_check_astropy_available():
    """Test that check_astropy_available returns a boolean."""
    result = check_astropy_available()
    assert isinstance(result, bool)


def test_evaluate_integration_potential():
    """Test the comprehensive evaluation function."""
    from libephemeris.astropy_integration import evaluate_integration_potential

    result = evaluate_integration_potential()

    # Check structure
    assert isinstance(result, dict)
    assert "astropy_available" in result
    assert "time_features" in result
    assert "coordinate_features" in result
    assert "recommended_integrations" in result
    assert "compatibility_notes" in result
    assert "implementation_approach" in result

    # Check types
    assert isinstance(result["astropy_available"], bool)
    assert isinstance(result["time_features"], list)
    assert len(result["time_features"]) > 0
    assert isinstance(result["coordinate_features"], list)
    assert len(result["coordinate_features"]) > 0


# =============================================================================
# TIME SCALE TESTS
# =============================================================================


class TestTimeConversions:
    """Tests for time scale conversions and comparisons."""

    def test_compare_time_conversions_j2000(self):
        """Test time conversion comparison at J2000.0 epoch."""
        from libephemeris.astropy_integration import compare_time_conversions

        jd_j2000 = 2451545.0  # Jan 1, 2000 12:00 TT
        result = compare_time_conversions(jd_j2000)

        # Check structure
        assert "jd_utc" in result
        assert "skyfield_delta_t_seconds" in result
        assert "astropy_delta_t_seconds" in result
        assert "delta_t_difference_ms" in result
        assert "skyfield_tt" in result
        assert "astropy_tt" in result
        assert "tt_difference_seconds" in result
        assert "astropy_tai" in result
        assert "astropy_tdb" in result
        assert "astropy_tcg" in result
        assert "astropy_tcb" in result

        # Delta T at J2000 should be around 63-64 seconds
        assert 60 < result["skyfield_delta_t_seconds"] < 70
        assert 60 < result["astropy_delta_t_seconds"] < 70

        # Difference between implementations can be up to ~1 second due to
        # different Delta T models (Skyfield uses Stephenson et al. 2016,
        # Astropy uses IERS data). This is documented and expected.
        assert abs(result["delta_t_difference_ms"]) < 1000  # Within 1 second

    def test_compare_time_conversions_modern(self):
        """Test time conversion comparison for a modern date."""
        from libephemeris.astropy_integration import compare_time_conversions

        jd_2020 = 2458849.5  # Jan 1, 2020
        result = compare_time_conversions(jd_2020)

        # Delta T in 2020 should be around 69 seconds
        assert 65 < result["skyfield_delta_t_seconds"] < 75
        assert 65 < result["astropy_delta_t_seconds"] < 75

    def test_get_extended_time_scales(self):
        """Test extended time scale conversion."""
        from libephemeris.astropy_integration import get_extended_time_scales

        jd = 2451545.0
        scales = get_extended_time_scales(jd)

        # Check all scales are present
        expected_scales = ["utc", "ut1", "tai", "tt", "tdb", "tcg", "tcb", "gps"]
        for scale in expected_scales:
            assert scale in scales
            assert isinstance(scales[scale], float)
            # All JD values should be close to the input
            assert abs(scales[scale] - jd) < 0.01  # Within ~14 minutes

        # TAI should be ahead of UTC by leap seconds (~32 sec in 2000)
        tai_utc_diff_seconds = (scales["tai"] - scales["utc"]) * 86400
        assert 30 < tai_utc_diff_seconds < 40

        # TT should be ahead of TAI by exactly 32.184 seconds
        tt_tai_diff_seconds = (scales["tt"] - scales["tai"]) * 86400
        assert abs(tt_tai_diff_seconds - 32.184) < 0.001

    def test_parse_time_string_iso(self):
        """Test ISO time string parsing."""
        from libephemeris.astropy_integration import parse_time_string

        # Parse J2000.0 epoch
        jd = parse_time_string("2000-01-01T12:00:00")

        # Should be very close to J2000.0 (2451545.0)
        # Note: There's a small difference because J2000 is defined in TT, not UTC
        assert abs(jd - 2451545.0) < 0.001

    def test_parse_time_string_various_formats(self):
        """Test parsing various time string formats."""
        from libephemeris.astropy_integration import parse_time_string

        # Test different format strings
        jd1 = parse_time_string("2000-01-01 12:00:00")
        jd2 = parse_time_string("2000-01-01T12:00:00.000")

        # Both should give approximately the same JD
        assert abs(jd1 - jd2) < 0.0001


# =============================================================================
# COORDINATE TRANSFORMATION TESTS
# =============================================================================


class TestCoordinateTransforms:
    """Tests for coordinate transformation comparisons."""

    def test_compare_coordinate_transforms(self):
        """Test Alt/Az comparison between Skyfield and Astropy."""
        from libephemeris.astropy_integration import compare_coordinate_transforms

        # Test with Sirius coordinates at Rome
        result = compare_coordinate_transforms(
            ra_deg=101.2875,  # Sirius RA
            dec_deg=-16.7161,  # Sirius Dec
            jd_utc=2451545.0,
            lon_deg=12.5,  # Rome longitude
            lat_deg=41.9,  # Rome latitude
        )

        # Check structure
        assert "icrs_ra" in result
        assert "icrs_dec" in result
        assert "skyfield_alt" in result
        assert "skyfield_az" in result
        assert "astropy_alt" in result
        assert "astropy_az" in result
        assert "alt_difference_arcsec" in result
        assert "az_difference_arcsec" in result
        assert "astropy_galactic_l" in result
        assert "astropy_galactic_b" in result
        assert "astropy_ecliptic_lon" in result
        assert "astropy_ecliptic_lat" in result

        # Altitude/Azimuth should be similar between implementations
        # Allowing for small differences due to different nutation models
        assert abs(result["alt_difference_arcsec"]) < 60  # Within 1 arcminute
        assert abs(result["az_difference_arcsec"]) < 60  # Within 1 arcminute

        # Galactic coordinates for Sirius (approximately l=227, b=-9)
        assert 220 < result["astropy_galactic_l"] < 235
        assert -15 < result["astropy_galactic_b"] < -5

    def test_icrs_to_galactic_center(self):
        """Test ICRS to Galactic conversion for Galactic center."""
        from libephemeris.astropy_integration import icrs_to_galactic

        # Galactic center in ICRS (approximately)
        # RA = 266.417°, Dec = -29.008° should give l~0, b~0
        gal_l, gal_b = icrs_to_galactic(266.417, -29.008)

        # Should be very close to Galactic center
        # Allow for small offset due to precise definition
        assert abs(gal_l) < 1 or abs(gal_l - 360) < 1  # l near 0° or 360°
        assert abs(gal_b) < 2  # b near 0°

    def test_icrs_to_galactic_pole(self):
        """Test ICRS to Galactic conversion for North Galactic Pole."""
        from libephemeris.astropy_integration import icrs_to_galactic

        # North Galactic Pole in ICRS (approximately RA=192.86°, Dec=27.13°)
        gal_l, gal_b = icrs_to_galactic(192.86, 27.13)

        # Should have b ~ 90° (North Galactic Pole)
        assert gal_b > 85

    def test_galactic_to_icrs_roundtrip(self):
        """Test Galactic to ICRS conversion and roundtrip."""
        from libephemeris.astropy_integration import icrs_to_galactic, galactic_to_icrs

        # Start with known ICRS coordinates
        ra_orig, dec_orig = 100.0, 30.0

        # Convert to Galactic and back
        gal_l, gal_b = icrs_to_galactic(ra_orig, dec_orig)
        ra_back, dec_back = galactic_to_icrs(gal_l, gal_b)

        # Should get back the same coordinates
        assert abs(ra_back - ra_orig) < 0.0001
        assert abs(dec_back - dec_orig) < 0.0001


# =============================================================================
# GEODETIC UTILITY TESTS
# =============================================================================


class TestGeodeticUtilities:
    """Tests for geodetic coordinate conversions."""

    def test_geodetic_to_geocentric_equator(self):
        """Test conversion at the equator."""
        from libephemeris.astropy_integration import geodetic_to_geocentric

        # Point on equator at prime meridian, on ellipsoid surface
        x, y, z = geodetic_to_geocentric(0, 0, 0)

        # Earth's equatorial radius is ~6378 km = 6.378e6 m
        assert 6.37e6 < x < 6.39e6
        assert abs(y) < 1000  # Should be very close to 0
        assert abs(z) < 1000  # Should be very close to 0

    def test_geodetic_to_geocentric_pole(self):
        """Test conversion at the North Pole."""
        from libephemeris.astropy_integration import geodetic_to_geocentric

        # North Pole on ellipsoid surface
        x, y, z = geodetic_to_geocentric(0, 90, 0)

        # Earth's polar radius is ~6357 km = 6.357e6 m
        assert abs(x) < 1000  # Should be very close to 0
        assert abs(y) < 1000  # Should be very close to 0
        assert 6.35e6 < z < 6.36e6

    def test_geocentric_to_geodetic_roundtrip(self):
        """Test roundtrip conversion."""
        from libephemeris.astropy_integration import (
            geodetic_to_geocentric,
            geocentric_to_geodetic,
        )

        # Original geodetic coordinates
        lon_orig, lat_orig, h_orig = 12.5, 41.9, 100.0  # Rome, 100m altitude

        # Convert to geocentric and back
        x, y, z = geodetic_to_geocentric(lon_orig, lat_orig, h_orig)
        lon_back, lat_back, h_back = geocentric_to_geodetic(x, y, z)

        # Should get back the same coordinates
        assert abs(lon_back - lon_orig) < 0.0001
        assert abs(lat_back - lat_orig) < 0.0001
        assert abs(h_back - h_orig) < 1  # Within 1 meter


# =============================================================================
# BARYCENTRIC CORRECTION TESTS
# =============================================================================


class TestBarycentricCorrections:
    """Tests for barycentric correction calculations."""

    def test_get_barycentric_correction(self):
        """Test barycentric correction calculation."""
        from libephemeris.astropy_integration import get_barycentric_correction

        result = get_barycentric_correction(
            ra_deg=180.0,  # Arbitrary direction
            dec_deg=45.0,
            jd_utc=2451545.0,
            lon_deg=0.0,  # Greenwich
            lat_deg=51.5,  # London latitude
        )

        # Check structure
        assert "bjd_tdb" in result
        assert "light_time_seconds" in result
        assert "rv_correction_km_s" in result

        # BJD should be close to input JD (within ~500 seconds = ~0.006 days)
        assert abs(result["bjd_tdb"] - 2451545.0) < 0.01

        # Light travel time should be less than ~500 seconds (8+ min to Sun)
        assert abs(result["light_time_seconds"]) < 600

        # RV correction should be less than ~30 km/s (Earth orbital velocity)
        assert abs(result["rv_correction_km_s"]) < 35

    def test_barycentric_correction_varies_with_direction(self):
        """Test that barycentric correction varies with viewing direction."""
        from libephemeris.astropy_integration import get_barycentric_correction

        jd = 2451545.0
        lon, lat = 0.0, 50.0

        # Two opposite directions
        result1 = get_barycentric_correction(0, 0, jd, lon, lat)
        result2 = get_barycentric_correction(180, 0, jd, lon, lat)

        # RV corrections should have opposite signs (roughly)
        # or at least be significantly different
        assert abs(result1["rv_correction_km_s"] - result2["rv_correction_km_s"]) > 1


# =============================================================================
# INTEGRATION EVALUATION TESTS
# =============================================================================


class TestIntegrationEvaluation:
    """Tests for the overall integration evaluation."""

    def test_evaluation_provides_recommendations(self):
        """Test that evaluation provides actionable recommendations."""
        result = evaluate_integration_potential()

        # Should have recommendations
        assert len(result["recommended_integrations"]) >= 3

        # Each recommendation should be a non-empty string
        for rec in result["recommended_integrations"]:
            assert isinstance(rec, str)
            assert len(rec) > 10

    def test_evaluation_lists_compatibility_notes(self):
        """Test that evaluation notes Skyfield/Astropy compatibility."""
        result = evaluate_integration_potential()

        # Should have compatibility notes
        assert len(result["compatibility_notes"]) >= 3

        # Notes should mention both Skyfield and ICRS
        notes_text = " ".join(result["compatibility_notes"]).lower()
        assert "skyfield" in notes_text or "icrs" in notes_text


# =============================================================================
# ERROR HANDLING TESTS
# =============================================================================


class TestErrorHandling:
    """Tests for proper error handling."""

    def test_invalid_time_string_raises_error(self):
        """Test that invalid time strings raise appropriate errors."""
        from libephemeris.astropy_integration import parse_time_string

        with pytest.raises((ValueError, Exception)):
            parse_time_string("not a valid time")

    def test_coordinate_transform_with_extreme_values(self):
        """Test coordinate transforms with edge case values."""
        from libephemeris.astropy_integration import icrs_to_galactic

        # Test at poles
        l1, b1 = icrs_to_galactic(0, 90)  # North celestial pole
        l2, b2 = icrs_to_galactic(0, -90)  # South celestial pole

        # Both should give valid results
        assert -90 <= b1 <= 90
        assert -90 <= b2 <= 90
        assert 0 <= l1 < 360 or abs(l1) < 0.001
        assert 0 <= l2 < 360 or abs(l2) < 0.001


# =============================================================================
# DOCUMENTATION VALIDATION TESTS
# =============================================================================


class TestDocumentationValidation:
    """Tests that validate claims made in ASTROPY_BENEFITS.md documentation."""

    def test_documented_time_scales_available(self):
        """Verify all time scales documented in ASTROPY_BENEFITS.md are available."""
        from libephemeris.astropy_integration import get_extended_time_scales

        scales = get_extended_time_scales(2451545.0)

        # Documented time scales: UTC, UT1, TAI, TT, TDB, TCG, TCB, GPS
        documented_scales = ["utc", "ut1", "tai", "tt", "tdb", "tcg", "tcb", "gps"]
        for scale in documented_scales:
            assert scale in scales, f"Documented time scale '{scale}' not available"
            assert isinstance(scales[scale], float)

    def test_documented_galactic_conversion_works(self):
        """Verify Galactic coordinate conversion documented in ASTROPY_BENEFITS.md works."""
        from libephemeris.astropy_integration import icrs_to_galactic, galactic_to_icrs

        # Documentation claims: Galactic coords not in Skyfield, available via astropy
        ra, dec = 180.0, 45.0
        gal_l, gal_b = icrs_to_galactic(ra, dec)

        # Should return valid Galactic coordinates
        assert 0 <= gal_l < 360 or abs(gal_l) < 0.001 or abs(gal_l - 360) < 0.001
        assert -90 <= gal_b <= 90

        # Roundtrip should work
        ra_back, dec_back = galactic_to_icrs(gal_l, gal_b)
        assert abs(ra_back - ra) < 0.001
        assert abs(dec_back - dec) < 0.001

    def test_documented_barycentric_correction_structure(self):
        """Verify barycentric correction returns documented fields."""
        from libephemeris.astropy_integration import get_barycentric_correction

        result = get_barycentric_correction(
            ra_deg=180.0,
            dec_deg=45.0,
            jd_utc=2451545.0,
            lon_deg=0.0,
            lat_deg=50.0,
        )

        # Documentation claims these specific fields
        assert "bjd_tdb" in result, "Missing documented field 'bjd_tdb'"
        assert "light_time_seconds" in result, (
            "Missing documented field 'light_time_seconds'"
        )
        assert "rv_correction_km_s" in result, (
            "Missing documented field 'rv_correction_km_s'"
        )

    def test_documented_geodetic_functions_exist(self):
        """Verify geodetic functions documented in ASTROPY_BENEFITS.md exist."""
        from libephemeris.astropy_integration import (
            geodetic_to_geocentric,
            geocentric_to_geodetic,
        )

        # Documentation claims WGS84 ellipsoid conversions
        x, y, z = geodetic_to_geocentric(12.5, 41.9, 100.0)
        assert isinstance(x, float)
        assert isinstance(y, float)
        assert isinstance(z, float)

        # Should be able to convert back
        lon, lat, h = geocentric_to_geodetic(x, y, z)
        assert abs(lon - 12.5) < 0.001
        assert abs(lat - 41.9) < 0.001
        assert abs(h - 100.0) < 1.0

    def test_documented_evaluation_function_keys(self):
        """Verify evaluate_integration_potential returns documented keys."""
        from libephemeris.astropy_integration import evaluate_integration_potential

        result = evaluate_integration_potential()

        # Documentation shows these specific keys
        documented_keys = [
            "astropy_available",
            "time_features",
            "coordinate_features",
            "recommended_integrations",
            "compatibility_notes",
            "implementation_approach",
        ]
        for key in documented_keys:
            assert key in result, f"Missing documented key '{key}'"

    def test_time_string_parsing_works(self):
        """Verify time string parsing documented in ASTROPY_BENEFITS.md works."""
        from libephemeris.astropy_integration import parse_time_string

        # Documentation shows parsing ISO format strings
        jd = parse_time_string("2000-01-01T12:00:00")
        assert isinstance(jd, float)
        # Should be close to J2000.0
        assert abs(jd - 2451545.0) < 0.01
