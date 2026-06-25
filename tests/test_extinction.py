"""
Tests for atmospheric extinction model in libephemeris.

Tests the calculation of atmospheric extinction for heliacal visibility
calculations based on Schaefer and Green formulas.

The basic extinction formula is:
    extinction (mag) = k * sec(z) ≈ 0.28 * sec(z)
where z is zenith angle and k is the extinction coefficient (~0.28 for V-band).

References:
    - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
    - Green, D.W.E. (1992) "Magnitude Corrections for Atmospheric Extinction"
    - Kasten, F. & Young, A.T. (1989) Applied Optics 28, 4735-4738
"""

import math

from libephemeris import (
    calc_airmass,
    calc_extinction_coefficient,
    calc_extinction_magnitude,
    calc_simple_extinction,
    apparent_magnitude_with_extinction,
    get_extinction_for_heliacal,
    calc_rayleigh_coefficient,
    calc_aerosol_coefficient,
    calc_ozone_coefficient,
    calc_water_vapor_coefficient,
    ExtinctionCoefficients,
    WAVELENGTH_U,
    WAVELENGTH_B,
    WAVELENGTH_V,
    WAVELENGTH_R,
    WAVELENGTH_I,
)


class TestAirmassCalculation:
    """Tests for airmass calculation."""

    def test_airmass_at_zenith(self):
        """Airmass should be 1.0 at zenith (90 degrees altitude)."""
        airmass = calc_airmass(90.0)
        assert abs(airmass - 1.0) < 0.001

    def test_airmass_at_45_degrees(self):
        """Airmass should be ~1.41 at 45 degrees (sec(45) = sqrt(2))."""
        airmass = calc_airmass(45.0)
        expected = 1.0 / math.cos(math.radians(45))  # sqrt(2) ≈ 1.414
        assert (
            abs(airmass - expected) < 0.05
        )  # Allow small deviation from simple formula

    def test_airmass_at_30_degrees(self):
        """Airmass should be ~2.0 at 30 degrees altitude (zenith angle 60)."""
        airmass = calc_airmass(30.0)
        # Simple formula: sec(60) = 2.0
        # Kasten-Young gives slightly different value
        assert 1.9 < airmass < 2.1

    def test_airmass_increases_toward_horizon(self):
        """Airmass should increase as altitude decreases."""
        altitudes = [90, 60, 45, 30, 15, 10, 5]
        airmasses = [calc_airmass(alt) for alt in altitudes]

        # Each value should be greater than the previous
        for i in range(1, len(airmasses)):
            assert airmasses[i] > airmasses[i - 1]

    def test_airmass_at_horizon(self):
        """Airmass should be high (but finite) at the horizon."""
        airmass = calc_airmass(0.0)
        # At horizon, airmass is typically 38-40
        assert 35 < airmass <= 40

    def test_airmass_below_horizon(self):
        """Airmass should return maximum value for objects below horizon."""
        airmass = calc_airmass(-5.0)
        assert airmass == 40.0

    def test_airmass_near_horizon(self):
        """Test airmass at very low altitudes (1-5 degrees)."""
        # At 1 degree altitude
        airmass_1deg = calc_airmass(1.0)
        assert 20 < airmass_1deg < 40

        # At 5 degrees altitude
        airmass_5deg = calc_airmass(5.0)
        assert 8 < airmass_5deg < 15

    def test_airmass_secant_method(self):
        """Test simple secant method for comparison."""
        # At 60 degrees, secant method should give exact sec(30) = 1.155
        airmass = calc_airmass(60.0, method="secant")
        expected = 1.0 / math.cos(math.radians(30))
        assert abs(airmass - expected) < 0.01

    def test_airmass_rozenberg_method(self):
        """Test Rozenberg method."""
        airmass = calc_airmass(30.0, method="rozenberg")
        assert 1.9 < airmass < 2.2


class TestExtinctionCoefficients:
    """Tests for individual extinction coefficient components."""

    def test_rayleigh_standard_conditions(self):
        """Test Rayleigh coefficient at standard conditions."""
        k_rayleigh = calc_rayleigh_coefficient(
            wavelength_nm=550.0,
            pressure_mbar=1013.25,
            temperature_c=15.0,
            altitude_m=0.0,
        )
        # Standard V-band Rayleigh coefficient is ~0.145
        assert 0.14 < k_rayleigh < 0.16

    def test_rayleigh_wavelength_dependence(self):
        """Rayleigh scattering should be stronger at shorter wavelengths."""
        k_blue = calc_rayleigh_coefficient(wavelength_nm=440.0)
        k_visual = calc_rayleigh_coefficient(wavelength_nm=550.0)
        k_red = calc_rayleigh_coefficient(wavelength_nm=640.0)

        # Blue should have strongest scattering (lambda^-4 dependence)
        assert k_blue > k_visual > k_red

    def test_rayleigh_pressure_dependence(self):
        """Higher pressure should increase Rayleigh scattering."""
        k_low_pressure = calc_rayleigh_coefficient(pressure_mbar=800.0)
        k_standard = calc_rayleigh_coefficient(pressure_mbar=1013.25)
        k_high_pressure = calc_rayleigh_coefficient(pressure_mbar=1050.0)

        assert k_low_pressure < k_standard < k_high_pressure

    def test_rayleigh_altitude_dependence(self):
        """Higher altitude should decrease Rayleigh scattering."""
        k_sea_level = calc_rayleigh_coefficient(altitude_m=0.0)
        k_mountain = calc_rayleigh_coefficient(altitude_m=2000.0)
        k_high_mountain = calc_rayleigh_coefficient(altitude_m=4000.0)

        assert k_sea_level > k_mountain > k_high_mountain

    def test_aerosol_standard_conditions(self):
        """Test aerosol coefficient at standard conditions."""
        k_aerosol = calc_aerosol_coefficient(humidity_percent=50.0, altitude_m=0.0)
        # Typical aerosol coefficient is 0.05-0.20
        assert 0.05 < k_aerosol < 0.25

    def test_aerosol_humidity_dependence(self):
        """Higher humidity should increase aerosol scattering."""
        k_dry = calc_aerosol_coefficient(humidity_percent=20.0)
        k_moderate = calc_aerosol_coefficient(humidity_percent=50.0)
        k_humid = calc_aerosol_coefficient(humidity_percent=90.0)

        assert k_dry < k_moderate < k_humid

    def test_aerosol_altitude_dependence(self):
        """Higher altitude should decrease aerosol scattering."""
        k_sea_level = calc_aerosol_coefficient(altitude_m=0.0)
        k_mountain = calc_aerosol_coefficient(altitude_m=2000.0)

        assert k_sea_level > k_mountain

    def test_aerosol_from_visibility(self):
        """Test aerosol coefficient calculation from visibility."""
        # High visibility = low aerosol
        k_clear = calc_aerosol_coefficient(visibility_km=50.0)
        k_hazy = calc_aerosol_coefficient(visibility_km=10.0)

        assert k_clear < k_hazy

    def test_ozone_coefficient(self):
        """Test ozone absorption coefficient."""
        k_ozone_v = calc_ozone_coefficient(wavelength_nm=550.0)
        k_ozone_uv = calc_ozone_coefficient(wavelength_nm=350.0)

        # Ozone absorption should be higher in UV
        assert k_ozone_uv > k_ozone_v
        # V-band ozone should be small
        assert k_ozone_v < 0.05

    def test_water_vapor_coefficient(self):
        """Test water vapor absorption coefficient."""
        k_water_v = calc_water_vapor_coefficient(
            humidity_percent=50.0, wavelength_nm=550.0
        )
        k_water_ir = calc_water_vapor_coefficient(
            humidity_percent=50.0, wavelength_nm=850.0
        )

        # Water vapor absorption higher in IR
        assert k_water_ir > k_water_v
        # V-band water vapor should be very small
        assert k_water_v < 0.01


class TestTotalExtinctionCoefficient:
    """Tests for total extinction coefficient calculation."""

    def test_standard_extinction_coefficient(self):
        """Test total extinction at standard conditions."""
        coeff = calc_extinction_coefficient()

        assert isinstance(coeff, ExtinctionCoefficients)
        # Total should be around 0.25-0.35 for typical conditions
        assert 0.20 < coeff.k_total < 0.40

    def test_extinction_components_sum(self):
        """Components should sum to total."""
        coeff = calc_extinction_coefficient()

        calculated_total = (
            coeff.k_rayleigh + coeff.k_aerosol + coeff.k_ozone + coeff.k_water
        )
        assert abs(calculated_total - coeff.k_total) < 0.001

    def test_extinction_excellent_conditions(self):
        """Test extinction at excellent observing site."""
        coeff = calc_extinction_coefficient(
            humidity_percent=20.0, altitude_m=2500.0, pressure_mbar=750.0
        )
        # Should be lower than standard
        assert coeff.k_total < 0.25

    def test_extinction_poor_conditions(self):
        """Test extinction at poor conditions."""
        coeff = calc_extinction_coefficient(
            humidity_percent=90.0, altitude_m=0.0, pressure_mbar=1020.0
        )
        # Should be higher than standard
        assert coeff.k_total > 0.30


class TestSimpleExtinction:
    """Tests for the simple Schaefer/Green extinction formula."""

    def test_simple_extinction_at_zenith(self):
        """At zenith, extinction should equal k."""
        extinction = calc_simple_extinction(90.0, k=0.28)
        assert abs(extinction - 0.28) < 0.01

    def test_simple_extinction_at_30_degrees(self):
        """At 30 degrees (airmass ~2), extinction should be ~2k."""
        extinction = calc_simple_extinction(30.0, k=0.28)
        # Should be approximately 0.28 * 2 = 0.56
        assert 0.50 < extinction < 0.65

    def test_simple_extinction_increases_toward_horizon(self):
        """Extinction should increase as altitude decreases."""
        altitudes = [90, 60, 45, 30, 15, 10, 5]
        extinctions = [calc_simple_extinction(alt) for alt in altitudes]

        for i in range(1, len(extinctions)):
            assert extinctions[i] > extinctions[i - 1]

    def test_simple_extinction_below_horizon(self):
        """Extinction should be very high for objects below horizon."""
        extinction = calc_simple_extinction(-5.0)
        assert extinction >= 10.0

    def test_simple_extinction_default_k(self):
        """Default k should be ~0.28."""
        extinction_45 = calc_simple_extinction(45.0)
        # At 45 degrees, airmass ~ 1.41
        # Extinction should be ~0.28 * 1.41 = 0.40
        assert 0.35 < extinction_45 < 0.50


class TestExtinctionMagnitude:
    """Tests for full extinction magnitude calculation."""

    def test_extinction_magnitude_at_zenith(self):
        """Extinction at zenith should be minimal (~0.28 mag)."""
        extinction = calc_extinction_magnitude(90.0)
        assert 0.20 < extinction < 0.40

    def test_extinction_magnitude_at_30_degrees(self):
        """Extinction at 30 degrees should be moderate."""
        extinction = calc_extinction_magnitude(30.0)
        # Airmass ~2, so extinction ~0.5-0.7 mag
        assert 0.45 < extinction < 0.80

    def test_extinction_magnitude_near_horizon(self):
        """Extinction near horizon should be several magnitudes."""
        extinction = calc_extinction_magnitude(5.0)
        # At 5 degrees, airmass ~10-12, so extinction ~3-4 mag
        assert 2.0 < extinction < 5.0

    def test_extinction_magnitude_below_horizon(self):
        """Objects below horizon should have very high extinction."""
        extinction = calc_extinction_magnitude(-5.0)
        assert extinction >= 90.0  # Essentially infinite

    def test_extinction_humidity_effect(self):
        """Higher humidity should increase extinction."""
        extinction_dry = calc_extinction_magnitude(30.0, humidity_percent=20.0)
        extinction_humid = calc_extinction_magnitude(30.0, humidity_percent=90.0)

        assert extinction_humid > extinction_dry

    def test_extinction_altitude_effect(self):
        """Higher observer altitude should decrease extinction."""
        extinction_sea = calc_extinction_magnitude(30.0, observer_altitude_m=0.0)
        extinction_mountain = calc_extinction_magnitude(
            30.0, observer_altitude_m=2500.0
        )

        assert extinction_mountain < extinction_sea


class TestApparentMagnitude:
    """Tests for apparent magnitude with extinction."""

    def test_apparent_magnitude_at_zenith(self):
        """Object at zenith should be slightly fainter than catalog."""
        catalog_mag = 0.0
        apparent = apparent_magnitude_with_extinction(catalog_mag, 90.0)
        # Should be ~0.28 mag fainter
        assert 0.2 < apparent < 0.4

    def test_apparent_magnitude_at_30_degrees(self):
        """Object at 30 degrees should be noticeably fainter."""
        catalog_mag = 0.0
        apparent = apparent_magnitude_with_extinction(catalog_mag, 30.0)
        # Should be ~0.5-0.7 mag fainter
        assert 0.4 < apparent < 0.8

    def test_venus_near_horizon(self):
        """Venus (mag -4) near horizon should still be visible but dimmer."""
        catalog_mag = -4.0
        apparent_5deg = apparent_magnitude_with_extinction(catalog_mag, 5.0)
        # Venus should still be bright but dimmed by 2-4 mag
        assert -2.0 < apparent_5deg < 1.0

    def test_sirius_at_various_altitudes(self):
        """Sirius (mag -1.46) visibility at various altitudes."""
        catalog_mag = -1.46

        apparent_zenith = apparent_magnitude_with_extinction(catalog_mag, 90.0)
        apparent_30deg = apparent_magnitude_with_extinction(catalog_mag, 30.0)
        apparent_10deg = apparent_magnitude_with_extinction(catalog_mag, 10.0)

        # Should get progressively fainter
        assert apparent_zenith < apparent_30deg < apparent_10deg

    def test_faint_star_visibility(self):
        """A 6th magnitude star should be visible at zenith, not at horizon."""
        catalog_mag = 6.0

        apparent_zenith = apparent_magnitude_with_extinction(catalog_mag, 90.0)
        apparent_5deg = apparent_magnitude_with_extinction(catalog_mag, 5.0)

        # At zenith, still close to 6th magnitude (barely visible)
        assert apparent_zenith < 7.0
        # Near horizon, should be much fainter than visible limit
        assert apparent_5deg > 8.0


class TestExtinctionForHeliacal:
    """Tests for heliacal visibility-specific extinction function."""

    def test_heliacal_extinction_standard(self):
        """Test standard extinction for heliacal calculations."""
        k = get_extinction_for_heliacal(zenith_angle_deg=60.0)
        # Should return reasonable k value
        assert 0.20 < k < 0.40

    def test_heliacal_extinction_humidity_effect(self):
        """Humidity should affect heliacal extinction."""
        k_dry = get_extinction_for_heliacal(60.0, humidity_percent=20.0)
        k_humid = get_extinction_for_heliacal(60.0, humidity_percent=90.0)

        assert k_humid > k_dry


class TestWavelengthConstants:
    """Tests for wavelength constants."""

    def test_wavelength_constants_defined(self):
        """All wavelength constants should be defined."""
        assert WAVELENGTH_U == 365.0
        assert WAVELENGTH_B == 440.0
        assert WAVELENGTH_V == 550.0
        assert WAVELENGTH_R == 640.0
        assert WAVELENGTH_I == 790.0

    def test_wavelength_order(self):
        """Wavelengths should be in order from UV to IR."""
        assert WAVELENGTH_U < WAVELENGTH_B < WAVELENGTH_V < WAVELENGTH_R < WAVELENGTH_I


class TestExtinctionPhysicalReasonability:
    """Tests to verify physically reasonable behavior."""

    def test_extinction_never_negative(self):
        """Extinction should never be negative."""
        for altitude in [90, 60, 45, 30, 15, 10, 5, 1]:
            extinction = calc_extinction_magnitude(altitude)
            assert extinction >= 0

    def test_extinction_smooth_variation(self):
        """Extinction should vary smoothly with altitude."""
        altitudes = list(range(5, 91, 5))
        extinctions = [calc_extinction_magnitude(alt) for alt in altitudes]

        # Check for monotonic decrease as altitude increases
        for i in range(1, len(extinctions)):
            assert extinctions[i] < extinctions[i - 1]

    def test_extreme_conditions(self):
        """Model should handle extreme conditions gracefully."""
        # Very high altitude
        coeff_high = calc_extinction_coefficient(altitude_m=5000.0)
        assert coeff_high.k_total > 0

        # Very low pressure
        coeff_low_p = calc_extinction_coefficient(pressure_mbar=500.0)
        assert coeff_low_p.k_total > 0

        # 100% humidity
        coeff_humid = calc_extinction_coefficient(humidity_percent=100.0)
        assert coeff_humid.k_total > 0
        assert coeff_humid.k_total < 1.0  # Should still be reasonable


class TestImportability:
    """Tests to verify proper import from libephemeris."""

    def test_all_functions_importable(self):
        """All extinction functions should be importable from libephemeris."""
        from libephemeris import (
            calc_airmass,
            calc_extinction_coefficient,
            calc_simple_extinction,
            WAVELENGTH_V,
        )

        # Verify they are callable/usable
        assert callable(calc_airmass)
        assert callable(calc_extinction_coefficient)
        assert callable(calc_simple_extinction)
        assert isinstance(WAVELENGTH_V, float)
