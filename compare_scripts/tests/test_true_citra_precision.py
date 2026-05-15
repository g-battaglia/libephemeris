"""
Tests for improved True Citra ayanamsha precision.

These tests verify that the True Citra ayanamsha calculation using high-precision
Hipparcos coordinates for Spica (HIP 65474) with full proper motion correction
(including parallax and radial velocity) achieves better precision than before.

The improvement targets:
- Previous precision: ±0.06° vs Swiss Ephemeris
- Target precision: ±0.02° vs Swiss Ephemeris (3x improvement)
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import SIDM_TRUE_CITRA, SUN, FLG_SIDEREAL


class TestTrueCitraPrecision:
    """Test improved True Citra ayanamsha precision."""

    @pytest.fixture(autouse=True)
    def reset_sidereal_mode(self):
        """Reset sidereal mode after each test."""
        yield
        swe.set_sid_mode(0)
        ephem.set_sid_mode(0)

    def test_true_citra_at_j2000(self):
        """Test True Citra ayanamsha at J2000.0 epoch."""
        jd = 2451545.0  # J2000.0

        # Swiss Ephemeris reference
        swe.set_sid_mode(SIDM_TRUE_CITRA)
        aya_swe = swe.get_ayanamsa_ut(jd)

        # LibEphemeris with improved Spica coordinates
        ephem.set_sid_mode(SIDM_TRUE_CITRA)
        aya_lib = ephem.get_ayanamsa_ut(jd)

        diff = abs(aya_swe - aya_lib)

        # Target: <0.02° precision (previous was ~0.06°)
        assert diff < 0.02, (
            f"True Citra at J2000.0: SWE={aya_swe:.6f}°, LIB={aya_lib:.6f}°, "
            f"Diff={diff:.6f}° (target <0.02°)"
        )

    def test_true_citra_at_multiple_dates(self):
        """Test True Citra ayanamsha at multiple dates across DE421 range."""
        test_dates = [
            (1900, 1, 1, 12.0, "1900"),
            (1950, 1, 1, 12.0, "1950"),
            (2000, 1, 1, 12.0, "2000"),
            (2024, 1, 1, 12.0, "2024"),
            (2050, 1, 1, 12.0, "2050"),
        ]

        max_diff = 0.0
        for year, month, day, hour, label in test_dates:
            jd = swe.julday(year, month, day, hour)

            swe.set_sid_mode(SIDM_TRUE_CITRA)
            aya_swe = swe.get_ayanamsa_ut(jd)

            ephem.set_sid_mode(SIDM_TRUE_CITRA)
            aya_lib = ephem.get_ayanamsa_ut(jd)

            diff = abs(aya_swe - aya_lib)
            max_diff = max(max_diff, diff)

            # Should be <0.02° at all dates
            assert diff < 0.02, (
                f"True Citra at {label}: SWE={aya_swe:.6f}°, LIB={aya_lib:.6f}°, "
                f"Diff={diff:.6f}° (target <0.02°)"
            )

        print(f"\nMax True Citra difference across all dates: {max_diff:.6f}°")

    def test_true_citra_proper_motion_effect(self):
        """Test that proper motion is correctly applied over time."""
        # Ayanamsha should change smoothly over time due to precession
        # and star proper motion

        jd_1950 = swe.julday(1950, 1, 1, 12.0)
        jd_2050 = swe.julday(2050, 1, 1, 12.0)

        ephem.set_sid_mode(SIDM_TRUE_CITRA)
        aya_1950 = ephem.get_ayanamsa_ut(jd_1950)
        aya_2050 = ephem.get_ayanamsa_ut(jd_2050)

        # Ayanamsha increases due to precession (~50.3"/year ≈ 1.397°/century)
        # Over 100 years: expect ~1.4° increase
        change = aya_2050 - aya_1950
        expected_change = 1.397  # degrees per century (approximate)

        # The change should be close to the precession rate
        assert abs(change - expected_change) < 0.1, (
            f"Ayanamsha change 1950-2050: {change:.4f}° "
            f"(expected ~{expected_change:.3f}°)"
        )

    def test_sidereal_sun_position_true_citra(self):
        """Test sidereal Sun position using True Citra matches Swiss Ephemeris."""
        jd = swe.julday(2024, 4, 14, 12.0)  # Sun in Aries (tropical)

        # Swiss Ephemeris
        swe.set_sid_mode(SIDM_TRUE_CITRA)
        pos_swe, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)
        sun_swe = pos_swe[0]

        # LibEphemeris
        ephem.set_sid_mode(SIDM_TRUE_CITRA)
        pos_lib, _ = ephem.calc_ut(jd, SUN, FLG_SIDEREAL)
        sun_lib = pos_lib[0]

        diff = abs(sun_swe - sun_lib)
        if diff > 180:
            diff = 360 - diff

        # Sidereal position should match closely
        # Planet position accuracy + ayanamsha accuracy = total accuracy
        assert diff < 0.03, (
            f"Sidereal Sun (True Citra): SWE={sun_swe:.6f}°, "
            f"LIB={sun_lib:.6f}°, Diff={diff:.6f}°"
        )


class TestSpicaStarData:
    """Test the improved Spica star data used for True Citra."""

    def test_spica_coordinates_precision(self):
        """Test that Spica coordinates match Hipparcos catalog."""
        from libephemeris.planets import STARS

        spica = STARS["SPICA"]

        # Hipparcos values for HIP 65474 (Spica)
        # RA: 13h 25m 11.5794s = 201.2982475°
        # Dec: -11° 09' 40.759" = -11.1613219°
        expected_ra = 201.2982475
        expected_dec = -11.1613219

        # Check position matches Hipparcos within 0.0001°
        assert abs(spica.ra_j2000 - expected_ra) < 0.0001, (
            f"Spica RA: {spica.ra_j2000} vs expected {expected_ra}"
        )
        assert abs(spica.dec_j2000 - expected_dec) < 0.0001, (
            f"Spica Dec: {spica.dec_j2000} vs expected {expected_dec}"
        )

    def test_spica_proper_motion_values(self):
        """Test that Spica proper motion matches Hipparcos catalog."""
        from libephemeris.planets import STARS

        spica = STARS["SPICA"]

        # Hipparcos proper motion for HIP 65474
        # pm_ra (mu_alpha*): -42.50 mas/yr = -0.04250 arcsec/yr
        # pm_dec (mu_delta): -31.73 mas/yr = -0.03173 arcsec/yr
        expected_pm_ra = -0.04250
        expected_pm_dec = -0.03173

        # Check proper motion matches Hipparcos
        assert abs(spica.pm_ra - expected_pm_ra) < 0.0001, (
            f"Spica pm_ra: {spica.pm_ra} vs expected {expected_pm_ra}"
        )
        assert abs(spica.pm_dec - expected_pm_dec) < 0.0001, (
            f"Spica pm_dec: {spica.pm_dec} vs expected {expected_pm_dec}"
        )

    def test_spica_parallax_and_radial_velocity(self):
        """Test that Spica parallax and radial velocity are set."""
        from libephemeris.planets import STARS

        spica = STARS["SPICA"]

        # Hipparcos parallax: 12.44 mas = 0.01244 arcsec
        # Radial velocity: +1.0 km/s
        assert spica.parallax > 0, "Spica should have parallax set"
        assert abs(spica.parallax - 0.01244) < 0.001, (
            f"Spica parallax: {spica.parallax} vs expected 0.01244"
        )
        assert spica.radial_velocity == 1.0, (
            f"Spica radial_velocity: {spica.radial_velocity} vs expected 1.0"
        )


class TestStarDataClass:
    """Test the StarData dataclass with optional fields."""

    def test_star_data_default_values(self):
        """Test that StarData defaults parallax and radial_velocity to 0."""
        from libephemeris.planets import StarData

        # Create StarData without optional fields
        star = StarData(
            ra_j2000=100.0,
            dec_j2000=20.0,
            pm_ra=0.01,
            pm_dec=-0.02,
        )

        assert star.parallax == 0.0, "Default parallax should be 0.0"
        assert star.radial_velocity == 0.0, "Default radial_velocity should be 0.0"

    def test_star_data_with_parallax(self):
        """Test StarData with parallax and radial velocity specified."""
        from libephemeris.planets import StarData

        star = StarData(
            ra_j2000=100.0,
            dec_j2000=20.0,
            pm_ra=0.01,
            pm_dec=-0.02,
            parallax=0.05,
            radial_velocity=-10.0,
        )

        assert star.parallax == 0.05
        assert star.radial_velocity == -10.0


if __name__ == "__main__":
    # Run tests with verbose output
    import sys

    pytest.main([__file__, "-v", "--tb=short"] + sys.argv[1:])
