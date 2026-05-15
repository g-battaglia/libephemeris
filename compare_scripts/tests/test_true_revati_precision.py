"""
Tests for improved True Revati ayanamsha precision.

These tests verify that the True Revati ayanamsha calculation using high-precision
Gaia DR3 coordinates for Zeta Piscium (HIP 5737) with full proper motion correction
(including parallax and radial velocity) achieves better precision than before.

The improvement targets:
- Previous precision: ±0.06° vs Swiss Ephemeris
- Current precision: <0.006° vs Swiss Ephemeris (10x improvement)
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import SIDM_TRUE_REVATI, SUN, FLG_SIDEREAL


class TestTrueRevatiPrecision:
    """Test improved True Revati ayanamsha precision."""

    @pytest.fixture(autouse=True)
    def reset_sidereal_mode(self):
        """Reset sidereal mode after each test."""
        yield
        swe.set_sid_mode(0)
        ephem.set_sid_mode(0)

    def test_true_revati_at_j2000(self):
        """Test True Revati ayanamsha at J2000.0 epoch."""
        jd = 2451545.0  # J2000.0

        # Swiss Ephemeris reference
        swe.set_sid_mode(SIDM_TRUE_REVATI)
        aya_swe = swe.get_ayanamsa_ut(jd)

        # LibEphemeris with improved Zeta Piscium coordinates
        ephem.set_sid_mode(SIDM_TRUE_REVATI)
        aya_lib = ephem.get_ayanamsa_ut(jd)

        diff = abs(aya_swe - aya_lib)

        # Target: <0.006° precision (calibrated offset at J2000)
        assert diff < 0.006, (
            f"True Revati at J2000.0: SWE={aya_swe:.6f}°, LIB={aya_lib:.6f}°, "
            f"Diff={diff:.6f}° (target <0.006°)"
        )

    def test_true_revati_at_multiple_dates(self):
        """Test True Revati ayanamsha at multiple dates across DE421 range."""
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

            swe.set_sid_mode(SIDM_TRUE_REVATI)
            aya_swe = swe.get_ayanamsa_ut(jd)

            ephem.set_sid_mode(SIDM_TRUE_REVATI)
            aya_lib = ephem.get_ayanamsa_ut(jd)

            diff = abs(aya_swe - aya_lib)
            max_diff = max(max_diff, diff)

            # Should be <0.006° at all dates
            assert diff < 0.006, (
                f"True Revati at {label}: SWE={aya_swe:.6f}°, LIB={aya_lib:.6f}°, "
                f"Diff={diff:.6f}° (target <0.006°)"
            )

        print(f"\nMax True Revati difference across all dates: {max_diff:.6f}°")

    def test_true_revati_proper_motion_effect(self):
        """Test that proper motion is correctly applied over time."""
        # Ayanamsha should change smoothly over time due to precession
        # and star proper motion

        jd_1950 = swe.julday(1950, 1, 1, 12.0)
        jd_2050 = swe.julday(2050, 1, 1, 12.0)

        ephem.set_sid_mode(SIDM_TRUE_REVATI)
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

    def test_sidereal_sun_position_true_revati(self):
        """Test sidereal Sun position using True Revati matches Swiss Ephemeris."""
        jd = swe.julday(2024, 4, 14, 12.0)  # Sun in Aries (tropical)

        # Swiss Ephemeris
        swe.set_sid_mode(SIDM_TRUE_REVATI)
        pos_swe, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)
        sun_swe = pos_swe[0]

        # LibEphemeris
        ephem.set_sid_mode(SIDM_TRUE_REVATI)
        pos_lib, _ = ephem.calc_ut(jd, SUN, FLG_SIDEREAL)
        sun_lib = pos_lib[0]

        diff = abs(sun_swe - sun_lib)
        if diff > 180:
            diff = 360 - diff

        # Sidereal position should match closely
        # Planet position accuracy + ayanamsha accuracy = total accuracy
        assert diff < 0.03, (
            f"Sidereal Sun (True Revati): SWE={sun_swe:.6f}°, "
            f"LIB={sun_lib:.6f}°, Diff={diff:.6f}°"
        )


class TestRevatiStarData:
    """Test the improved Revati (Zeta Piscium) star data used for True Revati."""

    def test_revati_coordinates_precision(self):
        """Test that Revati coordinates match Gaia DR3 catalog."""
        from libephemeris.planets import STARS

        revati = STARS["REVATI"]

        # Gaia DR3 values for HIP 5737 (Zeta Piscium A)
        # RA: 01h 13m 43.8860s = 18.4328583349°
        # Dec: +07° 34' 31.296" = 7.5753601597°
        expected_ra = 18.4328583349
        expected_dec = 7.5753601597

        # Check position matches Gaia DR3 within 0.0001°
        assert abs(revati.ra_j2000 - expected_ra) < 0.0001, (
            f"Revati RA: {revati.ra_j2000} vs expected {expected_ra}"
        )
        assert abs(revati.dec_j2000 - expected_dec) < 0.0001, (
            f"Revati Dec: {revati.dec_j2000} vs expected {expected_dec}"
        )

    def test_revati_proper_motion_values(self):
        """Test that Revati proper motion matches Gaia DR3 catalog."""
        from libephemeris.planets import STARS

        revati = STARS["REVATI"]

        # Gaia DR3 proper motion for HIP 5737
        # pm_ra (mu_alpha*): 142.693 mas/yr = 0.142693 arcsec/yr
        # pm_dec (mu_delta): -53.051 mas/yr = -0.053051 arcsec/yr
        expected_pm_ra = 0.142693
        expected_pm_dec = -0.053051

        # Check proper motion matches Gaia DR3
        assert abs(revati.pm_ra - expected_pm_ra) < 0.0001, (
            f"Revati pm_ra: {revati.pm_ra} vs expected {expected_pm_ra}"
        )
        assert abs(revati.pm_dec - expected_pm_dec) < 0.0001, (
            f"Revati pm_dec: {revati.pm_dec} vs expected {expected_pm_dec}"
        )

    def test_revati_parallax_and_radial_velocity(self):
        """Test that Revati parallax and radial velocity are set."""
        from libephemeris.planets import STARS

        revati = STARS["REVATI"]

        # Gaia DR3 parallax: 24.4595 mas = 0.0244595 arcsec
        # Radial velocity: +15.0 km/s
        assert revati.parallax > 0, "Revati should have parallax set"
        assert abs(revati.parallax - 0.0244595) < 0.001, (
            f"Revati parallax: {revati.parallax} vs expected 0.0244595"
        )
        assert revati.radial_velocity == 15.0, (
            f"Revati radial_velocity: {revati.radial_velocity} vs expected 15.0"
        )

    def test_revati_hip_number_correct(self):
        """Test that Revati uses the correct HIP number (5737, not 7097)."""
        # HIP 5737 = Zeta Piscium (Revati nakshatra star)
        # HIP 7097 = Alpha Piscium (Alrescha) - this was incorrectly used before
        from libephemeris.fixed_stars import STAR_NAME_TO_HIP

        # Verify STAR_NAME_TO_HIP has the correct HIP number
        assert STAR_NAME_TO_HIP.get("REVATI") == 5737, (
            f"REVATI should map to HIP 5737, got {STAR_NAME_TO_HIP.get('REVATI')}"
        )


if __name__ == "__main__":
    # Run tests with verbose output
    import sys

    pytest.main([__file__, "-v", "--tb=short"] + sys.argv[1:])
