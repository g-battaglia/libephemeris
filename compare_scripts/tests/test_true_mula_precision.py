"""
Tests for improved True Mula ayanamsha precision.

These tests verify that the True Mula ayanamsha calculation using high-precision
Hipparcos 2007 coordinates for Lambda Scorpii (HIP 85927) with full proper motion
correction (including parallax and radial velocity) achieves better precision.

Note: Lambda Scorpii (Shaula, V=1.63) is too bright for Gaia DR3, so we use
Hipparcos 2007 reprocessed data instead, following the pattern of True Citra (Spica).

The improvement targets:
- Previous precision: ~0.06 deg vs Swiss Ephemeris (without parallax/RV)
- Current precision: <0.006 deg vs Swiss Ephemeris (10x improvement)
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import SIDM_TRUE_MULA, SUN, FLG_SIDEREAL


class TestTrueMulaPrecision:
    """Test improved True Mula ayanamsha precision."""

    @pytest.fixture(autouse=True)
    def reset_sidereal_mode(self):
        """Reset sidereal mode after each test."""
        yield
        swe.set_sid_mode(0)
        ephem.set_sid_mode(0)

    def test_true_mula_at_j2000(self):
        """Test True Mula ayanamsha at J2000.0 epoch."""
        jd = 2451545.0  # J2000.0

        # Swiss Ephemeris reference
        swe.set_sid_mode(SIDM_TRUE_MULA)
        aya_swe = swe.get_ayanamsa_ut(jd)

        # LibEphemeris with improved Lambda Scorpii coordinates
        ephem.set_sid_mode(SIDM_TRUE_MULA)
        aya_lib = ephem.get_ayanamsa_ut(jd)

        diff = abs(aya_swe - aya_lib)

        # Target: <0.006 deg precision (calibrated offset at J2000)
        assert diff < 0.006, (
            f"True Mula at J2000.0: SWE={aya_swe:.6f} deg, LIB={aya_lib:.6f} deg, "
            f"Diff={diff:.6f} deg (target <0.006 deg)"
        )

    def test_true_mula_at_multiple_dates(self):
        """Test True Mula ayanamsha at multiple dates across DE421 range."""
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

            swe.set_sid_mode(SIDM_TRUE_MULA)
            aya_swe = swe.get_ayanamsa_ut(jd)

            ephem.set_sid_mode(SIDM_TRUE_MULA)
            aya_lib = ephem.get_ayanamsa_ut(jd)

            diff = abs(aya_swe - aya_lib)
            max_diff = max(max_diff, diff)

            # Should be <0.006 deg at all dates
            assert diff < 0.006, (
                f"True Mula at {label}: SWE={aya_swe:.6f} deg, LIB={aya_lib:.6f} deg, "
                f"Diff={diff:.6f} deg (target <0.006 deg)"
            )

        print(f"\nMax True Mula difference across all dates: {max_diff:.6f} deg")

    def test_true_mula_proper_motion_effect(self):
        """Test that proper motion is correctly applied over time."""
        # Ayanamsha should change smoothly over time due to precession
        # and star proper motion

        jd_1950 = swe.julday(1950, 1, 1, 12.0)
        jd_2050 = swe.julday(2050, 1, 1, 12.0)

        ephem.set_sid_mode(SIDM_TRUE_MULA)
        aya_1950 = ephem.get_ayanamsa_ut(jd_1950)
        aya_2050 = ephem.get_ayanamsa_ut(jd_2050)

        # Ayanamsha increases due to precession (~50.3"/year = ~1.397 deg/century)
        # Over 100 years: expect ~1.4 deg increase
        change = aya_2050 - aya_1950
        expected_change = 1.397  # degrees per century (approximate)

        # The change should be close to the precession rate
        assert abs(change - expected_change) < 0.1, (
            f"Ayanamsha change 1950-2050: {change:.4f} deg "
            f"(expected ~{expected_change:.3f} deg)"
        )

    def test_sidereal_sun_position_true_mula(self):
        """Test sidereal Sun position using True Mula matches Swiss Ephemeris."""
        jd = swe.julday(2024, 4, 14, 12.0)  # Sun in Aries (tropical)

        # Swiss Ephemeris
        swe.set_sid_mode(SIDM_TRUE_MULA)
        pos_swe, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)
        sun_swe = pos_swe[0]

        # LibEphemeris
        ephem.set_sid_mode(SIDM_TRUE_MULA)
        pos_lib, _ = ephem.calc_ut(jd, SUN, FLG_SIDEREAL)
        sun_lib = pos_lib[0]

        diff = abs(sun_swe - sun_lib)
        if diff > 180:
            diff = 360 - diff

        # Sidereal position should match closely
        # Planet position accuracy + ayanamsha accuracy = total accuracy
        assert diff < 0.02, (
            f"Sidereal Sun (True Mula): SWE={sun_swe:.6f} deg, "
            f"LIB={sun_lib:.6f} deg, Diff={diff:.6f} deg"
        )


class TestMulaStarData:
    """Test the improved Mula (Lambda Scorpii) star data used for True Mula."""

    def test_mula_coordinates_precision(self):
        """Test that Mula coordinates match Hipparcos 2007 catalog."""
        from libephemeris.planets import STARS

        mula = STARS["MULA"]

        # Hipparcos 2007 values for HIP 85927 (Lambda Scorpii)
        # RA: 17h 33m 36.52012s = 263.40216717 deg
        # Dec: -37 deg 06' 13.7648" = -37.10382356 deg
        expected_ra = 263.40216717
        expected_dec = -37.10382356

        # Check position matches Hipparcos 2007 within 0.0001 deg
        assert abs(mula.ra_j2000 - expected_ra) < 0.0001, (
            f"Mula RA: {mula.ra_j2000} vs expected {expected_ra}"
        )
        assert abs(mula.dec_j2000 - expected_dec) < 0.0001, (
            f"Mula Dec: {mula.dec_j2000} vs expected {expected_dec}"
        )

    def test_mula_proper_motion_values(self):
        """Test that Mula proper motion matches Hipparcos 2007 catalog."""
        from libephemeris.planets import STARS

        mula = STARS["MULA"]

        # Hipparcos 2007 proper motion for HIP 85927
        # pm_ra (mu_alpha*): -8.53 mas/yr = -0.00853 arcsec/yr
        # pm_dec (mu_delta): -30.80 mas/yr = -0.03080 arcsec/yr
        expected_pm_ra = -0.00853
        expected_pm_dec = -0.03080

        # Check proper motion matches Hipparcos 2007
        assert abs(mula.pm_ra - expected_pm_ra) < 0.0001, (
            f"Mula pm_ra: {mula.pm_ra} vs expected {expected_pm_ra}"
        )
        assert abs(mula.pm_dec - expected_pm_dec) < 0.0001, (
            f"Mula pm_dec: {mula.pm_dec} vs expected {expected_pm_dec}"
        )

    def test_mula_parallax_and_radial_velocity(self):
        """Test that Mula parallax and radial velocity are set."""
        from libephemeris.planets import STARS

        mula = STARS["MULA"]

        # Hipparcos 2007 parallax: 5.71 mas = 0.00571 arcsec
        # Radial velocity: -3.00 km/s (approaching)
        assert mula.parallax > 0, "Mula should have parallax set"
        assert abs(mula.parallax - 0.00571) < 0.001, (
            f"Mula parallax: {mula.parallax} vs expected 0.00571"
        )
        assert abs(mula.radial_velocity - (-3.00)) < 0.01, (
            f"Mula radial_velocity: {mula.radial_velocity} vs expected -3.00"
        )

    def test_mula_is_distant_star(self):
        """Test that Lambda Scorpii is a distant star with small parallax."""
        from libephemeris.planets import STARS

        mula = STARS["MULA"]

        # Lambda Scorpii has a parallax of ~5.71 mas, corresponding to
        # a distance of about 175 parsecs (570 light years)
        # This is a relatively distant star compared to nearby stars

        # Verify small parallax value
        assert mula.parallax < 0.01, (
            f"Lambda Scorpii should have small parallax (<0.01 arcsec), "
            f"got {mula.parallax}"
        )

        # Distance in parsecs = 1 / parallax (in arcsec)
        distance_pc = 1.0 / mula.parallax
        assert distance_pc > 150, (
            f"Lambda Scorpii should be distant (>150 pc), got {distance_pc:.1f} pc"
        )


class TestMulaStarIdentification:
    """Test that Lambda Scorpii (Shaula) is correctly identified."""

    def test_lambda_scorpii_is_hip_85927(self):
        """Test that Lambda Scorpii (Mula) is correctly identified as HIP 85927.

        Lambda Scorpii (Shaula) is the bright star at the tail of Scorpius,
        and is the primary nakshatra star for Mula.
        """
        from libephemeris.planets import STARS

        mula = STARS["MULA"]

        # Verify coordinates match HIP 85927 (Lambda Scorpii)
        # RA should be around 263.4 deg (17h 33m)
        # Dec should be around -37.1 deg

        assert 263.0 < mula.ra_j2000 < 264.0, (
            f"Mula RA should be ~263.4 deg (HIP 85927), got {mula.ra_j2000}"
        )

        assert -38.0 < mula.dec_j2000 < -36.0, (
            f"Mula Dec should be ~-37.1 deg (HIP 85927), got {mula.dec_j2000}"
        )

    def test_mula_star_not_available_in_gaia(self):
        """Document that Lambda Scorpii is too bright for Gaia DR3.

        Lambda Scorpii (V=1.63) is one of the brightest stars in the sky
        and saturates the Gaia detectors, so no Gaia DR3 data is available.
        We use Hipparcos 2007 reprocessed data instead.
        """
        # This test documents the reason for using Hipparcos instead of Gaia
        # Lambda Scorpii (V=1.63) saturates Gaia detectors
        # The same situation applies to Spica (True Citra)
        pass


if __name__ == "__main__":
    # Run tests with verbose output
    import sys

    pytest.main([__file__, "-v", "--tb=short"] + sys.argv[1:])
