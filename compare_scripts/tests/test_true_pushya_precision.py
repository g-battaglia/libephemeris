"""
Tests for improved True Pushya ayanamsha precision.

These tests verify that the True Pushya ayanamsha calculation using high-precision
Gaia DR3 coordinates for Delta Cancri (HIP 42911) with full proper motion correction
(including parallax and radial velocity) achieves better precision than before.

The improvement targets:
- Previous precision: ~0.06 deg vs Swiss Ephemeris (without parallax/RV)
- Current precision: <0.006 deg vs Swiss Ephemeris (10x improvement)
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import SIDM_TRUE_PUSHYA, SUN, FLG_SIDEREAL


class TestTruePushyaPrecision:
    """Test improved True Pushya ayanamsha precision."""

    @pytest.fixture(autouse=True)
    def reset_sidereal_mode(self):
        """Reset sidereal mode after each test."""
        yield
        swe.set_sid_mode(0)
        ephem.set_sid_mode(0)

    def test_true_pushya_at_j2000(self):
        """Test True Pushya ayanamsha at J2000.0 epoch."""
        jd = 2451545.0  # J2000.0

        # Swiss Ephemeris reference
        swe.set_sid_mode(SIDM_TRUE_PUSHYA)
        aya_swe = swe.get_ayanamsa_ut(jd)

        # LibEphemeris with improved Delta Cancri coordinates
        ephem.set_sid_mode(SIDM_TRUE_PUSHYA)
        aya_lib = ephem.get_ayanamsa_ut(jd)

        diff = abs(aya_swe - aya_lib)

        # Target: <0.006 deg precision (calibrated offset at J2000)
        assert diff < 0.006, (
            f"True Pushya at J2000.0: SWE={aya_swe:.6f} deg, LIB={aya_lib:.6f} deg, "
            f"Diff={diff:.6f} deg (target <0.006 deg)"
        )

    def test_true_pushya_at_multiple_dates(self):
        """Test True Pushya ayanamsha at multiple dates across DE421 range."""
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

            swe.set_sid_mode(SIDM_TRUE_PUSHYA)
            aya_swe = swe.get_ayanamsa_ut(jd)

            ephem.set_sid_mode(SIDM_TRUE_PUSHYA)
            aya_lib = ephem.get_ayanamsa_ut(jd)

            diff = abs(aya_swe - aya_lib)
            max_diff = max(max_diff, diff)

            # Should be <0.006 deg at all dates
            assert diff < 0.006, (
                f"True Pushya at {label}: SWE={aya_swe:.6f}, LIB={aya_lib:.6f}, "
                f"Diff={diff:.6f} deg (target <0.006 deg)"
            )

        print(f"\nMax True Pushya difference across all dates: {max_diff:.6f} deg")

    def test_true_pushya_proper_motion_effect(self):
        """Test that proper motion is correctly applied over time."""
        # Ayanamsha should change smoothly over time due to precession
        # and star proper motion

        jd_1950 = swe.julday(1950, 1, 1, 12.0)
        jd_2050 = swe.julday(2050, 1, 1, 12.0)

        ephem.set_sid_mode(SIDM_TRUE_PUSHYA)
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

    def test_sidereal_sun_position_true_pushya(self):
        """Test sidereal Sun position using True Pushya matches Swiss Ephemeris."""
        jd = swe.julday(2024, 4, 14, 12.0)  # Sun in Aries (tropical)

        # Swiss Ephemeris
        swe.set_sid_mode(SIDM_TRUE_PUSHYA)
        pos_swe, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)
        sun_swe = pos_swe[0]

        # LibEphemeris
        ephem.set_sid_mode(SIDM_TRUE_PUSHYA)
        pos_lib, _ = ephem.calc_ut(jd, SUN, FLG_SIDEREAL)
        sun_lib = pos_lib[0]

        diff = abs(sun_swe - sun_lib)
        if diff > 180:
            diff = 360 - diff

        # Sidereal position should match closely
        # Planet position accuracy + ayanamsha accuracy = total accuracy
        assert diff < 0.02, (
            f"Sidereal Sun (True Pushya): SWE={sun_swe:.6f} deg, "
            f"LIB={sun_lib:.6f} deg, Diff={diff:.6f} deg"
        )


class TestPushyaStarData:
    """Test the improved Pushya (Delta Cancri) star data used for True Pushya."""

    def test_pushya_coordinates_precision(self):
        """Test that Pushya coordinates match Gaia DR3 catalog."""
        from libephemeris.planets import STARS

        pushya = STARS["PUSHYA"]

        # Gaia DR3 values for HIP 42911 (Delta Cancri)
        # RA: 08h 44m 41.0991810454s = 131.1712460977 deg
        # Dec: +18 deg 09' 15.509048595" = 18.1543080691 deg
        expected_ra = 131.1712460977
        expected_dec = 18.1543080691

        # Check position matches Gaia DR3 within 0.0001 deg
        assert abs(pushya.ra_j2000 - expected_ra) < 0.0001, (
            f"Pushya RA: {pushya.ra_j2000} vs expected {expected_ra}"
        )
        assert abs(pushya.dec_j2000 - expected_dec) < 0.0001, (
            f"Pushya Dec: {pushya.dec_j2000} vs expected {expected_dec}"
        )

    def test_pushya_proper_motion_values(self):
        """Test that Pushya proper motion matches Gaia DR3 catalog."""
        from libephemeris.planets import STARS

        pushya = STARS["PUSHYA"]

        # Gaia DR3 proper motion for HIP 42911
        # pm_ra (mu_alpha*): -18.435 mas/yr = -0.018435 arcsec/yr
        # pm_dec (mu_delta): -227.813 mas/yr = -0.227813 arcsec/yr
        expected_pm_ra = -0.018435
        expected_pm_dec = -0.227813

        # Check proper motion matches Gaia DR3
        assert abs(pushya.pm_ra - expected_pm_ra) < 0.0001, (
            f"Pushya pm_ra: {pushya.pm_ra} vs expected {expected_pm_ra}"
        )
        assert abs(pushya.pm_dec - expected_pm_dec) < 0.0001, (
            f"Pushya pm_dec: {pushya.pm_dec} vs expected {expected_pm_dec}"
        )

    def test_pushya_parallax_and_radial_velocity(self):
        """Test that Pushya parallax and radial velocity are set."""
        from libephemeris.planets import STARS

        pushya = STARS["PUSHYA"]

        # Gaia DR3 parallax: 23.8271 mas = 0.0238271 arcsec
        # Radial velocity: +17.14 km/s
        assert pushya.parallax > 0, "Pushya should have parallax set"
        assert abs(pushya.parallax - 0.0238271) < 0.001, (
            f"Pushya parallax: {pushya.parallax} vs expected 0.0238271"
        )
        assert abs(pushya.radial_velocity - 17.14) < 0.01, (
            f"Pushya radial_velocity: {pushya.radial_velocity} vs expected 17.14"
        )

    def test_pushya_has_high_proper_motion(self):
        """Test that Delta Cancri has notably high proper motion in declination."""
        from libephemeris.planets import STARS

        pushya = STARS["PUSHYA"]

        # Delta Cancri has a notable proper motion of -227.813 mas/yr in declination
        # This is ~0.23 arcsec/year, which adds up over centuries
        # This makes proper motion correction particularly important for this star

        # Verify the high proper motion value
        assert abs(pushya.pm_dec) > 0.2, (
            f"Delta Cancri should have high Dec proper motion (>0.2 arcsec/yr), "
            f"got {pushya.pm_dec}"
        )

        # Verify proper motion correction is significant over 100 years
        # 0.227813 arcsec/yr * 100 years = 22.8 arcsec = 0.0063 deg
        pm_100yr = abs(pushya.pm_dec) * 100 / 3600  # convert to degrees
        assert pm_100yr > 0.005, (
            f"Proper motion over 100 years should be significant (>0.005 deg), "
            f"got {pm_100yr:.6f} deg"
        )


class TestPushyaHIPCorrection:
    """Test that Delta Cancri uses the correct HIP number."""

    def test_delta_cancri_is_hip_42911(self):
        """Test that Delta Cancri (Pushya) is correctly identified as HIP 42911.

        The code previously had a comment indicating HIP 43834, which is actually
        rho-2 Cancri (55 Cancri), not Delta Cancri (Asellus Australis).
        Delta Cancri is HIP 42911.
        """
        from libephemeris.planets import STARS

        pushya = STARS["PUSHYA"]

        # Verify coordinates match HIP 42911 (Delta Cancri), not HIP 43834
        # HIP 42911: RA = 08h 44m 41.1s, Dec = +18 deg 09' 15.5"
        # HIP 43834: RA = 08h 55m 39.7s, Dec = +27 deg 55' 38.9" (wrong star!)

        # RA should be around 131.17 deg (08h 44m), not 133.9 deg (08h 55m)
        assert 131.0 < pushya.ra_j2000 < 132.0, (
            f"Pushya RA should be ~131.17 deg (HIP 42911), got {pushya.ra_j2000}"
        )

        # Dec should be around +18.15 deg, not +27.93 deg
        assert 18.0 < pushya.dec_j2000 < 19.0, (
            f"Pushya Dec should be ~18.15 deg (HIP 42911), got {pushya.dec_j2000}"
        )


if __name__ == "__main__":
    # Run tests with verbose output
    import sys

    pytest.main([__file__, "-v", "--tb=short"] + sys.argv[1:])
