"""
Tests for proper motion propagation from Hipparcos J1991.25 reference epoch.

This module tests the propagate_proper_motion() function which calculates
star positions at any epoch by applying proper motion from a reference epoch.
The function implements the standard astrometric formula including the
cos(dec) correction for RA proper motion.

Tests verify:
- Basic proper motion propagation from J1991.25 to J2000.0
- Forward and backward propagation in time
- High proper motion stars (e.g., Sirius, Barnard's Star)
- Low proper motion stars
- Edge cases (polar declinations, zero proper motion)
- Consistency with catalog data
"""

import math

from libephemeris.fixed_stars import propagate_proper_motion
from libephemeris.constants import (
    J2000,
    J1991_25,
    DAYS_PER_JULIAN_YEAR,
)


class TestProperMotionBasics:
    """Basic tests for propagate_proper_motion function."""

    def test_no_proper_motion_same_position(self):
        """Star with zero proper motion should stay in place."""
        ra_in, dec_in = 100.0, 45.0
        ra_out, dec_out = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=0.0,
            pm_dec=0.0,
            from_jd=J1991_25,
            to_jd=J2000,
        )
        assert abs(ra_out - ra_in) < 1e-10, (
            "RA should be unchanged with zero proper motion"
        )
        assert abs(dec_out - dec_in) < 1e-10, (
            "Dec should be unchanged with zero proper motion"
        )

    def test_same_epoch_same_position(self):
        """Propagating to the same epoch should return the same position."""
        ra_in, dec_in = 180.0, -30.0
        pm_ra, pm_dec = 0.5, -0.3  # arcsec/year

        ra_out, dec_out = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=pm_ra,
            pm_dec=pm_dec,
            from_jd=J2000,
            to_jd=J2000,
        )
        assert abs(ra_out - ra_in) < 1e-10, "RA should be unchanged at same epoch"
        assert abs(dec_out - dec_in) < 1e-10, "Dec should be unchanged at same epoch"

    def test_returns_tuple_of_floats(self):
        """Function should return a tuple of two floats."""
        result = propagate_proper_motion(
            ra_epoch=100.0,
            dec_epoch=45.0,
            pm_ra_cosdec=0.1,
            pm_dec=0.1,
            from_jd=J1991_25,
            to_jd=J2000,
        )
        assert isinstance(result, tuple), "Should return a tuple"
        assert len(result) == 2, "Should return two values"
        assert isinstance(result[0], float), "RA should be float"
        assert isinstance(result[1], float), "Dec should be float"


class TestProperMotionDirection:
    """Test that proper motion moves stars in the expected direction."""

    def test_positive_pm_ra_increases_ra(self):
        """Positive pm_ra_cosdec should increase RA over time."""
        ra_in, dec_in = 100.0, 45.0
        pm_ra, pm_dec = 0.5, 0.0  # arcsec/year, eastward motion only

        ra_out, dec_out = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=pm_ra,
            pm_dec=pm_dec,
            from_jd=J1991_25,
            to_jd=J2000,  # ~8.75 years later
        )
        assert ra_out > ra_in, "RA should increase with positive pm_ra_cosdec"
        assert abs(dec_out - dec_in) < 1e-10, "Dec should be unchanged with zero pm_dec"

    def test_negative_pm_ra_decreases_ra(self):
        """Negative pm_ra_cosdec should decrease RA over time."""
        ra_in, dec_in = 200.0, 45.0
        pm_ra, pm_dec = -0.5, 0.0  # arcsec/year, westward motion

        ra_out, dec_out = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=pm_ra,
            pm_dec=pm_dec,
            from_jd=J1991_25,
            to_jd=J2000,
        )
        assert ra_out < ra_in, "RA should decrease with negative pm_ra_cosdec"

    def test_positive_pm_dec_increases_dec(self):
        """Positive pm_dec should increase Dec (move north) over time."""
        ra_in, dec_in = 100.0, 30.0
        pm_ra, pm_dec = 0.0, 0.5  # arcsec/year, northward motion

        ra_out, dec_out = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=pm_ra,
            pm_dec=pm_dec,
            from_jd=J1991_25,
            to_jd=J2000,
        )
        assert abs(ra_out - ra_in) < 1e-10, "RA should be unchanged with zero pm_ra"
        assert dec_out > dec_in, "Dec should increase with positive pm_dec"

    def test_negative_pm_dec_decreases_dec(self):
        """Negative pm_dec should decrease Dec (move south) over time."""
        ra_in, dec_in = 100.0, 30.0
        pm_ra, pm_dec = 0.0, -0.5  # arcsec/year, southward motion

        ra_out, dec_out = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=pm_ra,
            pm_dec=pm_dec,
            from_jd=J1991_25,
            to_jd=J2000,
        )
        assert dec_out < dec_in, "Dec should decrease with negative pm_dec"

    def test_backward_propagation_reverses_direction(self):
        """Propagating backward in time should reverse the effect."""
        ra_in, dec_in = 150.0, 20.0
        pm_ra, pm_dec = 0.3, 0.2  # arcsec/year

        # Propagate forward to J2000
        ra_j2000, dec_j2000 = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=pm_ra,
            pm_dec=pm_dec,
            from_jd=J1991_25,
            to_jd=J2000,
        )

        # Propagate backward to J1991.25
        ra_back, dec_back = propagate_proper_motion(
            ra_epoch=ra_j2000,
            dec_epoch=dec_j2000,
            pm_ra_cosdec=pm_ra,
            pm_dec=pm_dec,
            from_jd=J2000,
            to_jd=J1991_25,
        )

        # Should get back close to original (small error due to linear approx)
        assert abs(ra_back - ra_in) < 0.001, "RA should return to original value"
        assert abs(dec_back - dec_in) < 0.001, "Dec should return to original value"


class TestProperMotionMagnitude:
    """Test quantitative accuracy of proper motion calculations."""

    def test_proper_motion_magnitude_one_year(self):
        """Test that proper motion magnitude is correct for one year."""
        ra_in, dec_in = 180.0, 0.0  # On celestial equator, cos(dec)=1
        pm_ra_cosdec = 1.0  # 1 arcsec/year
        pm_dec = 1.0  # 1 arcsec/year

        # One Julian year later
        from_jd = J2000
        to_jd = J2000 + DAYS_PER_JULIAN_YEAR

        ra_out, dec_out = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=pm_ra_cosdec,
            pm_dec=pm_dec,
            from_jd=from_jd,
            to_jd=to_jd,
        )

        # Expected change: 1 arcsec = 1/3600 degrees
        expected_change_deg = 1.0 / 3600.0

        assert abs((ra_out - ra_in) - expected_change_deg) < 1e-10, (
            f"RA change should be {expected_change_deg} degrees"
        )
        assert abs((dec_out - dec_in) - expected_change_deg) < 1e-10, (
            f"Dec change should be {expected_change_deg} degrees"
        )

    def test_cos_dec_correction_at_equator(self):
        """At equator, cos(dec)=1, so pm_ra_cosdec equals actual RA change."""
        ra_in, dec_in = 180.0, 0.0
        pm_ra_cosdec = 3600.0  # 3600 arcsec/year = 1 degree/year

        from_jd = J2000
        to_jd = J2000 + DAYS_PER_JULIAN_YEAR

        ra_out, _ = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=pm_ra_cosdec,
            pm_dec=0.0,
            from_jd=from_jd,
            to_jd=to_jd,
        )

        # At equator, the RA change should be exactly 1 degree
        assert abs((ra_out - ra_in) - 1.0) < 1e-10, (
            "RA change at equator should be 1 degree"
        )

    def test_cos_dec_correction_at_60_degrees(self):
        """At dec=60, cos(60)=0.5, so actual RA change is 2x pm_ra_cosdec."""
        ra_in, dec_in = 180.0, 60.0
        pm_ra_cosdec = 3600.0  # 3600 arcsec/year

        from_jd = J2000
        to_jd = J2000 + DAYS_PER_JULIAN_YEAR

        ra_out, _ = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=pm_ra_cosdec,
            pm_dec=0.0,
            from_jd=from_jd,
            to_jd=to_jd,
        )

        # cos(60) = 0.5, so actual RA change = 1 / 0.5 = 2 degrees
        expected_ra_change = 1.0 / math.cos(math.radians(60.0))
        assert abs((ra_out - ra_in) - expected_ra_change) < 1e-8, (
            f"RA change at dec=60 should be {expected_ra_change} degrees"
        )

    def test_j1991_25_to_j2000_time_interval(self):
        """Test that J1991.25 to J2000.0 is approximately 8.75 years."""
        dt_years = (J2000 - J1991_25) / DAYS_PER_JULIAN_YEAR
        expected_years = 8.75
        assert abs(dt_years - expected_years) < 0.001, (
            f"Time from J1991.25 to J2000 should be ~{expected_years} years, got {dt_years}"
        )


class TestHighProperMotionStars:
    """Test with high proper motion stars to verify accuracy."""

    def test_sirius_proper_motion(self):
        """Test with Sirius - one of the highest proper motion bright stars.

        Sirius (HIP 32349) Hipparcos data:
        - pm_ra* = -546.01 mas/yr = -0.54601 arcsec/yr
        - pm_dec = -1223.07 mas/yr = -1.22307 arcsec/yr

        Over 8.75 years (J1991.25 to J2000), Sirius should move:
        - RA: ~-4.78 arcsec / cos(dec) = ~-4.98 arcsec in RA angle
        - Dec: ~-10.7 arcsec
        """
        # Approximate Sirius position at J1991.25 (from Hipparcos)
        # For testing, we use values that let us verify the calculation
        ra_1991 = 101.29  # degrees
        dec_1991 = -16.72  # degrees
        pm_ra_cosdec = -0.54601  # arcsec/yr
        pm_dec = -1.22307  # arcsec/yr

        ra_2000, dec_2000 = propagate_proper_motion(
            ra_epoch=ra_1991,
            dec_epoch=dec_1991,
            pm_ra_cosdec=pm_ra_cosdec,
            pm_dec=pm_dec,
            from_jd=J1991_25,
            to_jd=J2000,
        )

        # Verify Dec change: ~8.75 years * -1.22307 arcsec/yr = ~-10.7 arcsec
        dt_years = (J2000 - J1991_25) / DAYS_PER_JULIAN_YEAR
        expected_dec_change_arcsec = pm_dec * dt_years
        expected_dec_change_deg = expected_dec_change_arcsec / 3600.0
        actual_dec_change = dec_2000 - dec_1991

        assert abs(actual_dec_change - expected_dec_change_deg) < 1e-8, (
            f"Dec change should be {expected_dec_change_deg} deg"
        )

        # Verify RA change with cos(dec) correction
        cos_dec = math.cos(math.radians(dec_1991))
        expected_ra_change_deg = (pm_ra_cosdec * dt_years / 3600.0) / cos_dec
        actual_ra_change = ra_2000 - ra_1991

        assert abs(actual_ra_change - expected_ra_change_deg) < 1e-8, (
            f"RA change should be {expected_ra_change_deg} deg"
        )


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_ra_wraparound_at_360(self):
        """RA should wrap around from 360 to 0."""
        ra_in, dec_in = 359.5, 0.0
        pm_ra_cosdec = 3600.0  # 1 degree/year

        from_jd = J2000
        to_jd = J2000 + DAYS_PER_JULIAN_YEAR

        ra_out, _ = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=pm_ra_cosdec,
            pm_dec=0.0,
            from_jd=from_jd,
            to_jd=to_jd,
        )

        # RA should wrap: 359.5 + 1.0 = 360.5 -> 0.5
        expected_ra = 0.5
        assert abs(ra_out - expected_ra) < 1e-8, (
            f"RA should wrap to {expected_ra}, got {ra_out}"
        )

    def test_ra_wraparound_negative(self):
        """RA going negative should wrap to near 360."""
        ra_in, dec_in = 0.5, 0.0
        pm_ra_cosdec = -3600.0  # -1 degree/year

        from_jd = J2000
        to_jd = J2000 + DAYS_PER_JULIAN_YEAR

        ra_out, _ = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=pm_ra_cosdec,
            pm_dec=0.0,
            from_jd=from_jd,
            to_jd=to_jd,
        )

        # RA should wrap: 0.5 - 1.0 = -0.5 -> 359.5
        expected_ra = 359.5
        assert abs(ra_out - expected_ra) < 1e-8, (
            f"RA should wrap to {expected_ra}, got {ra_out}"
        )

    def test_dec_clamp_at_north_pole(self):
        """Dec should be clamped at +90 degrees."""
        ra_in, dec_in = 100.0, 89.9
        pm_dec = 3600.0  # Very high northward motion

        from_jd = J2000
        to_jd = J2000 + DAYS_PER_JULIAN_YEAR

        _, dec_out = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=0.0,
            pm_dec=pm_dec,
            from_jd=from_jd,
            to_jd=to_jd,
        )

        assert dec_out <= 90.0, f"Dec should not exceed 90, got {dec_out}"

    def test_dec_clamp_at_south_pole(self):
        """Dec should be clamped at -90 degrees."""
        ra_in, dec_in = 100.0, -89.9
        pm_dec = -3600.0  # Very high southward motion

        from_jd = J2000
        to_jd = J2000 + DAYS_PER_JULIAN_YEAR

        _, dec_out = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=0.0,
            pm_dec=pm_dec,
            from_jd=from_jd,
            to_jd=to_jd,
        )

        assert dec_out >= -90.0, f"Dec should not be less than -90, got {dec_out}"

    def test_near_pole_ra_stability(self):
        """Near poles, RA changes should be handled without numerical issues."""
        ra_in, dec_in = 100.0, 89.99  # Very close to north pole
        pm_ra_cosdec = 0.1  # Small but non-zero

        from_jd = J2000
        to_jd = J2000 + DAYS_PER_JULIAN_YEAR

        ra_out, dec_out = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=pm_ra_cosdec,
            pm_dec=0.0,
            from_jd=from_jd,
            to_jd=to_jd,
        )

        # Should not produce NaN or Inf
        assert math.isfinite(ra_out), "RA should be finite near pole"
        assert math.isfinite(dec_out), "Dec should be finite near pole"
        assert 0 <= ra_out < 360, "RA should be in valid range"

    def test_exact_pole_ra_unchanged(self):
        """At exactly the pole, RA should remain unchanged."""
        ra_in, dec_in = 100.0, 90.0
        pm_ra_cosdec = 1.0

        from_jd = J2000
        to_jd = J2000 + DAYS_PER_JULIAN_YEAR

        ra_out, dec_out = propagate_proper_motion(
            ra_epoch=ra_in,
            dec_epoch=dec_in,
            pm_ra_cosdec=pm_ra_cosdec,
            pm_dec=0.0,
            from_jd=from_jd,
            to_jd=to_jd,
        )

        # At exactly 90 dec, RA is undefined and should stay the same
        assert abs(ra_out - ra_in) < 1e-6, "RA should be unchanged at pole"


class TestEpochConstants:
    """Test that epoch constants are correct."""

    def test_j2000_value(self):
        """J2000.0 should be JD 2451545.0."""
        assert J2000 == 2451545.0, f"J2000 should be 2451545.0, got {J2000}"

    def test_j1991_25_value(self):
        """J1991.25 should be JD 2448349.0625."""
        assert J1991_25 == 2448349.0625, (
            f"J1991_25 should be 2448349.0625, got {J1991_25}"
        )

    def test_days_per_julian_year(self):
        """DAYS_PER_JULIAN_YEAR should be exactly 365.25."""
        assert DAYS_PER_JULIAN_YEAR == 365.25, (
            f"DAYS_PER_JULIAN_YEAR should be 365.25, got {DAYS_PER_JULIAN_YEAR}"
        )


class TestCatalogConsistency:
    """Test consistency with catalog data."""

    def test_propagation_from_j2000_to_future(self):
        """Test propagating a known star position forward in time."""
        # Use Regulus as example
        # At J2000: RA = 152.092958, Dec = 11.967208
        # pm_ra* = -249 mas/yr = -0.000249 arcsec/yr (very small)
        # pm_dec = +152 mas/yr = +0.00152 arcsec/yr

        ra_j2000, dec_j2000 = 152.092958, 11.967208
        pm_ra_cosdec = -0.00249  # arcsec/yr
        pm_dec = 0.00152  # arcsec/yr

        # Propagate to year 2100 (100 years forward)
        jd_2100 = J2000 + 100 * DAYS_PER_JULIAN_YEAR

        ra_2100, dec_2100 = propagate_proper_motion(
            ra_epoch=ra_j2000,
            dec_epoch=dec_j2000,
            pm_ra_cosdec=pm_ra_cosdec,
            pm_dec=pm_dec,
            from_jd=J2000,
            to_jd=jd_2100,
        )

        # Over 100 years:
        # Dec change: 100 * 0.00152 = 0.152 arcsec = 0.0000422 degrees
        expected_dec_change = 100 * pm_dec / 3600.0
        actual_dec_change = dec_2100 - dec_j2000

        assert abs(actual_dec_change - expected_dec_change) < 1e-10, (
            f"Dec change for Regulus should be {expected_dec_change}"
        )

    def test_roundtrip_propagation(self):
        """Propagating forward then backward should return to start."""
        ra_start, dec_start = 200.0, 35.0
        pm_ra_cosdec = 0.15
        pm_dec = -0.08

        # Forward 50 years
        jd_future = J2000 + 50 * DAYS_PER_JULIAN_YEAR
        ra_future, dec_future = propagate_proper_motion(
            ra_epoch=ra_start,
            dec_epoch=dec_start,
            pm_ra_cosdec=pm_ra_cosdec,
            pm_dec=pm_dec,
            from_jd=J2000,
            to_jd=jd_future,
        )

        # Backward 50 years
        ra_back, dec_back = propagate_proper_motion(
            ra_epoch=ra_future,
            dec_epoch=dec_future,
            pm_ra_cosdec=pm_ra_cosdec,
            pm_dec=pm_dec,
            from_jd=jd_future,
            to_jd=J2000,
        )

        assert abs(ra_back - ra_start) < 1e-6, "RA should return to start"
        assert abs(dec_back - dec_start) < 1e-6, "Dec should return to start"
