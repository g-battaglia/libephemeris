"""
Tests for SIDM_USER custom ayanamsha functionality.

Tests user-defined ayanamsha with custom initial value and precession rate
using set_sid_mode with SIDM_USER option.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import SIDM_USER, SIDM_LAHIRI, FLG_SIDEREAL, SUN


class TestUserDefinedAyanamsha:
    """Test SIDM_USER custom ayanamsha functionality."""

    @pytest.mark.unit
    def test_user_ayanamsha_at_reference_epoch(self):
        """User ayanamsha should return ayan_t0 at the reference epoch t0."""
        jd = 2451545.0  # J2000
        ayan_t0 = 24.0  # degrees

        ephem.set_sid_mode(SIDM_USER, t0=jd, ayan_t0=ayan_t0)
        ayan = ephem.get_ayanamsa_ut(jd)

        # At t0, ayanamsha should be exactly ayan_t0
        assert abs(ayan - ayan_t0) < 0.001, f"Expected {ayan_t0}, got {ayan}"

    @pytest.mark.unit
    def test_user_ayanamsha_precession_forward(self):
        """User ayanamsha should precess forward in time."""
        t0 = 2451545.0  # J2000
        ayan_t0 = 24.0  # degrees

        ephem.set_sid_mode(SIDM_USER, t0=t0, ayan_t0=ayan_t0)

        # One century later
        jd_future = t0 + 36525.0  # J2100
        ayan = ephem.get_ayanamsa_ut(jd_future)

        # Precession rate ~5027.8 arcsec/century = ~1.3966 deg/century
        expected = ayan_t0 + (5027.8 / 3600.0)  # ~25.396 degrees
        assert abs(ayan - expected) < 0.01, f"Expected ~{expected}, got {ayan}"

    @pytest.mark.unit
    def test_user_ayanamsha_precession_backward(self):
        """User ayanamsha should precess backward in time."""
        t0 = 2451545.0  # J2000
        ayan_t0 = 24.0  # degrees

        ephem.set_sid_mode(SIDM_USER, t0=t0, ayan_t0=ayan_t0)

        # One century earlier
        jd_past = t0 - 36525.0  # J1900
        ayan = ephem.get_ayanamsa_ut(jd_past)

        # Precession rate ~5027.8 arcsec/century = ~1.3966 deg/century
        expected = ayan_t0 - (5027.8 / 3600.0)  # ~22.603 degrees
        assert abs(ayan - expected) < 0.01, f"Expected ~{expected}, got {ayan}"

    @pytest.mark.unit
    def test_user_ayanamsha_different_t0(self):
        """User ayanamsha should work with different reference epochs."""
        # Use J1900 as reference
        t0 = 2415020.0  # J1900
        ayan_t0 = 22.5  # degrees

        ephem.set_sid_mode(SIDM_USER, t0=t0, ayan_t0=ayan_t0)

        # At J2000 (100 years later)
        jd = 2451545.0
        ayan = ephem.get_ayanamsa_ut(jd)

        # Should have precessed by ~1.4 degrees
        expected = ayan_t0 + (5027.8 / 3600.0)  # ~23.9 degrees
        assert abs(ayan - expected) < 0.05, f"Expected ~{expected}, got {ayan}"

    @pytest.mark.unit
    def test_user_ayanamsha_zero_value(self):
        """User ayanamsha of 0 degrees should work correctly."""
        t0 = 2451545.0  # J2000
        ayan_t0 = 0.0  # No ayanamsha at reference

        ephem.set_sid_mode(SIDM_USER, t0=t0, ayan_t0=ayan_t0)
        ayan = ephem.get_ayanamsa_ut(t0)

        assert abs(ayan) < 0.001, f"Expected 0, got {ayan}"

    @pytest.mark.unit
    def test_user_ayanamsha_large_value(self):
        """User ayanamsha with large value should work correctly."""
        t0 = 2451545.0  # J2000
        ayan_t0 = 359.0  # Near wrap-around

        ephem.set_sid_mode(SIDM_USER, t0=t0, ayan_t0=ayan_t0)
        ayan = ephem.get_ayanamsa_ut(t0)

        assert abs(ayan - ayan_t0) < 0.001, f"Expected {ayan_t0}, got {ayan}"

    @pytest.mark.unit
    def test_user_ayanamsha_name(self):
        """get_ayanamsa_name should return 'User Defined' for SIDM_USER."""
        name = ephem.get_ayanamsa_name(SIDM_USER)
        assert name == "User Defined", f"Expected 'User Defined', got '{name}'"

    @pytest.mark.comparison
    def test_user_ayanamsha_vs_pyswisseph(self):
        """User-defined ayanamsha should match pyswisseph for same parameters."""
        t0 = 2451545.0  # J2000
        ayan_t0 = 24.0

        ephem.set_sid_mode(SIDM_USER, t0=t0, ayan_t0=ayan_t0)
        swe.set_sid_mode(SIDM_USER, t0, ayan_t0)

        # Test at several dates
        test_dates = [
            2451545.0,  # J2000
            2451545.0 + 3652.5,  # 10 years later
            2451545.0 - 3652.5,  # 10 years earlier
            2451545.0 + 36525.0,  # 100 years later
        ]

        for jd in test_dates:
            ayan_lib = ephem.get_ayanamsa_ut(jd)
            ayan_swe = swe.get_ayanamsa_ut(jd)
            diff = abs(ayan_lib - ayan_swe)

            assert diff < 0.01, (
                f"At JD {jd}: lib={ayan_lib}, swe={ayan_swe}, diff={diff}"
            )


class TestUserAyanamshaWithSiderealPositions:
    """Test that user-defined ayanamsha correctly affects sidereal positions."""

    @pytest.mark.unit
    def test_sidereal_position_uses_user_ayanamsha(self):
        """Sidereal position should be tropical - user ayanamsha."""
        jd = 2451545.0
        ayan_t0 = 24.0

        # Get tropical position
        ephem.set_sid_mode(SIDM_LAHIRI)  # Reset first
        pos_tropical, _ = ephem.calc_ut(jd, SUN, 0)

        # Set user ayanamsha
        ephem.set_sid_mode(SIDM_USER, t0=jd, ayan_t0=ayan_t0)

        # Get sidereal position
        pos_sidereal, _ = ephem.calc_ut(jd, SUN, FLG_SIDEREAL)

        # Sidereal = Tropical - Ayanamsha
        expected_sidereal = (pos_tropical[0] - ayan_t0) % 360.0
        diff = abs(pos_sidereal[0] - expected_sidereal)
        # Handle wrap-around
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Expected ~{expected_sidereal}, got {pos_sidereal[0]}"

    @pytest.mark.comparison
    def test_sidereal_position_vs_pyswisseph(self):
        """Sidereal position with user ayanamsha should match pyswisseph."""
        jd = 2451545.0
        t0 = jd
        ayan_t0 = 23.5

        ephem.set_sid_mode(SIDM_USER, t0=t0, ayan_t0=ayan_t0)
        swe.set_sid_mode(SIDM_USER, t0, ayan_t0)

        pos_lib, _ = ephem.calc_ut(jd, SUN, FLG_SIDEREAL)
        pos_swe = swe.calc_ut(jd, SUN, FLG_SIDEREAL)

        # Compare longitude
        diff = abs(pos_lib[0] - pos_swe[0][0])
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Longitude diff {diff} >= 0.01"


class TestUserAyanamshaEdgeCases:
    """Test edge cases for user-defined ayanamsha."""

    @pytest.mark.unit
    def test_switching_from_user_to_standard(self):
        """Switching from user to standard ayanamsha should work."""
        jd = 2451545.0

        # Set user ayanamsha
        ephem.set_sid_mode(SIDM_USER, t0=jd, ayan_t0=50.0)
        ayan_user = ephem.get_ayanamsa_ut(jd)
        assert abs(ayan_user - 50.0) < 0.01

        # Switch to Lahiri
        ephem.set_sid_mode(SIDM_LAHIRI)
        ayan_lahiri = ephem.get_ayanamsa_ut(jd)

        # Should be Lahiri value (~23.9), not user value
        assert 23.0 < ayan_lahiri < 25.0, f"Expected Lahiri value, got {ayan_lahiri}"
        assert abs(ayan_lahiri - 50.0) > 20, "Still using user ayanamsha"

    @pytest.mark.unit
    def test_default_t0_uses_j2000(self):
        """When t0=0 is passed, should default to J2000."""
        jd = 2451545.0  # J2000
        ayan_t0 = 24.0

        # Pass t0=0 (should default to J2000)
        ephem.set_sid_mode(SIDM_USER, t0=0, ayan_t0=ayan_t0)
        ayan = ephem.get_ayanamsa_ut(jd)

        # At J2000, should get ayan_t0
        assert abs(ayan - ayan_t0) < 0.001, f"Expected {ayan_t0}, got {ayan}"

    @pytest.mark.unit
    def test_negative_ayanamsha(self):
        """Negative ayanamsha values should work (wraps to 360-x)."""
        jd = 2451545.0
        ayan_t0 = -1.0  # Should be equivalent to 359.0

        ephem.set_sid_mode(SIDM_USER, t0=jd, ayan_t0=ayan_t0)
        ayan = ephem.get_ayanamsa_ut(jd)

        # Result should be normalized to 0-360
        expected = 359.0
        assert abs(ayan - expected) < 0.01, f"Expected ~{expected}, got {ayan}"

    @pytest.mark.unit
    def test_mimicking_lahiri(self):
        """User ayanamsha can mimic Lahiri by using same parameters."""
        jd = 2451545.0

        # Get Lahiri value at J2000
        ephem.set_sid_mode(SIDM_LAHIRI)
        lahiri_at_j2000 = ephem.get_ayanamsa_ut(jd)

        # Set user ayanamsha to mimic Lahiri
        ephem.set_sid_mode(SIDM_USER, t0=jd, ayan_t0=lahiri_at_j2000)

        # Test at a different date
        test_jd = jd + 3652.5  # 10 years later

        ephem.set_sid_mode(SIDM_LAHIRI)
        lahiri_later = ephem.get_ayanamsa_ut(test_jd)

        ephem.set_sid_mode(SIDM_USER, t0=jd, ayan_t0=lahiri_at_j2000)
        user_later = ephem.get_ayanamsa_ut(test_jd)

        # Should be very close (both using same precession rate)
        diff = abs(lahiri_later - user_later)
        assert diff < 0.01, f"Expected values to match, diff={diff}"
