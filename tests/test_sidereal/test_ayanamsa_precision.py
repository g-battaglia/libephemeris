"""
Tests for improved ayanamsa precision using IAU 2006 precession model.

This test verifies that the precession rate used in ayanamsa calculations
matches the reference ephemeris closely (using IAU 2006 precession model with
quadratic term: 5028.796273"/century linear + 1.105608"/century² quadratic).

The expected precision is <0.01" (arcseconds) at all epochs from 1900-2100.
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
    SIDM_RAMAN,
    SIDM_KRISHNAMURTI,
)


class TestAyanamsaPrecession:
    """Test that ayanamsa precession rate matches the reference ephemeris."""

    # Test epochs (Julian Day)
    JD_J1900 = 2415020.0  # 1900-01-01 12:00 TT
    JD_J2000 = 2451545.0  # 2000-01-01 12:00 TT
    JD_J2100 = 2488070.0  # 2100-01-01 12:00 TT

    # Reference ephemeris values (computed with the reference ephemeris)
    # Format: (JD, sid_mode, expected_ayanamsa)
    REFERENCE_VALUES = [
        # Lahiri
        (JD_J1900, SIDM_LAHIRI, 22.46051148),
        (JD_J2000, SIDM_LAHIRI, 23.85709235),
        (JD_J2100, SIDM_LAHIRI, 25.25428740),
        # Fagan-Bradley
        (JD_J1900, SIDM_FAGAN_BRADLEY, 23.34371910),
        (JD_J2000, SIDM_FAGAN_BRADLEY, 24.74029999),
        (JD_J2100, SIDM_FAGAN_BRADLEY, 26.13749505),
        # Raman - actual reference ephemeris values
        (JD_J1900, SIDM_RAMAN, 21.01421001),
        (JD_J2000, SIDM_RAMAN, 22.41079104),
        (JD_J2100, SIDM_RAMAN, 23.80798618),
        # Krishnamurti - actual reference ephemeris values
        (JD_J1900, SIDM_KRISHNAMURTI, 22.36365901),
        (JD_J2000, SIDM_KRISHNAMURTI, 23.76024004),
        (JD_J2100, SIDM_KRISHNAMURTI, 25.15743518),
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("jd,sid_mode,expected", REFERENCE_VALUES)
    def test_ayanamsa_matches_reference(self, jd, sid_mode, expected):
        """
        Ayanamsa should match the reference ephemeris within 0.01 arcseconds.

        This tests the IAU 2006 precession model implementation with both
        linear rate (5028.796273"/century) and quadratic term (1.105608"/century²).
        """
        ephem.set_sid_mode(sid_mode)
        ayan = ephem.get_ayanamsa_ut(jd)

        # Tolerance: 0.01 arcseconds = 0.01/3600 degrees ≈ 0.0000028 degrees
        tolerance_arcsec = 0.01
        tolerance_deg = tolerance_arcsec / 3600

        diff = abs(ayan - expected)
        diff_arcsec = diff * 3600

        assert diff < tolerance_deg, (
            f"Ayanamsa mismatch at JD {jd}: "
            f"got {ayan:.8f}°, expected {expected:.8f}°, "
            f'diff {diff_arcsec:.4f}" exceeds {tolerance_arcsec}" tolerance'
        )

    @pytest.mark.unit
    def test_precession_rate_over_century(self):
        """
        Precession rate should be approximately 1.397°/century at J2000.

        The IAU 2006 precession model gives ~5028.8"/century at J2000,
        which is 1.39689°/century.
        """
        ephem.set_sid_mode(SIDM_LAHIRI)

        ayan_2000 = ephem.get_ayanamsa_ut(self.JD_J2000)
        ayan_2100 = ephem.get_ayanamsa_ut(self.JD_J2100)

        rate = ayan_2100 - ayan_2000  # degrees per century

        # Expected rate: ~1.397°/century
        expected_rate = 1.397
        tolerance = 0.002  # Allow 0.002° tolerance

        assert abs(rate - expected_rate) < tolerance, (
            f"Precession rate {rate:.4f}°/century differs from "
            f"expected {expected_rate:.4f}°/century by more than {tolerance}°"
        )

    @pytest.mark.unit
    def test_quadratic_term_effect(self):
        """
        The quadratic precession term should cause rate to increase over time.

        From 1900-2000, rate should be slightly less than from 2000-2100.
        """
        ephem.set_sid_mode(SIDM_LAHIRI)

        ayan_1900 = ephem.get_ayanamsa_ut(self.JD_J1900)
        ayan_2000 = ephem.get_ayanamsa_ut(self.JD_J2000)
        ayan_2100 = ephem.get_ayanamsa_ut(self.JD_J2100)

        rate_1900_2000 = ayan_2000 - ayan_1900
        rate_2000_2100 = ayan_2100 - ayan_2000

        # The rate should increase due to quadratic term
        # Expected increase: ~2.2"/century ≈ 0.0006°/century
        rate_increase = rate_2000_2100 - rate_1900_2000

        # Rate should increase (positive quadratic term)
        assert rate_increase > 0, (
            f"Precession rate should increase over time "
            f"(got rate_1900_2000={rate_1900_2000:.6f}, "
            f"rate_2000_2100={rate_2000_2100:.6f})"
        )

        # The increase should be approximately 0.0006°/century
        expected_increase = 0.0006
        assert abs(rate_increase - expected_increase) < 0.0002, (
            f"Rate increase {rate_increase:.6f}° differs from "
            f"expected ~{expected_increase:.6f}°"
        )
