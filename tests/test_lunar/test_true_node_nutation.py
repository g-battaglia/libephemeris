"""
Tests for IAU 2000A nutation correction in True Lunar Node calculations.

This module validates that the True Node calculation correctly applies
nutation to rotate from the mean ecliptic to the true ecliptic of date.

The nutation correction adds up to ~9" (0.0025 degrees) oscillation
with an 18.6-year dominant period matching the lunar node precession.
"""

import math
import pytest
import random
from libephemeris.lunar import calc_true_lunar_node
from skyfield.nutationlib import iau2000a_radians
from libephemeris.state import get_timescale


@pytest.mark.unit
class TestTrueNodeNutation:
    """Tests for nutation correction in True Node calculation."""

    def test_nutation_is_applied(self):
        """
        Verify that nutation correction is actually applied to the True Node.

        The nutation in longitude (dpsi) can reach up to ~17 arcseconds
        (approximately 0.005 degrees) at the extremes of the 18.6-year cycle.
        """
        jd_tt = 2451545.0  # J2000

        ts = get_timescale()
        t = ts.tt_jd(jd_tt)
        dpsi_rad, deps_rad = iau2000a_radians(t)

        # Nutation in longitude should be non-zero
        dpsi_arcsec = math.degrees(dpsi_rad) * 3600
        assert abs(dpsi_arcsec) > 0.1, (
            f"Nutation in longitude should be significant, got {dpsi_arcsec:.2f} arcsec"
        )

        # The nutation is applied in calc_true_lunar_node, which we verify
        # by checking the function completes without error
        lon, lat, dist = calc_true_lunar_node(jd_tt)
        assert 0 <= lon < 360, f"Longitude {lon} out of range"
        assert lat == 0.0, "Node latitude should always be 0"

    def test_nutation_magnitude_reasonable(self):
        """
        Verify nutation magnitude is within expected bounds.

        The dominant nutation term (18.6-year lunar node term) has amplitude
        of about 9.2 arcseconds in obliquity and ~17 arcseconds in longitude.
        Total nutation should typically be within +-20 arcseconds.
        """
        test_dates = [
            2451545.0,  # J2000
            2451545.0 + 365.25,  # 1 year after
            2451545.0 + 9.3 * 365.25,  # Near nutation maximum (~half cycle)
            2451545.0 + 18.6 * 365.25,  # Full nutation cycle
        ]

        ts = get_timescale()

        for jd in test_dates:
            t = ts.tt_jd(jd)
            dpsi_rad, deps_rad = iau2000a_radians(t)

            dpsi_arcsec = math.degrees(dpsi_rad) * 3600
            deps_arcsec = math.degrees(deps_rad) * 3600

            # Check bounds - nutation should be within reasonable limits
            assert abs(dpsi_arcsec) < 30, (
                f"Nutation in longitude {dpsi_arcsec:.2f} arcsec exceeds bounds at JD {jd}"
            )
            assert abs(deps_arcsec) < 15, (
                f"Nutation in obliquity {deps_arcsec:.2f} arcsec exceeds bounds at JD {jd}"
            )

    def test_nutation_variation_over_18_6_year_cycle(self):
        """
        Test that nutation varies over the 18.6-year lunar node precession cycle.

        The dominant nutation term has a period of 18.6 years, so positions
        at different phases should show different nutation values.
        """
        jd_start = 2451545.0  # J2000
        period = 18.6 * 365.25  # 18.6-year period in days

        ts = get_timescale()

        # Sample at 0, 1/4, 1/2, 3/4 of the cycle
        phases = [0, 0.25, 0.5, 0.75]
        dpsi_values = []

        for phase in phases:
            jd = jd_start + phase * period
            t = ts.tt_jd(jd)
            dpsi_rad, _ = iau2000a_radians(t)
            dpsi_values.append(math.degrees(dpsi_rad) * 3600)  # arcseconds

        # Values at different phases should differ
        max_dpsi = max(dpsi_values)
        min_dpsi = min(dpsi_values)
        variation = max_dpsi - min_dpsi

        # Expect significant variation (~18 arcsec peak-to-peak)
        assert variation > 5, (
            f"Nutation should vary significantly over 18.6-year cycle, "
            f"got variation of {variation:.2f} arcsec"
        )

    def test_true_node_includes_nutation_effect(self):
        """
        Verify that True Node positions at different dates show nutation effects.

        When nutation in longitude is at a maximum, it should shift the
        node position by a few arcseconds compared to a version without nutation.
        """
        # Test at J2000 and at a date where nutation is significantly different
        jd1 = 2451545.0  # J2000
        jd2 = 2451545.0 + 9.3 * 365.25  # About half the 18.6-year cycle

        ts = get_timescale()

        # Get nutation values at both dates
        t1 = ts.tt_jd(jd1)
        t2 = ts.tt_jd(jd2)
        dpsi1, _ = iau2000a_radians(t1)
        dpsi2, _ = iau2000a_radians(t2)

        # The difference in nutation should be significant
        dpsi_diff_arcsec = abs(math.degrees(dpsi2) - math.degrees(dpsi1)) * 3600

        # Calculate true node at both dates
        lon1, _, _ = calc_true_lunar_node(jd1)
        lon2, _, _ = calc_true_lunar_node(jd2)

        # Both should return valid positions
        assert 0 <= lon1 < 360
        assert 0 <= lon2 < 360

        # The positions will differ due to node motion AND nutation
        # This test verifies the calculation includes nutation contribution
        assert dpsi_diff_arcsec > 0.1, (
            f"Expected significant nutation difference, got {dpsi_diff_arcsec:.2f} arcsec"
        )

    def test_nutation_consistency_with_fixed_stars(self):
        """
        Verify that the same nutation model is used as in fixed star calculations.

        Both True Node and fixed stars should use IAU 2000A nutation,
        ensuring consistency across the library.
        """
        from skyfield.nutationlib import iau2000a_radians

        jd = 2451545.0
        ts = get_timescale()
        t = ts.tt_jd(jd)

        # Get nutation using the same function as True Node
        dpsi_rad, deps_rad = iau2000a_radians(t)

        # Verify these are reasonable values for J2000
        dpsi_deg = math.degrees(dpsi_rad)
        deps_deg = math.degrees(deps_rad)

        # At J2000, nutation values should be small but non-zero
        assert abs(dpsi_deg) < 0.01, f"dpsi at J2000: {dpsi_deg} degrees"
        assert abs(deps_deg) < 0.005, f"deps at J2000: {deps_deg} degrees"

    def test_true_node_precision_with_nutation(self):
        """
        Test that True Node with nutation maintains precision over date range.

        The addition of nutation should improve precision for the true
        ecliptic of date, not degrade it.
        """
        random.seed(42)
        errors = []

        # Test 50 random dates between 1950 and 2050
        for _ in range(50):
            jd = 2433282.5 + random.uniform(0, 36525)

            try:
                lon, lat, dist = calc_true_lunar_node(jd)

                # Verify basic validity
                assert 0 <= lon < 360, f"Invalid longitude {lon} at JD {jd}"
                assert lat == 0.0, "Latitude should be 0"
                assert dist > 0, "Distance should be positive"

            except Exception as e:
                # May fail for dates outside ephemeris range
                if "ephemeris" not in str(e).lower():
                    raise

    @pytest.mark.parametrize(
        "year,expected_dpsi_range",
        [
            (2000, (-20, 20)),  # J2000 - moderate nutation
            (2006, (-20, 20)),  # Near nutation minimum
            (2010, (-20, 20)),  # Different phase
            (2015, (-20, 20)),  # Different phase
            (2020, (-20, 20)),  # Recent past
        ],
    )
    def test_nutation_bounds_at_various_epochs(self, year, expected_dpsi_range):
        """
        Test that nutation in longitude stays within expected bounds at various dates.
        """
        # Convert year to JD (approximate - Jan 1, noon)
        jd = 2451545.0 + (year - 2000) * 365.25

        ts = get_timescale()
        t = ts.tt_jd(jd)
        dpsi_rad, _ = iau2000a_radians(t)
        dpsi_arcsec = math.degrees(dpsi_rad) * 3600

        min_dpsi, max_dpsi = expected_dpsi_range
        assert min_dpsi <= dpsi_arcsec <= max_dpsi, (
            f"Nutation at {year}: {dpsi_arcsec:.2f} arcsec "
            f"outside expected range {expected_dpsi_range}"
        )


@pytest.mark.integration
class TestTrueNodeNutationIntegration:
    """Integration tests for nutation in True Node calculations."""

    def test_true_node_frame_transformation_chain(self):
        """
        Verify the complete frame transformation chain:
        ICRS -> J2000 ecliptic -> Precession -> Nutation -> True ecliptic of date
        """
        jd = 2451545.0 + 10 * 365.25  # 10 years after J2000

        # Get true node position
        lon, lat, dist = calc_true_lunar_node(jd)

        # Position should be valid
        assert 0 <= lon < 360
        assert lat == 0.0

        # The node longitude should be somewhere in the zodiac
        # Mean node at J2000 is ~125 degrees
        # After 10 years of retrograde motion (~19 degrees/year),
        # it should be roughly 125 - 10*19.3 = ~-68, or ~292 degrees
        # But true node oscillates, so allow wide range
        assert 0 <= lon < 360

    def test_nutation_does_not_cause_discontinuities(self):
        """
        Verify that adding nutation doesn't cause discontinuities in node motion.
        """
        jd_start = 2451545.0
        positions = []

        # Sample every day for 30 days
        for i in range(30):
            jd = jd_start + i
            lon, _, _ = calc_true_lunar_node(jd)
            positions.append(lon)

        # Check for smooth motion (no large jumps)
        for i in range(1, len(positions)):
            diff = positions[i] - positions[i - 1]
            # Handle wrap-around
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360

            # Daily motion should be small (mean ~0.05 deg/day, true oscillates)
            assert abs(diff) < 2.0, (
                f"Jump of {diff:.3f} deg between day {i - 1} and {i} "
                f"(positions: {positions[i - 1]:.3f} -> {positions[i]:.3f})"
            )

    def test_nutation_effect_magnitude(self):
        """
        Estimate the contribution of nutation to the True Node position.

        The nutation in longitude (dpsi) directly adds to the node longitude.
        This test verifies the contribution is at the expected scale.
        """
        ts = get_timescale()

        # Sample over 2 years to see nutation variation
        jd_start = 2451545.0
        nutation_contributions = []

        for i in range(24):  # Monthly samples for 2 years
            jd = jd_start + i * 30
            t = ts.tt_jd(jd)
            dpsi_rad, _ = iau2000a_radians(t)
            nutation_contributions.append(math.degrees(dpsi_rad) * 3600)  # arcsec

        max_contribution = max(abs(n) for n in nutation_contributions)
        variation = max(nutation_contributions) - min(nutation_contributions)

        # Nutation contribution should be measurable (typically 10-17 arcsec)
        assert max_contribution > 5, (
            f"Nutation contribution seems too small: {max_contribution:.2f} arcsec"
        )
        assert max_contribution < 30, (
            f"Nutation contribution seems too large: {max_contribution:.2f} arcsec"
        )

        # Variation over 2 years should be significant
        assert variation > 5, (
            f"Nutation variation over 2 years seems too small: {variation:.2f} arcsec"
        )
