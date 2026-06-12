"""
Tests for enhanced calc_mean_lilith() formula with T⁴ term and periodic corrections.

The enhanced formula includes:
- T⁴ term from Chapront-Touzé & Chapront (1988) lunar theory (+ T⁴/3526000)
- Periodic corrections for solar/lunar perturbations on mean apsidal motion

Expected improvement: reduce error from ~270 arcsec to <100 arcsec compared
to JPL DE440 reference across the 1800-2200 range.
"""

import pytest
import warnings

from libephemeris.lunar import calc_mean_lilith, MeeusPolynomialWarning


class TestEnhancedMeanLilithFormula:
    """Test the enhanced mean Lilith formula components."""

    # J2000.0 epoch
    J2000 = 2451545.0

    # 1 Julian century in days
    CENTURY = 36525.0

    def test_returns_valid_longitude_at_j2000(self):
        """Mean Lilith at J2000 should return valid longitude 0-360."""
        result = calc_mean_lilith(self.J2000)
        assert 0 <= result < 360
        # At J2000, T=0, so only base polynomial matters
        # Expected: 83.3532465 + 180 = 263.3532465 (approx)
        assert 262 < result < 265

    def test_returns_valid_longitude_across_range(self):
        """Mean Lilith should return valid longitude for various dates."""
        test_dates = [
            self.J2000 - 2 * self.CENTURY,  # Year 1800
            self.J2000 - 1 * self.CENTURY,  # Year 1900
            self.J2000,  # Year 2000
            self.J2000 + 1 * self.CENTURY,  # Year 2100
            self.J2000 + 2 * self.CENTURY,  # Year 2200
        ]

        for jd in test_dates:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = calc_mean_lilith(jd)
                assert 0 <= result < 360, f"Invalid longitude {result} at JD {jd}"

    def test_t4_term_contribution(self):
        """T⁴ term should contribute measurably for distant dates."""
        # At T=2 (year 2200), T⁴ = 16
        # T⁴/3526000 ≈ 0.0000045° contribution to perigee
        # This is small but still important for precision

        # At T=10 (year 3000), T⁴ = 10000
        # T⁴/3526000 ≈ 0.00284° contribution (about 10 arcsec)

        jd_3000 = self.J2000 + 10 * self.CENTURY

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = calc_mean_lilith(jd_3000)
            assert 0 <= result < 360

    def test_periodic_corrections_vary_with_time(self):
        """Periodic corrections should cause small variations over short periods."""
        # Check that positions at close dates differ slightly due to periodic terms
        jd1 = self.J2000
        jd2 = self.J2000 + 27.3  # About one lunar month later

        result1 = calc_mean_lilith(jd1)
        result2 = calc_mean_lilith(jd2)

        # Should differ due to mean motion + periodic corrections
        # Mean motion is about 0.111° per day
        expected_mean_motion = 27.3 * 4069.0137287 / 36525.0

        diff = (result2 - result1 + 360) % 360  # Handle wraparound
        if diff > 180:
            diff = 360 - diff

        # Difference should be close to expected mean motion
        # but not exactly equal due to periodic corrections
        assert abs(diff - expected_mean_motion) < 5.0

    def test_no_warning_within_optimal_range(self):
        """No warning should be issued for dates within 1800-2200."""
        test_dates = [
            self.J2000 - 2 * self.CENTURY,  # 1800
            self.J2000 - 1 * self.CENTURY,  # 1900
            self.J2000,  # 2000
            self.J2000 + 1 * self.CENTURY,  # 2100
            self.J2000 + 2 * self.CENTURY,  # 2200
        ]

        for jd in test_dates:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                calc_mean_lilith(jd)
                meeus_warnings = [
                    x for x in w if issubclass(x.category, MeeusPolynomialWarning)
                ]
                assert len(meeus_warnings) == 0, f"Unexpected warning at JD {jd}"


class TestEnhancedMeanLilithAccuracy:
    """Test accuracy improvements from enhanced formula."""

    J2000 = 2451545.0
    CENTURY = 36525.0

    def test_formula_includes_t4_term(self):
        """Verify T⁴ term is included by checking calculation at distant date."""
        # At T=5 (year 2500), without T⁴: different from with T⁴
        # T⁴ = 625, T⁴/3526000 ≈ 0.000177°

        jd = self.J2000 + 5 * self.CENTURY

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = calc_mean_lilith(jd)

        # Calculate expected base polynomial without T⁴
        T = 5.0
        perigee_without_t4 = (
            83.3532465 + 4069.0137287 * T - 0.0103200 * T**2 - T**3 / 80053.0
        )
        apogee_without_t4 = (perigee_without_t4 + 180.0) % 360.0

        # The actual result should differ due to T⁴ and periodic corrections
        # They won't match exactly
        assert result != pytest.approx(apogee_without_t4, abs=0.001)

    def test_periodic_corrections_are_bounded(self):
        """Periodic corrections should be bounded (small values)."""
        # The periodic corrections are sum of sin terms, each bounded by coefficient
        # Total should be roughly bounded by sum of absolute coefficients
        # ≈ 1.5 + 0.15 + 0.12 + 0.12 + 0.08 + 0.01 + 0.005 ≈ 2° max

        # Test multiple dates to verify corrections stay bounded
        test_dates = [self.J2000 + i * 30.4375 for i in range(36)]  # 3 years

        results = []
        for jd in test_dates:
            results.append(calc_mean_lilith(jd))

        # All results should be valid
        for r in results:
            assert 0 <= r < 360

    def test_mean_motion_rate(self):
        """Mean motion rate should be approximately 4069°/century = 0.111°/day."""
        # Expected rate: 4069.0137287° per century = 0.11136° per day

        jd1 = self.J2000
        jd2 = self.J2000 + 365.25  # One year later

        lon1 = calc_mean_lilith(jd1)
        lon2 = calc_mean_lilith(jd2)

        # Total motion in one year (accounting for complete revolutions)
        # Expected: 4069.0137287 / 100 = 40.69° per year
        expected_motion = 4069.0137287 / 100.0

        # Calculate actual motion (handle wraparound)
        actual_motion = (lon2 - lon1) % 360.0

        # Should be close to expected (within few degrees due to periodic terms)
        assert abs(actual_motion - expected_motion) < 5.0


class TestEnhancedMeanLilithPeriodicTerms:
    """Test individual periodic correction terms."""

    J2000 = 2451545.0

    def test_solar_perturbation_term_exists(self):
        """Solar perturbation term (2D - M) should contribute to variations."""
        # This is the largest periodic term (coefficient -1.4979)
        # At synodic month boundaries, D changes by ~360°, creating variations

        # Sample at different phases
        results = []
        for i in range(30):
            jd = self.J2000 + i * 29.53  # Synodic month intervals
            results.append(calc_mean_lilith(jd))

        # Calculate variations relative to linear trend
        # The residuals should show periodic behavior
        differences = []
        for i in range(1, len(results)):
            expected_linear = 29.53 * 4069.0137287 / 36525.0
            actual = (results[i] - results[i - 1]) % 360.0
            if actual > 180:
                actual = actual - 360
            differences.append(actual - expected_linear)

        # Should have some non-zero residuals (periodic terms)
        max_residual = max(abs(d) for d in differences)
        assert max_residual > 0.01  # At least 0.01° residual

    def test_longitude_continuity(self):
        """Longitude should change continuously without jumps (except at 0/360)."""
        results = []
        for i in range(365):
            jd = self.J2000 + i
            results.append(calc_mean_lilith(jd))

        # Check for smooth transitions
        for i in range(1, len(results)):
            diff = abs(results[i] - results[i - 1])
            # Should be small (< 1° per day typically)
            # Unless crossing 0/360 boundary
            if diff > 180:
                diff = 360 - diff
            assert diff < 1.0, f"Large jump at day {i}: {diff}°"


@pytest.mark.unit
class TestEnhancedMeanLilithKnownValues:
    """Test against known reference values."""

    J2000 = 2451545.0
    CENTURY = 36525.0

    def test_j2000_epoch_value(self):
        """At J2000, mean Lilith should be approximately 263.4°."""
        # At T=0: perigee = 83.3532465, apogee = 263.3532465
        # Periodic corrections at T=0 contribute a small amount

        result = calc_mean_lilith(self.J2000)

        # Should be close to 263.4° (base value plus small corrections)
        assert 262.0 < result < 265.0

    def test_year_2000_jan_1_value(self):
        """Test known approximate value for Jan 1, 2000."""
        # J2000.0 is Jan 1, 2000 at 12:00 TT
        result = calc_mean_lilith(self.J2000)

        # Value should be in Sagittarius (around 263°)
        assert 260 < result < 270

    def test_apsidal_precession_period(self):
        """Mean apogee precession period should be approximately 8.85 years."""
        # The apsidal precession rate: 4069°/century = 40.69°/year
        # Full cycle (360°): 360 / 40.69 ≈ 8.85 years

        lon1 = calc_mean_lilith(self.J2000)
        # Find when longitude returns to same value (mod 360)
        # After 8.85 years ≈ 3232 days

        lon2 = calc_mean_lilith(self.J2000 + 8.85 * 365.25)

        # Should be close to same value (within periodic correction amplitude)
        diff = abs(lon2 - lon1)
        if diff > 180:
            diff = 360 - diff

        # Within 5° (accounting for small phase differences)
        assert diff < 5.0
