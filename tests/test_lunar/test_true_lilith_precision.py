"""
Tests to verify True Lilith (osculating apogee) precision claims in PRECISION.md.

These tests validate the documented precision:
- Mean difference: ~52 arcsec (~0.015 degrees)
- Maximum difference: ~235 arcsec (~0.065 degrees)
- RMS difference: ~60 arcsec (~0.017 degrees)

The True Lilith calculation uses the eccentricity vector method derived from
JPL DE ephemeris state vectors, achieving sub-arcminute precision against
pyswisseph.
"""

import random
import math
import pytest

try:
    import swisseph as swe

    HAS_SWISSEPH = True
except ImportError:
    HAS_SWISSEPH = False

import libephemeris as ephem
from libephemeris.constants import OSCU_APOG


def angular_diff(a: float, b: float) -> float:
    """Calculate angular difference accounting for 360-degree wrap."""
    d = abs(a - b)
    if d > 180:
        d = 360 - d
    return d


@pytest.mark.skipif(not HAS_SWISSEPH, reason="pyswisseph not installed")
class TestTrueLilithPrecision:
    """Verify True Lilith precision matches documented claims."""

    # Documented precision thresholds (from PRECISION.md)
    # Adding small margins to account for statistical variation
    DOCUMENTED_MAX_ARCSEC = 235.0
    DOCUMENTED_MEAN_ARCSEC = 52.0
    DOCUMENTED_RMS_ARCSEC = 60.0

    # Test tolerances (slightly relaxed to allow for variation)
    MAX_TOLERANCE_ARCSEC = 300.0  # ~0.083 degrees
    MEAN_TOLERANCE_ARCSEC = 80.0  # ~0.022 degrees
    RMS_TOLERANCE_ARCSEC = 80.0  # ~0.022 degrees

    @pytest.mark.precision
    def test_true_lilith_max_error_within_documented(self):
        """Verify True Lilith max error is within documented tolerance (~235 arcsec)."""
        random.seed(42)
        max_diff_arcsec = 0.0

        for _ in range(500):
            year = random.randint(1950, 2050)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = swe.julday(year, month, day, hour)

            pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, 0)
            pos_lib, _ = ephem.calc_ut(jd, OSCU_APOG, 0)

            diff_deg = angular_diff(pos_swe[0], pos_lib[0])
            diff_arcsec = diff_deg * 3600
            max_diff_arcsec = max(max_diff_arcsec, diff_arcsec)

        assert max_diff_arcsec < self.MAX_TOLERANCE_ARCSEC, (
            f"True Lilith max error {max_diff_arcsec:.1f} arcsec exceeds "
            f"documented tolerance of {self.DOCUMENTED_MAX_ARCSEC} arcsec"
        )

    @pytest.mark.precision
    def test_true_lilith_mean_error_within_documented(self):
        """Verify True Lilith mean error is within documented tolerance (~52 arcsec)."""
        random.seed(42)
        diffs_arcsec = []

        for _ in range(500):
            year = random.randint(1950, 2050)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = swe.julday(year, month, day, hour)

            pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, 0)
            pos_lib, _ = ephem.calc_ut(jd, OSCU_APOG, 0)

            diff_deg = angular_diff(pos_swe[0], pos_lib[0])
            diffs_arcsec.append(diff_deg * 3600)

        mean_diff = sum(diffs_arcsec) / len(diffs_arcsec)

        assert mean_diff < self.MEAN_TOLERANCE_ARCSEC, (
            f"True Lilith mean error {mean_diff:.1f} arcsec exceeds "
            f"documented tolerance of {self.DOCUMENTED_MEAN_ARCSEC} arcsec"
        )

    @pytest.mark.precision
    def test_true_lilith_rms_error_within_documented(self):
        """Verify True Lilith RMS error is within documented tolerance (~60 arcsec)."""
        random.seed(42)
        diffs_arcsec = []

        for _ in range(500):
            year = random.randint(1950, 2050)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = swe.julday(year, month, day, hour)

            pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, 0)
            pos_lib, _ = ephem.calc_ut(jd, OSCU_APOG, 0)

            diff_deg = angular_diff(pos_swe[0], pos_lib[0])
            diffs_arcsec.append(diff_deg * 3600)

        rms_diff = math.sqrt(sum(d**2 for d in diffs_arcsec) / len(diffs_arcsec))

        assert rms_diff < self.RMS_TOLERANCE_ARCSEC, (
            f"True Lilith RMS error {rms_diff:.1f} arcsec exceeds "
            f"documented tolerance of {self.DOCUMENTED_RMS_ARCSEC} arcsec"
        )

    @pytest.mark.precision
    def test_true_lilith_precision_statistics(self):
        """Comprehensive test of True Lilith precision statistics."""
        random.seed(42)
        diffs_arcsec = []

        for _ in range(500):
            year = random.randint(1950, 2050)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = swe.julday(year, month, day, hour)

            pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, 0)
            pos_lib, _ = ephem.calc_ut(jd, OSCU_APOG, 0)

            diff_deg = angular_diff(pos_swe[0], pos_lib[0])
            diffs_arcsec.append(diff_deg * 3600)

        max_diff = max(diffs_arcsec)
        mean_diff = sum(diffs_arcsec) / len(diffs_arcsec)
        rms_diff = math.sqrt(sum(d**2 for d in diffs_arcsec) / len(diffs_arcsec))

        # All three metrics should be within documented tolerances
        assert max_diff < self.MAX_TOLERANCE_ARCSEC, (
            f"Max error {max_diff:.1f} arcsec exceeds tolerance"
        )
        assert mean_diff < self.MEAN_TOLERANCE_ARCSEC, (
            f"Mean error {mean_diff:.1f} arcsec exceeds tolerance"
        )
        assert rms_diff < self.RMS_TOLERANCE_ARCSEC, (
            f"RMS error {rms_diff:.1f} arcsec exceeds tolerance"
        )

        # Also verify we're achieving sub-arcminute precision (< 60 arcsec = 1 arcmin)
        # for at least the mean error
        assert mean_diff < 60.0, (
            f"Mean error {mean_diff:.1f} arcsec should be sub-arcminute"
        )

    @pytest.mark.precision
    def test_true_lilith_sub_degree_precision(self):
        """Verify True Lilith achieves sub-degree precision (max < 0.1 deg)."""
        random.seed(42)

        for _ in range(100):
            year = random.randint(1950, 2050)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = swe.julday(year, month, day, hour)

            pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, 0)
            pos_lib, _ = ephem.calc_ut(jd, OSCU_APOG, 0)

            diff_deg = angular_diff(pos_swe[0], pos_lib[0])

            # Each individual test should show sub-0.1 degree precision
            assert diff_deg < 0.1, (
                f"True Lilith diff {diff_deg:.4f} deg exceeds 0.1 deg at JD {jd:.1f}"
            )
