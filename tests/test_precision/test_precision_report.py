"""
Tests for the precision report generator.

Validates the generate_precision_report.py script functionality including:
- PrecisionStats calculations (max, mean, std dev)
- Report generation for all calculation types
- JSON and CSV output formats
- Tolerance checking
"""

import sys
from pathlib import Path
import json
import pytest

# Add compare_scripts to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "compare_scripts"))

from generate_precision_report import (
    PrecisionStats,
    PrecisionReport,
    angular_diff,
    generate_random_jds,
    measure_planetary_precision,
    measure_velocity_precision,
    measure_house_precision,
    measure_ayanamsha_precision,
    measure_lunar_node_precision,
    measure_heliocentric_precision,
    measure_time_precision,
    generate_report,
)


# =============================================================================
# TEST PRECISION STATS
# =============================================================================


class TestPrecisionStats:
    """Test the PrecisionStats dataclass functionality."""

    def test_add_diff_basic(self):
        """Test adding differences to stats."""
        stats = PrecisionStats(
            calculation_type="Test",
            component="Test Component",
            unit="arcsec",
            tolerance=1.0,
        )

        stats.add_diff(0.1)
        stats.add_diff(0.2)
        stats.add_diff(0.3)

        assert stats.n_samples == 3
        assert len(stats._diffs) == 3

    def test_calculate_mean(self):
        """Test mean calculation."""
        stats = PrecisionStats(
            calculation_type="Test",
            component="Test Component",
            unit="arcsec",
            tolerance=1.0,
        )

        stats.add_diff(1.0)
        stats.add_diff(2.0)
        stats.add_diff(3.0)
        stats.calculate()

        assert stats.mean_diff == 2.0

    def test_calculate_max_min(self):
        """Test max and min calculation."""
        stats = PrecisionStats(
            calculation_type="Test",
            component="Test Component",
            unit="arcsec",
            tolerance=10.0,
        )

        stats.add_diff(1.0)
        stats.add_diff(5.0)
        stats.add_diff(3.0)
        stats.calculate()

        assert stats.max_diff == 5.0
        assert stats.min_diff == 1.0

    def test_calculate_std_dev(self):
        """Test standard deviation calculation."""
        stats = PrecisionStats(
            calculation_type="Test",
            component="Test Component",
            unit="arcsec",
            tolerance=10.0,
        )

        # Use values with known std dev
        # [2, 4, 4, 4, 5, 5, 7, 9] has mean=5
        # Sample std dev = sqrt(sum((x-mean)^2)/(n-1)) ~= 2.14
        for val in [2, 4, 4, 4, 5, 5, 7, 9]:
            stats.add_diff(float(val))
        stats.calculate()

        assert stats.mean_diff == 5.0
        # Sample std dev is ~2.14 (not population std dev of 2.0)
        assert abs(stats.std_dev - 2.14) < 0.1

    def test_calculate_percentiles(self):
        """Test percentile calculations."""
        stats = PrecisionStats(
            calculation_type="Test",
            component="Test Component",
            unit="arcsec",
            tolerance=1000.0,
        )

        # Add 100 values from 1 to 100
        for i in range(1, 101):
            stats.add_diff(float(i))
        stats.calculate()

        # P95 should be 95 or 96 (depends on exact index calculation)
        assert 95.0 <= stats.p95_diff <= 96.0
        # P99 should be 99 or 100
        assert 99.0 <= stats.p99_diff <= 100.0

    def test_within_tolerance_pass(self):
        """Test tolerance check when within limits."""
        stats = PrecisionStats(
            calculation_type="Test",
            component="Test Component",
            unit="arcsec",
            tolerance=10.0,
        )

        stats.add_diff(1.0)
        stats.add_diff(5.0)
        stats.add_diff(9.0)  # Below tolerance of 10
        stats.calculate()

        assert stats.within_tolerance is True

    def test_within_tolerance_fail(self):
        """Test tolerance check when exceeding limits."""
        stats = PrecisionStats(
            calculation_type="Test",
            component="Test Component",
            unit="arcsec",
            tolerance=5.0,
        )

        stats.add_diff(1.0)
        stats.add_diff(3.0)
        stats.add_diff(10.0)  # Above tolerance of 5
        stats.calculate()

        assert stats.within_tolerance is False

    def test_to_dict(self):
        """Test dictionary conversion."""
        stats = PrecisionStats(
            calculation_type="Test Type",
            component="Test Component",
            unit="arcsec",
            tolerance=5.0,
        )

        stats.add_diff(1.0)
        stats.add_diff(2.0)
        stats.calculate()

        d = stats.to_dict()

        assert d["calculation_type"] == "Test Type"
        assert d["component"] == "Test Component"
        assert d["unit"] == "arcsec"
        assert d["n_samples"] == 2
        assert isinstance(d["max_diff"], float)
        assert isinstance(d["within_tolerance"], bool)
        assert "_diffs" not in d  # Internal field excluded


# =============================================================================
# TEST PRECISION REPORT
# =============================================================================


class TestPrecisionReport:
    """Test the PrecisionReport dataclass functionality."""

    def test_to_json(self):
        """Test JSON serialization."""
        stats1 = PrecisionStats(
            calculation_type="Test",
            component="Component 1",
            unit="arcsec",
            tolerance=1.0,
        )
        stats1.add_diff(0.5)
        stats1.calculate()

        report = PrecisionReport(
            generated_at="2024-01-01T00:00:00",
            num_test_dates=10,
            date_range="1900-2050",
            stats=[stats1],
        )

        json_str = report.to_json()
        parsed = json.loads(json_str)

        assert parsed["generated_at"] == "2024-01-01T00:00:00"
        assert parsed["num_test_dates"] == 10
        assert len(parsed["stats"]) == 1
        assert parsed["stats"][0]["component"] == "Component 1"

    def test_to_csv(self):
        """Test CSV serialization."""
        stats1 = PrecisionStats(
            calculation_type="Test",
            component="Component 1",
            unit="arcsec",
            tolerance=1.0,
        )
        stats1.add_diff(0.5)
        stats1.calculate()

        report = PrecisionReport(
            generated_at="2024-01-01T00:00:00",
            num_test_dates=10,
            date_range="1900-2050",
            stats=[stats1],
        )

        csv_str = report.to_csv()
        lines = csv_str.strip().split("\n")

        # Check header
        assert "calculation_type" in lines[0]
        assert "component" in lines[0]
        assert "max_diff" in lines[0]

        # Check data row
        assert "Test" in lines[1]
        assert "Component 1" in lines[1]

    def test_to_dict(self):
        """Test dictionary conversion."""
        stats1 = PrecisionStats(
            calculation_type="Test",
            component="Component 1",
            unit="arcsec",
            tolerance=1.0,
        )
        stats1.add_diff(0.5)
        stats1.calculate()

        report = PrecisionReport(
            generated_at="2024-01-01T00:00:00",
            num_test_dates=10,
            date_range="1900-2050",
            stats=[stats1],
        )

        d = report.to_dict()

        assert d["generated_at"] == "2024-01-01T00:00:00"
        assert d["num_test_dates"] == 10
        assert isinstance(d["stats"], list)


# =============================================================================
# TEST HELPER FUNCTIONS
# =============================================================================


class TestHelperFunctions:
    """Test helper functions."""

    def test_angular_diff_simple(self):
        """Test simple angular difference."""
        assert angular_diff(10.0, 15.0) == 5.0
        assert angular_diff(15.0, 10.0) == 5.0

    def test_angular_diff_wrap(self):
        """Test angular difference with wrap-around."""
        # 359 to 1 should be 2 degrees, not 358
        assert angular_diff(359.0, 1.0) == 2.0
        assert angular_diff(1.0, 359.0) == 2.0

        # 350 to 10 should be 20 degrees
        assert angular_diff(350.0, 10.0) == 20.0

    def test_angular_diff_zero(self):
        """Test angular difference of zero."""
        assert angular_diff(45.0, 45.0) == 0.0

    def test_generate_random_jds(self):
        """Test Julian Day generation."""
        jds = generate_random_jds(10, seed=42)

        assert len(jds) == 10
        # All JDs should be within DE421 range (approx JD 2415021 to 2469807)
        for jd in jds:
            assert 2415000 < jd < 2470000

    def test_generate_random_jds_reproducible(self):
        """Test that same seed produces same JDs."""
        jds1 = generate_random_jds(5, seed=123)
        jds2 = generate_random_jds(5, seed=123)

        assert jds1 == jds2

    def test_generate_random_jds_different_seeds(self):
        """Test that different seeds produce different JDs."""
        jds1 = generate_random_jds(5, seed=1)
        jds2 = generate_random_jds(5, seed=2)

        assert jds1 != jds2


# =============================================================================
# TEST MEASUREMENT FUNCTIONS
# =============================================================================


class TestMeasurementFunctions:
    """Test the precision measurement functions."""

    @pytest.fixture
    def small_jd_set(self):
        """Fixture for a small set of test JDs."""
        return generate_random_jds(5, seed=42)

    def test_measure_planetary_precision(self, small_jd_set):
        """Test planetary precision measurement."""
        results = measure_planetary_precision(small_jd_set, verbose=False)

        # Should have 20 stats (10 planets x 2 components: lon + lat)
        assert len(results) == 20

        # Check structure
        for stat in results:
            assert stat.calculation_type == "Planetary Position"
            assert stat.unit == "arcsec"
            assert stat.n_samples == 5
            assert stat.max_diff >= 0
            assert stat.mean_diff >= 0

    def test_measure_velocity_precision(self, small_jd_set):
        """Test velocity precision measurement."""
        results = measure_velocity_precision(small_jd_set, verbose=False)

        # Should have 4 stats (Sun, Moon, Mars, Jupiter)
        assert len(results) == 4

        for stat in results:
            assert stat.calculation_type == "Planetary Velocity"
            assert stat.unit == "deg/day"

    def test_measure_house_precision(self, small_jd_set):
        """Test house cusp precision measurement."""
        results = measure_house_precision(small_jd_set, verbose=False)

        # Should have 16 stats (8 house systems x 2: cusps + ascendant)
        assert len(results) == 16

        # Check calculation types
        cusp_stats = [s for s in results if s.calculation_type == "House Cusps"]
        angle_stats = [s for s in results if s.calculation_type == "House Angles"]

        assert len(cusp_stats) == 8
        assert len(angle_stats) == 8

    def test_measure_ayanamsha_precision(self, small_jd_set):
        """Test ayanamsha precision measurement."""
        results = measure_ayanamsha_precision(small_jd_set, verbose=False)

        # Should have 4 stats (Lahiri, Fagan-Bradley, Raman, True Citra)
        assert len(results) == 4

        for stat in results:
            assert stat.calculation_type == "Ayanamsha"
            assert stat.unit == "deg"

    def test_measure_lunar_node_precision(self, small_jd_set):
        """Test lunar node precision measurement."""
        results = measure_lunar_node_precision(small_jd_set, verbose=False)

        # Should have 3 stats (Mean Node, True Node, Mean Lilith)
        assert len(results) == 3

        for stat in results:
            assert stat.calculation_type == "Lunar Points"
            assert stat.unit == "arcsec"

    def test_measure_heliocentric_precision(self, small_jd_set):
        """Test heliocentric precision measurement."""
        results = measure_heliocentric_precision(small_jd_set, verbose=False)

        # Should have 4 stats (Mercury, Venus, Mars, Jupiter)
        assert len(results) == 4

        for stat in results:
            assert stat.calculation_type == "Heliocentric"
            assert "Longitude" in stat.component

    def test_measure_time_precision(self, small_jd_set):
        """Test time function precision measurement."""
        results = measure_time_precision(small_jd_set, verbose=False)

        # Should have 2 stats (Delta T, Julian Day Conversion)
        assert len(results) == 2

        for stat in results:
            assert stat.calculation_type == "Time Functions"


# =============================================================================
# TEST FULL REPORT GENERATION
# =============================================================================


class TestReportGeneration:
    """Test full report generation."""

    def test_generate_report_basic(self):
        """Test basic report generation with minimal samples."""
        report = generate_report(num_dates=3, verbose=False, seed=42)

        assert report.num_test_dates == 3
        assert report.date_range == "1900-2050 (DE421 range)"
        assert len(report.stats) > 0
        assert report.generated_at is not None

    def test_generate_report_has_all_categories(self):
        """Test that report includes all calculation categories."""
        report = generate_report(num_dates=3, verbose=False, seed=42)

        calculation_types = set(s.calculation_type for s in report.stats)

        expected_types = {
            "Planetary Position",
            "Planetary Velocity",
            "House Cusps",
            "House Angles",
            "Ayanamsha",
            "Lunar Points",
            "Heliocentric",
            "Time Functions",
        }

        assert calculation_types == expected_types

    def test_generate_report_stats_count(self):
        """Test total number of stats in report."""
        report = generate_report(num_dates=3, verbose=False, seed=42)

        # Expected: 20 planetary + 4 velocity + 16 houses + 4 ayanamsha +
        #           3 lunar + 4 heliocentric + 2 time = 53
        assert len(report.stats) == 53

    def test_generate_report_reproducible(self):
        """Test that same seed produces same results."""
        report1 = generate_report(num_dates=5, verbose=False, seed=123)
        report2 = generate_report(num_dates=5, verbose=False, seed=123)

        # Compare stats
        for s1, s2 in zip(report1.stats, report2.stats):
            assert s1.max_diff == s2.max_diff
            assert s1.mean_diff == s2.mean_diff
            assert s1.component == s2.component

    def test_generate_report_json_serializable(self):
        """Test that report can be serialized to valid JSON."""
        report = generate_report(num_dates=3, verbose=False, seed=42)

        json_str = report.to_json()
        parsed = json.loads(json_str)

        assert "generated_at" in parsed
        assert "stats" in parsed
        assert len(parsed["stats"]) == 53

    def test_generate_report_csv_serializable(self):
        """Test that report can be serialized to valid CSV."""
        report = generate_report(num_dates=3, verbose=False, seed=42)

        csv_str = report.to_csv()
        lines = csv_str.strip().split("\n")

        # Header + 53 data rows
        assert len(lines) == 54

        # Check header has expected columns
        header = lines[0].split(",")
        assert "calculation_type" in header
        assert "max_diff" in header
        assert "within_tolerance" in header


# =============================================================================
# TEST TOLERANCE VALIDATION
# =============================================================================


class TestToleranceValidation:
    """Test tolerance checking in the report."""

    def test_planetary_tolerances_within_limits(self):
        """Test that planetary precision is within documented tolerances."""
        report = generate_report(num_dates=50, verbose=False, seed=42)

        planetary_stats = [
            s for s in report.stats if s.calculation_type == "Planetary Position"
        ]

        # Most should pass tolerance
        passed = sum(1 for s in planetary_stats if s.within_tolerance)
        total = len(planetary_stats)

        # At least 90% should pass
        assert passed / total >= 0.9, f"Only {passed}/{total} planetary stats passed"

    def test_house_tolerances_within_limits(self):
        """Test that house cusp precision is within documented tolerances."""
        report = generate_report(num_dates=50, verbose=False, seed=42)

        house_stats = [
            s
            for s in report.stats
            if s.calculation_type in ("House Cusps", "House Angles")
        ]

        # All should pass tolerance for house cusps
        for stat in house_stats:
            assert stat.within_tolerance, (
                f"{stat.component} exceeded tolerance: "
                f"{stat.max_diff:.4f} > {stat.tolerance:.4f}"
            )

    def test_ayanamsha_tolerances_within_limits(self):
        """Test that ayanamsha precision is within documented tolerances."""
        report = generate_report(num_dates=50, verbose=False, seed=42)

        aya_stats = [s for s in report.stats if s.calculation_type == "Ayanamsha"]

        # All should pass tolerance
        for stat in aya_stats:
            assert stat.within_tolerance, (
                f"{stat.component} exceeded tolerance: "
                f"{stat.max_diff:.6f} > {stat.tolerance:.6f}"
            )
