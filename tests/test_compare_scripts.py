"""
Tests for compare_scripts precision validation suite.

Tests verify that:
1. Comparison scripts run without errors
2. Comparison modules are properly registered in run_all_compare
3. Hypothetical planet comparisons execute
4. Minor body comparisons execute with proper category handling
"""

import sys
from pathlib import Path
import subprocess

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "compare_scripts"))
sys.path.insert(0, str(Path(__file__).parent.parent / "compare_scripts" / "tests"))


class TestCompareScriptsInfrastructure:
    """Test that comparison infrastructure is properly set up."""

    def test_comparison_modules_list(self):
        """Test that all comparison modules are registered."""
        from run_all_compare import COMPARISON_MODULES

        # Check we have a reasonable number of modules
        assert len(COMPARISON_MODULES) >= 19, (
            "Should have at least 19 comparison modules"
        )

        # Check each module has correct structure
        for module in COMPARISON_MODULES:
            assert len(module) == 4, f"Module {module} should have 4 elements"
            name, script, desc, is_quick = module
            assert isinstance(name, str)
            assert isinstance(script, str)
            assert script.endswith(".py")
            assert isinstance(desc, str)
            assert isinstance(is_quick, bool)

    def test_hypothetical_module_registered(self):
        """Test that hypothetical module is in the comparison list."""
        from run_all_compare import COMPARISON_MODULES

        module_names = [m[0] for m in COMPARISON_MODULES]
        assert "hypothetical" in module_names, (
            "hypothetical module should be registered"
        )

    def test_minor_bodies_module_registered(self):
        """Test that minor_bodies module is in the comparison list."""
        from run_all_compare import COMPARISON_MODULES

        module_names = [m[0] for m in COMPARISON_MODULES]
        assert "minor_bodies" in module_names, (
            "minor_bodies module should be registered"
        )


class TestCompareHypothetical:
    """Test hypothetical planet comparison script."""

    def test_import_compare_hypothetical(self):
        """Test that test_compare_hypothetical module exists and has expected content."""
        # The comparison tests have been converted to pytest-style tests
        # in compare_scripts/tests/test_compare_hypothetical.py
        from test_compare_hypothetical import (
            URANIAN_PLANETS,
            OTHER_HYPOTHETICAL,
        )

        assert len(URANIAN_PLANETS) == 8, "Should have 8 Uranian planets"
        assert len(OTHER_HYPOTHETICAL) >= 1, "Should have at least Transpluto"

    def test_uranian_planets_defined(self):
        """Test all 8 Uranian planets are defined."""
        from test_compare_hypothetical import URANIAN_PLANETS

        names = [p[2] for p in URANIAN_PLANETS]
        expected = [
            "Cupido",
            "Hades",
            "Zeus",
            "Kronos",
            "Apollon",
            "Admetos",
            "Vulkanus",
            "Poseidon",
        ]
        assert names == expected, f"Expected {expected}, got {names}"


class TestCompareMinorBodies:
    """Test minor bodies comparison script."""

    def test_import_compare_minor_bodies(self):
        """Test that test_compare_minor_bodies module exists and has expected content."""
        # The comparison tests have been converted to pytest-style tests
        # in compare_scripts/tests/test_compare_minor_bodies.py
        from test_compare_minor_bodies import (
            MAIN_ASTEROIDS,
            CENTAURS,
        )

        assert len(MAIN_ASTEROIDS) >= 4, (
            "Should have at least Ceres, Pallas, Juno, Vesta"
        )
        assert len(CENTAURS) >= 2, "Should have at least Chiron, Pholus"

    def test_tolerances_defined(self):
        """Test tolerances are defined in the module."""
        from test_compare_minor_bodies import (
            MAIN_ASTEROID_TOL,
            CENTAUR_TOL,
            TNO_TOL,
            NEA_TOL,
        )

        assert MAIN_ASTEROID_TOL > 0
        assert CENTAUR_TOL > 0
        assert TNO_TOL > 0
        assert NEA_TOL > 0

    def test_tnos_list(self):
        """Test TNO list includes major bodies."""
        from test_compare_minor_bodies import TNOS_STANDARD, TNOS_DISTANT

        all_tnos = TNOS_STANDARD + TNOS_DISTANT
        names = [t[1] for t in all_tnos]
        expected_in_list = ["Eris", "Makemake", "Sedna", "Orcus"]
        for expected in expected_in_list:
            assert expected in names, f"{expected} should be in TNO list"


class TestCompareLunar:
    """Test lunar comparison script."""

    def test_import_compare_lunar(self):
        """Test that test_compare_lunar module exists and has expected content."""
        # The comparison tests have been converted to pytest-style tests
        # in compare_scripts/tests/test_compare_lunar.py
        from test_compare_lunar import (
            MEAN_NODE_TOL,
            TRUE_NODE_TOL,
            MEAN_LILITH_TOL,
        )

        assert MEAN_NODE_TOL > 0
        assert TRUE_NODE_TOL > 0
        assert MEAN_LILITH_TOL > 0

    def test_threshold_constants(self):
        """Test threshold constants are defined."""
        from test_compare_lunar import (
            TRUE_NODE_TOL,
            TRUE_LILITH_TOL,
        )

        assert TRUE_NODE_TOL > 0
        assert TRUE_LILITH_TOL > 0


class TestComparisonUtils:
    """Test shared comparison utilities."""

    def test_import_comparison_utils(self):
        """Test that comparison_utils can be imported."""
        from comparison_utils import (
            Tolerances,
            angular_diff,
            format_coord,
            format_diff,
            format_status,
            TestStatistics,
            print_header,
            print_section,
        )

    def test_angular_diff(self):
        """Test angular difference calculation."""
        from comparison_utils import angular_diff

        # Simple case
        assert angular_diff(10.0, 15.0) == 5.0
        assert angular_diff(15.0, 10.0) == 5.0

        # Wrap around 360
        assert angular_diff(5.0, 355.0) == 10.0
        assert angular_diff(355.0, 5.0) == 10.0

    def test_tolerances_defined(self):
        """Test that tolerances are defined."""
        from comparison_utils import Tolerances

        assert Tolerances.LONGITUDE_STRICT > 0
        assert Tolerances.LONGITUDE_RELAXED > Tolerances.LONGITUDE_STRICT
        assert Tolerances.HOUSE_CUSP > 0
        assert Tolerances.AYANAMSHA > 0

    def test_test_statistics(self):
        """Test TestStatistics class."""
        from comparison_utils import TestStatistics

        stats = TestStatistics()
        assert stats.total == 0
        assert stats.passed == 0

        stats.add_result(True, 0.001)
        assert stats.total == 1
        assert stats.passed == 1

        stats.add_result(False, 0.1)
        assert stats.total == 2
        assert stats.passed == 1
        assert stats.failed == 1


class TestRunAllCompare:
    """Test run_all_compare.py master runner."""

    def test_module_result_dataclass(self):
        """Test ModuleResult dataclass."""
        from run_all_compare import ModuleResult

        result = ModuleResult(
            name="test",
            passed=8,
            total=10,
            duration=1.5,
        )

        assert result.pass_rate == 80.0
        assert result.status == "FAIL"

        result_pass = ModuleResult(name="test2", passed=10, total=10, duration=0.5)
        assert result_pass.status == "PASS"

    def test_print_summary_report(self):
        """Test summary report function exists."""
        from run_all_compare import print_summary_report

        # Just verify it can be called (output goes to stdout)
        from run_all_compare import ModuleResult

        results = [
            ModuleResult(name="test", passed=10, total=10, duration=0.1),
        ]
        # This returns a boolean
        result = print_summary_report(results)
        assert isinstance(result, bool)
