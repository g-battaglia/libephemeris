"""
Tests for the update_orbital_elements.py script.

These tests verify the script's functionality including:
- Body definitions and SBDB ID mapping
- Element comparison logic
- Report generation
- File update logic (mocked)
- API response parsing (mocked)
"""

import json
import os
import sys
import pytest
from unittest.mock import patch, MagicMock
from io import StringIO

# Add scripts directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))


class TestBodySbdbIds:
    """Test body SBDB ID definitions."""

    def test_all_minor_bodies_have_sbdb_ids(self):
        """Test that all bodies in MINOR_BODY_ELEMENTS have SBDB mappings."""
        from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS
        from scripts.update_orbital_elements import BODY_SBDB_IDS

        # Get body names from MINOR_BODY_ELEMENTS
        body_names = {elem.name for elem in MINOR_BODY_ELEMENTS.values()}

        # Check all are in SBDB mappings
        for name in body_names:
            assert name in BODY_SBDB_IDS, f"Body {name} missing from BODY_SBDB_IDS"

    def test_sbdb_ids_are_valid_strings(self):
        """Test that all SBDB IDs are non-empty strings."""
        from scripts.update_orbital_elements import BODY_SBDB_IDS

        for name, sbdb_id in BODY_SBDB_IDS.items():
            assert isinstance(sbdb_id, str), f"{name} SBDB ID should be string"
            assert len(sbdb_id) > 0, f"{name} SBDB ID should not be empty"
            # Most are numeric IDs
            assert sbdb_id.replace("-", "").isdigit() or sbdb_id.isalnum(), (
                f"{name} SBDB ID '{sbdb_id}' has unexpected format"
            )

    def test_expected_bodies_present(self):
        """Test that expected bodies are in the mapping."""
        from scripts.update_orbital_elements import BODY_SBDB_IDS

        expected_bodies = [
            "Chiron",
            "Pholus",
            "Ceres",
            "Pallas",
            "Juno",
            "Vesta",
            "Eris",
            "Sedna",
            "Haumea",
            "Makemake",
            "Ixion",
            "Orcus",
            "Quaoar",
        ]

        for body in expected_bodies:
            assert body in BODY_SBDB_IDS, f"Expected body {body} not in BODY_SBDB_IDS"

    def test_body_count_matches(self):
        """Test that SBDB ID count matches MINOR_BODY_ELEMENTS count."""
        from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS
        from scripts.update_orbital_elements import BODY_SBDB_IDS

        assert len(BODY_SBDB_IDS) == len(MINOR_BODY_ELEMENTS), (
            f"SBDB ID count ({len(BODY_SBDB_IDS)}) doesn't match "
            f"MINOR_BODY_ELEMENTS count ({len(MINOR_BODY_ELEMENTS)})"
        )


class TestFetchedElements:
    """Test FetchedElements dataclass."""

    def test_fetched_elements_creation(self):
        """Test creating FetchedElements with valid data."""
        from scripts.update_orbital_elements import FetchedElements

        elem = FetchedElements(
            name="Test",
            epoch=2460000.5,
            a=2.77,
            e=0.08,
            i=10.5,
            omega=73.3,
            Omega=80.2,
            M0=231.5,
            n=0.214,
        )

        assert elem.name == "Test"
        assert elem.epoch == 2460000.5
        assert elem.a == 2.77
        assert elem.e == 0.08
        assert elem.i == 10.5
        assert elem.omega == 73.3
        assert elem.Omega == 80.2
        assert elem.M0 == 231.5
        assert elem.n == 0.214


class TestCompareElements:
    """Test element comparison logic."""

    def test_compare_identical_elements(self):
        """Test comparing identical elements returns no differences."""
        from scripts.update_orbital_elements import compare_elements, FetchedElements
        from libephemeris.minor_bodies import OrbitalElements

        current = OrbitalElements(
            name="Test",
            epoch=2460000.5,
            a=2.77,
            e=0.08,
            i=10.5,
            omega=73.3,
            Omega=80.2,
            M0=231.5,
            n=0.214,
        )

        fetched = FetchedElements(
            name="Test",
            epoch=2460000.5,
            a=2.77,
            e=0.08,
            i=10.5,
            omega=73.3,
            Omega=80.2,
            M0=231.5,
            n=0.214,
        )

        diffs = compare_elements(current, fetched, threshold=0.001)
        assert len(diffs) == 0, "Identical elements should have no differences"

    def test_compare_different_elements(self):
        """Test comparing different elements returns differences."""
        from scripts.update_orbital_elements import compare_elements, FetchedElements
        from libephemeris.minor_bodies import OrbitalElements

        current = OrbitalElements(
            name="Test",
            epoch=2460000.5,
            a=2.77,
            e=0.08,
            i=10.5,
            omega=73.3,
            Omega=80.2,
            M0=231.5,
            n=0.214,
        )

        # Changed semi-major axis by 1%
        fetched = FetchedElements(
            name="Test",
            epoch=2460000.5,
            a=2.7977,  # ~1% different
            e=0.08,
            i=10.5,
            omega=73.3,
            Omega=80.2,
            M0=231.5,
            n=0.214,
        )

        diffs = compare_elements(current, fetched, threshold=0.001)
        assert "a" in diffs, "Should detect semi-major axis difference"
        assert len(diffs) == 1, "Should only have one difference"

    def test_compare_with_threshold(self):
        """Test that threshold correctly filters differences."""
        from scripts.update_orbital_elements import compare_elements, FetchedElements
        from libephemeris.minor_bodies import OrbitalElements

        current = OrbitalElements(
            name="Test",
            epoch=2460000.5,
            a=2.77,
            e=0.08,
            i=10.5,
            omega=73.3,
            Omega=80.2,
            M0=231.5,
            n=0.214,
        )

        # Changed by 0.05% (below 0.1% threshold)
        fetched = FetchedElements(
            name="Test",
            epoch=2460000.5,
            a=2.771385,  # 0.05% different
            e=0.08,
            i=10.5,
            omega=73.3,
            Omega=80.2,
            M0=231.5,
            n=0.214,
        )

        # With 0.1% threshold, should not report
        diffs = compare_elements(current, fetched, threshold=0.001)
        assert "a" not in diffs, "Small difference should be filtered by threshold"

        # With 0.01% threshold, should report
        diffs = compare_elements(current, fetched, threshold=0.0001)
        assert "a" in diffs, "Small difference should be reported with lower threshold"

    def test_compare_zero_current_value(self):
        """Test comparing when current value is zero."""
        from scripts.update_orbital_elements import compare_elements, FetchedElements
        from libephemeris.minor_bodies import OrbitalElements

        current = OrbitalElements(
            name="Test",
            epoch=2460000.5,
            a=2.77,
            e=0.0,  # Zero eccentricity
            i=10.5,
            omega=73.3,
            Omega=80.2,
            M0=231.5,
            n=0.214,
        )

        fetched = FetchedElements(
            name="Test",
            epoch=2460000.5,
            a=2.77,
            e=0.001,  # Non-zero
            i=10.5,
            omega=73.3,
            Omega=80.2,
            M0=231.5,
            n=0.214,
        )

        diffs = compare_elements(current, fetched, threshold=0.001)
        assert "e" in diffs, "Should handle zero current value"


class TestGenerateReport:
    """Test report generation."""

    def test_generate_empty_report(self):
        """Test report when no changes needed."""
        from scripts.update_orbital_elements import generate_report

        all_differences: dict[str, dict[str, tuple[float, float, float]]] = {
            "Ceres": {},
            "Chiron": {},
        }
        fetched_epochs: dict[str, float] = {}

        report = generate_report(all_differences, fetched_epochs)

        assert "ORBITAL ELEMENTS UPDATE REPORT" in report
        assert "All orbital elements are up to date" in report
        assert "No changes needed" in report

    def test_generate_report_with_changes(self):
        """Test report with changes."""
        from scripts.update_orbital_elements import generate_report

        all_differences: dict[str, dict[str, tuple[float, float, float]]] = {
            "Ceres": {
                "a": (2.77, 2.78, 0.36),
                "e": (0.08, 0.082, 2.5),
            },
        }
        fetched_epochs: dict[str, float] = {"Ceres": 2461000.5}

        report = generate_report(all_differences, fetched_epochs)

        assert "ORBITAL ELEMENTS UPDATE REPORT" in report
        assert "Ceres" in report
        assert "2.77" in report or "2.78" in report
        assert "python scripts/update_orbital_elements.py --update" in report

    def test_generate_report_shows_body_count(self):
        """Test that report shows number of bodies needing updates."""
        from scripts.update_orbital_elements import generate_report

        all_differences: dict[str, dict[str, tuple[float, float, float]]] = {
            "Ceres": {"a": (2.77, 2.78, 0.36)},
            "Chiron": {},  # No changes
            "Eris": {"e": (0.44, 0.45, 2.27)},
        }
        fetched_epochs: dict[str, float] = {}

        report = generate_report(all_differences, fetched_epochs)

        assert "2 bodies" in report  # 2 bodies with outdated elements


class TestCheckRequests:
    """Test requests library checking."""

    def test_check_requests_returns_boolean(self):
        """Test check returns a boolean value."""
        from scripts.update_orbital_elements import check_requests

        result = check_requests()
        assert isinstance(result, bool)

    def test_check_requests_when_not_available(self):
        """Test check when requests is not available."""
        import scripts.update_orbital_elements as update_module

        # Temporarily set requests to None
        original = update_module.requests
        update_module.requests = None

        try:
            result = update_module.check_requests()
            assert result is False
        finally:
            update_module.requests = original


# Check if requests is available for network tests
try:
    import requests as requests_module

    REQUESTS_AVAILABLE = True
except ImportError:
    REQUESTS_AVAILABLE = False
    requests_module = None  # type: ignore


class TestFetchOrbitalElements:
    """Test fetching orbital elements from API."""

    def test_fetch_unknown_body_returns_none(self):
        """Test that fetching unknown body returns None."""
        from scripts.update_orbital_elements import fetch_orbital_elements

        result = fetch_orbital_elements("UnknownBody123", verbose=False)
        assert result is None

    @pytest.mark.skipif(not REQUESTS_AVAILABLE, reason="requests not installed")
    def test_fetch_with_mocked_api_success(self):
        """Test fetching with mocked successful API response."""
        from scripts.update_orbital_elements import fetch_orbital_elements

        mock_response = {
            "orbit": {
                "epoch": 2461000.5,
                "elements": [
                    {"name": "a", "value": "2.77"},
                    {"name": "e", "value": "0.08"},
                    {"name": "i", "value": "10.5"},
                    {"name": "w", "value": "73.3"},
                    {"name": "om", "value": "80.2"},
                    {"name": "ma", "value": "231.5"},
                    {"name": "n", "value": "0.214"},
                ],
            }
        }

        with patch("requests.get") as mock_get:
            mock_get.return_value.json.return_value = mock_response
            mock_get.return_value.raise_for_status = MagicMock()

            result = fetch_orbital_elements("Ceres", verbose=False)

            assert result is not None
            assert result.name == "Ceres"
            assert result.epoch == 2461000.5
            assert result.a == 2.77
            assert result.e == 0.08

    @pytest.mark.skipif(not REQUESTS_AVAILABLE, reason="requests not installed")
    def test_fetch_with_api_error(self):
        """Test handling API error response."""
        from scripts.update_orbital_elements import fetch_orbital_elements

        mock_response = {"error": "Body not found"}

        with patch("requests.get") as mock_get:
            mock_get.return_value.json.return_value = mock_response
            mock_get.return_value.raise_for_status = MagicMock()

            result = fetch_orbital_elements("Ceres", verbose=False)
            assert result is None

    @pytest.mark.skipif(not REQUESTS_AVAILABLE, reason="requests not installed")
    def test_fetch_with_network_error(self):
        """Test handling network errors."""
        from scripts.update_orbital_elements import fetch_orbital_elements

        with patch("requests.get") as mock_get:
            mock_get.side_effect = requests_module.exceptions.RequestException(  # type: ignore[union-attr]
                "Network error"
            )

            result = fetch_orbital_elements("Ceres", verbose=False)
            assert result is None


class TestGeneratePythonCode:
    """Test Python code generation."""

    def test_generate_python_code_format(self):
        """Test that generated code has correct format."""
        from scripts.update_orbital_elements import (
            generate_python_code,
            FetchedElements,
        )

        fetched = FetchedElements(
            name="Ceres",
            epoch=2461000.5,
            a=2.765615651508659,
            e=0.07957631994408416,
            i=10.58788658206854,
            omega=73.29975464616518,
            Omega=80.24963090816965,
            M0=231.5397330043706,
            n=0.2142971214271186,
        )

        code = generate_python_code(fetched, "CERES")

        assert "CERES: OrbitalElements(" in code
        assert 'name="Ceres"' in code
        assert "epoch=2461000.5" in code
        assert "a=2.765615651508659" in code


class TestScriptIntegration:
    """Integration tests for the script."""

    def test_script_imports_correctly(self):
        """Test that the script can be imported without errors."""
        import scripts.update_orbital_elements as update_module

        assert hasattr(update_module, "main")
        assert hasattr(update_module, "fetch_orbital_elements")
        assert hasattr(update_module, "compare_elements")
        assert hasattr(update_module, "generate_report")
        assert hasattr(update_module, "BODY_SBDB_IDS")

    def test_main_with_help_flag(self):
        """Test that --help works."""
        import scripts.update_orbital_elements as update_module

        with patch.object(sys, "argv", ["update_orbital_elements.py", "--help"]):
            with pytest.raises(SystemExit) as exc_info:
                update_module.main()
            assert exc_info.value.code == 0

    def test_main_with_unknown_body(self):
        """Test main with unknown body name."""
        import scripts.update_orbital_elements as update_module

        with patch.object(
            sys, "argv", ["update_orbital_elements.py", "--body", "unknownbody123"]
        ):
            result = update_module.main()
            assert result == 1  # Should fail

    @pytest.mark.skipif(not REQUESTS_AVAILABLE, reason="requests not installed")
    def test_main_with_valid_body_mocked(self):
        """Test main with valid body and mocked API."""
        import scripts.update_orbital_elements as update_module

        mock_response = {
            "orbit": {
                "epoch": 2461000.5,
                "elements": [
                    {"name": "a", "value": "2.77"},
                    {"name": "e", "value": "0.08"},
                    {"name": "i", "value": "10.5"},
                    {"name": "w", "value": "73.3"},
                    {"name": "om", "value": "80.2"},
                    {"name": "ma", "value": "231.5"},
                    {"name": "n", "value": "0.214"},
                ],
            }
        }

        with patch("requests.get") as mock_get:
            mock_get.return_value.json.return_value = mock_response
            mock_get.return_value.raise_for_status = MagicMock()

            with patch.object(
                sys,
                "argv",
                ["update_orbital_elements.py", "--body", "ceres", "--quiet"],
            ):
                # Capture stdout
                captured = StringIO()
                with patch.object(sys, "stdout", captured):
                    result = update_module.main()

                # Should succeed
                assert result == 0


class TestSbdbApiEndpoint:
    """Test SBDB API endpoint configuration."""

    def test_api_endpoint_is_valid_url(self):
        """Test that API endpoint is a valid URL."""
        from scripts.update_orbital_elements import SBDB_API_URL

        assert SBDB_API_URL.startswith("https://")
        assert "jpl.nasa.gov" in SBDB_API_URL

    def test_api_endpoint_is_sbdb(self):
        """Test that endpoint is the SBDB API."""
        from scripts.update_orbital_elements import SBDB_API_URL

        assert "sbdb" in SBDB_API_URL.lower()
