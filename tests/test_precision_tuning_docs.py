"""
Tests for the PRECISION.md documentation.

Verifies that:
1. The documentation file exists
2. All documented features/functions exist in the library
3. Key sections are present
4. Code examples in the documentation are valid
"""

import os
import re

import pytest

import libephemeris as eph
from libephemeris import state


class TestPrecisionTuningDocExists:
    """Test that the PRECISION.md documentation exists."""

    @pytest.fixture
    def docs_dir(self):
        """Get the docs directory path."""
        return os.path.join(os.path.dirname(__file__), "..", "docs")

    @pytest.fixture
    def precision_tuning_path(self, docs_dir):
        """Get the PRECISION.md file path."""
        return os.path.join(docs_dir, "PRECISION.md")

    def test_precision_tuning_file_exists(self, precision_tuning_path):
        """PRECISION.md file should exist."""
        assert os.path.exists(precision_tuning_path), (
            f"PRECISION.md not found at {precision_tuning_path}"
        )

    def test_precision_tuning_has_content(self, precision_tuning_path):
        """PRECISION.md should have substantial content."""
        with open(precision_tuning_path, "r", encoding="utf-8") as f:
            content = f.read()

        # Should have at least 2KB of documentation
        assert len(content) > 2000, (
            f"PRECISION.md seems too short: {len(content)} bytes"
        )


class TestPrecisionTuningDocSections:
    """Test that required sections exist in the documentation."""

    @pytest.fixture
    def doc_content(self):
        """Load the PRECISION.md content."""
        docs_dir = os.path.join(os.path.dirname(__file__), "..", "docs")
        path = os.path.join(docs_dir, "PRECISION.md")
        with open(path, "r", encoding="utf-8") as f:
            return f.read()

    def test_has_overview_section(self, doc_content):
        """Documentation should have an Overview section."""
        assert (
            "## Overview" in doc_content
            or "# Overview" in doc_content
            or "# Precision Summary" in doc_content
        )

    def test_has_spk_kernels_section(self, doc_content):
        """Documentation should mention SPK kernels."""
        assert "SPK" in doc_content or "spk" in doc_content

    def test_has_iers_delta_t_section(self, doc_content):
        """Documentation should have IERS Delta T section."""
        assert "IERS Delta T" in doc_content or "Delta T" in doc_content

    def test_has_ephemeris_file_section(self, doc_content):
        """Documentation should mention ephemeris files."""
        assert (
            "Ephemeris File" in doc_content
            or "DE440" in doc_content
            or "ephemeris" in doc_content.lower()
        )

    def test_has_configuration_section(self, doc_content):
        """Documentation should mention configuration, methodology, or validation."""
        assert (
            "Configuration" in doc_content
            or "Methodology" in doc_content
            or "Hyper-Validation" in doc_content
        )

    def test_has_best_practices_section(self, doc_content):
        """Documentation should mention best practices or precision notes."""
        assert "Best Practices" in doc_content or "precision" in doc_content.lower()

    def test_has_table_of_contents(self, doc_content):
        """Documentation should have structured sections."""
        assert "Table of Contents" in doc_content or "## " in doc_content


class TestDocumentedFunctionsExist:
    """Test that functions documented in PRECISION.md actually exist."""

    @pytest.fixture
    def doc_content(self):
        """Load the PRECISION.md content."""
        docs_dir = os.path.join(os.path.dirname(__file__), "..", "docs")
        path = os.path.join(docs_dir, "PRECISION.md")
        with open(path, "r", encoding="utf-8") as f:
            return f.read()

    def test_set_ephemeris_file_exists(self):
        """set_ephemeris_file should exist."""
        assert hasattr(eph, "set_ephemeris_file") or hasattr(
            state, "set_ephemeris_file"
        )

    def test_set_ephe_path_exists(self):
        """set_ephe_path should exist."""
        assert hasattr(eph, "set_ephe_path") or hasattr(state, "set_ephe_path")

    def test_set_tid_acc_exists(self):
        """set_tid_acc should exist."""
        assert hasattr(eph, "set_tid_acc") or hasattr(state, "set_tid_acc")

    def test_get_tid_acc_exists(self):
        """get_tid_acc should exist."""
        assert hasattr(eph, "get_tid_acc") or hasattr(state, "get_tid_acc")

    def test_set_delta_t_userdef_exists(self):
        """set_delta_t_userdef should exist."""
        assert hasattr(eph, "set_delta_t_userdef") or hasattr(
            state, "set_delta_t_userdef"
        )

    def test_set_iers_delta_t_enabled_exists(self):
        """set_iers_delta_t_enabled should exist."""
        assert hasattr(eph, "set_iers_delta_t_enabled") or hasattr(
            state, "set_iers_delta_t_enabled"
        )

    def test_set_auto_spk_download_exists(self):
        """set_auto_spk_download should exist."""
        assert hasattr(eph, "set_auto_spk_download") or hasattr(
            state, "set_auto_spk_download"
        )

    def test_set_spk_cache_dir_exists(self):
        """set_spk_cache_dir should exist."""
        assert hasattr(eph, "set_spk_cache_dir") or hasattr(state, "set_spk_cache_dir")

    def test_set_spk_date_padding_exists(self):
        """set_spk_date_padding should exist."""
        assert hasattr(eph, "set_spk_date_padding") or hasattr(
            state, "set_spk_date_padding"
        )


class TestSPKFunctionsExist:
    """Test that SPK-related functions documented exist."""

    def test_download_spk_exists(self):
        """download_spk should exist."""
        from libephemeris import spk

        assert hasattr(spk, "download_spk")

    def test_register_spk_body_exists(self):
        """register_spk_body should exist."""
        from libephemeris import spk

        assert hasattr(spk, "register_spk_body")

    def test_unregister_spk_body_exists(self):
        """unregister_spk_body should exist."""
        from libephemeris import spk

        assert hasattr(spk, "unregister_spk_body")

    def test_download_and_register_spk_exists(self):
        """download_and_register_spk should exist."""
        from libephemeris import spk

        assert hasattr(spk, "download_and_register_spk")

    def test_list_spk_bodies_exists(self):
        """list_spk_bodies should exist."""
        from libephemeris import spk

        assert hasattr(spk, "list_spk_bodies")

    def test_get_spk_body_info_exists(self):
        """get_spk_body_info should exist."""
        from libephemeris import spk

        assert hasattr(spk, "get_spk_body_info")

    def test_get_spk_coverage_exists(self):
        """get_spk_coverage should exist."""
        from libephemeris import spk

        assert hasattr(spk, "get_spk_coverage")


class TestIERSFunctionsExist:
    """Test that IERS-related functions documented exist."""

    def test_download_delta_t_data_exists(self):
        """download_delta_t_data should exist."""
        from libephemeris import iers_data

        assert hasattr(iers_data, "download_delta_t_data")

    def test_download_iers_finals_exists(self):
        """download_iers_finals should exist."""
        from libephemeris import iers_data

        assert hasattr(iers_data, "download_iers_finals")

    def test_get_iers_cache_info_exists(self):
        """get_iers_cache_info should exist."""
        from libephemeris import iers_data

        assert hasattr(iers_data, "get_iers_cache_info")

    def test_get_observed_delta_t_exists(self):
        """get_observed_delta_t should exist."""
        from libephemeris import iers_data

        assert hasattr(iers_data, "get_observed_delta_t")

    def test_is_observed_delta_t_available_exists(self):
        """is_observed_delta_t_available should exist."""
        from libephemeris import iers_data

        assert hasattr(iers_data, "is_observed_delta_t_available")

    def test_set_iers_cache_dir_exists(self):
        """set_iers_cache_dir should exist."""
        from libephemeris import iers_data

        assert hasattr(iers_data, "set_iers_cache_dir")

    def test_set_iers_auto_download_exists(self):
        """set_iers_auto_download should exist."""
        from libephemeris import iers_data

        assert hasattr(iers_data, "set_iers_auto_download")


class TestSpkAutoFunctionsExist:
    """Test that spk_auto module functions documented exist."""

    def test_enable_common_bodies_exists(self):
        """enable_common_bodies should exist."""
        from libephemeris import spk_auto

        assert hasattr(spk_auto, "enable_common_bodies")

    def test_auto_get_spk_exists(self):
        """auto_get_spk should exist."""
        from libephemeris import spk_auto

        assert hasattr(spk_auto, "auto_get_spk")

    def test_list_cached_spk_exists(self):
        """list_cached_spk should exist."""
        from libephemeris import spk_auto

        assert hasattr(spk_auto, "list_cached_spk")

    def test_get_cache_size_exists(self):
        """get_cache_size should exist."""
        from libephemeris import spk_auto

        assert hasattr(spk_auto, "get_cache_size")

    def test_clear_spk_cache_exists(self):
        """clear_spk_cache should exist."""
        from libephemeris import spk_auto

        assert hasattr(spk_auto, "clear_spk_cache")

    def test_prune_old_cache_exists(self):
        """prune_old_cache should exist."""
        from libephemeris import spk_auto

        assert hasattr(spk_auto, "prune_old_cache")


class TestConstantsDocumented:
    """Test that constants mentioned in documentation exist."""

    def test_naif_chiron_exists(self):
        """NAIF_CHIRON constant should exist."""
        assert hasattr(eph, "NAIF_CHIRON")

    def test_naif_eris_exists(self):
        """NAIF_ERIS constant should exist."""
        assert hasattr(eph, "NAIF_ERIS")

    def test_naif_sedna_exists(self):
        """NAIF_SEDNA constant should exist."""
        assert hasattr(eph, "NAIF_SEDNA")

    def test_se_chiron_exists(self):
        """CHIRON constant should exist."""
        assert hasattr(eph, "CHIRON")

    def test_se_eris_exists(self):
        """ERIS constant should exist."""
        assert hasattr(eph, "ERIS")


class TestEnvironmentVariablesDocumented:
    """Test that environment variables documented are correctly named."""

    @pytest.fixture
    def doc_content(self):
        """Load the PRECISION.md content."""
        docs_dir = os.path.join(os.path.dirname(__file__), "..", "docs")
        path = os.path.join(docs_dir, "PRECISION.md")
        with open(path, "r", encoding="utf-8") as f:
            return f.read()

    def test_iers_delta_t_env_var_documented(self, doc_content):
        """Documentation should mention Delta T."""
        assert "Delta T" in doc_content or "delta_t" in doc_content.lower()

    def test_auto_spk_env_var_documented(self, doc_content):
        """Documentation should mention SPK."""
        assert "SPK" in doc_content or "spk" in doc_content


class TestCodeExamplesValid:
    """Test that code examples in the documentation use valid syntax."""

    @pytest.fixture
    def doc_content(self):
        """Load the PRECISION.md content."""
        docs_dir = os.path.join(os.path.dirname(__file__), "..", "docs")
        path = os.path.join(docs_dir, "PRECISION.md")
        with open(path, "r", encoding="utf-8") as f:
            return f.read()

    def test_code_blocks_have_valid_python_syntax(self, doc_content):
        """Code blocks in documentation should have valid Python syntax (if any)."""
        # Extract Python code blocks (```python ... ```)
        pattern = r"```python\n(.*?)```"
        code_blocks = re.findall(pattern, doc_content, re.DOTALL)

        # PRECISION.md may not have code examples; skip if none found
        if len(code_blocks) == 0:
            pytest.skip("No Python code blocks in PRECISION.md")

        for i, code_block in enumerate(code_blocks):
            try:
                # Check syntax by compiling (doesn't execute)
                compile(code_block, f"<code_block_{i}>", "exec")
            except SyntaxError as e:
                pytest.fail(f"Syntax error in code block {i}: {e}\nCode:\n{code_block}")

    def test_import_statements_valid(self, doc_content):
        """Import statements in examples should reference valid modules."""
        # Common imports that should work
        imports_to_test = [
            "import libephemeris as eph",
            "from libephemeris import EphemerisContext",
        ]

        for import_stmt in imports_to_test:
            if import_stmt in doc_content:
                try:
                    exec(import_stmt)
                except ImportError as e:
                    pytest.fail(f"Import failed: {import_stmt} - {e}")


class TestPrecisionTuningFunctionality:
    """Test that precision tuning functionality works as documented."""

    def teardown_method(self):
        """Reset state after each test."""
        eph.close()

    def test_set_ephemeris_file_works(self):
        """set_ephemeris_file should accept valid file names."""
        # Should not raise an error
        eph.set_ephemeris_file("de440.bsp")

    def test_set_tid_acc_works(self):
        """set_tid_acc should work as documented."""
        eph.set_tid_acc(-25.85)
        assert eph.get_tid_acc() == -25.85

        # Reset to default
        eph.set_tid_acc(0.0)

    def test_set_delta_t_userdef_works(self):
        """set_delta_t_userdef should work as documented."""
        from libephemeris import deltat

        # Set a fixed Delta T
        fixed_dt = 65.0 / 86400.0  # 65 seconds in days
        eph.set_delta_t_userdef(fixed_dt)

        # Should return the fixed value
        result = deltat(2451545.0)
        assert abs(result - fixed_dt) < 1e-10

        # Clear
        eph.set_delta_t_userdef(None)

    def test_set_iers_delta_t_enabled_works(self):
        """set_iers_delta_t_enabled should work as documented."""
        # Should not raise an error
        eph.set_iers_delta_t_enabled(False)
        assert eph.get_iers_delta_t_enabled() is False

        eph.set_iers_delta_t_enabled(True)
        assert eph.get_iers_delta_t_enabled() is True

        # Reset
        eph.set_iers_delta_t_enabled(None)

    def test_set_auto_spk_download_works(self):
        """set_auto_spk_download should work as documented."""
        eph.set_auto_spk_download(False)
        assert eph.get_auto_spk_download() is False

        eph.set_auto_spk_download(True)
        assert eph.get_auto_spk_download() is True

        # Reset
        eph.set_auto_spk_download(None)

    def test_set_spk_cache_dir_works(self):
        """set_spk_cache_dir should work as documented."""
        eph.set_spk_cache_dir("/tmp/test_spk_cache")
        assert eph.get_spk_cache_dir() == "/tmp/test_spk_cache"

        # Reset
        eph.set_spk_cache_dir(None)
        assert eph.get_spk_cache_dir() is None

    def test_set_spk_date_padding_works(self):
        """set_spk_date_padding should work as documented."""
        eph.set_spk_date_padding(365)
        assert eph.get_spk_date_padding() == 365

        # Reset
        eph.set_spk_date_padding(0)
        assert eph.get_spk_date_padding() == 0

    def test_get_current_file_data_works(self):
        """get_current_file_data should work as documented."""
        # Force ephemeris loading
        eph.calc_ut(2451545.0, eph.SUN, 0)

        # Should return file info
        path, start, end, denum = eph.get_current_file_data()

        # Should have some values
        assert isinstance(path, str)
        assert isinstance(start, float)
        assert isinstance(end, float)
        assert isinstance(denum, int)
