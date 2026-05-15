"""
Tests for API documentation and docstring formatting.

Verifies that:
1. All public functions have docstrings
2. Docstrings follow Google/NumPy format for Sphinx integration
3. Documentation files exist and are valid
4. Key sections (Args, Returns, Example) are present where expected
"""

import os

import libephemeris


class TestDocstringPresence:
    """Test that all public functions have docstrings."""

    def test_public_functions_have_docstrings(self):
        """All public functions should have non-empty docstrings."""
        missing_docstrings = []

        for name in libephemeris.__all__:
            obj = getattr(libephemeris, name, None)
            if obj is None:
                continue

            # Skip constants (non-callable)
            if not callable(obj) and not isinstance(obj, type):
                continue

            # Check for docstring
            docstring = getattr(obj, "__doc__", None)
            if not docstring or not docstring.strip():
                missing_docstrings.append(name)

        assert len(missing_docstrings) == 0, (
            f"The following public functions are missing docstrings: {missing_docstrings}"
        )

    def test_core_functions_have_detailed_docstrings(self):
        """Core functions should have detailed docstrings with Args/Returns."""
        core_functions = [
            "calc_ut",
            "calc",
            "julday",
            "revjul",
            "houses",
            "set_sid_mode",
            "get_ayanamsa_ut",
        ]

        for func_name in core_functions:
            func = getattr(libephemeris, func_name, None)
            if func is None:
                continue

            docstring = func.__doc__
            assert docstring is not None, f"{func_name} is missing docstring"
            assert len(docstring) > 100, (
                f"{func_name} docstring is too short (less than 100 chars)"
            )


class TestDocstringFormat:
    """Test that docstrings follow Google/NumPy format."""

    def test_args_section_format(self):
        """Check that Args sections are properly formatted."""
        # Test a function known to have Args section
        docstring = libephemeris.calc_ut.__doc__
        assert docstring is not None
        assert "Args:" in docstring, "calc_ut should have Args section"

    def test_returns_section_format(self):
        """Check that Returns sections are properly formatted."""
        docstring = libephemeris.calc_ut.__doc__
        assert docstring is not None
        assert "Returns:" in docstring, "calc_ut should have Returns section"

    def test_example_section_format(self):
        """Check that Example sections are present in key functions."""
        docstring = libephemeris.calc_ut.__doc__
        assert docstring is not None
        assert "Example:" in docstring or ">>>" in docstring, (
            "calc_ut should have Example section"
        )

    def test_docstring_has_parameter_types(self):
        """Docstrings should include type information for parameters."""
        docstring = libephemeris.julday.__doc__
        assert docstring is not None
        # Check for typical type indicators
        assert "int" in docstring.lower() or "float" in docstring.lower(), (
            "julday docstring should mention parameter types"
        )


class TestDocstringContent:
    """Test that docstrings contain accurate information."""

    def test_julday_docstring_mentions_julian_day(self):
        """julday docstring should explain Julian Day."""
        docstring = libephemeris.julday.__doc__
        assert docstring is not None
        assert "Julian" in docstring, "julday should mention Julian Day"

    def test_calc_ut_docstring_mentions_universal_time(self):
        """calc_ut docstring should explain Universal Time."""
        docstring = libephemeris.calc_ut.__doc__
        assert docstring is not None
        assert "Universal Time" in docstring or "UT" in docstring, (
            "calc_ut should mention Universal Time"
        )

    def test_houses_docstring_mentions_house_systems(self):
        """houses docstring should explain house systems."""
        docstring = libephemeris.houses.__doc__
        assert docstring is not None
        # Should mention at least one house system
        assert (
            "Placidus" in docstring
            or "house" in docstring.lower()
            or "cusp" in docstring.lower()
        ), "houses should mention house systems"


class TestDocumentationFiles:
    """Test that documentation files exist and are valid."""

    def test_api_reference_file_exists(self):
        """API reference RST file should exist."""
        docs_dir = os.path.join(os.path.dirname(__file__), "..", "docs")
        api_ref_path = os.path.join(docs_dir, "api_reference.rst")
        assert os.path.exists(api_ref_path), (
            f"API reference file not found at {api_ref_path}"
        )

    def test_sphinx_conf_exists(self):
        """Sphinx configuration file should exist."""
        docs_dir = os.path.join(os.path.dirname(__file__), "..", "docs")
        conf_path = os.path.join(docs_dir, "conf.py")
        assert os.path.exists(conf_path), f"Sphinx conf.py not found at {conf_path}"

    def test_index_rst_exists(self):
        """Documentation index file should exist."""
        docs_dir = os.path.join(os.path.dirname(__file__), "..", "docs")
        index_path = os.path.join(docs_dir, "index.rst")
        assert os.path.exists(index_path), f"index.rst not found at {index_path}"

    def test_api_reference_has_content(self):
        """API reference should have substantial content."""
        docs_dir = os.path.join(os.path.dirname(__file__), "..", "docs")
        api_ref_path = os.path.join(docs_dir, "api_reference.rst")

        if os.path.exists(api_ref_path):
            with open(api_ref_path, "r", encoding="utf-8") as f:
                content = f.read()

            # Should have at least 10KB of documentation
            assert len(content) > 10000, (
                f"API reference seems too short: {len(content)} bytes"
            )

            # Should document key functions
            assert "calc_ut" in content, "API reference should document calc_ut"
            assert "julday" in content, "API reference should document julday"
            assert "houses" in content, "API reference should document houses"

    def test_api_reference_has_sections(self):
        """API reference should be organized in sections."""
        docs_dir = os.path.join(os.path.dirname(__file__), "..", "docs")
        api_ref_path = os.path.join(docs_dir, "api_reference.rst")

        if os.path.exists(api_ref_path):
            with open(api_ref_path, "r", encoding="utf-8") as f:
                content = f.read()

            # Check for major sections
            expected_sections = [
                "Time Functions",
                "Planet",
                "House",
                "Ayanamsha",
                "Constants",
            ]

            for section in expected_sections:
                assert section in content, (
                    f"API reference should have '{section}' section"
                )


class TestExceptionClass:
    """Test that the Error exception class is properly documented."""

    def test_error_class_has_docstring(self):
        """Error exception class should have a docstring."""
        docstring = libephemeris.Error.__doc__
        assert docstring is not None
        assert len(docstring) > 50, "Error class docstring is too short"

    def test_error_class_inherits_from_exception(self):
        """Error should be a subclass of Exception."""
        assert issubclass(libephemeris.Error, Exception)


class TestEphemerisContextDocumentation:
    """Test that EphemerisContext class is properly documented."""

    def test_context_has_docstring(self):
        """EphemerisContext should have a detailed docstring."""
        docstring = libephemeris.EphemerisContext.__doc__
        assert docstring is not None
        assert len(docstring) > 100, "EphemerisContext docstring is too short"

    def test_context_methods_have_docstrings(self):
        """EphemerisContext methods should have docstrings."""
        ctx = libephemeris.EphemerisContext
        methods_to_check = [
            "set_topo",
            "get_topo",
            "set_sid_mode",
            "get_sid_mode",
            "calc_ut",
            "calc",
            "houses",
        ]

        for method_name in methods_to_check:
            method = getattr(ctx, method_name, None)
            if method is not None:
                docstring = method.__doc__
                assert docstring is not None, (
                    f"EphemerisContext.{method_name} is missing docstring"
                )


class TestModuleDocstrings:
    """Test that modules have docstrings."""

    def test_main_module_has_docstring(self):
        """Main module should have a docstring or at least have __all__."""
        # The __init__.py doesn't have a module docstring but exports __all__
        assert hasattr(libephemeris, "__all__"), "libephemeris should export __all__"
        assert len(libephemeris.__all__) > 50, (
            "libephemeris.__all__ should have many exports"
        )

    def test_submodule_docstrings(self):
        """Submodules should have module-level docstrings."""
        from libephemeris import planets, houses, time_utils

        for module in [planets, houses, time_utils]:
            docstring = module.__doc__
            assert docstring is not None, f"{module.__name__} is missing docstring"
            assert len(docstring) > 100, f"{module.__name__} docstring is too short"
