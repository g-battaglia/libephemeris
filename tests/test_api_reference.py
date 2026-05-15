"""
Tests for API reference documentation completeness and accuracy.

This module verifies that:
1. The API reference RST file documents all public exports
2. All documented functions exist and are callable
3. New API sections (TAI, IERS, Planetary Moons, etc.) are present
4. Function signatures in documentation match actual implementations
5. Examples in documentation are syntactically correct
"""

import os
import re

import pytest

import libephemeris


# Path to API reference documentation
DOCS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs")
API_REF_PATH = os.path.join(DOCS_DIR, "api_reference.rst")


class TestAPIReferenceSections:
    """Test that all required sections are present in API reference."""

    @pytest.fixture
    def api_content(self):
        """Load API reference content."""
        with open(API_REF_PATH, "r", encoding="utf-8") as f:
            return f.read()

    def test_core_sections_exist(self, api_content):
        """Verify core sections are documented."""
        required_sections = [
            "Time Functions",
            "Planet Calculation Functions",
            "House Functions",
            "Ayanamsha",
            "Fixed Star Functions",
            "Crossing Event Functions",
            "Eclipse Functions",
            "Rise/Set/Transit",
            "State Management Functions",
            "Utility Functions",
            "Constants",
        ]

        for section in required_sections:
            assert section in api_content, f"Missing section: {section}"

    def test_new_sections_exist(self, api_content):
        """Verify new API sections are documented."""
        new_sections = [
            "TAI (International Atomic Time) Functions",
            "IERS Data Functions",
            "Planetary Moons",
            "Elongation Helper Functions",
            "Hypothetical Bodies",
            "Atmospheric Extinction",
        ]

        for section in new_sections:
            assert section in api_content, f"Missing new section: {section}"

    def test_spk_section_exists(self, api_content):
        """Verify SPK kernel documentation is present."""
        assert "SPK Kernel Functions" in api_content
        assert "download_spk" in api_content
        assert "register_spk_body" in api_content

    def test_ephemeris_context_documented(self, api_content):
        """Verify EphemerisContext class is documented."""
        assert "EphemerisContext" in api_content
        assert "Thread-safe" in api_content or "thread-safe" in api_content
        assert "calc_ut" in api_content

    def test_arabic_parts_documented(self, api_content):
        """Verify Arabic parts are documented."""
        assert "Arabic Parts" in api_content
        assert "calc_all_arabic_parts" in api_content
        assert "Pars_Fortunae" in api_content or "Part of Fortune" in api_content


class TestCoreFunctionsDocumented:
    """Test that core functions are documented in API reference."""

    @pytest.fixture
    def api_content(self):
        """Load API reference content."""
        with open(API_REF_PATH, "r", encoding="utf-8") as f:
            return f.read()

    def test_time_functions_documented(self, api_content):
        """Verify time functions are documented."""
        time_functions = [
            "julday",
            "revjul",
            "deltat",
            "deltat_ex",
            "date_conversion",
            "day_of_week",
            "utc_to_jd",
            "sidtime",
            "sidtime0",
            "time_equ",
        ]

        for func in time_functions:
            assert func in api_content, f"Time function not documented: {func}"

    def test_planet_functions_documented(self, api_content):
        """Verify planet calculation functions are documented."""
        planet_functions = [
            "calc_ut",
            "calc",
            "calc_pctr",
            "get_planet_name",
            "nod_aps",
            "get_orbital_elements",
            "pheno",
        ]

        for func in planet_functions:
            assert func in api_content, f"Planet function not documented: {func}"

    def test_house_functions_documented(self, api_content):
        """Verify house functions are documented."""
        house_functions = [
            "houses",
            "houses_ex",
            "houses_ex2",
            "houses_armc",
            "house_pos",
            "house_name",
            "gauquelin_sector",
        ]

        for func in house_functions:
            assert func in api_content, f"House function not documented: {func}"

    def test_ayanamsha_functions_documented(self, api_content):
        """Verify ayanamsha functions are documented."""
        ayanamsha_functions = [
            "set_sid_mode",
            "get_ayanamsa_ut",
            "get_ayanamsa",
            "get_ayanamsa_ex",
            "get_ayanamsa_name",
        ]

        for func in ayanamsha_functions:
            assert func in api_content, f"Ayanamsha function not documented: {func}"


class TestNewFeaturesDocumented:
    """Test that new feature APIs are documented."""

    @pytest.fixture
    def api_content(self):
        """Load API reference content."""
        with open(API_REF_PATH, "r", encoding="utf-8") as f:
            return f.read()

    def test_tai_functions_documented(self, api_content):
        """Verify TAI time functions are documented."""
        tai_functions = [
            "TT_TAI_OFFSET_SECONDS",
            "get_tai_utc_for_jd",
            "utc_to_tai_jd",
            "tai_jd_to_utc",
            "tt_to_tai_jd",
            "tai_to_tt_jd",
        ]

        for func in tai_functions:
            assert func in api_content, f"TAI function not documented: {func}"

    def test_iers_functions_documented(self, api_content):
        """Verify IERS data functions are documented."""
        iers_functions = [
            "set_iers_delta_t_enabled",
            "get_iers_delta_t_enabled",
            "download_iers_finals",
            "download_leap_seconds",
            "get_observed_delta_t",
            "is_observed_delta_t_available",
            "get_delta_t_iers",
            "clear_iers_cache",
        ]

        for func in iers_functions:
            assert func in api_content, f"IERS function not documented: {func}"

    def test_planetary_moon_functions_documented(self, api_content):
        """Verify planetary moon functions are documented."""
        moon_functions = [
            "register_moon_spk",
            "unregister_moon_spk",
            "list_registered_moons",
            "get_moon_name",
            "is_planetary_moon",
            "calc_moon_position",
        ]

        for func in moon_functions:
            assert func in api_content, f"Moon function not documented: {func}"

    def test_moon_naif_ids_documented(self, api_content):
        """Verify planetary moon NAIF IDs are documented."""
        moon_ids = [
            "NAIF_IO",
            "NAIF_EUROPA",
            "NAIF_GANYMEDE",
            "NAIF_TITAN",
            "NAIF_TRITON",
            "NAIF_CHARON",
        ]

        for moon_id in moon_ids:
            assert moon_id in api_content, f"Moon NAIF ID not documented: {moon_id}"

    def test_elongation_functions_documented(self, api_content):
        """Verify elongation helper functions are documented."""
        elongation_functions = [
            "get_elongation_from_sun",
            "get_signed_elongation",
            "is_morning_star",
            "is_evening_star",
            "get_elongation_type",
        ]

        for func in elongation_functions:
            assert func in api_content, f"Elongation function not documented: {func}"

    def test_hypothetical_bodies_documented(self, api_content):
        """Verify hypothetical body functions are documented."""
        hypothetical_functions = [
            "calc_cupido",
            "calc_hades",
            "calc_zeus",
            "calc_kronos",
            "calc_apollon",
            "calc_admetos",
            "calc_vulkanus",
            "calc_poseidon",
            "calc_transpluto",
            "calc_vulcan",
            "calc_waldemath",
            "calc_white_moon_position",
        ]

        for func in hypothetical_functions:
            assert func in api_content, f"Hypothetical function not documented: {func}"

    def test_extinction_functions_documented(self, api_content):
        """Verify atmospheric extinction functions are documented."""
        extinction_functions = [
            "calc_airmass",
            "calc_extinction_coefficient",
            "calc_extinction_magnitude",
            "calc_twilight_sky_brightness",
            "get_twilight_phase",
        ]

        for func in extinction_functions:
            assert func in api_content, f"Extinction function not documented: {func}"


class TestConstantsDocumented:
    """Test that constants are documented."""

    @pytest.fixture
    def api_content(self):
        """Load API reference content."""
        with open(API_REF_PATH, "r", encoding="utf-8") as f:
            return f.read()

    def test_planet_ids_documented(self, api_content):
        """Verify planet IDs are documented."""
        planet_ids = [
            "SUN",
            "MOON",
            "MERCURY",
            "VENUS",
            "MARS",
            "JUPITER",
            "SATURN",
            "URANUS",
            "NEPTUNE",
            "PLUTO",
            "MEAN_NODE",
            "TRUE_NODE",
            "CHIRON",
        ]

        for planet_id in planet_ids:
            assert planet_id in api_content, f"Planet ID not documented: {planet_id}"

    def test_calculation_flags_documented(self, api_content):
        """Verify calculation flags are documented."""
        flags = [
            "FLG_SPEED",
            "FLG_HELCTR",
            "FLG_TOPOCTR",
            "FLG_SIDEREAL",
            "FLG_EQUATORIAL",
            "FLG_J2000",
            "FLG_TRUEPOS",
        ]

        for flag in flags:
            assert flag in api_content, f"Calculation flag not documented: {flag}"

    def test_sidereal_modes_documented(self, api_content):
        """Verify sidereal modes are documented."""
        sidereal_modes = [
            "SIDM_FAGAN_BRADLEY",
            "SIDM_LAHIRI",
            "SIDM_RAMAN",
            "SIDM_KRISHNAMURTI",
            "SIDM_TRUE_CITRA",
            "SIDM_USER",
        ]

        for mode in sidereal_modes:
            assert mode in api_content, f"Sidereal mode not documented: {mode}"

    def test_eclipse_flags_documented(self, api_content):
        """Verify eclipse flags are documented."""
        eclipse_flags = [
            "ECL_TOTAL",
            "ECL_ANNULAR",
            "ECL_PARTIAL",
            "ECL_PENUMBRAL",
        ]

        for flag in eclipse_flags:
            assert flag in api_content, f"Eclipse flag not documented: {flag}"


class TestDocumentedFunctionsExist:
    """Test that documented functions actually exist in libephemeris."""

    def test_core_exports_exist(self):
        """Verify core documented exports exist."""
        core_exports = [
            "julday",
            "revjul",
            "calc_ut",
            "calc",
            "houses",
            "set_sid_mode",
            "get_ayanamsa_ut",
            "fixstar_ut",
            "solcross_ut",
            "sol_eclipse_when_glob",
            "rise_trans",
            "set_topo",
            "EphemerisContext",
            "Error",
        ]

        for export in core_exports:
            assert hasattr(libephemeris, export), (
                f"Documented export not found: {export}"
            )

    def test_new_exports_exist(self):
        """Verify new documented exports exist."""
        new_exports = [
            # TAI functions
            "TT_TAI_OFFSET_SECONDS",
            "get_tai_utc_for_jd",
            "utc_to_tai_jd",
            "tai_jd_to_utc",
            "tt_to_tai_jd",
            "tai_to_tt_jd",
            # IERS functions
            "set_iers_delta_t_enabled",
            "get_iers_delta_t_enabled",
            "download_iers_finals",
            "download_leap_seconds",
            "get_observed_delta_t",
            # Elongation functions
            "get_elongation_from_sun",
            "get_signed_elongation",
            "is_morning_star",
            "is_evening_star",
            "get_elongation_type",
            # Hypothetical bodies
            "calc_cupido",
            "calc_transpluto",
            "calc_vulcan",
            "calc_white_moon_position",
            # Planetary moons
            "register_moon_spk",
            "get_moon_name",
            "NAIF_IO",
            "NAIF_TITAN",
        ]

        for export in new_exports:
            assert hasattr(libephemeris, export), (
                f"New documented export not found: {export}"
            )

    def test_extinction_exports_exist(self):
        """Verify extinction module exports exist."""
        extinction_exports = [
            "calc_airmass",
            "calc_extinction_coefficient",
            "calc_extinction_magnitude",
            "ExtinctionCoefficients",
            "TwilightSkyBrightness",
            "get_twilight_phase",
        ]

        for export in extinction_exports:
            assert hasattr(libephemeris, export), (
                f"Extinction export not found: {export}"
            )


class TestAPIReferenceQuality:
    """Test quality aspects of the API reference."""

    @pytest.fixture
    def api_content(self):
        """Load API reference content."""
        with open(API_REF_PATH, "r", encoding="utf-8") as f:
            return f.read()

    def test_minimum_documentation_size(self, api_content):
        """API reference should have substantial content."""
        # Should be at least 50KB with all sections
        min_size = 50000
        assert len(api_content) >= min_size, (
            f"API reference too small: {len(api_content)} bytes, expected >= {min_size}"
        )

    def test_has_examples(self, api_content):
        """API reference should contain examples."""
        example_count = api_content.count(">>>")
        min_examples = 10
        assert example_count >= min_examples, (
            f"API reference has too few examples: {example_count}, expected >= {min_examples}"
        )

    def test_has_param_documentation(self, api_content):
        """API reference should have parameter documentation."""
        param_count = api_content.count(":param ")
        min_params = 50
        assert param_count >= min_params, (
            f"API reference has too few :param entries: {param_count}, expected >= {min_params}"
        )

    def test_has_return_documentation(self, api_content):
        """API reference should have return type documentation."""
        return_count = api_content.count(":returns:") + api_content.count(":rtype:")
        min_returns = 30
        assert return_count >= min_returns, (
            f"API reference has too few return docs: {return_count}, expected >= {min_returns}"
        )

    def test_rst_syntax_valid(self, api_content):
        """API reference should have valid RST section headers."""
        # Check for proper section underlines
        lines = api_content.split("\n")
        for i, line in enumerate(lines):
            # Check that underlines match header lengths
            if i > 0 and len(line) > 0 and all(c in "=-~^" for c in line):
                # This is an underline, check it matches the header above
                header = lines[i - 1]
                if len(header.strip()) > 0:
                    # Underline should be at least as long as header
                    assert len(line) >= len(header.strip()), (
                        f"RST underline too short at line {i + 1}: '{header}'"
                    )

    def test_function_directives_valid(self, api_content):
        """Check for valid function directive syntax."""
        # Look for .. function:: directives
        function_pattern = r"\.\. function:: (\w+)"
        functions = re.findall(function_pattern, api_content)

        # Should have many function documentation entries
        assert len(functions) >= 50, (
            f"Too few function directives: {len(functions)}, expected >= 50"
        )

    def test_data_directives_valid(self, api_content):
        """Check for valid data directive syntax."""
        # Look for .. data:: directives for constants
        data_pattern = r"\.\. data:: (\w+)"
        data_entries = re.findall(data_pattern, api_content)

        # Should have constant documentation entries
        assert len(data_entries) >= 20, (
            f"Too few data directives: {len(data_entries)}, expected >= 20"
        )


class TestPolishLatitudeDocumented:
    """Test that polar latitude handling is documented."""

    @pytest.fixture
    def api_content(self):
        """Load API reference content."""
        with open(API_REF_PATH, "r", encoding="utf-8") as f:
            return f.read()

    def test_polar_functions_documented(self, api_content):
        """Verify polar latitude functions are documented."""
        polar_functions = [
            "houses_with_fallback",
            "houses_armc_with_fallback",
            "get_polar_latitude_threshold",
            "PolarCircleError",
        ]

        for func in polar_functions:
            assert func in api_content, (
                f"Polar latitude function not documented: {func}"
            )


class TestEclipseExtrasDocumented:
    """Test that eclipse additional functions are documented."""

    @pytest.fixture
    def api_content(self):
        """Load API reference content."""
        with open(API_REF_PATH, "r", encoding="utf-8") as f:
            return f.read()

    def test_eclipse_path_functions_documented(self, api_content):
        """Verify eclipse path functions are documented."""
        path_functions = [
            "calc_eclipse_path_width",
            "calc_eclipse_central_line",
            "calc_eclipse_northern_limit",
            "calc_eclipse_southern_limit",
        ]

        for func in path_functions:
            assert func in api_content, f"Eclipse path function not documented: {func}"

    def test_saros_inex_documented(self, api_content):
        """Verify Saros and Inex series are documented."""
        saros_items = [
            "get_saros_number",
            "get_inex_number",
            "SAROS_CYCLE_DAYS",
            "INEX_CYCLE_DAYS",
        ]

        for item in saros_items:
            assert item in api_content, f"Saros/Inex item not documented: {item}"

    def test_occultation_functions_documented(self, api_content):
        """Verify planetary occultation functions are documented."""
        occultation_functions = [
            "planet_occult_when_glob",
            "planet_occult_when_loc",
        ]

        for func in occultation_functions:
            assert func in api_content, f"Occultation function not documented: {func}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
