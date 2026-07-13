"""Tests for the live-generated API reference contract.

The API page is intentionally generated from ``libephemeris.__all__`` at
Sphinx build time.  These tests validate that generation contract and the
runtime documentation it consumes; the strict Sphinx gate validates the
rendered artifact itself.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

import libephemeris
from libephemeris.constants import __all__ as CONSTANT_EXPORTS


PROJECT_ROOT = Path(__file__).resolve().parents[1]
API_REF_PATH = PROJECT_ROOT / "docs" / "api_reference.rst"


EXPECTED_GROUPS = {
    "time": {
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
    },
    "planets": {
        "calc_ut",
        "calc",
        "calc_pctr",
        "get_planet_name",
        "nod_aps",
        "get_orbital_elements",
        "pheno",
    },
    "houses": {
        "houses",
        "houses_ex",
        "houses_ex2",
        "houses_armc",
        "house_pos",
        "house_name",
        "gauquelin_sector",
    },
    "ayanamsha": {
        "set_sid_mode",
        "get_ayanamsa_ut",
        "get_ayanamsa",
        "get_ayanamsa_ex",
        "get_ayanamsa_name",
    },
    "tai": {
        "TT_TAI_OFFSET_SECONDS",
        "get_tai_utc_for_jd",
        "utc_to_tai_jd",
        "tai_jd_to_utc",
        "tt_to_tai_jd",
        "tai_to_tt_jd",
    },
    "iers": {
        "set_iers_delta_t_enabled",
        "get_iers_delta_t_enabled",
        "download_iers_finals",
        "download_leap_seconds",
        "get_observed_delta_t",
        "is_observed_delta_t_available",
        "get_delta_t_iers",
        "clear_iers_cache",
    },
    "moons": {
        "register_moon_spk",
        "unregister_moon_spk",
        "list_registered_moons",
        "get_moon_name",
        "is_planetary_moon",
        "calc_moon_position",
        "NAIF_IO",
        "NAIF_EUROPA",
        "NAIF_GANYMEDE",
        "NAIF_TITAN",
        "NAIF_TRITON",
        "NAIF_CHARON",
    },
    "extensions": {
        "EphemerisContext",
        "reset_session",
        "get_elongation_from_sun",
        "get_signed_elongation",
        "is_morning_star",
        "is_evening_star",
        "get_elongation_type",
        "calc_all_arabic_parts",
        "detect_mean_motion_resonance",
        "calc_cupido",
        "calc_transpluto",
        "calc_vulcan",
        "calc_white_moon_position",
        "calc_airmass",
        "calc_extinction_coefficient",
        "calc_extinction_magnitude",
        "calc_twilight_sky_brightness",
        "get_twilight_phase",
        "houses_with_fallback",
        "houses_armc_with_fallback",
        "get_polar_latitude_threshold",
        "PolarCircleError",
        "calc_eclipse_path_width",
        "calc_eclipse_central_line",
        "calc_eclipse_northern_limit",
        "calc_eclipse_southern_limit",
        "get_saros_number",
        "get_inex_number",
        "planet_occult_when_glob",
        "planet_occult_when_loc",
    },
    "performance": {
        "reset_session",
        "set_leb_file",
        "get_leb_reader",
        "set_calc_mode",
        "get_calc_mode",
        "LIBEPHEMERIS_LOG_LEVEL_ENV",
    },
    "constants": {
        "SUN",
        "MOON",
        "MARS",
        "CHIRON",
        "FLG_SPEED",
        "FLG_HELCTR",
        "FLG_TOPOCTR",
        "FLG_SIDEREAL",
        "SIDM_FAGAN_BRADLEY",
        "SIDM_LAHIRI",
        "SIDM_USER",
        "ECL_TOTAL",
        "ECL_ANNULAR",
        "ECL_PARTIAL",
        "ECL_PENUMBRAL",
        "SAROS_CYCLE_DAYS",
        "INEX_CYCLE_DAYS",
    },
}


@pytest.fixture(scope="module")
def api_source() -> str:
    """Return the small source template used to generate the API page."""
    return API_REF_PATH.read_text(encoding="utf-8")


class TestGenerationContract:
    """Validate the source-of-truth relationship between Sphinx and __all__."""

    def test_api_page_uses_live_automodule(self, api_source: str) -> None:
        assert ".. automodule:: libephemeris" in api_source
        assert ":members:" in api_source
        assert ":imported-members:" in api_source
        assert ":undoc-members:" in api_source

    def test_api_page_generates_public_data(self, api_source: str) -> None:
        assert ".. public-data::" in api_source
        assert "libephemeris.__all__" in api_source

    def test_public_export_list_is_unique_and_resolvable(self) -> None:
        exports = libephemeris.__all__
        assert len(exports) == len(set(exports))
        assert len(exports) >= 680
        assert all(hasattr(libephemeris, name) for name in exports)

    def test_all_public_constants_are_package_exports(self) -> None:
        missing = set(CONSTANT_EXPORTS) - set(libephemeris.__all__)
        assert not missing

    @pytest.mark.parametrize(
        "group,names",
        sorted(EXPECTED_GROUPS.items()),
        ids=lambda value: value if isinstance(value, str) else None,
    )
    def test_expected_api_groups_are_exported(
        self, group: str, names: set[str]
    ) -> None:
        del group
        missing = names - set(libephemeris.__all__)
        assert not missing


class TestRuntimeDocumentationQuality:
    """Ensure live autodoc has substantive, inspectable material to render."""

    @staticmethod
    def _callable_docs() -> list[str]:
        return [
            inspect.getdoc(getattr(libephemeris, name)) or ""
            for name in libephemeris.__all__
            if callable(getattr(libephemeris, name))
        ]

    def test_every_public_callable_has_a_docstring(self) -> None:
        missing = [
            name
            for name in libephemeris.__all__
            if callable(getattr(libephemeris, name))
            and not inspect.getdoc(getattr(libephemeris, name))
        ]
        assert not missing

    def test_public_functions_have_runtime_signatures(self) -> None:
        missing = []
        for name in libephemeris.__all__:
            obj = getattr(libephemeris, name)
            if not inspect.isfunction(obj):
                continue
            try:
                inspect.signature(obj)
            except (TypeError, ValueError):
                missing.append(name)
        assert not missing

    def test_public_examples_do_not_import_private_helpers(self) -> None:
        offenders = []
        for name in libephemeris.__all__:
            obj = getattr(libephemeris, name)
            if not callable(obj):
                continue
            doc = inspect.getdoc(obj) or ""
            if any(
                "from libephemeris import" in line
                and any(part.strip().startswith("_") for part in line.split(","))
                for line in doc.splitlines()
            ):
                offenders.append(name)

        assert not offenders

    def test_docstrings_supply_examples_parameters_and_returns(self) -> None:
        docs = self._callable_docs()
        examples = sum(
            ">>>" in doc or "Example:" in doc or "Examples:" in doc for doc in docs
        )
        args = sum("Args:" in doc for doc in docs)
        returns = sum("Returns:" in doc for doc in docs)

        assert examples >= 150
        assert args >= 200
        assert returns >= 200

    def test_context_and_arabic_parts_have_conceptual_docs(self) -> None:
        context_doc = inspect.getdoc(libephemeris.EphemerisContext) or ""
        arabic_doc = inspect.getdoc(libephemeris.calc_all_arabic_parts) or ""

        assert "thread" in context_doc.lower()
        assert "Part of Fortune" in arabic_doc or "Pars_Fortunae" in arabic_doc


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
