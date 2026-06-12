"""Test commands: unit tests, comparison tests, lunar, LEB, LEB2, Horizons, coverage.

Replaces 48 poe test tasks with clearly named hierarchical commands.
Every command name makes it explicit WHAT is being tested and HOW (which backend).

Naming convention:
  - ``skyfield``      = unit tests using Skyfield/DE440 backend
  - ``leb-backend``   = unit tests using LEB precomputed backend
  - ``compare``       = libephemeris vs pyswisseph comparison
  - ``lunar``         = lunar module (nodes, apsides, Lilith)
  - ``leb-format``    = LEB binary format (reader/writer/precision)
  - ``leb2-format``   = LEB2 compressed format
  - ``horizons``      = Horizons API precision tests
  - ``coverage``      = test coverage reports
"""

from __future__ import annotations

import os
import subprocess
import sys

import click


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _pytest(args: list[str], env: dict[str, str] | None = None) -> None:
    """Run pytest with the given arguments, optionally setting env vars."""
    run_env: dict[str, str] | None = None
    if env:
        run_env = {**os.environ, **env}
    # sys.executable -m pytest: always the venv's pytest, not whatever
    # 'pytest' resolves to on PATH
    sys.exit(subprocess.call([sys.executable, "-m", "pytest", *args], env=run_env))


def _python(args: list[str], env: dict[str, str] | None = None) -> None:
    """Run a python script with the given arguments."""
    run_env: dict[str, str] | None = None
    if env:
        run_env = {**os.environ, **env}
    sys.exit(subprocess.call([sys.executable, *args], env=run_env))


# ---------------------------------------------------------------------------
# Root test group
# ---------------------------------------------------------------------------


@click.group(
    "test",
    short_help="Run tests: unit, comparison, lunar, LEB, Horizons, coverage.",
    help="Run tests: unit, comparison, lunar, LEB format, Horizons, coverage.\n\nEach subgroup specifies WHAT is tested and WHICH backend is used.\nUse TAB completion to explore: leph test <TAB>\n\nRecommended for daily development:\n\n  leph test skyfield essential       # ~490 tests, ~20s\n  leph test leb-backend unit-fast    # ~5890 tests, ~1 min\n\nNever run the full test suite -- always pick a targeted subcommand.",
)
def test_group() -> None:
    """Test runner with clear backend/target naming."""


# ===========================================================================
# leph test skyfield — Unit tests via Skyfield/DE440 backend
# ===========================================================================

# Test file lists for essential and smoke subsets
_ESSENTIAL_FILES = [
    "tests/test_planets/test_central_difference_velocity.py",
    "tests/test_planets/test_nod_aps.py",
    "tests/test_houses/test_house_algorithms.py",
    "tests/test_lunar_eclipse.py",
    "tests/test_lunar/test_interpolated_apogee.py",
    "tests/test_fixed_stars/test_fixstar2.py",
    "tests/test_asteroid_by_number.py",
    "tests/test_hypothetical.py",
    "tests/test_heliacal.py",
    "tests/test_rise_trans.py",
    "tests/test_arabic_parts.py",
    "tests/test_planetary_moons.py",
    "tests/test_time/test_date_conversion.py",
    "tests/test_exception_hierarchy.py",
    "tests/test_close.py",
    "tests/test_sidereal/test_ayanamsa_precision.py",
    "tests/test_deg_midp.py",
    "tests/test_get_planet_name.py",
]

_SMOKE_FILES = _ESSENTIAL_FILES + [
    "tests/test_cs2degstr.py",
    "tests/test_cache.py",
    "tests/test_spk.py",
    "tests/test_spk_auto.py",
    "tests/test_ephemeris_config.py",
    "tests/test_extinction.py",
    "tests/test_logging_config.py",
    "tests/test_dotenv.py",
    "tests/test_type_safety.py",
    "tests/test_golden_regression.py",
    "tests/test_coordinate_validation.py",
    "tests/test_profiling.py",
    "tests/test_retrograde_stations.py",
    "tests/test_edge_cases/test_boundaries.py",
    "tests/test_precision/test_erfa_nutation.py",
    "tests/test_leb/test_leb_format.py",
    "tests/test_leb/test_leb_reader.py",
    "tests/test_leb/test_fast_calc.py",
    "tests/test_lunar/test_mean_lilith_enhanced.py",
    "tests/test_lunar/test_true_node_terms.py",
    "tests/test_time/test_utc_to_jd.py",
    "tests/test_time/test_iers_delta_t.py",
    "tests/test_time/test_deltat_smh2016.py",
]


@click.group(
    "skyfield",
    short_help="Unit tests via Skyfield/DE440 backend (essential, smoke, unit, all).",
    help="Unit tests that calculate positions via the Skyfield/DE440 backend.\n\n"
    "This is the default calculation engine: it loads NASA JPL DE440 binary\n"
    "kernels via the Skyfield library and computes positions in real time.\n"
    "Requires DE440 data to be downloaded first (leph download spk-medium).\n\n"
    "Start here:\n\n"
    "  leph test skyfield essential   # ~490 tests, ~20s, parallel\n"
    "  leph test skyfield smoke       # ~1460 tests, ~30s, parallel",
)
def skyfield_group() -> None:
    """Unit tests using the Skyfield/DE440 backend."""


@skyfield_group.command(
    "all",
    short_help="Run all tests (unit + compare) sequentially, no @slow.",
)
def skyfield_all() -> None:
    """Run all unit tests sequentially, excluding @slow markers.

    Runs the tests/ tree (pytest.ini testpaths) — compare_scripts/tests/
    is NOT collected here; use the 'compare' commands for those.
    Sequential execution means slower but easier to debug failures.
    """
    _pytest(["-m", "not slow"])


@skyfield_group.command(
    "all-full",
    short_help="Run ALL tests including @slow, sequentially.",
)
def skyfield_all_full() -> None:
    """Run ALL tests including @slow-marked tests, sequentially.

    Includes iterative search tests, hypothesis property tests, and other
    slow tests that are normally skipped. Very slow (~10+ min).
    """
    _pytest([])


@skyfield_group.command(
    "all-fast",
    short_help="Run all tests excluding @slow in parallel (-n auto).",
)
def skyfield_all_fast() -> None:
    """Run all tests excluding @slow in parallel across all CPU cores (-n auto).

    Same scope as 'all' but uses pytest-xdist for parallel execution.
    Much faster (~2 min) but harder to read output on failures.
    """
    _pytest(["-n", "auto", "-m", "not slow"])


@skyfield_group.command(
    "all-full-fast",
    short_help="Run ALL tests including @slow in parallel.",
)
def skyfield_all_full_fast() -> None:
    """Run ALL tests including @slow in parallel. Still slow due to iterative tests."""
    _pytest(["-n", "auto"])


@skyfield_group.command(
    short_help="Run all tests with compact dots-only output (CI mode).",
)
def progress() -> None:
    """Run all tests with compact dots-only output (good for CI pipelines).

    Same scope as 'all' (no @slow) but uses --no-header -q for minimal output.
    """
    _pytest(["--no-header", "-q", "-m", "not slow"])


@skyfield_group.command(
    short_help="Quick sanity check: 1 test per module, parallel (~490 tests, ~20s).",
)
def essential() -> None:
    """Quick sanity check: 1 test file per major module, in parallel (~490 tests, ~20s).

    Covers planets, houses, eclipses, fixed stars, asteroids, hypothetical
    bodies, heliacal events, rise/transit, Arabic parts, planetary moons,
    date conversion, sidereal, and more. Best for a fast "did I break anything?" check.
    """
    _pytest(["-n", "auto", "-m", "not slow", *_ESSENTIAL_FILES])


@skyfield_group.command(
    short_help="Broader sanity check: 1+ test files per module, parallel (~1460 tests).",
)
def smoke() -> None:
    """Broader sanity check: 1+ test files per module, in parallel (~1460 tests, ~30s).

    Everything in 'essential' plus cache, SPK loader, ephemeris config,
    extinction, logging, packaging, dotenv, type safety, golden regression,
    coordinate validation, profiling, retrograde stations, edge cases,
    ERFA nutation, LEB format/reader, Lilith, true node, UTC, IERS delta-T.
    """
    _pytest(["-n", "auto", "-m", "not slow", *_SMOKE_FILES])


@skyfield_group.command(
    short_help="Run unit tests from tests/ sequentially, verbose, no @slow.",
)
def unit() -> None:
    """Run all unit tests from the tests/ directory sequentially with verbose output.

    Excludes @slow-marked tests. ~5890 tests, takes ~2 min.
    Does NOT include compare_scripts/ tests (use 'all' for those).
    """
    _pytest(["tests/", "-v", "-m", "not slow"])


@skyfield_group.command(
    "unit-full",
    short_help="Run ALL unit tests including @slow, sequentially.",
)
def unit_full() -> None:
    """Run ALL unit tests including @slow markers, sequentially (~10+ min).

    Includes iterative searches, hypothesis property-based tests, and other
    long-running tests that are excluded by default.
    """
    _pytest(["tests/", "-v"])


@skyfield_group.command(
    "unit-fast",
    short_help="Run unit tests in parallel (-n auto), no @slow (~1 min).",
)
def unit_fast() -> None:
    """Run unit tests in parallel across all CPU cores (-n auto), no @slow (~1 min).

    Same scope as 'unit' but parallelized for speed.
    """
    _pytest(["tests/", "-n", "auto", "-m", "not slow"])


test_group.add_command(skyfield_group)


# ===========================================================================
# leph test leb-backend — Unit tests using LEB precomputed backend
# ===========================================================================

_LEB_ENV = {"LIBEPHEMERIS_LEB": "data/leb/ephemeris_medium.leb"}


@click.group(
    "leb-backend",
    short_help="Unit tests via LEB precomputed Chebyshev backend (~14x faster).",
    help="Unit tests that calculate positions via the LEB precomputed Chebyshev backend.\n\n"
    "Instead of computing positions from DE440 kernels in real time (Skyfield),\n"
    "these tests use precomputed Chebyshev polynomial approximations stored in\n"
    ".leb files. ~14x faster than Skyfield at runtime.\n\n"
    "Requires: data/leb/ephemeris_medium.leb (generate with 'leph leb generate medium groups').\n\n"
    "Recommended for everyday development:\n\n"
    "  leph test leb-backend unit-fast   # ~5890 tests, ~1 min, parallel",
)
def leb_backend_group() -> None:
    """Unit tests using LEB as the calculation backend."""


@leb_backend_group.command(
    "essential",
    short_help="Quick sanity check in LEB mode: ~490 tests, parallel (~20s).",
)
def leb_essential() -> None:
    """Quick sanity check in LEB mode: 1 test file per major module, parallel (~490 tests, ~20s).

    Same file set as 'leph test skyfield essential' but positions are read from
    precomputed Chebyshev polynomials instead of computed from DE440.
    Requires data/leb/ephemeris_medium.leb.
    """
    _pytest(
        ["-n", "auto", "-m", "not slow", "--calc-mode", "leb", *_ESSENTIAL_FILES],
        env=_LEB_ENV,
    )


@leb_backend_group.command(
    "unit",
    short_help="Run unit tests sequentially in LEB mode, verbose output.",
)
def leb_unit() -> None:
    """Run unit tests sequentially in LEB mode with verbose output.

    Same test suite as 'leph test skyfield unit' but positions are read from
    precomputed Chebyshev polynomials instead of computed from DE440.
    Requires data/leb/ephemeris_medium.leb to exist.
    """
    _pytest(["tests/", "-v", "-m", "not slow", "--calc-mode", "leb"], env=_LEB_ENV)


@leb_backend_group.command(
    "unit-full",
    short_help="Run ALL unit tests in LEB mode including @slow, sequentially.",
)
def leb_unit_full() -> None:
    """Run ALL unit tests in LEB mode including @slow markers, sequentially.

    Very slow. Use for thorough pre-release validation.
    """
    _pytest(["tests/", "-v", "--calc-mode", "leb"], env=_LEB_ENV)


@leb_backend_group.command(
    "unit-fast",
    short_help="RECOMMENDED: run unit tests in LEB mode, parallel (~1 min).",
)
def leb_unit_fast() -> None:
    """RECOMMENDED: run unit tests in LEB mode, parallel across CPU cores (~1 min).

    Best balance of speed and coverage for everyday development.
    Runs ~5890 tests in ~1 minute using pytest-xdist parallelism.
    Requires data/leb/ephemeris_medium.leb (generate with 'leph leb generate medium groups').
    """
    _pytest(
        ["tests/", "-n", "auto", "-m", "not slow", "--calc-mode", "leb"], env=_LEB_ENV
    )


test_group.add_command(leb_backend_group)


# ===========================================================================
# leph test compare — libephemeris vs pyswisseph comparison
# ===========================================================================


@click.group(
    "compare",
    short_help="Compare libephemeris vs pyswisseph (Skyfield, LEB, or Horizons backend).",
    help="Compare libephemeris output against pyswisseph (the C reference implementation).\n\n"
    "These tests call the same function with the same inputs on both libraries\n"
    "and assert that the results match within tolerance. This is the primary\n"
    "way to verify that libephemeris is a faithful reimplementation.\n\n"
    "Requires: pyswisseph to be installed (pip install pyswisseph).\n"
    "Tests live in compare_scripts/tests/.\n\n"
    "You can compare via different backends to validate each calculation path:\n\n"
    "  leph test compare skyfield          # via Skyfield/DE440\n"
    "  leph test compare leb-backend       # via LEB precomputed\n"
    "  leph test compare horizons-backend  # via NASA Horizons API",
)
def compare_group() -> None:
    """Comparison tests against pyswisseph."""


@compare_group.command(
    short_help="Compare via Skyfield/DE440, excluding @slow.",
)
def skyfield() -> None:
    """Compare via Skyfield/DE440 backend, excluding @slow tests.

    Standard comparison: libephemeris computes positions using Skyfield,
    pyswisseph computes using its built-in Swiss Ephemeris. Results are
    compared for all supported bodies and flags.
    """
    _pytest(["compare_scripts/tests/", "-v", "-m", "not slow"])


@compare_group.command(
    "skyfield-full",
    short_help="Compare via Skyfield/DE440, ALL tests including @slow.",
)
def skyfield_full() -> None:
    """Compare via Skyfield/DE440 backend, ALL tests including @slow."""
    _pytest(["compare_scripts/tests/", "-v"])


@compare_group.command(
    "skyfield-fast",
    short_help="Compare via Skyfield/DE440, parallel (-n auto).",
)
def skyfield_fast() -> None:
    """Compare via Skyfield/DE440 backend, parallel across CPU cores (-n auto)."""
    _pytest(["compare_scripts/tests/", "-n", "auto", "-m", "not slow"])


@compare_group.command(
    "skyfield-jpl",
    short_help="Compare via Skyfield with explicit env var, no @slow.",
)
def skyfield_jpl() -> None:
    """Compare via Skyfield with explicit LIBEPHEMERIS_COMPARE_MODE=skyfield env var.

    Same calculation path as 'skyfield' but forces the env var for explicitness.
    Useful when testing compare mode dispatch logic.
    """
    _pytest(
        ["compare_scripts/tests/", "-v", "-m", "not slow"],
        env={"LIBEPHEMERIS_COMPARE_MODE": "skyfield"},
    )


@compare_group.command(
    "skyfield-jpl-full",
    short_help="Compare via Skyfield with explicit env var, ALL tests.",
)
def skyfield_jpl_full() -> None:
    """Compare via Skyfield with explicit env var, ALL tests including @slow."""
    _pytest(
        ["compare_scripts/tests/", "-v"],
        env={"LIBEPHEMERIS_COMPARE_MODE": "skyfield"},
    )


@compare_group.command(
    "leb-backend",
    short_help="Compare via LEB backend vs pyswisseph, no @slow.",
)
def leb_backend() -> None:
    """Compare via LEB precomputed backend vs pyswisseph, excluding @slow.

    Validates that LEB Chebyshev approximations produce the same results
    as pyswisseph. This tests end-to-end accuracy of the LEB pipeline.
    Requires data/leb/ephemeris_medium.leb.
    """
    _pytest(
        ["compare_scripts/tests/", "-v", "-m", "not slow"],
        env={
            "LIBEPHEMERIS_COMPARE_MODE": "leb",
            "LIBEPHEMERIS_LEB": "data/leb/ephemeris_medium.leb",
        },
    )


@compare_group.command(
    "leb-backend-full",
    short_help="Compare via LEB backend vs pyswisseph, ALL tests.",
)
def leb_backend_full() -> None:
    """Compare via LEB backend vs pyswisseph, ALL tests including @slow."""
    _pytest(
        ["compare_scripts/tests/", "-v"],
        env={
            "LIBEPHEMERIS_COMPARE_MODE": "leb",
            "LIBEPHEMERIS_LEB": "data/leb/ephemeris_medium.leb",
        },
    )


@compare_group.command(
    "horizons-backend",
    short_help="Compare via Horizons API vs pyswisseph (requires internet).",
)
def horizons_backend() -> None:
    """Compare via NASA JPL Horizons API vs pyswisseph (requires internet).

    Validates that the Horizons live-query backend produces results
    consistent with pyswisseph. Slower due to HTTP requests.
    """
    _pytest(
        ["compare_scripts/tests/", "-v", "-m", "not slow"],
        env={"LIBEPHEMERIS_COMPARE_MODE": "horizons"},
    )


@compare_group.command(
    "horizons-backend-full",
    short_help="Compare via Horizons API vs pyswisseph, ALL tests.",
)
def horizons_backend_full() -> None:
    """Compare via Horizons API vs pyswisseph, ALL tests including @slow.

    Full comparison including slow iterative tests. Requires internet.
    """
    _pytest(
        ["compare_scripts/tests/", "-v"],
        env={"LIBEPHEMERIS_COMPARE_MODE": "horizons"},
    )


@compare_group.command(
    "houses",
    short_help="Compare all 25 house systems (26 codes) vs pyswisseph with precision report.",
)
def compare_houses() -> None:
    """Compare all 25 house systems (26 codes) against pyswisseph.

    Tests cusps, all 8 ASCMC values, houses_armc, house_pos, and sidereal
    houses across 8 locations and 12 dates. Prints per-system precision report
    showing max error for each value.

    \b
    Coverage:
      - houses: 24 systems × 8 locations × 12 dates (cusps + 8 ASCMC)
      - houses_armc: 24 systems × 12 ARMC × 4 latitudes
      - house_pos: 12 systems × 18 longitudes × 3 body latitudes
      - houses_ex sidereal: 8 systems × 3 ayanamshas × 6 dates
    """
    _pytest([
        "compare_scripts/tests/test_houses_comprehensive.py",
        "-v", "-s",
    ])


test_group.add_command(compare_group)


# ===========================================================================
# leph test lunar — Lunar module tests
# ===========================================================================

_LILITH_FILES = [
    "tests/test_lunar/test_mean_lilith_enhanced.py",
    "tests/test_lunar/test_true_lilith_precision.py",
    "tests/test_lunar/test_true_lilith_annual_equation.py",
    "tests/test_lunar/test_true_lilith_evection_secondary.py",
    "tests/test_lunar/test_true_lilith_parallactic_inequality.py",
    "tests/test_lunar/test_true_lilith_reduction_to_ecliptic.py",
    "tests/test_lunar/test_true_lilith_solar_perturbation.py",
    "tests/test_lunar/test_true_lilith_variation.py",
]


@click.group(
    "lunar",
    short_help="Lunar module: nodes, apsides, Lilith, ELP2000 perturbation coefficients.",
    help="Tests for the lunar module: mean/true nodes, perigee/apogee apsides,\nmean and true Lilith (Black Moon), ELP2000 perturbation coefficients.\n\nThese tests validate the analytical lunar element calculations and\nthe interpolated perigee/apogee pipeline.\n\n  leph test lunar perigee   # ELP2000 coefficients + interpolated perigee\n  leph test lunar lilith    # Mean + true Lilith precision (8 test files)",
)
def lunar_group() -> None:
    """Lunar-specific test suites."""


@lunar_group.command(
    "all",
    short_help="Run all lunar tests (nodes, Lilith, perigee, apogee), no @slow.",
)
def lunar_all() -> None:
    """Run all lunar module tests (nodes, Lilith, perigee, apogee), excluding @slow."""
    _pytest(["tests/test_lunar/", "-v", "-m", "not slow"])


@lunar_group.command(
    short_help="Test ELP2000 perigee coefficients + interpolated perigee.",
)
def perigee() -> None:
    """Test ELP2000 perigee perturbation coefficients, interpolated and osculating perigee.

    Runs 3 test files: elp2000_perigee_perturbations, interpolated_perigee,
    osculating_perigee. Run this after calibrating perigee coefficients.
    """
    _pytest(
        [
            "-v",
            "tests/test_lunar/test_elp2000_perigee_perturbations.py",
            "tests/test_lunar/test_interpolated_perigee.py",
            "tests/test_lunar/test_osculating_perigee.py",
        ]
    )


@lunar_group.command(
    short_help="Test ELP2000 apogee coefficients and interpolated apogee.",
)
def apogee() -> None:
    """Test ELP2000 apogee perturbation coefficients and interpolated apogee.

    Runs 2 test files: elp2000_apogee_perturbations, interpolated_apogee.
    """
    _pytest(
        [
            "-v",
            "tests/test_lunar/test_elp2000_apogee_perturbations.py",
            "tests/test_lunar/test_interpolated_apogee.py",
        ]
    )


@lunar_group.command(
    short_help="Test mean + true Lilith precision across 8 test files.",
)
def lilith() -> None:
    """Test mean Lilith and true Lilith precision across 8 test files.

    Covers: mean_lilith_enhanced, true_lilith_precision, annual_equation,
    evection_secondary, parallactic_inequality, reduction_to_ecliptic,
    solar_perturbation, variation.
    """
    _pytest(["-v", *_LILITH_FILES])


test_group.add_command(lunar_group)


# ===========================================================================
# leph test leb-format — LEB binary format tests
# ===========================================================================

_LEB_COMPARE_IGNORES = [
    "--ignore=tests/test_leb/compare/base",
    "--ignore=tests/test_leb/compare/medium",
    "--ignore=tests/test_leb/compare/extended",
    "--ignore=tests/test_leb/compare/crosstier",
]


@click.group(
    "leb-format",
    short_help="LEB binary format internals: reader, writer, Chebyshev precision, vs Skyfield.",
    help="Tests for the LEB binary ephemeris FORMAT itself (reader, writer, precision).\n\n"
    "Unlike 'leb-backend' (which uses LEB as a calculation backend to run the\n"
    "standard test suite), these tests verify the LEB file format internals:\n"
    "binary read/write, Chebyshev polynomial precision, and accuracy compared\n"
    "to live Skyfield calculations across all tiers.\n\n"
    "  leph test leb-format all                    # Reader/writer/fast_calc tests\n"
    "  leph test leb-format precision              # Chebyshev error sweep, all tiers\n"
    "  leph test leb-format vs-skyfield medium     # LEB vs Skyfield accuracy, medium tier",
)
def leb_format_group() -> None:
    """Tests for the LEB binary ephemeris format itself (not using it as a backend)."""


@leb_format_group.command(
    "all", short_help="Reader, writer, Chebyshev reader, fast_calc pipeline (no @slow)."
)
def leb_format_all() -> None:
    """Run all LEB module tests: binary format read/write, Chebyshev reader, fast_calc pipeline.

    Excludes @slow. Tests the LEB file format internals without comparing
    against Skyfield (use 'vs-skyfield' subgroup for that).
    """
    _pytest(["tests/test_leb/", "-v", "-m", "not slow"])


@leb_format_group.command(
    "precision",
    short_help="Chebyshev approximation error sweep for all bodies, all tiers.",
)
def leb_precision() -> None:
    """Verify Chebyshev polynomial approximation error for all bodies across all tiers.

    Sweeps through every body in each tier and measures the maximum error
    between the LEB polynomial evaluation and the Skyfield reference value.
    Runs all tiers (base, medium, extended).
    """
    _pytest(["tests/test_leb/test_leb_precision.py", "-v"])


@leb_format_group.command(
    "precision-quick", short_help="Chebyshev error sweep for medium tier only (faster)."
)
def leb_precision_quick() -> None:
    """Chebyshev precision sweep for medium tier only (faster subset of 'precision')."""
    _pytest(["tests/test_leb/test_leb_precision.py", "-v", "-k", "medium"])


# --- vs-skyfield subgroup ---


@click.group(
    "vs-skyfield",
    short_help="Compare LEB Chebyshev output vs live Skyfield, per tier.",
    help="Compare LEB Chebyshev output against live Skyfield calculations.\n\n"
    "These tests evaluate the same body at the same Julian date using both\n"
    "the LEB reader and Skyfield, then assert the difference is within tolerance.\n"
    "Available per-tier: base, medium, extended, plus cross-tier consistency checks.",
)
def vs_skyfield_group() -> None:
    """Verify LEB output matches Skyfield reference."""


@vs_skyfield_group.command(
    short_help="Medium tier legacy tests from compare/ root (original layout).",
)
def legacy() -> None:
    """Medium tier legacy tests from the compare/ root directory (original flat layout)."""
    _pytest(
        ["tests/test_leb/compare/", "-v", "-m", "leb_compare", *_LEB_COMPARE_IGNORES]
    )


@vs_skyfield_group.command(
    "legacy-quick",
    short_help="Medium tier legacy tests, excluding @slow.",
)
def legacy_quick() -> None:
    """Medium tier legacy tests, excluding @slow markers."""
    _pytest(
        [
            "tests/test_leb/compare/",
            "-v",
            "-m",
            "leb_compare and not slow",
            *_LEB_COMPARE_IGNORES,
        ]
    )


@vs_skyfield_group.command(
    short_help="Base tier (1850-2150): positions, speeds, sidereal, equatorial.",
)
def base() -> None:
    """Base tier (de440s, 1850-2150): compare positions, speeds, sidereal, equatorial."""
    _pytest(["tests/test_leb/compare/base/", "-v", "-m", "leb_compare_base"])


@vs_skyfield_group.command(
    "base-quick",
    short_help="Base tier, excluding @slow markers.",
)
def base_quick() -> None:
    """Base tier, excluding @slow markers."""
    _pytest(
        ["tests/test_leb/compare/base/", "-v", "-m", "leb_compare_base and not slow"]
    )


@vs_skyfield_group.command(
    short_help="Medium tier (1550-2650): positions, speeds, sidereal, equatorial.",
)
def medium() -> None:
    """Medium tier (de440, 1550-2650): compare positions, speeds, sidereal, equatorial."""
    _pytest(["tests/test_leb/compare/medium/", "-v", "-m", "leb_compare_medium"])


@vs_skyfield_group.command(
    "medium-quick",
    short_help="Medium tier, excluding @slow markers (~2 min).",
)
def medium_quick() -> None:
    """Medium tier, excluding @slow markers (~2 min)."""
    _pytest(
        [
            "tests/test_leb/compare/medium/",
            "-v",
            "-m",
            "leb_compare_medium and not slow",
        ]
    )


@vs_skyfield_group.command(
    short_help="Extended tier (-5000 to 5000): positions, speeds, sidereal.",
)
def extended() -> None:
    """Extended tier (de441, -5000 to 5000): compare positions, speeds, sidereal."""
    _pytest(["tests/test_leb/compare/extended/", "-v", "-m", "leb_compare_extended"])


@vs_skyfield_group.command(
    "extended-quick",
    short_help="Extended tier, excluding @slow markers.",
)
def extended_quick() -> None:
    """Extended tier, excluding @slow markers."""
    _pytest(
        [
            "tests/test_leb/compare/extended/",
            "-v",
            "-m",
            "leb_compare_extended and not slow",
        ]
    )


@vs_skyfield_group.command(
    short_help="Cross-tier consistency: base/medium/extended agree on overlapping dates.",
)
def crosstier() -> None:
    """Cross-tier consistency: verify that base, medium, and extended agree on overlapping date ranges."""
    _pytest(["tests/test_leb/compare/crosstier/", "-v", "-m", "leb_compare_crosstier"])


@vs_skyfield_group.command(
    "all",
    short_help="Run ALL tier comparisons + cross-tier checks (comprehensive).",
)
def vs_skyfield_all() -> None:
    """Run ALL tier comparisons + cross-tier checks combined (comprehensive, slow)."""
    _pytest(
        [
            "tests/test_leb/compare/",
            "-v",
            "-m",
            "leb_compare or leb_compare_base or leb_compare_medium "
            "or leb_compare_extended or leb_compare_crosstier",
        ]
    )


leb_format_group.add_command(vs_skyfield_group)
test_group.add_command(leb_format_group)


# ===========================================================================
# leph test leb2-format — LEB2 compressed format tests
# ===========================================================================


@click.group(
    "leb2-format",
    short_help="LEB2 compressed format: compression roundtrip, precision vs LEB1.",
    help="Tests for the LEB2 compressed ephemeris format.\n\n"
    "LEB2 is a lossy-compressed version of LEB1 that achieves 4-10x smaller\n"
    'files while maintaining <0.001" precision. These tests verify the\n'
    "compression/decompression roundtrip and measure precision loss vs LEB1.\n\n"
    "  leph test leb2-format all             # Compression + reader unit tests (27)\n"
    "  leph test leb2-format precision-base  # 31 bodies x 6 flags x 200 dates (~15s)",
)
def leb2_format_group() -> None:
    """Tests for the LEB2 compressed format."""


@leb2_format_group.command(
    "all",
    short_help="Run LEB2 compression and reader unit tests (27 tests).",
)
def leb2_format_all() -> None:
    """Run LEB2 compression and reader unit tests (27 tests).

    Tests compress/decompress roundtrip, mantissa bit computation,
    LEB2Reader lazy decompression, and format auto-detection.
    """
    _pytest(
        [
            "tests/test_leb/test_leb_compression.py",
            "tests/test_leb/test_leb2_reader.py",
            "-v",
        ]
    )


@leb2_format_group.command(
    "precision-base",
    short_help="Measure LEB2 vs LEB1 precision for base tier (~15s).",
)
def leb2_precision_base() -> None:
    """Measure LEB2 vs LEB1 precision for base tier (~15s).

    Evaluates 31 bodies x 6 calculation flags x 200 random dates and
    reports the maximum arcsecond difference between LEB2 and LEB1.
    """
    _python(["scripts/test_leb2_precision.py", "base"])


@leb2_format_group.command(
    "precision-medium",
    short_help="Measure LEB2 vs LEB1 precision for medium tier (~15s).",
)
def leb2_precision_medium() -> None:
    """Measure LEB2 vs LEB1 precision for medium tier (~15s)."""
    _python(["scripts/test_leb2_precision.py", "medium"])


@leb2_format_group.command(
    "precision-extended",
    short_help="Measure LEB2 vs LEB1 precision for extended tier (~15s).",
)
def leb2_precision_extended() -> None:
    """Measure LEB2 vs LEB1 precision for extended tier (~15s)."""
    _python(["scripts/test_leb2_precision.py", "extended"])


@leb2_format_group.command(
    "precision-all",
    short_help="Measure LEB2 vs LEB1 precision for ALL tiers (~45s).",
)
def leb2_precision_all() -> None:
    """Measure LEB2 vs LEB1 precision for ALL tiers (base + medium + extended, ~45s)."""
    _python(["scripts/test_leb2_precision.py", "all"])


test_group.add_command(leb2_format_group)


# ===========================================================================
# leph test horizons — Horizons API precision tests
# ===========================================================================


@click.group(
    "horizons",
    short_help="NASA JPL Horizons API precision tests (requires internet).",
    help="Precision tests for the NASA JPL Horizons REST API backend.\n\n"
    "These tests query the live Horizons API and compare results against\n"
    "Skyfield and/or LEB to verify the Horizons backend produces accurate\n"
    "positions. Requires an internet connection.\n\n"
    "  leph test horizons precision-quick   # 50 dates, ~15s\n"
    "  leph test horizons precision         # 200 dates, ~45s\n"
    "  leph test horizons vs-leb            # Cross-validate against LEB2",
)
def horizons_group() -> None:
    """Horizons API tests."""


@horizons_group.command(
    "precision",
    short_help="Horizons vs Skyfield: 13 bodies, 200 dates, 6 flags (~45s).",
)
def horizons_precision() -> None:
    """Horizons vs Skyfield precision: 13 bodies, 200 random dates, 6 flags (~45s).

    Measures the maximum arcsecond difference between positions computed
    via the Horizons API and positions computed via Skyfield/DE440.
    Requires internet connection.
    """
    _python(["scripts/test_horizons_precision.py"])


@horizons_group.command(
    "precision-quick",
    short_help="Quick Horizons vs Skyfield check with 50 dates (~15s).",
)
def horizons_precision_quick() -> None:
    """Quick Horizons vs Skyfield precision check with 50 dates (~15s).

    Faster subset of 'precision' for quick validation. Requires internet.
    """
    _python(["scripts/test_horizons_precision.py", "--dates", "50"])


@horizons_group.command(
    "vs-leb",
    short_help="Cross-validate Horizons API against LEB2: 15 bodies, 100 dates.",
)
def vs_leb() -> None:
    """Cross-validate Horizons API against LEB2: 15 bodies, 100 dates.

    Verifies that the Horizons backend and the LEB2 precomputed backend
    produce consistent results. Requires internet and LEB2 files.
    """
    _python(["scripts/test_horizons_vs_leb.py"])


test_group.add_command(horizons_group)


# ===========================================================================
# leph test coverage — Coverage reports
# ===========================================================================


@click.group(
    "coverage",
    short_help="Run tests with code coverage reports (terminal + XML + HTML).",
    help="Run the test suite with code coverage measurement.\n\n"
    "Generates three report formats: terminal (with missing lines highlighted),\n"
    "XML (for CI integration), and HTML (for browsing in a browser).\n"
    "Reports are written to the project root.\n\n"
    "  leph test coverage run    # Standard run, no @slow\n"
    "  leph test coverage full   # Including @slow tests",
)
def coverage_group() -> None:
    """Coverage report generation."""


_COV_ARGS = [
    "--cov=libephemeris",
    "--cov-report=term-missing",
    "--cov-report=xml",
    "--cov-report=html",
]


@coverage_group.command(
    short_help="Run all tests (no @slow) with coverage reports.",
)
def run() -> None:
    """Run all tests (excluding @slow) with coverage. Generates term + XML + HTML reports."""
    _pytest([*_COV_ARGS, "-m", "not slow"])


@coverage_group.command(
    short_help="Run ALL tests (including @slow) with coverage reports.",
)
def full() -> None:
    """Run ALL tests (including @slow) with coverage. Very slow but most thorough."""
    _pytest(_COV_ARGS)


test_group.add_command(coverage_group)
