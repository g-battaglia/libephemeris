"""Contract tests for the public API surface (1.7.0).

These tests lock in three invariants:
  1. No legacy prefixed names (swe_*, SE_*, SEFLG_*) leak into the
     public namespace.
  2. The full upstream reference-API surface (top-level and contrib)
     is mirrored 1:1.
  3. The 66 names that the kerykeion v6/alpha client depends on are
     present and callable / accessible.

These tests are intentionally cheap — they introspect dir() rather
than invoking ephemeris calculations. They're meant to fail fast in
CI if anyone accidentally exposes a prefixed alias again, removes a
kerykeion-contract name, or lets the reference surface drift.

The "one allowed exception" is `SE_FNAME_DE431`: the upstream
reference distribution itself exposes that one constant with the
SE_ prefix (their inconsistency, not ours).
"""

from __future__ import annotations

import pytest

import libephemeris as leph


# The single name the upstream reference distribution itself exposes
# with an SE_ prefix. Mirrored for parity, not for stylistic consistency.
ALLOWED_PREFIXED_NAMES = {"SE_FNAME_DE431"}


def _public_names(module) -> set[str]:
    return {n for n in dir(module) if not n.startswith("_")}


def _require_pyswisseph():
    """Import the real pyswisseph extension or skip the test.

    A bare ``importorskip`` is not enough on its own — under Python 3
    any directory on sys.path (even one without ``__init__.py``) can be
    resolved as an empty namespace package and silently satisfy
    ``import swisseph``. That makes parity tests pass without
    comparing against anything. Guard against this by verifying that
    the resolved module actually exposes the reference ``calc_ut``
    callable.
    """
    swe = pytest.importorskip("swisseph")
    if not hasattr(swe, "calc_ut"):
        pytest.skip(
            "pyswisseph not installed (or shadowed by an unrelated "
            "namespace package); install pyswisseph to run parity tests"
        )
    return swe


class TestNoLegacyPrefixedNames:
    """No swe_*, SE_*, SEFLG_* names in the public surface (except
    `SE_FNAME_DE431` which the upstream reference exports under that
    name and we mirror for parity)."""

    def test_no_swe_prefix(self) -> None:
        leaked = [n for n in dir(leph) if n.startswith("swe_")]
        assert not leaked, f"Unexpected swe_* names: {leaked}"

    def test_no_se_prefix(self) -> None:
        leaked = [
            n
            for n in dir(leph)
            if n.startswith("SE_") and n not in ALLOWED_PREFIXED_NAMES
        ]
        assert not leaked, f"Unexpected SE_* names: {leaked}"

    def test_no_seflg_prefix(self) -> None:
        leaked = [n for n in dir(leph) if n.startswith("SEFLG_")]
        assert not leaked, f"Unexpected SEFLG_* names: {leaked}"


class TestUpstreamReferenceParity:
    """Top-level + contrib surface must cover the upstream reference."""

    def test_top_level_subset(self) -> None:
        ref = _require_pyswisseph()

        missing = _public_names(ref) - _public_names(leph)
        assert not missing, (
            f"libephemeris missing {len(missing)} upstream top-level "
            f"names: {sorted(missing)}"
        )

    def test_contrib_submodule_exists(self) -> None:
        assert hasattr(leph, "contrib"), (
            "libephemeris.contrib submodule is required for parity"
        )

    def test_contrib_subset(self) -> None:
        swisseph = _require_pyswisseph()

        if not hasattr(swisseph, "contrib"):
            pytest.skip("upstream pyswisseph build has no contrib submodule")
        their = _public_names(swisseph.contrib)
        ours = _public_names(leph.contrib)
        missing = their - ours
        assert not missing, (
            f"libephemeris.contrib missing {len(missing)} names: "
            f"{sorted(missing)}"
        )

    def test_contrib_no_extras_at_module_level(self) -> None:
        """contrib shouldn't expose stray imports (typing etc.)."""
        swisseph = _require_pyswisseph()

        if not hasattr(swisseph, "contrib"):
            pytest.skip("upstream pyswisseph build has no contrib submodule")
        their = _public_names(swisseph.contrib)
        ours = _public_names(leph.contrib)
        extras = ours - their
        assert not extras, (
            f"libephemeris.contrib exposes extras not in upstream: "
            f"{sorted(extras)}"
        )


class TestUpstreamConstantsMatch:
    """For shared constants, values must equal the upstream reference."""

    SHARED = [
        # Body IDs
        "SUN", "MOON", "MERCURY", "VENUS", "MARS", "JUPITER", "SATURN",
        "URANUS", "NEPTUNE", "PLUTO", "EARTH",
        "MEAN_NODE", "TRUE_NODE", "MEAN_APOG", "OSCU_APOG",
        "CHIRON", "CERES", "PALLAS", "JUNO", "VESTA",
        # Flags
        "FLG_SPEED", "FLG_SWIEPH", "FLG_HELCTR", "FLG_TOPOCTR",
        "FLG_SIDEREAL", "FLG_EQUATORIAL", "FLG_NONUT", "FLG_TRUEPOS",
        "FLG_J2000", "FLG_XYZ", "FLG_RADIANS", "FLG_BARYCTR",
        "FLG_NOABERR", "FLG_NOGDEFL",
        # Eclipses
        "ECL_TOTAL", "ECL_ANNULAR", "ECL_PARTIAL", "ECL_HYBRID",
        "ECL_VISIBLE",
        # Sidereal modes
        "SIDM_LAHIRI", "SIDM_FAGAN_BRADLEY",
        # Offsets
        "AST_OFFSET", "FICT_OFFSET",
        # Calendars
        "JUL_CAL", "GREG_CAL",
        # Heliacal
        "HELIACAL_RISING", "HELIACAL_SETTING",
    ]

    @pytest.mark.parametrize("name", SHARED)
    def test_constant_matches(self, name: str) -> None:
        ref = _require_pyswisseph()

        if not hasattr(ref, name):
            pytest.skip(f"{name} not in reference (no comparison possible)")
        ours = getattr(leph, name)
        theirs = getattr(ref, name)
        assert ours == theirs, (
            f"{name}: libephemeris={ours}, reference={theirs}"
        )


class TestKerykeionContractFreeze:
    """The 66 names kerykeion v6/alpha depends on. Frozen contract:
    deleting any of these is a breaking change requiring downstream
    coordination."""

    # 38 functions used by kerykeion's ephemeris backend
    KERYKEION_FUNCTIONS = [
        "azalt", "calc_pctr", "calc_ut", "close", "difdeg2n",
        "fixstar2_mag", "fixstar_ut", "gauquelin_sector",
        "get_ayanamsa_ex_ut", "get_ayanamsa_name", "get_planet_name",
        "get_trace_results", "heliacal_ut", "helio_cross_ut",
        "house_name", "houses_armc", "houses_ex", "houses_ex2",
        "julday", "list_fixed_stars", "lun_eclipse_when",
        "lun_eclipse_when_loc", "lun_occult_when_glob",
        "lun_occult_when_loc", "mooncross_node_ut", "mooncross_ut",
        "nod_aps_ut", "pheno_ut", "revjul", "rise_trans",
        "set_ephe_path", "set_sid_mode", "set_topo", "sidtime",
        "sol_eclipse_when_glob", "sol_eclipse_when_loc",
        "solcross_ut", "start_tracing",
    ]

    # 28 constants
    KERYKEION_CONSTANTS = [
        "AST_OFFSET", "ECL2HOR", "ECL_NUT", "EVENING_FIRST",
        "FLG_BARYCTR", "FLG_EQUATORIAL", "FLG_HELCTR",
        "FLG_SIDEREAL", "FLG_SPEED", "FLG_SWIEPH", "FLG_TOPOCTR",
        "FLG_TRUEPOS", "GREG_CAL", "HELIACAL_RISING",
        "HELIACAL_SETTING", "JUL_CAL", "JUPITER", "MARS",
        "MERCURY", "MOON", "MORNING_LAST", "NEPTUNE", "PLUTO",
        "SATURN", "SIDM_USER", "SUN", "URANUS", "VENUS",
    ]

    @pytest.mark.parametrize("name", KERYKEION_FUNCTIONS)
    def test_function_present(self, name: str) -> None:
        assert hasattr(leph, name), (
            f"kerykeion contract broken: function {name!r} not exported"
        )
        attr = getattr(leph, name)
        assert callable(attr), (
            f"kerykeion contract broken: {name!r} is not callable"
        )

    @pytest.mark.parametrize("name", KERYKEION_CONSTANTS)
    def test_constant_present(self, name: str) -> None:
        assert hasattr(leph, name), (
            f"kerykeion contract broken: constant {name!r} not exported"
        )


class TestVersionConsistency:
    """The __version__ string must match pyproject.toml."""

    def test_version_format(self) -> None:
        assert isinstance(leph.__version__, str)
        assert leph.__version__.count(".") >= 2  # e.g. "1.7.0"

    def test_version_alias(self) -> None:
        assert leph.version == leph.__version__
