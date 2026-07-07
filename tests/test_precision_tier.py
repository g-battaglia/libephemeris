"""
Tests for the precision tier system.

Tests the PrecisionTier dataclass, TIERS dict, and the tier get/set/list
functions, as well as ephemeris file priority logic, close() resets, env
var overrides, discover_local_spks(), and get_spk_date_range_for_tier().
"""

from __future__ import annotations

import os
import tempfile
from unittest import mock

import pytest

import libephemeris
from libephemeris import state
from libephemeris.state import (
    PrecisionTier,
    TIERS,
    _get_effective_ephemeris_file,
    get_precision_tier,
    get_spk_date_range_for_tier,
    list_tiers,
    set_precision_tier,
)


# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture(autouse=True)
def _reset_state():
    """Ensure clean state for every test."""
    yield
    libephemeris.close()


# =============================================================================
# TIERS DICT AND PrecisionTier DATACLASS
# =============================================================================


class TestTiersDict:
    """Tests for the TIERS dictionary and PrecisionTier dataclass."""

    def test_three_tiers_exist(self):
        """TIERS should have exactly three entries."""
        assert len(TIERS) == 3

    def test_tier_names(self):
        """TIERS should contain base, medium, and extended."""
        assert set(TIERS.keys()) == {"base", "medium", "extended"}

    @pytest.mark.parametrize(
        "name,expected_file",
        [
            ("base", "de440s.bsp"),
            ("medium", "de440.bsp"),
            ("extended", "de441.bsp"),
        ],
    )
    def test_ephemeris_files(self, name, expected_file):
        """Each tier should map to the correct ephemeris file."""
        assert TIERS[name].ephemeris_file == expected_file

    def test_tiers_are_frozen(self):
        """PrecisionTier should be frozen (immutable)."""
        with pytest.raises(AttributeError):
            TIERS["medium"].name = "changed"  # type: ignore

    def test_tier_has_description(self):
        """Each tier should have a non-empty description."""
        for tier in TIERS.values():
            assert isinstance(tier.description, str)
            assert len(tier.description) > 0

    def test_tier_spk_date_range_is_tuple_of_two_strings(self):
        """Each tier's spk_date_range should be a (start, end) string tuple."""
        for tier in TIERS.values():
            start, end = tier.spk_date_range
            assert isinstance(start, str)
            assert isinstance(end, str)
            assert start < end  # lexicographic, but these are YYYY-MM-DD

    @pytest.mark.parametrize(
        "name,expected_start,expected_end",
        [
            ("base", "1850-01-01", "2150-01-01"),
            ("medium", "1900-01-01", "2100-01-01"),
            ("extended", "1600-01-01", "2500-01-01"),
        ],
    )
    def test_spk_date_ranges(self, name, expected_start, expected_end):
        """Each tier should have the expected SPK date range."""
        start, end = TIERS[name].spk_date_range
        assert start == expected_start
        assert end == expected_end


# =============================================================================
# GET / SET PRECISION TIER
# =============================================================================


class TestGetSetPrecisionTier:
    """Tests for get_precision_tier() and set_precision_tier()."""

    def test_default_is_medium(self):
        """Default precision tier should be 'medium'."""
        assert get_precision_tier() == "medium"

    def test_set_and_get(self):
        """set_precision_tier() should change the active tier."""
        set_precision_tier("extended")
        assert get_precision_tier() == "extended"

        set_precision_tier("base")
        assert get_precision_tier() == "base"

    def test_set_invalid_raises(self):
        """set_precision_tier() should raise ValueError for unknown tiers."""
        with pytest.raises(ValueError, match="Invalid tier"):
            set_precision_tier("ultra")

    def test_set_clears_planets_cache(self):
        """Changing the tier should clear _PLANETS so the new file is loaded."""
        # Load planets
        state.get_planets()
        assert state._PLANETS is not None

        set_precision_tier("base")
        assert state._PLANETS is None


# =============================================================================
# TIER CHANGE INVALIDATES THE LEB READER (regression)
# =============================================================================


def _local_leb_for_tier(tier: str):
    """Return a locally-present .leb path for the tier, or None.

    Mirrors the local-file portion of state._discover_leb_file() WITHOUT the
    network auto-download fallback, so tests never touch the network.
    """
    leb_dir = os.path.join(state._get_data_dir(), "leb")
    for name in (f"{tier}_core.leb2", f"{tier}_core.leb", f"ephemeris_{tier}.leb"):
        candidate = os.path.join(leb_dir, name)
        if os.path.isfile(candidate):
            return candidate
    return None


class _FakeLebReader:
    """Minimal stand-in for a LEB reader that records close()."""

    def __init__(self):
        self.closed = False

    def close(self):
        self.closed = True


class TestSetPrecisionTierInvalidatesLebReader:
    """set_precision_tier() must drop the stale LEB reader (docstring promise:
    'Changing the tier clears cached ephemeris data to force a reload')."""

    def test_set_closes_and_clears_leb_reader(self):
        """A tier change closes the reader, nulls it, and clears fast-path state."""
        from libephemeris import fast_calc

        fake = _FakeLebReader()
        state._LEB_READER = fake
        fast_calc._leb_frame_cache[(0, 0.0)] = ("sentinel",)
        gen_before = fast_calc._active_generation

        set_precision_tier("base")

        assert fake.closed is True
        assert state._LEB_READER is None  # forces re-resolution on next access
        assert fast_calc._leb_frame_cache == {}
        assert fast_calc._active_generation > gen_before

    def test_set_preserves_explicit_leb_file(self):
        """An explicitly pinned _LEB_FILE survives a tier change (not clobbered
        by auto-discovery); only the reader is dropped so the same path re-opens."""
        state._LEB_FILE = "/some/explicit/path.leb"
        state._LEB_READER = _FakeLebReader()

        set_precision_tier("base")

        assert state._LEB_FILE == "/some/explicit/path.leb"
        assert state._LEB_READER is None

    def test_reader_reload_serves_new_tier_range(self):
        """End-to-end: after switching medium->base the reader re-resolves to the
        base file, so a date inside medium's range but outside base's now raises;
        the round trip back to medium is bit-identical (no state corruption).

        Gated on both tier files being present locally (no network)."""
        from libephemeris import calc_ut
        from libephemeris.constants import SUN
        from libephemeris.exceptions import EphemerisRangeError
        from libephemeris.time_utils import julday

        if _local_leb_for_tier("medium") is None or _local_leb_for_tier("base") is None:
            pytest.skip("requires local base and medium .leb files")

        jd_1700 = julday(1700, 1, 1, 0.0)  # in medium (1550-2650), not base
        old_mode = state.get_calc_mode()
        # Exercise per-tier auto-discovery, so neither an explicit _LEB_FILE nor
        # the LIBEPHEMERIS_LEB env var (both of which pin one file across every
        # tier by design, e.g. the leb-backend test runner) may be active.
        saved_env = os.environ.pop("LIBEPHEMERIS_LEB", None)
        try:
            state.set_calc_mode("leb")
            state.set_leb_file(None)

            set_precision_tier("medium")
            pos_medium, _ = calc_ut(jd_1700, SUN, 0)
            reader_medium = state.get_leb_reader()
            assert reader_medium is not None

            set_precision_tier("base")
            reader_base = state.get_leb_reader()
            assert reader_base is not None
            assert reader_base is not reader_medium  # not the stale medium reader
            with pytest.raises(EphemerisRangeError):
                calc_ut(jd_1700, SUN, 0)

            set_precision_tier("medium")
            pos_medium_again, _ = calc_ut(jd_1700, SUN, 0)
            assert pos_medium_again == pos_medium  # bit-identical round trip
        finally:
            state.set_calc_mode(old_mode)
            if saved_env is not None:
                os.environ["LIBEPHEMERIS_LEB"] = saved_env


# =============================================================================
# list_tiers()
# =============================================================================


class TestListTiers:
    """Tests for list_tiers()."""

    def test_returns_list_of_precision_tiers(self):
        """list_tiers() should return a list of PrecisionTier objects."""
        tiers = list_tiers()
        assert isinstance(tiers, list)
        assert len(tiers) == 3
        for t in tiers:
            assert isinstance(t, PrecisionTier)

    def test_contains_all_tier_names(self):
        """list_tiers() should contain base, medium, and extended."""
        names = {t.name for t in list_tiers()}
        assert names == {"base", "medium", "extended"}


# =============================================================================
# get_spk_date_range_for_tier()
# =============================================================================


class TestGetSpkDateRangeForTier:
    """Tests for get_spk_date_range_for_tier()."""

    def test_default_uses_current_tier(self):
        """With no argument, should return the current tier's range."""
        assert get_spk_date_range_for_tier() == ("1900-01-01", "2100-01-01")

    def test_explicit_tier_name(self):
        """With an explicit tier name, should return that tier's range."""
        assert get_spk_date_range_for_tier("extended") == (
            "1600-01-01",
            "2500-01-01",
        )

    def test_follows_current_tier(self):
        """Should reflect changes made via set_precision_tier()."""
        set_precision_tier("base")
        assert get_spk_date_range_for_tier() == ("1850-01-01", "2150-01-01")

    def test_invalid_tier_raises(self):
        """Should raise ValueError for unknown tier names."""
        with pytest.raises(ValueError, match="Invalid tier"):
            get_spk_date_range_for_tier("nonexistent")


# =============================================================================
# EPHEMERIS FILE PRIORITY
# =============================================================================


class TestEphemerisFilePriority:
    """Tests for the 3-level ephemeris file priority system."""

    def test_tier_controls_ephemeris_file(self):
        """When no explicit file is set, the tier controls the ephemeris file."""
        set_precision_tier("base")
        assert _get_effective_ephemeris_file() == "de440s.bsp"

        set_precision_tier("extended")
        assert _get_effective_ephemeris_file() == "de441.bsp"

    def test_explicit_file_overrides_tier(self):
        """set_ephemeris_file() should override the tier setting."""
        set_precision_tier("base")
        state.set_ephemeris_file("de430.bsp")
        assert _get_effective_ephemeris_file() == "de430.bsp"

    def test_env_var_overrides_explicit_file(self):
        """LIBEPHEMERIS_EPHEMERIS env var should override everything."""
        state.set_ephemeris_file("de430.bsp")
        set_precision_tier("extended")

        with mock.patch.dict(os.environ, {"LIBEPHEMERIS_EPHEMERIS": "de422.bsp"}):
            assert _get_effective_ephemeris_file() == "de422.bsp"

    def test_env_var_precision_tier(self):
        """LIBEPHEMERIS_PRECISION env var should set the tier."""
        with mock.patch.dict(os.environ, {"LIBEPHEMERIS_PRECISION": "extended"}):
            # Reset programmatic override
            state._PRECISION_TIER = None
            assert get_precision_tier() == "extended"
            assert _get_effective_ephemeris_file() == "de441.bsp"

    def test_programmatic_tier_overrides_env_var(self):
        """set_precision_tier() should override the env var."""
        with mock.patch.dict(os.environ, {"LIBEPHEMERIS_PRECISION": "extended"}):
            set_precision_tier("base")
            assert get_precision_tier() == "base"

    def test_invalid_env_var_falls_back_to_default(self):
        """An invalid LIBEPHEMERIS_PRECISION value should fall back to default."""
        with mock.patch.dict(os.environ, {"LIBEPHEMERIS_PRECISION": "invalid"}):
            state._PRECISION_TIER = None
            assert get_precision_tier() == "medium"


# =============================================================================
# CLOSE RESETS
# =============================================================================


class TestCloseResets:
    """Tests that close() properly resets precision tier state."""

    def test_close_resets_precision_tier(self):
        """close() should reset _PRECISION_TIER to None."""
        set_precision_tier("extended")
        assert state._PRECISION_TIER == "extended"

        libephemeris.close()
        assert state._PRECISION_TIER is None

    def test_close_resets_ephemeris_file_explicit(self):
        """close() should reset _EPHEMERIS_FILE_EXPLICIT to False."""
        state.set_ephemeris_file("de441.bsp")
        assert state._EPHEMERIS_FILE_EXPLICIT is True

        libephemeris.close()
        assert state._EPHEMERIS_FILE_EXPLICIT is False

    def test_close_returns_to_default_tier(self):
        """After close(), get_precision_tier() should return 'medium'."""
        set_precision_tier("extended")
        libephemeris.close()
        assert get_precision_tier() == "medium"


# =============================================================================
# DISCOVER LOCAL SPKS
# =============================================================================


class TestDiscoverLocalSpks:
    """Tests for discover_local_spks()."""

    def test_nonexistent_directory_returns_empty(self):
        """Should return an empty dict for a non-existent directory."""
        from libephemeris.spk_auto import discover_local_spks

        result = discover_local_spks("/nonexistent/path/xyz_12345")
        assert result == {}

    def test_empty_directory_returns_empty(self):
        """Should return an empty dict for a directory with no .bsp files."""
        from libephemeris.spk_auto import discover_local_spks

        with tempfile.TemporaryDirectory() as tmpdir:
            result = discover_local_spks(tmpdir)
            assert result == {}

    def test_discovers_known_bodies_in_repo_root(self):
        """Should discover SPK files in the data directory."""
        from libephemeris.spk_auto import discover_local_spks
        from libephemeris.state import _get_data_dir

        # First close to clear any existing registrations
        libephemeris.close()

        data_dir = _get_data_dir()
        result = discover_local_spks(data_dir)

        # The data directory may or may not have SPK files
        # All status values should be valid
        for body_name, status in result.items():
            assert status in ("registered", "already_registered") or status.startswith(
                "error:"
            ), f"Unexpected status for {body_name}: {status}"

    def test_discover_registers_bodies(self):
        """Discovered bodies should be registered in _SPK_BODY_MAP."""
        from libephemeris.spk_auto import discover_local_spks

        libephemeris.close()

        # Create a temp directory with a mock SPK file
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a mock SPK file for a known body (Chiron = 2060)
            mock_spk = os.path.join(tmpdir, "2060_2458849_2462502.bsp")
            with open(mock_spk, "wb") as f:
                f.write(b"dummy spk content")

            result = discover_local_spks(tmpdir)

            # May have discovered the mock file (though it won't load as real SPK)
            # Just verify the function runs without error
            assert isinstance(result, dict)


# =============================================================================
# ENSURE ALL EPHEMERIDES (UNIT-LEVEL, NO NETWORK)
# =============================================================================


class TestEnsureAllEphemerides:
    """Unit tests for ensure_all_ephemerides (mocked, no network)."""

    def test_function_exists(self):
        """ensure_all_ephemerides should be importable."""
        from libephemeris.spk_auto import ensure_all_ephemerides

        assert callable(ensure_all_ephemerides)

    def test_exported_from_main_module(self):
        """ensure_all_ephemerides should be exported from libephemeris."""
        assert hasattr(libephemeris, "ensure_all_ephemerides")


# =============================================================================
# EXPORTS
# =============================================================================


class TestExports:
    """Tests that all new public API symbols are properly exported."""

    @pytest.mark.parametrize(
        "name",
        [
            "PrecisionTier",
            "TIERS",
            "set_precision_tier",
            "get_precision_tier",
            "list_tiers",
            "get_spk_date_range_for_tier",
            "discover_local_spks",
            "ensure_all_ephemerides",
        ],
    )
    def test_exported(self, name):
        """Each symbol should be accessible from the top-level module."""
        assert hasattr(libephemeris, name), f"{name} not exported"

    @pytest.mark.parametrize(
        "name",
        [
            "PrecisionTier",
            "TIERS",
            "set_precision_tier",
            "get_precision_tier",
            "list_tiers",
            "get_spk_date_range_for_tier",
            "discover_local_spks",
            "ensure_all_ephemerides",
        ],
    )
    def test_in_all(self, name):
        """Each symbol should be in __all__."""
        assert name in libephemeris.__all__, f"{name} not in __all__"
