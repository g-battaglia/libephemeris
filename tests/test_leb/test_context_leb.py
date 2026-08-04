"""
Tests for EphemerisContext LEB integration.

Validates that EphemerisContext.calc_ut() / calc() correctly use the LEB
fast path when a .leb file is configured, and falls back to Skyfield
when needed.
"""

from __future__ import annotations


import pytest

import libephemeris as ephem
from libephemeris.constants import (
    MEAN_NODE,
    MOON,
    SUN,
    FLG_EQUATORIAL,
    FLG_SPEED,
    FLG_TOPOCTR,
)
from libephemeris.context import EphemerisContext


class TestContextLEBFastPath:
    """Verify that EphemerisContext uses the LEB fast path."""

    @pytest.fixture
    def ctx_with_leb(self, test_leb_file):
        """Create an EphemerisContext with .leb file configured."""
        ctx = EphemerisContext()
        ctx.set_leb_file(test_leb_file)
        yield ctx
        ctx.set_leb_file(None)

    @pytest.fixture
    def jd_mid(self, test_leb_file):
        """Return the midpoint JD of the test LEB file."""
        from libephemeris.leb_reader import LEBReader

        reader = LEBReader(test_leb_file)
        jd_start, jd_end = reader.jd_range
        reader.close()
        return (jd_start + jd_end) / 2.0

    @pytest.mark.integration
    def test_calc_ut_with_leb(self, ctx_with_leb, jd_mid):
        """calc_ut() should produce valid results when .leb is set."""
        result, retflag = ctx_with_leb.calc_ut(jd_mid, SUN, FLG_SPEED)

        assert len(result) == 6
        assert 0.0 <= result[0] < 360.0
        assert -90.0 <= result[1] <= 90.0
        # Sun speed should be ~1 deg/day
        assert abs(result[3]) > 0.5

    @pytest.mark.integration
    def test_calc_ut_matches_global(self, ctx_with_leb, jd_mid, test_leb_file):
        """Context calc_ut() with LEB should match global calc_ut() with LEB."""
        # Set up global LEB too
        old_leb = None
        old_calc_mode = None
        try:
            from libephemeris import state

            old_leb = state._LEB_FILE
            old_calc_mode = state._CALC_MODE

            # Ensure calc mode allows LEB usage.  A preceding test on this
            # xdist worker may have set _CALC_MODE to "skyfield" or popped
            # the LIBEPHEMERIS_MODE env var, causing get_leb_reader() to
            # return None even though set_leb_file() was called.
            state.set_calc_mode("auto")
            ephem.set_leb_file(test_leb_file)

            ctx_result, _ = ctx_with_leb.calc_ut(jd_mid, SUN, FLG_SPEED)
            global_result, _ = ephem.calc_ut(jd_mid, SUN, FLG_SPEED)

            # Should match exactly (same LEB reader, same code path)
            for i in range(6):
                assert abs(ctx_result[i] - global_result[i]) < 1e-10, (
                    f"Component {i}: ctx={ctx_result[i]}, global={global_result[i]}"
                )
        finally:
            state.set_calc_mode(old_calc_mode)
            ephem.set_leb_file(old_leb)

    @pytest.mark.integration
    def test_calc_tt_with_leb(self, ctx_with_leb, jd_mid):
        """calc() (TT) should produce valid results when .leb is set."""
        result, retflag = ctx_with_leb.calc(jd_mid, SUN, FLG_SPEED)

        assert len(result) == 6
        assert 0.0 <= result[0] < 360.0

    @pytest.mark.integration
    def test_calc_ut_moon(self, ctx_with_leb, jd_mid):
        """Moon position via context + LEB should be valid."""
        result, _ = ctx_with_leb.calc_ut(jd_mid, MOON, FLG_SPEED)

        assert 0.0 <= result[0] < 360.0
        assert -10.0 <= result[1] <= 10.0  # Moon latitude within ~5 deg

    @pytest.mark.integration
    def test_calc_ut_mean_node(self, ctx_with_leb, jd_mid):
        """Mean node via context + LEB (ecliptic direct pipeline)."""
        result, _ = ctx_with_leb.calc_ut(jd_mid, MEAN_NODE, FLG_SPEED)

        assert 0.0 <= result[0] < 360.0
        # Mean node retrograde speed: ~-0.053 deg/day
        assert result[3] < 0.0, f"Mean node speed should be negative: {result[3]}"

    @pytest.mark.integration
    def test_calc_ut_equatorial(self, ctx_with_leb, jd_mid):
        """Equatorial coordinates via context + LEB."""
        result, _ = ctx_with_leb.calc_ut(jd_mid, SUN, FLG_EQUATORIAL | FLG_SPEED)

        assert 0.0 <= result[0] < 360.0
        assert -90.0 <= result[1] <= 90.0


class TestContextLEBFallback:
    """Verify that context falls back correctly when LEB can't serve a request."""

    @pytest.fixture
    def ctx_with_leb(self, test_leb_file):
        ctx = EphemerisContext()
        ctx.set_leb_file(test_leb_file)
        yield ctx
        ctx.set_leb_file(None)

    @pytest.fixture
    def jd_mid(self, test_leb_file):
        from libephemeris.leb_reader import LEBReader

        reader = LEBReader(test_leb_file)
        jd_start, jd_end = reader.jd_range
        reader.close()
        return (jd_start + jd_end) / 2.0

    @pytest.mark.integration
    def test_fallback_for_topoctr(self, ctx_with_leb, jd_mid):
        """FLG_TOPOCTR should fall back to Skyfield (not crash)."""
        ctx_with_leb.set_topo(12.5, 41.9, 0)
        result, _ = ctx_with_leb.calc_ut(jd_mid, SUN, FLG_TOPOCTR | FLG_SPEED)

        # Should get valid results from Skyfield fallback
        assert 0.0 <= result[0] < 360.0

    @pytest.mark.integration
    def test_fallback_for_unknown_body(self, ctx_with_leb, jd_mid):
        """Body not in .leb should fall back to Skyfield."""
        # MARS (4) is in the test LEB, but try a body that isn't
        # Jupiter (5) is not in the test fixture
        result, _ = ctx_with_leb.calc_ut(jd_mid, 5, FLG_SPEED)

        # Should get valid results from Skyfield fallback
        assert 0.0 <= result[0] < 360.0


class TestContextLEBGracefulError:
    """Test graceful error handling for invalid .leb files."""

    @pytest.mark.integration
    def test_invalid_leb_path_falls_back(self):
        """An absent .leb falls back in ``auto`` mode and raises when sealed.

        Sealed ``leb`` mode treats a missing or unreadable artifact as a
        provisioning error rather than silently substituting another source,
        so the mode has to be pinned explicitly here instead of inherited
        from LIBEPHEMERIS_MODE.
        """
        import pytest

        from libephemeris.state import get_calc_mode, set_calc_mode

        saved_mode = get_calc_mode()
        ctx = EphemerisContext()
        ctx.set_leb_file("/nonexistent/path/fake.leb")
        jd = 2460000.0  # ~2023
        try:
            set_calc_mode("auto")
            # get_leb_reader() should return None (not raise)
            assert ctx.get_leb_reader() is None
            # calc_ut should still work via the Skyfield fallback
            result, _ = ctx.calc_ut(jd, SUN, FLG_SPEED)
            assert 0.0 <= result[0] < 360.0

            set_calc_mode("leb")
            with pytest.raises((RuntimeError, OSError)):
                ctx.get_leb_reader()
        finally:
            set_calc_mode(saved_mode)

    @pytest.mark.integration
    def test_global_invalid_leb_path_falls_back(self):
        """Global set_leb_file with invalid path should fall back in auto mode."""
        import os

        from libephemeris import state

        old_leb = state._LEB_FILE
        old_mode_env = os.environ.pop("LIBEPHEMERIS_MODE", None)
        old_leb_env = os.environ.pop("LIBEPHEMERIS_LEB", None)
        try:
            state.set_calc_mode("auto")
            state._LEB_READER = None
            ephem.set_leb_file("/nonexistent/path/fake.leb")
            reader = state.get_leb_reader()
            assert reader is None

            # calc_ut should still work via Skyfield fallback
            jd = 2460000.0
            result, _ = ephem.calc_ut(jd, SUN, FLG_SPEED)
            assert 0.0 <= result[0] < 360.0
        finally:
            state.set_calc_mode(None)
            ephem.set_leb_file(old_leb)
            if old_mode_env is not None:
                os.environ["LIBEPHEMERIS_MODE"] = old_mode_env
            if old_leb_env is not None:
                os.environ["LIBEPHEMERIS_LEB"] = old_leb_env


class TestContextLEBGlobalFallthrough:
    """Test that context picks up the global LEB reader when no context-local one."""

    @pytest.mark.integration
    def test_context_uses_global_leb(self, test_leb_file):
        """Context without own .leb should use global reader if set."""
        from libephemeris import state
        from libephemeris.leb_reader import LEBReader

        old_leb = state._LEB_FILE
        try:
            # Set global LEB
            ephem.set_leb_file(test_leb_file)

            # Create context WITHOUT setting .leb
            ctx = EphemerisContext()
            assert ctx._leb_file is None

            reader_tmp = LEBReader(test_leb_file)
            jd_start, jd_end = reader_tmp.jd_range
            jd_mid = (jd_start + jd_end) / 2.0
            reader_tmp.close()

            # calc_ut should use the global LEB reader
            result, _ = ctx.calc_ut(jd_mid, SUN, FLG_SPEED)
            assert 0.0 <= result[0] < 360.0
            assert abs(result[3]) > 0.5  # Sun speed confirms it worked
        finally:
            ephem.set_leb_file(old_leb)


class TestCalcMode:
    """Tests for set_calc_mode() / get_calc_mode() and LIBEPHEMERIS_MODE."""

    def setup_method(self):
        """Reset calc mode before each test."""
        import os

        from libephemeris.state import set_calc_mode

        set_calc_mode(None)
        # Clear env var that may have been loaded from .env by _load_dotenv()
        os.environ.pop("LIBEPHEMERIS_MODE", None)

    def teardown_method(self):
        """Clean up after each test."""
        import os

        from libephemeris.state import set_calc_mode

        set_calc_mode(None)
        os.environ.pop("LIBEPHEMERIS_MODE", None)

    def test_default_mode_is_auto(self):
        """Default calc mode should be 'auto' when env var is unset."""
        from libephemeris.state import get_calc_mode

        assert get_calc_mode() == "auto"

    def test_set_mode_skyfield(self):
        """set_calc_mode('skyfield') should be respected."""
        from libephemeris.state import get_calc_mode, set_calc_mode

        set_calc_mode("skyfield")
        assert get_calc_mode() == "skyfield"

    def test_set_mode_leb(self):
        """set_calc_mode('leb') should be respected."""
        from libephemeris.state import get_calc_mode, set_calc_mode

        set_calc_mode("leb")
        assert get_calc_mode() == "leb"

    def test_set_mode_auto(self):
        """set_calc_mode('auto') should be respected."""
        from libephemeris.state import get_calc_mode, set_calc_mode

        set_calc_mode("auto")
        assert get_calc_mode() == "auto"

    def test_set_mode_none_resets(self):
        """set_calc_mode(None) should reset to env var / default."""
        from libephemeris.state import get_calc_mode, set_calc_mode

        set_calc_mode("skyfield")
        set_calc_mode(None)
        assert get_calc_mode() == "auto"

    def test_invalid_mode_raises(self):
        """set_calc_mode with invalid mode should raise ValueError."""
        from libephemeris.state import set_calc_mode

        with pytest.raises(ValueError, match="Invalid mode"):
            set_calc_mode("invalid")

    def test_mode_case_insensitive(self):
        """set_calc_mode should be case-insensitive."""
        from libephemeris.state import get_calc_mode, set_calc_mode

        set_calc_mode("SKYFIELD")
        assert get_calc_mode() == "skyfield"

    def test_env_var_mode(self):
        """LIBEPHEMERIS_MODE env var should set the mode."""
        import os

        from libephemeris.state import get_calc_mode

        os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
        assert get_calc_mode() == "skyfield"

    def test_env_var_case_insensitive(self):
        """LIBEPHEMERIS_MODE env var should be case-insensitive."""
        import os

        from libephemeris.state import get_calc_mode

        os.environ["LIBEPHEMERIS_MODE"] = "LEB"
        assert get_calc_mode() == "leb"

    def test_env_var_invalid_rejected(self):
        """An invalid LIBEPHEMERIS_MODE must fail loudly, not default to auto.

        Silently widening a typo'd mode back to ``auto`` would re-enable
        source fallback and the network in a deployment that asked to be
        sealed (see state._validated_calc_mode).
        """
        import os

        import pytest

        from libephemeris.state import get_calc_mode

        saved = os.environ.get("LIBEPHEMERIS_MODE")
        os.environ["LIBEPHEMERIS_MODE"] = "bogus"
        try:
            with pytest.raises(ValueError, match="Invalid calculation mode"):
                get_calc_mode()
        finally:
            if saved is None:
                os.environ.pop("LIBEPHEMERIS_MODE", None)
            else:
                os.environ["LIBEPHEMERIS_MODE"] = saved

    def test_programmatic_overrides_env_var(self):
        """Programmatic set_calc_mode should override env var."""
        import os

        from libephemeris.state import get_calc_mode, set_calc_mode

        os.environ["LIBEPHEMERIS_MODE"] = "leb"
        set_calc_mode("skyfield")
        assert get_calc_mode() == "skyfield"

    def test_skyfield_mode_returns_none_reader(self, test_leb_file):
        """In skyfield mode, get_leb_reader() should return None."""
        from libephemeris.state import get_leb_reader, set_calc_mode

        old_leb = ephem.state._LEB_FILE
        try:
            ephem.set_leb_file(test_leb_file)
            set_calc_mode("skyfield")
            assert get_leb_reader() is None
        finally:
            set_calc_mode(None)
            ephem.set_leb_file(old_leb)

    def test_leb_mode_without_file_raises(self):
        """In leb mode, get_leb_reader() should raise if no file configured."""
        import os
        from unittest.mock import patch

        from libephemeris.state import get_leb_reader, set_calc_mode

        old_leb = ephem.state._LEB_FILE
        old_leb_env = os.environ.pop("LIBEPHEMERIS_LEB", None)
        try:
            ephem.set_leb_file(None)
            ephem.state._LEB_READER = None
            set_calc_mode("leb")
            # Block both auto-discovery channels: a discovered reviewed tier
            # core would legitimately serve instead of raising.
            with (
                patch("libephemeris.state._discover_leb_file", return_value=None),
                patch(
                    "libephemeris.state._discover_reviewed_leb_tier_cores",
                    return_value={},
                ),
            ):
                with pytest.raises(RuntimeError, match="no .leb file configured"):
                    get_leb_reader()
        finally:
            set_calc_mode(None)
            ephem.set_leb_file(old_leb)
            if old_leb_env is not None:
                os.environ["LIBEPHEMERIS_LEB"] = old_leb_env

    def test_leb_mode_with_invalid_file_raises(self):
        """In leb mode, invalid .leb file should raise RuntimeError."""
        from libephemeris.state import get_leb_reader, set_calc_mode

        old_leb = ephem.state._LEB_FILE
        try:
            ephem.set_leb_file("/nonexistent/path/fake.leb")
            set_calc_mode("leb")
            with pytest.raises(RuntimeError, match="failed to open LEB file"):
                get_leb_reader()
        finally:
            set_calc_mode(None)
            ephem.set_leb_file(old_leb)

    @pytest.mark.integration
    def test_skyfield_mode_calc_works(self, test_leb_file):
        """In skyfield mode, calc_ut should use Skyfield (not LEB)."""
        from libephemeris.state import set_calc_mode

        old_leb = ephem.state._LEB_FILE
        try:
            ephem.set_leb_file(test_leb_file)
            set_calc_mode("skyfield")

            jd = 2451545.0  # J2000
            result, _ = ephem.calc_ut(jd, SUN, FLG_SPEED)
            assert 0.0 <= result[0] < 360.0
            assert abs(result[3]) > 0.5
        finally:
            set_calc_mode(None)
            ephem.set_leb_file(old_leb)

    @pytest.mark.integration
    def test_auto_mode_uses_leb(self, test_leb_file):
        """In auto mode with LEB configured, should use LEB fast path."""
        from libephemeris.state import get_leb_reader, set_calc_mode

        old_leb = ephem.state._LEB_FILE
        try:
            ephem.set_leb_file(test_leb_file)
            set_calc_mode("auto")

            reader = get_leb_reader()
            assert reader is not None
        finally:
            set_calc_mode(None)
            ephem.set_leb_file(old_leb)

    @pytest.mark.integration
    def test_close_preserves_calc_mode(self):
        """close() preserves the explicit calc mode (issue #30).

        The mode set via set_calc_mode() is process-level configuration;
        resetting it on close() silently re-enabled the auto fallback
        chain. Clearing requires an explicit set_calc_mode(None).
        """
        from libephemeris.state import get_calc_mode, set_calc_mode

        set_calc_mode("skyfield")
        assert get_calc_mode() == "skyfield"
        ephem.close()
        assert get_calc_mode() == "skyfield"
        set_calc_mode(None)

    @pytest.mark.integration
    def test_context_respects_global_skyfield_mode(self, test_leb_file):
        """EphemerisContext should respect global skyfield mode."""
        from libephemeris.state import set_calc_mode

        old_leb = ephem.state._LEB_FILE
        try:
            ephem.set_leb_file(test_leb_file)
            set_calc_mode("skyfield")

            # Context without own LEB should pick up global mode
            ctx = EphemerisContext()
            jd = 2451545.0
            result, _ = ctx.calc_ut(jd, SUN, FLG_SPEED)
            assert 0.0 <= result[0] < 360.0
        finally:
            set_calc_mode(None)
            ephem.set_leb_file(old_leb)
