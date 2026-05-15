"""
Tests for the close() function that releases ephemeris resources.

Tests the module-level close() function.
Note: EphemerisContext.close() tests are skipped due to pre-existing
context loading issues in the codebase.
"""

import libephemeris
from libephemeris import (
    close,
    calc_ut,
    SUN,
)
from libephemeris import state


class TestCloseFunction:
    """Tests for the module-level close() function."""

    def test_close_resets_planets(self):
        """close() should reset the _PLANETS global to None."""
        # Explicitly load planets (calc_ut may use LEB fast path)
        state.get_planets()
        assert state._PLANETS is not None

        # Call close
        close()

        # Verify _PLANETS is reset
        assert state._PLANETS is None

    def test_close_resets_timescale(self):
        """close() should reset the _TS global to None."""
        # Ensure timescale is loaded
        calc_ut(2451545.0, SUN, 0)
        assert state._TS is not None

        # Call close
        close()

        # Verify _TS is reset
        assert state._TS is None

    def test_close_resets_loader(self):
        """close() should reset the _LOADER global to None."""
        # Explicitly load planets to ensure loader is created
        state.get_planets()
        assert state._LOADER is not None

        # Call close
        close()

        # Verify _LOADER is reset
        assert state._LOADER is None

    def test_close_resets_ephemeris_path(self):
        """close() should reset the _EPHEMERIS_PATH to None."""
        libephemeris.set_ephe_path("/some/path")
        assert state._EPHEMERIS_PATH == "/some/path"

        close()

        assert state._EPHEMERIS_PATH is None

    def test_close_resets_ephemeris_file(self):
        """close() should reset the _EPHEMERIS_FILE to default."""
        libephemeris.set_ephemeris_file("de441.bsp")
        assert state._EPHEMERIS_FILE == "de441.bsp"

        close()

        assert state._EPHEMERIS_FILE == "de440.bsp"

    def test_close_resets_topo(self):
        """close() should reset the _TOPO to None."""
        libephemeris.set_topo(12.5, 41.9, 0)
        assert state._TOPO is not None

        close()

        assert state._TOPO is None

    def test_close_resets_sidereal_mode(self):
        """close() should reset sidereal mode settings."""
        libephemeris.set_sid_mode(3, 2451545.0, 23.5)
        assert state._SIDEREAL_MODE == 3
        assert state._SIDEREAL_AYAN_T0 == 23.5

        close()

        assert state._SIDEREAL_MODE is None
        assert state._SIDEREAL_AYAN_T0 == 0.0
        assert state._SIDEREAL_T0 == 0.0

    def test_close_clears_angles_cache(self):
        """close() should clear the angles cache."""
        state.set_angles_cache({"Sun": 120.0, "Moon": 240.0})
        assert state.get_angles_cache() == {"Sun": 120.0, "Moon": 240.0}

        close()

        assert state.get_angles_cache() == {}

    def test_close_resets_tidal_acceleration(self):
        """close() should reset tidal acceleration to None."""
        libephemeris.set_tid_acc(-25.85)
        assert state._TIDAL_ACCELERATION == -25.85

        close()

        assert state._TIDAL_ACCELERATION is None

    def test_close_resets_delta_t_userdef(self):
        """close() should reset user-defined delta T to None."""
        libephemeris.set_delta_t_userdef(0.001)
        assert state._DELTA_T_USERDEF == 0.001

        close()

        assert state._DELTA_T_USERDEF is None

    def test_close_resets_lapse_rate(self):
        """close() should reset lapse rate to None (default)."""
        libephemeris.set_lapse_rate(0.008)
        assert state._LAPSE_RATE == 0.008

        close()

        assert state._LAPSE_RATE is None

    def test_calculations_work_after_close(self):
        """Calculations should work normally after calling close()."""
        # First calculation loads ephemeris
        pos1, _ = calc_ut(2451545.0, SUN, 0)

        # Close releases resources
        close()

        # Second calculation should reload and produce same results
        pos2, _ = calc_ut(2451545.0, SUN, 0)

        # Compare longitude to high precision (LEB vs Skyfield may differ
        # at sub-arcsecond level after close() clears cached state)
        assert abs(pos1[0] - pos2[0]) < 1e-6, (
            f"Longitude differs after close(): {pos1[0]} vs {pos2[0]}"
        )
        assert abs(pos1[1] - pos2[1]) < 1e-6, (
            f"Latitude differs after close(): {pos1[1]} vs {pos2[1]}"
        )
        assert abs(pos1[2] - pos2[2]) < 1e-6, (
            f"Distance differs after close(): {pos1[2]} vs {pos2[2]}"
        )

