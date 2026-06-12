"""
Tests for state management and calc_mode switching.

Verifies that set_calc_mode, get_calc_mode, close,
and related state functions behave correctly.
"""

from __future__ import annotations

import pytest

import libephemeris as swe
from libephemeris.state import set_calc_mode, get_calc_mode
from libephemeris.constants import (
    SUN,
    MARS,
    FLG_SPEED,
    FLG_SIDEREAL,
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
    SIDM_RAMAN,
)


class TestCalcModeBasic:
    """Basic calc_mode getter/setter tests."""

    @pytest.mark.unit
    def test_get_calc_mode_returns_string(self):
        """get_calc_mode returns a string."""
        mode = get_calc_mode()
        assert isinstance(mode, str)

    @pytest.mark.unit
    def test_set_calc_mode_skyfield(self):
        """Can set mode to skyfield."""
        set_calc_mode("skyfield")
        assert get_calc_mode() == "skyfield"

    @pytest.mark.unit
    def test_set_calc_mode_auto(self):
        """Can set mode to auto."""
        set_calc_mode("auto")
        assert get_calc_mode() == "auto"

    @pytest.mark.unit
    def test_set_calc_mode_none_resets(self):
        """Setting mode to None resets to default."""
        set_calc_mode("skyfield")
        set_calc_mode(None)
        mode = get_calc_mode()
        # After reset, resolves to env var or auto/leb
        assert mode in ("auto", "leb", "skyfield", "horizons")

    @pytest.mark.unit
    def test_set_calc_mode_invalid_raises(self):
        """Setting an invalid mode raises ValueError."""
        with pytest.raises(ValueError):
            set_calc_mode("invalid_mode")

    @pytest.mark.unit
    def test_set_calc_mode_case_insensitive(self):
        """Mode setting is case-insensitive."""
        set_calc_mode("SKYFIELD")
        assert get_calc_mode() == "skyfield"

    @pytest.mark.unit
    def test_set_calc_mode_with_whitespace(self):
        """Mode setting handles whitespace."""
        set_calc_mode("  skyfield  ")
        assert get_calc_mode() == "skyfield"


class TestCalcModeCalculations:
    """Test that calc_mode affects calculations."""

    @pytest.mark.unit
    def test_skyfield_mode_produces_results(self):
        """Skyfield mode produces valid results."""
        set_calc_mode("skyfield")
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, SUN, FLG_SPEED)
        assert 0 <= result[0] < 360

    @pytest.mark.unit
    def test_auto_mode_produces_results(self):
        """Auto mode produces valid results."""
        set_calc_mode("auto")
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, SUN, FLG_SPEED)
        assert 0 <= result[0] < 360

    @pytest.mark.unit
    def test_mode_switch_consistent_results(self):
        """Switching modes should produce reasonably consistent results."""
        jd = 2451545.0

        set_calc_mode("skyfield")
        r_sky, _ = swe.calc_ut(jd, MARS, FLG_SPEED)

        set_calc_mode("auto")
        r_auto, _ = swe.calc_ut(jd, MARS, FLG_SPEED)

        # Results should be very close regardless of mode
        lon_diff = abs(r_sky[0] - r_auto[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        assert lon_diff < 0.01, f"Skyfield vs Auto Mars: {lon_diff:.6f}° difference"


class TestSweClose:
    """Test close behavior."""

    @pytest.mark.unit
    def test_swe_close_does_not_crash(self):
        """close should not crash."""
        swe.close()

    @pytest.mark.unit
    def test_calc_works_after_close(self):
        """Calculations should work after close."""
        swe.close()
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, SUN, FLG_SPEED)
        assert 0 <= result[0] < 360

    @pytest.mark.unit
    def test_multiple_close_calls(self):
        """Multiple close calls should be safe."""
        for _ in range(5):
            swe.close()
        result, _ = swe.calc_ut(2451545.0, SUN, 0)
        assert 0 <= result[0] < 360


class TestSiderealModeState:
    """Test sidereal mode state management."""

    @pytest.mark.unit
    def test_set_sid_mode_lahiri(self):
        """Can set Lahiri sidereal mode."""
        swe.set_sid_mode(SIDM_LAHIRI)
        # Should not raise

    @pytest.mark.unit
    def test_set_sid_mode_fagan(self):
        """Can set Fagan-Bradley sidereal mode."""
        swe.set_sid_mode(SIDM_FAGAN_BRADLEY)
        # Should not raise

    @pytest.mark.unit
    def test_sid_mode_affects_results(self):
        """Different sidereal modes should produce different longitudes."""
        jd = 2451545.0

        swe.set_sid_mode(SIDM_LAHIRI)
        r1, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)

        swe.set_sid_mode(SIDM_FAGAN_BRADLEY)
        r2, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)

        diff = abs(r1[0] - r2[0])
        if diff > 180:
            diff = 360 - diff
        assert diff > 0.01, f"Lahiri vs Fagan same result: diff={diff:.4f}°"

    @pytest.mark.unit
    def test_three_sid_modes_different(self):
        """Three sidereal modes should all differ."""
        jd = 2451545.0
        results = []
        for mode in [SIDM_LAHIRI, SIDM_FAGAN_BRADLEY, SIDM_RAMAN]:
            swe.set_sid_mode(mode)
            r, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)
            results.append(r[0])

        # All three should be distinct
        for i in range(3):
            for j in range(i + 1, 3):
                diff = abs(results[i] - results[j])
                if diff > 180:
                    diff = 360 - diff
                assert diff > 0.01, f"Modes {i} and {j} too similar: diff={diff:.4f}°"


class TestAyanamsaState:
    """Test ayanamsa-related state functions."""

    @pytest.mark.unit
    def test_get_ayanamsa_ut_returns_float(self):
        """get_ayanamsa_ut returns a float."""
        swe.set_sid_mode(SIDM_LAHIRI)
        ayan = swe.get_ayanamsa_ut(2451545.0)
        assert type(ayan) is float

    @pytest.mark.unit
    def test_ayanamsa_positive(self):
        """Ayanamsa should be positive for modern dates."""
        swe.set_sid_mode(SIDM_LAHIRI)
        ayan = swe.get_ayanamsa_ut(2451545.0)
        assert ayan > 0, f"Ayanamsa {ayan} not positive"

    @pytest.mark.unit
    def test_ayanamsa_lahiri_approx(self):
        """Lahiri ayanamsa at J2000 should be ~23.85°."""
        swe.set_sid_mode(SIDM_LAHIRI)
        ayan = swe.get_ayanamsa_ut(2451545.0)
        assert 23 < ayan < 25, f"Lahiri ayanamsa {ayan}° (expected ~23.85)"

    @pytest.mark.unit
    def test_ayanamsa_increases_over_time(self):
        """Ayanamsa should increase over centuries (precession)."""
        swe.set_sid_mode(SIDM_LAHIRI)
        ayan_2000 = swe.get_ayanamsa_ut(2451545.0)
        ayan_2100 = swe.get_ayanamsa_ut(2451545.0 + 36525.0)
        assert ayan_2100 > ayan_2000, (
            f"Ayanamsa not increasing: {ayan_2000} -> {ayan_2100}"
        )

    @pytest.mark.unit
    def test_get_ayanamsa_ex_ut_returns_tuple(self):
        """get_ayanamsa_ex_ut returns (retflag, ayanamsa)."""
        swe.set_sid_mode(SIDM_LAHIRI)
        result = swe.get_ayanamsa_ex_ut(2451545.0, 0)
        assert len(result) == 2
        retflag, ayan = result
        assert isinstance(retflag, int)
        assert type(ayan) is float


class TestPlanetName:
    """Test get_planet_name."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,expected",
        [
            (0, "Sun"),
            (1, "Moon"),
            (2, "Mercury"),
            (3, "Venus"),
            (4, "Mars"),
            (5, "Jupiter"),
            (6, "Saturn"),
            (7, "Uranus"),
            (8, "Neptune"),
            (9, "Pluto"),
        ],
    )
    def test_planet_name(self, body_id: int, expected: str):
        """get_planet_name returns correct name."""
        name = swe.get_planet_name(body_id)
        assert expected.lower() in name.lower(), (
            f"Body {body_id}: got '{name}', expected '{expected}'"
        )

    @pytest.mark.unit
    def test_planet_name_returns_string(self):
        """get_planet_name always returns a string."""
        for body_id in range(10):
            name = swe.get_planet_name(body_id)
            assert isinstance(name, str)
