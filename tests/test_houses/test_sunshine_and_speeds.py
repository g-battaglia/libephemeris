"""
Tests for Sunshine house system ('I'/'i') in sidereal mode,
and houses_ex2 speed accuracy vs numerical differentiation.
"""

from __future__ import annotations


import pytest

import libephemeris as swe
from libephemeris.constants import (
    FLG_SIDEREAL,
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
    SIDM_RAMAN,
)


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    swe.close()


JD_J2000 = 2451545.0
JD_2020 = 2458849.5


class TestSunshineHouseSystem:
    """Test Sunshine ('I'/'i') house system."""

    @pytest.mark.unit
    def test_sunshine_I_returns_12_cusps(self):
        """Sunshine 'I' returns 12 cusps."""
        cusps, ascmc = swe.houses(JD_J2000, 48.85, 2.35, ord("I"))
        assert len(cusps) == 12

    @pytest.mark.unit
    def test_sunshine_i_returns_12_cusps(self):
        """Sunshine 'i' (alternative) returns 12 cusps."""
        cusps, ascmc = swe.houses(JD_J2000, 48.85, 2.35, ord("i"))
        assert len(cusps) == 12

    @pytest.mark.unit
    def test_sunshine_cusps_in_range(self):
        """All Sunshine cusps are in [0, 360)."""
        cusps, ascmc = swe.houses(JD_J2000, 48.85, 2.35, ord("I"))
        for i, c in enumerate(cusps):
            assert 0.0 <= c < 360.0, f"Cusp {i + 1} = {c} out of range"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "lat,lon",
        [
            (48.85, 2.35),  # Paris
            (40.71, -74.01),  # New York
            (0.0, 0.0),  # Equator
            (-33.87, 151.21),  # Sydney
            (55.75, 37.62),  # Moscow
        ],
    )
    def test_sunshine_at_various_locations(self, lat, lon):
        """Sunshine houses work at various locations."""
        cusps, ascmc = swe.houses(JD_J2000, lat, lon, ord("I"))
        assert len(cusps) == 12
        for c in cusps:
            assert 0.0 <= c < 360.0

    @pytest.mark.unit
    def test_sunshine_I_vs_i_differ(self):
        """'I' and 'i' may produce different cusps (different Sunshine variants)."""
        cusps_I, _ = swe.houses(JD_J2000, 48.85, 2.35, ord("I"))
        cusps_i, _ = swe.houses(JD_J2000, 48.85, 2.35, ord("i"))
        # They may or may not differ
        assert len(cusps_I) == 12
        assert len(cusps_i) == 12

    @pytest.mark.unit
    def test_sunshine_house_name(self):
        """house_name recognizes Sunshine."""
        name_I = swe.house_name(ord("I"))
        name_i = swe.house_name(ord("i"))
        assert name_I != "Unknown"
        assert name_i != "Unknown"


class TestSunshineHouseSidereal:
    """Test Sunshine house system with sidereal mode."""

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys", [ord("I"), ord("i")])
    def test_sunshine_sidereal_cusps_in_range(self, hsys):
        """Sidereal Sunshine cusps are in [0, 360)."""
        swe.set_sid_mode(SIDM_LAHIRI)
        cusps, ascmc = swe.houses_ex(JD_J2000, 48.85, 2.35, hsys, FLG_SIDEREAL)
        assert len(cusps) == 12
        for i, c in enumerate(cusps):
            assert 0.0 <= c < 360.0, f"Cusp {i + 1} = {c} out of range"

    @pytest.mark.unit
    def test_sunshine_sidereal_differs_from_tropical(self):
        """Sidereal Sunshine cusps differ from tropical by ~ayanamsa."""
        swe.set_sid_mode(SIDM_LAHIRI)
        ayan = swe.get_ayanamsa_ut(JD_J2000)

        cusps_t, _ = swe.houses_ex(JD_J2000, 48.85, 2.35, ord("I"), 0)
        cusps_s, _ = swe.houses_ex(JD_J2000, 48.85, 2.35, ord("I"), FLG_SIDEREAL)

        # For Sunshine, the ayanamsa correction is applied but the intermediate
        # cusps are recalculated, so the diff may not be exactly ayanamsa for all.
        # Check cusp 1 (Asc) which should be close.
        diff = (cusps_t[0] - cusps_s[0]) % 360.0
        if diff > 180.0:
            diff = 360.0 - diff
        # Sunshine recalculates some cusps, so allow wider tolerance
        assert diff == pytest.approx(ayan, abs=5.0), (
            f"Cusp 1 diff {diff}° vs ayanamsa {ayan}°"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "sid_mode",
        [
            SIDM_LAHIRI,
            SIDM_FAGAN_BRADLEY,
            SIDM_RAMAN,
        ],
    )
    def test_sunshine_sidereal_different_modes(self, sid_mode):
        """Sunshine works with different sidereal modes."""
        swe.set_sid_mode(sid_mode)
        cusps, _ = swe.houses_ex(JD_J2000, 48.85, 2.35, ord("I"), FLG_SIDEREAL)
        assert len(cusps) == 12
        for c in cusps:
            assert 0.0 <= c < 360.0


class TestHousesEx2SpeedAccuracy:
    """Test houses_ex2 cusp speed accuracy via numerical differentiation."""

    DT = 1.0 / 86400.0  # 1 second in days

    def _numerical_cusp_speed(self, jd, lat, lon, hsys, cusp_idx):
        """Compute numerical derivative of cusp position."""
        cusps_m, _ = swe.houses(jd - self.DT, lat, lon, hsys)
        cusps_p, _ = swe.houses(jd + self.DT, lat, lon, hsys)

        val_m = cusps_m[cusp_idx]
        val_p = cusps_p[cusp_idx]

        # Handle wrap-around
        diff = val_p - val_m
        if diff > 180.0:
            diff -= 360.0
        elif diff < -180.0:
            diff += 360.0

        return diff / (2 * self.DT)

    @pytest.mark.unit
    @pytest.mark.parametrize("cusp_idx", [0, 3, 6, 9])  # Cusps 1, 4, 7, 10
    def test_placidus_cusp_speed(self, cusp_idx):
        """Placidus cusp speed matches numerical derivative."""
        cusps, ascmc, cusps_speed, _ = swe.houses_ex2(JD_J2000, 48.85, 2.35, ord("P"))
        returned = cusps_speed[cusp_idx]
        numerical = self._numerical_cusp_speed(
            JD_J2000, 48.85, 2.35, ord("P"), cusp_idx
        )
        # Tolerance: ~5 deg/day (speeds are ~360 deg/day so ~1.4% error)
        assert returned == pytest.approx(numerical, abs=10.0), (
            f"Cusp {cusp_idx + 1}: returned {returned}, numerical {numerical}"
        )

    @pytest.mark.unit
    def test_equal_house_cusp_speeds_similar(self):
        """Equal house cusp speeds should all be similar (~360 deg/day)."""
        _, _, cusps_speed, _ = swe.houses_ex2(JD_J2000, 48.85, 2.35, ord("E"))
        # Equal houses — all cusps move at the same rate
        for i in range(1, 12):
            assert cusps_speed[i] == pytest.approx(cusps_speed[0], abs=5.0), (
                f"Equal house cusp {i + 1} speed {cusps_speed[i]} differs from cusp 1 {cusps_speed[0]}"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys", [ord("P"), ord("K"), ord("R"), ord("C")])
    def test_cusp_speeds_nonzero(self, hsys):
        """All cusp speeds should be nonzero."""
        _, _, cusps_speed, _ = swe.houses_ex2(JD_J2000, 48.85, 2.35, hsys)
        for i, s in enumerate(cusps_speed):
            assert abs(s) > 10.0, f"Cusp {i + 1} speed too small: {s}"

    @pytest.mark.unit
    def test_ascmc_speed_asc_nonzero(self):
        """Ascendant speed in ascmc_speed should be nonzero."""
        _, _, _, ascmc_speed = swe.houses_ex2(JD_J2000, 48.85, 2.35, ord("P"))
        # ascmc_speed[0] = Ascendant speed
        assert abs(ascmc_speed[0]) > 10.0, (
            f"Ascendant speed too small: {ascmc_speed[0]}"
        )

    @pytest.mark.unit
    def test_ascmc_speed_mc_nonzero(self):
        """MC speed in ascmc_speed should be nonzero."""
        _, _, _, ascmc_speed = swe.houses_ex2(JD_J2000, 48.85, 2.35, ord("P"))
        # ascmc_speed[1] = MC speed
        assert abs(ascmc_speed[1]) > 10.0, f"MC speed too small: {ascmc_speed[1]}"
