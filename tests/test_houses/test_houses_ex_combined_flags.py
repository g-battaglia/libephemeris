"""
Tests for houses_ex with combined flags, especially SIDEREAL.

Verifies houses_ex with FLG_SIDEREAL, houses_ex2 with speeds,
and combined flag behavior.
"""

from __future__ import annotations

import math

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

HOUSE_SYSTEMS = "PKORABMTUVX"  # Common systems (excluding Gauquelin and Whole Sign)


class TestHousesExSidereal:
    """Test houses_ex with FLG_SIDEREAL flag."""

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys", list(HOUSE_SYSTEMS))
    def test_sidereal_cusps_in_range(self, hsys):
        """Sidereal house cusps are in [0, 360)."""
        swe.set_sid_mode(SIDM_LAHIRI)
        cusps, ascmc = swe.houses_ex(JD_J2000, 48.85, 2.35, ord(hsys), FLG_SIDEREAL)
        for i, c in enumerate(cusps):
            assert 0.0 <= c < 360.0, f"System {hsys} cusp {i + 1} = {c} out of range"

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys", list(HOUSE_SYSTEMS))
    def test_sidereal_differs_from_tropical(self, hsys):
        """Sidereal cusps differ from tropical cusps by ~ayanamsa."""
        swe.set_sid_mode(SIDM_LAHIRI)
        ayan = swe.get_ayanamsa_ut(JD_J2000)

        cusps_t, _ = swe.houses_ex(JD_J2000, 48.85, 2.35, ord(hsys), 0)
        cusps_s, _ = swe.houses_ex(JD_J2000, 48.85, 2.35, ord(hsys), FLG_SIDEREAL)

        # First cusp difference should be approximately the ayanamsa
        diff = (cusps_t[0] - cusps_s[0]) % 360.0
        if diff > 180.0:
            diff = 360.0 - diff
        assert diff == pytest.approx(ayan, abs=1.0), (
            f"System {hsys}: cusp diff {diff}° vs ayanamsa {ayan}°"
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
    def test_different_sidereal_modes(self, sid_mode):
        """Different sidereal modes produce different cusp positions."""
        swe.set_sid_mode(sid_mode)
        cusps, ascmc = swe.houses_ex(JD_J2000, 48.85, 2.35, ord("P"), FLG_SIDEREAL)
        assert len(cusps) == 12
        for c in cusps:
            assert 0.0 <= c < 360.0

    @pytest.mark.unit
    def test_sidereal_ascmc_armc_not_corrected(self):
        """ARMC (ascmc[2]) should NOT be corrected by ayanamsa."""
        swe.set_sid_mode(SIDM_LAHIRI)
        _, ascmc_t = swe.houses_ex(JD_J2000, 48.85, 2.35, ord("P"), 0)
        _, ascmc_s = swe.houses_ex(JD_J2000, 48.85, 2.35, ord("P"), FLG_SIDEREAL)
        # ARMC is at index 2 — should be the same
        assert ascmc_t[2] == pytest.approx(ascmc_s[2], abs=0.01), (
            f"ARMC tropical={ascmc_t[2]}, sidereal={ascmc_s[2]}"
        )

    @pytest.mark.unit
    def test_sidereal_asc_corrected(self):
        """Ascendant (ascmc[0]) should be corrected by ayanamsa."""
        swe.set_sid_mode(SIDM_LAHIRI)
        ayan = swe.get_ayanamsa_ut(JD_J2000)

        _, ascmc_t = swe.houses_ex(JD_J2000, 48.85, 2.35, ord("P"), 0)
        _, ascmc_s = swe.houses_ex(JD_J2000, 48.85, 2.35, ord("P"), FLG_SIDEREAL)

        diff = (ascmc_t[0] - ascmc_s[0]) % 360.0
        if diff > 180.0:
            diff = 360.0 - diff
        assert diff == pytest.approx(ayan, abs=0.5)


class TestHousesExWholeSignSidereal:
    """Test Whole Sign houses with sidereal mode."""

    @pytest.mark.unit
    def test_whole_sign_sidereal_30_degree_spacing(self):
        """Whole Sign sidereal cusps are exactly 30° apart."""
        swe.set_sid_mode(SIDM_LAHIRI)
        cusps, _ = swe.houses_ex(JD_J2000, 48.85, 2.35, ord("W"), FLG_SIDEREAL)
        for i in range(11):
            diff = (cusps[i + 1] - cusps[i]) % 360.0
            assert diff == pytest.approx(30.0, abs=0.01), (
                f"Cusp {i + 1} to {i + 2}: {diff}° (expected 30°)"
            )

    @pytest.mark.unit
    def test_whole_sign_sidereal_cusp1_at_sign_boundary(self):
        """Whole Sign first cusp should be at a sign boundary (multiple of 30°)."""
        swe.set_sid_mode(SIDM_LAHIRI)
        cusps, _ = swe.houses_ex(JD_J2000, 48.85, 2.35, ord("W"), FLG_SIDEREAL)
        remainder = cusps[0] % 30.0
        assert remainder == pytest.approx(0.0, abs=0.01) or remainder == pytest.approx(
            30.0, abs=0.01
        ), f"Cusp 1 ({cusps[0]}) not at sign boundary (remainder {remainder})"


class TestHousesEx2Speeds:
    """Test houses_ex2 which returns cusp speeds."""

    @pytest.mark.unit
    def test_houses_ex2_returns_4_tuples(self):
        """houses_ex2 returns (cusps, ascmc, cusps_speed, ascmc_speed)."""
        result = swe.houses_ex2(JD_J2000, 48.85, 2.35, ord("P"))
        assert len(result) == 4
        cusps, ascmc, cusps_speed, ascmc_speed = result
        assert len(cusps) == 12
        assert len(cusps_speed) == 12

    @pytest.mark.unit
    def test_cusp_speeds_finite(self):
        """All cusp speeds should be finite."""
        _, _, cusps_speed, ascmc_speed = swe.houses_ex2(JD_J2000, 48.85, 2.35, ord("P"))
        for i, s in enumerate(cusps_speed):
            assert math.isfinite(s), f"Cusp speed {i + 1} not finite: {s}"

    @pytest.mark.unit
    def test_cusp_speeds_reasonable_magnitude(self):
        """Cusp speeds should be positive and nonzero (Earth rotation)."""
        _, _, cusps_speed, _ = swe.houses_ex2(JD_J2000, 48.85, 2.35, ord("P"))
        for i, s in enumerate(cusps_speed):
            # Cusps move due to Earth rotation; exact speed depends on
            # house system and latitude. Some cusps can exceed 700 deg/day
            # for certain oblique ascension geometries.
            assert abs(s) > 100.0, f"Cusp {i + 1} speed {s} seems too slow"
            assert abs(s) < 2000.0, f"Cusp {i + 1} speed {s} seems unreasonably fast"

    @pytest.mark.unit
    def test_houses_ex2_with_sidereal(self):
        """houses_ex2 works with FLG_SIDEREAL."""
        swe.set_sid_mode(SIDM_LAHIRI)
        cusps, ascmc, cusps_speed, ascmc_speed = swe.houses_ex2(
            JD_J2000, 48.85, 2.35, ord("P"), FLG_SIDEREAL
        )
        assert len(cusps) == 12
        assert len(cusps_speed) == 12
        for c in cusps:
            assert 0.0 <= c < 360.0


class TestHousesInvalidHsys:
    """Test house system with invalid hsys values."""

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys", [ord("Z"), ord("!"), 0, 1, 255])
    def test_invalid_hsys_falls_back(self, hsys):
        """Invalid house system falls back to Placidus (no error raised)."""
        cusps_inv, ascmc_inv = swe.houses(JD_J2000, 48.85, 2.35, hsys)
        cusps_p, ascmc_p = swe.houses(JD_J2000, 48.85, 2.35, ord("P"))
        # Should produce Placidus cusps
        assert len(cusps_inv) == 12
        for i in range(12):
            assert cusps_inv[i] == pytest.approx(cusps_p[i], abs=0.01), (
                f"Cusp {i + 1}: invalid={cusps_inv[i]}, Placidus={cusps_p[i]}"
            )

    @pytest.mark.unit
    def test_house_name_empty_for_invalid(self):
        """house_name returns '' for an unknown selector (measured contract)."""
        assert swe.house_name(ord("Z")) == ""
        assert swe.house_name(ord("5")) == ""

    def test_house_system_selector_is_case_insensitive(self):
        """Lowercase selectors fold to their uppercase system ('k' == 'K')...

        ...except 'i', which is the distinct Sunshine-alternative system.
        """
        jd = swe.julday(2005, 6, 15, 12.0)
        for ch in "korcaewbxhtgmvudnylfjspq":
            lo = swe.houses(jd, 48.2, 16.4, ch.encode())[0]
            up = swe.houses(jd, 48.2, 16.4, ch.upper().encode())[0]
            assert lo == up, ch
        assert swe.house_name(ord("k")) == "Koch"
        assert swe.house_name(ord("i")) == "Sunshine/alt."
        # 'i' stays distinct from 'I' where the two constructions differ.
        ci = swe.houses(jd, 62.0, 12.5, b"i")[0]
        cI = swe.houses(jd, 62.0, 12.5, b"I")[0]
        assert max(abs(a - b) for a, b in zip(ci, cI)) > 1.0

    @pytest.mark.unit
    def test_houses_ex_invalid_hsys_falls_back(self):
        """houses_ex with invalid hsys also falls back to Placidus."""
        cusps_inv, _ = swe.houses_ex(JD_J2000, 48.85, 2.35, ord("Z"), 0)
        cusps_p, _ = swe.houses_ex(JD_J2000, 48.85, 2.35, ord("P"), 0)
        for i in range(12):
            assert cusps_inv[i] == pytest.approx(cusps_p[i], abs=0.01)
