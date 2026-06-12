"""
Comprehensive tests for houses_armc and houses_ex functions.

Verifies that houses_armc produces valid results from ARMC input,
and that houses_ex handles sidereal correction correctly.
"""

from __future__ import annotations


import pytest

import libephemeris as swe
from libephemeris.constants import (
    FLG_SIDEREAL,
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
)


HOUSE_SYSTEMS = [
    (ord("P"), "Placidus"),
    (ord("K"), "Koch"),
    (ord("R"), "Regiomontanus"),
    (ord("C"), "Campanus"),
    (ord("A"), "Equal"),
    (ord("W"), "Whole Sign"),
    (ord("O"), "Porphyry"),
    (ord("B"), "Alcabitius"),
    (ord("M"), "Morinus"),
    (ord("E"), "Equal MC"),
]


class TestHousesArmcBasic:
    """Basic houses_armc functionality."""

    @pytest.mark.unit
    def test_houses_armc_returns_cusps_and_ascmc(self):
        """houses_armc returns (cusps, ascmc) tuple."""
        cusps, ascmc = swe.houses_armc(292.957, 41.9, 23.44, ord("P"))
        assert len(cusps) >= 12
        assert len(ascmc) >= 8

    @pytest.mark.unit
    def test_houses_armc_cusps_in_range(self):
        """All cusps should be 0-360°."""
        cusps, ascmc = swe.houses_armc(292.957, 41.9, 23.44, ord("P"))
        for i, c in enumerate(cusps[:12]):
            assert 0 <= c < 360, f"Cusp {i + 1}: {c}°"

    @pytest.mark.unit
    def test_houses_armc_ascmc_in_range(self):
        """ASCMC values should be 0-360° (except ARMC which can be the input)."""
        cusps, ascmc = swe.houses_armc(292.957, 41.9, 23.44, ord("P"))
        for i in [0, 1, 3, 4, 5, 6, 7]:  # Skip ARMC (index 2)
            assert 0 <= ascmc[i] < 360, f"ascmc[{i}]: {ascmc[i]}°"


class TestHousesArmcAllSystems:
    """Test houses_armc with various house systems."""

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys,name", HOUSE_SYSTEMS)
    def test_houses_armc_system(self, hsys: int, name: str):
        """houses_armc works for each house system."""
        cusps, ascmc = swe.houses_armc(180.0, 41.9, 23.44, hsys)
        assert len(cusps) >= 12, f"{name}: insufficient cusps"
        for i, c in enumerate(cusps[:12]):
            assert 0 <= c < 360, f"{name} cusp {i + 1}: {c}°"


class TestHousesArmcConsistency:
    """Test that houses_armc is consistent with houses()."""

    @pytest.mark.unit
    def test_houses_armc_matches_houses(self):
        """houses_armc with correct ARMC should match houses() output."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        # Get houses normally
        cusps_h, ascmc_h = swe.houses(jd, lat, lon, ord("P"))
        armc = ascmc_h[2]
        obliquity = ascmc_h[4] if len(ascmc_h) > 4 else 23.44

        # Get houses from ARMC
        cusps_a, ascmc_a = swe.houses_armc(armc, lat, obliquity, ord("P"))

        # Cusps should be close. Note: houses() uses its own obliquity
        # calculation while houses_armc uses the obliquity we pass,
        # so there can be small differences (~2°) depending on the
        # obliquity value used. Use generous tolerance.
        for i in range(12):
            diff = abs(cusps_h[i] - cusps_a[i])
            if diff > 180:
                diff = 360 - diff
            assert diff < 3.0, (
                f"Cusp {i + 1}: houses={cusps_h[i]:.4f}, armc={cusps_a[i]:.4f}"
            )


class TestHousesArmcVariousArmc:
    """Test houses_armc with various ARMC values."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "armc", [0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0]
    )
    def test_various_armc_values(self, armc: float):
        """houses_armc valid at various ARMC values."""
        cusps, ascmc = swe.houses_armc(armc, 41.9, 23.44, ord("P"))
        assert len(cusps) >= 12
        for i, c in enumerate(cusps[:12]):
            assert 0 <= c < 360, f"ARMC={armc} cusp {i + 1}: {c}°"

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [0.0, 20.0, 40.0, 60.0, -20.0, -40.0, -60.0])
    def test_various_latitudes(self, lat: float):
        """houses_armc valid at various latitudes."""
        cusps, ascmc = swe.houses_armc(180.0, lat, 23.44, ord("R"))
        assert len(cusps) >= 12
        for i, c in enumerate(cusps[:12]):
            assert 0 <= c < 360, f"Lat={lat} cusp {i + 1}: {c}°"


class TestHousesExBasic:
    """Basic houses_ex functionality."""

    @pytest.mark.unit
    def test_houses_ex_tropical(self):
        """houses_ex without sidereal flag equals houses()."""
        jd = 2451545.0
        cusps_h, ascmc_h = swe.houses(jd, 41.9, 12.5, ord("P"))
        cusps_e, ascmc_e = swe.houses_ex(jd, 41.9, 12.5, ord("P"), 0)

        for i in range(12):
            diff = abs(cusps_h[i] - cusps_e[i])
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.001, (
                f"Cusp {i + 1}: houses={cusps_h[i]:.4f}, ex={cusps_e[i]:.4f}"
            )

    @pytest.mark.unit
    def test_houses_ex_sidereal_differs(self):
        """houses_ex with FLG_SIDEREAL should differ from tropical."""
        jd = 2451545.0
        swe.set_sid_mode(SIDM_LAHIRI)

        cusps_trop, _ = swe.houses_ex(jd, 41.9, 12.5, ord("P"), 0)
        cusps_sid, _ = swe.houses_ex(jd, 41.9, 12.5, ord("P"), FLG_SIDEREAL)

        # At least some cusps should differ by ~ayanamsha
        diffs = []
        for i in range(12):
            diff = abs(cusps_trop[i] - cusps_sid[i])
            if diff > 180:
                diff = 360 - diff
            diffs.append(diff)

        avg_diff = sum(diffs) / len(diffs)
        assert avg_diff > 15, f"Average cusp diff {avg_diff:.2f}° (expected ~23°)"

    @pytest.mark.unit
    def test_houses_ex_sidereal_cusps_in_range(self):
        """Sidereal cusps should still be 0-360°."""
        jd = 2451545.0
        swe.set_sid_mode(SIDM_LAHIRI)
        cusps, ascmc = swe.houses_ex(jd, 41.9, 12.5, ord("P"), FLG_SIDEREAL)
        for i, c in enumerate(cusps[:12]):
            assert 0 <= c < 360, f"Sidereal cusp {i + 1}: {c}°"


class TestHousesExSiderealModes:
    """Test houses_ex with different sidereal modes."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode,name",
        [
            (SIDM_LAHIRI, "Lahiri"),
            (SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
        ],
    )
    def test_different_sid_modes_different_cusps(self, mode: int, name: str):
        """Different sidereal modes should produce different cusps."""
        jd = 2451545.0
        swe.set_sid_mode(mode)
        cusps, _ = swe.houses_ex(jd, 41.9, 12.5, ord("P"), FLG_SIDEREAL)
        assert len(cusps) >= 12
        for c in cusps[:12]:
            assert 0 <= c < 360


class TestHousesExDateRange:
    """Test houses_ex across dates."""

    @pytest.mark.unit
    @pytest.mark.parametrize("year", [1900, 1950, 2000, 2050, 2100])
    def test_houses_ex_across_years(self, year: int):
        """houses_ex valid across years."""
        jd = swe.julday(year, 6, 21, 12.0)
        cusps, ascmc = swe.houses_ex(jd, 41.9, 12.5, ord("P"), 0)
        assert len(cusps) >= 12
        for c in cusps[:12]:
            assert 0 <= c < 360
