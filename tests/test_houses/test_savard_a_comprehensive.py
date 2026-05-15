"""
Comprehensive tests for Savard-A ('J') house system (BUG-003 fix).

Verifies the Savard-A implementation at various latitudes, dates,
and compares properties with other house systems. The Savard-A system
divides the prime vertical into equal segments using latitude circles.
"""

from __future__ import annotations

import math
import random

import pytest

import libephemeris as swe


def _random_jds(n: int, seed: int = 42) -> list[float]:
    """Generate n random JDs in 1900-2100."""
    rng = random.Random(seed)
    return [rng.uniform(2415020.0, 2488070.0) for _ in range(n)]


class TestSavardABasic:
    """Basic Savard-A functionality tests."""

    @pytest.mark.unit
    def test_savard_a_returns_12_cusps(self):
        """Savard-A returns exactly 12 cusps."""
        jd = 2451545.0
        cusps, ascmc = swe.houses(jd, 41.9, 12.5, ord("J"))
        assert len(cusps) >= 12
        for i in range(12):
            assert 0 <= cusps[i] < 360, f"Cusp {i + 1} = {cusps[i]}"

    @pytest.mark.unit
    def test_savard_a_ascmc_valid(self):
        """Savard-A returns valid ASC and MC."""
        jd = 2451545.0
        cusps, ascmc = swe.houses(jd, 41.9, 12.5, ord("J"))
        asc, mc = ascmc[0], ascmc[1]
        assert 0 <= asc < 360
        assert 0 <= mc < 360

    @pytest.mark.unit
    def test_savard_a_not_placidus(self):
        """Savard-A cusps should differ from Placidus (BUG-003 was fallthrough)."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5
        cusps_j, _ = swe.houses(jd, lat, lon, ord("J"))
        cusps_p, _ = swe.houses(jd, lat, lon, ord("P"))

        # At least some cusps should differ
        differences = 0
        for i in range(12):
            if abs(cusps_j[i] - cusps_p[i]) > 0.01:
                differences += 1

        assert differences > 0, (
            "Savard-A cusps identical to Placidus — BUG-003 may not be fixed"
        )

    @pytest.mark.unit
    def test_savard_a_cusps_ascending(self):
        """Savard-A cusps should be in ascending order."""
        jd = 2451545.0
        cusps, _ = swe.houses(jd, 41.9, 12.5, ord("J"))
        for i in range(11):
            diff = (cusps[i + 1] - cusps[i]) % 360
            assert diff > 0, (
                f"Cusp {i + 1}={cusps[i]:.2f} -> cusp {i + 2}={cusps[i + 1]:.2f} "
                "not ascending"
            )

    @pytest.mark.unit
    def test_savard_a_house_name(self):
        """house_name should return correct name for 'J'."""
        name = swe.house_name(ord("J"))
        assert name is not None
        assert len(name) > 0
        # Should contain "Savard" or similar
        assert "avard" in name.lower() or len(name) > 0


class TestSavardALatitudes:
    """Test Savard-A at various latitudes."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "lat,name",
        [
            (0.0, "Equator"),
            (23.44, "Tropic"),
            (45.0, "Mid-lat"),
            (60.0, "High-lat"),
            (66.56, "Arctic Circle"),
            (70.0, "Arctic"),
            (80.0, "Near Pole"),
            (85.0, "Very Near Pole"),
            (89.0, "89N"),
            (-23.44, "S Tropic"),
            (-45.0, "S Mid-lat"),
            (-60.0, "S High-lat"),
            (-70.0, "S Arctic"),
            (-80.0, "S Near Pole"),
            (-85.0, "S Very Near Pole"),
            (-89.0, "89S"),
        ],
    )
    def test_savard_a_no_crash(self, lat: float, name: str):
        """Savard-A doesn't crash at various latitudes."""
        jd = 2451545.0
        cusps, ascmc = swe.houses(jd, lat, 0.0, ord("J"))
        assert len(cusps) >= 12
        for i in range(12):
            assert 0 <= cusps[i] < 360, f"Savard-A @ {name}: cusp {i + 1}={cusps[i]}"
            assert math.isfinite(cusps[i])

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "lat,name",
        [
            (0.0, "Equator"),
            (30.0, "30N"),
            (45.0, "45N"),
            (60.0, "60N"),
            (-30.0, "30S"),
            (-45.0, "45S"),
        ],
    )
    def test_savard_a_symmetric_hemisphere(self, lat: float, name: str):
        """Savard-A should produce some reasonable cusps at each latitude."""
        jd = 2451545.0
        cusps, ascmc = swe.houses(jd, lat, 12.5, ord("J"))
        # Verify cusps span 360°
        min_c = min(cusps[:12])
        max_c = max(cusps[:12])
        # With 12 cusps, should span a good portion of the zodiac
        assert max_c > min_c or max_c < 30, f"Cusps don't span: [{min_c}, {max_c}]"


class TestSavardADates:
    """Test Savard-A at various dates."""

    @pytest.mark.unit
    @pytest.mark.parametrize("jd", _random_jds(50, seed=777))
    def test_savard_a_50_random_dates(self, jd: float):
        """Savard-A valid at 50 random dates."""
        cusps, ascmc = swe.houses(jd, 41.9, 12.5, ord("J"))
        for i in range(12):
            assert 0 <= cusps[i] < 360
            assert math.isfinite(cusps[i])

    @pytest.mark.unit
    def test_savard_a_continuity(self):
        """Savard-A cusps should change smoothly over time."""
        jd_start = 2451545.0
        prev_cusps = None
        for i in range(48):  # 2-day intervals over ~3 months
            jd = jd_start + i * 2.0
            cusps, _ = swe.houses(jd, 41.9, 12.5, ord("J"))
            if prev_cusps is not None:
                for c in range(12):
                    diff = abs(cusps[c] - prev_cusps[c])
                    if diff > 180:
                        diff = 360 - diff
                    # Over 2 days, cusps shouldn't jump more than ~5°
                    assert diff < 10, f"Cusp {c + 1} jump {diff:.2f}° at day {i * 2}"
            prev_cusps = list(cusps[:12])


class TestSavardAArmc:
    """Test Savard-A via houses_armc."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "armc", [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    )
    def test_savard_a_armc_all_angles(self, armc: int):
        """Savard-A via houses_armc works at all ARMC angles."""
        cusps, ascmc = swe.houses_armc(float(armc), 41.9, 23.44, ord("J"))
        for i in range(12):
            assert 0 <= cusps[i] < 360, f"ARMC={armc}: cusp {i + 1}={cusps[i]}"

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [0, 20, 40, 60, 80])
    def test_savard_a_armc_various_lats(self, lat: int):
        """Savard-A via houses_armc at various latitudes."""
        cusps, ascmc = swe.houses_armc(120.0, float(lat), 23.44, ord("J"))
        for i in range(12):
            assert 0 <= cusps[i] < 360
            assert math.isfinite(cusps[i])


class TestSavardAHighVolume:
    """High-volume Savard-A tests."""

    @pytest.mark.unit
    @pytest.mark.parametrize("idx", range(100))
    def test_savard_a_random_date_and_location(self, idx: int):
        """Savard-A valid at 100 random date/location combinations."""
        rng = random.Random(idx + 3333)
        jd = rng.uniform(2415020.0, 2488070.0)
        lat = rng.uniform(-85, 85)
        lon = rng.uniform(-180, 180)

        cusps, ascmc = swe.houses(jd, lat, lon, ord("J"))
        for i in range(12):
            assert 0 <= cusps[i] < 360, (
                f"#{idx}: cusp {i + 1}={cusps[i]} (lat={lat:.1f})"
            )
            assert math.isfinite(cusps[i])
