"""
Comprehensive house system edge case tests.

Tests all 25 supported house systems at extreme latitudes, various
longitudes, and boundary conditions. Verifies no crashes, valid
output ranges, and correct cusp ordering.
"""

from __future__ import annotations

import math
import random

import pytest

import libephemeris as swe
from libephemeris.constants import FLG_SPEED
from libephemeris.exceptions import PolarCircleError


def _random_jds(n: int, seed: int = 42) -> list[float]:
    """Generate n random JDs in 1900-2100."""
    rng = random.Random(seed)
    return [rng.uniform(2415020.0, 2488070.0) for _ in range(n)]


# All 25 house systems (char -> name)
ALL_HOUSE_SYSTEMS = [
    ("P", "Placidus"),
    ("K", "Koch"),
    ("O", "Porphyry"),
    ("R", "Regiomontanus"),
    ("C", "Campanus"),
    ("E", "Equal (Asc)"),
    ("A", "Equal (Asc) alt"),
    ("W", "Whole Sign"),
    ("M", "Morinus"),
    ("B", "Alcabitius"),
    ("T", "Topocentric"),
    ("X", "Meridian"),
    ("V", "Vehlow Equal"),
    ("H", "Horizontal"),
    ("G", "Gauquelin"),
    ("U", "Krusinski-Pisa"),
    ("F", "Carter Poli-Equatorial"),
    ("Y", "APC Houses"),
    ("N", "Natural Graduation"),
    ("D", "Equal MC"),
    ("L", "Pullen SD"),
    ("Q", "Pullen SR"),
    ("S", "Sripati"),
    ("I", "Sunshine (Treindl)"),
    ("J", "Savard-A"),
]

EXTREME_LATITUDES = [
    (0.0, "Equator"),
    (23.44, "Tropic of Cancer"),
    (-23.44, "Tropic of Capricorn"),
    (45.0, "Mid-northern"),
    (-45.0, "Mid-southern"),
    (60.0, "Subarctic"),
    (-60.0, "Subantarctic"),
    (66.56, "Arctic Circle"),
    (-66.56, "Antarctic Circle"),
    (70.0, "Arctic 70N"),
    (-70.0, "Antarctic 70S"),
    (80.0, "Arctic 80N"),
    (-80.0, "Antarctic 80S"),
    (85.0, "Near North Pole"),
    (-85.0, "Near South Pole"),
    (89.0, "89N"),
    (-89.0, "89S"),
]


class TestAllHouseSystemsBasic:
    """Basic tests for all 25 house systems."""

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys,name", ALL_HOUSE_SYSTEMS)
    def test_house_system_returns_valid_cusps(self, hsys: str, name: str):
        """Each house system returns 12 valid cusps at a normal latitude."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5  # Rome
        cusps, ascmc = swe.houses(jd, lat, lon, ord(hsys))

        # Should have at least 12 cusps (Gauquelin has 36)
        n_cusps = 36 if hsys == "G" else 12
        assert len(cusps) >= n_cusps, (
            f"{name}: expected >={n_cusps} cusps, got {len(cusps)}"
        )

        for i in range(n_cusps):
            assert 0 <= cusps[i] < 360, (
                f"{name}: cusp {i + 1} = {cusps[i]} out of range"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys,name", ALL_HOUSE_SYSTEMS)
    def test_house_system_ascmc_valid(self, hsys: str, name: str):
        """Each house system returns valid ASCMC values."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5
        cusps, ascmc = swe.houses(jd, lat, lon, ord(hsys))

        # ascmc[0] = ASC, ascmc[1] = MC
        asc, mc = ascmc[0], ascmc[1]
        assert 0 <= asc < 360, f"{name}: ASC {asc} out of range"
        assert 0 <= mc < 360, f"{name}: MC {mc} out of range"


class TestExtremeLatitudes:
    """Test house systems at extreme latitudes."""

    @pytest.mark.unit
    @pytest.mark.parametrize("lat,lat_name", EXTREME_LATITUDES)
    @pytest.mark.parametrize("hsys,hsys_name", ALL_HOUSE_SYSTEMS)
    def test_no_crash_at_extreme_latitude(
        self, hsys: str, hsys_name: str, lat: float, lat_name: str
    ):
        """House system doesn't crash at extreme latitudes.

        Some systems (Placidus, Koch, Gauquelin) raise PolarCircleError
        at extreme latitudes — that's expected behaviour, not a crash.
        """
        jd = 2451545.0
        lon = 0.0
        try:
            cusps, ascmc = swe.houses(jd, lat, lon, ord(hsys))
        except PolarCircleError:
            # Expected for Placidus/Koch/Gauquelin at polar latitudes
            return

        n_cusps = 36 if hsys == "G" else 12
        assert len(cusps) >= n_cusps, (
            f"{hsys_name} @ {lat_name}: expected >={n_cusps} cusps"
        )

        for i in range(n_cusps):
            assert math.isfinite(cusps[i]), (
                f"{hsys_name} @ {lat_name}: cusp {i + 1} not finite"
            )
            assert 0 <= cusps[i] < 360, (
                f"{hsys_name} @ {lat_name}: cusp {i + 1} = {cusps[i]}"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,name",
        [
            ("O", "Porphyry"),
            ("R", "Regiomontanus"),
            ("E", "Equal"),
            ("W", "Whole Sign"),
        ],
    )
    def test_polar_circle_non_failing_systems(self, hsys: str, name: str):
        """Non-quadrant house systems work at all polar latitudes."""
        jd = 2451545.0
        for lat in [70, 75, 80, 85, 89]:
            cusps, ascmc = swe.houses(jd, lat, 0.0, ord(hsys))
            for i in range(12):
                assert 0 <= cusps[i] < 360, (
                    f"{name} @ {lat}°N: cusp {i + 1} = {cusps[i]}"
                )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,name",
        [
            ("P", "Placidus"),
            ("K", "Koch"),
        ],
    )
    def test_polar_circle_raises_for_quadrant_systems(self, hsys: str, name: str):
        """Placidus/Koch correctly raise PolarCircleError at extreme latitudes."""
        jd = 2451545.0
        with pytest.raises(PolarCircleError):
            swe.houses(jd, 89.0, 0.0, ord(hsys))


class TestHouseSystemConsistency:
    """Test consistency properties of house systems."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,name",
        [
            ("E", "Equal"),
            ("A", "Equal alt"),
            ("W", "Whole Sign"),
            ("V", "Vehlow"),
            ("D", "Equal MC"),
        ],
    )
    def test_equal_house_30_degree_spacing(self, hsys: str, name: str):
        """Equal-type houses should have exactly 30° spacing."""
        jd = 2451545.0
        cusps, _ = swe.houses(jd, 41.9, 12.5, ord(hsys))

        for i in range(11):
            spacing = (cusps[i + 1] - cusps[i]) % 360
            assert abs(spacing - 30.0) < 0.001, (
                f"{name}: cusps {i + 1}->{i + 2} spacing {spacing}° != 30°"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys,name", ALL_HOUSE_SYSTEMS)
    def test_cusps_ascending_order(self, hsys: str, name: str):
        """Cusps should be in ascending order (mod 360)."""
        jd = 2451545.0
        cusps, _ = swe.houses(jd, 41.9, 12.5, ord(hsys))
        n = 36 if hsys == "G" else 12

        for i in range(n - 1):
            diff = (cusps[i + 1] - cusps[i]) % 360
            assert diff > 0, (
                f"{name}: cusp {i + 1}={cusps[i]:.2f} -> "
                f"cusp {i + 2}={cusps[i + 1]:.2f} not ascending"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,name",
        [
            ("P", "Placidus"),
            ("K", "Koch"),
            ("O", "Porphyry"),
            ("R", "Regiomontanus"),
            ("C", "Campanus"),
            ("B", "Alcabitius"),
            ("T", "Topocentric"),
            ("U", "Krusinski-Pisa"),
        ],
    )
    def test_cusp1_equals_asc(self, hsys: str, name: str):
        """For quadrant house systems, cusp 1 should equal ASC."""
        jd = 2451545.0
        cusps, ascmc = swe.houses(jd, 41.9, 12.5, ord(hsys))
        asc = ascmc[0]
        assert abs(cusps[0] - asc) < 0.001, (
            f"{name}: cusp 1 {cusps[0]:.4f} != ASC {asc:.4f}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,name",
        [
            ("P", "Placidus"),
            ("K", "Koch"),
            ("O", "Porphyry"),
            ("R", "Regiomontanus"),
            ("C", "Campanus"),
            ("B", "Alcabitius"),
            ("T", "Topocentric"),
        ],
    )
    def test_cusp10_equals_mc(self, hsys: str, name: str):
        """For quadrant house systems, cusp 10 should equal MC."""
        jd = 2451545.0
        cusps, ascmc = swe.houses(jd, 41.9, 12.5, ord(hsys))
        mc = ascmc[1]
        assert abs(cusps[9] - mc) < 0.001, (
            f"{name}: cusp 10 {cusps[9]:.4f} != MC {mc:.4f}"
        )


class TestMultipleDatesAndLocations:
    """Test house systems across multiple dates and locations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,name",
        [
            ("P", "Placidus"),
            ("E", "Equal"),
            ("W", "Whole Sign"),
            ("O", "Porphyry"),
            ("K", "Koch"),
        ],
    )
    def test_common_systems_50_dates(self, hsys: str, name: str):
        """Common systems produce valid results at 50 random dates."""
        jds = _random_jds(50, seed=ord(hsys) * 7)
        for jd in jds:
            cusps, ascmc = swe.houses(jd, 41.9, 12.5, ord(hsys))
            for i in range(12):
                assert 0 <= cusps[i] < 360, (
                    f"{name} @ JD {jd:.1f}: cusp {i + 1}={cusps[i]}"
                )

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys,name", ALL_HOUSE_SYSTEMS)
    def test_southern_hemisphere(self, hsys: str, name: str):
        """All systems work in the southern hemisphere."""
        jd = 2451545.0
        # Sydney, Australia
        cusps, ascmc = swe.houses(jd, -33.87, 151.21, ord(hsys))
        n = 36 if hsys == "G" else 12
        for i in range(n):
            assert 0 <= cusps[i] < 360, f"{name} @ Sydney: cusp {i + 1}={cusps[i]}"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,name",
        [
            ("P", "Placidus"),
            ("K", "Koch"),
            ("O", "Porphyry"),
            ("R", "Regiomontanus"),
            ("E", "Equal"),
            ("W", "Whole Sign"),
        ],
    )
    def test_polar_circle_common_systems(self, hsys: str, name: str):
        """Common house systems work at polar latitudes (70-89°).

        Placidus and Koch raise PolarCircleError at extreme latitudes —
        that's correct behaviour. Other systems should work fine.
        """
        jd = 2451545.0
        for lat in [70, 75, 80, 85, 89]:
            try:
                cusps, ascmc = swe.houses(jd, lat, 0.0, ord(hsys))
            except PolarCircleError:
                # Expected for Placidus/Koch at high latitudes
                continue
            for i in range(12):
                assert 0 <= cusps[i] < 360, (
                    f"{name} @ {lat}°N: cusp {i + 1} = {cusps[i]}"
                )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "lon",
        [0, 30, 60, 90, 120, 150, 180, -30, -60, -90, -120, -150],
    )
    def test_placidus_various_longitudes(self, lon: int):
        """Placidus works at various geographic longitudes."""
        jd = 2451545.0
        cusps, ascmc = swe.houses(jd, 45.0, float(lon), ord("P"))
        for i in range(12):
            assert 0 <= cusps[i] < 360, f"Placidus @ lon={lon}: cusp {i + 1}={cusps[i]}"


class TestHousesArmc:
    """Tests for houses_armc function."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,name",
        [
            ("P", "Placidus"),
            ("K", "Koch"),
            ("O", "Porphyry"),
            ("E", "Equal"),
            ("W", "Whole Sign"),
            ("R", "Regiomontanus"),
            ("C", "Campanus"),
            ("J", "Savard-A"),
        ],
    )
    def test_houses_armc_valid(self, hsys: str, name: str):
        """houses_armc returns valid cusps for common systems."""
        armc = 120.0
        eps = 23.44
        lat = 41.9
        cusps, ascmc = swe.houses_armc(armc, lat, eps, ord(hsys))
        n = 36 if hsys == "G" else 12
        for i in range(n):
            assert 0 <= cusps[i] < 360, f"{name}: cusp {i + 1}={cusps[i]}"
