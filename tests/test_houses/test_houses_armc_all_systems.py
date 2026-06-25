"""
Tests for houses_armc with all house systems.

Verifies houses_armc produces valid cusps for all 25 systems,
consistency with houses() at equivalent ARMC, and edge cases.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.exceptions import PolarCircleError


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    swe.close()


JD_J2000 = 2451545.0

# All house system characters
ALL_SYSTEMS = "PKORECABMTUWVXHFSLQNYDIJ"
# Systems that may raise PolarCircleError at high latitudes
POLAR_SENSITIVE = "PKG"


class TestHousesArmcAllSystems:
    """Test houses_armc with all house systems."""

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys", list(ALL_SYSTEMS))
    def test_houses_armc_returns_valid(self, hsys):
        """houses_armc returns valid cusps for each system."""
        # Typical ARMC, latitude, obliquity
        armc = 30.0
        lat = 48.85
        eps = 23.44
        try:
            cusps, ascmc = swe.houses_armc(armc, lat, eps, ord(hsys))
            if hsys == "G":
                assert len(cusps) == 36
            else:
                assert len(cusps) == 12
            for i, c in enumerate(cusps):
                assert 0.0 <= c < 360.0, (
                    f"System {hsys} cusp {i + 1} = {c} out of range"
                )
        except PolarCircleError:
            pytest.skip(f"System {hsys} raised PolarCircleError at lat {lat}")

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys", list("PKORECABMTUWVX"))
    def test_houses_armc_at_equator(self, hsys):
        """houses_armc works at equator (lat=0)."""
        cusps, ascmc = swe.houses_armc(30.0, 0.0, 23.44, ord(hsys))
        assert len(cusps) == 12
        for c in cusps:
            assert 0.0 <= c < 360.0

    @pytest.mark.unit
    @pytest.mark.parametrize("armc", [0.0, 90.0, 180.0, 270.0, 359.99])
    def test_houses_armc_various_armc(self, armc):
        """houses_armc works at various ARMC values."""
        cusps, ascmc = swe.houses_armc(armc, 48.85, 23.44, ord("P"))
        assert len(cusps) == 12
        for c in cusps:
            assert 0.0 <= c < 360.0

    @pytest.mark.unit
    def test_houses_armc_ascmc_valid(self):
        """houses_armc ascmc tuple has valid angles."""
        cusps, ascmc = swe.houses_armc(30.0, 48.85, 23.44, ord("P"))
        # ascmc[0] = Asc, ascmc[1] = MC
        assert 0.0 <= ascmc[0] < 360.0
        assert 0.0 <= ascmc[1] < 360.0

    @pytest.mark.unit
    def test_mc_near_armc(self):
        """MC should be close to ARMC (modulo obliquity effects)."""
        armc = 120.0
        cusps, ascmc = swe.houses_armc(armc, 48.85, 23.44, ord("P"))
        mc = ascmc[1]
        # MC = atan(tan(ARMC) / cos(eps)) — close to ARMC
        diff = abs(mc - armc) % 360.0
        if diff > 180:
            diff = 360 - diff
        assert diff < 30.0, f"MC {mc} too far from ARMC {armc}"


class TestHousesArmcConsistency:
    """Test consistency between houses_armc and houses."""

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys", list("PKORCABMT"))
    def test_houses_armc_matches_houses(self, hsys):
        """houses_armc should match houses() when given the same ARMC."""
        # Get cusps from houses() to extract ARMC
        cusps_h, ascmc_h = swe.houses(JD_J2000, 48.85, 2.35, ord(hsys))
        armc = ascmc_h[2]  # ARMC is at index 2

        # Get obliquity
        eps = 23.44  # Approximate

        try:
            cusps_a, ascmc_a = swe.houses_armc(armc, 48.85, eps, ord(hsys))
            # Cusps should be close (not exact due to eps approximation)
            for i in range(min(len(cusps_h), len(cusps_a))):
                diff = abs(cusps_h[i] - cusps_a[i])
                if diff > 180:
                    diff = 360 - diff
                # Allow 1° tolerance due to obliquity approximation
                assert diff < 1.0, (
                    f"System {hsys} cusp {i + 1}: houses={cusps_h[i]}, armc={cusps_a[i]}"
                )
        except PolarCircleError:
            pass


class TestHousesArmcEx2:
    """Test houses_armc_ex2 with speeds."""

    @pytest.mark.unit
    def test_houses_armc_ex2_returns_4_tuples(self):
        """houses_armc_ex2 returns (cusps, ascmc, cusps_speed, ascmc_speed)."""
        result = swe.houses_armc_ex2(30.0, 48.85, 23.44, ord("P"))
        assert len(result) == 4
        cusps, ascmc, cusps_speed, ascmc_speed = result
        assert len(cusps) == 12
        assert len(cusps_speed) == 12

    @pytest.mark.unit
    def test_houses_armc_ex2_speeds_finite(self):
        """All speeds from houses_armc_ex2 are finite."""
        _, _, cusps_speed, ascmc_speed = swe.houses_armc_ex2(
            30.0, 48.85, 23.44, ord("P")
        )
        for i, s in enumerate(cusps_speed):
            assert math.isfinite(s), f"Cusp speed {i + 1} not finite: {s}"

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys", list("PKORCABMT"))
    def test_houses_armc_ex2_all_systems(self, hsys):
        """houses_armc_ex2 works with various systems."""
        try:
            result = swe.houses_armc_ex2(30.0, 48.85, 23.44, ord(hsys))
            assert len(result) == 4
        except PolarCircleError:
            pass


class TestHousesExtremeLat:
    """Test houses at equator and extreme latitudes."""

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys", list("EORABMUWVX"))
    def test_equator(self, hsys):
        """Houses at equator (lat=0)."""
        cusps, ascmc = swe.houses(JD_J2000, 0.0, 0.0, ord(hsys))
        assert len(cusps) == 12

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [60.0, 64.0, 65.0, 66.0])
    def test_high_latitude_placidus(self, lat):
        """Placidus at high latitudes near polar circle."""
        try:
            cusps, ascmc = swe.houses(JD_J2000, lat, 0.0, ord("P"))
            assert len(cusps) == 12
        except PolarCircleError:
            pass  # Expected near polar circle

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys", list("EORABMUWVX"))
    def test_lat_65(self, hsys):
        """Houses at lat=65 (near but below polar circle)."""
        cusps, ascmc = swe.houses(JD_J2000, 65.0, 25.0, ord(hsys))
        assert len(cusps) == 12

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys", list("EOWMVXN"))
    def test_polar_latitude(self, hsys):
        """Some house systems should work at polar latitudes."""
        # Equal, Whole Sign, Morinus, etc. should not need semi-arc
        try:
            cusps, ascmc = swe.houses(JD_J2000, 80.0, 25.0, ord(hsys))
            assert len(cusps) == 12
        except PolarCircleError:
            pass  # Some systems may still fail

    @pytest.mark.unit
    def test_southern_hemisphere(self):
        """Houses work in southern hemisphere."""
        cusps, ascmc = swe.houses(JD_J2000, -33.87, 151.21, ord("P"))
        assert len(cusps) == 12
        for c in cusps:
            assert 0.0 <= c < 360.0
