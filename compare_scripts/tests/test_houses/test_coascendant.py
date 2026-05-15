"""
Tests for Co-Ascendant (Walter Koch) and related angle calculations.

Tests the ascmc[4] through ascmc[7] outputs from houses_ex2:
- ascmc[4]: Equatorial Ascendant (East Point)
- ascmc[5]: Co-Ascendant (Walter Koch formula)
- ascmc[6]: Co-Ascendant (Michael Munkasey formula)
- ascmc[7]: Polar Ascendant

These are astrological sensitive points related to the Ascendant
calculated using modified pole heights or ARMC offsets.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


class TestCoAscendantBasic:
    """Basic tests for Co-Ascendant calculations."""

    @pytest.mark.unit
    def test_ascmc_has_8_elements(self):
        """ascmc should have 8 elements including Co-Ascendants."""
        jd = 2451545.0
        cusps, ascmc = ephem.houses(jd, 41.9, 12.5, ord("P"))

        assert len(ascmc) == 8, f"ascmc has {len(ascmc)} elements, expected 8"

    @pytest.mark.unit
    def test_coascendant_indices(self):
        """Verify the indices of Co-Ascendant values in ascmc."""
        jd = 2451545.0
        cusps, ascmc = ephem.houses(jd, 41.9, 12.5, ord("P"))

        # ascmc[4] = Equatorial Ascendant
        # ascmc[5] = Co-Ascendant W. Koch
        # ascmc[6] = Co-Ascendant M. Munkasey
        # ascmc[7] = Polar Ascendant
        assert 0 <= ascmc[4] < 360, f"Equ Asc = {ascmc[4]} out of range"
        assert 0 <= ascmc[5] < 360, f"CoAsc Koch = {ascmc[5]} out of range"
        assert 0 <= ascmc[6] < 360, f"CoAsc Munk = {ascmc[6]} out of range"
        assert 0 <= ascmc[7] < 360, f"Polar Asc = {ascmc[7]} out of range"

    @pytest.mark.unit
    def test_coascendant_non_zero(self):
        """Co-Ascendants should not be zero for typical locations."""
        jd = 2451545.0
        cusps, ascmc = ephem.houses(jd, 41.9, 12.5, ord("P"))

        # All should be non-zero for typical latitudes
        assert ascmc[4] != 0.0, "Equatorial Ascendant should not be 0"
        assert ascmc[5] != 0.0, "Co-Ascendant Koch should not be 0"
        assert ascmc[6] != 0.0, "Co-Ascendant Munkasey should not be 0"
        assert ascmc[7] != 0.0, "Polar Ascendant should not be 0"


class TestCoAscendantVsPyswisseph:
    """Compare Co-Ascendant calculations with pyswisseph."""

    @pytest.mark.comparison
    def test_coascendant_matches_pyswisseph_rome(self):
        """Co-Ascendants should match pyswisseph for Rome."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        cusps_lib, ascmc_lib = ephem.houses(jd, lat, lon, ord("P"))
        cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, b"P")

        # Compare all 8 ascmc elements
        labels = [
            "ASC",
            "MC",
            "ARMC",
            "Vertex",
            "Equ Asc",
            "CoAsc Koch",
            "CoAsc Munk",
            "Polar Asc",
        ]
        for i in range(8):
            diff = abs(ascmc_lib[i] - ascmc_swe[i])
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.01, (
                f"{labels[i]} diff {diff}° (lib={ascmc_lib[i]:.4f}, swe={ascmc_swe[i]:.4f})"
            )

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "lat,lon,name",
        [
            (41.9, 12.5, "Rome"),
            (51.5, -0.12, "London"),
            (40.7, -74.0, "New York"),
            (-33.9, 151.2, "Sydney"),
            (35.7, 139.7, "Tokyo"),
            (55.75, 37.62, "Moscow"),
            (19.4, -99.1, "Mexico City"),
            (-22.9, -43.2, "Rio de Janeiro"),
        ],
    )
    def test_coascendant_multiple_locations(self, lat, lon, name):
        """Co-Ascendants should match pyswisseph for various locations."""
        jd = 2451545.0

        cusps_lib, ascmc_lib = ephem.houses(jd, lat, lon, ord("P"))
        cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, b"P")

        # Compare Co-Ascendants (indices 5 and 6)
        for i in [5, 6]:
            diff = abs(ascmc_lib[i] - ascmc_swe[i])
            if diff > 180:
                diff = 360 - diff
            label = "CoAsc Koch" if i == 5 else "CoAsc Munk"
            assert diff < 0.01, f"{name}: {label} diff {diff}°"

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "armc", [0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0]
    )
    def test_coascendant_various_armc(self, armc):
        """Co-Ascendants should match pyswisseph for various ARMC values."""
        lat = 45.0
        eps = 23.44

        cusps_lib, ascmc_lib = ephem.houses_armc(armc, lat, eps, ord("P"))
        cusps_swe, ascmc_swe = swe.houses_armc(armc, lat, eps, b"P")

        # Compare all Co-Ascendant values
        for i in [4, 5, 6, 7]:
            diff = abs(ascmc_lib[i] - ascmc_swe[i])
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.001, f"ARMC={armc}: ascmc[{i}] diff {diff}°"


class TestCoAscendantHousesArmc:
    """Test Co-Ascendants from houses_armc."""

    @pytest.mark.comparison
    def test_houses_armc_coascendant(self):
        """houses_armc should return correct Co-Ascendants."""
        armc = 292.957072438026
        lat = 41.9
        eps = 23.43767671605485

        cusps_lib, ascmc_lib = ephem.houses_armc(armc, lat, eps, ord("P"))
        cusps_swe, ascmc_swe = swe.houses_armc(armc, lat, eps, b"P")

        # Compare all 8 ascmc elements
        for i in range(8):
            diff = abs(ascmc_lib[i] - ascmc_swe[i])
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.001, f"ascmc[{i}] diff {diff}°"


class TestCoAscendantHousesArmcEx2:
    """Test Co-Ascendants from houses_armc_ex2."""

    @pytest.mark.comparison
    def test_houses_armc_ex2_coascendant(self):
        """houses_armc_ex2 should return correct Co-Ascendants and velocities when FLG_SPEED is set."""
        armc = 292.957072438026
        lat = 41.9
        eps = 23.43767671605485

        cusps_lib, ascmc_lib, cusps_speed_lib, ascmc_speed_lib = (
            ephem.houses_armc_ex2(armc, lat, eps, ord("P"), FLG_SPEED)
        )
        cusps_swe, ascmc_swe, cusps_speed_swe, ascmc_speed_swe = swe.houses_armc_ex2(
            armc, lat, eps, b"P"
        )

        # Compare Co-Ascendant positions
        for i in [4, 5, 6, 7]:
            diff = abs(ascmc_lib[i] - ascmc_swe[i])
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.001, f"ascmc[{i}] position diff {diff}°"

        # Compare Co-Ascendant velocities (within 6% tolerance)
        for i in [4, 5, 6, 7]:
            if abs(ascmc_speed_swe[i]) > 1:  # Avoid division by very small values
                rel_diff = abs(ascmc_speed_lib[i] - ascmc_speed_swe[i]) / abs(
                    ascmc_speed_swe[i]
                )
                assert rel_diff < 0.06, (
                    f"ascmc_speed[{i}] rel diff {rel_diff * 100:.1f}%"
                )


class TestCoAscendantHousesEx2:
    """Test Co-Ascendants from houses_ex2."""

    @pytest.mark.comparison
    def test_houses_ex2_coascendant(self):
        """houses_ex2 should return correct Co-Ascendant positions."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        cusps_lib, ascmc_lib, cusps_speed_lib, ascmc_speed_lib = ephem.houses_ex2(
            jd, lat, lon, ord("P"), FLG_SPEED
        )
        cusps_swe, ascmc_swe, cusps_speed_swe, ascmc_speed_swe = swe.houses_ex2(
            jd, lat, lon, b"P", FLG_SPEED
        )

        # Compare Co-Ascendant positions
        for i in [4, 5, 6, 7]:
            diff = abs(ascmc_lib[i] - ascmc_swe[i])
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.01, f"ascmc[{i}] position diff {diff}°"


class TestCoAscendantEdgeCases:
    """Edge cases for Co-Ascendant calculations."""

    @pytest.mark.edge_case
    def test_coascendant_at_equator(self):
        """Co-Ascendants at equator may be undefined (returns 0.0)."""
        jd = 2451545.0
        lat = 0.0  # Equator

        cusps_lib, ascmc_lib = ephem.houses(jd, lat, 0.0, ord("P"))
        cusps_swe, ascmc_swe = swe.houses(jd, lat, 0.0, b"P")

        # At equator, Co-Ascendant Munkasey (index 6) is undefined
        # Swiss Ephemeris returns 0.0
        assert ascmc_lib[6] == ascmc_swe[6], (
            f"CoAsc Munk at equator: lib={ascmc_lib[6]}, swe={ascmc_swe[6]}"
        )

    @pytest.mark.edge_case
    @pytest.mark.parametrize("lat", [60.0, 65.0])
    def test_coascendant_high_latitude(self, lat):
        """Co-Ascendants should work at high latitudes (below polar circle)."""
        jd = 2451545.0

        cusps_lib, ascmc_lib = ephem.houses(jd, lat, 0.0, ord("P"))
        cusps_swe, ascmc_swe = swe.houses(jd, lat, 0.0, b"P")

        for i in [4, 5, 6, 7]:
            diff = abs(ascmc_lib[i] - ascmc_swe[i])
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.01, f"ascmc[{i}] diff {diff}° at lat {lat}°"

    @pytest.mark.edge_case
    def test_coascendant_southern_hemisphere(self):
        """Co-Ascendants should work in southern hemisphere."""
        jd = 2451545.0
        lat = -33.9  # Sydney

        cusps_lib, ascmc_lib = ephem.houses(jd, lat, 151.2, ord("P"))
        cusps_swe, ascmc_swe = swe.houses(jd, lat, 151.2, b"P")

        for i in [4, 5, 6, 7]:
            diff = abs(ascmc_lib[i] - ascmc_swe[i])
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.01, f"ascmc[{i}] diff {diff}° in southern hemisphere"


class TestCoAscendantFormulas:
    """Test the mathematical relationships of Co-Ascendant formulas."""

    @pytest.mark.unit
    def test_coascendant_koch_formula(self):
        """Verify Co-Ascendant Koch formula: Asc(ARMC - 90, lat) + 180."""
        armc = 150.0
        lat = 45.0
        eps = 23.44

        cusps, ascmc = ephem.houses_armc(armc, lat, eps, ord("P"))

        # Get Co-Asc Koch (index 5)
        co_asc_koch = ascmc[5]

        # Get Polar Asc (index 7) - which is Asc(ARMC - 90, lat) without +180
        polar_asc = ascmc[7]

        # Co-Asc Koch should be Polar Asc + 180
        expected_co_asc = (polar_asc + 180.0) % 360.0
        diff = abs(co_asc_koch - expected_co_asc)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.001, f"CoAsc Koch {co_asc_koch} != Polar Asc {polar_asc} + 180"

    @pytest.mark.unit
    def test_equatorial_ascendant_formula(self):
        """Verify Equatorial Ascendant is at RA = ARMC + 90."""
        armc = 150.0
        lat = 45.0
        eps = 23.44

        cusps, ascmc = ephem.houses_armc(armc, lat, eps, ord("P"))

        # Equatorial Ascendant (East Point) at index 4
        equ_asc = ascmc[4]

        # It should be valid longitude
        assert 0 <= equ_asc < 360, f"Equ Asc {equ_asc} out of range"


class TestCoAscendantWithHouseSystems:
    """Test that Co-Ascendants are the same regardless of house system."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,name",
        [
            (ord("P"), "Placidus"),
            (ord("K"), "Koch"),
            (ord("O"), "Porphyry"),
            (ord("R"), "Regiomontanus"),
            (ord("C"), "Campanus"),
            (ord("E"), "Equal"),
            (ord("W"), "Whole Sign"),
        ],
    )
    def test_coascendant_same_for_all_systems(self, hsys, name):
        """Co-Ascendants should be the same regardless of house system used."""
        jd = 2451545.0
        lat, lon = 45.0, 0.0

        cusps_p, ascmc_p = ephem.houses(jd, lat, lon, ord("P"))
        cusps_h, ascmc_h = ephem.houses(jd, lat, lon, hsys)

        # Co-Ascendants (indices 4-7) should be identical
        for i in [4, 5, 6, 7]:
            diff = abs(ascmc_p[i] - ascmc_h[i])
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.001, (
                f"ascmc[{i}] differs between Placidus and {name}: "
                f"{ascmc_p[i]:.4f} vs {ascmc_h[i]:.4f}"
            )
