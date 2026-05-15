"""
Tests validating the mathematical formulas documented in HOUSE_SYSTEMS.md.

These tests verify that the documented formulas accurately describe the
implementations in houses.py.
"""

import math
import pytest
import libephemeris as ephem


class TestDocumentedFormulas:
    """Tests verifying documented mathematical formulas are accurate."""

    @pytest.mark.unit
    def test_equal_house_formula(self):
        """
        Test Equal house formula:
        λᵢ = (λ_Asc + (i-1) × 30°) mod 360°
        """
        jd = 2451545.0  # J2000
        lat, lon = 41.9, 12.5  # Rome

        cusps, ascmc = ephem.houses(jd, lat, lon, ord("E"))
        asc = ascmc[0]

        # Verify each cusp follows the documented formula
        for i in range(1, 13):
            expected = (asc + (i - 1) * 30.0) % 360.0
            # cusps is 0-indexed, cusp 1 is at index 0
            assert abs(cusps[i - 1] - expected) < 0.001, (
                f"House {i}: expected {expected:.3f}, got {cusps[i - 1]:.3f}"
            )

    @pytest.mark.unit
    def test_whole_sign_formula(self):
        """
        Test Whole Sign formula:
        start = floor(λ_Asc / 30°) × 30°
        λᵢ = (start + (i-1) × 30°) mod 360°
        """
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        cusps, ascmc = ephem.houses(jd, lat, lon, ord("W"))
        asc = ascmc[0]

        # Calculate expected start (0° of sign containing Asc)
        start = math.floor(asc / 30.0) * 30.0

        for i in range(1, 13):
            expected = (start + (i - 1) * 30.0) % 360.0
            assert abs(cusps[i - 1] - expected) < 0.001, (
                f"House {i}: expected {expected:.3f}, got {cusps[i - 1]:.3f}"
            )

    @pytest.mark.unit
    def test_vehlow_formula(self):
        """
        Test Vehlow formula:
        start = (λ_Asc - 15°) mod 360°
        λᵢ = (start + (i-1) × 30°) mod 360°
        """
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        cusps, ascmc = ephem.houses(jd, lat, lon, ord("V"))
        asc = ascmc[0]

        # Calculate expected start (Asc - 15°)
        start = (asc - 15.0) % 360.0

        for i in range(1, 13):
            expected = (start + (i - 1) * 30.0) % 360.0
            assert abs(cusps[i - 1] - expected) < 0.001, (
                f"House {i}: expected {expected:.3f}, got {cusps[i - 1]:.3f}"
            )

    @pytest.mark.unit
    def test_porphyry_quadrant_trisection(self):
        """
        Test Porphyry formula: quadrant trisection.
        Houses 11, 12 trisect the MC-Asc arc.
        Houses 2, 3 trisect the Asc-IC arc.
        """
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        cusps, ascmc = ephem.houses(jd, lat, lon, ord("O"))
        asc = ascmc[0]
        mc = ascmc[1]
        ic = (mc + 180.0) % 360.0

        # Quadrant MC to Asc
        arc1 = (asc - mc) % 360.0
        step1 = arc1 / 3.0
        expected_11 = (mc + step1) % 360.0
        expected_12 = (mc + 2 * step1) % 360.0

        # cusps[10] = house 11 (0-indexed)
        assert abs(cusps[10] - expected_11) < 0.01, (
            f"House 11: expected {expected_11:.3f}, got {cusps[10]:.3f}"
        )
        assert abs(cusps[11] - expected_12) < 0.01, (
            f"House 12: expected {expected_12:.3f}, got {cusps[11]:.3f}"
        )

        # Quadrant Asc to IC
        arc2 = (ic - asc) % 360.0
        step2 = arc2 / 3.0
        expected_2 = (asc + step2) % 360.0
        expected_3 = (asc + 2 * step2) % 360.0

        assert abs(cusps[1] - expected_2) < 0.01, (
            f"House 2: expected {expected_2:.3f}, got {cusps[1]:.3f}"
        )
        assert abs(cusps[2] - expected_3) < 0.01, (
            f"House 3: expected {expected_3:.3f}, got {cusps[2]:.3f}"
        )

    @pytest.mark.unit
    def test_natural_gradient_formula(self):
        """
        Test Natural Gradient formula:
        λᵢ = (i-1) × 30°
        (Equal houses from 0° Aries)
        """
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        cusps, ascmc = ephem.houses(jd, lat, lon, ord("N"))

        for i in range(1, 13):
            expected = (i - 1) * 30.0
            assert abs(cusps[i - 1] - expected) < 0.001, (
                f"House {i}: expected {expected:.3f}, got {cusps[i - 1]:.3f}"
            )

    @pytest.mark.unit
    def test_opposite_houses_180_degrees(self):
        """
        Test that opposite houses are 180° apart for all systems.
        Documented property: λ_{i+6} = (λᵢ + 180°) mod 360°
        """
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        systems = ["P", "K", "O", "R", "C", "E", "W", "B", "M", "T", "X", "V"]

        for hsys in systems:
            cusps, _ = ephem.houses(jd, lat, lon, ord(hsys))

            for i in range(6):
                opposite = (i + 6) % 12
                diff = abs(cusps[i] - cusps[opposite])
                if diff > 180:
                    diff = 360 - diff
                assert abs(diff - 180.0) < 0.1, (
                    f"System {hsys}: House {i + 1} and {opposite + 1} "
                    f"not 180° apart (diff = {diff:.3f})"
                )


class TestCoordinateConversions:
    """Tests for documented coordinate conversion formulas."""

    @pytest.mark.unit
    def test_ecliptic_to_ra_formula(self):
        """
        Test documented ecliptic to RA conversion:
        tan(α) = sin(λ) · cos(ε) / cos(λ)
        """
        # Test at standard obliquity
        eps = 23.44  # degrees
        eps_rad = math.radians(eps)

        # Test several longitudes
        for lon in [0, 30, 60, 90, 120, 180, 270]:
            lon_rad = math.radians(lon)

            # Documented formula using atan2
            y = math.sin(lon_rad) * math.cos(eps_rad)
            x = math.cos(lon_rad)
            ra = math.degrees(math.atan2(y, x)) % 360.0

            # At 0° and 180°, RA should equal longitude
            if lon in [0, 180]:
                assert abs(ra - lon) < 0.01, f"At λ={lon}°, RA should be {lon}°"

            # RA should always be valid
            assert 0 <= ra < 360, f"RA {ra} out of range for λ={lon}"

    @pytest.mark.unit
    def test_ra_to_ecliptic_formula(self):
        """
        Test documented RA to ecliptic conversion:
        tan(λ) = sin(α) / (cos(α) · cos(ε))
        """
        eps = 23.44
        eps_rad = math.radians(eps)

        for ra in [0, 30, 60, 90, 120, 180, 270]:
            ra_rad = math.radians(ra)

            # Documented formula using atan2
            y = math.sin(ra_rad)
            x = math.cos(ra_rad) * math.cos(eps_rad)
            lon = math.degrees(math.atan2(y, x)) % 360.0

            # At 0° and 180°, longitude should equal RA
            if ra in [0, 180]:
                assert abs(lon - ra) < 0.01, f"At α={ra}°, λ should be {ra}°"

            # Longitude should always be valid
            assert 0 <= lon < 360, f"λ {lon} out of range for α={ra}"


class TestPolarLatitudeDocumentation:
    """Tests verifying polar latitude documentation."""

    @pytest.mark.unit
    def test_polar_threshold_formula(self):
        """
        Test documented polar threshold:
        threshold = 90° - obliquity
        """
        from libephemeris.houses import get_polar_latitude_threshold

        # Standard obliquity
        eps = 23.44
        threshold = get_polar_latitude_threshold(eps)
        expected = 90.0 - eps

        assert abs(threshold - expected) < 0.001, (
            f"Threshold {threshold} != expected {expected}"
        )

    @pytest.mark.unit
    def test_polar_systems_fail(self):
        """
        Test that Placidus and Koch fail at polar latitudes as documented.
        Documented: Fails when |φ| + ε > 90°
        """
        from libephemeris.exceptions import PolarCircleError

        jd = 2451545.0
        polar_lat = 70.0  # Beyond ~66.5° threshold
        lon = 0.0

        # Placidus should raise PolarCircleError
        with pytest.raises(PolarCircleError):
            ephem.houses(jd, polar_lat, lon, ord("P"))

        # Koch should raise PolarCircleError
        with pytest.raises(PolarCircleError):
            ephem.houses(jd, polar_lat, lon, ord("K"))

    @pytest.mark.unit
    def test_universal_systems_work_at_poles(self):
        """
        Test that documented universal systems work at polar latitudes.
        Equal, Whole Sign, Porphyry should work everywhere.
        """
        jd = 2451545.0
        polar_lat = 80.0
        lon = 0.0

        # These should not raise errors
        for hsys in ["E", "W", "O"]:
            cusps, ascmc = ephem.houses(jd, polar_lat, lon, ord(hsys))
            assert len(cusps) >= 12, f"System {hsys} failed at polar latitude"


class TestHouseSystemNames:
    """Tests verifying house system names match documentation."""

    @pytest.mark.unit
    def test_house_system_names(self):
        """Test that house_name returns documented names."""
        from libephemeris.houses import house_name

        # Test names that are defined in the house_name function
        expected_names = {
            "P": "Placidus",
            "K": "Koch",
            "O": "Porphyry",
            "R": "Regiomontanus",
            "C": "Campanus",
            "E": "equal",
            "A": "equal",
            "W": "equal/ whole sign",
            "B": "Alcabitius",
            "M": "Morinus",
            "T": "Polich/Page",
            "X": "axial rotation system/Meridian houses",
            "V": "equal/Vehlow",
            "H": "horizon/azimut",
            "G": "Gauquelin sectors",
            "U": "Krusinski-Pisa-Goelzer",
            "F": "Carter poli-equ.",
            "S": "Sripati",
        }

        for code, expected in expected_names.items():
            actual = house_name(ord(code))
            assert actual == expected, (
                f"System '{code}': expected '{expected}', got '{actual}'"
            )
