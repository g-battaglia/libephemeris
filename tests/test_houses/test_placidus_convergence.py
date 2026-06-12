"""
Tests for Placidus house convergence accuracy.

These tests verify that the Placidus house system iteration converges
to match pyswisseph within 0.0001° (0.36 arcseconds).
"""

import pytest

swe = pytest.importorskip("swisseph")
import libephemeris as ephem


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# High precision tolerance: 0.0001° = 0.36 arcseconds
HIGH_PRECISION_TOL = 0.0001


class TestPlacidusConvergence:
    """Tests for Placidus house cusp convergence to pyswisseph precision."""

    @pytest.mark.parametrize(
        "name,jd,lat,lon",
        [
            # john_lennon (1940-10-09 07:00:00, Liverpool)
            ("john_lennon", 2428867.7916666665, 53.4084, -2.9916),
            # paul_mccartney (1942-06-18 14:30:00, Liverpool)
            ("paul_mccartney", 2431259.4375, 53.4084, -2.9916),
            # johnny_depp (1963-06-09 03:45:00, San Francisco)
            ("johnny_depp", 2439606.15625, 37.7749, -122.4194),
        ],
    )
    def test_placidus_high_precision_match(self, name, jd, lat, lon):
        """Test that Placidus cusps match pyswisseph within 0.0001°."""
        cusps_swe, _ = swe.houses(jd, lat, lon, b"P")
        cusps_py, _ = ephem.houses(jd, lat, lon, ord("P"))

        max_diff = 0.0
        worst_cusp = 0
        for i in range(12):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            if diff > max_diff:
                max_diff = diff
                worst_cusp = i + 1

        assert max_diff < HIGH_PRECISION_TOL, (
            f"{name}: Placidus cusp {worst_cusp} differs by {max_diff:.6f}° "
            f"(target: <{HIGH_PRECISION_TOL}°)"
        )

    @pytest.mark.parametrize(
        "name,lat,lon",
        [
            ("Rome", 41.9028, 12.4964),
            ("New York", 40.7128, -74.0060),
            ("Sydney", -33.8688, 151.2093),
            ("Equator", 0.0, 0.0),
            ("London", 51.5074, -0.1278),
            ("Tokyo", 35.6762, 139.6503),
        ],
    )
    @pytest.mark.parametrize(
        "year,month,day,hour",
        [
            (2000, 1, 1, 12.0),  # J2000
            (2024, 6, 15, 0.0),  # Summer 2024
            (1980, 5, 20, 14.5),  # Past
            (1940, 10, 9, 7.0),  # john_lennon era
        ],
    )
    def test_placidus_convergence_across_dates(
        self, name, lat, lon, year, month, day, hour
    ):
        """Test Placidus convergence across various dates and locations."""
        jd = swe.julday(year, month, day, hour)

        cusps_swe, _ = swe.houses(jd, lat, lon, b"P")
        cusps_py, _ = ephem.houses(jd, lat, lon, ord("P"))

        max_diff = 0.0
        for i in range(12):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            max_diff = max(max_diff, diff)

        assert max_diff < HIGH_PRECISION_TOL, (
            f"Placidus at {name} ({year}-{month:02d}-{day:02d}): "
            f"max diff {max_diff:.6f}° exceeds {HIGH_PRECISION_TOL}°"
        )


class TestPlacidusIntermediateCusps:
    """
    Test intermediate cusps (2, 3, 5, 6, 8, 9, 11, 12) specifically.

    These cusps require iterative calculation and are most affected by
    the convergence algorithm's precision.
    """

    @pytest.mark.parametrize(
        "name,jd,lat,lon",
        [
            ("john_lennon", 2428867.7916666665, 53.4084, -2.9916),
            ("paul_mccartney", 2431259.4375, 53.4084, -2.9916),
        ],
    )
    def test_intermediate_cusps_precision(self, name, jd, lat, lon):
        """Test that intermediate cusps (calculated iteratively) are precise."""
        cusps_swe, _ = swe.houses(jd, lat, lon, b"P")
        cusps_py, _ = ephem.houses(jd, lat, lon, ord("P"))

        # Intermediate cusps (not Asc/IC/Dsc/MC)
        intermediate_cusps = [2, 3, 5, 6, 8, 9, 11, 12]

        for cusp_num in intermediate_cusps:
            i = cusp_num - 1  # 0-indexed
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            assert diff < HIGH_PRECISION_TOL, (
                f"{name}: Cusp {cusp_num} diff {diff:.6f}° exceeds {HIGH_PRECISION_TOL}°"
            )


# Tolerance for polar-adjacent latitudes (slightly relaxed due to convergence complexity)
POLAR_ADJACENT_TOL = 0.01  # 0.01° = 36 arcseconds


class TestPolarAdjacentLatitudes:
    """
    Tests for Placidus and Koch at polar-adjacent latitudes (64-66°).

    These latitudes are challenging because they are close to the polar threshold
    (~66.56°) where Placidus/Koch become undefined. The iteration may take longer
    to converge or may oscillate.

    The tests verify that cusps match pyswisseph within 0.01° (36 arcseconds).
    """

    @pytest.mark.parametrize(
        "name,jd,lat,lon",
        [
            # Reykjavik, Iceland (64.1466°N)
            ("reykjavik_64n", 2451545.0, 64.1466, -21.9426),
            # Additional Reykjavik times to test stability
            ("reykjavik_64n_summer", 2451728.0, 64.1466, -21.9426),
            ("reykjavik_64n_winter", 2451545.0 + 180, 64.1466, -21.9426),
            # Arctic circle edge cases (66°N)
            ("arctic_circle_66n", 2451545.0, 66.0, 0.0),
            ("arctic_circle_66n_summer", 2451728.0, 66.0, 0.0),
            # Southern polar-adjacent
            ("antarctic_circle_66s", 2451545.0, -66.0, 0.0),
            ("antarctic_circle_66s_summer", 2451728.0, -66.0, 0.0),
            # High northern latitudes - Scandinavia
            ("tromso_69_4n", 2451545.0, 65.5, 19.0),  # Near Tromsø
            ("fairbanks_64_8n", 2451545.0, 64.8, -147.7),  # Fairbanks, Alaska
        ],
    )
    def test_placidus_polar_adjacent_match(self, name, jd, lat, lon):
        """Test that Placidus cusps at polar-adjacent latitudes match pyswisseph."""
        cusps_swe, _ = swe.houses(jd, lat, lon, b"P")
        cusps_py, _ = ephem.houses(jd, lat, lon, ord("P"))

        max_diff = 0.0
        worst_cusp = 0
        for i in range(12):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            if diff > max_diff:
                max_diff = diff
                worst_cusp = i + 1

        assert max_diff < POLAR_ADJACENT_TOL, (
            f"{name}: Placidus cusp {worst_cusp} differs by {max_diff:.6f}° "
            f"(target: <{POLAR_ADJACENT_TOL}°)"
        )

    @pytest.mark.parametrize(
        "name,jd,lat,lon",
        [
            # Reykjavik, Iceland (64.1466°N)
            ("reykjavik_64n", 2451545.0, 64.1466, -21.9426),
            # Arctic circle edge cases (66°N)
            ("arctic_circle_66n", 2451545.0, 66.0, 0.0),
            # Southern polar-adjacent
            ("antarctic_circle_66s", 2451545.0, -66.0, 0.0),
        ],
    )
    def test_koch_polar_adjacent_match(self, name, jd, lat, lon):
        """Test that Koch cusps at polar-adjacent latitudes match pyswisseph."""
        cusps_swe, _ = swe.houses(jd, lat, lon, b"K")
        cusps_py, _ = ephem.houses(jd, lat, lon, ord("K"))

        max_diff = 0.0
        worst_cusp = 0
        for i in range(12):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            if diff > max_diff:
                max_diff = diff
                worst_cusp = i + 1

        assert max_diff < POLAR_ADJACENT_TOL, (
            f"{name}: Koch cusp {worst_cusp} differs by {max_diff:.6f}° "
            f"(target: <{POLAR_ADJACENT_TOL}°)"
        )


class TestPlacidusConvergenceStability:
    """
    Test Placidus iteration stability near the polar threshold.

    At latitudes near 66° (the polar threshold), the iteration may become
    slow to converge or oscillate. These tests ensure the algorithm detects
    non-convergence and falls back gracefully.
    """

    @pytest.mark.parametrize(
        "lat",
        [64.0, 64.5, 65.0, 65.5, 66.0, 66.2, 66.4],
    )
    def test_placidus_stable_at_various_polar_adjacent_latitudes(self, lat):
        """Test Placidus stability at various latitudes approaching polar threshold."""
        jd = 2451545.0  # J2000

        cusps_swe, _ = swe.houses(jd, lat, 0.0, b"P")
        cusps_py, _ = ephem.houses(jd, lat, 0.0, ord("P"))

        max_diff = 0.0
        for i in range(12):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            max_diff = max(max_diff, diff)

        # At these latitudes, we should match pyswisseph within tolerance
        assert max_diff < POLAR_ADJACENT_TOL, (
            f"Placidus at lat={lat}°: max diff {max_diff:.6f}° exceeds tolerance"
        )

    @pytest.mark.parametrize(
        "jd_offset",
        [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330],  # Monthly through year
    )
    def test_placidus_seasonal_stability_at_66n(self, jd_offset):
        """Test Placidus stability throughout the year at 66°N."""
        jd = 2451545.0 + jd_offset  # J2000 + offset days
        lat = 66.0

        cusps_swe, _ = swe.houses(jd, lat, 0.0, b"P")
        cusps_py, _ = ephem.houses(jd, lat, 0.0, ord("P"))

        max_diff = 0.0
        for i in range(12):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            max_diff = max(max_diff, diff)

        assert max_diff < POLAR_ADJACENT_TOL, (
            f"Placidus at 66°N, jd_offset={jd_offset}: "
            f"max diff {max_diff:.6f}° exceeds tolerance"
        )
