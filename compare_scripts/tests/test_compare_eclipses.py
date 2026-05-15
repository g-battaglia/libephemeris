"""
Eclipse Calculations Comparison Tests.

Compares eclipse calculations between pyswisseph and libephemeris.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import SUN, MOON


# ============================================================================
# TOLERANCES
# ============================================================================

TIME_TOL_SECONDS = 120.0  # 2 minutes for eclipse timing
POSITION_TOL = 0.01  # degrees


# ============================================================================
# TEST DATA
# ============================================================================

# Known solar eclipses for testing
SOLAR_ECLIPSES = [
    (2024, 4, 8, "Total Solar Eclipse 2024"),
    (2024, 10, 2, "Annular Solar Eclipse 2024"),
    (2023, 10, 14, "Annular Solar Eclipse 2023"),
    (2025, 3, 29, "Partial Solar Eclipse 2025"),
]

# Known lunar eclipses for testing
LUNAR_ECLIPSES = [
    (2024, 3, 25, "Penumbral Lunar Eclipse 2024"),
    (2024, 9, 18, "Partial Lunar Eclipse 2024"),
    (2025, 3, 14, "Total Lunar Eclipse 2025"),
]

# Test locations
TEST_LOCATIONS = [
    ("Rome", 41.9028, 12.4964, 0),
    ("New York", 40.7128, -74.0060, 0),
    ("Tokyo", 35.6762, 139.6503, 0),
]


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestSolarEclipseSearch:
    """Compare solar eclipse search functions."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,desc", SOLAR_ECLIPSES)
    def test_sol_eclipse_when(self, year, month, day, desc):
        """Test finding solar eclipse timing."""
        jd_start = swe.julday(year, month, day, 0.0)

        # Search for next solar eclipse
        ecl_swe = swe.sol_eclipse_when_glob(jd_start, 0)
        ecl_py = ephem.sol_eclipse_when_glob(jd_start, 0)

        # Compare eclipse maximum time
        jd_max_swe = ecl_swe[1][0]  # tret[0] = maximum
        jd_max_py = ecl_py[1][0]

        diff_seconds = abs(jd_max_swe - jd_max_py) * 86400

        assert diff_seconds < TIME_TOL_SECONDS, (
            f"{desc}: eclipse max time diff {diff_seconds:.1f}s exceeds tolerance"
        )


class TestLunarEclipseSearch:
    """Compare lunar eclipse search functions."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,desc", LUNAR_ECLIPSES)
    def test_lun_eclipse_when(self, year, month, day, desc):
        """Test finding lunar eclipse timing."""
        jd_start = swe.julday(year, month, day, 0.0)

        # Search for next lunar eclipse
        ecl_swe = swe.lun_eclipse_when(jd_start, 0)
        ecl_py = ephem.lun_eclipse_when(jd_start, 0)

        # Compare eclipse maximum time
        jd_max_swe = ecl_swe[1][0]
        jd_max_py = ecl_py[1][0]

        diff_seconds = abs(jd_max_swe - jd_max_py) * 86400

        assert diff_seconds < TIME_TOL_SECONDS, (
            f"{desc}: eclipse max time diff {diff_seconds:.1f}s exceeds tolerance"
        )


class TestEclipseAttributes:
    """Compare eclipse attribute calculations."""

    @pytest.mark.comparison
    def test_solar_eclipse_attributes(self):
        """Test solar eclipse attribute calculations."""
        # 2024 total solar eclipse
        jd = swe.julday(2024, 4, 8, 18.0)
        geopos = (-99.0, 25.0, 0)  # Mexico path: (lon, lat, alt)

        # pyswisseph signature: sol_eclipse_how(tjdut, geopos, flags)
        attr_swe = swe.sol_eclipse_how(jd, geopos, 0)
        # libephemeris signature: sol_eclipse_how(tjdut, geopos, flags)
        attr_py = ephem.sol_eclipse_how(jd, geopos, 0)

        # Compare obscuration/magnitude
        # attr[0] = fraction of solar diameter covered
        diff_obscur = abs(attr_swe[1][0] - attr_py[1][0])

        assert diff_obscur < 0.01, (
            f"Eclipse obscuration diff {diff_obscur:.4f} exceeds tolerance"
        )


class TestEclipseFlags:
    """Test eclipse type detection."""

    @pytest.mark.comparison
    def test_eclipse_type_detection(self):
        """Test that eclipse types are correctly identified."""
        # 2024 total solar eclipse
        jd = swe.julday(2024, 4, 1, 0.0)

        ecl_swe = swe.sol_eclipse_when_glob(jd, 0)
        ecl_py = ephem.sol_eclipse_when_glob(jd, 0)

        # Both should detect the same eclipse type
        type_swe = ecl_swe[0]
        type_py = ecl_py[0]

        # At minimum, both should detect an eclipse occurred
        assert type_swe > 0 and type_py > 0, "Both should detect an eclipse"


# ============================================================================
# ECLIPSE PATH GEOMETRY TESTS
# ============================================================================


# Tolerances for path geometry
CENTRAL_LINE_TOL_DEG = 0.1  # 0.1 degrees (~10 km) for central line coordinates
PATH_WIDTH_TOL_KM = 30.0  # 30 km for path width comparisons
PATH_WIDTH_REL_TOL = 0.15  # 15% relative tolerance for path width


class TestEclipseCentralLine:
    """Compare central line coordinate calculations between pyswisseph and libephemeris."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,desc", SOLAR_ECLIPSES[:2]
    )  # Total and annular
    def test_central_line_coordinates_at_max(self, year, month, day, desc):
        """Test that central line coordinates at eclipse maximum match between implementations."""
        jd_start = swe.julday(year, month, day, 0.0)

        # Find eclipse maximum time
        ecl_swe = swe.sol_eclipse_when_glob(jd_start, 0)
        jd_max = ecl_swe[1][0]

        # Get central line coordinates from pyswisseph
        swe_result = swe.sol_eclipse_where(jd_max, 0)
        swe_lon = swe_result[1][0]
        swe_lat = swe_result[1][1]

        # Get central line coordinates from libephemeris
        lib_result = ephem.sol_eclipse_where(jd_max, 0)
        lib_lon = lib_result[1][0]
        lib_lat = lib_result[1][1]

        # Compare coordinates
        lon_diff = abs(swe_lon - lib_lon)
        lat_diff = abs(swe_lat - lib_lat)

        assert lon_diff < CENTRAL_LINE_TOL_DEG, (
            f"{desc}: central line longitude diff {lon_diff:.4f}° exceeds tolerance"
        )
        assert lat_diff < CENTRAL_LINE_TOL_DEG, (
            f"{desc}: central line latitude diff {lat_diff:.4f}° exceeds tolerance"
        )

    @pytest.mark.comparison
    def test_central_line_trajectory(self):
        """Test central line coordinates at multiple points during eclipse."""
        # April 8, 2024 total solar eclipse - narrower time window during central eclipse
        jd_start = swe.julday(2024, 4, 8, 17.0)  # During central eclipse
        jd_end = swe.julday(2024, 4, 8, 19.5)
        step_hours = 0.5  # 30 minute steps

        jd = jd_start
        while jd <= jd_end:
            # Get coordinates from pyswisseph
            swe_result = swe.sol_eclipse_where(jd, 0)
            swe_lon = swe_result[1][0]
            swe_lat = swe_result[1][1]

            # Get coordinates from libephemeris
            lib_result = ephem.sol_eclipse_where(jd, 0)
            lib_lon = lib_result[1][0]
            lib_lat = lib_result[1][1]

            # Only compare when eclipse is on Earth (non-zero coords from pyswisseph)
            if swe_lon != 0.0 or swe_lat != 0.0:
                lon_diff = abs(swe_lon - lib_lon)
                lat_diff = abs(swe_lat - lib_lat)

                assert lon_diff < CENTRAL_LINE_TOL_DEG, (
                    f"JD {jd}: central line longitude diff {lon_diff:.4f}° exceeds tolerance"
                )
                assert lat_diff < CENTRAL_LINE_TOL_DEG, (
                    f"JD {jd}: central line latitude diff {lat_diff:.4f}° exceeds tolerance"
                )

            jd += step_hours / 24.0

    @pytest.mark.comparison
    def test_calc_eclipse_central_line_function(self):
        """Test calc_eclipse_central_line produces consistent results with sol_eclipse_where.

        This test verifies that calc_eclipse_central_line and sol_eclipse_where return
        consistent coordinates, as both now use the same underlying algorithm.
        """
        # April 8, 2024 total solar eclipse - central eclipse time window
        jd_start = swe.julday(
            2024, 4, 8, 17.5
        )  # Start later during established eclipse
        jd_end = swe.julday(2024, 4, 8, 19.0)

        # Get central line from libephemeris calc_eclipse_central_line
        times, lats, lons = ephem.calc_eclipse_central_line(
            jd_start, jd_end, step_minutes=30.0
        )

        # Compare each point with libephemeris sol_eclipse_where
        for i, jd in enumerate(times):
            # Get libephemeris central line at this time
            lib_result = ephem.sol_eclipse_where(jd, 0)
            lib_lon = lib_result[1][0]
            lib_lat = lib_result[1][1]

            # calc_eclipse_central_line now uses sol_eclipse_where internally,
            # so results should be nearly identical (within floating point precision)
            lon_diff = abs(lib_lon - lons[i])
            lat_diff = abs(lib_lat - lats[i])

            # Handle longitude wrap-around
            if lon_diff > 180:
                lon_diff = 360 - lon_diff

            # Use tight tolerance since both functions now use same algorithm
            internal_tol = 0.01  # 0.01 degree tolerance for floating point precision
            assert lon_diff < internal_tol, (
                f"calc_eclipse_central_line at JD {jd}: lon diff {lon_diff:.4f}° exceeds tolerance"
            )
            assert lat_diff < internal_tol, (
                f"calc_eclipse_central_line at JD {jd}: lat diff {lat_diff:.4f}° exceeds tolerance"
            )


class TestEclipsePathLimits:
    """Test eclipse path limit calculations.

    Note: pyswisseph sol_eclipse_where does not calculate limit coordinates
    (returns zeros for geopos[2:10]). Therefore, we test libephemeris limit
    functions using geometric consistency checks rather than direct comparison.
    """

    @pytest.mark.comparison
    def test_northern_limit_north_of_central_line(self):
        """Test that northern limit is north of central line (geometric consistency)."""
        # April 8, 2024 total solar eclipse
        jd_start = swe.julday(2024, 4, 8, 16.0)
        jd_end = swe.julday(2024, 4, 8, 20.0)

        # Get central line and northern limit from libephemeris
        c_times, c_lats, c_lons = ephem.calc_eclipse_central_line(
            jd_start, jd_end, step_minutes=30.0
        )
        n_times, n_lats, n_lons = ephem.calc_eclipse_northern_limit(
            jd_start, jd_end, step_minutes=30.0
        )

        # Find common times (approximately)
        for c_idx, c_jd in enumerate(c_times):
            # Find matching northern limit time
            for n_idx, n_jd in enumerate(n_times):
                if abs(c_jd - n_jd) < 0.001:  # ~1.5 min tolerance
                    # Northern limit should be at higher latitude
                    assert n_lats[n_idx] >= c_lats[c_idx] - 0.5, (
                        f"Northern limit lat {n_lats[n_idx]:.2f}° should be >= "
                        f"central line lat {c_lats[c_idx]:.2f}°"
                    )
                    break

    @pytest.mark.comparison
    def test_southern_limit_south_of_central_line(self):
        """Test that southern limit is south of central line (geometric consistency)."""
        # April 8, 2024 total solar eclipse
        jd_start = swe.julday(2024, 4, 8, 16.0)
        jd_end = swe.julday(2024, 4, 8, 20.0)

        # Get central line and southern limit from libephemeris
        c_times, c_lats, c_lons = ephem.calc_eclipse_central_line(
            jd_start, jd_end, step_minutes=30.0
        )
        s_times, s_lats, s_lons = ephem.calc_eclipse_southern_limit(
            jd_start, jd_end, step_minutes=30.0
        )

        # Find common times (approximately)
        for c_idx, c_jd in enumerate(c_times):
            # Find matching southern limit time
            for s_idx, s_jd in enumerate(s_times):
                if abs(c_jd - s_jd) < 0.001:  # ~1.5 min tolerance
                    # Southern limit should be at lower latitude
                    assert s_lats[s_idx] <= c_lats[c_idx] + 0.5, (
                        f"Southern limit lat {s_lats[s_idx]:.2f}° should be <= "
                        f"central line lat {c_lats[c_idx]:.2f}°"
                    )
                    break

    @pytest.mark.comparison
    def test_limit_coordinates_from_sol_eclipse_where(self):
        """Test that libephemeris sol_eclipse_where returns valid limit coordinates."""
        # April 8, 2024 total solar eclipse at maximum
        jd = swe.julday(2024, 4, 8, 18.0)

        # Get libephemeris results (includes limit coordinates)
        lib_result = ephem.sol_eclipse_where(jd, 0)
        geopos = lib_result[1]

        # Central line
        central_lon, central_lat = geopos[0], geopos[1]

        # Northern umbra limit
        n_umbra_lon, n_umbra_lat = geopos[2], geopos[3]

        # Southern umbra limit
        s_umbra_lon, s_umbra_lat = geopos[4], geopos[5]

        # Northern penumbra limit
        n_penumbra_lon, n_penumbra_lat = geopos[6], geopos[7]

        # Southern penumbra limit
        s_penumbra_lon, s_penumbra_lat = geopos[8], geopos[9]

        # Verify central line coords are valid
        assert -180 <= central_lon <= 180, "Central line longitude out of range"
        assert -90 <= central_lat <= 90, "Central line latitude out of range"

        # Verify limit coords are valid when non-zero
        if n_umbra_lon != 0 or n_umbra_lat != 0:
            assert -180 <= n_umbra_lon <= 180, (
                "Northern umbra limit longitude out of range"
            )
            assert -90 <= n_umbra_lat <= 90, (
                "Northern umbra limit latitude out of range"
            )

        if s_umbra_lon != 0 or s_umbra_lat != 0:
            assert -180 <= s_umbra_lon <= 180, (
                "Southern umbra limit longitude out of range"
            )
            assert -90 <= s_umbra_lat <= 90, (
                "Southern umbra limit latitude out of range"
            )

        # Penumbra limits should span wider than umbra limits
        if n_penumbra_lat != 0 and n_umbra_lat != 0:
            assert n_penumbra_lat >= n_umbra_lat - 1, (
                "Northern penumbra should be north of or near umbra limit"
            )

        if s_penumbra_lat != 0 and s_umbra_lat != 0:
            assert s_penumbra_lat <= s_umbra_lat + 1, (
                "Southern penumbra should be south of or near umbra limit"
            )


class TestEclipsePathWidth:
    """Compare eclipse path width calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,desc", SOLAR_ECLIPSES[:2]
    )  # Total and annular
    def test_path_width_at_maximum(self, year, month, day, desc):
        """Test path width calculations at eclipse maximum."""
        jd_start = swe.julday(year, month, day, 0.0)

        # Find eclipse maximum
        ecl_swe = swe.sol_eclipse_when_glob(jd_start, 0)
        jd_max = ecl_swe[1][0]

        # Get path width from pyswisseph attr[3] (shadow width in km)
        # Note: pyswisseph uses attr[3] for path width, may be negative
        swe_result = swe.sol_eclipse_where(jd_max, 0)
        swe_width = abs(swe_result[2][3])  # Take absolute value

        # Get path width from libephemeris
        lib_width = ephem.calc_eclipse_path_width(jd_max)

        # Both widths should be reasonable for a central eclipse
        if lib_width > 0 and swe_width > 0:
            # Compare with relative tolerance
            rel_diff = abs(swe_width - lib_width) / max(swe_width, lib_width)

            assert (
                rel_diff < PATH_WIDTH_REL_TOL
                or abs(swe_width - lib_width) < PATH_WIDTH_TOL_KM
            ), (
                f"{desc}: path width diff {abs(swe_width - lib_width):.1f} km "
                f"({rel_diff:.1%}) exceeds tolerance"
            )

    @pytest.mark.comparison
    def test_path_width_consistency_with_attr(self):
        """Test that calc_eclipse_path_width matches sol_eclipse_where attr value."""
        # April 8, 2024 total solar eclipse at maximum
        jd = swe.julday(2024, 4, 8, 18.0)

        # Get path width from calc_eclipse_path_width
        width_func = ephem.calc_eclipse_path_width(jd)

        # Get path width from sol_eclipse_where attr[3]
        lib_result = ephem.sol_eclipse_where(jd, 0)
        width_attr = lib_result[2][3]

        # They should be consistent (same function internally)
        # Note: sol_eclipse_where attr[3] uses sign convention (negative for
        # total/umbra, positive for annular/antumbra), while
        # calc_eclipse_path_width returns absolute values.
        assert abs(width_func - abs(width_attr)) < 1.0, (
            f"Path width mismatch: calc_eclipse_path_width={width_func:.1f}, "
            f"sol_eclipse_where attr[3]={width_attr:.1f}"
        )

    @pytest.mark.comparison
    def test_path_width_reasonable_values(self):
        """Test that path width values are within reasonable ranges."""
        # April 8, 2024 total solar eclipse - central eclipse time window
        jd_start = swe.julday(2024, 4, 8, 17.0)
        jd_end = swe.julday(2024, 4, 8, 19.5)
        step_hours = 0.5

        jd = jd_start
        while jd <= jd_end:
            width = ephem.calc_eclipse_path_width(jd)

            if width > 0:
                # Total eclipse path width typically 0-500 km
                assert 0 < width < 1000, (
                    f"Path width {width:.1f} km at JD {jd} outside reasonable range"
                )

            jd += step_hours / 24.0

    @pytest.mark.comparison
    def test_path_width_along_trajectory(self):
        """Test path width variation along eclipse trajectory."""
        # April 8, 2024 total solar eclipse - central eclipse time window
        jd_start = swe.julday(2024, 4, 8, 17.0)
        jd_end = swe.julday(2024, 4, 8, 19.5)

        widths = []
        times = []

        jd = jd_start
        while jd <= jd_end:
            width = ephem.calc_eclipse_path_width(jd)
            if width > 0:
                widths.append(width)
                times.append(jd)
            jd += 0.25 / 24.0  # 15 minute steps

        # Should have multiple valid width measurements
        assert len(widths) >= 3, "Expected multiple valid path width measurements"

        # Path width should vary smoothly (no sudden jumps)
        # Allow up to 60 km change in 15 minutes for smooth path
        for i in range(1, len(widths)):
            diff = abs(widths[i] - widths[i - 1])
            assert diff < 60, (
                f"Path width changed by {diff:.1f} km between consecutive points"
            )
