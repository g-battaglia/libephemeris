"""
Coordinate Transformation Comparison Tests.

Compares coordinate transformation functions between pyswisseph and libephemeris:
- cotrans / cotrans_sp - ecliptic/equatorial coordinate transforms
- azalt / azalt_rev - azimuth/altitude calculations
- refrac / refrac_extended - atmospheric refraction
"""

import pytest
import swisseph as swe
import libephemeris as ephem


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TOLERANCES
# ============================================================================

ANGLE_TOL = 0.0001  # degrees (0.36 arcsec)
REFRACTION_TOL = 0.004  # degrees (~14 arcsec, low-altitude refraction precision limit)


# ============================================================================
# TEST DATA
# ============================================================================

TEST_COORDINATES = [
    (0.0, 0.0, "Vernal Equinox"),
    (90.0, 0.0, "Summer Solstice"),
    (180.0, 0.0, "Autumnal Equinox"),
    (270.0, 0.0, "Winter Solstice"),
    (45.0, 23.5, "Mid ecliptic with lat"),
    (120.0, -15.0, "Negative latitude"),
    (315.0, 45.0, "High latitude"),
]

OBLIQUITY_VALUES = [
    (23.4393, "Modern obliquity"),
    (23.5, "Round value"),
    (23.0, "Lower obliquity"),
    (24.0, "Higher obliquity"),
]

AZALT_LOCATIONS = [
    ("Equator", 0.0, 0.0, 0),
    ("Mid-lat North", 45.0, 0.0, 0),
    ("High-lat North", 70.0, 0.0, 0),
    ("Mid-lat South", -35.0, 0.0, 0),
]

ALTITUDE_VALUES = [0.0, 5.0, 10.0, 20.0, 45.0, 70.0, 85.0, 90.0]


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestCotrans:
    """Compare cotrans coordinate transformation."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("lon,lat,desc", TEST_COORDINATES)
    @pytest.mark.parametrize("eps,eps_desc", OBLIQUITY_VALUES)
    def test_cotrans(self, lon, lat, desc, eps, eps_desc):
        """Test cotrans ecliptic to equatorial transformation."""
        coords = (lon, lat, 1.0)

        result_swe = swe.cotrans(coords, -eps)
        result_py = ephem.cotrans(coords, -eps)

        diff_lon = angular_diff(result_swe[0], result_py[0])
        diff_lat = abs(result_swe[1] - result_py[1])
        max_diff = max(diff_lon, diff_lat)

        assert max_diff < ANGLE_TOL, (
            f"cotrans {desc} with {eps_desc}: max diff {max_diff:.8f}°"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("lon,lat,desc", TEST_COORDINATES[:4])
    def test_cotrans_roundtrip(self, lon, lat, desc):
        """Test cotrans ecliptic -> equatorial -> ecliptic roundtrip."""
        eps = 23.4393
        coords = (lon, lat, 1.0)

        # Ecliptic to equatorial
        equatorial = ephem.cotrans(coords, -eps)

        # Equatorial back to ecliptic
        ecliptic = ephem.cotrans((equatorial[0], equatorial[1], 1.0), eps)

        diff_lon = angular_diff(lon, ecliptic[0])
        diff_lat = abs(lat - ecliptic[1])

        assert diff_lon < 1e-8, f"Roundtrip lon diff {diff_lon:.12f}° for {desc}"
        assert diff_lat < 1e-8, f"Roundtrip lat diff {diff_lat:.12f}° for {desc}"


VELOCITY_VALUES = [
    ((1.0, 0.0, 0.0), "Pure lon speed"),
    ((0.0, 0.5, 0.0), "Pure lat speed"),
    ((1.0, 0.1, 0.0), "Typical planet speed"),
    ((0.5, -0.2, 0.0), "Mixed with negative lat speed"),
    ((13.2, 0.0, 0.0), "High speed (Moon-like)"),
    ((-0.5, 0.05, 0.0), "Retrograde motion"),
]

SPEED_TOL = 0.001  # degrees/day tolerance for velocity


class TestCotransSp:
    """Compare cotrans_sp with velocity transformation."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("lon,lat,desc", TEST_COORDINATES)
    @pytest.mark.parametrize("eps,eps_desc", OBLIQUITY_VALUES)
    @pytest.mark.parametrize("speed,speed_desc", VELOCITY_VALUES)
    def test_cotrans_sp_ecl_to_eq(
        self, lon, lat, desc, eps, eps_desc, speed, speed_desc
    ):
        """Test cotrans_sp ecliptic to equatorial with various speeds."""
        coord = (lon, lat, 1.0)

        # libephemeris uses separate coord and speed tuples
        result_coord, result_speed = ephem.cotrans_sp(coord, speed, -eps)

        # pyswisseph uses a 6-tuple (lon, lat, dist, lon_sp, lat_sp, dist_sp)
        coords_6 = coord + speed
        result_swe = swe.cotrans_sp(coords_6, -eps)

        # Compare positions
        diff_lon = angular_diff(result_swe[0], result_coord[0])
        diff_lat = abs(result_swe[1] - result_coord[1])

        assert diff_lon < ANGLE_TOL, (
            f"cotrans_sp lon diff {diff_lon:.8f}° for {desc}, {eps_desc}, {speed_desc}"
        )
        assert diff_lat < ANGLE_TOL, (
            f"cotrans_sp lat diff {diff_lat:.8f}° for {desc}, {eps_desc}, {speed_desc}"
        )

        # Compare velocities
        diff_lon_speed = abs(result_swe[3] - result_speed[0])
        diff_lat_speed = abs(result_swe[4] - result_speed[1])

        assert diff_lon_speed < SPEED_TOL, (
            f"cotrans_sp lon_speed diff {diff_lon_speed:.8f} for {desc}, {speed_desc}"
        )
        assert diff_lat_speed < SPEED_TOL, (
            f"cotrans_sp lat_speed diff {diff_lat_speed:.8f} for {desc}, {speed_desc}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("lon,lat,desc", TEST_COORDINATES)
    @pytest.mark.parametrize("speed,speed_desc", VELOCITY_VALUES[:3])
    def test_cotrans_sp_eq_to_ecl(self, lon, lat, desc, speed, speed_desc):
        """Test cotrans_sp equatorial to ecliptic with various speeds."""
        eps = 23.4393
        coord = (lon, lat, 1.0)

        # Positive obliquity for equatorial to ecliptic
        result_coord, result_speed = ephem.cotrans_sp(coord, speed, eps)

        coords_6 = coord + speed
        result_swe = swe.cotrans_sp(coords_6, eps)

        diff_lon = angular_diff(result_swe[0], result_coord[0])
        diff_lat = abs(result_swe[1] - result_coord[1])

        assert diff_lon < ANGLE_TOL, (
            f"cotrans_sp eq->ecl lon diff {diff_lon:.8f}° for {desc}"
        )
        assert diff_lat < ANGLE_TOL, (
            f"cotrans_sp eq->ecl lat diff {diff_lat:.8f}° for {desc}"
        )

        diff_lon_speed = abs(result_swe[3] - result_speed[0])
        diff_lat_speed = abs(result_swe[4] - result_speed[1])

        assert diff_lon_speed < SPEED_TOL, (
            f"cotrans_sp eq->ecl lon_speed diff {diff_lon_speed:.8f} for {desc}"
        )
        assert diff_lat_speed < SPEED_TOL, (
            f"cotrans_sp eq->ecl lat_speed diff {diff_lat_speed:.8f} for {desc}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("lon,lat,desc", TEST_COORDINATES)
    @pytest.mark.parametrize("speed,speed_desc", VELOCITY_VALUES[:4])
    def test_cotrans_sp_roundtrip_velocity(self, lon, lat, desc, speed, speed_desc):
        """Test cotrans_sp roundtrip preserves velocity coherence.

        Velocity should be preserved when transforming:
        ecliptic -> equatorial -> ecliptic
        """
        eps = 23.4393
        coord = (lon, lat, 1.0)

        # Ecliptic to equatorial
        eq_coord, eq_speed = ephem.cotrans_sp(coord, speed, -eps)

        # Equatorial back to ecliptic
        result_coord, result_speed = ephem.cotrans_sp(eq_coord, eq_speed, eps)

        # Position roundtrip
        diff_lon = angular_diff(lon, result_coord[0])
        diff_lat = abs(lat - result_coord[1])

        assert diff_lon < 1e-8, (
            f"Roundtrip lon diff {diff_lon:.12f}° for {desc}, {speed_desc}"
        )
        assert diff_lat < 1e-8, (
            f"Roundtrip lat diff {diff_lat:.12f}° for {desc}, {speed_desc}"
        )

        # Velocity roundtrip - velocities should be recovered
        diff_lon_speed = abs(speed[0] - result_speed[0])
        diff_lat_speed = abs(speed[1] - result_speed[1])
        diff_dist_speed = abs(speed[2] - result_speed[2])

        assert diff_lon_speed < 1e-8, (
            f"Roundtrip lon_speed diff {diff_lon_speed:.12f} for {desc}, {speed_desc}"
        )
        assert diff_lat_speed < 1e-8, (
            f"Roundtrip lat_speed diff {diff_lat_speed:.12f} for {desc}, {speed_desc}"
        )
        assert diff_dist_speed < 1e-10, (
            f"Roundtrip dist_speed diff {diff_dist_speed:.12f} for {desc}, {speed_desc}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("lon,lat,desc", TEST_COORDINATES[:4])
    def test_cotrans_sp_position_matches_cotrans(self, lon, lat, desc):
        """Test that cotrans_sp position matches cotrans (velocity-free)."""
        eps = 23.4393
        coord = (lon, lat, 1.0)
        speed = (1.0, 0.1, 0.0)

        # Get position from cotrans_sp
        result_coord, _ = ephem.cotrans_sp(coord, speed, -eps)

        # Get position from cotrans
        result_cotrans = ephem.cotrans(coord, -eps)

        diff_lon = angular_diff(result_coord[0], result_cotrans[0])
        diff_lat = abs(result_coord[1] - result_cotrans[1])

        assert diff_lon < 1e-10, (
            f"cotrans_sp/cotrans lon mismatch {diff_lon:.12f}° for {desc}"
        )
        assert diff_lat < 1e-10, (
            f"cotrans_sp/cotrans lat mismatch {diff_lat:.12f}° for {desc}"
        )

    @pytest.mark.comparison
    def test_cotrans_sp_zero_velocity(self):
        """Test cotrans_sp with zero velocity returns zero velocity."""
        eps = 23.4393
        coord = (90.0, 23.0, 1.0)
        speed = (0.0, 0.0, 0.0)

        result_coord, result_speed = ephem.cotrans_sp(coord, speed, -eps)

        # pyswisseph reference
        coords_6 = coord + speed
        result_swe = swe.cotrans_sp(coords_6, -eps)

        # Transformed speed should be zero (or very close)
        assert abs(result_speed[0]) < 1e-10, (
            f"Zero input should give zero lon_speed, got {result_speed[0]}"
        )
        assert abs(result_speed[1]) < 1e-10, (
            f"Zero input should give zero lat_speed, got {result_speed[1]}"
        )
        assert abs(result_speed[2]) < 1e-10, (
            f"Zero input should give zero dist_speed, got {result_speed[2]}"
        )

        # Should match pyswisseph
        assert abs(result_swe[3] - result_speed[0]) < 1e-10
        assert abs(result_swe[4] - result_speed[1]) < 1e-10

    @pytest.mark.comparison
    def test_cotrans_sp_at_poles(self):
        """Test cotrans_sp at ecliptic poles with velocity."""
        eps = 23.4393
        speed = (1.0, 0.1, 0.0)

        # North ecliptic pole
        coord_n = (0.0, 89.0, 1.0)  # Near pole to avoid singularity
        result_n, speed_n = ephem.cotrans_sp(coord_n, speed, -eps)
        swe_n = swe.cotrans_sp(coord_n + speed, -eps)

        diff_lat = abs(result_n[1] - swe_n[1])
        assert diff_lat < ANGLE_TOL, f"Pole lat diff {diff_lat:.10f}°"

        diff_lat_speed = abs(speed_n[1] - swe_n[4])
        assert diff_lat_speed < SPEED_TOL, f"Pole lat_speed diff {diff_lat_speed:.10f}"

        # South ecliptic pole
        coord_s = (0.0, -89.0, 1.0)
        result_s, speed_s = ephem.cotrans_sp(coord_s, speed, -eps)
        swe_s = swe.cotrans_sp(coord_s + speed, -eps)

        diff_lat = abs(result_s[1] - swe_s[1])
        assert diff_lat < ANGLE_TOL, f"S Pole lat diff {diff_lat:.10f}°"


class TestAzalt:
    """Compare azalt azimuth/altitude calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("name,lat,lon,alt", AZALT_LOCATIONS)
    def test_azalt(self, name, lat, lon, alt):
        """Test azalt calculations."""
        jd = swe.julday(2000, 1, 1, 12.0)
        geopos = (lon, lat, alt)
        atpress = 1013.25  # Standard pressure
        attemp = 15.0  # Standard temperature

        # Sun's equatorial coordinates at J2000
        xin = (280.0, -23.0, 1.0)  # Approximate RA, Dec, distance

        result_swe = swe.azalt(jd, swe.ECL2HOR, geopos, atpress, attemp, xin)
        result_py = ephem.azalt(jd, 0, geopos, atpress, attemp, xin)  # ECL2HOR = 0

        diff_az = angular_diff(result_swe[0], result_py[0])
        diff_alt = abs(result_swe[1] - result_py[1])

        assert diff_az < 0.01, f"azalt azimuth diff {diff_az:.6f}° at {name}"
        assert diff_alt < 0.01, f"azalt altitude diff {diff_alt:.6f}° at {name}"


class TestAzaltRev:
    """Compare azalt_rev reverse transformation."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("name,lat,lon,alt", AZALT_LOCATIONS)
    def test_azalt_rev(self, name, lat, lon, alt):
        """Test azalt_rev reverse calculations."""
        jd = swe.julday(2000, 1, 1, 12.0)
        geopos = (lon, lat, alt)

        # Test with azimuth/altitude values (separate args, not tuple)
        azimuth = 180.0
        altitude = 45.0

        result_swe = swe.azalt_rev(jd, swe.HOR2ECL, geopos, azimuth, altitude)
        result_py = ephem.azalt_rev(jd, ephem.HOR2ECL, geopos, azimuth, altitude)

        diff_lon = angular_diff(result_swe[0], result_py[0])
        diff_lat = abs(result_swe[1] - result_py[1])

        assert diff_lon < 0.01, f"azalt_rev lon diff {diff_lon:.6f}° at {name}"
        assert diff_lat < 0.01, f"azalt_rev lat diff {diff_lat:.6f}° at {name}"


class TestRefrac:
    """Compare atmospheric refraction calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("alt_true", ALTITUDE_VALUES)
    def test_refrac(self, alt_true):
        """Test refrac atmospheric refraction."""
        pressure = 1013.25
        temperature = 15.0

        refrac_swe = swe.refrac(alt_true, pressure, temperature, swe.TRUE_TO_APP)
        refrac_py = ephem.refrac(
            alt_true, pressure, temperature, 0
        )  # TRUE_TO_APP = 0

        diff = abs(refrac_swe - refrac_py)

        assert diff < REFRACTION_TOL, f"refrac diff {diff:.6f}° at altitude {alt_true}°"

    @pytest.mark.comparison
    @pytest.mark.parametrize("alt_app", [0.5, 5.0, 10.0, 30.0, 60.0, 90.0])
    def test_refrac_reverse(self, alt_app):
        """Test refrac apparent to true direction."""
        pressure = 1013.25
        temperature = 15.0

        refrac_swe = swe.refrac(alt_app, pressure, temperature, swe.APP_TO_TRUE)
        refrac_py = ephem.refrac(
            alt_app, pressure, temperature, 1
        )  # APP_TO_TRUE = 1

        diff = abs(refrac_swe - refrac_py)

        assert diff < REFRACTION_TOL, (
            f"refrac reverse diff {diff:.6f}° at altitude {alt_app}°"
        )


class TestRefracExtended:
    """Compare extended refraction calculations."""

    @pytest.mark.comparison
    def test_refrac_extended(self):
        """Test refrac_extended with various parameters."""
        alt_true = 10.0
        geoalt = 0.0  # Sea level
        pressure = 1013.25
        temperature = 15.0
        lapse_rate = 0.0065  # Standard lapse rate
        calc_flag = 0  # TRUE_TO_APP

        result_swe = swe.refrac_extended(
            alt_true, geoalt, pressure, temperature, lapse_rate, calc_flag
        )
        result_py = ephem.refrac_extended(
            alt_true, geoalt, pressure, temperature, lapse_rate, calc_flag
        )

        # refrac_extended returns tuple (refracted_alt, dret)
        diff = abs(result_swe[0] - result_py[0])

        assert diff < REFRACTION_TOL, f"refrac_extended diff {diff:.6f}°"


class TestSpecialCases:
    """Test special cases and edge conditions."""

    @pytest.mark.comparison
    def test_cotrans_at_poles(self):
        """Test cotrans at ecliptic poles."""
        eps = 23.4393

        # North ecliptic pole
        coords_n = (0.0, 90.0, 1.0)
        result_swe = swe.cotrans(coords_n, -eps)
        result_py = ephem.cotrans(coords_n, -eps)

        diff_lat = abs(result_swe[1] - result_py[1])
        assert diff_lat < ANGLE_TOL, f"Cotrans at N pole: lat diff {diff_lat:.10f}°"

        # South ecliptic pole
        coords_s = (0.0, -90.0, 1.0)
        result_swe = swe.cotrans(coords_s, -eps)
        result_py = ephem.cotrans(coords_s, -eps)

        diff_lat = abs(result_swe[1] - result_py[1])
        assert diff_lat < ANGLE_TOL, f"Cotrans at S pole: lat diff {diff_lat:.10f}°"

    @pytest.mark.comparison
    def test_refrac_at_horizon(self):
        """Test refraction at horizon (maximum effect)."""
        pressure = 1013.25
        temperature = 15.0

        # At horizon, refraction is maximum (~0.57°)
        refrac_swe = swe.refrac(0.0, pressure, temperature, swe.TRUE_TO_APP)
        refrac_py = ephem.refrac(0.0, pressure, temperature, 0)

        diff = abs(refrac_swe - refrac_py)

        assert diff < REFRACTION_TOL, f"Horizon refraction diff {diff:.6f}°"

        # Both should give significant refraction at horizon (typically ~0.5°)
        assert abs(refrac_swe) > 0.4, "Expected significant refraction at horizon"
        assert abs(refrac_py) > 0.4, "Expected significant refraction at horizon"
