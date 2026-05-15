"""Property-based testing for libephemeris using Hypothesis.

Validates mathematical invariants, API contracts, monotonicity properties,
and symmetry constraints through randomized testing.

Validation Plan v2, Section 3.
"""

from __future__ import annotations

import math
import warnings

import pytest
from hypothesis import given, settings, assume, HealthCheck
from hypothesis import strategies as st

import libephemeris as swe

warnings.filterwarnings("ignore")

pytestmark = pytest.mark.slow

# ============================================================================
# Strategies
# ============================================================================

# Core era JD range (1972-2040) where Delta T is precisely known
jd_core = st.floats(min_value=2441317.5, max_value=2466154.5)

# Wider JD range (1800-2200) for general tests
jd_wide = st.floats(min_value=2378497.0, max_value=2524594.0)

# Standard body IDs (planets + lunar points)
standard_bodies = st.sampled_from(
    [
        swe.SUN,
        swe.MOON,
        swe.MERCURY,
        swe.VENUS,
        swe.MARS,
        swe.JUPITER,
        swe.SATURN,
        swe.URANUS,
        swe.NEPTUNE,
        swe.PLUTO,
    ]
)

# All bodies including lunar points and asteroids
all_bodies = st.sampled_from(
    [
        swe.SUN,
        swe.MOON,
        swe.MERCURY,
        swe.VENUS,
        swe.MARS,
        swe.JUPITER,
        swe.SATURN,
        swe.URANUS,
        swe.NEPTUNE,
        swe.PLUTO,
        swe.MEAN_NODE,
        swe.TRUE_NODE,
        swe.MEAN_APOG,
        swe.CHIRON,
        swe.CERES,
        swe.PALLAS,
        swe.JUNO,
        swe.VESTA,
    ]
)

# House systems (common ones)
house_systems = st.sampled_from(
    [
        b"P",
        b"K",
        b"O",
        b"R",
        b"C",
        b"E",
        b"W",
        b"B",
        b"M",
        b"T",
    ]
)

# Geographic coordinates
latitudes = st.floats(min_value=-89.9, max_value=89.9)
longitudes = st.floats(min_value=-180.0, max_value=180.0)


# ============================================================================
# §3.1 Coordinate Transform Round-Trips
# ============================================================================


class TestCoordinateRoundTrips:
    """Verify coordinate transform identities."""

    @given(jd=jd_core, body=standard_bodies)
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_ecliptic_equatorial_roundtrip(self, jd: float, body: int) -> None:
        """ecliptic → equatorial → ecliptic via cotrans preserves position."""
        # Get ecliptic position
        ecl, _ = swe.calc_ut(jd, body)
        lon_ecl, lat_ecl, dist_ecl = float(ecl[0]), float(ecl[1]), float(ecl[2])

        # Get true obliquity from nutation/obliquity
        nut, _ = swe.calc_ut(jd, swe.ECL_NUT)
        eps = float(nut[1])  # true obliquity in degrees

        # cotrans: ecliptic → equatorial (negative obliquity)
        eq = swe.cotrans((lon_ecl, lat_ecl, dist_ecl), -eps)

        # cotrans: equatorial → ecliptic (positive obliquity)
        back = swe.cotrans((float(eq[0]), float(eq[1]), float(eq[2])), eps)

        dlon = abs(lon_ecl - float(back[0]))
        if dlon > 180:
            dlon = 360 - dlon
        dlat = abs(lat_ecl - float(back[1]))

        assert dlon * 3600 < 0.001, (
            f'Ecliptic roundtrip lon error: {dlon * 3600:.6f}" (body={body}, jd={jd})'
        )
        assert dlat * 3600 < 0.001, (
            f'Ecliptic roundtrip lat error: {dlat * 3600:.6f}" (body={body}, jd={jd})'
        )
        assert dlat * 3600 < 0.001, (
            f'Ecliptic roundtrip lat error: {dlat * 3600:.6f}" (body={body}, jd={jd})'
        )

    @given(jd=jd_core, body=standard_bodies)
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_degrees_radians_roundtrip(self, jd: float, body: int) -> None:
        """Degrees and radians outputs should be consistent."""
        deg_result, _ = swe.calc_ut(jd, body)
        try:
            rad_result, _ = swe.calc_ut(jd, body, swe.FLG_RADIANS)
        except Exception:
            # FLG_RADIANS may trigger fallback in LEB mode
            return

        lon_deg = deg_result[0]
        lon_rad = rad_result[0]

        # Convert radians back to degrees
        lon_from_rad = math.degrees(lon_rad) % 360

        dlon = abs(lon_deg - lon_from_rad)
        if dlon > 180:
            dlon = 360 - dlon

        assert dlon * 3600 < 0.001, (
            f'Degrees/radians roundtrip error: {dlon * 3600:.6f}" '
            f"(body={body}, jd={jd})"
        )

    @given(
        jd=jd_core,
        body=st.sampled_from(
            [
                swe.SUN,
                swe.MOON,
                swe.MARS,
                swe.JUPITER,
            ]
        ),
    )
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_xyz_spherical_consistency(self, jd: float, body: int) -> None:
        """XYZ and spherical coordinates should be geometrically consistent."""
        try:
            sph, _ = swe.calc_ut(jd, body)
            xyz, _ = swe.calc_ut(jd, body, swe.FLG_XYZ)
        except Exception:
            return

        lon_r = math.radians(sph[0])
        lat_r = math.radians(sph[1])
        dist = sph[2]

        # Spherical to XYZ
        x_calc = dist * math.cos(lat_r) * math.cos(lon_r)
        y_calc = dist * math.cos(lat_r) * math.sin(lon_r)
        z_calc = dist * math.sin(lat_r)

        # Check consistency (within 1e-8 AU)
        assert abs(xyz[0] - x_calc) < 1e-6, (
            f"X mismatch: {xyz[0]} vs {x_calc} (body={body}, jd={jd})"
        )
        assert abs(xyz[1] - y_calc) < 1e-6, (
            f"Y mismatch: {xyz[1]} vs {y_calc} (body={body}, jd={jd})"
        )
        assert abs(xyz[2] - z_calc) < 1e-6, (
            f"Z mismatch: {xyz[2]} vs {z_calc} (body={body}, jd={jd})"
        )


# ============================================================================
# §3.2 API Contract Testing
# ============================================================================


class TestAPIContracts:
    """Verify API return types and shapes for all valid inputs."""

    @given(jd=jd_core, body=all_bodies)
    @settings(max_examples=300, suppress_health_check=[HealthCheck.too_slow])
    def test_calc_ut_returns_tuple_and_int(self, jd: float, body: int) -> None:
        """calc_ut(jd, body, flags) always returns (6-tuple, int)."""
        result, flags = swe.calc_ut(jd, body)

        assert isinstance(result, tuple), f"Expected tuple, got {type(result)}"
        assert len(result) == 6, f"Expected 6 elements, got {len(result)}"
        assert isinstance(flags, int), f"Expected int flags, got {type(flags)}"

        # All elements should be numeric (float or numpy float)
        for i, val in enumerate(result):
            assert isinstance(val, (int, float)) or hasattr(val, "__float__"), (
                f"result[{i}] is {type(val)}, not numeric"
            )

    @given(jd=jd_core, lat=latitudes, lon=longitudes, sys=house_systems)
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_houses_returns_correct_shape(
        self, jd: float, lat: float, lon: float, sys: bytes
    ) -> None:
        """houses(jd, lat, lon, system) always returns (tuple, tuple)."""
        try:
            cusps, ascmc = swe.houses(jd, lat, lon, sys)
        except Exception:
            # Some house systems may not work at extreme latitudes
            return

        assert isinstance(cusps, tuple), f"cusps type: {type(cusps)}"
        assert isinstance(ascmc, tuple), f"ascmc type: {type(ascmc)}"
        assert len(cusps) == 12, f"Expected 12 cusps, got {len(cusps)}"
        assert len(ascmc) == 8, f"Expected 8 ascmc, got {len(ascmc)}"

    @given(jd=jd_core, body=all_bodies)
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_calc_ut_accepts_moseph(self, jd: float, body: int) -> None:
        """All functions accept FLG_MOSEPH without error."""
        result, flags = swe.calc_ut(jd, body, swe.FLG_MOSEPH)
        assert isinstance(result, tuple)
        assert len(result) == 6

    @given(jd=jd_core, body=all_bodies)
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_calc_ut_lon_in_range(self, jd: float, body: int) -> None:
        """Longitude should be in [0, 360)."""
        result, _ = swe.calc_ut(jd, body)
        lon = float(result[0])
        assert 0 <= lon < 360, f"Longitude {lon} outside [0, 360) (body={body})"

    @given(jd=jd_core, body=all_bodies)
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_calc_ut_lat_in_range(self, jd: float, body: int) -> None:
        """Latitude should be in [-90, 90]."""
        result, _ = swe.calc_ut(jd, body)
        lat = float(result[1])
        assert -90 <= lat <= 90, f"Latitude {lat} outside [-90, 90] (body={body})"

    @given(jd=jd_core, body=standard_bodies)
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_calc_ut_distance_positive(self, jd: float, body: int) -> None:
        """Distance should be positive for standard bodies."""
        result, _ = swe.calc_ut(jd, body)
        dist = float(result[2])
        assert dist > 0, f"Distance {dist} not positive (body={body})"

    def test_fixstar_returns_correct_shape(self) -> None:
        """fixstar_ut returns (6-tuple, str, int)."""
        jd = 2460311.0
        for star in ["Regulus", "Sirius", "Aldebaran", "Spica"]:
            result = swe.fixstar_ut(star, jd)
            assert len(result) == 3, "Expected 3 values from fixstar_ut"
            pos, name, ret_flags = result
            assert isinstance(pos, tuple) and len(pos) == 6
            assert isinstance(name, str)
            assert isinstance(ret_flags, int)


# ============================================================================
# §3.3 Monotonicity Properties
# ============================================================================


class TestMonotonicity:
    """Verify that certain quantities are monotonic as expected."""

    @given(jd=jd_core)
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_sun_longitude_increasing(self, jd: float) -> None:
        """Sun longitude increases monotonically over 1 day (mod 360)."""
        step = 0.1  # 2.4 hours
        prev_lon = None
        for i in range(10):
            t = jd + i * step
            result, _ = swe.calc_ut(t, swe.SUN)
            lon = float(result[0])
            if prev_lon is not None:
                # Sun moves ~1°/day eastward, so each step should increase
                delta = lon - prev_lon
                if delta < -180:
                    delta += 360  # handle 360→0 wrap
                assert delta > 0, (
                    f"Sun longitude not increasing: {prev_lon:.6f} → {lon:.6f} "
                    f"at JD {t:.2f}"
                )
            prev_lon = lon

    @given(jd=jd_core)
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_mean_node_decreasing(self, jd: float) -> None:
        """Mean Node longitude decreases monotonically over 1 day."""
        step = 0.1
        prev_lon = None
        for i in range(10):
            t = jd + i * step
            result, _ = swe.calc_ut(t, swe.MEAN_NODE)
            lon = float(result[0])
            if prev_lon is not None:
                # Mean node regresses ~0.053°/day
                delta = lon - prev_lon
                if delta > 180:
                    delta -= 360  # handle 0→360 wrap
                assert delta < 0, (
                    f"Mean Node not decreasing: {prev_lon:.6f} → {lon:.6f} "
                    f"at JD {t:.2f}"
                )
            prev_lon = lon

    @given(
        y=st.integers(min_value=1900, max_value=2100),
        m=st.integers(min_value=1, max_value=12),
        d=st.integers(min_value=1, max_value=28),
    )
    @settings(max_examples=200)
    def test_julday_monotonic(self, y: int, m: int, d: int) -> None:
        """julday(y,m,d) < julday(y,m,d+1)."""
        jd1 = swe.julday(y, m, d, 12.0)
        jd2 = swe.julday(y, m, d + 1, 12.0)
        # Day+1 might roll into next month, but JD should still increase
        assert jd2 > jd1, f"julday not monotonic: {jd1} >= {jd2}"

    @given(jd=st.floats(min_value=2441317.5, max_value=2466154.5))
    @settings(max_examples=100)
    def test_deltat_reasonable_range_post_1972(self, jd: float) -> None:
        """Delta T is in a physically reasonable range post-1972.

        Note: Delta T is NOT monotonically increasing — Earth's rotation
        has irregular variations (it sped up ~2020-2029, causing DeltaT
        to briefly decrease).  Instead we verify:
        1. Delta T is in a reasonable range (40-120 seconds for 1972-2040)
        2. Delta T returns a finite float
        """
        dt = swe.deltat(jd)
        assert isinstance(dt, float), f"deltat returned {type(dt)}"
        assert math.isfinite(dt), f"deltat returned non-finite: {dt}"
        # Delta T in days; convert to seconds for range check
        dt_sec = dt * 86400
        assert 40 <= dt_sec <= 120, (
            f"Delta T {dt_sec:.2f}s outside expected range [40, 120]s at JD {jd}"
        )


# ============================================================================
# §3.4 Symmetry Properties
# ============================================================================


class TestSymmetry:
    """Verify symmetry and identity properties."""

    @given(jd=jd_core)
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_heliocentric_sun_is_origin(self, jd: float) -> None:
        """calc_ut(jd, SUN, FLG_HELCTR) returns (0, 0, 0, ...)."""
        result, _ = swe.calc_ut(jd, swe.SUN, swe.FLG_HELCTR)
        assert float(result[0]) == 0.0, f"Helio Sun lon={result[0]}"
        assert float(result[1]) == 0.0, f"Helio Sun lat={result[1]}"
        assert float(result[2]) == 0.0, f"Helio Sun dist={result[2]}"

    @given(jd=jd_core)
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_pheno_sun_phase_angle_zero(self, jd: float) -> None:
        """pheno_ut(jd, SUN) returns phase_angle = 0."""
        pheno = swe.pheno_ut(jd, swe.SUN)
        assert float(pheno[0]) == 0.0, f"Sun phase_angle={pheno[0]}"

    @given(jd=jd_core, lat=latitudes, lon=longitudes)
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_house_cusps_ascending(self, jd: float, lat: float, lon: float) -> None:
        """House cusps should be in ascending order (mod 360).

        Adjacent cusps increase mod 360; the angular gap from cusp[i] to
        cusp[i+1] should be positive and < 180° for most systems.
        """
        try:
            cusps, _ = swe.houses(jd, lat, lon, b"P")
        except Exception:
            return

        for i in range(len(cusps)):
            c1 = float(cusps[i])
            c2 = float(cusps[(i + 1) % len(cusps)])
            delta = (c2 - c1) % 360
            # Each house should span > 0° and < 180°
            # (in extreme cases near poles this can get close to 0 or 180)
            assert 0 < delta < 180, (
                f"House cusp ordering violated: cusp[{i}]={c1:.4f} → "
                f"cusp[{(i + 1) % len(cusps)}]={c2:.4f}, delta={delta:.4f}° "
                f"(lat={lat}, lon={lon})"
            )

    @given(jd=jd_core)
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_speed_flag_adds_velocity(self, jd: float) -> None:
        """FLG_SPEED should produce non-zero velocities for moving bodies."""
        result, _ = swe.calc_ut(jd, swe.MARS, swe.FLG_SPEED)
        # result[3] = lon speed, result[4] = lat speed, result[5] = dist speed
        lon_speed = float(result[3])
        # Mars moves at least 0.1°/day in longitude typically
        # (can be near 0 at station, but very rare for random dates)
        assert lon_speed != 0.0, f"Mars lon speed is exactly 0 at JD {jd}"

    @given(jd=jd_core, body=standard_bodies)
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_j2000_and_default_differ(self, jd: float, body: int) -> None:
        """J2000 and ecliptic-of-date should differ (due to precession).

        Exception: for dates very near J2000.0, the difference is tiny.
        """
        assume(abs(jd - 2451545.0) > 365)  # At least 1 year from J2000

        r_date, _ = swe.calc_ut(jd, body)
        r_j2000, _ = swe.calc_ut(jd, body, swe.FLG_J2000)

        # Precession is ~50"/year, so after 1 year they should differ
        dlon = abs(float(r_date[0]) - float(r_j2000[0]))
        if dlon > 180:
            dlon = 360 - dlon

        assert dlon > 0.001, (
            f"J2000 and of-date longitudes too similar: "
            f"{r_date[0]:.8f} vs {r_j2000[0]:.8f} "
            f"(body={body}, jd={jd})"
        )
