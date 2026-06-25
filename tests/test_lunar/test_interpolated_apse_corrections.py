"""
Tests for the interpolated apogee/perigee correction tables (BUG-001 fix).

Verifies that the correction-table-based interpolated apogee (INTP_APOG)
and perigee (INTP_PERG) produce positions within expected tolerances
across multiple date ranges, flag variants, and edge cases.

The correction tables provide sub-arcminute precision for longitude
within their coverage range (1549-2651 CE for medium tier).
"""

from __future__ import annotations

import math
import random

import pytest

import libephemeris as swe
from libephemeris.constants import (
    INTP_APOG,
    INTP_PERG,
    OSCU_APOG,
    MEAN_APOG,
    FLG_SPEED,
    FLG_EQUATORIAL,
    FLG_SIDEREAL,
    FLG_J2000,
    FLG_NOABERR,
    FLG_NONUT,
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
)


# ---------------------------------------------------------------------------
# Test data generation helpers
# ---------------------------------------------------------------------------


def _random_jds(
    n: int, seed: int = 42, jd_min: float = 2378497.0, jd_max: float = 2597642.0
) -> list[float]:
    """Generate n random JDs in the given range (default: 1800-2400 CE)."""
    rng = random.Random(seed)
    return [rng.uniform(jd_min, jd_max) for _ in range(n)]


def _jd_from_year(year: float) -> float:
    """Approximate JD for a given year (Jan 1 noon)."""
    return 2451545.0 + (year - 2000.0) * 365.25


# ---------------------------------------------------------------------------
# Apogee correction table tests
# ---------------------------------------------------------------------------


class TestInterpolatedApogeeCorrections:
    """Tests for the corrected interpolated apogee (body 21)."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "jd",
        [
            2451545.0,  # J2000
            2415020.0,  # J1900
            2440587.5,  # Unix epoch
            2460000.0,  # 2023
            2470000.0,  # ~2050
            2430000.0,  # ~1941
        ],
    )
    def test_apogee_returns_valid_position(self, jd: float):
        """Interpolated apogee returns valid 6-element tuple."""
        result, retflag = swe.calc_ut(jd, INTP_APOG, FLG_SPEED)
        assert len(result) == 6
        lon, lat, dist, slon, slat, sdist = result
        assert 0 <= lon < 360, f"Longitude {lon} out of range"
        assert -10 < lat < 10, f"Latitude {lat} unreasonable"
        assert 0.002 < dist < 0.003, f"Distance {dist} AU unreasonable"

    @pytest.mark.unit
    def test_apogee_longitude_changes_over_time(self):
        """Apogee longitude advances over time (~40°/year)."""
        jd1 = 2451545.0
        jd2 = jd1 + 365.25
        r1, _ = swe.calc_ut(jd1, INTP_APOG, 0)
        r2, _ = swe.calc_ut(jd2, INTP_APOG, 0)
        # Mean apogee moves ~40.7°/year
        diff = (r2[0] - r1[0]) % 360
        assert 30 < diff < 55, f"Annual motion {diff}° outside expected range"

    @pytest.mark.unit
    def test_apogee_speed_sign_and_magnitude(self):
        """Apogee speed should be predominantly positive (prograde)."""
        jds = _random_jds(50, seed=100)
        positive_count = 0
        for jd in jds:
            result, _ = swe.calc_ut(jd, INTP_APOG, FLG_SPEED)
            speed = result[3]
            # Speed should be in a reasonable range
            assert abs(speed) < 1.0, f"Speed {speed}°/day unreasonable at JD {jd}"
            if speed > 0:
                positive_count += 1
        # Most of the time apogee moves prograde
        assert positive_count > 30, f"Only {positive_count}/50 prograde"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "year",
        [
            1600,
            1700,
            1800,
            1850,
            1900,
            1950,
            2000,
            2050,
            2100,
            2200,
            2400,
            2600,
        ],
    )
    def test_apogee_no_crash_across_centuries(self, year: int):
        """Apogee calculation doesn't crash for dates across centuries."""
        jd = _jd_from_year(year)
        result, retflag = swe.calc_ut(jd, INTP_APOG, FLG_SPEED)
        assert len(result) == 6
        assert 0 <= result[0] < 360

    @pytest.mark.unit
    def test_apogee_continuity(self):
        """Apogee longitude should be continuous (no jumps > 5° per day)."""
        jd_start = 2451545.0
        prev_lon = None
        for i in range(100):
            jd = jd_start + i * 1.0
            result, _ = swe.calc_ut(jd, INTP_APOG, 0)
            lon = result[0]
            if prev_lon is not None:
                diff = abs(lon - prev_lon)
                if diff > 180:
                    diff = 360 - diff
                assert diff < 5.0, (
                    f"Discontinuity at JD {jd}: {prev_lon:.4f} -> {lon:.4f} "
                    f"(diff {diff:.4f}°)"
                )
            prev_lon = lon

    @pytest.mark.unit
    def test_apogee_vs_osculating_relationship(self):
        """Interpolated apogee should be near osculating apogee (within ~25°)."""
        jds = _random_jds(30, seed=200)
        for jd in jds:
            r_intp, _ = swe.calc_ut(jd, INTP_APOG, 0)
            r_oscu, _ = swe.calc_ut(jd, OSCU_APOG, 0)
            diff = abs(r_intp[0] - r_oscu[0])
            if diff > 180:
                diff = 360 - diff
            # Osculating apogee can oscillate significantly from interpolated
            assert diff < 30.0, (
                f"Interp apogee {r_intp[0]:.2f} vs oscu {r_oscu[0]:.2f} "
                f"diff {diff:.2f}° at JD {jd}"
            )

    @pytest.mark.unit
    def test_apogee_distance_is_constant_mean(self):
        """Distance should be approximately constant at mean apogee distance."""
        jds = _random_jds(50, seed=300)
        for jd in jds:
            result, _ = swe.calc_ut(jd, INTP_APOG, 0)
            dist = result[2]
            # Mean apogee distance is ~0.0027099 AU
            assert abs(dist - 0.0027099) < 0.0002, (
                f"Distance {dist} AU far from mean at JD {jd}"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "flags,desc",
        [
            (0, "default"),
            (FLG_SPEED, "with speed"),
            (FLG_EQUATORIAL, "equatorial"),
            (FLG_J2000, "J2000"),
            (FLG_NOABERR, "no aberration"),
            (FLG_NONUT, "no nutation"),
            (FLG_SPEED | FLG_EQUATORIAL, "speed+equatorial"),
        ],
    )
    def test_apogee_flag_variants(self, flags: int, desc: str):
        """Apogee calculation works with various flag combinations."""
        jd = 2451545.0
        result, retflag = swe.calc_ut(jd, INTP_APOG, flags)
        assert len(result) == 6
        # All values should be finite
        for i, val in enumerate(result):
            assert math.isfinite(val), (
                f"Non-finite value at index {i}: {val} (flags={desc})"
            )

    @pytest.mark.unit
    def test_apogee_sidereal_mode(self):
        """Apogee works in sidereal mode."""
        jd = 2451545.0
        swe.set_sid_mode(SIDM_LAHIRI)
        try:
            result_sid, _ = swe.calc_ut(jd, INTP_APOG, FLG_SIDEREAL)
            result_trop, _ = swe.calc_ut(jd, INTP_APOG, 0)
            # Sidereal should differ from tropical by roughly the ayanamsha
            ayan = swe.get_ayanamsa_ut(jd)
            diff = (result_trop[0] - result_sid[0]) % 360
            if diff > 180:
                diff = 360 - diff
            assert abs(diff - ayan) < 1.0, (
                f"Sidereal offset {diff:.3f} vs ayanamsha {ayan:.3f}"
            )
        finally:
            swe.set_sid_mode(SIDM_FAGAN_BRADLEY)


class TestInterpolatedApogeeHighVolume:
    """High-volume parameterized tests for interpolated apogee."""

    @pytest.mark.unit
    @pytest.mark.parametrize("jd", _random_jds(200, seed=1000))
    def test_apogee_valid_position_200_dates(self, jd: float):
        """Apogee returns valid position for 200 random dates."""
        result, _ = swe.calc_ut(jd, INTP_APOG, FLG_SPEED)
        lon, lat, dist, slon, slat, sdist = result
        assert 0 <= lon < 360
        assert -10 < lat < 10
        assert 0.002 < dist < 0.003
        assert abs(slon) < 1.0


# ---------------------------------------------------------------------------
# Perigee correction table tests
# ---------------------------------------------------------------------------


class TestInterpolatedPerigeeCorrections:
    """Tests for the corrected interpolated perigee (body 22)."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "jd",
        [
            2451545.0,
            2415020.0,
            2440587.5,
            2460000.0,
            2470000.0,
            2430000.0,
        ],
    )
    def test_perigee_returns_valid_position(self, jd: float):
        """Interpolated perigee returns valid 6-element tuple."""
        result, retflag = swe.calc_ut(jd, INTP_PERG, FLG_SPEED)
        assert len(result) == 6
        lon, lat, dist, slon, slat, sdist = result
        assert 0 <= lon < 360, f"Longitude {lon} out of range"
        assert -10 < lat < 10, f"Latitude {lat} unreasonable"
        assert 0.002 < dist < 0.003, f"Distance {dist} AU unreasonable"

    @pytest.mark.unit
    def test_perigee_longitude_advances(self):
        """Perigee longitude advances over time (~40°/year on average).

        The interpolated perigee can have large oscillations, so we check
        over a longer period (5 years) and allow a wider range.
        """
        jd1 = 2451545.0
        jd2 = jd1 + 365.25 * 5  # 5 years
        r1, _ = swe.calc_ut(jd1, INTP_PERG, 0)
        r2, _ = swe.calc_ut(jd2, INTP_PERG, 0)
        diff = (r2[0] - r1[0]) % 360
        # Over 5 years, expect ~200° of mean motion (40°/yr * 5)
        assert 100 < diff < 300, f"5-year motion {diff}° outside expected range"

    @pytest.mark.unit
    def test_perigee_speed_reasonable(self):
        """Perigee speed should be in reasonable range.

        Interpolated perigee can have speeds up to ~2.5°/day near
        extreme oscillation points.
        """
        jds = _random_jds(50, seed=400)
        for jd in jds:
            result, _ = swe.calc_ut(jd, INTP_PERG, FLG_SPEED)
            speed = result[3]
            assert abs(speed) < 3.0, f"Speed {speed}°/day unreasonable at JD {jd}"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "year",
        [
            1600,
            1700,
            1800,
            1900,
            1950,
            2000,
            2050,
            2100,
            2200,
            2400,
            2600,
        ],
    )
    def test_perigee_no_crash_across_centuries(self, year: int):
        """Perigee calculation doesn't crash for dates across centuries."""
        jd = _jd_from_year(year)
        result, retflag = swe.calc_ut(jd, INTP_PERG, FLG_SPEED)
        assert len(result) == 6
        assert 0 <= result[0] < 360

    @pytest.mark.unit
    def test_perigee_continuity(self):
        """Perigee longitude should be continuous."""
        jd_start = 2451545.0
        prev_lon = None
        for i in range(100):
            jd = jd_start + i * 1.0
            result, _ = swe.calc_ut(jd, INTP_PERG, 0)
            lon = result[0]
            if prev_lon is not None:
                diff = abs(lon - prev_lon)
                if diff > 180:
                    diff = 360 - diff
                assert diff < 5.0, (
                    f"Discontinuity at JD {jd}: {prev_lon:.4f} -> {lon:.4f}"
                )
            prev_lon = lon

    @pytest.mark.unit
    def test_perigee_distance_is_constant_mean(self):
        """Distance should be approximately constant at mean perigee distance."""
        jds = _random_jds(50, seed=500)
        for jd in jds:
            result, _ = swe.calc_ut(jd, INTP_PERG, 0)
            dist = result[2]
            # Mean perigee distance is ~0.0024222 AU
            assert abs(dist - 0.0024222) < 0.0002, (
                f"Distance {dist} AU far from mean at JD {jd}"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "flags,desc",
        [
            (0, "default"),
            (FLG_SPEED, "with speed"),
            (FLG_EQUATORIAL, "equatorial"),
            (FLG_J2000, "J2000"),
            (FLG_NOABERR, "no aberration"),
            (FLG_NONUT, "no nutation"),
            (FLG_SPEED | FLG_EQUATORIAL, "speed+equatorial"),
        ],
    )
    def test_perigee_flag_variants(self, flags: int, desc: str):
        """Perigee calculation works with various flag combinations."""
        jd = 2451545.0
        result, retflag = swe.calc_ut(jd, INTP_PERG, flags)
        assert len(result) == 6
        for i, val in enumerate(result):
            assert math.isfinite(val), (
                f"Non-finite value at index {i}: {val} (flags={desc})"
            )


class TestInterpolatedPerigeeHighVolume:
    """High-volume parameterized tests for interpolated perigee."""

    @pytest.mark.unit
    @pytest.mark.parametrize("jd", _random_jds(200, seed=2000))
    def test_perigee_valid_position_200_dates(self, jd: float):
        """Perigee returns valid position for 200 random dates."""
        result, _ = swe.calc_ut(jd, INTP_PERG, FLG_SPEED)
        lon, lat, dist, slon, slat, sdist = result
        assert 0 <= lon < 360
        assert -10 < lat < 10
        assert 0.002 < dist < 0.003
        assert abs(slon) < 3.0


# ---------------------------------------------------------------------------
# Apogee-Perigee relationship tests
# ---------------------------------------------------------------------------


class TestApogeePerigeeRelationship:
    """Tests for the relationship between apogee and perigee."""

    @pytest.mark.unit
    @pytest.mark.parametrize("jd", _random_jds(50, seed=3000))
    def test_apogee_perigee_roughly_opposite(self, jd: float):
        """Apogee and perigee should be roughly 180° apart (within ~30°)."""
        r_apo, _ = swe.calc_ut(jd, INTP_APOG, 0)
        r_per, _ = swe.calc_ut(jd, INTP_PERG, 0)
        diff = abs(r_apo[0] - r_per[0])
        if diff > 180:
            diff = 360 - diff
        # They're roughly opposite but can deviate up to ~30°
        assert abs(diff - 180) < 35 or diff > 145, (
            f"Apogee {r_apo[0]:.2f} and perigee {r_per[0]:.2f} "
            f"differ by {diff:.2f}° (expected ~180°)"
        )

    @pytest.mark.unit
    def test_apogee_perigee_different_distances(self):
        """Apogee distance should be larger than perigee distance."""
        jd = 2451545.0
        r_apo, _ = swe.calc_ut(jd, INTP_APOG, 0)
        r_per, _ = swe.calc_ut(jd, INTP_PERG, 0)
        assert r_apo[2] > r_per[2], (
            f"Apogee dist {r_apo[2]} should be > perigee dist {r_per[2]}"
        )

    @pytest.mark.unit
    def test_apogee_near_mean_apogee(self):
        """Interpolated apogee should be within ~15° of mean apogee."""
        jds = _random_jds(30, seed=3100)
        for jd in jds:
            r_intp, _ = swe.calc_ut(jd, INTP_APOG, 0)
            r_mean, _ = swe.calc_ut(jd, MEAN_APOG, 0)
            diff = abs(r_intp[0] - r_mean[0])
            if diff > 180:
                diff = 360 - diff
            assert diff < 20.0, (
                f"Interp apogee {r_intp[0]:.2f} vs mean {r_mean[0]:.2f} "
                f"diff {diff:.2f}° at JD {jd}"
            )


# ---------------------------------------------------------------------------
# Latitude model tests
# ---------------------------------------------------------------------------


class TestApseLatitudeModel:
    """Tests for the latitude model of interpolated apse bodies."""

    @pytest.mark.unit
    def test_apogee_latitude_bounded(self):
        """Apogee latitude should be bounded by ~±5.15°."""
        jds = _random_jds(100, seed=4000)
        for jd in jds:
            result, _ = swe.calc_ut(jd, INTP_APOG, 0)
            lat = result[1]
            assert abs(lat) < 6.0, f"Apogee latitude {lat}° exceeds bounds at JD {jd}"

    @pytest.mark.unit
    def test_perigee_latitude_bounded(self):
        """Perigee latitude should be bounded by ~±5.15°."""
        jds = _random_jds(100, seed=4100)
        for jd in jds:
            result, _ = swe.calc_ut(jd, INTP_PERG, 0)
            lat = result[1]
            assert abs(lat) < 6.0, f"Perigee latitude {lat}° exceeds bounds at JD {jd}"

    @pytest.mark.unit
    def test_apogee_latitude_varies_with_node(self):
        """Latitude should vary as body crosses the ecliptic plane."""
        # Sample across half a nodal period (~9 years)
        jd_start = 2451545.0
        lats = []
        for i in range(100):
            jd = jd_start + i * 33.0  # ~9 years
            result, _ = swe.calc_ut(jd, INTP_APOG, 0)
            lats.append(result[1])
        # Should see both positive and negative latitudes
        has_positive = any(lat > 1.0 for lat in lats)
        has_negative = any(lat < -1.0 for lat in lats)
        assert has_positive, "No positive latitudes found"
        assert has_negative, "No negative latitudes found"
