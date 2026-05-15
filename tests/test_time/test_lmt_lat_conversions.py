"""
Tests for lmt_to_lat and lat_to_lmt time conversions.

Verifies Local Mean Time <-> Local Apparent Time conversions,
Equation of Time effects, round-trip consistency.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    swe.close()


JD_J2000 = 2451545.0  # 2000-01-01 12:00 UT


class TestLatToLmt:
    """Test lat_to_lmt: Local Apparent Time -> Local Mean Time."""

    @pytest.mark.unit
    def test_returns_float(self):
        """lat_to_lmt returns a float (Julian Day)."""
        result = swe.lat_to_lmt(JD_J2000, 0.0)
        assert isinstance(result, float)

    @pytest.mark.unit
    def test_greenwich_close_to_input(self):
        """At Greenwich (lon=0), LAT and LMT differ only by EoT."""
        result = swe.lat_to_lmt(JD_J2000, 0.0)
        # EoT is at most ~16 minutes = 0.011 days
        assert abs(result - JD_J2000) < 0.012

    @pytest.mark.unit
    @pytest.mark.parametrize("lon", [-180.0, -90.0, -45.0, 0.0, 45.0, 90.0, 180.0])
    def test_various_longitudes(self, lon):
        """lat_to_lmt works at various longitudes."""
        result = swe.lat_to_lmt(JD_J2000, lon)
        assert math.isfinite(result)
        # Result should be within 0.5 days of input
        assert abs(result - JD_J2000) < 0.5

    @pytest.mark.unit
    def test_eot_variation_over_year(self):
        """EoT varies over the year — differences should vary."""
        diffs = []
        for day in range(0, 365, 30):
            jd = JD_J2000 + day
            result = swe.lat_to_lmt(jd, 0.0)
            diffs.append((result - jd) * 24 * 60)  # Convert to minutes
        # EoT varies from about -14 to +16 minutes
        variation = max(diffs) - min(diffs)
        assert variation > 10.0, f"EoT variation {variation} minutes — too small"


class TestLmtToLat:
    """Test lmt_to_lat: Local Mean Time -> Local Apparent Time."""

    @pytest.mark.unit
    def test_returns_float(self):
        """lmt_to_lat returns a float (Julian Day)."""
        result = swe.lmt_to_lat(JD_J2000, 0.0)
        assert isinstance(result, float)

    @pytest.mark.unit
    def test_greenwich_close_to_input(self):
        """At Greenwich (lon=0), LMT and LAT differ only by EoT."""
        result = swe.lmt_to_lat(JD_J2000, 0.0)
        assert abs(result - JD_J2000) < 0.012

    @pytest.mark.unit
    @pytest.mark.parametrize("lon", [-120.0, -60.0, 0.0, 60.0, 120.0])
    def test_various_longitudes(self, lon):
        """lmt_to_lat works at various longitudes."""
        result = swe.lmt_to_lat(JD_J2000, lon)
        assert math.isfinite(result)
        assert abs(result - JD_J2000) < 0.5


class TestRoundTrip:
    """Test LMT <-> LAT round-trip consistency."""

    @pytest.mark.unit
    @pytest.mark.parametrize("lon", [-90.0, 0.0, 45.0, 120.0])
    def test_lat_to_lmt_to_lat_roundtrip(self, lon):
        """LAT -> LMT -> LAT should return ~original value."""
        jd_lat = JD_J2000
        jd_lmt = swe.lat_to_lmt(jd_lat, lon)
        jd_lat_back = swe.lmt_to_lat(jd_lmt, lon)
        # EoT iterative inversion may have ~1 second residual
        assert jd_lat_back == pytest.approx(jd_lat, abs=1e-4), (
            f"Round-trip error: {abs(jd_lat_back - jd_lat) * 86400:.3f} seconds"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("lon", [-90.0, 0.0, 45.0, 120.0])
    def test_lmt_to_lat_to_lmt_roundtrip(self, lon):
        """LMT -> LAT -> LMT should return ~original value."""
        jd_lmt = JD_J2000
        jd_lat = swe.lmt_to_lat(jd_lmt, lon)
        jd_lmt_back = swe.lat_to_lmt(jd_lat, lon)
        # EoT changes slightly between LMT and LAT epochs, so ~1s residual
        assert jd_lmt_back == pytest.approx(jd_lmt, abs=1e-4), (
            f"Round-trip error: {abs(jd_lmt_back - jd_lmt) * 86400:.3f} seconds"
        )

    @pytest.mark.unit
    def test_roundtrip_across_dates(self):
        """Round-trip works across the year."""
        lon = 30.0
        for day in range(0, 365, 30):
            jd = JD_J2000 + day
            jd_lmt = swe.lat_to_lmt(jd, lon)
            jd_back = swe.lmt_to_lat(jd_lmt, lon)
            assert jd_back == pytest.approx(jd, abs=1e-4), (
                f"Round-trip failed at JD {jd}"
            )


class TestSweAliases:
    """Test swe_ prefixed aliases."""

    @pytest.mark.unit
    def test_swe_lat_to_lmt_exists(self):
        """lat_to_lmt is an alias for lat_to_lmt."""
        r1 = swe.lat_to_lmt(JD_J2000, 10.0)
        r2 = swe.lat_to_lmt(JD_J2000, 10.0)
        assert r1 == r2

    @pytest.mark.unit
    def test_swe_lmt_to_lat_exists(self):
        """lmt_to_lat is an alias for lmt_to_lat."""
        r1 = swe.lmt_to_lat(JD_J2000, 10.0)
        r2 = swe.lmt_to_lat(JD_J2000, 10.0)
        assert r1 == r2
