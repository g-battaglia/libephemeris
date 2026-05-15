"""Tests for azalt, azalt_rev, refrac, and refrac_extended coordinate transforms."""

from __future__ import annotations

import math
import pytest
import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    JUPITER,
    SATURN,
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_EQUATORIAL,
)
from libephemeris.utils import (
    ECL2HOR,
    EQU2HOR,
    HOR2ECL,
    HOR2EQU,
    TRUE_TO_APP,
    APP_TO_TRUE,
)

JD_J2000 = 2451545.0  # 2000-01-01 12:00 TT
ROME = (12.5, 41.9, 50.0)  # lon, lat, alt
LONDON = (-0.1, 51.5, 11.0)
EQUATOR = (0.0, 0.0, 0.0)
TOKYO = (139.7, 35.7, 40.0)
SYDNEY = (151.2, -33.9, 58.0)


@pytest.mark.unit
class TestAzalt:
    """Test azalt horizontal coordinate conversion."""

    def test_sun_azalt_returns_3_values(self):
        """azalt should return a 3-element tuple."""
        sun, _ = swe.calc_ut(JD_J2000, SUN, FLG_SWIEPH)
        result = swe.azalt(
            JD_J2000, ECL2HOR, ROME, 1013.25, 15.0, (sun[0], sun[1], sun[2])
        )
        assert len(result) == 3

    def test_sun_azalt_values_finite(self):
        """All azalt output values should be finite."""
        sun, _ = swe.calc_ut(JD_J2000, SUN, FLG_SWIEPH)
        az, true_alt, app_alt = swe.azalt(
            JD_J2000, ECL2HOR, ROME, 1013.25, 15.0, (sun[0], sun[1], sun[2])
        )
        assert math.isfinite(az)
        assert math.isfinite(true_alt)
        assert math.isfinite(app_alt)

    def test_azimuth_range(self):
        """Azimuth should be in [0, 360)."""
        sun, _ = swe.calc_ut(JD_J2000, SUN, FLG_SWIEPH)
        az, _, _ = swe.azalt(
            JD_J2000, ECL2HOR, ROME, 1013.25, 15.0, (sun[0], sun[1], sun[2])
        )
        assert 0.0 <= az < 360.0, f"Azimuth {az} out of range"

    def test_altitude_range(self):
        """True altitude should be in [-90, 90]."""
        sun, _ = swe.calc_ut(JD_J2000, SUN, FLG_SWIEPH)
        _, true_alt, _ = swe.azalt(
            JD_J2000, ECL2HOR, ROME, 1013.25, 15.0, (sun[0], sun[1], sun[2])
        )
        assert -90.0 <= true_alt <= 90.0, f"True alt {true_alt} out of range"

    def test_apparent_geq_true_altitude(self):
        """Apparent altitude should be >= true altitude (refraction lifts)."""
        sun, _ = swe.calc_ut(JD_J2000, SUN, FLG_SWIEPH)
        _, true_alt, app_alt = swe.azalt(
            JD_J2000, ECL2HOR, ROME, 1013.25, 15.0, (sun[0], sun[1], sun[2])
        )
        # Only valid when object is above horizon (refraction adds to alt)
        if true_alt > -1.0:
            assert app_alt >= true_alt - 0.001, (
                f"Apparent alt {app_alt} < true alt {true_alt}"
            )

    def test_equatorial_input_flag(self):
        """azalt with EQU2HOR and equatorial coords should work."""
        sun, _ = swe.calc_ut(JD_J2000, SUN, FLG_SWIEPH | FLG_EQUATORIAL)
        result = swe.azalt(
            JD_J2000, EQU2HOR, ROME, 1013.25, 15.0, (sun[0], sun[1], sun[2])
        )
        assert len(result) == 3
        assert 0.0 <= result[0] < 360.0

    @pytest.mark.parametrize(
        "loc,name",
        [
            (ROME, "Rome"),
            (LONDON, "London"),
            (EQUATOR, "Equator"),
            (TOKYO, "Tokyo"),
            (SYDNEY, "Sydney"),
        ],
    )
    def test_azalt_different_locations(self, loc, name):
        """azalt should produce valid results at various locations."""
        sun, _ = swe.calc_ut(JD_J2000, SUN, FLG_SWIEPH)
        az, true_alt, app_alt = swe.azalt(
            JD_J2000, ECL2HOR, loc, 1013.25, 15.0, (sun[0], sun[1], sun[2])
        )
        assert 0.0 <= az < 360.0, f"Azimuth at {name}: {az}"
        assert -90.0 <= true_alt <= 90.0, f"True alt at {name}: {true_alt}"

    @pytest.mark.parametrize(
        "body,name",
        [
            (MOON, "Moon"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_azalt_different_bodies(self, body, name):
        """azalt should work for various bodies."""
        pos, _ = swe.calc_ut(JD_J2000, body, FLG_SWIEPH)
        az, true_alt, app_alt = swe.azalt(
            JD_J2000, ECL2HOR, ROME, 1013.25, 15.0, (pos[0], pos[1], pos[2])
        )
        assert math.isfinite(az), f"{name} azimuth not finite"
        assert math.isfinite(true_alt), f"{name} true alt not finite"


@pytest.mark.unit
class TestAzaltRev:
    """Test azalt_rev inverse horizontal conversion."""

    def test_azalt_rev_returns_2_values(self):
        """azalt_rev should return a 2-element tuple."""
        result = swe.azalt_rev(JD_J2000, HOR2EQU, ROME, 180.0, 45.0)
        assert len(result) == 2

    def test_azalt_rev_finite(self):
        """azalt_rev output should be finite."""
        ra, dec = swe.azalt_rev(JD_J2000, HOR2EQU, ROME, 180.0, 45.0)
        assert math.isfinite(ra)
        assert math.isfinite(dec)

    def test_azalt_rev_equatorial_range(self):
        """RA should be in [0, 360), Dec in [-90, 90]."""
        ra, dec = swe.azalt_rev(JD_J2000, HOR2EQU, ROME, 180.0, 45.0)
        assert 0.0 <= ra < 360.0, f"RA {ra} out of range"
        assert -90.0 <= dec <= 90.0, f"Dec {dec} out of range"

    def test_azalt_rev_ecliptic_mode(self):
        """azalt_rev with HOR2ECL should return ecliptic coords."""
        lon, lat = swe.azalt_rev(JD_J2000, HOR2ECL, ROME, 180.0, 45.0)
        assert math.isfinite(lon)
        assert math.isfinite(lat)

    def test_azalt_round_trip_equatorial(self):
        """azalt -> azalt_rev should approximately recover the input coords."""
        # Use equatorial coords for a known body
        sun, _ = swe.calc_ut(JD_J2000, SUN, FLG_SWIEPH | FLG_EQUATORIAL)
        ra_in, dec_in = sun[0], sun[1]

        # Forward: equatorial -> horizontal
        az, true_alt, _ = swe.azalt(
            JD_J2000, EQU2HOR, ROME, 1013.25, 15.0, (ra_in, dec_in, sun[2])
        )

        # Reverse: horizontal -> equatorial
        ra_out, dec_out = swe.azalt_rev(JD_J2000, HOR2EQU, ROME, az, true_alt)

        # Should recover within ~0.01° (refraction and precision limits)
        ra_diff = abs(ra_out - ra_in)
        if ra_diff > 180:
            ra_diff = 360 - ra_diff
        assert ra_diff < 0.05, f"RA round-trip error: {ra_diff}°"
        assert abs(dec_out - dec_in) < 0.05, (
            f"Dec round-trip error: {abs(dec_out - dec_in)}°"
        )

    @pytest.mark.parametrize("azimuth", [0.0, 90.0, 180.0, 270.0, 359.9])
    def test_azalt_rev_cardinal_azimuths(self, azimuth):
        """azalt_rev should handle cardinal azimuths without issues."""
        ra, dec = swe.azalt_rev(JD_J2000, HOR2EQU, ROME, azimuth, 30.0)
        assert math.isfinite(ra)
        assert math.isfinite(dec)


@pytest.mark.unit
class TestRefrac:
    """Test refrac atmospheric refraction."""

    def test_horizon_refraction(self):
        """Refraction at the horizon should be ~34 arcminutes."""
        app_alt = swe.refrac(0.0, 1013.25, 15.0, TRUE_TO_APP)
        refraction = app_alt  # true alt = 0, so refraction = apparent alt
        assert 0.4 < refraction < 0.7, (
            f"Horizon refraction {refraction}° (expected ~0.57° = 34')"
        )

    def test_refrac_true_to_app_positive(self):
        """TRUE_TO_APP should increase altitude (for positive alt)."""
        true_alt = 10.0
        app_alt = swe.refrac(true_alt, 1013.25, 15.0, TRUE_TO_APP)
        assert app_alt > true_alt, f"Apparent alt {app_alt} not > true alt {true_alt}"

    def test_refrac_decreases_with_altitude(self):
        """Refraction should decrease as altitude increases."""
        ref_10 = swe.refrac(10.0, 1013.25, 15.0, TRUE_TO_APP) - 10.0
        ref_45 = swe.refrac(45.0, 1013.25, 15.0, TRUE_TO_APP) - 45.0
        ref_80 = swe.refrac(80.0, 1013.25, 15.0, TRUE_TO_APP) - 80.0
        assert ref_10 > ref_45 > ref_80, (
            f"Refraction not decreasing: {ref_10}, {ref_45}, {ref_80}"
        )

    def test_refrac_zenith_tiny(self):
        """Refraction at the zenith should be very small."""
        ref_90 = swe.refrac(89.0, 1013.25, 15.0, TRUE_TO_APP) - 89.0
        assert ref_90 < 0.01, f"Zenith refraction {ref_90}° too large"

    def test_refrac_round_trip(self):
        """TRUE_TO_APP then APP_TO_TRUE should recover the original altitude."""
        true_alt = 15.0
        app_alt = swe.refrac(true_alt, 1013.25, 15.0, TRUE_TO_APP)
        recovered = swe.refrac(app_alt, 1013.25, 15.0, APP_TO_TRUE)
        assert abs(recovered - true_alt) < 0.01, (
            f"Round-trip error: {true_alt} -> {app_alt} -> {recovered}"
        )

    @pytest.mark.parametrize("true_alt", [0.0, 5.0, 15.0, 30.0, 60.0, 85.0])
    def test_refrac_round_trip_various(self, true_alt):
        """Round-trip should work for various altitudes."""
        app_alt = swe.refrac(true_alt, 1013.25, 15.0, TRUE_TO_APP)
        recovered = swe.refrac(app_alt, 1013.25, 15.0, APP_TO_TRUE)
        assert abs(recovered - true_alt) < 0.02, (
            f"Round-trip at {true_alt}°: error {abs(recovered - true_alt)}°"
        )

    def test_refrac_zero_pressure_no_refraction(self):
        """With zero atmospheric pressure, there should be no refraction."""
        app_alt = swe.refrac(10.0, 0.0, 15.0, TRUE_TO_APP)
        assert abs(app_alt - 10.0) < 0.001, (
            f"Zero pressure: expected 10.0, got {app_alt}"
        )

    def test_refrac_temperature_effect(self):
        """Cold air refracts more than warm air."""
        ref_cold = swe.refrac(10.0, 1013.25, -10.0, TRUE_TO_APP) - 10.0
        ref_warm = swe.refrac(10.0, 1013.25, 30.0, TRUE_TO_APP) - 10.0
        assert ref_cold > ref_warm, f"Cold refraction {ref_cold} not > warm {ref_warm}"


@pytest.mark.unit
class TestRefracExtended:
    """Test refrac_extended with observer altitude."""

    def test_returns_tuple_of_two(self):
        """refrac_extended should return (alt, details_tuple)."""
        result = swe.refrac_extended(10.0, 0.0, 1013.25, 15.0, 0.0065, TRUE_TO_APP)
        assert len(result) == 2
        alt, details = result
        assert math.isfinite(alt)
        assert len(details) == 4

    def test_details_all_finite(self):
        """All detail values should be finite."""
        _, details = swe.refrac_extended(
            10.0, 0.0, 1013.25, 15.0, 0.0065, TRUE_TO_APP
        )
        for i, val in enumerate(details):
            assert math.isfinite(val), f"Detail[{i}] = {val} not finite"

    def test_elevated_observer_dip(self):
        """Elevated observer should have negative dip (horizon below horizontal)."""
        _, details = swe.refrac_extended(
            0.0, 2000.0, 1013.25, 15.0, 0.0065, TRUE_TO_APP
        )
        dip = details[3]
        assert dip < 0.0, f"Dip at 2000m altitude: {dip}° (expected negative)"

    def test_sea_level_minimal_dip(self):
        """At sea level, dip should be ~0."""
        _, details = swe.refrac_extended(
            10.0, 0.0, 1013.25, 15.0, 0.0065, TRUE_TO_APP
        )
        dip = details[3]
        assert abs(dip) < 0.1, f"Sea level dip: {dip}° (expected ~0)"

    def test_refrac_extended_round_trip(self):
        """Extended refrac round-trip should recover original altitude."""
        true_alt = 15.0
        app_alt, _ = swe.refrac_extended(
            true_alt, 100.0, 1013.25, 15.0, 0.0065, TRUE_TO_APP
        )
        recovered, _ = swe.refrac_extended(
            app_alt, 100.0, 1013.25, 15.0, 0.0065, APP_TO_TRUE
        )
        assert abs(recovered - true_alt) < 0.02, (
            f"Extended round-trip: {true_alt} -> {app_alt} -> {recovered}"
        )
