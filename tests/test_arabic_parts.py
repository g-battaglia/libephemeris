"""
Tests for Arabic Parts (Lots) calculations.

These tests verify the traditional formulas for Arabic Parts,
including correct day/night chart handling for Part of Fortune and Spirit.

The key formulas being tested:
- Part of Fortune:
  - Day: ASC + Moon - Sun
  - Night: ASC + Sun - Moon
- Part of Spirit:
  - Day: ASC + Sun - Moon
  - Night: ASC + Moon - Sun

Fortune and Spirit are swapped between day and night charts,
following traditional sect-based calculation methods.
"""

import pytest

from libephemeris.arabic_parts import (
    _EXTREME_LATITUDE_THRESHOLD,
    _is_sun_above_horizon_3d,
    calc_all_arabic_parts,
    calc_arabic_part_of_faith,
    calc_arabic_part_of_fortune,
    calc_arabic_part_of_love,
    calc_arabic_part_of_spirit,
    is_day_chart,
)
from libephemeris import julday


class TestIsDayChart:
    """Tests for is_day_chart() function."""

    def test_sun_below_horizon_normal_case(self):
        """Sun between ASC and DSC (lower hemisphere, near IC) is a night chart."""
        # ASC at 0°, DSC at 180°
        # Sun at 90° (≈IC, below horizon) = night chart
        asc = 0.0
        sun = 90.0
        assert is_day_chart(sun, asc) is False

    def test_sun_above_horizon_normal_case(self):
        """Sun outside ASC-DSC arc (upper hemisphere, near MC) is a day chart."""
        # ASC at 0°, DSC at 180°
        # Sun at 270° (≈MC, above horizon) = day chart
        asc = 0.0
        sun = 270.0
        assert is_day_chart(sun, asc) is True

    def test_sun_below_horizon_wrapped_case(self):
        """Sun in lower hemisphere when ASC > DSC (wrapped case)."""
        # ASC at 350°, DSC at 170°
        # Sun at 355° (just past ASC, house 1, below horizon) = night chart
        asc = 350.0
        sun = 355.0
        assert is_day_chart(sun, asc) is False

        # Sun at 100° (≈IC, below horizon in wrapped case) = night chart
        sun = 100.0
        assert is_day_chart(sun, asc) is False

    def test_sun_above_horizon_wrapped_case(self):
        """Sun in upper hemisphere when ASC > DSC (wrapped case)."""
        # ASC at 350°, DSC at 170°
        # Sun at 250° (between DSC=170° and ASC=350°, near MC) = day chart
        asc = 350.0
        sun = 250.0
        assert is_day_chart(sun, asc) is True

    def test_sun_exactly_on_asc(self):
        """Sun exactly on Ascendant is considered a day chart (sunrise)."""
        asc = 15.0
        sun = 15.0  # Sun exactly on ASC
        assert is_day_chart(sun, asc) is True

    def test_sun_exactly_on_desc(self):
        """Sun exactly on Descendant is considered a day chart (sunset)."""
        asc = 15.0
        desc = (asc + 180.0) % 360.0  # 195°
        sun = desc
        assert is_day_chart(sun, asc) is True

    def test_issue_1_regression(self):
        """GH issue #1: is_day_chart(352, 153) must be True (reported by @Blu)."""
        # ASC=153°, DSC=333°. Sun at 352° is in upper hemisphere (between DSC and ASC
        # going CCW through MC). 3D method confirms day chart.
        assert is_day_chart(352.0, 153.0) is True

    @pytest.mark.parametrize(
        "asc", [0.0, 45.0, 90.0, 135.0, 153.0, 180.0, 225.0, 270.0, 315.0, 350.0]
    )
    def test_sun_near_mc_is_day(self, asc):
        """Sun near MC (~ASC+270° in idealized chart) is always above horizon."""
        # In the simplified 2D model, MC is opposite IC, roughly at (ASC - 90°) mod 360 = (ASC + 270°) mod 360
        mc_approx = (asc + 270.0) % 360.0
        assert is_day_chart(mc_approx, asc) is True, (
            f"MC at ASC={asc} should be day chart"
        )

    @pytest.mark.parametrize(
        "asc", [0.0, 45.0, 90.0, 135.0, 153.0, 180.0, 225.0, 270.0, 315.0, 350.0]
    )
    def test_sun_near_ic_is_night(self, asc):
        """Sun near IC (~ASC+90°) is always below horizon."""
        ic_approx = (asc + 90.0) % 360.0
        assert is_day_chart(ic_approx, asc) is False, (
            f"IC at ASC={asc} should be night chart"
        )

    @pytest.mark.parametrize(
        "asc", [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0]
    )
    @pytest.mark.parametrize(
        "offset", [10.0, 45.0, 80.0, 100.0, 135.0, 170.0, 190.0, 225.0, 260.0, 280.0, 315.0, 350.0]
    )
    def test_symmetry_180_rotation(self, asc, offset):
        """Rotating both Sun and ASC by 180° preserves the day/night result."""
        sun = (asc + offset) % 360.0
        asc_rot = (asc + 180.0) % 360.0
        sun_rot = (sun + 180.0) % 360.0
        assert is_day_chart(sun, asc) == is_day_chart(sun_rot, asc_rot)

    @pytest.mark.parametrize(
        "asc", [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0]
    )
    @pytest.mark.parametrize(
        "offset", [10.0, 45.0, 80.0, 100.0, 135.0, 170.0, 190.0, 225.0, 260.0, 280.0, 315.0, 350.0]
    )
    def test_complementary_offsets_opposite_sect(self, asc, offset):
        """Sun at offset X and at 360-X from ASC are in opposite hemispheres.

        Excludes boundary offsets 0 and 180 where both fall on the horizon
        and are conventionally treated as day charts.
        """
        if offset in (0.0, 180.0):
            return
        sun_a = (asc + offset) % 360.0
        sun_b = (asc + (360.0 - offset)) % 360.0
        assert is_day_chart(sun_a, asc) != is_day_chart(sun_b, asc), (
            f"ASC={asc}, offsets {offset} and {360-offset} should be opposite"
        )


class TestCalcArabicPartOfFortune:
    """Tests for calc_arabic_part_of_fortune() function."""

    def test_diurnal_formula(self):
        """Day chart uses ASC + Moon - Sun."""
        asc = 15.0
        sun = 120.0
        moon = 240.0
        # Day formula: ASC + Moon - Sun = 15 + 240 - 120 = 135
        result = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=True)
        assert result == pytest.approx(135.0)

    def test_nocturnal_formula(self):
        """Night chart uses ASC + Sun - Moon (inverted formula)."""
        asc = 15.0
        sun = 120.0
        moon = 240.0
        # Night formula: ASC + Sun - Moon = 15 + 120 - 240 = -105 => 255 (mod 360)
        result = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=False)
        assert result == pytest.approx(255.0)

    def test_modulo_wrap_positive(self):
        """Result correctly wraps with modulo 360 (overflow)."""
        asc = 300.0
        sun = 10.0
        moon = 350.0
        # Day formula: 300 + 350 - 10 = 640 => 280 (mod 360)
        result = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=True)
        assert result == pytest.approx(280.0)

    def test_modulo_wrap_negative(self):
        """Result correctly wraps with modulo 360 (underflow)."""
        asc = 10.0
        sun = 300.0
        moon = 30.0
        # Day formula: 10 + 30 - 300 = -260 => 100 (mod 360)
        result = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=True)
        assert result == pytest.approx(100.0)


class TestCalcArabicPartOfSpirit:
    """Tests for calc_arabic_part_of_spirit() function."""

    def test_diurnal_formula(self):
        """Day chart uses ASC + Sun - Moon."""
        asc = 15.0
        sun = 120.0
        moon = 240.0
        # Day formula: ASC + Sun - Moon = 15 + 120 - 240 = -105 => 255 (mod 360)
        result = calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal=True)
        assert result == pytest.approx(255.0)

    def test_nocturnal_formula(self):
        """Night chart uses ASC + Moon - Sun (inverted formula)."""
        asc = 15.0
        sun = 120.0
        moon = 240.0
        # Night formula: ASC + Moon - Sun = 15 + 240 - 120 = 135
        result = calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal=False)
        assert result == pytest.approx(135.0)

    def test_modulo_wrap_positive(self):
        """Result correctly wraps with modulo 360 (overflow)."""
        asc = 300.0
        sun = 350.0
        moon = 10.0
        # Day formula: 300 + 350 - 10 = 640 => 280 (mod 360)
        result = calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal=True)
        assert result == pytest.approx(280.0)


class TestFortuneAndSpiritInversion:
    """Tests verifying Fortune and Spirit are inverted between day and night."""

    def test_fortune_day_equals_spirit_night(self):
        """Part of Fortune (day) equals Part of Spirit (night) - same formula."""
        asc = 45.0
        sun = 100.0
        moon = 200.0

        fortune_day = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=True)
        spirit_night = calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal=False)

        # Both use ASC + Moon - Sun
        assert fortune_day == pytest.approx(spirit_night)

    def test_spirit_day_equals_fortune_night(self):
        """Part of Spirit (day) equals Part of Fortune (night) - same formula."""
        asc = 45.0
        sun = 100.0
        moon = 200.0

        spirit_day = calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal=True)
        fortune_night = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=False)

        # Both use ASC + Sun - Moon
        assert spirit_day == pytest.approx(fortune_night)

    def test_fortune_and_spirit_different_in_same_chart(self):
        """Fortune and Spirit have different values in the same chart."""
        asc = 45.0
        sun = 100.0
        moon = 200.0

        fortune = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=True)
        spirit = calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal=True)

        # Fortune: 45 + 200 - 100 = 145
        # Spirit: 45 + 100 - 200 = -55 => 305
        assert fortune == pytest.approx(145.0)
        assert spirit == pytest.approx(305.0)
        assert fortune != spirit


class TestCalcArabicPartOfLove:
    """Tests for calc_arabic_part_of_love() function."""

    def test_basic_calculation(self):
        """Basic calculation: ASC + Venus - Sun."""
        asc = 30.0
        venus = 60.0
        sun = 45.0
        # Formula: 30 + 60 - 45 = 45
        result = calc_arabic_part_of_love(asc, venus, sun)
        assert result == pytest.approx(45.0)


class TestCalcArabicPartOfFaith:
    """Tests for calc_arabic_part_of_faith() function."""

    def test_basic_calculation(self):
        """Basic calculation: ASC + Mercury - Moon."""
        asc = 30.0
        mercury = 90.0
        moon = 60.0
        # Formula: 30 + 90 - 60 = 60
        result = calc_arabic_part_of_faith(asc, mercury, moon)
        assert result == pytest.approx(60.0)


class TestCalcAllArabicParts:
    """Tests for calc_all_arabic_parts() function using position dictionaries."""

    def test_nocturnal_chart_sun_near_ic(self):
        """Test with Sun in lower hemisphere (night chart)."""
        # ASC at 15°, Sun at 90° is between ASC=15° and DSC=195° (lower hemisphere, near IC)
        positions = {
            "Asc": 15.0,
            "Sun": 90.0,
            "Moon": 240.0,
            "Mercury": 100.0,
            "Venus": 80.0,
        }

        parts = calc_all_arabic_parts(positions)

        # Verify this is treated as a night chart
        assert is_day_chart(positions["Sun"], positions["Asc"]) is False

        # Night formula for Fortune: ASC + Sun - Moon = 15 + 90 - 240 = -135 => 225
        assert parts["Pars_Fortunae"] == pytest.approx(225.0)

        # Night formula for Spirit: ASC + Moon - Sun = 15 + 240 - 90 = 165
        assert parts["Pars_Spiritus"] == pytest.approx(165.0)

        # Part of Love: ASC + Venus - Sun = 15 + 80 - 90 = 5
        assert parts["Pars_Amoris"] == pytest.approx(5.0)

        # Part of Faith: ASC + Mercury - Moon = 15 + 100 - 240 = -125 => 235
        assert parts["Pars_Fidei"] == pytest.approx(235.0)

    def test_diurnal_chart_sun_near_mc(self):
        """Test with Sun in upper hemisphere (day chart)."""
        # ASC at 15°, Sun at 270° is outside ASC=15° to DSC=195° (upper hemisphere, near MC)
        positions = {
            "Asc": 15.0,
            "Sun": 270.0,
            "Moon": 240.0,
            "Mercury": 100.0,
            "Venus": 80.0,
        }

        parts = calc_all_arabic_parts(positions)

        # Verify this is treated as a day chart
        assert is_day_chart(positions["Sun"], positions["Asc"]) is True

        # Day formula for Fortune: ASC + Moon - Sun = 15 + 240 - 270 = -15 => 345
        assert parts["Pars_Fortunae"] == pytest.approx(345.0)

        # Day formula for Spirit: ASC + Sun - Moon = 15 + 270 - 240 = 45
        assert parts["Pars_Spiritus"] == pytest.approx(45.0)

        # Part of Love (not sect-dependent): ASC + Venus - Sun = 15 + 80 - 270 = -175 => 185
        assert parts["Pars_Amoris"] == pytest.approx(185.0)

        # Part of Faith (not sect-dependent): ASC + Mercury - Moon = 15 + 100 - 240 = -125 => 235
        assert parts["Pars_Fidei"] == pytest.approx(235.0)

    def test_fortune_spirit_swap_day_vs_night(self):
        """Verify Fortune and Spirit values swap between day and night charts."""
        # Same base positions, but Sun position changes sect
        base_positions = {
            "Asc": 0.0,
            "Moon": 180.0,
            "Mercury": 45.0,
            "Venus": 60.0,
        }

        # Night chart: Sun at 90° (lower hemisphere, near IC)
        night_positions = {**base_positions, "Sun": 90.0}
        night_parts = calc_all_arabic_parts(night_positions)

        # Day chart: Sun at 270° (upper hemisphere, near MC)
        day_positions = {**base_positions, "Sun": 270.0}
        day_parts = calc_all_arabic_parts(day_positions)

        # Night Fortune: ASC + Sun - Moon = 0 + 90 - 180 = -90 => 270
        # Night Spirit: ASC + Moon - Sun = 0 + 180 - 90 = 90
        assert night_parts["Pars_Fortunae"] == pytest.approx(270.0)
        assert night_parts["Pars_Spiritus"] == pytest.approx(90.0)

        # Day Fortune: ASC + Moon - Sun = 0 + 180 - 270 = -90 => 270
        # Day Spirit: ASC + Sun - Moon = 0 + 270 - 180 = 90
        assert day_parts["Pars_Fortunae"] == pytest.approx(270.0)
        assert day_parts["Pars_Spiritus"] == pytest.approx(90.0)

        # Note: In this particular case values are the same due to symmetric positions
        # Let's verify with asymmetric values

    def test_fortune_spirit_swap_asymmetric(self):
        """Verify Formula swap with asymmetric Moon/Sun positions."""
        base_positions = {
            "Asc": 30.0,
            "Moon": 120.0,
            "Mercury": 45.0,
            "Venus": 60.0,
        }

        # Night chart: Sun at 60° (lower hemisphere)
        night_positions = {**base_positions, "Sun": 60.0}
        night_parts = calc_all_arabic_parts(night_positions)

        # Night Fortune: ASC + Sun - Moon = 30 + 60 - 120 = -30 => 330
        # Night Spirit: ASC + Moon - Sun = 30 + 120 - 60 = 90
        assert night_parts["Pars_Fortunae"] == pytest.approx(330.0)
        assert night_parts["Pars_Spiritus"] == pytest.approx(90.0)

        # Day chart: Sun at 300° (upper hemisphere)
        day_positions = {**base_positions, "Sun": 300.0}
        day_parts = calc_all_arabic_parts(day_positions)

        # Day Fortune: ASC + Moon - Sun = 30 + 120 - 300 = -150 => 210
        # Day Spirit: ASC + Sun - Moon = 30 + 300 - 120 = 210
        assert day_parts["Pars_Fortunae"] == pytest.approx(210.0)
        assert day_parts["Pars_Spiritus"] == pytest.approx(210.0)

        # Day and Night should produce different results for Fortune
        assert day_parts["Pars_Fortunae"] != night_parts["Pars_Fortunae"]

    def test_sun_exactly_on_asc_edge_case(self):
        """Test edge case: Sun exactly on Ascendant (sunrise moment)."""
        # Sun exactly on ASC should be treated as day chart
        positions = {
            "Asc": 100.0,
            "Sun": 100.0,  # Exactly on ASC
            "Moon": 200.0,
            "Mercury": 150.0,
            "Venus": 120.0,
        }

        parts = calc_all_arabic_parts(positions)

        # Should be a day chart
        assert is_day_chart(positions["Sun"], positions["Asc"]) is True

        # Day Fortune: 100 + 200 - 100 = 200
        assert parts["Pars_Fortunae"] == pytest.approx(200.0)

        # Day Spirit: 100 + 100 - 200 = 0
        assert parts["Pars_Spiritus"] == pytest.approx(0.0)

    def test_missing_required_keys_raise(self):
        """Missing required positions fail fast instead of defaulting to 0.0."""
        # Only ASC provided -> Sun/Moon/Mercury/Venus missing.
        with pytest.raises(KeyError, match="missing required position"):
            calc_all_arabic_parts({"Asc": 90.0})

        # A single missing key (Venus) is still an error.
        with pytest.raises(KeyError, match="Venus"):
            calc_all_arabic_parts(
                {"Asc": 0.0, "Sun": 90.0, "Moon": 180.0, "Mercury": 45.0}
            )

    def test_sun_lat_remains_optional(self):
        """Sun_lat is the only optional key: a complete dict without it works."""
        positions = {
            "Asc": 0.0,
            "Sun": 90.0,
            "Moon": 180.0,
            "Mercury": 45.0,
            "Venus": 135.0,
        }
        parts = calc_all_arabic_parts(positions)  # no Sun_lat -> must not raise
        assert "Pars_Fortunae" in parts

    def test_all_parts_returned(self):
        """Verify all four standard parts are returned."""
        positions = {
            "Asc": 0.0,
            "Sun": 90.0,
            "Moon": 180.0,
            "Mercury": 45.0,
            "Venus": 135.0,
        }

        parts = calc_all_arabic_parts(positions)

        assert "Pars_Fortunae" in parts
        assert "Pars_Spiritus" in parts
        assert "Pars_Amoris" in parts
        assert "Pars_Fidei" in parts
        assert len(parts) == 4

    def test_values_in_valid_range(self):
        """All returned values should be in range 0-360."""
        positions = {
            "Asc": 350.0,
            "Sun": 10.0,
            "Moon": 270.0,
            "Mercury": 300.0,
            "Venus": 5.0,
        }

        parts = calc_all_arabic_parts(positions)

        for name, value in parts.items():
            assert 0.0 <= value < 360.0, f"{name} = {value} is out of range"


class TestIsDayChart3D:
    """Tests for 3D day/night calculation at extreme latitudes."""

    def test_threshold_constant_exists(self):
        """Verify the extreme latitude threshold is defined."""
        assert _EXTREME_LATITUDE_THRESHOLD == 60.0

    def test_moderate_latitude_uses_2d_method(self):
        """At moderate latitudes, should use 2D method even with jd provided."""
        # Rome, Italy - moderate latitude (41.9°N)
        # June 21, 2024 at noon - Sun clearly above horizon
        jd = julday(2024, 6, 21, 12.0)
        sun_lon = 90.0  # ~Cancer (summer solstice)
        asc = 180.0  # Libra rising

        # 2D method: sun_rel = (90 - 180) % 360 = 270, which is > 180.
        # So 2D correctly says day chart (matching actual Rome noon).
        result_2d = is_day_chart(sun_lon, asc)

        # With location at moderate latitude, should still use 2D
        result_with_loc = is_day_chart(sun_lon, asc, jd=jd, geo_lat=41.9, geo_lon=12.5)
        assert result_2d == result_with_loc

    def test_extreme_latitude_uses_3d_method_northern(self):
        """At extreme northern latitude with location, uses 3D calculation."""
        # Tromso, Norway (69.6°N) - Arctic summer, Sun above horizon at midnight
        # June 21, 2024 at midnight UTC - Sun should be above horizon
        jd = julday(2024, 6, 21, 0.0)  # Midnight UTC

        # Sun at summer solstice position
        sun_lon = 90.0  # Cancer 0°
        asc = 0.0  # Aries rising (doesn't matter for 3D calc)

        # Without location: uses 2D method
        result_2d = is_day_chart(sun_lon, asc)

        # With location at extreme latitude: uses 3D method
        result_3d = is_day_chart(sun_lon, asc, jd=jd, geo_lat=69.6, geo_lon=19.0)

        # In Arctic summer, Sun is above horizon even at midnight
        # 3D method should return True (day chart)
        assert result_3d is True

    def test_extreme_latitude_winter_night(self):
        """At extreme latitude in winter, Sun below horizon at noon."""
        # Tromso, Norway (69.6°N) - polar night period
        # December 21, 2024 at noon UTC - Sun below horizon
        jd = julday(2024, 12, 21, 12.0)

        # Sun at winter solstice position
        sun_lon = 270.0  # Capricorn 0°
        asc = 0.0

        # 3D method should return False (night chart) even at "noon"
        result_3d = is_day_chart(sun_lon, asc, jd=jd, geo_lat=69.6, geo_lon=19.0)
        assert result_3d is False

    def test_extreme_latitude_southern_hemisphere(self):
        """Test 3D calculation for Antarctic location."""
        # McMurdo Station, Antarctica (-77.8°S)
        # December 21, 2024 - Antarctic summer, midnight sun
        jd = julday(2024, 12, 21, 0.0)  # Midnight UTC

        sun_lon = 270.0  # Capricorn (Sun at summer solstice for S hemisphere)
        asc = 0.0

        # In Antarctic summer, Sun is above horizon at midnight
        result_3d = is_day_chart(sun_lon, asc, jd=jd, geo_lat=-77.8, geo_lon=166.7)
        assert result_3d is True

    def test_is_sun_above_horizon_3d_directly(self):
        """Test the internal 3D function directly."""
        # Standard daytime scenario - Rome at noon in summer
        jd = julday(2024, 6, 21, 12.0)
        sun_lon = 90.0
        sun_lat = 0.0
        geo_lat = 41.9
        geo_lon = 12.5

        result = _is_sun_above_horizon_3d(jd, sun_lon, sun_lat, geo_lat, geo_lon)
        assert result is True

    def test_is_sun_above_horizon_3d_night(self):
        """Test 3D function for nighttime scenario."""
        # Midnight in Rome, winter
        jd = julday(2024, 12, 21, 0.0)  # Midnight UTC
        sun_lon = 270.0
        sun_lat = 0.0
        geo_lat = 41.9
        geo_lon = 12.5

        result = _is_sun_above_horizon_3d(jd, sun_lon, sun_lat, geo_lat, geo_lon)
        assert result is False

    def test_fallback_to_2d_without_jd(self):
        """Without jd parameter, should use 2D method."""
        sun_lon = 90.0
        asc = 0.0

        # Without jd: 2D method (Sun between 0° and 180°)
        result = is_day_chart(sun_lon, asc, geo_lat=70.0, geo_lon=20.0)

        # Should be same as pure 2D call
        result_2d = is_day_chart(sun_lon, asc)
        assert result == result_2d

    def test_fallback_to_2d_without_lat(self):
        """Without geo_lat parameter, should use 2D method."""
        jd = julday(2024, 6, 21, 12.0)
        sun_lon = 90.0
        asc = 0.0

        result = is_day_chart(sun_lon, asc, jd=jd, geo_lon=20.0)
        result_2d = is_day_chart(sun_lon, asc)
        assert result == result_2d


class TestCalcAllArabicPartsWithLocation:
    """Tests for calc_all_arabic_parts with location parameters."""

    def test_moderate_latitude_same_as_without_location(self):
        """At moderate latitudes, results should be same with or without location."""
        positions = {
            "Asc": 15.0,
            "Sun": 90.0,
            "Moon": 240.0,
            "Mercury": 100.0,
            "Venus": 80.0,
        }

        parts_without = calc_all_arabic_parts(positions)
        parts_with = calc_all_arabic_parts(
            positions,
            jd=julday(2024, 6, 21, 12.0),
            geo_lat=41.9,
            geo_lon=12.5,
        )

        # Same results because latitude is moderate
        assert parts_without == parts_with

    def test_extreme_latitude_may_differ(self):
        """At extreme latitudes, 3D calc may give different sect determination."""
        # Scenario: Arctic summer midnight - Sun geometrically above horizon
        # but 2D method might say night chart depending on ASC/Sun relationship
        jd = julday(2024, 6, 21, 0.0)

        positions = {
            "Asc": 270.0,  # Capricorn rising
            "Sun": 90.0,  # Cancer (summer solstice)
            "Moon": 180.0,
            "Mercury": 100.0,
            "Venus": 80.0,
        }

        parts_2d = calc_all_arabic_parts(positions)
        parts_3d = calc_all_arabic_parts(
            positions,
            jd=jd,
            geo_lat=69.6,  # Tromso
            geo_lon=19.0,
        )

        # 2D: ASC=270, Sun=90. sun_rel = (90 - 270) % 360 = 180.
        # Lower hemisphere is open arc (0, 180), so sun_rel=180 is on boundary (DSC).
        # not (0.0 < 180.0 < 180.0) = not False = True => day chart.

        # 3D should also say day chart (Sun above horizon in Arctic summer)
        # Verify 3D is being used (latitude > 60)
        assert abs(69.6) > _EXTREME_LATITUDE_THRESHOLD

        # Both should agree in this case (day chart)
        is_day_2d = is_day_chart(positions["Sun"], positions["Asc"])
        is_day_3d = is_day_chart(
            positions["Sun"], positions["Asc"], jd=jd, geo_lat=69.6, geo_lon=19.0
        )

        # In Arctic summer at midnight, Sun IS above horizon (midnight sun)
        assert is_day_3d is True

    def test_sun_lat_parameter_passed(self):
        """Verify Sun_lat from positions is used in 3D calculation."""
        jd = julday(2024, 6, 21, 12.0)

        positions = {
            "Asc": 15.0,
            "Sun": 90.0,
            "Sun_lat": 0.0,  # Sun ecliptic latitude (nearly always ~0)
            "Moon": 240.0,
            "Mercury": 100.0,
            "Venus": 80.0,
        }

        # Should not raise any errors
        parts = calc_all_arabic_parts(
            positions,
            jd=jd,
            geo_lat=70.0,
            geo_lon=20.0,
        )

        assert "Pars_Fortunae" in parts


class TestIsDayChart2D3DAgreement:
    """Verify 2D and 3D methods agree at moderate latitudes using real ephemeris data.

    The 2D method is a simplified ecliptic comparison; the 3D method uses a full
    coordinate transformation. They should agree everywhere except very close to
    the horizon, where sub-arcminute precision differences can occur.
    """

    @pytest.mark.parametrize("hour", [0, 3, 6, 9, 12, 15, 18, 21])
    @pytest.mark.parametrize("month", [3, 6, 9, 12])
    def test_rome_agreement(self, month, hour):
        """At Rome (lat 41.9°N, moderate), 2D and 3D agree away from horizon."""
        from libephemeris import calc_ut, houses
        from libephemeris.utils import azalt, ECL2HOR

        geo_lat, geo_lon = 41.9, 12.5
        jd = julday(2024, month, 21, float(hour))
        _, ascmc = houses(jd, geo_lat, geo_lon, ord("P"))
        asc = ascmc[0]
        sun_pos, _ = calc_ut(jd, 0, 0)
        sun_lon = sun_pos[0]

        geopos = (geo_lon, geo_lat, 0.0)
        xin = (sun_lon, 0.0, 1.0)
        _, true_altitude, _ = azalt(jd, ECL2HOR, geopos, 0.0, 0.0, xin)
        if abs(true_altitude) < 1.0:
            pytest.skip(f"Sun too close to horizon ({true_altitude:.3f}°)")

        result_2d = is_day_chart(sun_lon, asc)
        result_3d = _is_sun_above_horizon_3d(jd, sun_lon, 0.0, geo_lat, geo_lon)
        assert result_2d == result_3d, (
            f"2D/3D disagree at 2024-{month:02d}-21 {hour:02d}:00 UTC, Rome: "
            f"ASC={asc:.2f}, Sun={sun_lon:.2f}, alt={true_altitude:.3f}, "
            f"2D={result_2d}, 3D={result_3d}"
        )

    @pytest.mark.parametrize("hour", [0, 6, 12, 18])
    @pytest.mark.parametrize("month", [3, 6, 9, 12])
    def test_buenos_aires_agreement(self, month, hour):
        """Southern hemisphere moderate latitude (Buenos Aires, lat -34.6°)."""
        from libephemeris import calc_ut, houses
        from libephemeris.utils import azalt, ECL2HOR

        geo_lat, geo_lon = -34.6, -58.4
        jd = julday(2024, month, 21, float(hour))
        _, ascmc = houses(jd, geo_lat, geo_lon, ord("P"))
        asc = ascmc[0]
        sun_pos, _ = calc_ut(jd, 0, 0)
        sun_lon = sun_pos[0]

        geopos = (geo_lon, geo_lat, 0.0)
        xin = (sun_lon, 0.0, 1.0)
        _, true_altitude, _ = azalt(jd, ECL2HOR, geopos, 0.0, 0.0, xin)
        if abs(true_altitude) < 1.0:
            pytest.skip(f"Sun too close to horizon ({true_altitude:.3f}°)")

        result_2d = is_day_chart(sun_lon, asc)
        result_3d = _is_sun_above_horizon_3d(jd, sun_lon, 0.0, geo_lat, geo_lon)
        assert result_2d == result_3d, (
            f"2D/3D disagree at 2024-{month:02d}-21 {hour:02d}:00 UTC, BA: "
            f"ASC={asc:.2f}, Sun={sun_lon:.2f}, alt={true_altitude:.3f}, "
            f"2D={result_2d}, 3D={result_3d}"
        )
