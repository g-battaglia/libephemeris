"""
Unit tests for retrograde station handling functions.

Tests the stability of calculations when planets are near stationary points
with near-zero velocity.
"""

import pytest
import libephemeris as ephem
from libephemeris.crossing import (
    is_retrograde,
    get_station_type,
    find_station_ut,
    next_retrograde_ut,
    calc_velocity_at_station,
    get_station_info,
    _is_near_station,
    STATION_SPEED_THRESHOLD,
)
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    FLG_SPEED,
)

pytestmark = pytest.mark.slow


@pytest.mark.unit
class TestIsNearStation:
    """Tests for _is_near_station helper function."""

    def test_near_station_threshold(self):
        """Test station detection threshold."""
        # Speed below threshold is near station
        assert _is_near_station(0.0001) is True
        assert _is_near_station(-0.0001) is True
        assert _is_near_station(0.0) is True

        # Speed above threshold is not near station
        assert _is_near_station(0.01) is False
        assert _is_near_station(-0.01) is False
        assert _is_near_station(1.0) is False

    def test_threshold_boundary(self):
        """Test exact threshold value."""
        # Exactly at threshold
        assert _is_near_station(STATION_SPEED_THRESHOLD) is False
        # Just below threshold
        assert _is_near_station(STATION_SPEED_THRESHOLD * 0.99) is True


@pytest.mark.unit
class TestIsRetrograde:
    """Tests for is_retrograde function."""

    def test_sun_never_retrograde(self):
        """Sun should never appear retrograde."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        assert not is_retrograde(SUN, jd)

    def test_moon_never_retrograde(self):
        """Moon should never appear retrograde."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        assert not is_retrograde(MOON, jd)

    def test_mercury_retrograde_detection(self):
        """Test Mercury retrograde detection at known retrograde period."""
        # Mercury was retrograde around April 1-25, 2024
        # During retrograde, velocity should be negative
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Find a station
        jd_station, stype = find_station_ut(MERCURY, jd_start)

        # Check retrograde status a few days after station
        if stype == "SR":
            # After SR station, should be retrograde
            jd_check = jd_station + 10.0
            assert is_retrograde(MERCURY, jd_check)
        else:
            # After SD station, should not be retrograde
            jd_check = jd_station + 10.0
            assert not is_retrograde(MERCURY, jd_check)

    def test_direct_motion_detection(self):
        """Test that direct motion is properly detected."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        # Find next direct station
        jd_sd, _ = find_station_ut(MARS, jd, "direct")

        # 30 days after SD, should be in direct motion
        jd_direct = jd_sd + 30.0
        assert not is_retrograde(MARS, jd_direct)


@pytest.mark.unit
class TestGetStationType:
    """Tests for get_station_type function."""

    def test_sun_always_direct(self):
        """Sun should always return 'direct'."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        assert get_station_type(SUN, jd) == "direct"

    def test_moon_always_direct(self):
        """Moon should always return 'direct'."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        assert get_station_type(MOON, jd) == "direct"

    def test_retrograde_type_detection(self):
        """Test retrograde state detection."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        # Find SR station
        jd_sr, _ = find_station_ut(MERCURY, jd, "retrograde")

        # During retrograde (10 days after SR)
        station_type = get_station_type(MERCURY, jd_sr + 10.0)
        assert station_type == "retrograde"

    def test_direct_type_detection(self):
        """Test direct state detection."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        # Find SD station
        jd_sd, _ = find_station_ut(MERCURY, jd, "direct")

        # During direct motion (20 days after SD)
        station_type = get_station_type(MERCURY, jd_sd + 20.0)
        assert station_type == "direct"


@pytest.mark.unit
class TestSweFindStationUt:
    """Tests for find_station_ut function."""

    def test_sun_raises_error(self):
        """Sun should raise ValueError."""
        jd = ephem.julday(2024, 1, 1, 0.0)
        with pytest.raises(ValueError, match="Sun and Moon"):
            find_station_ut(SUN, jd)

    def test_moon_raises_error(self):
        """Moon should raise ValueError."""
        jd = ephem.julday(2024, 1, 1, 0.0)
        with pytest.raises(ValueError, match="Sun and Moon"):
            find_station_ut(MOON, jd)

    def test_find_any_station_mercury(self):
        """Find any Mercury station."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        jd_station, stype = find_station_ut(MERCURY, jd)

        # Should find a station
        assert jd_station > jd
        assert stype in ("SR", "SD")

        # Station should be within Mercury's synodic period (~116 days)
        assert jd_station - jd < 116

    def test_find_retrograde_station(self):
        """Find next retrograde (SR) station."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        jd_sr, stype = find_station_ut(MERCURY, jd, "retrograde")

        assert stype == "SR"
        assert jd_sr > jd

    def test_find_direct_station(self):
        """Find next direct (SD) station."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        jd_sd, stype = find_station_ut(MERCURY, jd, "direct")

        assert stype == "SD"
        assert jd_sd > jd

    def test_station_velocity_near_zero(self):
        """At station, velocity should be near zero."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        jd_station, _ = find_station_ut(MARS, jd)

        # Check velocity at station
        pos, _ = ephem.calc_ut(jd_station, MARS, FLG_SPEED)
        velocity = pos[3]

        # Velocity should be very small at station
        assert abs(velocity) < 0.01, f"Velocity at station: {velocity}"

    @pytest.mark.parametrize("planet", [MERCURY, VENUS, MARS, JUPITER, SATURN])
    def test_find_station_various_planets(self, planet):
        """Test station finding for various planets."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        jd_station, stype = find_station_ut(planet, jd)

        assert jd_station > jd
        assert stype in ("SR", "SD")


@pytest.mark.unit
class TestSweNextRetrogradeUt:
    """Tests for next_retrograde_ut function."""

    def test_find_next_retrograde_period(self):
        """Find next retrograde period."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        jd_sr, jd_sd = next_retrograde_ut(MERCURY, jd)

        # SR should come before SD
        assert jd_sr < jd_sd
        # Both should be in the future
        assert jd_sr > jd

        # Retrograde period should be reasonable (Mercury ~3 weeks)
        retrograde_duration = jd_sd - jd_sr
        assert 15 < retrograde_duration < 30, f"Duration: {retrograde_duration}"

    def test_retrograde_period_mars(self):
        """Test Mars retrograde period."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        jd_sr, jd_sd = next_retrograde_ut(MARS, jd)

        # Mars retrograde lasts ~60-80 days
        retrograde_duration = jd_sd - jd_sr
        assert 50 < retrograde_duration < 90, f"Mars Rx duration: {retrograde_duration}"


@pytest.mark.unit
class TestCalcVelocityAtStation:
    """Tests for calc_velocity_at_station function."""

    def test_velocity_near_zero_at_station(self):
        """Velocity should be near zero at station point."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        # Find a station
        jd_station, _ = find_station_ut(MERCURY, jd)

        # Calculate velocity with wider timestep
        v_lon, v_lat, v_dist = calc_velocity_at_station(MERCURY, jd_station)

        # Longitude velocity should be very small
        assert abs(v_lon) < 0.1, f"Velocity at station: {v_lon}"

    def test_velocity_comparison_with_swe_calc(self):
        """Compare with calc_ut velocity."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        jd_station, _ = find_station_ut(MARS, jd)

        # Standard velocity
        pos, _ = ephem.calc_ut(jd_station, MARS, FLG_SPEED)
        v_standard = pos[3]

        # Stable velocity
        v_stable, _, _ = calc_velocity_at_station(MARS, jd_station)

        # Both should be near zero at station
        assert abs(v_standard) < 0.01
        assert abs(v_stable) < 0.01

    def test_velocity_away_from_station(self):
        """Velocity should be non-zero away from station."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        jd_station, stype = find_station_ut(MERCURY, jd)

        # 30 days after station
        jd_away = jd_station + 30.0
        v_lon, _, _ = calc_velocity_at_station(MERCURY, jd_away)

        # Should have noticeable velocity away from station
        assert abs(v_lon) > 0.5


@pytest.mark.unit
class TestGetStationInfo:
    """Tests for get_station_info function."""

    def test_sun_raises_error(self):
        """Sun should raise ValueError."""
        jd = ephem.julday(2024, 1, 1, 0.0)
        with pytest.raises(ValueError, match="Sun and Moon"):
            get_station_info(SUN, jd)

    def test_station_info_structure(self):
        """Test returned dictionary structure."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        info = get_station_info(MERCURY, jd)

        # Check required keys
        assert "jd_station" in info
        assert "station_type" in info
        assert "days_to_station" in info
        assert "longitude_at_station" in info
        assert "is_currently_retrograde" in info
        assert "velocity" in info

    def test_station_info_values(self):
        """Test station info values are reasonable."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        info = get_station_info(MERCURY, jd)

        # Station should be in the future
        assert info["days_to_station"] > 0

        # Station type should be valid
        assert info["station_type"] in ("SR", "SD")

        # Longitude should be valid
        assert 0 <= info["longitude_at_station"] < 360

        # Velocity should be a number
        assert isinstance(info["velocity"], float)


@pytest.mark.integration
class TestRetrogradeStationStability:
    """Integration tests for stability near retrograde stations."""

    def test_crossing_near_station(self):
        """Test that crossing calculations work near stations."""
        from libephemeris.crossing import cross_ut

        jd = ephem.julday(2024, 1, 1, 0.0)

        # Find a Mercury station
        jd_station, _ = find_station_ut(MERCURY, jd)

        # Get Mercury's position at station
        pos, _ = ephem.calc_ut(jd_station, MERCURY, FLG_SPEED)
        lon_at_station = pos[0]

        # Try to find crossing near station (this tests the Brent's method fallback)
        # Search for crossing just after station
        target_lon = (lon_at_station + 5) % 360
        jd_cross = cross_ut(MERCURY, target_lon, jd_station - 5)

        # Should find a valid crossing
        assert isinstance(jd_cross, float)

    def test_position_stability_at_station(self):
        """Test that position calculations remain stable at station."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        jd_station, _ = find_station_ut(MARS, jd)

        # Calculate positions at very small intervals around station
        positions = []
        for offset in [-0.001, -0.0001, 0, 0.0001, 0.001]:
            pos, _ = ephem.calc_ut(jd_station + offset, MARS, FLG_SPEED)
            positions.append(pos[0])

        # Position changes should be very small (< 1 arcsecond = 1/3600 deg)
        for i in range(len(positions) - 1):
            diff = abs(positions[i + 1] - positions[i])
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.01, f"Position jump at station: {diff} degrees"

    def test_velocity_sign_change_at_station(self):
        """Verify velocity sign change occurs at station."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        jd_sr, stype = find_station_ut(MERCURY, jd, "retrograde")

        # Get velocities before and after station
        pos_before, _ = ephem.calc_ut(jd_sr - 1.0, MERCURY, FLG_SPEED)
        pos_after, _ = ephem.calc_ut(jd_sr + 1.0, MERCURY, FLG_SPEED)

        # Velocity should change sign (positive before, negative after for SR)
        assert pos_before[3] > 0, f"Before SR: velocity={pos_before[3]}"
        assert pos_after[3] < 0, f"After SR: velocity={pos_after[3]}"

    def test_numerical_stability_slow_planet(self):
        """Test stability for slow outer planets near station."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        # Saturn is a slow planet
        jd_station, _ = find_station_ut(SATURN, jd)

        # Calculate velocity at station using stable method
        v_lon, v_lat, v_dist = calc_velocity_at_station(SATURN, jd_station)

        # Velocity should be extremely small
        assert abs(v_lon) < 0.005, f"Saturn velocity at station: {v_lon}"

        # Standard calculation should also work
        pos, _ = ephem.calc_ut(jd_station, SATURN, FLG_SPEED)
        assert abs(pos[3]) < 0.005, f"Saturn calc velocity: {pos[3]}"
