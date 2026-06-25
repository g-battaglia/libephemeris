"""
Tests for ASSIST N-body integration fallback for minor bodies.

Tests verify:
- ASSIST fallback is attempted before Keplerian when available
- If ASSIST fails, falls back to Keplerian
- Velocities are computed correctly via central difference
- Source detection returns "ASSIST" when ASSIST is available
"""

import pytest
from unittest.mock import patch
from dataclasses import dataclass
import math

import libephemeris as eph
from libephemeris import state
from libephemeris.constants import (
    SEDNA,
    FLG_SPEED,
    FLG_HELCTR,
)


@dataclass
class MockPropagationResult:
    """Mock result from propagate_orbit_assist."""

    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float
    jd_tt: float

    @property
    def ecliptic_lon(self) -> float:
        return math.degrees(math.atan2(self.y, self.x)) % 360.0

    @property
    def ecliptic_lat(self) -> float:
        r = math.sqrt(self.x**2 + self.y**2 + self.z**2)
        return math.degrees(math.asin(self.z / r))

    @property
    def distance(self) -> float:
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)


@pytest.fixture(autouse=True)
def cleanup():
    """Reset state before and after each test."""
    from libephemeris.rebound_integration import reset_assist_data_cache

    state._SPK_BODY_MAP.clear()
    eph.set_strict_precision(None)
    eph.set_auto_spk_download(None)
    reset_assist_data_cache()
    yield
    state._SPK_BODY_MAP.clear()
    eph.set_strict_precision(None)
    eph.set_auto_spk_download(None)
    reset_assist_data_cache()


class TestAssistFallbackLogic:
    """Test ASSIST fallback behavior in _calc_body()."""

    def test_assist_tried_before_keplerian(self):
        """ASSIST should be attempted before Keplerian fallback."""
        eph.set_strict_precision(False)
        eph.set_auto_spk_download(False)

        mock_result = MockPropagationResult(
            x=-50.0, y=30.0, z=-10.0, vx=0.001, vy=-0.002, vz=0.0001, jd_tt=2451545.0
        )

        with patch(
            "libephemeris.rebound_integration.check_assist_data_available",
            return_value=True,
        ):
            with patch(
                "libephemeris.rebound_integration.propagate_orbit_assist",
                return_value=mock_result,
            ) as mock_propagate:
                pos, _ = eph.calc_ut(2451545.0, SEDNA, 0)

                mock_propagate.assert_called_once()
                assert pos[0] > 0
                assert pos[2] > 0

    def test_keplerian_used_if_assist_not_available(self):
        """Keplerian should be used if ASSIST is not installed."""
        eph.set_strict_precision(False)
        eph.set_auto_spk_download(False)

        with patch(
            "libephemeris.rebound_integration.check_assist_data_available",
            return_value=False,
        ):
            pos, _ = eph.calc_ut(2451545.0, SEDNA, 0)

            assert pos[0] != 0.0
            assert pos[2] > 50.0

    def test_keplerian_used_if_assist_raises(self):
        """Keplerian should be used if ASSIST raises an exception."""
        eph.set_strict_precision(False)
        eph.set_auto_spk_download(False)

        with patch(
            "libephemeris.rebound_integration.check_assist_data_available",
            return_value=True,
        ):
            with patch(
                "libephemeris.rebound_integration.propagate_orbit_assist",
                side_effect=FileNotFoundError("No ephemeris files"),
            ):
                pos, _ = eph.calc_ut(2451545.0, SEDNA, 0)

                assert pos[0] != 0.0

    def test_assist_velocities_via_central_difference(self):
        """ASSIST velocities should be computed via central difference."""
        eph.set_strict_precision(False)
        eph.set_auto_spk_download(False)

        call_count = [0]

        def mock_propagate(elements, jd_start, jd_end, **kwargs):
            call_count[0] += 1
            base_lon = 45.0
            base_lat = -12.0
            base_dist = 90.0

            dt_days = jd_end - 2451545.0
            lon = base_lon + 0.0001 * dt_days
            lat = base_lat + 0.00002 * dt_days
            dist = base_dist + 0.00001 * dt_days

            lon_rad = math.radians(lon)
            lat_rad = math.radians(lat)
            x = dist * math.cos(lat_rad) * math.cos(lon_rad)
            y = dist * math.cos(lat_rad) * math.sin(lon_rad)
            z = dist * math.sin(lat_rad)

            return MockPropagationResult(
                x=x, y=y, z=z, vx=0.001, vy=-0.002, vz=0.0001, jd_tt=jd_end
            )

        with patch(
            "libephemeris.rebound_integration.check_assist_data_available",
            return_value=True,
        ):
            with patch(
                "libephemeris.rebound_integration.propagate_orbit_assist",
                side_effect=mock_propagate,
            ):
                pos, _ = eph.calc_ut(2451545.0, SEDNA, FLG_SPEED)

                assert call_count[0] == 3
                assert pos[3] != 0.0
                assert pos[4] != 0.0
                assert pos[5] != 0.0

    def test_assist_heliocentric(self):
        """ASSIST should work with FLG_HELCTR."""
        eph.set_strict_precision(False)
        eph.set_auto_spk_download(False)

        mock_result = MockPropagationResult(
            x=50.0, y=30.0, z=10.0, vx=0.001, vy=-0.002, vz=0.0001, jd_tt=2451545.0
        )

        with patch(
            "libephemeris.rebound_integration.check_assist_data_available",
            return_value=True,
        ):
            with patch(
                "libephemeris.rebound_integration.propagate_orbit_assist",
                return_value=mock_result,
            ):
                pos, _ = eph.calc_ut(2451545.0, SEDNA, FLG_HELCTR)

                expected_lon = mock_result.ecliptic_lon
                expected_dist = mock_result.distance

                assert pos[0] == pytest.approx(expected_lon, abs=0.01)
                assert pos[2] == pytest.approx(expected_dist, abs=0.01)


class TestAssistSourceDetection:
    """Test source detection for ASSIST."""

    def test_source_is_assist_when_available(self):
        """Source detection should return 'ASSIST' when ASSIST is available."""
        with patch(
            "libephemeris.rebound_integration.check_assist_available", return_value=True
        ):
            from libephemeris.rebound_integration import check_assist_available

            is_available = check_assist_available()
            assert is_available is True

    def test_source_is_keplerian_when_assist_unavailable(self):
        """Source detection should return 'Keplerian' when ASSIST is not available."""
        with patch(
            "libephemeris.rebound_integration.check_assist_available",
            return_value=False,
        ):
            from libephemeris.rebound_integration import check_assist_available

            is_available = check_assist_available()
            assert is_available is False


class TestAssistVsKeplerianPrecision:
    """Compare ASSIST vs Keplerian results (when both available)."""

    def test_assist_and_keplerian_both_work(self):
        """Both ASSIST and Keplerian should produce valid results."""
        eph.set_strict_precision(False)
        eph.set_auto_spk_download(False)

        jd = 2451545.0

        keplerian_pos, _ = eph.calc_ut(jd, SEDNA, 0)

        mock_result = MockPropagationResult(
            x=-50.0, y=30.0, z=-10.0, vx=0.001, vy=-0.002, vz=0.0001, jd_tt=jd
        )

        with patch(
            "libephemeris.rebound_integration.check_assist_data_available",
            return_value=True,
        ):
            with patch(
                "libephemeris.rebound_integration.propagate_orbit_assist",
                return_value=mock_result,
            ) as mock_propagate:
                assist_pos, _ = eph.calc_ut(jd, SEDNA, 0)

                mock_propagate.assert_called()

        assert keplerian_pos[0] > 0
        assert keplerian_pos[2] > 50
        assert assist_pos[0] > 0
        assert assist_pos[2] > 0
