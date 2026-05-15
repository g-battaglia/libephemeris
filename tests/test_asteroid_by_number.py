"""
Tests for generic asteroid lookup by number.

Tests the calc_asteroid_by_number function which allows calculating
positions for any numbered asteroid by fetching orbital elements from
JPL SBDB API on demand.
"""

import math
import pytest
from unittest.mock import patch, MagicMock

from libephemeris.minor_bodies import (
    calc_asteroid_by_number,
    fetch_orbital_elements_from_sbdb,
    clear_asteroid_elements_cache,
    OrbitalElements,
    _ASTEROID_ELEMENTS_CACHE,
)
from libephemeris.constants import AST_OFFSET


# Sample SBDB API response for Eros (433)
MOCK_SBDB_RESPONSE_EROS = {
    "object": {
        "fullname": "433 Eros (A898 PA)",
    },
    "orbit": {
        "epoch": 2461000.5,
        "elements": [
            {"name": "e", "value": "0.2228359407071628"},
            {"name": "a", "value": "1.458120998474684"},
            {"name": "i", "value": "10.82846651399785"},
            {"name": "om", "value": "304.2701025753316"},
            {"name": "w", "value": "178.9297536744151"},
            {"name": "ma", "value": "310.5543277370992"},
            {"name": "n", "value": "0.5597752949285997"},
        ],
    },
}

# Sample SBDB API response for a fictional small asteroid
MOCK_SBDB_RESPONSE_SMALL = {
    "object": {
        "fullname": "12345 TestAsteroid",
    },
    "orbit": {
        "epoch": 2461000.5,
        "elements": [
            {"name": "e", "value": "0.15"},
            {"name": "a", "value": "2.5"},
            {"name": "i", "value": "5.0"},
            {"name": "om", "value": "100.0"},
            {"name": "w", "value": "200.0"},
            {"name": "ma", "value": "50.0"},
            {"name": "n", "value": "0.25"},
        ],
    },
}


class TestFetchOrbitalElementsFromSBDB:
    """Tests for fetch_orbital_elements_from_sbdb function."""

    def setup_method(self):
        """Clear cache before each test."""
        clear_asteroid_elements_cache()

    def test_fetch_returns_orbital_elements(self):
        """Test that fetching returns valid OrbitalElements object."""
        import json
        import urllib.request

        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps(MOCK_SBDB_RESPONSE_EROS).encode(
            "utf-8"
        )
        mock_response.__enter__ = lambda x: mock_response
        mock_response.__exit__ = lambda *args: None

        with patch.object(urllib.request, "urlopen", return_value=mock_response):
            elements = fetch_orbital_elements_from_sbdb(433)

        assert elements is not None
        assert isinstance(elements, OrbitalElements)
        assert elements.name == "433 Eros (A898 PA)"
        assert abs(elements.a - 1.458120998474684) < 1e-10
        assert abs(elements.e - 0.2228359407071628) < 1e-10
        assert abs(elements.i - 10.82846651399785) < 1e-10

    def test_fetch_caches_result(self):
        """Test that results are cached to avoid redundant API calls."""
        import json
        import urllib.request

        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps(MOCK_SBDB_RESPONSE_EROS).encode(
            "utf-8"
        )
        mock_response.__enter__ = lambda x: mock_response
        mock_response.__exit__ = lambda *args: None

        with patch.object(
            urllib.request, "urlopen", return_value=mock_response
        ) as mock_urlopen:
            # First call
            elements1 = fetch_orbital_elements_from_sbdb(433)
            # Second call should use cache
            elements2 = fetch_orbital_elements_from_sbdb(433)

        # urlopen should only be called once
        assert mock_urlopen.call_count == 1
        assert elements1 is elements2

    def test_fetch_returns_none_on_error(self):
        """Test that errors return None instead of raising."""
        import urllib.request
        import urllib.error

        with patch.object(
            urllib.request,
            "urlopen",
            side_effect=urllib.error.URLError("Network error"),
        ):
            elements = fetch_orbital_elements_from_sbdb(99999999)

        assert elements is None

    def test_fetch_returns_none_for_api_error(self):
        """Test that API error responses return None."""
        import json
        import urllib.request

        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps(
            {"error": "specified object was not found"}
        ).encode("utf-8")
        mock_response.__enter__ = lambda x: mock_response
        mock_response.__exit__ = lambda *args: None

        with patch.object(urllib.request, "urlopen", return_value=mock_response):
            elements = fetch_orbital_elements_from_sbdb(99999999)

        assert elements is None


class TestClearAsteroidElementsCache:
    """Tests for clear_asteroid_elements_cache function."""

    def test_clear_cache_returns_count(self):
        """Test that clearing cache returns the number of entries cleared."""
        import json
        import urllib.request

        # Add some entries to cache
        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps(MOCK_SBDB_RESPONSE_EROS).encode(
            "utf-8"
        )
        mock_response.__enter__ = lambda x: mock_response
        mock_response.__exit__ = lambda *args: None

        with patch.object(urllib.request, "urlopen", return_value=mock_response):
            fetch_orbital_elements_from_sbdb(433)

        # Cache should have 1 entry
        assert len(_ASTEROID_ELEMENTS_CACHE) >= 1

        cleared = clear_asteroid_elements_cache()
        assert cleared >= 1
        assert len(_ASTEROID_ELEMENTS_CACHE) == 0

    def test_clear_empty_cache_returns_zero(self):
        """Test that clearing empty cache returns 0."""
        clear_asteroid_elements_cache()  # Ensure empty
        cleared = clear_asteroid_elements_cache()
        assert cleared == 0


class TestCalcAsteroidByNumber:
    """Tests for calc_asteroid_by_number function."""

    def setup_method(self):
        """Clear cache before each test."""
        clear_asteroid_elements_cache()

    def test_calc_known_asteroid_uses_predefined_elements(self):
        """Test that known asteroids use predefined MINOR_BODY_ELEMENTS."""
        # Ceres (1) is in MINOR_BODY_ELEMENTS
        lon, lat, dist = calc_asteroid_by_number(1, 2451545.0)

        # Should return valid coordinates
        assert 0 <= lon < 360
        assert -90 <= lat <= 90
        assert dist > 0

    def test_calc_special_body_chiron(self):
        """Test that Chiron (2060) uses its predefined elements."""
        lon, lat, dist = calc_asteroid_by_number(2060, 2451545.0)

        # Should return valid coordinates
        assert 0 <= lon < 360
        assert -90 <= lat <= 90
        assert dist > 0

    def test_calc_unknown_asteroid_fetches_from_sbdb(self):
        """Test that unknown asteroids are fetched from SBDB."""
        import json
        import urllib.request

        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps(MOCK_SBDB_RESPONSE_SMALL).encode(
            "utf-8"
        )
        mock_response.__enter__ = lambda x: mock_response
        mock_response.__exit__ = lambda *args: None

        with patch.object(urllib.request, "urlopen", return_value=mock_response):
            lon, lat, dist = calc_asteroid_by_number(12345, 2451545.0)

        # Should return valid coordinates
        assert 0 <= lon < 360
        assert -90 <= lat <= 90
        assert dist > 0

    def test_calc_invalid_asteroid_number_raises(self):
        """Test that invalid asteroid numbers raise ValueError."""
        with pytest.raises(ValueError, match="positive integer"):
            calc_asteroid_by_number(0, 2451545.0)

        with pytest.raises(ValueError, match="positive integer"):
            calc_asteroid_by_number(-1, 2451545.0)

    def test_calc_not_found_asteroid_raises(self):
        """Test that asteroids not found in SBDB raise ValueError."""
        import urllib.request
        import urllib.error

        with patch.object(
            urllib.request,
            "urlopen",
            side_effect=urllib.error.URLError("Not found"),
        ):
            with pytest.raises(ValueError, match="not found in JPL"):
                calc_asteroid_by_number(99999999, 2451545.0)

    def test_calc_returns_tuple_of_three_floats(self):
        """Test that result is a tuple of three floats."""
        lon, lat, dist = calc_asteroid_by_number(1, 2451545.0)

        assert isinstance(lon, float)
        assert isinstance(lat, float)
        assert isinstance(dist, float)

    def test_calc_with_use_spk_false(self):
        """Test that use_spk=False forces Keplerian calculation."""
        # Even for known asteroids, should still work
        lon, lat, dist = calc_asteroid_by_number(1, 2451545.0, use_spk=False)

        assert 0 <= lon < 360
        assert -90 <= lat <= 90
        assert dist > 0


class TestCalcAsteroidConsistency:
    """Tests for consistency between calc methods."""

    def test_asteroid_position_consistency_ceres(self):
        """Test that calc_asteroid_by_number gives same result as calc_minor_body_heliocentric for Ceres."""
        from libephemeris.minor_bodies import calc_minor_body_heliocentric
        from libephemeris.constants import CERES

        jd = 2451545.0  # J2000.0

        lon1, lat1, dist1 = calc_asteroid_by_number(1, jd)
        lon2, lat2, dist2 = calc_minor_body_heliocentric(CERES, jd)

        # Should be identical
        assert abs(lon1 - lon2) < 1e-10
        assert abs(lat1 - lat2) < 1e-10
        assert abs(dist1 - dist2) < 1e-10

    def test_asteroid_position_consistency_eros(self):
        """Test consistency for Eros (433) which is in MINOR_BODY_ELEMENTS."""
        from libephemeris.minor_bodies import calc_minor_body_heliocentric
        from libephemeris.constants import EROS

        jd = 2451545.0

        lon1, lat1, dist1 = calc_asteroid_by_number(433, jd)
        lon2, lat2, dist2 = calc_minor_body_heliocentric(EROS, jd)

        # Should be identical
        assert abs(lon1 - lon2) < 1e-10
        assert abs(lat1 - lat2) < 1e-10
        assert abs(dist1 - dist2) < 1e-10


class TestCalcAsteroidPrecision:
    """Tests for precision of asteroid position calculations."""

    def test_main_belt_asteroid_reasonable_position(self):
        """Test that main belt asteroid positions are reasonable."""
        # Using Ceres which is well known
        lon, lat, dist = calc_asteroid_by_number(1, 2451545.0)

        # Ceres should be in main belt (distance ~2-3 AU from Sun)
        assert 2.0 < dist < 4.0

        # Latitude should be moderate (Ceres has ~10 deg inclination)
        assert abs(lat) < 20

    def test_near_earth_asteroid_distance(self):
        """Test that near-Earth asteroids have appropriate distances."""
        # Eros comes close to Earth
        lon, lat, dist = calc_asteroid_by_number(433, 2451545.0)

        # Eros should be within inner solar system
        assert 0.5 < dist < 3.0


# Integration test that requires network (skipped by default)
@pytest.mark.skipif(
    True,  # Skip by default since it requires network
    reason="Network test - requires internet access to JPL SBDB",
)
class TestCalcAsteroidNetworkIntegration:
    """Integration tests that actually call JPL SBDB API."""

    def test_fetch_real_asteroid_from_sbdb(self):
        """Test fetching a real asteroid from JPL SBDB."""
        clear_asteroid_elements_cache()

        # Fetch Eros (433) - well-known asteroid
        elements = fetch_orbital_elements_from_sbdb(433)

        assert elements is not None
        assert "Eros" in elements.name
        # Eros semi-major axis is about 1.458 AU
        assert 1.4 < elements.a < 1.6

    def test_calc_real_asteroid_position(self):
        """Test calculating position for a real asteroid."""
        clear_asteroid_elements_cache()

        # This would actually fetch from SBDB for a non-cached asteroid
        # Using Vesta (4) which is in predefined elements
        lon, lat, dist = calc_asteroid_by_number(4, 2451545.0)

        assert 0 <= lon < 360
        assert -90 <= lat <= 90
        assert dist > 0
