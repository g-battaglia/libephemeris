"""
Unit tests for asteroid name lookup functionality.

Tests verify:
- Local lookup for known asteroids in MINOR_BODY_ELEMENTS
- Case-insensitive name matching
- Cache clearing functionality
- Integration with JPL SBDB API (network tests marked)
"""

import pytest
from unittest.mock import patch, MagicMock
import json

from libephemeris.minor_bodies import (
    get_asteroid_number,
    clear_asteroid_name_cache,
)


@pytest.fixture(autouse=True)
def clear_cache():
    """Clear the name cache before and after each test."""
    clear_asteroid_name_cache()
    yield
    clear_asteroid_name_cache()


@pytest.mark.unit
class TestLocalAsteroidNameLookup:
    """Tests for local asteroid name lookup (no network required)."""

    def test_lookup_ceres(self):
        """Ceres (1) should be found locally."""
        result = get_asteroid_number("Ceres")
        assert result == 1

    def test_lookup_chiron(self):
        """Chiron (2060) should be found locally."""
        result = get_asteroid_number("Chiron")
        assert result == 2060

    def test_lookup_pholus(self):
        """Pholus (5145) should be found locally."""
        result = get_asteroid_number("Pholus")
        assert result == 5145

    def test_lookup_pallas(self):
        """Pallas (2) should be found locally."""
        result = get_asteroid_number("Pallas")
        assert result == 2

    def test_lookup_juno(self):
        """Juno (3) should be found locally."""
        result = get_asteroid_number("Juno")
        assert result == 3

    def test_lookup_vesta(self):
        """Vesta (4) should be found locally."""
        result = get_asteroid_number("Vesta")
        assert result == 4

    def test_lookup_eros(self):
        """Eros (433) should be found locally."""
        result = get_asteroid_number("Eros")
        assert result == 433

    def test_lookup_apophis(self):
        """Apophis (99942) should be found locally."""
        result = get_asteroid_number("Apophis")
        assert result == 99942

    def test_lookup_eris(self):
        """Eris (136199) should be found locally."""
        result = get_asteroid_number("Eris")
        assert result == 136199

    def test_lookup_sedna(self):
        """Sedna (90377) should be found locally."""
        result = get_asteroid_number("Sedna")
        assert result == 90377

    def test_lookup_haumea(self):
        """Haumea (136108) should be found locally."""
        result = get_asteroid_number("Haumea")
        assert result == 136108

    def test_lookup_makemake(self):
        """Makemake (136472) should be found locally."""
        result = get_asteroid_number("Makemake")
        assert result == 136472

    def test_lookup_ixion(self):
        """Ixion (28978) should be found locally."""
        result = get_asteroid_number("Ixion")
        assert result == 28978

    def test_lookup_orcus(self):
        """Orcus (90482) should be found locally."""
        result = get_asteroid_number("Orcus")
        assert result == 90482

    def test_lookup_quaoar(self):
        """Quaoar (50000) should be found locally."""
        result = get_asteroid_number("Quaoar")
        assert result == 50000

    def test_lookup_nessus(self):
        """Nessus (7066) should be found locally."""
        result = get_asteroid_number("Nessus")
        assert result == 7066

    def test_lookup_asbolus(self):
        """Asbolus (8405) should be found locally."""
        result = get_asteroid_number("Asbolus")
        assert result == 8405

    def test_lookup_chariklo(self):
        """Chariklo (10199) should be found locally."""
        result = get_asteroid_number("Chariklo")
        assert result == 10199

    def test_lookup_gonggong(self):
        """Gonggong (225088) should be found locally."""
        result = get_asteroid_number("Gonggong")
        assert result == 225088

    def test_lookup_varuna(self):
        """Varuna (20000) should be found locally."""
        result = get_asteroid_number("Varuna")
        assert result == 20000

    def test_lookup_bennu(self):
        """Bennu (101955) should be found locally."""
        result = get_asteroid_number("Bennu")
        assert result == 101955

    def test_lookup_ryugu(self):
        """Ryugu (162173) should be found locally."""
        result = get_asteroid_number("Ryugu")
        assert result == 162173

    def test_lookup_hygiea(self):
        """Hygiea (10) should be found locally."""
        result = get_asteroid_number("Hygiea")
        assert result == 10

    def test_lookup_psyche(self):
        """Psyche (16) should be found locally."""
        result = get_asteroid_number("Psyche")
        assert result == 16

    def test_lookup_itokawa(self):
        """Itokawa (25143) should be found locally."""
        result = get_asteroid_number("Itokawa")
        assert result == 25143

    def test_lookup_toutatis(self):
        """Toutatis (4179) should be found locally."""
        result = get_asteroid_number("Toutatis")
        assert result == 4179

    def test_lookup_hidalgo(self):
        """Hidalgo (944) should be found locally."""
        result = get_asteroid_number("Hidalgo")
        assert result == 944

    def test_lookup_sappho(self):
        """Sappho (80) should be found locally."""
        result = get_asteroid_number("Sappho")
        assert result == 80


@pytest.mark.unit
class TestCaseInsensitiveLookup:
    """Tests for case-insensitive name matching."""

    def test_lowercase(self):
        """Lowercase name should work."""
        result = get_asteroid_number("ceres")
        assert result == 1

    def test_uppercase(self):
        """Uppercase name should work."""
        result = get_asteroid_number("CERES")
        assert result == 1

    def test_mixed_case(self):
        """Mixed case name should work."""
        result = get_asteroid_number("CeReS")
        assert result == 1

    def test_with_whitespace(self):
        """Name with leading/trailing whitespace should work."""
        result = get_asteroid_number("  Ceres  ")
        assert result == 1

    def test_chiron_lowercase(self):
        """Chiron in lowercase should work."""
        result = get_asteroid_number("chiron")
        assert result == 2060


@pytest.mark.unit
class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_empty_string(self):
        """Empty string should return None."""
        result = get_asteroid_number("")
        assert result is None

    def test_whitespace_only(self):
        """Whitespace-only string should return None."""
        result = get_asteroid_number("   ")
        assert result is None


@pytest.mark.unit
class TestCacheClearing:
    """Tests for cache clearing functionality."""

    def test_clear_empty_cache(self):
        """Clearing empty cache should return 0."""
        clear_asteroid_name_cache()  # Ensure empty
        result = clear_asteroid_name_cache()
        assert result == 0

    def test_clear_populated_cache(self):
        """Clearing populated cache should return count."""
        # Trigger local cache population
        get_asteroid_number("Ceres")
        get_asteroid_number("Vesta")

        # The local cache is populated on first lookup
        result = clear_asteroid_name_cache()
        # Should have at least the built-in entries
        assert result > 0


@pytest.mark.unit
class TestMockedSBDBLookup:
    """Tests for JPL SBDB API lookup with mocked responses."""

    @patch("urllib.request.urlopen")
    def test_sbdb_lookup_success(self, mock_urlopen):
        """Test successful SBDB API lookup for unknown asteroid."""
        # Mock response for asteroid Flora (8)
        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps(
            {
                "object": {
                    "spkid": "2000008",
                    "fullname": "8 Flora",
                    "des": "8",
                },
                "orbit": {
                    "epoch": "2461000.5",
                    "elements": [],
                },
            }
        ).encode("utf-8")
        mock_response.__enter__ = lambda s: s
        mock_response.__exit__ = lambda s, *args: None
        mock_urlopen.return_value = mock_response

        # Clear cache first
        clear_asteroid_name_cache()

        # Use a name not in local cache
        result = get_asteroid_number("Flora")
        assert result == 8

    @patch("urllib.request.urlopen")
    def test_sbdb_lookup_not_found(self, mock_urlopen):
        """Test SBDB API lookup for non-existent asteroid."""
        # Mock error response
        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps(
            {
                "error": "specified object was not found",
            }
        ).encode("utf-8")
        mock_response.__enter__ = lambda s: s
        mock_response.__exit__ = lambda s, *args: None
        mock_urlopen.return_value = mock_response

        # Clear cache first
        clear_asteroid_name_cache()

        result = get_asteroid_number("NonExistentAsteroid12345")
        assert result is None

    @patch("urllib.request.urlopen")
    def test_sbdb_lookup_with_fullname_parsing(self, mock_urlopen):
        """Test SBDB API lookup parsing fullname format."""
        # Mock response with number in fullname
        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps(
            {
                "object": {
                    "fullname": "15 Eunomia",
                },
                "orbit": {
                    "epoch": "2461000.5",
                    "elements": [],
                },
            }
        ).encode("utf-8")
        mock_response.__enter__ = lambda s: s
        mock_response.__exit__ = lambda s, *args: None
        mock_urlopen.return_value = mock_response

        # Clear cache first
        clear_asteroid_name_cache()

        result = get_asteroid_number("Eunomia")
        assert result == 15

    @patch("urllib.request.urlopen")
    def test_sbdb_network_error(self, mock_urlopen):
        """Test handling of network errors."""
        import urllib.error

        mock_urlopen.side_effect = urllib.error.URLError("Network error")

        # Clear cache first
        clear_asteroid_name_cache()

        result = get_asteroid_number("UnknownAsteroid")
        assert result is None


@pytest.mark.unit
class TestLocalCachePersistence:
    """Tests for local cache persistence across lookups."""

    def test_cached_result_returned(self):
        """Verify cached result is returned on second lookup."""
        # First lookup
        result1 = get_asteroid_number("Ceres")
        # Second lookup (should use cache)
        result2 = get_asteroid_number("ceres")

        assert result1 == result2 == 1

    def test_different_asteroids_cached(self):
        """Verify different asteroids are cached correctly."""
        result_ceres = get_asteroid_number("Ceres")
        result_vesta = get_asteroid_number("Vesta")
        result_chiron = get_asteroid_number("Chiron")

        assert result_ceres == 1
        assert result_vesta == 4
        assert result_chiron == 2060


@pytest.mark.unit
class TestIntegrationWithCalcAsteroidByNumber:
    """Tests for integration between name lookup and position calculation."""

    def test_lookup_and_calculate_position(self):
        """Verify we can look up a name and calculate its position."""
        from libephemeris.minor_bodies import calc_asteroid_by_number

        # Look up Ceres
        asteroid_num = get_asteroid_number("Ceres")
        assert asteroid_num == 1

        # Calculate position at J2000.0
        jd = 2451545.0
        lon, lat, dist = calc_asteroid_by_number(asteroid_num, jd)

        # Verify we get valid results
        assert 0 <= lon < 360
        assert -90 <= lat <= 90
        assert dist > 0


@pytest.mark.unit
class TestAllKnownAsteroids:
    """Test all known asteroids in the local database."""

    @pytest.mark.parametrize(
        "name,expected_number",
        [
            ("Chiron", 2060),
            ("Pholus", 5145),
            ("Ceres", 1),
            ("Pallas", 2),
            ("Juno", 3),
            ("Vesta", 4),
            ("Eris", 136199),
            ("Sedna", 90377),
            ("Haumea", 136108),
            ("Makemake", 136472),
            ("Ixion", 28978),
            ("Orcus", 90482),
            ("Quaoar", 50000),
            ("Nessus", 7066),
            ("Asbolus", 8405),
            ("Chariklo", 10199),
            ("Gonggong", 225088),
            ("Varuna", 20000),
            ("Apophis", 99942),
            ("Hygiea", 10),
            ("Interamnia", 704),
            ("Davida", 511),
            ("Europa", 52),
            ("Sylvia", 87),
            ("Psyche", 16),
            ("Eros", 433),
            ("Amor", 1221),
            ("Icarus", 1566),
            ("Toro", 1685),
            ("Sappho", 80),
            ("Pandora", 55),
            ("Lilith", 1181),
            ("Hidalgo", 944),
            ("Toutatis", 4179),
            ("Itokawa", 25143),
            ("Bennu", 101955),
            ("Ryugu", 162173),
        ],
    )
    def test_known_asteroid_lookup(self, name, expected_number):
        """Verify all known asteroids return correct catalog numbers."""
        result = get_asteroid_number(name)
        assert result == expected_number, (
            f"{name}: expected {expected_number}, got {result}"
        )


@pytest.mark.network
class TestNetworkSBDBLookup:
    """
    Tests that require network access to JPL SBDB.

    These tests are marked with @pytest.mark.network and should only be run
    when network access is available. Run with: pytest -m network
    """

    @pytest.mark.slow
    def test_real_sbdb_lookup_flora(self):
        """Test real SBDB API lookup for Flora (8)."""
        # Clear cache to force network lookup
        clear_asteroid_name_cache()

        result = get_asteroid_number("Flora")
        assert result == 8

    @pytest.mark.slow
    def test_real_sbdb_lookup_eunomia(self):
        """Test real SBDB API lookup for Eunomia (15)."""
        # Clear cache to force network lookup
        clear_asteroid_name_cache()

        result = get_asteroid_number("Eunomia")
        assert result == 15
