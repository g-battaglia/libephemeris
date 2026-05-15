"""
Unit tests for Crux (Southern Cross) constellation stars.

The Southern Cross is one of the most famous asterisms in the southern sky,
consisting of 4 bright stars that form a distinctive cross shape:

- Acrux (Alpha Crucis) - brightest star, southernmost point of the cross
- Mimosa/Becrux (Beta Crucis) - second brightest, eastern arm
- Gacrux (Gamma Crucis) - red giant, northern point of the cross
- Delta Crucis - western arm, completes the cross shape

The Southern Cross is used for navigation in the southern hemisphere,
as it helps locate the South Celestial Pole. The two brightest stars
(Acrux and Gacrux) point toward the pole.
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    ACRUX,
    MIMOSA,
    GACRUX,
    DELTA_CRUCIS,
)
from libephemeris.fixed_stars import (
    STAR_CATALOG,
    resolve_star_name,
    get_canonical_star_name,
)


# The 4 Southern Cross stars
# (constant, name, hip_number, magnitude)
CRUX_STARS = [
    (ACRUX, "Acrux", 60718, 0.76),  # Alpha Cru - brightest, south point
    (MIMOSA, "Mimosa", 62434, 1.25),  # Beta Cru - east arm
    (GACRUX, "Gacrux", 61084, 1.64),  # Gamma Cru - north point
    (DELTA_CRUCIS, "Delta Crucis", 59747, 2.80),  # Delta Cru - west arm
]


@pytest.mark.unit
class TestCruxStarsCatalog:
    """Test that all Crux stars are in the catalog."""

    def test_all_crux_stars_present(self):
        """Verify all 4 Crux stars are in the STAR_CATALOG."""
        catalog_ids = {entry.id for entry in STAR_CATALOG}

        for star_id, name, _, _ in CRUX_STARS:
            assert star_id in catalog_ids, (
                f"Crux star {name} (ID={star_id}) not in catalog"
            )

    def test_crux_star_count(self):
        """Verify we have exactly 4 Crux stars defined."""
        assert len(CRUX_STARS) == 4, "Should have exactly 4 Crux stars"

    def test_acrux_is_brightest(self):
        """Verify Acrux is the brightest Crux star (lowest magnitude)."""
        acrux_mag = None
        other_mags = []

        for star_id, name, _, mag in CRUX_STARS:
            if star_id == ACRUX:
                acrux_mag = mag
            else:
                other_mags.append((name, mag))

        assert acrux_mag is not None, "Acrux not found"
        for name, mag in other_mags:
            assert acrux_mag < mag, (
                f"Acrux ({acrux_mag}) should be brighter than {name} ({mag})"
            )

    @pytest.mark.parametrize("star_id,name,hip,mag", CRUX_STARS)
    def test_star_has_proper_motion(self, star_id, name, hip, mag):
        """Verify each Crux star has proper motion data."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        # Crux stars should have proper motion data
        assert entry.data.pm_ra != 0 or entry.data.pm_dec != 0, (
            f"Star {name} should have proper motion data"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", CRUX_STARS)
    def test_star_has_correct_hip_number(self, star_id, name, hip, mag):
        """Verify each Crux star has correct Hipparcos number."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert entry.hip_number == hip, (
            f"Star {name} should have HIP {hip}, got {entry.hip_number}"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", CRUX_STARS)
    def test_star_nomenclature_is_crux(self, star_id, name, hip, mag):
        """Verify each Crux star has Cru in nomenclature."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert "Cru" in entry.nomenclature, (
            f"Star {name} nomenclature should contain 'Cru', got {entry.nomenclature}"
        )


@pytest.mark.unit
class TestCruxStarsCalculation:
    """Test position calculations for Crux stars."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch."""
        return 2451545.0

    @pytest.mark.parametrize("star_id,name,hip,mag", CRUX_STARS)
    def test_star_position_reasonable(self, standard_jd, star_id, name, hip, mag):
        """Test each Crux star returns a reasonable position."""
        pos, _ = ephem.calc_ut(standard_jd, star_id, 0)

        # Longitude should be 0-360
        assert 0 <= pos[0] < 360, f"{name} longitude {pos[0]}deg out of range"

        # Latitude should be reasonable (-90 to 90)
        assert -90 <= pos[1] <= 90, f"{name} latitude {pos[1]}deg out of range"

        # Fixed stars should be very distant
        assert pos[2] > 1000, f"{name} should have large distance, got {pos[2]}"

    def test_crux_stars_in_southern_sky(self, standard_jd):
        """Test that all Crux stars have negative ecliptic latitude (southern sky)."""
        for star_id, name, _, _ in CRUX_STARS:
            pos, _ = ephem.calc_ut(standard_jd, star_id, 0)
            # Crux stars are all in the southern sky (negative ecliptic latitude)
            assert pos[1] < 0, (
                f"{name} should be in southern sky, got latitude {pos[1]:.2f}"
            )

    def test_crux_stars_form_cross(self, standard_jd):
        """Test that the Crux stars form a rough cross pattern."""
        positions = {}
        for star_id, name, _, _ in CRUX_STARS:
            pos, _ = ephem.calc_ut(standard_jd, star_id, 0)
            positions[name] = (pos[0], pos[1])

        # All four stars should be in a compact region
        lons = [p[0] for p in positions.values()]
        lon_range = max(lons) - min(lons)
        # Cross spans about 5-7 degrees
        assert lon_range < 15, f"Crux longitude spread {lon_range}deg too large"

        lats = [p[1] for p in positions.values()]
        lat_range = max(lats) - min(lats)
        assert lat_range < 15, f"Crux latitude spread {lat_range}deg too large"

    def test_gacrux_northernmost(self, standard_jd):
        """Test that Gacrux is the northernmost Crux star (largest latitude)."""
        positions = {}
        for star_id, name, _, _ in CRUX_STARS:
            pos, _ = ephem.calc_ut(standard_jd, star_id, 0)
            positions[name] = (pos[0], pos[1])

        gacrux_lat = positions["Gacrux"][1]
        for name, (lon, lat) in positions.items():
            if name != "Gacrux":
                # Gacrux should have the highest (least negative) latitude
                assert gacrux_lat > lat, (
                    f"Gacrux ({gacrux_lat:.2f}) should be north of {name} ({lat:.2f})"
                )

    def test_acrux_southernmost(self, standard_jd):
        """Test that Acrux is the southernmost Crux star (smallest latitude)."""
        positions = {}
        for star_id, name, _, _ in CRUX_STARS:
            pos, _ = ephem.calc_ut(standard_jd, star_id, 0)
            positions[name] = (pos[0], pos[1])

        acrux_lat = positions["Acrux"][1]
        for name, (lon, lat) in positions.items():
            if name != "Acrux":
                # Acrux should have the lowest (most negative) latitude
                assert acrux_lat < lat, (
                    f"Acrux ({acrux_lat:.2f}) should be south of {name} ({lat:.2f})"
                )

    def test_delta_crucis_position(self, standard_jd):
        """Test that Delta Crucis is positioned correctly in the cross."""
        pos, _ = ephem.calc_ut(standard_jd, DELTA_CRUCIS, 0)

        # Delta Crucis should be in the Virgo/Libra ecliptic region
        assert 170 < pos[0] < 220, (
            f"Delta Crucis longitude {pos[0]:.1f} out of expected range"
        )

        # Delta Crucis has negative ecliptic latitude (far south)
        assert -60 < pos[1] < -40, (
            f"Delta Crucis latitude {pos[1]:.2f} out of expected range"
        )


@pytest.mark.unit
class TestCruxStarsNameResolution:
    """Test name resolution for Crux stars."""

    def test_resolve_acrux(self):
        """Test Acrux name resolution."""
        assert resolve_star_name("Acrux") == ACRUX
        assert resolve_star_name("Alpha Crucis") == ACRUX
        assert resolve_star_name("Alpha Cru") == ACRUX

    def test_resolve_mimosa(self):
        """Test Mimosa name resolution."""
        assert resolve_star_name("Mimosa") == MIMOSA
        assert resolve_star_name("Beta Crucis") == MIMOSA
        assert resolve_star_name("Beta Cru") == MIMOSA
        assert resolve_star_name("Becrux") == MIMOSA

    def test_resolve_gacrux(self):
        """Test Gacrux name resolution."""
        assert resolve_star_name("Gacrux") == GACRUX
        assert resolve_star_name("Gamma Crucis") == GACRUX
        assert resolve_star_name("Gamma Cru") == GACRUX
        assert resolve_star_name("Rubidea") == GACRUX

    def test_resolve_delta_crucis(self):
        """Test Delta Crucis name resolution."""
        assert resolve_star_name("Delta Crucis") == DELTA_CRUCIS
        assert resolve_star_name("Delta Cru") == DELTA_CRUCIS
        assert resolve_star_name("Decrux") == DELTA_CRUCIS

    def test_canonical_names(self):
        """Test canonical name retrieval for Crux stars."""
        assert get_canonical_star_name(ACRUX) == "Acrux"
        assert get_canonical_star_name(MIMOSA) == "Mimosa"
        assert get_canonical_star_name(GACRUX) == "Gacrux"
        assert get_canonical_star_name(DELTA_CRUCIS) == "Delta Crucis"
