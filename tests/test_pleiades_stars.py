"""
Unit tests for Pleiades cluster stars.

The Pleiades (M45) is an open star cluster in Taurus, also known as the Seven Sisters.
These tests verify the 9 visible Pleiades stars are correctly defined and calculable.
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    ALCYONE,
    ASTEROPE,
    CELAENO,
    ELECTRA,
    MAIA,
    MEROPE,
    TAYGETA,
    ATLAS,
    PLEIONE,
)
from libephemeris.fixed_stars import (
    STAR_CATALOG,
    resolve_star_name,
    get_canonical_star_name,
)


# The visible Pleiades stars with their approximate ecliptic longitudes
# All Pleiades stars are located near ~60° ecliptic longitude (near 0° Gemini)
PLEIADES_STARS = [
    # (constant, name, hip_number, magnitude)
    (ALCYONE, "Alcyone", 17702, 2.87),  # Brightest, Eta Tauri
    (ASTEROPE, "Asterope", 17579, 5.76),  # 21 Tauri
    (CELAENO, "Celaeno", 17489, 5.45),  # 16 Tauri
    (ELECTRA, "Electra", 17499, 3.70),  # 17 Tauri
    (MAIA, "Maia", 17573, 3.87),  # 20 Tauri
    (MEROPE, "Merope", 17608, 4.14),  # 23 Tauri
    (TAYGETA, "Taygeta", 17531, 4.30),  # 19 Tauri
    (ATLAS, "Atlas", 17847, 3.62),  # 27 Tauri
    (PLEIONE, "Pleione", 17851, 5.09),  # 28 Tauri
]


@pytest.mark.unit
class TestPleiadesStarsCatalog:
    """Test that all 9 visible Pleiades stars are in the catalog."""

    def test_all_pleiades_stars_present(self):
        """Verify all 9 Pleiades stars are in the STAR_CATALOG."""
        catalog_ids = {entry.id for entry in STAR_CATALOG}

        for star_id, name, _, _ in PLEIADES_STARS:
            assert star_id in catalog_ids, (
                f"Pleiades star {name} (ID={star_id}) not in catalog"
            )

    def test_pleiades_count(self):
        """Verify we have exactly 9 Pleiades stars defined."""
        assert len(PLEIADES_STARS) == 9, "Should have exactly 9 Pleiades stars"

    def test_alcyone_is_brightest(self):
        """Verify Alcyone is the brightest Pleiades star."""
        alcyone_mag = None
        other_mags = []

        for star_id, name, _, mag in PLEIADES_STARS:
            if star_id == ALCYONE:
                alcyone_mag = mag
            else:
                other_mags.append((name, mag))

        assert alcyone_mag is not None, "Alcyone not found"
        for name, mag in other_mags:
            assert alcyone_mag < mag, (
                f"Alcyone ({alcyone_mag}) should be brighter than {name} ({mag})"
            )

    @pytest.mark.parametrize("star_id,name,hip,mag", PLEIADES_STARS)
    def test_star_has_proper_motion(self, star_id, name, hip, mag):
        """Verify each Pleiades star has proper motion data."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        # Pleiades cluster stars have similar proper motions
        assert entry.data.pm_ra != 0 or entry.data.pm_dec != 0, (
            f"Star {name} should have proper motion data"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", PLEIADES_STARS)
    def test_star_has_correct_hip_number(self, star_id, name, hip, mag):
        """Verify each Pleiades star has correct Hipparcos number."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert entry.hip_number == hip, (
            f"Star {name} should have HIP {hip}, got {entry.hip_number}"
        )


@pytest.mark.unit
class TestPleiadesStarsCalculation:
    """Test position calculations for Pleiades stars."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch."""
        return 2451545.0

    @pytest.mark.parametrize("star_id,name,hip,mag", PLEIADES_STARS)
    def test_star_position_reasonable(self, standard_jd, star_id, name, hip, mag):
        """Test each Pleiades star returns a reasonable position."""
        pos, _ = ephem.calc_ut(standard_jd, star_id, 0)

        # Longitude should be 0-360
        assert 0 <= pos[0] < 360, f"{name} longitude {pos[0]}deg out of range"

        # Latitude should be reasonable (-90 to 90)
        assert -90 <= pos[1] <= 90, f"{name} latitude {pos[1]}deg out of range"

        # Fixed stars should be very distant
        assert pos[2] > 1000, f"{name} should have large distance, got {pos[2]}"

    def test_pleiades_stars_clustered(self, standard_jd):
        """Test that all Pleiades stars are clustered together in the sky."""
        positions = []
        for star_id, name, _, _ in PLEIADES_STARS:
            pos, _ = ephem.calc_ut(standard_jd, star_id, 0)
            positions.append((name, pos[0], pos[1]))

        # All Pleiades stars should be within ~5 degrees of each other
        lons = [p[1] for p in positions]
        lats = [p[2] for p in positions]

        lon_range = max(lons) - min(lons)
        lat_range = max(lats) - min(lats)

        # The cluster spans about 2 degrees
        assert lon_range < 5, f"Pleiades longitude spread {lon_range}deg too large"
        assert lat_range < 5, f"Pleiades latitude spread {lat_range}deg too large"

    def test_pleiades_near_taurus(self, standard_jd):
        """Test that Pleiades are located in Taurus region (~60 deg ecliptic)."""
        # Alcyone is the reference point for the Pleiades
        pos, _ = ephem.calc_ut(standard_jd, ALCYONE, 0)

        # Should be near 60 degrees ecliptic (late Taurus / early Gemini)
        # Allow some tolerance for precession
        assert 55 < pos[0] < 65, (
            f"Alcyone should be near 60 deg ecliptic, got {pos[0]:.1f}"
        )


@pytest.mark.unit
class TestPleiadesStarsNameResolution:
    """Test name resolution for Pleiades stars."""

    def test_resolve_alcyone(self):
        """Test Alcyone name resolution (already tested in behenian, but verify)."""
        assert resolve_star_name("Alcyone") == ALCYONE
        assert resolve_star_name("Eta Tauri") == ALCYONE
        assert resolve_star_name("Pleiades") == ALCYONE
        assert resolve_star_name("Seven Sisters") == ALCYONE

    def test_resolve_asterope(self):
        """Test Asterope name resolution."""
        assert resolve_star_name("Asterope") == ASTEROPE
        assert resolve_star_name("21 Tauri") == ASTEROPE
        assert resolve_star_name("21 Tau") == ASTEROPE
        assert resolve_star_name("Sterope") == ASTEROPE

    def test_resolve_celaeno(self):
        """Test Celaeno name resolution."""
        assert resolve_star_name("Celaeno") == CELAENO
        assert resolve_star_name("16 Tauri") == CELAENO
        assert resolve_star_name("16 Tau") == CELAENO
        assert resolve_star_name("Celeno") == CELAENO

    def test_resolve_electra(self):
        """Test Electra name resolution."""
        assert resolve_star_name("Electra") == ELECTRA
        assert resolve_star_name("17 Tauri") == ELECTRA
        assert resolve_star_name("17 Tau") == ELECTRA

    def test_resolve_maia(self):
        """Test Maia name resolution."""
        assert resolve_star_name("Maia") == MAIA
        assert resolve_star_name("20 Tauri") == MAIA
        assert resolve_star_name("20 Tau") == MAIA

    def test_resolve_merope(self):
        """Test Merope name resolution."""
        assert resolve_star_name("Merope") == MEROPE
        assert resolve_star_name("23 Tauri") == MEROPE
        assert resolve_star_name("23 Tau") == MEROPE

    def test_resolve_taygeta(self):
        """Test Taygeta name resolution."""
        assert resolve_star_name("Taygeta") == TAYGETA
        assert resolve_star_name("19 Tauri") == TAYGETA
        assert resolve_star_name("19 Tau") == TAYGETA
        assert resolve_star_name("Taygete") == TAYGETA

    def test_resolve_atlas(self):
        """Test Atlas name resolution."""
        assert resolve_star_name("Atlas") == ATLAS
        assert resolve_star_name("27 Tauri") == ATLAS
        assert resolve_star_name("27 Tau") == ATLAS

    def test_resolve_pleione(self):
        """Test Pleione name resolution."""
        assert resolve_star_name("Pleione") == PLEIONE
        assert resolve_star_name("28 Tauri") == PLEIONE
        assert resolve_star_name("28 Tau") == PLEIONE

    @pytest.mark.parametrize("star_id,name,hip,mag", PLEIADES_STARS)
    def test_canonical_name_retrieval(self, star_id, name, hip, mag):
        """Test canonical name retrieval for all Pleiades stars."""
        canonical = get_canonical_star_name(star_id)
        assert canonical == name, f"Expected {name}, got {canonical}"


@pytest.mark.unit
class TestPleiadesStarsData:
    """Test the Pleiades stars have correct catalog data."""

    def test_alcyone_data(self):
        """Test Alcyone catalog entry."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == ALCYONE:
                entry = e
                break

        assert entry is not None
        assert entry.name == "Alcyone"
        assert entry.nomenclature == "etTau"
        assert entry.hip_number == 17702
        assert 56 < entry.data.ra_j2000 < 57  # ~56.87deg
        assert 24 < entry.data.dec_j2000 < 25  # ~24.1deg
        assert 2.5 < entry.magnitude < 3.0  # ~2.87

    def test_electra_data(self):
        """Test Electra catalog entry (third brightest)."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == ELECTRA:
                entry = e
                break

        assert entry is not None
        assert entry.name == "Electra"
        assert entry.nomenclature == "17Tau"
        assert entry.hip_number == 17499
        assert 56 < entry.data.ra_j2000 < 57
        assert 24 < entry.data.dec_j2000 < 25
        assert 3.5 < entry.magnitude < 4.0  # ~3.70

    def test_atlas_data(self):
        """Test Atlas catalog entry (second brightest)."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == ATLAS:
                entry = e
                break

        assert entry is not None
        assert entry.name == "Atlas"
        assert entry.nomenclature == "27Tau"
        assert entry.hip_number == 17847
        assert 57 < entry.data.ra_j2000 < 58  # ~57.29deg
        assert 24 < entry.data.dec_j2000 < 25
        assert 3.5 < entry.magnitude < 3.8  # ~3.62

    def test_pleione_data(self):
        """Test Pleione catalog entry (mother of the Pleiades)."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == PLEIONE:
                entry = e
                break

        assert entry is not None
        assert entry.name == "Pleione"
        assert entry.nomenclature == "28Tau"
        assert entry.hip_number == 17851
        assert 57 < entry.data.ra_j2000 < 58  # ~57.30deg
        assert 24 < entry.data.dec_j2000 < 25
        assert 5.0 < entry.magnitude < 5.2  # ~5.09

    def test_pleiades_cluster_proper_motion(self):
        """Test that Pleiades stars have similar proper motions (cluster motion)."""
        pm_ras = []
        pm_decs = []

        for star_id, name, _, _ in PLEIADES_STARS:
            for entry in STAR_CATALOG:
                if entry.id == star_id:
                    pm_ras.append(entry.data.pm_ra)
                    pm_decs.append(entry.data.pm_dec)
                    break

        # Pleiades cluster has a common space motion
        # All proper motions should be similar (within ~0.01 arcsec/yr)
        assert max(pm_ras) - min(pm_ras) < 0.01, (
            "Pleiades PM(RA) should be similar across cluster"
        )
        assert max(pm_decs) - min(pm_decs) < 0.01, (
            "Pleiades PM(Dec) should be similar across cluster"
        )
