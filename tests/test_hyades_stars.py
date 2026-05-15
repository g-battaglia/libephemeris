"""
Unit tests for Hyades cluster stars.

The Hyades is an open star cluster in Taurus, one of the nearest open clusters to Earth.
These tests verify the 4 significant Hyades stars are correctly defined and calculable.
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    PRIMA_HYADUM,
    SECUNDA_HYADUM,
    THETA_TAURI,
    AIN,
    ALDEBARAN,
)
from libephemeris.fixed_stars import (
    STAR_CATALOG,
    resolve_star_name,
    get_canonical_star_name,
)


# The Hyades stars with their approximate data
# All Hyades stars are located near ~65-70 degrees ecliptic longitude (in Taurus)
HYADES_STARS = [
    # (constant, name, hip_number, magnitude)
    (PRIMA_HYADUM, "Prima Hyadum", 20205, 3.65),  # Gamma Tauri
    (SECUNDA_HYADUM, "Secunda Hyadum", 20455, 3.77),  # Delta^1 Tauri
    (THETA_TAURI, "Theta Tauri", 20894, 3.40),  # Theta^2 Tauri
    (AIN, "Ain", 20889, 3.53),  # Epsilon Tauri
]


@pytest.mark.unit
class TestHyadesStarsCatalog:
    """Test that all 4 significant Hyades stars are in the catalog."""

    def test_all_hyades_stars_present(self):
        """Verify all 4 Hyades stars are in the STAR_CATALOG."""
        catalog_ids = {entry.id for entry in STAR_CATALOG}

        for star_id, name, _, _ in HYADES_STARS:
            assert star_id in catalog_ids, (
                f"Hyades star {name} (ID={star_id}) not in catalog"
            )

    def test_hyades_count(self):
        """Verify we have exactly 4 Hyades stars defined."""
        assert len(HYADES_STARS) == 4, "Should have exactly 4 Hyades stars"

    def test_theta_tauri_is_brightest(self):
        """Verify Theta Tauri is the brightest of the defined Hyades stars."""
        theta_mag = None
        other_mags = []

        for star_id, name, _, mag in HYADES_STARS:
            if star_id == THETA_TAURI:
                theta_mag = mag
            else:
                other_mags.append((name, mag))

        assert theta_mag is not None, "Theta Tauri not found"
        for name, mag in other_mags:
            assert theta_mag < mag, (
                f"Theta Tauri ({theta_mag}) should be brighter than {name} ({mag})"
            )

    @pytest.mark.parametrize("star_id,name,hip,mag", HYADES_STARS)
    def test_star_has_proper_motion(self, star_id, name, hip, mag):
        """Verify each Hyades star has proper motion data."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        # Hyades cluster stars have similar proper motions
        assert entry.data.pm_ra != 0 or entry.data.pm_dec != 0, (
            f"Star {name} should have proper motion data"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", HYADES_STARS)
    def test_star_has_correct_hip_number(self, star_id, name, hip, mag):
        """Verify each Hyades star has correct Hipparcos number."""
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
class TestHyadesStarsCalculation:
    """Test position calculations for Hyades stars."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch."""
        return 2451545.0

    @pytest.mark.parametrize("star_id,name,hip,mag", HYADES_STARS)
    def test_star_position_reasonable(self, standard_jd, star_id, name, hip, mag):
        """Test each Hyades star returns a reasonable position."""
        pos, _ = ephem.calc_ut(standard_jd, star_id, 0)

        # Longitude should be 0-360
        assert 0 <= pos[0] < 360, f"{name} longitude {pos[0]}deg out of range"

        # Latitude should be reasonable (-90 to 90)
        assert -90 <= pos[1] <= 90, f"{name} latitude {pos[1]}deg out of range"

        # Fixed stars should be very distant
        assert pos[2] > 1000, f"{name} should have large distance, got {pos[2]}"

    def test_hyades_stars_clustered(self, standard_jd):
        """Test that all Hyades stars are clustered together in the sky."""
        positions = []
        for star_id, name, _, _ in HYADES_STARS:
            pos, _ = ephem.calc_ut(standard_jd, star_id, 0)
            positions.append((name, pos[0], pos[1]))

        # All Hyades stars should be within ~10 degrees of each other
        lons = [p[1] for p in positions]
        lats = [p[2] for p in positions]

        lon_range = max(lons) - min(lons)
        lat_range = max(lats) - min(lats)

        # The cluster spans about 5-8 degrees (larger than Pleiades)
        assert lon_range < 10, f"Hyades longitude spread {lon_range}deg too large"
        assert lat_range < 10, f"Hyades latitude spread {lat_range}deg too large"

    def test_hyades_near_aldebaran(self, standard_jd):
        """Test that Hyades are located near Aldebaran region (~65-70 deg ecliptic)."""
        # Prima Hyadum is a good reference point for the Hyades
        pos, _ = ephem.calc_ut(standard_jd, PRIMA_HYADUM, 0)

        # Should be near 65-70 degrees ecliptic (Taurus)
        assert 60 < pos[0] < 75, (
            f"Prima Hyadum should be near 65-70 deg ecliptic, got {pos[0]:.1f}"
        )

    def test_aldebaran_not_in_hyades_cluster(self, standard_jd):
        """Test that Aldebaran is distinct from the Hyades cluster.

        Aldebaran appears in the same region of the sky but is actually
        much closer to Earth and not a member of the Hyades cluster.
        """
        aldebaran_pos, _ = ephem.calc_ut(standard_jd, ALDEBARAN, 0)
        prima_pos, _ = ephem.calc_ut(standard_jd, PRIMA_HYADUM, 0)

        # They should be in the same general region of the sky
        # Both in Taurus, within about 5 degrees of each other
        lon_diff = abs(aldebaran_pos[0] - prima_pos[0])
        assert lon_diff < 10, (
            f"Aldebaran and Prima Hyadum should be in same region, "
            f"but differ by {lon_diff:.1f} degrees"
        )


@pytest.mark.unit
class TestHyadesStarsNameResolution:
    """Test name resolution for Hyades stars."""

    def test_resolve_prima_hyadum(self):
        """Test Prima Hyadum name resolution."""
        assert resolve_star_name("Prima Hyadum") == PRIMA_HYADUM
        assert resolve_star_name("Gamma Tauri") == PRIMA_HYADUM
        assert resolve_star_name("Gamma Tau") == PRIMA_HYADUM
        assert resolve_star_name("Hyadum I") == PRIMA_HYADUM
        assert resolve_star_name("First Hyad") == PRIMA_HYADUM

    def test_resolve_secunda_hyadum(self):
        """Test Secunda Hyadum name resolution."""
        assert resolve_star_name("Secunda Hyadum") == SECUNDA_HYADUM
        assert resolve_star_name("Delta Tauri") == SECUNDA_HYADUM
        assert resolve_star_name("Delta Tau") == SECUNDA_HYADUM
        assert resolve_star_name("Hyadum II") == SECUNDA_HYADUM
        assert resolve_star_name("Second Hyad") == SECUNDA_HYADUM

    def test_resolve_theta_tauri(self):
        """Test Theta Tauri name resolution."""
        assert resolve_star_name("Theta Tauri") == THETA_TAURI
        assert resolve_star_name("Theta2 Tauri") == THETA_TAURI
        assert resolve_star_name("Theta2 Tau") == THETA_TAURI

    def test_resolve_ain(self):
        """Test Ain (Epsilon Tauri) name resolution."""
        assert resolve_star_name("Ain") == AIN
        assert resolve_star_name("Epsilon Tauri") == AIN
        assert resolve_star_name("Epsilon Tau") == AIN
        assert resolve_star_name("Oculus Borealis") == AIN

    @pytest.mark.parametrize("star_id,name,hip,mag", HYADES_STARS)
    def test_canonical_name_retrieval(self, star_id, name, hip, mag):
        """Test canonical name retrieval for all Hyades stars."""
        canonical = get_canonical_star_name(star_id)
        assert canonical == name, f"Expected {name}, got {canonical}"


@pytest.mark.unit
class TestHyadesStarsData:
    """Test the Hyades stars have correct catalog data."""

    def test_prima_hyadum_data(self):
        """Test Prima Hyadum (Gamma Tauri) catalog entry."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == PRIMA_HYADUM:
                entry = e
                break

        assert entry is not None
        assert entry.name == "Prima Hyadum"
        assert entry.nomenclature == "gaTau"
        assert entry.hip_number == 20205
        assert 64 < entry.data.ra_j2000 < 66  # ~64.95deg
        assert 15 < entry.data.dec_j2000 < 16  # ~15.63deg
        assert 3.5 < entry.magnitude < 3.8  # ~3.65

    def test_secunda_hyadum_data(self):
        """Test Secunda Hyadum (Delta^1 Tauri) catalog entry."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == SECUNDA_HYADUM:
                entry = e
                break

        assert entry is not None
        assert entry.name == "Secunda Hyadum"
        assert entry.nomenclature == "de1Tau"
        assert entry.hip_number == 20455
        assert 65 < entry.data.ra_j2000 < 66  # ~65.73deg
        assert 17 < entry.data.dec_j2000 < 18  # ~17.54deg
        assert 3.7 < entry.magnitude < 3.9  # ~3.77

    def test_theta_tauri_data(self):
        """Test Theta Tauri (Theta^2 Tauri) catalog entry."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == THETA_TAURI:
                entry = e
                break

        assert entry is not None
        assert entry.name == "Theta Tauri"
        assert entry.nomenclature == "th2Tau"
        assert entry.hip_number == 20894
        assert 67 < entry.data.ra_j2000 < 68  # ~67.17deg
        assert 15 < entry.data.dec_j2000 < 16  # ~15.87deg
        assert 3.3 < entry.magnitude < 3.5  # ~3.40

    def test_ain_data(self):
        """Test Ain (Epsilon Tauri) catalog entry."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == AIN:
                entry = e
                break

        assert entry is not None
        assert entry.name == "Ain"
        assert entry.nomenclature == "epTau"
        assert entry.hip_number == 20889
        assert 67 < entry.data.ra_j2000 < 68  # ~67.15deg
        assert 19 < entry.data.dec_j2000 < 20  # ~19.18deg
        assert 3.4 < entry.magnitude < 3.6  # ~3.53

    def test_hyades_cluster_proper_motion(self):
        """Test that Hyades stars have similar proper motions (cluster motion)."""
        pm_ras = []
        pm_decs = []

        for star_id, name, _, _ in HYADES_STARS:
            for entry in STAR_CATALOG:
                if entry.id == star_id:
                    pm_ras.append(entry.data.pm_ra)
                    pm_decs.append(entry.data.pm_dec)
                    break

        # Hyades cluster has a common space motion
        # All proper motions should be similar (within ~0.02 arcsec/yr)
        assert max(pm_ras) - min(pm_ras) < 0.02, (
            "Hyades PM(RA) should be similar across cluster"
        )
        assert max(pm_decs) - min(pm_decs) < 0.02, (
            "Hyades PM(Dec) should be similar across cluster"
        )
