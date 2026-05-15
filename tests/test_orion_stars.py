"""
Unit tests for Orion constellation stars.

Orion is one of the most recognizable constellations, featuring 8 major stars:
- Betelgeuse (Alpha Ori) - Red supergiant, shoulder
- Rigel (Beta Ori) - Blue supergiant, foot
- Bellatrix (Gamma Ori) - Amazon Star, shoulder
- Alnilam (Epsilon Ori) - Belt center
- Alnitak (Zeta Ori) - Belt left
- Mintaka (Delta Ori) - Belt right
- Saiph (Kappa Ori) - Foot
- Meissa (Lambda Ori) - Head
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    BETELGEUSE,
    RIGEL,
    BELLATRIX,
    ALNILAM,
    ALNITAK,
    MINTAKA,
    SAIPH,
    MEISSA,
)
from libephemeris.fixed_stars import (
    STAR_CATALOG,
    resolve_star_name,
    get_canonical_star_name,
)


# The 8 major Orion stars with their properties
# (constant, name, hip_number, magnitude)
ORION_STARS = [
    (BETELGEUSE, "Betelgeuse", 27989, 0.42),  # Alpha Ori - red supergiant
    (RIGEL, "Rigel", 24436, 0.12),  # Beta Ori - brightest
    (BELLATRIX, "Bellatrix", 25336, 1.64),  # Gamma Ori
    (ALNILAM, "Alnilam", 26311, 1.69),  # Epsilon Ori - belt center
    (ALNITAK, "Alnitak", 26727, 1.74),  # Zeta Ori - belt left
    (MINTAKA, "Mintaka", 25930, 2.23),  # Delta Ori - belt right
    (SAIPH, "Saiph", 27366, 2.06),  # Kappa Ori
    (MEISSA, "Meissa", 26207, 3.33),  # Lambda Ori - head
]

# The three belt stars
ORION_BELT_STARS = [
    (ALNILAM, "Alnilam"),
    (ALNITAK, "Alnitak"),
    (MINTAKA, "Mintaka"),
]


@pytest.mark.unit
class TestOrionStarsCatalog:
    """Test that all 8 major Orion stars are in the catalog."""

    def test_all_orion_stars_present(self):
        """Verify all 8 Orion stars are in the STAR_CATALOG."""
        catalog_ids = {entry.id for entry in STAR_CATALOG}

        for star_id, name, _, _ in ORION_STARS:
            assert star_id in catalog_ids, (
                f"Orion star {name} (ID={star_id}) not in catalog"
            )

    def test_orion_star_count(self):
        """Verify we have exactly 8 major Orion stars defined."""
        assert len(ORION_STARS) == 8, "Should have exactly 8 major Orion stars"

    def test_rigel_is_brightest(self):
        """Verify Rigel is the brightest Orion star (lowest magnitude)."""
        rigel_mag = None
        other_mags = []

        for star_id, name, _, mag in ORION_STARS:
            if star_id == RIGEL:
                rigel_mag = mag
            else:
                other_mags.append((name, mag))

        assert rigel_mag is not None, "Rigel not found"
        for name, mag in other_mags:
            assert rigel_mag < mag, (
                f"Rigel ({rigel_mag}) should be brighter than {name} ({mag})"
            )

    @pytest.mark.parametrize("star_id,name,hip,mag", ORION_STARS)
    def test_star_has_proper_motion(self, star_id, name, hip, mag):
        """Verify each Orion star has proper motion data."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        # Orion stars should have proper motion data
        assert entry.data.pm_ra != 0 or entry.data.pm_dec != 0, (
            f"Star {name} should have proper motion data"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", ORION_STARS)
    def test_star_has_correct_hip_number(self, star_id, name, hip, mag):
        """Verify each Orion star has correct Hipparcos number."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert entry.hip_number == hip, (
            f"Star {name} should have HIP {hip}, got {entry.hip_number}"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", ORION_STARS)
    def test_star_nomenclature_is_orion(self, star_id, name, hip, mag):
        """Verify each Orion star has Ori in nomenclature."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert "Ori" in entry.nomenclature, (
            f"Star {name} nomenclature should contain 'Ori', got {entry.nomenclature}"
        )


@pytest.mark.unit
class TestOrionStarsCalculation:
    """Test position calculations for Orion stars."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch."""
        return 2451545.0

    @pytest.mark.parametrize("star_id,name,hip,mag", ORION_STARS)
    def test_star_position_reasonable(self, standard_jd, star_id, name, hip, mag):
        """Test each Orion star returns a reasonable position."""
        pos, _ = ephem.calc_ut(standard_jd, star_id, 0)

        # Longitude should be 0-360
        assert 0 <= pos[0] < 360, f"{name} longitude {pos[0]}deg out of range"

        # Latitude should be reasonable (-90 to 90)
        assert -90 <= pos[1] <= 90, f"{name} latitude {pos[1]}deg out of range"

        # Fixed stars should be very distant
        assert pos[2] > 1000, f"{name} should have large distance, got {pos[2]}"

    def test_orion_belt_stars_aligned(self, standard_jd):
        """Test that the three belt stars are roughly aligned."""
        positions = []
        for star_id, name in ORION_BELT_STARS:
            pos, _ = ephem.calc_ut(standard_jd, star_id, 0)
            positions.append((name, pos[0], pos[1]))

        # Belt stars should be within ~3 degrees of each other in longitude
        lons = [p[1] for p in positions]
        lon_range = max(lons) - min(lons)
        assert lon_range < 5, f"Belt longitude spread {lon_range}deg too large"

        # Belt stars should be close in latitude too
        lats = [p[2] for p in positions]
        lat_range = max(lats) - min(lats)
        assert lat_range < 3, f"Belt latitude spread {lat_range}deg too large"

    def test_orion_in_gemini_region(self, standard_jd):
        """Test that Orion stars are in the Gemini ecliptic region (~60-90 deg)."""
        # Betelgeuse is a good reference for Orion's position
        pos, _ = ephem.calc_ut(standard_jd, BETELGEUSE, 0)

        # Should be near 88 degrees ecliptic (Gemini)
        assert 80 < pos[0] < 95, (
            f"Betelgeuse should be near 88 deg ecliptic, got {pos[0]:.1f}"
        )

    def test_meissa_in_orion_region(self, standard_jd):
        """Test that Meissa is located in the correct region (Orion's head)."""
        pos, _ = ephem.calc_ut(standard_jd, MEISSA, 0)

        # Meissa should be near 83 degrees ecliptic longitude
        assert 75 < pos[0] < 90, f"Meissa longitude {pos[0]:.1f} out of expected range"

        # Meissa is below the ecliptic (southern ecliptic latitude)
        assert -20 < pos[1] < 0, f"Meissa latitude {pos[1]:.2f} out of expected range"


@pytest.mark.unit
class TestOrionStarsNameResolution:
    """Test name resolution for Orion stars."""

    def test_resolve_betelgeuse(self):
        """Test Betelgeuse name resolution."""
        assert resolve_star_name("Betelgeuse") == BETELGEUSE
        assert resolve_star_name("Alpha Orionis") == BETELGEUSE
        assert resolve_star_name("Alpha Ori") == BETELGEUSE

    def test_resolve_rigel(self):
        """Test Rigel name resolution."""
        assert resolve_star_name("Rigel") == RIGEL
        assert resolve_star_name("Beta Orionis") == RIGEL
        assert resolve_star_name("Beta Ori") == RIGEL

    def test_resolve_bellatrix(self):
        """Test Bellatrix name resolution."""
        assert resolve_star_name("Bellatrix") == BELLATRIX
        assert resolve_star_name("Gamma Orionis") == BELLATRIX

    def test_resolve_belt_stars(self):
        """Test belt stars name resolution."""
        assert resolve_star_name("Alnilam") == ALNILAM
        assert resolve_star_name("Epsilon Orionis") == ALNILAM

        assert resolve_star_name("Alnitak") == ALNITAK
        assert resolve_star_name("Zeta Orionis") == ALNITAK

        assert resolve_star_name("Mintaka") == MINTAKA
        assert resolve_star_name("Delta Orionis") == MINTAKA

    def test_resolve_saiph(self):
        """Test Saiph name resolution."""
        assert resolve_star_name("Saiph") == SAIPH
        assert resolve_star_name("Kappa Orionis") == SAIPH

    def test_resolve_meissa(self):
        """Test Meissa name resolution."""
        assert resolve_star_name("Meissa") == MEISSA
        assert resolve_star_name("Lambda Orionis") == MEISSA
        assert resolve_star_name("Lambda Ori") == MEISSA
        assert resolve_star_name("Heka") == MEISSA
        assert resolve_star_name("Head of Orion") == MEISSA

    def test_canonical_names(self):
        """Test canonical name retrieval for Orion stars."""
        assert get_canonical_star_name(BETELGEUSE) == "Betelgeuse"
        assert get_canonical_star_name(RIGEL) == "Rigel"
        assert get_canonical_star_name(BELLATRIX) == "Bellatrix"
        assert get_canonical_star_name(ALNILAM) == "Alnilam"
        assert get_canonical_star_name(ALNITAK) == "Alnitak"
        assert get_canonical_star_name(MINTAKA) == "Mintaka"
        assert get_canonical_star_name(SAIPH) == "Saiph"
        assert get_canonical_star_name(MEISSA) == "Meissa"
