"""
Unit tests for Leo constellation stars.

Leo is one of the zodiac constellations, featuring these 4 major bright stars:
- Regulus (Alpha Leo) - Brightest star, Royal Star of Persia ("Watcher of the North")
- Denebola (Beta Leo) - Lion's tail
- Algieba (Gamma Leo) - Lion's mane, beautiful double star
- Zosma (Delta Leo) - Lion's hip/back
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    REGULUS,
    DENEBOLA,
    ALGIEBA,
    ZOSMA,
)
from libephemeris.fixed_stars import (
    STAR_CATALOG,
    resolve_star_name,
    get_canonical_star_name,
)


# The 4 major Leo stars with their properties
# (constant, name, hip_number, magnitude)
LEO_STARS = [
    (REGULUS, "Regulus", 49669, 1.40),  # Alpha Leo - Royal Star
    (DENEBOLA, "Denebola", 57632, 2.14),  # Beta Leo - lion's tail
    (ALGIEBA, "Algieba", 50583, 2.08),  # Gamma Leo - lion's mane
    (ZOSMA, "Zosma", 54872, 2.56),  # Delta Leo - lion's hip
]


@pytest.mark.unit
class TestLeoStarsCatalog:
    """Test that all 4 major Leo stars are in the catalog."""

    def test_all_leo_stars_present(self):
        """Verify all 4 Leo stars are in the STAR_CATALOG."""
        catalog_ids = {entry.id for entry in STAR_CATALOG}

        for star_id, name, _, _ in LEO_STARS:
            assert star_id in catalog_ids, (
                f"Leo star {name} (ID={star_id}) not in catalog"
            )

    def test_leo_star_count(self):
        """Verify we have exactly 4 major Leo stars defined."""
        assert len(LEO_STARS) == 4, "Should have exactly 4 major Leo stars"

    def test_regulus_is_brightest(self):
        """Verify Regulus is the brightest Leo star (lowest magnitude)."""
        regulus_mag = None
        other_mags = []

        for star_id, name, _, mag in LEO_STARS:
            if star_id == REGULUS:
                regulus_mag = mag
            else:
                other_mags.append((name, mag))

        assert regulus_mag is not None, "Regulus not found"
        for name, mag in other_mags:
            assert regulus_mag < mag, (
                f"Regulus ({regulus_mag}) should be brighter than {name} ({mag})"
            )

    @pytest.mark.parametrize("star_id,name,hip,mag", LEO_STARS)
    def test_star_has_proper_motion(self, star_id, name, hip, mag):
        """Verify each Leo star has proper motion data."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        # Leo stars should have proper motion data
        assert entry.data.pm_ra != 0 or entry.data.pm_dec != 0, (
            f"Star {name} should have proper motion data"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", LEO_STARS)
    def test_star_has_correct_hip_number(self, star_id, name, hip, mag):
        """Verify each Leo star has correct Hipparcos number."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert entry.hip_number == hip, (
            f"Star {name} should have HIP {hip}, got {entry.hip_number}"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", LEO_STARS)
    def test_star_nomenclature_is_leo(self, star_id, name, hip, mag):
        """Verify each Leo star has Leo in nomenclature."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert "Leo" in entry.nomenclature, (
            f"Star {name} nomenclature should contain 'Leo', got {entry.nomenclature}"
        )


@pytest.mark.unit
class TestLeoStarsCalculation:
    """Test position calculations for Leo stars."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch."""
        return 2451545.0

    @pytest.mark.parametrize("star_id,name,hip,mag", LEO_STARS)
    def test_star_position_reasonable(self, standard_jd, star_id, name, hip, mag):
        """Test each Leo star returns a reasonable position."""
        pos, _ = ephem.calc_ut(standard_jd, star_id, 0)

        # Longitude should be 0-360
        assert 0 <= pos[0] < 360, f"{name} longitude {pos[0]}deg out of range"

        # Latitude should be reasonable (-90 to 90)
        assert -90 <= pos[1] <= 90, f"{name} latitude {pos[1]}deg out of range"

        # Fixed stars should be very distant
        assert pos[2] > 1000, f"{name} should have large distance, got {pos[2]}"

    def test_leo_stars_in_leo_virgo_region(self, standard_jd):
        """Test that Leo stars are in the Leo/Virgo ecliptic region (~120-180 deg)."""
        for star_id, name, _, _ in LEO_STARS:
            pos, _ = ephem.calc_ut(standard_jd, star_id, 0)

            # Leo stars should be between ~145 (late Leo) and ~185 (early Virgo)
            # This accounts for proper motion and precession
            assert 140 < pos[0] < 185, (
                f"{name} should be in Leo/Virgo region, got {pos[0]:.1f} degrees"
            )

    def test_regulus_near_ecliptic(self, standard_jd):
        """Test that Regulus is very close to the ecliptic."""
        pos, _ = ephem.calc_ut(standard_jd, REGULUS, 0)

        # Regulus is famous for being almost exactly on the ecliptic
        # Its ecliptic latitude should be less than 1 degree
        assert abs(pos[1]) < 1.0, (
            f"Regulus should be near ecliptic, got latitude {pos[1]:.2f}"
        )

    def test_regulus_position(self, standard_jd):
        """Test that Regulus is positioned correctly."""
        pos, _ = ephem.calc_ut(standard_jd, REGULUS, 0)

        # Regulus should be near 150 degrees ecliptic longitude (Leo)
        assert 145 < pos[0] < 155, (
            f"Regulus longitude {pos[0]:.1f} out of expected range"
        )

    def test_denebola_position(self, standard_jd):
        """Test that Denebola (lion's tail) is positioned after Regulus."""
        regulus_pos, _ = ephem.calc_ut(standard_jd, REGULUS, 0)
        denebola_pos, _ = ephem.calc_ut(standard_jd, DENEBOLA, 0)

        # Denebola should be further along the ecliptic than Regulus
        assert denebola_pos[0] > regulus_pos[0], (
            f"Denebola ({denebola_pos[0]:.1f}) should be after Regulus ({regulus_pos[0]:.1f})"
        )

    def test_zosma_position(self, standard_jd):
        """Test that Zosma is positioned correctly in Leo."""
        pos, _ = ephem.calc_ut(standard_jd, ZOSMA, 0)

        # Zosma should be near 155-165 degrees ecliptic longitude
        assert 150 < pos[0] < 170, f"Zosma longitude {pos[0]:.1f} out of expected range"

        # Zosma has positive ecliptic latitude (north of ecliptic)
        assert pos[1] > 0, f"Zosma should be north of ecliptic, got lat {pos[1]:.2f}"


@pytest.mark.unit
class TestLeoStarsNameResolution:
    """Test name resolution for Leo stars."""

    def test_resolve_regulus(self):
        """Test Regulus name resolution."""
        assert resolve_star_name("Regulus") == REGULUS
        assert resolve_star_name("Alpha Leonis") == REGULUS
        assert resolve_star_name("Alpha Leo") == REGULUS
        assert resolve_star_name("Cor Leonis") == REGULUS
        assert resolve_star_name("Watcher of the North") == REGULUS

    def test_resolve_denebola(self):
        """Test Denebola name resolution."""
        assert resolve_star_name("Denebola") == DENEBOLA
        assert resolve_star_name("Beta Leonis") == DENEBOLA
        assert resolve_star_name("Beta Leo") == DENEBOLA
        assert resolve_star_name("Lion's Tail") == DENEBOLA

    def test_resolve_algieba(self):
        """Test Algieba name resolution."""
        assert resolve_star_name("Algieba") == ALGIEBA
        assert resolve_star_name("Gamma Leonis") == ALGIEBA
        assert resolve_star_name("Gamma Leo") == ALGIEBA
        assert resolve_star_name("Lion's Mane") == ALGIEBA

    def test_resolve_zosma(self):
        """Test Zosma name resolution."""
        assert resolve_star_name("Zosma") == ZOSMA
        assert resolve_star_name("Delta Leonis") == ZOSMA
        assert resolve_star_name("Delta Leo") == ZOSMA
        assert resolve_star_name("Dhur") == ZOSMA
        assert resolve_star_name("Duhr") == ZOSMA
        assert resolve_star_name("Lion's Hip") == ZOSMA
        assert resolve_star_name("Lion's Back") == ZOSMA

    def test_canonical_names(self):
        """Test canonical name retrieval for Leo stars."""
        assert get_canonical_star_name(REGULUS) == "Regulus"
        assert get_canonical_star_name(DENEBOLA) == "Denebola"
        assert get_canonical_star_name(ALGIEBA) == "Algieba"
        assert get_canonical_star_name(ZOSMA) == "Zosma"
