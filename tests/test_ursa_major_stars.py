"""
Unit tests for Ursa Major (Big Dipper) constellation stars.

The Big Dipper is one of the most recognizable asterisms in the northern sky,
consisting of 7 bright stars that form a distinctive ladle/plough shape:

Bowl stars:
- Dubhe (Alpha UMa) - pointer star, upper lip of bowl
- Merak (Beta UMa) - pointer star, lower lip of bowl
- Phecda (Gamma UMa) - bottom corner of bowl
- Megrez (Delta UMa) - connects bowl to handle

Handle stars:
- Alioth (Epsilon UMa) - brightest of the seven
- Mizar (Zeta UMa) - famous double star, "horse"
- Alcor (80 UMa) - optical double with Mizar, "rider"
- Alkaid (Eta UMa) - end of handle

The two "pointer" stars (Dubhe and Merak) point to Polaris, the North Star.
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    DUBHE,
    MERAK,
    PHECDA,
    MEGREZ,
    ALIOTH,
    MIZAR,
    ALCOR,
    ALKAID,
)
from libephemeris.fixed_stars import (
    STAR_CATALOG,
    resolve_star_name,
    get_canonical_star_name,
)


# The 7 Big Dipper stars plus Alcor (visual companion to Mizar)
# (constant, name, hip_number, magnitude)
BIG_DIPPER_STARS = [
    (DUBHE, "Dubhe", 54061, 1.79),  # Alpha UMa - pointer, upper
    (MERAK, "Merak", 53910, 2.37),  # Beta UMa - pointer, lower
    (PHECDA, "Phecda", 58001, 2.44),  # Gamma UMa - bowl corner
    (MEGREZ, "Megrez", 59774, 3.31),  # Delta UMa - bowl-handle junction
    (ALIOTH, "Alioth", 62956, 1.77),  # Epsilon UMa - brightest
    (MIZAR, "Mizar", 65378, 2.23),  # Zeta UMa - double star
    (ALCOR, "Alcor", 65477, 3.99),  # 80 UMa - companion to Mizar
    (ALKAID, "Alkaid", 67301, 1.86),  # Eta UMa - end of handle
]

# The bowl stars (without handle)
BOWL_STARS = [
    (DUBHE, "Dubhe"),
    (MERAK, "Merak"),
    (PHECDA, "Phecda"),
    (MEGREZ, "Megrez"),
]

# The handle stars
HANDLE_STARS = [
    (MEGREZ, "Megrez"),  # junction point
    (ALIOTH, "Alioth"),
    (MIZAR, "Mizar"),
    (ALKAID, "Alkaid"),
]

# Pointer stars (point to Polaris)
POINTER_STARS = [
    (MERAK, "Merak"),
    (DUBHE, "Dubhe"),
]


@pytest.mark.unit
class TestBigDipperStarsCatalog:
    """Test that all Big Dipper stars are in the catalog."""

    def test_all_big_dipper_stars_present(self):
        """Verify all 8 Big Dipper stars are in the STAR_CATALOG."""
        catalog_ids = {entry.id for entry in STAR_CATALOG}

        for star_id, name, _, _ in BIG_DIPPER_STARS:
            assert star_id in catalog_ids, (
                f"Big Dipper star {name} (ID={star_id}) not in catalog"
            )

    def test_big_dipper_star_count(self):
        """Verify we have exactly 8 Big Dipper stars defined (7 + Alcor)."""
        assert len(BIG_DIPPER_STARS) == 8, (
            "Should have exactly 8 Big Dipper stars (including Alcor)"
        )

    def test_alioth_is_brightest(self):
        """Verify Alioth is the brightest Big Dipper star (lowest magnitude)."""
        alioth_mag = None
        other_mags = []

        for star_id, name, _, mag in BIG_DIPPER_STARS:
            if star_id == ALIOTH:
                alioth_mag = mag
            else:
                other_mags.append((name, mag))

        assert alioth_mag is not None, "Alioth not found"
        for name, mag in other_mags:
            assert alioth_mag < mag, (
                f"Alioth ({alioth_mag}) should be brighter than {name} ({mag})"
            )

    @pytest.mark.parametrize("star_id,name,hip,mag", BIG_DIPPER_STARS)
    def test_star_has_proper_motion(self, star_id, name, hip, mag):
        """Verify each Big Dipper star has proper motion data."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        # Big Dipper stars should have proper motion data
        assert entry.data.pm_ra != 0 or entry.data.pm_dec != 0, (
            f"Star {name} should have proper motion data"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", BIG_DIPPER_STARS)
    def test_star_has_correct_hip_number(self, star_id, name, hip, mag):
        """Verify each Big Dipper star has correct Hipparcos number."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert entry.hip_number == hip, (
            f"Star {name} should have HIP {hip}, got {entry.hip_number}"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", BIG_DIPPER_STARS)
    def test_star_nomenclature_is_ursa_major(self, star_id, name, hip, mag):
        """Verify each Big Dipper star has UMa in nomenclature."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert "UMa" in entry.nomenclature, (
            f"Star {name} nomenclature should contain 'UMa', got {entry.nomenclature}"
        )


@pytest.mark.unit
class TestBigDipperStarsCalculation:
    """Test position calculations for Big Dipper stars."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch."""
        return 2451545.0

    @pytest.mark.parametrize("star_id,name,hip,mag", BIG_DIPPER_STARS)
    def test_star_position_reasonable(self, standard_jd, star_id, name, hip, mag):
        """Test each Big Dipper star returns a reasonable position."""
        pos, _ = ephem.calc_ut(standard_jd, star_id, 0)

        # Longitude should be 0-360
        assert 0 <= pos[0] < 360, f"{name} longitude {pos[0]}deg out of range"

        # Latitude should be reasonable (-90 to 90)
        assert -90 <= pos[1] <= 90, f"{name} latitude {pos[1]}deg out of range"

        # Fixed stars should be very distant
        assert pos[2] > 1000, f"{name} should have large distance, got {pos[2]}"

    def test_bowl_stars_form_quadrilateral(self, standard_jd):
        """Test that the bowl stars form a rough quadrilateral."""
        positions = []
        for star_id, name in BOWL_STARS:
            pos, _ = ephem.calc_ut(standard_jd, star_id, 0)
            positions.append((name, pos[0], pos[1]))

        # Bowl stars should be within a reasonable region
        lons = [p[1] for p in positions]
        lon_range = max(lons) - min(lons)
        assert lon_range < 20, f"Bowl longitude spread {lon_range}deg too large"

        lats = [p[2] for p in positions]
        lat_range = max(lats) - min(lats)
        assert lat_range < 15, f"Bowl latitude spread {lat_range}deg too large"

    def test_handle_stars_roughly_aligned(self, standard_jd):
        """Test that handle stars form a rough line."""
        positions = []
        for star_id, name in HANDLE_STARS:
            pos, _ = ephem.calc_ut(standard_jd, star_id, 0)
            positions.append((name, pos[0], pos[1]))

        # Handle stars should span along an arc
        lons = [p[1] for p in positions]
        lon_range = max(lons) - min(lons)
        # Handle spans from Megrez to Alkaid, about 20-25 degrees
        assert 10 < lon_range < 35, (
            f"Handle longitude spread {lon_range}deg not in expected range"
        )

    def test_mizar_alcor_close_together(self, standard_jd):
        """Test that Mizar and Alcor are very close (famous naked-eye double)."""
        mizar_pos, _ = ephem.calc_ut(standard_jd, MIZAR, 0)
        alcor_pos, _ = ephem.calc_ut(standard_jd, ALCOR, 0)

        # Calculate angular separation (simplified for small angles)
        lon_diff = abs(mizar_pos[0] - alcor_pos[0])
        lat_diff = abs(mizar_pos[1] - alcor_pos[1])
        separation = (lon_diff**2 + lat_diff**2) ** 0.5

        # Mizar and Alcor are about 12 arcminutes apart = 0.2 degrees
        assert separation < 0.5, (
            f"Mizar-Alcor separation {separation:.3f}deg too large, should be ~0.2deg"
        )

    def test_big_dipper_in_leo_virgo_region(self, standard_jd):
        """Test that Big Dipper stars are in the Leo/Virgo ecliptic region (~120-180 deg)."""
        # Dubhe is a good reference for the Big Dipper's position
        pos, _ = ephem.calc_ut(standard_jd, DUBHE, 0)

        # Should be in the Leo/Virgo region
        assert 130 < pos[0] < 180, (
            f"Dubhe should be in Leo/Virgo region, got {pos[0]:.1f} degrees"
        )

    def test_phecda_position(self, standard_jd):
        """Test that Phecda (new star) is positioned correctly in the bowl."""
        pos, _ = ephem.calc_ut(standard_jd, PHECDA, 0)

        # Phecda should be near 147 degrees ecliptic longitude (Leo region)
        assert 140 < pos[0] < 160, (
            f"Phecda longitude {pos[0]:.1f} out of expected range"
        )

        # Phecda has positive ecliptic latitude (north of ecliptic)
        assert 40 < pos[1] < 60, f"Phecda latitude {pos[1]:.2f} out of expected range"

    def test_megrez_position(self, standard_jd):
        """Test that Megrez (new star) is positioned correctly at bowl-handle junction."""
        pos, _ = ephem.calc_ut(standard_jd, MEGREZ, 0)

        # Megrez should be near 152 degrees ecliptic longitude
        assert 145 < pos[0] < 165, (
            f"Megrez longitude {pos[0]:.1f} out of expected range"
        )

        # Megrez has positive ecliptic latitude
        assert 45 < pos[1] < 65, f"Megrez latitude {pos[1]:.2f} out of expected range"


@pytest.mark.unit
class TestBigDipperStarsNameResolution:
    """Test name resolution for Big Dipper stars."""

    def test_resolve_dubhe(self):
        """Test Dubhe name resolution."""
        assert resolve_star_name("Dubhe") == DUBHE
        assert resolve_star_name("Alpha Ursae Majoris") == DUBHE
        assert resolve_star_name("Alpha UMa") == DUBHE

    def test_resolve_merak(self):
        """Test Merak name resolution."""
        assert resolve_star_name("Merak") == MERAK
        assert resolve_star_name("Beta Ursae Majoris") == MERAK
        assert resolve_star_name("Beta UMa") == MERAK

    def test_resolve_phecda(self):
        """Test Phecda name resolution."""
        assert resolve_star_name("Phecda") == PHECDA
        assert resolve_star_name("Gamma Ursae Majoris") == PHECDA
        assert resolve_star_name("Gamma UMa") == PHECDA
        assert resolve_star_name("Phad") == PHECDA
        assert resolve_star_name("Phekda") == PHECDA

    def test_resolve_megrez(self):
        """Test Megrez name resolution."""
        assert resolve_star_name("Megrez") == MEGREZ
        assert resolve_star_name("Delta Ursae Majoris") == MEGREZ
        assert resolve_star_name("Delta UMa") == MEGREZ
        assert resolve_star_name("Kaffa") == MEGREZ

    def test_resolve_alioth(self):
        """Test Alioth name resolution."""
        assert resolve_star_name("Alioth") == ALIOTH
        assert resolve_star_name("Epsilon Ursae Majoris") == ALIOTH
        assert resolve_star_name("Epsilon UMa") == ALIOTH

    def test_resolve_mizar(self):
        """Test Mizar name resolution."""
        assert resolve_star_name("Mizar") == MIZAR
        assert resolve_star_name("Zeta Ursae Majoris") == MIZAR
        assert resolve_star_name("Zeta UMa") == MIZAR
        assert resolve_star_name("Horse and Rider") == MIZAR

    def test_resolve_alcor(self):
        """Test Alcor name resolution."""
        assert resolve_star_name("Alcor") == ALCOR
        assert resolve_star_name("80 Ursae Majoris") == ALCOR
        assert resolve_star_name("80 UMa") == ALCOR
        assert resolve_star_name("Suha") == ALCOR

    def test_resolve_alkaid(self):
        """Test Alkaid name resolution."""
        assert resolve_star_name("Alkaid") == ALKAID
        assert resolve_star_name("Eta Ursae Majoris") == ALKAID
        assert resolve_star_name("Eta UMa") == ALKAID
        assert resolve_star_name("Benetnash") == ALKAID

    def test_canonical_names(self):
        """Test canonical name retrieval for Big Dipper stars."""
        assert get_canonical_star_name(DUBHE) == "Dubhe"
        assert get_canonical_star_name(MERAK) == "Merak"
        assert get_canonical_star_name(PHECDA) == "Phecda"
        assert get_canonical_star_name(MEGREZ) == "Megrez"
        assert get_canonical_star_name(ALIOTH) == "Alioth"
        assert get_canonical_star_name(MIZAR) == "Mizar"
        assert get_canonical_star_name(ALCOR) == "Alcor"
        assert get_canonical_star_name(ALKAID) == "Alkaid"
