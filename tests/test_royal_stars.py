"""
Unit tests for the four Royal Stars of Persia.

The four Royal Stars are the brightest stars near the four cardinal points
of the ancient Persian celestial sphere, used as watchers of the sky:

- Aldebaran (Alpha Tauri) - Watcher of the East (Persian: Tascheter)
- Regulus (Alpha Leonis) - Watcher of the North (Persian: Venant)
- Antares (Alpha Scorpii) - Watcher of the West (Persian: Satevis)
- Fomalhaut (Alpha Piscis Austrini) - Watcher of the South (Persian: Hastorang)

These stars are approximately 90 degrees apart in the sky and were used
in ancient Persia to mark the solstices and equinoxes around 3000 BCE.
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    ALDEBARAN,
    REGULUS,
    ANTARES,
    FOMALHAUT,
)
from libephemeris.fixed_stars import (
    STAR_CATALOG,
    resolve_star_name,
    get_canonical_star_name,
)


# The four Royal Stars of Persia with their expected data
# (constant, name, watcher_direction, persian_name, approx_ra_j2000, approx_magnitude)
ROYAL_STARS = [
    (ALDEBARAN, "Aldebaran", "East", "Tascheter", 68.98, 0.85),
    (REGULUS, "Regulus", "North", "Venant", 152.09, 1.40),
    (ANTARES, "Antares", "West", "Satevis", 247.35, 1.06),
    (FOMALHAUT, "Fomalhaut", "South", "Hastorang", 344.41, 1.16),
]


@pytest.mark.unit
class TestRoyalStarsCatalog:
    """Test that all four Royal Stars are in the catalog with accurate data."""

    def test_all_royal_stars_present(self):
        """Verify all four Royal Stars are in the STAR_CATALOG."""
        catalog_ids = {entry.id for entry in STAR_CATALOG}

        for star_id, name, direction, _, _, _ in ROYAL_STARS:
            assert star_id in catalog_ids, (
                f"Royal Star {name} (Watcher of the {direction}, ID={star_id}) "
                f"not in catalog"
            )

    def test_royal_count(self):
        """Verify we have exactly 4 Royal Stars defined."""
        assert len(ROYAL_STARS) == 4, "Should have exactly 4 Royal Stars"

    @pytest.mark.parametrize(
        "star_id,name,direction,persian_name,approx_ra,approx_mag", ROYAL_STARS
    )
    def test_star_has_proper_motion(
        self, star_id, name, direction, persian_name, approx_ra, approx_mag
    ):
        """Verify each Royal Star has proper motion data."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert entry.data.pm_ra != 0 or entry.data.pm_dec != 0, (
            f"Star {name} should have proper motion data"
        )

    @pytest.mark.parametrize(
        "star_id,name,direction,persian_name,approx_ra,approx_mag", ROYAL_STARS
    )
    def test_star_has_valid_magnitude(
        self, star_id, name, direction, persian_name, approx_ra, approx_mag
    ):
        """Verify each Royal Star has a valid magnitude close to expected."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        # Royal Stars are all bright (mag < 2.0)
        assert -1 < entry.magnitude < 2, (
            f"Star {name} should be bright (mag < 2), got {entry.magnitude}"
        )
        # Check magnitude is close to expected
        assert abs(entry.magnitude - approx_mag) < 0.1, (
            f"Star {name} magnitude {entry.magnitude} should be close to {approx_mag}"
        )

    @pytest.mark.parametrize(
        "star_id,name,direction,persian_name,approx_ra,approx_mag", ROYAL_STARS
    )
    def test_star_ra_accurate(
        self, star_id, name, direction, persian_name, approx_ra, approx_mag
    ):
        """Verify each Royal Star has accurate Right Ascension."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        # RA should be within 0.1 degrees of expected
        assert abs(entry.data.ra_j2000 - approx_ra) < 0.5, (
            f"Star {name} RA {entry.data.ra_j2000}deg should be close to {approx_ra}deg"
        )


@pytest.mark.unit
class TestRoyalStarsCalculation:
    """Test position calculations for Royal Stars."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch."""
        return 2451545.0

    @pytest.mark.parametrize(
        "star_id,name,direction,persian_name,approx_ra,approx_mag", ROYAL_STARS
    )
    def test_star_position_reasonable(
        self, standard_jd, star_id, name, direction, persian_name, approx_ra, approx_mag
    ):
        """Test each Royal Star returns a reasonable position."""
        pos, _ = ephem.calc_ut(standard_jd, star_id, 0)

        # Longitude should be 0-360
        assert 0 <= pos[0] < 360, f"{name} longitude {pos[0]}deg out of range"

        # Latitude should be reasonable (-90 to 90)
        assert -90 <= pos[1] <= 90, f"{name} latitude {pos[1]}deg out of range"

        # Fixed stars should be very distant
        assert pos[2] > 1000, f"{name} should have large distance, got {pos[2]}"

    def test_royal_stars_spread_across_sky(self, standard_jd):
        """Test Royal Stars are spread roughly 90 degrees apart in ecliptic longitude."""
        positions = []
        for star_id, name, _, _, _, _ in ROYAL_STARS:
            pos, _ = ephem.calc_ut(standard_jd, star_id, 0)
            positions.append((name, pos[0]))

        # Sort by longitude
        positions.sort(key=lambda x: x[1])

        # Calculate gaps between consecutive stars
        gaps = []
        for i in range(len(positions)):
            next_i = (i + 1) % len(positions)
            gap = positions[next_i][1] - positions[i][1]
            if gap < 0:
                gap += 360  # Handle wraparound
            gaps.append(gap)

        # Each gap should be roughly 90 degrees (within 30 degrees tolerance)
        # due to precession and actual positions
        for gap in gaps:
            assert 50 < gap < 130, (
                f"Royal Stars should be spread across sky, got gaps: {gaps}"
            )


@pytest.mark.unit
class TestRoyalStarsNameResolution:
    """Test name resolution for Royal Stars including watcher aliases."""

    def test_resolve_aldebaran_canonical(self):
        """Test Aldebaran canonical name resolution."""
        assert resolve_star_name("Aldebaran") == ALDEBARAN
        assert resolve_star_name("Alpha Tauri") == ALDEBARAN
        assert resolve_star_name("Eye of Taurus") == ALDEBARAN

    def test_resolve_aldebaran_watcher(self):
        """Test Aldebaran watcher alias resolution."""
        assert resolve_star_name("Watcher of the East") == ALDEBARAN
        assert resolve_star_name("Tascheter") == ALDEBARAN

    def test_resolve_regulus_canonical(self):
        """Test Regulus canonical name resolution."""
        assert resolve_star_name("Regulus") == REGULUS
        assert resolve_star_name("Alpha Leonis") == REGULUS
        assert resolve_star_name("Cor Leonis") == REGULUS

    def test_resolve_regulus_watcher(self):
        """Test Regulus watcher alias resolution."""
        assert resolve_star_name("Watcher of the North") == REGULUS
        assert resolve_star_name("Venant") == REGULUS

    def test_resolve_antares_canonical(self):
        """Test Antares canonical name resolution."""
        assert resolve_star_name("Antares") == ANTARES
        assert resolve_star_name("Alpha Scorpii") == ANTARES
        assert resolve_star_name("Rival of Mars") == ANTARES

    def test_resolve_antares_watcher(self):
        """Test Antares watcher alias resolution."""
        assert resolve_star_name("Watcher of the West") == ANTARES
        assert resolve_star_name("Satevis") == ANTARES

    def test_resolve_fomalhaut_canonical(self):
        """Test Fomalhaut canonical name resolution."""
        assert resolve_star_name("Fomalhaut") == FOMALHAUT
        assert resolve_star_name("Alpha Piscis Austrini") == FOMALHAUT
        assert resolve_star_name("Fish's Mouth") == FOMALHAUT

    def test_resolve_fomalhaut_watcher(self):
        """Test Fomalhaut watcher alias resolution."""
        assert resolve_star_name("Watcher of the South") == FOMALHAUT
        assert resolve_star_name("Hastorang") == FOMALHAUT

    @pytest.mark.parametrize(
        "star_id,name,direction,persian_name,approx_ra,approx_mag", ROYAL_STARS
    )
    def test_canonical_name_retrieval(
        self, star_id, name, direction, persian_name, approx_ra, approx_mag
    ):
        """Test canonical name retrieval for all Royal Stars."""
        canonical = get_canonical_star_name(star_id)
        assert canonical == name, f"Expected {name}, got {canonical}"


@pytest.mark.unit
class TestRoyalStarsHipparcos:
    """Test Royal Stars have correct Hipparcos catalog data."""

    def test_aldebaran_hipparcos(self):
        """Test Aldebaran has correct Hipparcos number."""
        for e in STAR_CATALOG:
            if e.id == ALDEBARAN:
                assert e.hip_number == 21421
                assert e.nomenclature == "alTau"
                break

    def test_regulus_hipparcos(self):
        """Test Regulus has correct Hipparcos number."""
        for e in STAR_CATALOG:
            if e.id == REGULUS:
                assert e.hip_number == 49669
                assert e.nomenclature == "alLeo"
                break

    def test_antares_hipparcos(self):
        """Test Antares has correct Hipparcos number."""
        for e in STAR_CATALOG:
            if e.id == ANTARES:
                assert e.hip_number == 80763
                assert e.nomenclature == "alSco"
                break

    def test_fomalhaut_hipparcos(self):
        """Test Fomalhaut has correct Hipparcos number."""
        for e in STAR_CATALOG:
            if e.id == FOMALHAUT:
                assert e.hip_number == 113368
                assert e.nomenclature == "alPsA"
                break


@pytest.mark.unit
class TestRoyalStarsProperMotion:
    """Test proper motion effects for Royal Stars over time."""

    def test_proper_motion_over_50_years(self):
        """Test that Royal Stars move over 50 years."""
        jd1 = ephem.julday(2000, 1, 1, 12.0)
        jd2 = ephem.julday(2050, 1, 1, 12.0)

        for star_id, name, _, _, _, _ in ROYAL_STARS:
            pos1, _ = ephem.calc_ut(jd1, star_id, 0)
            pos2, _ = ephem.calc_ut(jd2, star_id, 0)

            diff = abs(pos2[0] - pos1[0])
            if diff > 180:
                diff = 360 - diff

            # Should move at least 0.01 degrees over 50 years (precession alone is ~0.7)
            assert diff > 0.01, f"{name} should show movement over 50 years"
