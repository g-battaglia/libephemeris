"""
Unit tests for Centaurus constellation stars.

Centaurus is a large southern constellation containing some of the
brightest stars in the sky, including:

- Rigil Kentaurus (Alpha Centauri) - closest star system to the Sun
- Hadar (Beta Centauri) - second brightest in the constellation
- Menkent (Theta Centauri) - third brightest
- Muhlifain (Gamma Centauri) - a double star
- Epsilon Centauri - a blue giant
- Eta Centauri - a Be star
- Zeta Centauri - a spectroscopic binary

The two brightest stars, Alpha and Beta Centauri, are often called
the "Pointer Stars" as they point toward the Southern Cross.
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    RIGIL_KENT,
    HADAR,
    MENKENT,
    MUHLIFAIN,
    EPSILON_CENTAURI,
    ETA_CENTAURI,
    ZETA_CENTAURI,
)
from libephemeris.fixed_stars import (
    STAR_CATALOG,
    resolve_star_name,
    get_canonical_star_name,
)


# The 7 bright Centaurus stars
# (constant, name, hip_number, magnitude)
CENTAURUS_STARS = [
    (RIGIL_KENT, "Rigil Kentaurus", 71683, -0.27),  # Alpha Cen - brightest
    (HADAR, "Hadar", 68702, 0.61),  # Beta Cen - second brightest
    (MENKENT, "Menkent", 68933, 2.06),  # Theta Cen - third brightest
    (MUHLIFAIN, "Muhlifain", 61932, 2.17),  # Gamma Cen
    (EPSILON_CENTAURI, "Epsilon Centauri", 66657, 2.30),  # Epsilon Cen
    (ETA_CENTAURI, "Eta Centauri", 71352, 2.31),  # Eta Cen
    (ZETA_CENTAURI, "Zeta Centauri", 68002, 2.55),  # Zeta Cen
]


@pytest.mark.unit
class TestCentaurusStarsCatalog:
    """Test that all Centaurus stars are in the catalog."""

    def test_all_centaurus_stars_present(self):
        """Verify all 7 Centaurus stars are in the STAR_CATALOG."""
        catalog_ids = {entry.id for entry in STAR_CATALOG}

        for star_id, name, _, _ in CENTAURUS_STARS:
            assert star_id in catalog_ids, (
                f"Centaurus star {name} (ID={star_id}) not in catalog"
            )

    def test_centaurus_star_count(self):
        """Verify we have exactly 7 Centaurus stars defined."""
        assert len(CENTAURUS_STARS) == 7, "Should have exactly 7 Centaurus stars"

    def test_rigil_kent_is_brightest(self):
        """Verify Rigil Kentaurus is the brightest Centaurus star (lowest magnitude)."""
        rigil_mag = None
        other_mags = []

        for star_id, name, _, mag in CENTAURUS_STARS:
            if star_id == RIGIL_KENT:
                rigil_mag = mag
            else:
                other_mags.append((name, mag))

        assert rigil_mag is not None, "Rigil Kentaurus not found"
        for name, mag in other_mags:
            assert rigil_mag < mag, (
                f"Rigil Kentaurus ({rigil_mag}) should be brighter than {name} ({mag})"
            )

    @pytest.mark.parametrize("star_id,name,hip,mag", CENTAURUS_STARS)
    def test_star_has_proper_motion(self, star_id, name, hip, mag):
        """Verify each Centaurus star has proper motion data."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        # Centaurus stars should have proper motion data
        assert entry.data.pm_ra != 0 or entry.data.pm_dec != 0, (
            f"Star {name} should have proper motion data"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", CENTAURUS_STARS)
    def test_star_has_correct_hip_number(self, star_id, name, hip, mag):
        """Verify each Centaurus star has correct Hipparcos number."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert entry.hip_number == hip, (
            f"Star {name} should have HIP {hip}, got {entry.hip_number}"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", CENTAURUS_STARS)
    def test_star_nomenclature_is_centaurus(self, star_id, name, hip, mag):
        """Verify each Centaurus star has Cen in nomenclature."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        assert "Cen" in entry.nomenclature, (
            f"Star {name} nomenclature should contain 'Cen', got {entry.nomenclature}"
        )


@pytest.mark.unit
class TestCentaurusStarsCalculation:
    """Test position calculations for Centaurus stars."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch."""
        return 2451545.0

    @pytest.mark.parametrize("star_id,name,hip,mag", CENTAURUS_STARS)
    def test_star_position_reasonable(self, standard_jd, star_id, name, hip, mag):
        """Test each Centaurus star returns a reasonable position."""
        pos, _ = ephem.calc_ut(standard_jd, star_id, 0)

        # Longitude should be 0-360
        assert 0 <= pos[0] < 360, f"{name} longitude {pos[0]}deg out of range"

        # Latitude should be reasonable (-90 to 90)
        assert -90 <= pos[1] <= 90, f"{name} latitude {pos[1]}deg out of range"

        # Fixed stars should be very distant
        assert pos[2] > 1000, f"{name} should have large distance, got {pos[2]}"

    def test_centaurus_stars_in_southern_sky(self, standard_jd):
        """Test that all Centaurus stars have negative ecliptic latitude (southern sky)."""
        for star_id, name, _, _ in CENTAURUS_STARS:
            pos, _ = ephem.calc_ut(standard_jd, star_id, 0)
            # Centaurus stars are all in the southern sky (negative ecliptic latitude)
            assert pos[1] < 0, (
                f"{name} should be in southern sky, got latitude {pos[1]:.2f}"
            )

    def test_centaurus_stars_in_expected_region(self, standard_jd):
        """Test that Centaurus stars are in expected ecliptic region."""
        for star_id, name, _, _ in CENTAURUS_STARS:
            pos, _ = ephem.calc_ut(standard_jd, star_id, 0)
            # Centaurus stars should be roughly in Libra/Scorpio ecliptic region
            assert 170 < pos[0] < 260, (
                f"{name} longitude {pos[0]:.1f} out of expected range"
            )

    def test_rigil_kent_position(self, standard_jd):
        """Test Alpha Centauri position is accurate."""
        pos, _ = ephem.calc_ut(standard_jd, RIGIL_KENT, 0)

        # Alpha Centauri should be around 29 Scorpio ecliptic longitude
        assert 205 < pos[0] < 240, (
            f"Rigil Kentaurus longitude {pos[0]:.1f} out of expected range"
        )

        # Deep southern ecliptic latitude
        assert -55 < pos[1] < -35, (
            f"Rigil Kentaurus latitude {pos[1]:.2f} out of expected range"
        )

    def test_hadar_position(self, standard_jd):
        """Test Beta Centauri position is accurate."""
        pos, _ = ephem.calc_ut(standard_jd, HADAR, 0)

        # Beta Centauri should be near Alpha Centauri in the sky
        assert 195 < pos[0] < 240, f"Hadar longitude {pos[0]:.1f} out of expected range"

        # Deep southern ecliptic latitude
        assert -55 < pos[1] < -35, f"Hadar latitude {pos[1]:.2f} out of expected range"

    def test_pointer_stars_proximity(self, standard_jd):
        """Test that Alpha and Beta Centauri (pointer stars) are close together."""
        pos_alpha, _ = ephem.calc_ut(standard_jd, RIGIL_KENT, 0)
        pos_beta, _ = ephem.calc_ut(standard_jd, HADAR, 0)

        # The pointer stars should be within about 10 degrees of each other
        lon_diff = abs(pos_alpha[0] - pos_beta[0])
        lat_diff = abs(pos_alpha[1] - pos_beta[1])

        assert lon_diff < 15, (
            f"Pointer stars longitude separation {lon_diff:.1f} too large"
        )
        assert lat_diff < 10, (
            f"Pointer stars latitude separation {lat_diff:.1f} too large"
        )


@pytest.mark.unit
class TestCentaurusStarsNameResolution:
    """Test name resolution for Centaurus stars."""

    def test_resolve_rigil_kent(self):
        """Test Rigil Kentaurus name resolution."""
        assert resolve_star_name("Rigil Kentaurus") == RIGIL_KENT
        assert resolve_star_name("Alpha Centauri") == RIGIL_KENT
        assert resolve_star_name("Alpha Cen") == RIGIL_KENT

    def test_resolve_toliman_is_alpha_cen_b(self):
        """Toliman is the IAU name of alpha Cen B (HIP 71681), a
        separate catalog entry from Rigil Kentaurus (alpha Cen A)."""
        from libephemeris.fixed_stars import STAR_CATALOG

        toliman_id = resolve_star_name("Toliman")
        entry = next(e for e in STAR_CATALOG if e.id == toliman_id)
        assert entry.hip_number == 71681
        assert toliman_id != RIGIL_KENT

    def test_resolve_hadar(self):
        """Test Hadar name resolution."""
        assert resolve_star_name("Hadar") == HADAR
        assert resolve_star_name("Beta Centauri") == HADAR
        assert resolve_star_name("Beta Cen") == HADAR
        assert resolve_star_name("Agena") == HADAR

    def test_resolve_menkent(self):
        """Test Menkent name resolution."""
        assert resolve_star_name("Menkent") == MENKENT
        assert resolve_star_name("Theta Centauri") == MENKENT
        assert resolve_star_name("Theta Cen") == MENKENT

    def test_resolve_muhlifain(self):
        """Test Muhlifain name resolution."""
        assert resolve_star_name("Muhlifain") == MUHLIFAIN
        assert resolve_star_name("Gamma Centauri") == MUHLIFAIN
        assert resolve_star_name("Gamma Cen") == MUHLIFAIN

    def test_resolve_epsilon_centauri(self):
        """Test Epsilon Centauri name resolution."""
        assert resolve_star_name("Epsilon Centauri") == EPSILON_CENTAURI
        assert resolve_star_name("Epsilon Cen") == EPSILON_CENTAURI

    def test_resolve_eta_centauri(self):
        """Test Eta Centauri name resolution."""
        assert resolve_star_name("Eta Centauri") == ETA_CENTAURI
        assert resolve_star_name("Eta Cen") == ETA_CENTAURI

    def test_resolve_zeta_centauri(self):
        """Test Zeta Centauri name resolution."""
        assert resolve_star_name("Zeta Centauri") == ZETA_CENTAURI
        assert resolve_star_name("Zeta Cen") == ZETA_CENTAURI

    def test_canonical_names(self):
        """Test canonical name retrieval for Centaurus stars."""
        assert get_canonical_star_name(RIGIL_KENT) == "Rigil Kentaurus"
        assert get_canonical_star_name(HADAR) == "Hadar"
        assert get_canonical_star_name(MENKENT) == "Menkent"
        assert get_canonical_star_name(MUHLIFAIN) == "Muhlifain"
        assert get_canonical_star_name(EPSILON_CENTAURI) == "Epsilon Centauri"
        assert get_canonical_star_name(ETA_CENTAURI) == "Eta Centauri"
        assert get_canonical_star_name(ZETA_CENTAURI) == "Zeta Centauri"
