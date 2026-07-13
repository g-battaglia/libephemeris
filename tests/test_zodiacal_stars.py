"""
Unit tests for zodiacal constellation bright stars.

This module tests the bright stars from each zodiacal constellation used in
astrological interpretation. The 12 zodiacal constellations lie along the
ecliptic and are associated with the astrological signs.

Zodiacal constellations and their key stars:
- Aries: Hamal, Sheratan, Mesarthim
- Taurus: Aldebaran, Elnath (+ Pleiades, Hyades clusters)
- Gemini: Pollux, Castor
- Cancer: Acubens, Tarf, Asellus Borealis, Asellus Australis
- Leo: Regulus, Denebola, Algieba, Zosma
- Virgo: Spica, Vindemiatrix
- Libra: Zubenelgenubi, Zubeneschamali
- Scorpius: Antares, Shaula, Sargas, Dschubba, Graffias, Lesath
- Sagittarius: Kaus Australis, Nunki, Kaus Media, Kaus Borealis, Ascella
- Capricornus: Deneb Algedi, Algedi, Dabih, Nashira
- Aquarius: Sadalsuud, Sadalmelik, Skat
- Pisces: Eta Piscium, Alrescha
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    # Aries
    HAMAL,
    SHERATAN,
    MESARTHIM,
    # Cancer
    ACUBENS,
    TARF,
    ASELLUS_BOREALIS,
    ASELLUS_AUSTRALIS,
    # Sagittarius
    KAUS_AUSTRALIS,
    NUNKI,
    KAUS_MEDIA,
    KAUS_BOREALIS,
    ASCELLA,
    # Capricornus
    DENEB_ALGEDI,
    ALGEDI,
    DABIH,
    NASHIRA,
    # Aquarius
    SADALSUUD,
    SADALMELIK,
    SKAT,
    # Pisces
    ETA_PISCIUM,
    ALRESCHA,
)
from libephemeris.fixed_stars import (
    STAR_CATALOG,
    resolve_star_name,
    get_canonical_star_name,
)


# ======== ARIES CONSTELLATION STARS ========
ARIES_STARS = [
    (HAMAL, "Hamal", 9884, 2.00),  # Alpha Ari - brightest
    (SHERATAN, "Sheratan", 8903, 2.64),  # Beta Ari
    (MESARTHIM, "Mesarthim", 8832, 3.88),  # Gamma Ari
]


# ======== CANCER CONSTELLATION STARS ========
CANCER_STARS = [
    (TARF, "Tarf", 40526, 3.52),  # Beta Cnc - brightest
    (ASELLUS_AUSTRALIS, "Asellus Australis", 42911, 3.94),  # Delta Cnc
    (ACUBENS, "Acubens", 44066, 4.25),  # Alpha Cnc
    (ASELLUS_BOREALIS, "Asellus Borealis", 42806, 4.66),  # Gamma Cnc
]


# ======== SAGITTARIUS CONSTELLATION STARS ========
SAGITTARIUS_STARS = [
    (KAUS_AUSTRALIS, "Kaus Australis", 90185, 1.85),  # Epsilon Sgr - brightest
    (NUNKI, "Nunki", 92855, 2.02),  # Sigma Sgr
    (ASCELLA, "Ascella", 93506, 2.59),  # Zeta Sgr
    (KAUS_MEDIA, "Kaus Media", 89931, 2.70),  # Delta Sgr
    (KAUS_BOREALIS, "Kaus Borealis", 90496, 2.81),  # Lambda Sgr
]


# ======== CAPRICORNUS CONSTELLATION STARS ========
CAPRICORNUS_STARS = [
    (DENEB_ALGEDI, "Deneb Algedi", 107556, 2.81),  # Delta Cap - brightest
    (DABIH, "Dabih", 100345, 3.08),  # Beta Cap
    (ALGEDI, "Algedi", 100064, 3.57),  # Alpha Cap
    (NASHIRA, "Nashira", 106985, 3.68),  # Gamma Cap
]


# ======== AQUARIUS CONSTELLATION STARS ========
AQUARIUS_STARS = [
    (SADALSUUD, "Sadalsuud", 106278, 2.87),  # Beta Aqr - brightest
    (SADALMELIK, "Sadalmelik", 109074, 2.96),  # Alpha Aqr
    (SKAT, "Skat", 113136, 3.27),  # Delta Aqr
]


# ======== PISCES CONSTELLATION STARS ========
PISCES_STARS = [
    (ETA_PISCIUM, "Eta Piscium", 7097, 3.62),  # Eta Psc - brightest
    (ALRESCHA, "Alrescha", 9487, 3.82),  # Alpha Psc
]


# Combine all new zodiacal stars for parametrized tests
ALL_NEW_ZODIACAL_STARS = (
    ARIES_STARS
    + CANCER_STARS
    + SAGITTARIUS_STARS
    + CAPRICORNUS_STARS
    + AQUARIUS_STARS
    + PISCES_STARS
)


@pytest.mark.unit
class TestZodiacalStarsCatalog:
    """Test that all zodiacal constellation stars are in the catalog."""

    def test_all_new_zodiacal_stars_present(self):
        """Verify all new zodiacal stars are in the STAR_CATALOG."""
        catalog_ids = {entry.id for entry in STAR_CATALOG}

        for star_id, name, _, _ in ALL_NEW_ZODIACAL_STARS:
            assert star_id in catalog_ids, (
                f"Zodiacal star {name} (ID={star_id}) not in catalog"
            )

    @pytest.mark.parametrize("star_id,name,hip,mag", ALL_NEW_ZODIACAL_STARS)
    def test_star_has_proper_motion(self, star_id, name, hip, mag):
        """Verify each zodiacal star has proper motion data."""
        entry = None
        for e in STAR_CATALOG:
            if e.id == star_id:
                entry = e
                break

        assert entry is not None, f"Star {name} not found in catalog"
        # Stars should have proper motion data
        assert entry.data.pm_ra != 0 or entry.data.pm_dec != 0, (
            f"Star {name} should have proper motion data"
        )

    @pytest.mark.parametrize("star_id,name,hip,mag", ALL_NEW_ZODIACAL_STARS)
    def test_star_has_correct_hip_number(self, star_id, name, hip, mag):
        """Verify each zodiacal star has correct Hipparcos number."""
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
class TestZodiacalStarsCalculation:
    """Test position calculations for zodiacal constellation stars."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch."""
        return 2451545.0

    @pytest.mark.parametrize("star_id,name,hip,mag", ALL_NEW_ZODIACAL_STARS)
    def test_star_position_reasonable(self, standard_jd, star_id, name, hip, mag):
        """Test each zodiacal star returns a reasonable position."""
        pos, _ = ephem.calc_ut(standard_jd, star_id, 0)

        # Longitude should be 0-360
        assert 0 <= pos[0] < 360, f"{name} longitude {pos[0]}deg out of range"

        # Latitude should be reasonable (-90 to 90)
        assert -90 <= pos[1] <= 90, f"{name} latitude {pos[1]}deg out of range"

        # Fixed stars should be very distant
        assert pos[2] > 1000, f"{name} should have large distance, got {pos[2]}"

    def test_aries_stars_in_aries_region(self, standard_jd):
        """Test that Aries stars are in the Aries/Taurus ecliptic region."""
        for star_id, name, _, _ in ARIES_STARS:
            pos, _ = ephem.calc_ut(standard_jd, star_id, 0)
            # Aries stars should be between ~0 and ~45 degrees ecliptic
            assert 0 < pos[0] < 50, (
                f"{name} should be in Aries region, got {pos[0]:.1f} degrees"
            )

    def test_sagittarius_stars_in_sagittarius_region(self, standard_jd):
        """Test that Sagittarius stars are in the Sagittarius ecliptic region."""
        for star_id, name, _, _ in SAGITTARIUS_STARS:
            pos, _ = ephem.calc_ut(standard_jd, star_id, 0)
            # Sagittarius stars should be between ~260 and ~300 degrees ecliptic
            assert 255 < pos[0] < 310, (
                f"{name} should be in Sagittarius region, got {pos[0]:.1f} degrees"
            )

    def test_kaus_australis_is_brightest_sagittarius(self, standard_jd):
        """Test that Kaus Australis is the brightest Sagittarius star."""
        kaus_entry = None
        for e in STAR_CATALOG:
            if e.id == KAUS_AUSTRALIS:
                kaus_entry = e
                break

        assert kaus_entry is not None, "Kaus Australis not found"

        for star_id, name, _, mag in SAGITTARIUS_STARS:
            if star_id != KAUS_AUSTRALIS:
                assert kaus_entry.magnitude < mag, (
                    f"Kaus Australis ({kaus_entry.magnitude}) "
                    f"should be brighter than {name} ({mag})"
                )

    def test_hamal_is_brightest_aries(self, standard_jd):
        """Test that Hamal is the brightest Aries star."""
        hamal_entry = None
        for e in STAR_CATALOG:
            if e.id == HAMAL:
                hamal_entry = e
                break

        assert hamal_entry is not None, "Hamal not found"

        for star_id, name, _, mag in ARIES_STARS:
            if star_id != HAMAL:
                assert hamal_entry.magnitude < mag, (
                    f"Hamal ({hamal_entry.magnitude}) "
                    f"should be brighter than {name} ({mag})"
                )


@pytest.mark.unit
class TestZodiacalStarsNameResolution:
    """Test name resolution for zodiacal constellation stars."""

    # Aries
    def test_resolve_hamal(self):
        """Test Hamal name resolution."""
        assert resolve_star_name("Hamal") == HAMAL
        assert resolve_star_name("Alpha Arietis") == HAMAL
        assert resolve_star_name("Alpha Ari") == HAMAL

    def test_resolve_sheratan(self):
        """Test Sheratan name resolution."""
        assert resolve_star_name("Sheratan") == SHERATAN
        assert resolve_star_name("Beta Arietis") == SHERATAN
        assert resolve_star_name("Beta Ari") == SHERATAN

    # Cancer
    def test_resolve_tarf(self):
        """Test Tarf name resolution."""
        assert resolve_star_name("Tarf") == TARF
        assert resolve_star_name("Beta Cancri") == TARF
        assert resolve_star_name("Beta Cnc") == TARF

    def test_resolve_asellus_stars(self):
        """Test Asellus (donkey) stars name resolution."""
        assert resolve_star_name("Asellus Borealis") == ASELLUS_BOREALIS
        assert resolve_star_name("Northern Donkey") == ASELLUS_BOREALIS
        assert resolve_star_name("Asellus Australis") == ASELLUS_AUSTRALIS
        assert resolve_star_name("Southern Donkey") == ASELLUS_AUSTRALIS

    # Sagittarius
    def test_resolve_kaus_australis(self):
        """Test Kaus Australis name resolution."""
        assert resolve_star_name("Kaus Australis") == KAUS_AUSTRALIS
        assert resolve_star_name("Epsilon Sagittarii") == KAUS_AUSTRALIS
        assert resolve_star_name("Epsilon Sgr") == KAUS_AUSTRALIS

    def test_resolve_nunki(self):
        """Test Nunki name resolution."""
        assert resolve_star_name("Nunki") == NUNKI
        assert resolve_star_name("Sigma Sagittarii") == NUNKI
        assert resolve_star_name("Sigma Sgr") == NUNKI

    # Capricornus
    def test_resolve_algedi(self):
        """Test Algedi name resolution."""
        assert resolve_star_name("Algedi") == ALGEDI
        assert resolve_star_name("Alpha Capricorni") == ALGEDI
        assert resolve_star_name("Alpha Cap") == ALGEDI
        assert resolve_star_name("Giedi") == ALGEDI

    def test_resolve_nashira(self):
        """Test Nashira name resolution."""
        assert resolve_star_name("Nashira") == NASHIRA
        assert resolve_star_name("Gamma Capricorni") == NASHIRA
        assert resolve_star_name("Fortunate One") == NASHIRA

    # Aquarius
    def test_resolve_sadalsuud(self):
        """Test Sadalsuud name resolution."""
        assert resolve_star_name("Sadalsuud") == SADALSUUD
        assert resolve_star_name("Beta Aquarii") == SADALSUUD
        assert resolve_star_name("Beta Aqr") == SADALSUUD

    def test_resolve_sadalmelik(self):
        """Test Sadalmelik name resolution."""
        assert resolve_star_name("Sadalmelik") == SADALMELIK
        assert resolve_star_name("Alpha Aquarii") == SADALMELIK
        assert resolve_star_name("Alpha Aqr") == SADALMELIK

    # Pisces
    def test_resolve_alrescha(self):
        """Test Alrescha name resolution."""
        assert resolve_star_name("Alrescha") == ALRESCHA
        assert resolve_star_name("Alpha Piscium") == ALRESCHA
        assert resolve_star_name("The Knot") == ALRESCHA

    def test_resolve_eta_piscium(self):
        """Test Eta Piscium name resolution."""
        assert resolve_star_name("Eta Piscium") == ETA_PISCIUM
        assert resolve_star_name("Eta Psc") == ETA_PISCIUM

    # Canonical names
    def test_canonical_names(self):
        """Test canonical name retrieval for new zodiacal stars."""
        assert get_canonical_star_name(HAMAL) == "Hamal"
        assert get_canonical_star_name(TARF) == "Tarf"
        assert get_canonical_star_name(KAUS_AUSTRALIS) == "Kaus Australis"
        assert get_canonical_star_name(SADALSUUD) == "Sadalsuud"
        assert get_canonical_star_name(ALRESCHA) == "Alrescha"


@pytest.mark.unit
class TestZodiacalConstellationCoverage:
    """Test that all 12 zodiacal constellations have stars."""

    def test_all_zodiacal_constellations_have_stars(self):
        """Verify each zodiacal constellation has at least one star."""
        # Map constellation abbreviations to expected nomenclature patterns
        zodiacal_constellations = {
            "Aries": "Ari",
            "Taurus": "Tau",
            "Gemini": "Gem",
            "Cancer": "Cnc",
            "Leo": "Leo",
            "Virgo": "Vir",
            "Libra": "Lib",
            "Scorpius": "Sco",
            "Sagittarius": "Sgr",
            "Capricornus": "Cap",
            "Aquarius": "Aqr",
            "Pisces": "Psc",
        }

        for constellation, abbrev in zodiacal_constellations.items():
            found = False
            for entry in STAR_CATALOG:
                if abbrev in entry.nomenclature:
                    found = True
                    break
            assert found, (
                f"Zodiacal constellation {constellation} ({abbrev}) "
                "has no stars in catalog"
            )

    def test_each_zodiacal_constellation_has_minimum_stars(self):
        """Verify each zodiacal constellation has at least 2 bright stars."""
        zodiacal_constellations = {
            "Aries": ("Ari", 3),  # Hamal, Sheratan, Mesarthim
            "Taurus": ("Tau", 2),  # Aldebaran, Elnath (+ clusters)
            "Gemini": ("Gem", 2),  # Pollux, Castor
            "Cancer": ("Cnc", 4),  # Tarf, Acubens, Asellus Borealis/Australis
            "Leo": ("Leo", 4),  # Regulus, Denebola, Algieba, Zosma
            "Virgo": ("Vir", 2),  # Spica, Vindemiatrix
            "Libra": ("Lib", 2),  # Zubenelgenubi, Zubeneschamali
            "Scorpius": ("Sco", 6),  # Antares, Shaula, etc.
            "Sagittarius": ("Sgr", 5),  # Kaus Australis, Nunki, etc.
            "Capricornus": ("Cap", 4),  # Deneb Algedi, Algedi, Dabih, Nashira
            "Aquarius": ("Aqr", 3),  # Sadalsuud, Sadalmelik, Skat
            "Pisces": ("Psc", 2),  # Eta Piscium, Alrescha
        }

        for constellation, (abbrev, min_count) in zodiacal_constellations.items():
            count = sum(1 for entry in STAR_CATALOG if abbrev in entry.nomenclature)
            assert count >= min_count, (
                f"Zodiacal constellation {constellation} should have at least "
                f"{min_count} stars, found {count}"
            )


# Regression pins for the three catalog entries that were corrupted in an
# earlier hand-maintained catalog and fixed in WS6c/WS9a (see
# REVIEW-2026-06-10.md item C9): Tarf carried delta-Cnc's HIP/position,
# Eta Piscium's RA was off by ~4.4°, and Alrescha mixed alpha-Psc with
# eta-Psc's proper motion and the wrong HIP. The catalog is now generated
# from van Leeuwen 2007 (CDS I/311) by scripts/build_star_catalog_v2.py;
# these pins keep a regeneration from silently re-introducing the swap.
# Reference values: van Leeuwen, F. (2007), A&A 474, 653 (CDS I/311/hip2),
# positions propagated J1991.25 -> J2000.0.
C9_CATALOG_PINS = [
    # name,           HIP,    RA_J2000_deg,   Dec_J2000_deg, nomenclature
    ("Tarf", 40526, 124.128837533, 9.185543869, "beCnc"),
    ("Eta Piscium", 7097, 22.870876074, 15.345824593, "etPsc"),
    ("Alrescha", 9487, 30.511748813, 2.763761377, "alPsc"),
]


class TestC9CatalogRegression:
    """Pin the three formerly-corrupt entries (REVIEW-2026-06-10 C9)."""

    @pytest.mark.parametrize("name,hip,ra,dec,nomen", C9_CATALOG_PINS)
    def test_c9_entry_astrometry(self, name, hip, ra, dec, nomen):
        matches = [e for e in STAR_CATALOG if e.name == name]
        assert len(matches) == 1, f"expected exactly one {name!r} entry"
        entry = matches[0]
        assert entry.hip_number == hip, (
            f"{name}: HIP {entry.hip_number} != {hip} (C9 swap regression)"
        )
        assert entry.nomenclature == nomen, (
            f"{name}: nomenclature {entry.nomenclature!r} != {nomen!r}"
        )
        assert entry.data.ra_j2000 == pytest.approx(ra, abs=1e-6), (
            f"{name}: RA {entry.data.ra_j2000} != {ra}"
        )
        assert entry.data.dec_j2000 == pytest.approx(dec, abs=1e-6), (
            f"{name}: Dec {entry.data.dec_j2000} != {dec}"
        )
