"""
Unit tests for Bayer designation star search.

Tests the ability to search for stars by their Bayer designations using
full Greek letter names, e.g., "Alpha Leonis", "Beta Persei", "Gamma Virginis".

This feature is required for fixstar2 compatibility with designation search.
"""

import pytest

from libephemeris.exceptions import Error
from libephemeris.fixed_stars import (
    _parse_bayer_designation,
    resolve_star_name,
    _resolve_star2,
    fixstar2_ut,
    fixstar2,
    GREEK_LETTER_ABBREV,
    CONSTELLATION_ABBREV,
)
from libephemeris.constants import (
    REGULUS,
    SPICA_STAR,
    ALGOL,
    SIRIUS,
    ALDEBARAN,
    ANTARES,
    VEGA,
    BETELGEUSE,
    RIGEL,
    BELLATRIX,
    DENEB,
    ALTAIR,
)


@pytest.mark.unit
class TestBayerDesignationParsing:
    """Tests for _parse_bayer_designation function."""

    def test_parse_alpha_leonis(self):
        """Alpha Leonis should parse to alLeo."""
        result = _parse_bayer_designation("Alpha Leonis")
        assert result == "alLeo"

    def test_parse_beta_persei(self):
        """Beta Persei should parse to bePer."""
        result = _parse_bayer_designation("Beta Persei")
        assert result == "bePer"

    def test_parse_gamma_virginis(self):
        """Gamma Virginis should parse to gaVir."""
        result = _parse_bayer_designation("Gamma Virginis")
        assert result == "gaVir"

    def test_parse_alpha_orionis(self):
        """Alpha Orionis should parse to alOri."""
        result = _parse_bayer_designation("Alpha Orionis")
        assert result == "alOri"

    def test_parse_epsilon_virginis(self):
        """Epsilon Virginis should parse to epVir."""
        result = _parse_bayer_designation("Epsilon Virginis")
        assert result == "epVir"

    def test_parse_case_insensitive(self):
        """Parsing should be case-insensitive."""
        result1 = _parse_bayer_designation("ALPHA LEONIS")
        result2 = _parse_bayer_designation("alpha leonis")
        result3 = _parse_bayer_designation("Alpha Leonis")
        assert result1 == result2 == result3 == "alLeo"

    def test_parse_nominative_form(self):
        """Should support nominative constellation names like 'Alpha Leo'."""
        result = _parse_bayer_designation("Alpha Leo")
        assert result == "alLeo"

    def test_parse_alpha_orion(self):
        """Should support nominative form 'Alpha Orion'."""
        result = _parse_bayer_designation("Alpha Orion")
        assert result == "alOri"

    def test_parse_invalid_greek_letter(self):
        """Invalid Greek letter should return None."""
        result = _parse_bayer_designation("Xyzzy Leonis")
        assert result is None

    def test_parse_invalid_constellation(self):
        """Invalid constellation should return None."""
        result = _parse_bayer_designation("Alpha Xyzzy")
        assert result is None

    def test_parse_empty_string(self):
        """Empty string should return None."""
        result = _parse_bayer_designation("")
        assert result is None

    def test_parse_single_word(self):
        """Single word (missing constellation) should return None."""
        result = _parse_bayer_designation("Alpha")
        assert result is None

    def test_parse_multi_word_constellation(self):
        """Should support multi-word constellations like 'Ursa Major'."""
        result = _parse_bayer_designation("Alpha Ursae Majoris")
        assert result == "alUMa"

    def test_parse_alpha_canis_majoris(self):
        """Alpha Canis Majoris (Sirius) should parse to alCMa."""
        result = _parse_bayer_designation("Alpha Canis Majoris")
        assert result == "alCMa"

    def test_parse_lambda_scorpii(self):
        """Lambda Scorpii (Shaula) should parse to laSco."""
        result = _parse_bayer_designation("Lambda Scorpii")
        assert result == "laSco"

    def test_parse_theta_centauri(self):
        """Theta Centauri should parse to thCen."""
        result = _parse_bayer_designation("Theta Centauri")
        assert result == "thCen"


@pytest.mark.unit
class TestGreekLetterMappings:
    """Tests for Greek letter abbreviation mappings."""

    def test_all_greek_letters_present(self):
        """All 24 Greek letters should be mapped."""
        expected_letters = {
            "ALPHA",
            "BETA",
            "GAMMA",
            "DELTA",
            "EPSILON",
            "ZETA",
            "ETA",
            "THETA",
            "IOTA",
            "KAPPA",
            "LAMBDA",
            "MU",
            "NU",
            "XI",
            "OMICRON",
            "PI",
            "RHO",
            "SIGMA",
            "TAU",
            "UPSILON",
            "PHI",
            "CHI",
            "PSI",
            "OMEGA",
        }
        assert expected_letters.issubset(set(GREEK_LETTER_ABBREV.keys()))

    def test_abbreviations_are_two_or_three_chars(self):
        """Abbreviations are 2 chars, except omicron's 3-char 'omi'.

        The reference catalog distinguishes omicron ('omi') from omega
        ('om'), so a single 2-char rule cannot hold for both.
        """
        for letter, abbrev in GREEK_LETTER_ABBREV.items():
            assert 2 <= len(abbrev) <= 3, (
                f"{letter} abbreviation '{abbrev}' is not 2-3 chars"
            )
        assert GREEK_LETTER_ABBREV["OMICRON"] == "omi"
        assert GREEK_LETTER_ABBREV["OMEGA"] == "om"


@pytest.mark.unit
class TestConstellationMappings:
    """Tests for constellation abbreviation mappings."""

    def test_zodiac_constellations_present(self):
        """All 12 zodiac constellations should be mapped (genitive forms)."""
        zodiac_genitive = [
            "ARIETIS",
            "TAURI",
            "GEMINORUM",
            "CANCRI",
            "LEONIS",
            "VIRGINIS",
            "LIBRAE",
            "SCORPII",
            "SAGITTARII",
            "CAPRICORNI",
            "AQUARII",
            "PISCIUM",
        ]
        for constellation in zodiac_genitive:
            assert constellation in CONSTELLATION_ABBREV, f"Missing: {constellation}"

    def test_common_constellations_present(self):
        """Common constellations should be mapped."""
        common = [
            "ORIONIS",
            "CYGNI",
            "LYRAE",
            "AQUILAE",
            "CENTAURI",
            "CRUCIS",
            "DRACONIS",
            "URSAE MAJORIS",
        ]
        for constellation in common:
            assert constellation in CONSTELLATION_ABBREV, f"Missing: {constellation}"


@pytest.mark.unit
class TestResolveStarNameWithBayer:
    """Tests for resolve_star_name with Bayer designations."""

    def test_resolve_alpha_leonis(self):
        """Alpha Leonis should resolve to REGULUS."""
        result = resolve_star_name("Alpha Leonis")
        assert result == REGULUS

    def test_resolve_beta_persei(self):
        """Beta Persei should resolve to ALGOL."""
        result = resolve_star_name("Beta Persei")
        assert result == ALGOL

    def test_resolve_alpha_virginis(self):
        """Alpha Virginis should resolve to SPICA_STAR."""
        result = resolve_star_name("Alpha Virginis")
        assert result == SPICA_STAR

    def test_resolve_alpha_scorpii(self):
        """Alpha Scorpii should resolve to ANTARES."""
        result = resolve_star_name("Alpha Scorpii")
        assert result == ANTARES

    def test_resolve_alpha_tauri(self):
        """Alpha Tauri should resolve to ALDEBARAN."""
        result = resolve_star_name("Alpha Tauri")
        assert result == ALDEBARAN

    def test_resolve_alpha_lyrae(self):
        """Alpha Lyrae should resolve to VEGA."""
        result = resolve_star_name("Alpha Lyrae")
        assert result == VEGA

    def test_resolve_alpha_orionis(self):
        """Alpha Orionis should resolve to BETELGEUSE."""
        result = resolve_star_name("Alpha Orionis")
        assert result == BETELGEUSE

    def test_resolve_beta_orionis(self):
        """Beta Orionis should resolve to RIGEL."""
        result = resolve_star_name("Beta Orionis")
        assert result == RIGEL

    def test_resolve_gamma_orionis(self):
        """Gamma Orionis should resolve to BELLATRIX."""
        result = resolve_star_name("Gamma Orionis")
        assert result == BELLATRIX

    def test_resolve_alpha_cygni(self):
        """Alpha Cygni should resolve to DENEB."""
        result = resolve_star_name("Alpha Cygni")
        assert result == DENEB

    def test_resolve_alpha_aquilae(self):
        """Alpha Aquilae should resolve to ALTAIR."""
        result = resolve_star_name("Alpha Aquilae")
        assert result == ALTAIR

    def test_resolve_alpha_canis_majoris(self):
        """Alpha Canis Majoris should resolve to SIRIUS."""
        result = resolve_star_name("Alpha Canis Majoris")
        assert result == SIRIUS


@pytest.mark.unit
class TestResolveStar2WithBayer:
    """Tests for _resolve_star2 with Bayer designations."""

    def test_resolve2_alpha_leonis(self):
        """Alpha Leonis should resolve to Regulus entry."""
        entry, err = _resolve_star2("Alpha Leonis")
        assert err is None
        assert entry is not None
        assert entry.name == "Regulus"
        assert entry.nomenclature == "alLeo"

    def test_resolve2_beta_persei(self):
        """Beta Persei should resolve to Algol entry."""
        entry, err = _resolve_star2("Beta Persei")
        assert err is None
        assert entry is not None
        assert entry.name == "Algol"
        assert entry.nomenclature == "bePer"

    def test_resolve2_gamma_orionis(self):
        """Gamma Orionis should resolve to Bellatrix entry."""
        entry, err = _resolve_star2("Gamma Orionis")
        assert err is None
        assert entry is not None
        assert entry.name == "Bellatrix"
        assert entry.nomenclature == "gaOri"

    def test_resolve2_epsilon_virginis(self):
        """Epsilon Virginis should resolve to Vindemiatrix entry."""
        entry, err = _resolve_star2("Epsilon Virginis")
        assert err is None
        assert entry is not None
        assert entry.name == "Vindemiatrix"
        assert entry.nomenclature == "epVir"

    def test_resolve2_star_not_in_catalog(self):
        """Bayer designation for a star not in the catalog returns an error.

        (Porrima/Gamma Virginis joined the catalog with the 1447-star
        expansion, so a constellation/letter pair with no star at all is
        used instead.)
        """
        entry, err = _resolve_star2("Omega Crucis")  # no such Bayer star
        assert entry is None
        assert "no fixed star matches" in err.lower()


@pytest.mark.unit
class TestSweFixstar2WithBayer:
    """Tests for fixstar2 and fixstar2_ut with Bayer designations."""

    def test_swe_fixstar2_ut_alpha_leonis(self):
        """fixstar2_ut should find Regulus via 'Alpha Leonis'."""
        jd = 2451545.0  # J2000
        pos, name, flag = fixstar2_ut("Alpha Leonis", jd, 0)

        assert "Regulus" in name
        assert "alLeo" in name
        # Regulus at ~149-150° at J2000
        assert 149 < pos[0] < 151

    def test_swe_fixstar2_ut_beta_persei(self):
        """fixstar2_ut should find Algol via 'Beta Persei'."""
        jd = 2451545.0
        pos, name, flag = fixstar2_ut("Beta Persei", jd, 0)

        assert "Algol" in name
        assert "bePer" in name
        # Algol at ~55-57° at J2000
        assert 55 < pos[0] < 58

    def test_swe_fixstar2_alpha_virginis(self):
        """fixstar2 (TT) should find Spica via 'Alpha Virginis'."""
        jd = 2451545.0
        pos, name, flag = fixstar2("Alpha Virginis", jd, 0)

        assert "Spica" in name
        assert "alVir" in name
        # Spica at ~203-204° at J2000
        assert 203 < pos[0] < 205

    def test_swe_fixstar2_gamma_orionis(self):
        """fixstar2 should find Bellatrix via 'Gamma Orionis'."""
        jd = 2451545.0
        pos, name, flag = fixstar2("Gamma Orionis", jd, 0)

        assert "Bellatrix" in name
        assert "gaOri" in name

    def test_swe_fixstar2_ut_star_not_found(self):
        """fixstar2_ut should raise error for star not in catalog."""
        jd = 2451545.0
        with pytest.raises(Error):
            fixstar2_ut("NoSuchStarXYZ", jd, 0)


@pytest.mark.unit
class TestBayerDesignationEdgeCases:
    """Edge case tests for Bayer designation parsing."""

    def test_extra_whitespace_handling(self):
        """Should handle extra whitespace in designation."""
        result = _parse_bayer_designation("  Alpha   Leonis  ")
        assert result == "alLeo"

    def test_mixed_case_handling(self):
        """Should handle mixed case in designation."""
        result = _parse_bayer_designation("aLpHa LeOnIs")
        assert result == "alLeo"

    def test_resolve_with_mixed_case(self):
        """resolve_star_name should work with mixed case Bayer designation."""
        result = resolve_star_name("ALPHA LEONIS")
        assert result == REGULUS

        result = resolve_star_name("alpha leonis")
        assert result == REGULUS

    def test_all_greek_letters_parse(self):
        """All Greek letters should parse correctly with a known constellation."""
        greek_letters = [
            "Alpha",
            "Beta",
            "Gamma",
            "Delta",
            "Epsilon",
            "Zeta",
            "Eta",
            "Theta",
            "Iota",
            "Kappa",
            "Lambda",
            "Mu",
            "Nu",
            "Xi",
            "Omicron",
            "Pi",
            "Rho",
            "Sigma",
            "Tau",
            "Upsilon",
            "Phi",
            "Chi",
            "Psi",
            "Omega",
        ]
        for letter in greek_letters:
            result = _parse_bayer_designation(f"{letter} Leonis")
            assert result is not None, f"Failed to parse '{letter} Leonis'"
            assert result.endswith("Leo"), (
                f"'{letter} Leonis' -> '{result}' should end with 'Leo'"
            )
