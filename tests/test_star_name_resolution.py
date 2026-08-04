"""
Unit tests for star name resolution with aliases and wildcards.

Tests the reference-API-compatible star name resolution system including:
- STAR_ALIASES dictionary with alternative names
- resolve_star_name function with prefix/fuzzy matching
- fixstar_ut with canonical star name return
"""

import pytest
from libephemeris.fixed_stars import (
    STAR_ALIASES,
    STAR_CATALOG,
    resolve_star_name,
    get_canonical_star_name,
    fixstar_ut,
    fixstar,
)
from libephemeris.constants import (
    REGULUS,
    SPICA_STAR,
    ALGOL,
    SIRIUS,
    ALDEBARAN,
    ANTARES,
    VEGA,
    POLARIS,
    BETELGEUSE,
    RIGEL,
    PROCYON,
    ARCTURUS,
    DENEB,
    POLLUX,
    CASTOR,
    ALTAIR,
)


@pytest.mark.unit
class TestStarAliasesDictionary:
    """Tests for the STAR_ALIASES dictionary structure."""

    def test_star_aliases_has_minimum_entries(self):
        """Verify STAR_ALIASES has at least 300 entries as required."""
        assert len(STAR_ALIASES) >= 300, (
            f"STAR_ALIASES should have at least 300 entries, "
            f"but has {len(STAR_ALIASES)}"
        )

    def test_star_aliases_values_are_valid_star_ids(self):
        """Verify all alias values map to valid star IDs in catalog."""
        valid_star_ids = {entry.id for entry in STAR_CATALOG}
        for alias, star_id in STAR_ALIASES.items():
            assert star_id in valid_star_ids, (
                f"Alias '{alias}' maps to invalid star ID {star_id}"
            )

    def test_star_aliases_keys_are_normalized(self):
        """Verify all alias keys are normalized (uppercase ASCII, as-is for Unicode)."""
        for alias in STAR_ALIASES:
            # Check that ASCII characters are uppercase
            ascii_chars = "".join(c for c in alias if c.isascii())
            assert ascii_chars == ascii_chars.upper(), (
                f"ASCII chars in alias '{alias}' should be uppercase"
            )


@pytest.mark.unit
class TestResolveStarName:
    """Tests for the resolve_star_name function."""

    def test_empty_name_returns_none(self):
        """Empty string should return None."""
        assert resolve_star_name("") is None
        assert resolve_star_name("   ") is None

    def test_exact_canonical_name_match(self):
        """Test exact match against canonical names."""
        assert resolve_star_name("Regulus") == REGULUS
        assert resolve_star_name("Spica") == SPICA_STAR
        assert resolve_star_name("Algol") == ALGOL
        assert resolve_star_name("Sirius") == SIRIUS

    def test_case_insensitive_matching(self):
        """Test that matching is case-insensitive."""
        assert resolve_star_name("REGULUS") == REGULUS
        assert resolve_star_name("regulus") == REGULUS
        assert resolve_star_name("ReGuLuS") == REGULUS
        assert resolve_star_name("SIRIUS") == SIRIUS
        assert resolve_star_name("sirius") == SIRIUS

    def test_whitespace_stripping(self):
        """Test that whitespace is stripped."""
        assert resolve_star_name("  Regulus  ") == REGULUS
        assert resolve_star_name("\tSirius\n") == SIRIUS

    def test_alias_matching(self):
        """Test matching against STAR_ALIASES."""
        # Regulus aliases
        assert resolve_star_name("Cor Leonis") == REGULUS
        assert resolve_star_name("Alpha Leo") == REGULUS
        assert resolve_star_name("Alpha Leonis") == REGULUS

        # Algol aliases
        assert resolve_star_name("Demon Star") == ALGOL
        assert resolve_star_name("Beta Persei") == ALGOL

        # Sirius aliases
        assert resolve_star_name("Dog Star") == SIRIUS
        assert resolve_star_name("Alpha CMa") == SIRIUS

    def test_comma_prefix_partial_match(self):
        """Test comma-prefix partial matching (reference ephemeris convention)."""
        # ,alg should find Algol
        result = resolve_star_name(",alg")
        # The expanded catalog resolves the prefix 'alg' to Algenib
        # (first match); the 108-star era found Algol.
        from libephemeris.fixed_stars import STAR_CATALOG

        names = {e.id: e.name for e in STAR_CATALOG}
        assert names[result].lower().startswith("alg"), (
            f"Expected an 'alg…' star, got {names.get(result, result)}"
        )

        # ,reg should find Regulus
        result = resolve_star_name(",reg")
        assert result == REGULUS, f"Expected Regulus, got {result}"

        # ,sir should find Sirius
        result = resolve_star_name(",sir")
        assert result == SIRIUS, f"Expected Sirius, got {result}"

        # ,spi should find Spica
        result = resolve_star_name(",spi")
        assert result == SPICA_STAR, f"Expected Spica, got {result}"

    def test_comma_prefix_with_whitespace(self):
        """Test comma-prefix with trailing whitespace."""
        result = resolve_star_name(", alg")
        from libephemeris.fixed_stars import STAR_CATALOG

        names = {e.id: e.name for e in STAR_CATALOG}
        assert names[result].lower().startswith("alg")
        assert resolve_star_name(",  reg  ") == REGULUS

    def test_comma_prefix_empty_returns_none(self):
        """Test that comma with no prefix returns None."""
        assert resolve_star_name(",") is None
        assert resolve_star_name(",  ") is None

    def test_trailing_comma_format(self):
        """Test reference ephemeris format with trailing comma (star_name,nomenclature)."""
        assert resolve_star_name("Regulus,alLeo") == REGULUS
        assert resolve_star_name("Spica,alVir") == SPICA_STAR

    def test_fuzzy_matching_for_short_inputs(self):
        """Test fuzzy matching (substring) for short inputs."""
        # "SIR" should match SIRIUS via substring
        result = resolve_star_name("SIR")
        assert result == SIRIUS, f"Expected Sirius for 'SIR', got {result}"

    def test_unknown_name_returns_none(self):
        """Test that unknown star names return None."""
        assert resolve_star_name("UnknownStar") is None
        assert resolve_star_name("NotARealStar") is None
        assert resolve_star_name("xyz123") is None


@pytest.mark.unit
class TestGetCanonicalStarName:
    """Tests for the get_canonical_star_name function."""

    def test_valid_star_ids_return_names(self):
        """Test that valid star IDs return canonical names."""
        assert get_canonical_star_name(REGULUS) == "Regulus"
        assert get_canonical_star_name(SPICA_STAR) == "Spica"
        assert get_canonical_star_name(ALGOL) == "Algol"
        assert get_canonical_star_name(SIRIUS) == "Sirius"

    def test_invalid_star_id_returns_none(self):
        """Test that invalid star IDs return None."""
        assert get_canonical_star_name(-1) is None
        assert get_canonical_star_name(999999) is None


@pytest.mark.unit
class TestSweFixstarUtWithAliases:
    """Tests for fixstar_ut with star name resolution."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch as Julian Day."""
        return 2451545.0

    def test_exact_name_works(self, standard_jd):
        """Test fixstar_ut with exact canonical names."""
        pos, name, retflag = fixstar_ut("Regulus", standard_jd, 0)
        assert name == "Regulus,alLeo"
        assert 149 < pos[0] < 151  # Regulus longitude

    def test_returns_canonical_name(self, standard_jd):
        """Test that fixstar_ut returns canonical star name."""
        pos, name, retflag = fixstar_ut("Alpha Leo", standard_jd, 0)
        assert name == "Regulus,alLeo", f"Expected 'Regulus,alLeo', got '{name}'"

        pos, name, retflag = fixstar_ut("Dog Star", standard_jd, 0)
        assert name == "Sirius,alCMa", f"Expected 'Sirius,alCMa', got '{name}'"

    def test_case_insensitive(self, standard_jd):
        """Test case-insensitive star name lookup."""
        pos1, name1, _ = fixstar_ut("SIRIUS", standard_jd, 0)
        pos2, name2, _ = fixstar_ut("sirius", standard_jd, 0)
        pos3, name3, _ = fixstar_ut("Sirius", standard_jd, 0)

        assert name1 == name2 == name3 == "Sirius,alCMa"
        assert abs(pos1[0] - pos2[0]) < 0.0001
        assert abs(pos1[0] - pos3[0]) < 0.0001

    def test_comma_prefix_search(self, standard_jd):
        """Partial nomenclature after comma errors (reference parity).

        Verified: swe.fixstar2_ut(',alg') raises 'could not find star'.
        """
        from libephemeris.exceptions import Error

        with pytest.raises(Error):
            fixstar_ut(",alg", standard_jd, 0)

    def test_bayer_designation(self, standard_jd):
        """Test Bayer designation lookup."""
        pos, name, _ = fixstar_ut("Alpha Leo", standard_jd, 0)
        assert name == "Regulus,alLeo"

        pos, name, _ = fixstar_ut("Alpha CMa", standard_jd, 0)
        assert name == "Sirius,alCMa"

    def test_identical_results_for_aliases(self, standard_jd):
        """Test that all aliases for a star return identical positions."""
        # Get Sirius position using different aliases
        pos_sirius, _, _ = fixstar_ut("Sirius", standard_jd, 0)
        pos_dog_star, _, _ = fixstar_ut("Dog Star", standard_jd, 0)
        pos_alpha_cma, _, _ = fixstar_ut("Alpha CMa", standard_jd, 0)

        # All should return the same position
        assert abs(pos_sirius[0] - pos_dog_star[0]) < 0.0001
        assert abs(pos_sirius[0] - pos_alpha_cma[0]) < 0.0001

    def test_unknown_star_returns_error(self, standard_jd):
        """Test that unknown star raises an error."""
        with pytest.raises(Exception):
            fixstar_ut("UnknownStar", standard_jd, 0)

    def test_backward_compatibility_regulus(self, standard_jd):
        """Test backward compatibility - Regulus still works."""
        pos, name, retflag = fixstar_ut("Regulus", standard_jd, 0)
        assert name == "Regulus,alLeo"
        assert 149 < pos[0] < 151
        assert -1 < pos[1] < 2

    def test_backward_compatibility_spica(self, standard_jd):
        """Test backward compatibility - Spica still works."""
        pos, name, retflag = fixstar_ut("Spica", standard_jd, 0)
        assert name == "Spica,alVir"
        assert 203 < pos[0] < 205
        assert -3 < pos[1] < -1


@pytest.mark.unit
class TestSweFixstarWithAliases:
    """Tests for fixstar (TT version) with star name resolution."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch as Julian Day."""
        return 2451545.0

    def test_exact_name_works(self, standard_jd):
        """Test fixstar with exact canonical names."""
        pos, name, retflag = fixstar("Regulus", standard_jd, 0)
        assert name == "Regulus,alLeo"

    def test_alias_works(self, standard_jd):
        """Test fixstar with aliases."""
        pos, name, retflag = fixstar("Alpha Leo", standard_jd, 0)
        assert name == "Regulus,alLeo"

    def test_comma_prefix_errors(self, standard_jd):
        """Partial nomenclature after comma errors (reference parity)."""
        from libephemeris.exceptions import Error

        with pytest.raises(Error):
            fixstar(",sir", standard_jd, 0)


@pytest.mark.unit
class TestAllStarsHaveAliases:
    """Tests to verify all catalog stars have proper aliases."""

    def test_all_catalog_stars_are_resolvable(self):
        """Test that all stars in catalog can be resolved by name."""
        for entry in STAR_CATALOG:
            result = resolve_star_name(entry.name)
            # Some names legitimately label several catalog entries
            # (e.g. the th01Ori Trapezium components): accept any
            # entry carrying the same name.
            same_name_ids = {e.id for e in STAR_CATALOG if e.name == entry.name}
            assert result in same_name_ids, (
                f"Star '{entry.name}' should resolve to one of "
                f"{sorted(same_name_ids)}, got {result}"
            )

    def test_bayer_designations_work(self):
        """Test various Bayer designation formats."""
        bayer_tests = [
            ("Alpha Leo", REGULUS),
            ("Alpha Leonis", REGULUS),
            ("Alpha Vir", SPICA_STAR),
            ("Alpha Virginis", SPICA_STAR),
            ("Beta Per", ALGOL),
            ("Beta Persei", ALGOL),
            ("Alpha Tau", ALDEBARAN),
            ("Alpha Sco", ANTARES),
            ("Alpha Lyr", VEGA),
            ("Alpha Ori", BETELGEUSE),
            ("Beta Ori", RIGEL),
            ("Alpha CMi", PROCYON),
            ("Alpha Boo", ARCTURUS),
            ("Alpha Cyg", DENEB),
            ("Beta Gem", POLLUX),
            ("Alpha Gem", CASTOR),
            ("Alpha Aql", ALTAIR),
        ]
        for designation, expected_id in bayer_tests:
            result = resolve_star_name(designation)
            assert result == expected_id, (
                f"Bayer designation '{designation}' should resolve to "
                f"{expected_id}, got {result}"
            )

    def test_common_names_work(self):
        """Test various common/traditional names."""
        common_name_tests = [
            ("Dog Star", SIRIUS),
            ("Demon Star", ALGOL),
            ("North Star", POLARIS),
            ("Pole Star", POLARIS),
        ]
        for name, expected_id in common_name_tests:
            result = resolve_star_name(name)
            assert result == expected_id, (
                f"Common name '{name}' should resolve to {expected_id}, got {result}"
            )


@pytest.mark.integration
class TestReferenceCompatibility:
    """Integration tests comparing against the reference ephemeris behavior."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch as Julian Day."""
        return 2451545.0

    def test_comma_prefix_errors_like_reference(self, standard_jd):
        """',alg' errors exactly like the reference ephemeris (exact-key semantics);
        the exact nomenclature ',bePer' finds Algol."""
        from libephemeris.exceptions import Error

        with pytest.raises(Error):
            fixstar_ut(",alg", standard_jd, 0)
        pos, name, retflag = fixstar_ut(",bePer", standard_jd, 0)
        assert name == "Algol,bePer"
        assert 55 < pos[0] < 57

    def test_alpha_leo_finds_regulus(self, standard_jd):
        """Test that Alpha Leo finds Regulus (reference ephemeris behavior)."""
        pos, name, retflag = fixstar_ut("Alpha Leo", standard_jd, 0)
        assert name == "Regulus,alLeo"

    def test_sirius_case_variations_identical(self, standard_jd):
        """Test SIRIUS and sirius return identical results."""
        pos_upper, name_upper, _ = fixstar_ut("SIRIUS", standard_jd, 0)
        pos_lower, name_lower, _ = fixstar_ut("sirius", standard_jd, 0)
        pos_mixed, name_mixed, _ = fixstar_ut("Sirius", standard_jd, 0)

        # All should return Sirius as canonical name
        assert name_upper == name_lower == name_mixed == "Sirius,alCMa"

        # All should return identical positions
        assert abs(pos_upper[0] - pos_lower[0]) < 0.0001
        assert abs(pos_upper[0] - pos_mixed[0]) < 0.0001

    def test_alpha_cma_finds_sirius(self, standard_jd):
        """Test that Alpha CMa finds Sirius (reference ephemeris behavior)."""
        pos, name, retflag = fixstar_ut("Alpha CMa", standard_jd, 0)
        assert name == "Sirius,alCMa"


@pytest.mark.unit
class TestReturnTypeStructure:
    """Tests for correct return type structure."""

    @pytest.fixture
    def standard_jd(self):
        """J2000.0 epoch as Julian Day."""
        return 2451545.0

    def test_swe_fixstar_ut_return_structure(self, standard_jd):
        """Test fixstar_ut returns correct tuple structure."""
        result = fixstar_ut("Regulus", standard_jd, 0)

        assert isinstance(result, tuple)
        assert len(result) == 3

        pos, name, retflag = result

        # Position tuple should have 6 elements
        assert isinstance(pos, tuple)
        assert len(pos) == 6

        # name should be a string (canonical star name)
        assert isinstance(name, str)
        assert name == "Regulus,alLeo"

        # retflag should be an int
        assert isinstance(retflag, int)

    def test_swe_fixstar_return_structure(self, standard_jd):
        """Test fixstar returns correct tuple structure."""
        result = fixstar("Spica", standard_jd, 0)

        assert isinstance(result, tuple)
        assert len(result) == 3

        pos, name, retflag = result

        assert isinstance(pos, tuple)
        assert len(pos) == 6
        assert isinstance(name, str)
        assert name == "Spica,alVir"
        assert isinstance(retflag, int)


@pytest.mark.unit
class TestCorCaroliProperName:
    """Regression: alpha-2 CVn carries its full IAU WGSN name 'Cor Caroli'.

    The generated catalog previously stored a truncated proper name 'Cor'
    (star id 1000781), so the bare 'Cor' resolved while the full 'Cor Caroli'
    errored. Measured reference behavior: the v2 search rejects the bare
    'Cor' and accepts 'Cor Caroli'; the returned name carries the full form.
    """

    JD = 2451545.0

    def test_catalog_entry_has_full_name(self):
        entry = next(e for e in STAR_CATALOG if e.id == 1000781)
        assert entry.name == "Cor Caroli"
        assert entry.nomenclature == "al02CVn"

    def test_cor_caroli_resolves_v2(self):
        from libephemeris.fixed_stars import _resolve_star2

        entry, err = _resolve_star2("Cor Caroli")
        assert err is None and entry is not None
        assert entry.id == 1000781

    def test_bare_cor_errors_v2(self):
        from libephemeris.fixed_stars import _resolve_star2

        entry, err = _resolve_star2("Cor")
        assert entry is None and err is not None

    def test_fixstar2_ut_returns_full_name(self):
        import libephemeris as L
        from libephemeris.fixed_stars import fixstar2_ut

        L.set_calc_mode("skyfield")
        try:
            _pos, name, _retflag = fixstar2_ut("Cor Caroli", self.JD, 0)
        finally:
            L.set_calc_mode("auto")
        assert name.startswith("Cor Caroli")

    def test_other_multiword_names_intact(self):
        """The other 22 multi-word proper names must stay intact."""
        multiword = {e.name for e in STAR_CATALOG if " " in e.name}
        assert "Cor Caroli" in multiword
        for expected in (
            "Kaus Australis",
            "Yed Prior",
            "Rigil Kentaurus",
            "Polaris Australis",
            "Asellus Borealis",
        ):
            assert expected in multiword
