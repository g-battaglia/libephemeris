"""Guards for the fixed-star name resolvers.

The recorded golden family exercises the two resolvers hard, but three of
their rules leave no trace in it: an implementation that deleted the alias
table, the curated name corrections and the v1 implicit prefix search would
still replay the whole family without a single difference.  These tests are
the gate those rules do not otherwise have.

They also pin, on the same string, the five places where the two families
answer differently.  Reading them as one function is the mistake this file
exists to catch:

* the comma form is a mirror — v1 keys on the half before the comma, v2 on
  the half after it;
* a trailing comma resolves in v1 and refuses in v2;
* the wildcard belongs to v2 alone;
* a bare number indexes two different orders of the catalogue;
* the last resort differs — v1 searches names by prefix, v2 consults the
  curated corrections — so three references reach a star in both families
  and reach a *different* star.

Named answers are quoted as "Name,Designation", which is what both public
families report.
"""

import pytest

from libephemeris.exceptions import StarNotFoundError
from libephemeris.fixed_stars import (
    _NAME_HIP_FIX,
    _format_star_name,
    _resolve_star2,
    _resolve_star_id,
    _resolve_star_ref,
    fixstar2_ut,
    fixstar_ut,
)

JD = 2451545.0


def v1(reference: str) -> str | None:
    """The star the v1 resolver selects, as "Name,Designation", or None."""
    star_id, error, name = _resolve_star_ref(reference)
    if error is not None:
        assert star_id == -1 and name is None
        return None
    assert star_id != -1
    return name


def v2(reference: str) -> str | None:
    """The star the v2 resolver selects, as "Name,Designation", or None."""
    entry, error = _resolve_star2(reference)
    if error is not None:
        assert entry is None
        return None
    assert entry is not None
    return _format_star_name(entry)


@pytest.mark.unit
class TestAliasTableIsReachable:
    """The alias table resolves names no catalogue name carries.

    Unwitnessed by the golden family: every one of these would refuse if the
    table were dropped, and no recorded case would notice.
    """

    @pytest.mark.parametrize(
        "reference,expected",
        [
            ("Dog Star", "Sirius,alCMa"),
            ("dogstar", "Sirius,alCMa"),
            ("DOG STAR", "Sirius,alCMa"),
            ("Cor Leonis", "Regulus,alLeo"),
            ("Betelgeux", "Betelgeuse,alOri"),
        ],
    )
    def test_alias_resolves_in_both_families(self, reference, expected):
        """An alias is a name, under the same fold, in both families."""
        assert v1(reference) == expected
        assert v2(reference) == expected

    @pytest.mark.parametrize("reference", ["α LEO", "Α LEO", "α leo"])
    def test_greek_alias_survives_the_fold(self, reference):
        """Lowercasing is the Unicode one, so the Greek keys stay reachable."""
        assert v1(reference) == "Regulus,alLeo"
        assert v2(reference) == "Regulus,alLeo"

    @pytest.mark.parametrize("reference", ["alLeo", "ALLEO", "alleo"])
    def test_designation_shaped_alias_is_not_a_bare_name(self, reference):
        """A bare designation must never resolve — only the comma form does.

        The alias table carries keys that merely mirror a catalogue
        designation; they stay unreachable as bare names, or the rule above
        would be defeated.
        """
        assert v1(reference) is None
        assert v2(reference) is None
        assert v1(",alLeo") == "Regulus,alLeo"
        assert v2(",alLeo") == "Regulus,alLeo"

    @pytest.mark.parametrize(
        "reference,expected",
        [("Suhail", "Suhail,laVel"), ("Alnair", "Alnair,alGru")],
    )
    def test_catalogue_name_outranks_a_colliding_alias(self, reference, expected):
        """Both collisions select the catalogue's own entry, not the alias."""
        assert v1(reference) == expected
        assert v2(reference) == expected

    def test_digit_leading_alias_is_read_as_a_number_first(self):
        """An alias whose key begins with a digit is read as a number first."""
        assert v1("32 Leo") == "134Tau,134Tau"
        assert v2("32 Leo") == "134Tau,134Tau"


@pytest.mark.unit
class TestCuratedCorrectionsAreReachable:
    """The five curated corrections are v2's last resort.

    Unwitnessed by the golden family. Each key is a traditional name the
    catalogue does not carry, so the table never competes with it.
    """

    @pytest.mark.parametrize(
        "reference,expected",
        [
            ("Alaraph", "Zavijava,beVir"),
            ("Gienah Corvi", "Gienah,gaCrv"),
            ("Atri", "Megrez,deUMa"),
            ("Nash", "Alnasl,gaSgr"),
            ("Deli", "etAqr,etAqr"),
        ],
    )
    def test_correction_resolves_in_v2(self, reference, expected):
        assert v2(reference) == expected

    def test_every_correction_key_is_reachable(self):
        """No key of the table is shadowed: all five resolve in v2."""
        for key in _NAME_HIP_FIX:
            assert v2(key) is not None, key

    def test_v1_does_not_consult_the_corrections(self):
        """v1 ends with a prefix search instead, and no name begins so."""
        assert v1("Alaraph") is None
        assert v1("Gienah Corvi") is None


@pytest.mark.unit
class TestImplicitPrefixSearchIsReachable:
    """v1's last resort: the first catalogue row whose name starts with the key.

    Unwitnessed by the golden family, and the only rule that answers these
    references at all.
    """

    def test_partial_name_resolves_in_v1_only(self):
        assert v1("aldeb") == "Aldebaran,alTau"
        assert v2("aldeb") is None

    @pytest.mark.parametrize(
        "reference,expected",
        [
            ("t", "thOct,thOct"),
            ("z", "zeScl,zeScl"),
            ("a", "Alpheratz,alAnd"),
            ("c", "Caph,beCas"),
        ],
    )
    def test_the_scan_follows_the_catalogue_row_order(self, reference, expected):
        """A one-letter key discriminates the row order from the name order."""
        assert v1(reference) == expected

    def test_the_prefix_search_overrides_the_correction_table(self):
        """Documented, not fixed: "Atri" is a prefix of "Atria" in v1.

        The correction table says "Atri" means Megrez, and v2 answers that;
        v1 never reaches the table, so its prefix search answers Atria.
        """
        assert v1("Atri") == "Atria,alTrA"
        assert v2("Atri") == "Megrez,deUMa"


@pytest.mark.unit
class TestBothFamiliesDifferentStars:
    """The three references that resolve in both families to different stars.

    They are the sharpest witness that the two resolvers are not one
    function: each reaches a star through v1's prefix search and another
    through v2's correction table.
    """

    @pytest.mark.parametrize(
        "reference,by_v1,by_v2",
        [
            ("Atri", "Atria,alTrA", "Megrez,deUMa"),
            ("Nash", "Nashira,gaCap", "Alnasl,gaSgr"),
            ("Deli", "deLib,deLib", "etAqr,etAqr"),
        ],
    )
    def test_same_string_two_stars(self, reference, by_v1, by_v2):
        assert v1(reference) == by_v1
        assert v2(reference) == by_v2
        assert by_v1 != by_v2


@pytest.mark.unit
class TestTheFiveDifferences:
    """One test per row of the measured difference table."""

    @pytest.mark.parametrize(
        "reference,by_v1,by_v2",
        [
            ("Regulus,alVir", "Regulus,alLeo", "Spica,alVir"),
            ("Spica,alLeo", "Spica,alVir", "Regulus,alLeo"),
        ],
    )
    def test_the_comma_form_is_a_mirror(self, reference, by_v1, by_v2):
        """v1 keys on the half before the comma, v2 on the half after it."""
        assert v1(reference) == by_v1
        assert v2(reference) == by_v2

    def test_only_the_first_comma_splits(self):
        assert v1("Regulus,alLeo,alVir") == "Regulus,alLeo"
        assert v2("Regulus,alLeo,alVir") is None

    def test_a_trailing_comma_resolves_in_v1_and_refuses_in_v2(self):
        assert v1("Regulus,") == "Regulus,alLeo"
        assert v2("Regulus,") is None

    def test_the_wildcard_belongs_to_v2(self):
        assert v1("Regu%") is None
        assert v2("Regu%") == "Regulus,alLeo"

    @pytest.mark.parametrize(
        "reference,by_v1,by_v2",
        [
            ("1", "101Her,101Her", "101Her,101Her"),
            ("1447", "Zubeneschamali,beLib", "zeVol,zeVol"),
        ],
    )
    def test_a_bare_number_indexes_two_orders(self, reference, by_v1, by_v2):
        """v1 indexes the name-key order, v2 the designation order."""
        assert v1(reference) == by_v1
        assert v2(reference) == by_v2

    def test_the_last_resort_differs(self):
        """v1 ends with the prefix search, v2 with the correction table."""
        assert v1("aldeb") == "Aldebaran,alTau" and v2("aldeb") is None
        assert v1("Alaraph") is None and v2("Alaraph") == "Zavijava,beVir"


@pytest.mark.unit
class TestFormPrecedence:
    """Which form a reference has, where the two families order it differently."""

    def test_a_percent_sign_is_an_ordinary_character_in_v1(self):
        """The key "12%" is a sequential index in v1 and a wildcard in v2."""
        assert v1("12%") == "10UMa,10UMa"
        assert v2("12%") == "125Tau,125Tau"

    def test_the_wildcard_outranks_the_comma_in_v2(self):
        assert v1("Regulus,alLeo%") == "Regulus,alLeo"
        assert v2("Regulus,alLeo%") is None

    def test_a_digit_leading_name_half_is_still_a_number_in_v1(self):
        assert v1("12,alLeo") == "10UMa,10UMa"
        assert v2("12,alLeo") == "Regulus,alLeo"


@pytest.mark.unit
class TestWildcardEdges:
    """The v2 wildcard, including the case the family never records."""

    def test_a_prefix_no_name_begins_with_refuses(self):
        entry, error = _resolve_star2("zzzz%")
        assert entry is None
        assert error is not None

    def test_the_bare_wildcard_selects_the_first_name(self):
        assert v2("%") == "101Her,101Her"

    @pytest.mark.parametrize("reference", ["Al%de", "Alde%%", "%Alde"])
    def test_a_percent_sign_elsewhere_refuses(self, reference):
        assert v2(reference) is None


@pytest.mark.unit
class TestSearchTypeClassification:
    """The word the raised error carries, for each of its five values.

    The harness records the type and the message only, so this attribute is
    contract the golden cannot police.
    """

    @pytest.mark.parametrize(
        "reference,expected",
        [
            ("", "empty"),
            ("   ", "empty"),
            (",ZZZZZ", "nomenclature"),
            ("999999", "sequential"),
            ("Nosuchstarname", "name"),
            ("Al%de", "name"),
        ],
    )
    def test_v1_classification(self, reference, expected):
        with pytest.raises(StarNotFoundError) as excinfo:
            fixstar_ut(reference, JD, 0)
        assert excinfo.value.search_type == expected

    @pytest.mark.parametrize(
        "reference,expected",
        [
            ("", "empty"),
            (",ZZZZZ", "nomenclature"),
            ("999999", "sequential"),
            ("Nosuchstarname", "name"),
            ("Al%de", "wildcard"),
        ],
    )
    def test_v2_classification(self, reference, expected):
        with pytest.raises(StarNotFoundError) as excinfo:
            fixstar2_ut(reference, JD, 0)
        assert excinfo.value.search_type == expected

    def test_the_same_string_is_classified_differently(self):
        """The key "Al%de" is a name to v1 and a wildcard to v2."""
        with pytest.raises(StarNotFoundError) as v1_error:
            fixstar_ut("Al%de", JD, 0)
        with pytest.raises(StarNotFoundError) as v2_error:
            fixstar2_ut("Al%de", JD, 0)
        assert v1_error.value.search_type == "name"
        assert v2_error.value.search_type == "wildcard"


@pytest.mark.unit
class TestResolverIdentityAndPurity:
    """The two names of the v1 resolver, and its independence of everything."""

    def test_resolve_star_id_is_the_same_answer(self):
        assert _resolve_star_id("Regulus") == _resolve_star_ref("Regulus")
        assert _resolve_star_id("Nosuchstar") == _resolve_star_ref("Nosuchstar")

    def test_neither_resolver_raises_on_a_refusal(self):
        star_id, error, name = _resolve_star_ref("Nosuchstarname")
        assert (star_id, name) == (-1, None) and error
        entry, error = _resolve_star2("Nosuchstarname")
        assert entry is None and error
