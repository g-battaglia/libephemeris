"""
Unit tests: HIP-number search strings do NOT resolve in the v2 star search.

Measured reference behavior: the v2 family (fixstar2/fixstar2_ut/
fixstar2_mag) resolves traditional names, comma-prefixed nomenclature,
sequential catalog numbers and trailing-'%' wildcards — never Hipparcos
numbers. "HIP 49669" fails as an unknown name, ",49669" fails as an
unknown nomenclature, and a bare "49669" reads as an (out-of-range)
sequential index.

HIP numbers remain available through the library-specific helpers
(get_hip_from_star_name / STAR_NAME_TO_HIP), not through the
reference-named fixstar functions.
"""

import pytest
from libephemeris.exceptions import Error
from libephemeris.fixed_stars import _resolve_star2, fixstar2_ut


@pytest.mark.unit
class TestHipPrefixRejected:
    """ "HIP NNNNN" strings are unknown names for _resolve_star2."""

    @pytest.mark.parametrize(
        "query",
        [
            "HIP 49669",
            "HIP 65474",
            "HIP65474",
            "hip 49669",
            "Hip 49669",
            "  HIP 49669  ",
            "HIP   49669",
            "HIP 999999999",
            "HIP abc",
        ],
    )
    def test_hip_prefix_never_resolves(self, query):
        """Any HIP-prefixed string fails as an unknown traditional name."""
        entry, err = _resolve_star2(query)
        assert entry is None
        assert err is not None
        assert "could not find star name" in err

    def test_hip_prefix_error_echoes_normalized_name(self):
        """The error echoes the lowercased, space-stripped search key."""
        entry, err = _resolve_star2("HIP 49669")
        assert entry is None
        assert err == "could not find star name hip49669"


@pytest.mark.unit
class TestHipFormsViaFixstar2:
    """The public fixstar2_ut rejects every HIP-style form."""

    JD = 2451545.0  # J2000.0

    def test_fixstar2_ut_hip_prefix_raises(self):
        with pytest.raises(Error, match="could not find star name hip49669"):
            fixstar2_ut("HIP 49669", self.JD, 0)

    def test_fixstar2_ut_comma_number_raises(self):
        """A comma form keys on the nomenclature, never on a HIP number."""
        with pytest.raises(Error, match="could not find star name ,65474"):
            fixstar2_ut(",65474", self.JD, 0)

    def test_fixstar2_ut_bare_number_is_sequential(self):
        """A bare number is a sequential index (here: out of range)."""
        with pytest.raises(Error, match="sequential fixed star number"):
            fixstar2_ut("49669", self.JD, 0)
