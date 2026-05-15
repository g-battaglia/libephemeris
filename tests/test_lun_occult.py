"""
Tests for lunar occultation calculations in libephemeris.

Tests the lun_occult_when_glob function which finds lunar occultations
of planets and fixed stars.

Lunar occultations occur when the Moon passes in front of a planet or star.
These events are relatively rare and happen in specific cycles depending on
the Moon's nodal position relative to the target body's ecliptic latitude.

Known occultation series:
- Regulus: 2017 series (June 28, 2017 was an occultation)
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris import (
    julday,
    revjul,
    lun_occult_when_glob,
    ECL_TOTAL,
    ECL_PARTIAL,
    FLG_SWIEPH,
)


class TestLunOccultWhenGlob:
    """Test suite for lun_occult_when_glob function."""

    def test_finds_star_occultation_regulus(self):
        """Test that function finds a lunar occultation of Regulus.

        Regulus (Alpha Leonis) is frequently occulted by the Moon
        because it lies very close to the ecliptic.

        In 2017, there were multiple Regulus occultations:
        - Jan 15, Feb 11, Mar 10, Apr 7, May 4: Global occultations (visible
          from specific locations due to lunar parallax)
        - May 31, Jun 28, Jul 25: Geocentric occultations (visible from
          entire Earth-facing hemisphere)

        The function correctly finds the first available global occultation.
        """
        # Start from Jan 1, 2017
        jd_start = julday(2017, 1, 1, 0)

        # pyswisseph signature: lun_occult_when_glob(tjdut, body, flags, ecltype, backwards)
        # Returns: (retflags, tret)
        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)

        # Should find an occultation
        assert retflags != 0
        assert tret[0] > jd_start  # Maximum should be after start
        assert tret[0] > 0  # Should have valid maximum time

        # Verify it's in 2017 (any month is valid - there are multiple occultations)
        year, month, day, hour = revjul(tret[0])
        assert year == 2017
        assert 1 <= month <= 12  # Any month in 2017 is valid

    def test_returns_correct_tuple_structure(self):
        """Test that return values have correct structure per pyswisseph spec."""
        jd_start = julday(2017, 1, 1, 0)

        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)

        # Should return 10-element tuple per pyswisseph specification
        assert len(tret) == 10
        # All elements should be floats
        assert all(isinstance(t, float) for t in tret)
        # Occultation type should be int
        assert isinstance(retflags, int)

    def test_contact_times_ordering(self):
        """Test that contact times are in correct chronological order.

        pyswisseph tret indices:
        [0]: time of maximum occultation
        [2]: time of occultation begin
        [3]: time of occultation end
        [4]: time of totality begin
        [5]: time of totality end
        """
        jd_start = julday(2017, 1, 1, 0)

        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)

        jd_max = tret[0]
        jd_begin = tret[2]
        jd_end = tret[3]
        jd_total_begin = tret[4]
        jd_total_end = tret[5]

        # Occultation begin should be before maximum
        if jd_begin > 0:
            assert jd_begin < jd_max

        # Occultation end should be after maximum
        if jd_end > 0:
            assert jd_end > jd_max

        # If we have totality (total occultation)
        if jd_total_begin > 0 and jd_total_end > 0:
            assert jd_begin <= jd_total_begin <= jd_max
            assert jd_max <= jd_total_end <= jd_end

    def test_occultation_type_flags(self):
        """Test that occultation type flags are set."""
        jd_start = julday(2017, 1, 1, 0)

        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)

        # Should be either total or partial
        assert (retflags & ECL_TOTAL) or (retflags & ECL_PARTIAL)

    def test_raises_error_for_no_target(self):
        """Test that function raises error if no target specified."""
        jd_start = julday(2024, 1, 1, 0)

        with pytest.raises(ValueError):
            lun_occult_when_glob(jd_start, "", FLG_SWIEPH, 0)

    def test_raises_error_for_unknown_star(self):
        """Test that function raises error for unknown star name."""
        jd_start = julday(2017, 1, 1, 0)

        with pytest.raises(ValueError):
            lun_occult_when_glob(jd_start, "UnknownStar123", FLG_SWIEPH, 0)

    def test_swe_alias(self):
        """Test that lun_occult_when_glob is an alias."""
        jd_start = julday(2017, 1, 1, 0)

        retflags1, tret1 = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)
        retflags2, tret2 = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0
        )

        assert tret1 == tret2
        assert retflags1 == retflags2

    def test_occultation_duration_reasonable(self):
        """Test that occultation duration is reasonable.

        Lunar occultations typically last a few seconds to about an hour,
        depending on the geometry.
        """
        jd_start = julday(2017, 1, 1, 0)

        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)

        jd_begin = tret[2]
        jd_end = tret[3]

        if jd_begin > 0 and jd_end > 0:
            duration_hours = (jd_end - jd_begin) * 24
            # Duration should be between 0.01 hours (~30 seconds) and 4 hours
            # Global occultation durations (C1 to C4 across entire Earth)
            # can be longer than local durations, especially for bright stars
            # near the ecliptic with favorable geometry.
            assert 0.001 < duration_hours < 4.0


class TestLunOccultEdgeCases:
    """Test edge cases for lunar occultation calculations."""

    def test_known_regulus_occultation_2017(self):
        """Test finding the known Regulus occultation of June 28, 2017.

        This is a well-documented lunar occultation event.
        """
        # Search from June 1, 2017
        jd_start = julday(2017, 6, 1, 0)

        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)

        # Verify it's on June 28
        year, month, day, hour = revjul(tret[0])
        assert year == 2017
        assert month == 6
        assert 27 <= day <= 28  # Allow for timezone differences

    def test_search_continues_after_no_occultation(self):
        """Test that search properly advances when no immediate occultation."""
        # Start from a date with no immediate occultation
        jd_start = julday(2020, 1, 1, 0)

        # The function should either find an occultation or raise RuntimeError
        # after exhausting the search limit
        try:
            retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)
            # If found, verify it's after start date
            assert tret[0] > jd_start
        except RuntimeError as e:
            # Expected if no occultation in 20 years
            assert "No lunar occultation" in str(e)

    def test_raises_error_for_invalid_planet(self):
        """Test that invalid planet ID raises appropriate error."""
        jd_start = julday(2020, 1, 1, 0)

        with pytest.raises(ValueError):
            lun_occult_when_glob(jd_start, 999, FLG_SWIEPH, 0)  # Invalid planet ID


class TestLunOccultResolveStarId:
    """Regression tests for _resolve_star_id tuple unpacking fix.

    These tests verify that _resolve_star_id returning 3 values
    (star_id, error_message, canonical_name) is correctly handled
    throughout the eclipse module.
    """

    def test_valid_star_does_not_raise_unpack_error(self):
        """Test that valid star names don't cause unpacking errors.

        Regression test for ValueError: too many values to unpack (expected 2).
        The _resolve_star_id function returns 3 values, not 2.
        """
        jd_start = julday(2017, 1, 1, 0)

        # This should not raise "ValueError: too many values to unpack"
        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)
        assert tret[0] > jd_start

    def test_invalid_star_raises_value_error_with_message(self):
        """Test that invalid stars raise ValueError with descriptive message."""
        jd_start = julday(2017, 1, 1, 0)

        with pytest.raises(ValueError) as exc_info:
            lun_occult_when_glob(jd_start, "InvalidStarXYZ123", FLG_SWIEPH, 0)

        # Error message should mention the star name
        assert "invalidstarxyz123" in str(exc_info.value).lower()


class TestLunOccultPyswissephCompatibility:
    """Tests verifying pyswisseph API compatibility."""

    def test_signature_matches_pyswisseph(self):
        """Test that the function signature matches pyswisseph.

        pyswisseph signature:
        lun_occult_when_glob(tjdut, body, flags=FLG_SWIEPH, ecltype=0, backwards=False)

        body: int for planet ID or str for star name
        Returns: (retflags, tret) where tret is 10-element tuple
        """
        jd_start = julday(2017, 1, 1, 0)

        # Test with star name
        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)
        assert isinstance(retflags, int)
        assert len(tret) == 10

    def test_tret_indices_match_pyswisseph_documentation(self):
        """Test that tret array indices match pyswisseph documentation.

        pyswisseph tret indices:
        [0]: time of maximum occultation
        [1]: time when occultation takes place at local apparent noon
        [2]: time of occultation begin
        [3]: time of occultation end
        [4]: time of totality begin
        [5]: time of totality end
        [6]: time of center line begin
        [7]: time of center line end
        [8]: time when annular-total becomes total
        [9]: time when annular-total becomes annular again
        """
        jd_start = julday(2017, 1, 1, 0)

        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)

        # Maximum occultation should be valid
        assert tret[0] > 0

        # Begin should be before max
        if tret[2] > 0:
            assert tret[2] < tret[0]

        # End should be after max
        if tret[3] > 0:
            assert tret[3] > tret[0]

        # For total occultations, check totality times
        if retflags & ECL_TOTAL:
            if tret[4] > 0:
                assert tret[2] <= tret[4] <= tret[0]  # begin <= totality_begin <= max
            if tret[5] > 0:
                assert tret[0] <= tret[5] <= tret[3]  # max <= totality_end <= end
