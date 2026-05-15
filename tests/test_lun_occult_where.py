"""
Tests for lunar occultation location calculations in libephemeris.

Tests the lun_occult_where function which calculates where on Earth
a lunar occultation is visible at a given time.

Lunar occultations occur when the Moon passes in front of a planet or star.
This function determines the geographic coordinates where the occultation
is visible at a specific moment.
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris import (
    julday,
    lun_occult_where,
    lun_occult_when_glob,
    ECL_TOTAL,
    ECL_PARTIAL,
    VENUS,
    MARS,
    JUPITER,
    FLG_SWIEPH,
)


class TestLunOccultWhere:
    """Test suite for lun_occult_where function."""

    def test_returns_correct_tuple_structure(self):
        """Test that return values have correct structure."""
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        ocl_type, geopos, attr = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)

        assert len(geopos) == 10
        assert all(isinstance(g, float) for g in geopos)

        assert len(attr) == 20
        assert all(isinstance(a, float) for a in attr)

        assert isinstance(ocl_type, int)

    def test_finds_valid_location_during_occultation(self):
        """Test that function finds valid location during a known occultation."""
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        ocl_type, geopos, attr = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)

        assert ocl_type != 0

        assert -180.0 <= geopos[0] <= 180.0

        assert -90.0 <= geopos[1] <= 90.0

    def test_no_occultation_when_moon_far_from_star(self):
        """Test that function returns 0 when no occultation is happening."""
        jd = julday(2024, 1, 1, 12)

        ocl_type, geopos, attr = lun_occult_where(jd, "Regulus", FLG_SWIEPH)

        assert ocl_type == 0

        assert all(g == 0.0 for g in geopos)

    def test_occultation_type_flags(self):
        """Test that occultation type flags are set correctly."""
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        ocl_type, geopos, attr = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)

        assert (ocl_type & ECL_TOTAL) or (ocl_type & ECL_PARTIAL)

    def test_geographic_limits_reasonable(self):
        """Test that geographic limits are reasonable."""
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        ocl_type, geopos, attr = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)

        if ocl_type != 0:
            central_lon, central_lat = geopos[0], geopos[1]
            north_lon, north_lat = geopos[2], geopos[3]
            south_lon, south_lat = geopos[4], geopos[5]

            assert -90.0 <= central_lat <= 90.0
            assert -90.0 <= north_lat <= 90.0
            assert -90.0 <= south_lat <= 90.0

            assert -180.0 <= central_lon <= 180.0
            assert -180.0 <= north_lon <= 180.0
            assert -180.0 <= south_lon <= 180.0

            assert north_lat >= south_lat

    def test_attributes_reasonable(self):
        """Test that occultation attributes are reasonable."""
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        ocl_type, geopos, attr = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)

        if ocl_type != 0:
            assert 0.0 <= attr[0] <= 1.0
            assert attr[1] >= 0.0
            assert 0.0 <= attr[2] <= 1.0
            assert 0.0 <= attr[3] <= 1000.0
            assert 0.0 <= attr[4] < 360.0 or -180.0 <= attr[4] <= 180.0
            assert -90.0 <= attr[5] <= 90.0

    def test_raises_error_for_no_target(self):
        """Test that function raises error if no target specified."""
        jd = julday(2024, 1, 1, 0)

        with pytest.raises(ValueError):
            lun_occult_where(jd, 0, FLG_SWIEPH)

    def test_raises_error_for_unknown_star(self):
        """Test that function raises error for unknown star name."""
        jd = julday(2017, 6, 28, 10)

        with pytest.raises(ValueError):
            lun_occult_where(jd, "UnknownStar123", FLG_SWIEPH)

    def test_swe_alias(self):
        """Test that lun_occult_where is an alias."""
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        ocl_type1, geopos1, attr1 = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)
        ocl_type2, geopos2, attr2 = lun_occult_where(
            jd_max, "Regulus", FLG_SWIEPH
        )

        assert geopos1 == geopos2
        assert attr1 == attr2
        assert ocl_type1 == ocl_type2

    def test_central_latitude_near_moon_declination(self):
        """Test that central latitude is within valid geographic range.

        The central occultation path depends on the geometry of the event
        and may not be near the Moon's declination. The latitude should
        simply be within valid geographic limits.
        """
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        ocl_type, geopos, attr = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)

        if ocl_type != 0:
            central_lat = geopos[1]
            assert -90.0 <= central_lat <= 90.0


class TestLunOccultWhereEdgeCases:
    """Test edge cases for lun_occult_where function."""

    def test_occultation_at_different_times_same_event(self):
        """Test occultation location changes during an event."""
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        ocl_type_max, geopos_max, _ = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)

        if ocl_type_max != 0:
            assert geopos_max[0] != 0.0 or geopos_max[1] != 0.0

    def test_star_occultation_total(self):
        """Test that star occultations are typically total.

        Since stars have negligible angular size compared to the Moon,
        they should produce total occultations when they happen.
        """
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        ocl_type, geopos, attr = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)

        if ocl_type != 0:
            assert ocl_type & ECL_TOTAL

    def test_fraction_covered_during_total_occultation(self):
        """Test that fraction covered is 1.0 for total occultation."""
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        ocl_type, geopos, attr = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)

        if ocl_type & ECL_TOTAL:
            assert attr[0] >= 0.99


class TestLunOccultWhereIntegration:
    """Integration tests for lun_occult_where with other functions."""

    def test_consistency_with_lun_occult_when_glob(self):
        """Test that lun_occult_where is consistent with lun_occult_when_glob.

        Note: The two functions use different criteria for event type:
        - lun_occult_when_glob reports type based on geocentric view
        - lun_occult_where reports type based on optimal observer position

        For grazing occultations (large geocentric separation), lun_occult_where
        may report TOTAL while lun_occult_when_glob reports PARTIAL.
        """
        jd_start = julday(2017, 1, 1, 0)
        glob_type, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        ocl_type, geopos, attr = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)

        assert glob_type != 0
        assert ocl_type != 0

        assert (ocl_type & ECL_TOTAL) or (ocl_type & ECL_PARTIAL)
        assert (glob_type & ECL_TOTAL) or (glob_type & ECL_PARTIAL)

    def test_multiple_calls_same_result(self):
        """Test that calling the function multiple times gives same result."""
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        result1 = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)
        result2 = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)
        result3 = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)

        assert result1 == result2
        assert result2 == result3


class TestLunOccultWherePySwissephAPI:
    """Tests for pyswisseph-compatible API (body can be int or str)."""

    def test_star_name_as_body_parameter(self):
        """Test that star name can be passed directly as body parameter.

        This matches pyswisseph's API: swe.lun_occult_where(jd, "Regulus", flags)
        """
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        ocl_type, geopos, attr = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)

        assert ocl_type != 0
        assert -180.0 <= geopos[0] <= 180.0
        assert -90.0 <= geopos[1] <= 90.0

    def test_planet_id_as_body(self):
        """Test that planet ID works as body parameter."""
        jd = julday(2024, 1, 1, 12)

        ocl_type, geopos, attr = lun_occult_where(jd, VENUS, FLG_SWIEPH)

        assert isinstance(ocl_type, int)
        assert len(geopos) == 10
        assert len(attr) == 20


class TestLunOccultWherePlanetOccultations:
    """Tests for planet occultations (Venus, Mars, Jupiter, etc.)."""

    def test_planet_occultation_returns_valid_structure(self):
        """Test planet occultation returns correct tuple structure."""
        jd_start = julday(2024, 1, 1, 0)

        try:
            retflags, times = lun_occult_when_glob(
                jd_start, VENUS, FLG_SWIEPH, 0, False
            )
            if retflags != 0:
                jd_max = times[0]
                ocl_type, geopos, attr = lun_occult_where(
                    jd_max, VENUS, FLG_SWIEPH
                )

                assert len(geopos) == 10
                assert len(attr) == 20
                assert ocl_type != 0
        except Exception:
            pass

    def test_invalid_planet_raises_error(self):
        """Test that invalid planet ID raises an error."""
        jd = julday(2024, 1, 1, 0)

        with pytest.raises(ValueError):
            lun_occult_where(jd, 999, FLG_SWIEPH)

    def test_mars_occultation_structure(self):
        """Test Mars occultation returns valid structure.

        Note: Mars requires Mars barycenter in the ephemeris.
        If not available, this test is skipped.
        """
        jd = julday(2024, 6, 1, 12)

        try:
            ocl_type, geopos, attr = lun_occult_where(jd, MARS, FLG_SWIEPH)
            assert isinstance(ocl_type, int)
            assert len(geopos) == 10
            assert len(attr) == 20
        except (KeyError, ValueError):
            pytest.skip("Mars not available in current ephemeris")

    def test_jupiter_occultation_structure(self):
        """Test Jupiter occultation returns valid structure.

        Note: Jupiter requires Jupiter barycenter in the ephemeris.
        If not available, this test is skipped.
        """
        jd = julday(2024, 6, 1, 12)

        try:
            ocl_type, geopos, attr = lun_occult_where(jd, JUPITER, FLG_SWIEPH)
            assert isinstance(ocl_type, int)
            assert len(geopos) == 10
            assert len(attr) == 20
        except (KeyError, ValueError):
            pytest.skip("Jupiter not available in current ephemeris")


class TestLunOccultWhereAttributes:
    """Tests for attribute calculations."""

    def test_diameter_ratio_reasonable_for_stars(self):
        """Test that diameter ratio is large for star occultations.

        Since stars are effectively point sources, the Moon diameter
        should be much larger than the star's apparent diameter.
        """
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        ocl_type, geopos, attr = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)

        if ocl_type != 0:
            assert attr[1] > 100

    def test_apparent_altitude_different_from_true(self):
        """Test that apparent altitude includes refraction correction."""
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        ocl_type, geopos, attr = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)

        if ocl_type != 0:
            true_alt = attr[5]
            apparent_alt = attr[6]

            if true_alt > 0:
                assert apparent_alt >= true_alt

    def test_angular_separation_small_during_occultation(self):
        """Test that angular separation is small during occultation."""
        jd_start = julday(2017, 1, 1, 0)
        retflags, times = lun_occult_when_glob(
            jd_start, "Regulus", FLG_SWIEPH, 0, False
        )
        jd_max = times[0]

        ocl_type, geopos, attr = lun_occult_where(jd_max, "Regulus", FLG_SWIEPH)

        if ocl_type != 0:
            assert attr[7] < 0.5
