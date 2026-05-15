"""
Tests for house_pos and house_pos functions.

Tests the calculation of house position for celestial bodies.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


class TestHousePosBasic:
    """Basic tests for house_pos function."""

    @pytest.mark.unit
    def test_house_pos_returns_float(self):
        """house_pos() should return a float."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393
        lon = 15.0

        result = ephem.house_pos(armc, lat, eps, ord("P"), lon, 0.0)

        assert isinstance(result, float)

    @pytest.mark.unit
    def test_house_pos_in_valid_range(self):
        """house_pos() should return value between 1.0 and 13.0."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393
        lon = 15.0

        result = ephem.house_pos(armc, lat, eps, ord("P"), lon, 0.0)

        assert 1.0 <= result < 13.0

    @pytest.mark.unit
    def test_house_pos_integer_part_is_house_number(self):
        """Integer part of result should be house number (1-12)."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393
        lon = 15.0

        result = ephem.house_pos(armc, lat, eps, ord("P"), lon, 0.0)
        house_num = int(result)

        assert 1 <= house_num <= 12

    @pytest.mark.unit
    def test_house_pos_alias_matches(self):
        """house_pos should give same result as house_pos."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393
        lon = 15.0

        result1 = ephem.house_pos(armc, lat, eps, ord("P"), lon, 0.0)
        result2 = ephem.house_pos(armc, lat, eps, ord("P"), lon, 0.0)

        assert result1 == result2


class TestHousePosConsistency:
    """Test house_pos consistency with house cusps."""

    @pytest.mark.unit
    def test_body_at_cusp_start_has_near_zero_fraction(self):
        """Body exactly at cusp start should have near-zero fractional part."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        # Get the house cusps first
        cusps, ascmc = ephem.houses_armc(armc, lat, eps, ord("P"))

        # Place body at the start of house 1 (the ASC)
        lon = cusps[0]

        result = ephem.house_pos(armc, lat, eps, ord("P"), lon, 0.0)
        house_num = int(result)
        fraction = result - house_num

        # Should be in house 1 with near-zero fraction
        assert house_num == 1
        assert fraction < 0.001, f"Expected fraction near 0, got {fraction}"

    @pytest.mark.unit
    def test_body_just_before_next_cusp(self):
        """Body just before next cusp should have fraction close to 1."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        # Get the house cusps
        cusps, ascmc = ephem.houses_armc(armc, lat, eps, ord("P"))

        # Place body just before house 2 cusp (still in house 1)
        # Calculate a position that's 99% through house 1
        cusp1 = cusps[0]
        cusp2 = cusps[1]
        house_size = (cusp2 - cusp1 + 360.0) % 360.0
        lon = (cusp1 + house_size * 0.99) % 360.0

        result = ephem.house_pos(armc, lat, eps, ord("P"), lon, 0.0)
        house_num = int(result)
        fraction = result - house_num

        # Should be in house 1 with fraction close to 1
        assert house_num == 1
        assert fraction > 0.95, f"Expected fraction near 1, got {fraction}"


class TestHousePosMultipleSystems:
    """Test house_pos with different house systems."""

    HOUSE_SYSTEMS = ["P", "K", "R", "C", "E", "W", "O", "B"]

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys", HOUSE_SYSTEMS)
    def test_house_pos_all_systems_return_valid(self, hsys):
        """house_pos should return valid values for all major house systems."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393
        lon = 45.0  # 15 Taurus

        result = ephem.house_pos(armc, lat, eps, ord(hsys), lon, 0.0)

        assert 1.0 <= result < 13.0, f"Invalid result {result} for system {hsys}"


class TestHousePosComparisonWithSwisseph:
    """Compare house_pos results with Swiss Ephemeris."""

    # Test cases: (armc, lat, eps, lon, lat_body, hsys)
    TEST_CASES = [
        # Rome, Placidus, Sun at 15 Aries
        (292.957, 41.9, 23.4393, 15.0, 0.0, "P"),
        # London, Koch, Sun at 45 degrees
        (180.0, 51.5, 23.4393, 45.0, 0.0, "K"),
        # New York, Regiomontanus, Sun at 120 degrees
        (0.0, 40.7, 23.4393, 120.0, 0.0, "R"),
        # Equator, Placidus, Sun at 270 degrees
        (90.0, 0.0, 23.4393, 270.0, 0.0, "P"),
        # Southern hemisphere, Campanus
        (45.0, -33.9, 23.4393, 180.0, 0.0, "C"),
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("armc,lat,eps,lon,lat_body,hsys", TEST_CASES)
    def test_house_pos_matches_swisseph(self, armc, lat, eps, lon, lat_body, hsys):
        """house_pos should match Swiss Ephemeris house_pos."""
        # Get result from libephemeris
        result_lib = ephem.house_pos(armc, lat, eps, ord(hsys), lon, lat_body)

        # Get result from Swiss Ephemeris
        # swe.house_pos(armc, geolat, eps, objcoord, hsys)
        result_swe = swe.house_pos(armc, lat, eps, (lon, lat_body), hsys.encode())

        # Compare house numbers
        house_lib = int(result_lib)
        house_swe = int(result_swe)

        assert house_lib == house_swe, (
            f"House mismatch: libephemeris={house_lib}, swisseph={house_swe}"
        )

        # Compare fractions with tolerance
        frac_lib = result_lib - house_lib
        frac_swe = result_swe - house_swe

        # Allow 6% tolerance for fractional position
        # (house calculations can have minor differences between implementations)
        assert abs(frac_lib - frac_swe) < 0.06, (
            f"Fraction mismatch: libephemeris={frac_lib:.4f}, swisseph={frac_swe:.4f}"
        )


class TestHousePosEdgeCases:
    """Test edge cases for house_pos."""

    @pytest.mark.unit
    def test_house_pos_longitude_at_360_normalized(self):
        """Longitude at 360 should be normalized to 0."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        result1 = ephem.house_pos(armc, lat, eps, ord("P"), 0.0, 0.0)
        result2 = ephem.house_pos(armc, lat, eps, ord("P"), 360.0, 0.0)

        assert result1 == result2

    @pytest.mark.unit
    def test_house_pos_negative_longitude(self):
        """Negative longitude should be normalized."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        # -30 should be same as 330
        result1 = ephem.house_pos(armc, lat, eps, ord("P"), -30.0, 0.0)
        result2 = ephem.house_pos(armc, lat, eps, ord("P"), 330.0, 0.0)

        assert result1 == result2

    @pytest.mark.unit
    def test_house_pos_with_ecliptic_latitude(self):
        """house_pos should accept non-zero ecliptic latitude."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393
        lon = 45.0

        # Should not raise exception for non-zero latitude
        result = ephem.house_pos(armc, lat, eps, ord("P"), lon, 5.0)

        assert 1.0 <= result < 13.0


class TestHousePosAllHouses:
    """Test that bodies can be found in all 12 houses."""

    @pytest.mark.unit
    def test_can_find_body_in_each_house(self):
        """Should be able to place a body in each of the 12 houses."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        # Get cusps
        cusps, ascmc = ephem.houses_armc(armc, lat, eps, ord("P"))

        houses_found = set()

        for i in range(12):
            # Place body in middle of each house
            cusp_start = cusps[i]
            cusp_end = cusps[(i + 1) % 12]
            house_size = (cusp_end - cusp_start + 360.0) % 360.0
            lon = (cusp_start + house_size * 0.5) % 360.0

            result = ephem.house_pos(armc, lat, eps, ord("P"), lon, 0.0)
            house_num = int(result)
            houses_found.add(house_num)

        # Should have found all 12 houses
        assert houses_found == set(range(1, 13)), (
            f"Missing houses: {set(range(1, 13)) - houses_found}"
        )


class TestHousePosTypeSignatures:
    """Test type signature overloads for house_pos and house_pos.

    Tests that all supported calling conventions work correctly:
    1. 5-arg pyswisseph form: house_pos(armc, lat, obliquity, objcoord, hsys)
       where objcoord is tuple (lon, lat_body) and hsys is bytes/str
    2. 6-arg extended form: house_pos(armc, lat, obliquity, hsys, lon, lat_body)
       where hsys is int/bytes/str
    """

    # Standard test parameters
    ARMC = 292.957
    LAT = 41.9
    EPS = 23.4393
    LON = 15.0
    LAT_BODY = 0.0

    @pytest.mark.unit
    def test_6arg_form_hsys_as_int(self):
        """6-arg form with hsys as int (ord('P')) should work."""
        result = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, ord("P"), self.LON, self.LAT_BODY
        )
        assert isinstance(result, float)
        assert 1.0 <= result < 13.0

    @pytest.mark.unit
    def test_6arg_form_hsys_as_bytes(self):
        """6-arg form with hsys as bytes should work."""
        result = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, b"P", self.LON, self.LAT_BODY
        )
        assert isinstance(result, float)
        assert 1.0 <= result < 13.0

    @pytest.mark.unit
    def test_6arg_form_hsys_as_str(self):
        """6-arg form with hsys as str should work."""
        result = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, "P", self.LON, self.LAT_BODY
        )
        assert isinstance(result, float)
        assert 1.0 <= result < 13.0

    @pytest.mark.unit
    def test_5arg_form_objcoord_tuple_hsys_bytes(self):
        """5-arg pyswisseph form with tuple objcoord and bytes hsys should work."""
        result = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, (self.LON, self.LAT_BODY), b"P"
        )
        assert isinstance(result, float)
        assert 1.0 <= result < 13.0

    @pytest.mark.unit
    def test_5arg_form_objcoord_tuple_hsys_str(self):
        """5-arg pyswisseph form with tuple objcoord and str hsys should work."""
        result = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, (self.LON, self.LAT_BODY), "P"
        )
        assert isinstance(result, float)
        assert 1.0 <= result < 13.0

    @pytest.mark.unit
    def test_4arg_form_tuple_default_hsys(self):
        """4-arg form with tuple objcoord and default hsys should work."""
        result = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, (self.LON, self.LAT_BODY)
        )
        assert isinstance(result, float)
        assert 1.0 <= result < 13.0

    @pytest.mark.unit
    def test_calling_conventions_produce_same_result(self):
        """All calling conventions should produce the same result."""
        # 6-arg form with int hsys
        result_6arg_int = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, ord("P"), self.LON, self.LAT_BODY
        )

        # 6-arg form with bytes hsys
        result_6arg_bytes = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, b"P", self.LON, self.LAT_BODY
        )

        # 6-arg form with str hsys
        result_6arg_str = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, "P", self.LON, self.LAT_BODY
        )

        # 5-arg form with tuple and bytes
        result_5arg_bytes = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, (self.LON, self.LAT_BODY), b"P"
        )

        # 5-arg form with tuple and str
        result_5arg_str = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, (self.LON, self.LAT_BODY), "P"
        )

        # All should be equal
        assert result_6arg_int == result_6arg_bytes
        assert result_6arg_bytes == result_6arg_str
        assert result_6arg_str == result_5arg_bytes
        assert result_5arg_bytes == result_5arg_str

    @pytest.mark.unit
    def test_swe_house_pos_same_as_house_pos_6arg_int(self):
        """house_pos should match house_pos for 6-arg int form."""
        result_hp = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, ord("P"), self.LON, self.LAT_BODY
        )
        result_shp = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, ord("P"), self.LON, self.LAT_BODY
        )
        assert result_hp == result_shp

    @pytest.mark.unit
    def test_swe_house_pos_same_as_house_pos_5arg_tuple(self):
        """house_pos should match house_pos for 5-arg tuple form."""
        result_hp = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, (self.LON, self.LAT_BODY), b"P"
        )
        result_shp = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, (self.LON, self.LAT_BODY), b"P"
        )
        assert result_hp == result_shp

    @pytest.mark.unit
    def test_pyswisseph_compatible_5arg_form(self):
        """5-arg form should match pyswisseph swe.house_pos signature."""
        # This is the exact signature pyswisseph uses:
        # swe.house_pos(armc, geolat, eps, objcoord, hsys)
        result_lib = ephem.house_pos(
            self.ARMC, self.LAT, self.EPS, (self.LON, self.LAT_BODY), b"P"
        )
        result_swe = swe.house_pos(
            self.ARMC, self.LAT, self.EPS, (self.LON, self.LAT_BODY), b"P"
        )

        # House numbers should match
        assert int(result_lib) == int(result_swe)

        # Fractions should be close
        frac_lib = result_lib - int(result_lib)
        frac_swe = result_swe - int(result_swe)
        assert abs(frac_lib - frac_swe) < 0.06

    @pytest.mark.unit
    def test_multiple_house_systems_all_forms(self):
        """All calling conventions should work with different house systems."""
        house_systems = ["P", "K", "R", "C", "E", "W"]

        for hsys in house_systems:
            # 6-arg int form
            result1 = ephem.house_pos(
                self.ARMC, self.LAT, self.EPS, ord(hsys), self.LON, self.LAT_BODY
            )
            assert 1.0 <= result1 < 13.0, f"Failed for {hsys} (int form)"

            # 5-arg tuple form with bytes
            result2 = ephem.house_pos(
                self.ARMC, self.LAT, self.EPS, (self.LON, self.LAT_BODY), hsys.encode()
            )
            assert 1.0 <= result2 < 13.0, f"Failed for {hsys} (tuple/bytes form)"

            # Results should match
            assert int(result1) == int(result2), (
                f"House mismatch for {hsys}: int form={int(result1)}, tuple form={int(result2)}"
            )
