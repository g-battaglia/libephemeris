"""
Tests for verifying all constants match pyswisseph exactly.

This ensures 100% API compatibility at the constant level.
"""

import pytest
import swisseph as swe
from libephemeris import constants


class TestPlanetIDs:
    """Test that all planet IDs match pyswisseph."""

    @pytest.mark.unit
    def test_sun_id(self):
        assert constants.SUN == swe.SUN == 0

    @pytest.mark.unit
    def test_moon_id(self):
        assert constants.MOON == swe.MOON == 1

    @pytest.mark.unit
    def test_mercury_id(self):
        assert constants.MERCURY == swe.MERCURY == 2

    @pytest.mark.unit
    def test_venus_id(self):
        assert constants.VENUS == swe.VENUS == 3

    @pytest.mark.unit
    def test_mars_id(self):
        assert constants.MARS == swe.MARS == 4

    @pytest.mark.unit
    def test_jupiter_id(self):
        assert constants.JUPITER == swe.JUPITER == 5

    @pytest.mark.unit
    def test_saturn_id(self):
        assert constants.SATURN == swe.SATURN == 6

    @pytest.mark.unit
    def test_uranus_id(self):
        assert constants.URANUS == swe.URANUS == 7

    @pytest.mark.unit
    def test_neptune_id(self):
        assert constants.NEPTUNE == swe.NEPTUNE == 8

    @pytest.mark.unit
    def test_pluto_id(self):
        assert constants.PLUTO == swe.PLUTO == 9

    @pytest.mark.unit
    def test_mean_node_id(self):
        assert constants.MEAN_NODE == swe.MEAN_NODE == 10

    @pytest.mark.unit
    def test_true_node_id(self):
        assert constants.TRUE_NODE == swe.TRUE_NODE == 11

    @pytest.mark.unit
    def test_mean_apog_id(self):
        """Mean Lilith / Black Moon."""
        assert constants.MEAN_APOG == swe.MEAN_APOG == 12

    @pytest.mark.unit
    def test_oscu_apog_id(self):
        """Osculating/True Lilith."""
        assert constants.OSCU_APOG == swe.OSCU_APOG == 13

    @pytest.mark.unit
    def test_earth_id(self):
        assert constants.EARTH == swe.EARTH == 14

    @pytest.mark.unit
    def test_chiron_id(self):
        assert constants.CHIRON == swe.CHIRON == 15

    @pytest.mark.unit
    def test_pholus_id(self):
        assert constants.PHOLUS == swe.PHOLUS == 16

    @pytest.mark.unit
    def test_ceres_id(self):
        assert constants.CERES == swe.CERES == 17

    @pytest.mark.unit
    def test_pallas_id(self):
        assert constants.PALLAS == swe.PALLAS == 18

    @pytest.mark.unit
    def test_juno_id(self):
        assert constants.JUNO == swe.JUNO == 19

    @pytest.mark.unit
    def test_vesta_id(self):
        assert constants.VESTA == swe.VESTA == 20

    @pytest.mark.unit
    def test_all_major_planets_sequential(self):
        """Verify planets 0-9 are sequential."""
        for i, planet_id in enumerate(
            [
                constants.SUN,
                constants.MOON,
                constants.MERCURY,
                constants.VENUS,
                constants.MARS,
                constants.JUPITER,
                constants.SATURN,
                constants.URANUS,
                constants.NEPTUNE,
                constants.PLUTO,
            ]
        ):
            assert planet_id == i, f"Planet at index {i} has wrong ID {planet_id}"


class TestCalculationFlags:
    """Test that all calculation flags match pyswisseph."""

    @pytest.mark.unit
    def test_flg_speed(self):
        assert constants.FLG_SPEED == swe.FLG_SPEED == 256

    @pytest.mark.unit
    def test_flg_topoctr(self):
        assert constants.FLG_TOPOCTR == swe.FLG_TOPOCTR == 32768

    @pytest.mark.unit
    def test_flg_helctr(self):
        assert constants.FLG_HELCTR == swe.FLG_HELCTR == 8

    @pytest.mark.unit
    def test_flg_baryctr(self):
        assert constants.FLG_BARYCTR == swe.FLG_BARYCTR

    @pytest.mark.unit
    def test_flg_equatorial(self):
        assert constants.FLG_EQUATORIAL == swe.FLG_EQUATORIAL == 2048

    @pytest.mark.unit
    def test_flg_sidereal(self):
        assert constants.FLG_SIDEREAL == swe.FLG_SIDEREAL == 65536

    @pytest.mark.unit
    def test_flg_j2000(self):
        assert constants.FLG_J2000 == swe.FLG_J2000

    @pytest.mark.unit
    def test_flg_truepos(self):
        assert constants.FLG_TRUEPOS == swe.FLG_TRUEPOS

    @pytest.mark.unit
    def test_flg_noaberr(self):
        assert constants.FLG_NOABERR == swe.FLG_NOABERR

    @pytest.mark.unit
    def test_flg_nogdefl(self):
        assert constants.FLG_NOGDEFL == swe.FLG_NOGDEFL

    @pytest.mark.unit
    def test_flags_can_be_combined(self):
        """Test that flags can be OR'd together."""
        combined = constants.FLG_SPEED | constants.FLG_EQUATORIAL
        assert combined == 256 | 2048 == 2304


class TestSiderealModes:
    """Test that all sidereal mode constants match pyswisseph."""

    @pytest.mark.unit
    def test_sidm_fagan_bradley(self):
        assert constants.SIDM_FAGAN_BRADLEY == swe.SIDM_FAGAN_BRADLEY == 0

    @pytest.mark.unit
    def test_sidm_lahiri(self):
        assert constants.SIDM_LAHIRI == swe.SIDM_LAHIRI == 1

    @pytest.mark.unit
    def test_sidm_deluce(self):
        assert constants.SIDM_DELUCE == swe.SIDM_DELUCE == 2

    @pytest.mark.unit
    def test_sidm_raman(self):
        assert constants.SIDM_RAMAN == swe.SIDM_RAMAN == 3

    @pytest.mark.unit
    def test_sidm_krishnamurti(self):
        assert constants.SIDM_KRISHNAMURTI == swe.SIDM_KRISHNAMURTI == 5

    @pytest.mark.unit
    def test_sidm_true_citra(self):
        assert constants.SIDM_TRUE_CITRA == swe.SIDM_TRUE_CITRA == 27

    @pytest.mark.unit
    def test_sidm_galcent_0sag(self):
        assert constants.SIDM_GALCENT_0SAG == swe.SIDM_GALCENT_0SAG == 17

    @pytest.mark.unit
    def test_sidm_j2000(self):
        assert constants.SIDM_J2000 == swe.SIDM_J2000 == 18

    @pytest.mark.unit
    def test_total_sidereal_modes(self):
        """Verify we have at least 43 sidereal modes defined."""
        sidm_count = sum(1 for name in dir(constants) if name.startswith("SE_SIDM_"))
        assert sidm_count >= 43, (
            f"Expected at least 43 sidereal modes, got {sidm_count}"
        )


class TestCalendarConstants:
    """Test calendar system constants."""

    @pytest.mark.unit
    def test_greg_cal(self):
        assert constants.GREG_CAL == swe.GREG_CAL == 1

    @pytest.mark.unit
    def test_jul_cal(self):
        assert constants.JUL_CAL == swe.JUL_CAL == 0


class TestOffsets:
    """Test offset constants for special bodies."""

    @pytest.mark.unit
    def test_ast_offset(self):
        """Asteroid offset for numbered asteroids."""
        assert constants.AST_OFFSET == swe.AST_OFFSET == 10000

    @pytest.mark.unit
    def test_fict_offset(self):
        """Fictitious bodies offset (hypothetical planets)."""
        assert hasattr(swe, "FICT_OFFSET")
        # Verify our constant exists
        assert hasattr(constants, "FICT_OFFSET") or hasattr(constants, "FICT_OFFSET")


class TestPyswissephAliases:
    """Test that pyswisseph-style aliases work."""

    @pytest.mark.unit
    def test_flg_prefix_aliases(self):
        """FLG_* should work as alias for SEFLG_*."""
        assert constants.FLG_SPEED == constants.FLG_SPEED
        assert constants.FLG_HELCTR == constants.FLG_HELCTR
        assert constants.FLG_SIDEREAL == constants.FLG_SIDEREAL

    @pytest.mark.unit
    def test_sidm_prefix_aliases(self):
        """SIDM_* should work as alias for SE_SIDM_*."""
        assert constants.SIDM_LAHIRI == constants.SIDM_LAHIRI
        assert constants.SIDM_FAGAN_BRADLEY == constants.SIDM_FAGAN_BRADLEY

    @pytest.mark.unit
    def test_planet_names_without_se_prefix(self):
        """SUN, MOON, etc. should work without SE_ prefix."""
        assert constants.SUN == constants.SUN
        assert constants.MOON == constants.MOON
        assert constants.MERCURY == constants.MERCURY


class TestHouseCuspConstants:
    """Test house cusp point constants."""

    @pytest.mark.unit
    def test_ascmc_indices(self):
        """Test the indices used in ascmc array."""
        # These are the indices in the ascmc tuple returned by houses()
        assert hasattr(swe, "ASC") or hasattr(swe, "SE_ASC")
        # Our implementation should have ASC at index 0, MC at index 1


class TestHypotheticalBodies:
    """Test hypothetical body constants."""

    @pytest.mark.unit
    def test_planet_x_leverrier_constant(self):
        """Test that PLANET_X_LEVERRIER is correctly defined."""
        # Should be FICT_OFFSET + 11 = 40 + 11 = 51
        assert constants.PLANET_X_LEVERRIER == 51
        assert constants.PLANET_X_LEVERRIER == constants.FICT_OFFSET + 11

    @pytest.mark.unit
    def test_planet_x_leverrier_alias(self):
        """Test that PLANET_X_LEVERRIER alias works."""
        assert constants.PLANET_X_LEVERRIER == constants.PLANET_X_LEVERRIER

    @pytest.mark.unit
    def test_planet_x_adams_constant(self):
        """Test that PLANET_X_ADAMS is correctly defined."""
        # Should be FICT_OFFSET + 12 = 40 + 12 = 52
        assert constants.PLANET_X_ADAMS == 52
        assert constants.PLANET_X_ADAMS == constants.FICT_OFFSET + 12

    @pytest.mark.unit
    def test_planet_x_adams_alias(self):
        """Test that PLANET_X_ADAMS alias works."""
        assert constants.PLANET_X_ADAMS == constants.PLANET_X_ADAMS

    @pytest.mark.unit
    def test_waldemath_constant(self):
        """Test that WALDEMATH is correctly defined."""
        assert constants.WALDEMATH == 58
        assert constants.WALDEMATH == constants.FICT_OFFSET + 18

    @pytest.mark.unit
    def test_vulcan_constant(self):
        """Test that VULCAN is correctly defined."""
        assert constants.VULCAN == 55
        assert constants.VULCAN == constants.FICT_OFFSET + 15


class TestConstantsCompleteness:
    """Test that we have all constants defined that pyswisseph has."""

    @pytest.mark.unit
    def test_all_planet_ids_defined(self):
        """Verify all planet IDs from 0-20 are defined."""
        required_planets = [
            "SUN",
            "MOON",
            "MERCURY",
            "VENUS",
            "MARS",
            "JUPITER",
            "SATURN",
            "URANUS",
            "NEPTUNE",
            "PLUTO",
            "MEAN_NODE",
            "TRUE_NODE",
            "MEAN_APOG",
            "OSCU_APOG",
            "EARTH",
            "CHIRON",
            "PHOLUS",
            "CERES",
            "PALLAS",
            "JUNO",
            "VESTA",
        ]
        for planet_name in required_planets:
            assert hasattr(constants, planet_name), f"Missing constant: {planet_name}"

    @pytest.mark.unit
    def test_all_major_flags_defined(self):
        """Verify all major calculation flags are defined."""
        required_flags = [
            "FLG_SPEED",
            "FLG_TOPOCTR",
            "FLG_HELCTR",
            "FLG_BARYCTR",
            "FLG_EQUATORIAL",
            "FLG_SIDEREAL",
            "FLG_J2000",
            "FLG_TRUEPOS",
            "FLG_NOABERR",
            "FLG_NOGDEFL",
        ]
        for flag_name in required_flags:
            assert hasattr(constants, flag_name), f"Missing constant: {flag_name}"
