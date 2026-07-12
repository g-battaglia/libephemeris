"""
Tests for planetary moon support in libephemeris.

These tests verify:
- Moon ID constants and NAIF ID mappings
- Moon SPK registration and unregistration
- Moon name lookup
- is_planetary_moon() function
- Position calculation (when SPK is available)
- Error handling for unregistered moons

Note: Integration tests that require satellite SPK files (jup365.bsp, etc.)
are skipped by default. Set LIBEPHEMERIS_TEST_MOON_SPK=1 to run them.
"""

import os
import pytest

import libephemeris as eph
from libephemeris import planetary_moons, state
from libephemeris.constants import (
    MOON_OFFSET,
    MOON_IO,
    MOON_EUROPA,
    MOON_GANYMEDE,
    MOON_CALLISTO,
    MOON_TITAN,
    MOON_ENCELADUS,
    MOON_TRITON,
    MOON_PHOBOS,
    MOON_DEIMOS,
    MOON_CHARON,
    SUN,
    MARS,
    FLG_SPEED,
    NAIF_IO,
    NAIF_EUROPA,
    NAIF_GANYMEDE,
    NAIF_CALLISTO,
    NAIF_TITAN,
    NAIF_TRITON,
    NAIF_PHOBOS,
    NAIF_CHARON,
)


class TestMoonConstants:
    """Test planetary moon constant definitions."""

    def test_moon_offset(self):
        """Verify MOON_OFFSET is correctly defined."""
        assert MOON_OFFSET == 9000

    def test_galilean_moon_ids(self):
        """Verify Galilean moon IDs are correctly defined."""
        assert MOON_IO == MOON_OFFSET + 1  # 9001
        assert MOON_EUROPA == MOON_OFFSET + 2  # 9002
        assert MOON_GANYMEDE == MOON_OFFSET + 3  # 9003
        assert MOON_CALLISTO == MOON_OFFSET + 4  # 9004

    def test_saturn_moon_ids(self):
        """Verify Saturn moon IDs are correctly defined."""
        assert MOON_TITAN == MOON_OFFSET + 16  # 9016
        assert MOON_ENCELADUS == MOON_OFFSET + 12  # 9012

    def test_other_moon_ids(self):
        """Verify other moon IDs are correctly defined."""
        assert MOON_TRITON == MOON_OFFSET + 31  # 9031
        assert MOON_PHOBOS == MOON_OFFSET + 41  # 9041
        assert MOON_DEIMOS == MOON_OFFSET + 42  # 9042
        assert MOON_CHARON == MOON_OFFSET + 51  # 9051


class TestCanonicalMoonIds:
    """Canonical scheme: ipl = PLMOON_OFFSET + NAIF id (reference parity)."""

    def test_canonical_constants(self):
        """PLMOON_* constants encode PLMOON_OFFSET + NAIF id."""
        assert eph.PLMOON_OFFSET == 9000
        assert eph.PLMOON_IO == 9501
        assert eph.PLMOON_TITAN == 9606
        assert eph.PLMOON_CHARON == 9901
        assert eph.PLMOON_PHOBOS == 9401

    def test_canonical_ids_are_moons(self):
        """is_planetary_moon accepts canonical ids, including unnamed ones."""
        assert planetary_moons.is_planetary_moon(9501) is True  # Io
        assert planetary_moons.is_planetary_moon(9606) is True  # Titan
        assert planetary_moons.is_planetary_moon(9505) is True  # Amalthea
        assert planetary_moons.is_planetary_moon(9599) is True  # Jupiter COB
        # Below the NAIF satellite range: not moons
        assert planetary_moons.is_planetary_moon(9000) is False
        assert planetary_moons.is_planetary_moon(9400) is False

    def test_legacy_and_canonical_resolve_same_naif(self):
        """Legacy aliases and canonical ids map to the same NAIF target."""
        assert planetary_moons._moon_naif_id(MOON_IO) == 501
        assert planetary_moons._moon_naif_id(9501) == 501
        assert planetary_moons._moon_naif_id(MOON_TITAN) == 606
        assert planetary_moons._moon_naif_id(9606) == 606

    def test_canonical_names(self):
        """get_moon_name and get_planet_name resolve canonical ids."""
        assert planetary_moons.get_moon_name(9501) == "Io"
        assert planetary_moons.get_moon_name(9606) == "Titan"
        assert eph.get_planet_name(9501) == "Io"
        assert eph.get_planet_name(MOON_IO) == "Io"


class TestNaifMapping:
    """Test NAIF ID mappings for planetary moons."""

    def test_naif_id_values(self):
        """Verify NAIF ID values are correct."""
        assert NAIF_IO == 501
        assert NAIF_EUROPA == 502
        assert NAIF_GANYMEDE == 503
        assert NAIF_CALLISTO == 504
        assert NAIF_TITAN == 606
        assert NAIF_TRITON == 801
        assert NAIF_PHOBOS == 401
        assert NAIF_CHARON == 901

    def test_moon_naif_map_contains_galilean(self):
        """Verify MOON_NAIF_MAP contains Galilean moons."""
        from libephemeris.planetary_moons import MOON_NAIF_MAP

        assert MOON_IO in MOON_NAIF_MAP
        assert MOON_EUROPA in MOON_NAIF_MAP
        assert MOON_GANYMEDE in MOON_NAIF_MAP
        assert MOON_CALLISTO in MOON_NAIF_MAP

        assert MOON_NAIF_MAP[MOON_IO] == NAIF_IO
        assert MOON_NAIF_MAP[MOON_EUROPA] == NAIF_EUROPA
        assert MOON_NAIF_MAP[MOON_GANYMEDE] == NAIF_GANYMEDE
        assert MOON_NAIF_MAP[MOON_CALLISTO] == NAIF_CALLISTO


class TestMoonNames:
    """Test moon name lookup functionality."""

    def test_get_moon_name_galilean(self):
        """Test getting names of Galilean moons."""
        assert planetary_moons.get_moon_name(MOON_IO) == "Io"
        assert planetary_moons.get_moon_name(MOON_EUROPA) == "Europa"
        assert planetary_moons.get_moon_name(MOON_GANYMEDE) == "Ganymede"
        assert planetary_moons.get_moon_name(MOON_CALLISTO) == "Callisto"

    def test_get_moon_name_saturn(self):
        """Test getting names of Saturn moons."""
        assert planetary_moons.get_moon_name(MOON_TITAN) == "Titan"
        assert planetary_moons.get_moon_name(MOON_ENCELADUS) == "Enceladus"

    def test_get_moon_name_other(self):
        """Test getting names of other moons."""
        assert planetary_moons.get_moon_name(MOON_TRITON) == "Triton"
        assert planetary_moons.get_moon_name(MOON_PHOBOS) == "Phobos"
        assert planetary_moons.get_moon_name(MOON_DEIMOS) == "Deimos"
        assert planetary_moons.get_moon_name(MOON_CHARON) == "Charon"

    def test_get_moon_name_unknown(self):
        """Test getting name for unknown moon ID."""
        unknown_id = 99999
        name = planetary_moons.get_moon_name(unknown_id)
        assert "Unknown" in name
        assert str(unknown_id) in name


class TestIsPlanetaryMoon:
    """Test is_planetary_moon() function."""

    def test_is_planetary_moon_true(self):
        """Test that planetary moon IDs return True."""
        assert planetary_moons.is_planetary_moon(MOON_IO) is True
        assert planetary_moons.is_planetary_moon(MOON_TITAN) is True
        assert planetary_moons.is_planetary_moon(MOON_TRITON) is True
        assert planetary_moons.is_planetary_moon(MOON_PHOBOS) is True

    def test_is_planetary_moon_false(self):
        """Test that non-moon IDs return False."""
        assert planetary_moons.is_planetary_moon(SUN) is False
        assert planetary_moons.is_planetary_moon(MARS) is False
        assert planetary_moons.is_planetary_moon(0) is False
        assert planetary_moons.is_planetary_moon(99999) is False


class TestMoonRegistration:
    """Test moon SPK registration without actual SPK files."""

    def setup_method(self):
        """Reset state before each test."""
        planetary_moons.close_moon_kernels()

    def teardown_method(self):
        """Clean up after each test."""
        planetary_moons.close_moon_kernels()

    def test_list_registered_moons_empty(self):
        """Test that no moons are registered initially."""
        moons = planetary_moons.list_registered_moons()
        assert len(moons) == 0

    def test_register_moon_spk_file_not_found(self):
        """Test error when SPK file doesn't exist."""
        with pytest.raises(FileNotFoundError):
            planetary_moons.register_moon_spk("nonexistent_file.bsp")

    def test_get_moon_coverage_unregistered(self):
        """Test coverage returns None for unregistered moon."""
        coverage = planetary_moons.get_moon_coverage(MOON_IO)
        assert coverage is None

    def test_calc_moon_position_unregistered(self):
        """Test position returns None for unregistered moon."""
        ts = state.get_timescale()
        t = ts.tt_jd(2451545.0)
        result = planetary_moons.calc_moon_position(t, MOON_IO, FLG_SPEED)
        assert result is None


class TestMoonPositionUnregistered:
    """Test moon position calculation via calc_ut when not registered."""

    def setup_method(self):
        """Reset state before each test."""
        planetary_moons.close_moon_kernels()
        eph.close()

    def teardown_method(self):
        """Clean up after each test."""
        planetary_moons.close_moon_kernels()

    def test_calc_ut_unregistered_moon_raises(self):
        """calc_ut raises Error for an unregistered moon (reference parity).

        Black-box reference calls also raise when satellite ephemeris coverage
        is unavailable; silent zeros would be indistinguishable from a real
        position. No reference-distribution file is needed by this test.
        """
        jd = 2451545.0  # J2000.0
        with pytest.raises(eph.Error):
            eph.calc_ut(jd, MOON_IO, FLG_SPEED)
        # Canonical id (PLMOON_OFFSET + NAIF, Io = 9501) behaves the same
        with pytest.raises(eph.Error):
            eph.calc_ut(jd, 9501, FLG_SPEED)


class TestMoonParentMapping:
    """Test parent planet mapping for moons."""

    def test_galilean_moons_parent_jupiter(self):
        """Test that Galilean moons have Jupiter as parent."""
        from libephemeris.planetary_moons import (
            MOON_PARENT_MAP,
            NAIF_JUPITER_BARYCENTER,
        )

        assert MOON_PARENT_MAP[MOON_IO] == NAIF_JUPITER_BARYCENTER
        assert MOON_PARENT_MAP[MOON_EUROPA] == NAIF_JUPITER_BARYCENTER
        assert MOON_PARENT_MAP[MOON_GANYMEDE] == NAIF_JUPITER_BARYCENTER
        assert MOON_PARENT_MAP[MOON_CALLISTO] == NAIF_JUPITER_BARYCENTER

    def test_titan_parent_saturn(self):
        """Test that Titan has Saturn as parent."""
        from libephemeris.planetary_moons import MOON_PARENT_MAP, NAIF_SATURN_BARYCENTER

        assert MOON_PARENT_MAP[MOON_TITAN] == NAIF_SATURN_BARYCENTER

    def test_triton_parent_neptune(self):
        """Test that Triton has Neptune as parent."""
        from libephemeris.planetary_moons import (
            MOON_PARENT_MAP,
            NAIF_NEPTUNE_BARYCENTER,
        )

        assert MOON_PARENT_MAP[MOON_TRITON] == NAIF_NEPTUNE_BARYCENTER

    def test_mars_moons_parent_mars(self):
        """Test that Mars moons have Mars as parent."""
        from libephemeris.planetary_moons import MOON_PARENT_MAP, NAIF_MARS_BARYCENTER

        assert MOON_PARENT_MAP[MOON_PHOBOS] == NAIF_MARS_BARYCENTER
        assert MOON_PARENT_MAP[MOON_DEIMOS] == NAIF_MARS_BARYCENTER


# =============================================================================
# INTEGRATION TESTS - Require satellite SPK files
# =============================================================================
# These tests require actual satellite SPK files to be available.
# Skip by default unless LIBEPHEMERIS_TEST_MOON_SPK=1 is set.


def has_moon_spk_testing():
    """Check if moon SPK testing is enabled."""
    return os.environ.get("LIBEPHEMERIS_TEST_MOON_SPK", "").lower() in (
        "1",
        "true",
        "yes",
    )


def get_jupiter_spk_path():
    """Get path to Jupiter satellite SPK file if available."""
    # Try common locations
    candidates = [
        "jup365.bsp",
        os.path.expanduser("~/.libephemeris/spk/jup365.bsp"),
        os.path.join(os.path.dirname(__file__), "..", "jup365.bsp"),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    return None


@pytest.mark.skipif(
    not has_moon_spk_testing(),
    reason="Moon SPK testing disabled. Set LIBEPHEMERIS_TEST_MOON_SPK=1 to enable.",
)
class TestMoonPositionIntegration:
    """Integration tests requiring satellite SPK files."""

    def setup_method(self):
        """Reset state before each test."""
        planetary_moons.close_moon_kernels()
        eph.close()

    def teardown_method(self):
        """Clean up after each test."""
        planetary_moons.close_moon_kernels()

    def test_register_jupiter_moons(self):
        """Test registering Jupiter satellite SPK."""
        spk_path = get_jupiter_spk_path()
        if spk_path is None:
            pytest.skip("Jupiter SPK file not found")

        planetary_moons.register_moon_spk(spk_path)
        moons = planetary_moons.list_registered_moons()

        # Should have at least the Galilean moons
        assert MOON_IO in moons or len(moons) > 0

    def test_calc_io_position(self):
        """Test calculating Io's position."""
        spk_path = get_jupiter_spk_path()
        if spk_path is None:
            pytest.skip("Jupiter SPK file not found")

        planetary_moons.register_moon_spk(spk_path)

        jd = 2451545.0  # J2000.0
        pos, flag = eph.calc_ut(jd, MOON_IO, FLG_SPEED)

        # Should have non-zero position
        assert pos[2] > 0  # Distance should be positive

        # Io should be roughly 4-6 AU from Earth (Jupiter's distance)
        assert 3.5 < pos[2] < 7.0, f"Io distance {pos[2]:.2f} AU unexpected"

    def test_calc_io_with_speed(self):
        """Test calculating Io's position with velocity."""
        spk_path = get_jupiter_spk_path()
        if spk_path is None:
            pytest.skip("Jupiter SPK file not found")

        planetary_moons.register_moon_spk(spk_path)

        jd = 2451545.0
        pos, flag = eph.calc_ut(jd, MOON_IO, FLG_SPEED)

        # Io has fast orbital motion (~17 hours period)
        # Should have significant daily motion
        assert pos[3] != 0.0, "Io should have non-zero velocity"


class TestMoonConstantsExport:
    """Test that moon constants are properly exported from libephemeris."""

    def test_moon_constants_accessible(self):
        """Test that moon constants are accessible from libephemeris."""
        assert hasattr(eph, "MOON_IO")
        assert hasattr(eph, "MOON_EUROPA")
        assert hasattr(eph, "MOON_GANYMEDE")
        assert hasattr(eph, "MOON_CALLISTO")
        assert hasattr(eph, "MOON_TITAN")
        assert hasattr(eph, "MOON_TRITON")

    def test_moon_functions_accessible(self):
        """Test that moon functions are accessible from libephemeris."""
        assert hasattr(eph, "register_moon_spk")
        assert hasattr(eph, "unregister_moon_spk")
        assert hasattr(eph, "list_registered_moons")
        assert hasattr(eph, "get_moon_name")
        assert hasattr(eph, "is_planetary_moon")

    def test_naif_constants_accessible(self):
        """Test that NAIF ID constants are accessible from libephemeris."""
        assert hasattr(eph, "NAIF_IO")
        assert hasattr(eph, "NAIF_TITAN")
        assert hasattr(eph, "NAIF_TRITON")


class TestCloseFunction:
    """Test that close() properly cleans up moon resources."""

    def setup_method(self):
        """Reset state before each test."""
        planetary_moons.close_moon_kernels()

    def test_close_clears_moon_registrations(self):
        """Test that eph.close() clears moon registrations."""
        # Add a mock registration
        planetary_moons._MOON_SPK_BY_BODY[MOON_IO] = "test.bsp"

        # Close should clear it
        planetary_moons.close_moon_kernels()

        moons = planetary_moons.list_registered_moons()
        assert len(moons) == 0
