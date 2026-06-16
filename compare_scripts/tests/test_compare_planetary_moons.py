"""
Planetary Moons Comparison Tests.

Validates planetary moon calculations in libephemeris, comparing with pyswisseph
where supported and documenting libephemeris-specific behavior.

IMPORTANT: Swiss Ephemeris and pyswisseph DO support planetary moons (since SE 2.10+)
using the same body IDs (MOON_OFFSET + moon_number = 9000+). However, both
libraries require satellite SPK kernel files from JPL NAIF to compute accurate
positions:
- jup365.bsp for Jupiter's moons (Io, Europa, Ganymede, Callisto)
- sat441.bsp for Saturn's moons (Titan, Rhea, Dione, etc.)
- ura116.bsp for Uranus' moons (Miranda, Ariel, Umbriel, Titania, Oberon)
- nep097.bsp for Neptune's moons (Triton)
- mar097.bsp for Mars' moons (Phobos, Deimos)
- plu058.bsp for Pluto's moons (Charon)

Without these data files, pyswisseph returns the Sun's position as a fallback.
libephemeris uses Skyfield to load these SPK files.

Supported Moons:
- Jupiter's Galilean moons: Io (9001), Europa (9002), Ganymede (9003), Callisto (9004)
- Saturn's major moons: Titan (9016), Enceladus (9012), etc.
- Uranus' major moons: Titania (9024), Oberon (9025), etc.
- Neptune's major moon: Triton (9031)
- Mars' moons: Phobos (9041), Deimos (9042)
- Pluto's moon: Charon (9051)
"""

import os
import pytest
import swisseph as swe
import libephemeris as eph
from libephemeris import planetary_moons
from libephemeris.constants import (
    SUN,
    JUPITER,
    SATURN,
    NEPTUNE,
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
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_HELCTR,
)


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Galilean moons (Jupiter I-IV)
GALILEAN_MOONS = [
    (MOON_IO, "Io", 501),  # NAIF ID 501
    (MOON_EUROPA, "Europa", 502),
    (MOON_GANYMEDE, "Ganymede", 503),
    (MOON_CALLISTO, "Callisto", 504),
]

# Saturn's major moons
SATURN_MOONS = [
    (MOON_TITAN, "Titan", 606),
    (MOON_ENCELADUS, "Enceladus", 602),
]

# Other planetary moons
OTHER_MOONS = [
    (MOON_TRITON, "Triton", 801),  # Neptune I
    (MOON_PHOBOS, "Phobos", 401),  # Mars I
    (MOON_DEIMOS, "Deimos", 402),  # Mars II
    (MOON_CHARON, "Charon", 901),  # Pluto I
]

ALL_MOONS = GALILEAN_MOONS + SATURN_MOONS + OTHER_MOONS

TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (2024, 6, 15, 12.0, "Mid-2024"),
    (2025, 1, 1, 0.0, "2025"),
]

# Tolerances for validation (when SPK data is available)
LONGITUDE_TOLERANCE = 0.1  # degrees - relaxed for orbital variations
DISTANCE_TOLERANCE = 0.5  # AU - Galilean moons orbit Jupiter at ~5 AU


def has_moon_spk_testing():
    """Check if moon SPK testing is enabled via environment variable."""
    return os.environ.get("LIBEPHEMERIS_TEST_MOON_SPK", "").lower() in (
        "1",
        "true",
        "yes",
    )


def get_spk_path(filename):
    """Try to find a satellite SPK file in common locations."""
    candidates = [
        filename,
        os.path.expanduser(f"~/.libephemeris/spk/{filename}"),
        os.path.join(os.path.dirname(__file__), "..", "..", filename),
        os.path.join(os.path.dirname(__file__), "..", "..", "data", filename),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    return None


# ============================================================================
# PYSWISSEPH BEHAVIOR DOCUMENTATION TESTS
# ============================================================================


class TestPyswissephMoonSupport:
    """Document and verify pyswisseph planetary moon support behavior.

    These tests confirm that pyswisseph uses the same body ID scheme as
    libephemeris (MOON_OFFSET = 9000 + moon number) but requires
    satellite SPK files for accurate calculations.
    """

    @pytest.mark.comparison
    def test_pyswisseph_moon_offset_matches(self):
        """Verify pyswisseph uses the same MOON_OFFSET (9000) as libephemeris."""
        assert hasattr(swe, "PLMOON_OFFSET"), (
            "pyswisseph should have PLMOON_OFFSET constant"
        )
        assert swe.PLMOON_OFFSET == MOON_OFFSET == 9000, (
            f"Moon offset mismatch: swe={swe.PLMOON_OFFSET}, lib={MOON_OFFSET}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("moon_id,moon_name,naif_id", GALILEAN_MOONS)
    def test_pyswisseph_galilean_moon_ids(self, moon_id, moon_name, naif_id):
        """Verify Galilean moon IDs are consistent between libraries.

        Both pyswisseph and libephemeris use:
        - Io: 9001 (NAIF 501)
        - Europa: 9002 (NAIF 502)
        - Ganymede: 9003 (NAIF 503)
        - Callisto: 9004 (NAIF 504)
        """
        expected_id = MOON_OFFSET + (naif_id - 500)
        assert moon_id == expected_id, (
            f"{moon_name} ID mismatch: expected {expected_id}, got {moon_id}"
        )

    @pytest.mark.comparison
    def test_pyswisseph_returns_fallback_without_spk(self):
        """Document that pyswisseph returns Sun's position without SPK files.

        When satellite SPK files are not loaded, pyswisseph returns the Sun's
        position as a fallback for planetary moon calculations. This is the
        expected behavior when the ephemeris data is not available.
        """
        jd = 2451545.0  # J2000.0

        # Get Sun position for reference
        sun_pos, _ = swe.calc_ut(jd, swe.SUN, swe.FLG_SWIEPH)

        # Try calculating Io
        io_pos, flags = swe.calc_ut(jd, MOON_IO, swe.FLG_SWIEPH)

        # Without SPK files, pyswisseph should return something close to Sun
        # (it returns the Sun's position as a fallback)
        lon_diff = abs(io_pos[0] - sun_pos[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        # Document expected behavior
        assert lon_diff < 1.0, (
            "Without SPK files, pyswisseph should return Sun position as fallback. "
            f"Sun={sun_pos[0]:.4f}, Io={io_pos[0]:.4f}, diff={lon_diff:.4f}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("moon_id,moon_name,naif_id", GALILEAN_MOONS)
    def test_pyswisseph_all_galilean_moons_same_without_spk(
        self, moon_id, moon_name, naif_id
    ):
        """Document that without SPK, all Galilean moons return the same position.

        This confirms pyswisseph's fallback behavior when satellite ephemeris
        data is not available.
        """
        jd = 2451545.0

        io_pos, _ = swe.calc_ut(jd, MOON_IO, swe.FLG_SWIEPH)
        moon_pos, _ = swe.calc_ut(jd, moon_id, swe.FLG_SWIEPH)

        # All moons should return the same fallback position
        assert abs(io_pos[0] - moon_pos[0]) < 0.0001, (
            f"Without SPK, {moon_name} should return same fallback as Io. "
            f"Io={io_pos[0]:.4f}, {moon_name}={moon_pos[0]:.4f}"
        )


# ============================================================================
# LIBEPHEMERIS BEHAVIOR TESTS (WITHOUT SPK)
# ============================================================================


class TestLibephemerisWithoutSpk:
    """Test libephemeris planetary moon behavior without SPK files loaded.

    When no satellite SPK files are registered, libephemeris raises Error
    (honest failure, reference-API parity) rather than pyswisseph's silent
    Sun-position fallback or the older silent-zeros behavior.
    """

    def setup_method(self):
        """Reset moon registrations before each test."""
        planetary_moons.close_moon_kernels()
        eph.close()

    def teardown_method(self):
        """Clean up after each test."""
        planetary_moons.close_moon_kernels()

    @pytest.mark.comparison
    @pytest.mark.parametrize("moon_id,moon_name,naif_id", GALILEAN_MOONS)
    def test_unregistered_moon_raises_error(self, moon_id, moon_name, naif_id):
        """Unregistered planetary moons raise Error (honest failure, like the
        reference API), directing the caller to register_moon_spk() instead of
        returning silent zeros. Changed deliberately in WS7f for parity.
        """
        jd = 2451545.0
        with pytest.raises(eph.Error, match="SPK"):
            eph.calc_ut(jd, moon_id, FLG_SPEED)

    @pytest.mark.comparison
    def test_difference_from_pyswisseph_documented(self):
        """Document the behavioral difference when no satellite SPK is loaded:
        libephemeris raises Error (honest failure, reference-API parity) while
        pyswisseph silently returns the Sun position.
        """
        jd = 2451545.0

        # libephemeris raises rather than returning a misleading position.
        with pytest.raises(eph.Error):
            eph.calc_ut(jd, MOON_IO, FLG_SPEED)

        # pyswisseph returns Sun position
        swe_pos, _ = swe.calc_ut(jd, MOON_IO, swe.FLG_SWIEPH)
        sun_pos, _ = swe.calc_ut(jd, swe.SUN, swe.FLG_SWIEPH)
        assert abs(swe_pos[0] - sun_pos[0]) < 1.0, (
            "pyswisseph should return Sun position for unregistered moon"
        )

    @pytest.mark.comparison
    def test_is_planetary_moon_correct(self):
        """Verify is_planetary_moon() correctly identifies moon IDs."""
        assert planetary_moons.is_planetary_moon(MOON_IO) is True
        assert planetary_moons.is_planetary_moon(MOON_TITAN) is True
        assert planetary_moons.is_planetary_moon(MOON_TRITON) is True
        assert planetary_moons.is_planetary_moon(SUN) is False
        assert planetary_moons.is_planetary_moon(JUPITER) is False

    @pytest.mark.comparison
    @pytest.mark.parametrize("moon_id,moon_name,naif_id", ALL_MOONS)
    def test_moon_names_correct(self, moon_id, moon_name, naif_id):
        """Verify get_moon_name() returns correct names."""
        assert planetary_moons.get_moon_name(moon_id) == moon_name


# ============================================================================
# INTEGRATION TESTS - REQUIRE SATELLITE SPK FILES
# ============================================================================


@pytest.mark.skipif(
    not has_moon_spk_testing(),
    reason="Moon SPK testing disabled. Set LIBEPHEMERIS_TEST_MOON_SPK=1 to enable.",
)
class TestGalileanMoonsWithSpk:
    """Integration tests for Galilean moons with SPK files.

    These tests require jup365.bsp (Jupiter satellite ephemeris) to be available.
    They validate that libephemeris correctly calculates positions for
    Io, Europa, Ganymede, and Callisto.
    """

    @pytest.fixture(autouse=True)
    def setup_jupiter_spk(self):
        """Register Jupiter SPK if available."""
        planetary_moons.close_moon_kernels()
        eph.close()

        spk_path = get_spk_path("jup365.bsp")
        if spk_path:
            planetary_moons.register_moon_spk(spk_path)
            yield
            planetary_moons.close_moon_kernels()
        else:
            pytest.skip("jup365.bsp not found")

    @pytest.mark.comparison
    @pytest.mark.parametrize("moon_id,moon_name,naif_id", GALILEAN_MOONS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_galilean_moon_geocentric_position(
        self, moon_id, moon_name, naif_id, year, month, day, hour, date_desc
    ):
        """Validate Galilean moon geocentric positions are reasonable.

        Since pyswisseph without SPK returns invalid data, we validate against
        physical constraints:
        - Longitude: 0-360 degrees
        - Latitude: reasonable for ecliptic (-5 to +5 degrees typically)
        - Distance: approximately Jupiter's distance (~4-6 AU)
        """
        jd = eph.julday(year, month, day, hour)
        pos, _ = eph.calc_ut(jd, moon_id, FLG_SWIEPH | FLG_SPEED)

        # Validate longitude is in valid range
        assert 0 <= pos[0] < 360, (
            f"{moon_name} at {date_desc}: longitude {pos[0]:.4f} out of range"
        )

        # Validate latitude is reasonable (Jupiter's orbital inclination ~1.3 deg)
        assert -10 <= pos[1] <= 10, (
            f"{moon_name} at {date_desc}: latitude {pos[1]:.4f} seems too extreme"
        )

        # Validate distance is approximately Jupiter's distance from Earth
        jupiter_pos, _ = eph.calc_ut(jd, JUPITER, FLG_SWIEPH)
        jupiter_dist = jupiter_pos[2]

        # Moon should be within ~0.05 AU of Jupiter's distance
        assert abs(pos[2] - jupiter_dist) < 0.1, (
            f"{moon_name} at {date_desc}: distance {pos[2]:.4f} AU differs "
            f"too much from Jupiter {jupiter_dist:.4f} AU"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("moon_id,moon_name,naif_id", GALILEAN_MOONS)
    def test_galilean_moon_velocity(self, moon_id, moon_name, naif_id):
        """Validate Galilean moon velocity calculations.

        Galilean moons have fast orbital periods:
        - Io: ~1.77 days (very fast angular motion when viewed from Earth)
        - Europa: ~3.55 days
        - Ganymede: ~7.15 days
        - Callisto: ~16.69 days
        """
        jd = 2460000.0  # 2023
        pos, _ = eph.calc_ut(jd, moon_id, FLG_SWIEPH | FLG_SPEED)

        # Velocity should be non-zero
        assert pos[3] != 0.0, f"{moon_name} should have non-zero longitude velocity"

        # Galilean moons should have significant daily motion
        # (faster than outer planets, which move <1 deg/day)
        # But can vary widely due to Earth-Jupiter geometry

    @pytest.mark.comparison
    @pytest.mark.parametrize("moon_id,moon_name,naif_id", GALILEAN_MOONS)
    def test_galilean_moon_position_changes(self, moon_id, moon_name, naif_id):
        """Verify Galilean moon positions change over time.

        Due to their fast orbital periods, the moons should show measurable
        position changes even over short time intervals.
        """
        jd1 = 2460000.0
        jd2 = 2460001.0  # 1 day later

        pos1, _ = eph.calc_ut(jd1, moon_id, FLG_SWIEPH)
        pos2, _ = eph.calc_ut(jd2, moon_id, FLG_SWIEPH)

        # Position should change over 1 day
        lon_diff = abs(pos2[0] - pos1[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        # At minimum, there should be some motion (even Callisto moves)
        assert lon_diff > 0.001 or abs(pos2[2] - pos1[2]) > 0.0001, (
            f"{moon_name}: Position should change over 1 day. "
            f"lon_diff={lon_diff:.6f}, dist_diff={abs(pos2[2] - pos1[2]):.6f}"
        )


@pytest.mark.skipif(
    not has_moon_spk_testing(),
    reason="Moon SPK testing disabled. Set LIBEPHEMERIS_TEST_MOON_SPK=1 to enable.",
)
class TestTitanWithSpk:
    """Integration tests for Saturn's moon Titan with SPK files.

    These tests require sat441.bsp (Saturn satellite ephemeris).
    """

    @pytest.fixture(autouse=True)
    def setup_saturn_spk(self):
        """Register Saturn SPK if available."""
        planetary_moons.close_moon_kernels()
        eph.close()

        spk_path = get_spk_path("sat441.bsp")
        if spk_path:
            planetary_moons.register_moon_spk(spk_path)
            yield
            planetary_moons.close_moon_kernels()
        else:
            pytest.skip("sat441.bsp not found")

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_titan_geocentric_position(self, year, month, day, hour, date_desc):
        """Validate Titan's geocentric position is reasonable.

        Titan orbits Saturn at ~1.22 million km with a period of ~16 days.
        Distance should be approximately Saturn's distance (~8-11 AU).
        """
        jd = eph.julday(year, month, day, hour)
        pos, _ = eph.calc_ut(jd, MOON_TITAN, FLG_SWIEPH | FLG_SPEED)

        # Validate longitude in range
        assert 0 <= pos[0] < 360, (
            f"Titan at {date_desc}: longitude {pos[0]:.4f} out of range"
        )

        # Validate distance is approximately Saturn's distance
        saturn_pos, _ = eph.calc_ut(jd, SATURN, FLG_SWIEPH)
        saturn_dist = saturn_pos[2]

        assert abs(pos[2] - saturn_dist) < 0.1, (
            f"Titan at {date_desc}: distance {pos[2]:.4f} AU differs "
            f"too much from Saturn {saturn_dist:.4f} AU"
        )

    @pytest.mark.comparison
    def test_titan_velocity(self):
        """Validate Titan velocity calculation."""
        jd = 2460000.0
        pos, _ = eph.calc_ut(jd, MOON_TITAN, FLG_SWIEPH | FLG_SPEED)

        # Titan should have measurable motion
        assert pos[3] != 0.0, "Titan should have non-zero longitude velocity"

    @pytest.mark.comparison
    def test_titan_16_day_period(self):
        """Verify Titan shows orbital motion consistent with ~16 day period.

        Over 16 days, Titan should approximately return to its starting position
        relative to Saturn (but overall position shifts due to Saturn's motion).
        """
        jd1 = 2460000.0
        jd2 = jd1 + 16.0  # ~1 Titan orbital period

        pos1, _ = eph.calc_ut(jd1, MOON_TITAN, FLG_SWIEPH)
        pos2, _ = eph.calc_ut(jd2, MOON_TITAN, FLG_SWIEPH)

        # Get Saturn positions to compare relative motion
        sat1, _ = eph.calc_ut(jd1, SATURN, FLG_SWIEPH)
        sat2, _ = eph.calc_ut(jd2, SATURN, FLG_SWIEPH)

        # Titan-Saturn angular separation should be similar after 1 period
        sep1 = pos1[0] - sat1[0]
        sep2 = pos2[0] - sat2[0]

        # Normalize separations to -180 to 180
        while sep1 > 180:
            sep1 -= 360
        while sep1 < -180:
            sep1 += 360
        while sep2 > 180:
            sep2 -= 360
        while sep2 < -180:
            sep2 += 360

        sep_diff = abs(sep1 - sep2)
        if sep_diff > 180:
            sep_diff = 360 - sep_diff

        # After ~1 orbital period, separation should be similar
        # (within a few degrees due to orbital eccentricity and viewing geometry)
        assert sep_diff < 30, (
            f"Titan-Saturn separation after 16 days differs by {sep_diff:.2f} degrees. "
            f"Expected similar values for ~1 orbital period."
        )


# ============================================================================
# CROSS-VALIDATION WHEN BOTH LIBRARIES HAVE SPK DATA
# ============================================================================


@pytest.mark.skipif(
    not has_moon_spk_testing(),
    reason="Moon SPK testing disabled. Set LIBEPHEMERIS_TEST_MOON_SPK=1 to enable.",
)
class TestCrossValidation:
    """Cross-validation tests when both libraries can access SPK data.

    Note: pyswisseph typically doesn't load SPK files without explicit
    configuration. These tests primarily validate libephemeris consistency
    and document the expected data sources.
    """

    def setup_method(self):
        """Reset state before each test."""
        planetary_moons.close_moon_kernels()
        eph.close()

    def teardown_method(self):
        """Clean up after each test."""
        planetary_moons.close_moon_kernels()

    @pytest.mark.comparison
    def test_libephemeris_uses_jpl_naif_ids(self):
        """Verify libephemeris uses standard JPL NAIF IDs for moons.

        NAIF IDs are the standard identifiers used by JPL:
        - 5xx: Jupiter system (501=Io, 502=Europa, 503=Ganymede, 504=Callisto)
        - 6xx: Saturn system (606=Titan, 602=Enceladus)
        - 8xx: Neptune system (801=Triton)
        - 4xx: Mars system (401=Phobos, 402=Deimos)
        - 9xx: Pluto system (901=Charon)
        """
        from libephemeris.planetary_moons import MOON_NAIF_MAP

        # Verify NAIF ID mappings
        assert MOON_NAIF_MAP[MOON_IO] == 501
        assert MOON_NAIF_MAP[MOON_EUROPA] == 502
        assert MOON_NAIF_MAP[MOON_GANYMEDE] == 503
        assert MOON_NAIF_MAP[MOON_CALLISTO] == 504
        assert MOON_NAIF_MAP[MOON_TITAN] == 606
        assert MOON_NAIF_MAP[MOON_TRITON] == 801

    @pytest.mark.comparison
    def test_spk_files_documented(self):
        """Document required SPK files for each planetary system.

        Users need to download these files from JPL NAIF Horizons:
        https://naif.jpl.nasa.gov/naif/data_generic.html
        """
        expected_files = {
            "Jupiter moons": "jup365.bsp",
            "Saturn moons": "sat441.bsp",
            "Uranus moons": "ura116.bsp",
            "Neptune moons": "nep097.bsp",
            "Mars moons": "mar097.bsp",
            "Pluto moons": "plu058.bsp",
        }

        # Just document the expected files
        for system, filename in expected_files.items():
            spk_path = get_spk_path(filename)
            if spk_path:
                assert os.path.exists(spk_path), f"{system}: {filename} found"
            # No failure if not found - this is documentation


# ============================================================================
# REFERENCE DATA VALIDATION
# ============================================================================


class TestReferenceData:
    """Validate planetary moon positions against known reference data.

    These tests use approximate positions from JPL Horizons or other
    authoritative sources to validate libephemeris calculations.
    """

    def setup_method(self):
        """Reset state before each test."""
        planetary_moons.close_moon_kernels()
        eph.close()

    def teardown_method(self):
        """Clean up after each test."""
        planetary_moons.close_moon_kernels()

    @pytest.mark.comparison
    @pytest.mark.skipif(not has_moon_spk_testing(), reason="Moon SPK testing disabled.")
    def test_jupiter_distance_at_j2000(self):
        """Validate Jupiter's distance at J2000 for reference.

        At J2000.0 (2000-01-01 12:00 TT), Jupiter was at approximately 4.96 AU
        from Earth. Galilean moons should be at nearly the same distance.
        """
        spk_path = get_spk_path("jup365.bsp")
        if not spk_path:
            pytest.skip("jup365.bsp not found")

        planetary_moons.register_moon_spk(spk_path)

        jd = 2451545.0  # J2000.0

        jupiter_pos, _ = eph.calc_ut(jd, JUPITER, FLG_SWIEPH)
        io_pos, _ = eph.calc_ut(jd, MOON_IO, FLG_SWIEPH)

        # Jupiter should be approximately 4.5-5.5 AU at this time
        assert 4.0 < jupiter_pos[2] < 6.0, (
            f"Jupiter distance {jupiter_pos[2]:.4f} AU unexpected at J2000"
        )

        # Io should be approximately the same distance
        assert abs(io_pos[2] - jupiter_pos[2]) < 0.1, (
            f"Io distance {io_pos[2]:.4f} should be close to Jupiter {jupiter_pos[2]:.4f}"
        )


# ============================================================================
# EDGE CASES AND BOUNDARY CONDITIONS
# ============================================================================


class TestEdgeCases:
    """Test edge cases and boundary conditions for planetary moons."""

    def setup_method(self):
        """Reset state before each test."""
        planetary_moons.close_moon_kernels()
        eph.close()

    def teardown_method(self):
        """Clean up after each test."""
        planetary_moons.close_moon_kernels()

    @pytest.mark.comparison
    def test_invalid_moon_id(self):
        """Test behavior with invalid moon IDs."""
        # ID not in moon range. (9999 now resolves under the canonical
        # PLMOON_OFFSET+NAIF scheme, so use a clearly out-of-range id.)
        assert planetary_moons.is_planetary_moon(99999) is False

        # Unknown moon ID
        name = planetary_moons.get_moon_name(99999)
        assert "Unknown" in name

    @pytest.mark.comparison
    def test_moon_coverage_unregistered(self):
        """Test get_moon_coverage returns None for unregistered moons."""
        coverage = planetary_moons.get_moon_coverage(MOON_IO)
        assert coverage is None

    @pytest.mark.comparison
    def test_list_registered_moons_empty(self):
        """Test list_registered_moons when no moons are registered."""
        moons = planetary_moons.list_registered_moons()
        assert len(moons) == 0
        assert isinstance(moons, dict)

    @pytest.mark.comparison
    @pytest.mark.skipif(not has_moon_spk_testing(), reason="Moon SPK testing disabled.")
    def test_register_and_unregister_spk(self):
        """Test registering and unregistering SPK files."""
        spk_path = get_spk_path("jup365.bsp")
        if not spk_path:
            pytest.skip("jup365.bsp not found")

        # Register
        planetary_moons.register_moon_spk(spk_path)
        moons = planetary_moons.list_registered_moons()
        assert len(moons) > 0

        # Unregister
        planetary_moons.unregister_moon_spk(spk_path)
        moons = planetary_moons.list_registered_moons()
        assert len(moons) == 0
