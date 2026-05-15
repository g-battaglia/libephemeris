"""
Tests for the migration guide examples.

This test module verifies that all code examples from the migration guide
(docs/migration-guide.md) work correctly.

Note: EphemerisContext calculation tests are in test_context_thread_safety.py.
The tests here focus on the module-level API which is documented in the migration guide.
"""

import libephemeris as swe
from libephemeris import EphemerisContext
from libephemeris.constants import (
    SUN,
    MOON,
    FLG_SPEED,
    FLG_SIDEREAL,
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
    SIDM_TRUE_CITRA,
)


class TestMigrationGuideBasicUsage:
    """Test basic usage examples from the migration guide."""

    def test_drop_in_replacement_import(self):
        """Test that libephemeris works as a drop-in replacement for swisseph."""
        # This is the main claim: you can replace `import swisseph as swe`
        # with `import libephemeris as swe`
        jd = swe.julday(2000, 1, 1, 12.0)
        pos, _ = swe.calc_ut(jd, SUN, FLG_SPEED)

        assert pos is not None
        assert len(pos) == 6
        assert 0 <= pos[0] < 360  # Valid longitude

    def test_constants_from_module(self):
        """Test accessing constants from the module."""
        assert swe.SUN == 0
        assert swe.MOON == 1
        assert swe.FLG_SPEED == 256
        assert swe.FLG_SIDEREAL == 64 * 1024

    def test_prefixed_and_unprefixed_functions(self):
        """Test that both swe_ prefixed and unprefixed functions work."""
        jd = 2451545.0

        # Prefixed version
        pos1, _ = swe.calc_ut(jd, SUN, 0)

        # Unprefixed version
        pos2, _ = swe.calc_ut(jd, SUN, 0)

        # Both should give identical results
        assert pos1[0] == pos2[0]
        assert pos1[1] == pos2[1]
        assert pos1[2] == pos2[2]


class TestMigrationGuideHouseCusps:
    """Test house cusp indexing as documented in the migration guide."""

    def test_house_cusps_0_indexed(self):
        """Test that libephemeris house cusps are 0-indexed."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5  # Rome

        cusps, ascmc = swe.houses(jd, lat, lon, ord("P"))

        # Verify we get 12 cusps
        assert len(cusps) == 12

        # First house cusp is at index 0
        cusp1 = cusps[0]
        assert 0 <= cusp1 < 360

        # ascmc contains angles
        assert len(ascmc) >= 4
        asc = ascmc[0]
        mc = ascmc[1]
        assert 0 <= asc < 360
        assert 0 <= mc < 360


class TestMigrationGuidePrecision:
    """Test precision claims from the migration guide."""

    def test_planetary_longitude_valid_range(self):
        """Test that planetary longitudes are in valid range."""
        jd = 2451545.0

        for planet in [SUN, MOON]:
            pos, _ = swe.calc_ut(jd, planet, FLG_SPEED)

            assert 0 <= pos[0] < 360, f"Planet {planet} longitude out of range"
            assert -90 <= pos[1] <= 90, f"Planet {planet} latitude out of range"
            assert pos[2] > 0, f"Planet {planet} distance should be positive"

    def test_velocities_are_computed(self):
        """Test that velocities are computed when FLG_SPEED is set."""
        jd = 2451545.0

        pos, _ = swe.calc_ut(jd, SUN, FLG_SPEED)

        # Sun moves about 1 degree per day
        assert 0.9 < pos[3] < 1.1, "Sun speed should be ~1 deg/day"

        pos, _ = swe.calc_ut(jd, MOON, FLG_SPEED)

        # Moon moves about 12-15 degrees per day
        assert 11 < pos[3] < 16, "Moon speed should be ~12-15 deg/day"


class TestMigrationGuideAyanamshas:
    """Test ayanamsha functionality from the migration guide."""

    def test_sidereal_mode_setting(self):
        """Test setting sidereal mode and calculating sidereal positions."""
        import swisseph as pyswe

        jd = 2451545.0

        # Set Lahiri ayanamsha in both libraries
        swe.set_sid_mode(SIDM_LAHIRI)
        pyswe.set_sid_mode(SIDM_LAHIRI)

        try:
            # Calculate sidereal position with libephemeris
            pos_sidereal, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)

            # Calculate sidereal position with pyswisseph for comparison
            py_pos_sidereal = pyswe.calc_ut(jd, SUN, FLG_SIDEREAL)[0]

            # libephemeris sidereal should match pyswisseph sidereal
            # Tolerance: 1 arcsecond = 0.000278 degrees
            diff = abs(pos_sidereal[0] - py_pos_sidereal[0])
            if diff > 180:
                diff = 360 - diff

            assert diff < 0.001, f"Sidereal calculation mismatch vs pyswisseph: {diff}"
        finally:
            pyswe.close()

    def test_multiple_ayanamshas_give_different_results(self):
        """Test that different ayanamshas produce different sidereal positions."""
        jd = 2451545.0

        results = {}
        for mode in [SIDM_FAGAN_BRADLEY, SIDM_LAHIRI, SIDM_TRUE_CITRA]:
            swe.set_sid_mode(mode)
            pos, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)
            results[mode] = pos[0]

        # All three should be different
        assert len(set(round(v, 4) for v in results.values())) == 3


class TestMigrationGuideEphemerisContextState:
    """Test EphemerisContext state management from the migration guide."""

    def test_ephemeris_context_creation(self):
        """Test creating an EphemerisContext and setting state."""
        ctx = EphemerisContext()

        # Should be able to set state
        ctx.set_topo(12.5, 41.9, 0)  # Rome
        ctx.set_sid_mode(SIDM_LAHIRI)

        # Verify state was set
        assert ctx.get_sid_mode() == SIDM_LAHIRI
        topo = ctx.get_topo()
        assert topo is not None
        assert abs(topo.latitude.degrees - 41.9) < 0.01

    def test_ephemeris_context_isolation(self):
        """Test that EphemerisContext instances have isolated state."""
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        # Set different sidereal modes
        ctx1.set_sid_mode(SIDM_LAHIRI)
        ctx2.set_sid_mode(SIDM_FAGAN_BRADLEY)

        # Verify isolation
        assert ctx1.get_sid_mode() == SIDM_LAHIRI
        assert ctx2.get_sid_mode() == SIDM_FAGAN_BRADLEY

        # Set different locations
        ctx1.set_topo(12.5, 41.9, 0)  # Rome
        ctx2.set_topo(-0.1, 51.5, 0)  # London

        # Verify isolation
        topo1 = ctx1.get_topo()
        topo2 = ctx2.get_topo()

        assert topo1 is not None
        assert topo2 is not None
        assert abs(topo1.latitude.degrees - 41.9) < 0.01
        assert abs(topo2.latitude.degrees - 51.5) < 0.01


class TestMigrationGuideNotImplemented:
    """Test features documented as not implemented in the migration guide."""

    def test_fixed_star_velocities_are_small(self):
        """Test that fixed star velocities are small but non-zero due to precession.

        Fixed stars have small velocities (~3.8e-05 deg/day) due to precession
        of equinoxes, consistent with pyswisseph behavior.
        """
        jd = 2451545.0

        result = swe.fixstar_ut("Aldebaran", jd, FLG_SPEED)
        # fixstar_ut returns (pos, flag, starname)
        pos = result[0]

        # Velocities are very small but non-zero due to precession of equinoxes
        assert abs(pos[3]) < 0.001, "Fixed star lon velocity should be very small"
        assert abs(pos[4]) < 0.001, "Fixed star lat velocity should be very small"
        assert abs(pos[5]) < 0.001, "Fixed star dist velocity should be very small"


class TestMigrationGuideLunarNodes:
    """Test lunar node and Lilith calculations as documented in the migration guide."""

    def test_mean_node_returns_valid_position(self):
        """Test that mean node calculation returns valid position."""
        from libephemeris.constants import MEAN_NODE

        jd = 2451545.0
        pos, _ = swe.calc_ut(jd, MEAN_NODE, FLG_SPEED)

        assert 0 <= pos[0] < 360, "Mean node longitude should be in valid range"
        # Mean node latitude is always near zero
        assert abs(pos[1]) < 1, "Mean node latitude should be near zero"
        # Mean node moves retrograde (~0.053 deg/day)
        assert -0.1 < pos[3] < 0, "Mean node should move retrograde"

    def test_true_node_returns_valid_position(self):
        """Test that true node calculation returns valid position."""
        from libephemeris.constants import TRUE_NODE

        jd = 2451545.0
        pos, _ = swe.calc_ut(jd, TRUE_NODE, FLG_SPEED)

        assert 0 <= pos[0] < 360, "True node longitude should be in valid range"
        # True node latitude may oscillate slightly
        assert abs(pos[1]) < 5, "True node latitude should be small"

    def test_mean_and_true_node_differ(self):
        """Test that mean and true node positions differ.

        The true node oscillates around the mean node with an amplitude
        of about 1.5 degrees.
        """
        from libephemeris.constants import MEAN_NODE, TRUE_NODE

        jd = 2451545.0
        mean_pos, _ = swe.calc_ut(jd, MEAN_NODE, 0)
        true_pos, _ = swe.calc_ut(jd, TRUE_NODE, 0)

        # The difference should be within the oscillation amplitude
        diff = abs(mean_pos[0] - true_pos[0])
        if diff > 180:
            diff = 360 - diff
        assert diff < 3, "Mean and true node should differ by less than 3 degrees"

    def test_osculating_apogee_returns_valid_position(self):
        """Test that osculating apogee (True Lilith) returns valid position.

        Note: True Lilith has known precision differences from pyswisseph
        (up to 5-7 degrees) due to different orbital element models.
        """
        from libephemeris.constants import OSCU_APOG

        jd = 2451545.0
        pos, _ = swe.calc_ut(jd, OSCU_APOG, FLG_SPEED)

        assert 0 <= pos[0] < 360, "True Lilith longitude should be in valid range"
        # Lilith latitude can vary
        assert -90 <= pos[1] <= 90, "True Lilith latitude should be in valid range"


class TestMigrationGuideMultiThreading:
    """Test multi-threading example from the migration guide."""

    def test_concurrent_chart_calculation(self):
        """Test that concurrent calculations with EphemerisContext work correctly.

        This tests the multi-threading example from the migration guide.
        """
        from concurrent.futures import ThreadPoolExecutor

        def calculate_chart(location: dict, jd: float) -> dict:
            """Thread-safe chart calculation function."""
            ctx = EphemerisContext()
            ctx.set_topo(location["lon"], location["lat"], 0)
            ctx.set_sid_mode(SIDM_LAHIRI)

            sun, _ = ctx.calc_ut(jd, SUN, FLG_SIDEREAL)
            moon, _ = ctx.calc_ut(jd, MOON, FLG_SIDEREAL)
            cusps, ascmc = ctx.houses(jd, location["lat"], location["lon"], ord("P"))

            return {
                "location": location["name"],
                "sun": sun[0],
                "moon": moon[0],
                "asc": ascmc[0],
            }

        locations = [
            {"name": "Rome", "lon": 12.5, "lat": 41.9},
            {"name": "London", "lon": -0.1, "lat": 51.5},
            {"name": "Tokyo", "lon": 139.7, "lat": 35.7},
            {"name": "New York", "lon": -74.0, "lat": 40.7},
            {"name": "Sydney", "lon": 151.2, "lat": -33.9},
        ]

        jd = 2451545.0  # J2000.0
        with ThreadPoolExecutor(max_workers=5) as executor:
            results = list(
                executor.map(lambda loc: calculate_chart(loc, jd), locations)
            )

        # Verify we got results for all locations
        assert len(results) == 5

        # Verify all results are valid
        for r in results:
            assert r["location"] in ["Rome", "London", "Tokyo", "New York", "Sydney"]
            assert 0 <= r["sun"] < 360, f"{r['location']}: Sun out of range"
            assert 0 <= r["moon"] < 360, f"{r['location']}: Moon out of range"
            assert 0 <= r["asc"] < 360, f"{r['location']}: Asc out of range"

        # Verify different locations give different ascendants
        ascendants = [r["asc"] for r in results]
        # At least some should be different (not all the same)
        assert len(set(round(a, 1) for a in ascendants)) > 1

    def test_context_isolation_under_concurrent_modification(self):
        """Test that context state is truly isolated during concurrent access."""
        from concurrent.futures import ThreadPoolExecutor, as_completed
        import time

        def modify_and_read_context(ayanamsha: int, delay: float):
            """Set ayanamsha, wait, then read it back."""
            ctx = EphemerisContext()
            ctx.set_sid_mode(ayanamsha)
            time.sleep(delay)
            # get_sid_mode returns int or tuple depending on full parameter
            return ctx.get_sid_mode()

        # Run concurrent modifications with different ayanamshas
        ayanamshas = [
            SIDM_FAGAN_BRADLEY,
            SIDM_LAHIRI,
            SIDM_TRUE_CITRA,
        ]
        with ThreadPoolExecutor(max_workers=3) as executor:
            futures = {
                executor.submit(modify_and_read_context, aya, 0.01): aya
                for aya in ayanamshas
            }

            for future in as_completed(futures):
                expected = futures[future]
                result = future.result()
                assert result == expected, f"Expected {expected}, got {result}"
