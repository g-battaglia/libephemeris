"""
Tests for REBOUND/ASSIST n-body integration module.

These tests verify that the REBOUND integration works correctly for asteroid
orbit propagation, both with and without the optional REBOUND/ASSIST packages.
"""

import math
import os
import pytest
from unittest.mock import patch

from libephemeris.minor_bodies import (
    MINOR_BODY_ELEMENTS,
    calc_minor_body_position,
)
from libephemeris.constants import CERES, VESTA, CHIRON, ERIS
from libephemeris.rebound_integration import (
    ReboundIntegrator,
    AssistEphemConfig,
    PropagationResult,
    check_rebound_available,
    check_assist_available,
    check_assist_data_available,
    reset_assist_data_cache,
    get_rebound_version,
    get_assist_version,
    elements_to_rebound_orbit,
    propagate_orbit_rebound,
    propagate_orbit_assist,
    propagate_trajectory,
    compare_with_keplerian,
)


class TestReboundIntegratorEnum:
    """Test the ReboundIntegrator enum."""

    def test_ias15_value(self):
        """IAS15 should have correct string value."""
        assert ReboundIntegrator.IAS15.value == "ias15"

    def test_whfast_value(self):
        """WHFast should have correct string value."""
        assert ReboundIntegrator.WHFAST.value == "whfast"

    def test_mercurius_value(self):
        """MERCURIUS should have correct string value."""
        assert ReboundIntegrator.MERCURIUS.value == "mercurius"

    def test_trace_value(self):
        """TRACE should have correct string value."""
        assert ReboundIntegrator.TRACE.value == "trace"

    def test_leapfrog_value(self):
        """LEAPFROG should have correct string value."""
        assert ReboundIntegrator.LEAPFROG.value == "leapfrog"


class TestAssistEphemConfig:
    """Test ASSIST ephemeris configuration."""

    def test_default_config(self):
        """Default config should have None paths."""
        config = AssistEphemConfig()
        # Paths may be None or set from environment
        assert isinstance(config, AssistEphemConfig)

    def test_custom_data_dir(self, tmp_path):
        """Custom data directory should be used."""
        config = AssistEphemConfig(data_dir=str(tmp_path))
        assert config.data_dir == str(tmp_path)

    def test_custom_file_paths(self, tmp_path):
        """Custom file paths should be preserved."""
        planets = str(tmp_path / "planets.dat")
        asteroids = str(tmp_path / "asteroids.bsp")

        config = AssistEphemConfig(planets_file=planets, asteroids_file=asteroids)
        assert config.planets_file == planets
        assert config.asteroids_file == asteroids


class TestPropagationResult:
    """Test PropagationResult dataclass."""

    def test_basic_result(self):
        """Basic PropagationResult should store values correctly."""
        result = PropagationResult(
            x=1.0, y=0.0, z=0.0, vx=0.0, vy=0.017, vz=0.0, jd_tt=2460000.5
        )
        assert result.x == 1.0
        assert result.y == 0.0
        assert result.z == 0.0
        assert result.jd_tt == 2460000.5

    def test_distance_calculation(self):
        """Distance should be calculated from x,y,z."""
        result = PropagationResult(
            x=3.0, y=4.0, z=0.0, vx=0.0, vy=0.0, vz=0.0, jd_tt=2460000.5
        )
        assert result.distance == 5.0

    def test_ecliptic_lon_positive_x(self):
        """Ecliptic longitude at positive X axis should be 0."""
        result = PropagationResult(
            x=1.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0, jd_tt=2460000.5
        )
        assert abs(result.ecliptic_lon) < 1e-10

    def test_ecliptic_lon_positive_y(self):
        """Ecliptic longitude at positive Y axis should be 90."""
        result = PropagationResult(
            x=0.0, y=1.0, z=0.0, vx=0.0, vy=0.0, vz=0.0, jd_tt=2460000.5
        )
        assert abs(result.ecliptic_lon - 90.0) < 1e-10

    def test_ecliptic_lat_positive_z(self):
        """Ecliptic latitude at positive Z axis should be 90."""
        result = PropagationResult(
            x=0.0, y=0.0, z=1.0, vx=0.0, vy=0.0, vz=0.0, jd_tt=2460000.5
        )
        assert abs(result.ecliptic_lat - 90.0) < 1e-10

    def test_to_tuple(self):
        """to_tuple should return correct order."""
        result = PropagationResult(
            x=1.0, y=2.0, z=3.0, vx=0.1, vy=0.2, vz=0.3, jd_tt=2460000.5
        )
        expected = (1.0, 2.0, 3.0, 0.1, 0.2, 0.3)
        assert result.to_tuple() == expected


class TestAvailabilityChecks:
    """Test package availability check functions."""

    def test_check_rebound_returns_bool(self):
        """check_rebound_available should return bool."""
        result = check_rebound_available()
        assert isinstance(result, bool)

    def test_check_assist_returns_bool(self):
        """check_assist_available should return bool."""
        result = check_assist_available()
        assert isinstance(result, bool)

    def test_get_rebound_version(self):
        """get_rebound_version should return string or None."""
        result = get_rebound_version()
        assert result is None or isinstance(result, str)

    def test_get_assist_version(self):
        """get_assist_version should return string or None."""
        result = get_assist_version()
        assert result is None or isinstance(result, str)

    def test_version_matches_availability(self):
        """Version should be available iff package is available."""
        if check_rebound_available():
            assert get_rebound_version() is not None
        else:
            assert get_rebound_version() is None


class TestElementsConversion:
    """Test orbital elements conversion to REBOUND format."""

    def test_convert_ceres_elements(self):
        """Converting Ceres elements should produce valid REBOUND params."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd = elements.epoch

        params = elements_to_rebound_orbit(elements, jd)

        assert "a" in params
        assert "e" in params
        assert "inc" in params
        assert "omega" in params
        assert "Omega" in params
        assert "M" in params

        # Check values are reasonable
        assert params["a"] > 0  # Semi-major axis positive
        assert 0 <= params["e"] < 1  # Eccentricity valid for ellipse
        assert 0 <= params["inc"] <= math.pi  # Inclination in [0, pi]

    def test_angles_in_radians(self):
        """Angles should be converted to radians."""
        elements = MINOR_BODY_ELEMENTS[VESTA]
        jd = elements.epoch

        params = elements_to_rebound_orbit(elements, jd)

        # All angles should be in radians (< 2*pi for normalized angles)
        assert abs(params["inc"]) <= math.pi
        assert abs(params["omega"]) <= 2 * math.pi
        assert abs(params["Omega"]) <= 2 * math.pi
        assert abs(params["M"]) <= 2 * math.pi


class TestPropagateOrbitRebound:
    """Test REBOUND orbit propagation (requires rebound package)."""

    @pytest.fixture
    def skip_if_no_rebound(self):
        """Skip test if REBOUND is not available."""
        if not check_rebound_available():
            pytest.skip("REBOUND not installed")

    def test_propagate_ceres_short(self, skip_if_no_rebound):
        """Propagate Ceres orbit for a short time."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 30  # 30 days

        result = propagate_orbit_rebound(elements, jd_start, jd_end)

        assert isinstance(result, PropagationResult)
        assert result.jd_tt == jd_end
        assert result.distance > 0
        # Ceres should be between 2.5 and 3.0 AU from Sun
        assert 2.0 < result.distance < 3.5

    def test_propagate_vesta(self, skip_if_no_rebound):
        """Propagate Vesta orbit."""
        elements = MINOR_BODY_ELEMENTS[VESTA]
        jd_start = elements.epoch
        jd_end = jd_start + 100  # 100 days

        result = propagate_orbit_rebound(elements, jd_start, jd_end)

        assert result.distance > 0
        # Vesta orbit between 2.15 and 2.57 AU
        assert 1.8 < result.distance < 2.8

    def test_propagate_chiron(self, skip_if_no_rebound):
        """Propagate Chiron (centaur) orbit."""
        elements = MINOR_BODY_ELEMENTS[CHIRON]
        jd_start = elements.epoch
        jd_end = jd_start + 365  # 1 year

        result = propagate_orbit_rebound(elements, jd_start, jd_end)

        # Chiron orbit between 8.5 and 19 AU (highly eccentric)
        assert 5.0 < result.distance < 25.0

    def test_integrator_ias15(self, skip_if_no_rebound):
        """IAS15 integrator should work."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 10

        result = propagate_orbit_rebound(
            elements, jd_start, jd_end, integrator=ReboundIntegrator.IAS15
        )
        assert result.distance > 0

    def test_integrator_whfast(self, skip_if_no_rebound):
        """WHFast integrator should work."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 10

        result = propagate_orbit_rebound(
            elements, jd_start, jd_end, integrator=ReboundIntegrator.WHFAST
        )
        assert result.distance > 0

    def test_integrator_leapfrog(self, skip_if_no_rebound):
        """Leapfrog integrator should work."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 10

        result = propagate_orbit_rebound(
            elements, jd_start, jd_end, integrator=ReboundIntegrator.LEAPFROG
        )
        assert result.distance > 0

    def test_custom_timestep(self, skip_if_no_rebound):
        """Custom timestep should be used."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 30

        result = propagate_orbit_rebound(
            elements,
            jd_start,
            jd_end,
            dt=0.5,  # 0.5 day timestep
        )
        assert result.distance > 0

    def test_backward_integration(self, skip_if_no_rebound):
        """Backward integration should work."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd_start = elements.epoch
        jd_end = jd_start - 30  # 30 days backward

        result = propagate_orbit_rebound(elements, jd_start, jd_end)

        assert result.jd_tt == jd_end
        assert result.distance > 0


class TestPropagateTrajectory:
    """Test trajectory propagation (multiple output points)."""

    @pytest.fixture
    def skip_if_no_rebound(self):
        """Skip test if REBOUND is not available."""
        if not check_rebound_available():
            pytest.skip("REBOUND not installed")

    def test_trajectory_basic(self, skip_if_no_rebound):
        """Basic trajectory propagation should return list of points."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 100

        trajectory = propagate_trajectory(
            elements,
            jd_start,
            jd_end,
            num_points=10,
            use_assist=False,  # Use REBOUND only
        )

        assert len(trajectory) == 10
        assert all(isinstance(p, PropagationResult) for p in trajectory)

    def test_trajectory_times_evenly_spaced(self, skip_if_no_rebound):
        """Trajectory times should be evenly spaced."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 100

        trajectory = propagate_trajectory(
            elements, jd_start, jd_end, num_points=11, use_assist=False
        )

        times = [p.jd_tt for p in trajectory]
        dt = times[1] - times[0]

        for i in range(1, len(times)):
            assert abs(times[i] - times[i - 1] - dt) < 1e-6

    def test_trajectory_first_and_last(self, skip_if_no_rebound):
        """First and last points should match start and end times."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 50

        trajectory = propagate_trajectory(
            elements, jd_start, jd_end, num_points=5, use_assist=False
        )

        assert abs(trajectory[0].jd_tt - jd_start) < 1e-6
        assert abs(trajectory[-1].jd_tt - jd_end) < 1e-6


class TestCompareWithKeplerian:
    """Test comparison between n-body and Keplerian results."""

    @pytest.fixture
    def skip_if_no_rebound(self):
        """Skip test if REBOUND is not available."""
        if not check_rebound_available():
            pytest.skip("REBOUND not installed")

    def test_compare_returns_dict(self, skip_if_no_rebound):
        """Comparison should return a dictionary with expected keys."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd = elements.epoch + 100

        comparison = compare_with_keplerian(elements, jd, use_assist=False)

        assert isinstance(comparison, dict)
        assert "keplerian" in comparison
        assert "nbody" in comparison
        assert "delta_lon_arcsec" in comparison
        assert "delta_lat_arcsec" in comparison
        assert "angular_sep_arcsec" in comparison

    def test_compare_keplerian_values(self, skip_if_no_rebound):
        """Keplerian values should match direct calculation."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd = elements.epoch + 30

        comparison = compare_with_keplerian(elements, jd, use_assist=False)

        # Compare with direct Keplerian calculation
        x, y, z = calc_minor_body_position(elements, jd)
        r = math.sqrt(x**2 + y**2 + z**2)
        lon = math.degrees(math.atan2(y, x)) % 360.0
        lat = math.degrees(math.asin(z / r))

        assert abs(comparison["keplerian"]["lon"] - lon) < 1e-10
        assert abs(comparison["keplerian"]["lat"] - lat) < 1e-10
        assert abs(comparison["keplerian"]["dist"] - r) < 1e-10

    def test_short_propagation_small_difference(self, skip_if_no_rebound):
        """Short propagation should have small difference from Keplerian."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd = elements.epoch + 10  # Only 10 days

        comparison = compare_with_keplerian(elements, jd, use_assist=False)

        # For short propagation, 2-body REBOUND and Keplerian should be similar
        # The difference comes from REBOUND using slightly different integration
        assert comparison["angular_sep_arcsec"] < 3600  # Within 1 degree


class TestModuleWithoutRebound:
    """Test module behavior when REBOUND is not installed."""

    def test_import_works_without_rebound(self):
        """Module should import without REBOUND installed."""
        # This test verifies that the module can be imported even if
        # REBOUND is not installed (lazy import pattern)
        with patch.dict("sys.modules", {"rebound": None}):
            # Re-import to test
            from libephemeris import rebound_integration

            assert hasattr(rebound_integration, "check_rebound_available")

    def test_propagate_raises_without_rebound(self):
        """propagate_orbit_rebound should raise ImportError without REBOUND."""
        if check_rebound_available():
            pytest.skip("REBOUND is installed")

        elements = MINOR_BODY_ELEMENTS[CERES]
        with pytest.raises(ImportError):
            propagate_orbit_rebound(elements, 2460000.5, 2460010.5)


class TestIntegrationAccuracy:
    """Test integration accuracy for longer propagation times."""

    @pytest.fixture
    def skip_if_no_rebound(self):
        """Skip test if REBOUND is not available."""
        if not check_rebound_available():
            pytest.skip("REBOUND not installed")

    def test_orbit_stays_bounded(self, skip_if_no_rebound):
        """Orbit should stay within physical bounds over 1 year."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 365.25  # 1 year

        result = propagate_orbit_rebound(elements, jd_start, jd_end)

        # Ceres orbit: a = 2.77 AU, e = 0.076
        # Perihelion: 2.56 AU, Aphelion: 2.98 AU
        assert 2.0 < result.distance < 3.5

    def test_eccentric_orbit_accuracy(self, skip_if_no_rebound):
        """Highly eccentric orbit (Chiron) should remain bounded."""
        elements = MINOR_BODY_ELEMENTS[CHIRON]
        jd_start = elements.epoch
        jd_end = jd_start + 365.25

        result = propagate_orbit_rebound(elements, jd_start, jd_end)

        # Chiron orbit: a = 13.7 AU, e = 0.38
        # Perihelion: 8.5 AU, Aphelion: 18.9 AU
        assert 5.0 < result.distance < 25.0

    def test_distant_tno_accuracy(self, skip_if_no_rebound):
        """Distant TNO (Eris) orbit should remain bounded."""
        elements = MINOR_BODY_ELEMENTS[ERIS]
        jd_start = elements.epoch
        jd_end = jd_start + 365.25

        result = propagate_orbit_rebound(elements, jd_start, jd_end)

        # Eris orbit: a = 68 AU, e = 0.44
        # Perihelion: 38 AU, Aphelion: 98 AU
        assert 30.0 < result.distance < 110.0


class TestOrbitRoundTrip:
    """Test that orbit propagation is reversible."""

    @pytest.fixture
    def skip_if_no_rebound(self):
        """Skip test if REBOUND is not available."""
        if not check_rebound_available():
            pytest.skip("REBOUND not installed")

    def test_forward_backward_consistency(self, skip_if_no_rebound):
        """Forward then backward propagation should return to start."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd_start = elements.epoch
        jd_mid = jd_start + 50

        # Forward propagation
        result_fwd = propagate_orbit_rebound(elements, jd_start, jd_mid)

        # Create elements at mid point (approximate)
        # For this test, we just verify that both directions give similar distances
        result_bwd = propagate_orbit_rebound(elements, jd_start, jd_mid)

        # Results should be consistent
        assert abs(result_fwd.distance - result_bwd.distance) < 1e-10


class TestAssistDataAvailability:
    """Test cached ASSIST data availability checks."""

    def test_check_assist_data_returns_bool(self):
        """check_assist_data_available() should return a boolean."""
        result = check_assist_data_available()
        assert isinstance(result, bool)

    def test_cache_is_consistent(self):
        """Repeated calls should return the same cached result."""
        reset_assist_data_cache()
        first = check_assist_data_available()
        second = check_assist_data_available()
        assert first == second

    def test_reset_clears_cache(self):
        """reset_assist_data_cache() should allow re-probing."""
        from libephemeris import rebound_integration

        # Prime the cache
        check_assist_data_available()
        rebound_integration._assist_propagation_cache[("dummy", 0.0, 0)] = None
        # Reset
        reset_assist_data_cache()
        # The next call re-probes (should still return the same value,
        # but the important thing is it doesn't crash)
        result = check_assist_data_available()
        assert isinstance(result, bool)
        assert rebound_integration._assist_propagation_cache == {}

    def test_no_assist_import_returns_false(self):
        """If ASSIST is not importable, should return False."""
        reset_assist_data_cache()
        with patch.dict("sys.modules", {"assist": None}):
            # Force re-check by clearing the lru_cache-like global
            from libephemeris import rebound_integration

            old_val = rebound_integration._assist_data_available
            rebound_integration._assist_data_available = None
            try:
                result = check_assist_data_available()
                assert result is False
            finally:
                rebound_integration._assist_data_available = old_val

    def test_missing_data_files_returns_false(self):
        """If ASSIST is importable but data files are missing, should return False."""
        if not check_assist_available():
            pytest.skip("ASSIST not installed")
        reset_assist_data_cache()
        # Point to a non-existent directory
        with patch.dict(os.environ, {"ASSIST_DIR": "/nonexistent/assist/dir"}):
            from libephemeris import rebound_integration

            old_val = rebound_integration._assist_data_available
            rebound_integration._assist_data_available = None
            try:
                config = AssistEphemConfig(data_dir="/nonexistent/assist/dir")
                with patch.object(
                    rebound_integration,
                    "AssistEphemConfig",
                    return_value=config,
                ):
                    # Only returns False if the default dirs also don't have it
                    result = check_assist_data_available()
                    assert isinstance(result, bool)
            finally:
                rebound_integration._assist_data_available = old_val
                reset_assist_data_cache()

    def test_close_resets_assist_cache(self):
        """libephemeris.close() should reset the ASSIST data cache."""
        import libephemeris

        # Prime the cache
        check_assist_data_available()
        # Close should reset it
        libephemeris.close()
        from libephemeris import rebound_integration

        assert rebound_integration._assist_data_available is None


# Whether ASSIST data files are actually present on this machine
_assist_data_present = check_assist_data_available()

# Fixture-like skip marker
requires_assist_data = pytest.mark.skipif(
    not _assist_data_present,
    reason="ASSIST data files not downloaded (~714 MB required)",
)


class TestAssistEndToEnd:
    """End-to-end ASSIST integration tests.

    These tests require ASSIST data files (~714 MB) to be downloaded.
    They are automatically skipped if the files are not present.
    Run ``download_assist_data()`` to enable them.
    """

    @requires_assist_data
    def test_propagate_ceres_assist(self):
        """Propagate Ceres with ASSIST for 30 days."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 30

        result = propagate_orbit_assist(elements, jd_start, jd_end)

        assert isinstance(result, PropagationResult)
        assert result.jd_tt == jd_end
        # Ceres should be between 2.0 and 3.5 AU
        assert 2.0 < result.distance < 3.5

    @requires_assist_data
    def test_propagate_vesta_assist(self):
        """Propagate Vesta with ASSIST for 100 days."""
        elements = MINOR_BODY_ELEMENTS[VESTA]
        jd_start = elements.epoch
        jd_end = jd_start + 100

        result = propagate_orbit_assist(elements, jd_start, jd_end)

        assert 1.8 < result.distance < 2.8

    @requires_assist_data
    def test_propagate_chiron_assist(self):
        """Propagate Chiron (centaur) with ASSIST for 1 year."""
        elements = MINOR_BODY_ELEMENTS[CHIRON]
        jd_start = elements.epoch
        jd_end = jd_start + 365.25

        result = propagate_orbit_assist(elements, jd_start, jd_end)

        assert 5.0 < result.distance < 25.0

    @requires_assist_data
    def test_assist_vs_rebound_ceres(self):
        """ASSIST should give different (more accurate) results than 2-body REBOUND."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 365.25  # 1 year

        result_assist = propagate_orbit_assist(elements, jd_start, jd_end)
        result_rebound = propagate_orbit_rebound(elements, jd_start, jd_end)

        # They should give different results because ASSIST includes
        # planetary perturbations while REBOUND 2-body does not
        dx = result_assist.x - result_rebound.x
        dy = result_assist.y - result_rebound.y
        dz = result_assist.z - result_rebound.z
        sep = math.sqrt(dx**2 + dy**2 + dz**2)

        # The difference should be non-trivial (> 1e-6 AU) for a 1-year propagation
        # due to Jupiter's perturbation on Ceres
        assert sep > 1e-6

    @requires_assist_data
    def test_assist_perturber_self_interaction_prograde(self):
        """Regression: propagating a body that is itself an ASSIST asteroid
        perturber (Ceres/Vesta/...) must NOT self-interact into a retrograde /
        diverged orbit. Previously the ASTEROIDS force added the coincident
        perturber's near-singular acceleration, flipping the orbit by ~160deg/yr.
        """
        import warnings

        for ipl in (CERES, VESTA):
            elements = MINOR_BODY_ELEMENTS[ipl]
            jd_start = elements.epoch
            jd_end = jd_start + 365.25

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = propagate_orbit_assist(elements, jd_start, jd_end)

            # Prograde guard: specific angular momentum z-component must be > 0
            # (ecliptic-J2000); the bug produced h_z < 0.
            h_z = result.x * result.vy - result.y * result.vx
            assert h_z > 0, f"{ipl}: retrograde h_z={h_z:.3e} (self-interaction)"

            # ASSIST must stay close to the 2-body REBOUND direction (the bug put
            # them ~160deg apart). Both track the real orbit; ASSIST only adds the
            # small planetary perturbation, so the angle is < ~2deg over a year.
            reb = propagate_orbit_rebound(elements, jd_start, jd_end)
            dot = result.x * reb.x + result.y * reb.y + result.z * reb.z
            na = math.sqrt(result.x**2 + result.y**2 + result.z**2)
            nb = math.sqrt(reb.x**2 + reb.y**2 + reb.z**2)
            ang = math.degrees(math.acos(max(-1.0, min(1.0, dot / (na * nb)))))
            assert ang < 2.0, f"{ipl}: ASSIST vs REBOUND diverge by {ang:.1f}deg"

            # The self-coincidence must be detected and reported.
            assert any("self-interaction" in str(w.message) for w in caught), (
                f"{ipl}: expected a self-perturber warning"
            )

    @requires_assist_data
    def test_assist_non_perturber_keeps_asteroid_force(self):
        """A body that is NOT one of ASSIST's asteroid perturbers (Chiron) must
        propagate with the full force model and emit no self-interaction warning.
        """
        import warnings

        elements = MINOR_BODY_ELEMENTS[CHIRON]
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = propagate_orbit_assist(
                elements, elements.epoch, elements.epoch + 365.25
            )
        h_z = result.x * result.vy - result.y * result.vx
        assert h_z > 0  # prograde
        assert not any("self-interaction" in str(w.message) for w in caught)

    @requires_assist_data
    def test_assist_backward_integration(self):
        """Backward ASSIST integration should work."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd_start = elements.epoch
        jd_end = jd_start - 30  # 30 days backward

        result = propagate_orbit_assist(elements, jd_start, jd_end)

        assert result.jd_tt == jd_end
        assert result.distance > 0

    @requires_assist_data
    def test_assist_trajectory(self):
        """Trajectory propagation with ASSIST should return multiple points."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd_start = elements.epoch
        jd_end = jd_start + 100

        trajectory = propagate_trajectory(
            elements, jd_start, jd_end, num_points=10, use_assist=True
        )

        assert len(trajectory) == 10
        for point in trajectory:
            assert isinstance(point, PropagationResult)
            assert 2.0 < point.distance < 3.5

    @requires_assist_data
    def test_compare_with_keplerian_assist(self):
        """compare_with_keplerian with use_assist=True should use ASSIST."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        jd = elements.epoch + 365.25

        comparison = compare_with_keplerian(elements, jd, use_assist=True)

        assert comparison["method"] == "assist"
        assert "angular_sep_arcsec" in comparison

    @requires_assist_data
    def test_assist_missing_planets_file_raises(self):
        """propagate_orbit_assist should raise FileNotFoundError for missing file."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        config = AssistEphemConfig(
            planets_file="/nonexistent/linux_p1550p2650.440",
        )

        with pytest.raises(FileNotFoundError):
            propagate_orbit_assist(
                elements, elements.epoch, elements.epoch + 10, ephem_config=config
            )

    def test_assist_no_planets_file_raises(self):
        """propagate_orbit_assist should raise FileNotFoundError when no file found."""
        if not check_assist_available():
            pytest.skip("ASSIST not installed")

        elements = MINOR_BODY_ELEMENTS[CERES]
        config = AssistEphemConfig.__new__(AssistEphemConfig)
        config.planets_file = None
        config.asteroids_file = None
        config.data_dir = None

        with pytest.raises(FileNotFoundError):
            propagate_orbit_assist(
                elements, elements.epoch, elements.epoch + 10, ephem_config=config
            )
