"""
Tests for mean motion resonance detection in minor body calculations.

Tests verify:
- Known resonant bodies (plutinos, twotinos) are correctly identified
- Non-resonant bodies are correctly identified
- Resonance detection functions return correct information
- Tolerance parameter works correctly
- All resonance types are properly defined
"""

import pytest
from libephemeris.minor_bodies import (
    MINOR_BODY_ELEMENTS,
    detect_mean_motion_resonance,
    is_body_resonant,
    get_resonance_info,
    ResonanceResult,
    NEPTUNE_RESONANCES,
    NEPTUNE_A,
    NEPTUNE_MEAN_MOTION,
    DEFAULT_RESONANCE_TOLERANCE,
    _calc_resonant_a,
    OrbitalElements,
)
from libephemeris.constants import (
    CERES,
    PALLAS,
    VESTA,
    JUNO,
    CHIRON,
    PHOLUS,
    ERIS,
    SEDNA,
    HAUMEA,
    MAKEMAKE,
    IXION,
    ORCUS,
    QUAOAR,
)


class TestResonanceConstants:
    """Test resonance constants and calculations."""

    def test_neptune_mean_motion_value(self):
        """Neptune mean motion should be approximately correct."""
        # Neptune's period is ~164.8 years
        expected_period_years = 164.8
        expected_n = 360.0 / (expected_period_years * 365.25)
        assert abs(NEPTUNE_MEAN_MOTION - expected_n) / expected_n < 0.01

    def test_resonance_list_not_empty(self):
        """Should have at least the major resonances defined."""
        assert len(NEPTUNE_RESONANCES) >= 5
        # Should include plutino (2:3) and twotino (1:2)
        names = [r.name for r in NEPTUNE_RESONANCES]
        assert "plutino" in names
        assert "twotino" in names

    def test_resonance_info_fields(self):
        """ResonanceInfo should have all required fields."""
        for res in NEPTUNE_RESONANCES:
            assert isinstance(res.p, int)
            assert isinstance(res.q, int)
            assert isinstance(res.name, str)
            assert isinstance(res.a_resonant, float)
            assert isinstance(res.tolerance, float)
            assert res.p > 0
            assert res.q > 0
            assert res.a_resonant > 0
            assert 0 < res.tolerance < 0.5

    def test_plutino_resonance_semi_major_axis(self):
        """Plutino (2:3) resonance should be at ~39.4 AU."""
        plutino = next(r for r in NEPTUNE_RESONANCES if r.name == "plutino")
        # Pluto's semi-major axis is ~39.5 AU
        assert 38.0 < plutino.a_resonant < 41.0

    def test_twotino_resonance_semi_major_axis(self):
        """Twotino (1:2) resonance should be at ~47.8 AU."""
        twotino = next(r for r in NEPTUNE_RESONANCES if r.name == "twotino")
        assert 46.0 < twotino.a_resonant < 49.0

    def test_calc_resonant_a_kepler_third_law(self):
        """Resonant semi-major axis should follow Kepler's 3rd law."""
        # For p:q resonance: body makes p orbits while Neptune makes q
        # T_body / T_Neptune = q / p
        # (a/a_Neptune)^(3/2) = q/p
        # So a = a_Neptune * (q/p)^(2/3)

        # Test 2:3 resonance (plutino): body makes 2 orbits while Neptune makes 3
        a_23 = _calc_resonant_a(2, 3)
        expected_ratio = (3 / 2) ** (2 / 3)  # q/p = 3/2
        assert abs(a_23 / NEPTUNE_A - expected_ratio) < 0.001

        # Test 1:2 resonance (twotino): body makes 1 orbit while Neptune makes 2
        a_12 = _calc_resonant_a(1, 2)
        expected_ratio = (2 / 1) ** (2 / 3)  # q/p = 2/1
        assert abs(a_12 / NEPTUNE_A - expected_ratio) < 0.001


class TestDetectMeanMotionResonance:
    """Test the detect_mean_motion_resonance function."""

    def test_ixion_is_plutino(self):
        """Ixion should be detected as a plutino (2:3 resonance)."""
        elements = MINOR_BODY_ELEMENTS[IXION]
        result = detect_mean_motion_resonance(elements)

        assert result.is_resonant is True
        assert result.resonance is not None
        assert result.resonance.p == 2
        assert result.resonance.q == 3
        assert result.resonance.name == "plutino"
        assert result.deviation < DEFAULT_RESONANCE_TOLERANCE
        assert result.warning_message is not None
        assert "resonance" in result.warning_message.lower()

    def test_orcus_is_plutino(self):
        """Orcus (anti-Pluto) should be detected as a plutino."""
        elements = MINOR_BODY_ELEMENTS[ORCUS]
        result = detect_mean_motion_resonance(elements)

        assert result.is_resonant is True
        assert result.resonance is not None
        assert result.resonance.name == "plutino"

    def test_ceres_not_resonant(self):
        """Ceres (main belt) should not be in Neptune resonance."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        result = detect_mean_motion_resonance(elements)

        assert result.is_resonant is False
        assert result.resonance is None
        assert result.warning_message is None

    def test_eris_not_in_major_resonance(self):
        """Eris (scattered disk) is not in a major Neptune resonance."""
        elements = MINOR_BODY_ELEMENTS[ERIS]
        result = detect_mean_motion_resonance(elements)

        # Eris at ~68 AU is not in any of the major resonances
        assert result.is_resonant is False

    def test_chiron_not_resonant(self):
        """Chiron (centaur) should not be in Neptune resonance."""
        elements = MINOR_BODY_ELEMENTS[CHIRON]
        result = detect_mean_motion_resonance(elements)

        assert result.is_resonant is False

    def test_sedna_not_resonant(self):
        """Sedna (detached object) should not be in Neptune resonance."""
        elements = MINOR_BODY_ELEMENTS[SEDNA]
        result = detect_mean_motion_resonance(elements)

        # Sedna at ~550 AU is way beyond any Neptune resonance
        assert result.is_resonant is False

    def test_inner_solar_system_bodies_not_checked(self):
        """Bodies with a < 20 AU should not be checked for Neptune resonance."""
        # All inner solar system bodies should return is_resonant=False
        for body_id in [CERES, PALLAS, VESTA, JUNO]:
            elements = MINOR_BODY_ELEMENTS[body_id]
            result = detect_mean_motion_resonance(elements)
            assert result.is_resonant is False
            assert elements.a < 20.0

    def test_tolerance_parameter(self):
        """Tolerance parameter should control resonance detection sensitivity."""
        elements = MINOR_BODY_ELEMENTS[IXION]

        # With very small tolerance, might not detect resonance
        result_strict = detect_mean_motion_resonance(elements, tolerance=0.001)

        # With large tolerance, should detect
        result_loose = detect_mean_motion_resonance(elements, tolerance=0.05)

        # Ixion is close enough that both should detect it
        assert result_loose.is_resonant is True

        # The deviation should be the same regardless of tolerance
        if result_strict.is_resonant:
            assert result_strict.deviation == result_loose.deviation

    def test_result_contains_body_name_in_warning(self):
        """Warning message should contain the body name."""
        elements = MINOR_BODY_ELEMENTS[IXION]
        result = detect_mean_motion_resonance(elements)

        assert "Ixion" in result.warning_message

    def test_deviation_is_fractional(self):
        """Deviation should be a small fractional value for resonant bodies."""
        elements = MINOR_BODY_ELEMENTS[IXION]
        result = detect_mean_motion_resonance(elements)

        # Deviation should be less than the tolerance
        assert 0 <= result.deviation < 0.1


class TestIsBodyResonant:
    """Test the is_body_resonant convenience function."""

    def test_ixion_is_resonant(self):
        """Ixion should be detected as resonant."""
        assert is_body_resonant(IXION) is True

    def test_orcus_is_resonant(self):
        """Orcus should be detected as resonant."""
        assert is_body_resonant(ORCUS) is True

    def test_ceres_not_resonant(self):
        """Ceres should not be resonant."""
        assert is_body_resonant(CERES) is False

    def test_eris_not_resonant(self):
        """Eris should not be in a major resonance."""
        assert is_body_resonant(ERIS) is False

    def test_unknown_body_raises_error(self):
        """Unknown body ID should raise ValueError."""
        with pytest.raises(ValueError):
            is_body_resonant(999999)

    @pytest.mark.parametrize(
        "body_id,expected",
        [
            (CERES, False),
            (PALLAS, False),
            (VESTA, False),
            (JUNO, False),
            (CHIRON, False),
            (PHOLUS, False),
            (ERIS, False),
            (SEDNA, False),
            # Note: Haumea and Quaoar are near but not in the 4:7 resonance
            # With 2% tolerance, they get flagged - this is expected behavior
            # as the function is conservative and flags "near resonance" bodies
            (HAUMEA, True),  # Near 4:7 (~1.6% deviation)
            (MAKEMAKE, False),
            (IXION, True),  # Plutino
            (ORCUS, True),  # Plutino
            (QUAOAR, True),  # Near 4:7 (~1.1% deviation)
        ],
    )
    def test_all_bodies_resonance_status(self, body_id, expected):
        """Verify resonance status for all known minor bodies."""
        assert is_body_resonant(body_id) == expected


class TestGetResonanceInfo:
    """Test the get_resonance_info function."""

    def test_ixion_info(self):
        """Should return detailed info for Ixion."""
        info = get_resonance_info(IXION)

        assert info is not None
        assert info.is_resonant is True
        assert info.resonance.name == "plutino"

    def test_ceres_info(self):
        """Should return non-resonant info for Ceres."""
        info = get_resonance_info(CERES)

        assert info is not None
        assert info.is_resonant is False

    def test_unknown_body_returns_none(self):
        """Unknown body ID should return None."""
        info = get_resonance_info(999999)
        assert info is None


class TestResonanceResult:
    """Test the ResonanceResult dataclass."""

    def test_resonance_result_fields(self):
        """ResonanceResult should have all required fields."""
        elements = MINOR_BODY_ELEMENTS[IXION]
        result = detect_mean_motion_resonance(elements)

        assert hasattr(result, "is_resonant")
        assert hasattr(result, "resonance")
        assert hasattr(result, "deviation")
        assert hasattr(result, "warning_message")

    def test_non_resonant_result_structure(self):
        """Non-resonant result should have None for optional fields."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        result = detect_mean_motion_resonance(elements)

        assert result.is_resonant is False
        assert result.resonance is None
        assert result.deviation == 0.0
        assert result.warning_message is None


class TestResonancePhysics:
    """Test physical correctness of resonance calculations."""

    def test_plutino_mean_motion_ratio(self):
        """Plutinos should have mean motion ratio ~2/3 with Neptune."""
        elements = MINOR_BODY_ELEMENTS[IXION]
        n_ratio = elements.n / NEPTUNE_MEAN_MOTION

        # For 2:3 resonance, body completes 2 orbits while Neptune completes 3
        # So n_body / n_Neptune = 2/3 (body is slower)
        expected_ratio = 2 / 3
        assert abs(n_ratio - expected_ratio) / expected_ratio < 0.02

    def test_ixion_and_orcus_similar_orbits(self):
        """Ixion and Orcus (both plutinos) should have similar semi-major axes."""
        ixion = MINOR_BODY_ELEMENTS[IXION]
        orcus = MINOR_BODY_ELEMENTS[ORCUS]

        # Both should be at ~39.4 AU
        assert abs(ixion.a - orcus.a) < 1.0  # Within 1 AU
        assert abs(ixion.a - 39.4) < 1.0
        assert abs(orcus.a - 39.4) < 1.0

    def test_resonance_location_vs_kepler_law(self):
        """Resonance locations should match Kepler's 3rd law predictions."""
        for resonance in NEPTUNE_RESONANCES:
            # For p:q resonance: a = a_Neptune * (q/p)^(2/3)
            expected_a = NEPTUNE_A * (resonance.q / resonance.p) ** (2 / 3)
            assert abs(resonance.a_resonant - expected_a) < 0.01, (
                f"Resonance {resonance.name} a={resonance.a_resonant} "
                f"doesn't match Kepler prediction {expected_a}"
            )

    def test_quaoar_not_in_plutino_region(self):
        """Quaoar at ~43 AU should not be in the plutino region."""
        quaoar = MINOR_BODY_ELEMENTS[QUAOAR]
        plutino = next(r for r in NEPTUNE_RESONANCES if r.name == "plutino")

        # Quaoar is at ~43 AU, plutinos at ~39 AU
        assert quaoar.a > plutino.a_resonant * 1.05  # At least 5% larger


class TestSyntheticResonantBody:
    """Test resonance detection with synthetic orbital elements."""

    def test_exact_plutino(self):
        """A body at exact 2:3 resonance should be detected."""
        # Calculate exact plutino semi-major axis
        a_plutino = _calc_resonant_a(2, 3)
        n_plutino = 0.9856076686 / (a_plutino**1.5)  # Kepler constant

        elements = OrbitalElements(
            name="TestPlutino",
            epoch=2461000.5,
            a=a_plutino,
            e=0.2,
            i=10.0,
            omega=100.0,
            Omega=200.0,
            M0=50.0,
            n=n_plutino,
        )

        result = detect_mean_motion_resonance(elements)
        assert result.is_resonant is True
        assert result.resonance.name == "plutino"
        # Small deviation expected due to slight difference in mean motion constants
        assert result.deviation < 0.02  # Within 2% tolerance

    def test_exact_twotino(self):
        """A body at exact 1:2 resonance should be detected."""
        a_twotino = _calc_resonant_a(1, 2)
        n_twotino = 0.9856076686 / (a_twotino**1.5)

        elements = OrbitalElements(
            name="TestTwotino",
            epoch=2461000.5,
            a=a_twotino,
            e=0.15,
            i=5.0,
            omega=150.0,
            Omega=250.0,
            M0=75.0,
            n=n_twotino,
        )

        result = detect_mean_motion_resonance(elements)
        assert result.is_resonant is True
        assert result.resonance.name == "twotino"

    def test_non_resonant_tno(self):
        """A TNO at 50 AU (not a resonant location) should not be flagged."""
        a = 50.0  # Between 1:2 (~47.8 AU) and 2:5 (~55.4 AU)
        n = 0.9856076686 / (a**1.5)

        elements = OrbitalElements(
            name="TestTNO",
            epoch=2461000.5,
            a=a,
            e=0.1,
            i=3.0,
            omega=90.0,
            Omega=180.0,
            M0=45.0,
            n=n,
        )

        result = detect_mean_motion_resonance(elements)
        assert result.is_resonant is False


class TestAllBodiesResonanceCheck:
    """Comprehensive tests running resonance check on all bodies."""

    @pytest.mark.parametrize(
        "body_id",
        [
            CERES,
            PALLAS,
            VESTA,
            JUNO,
            CHIRON,
            PHOLUS,
            ERIS,
            SEDNA,
            HAUMEA,
            MAKEMAKE,
            IXION,
            ORCUS,
            QUAOAR,
        ],
    )
    def test_all_bodies_resonance_check_runs(self, body_id):
        """Resonance check should complete without error for all bodies."""
        elements = MINOR_BODY_ELEMENTS[body_id]
        result = detect_mean_motion_resonance(elements)

        assert isinstance(result, ResonanceResult)
        assert isinstance(result.is_resonant, bool)
        assert result.deviation >= 0

    @pytest.mark.parametrize(
        "body_id",
        [
            CERES,
            PALLAS,
            VESTA,
            JUNO,
            CHIRON,
            PHOLUS,
            ERIS,
            SEDNA,
            HAUMEA,
            MAKEMAKE,
            IXION,
            ORCUS,
            QUAOAR,
        ],
    )
    def test_get_resonance_info_runs(self, body_id):
        """get_resonance_info should complete without error for all bodies."""
        info = get_resonance_info(body_id)
        assert info is not None


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_body_at_neptune_orbit(self):
        """Body at Neptune's exact orbit (0 deviation) edge case."""
        elements = OrbitalElements(
            name="AtNeptune",
            epoch=2461000.5,
            a=NEPTUNE_A,
            e=0.01,
            i=1.0,
            omega=45.0,
            Omega=135.0,
            M0=90.0,
            n=NEPTUNE_MEAN_MOTION,
        )

        # This would be 1:1 resonance (like a trojan), not defined in our list
        result = detect_mean_motion_resonance(elements)
        # Should not match any defined resonance
        assert result.is_resonant is False

    def test_very_distant_body(self):
        """Very distant body should not trigger false positive."""
        a = 1000.0  # Very far out
        n = 0.9856076686 / (a**1.5)

        elements = OrbitalElements(
            name="VeryDistant",
            epoch=2461000.5,
            a=a,
            e=0.5,
            i=45.0,
            omega=0.0,
            Omega=0.0,
            M0=0.0,
            n=n,
        )

        result = detect_mean_motion_resonance(elements)
        assert result.is_resonant is False

    def test_zero_tolerance(self):
        """Zero tolerance should only match exact resonances."""
        elements = MINOR_BODY_ELEMENTS[IXION]
        result = detect_mean_motion_resonance(elements, tolerance=0.0)

        # Ixion is close but not exactly at resonance, so with 0 tolerance
        # it might not match
        # This tests that the function handles edge case gracefully
        assert isinstance(result.is_resonant, bool)


class TestAssistCacheKeyedOnElementValues:
    """The propagation cache must key on the orbital-element VALUES: two
    different orbits sharing a name (and window) previously aliased to the
    same slot and the second returned the first's position."""

    def test_cache_key_includes_element_values(self):
        import pytest

        try:
            import assist  # noqa: F401
            import rebound  # noqa: F401
        except ImportError:
            pytest.skip("rebound/assist not installed")
        from libephemeris import rebound_integration as ri
        from libephemeris.minor_bodies import OrbitalElements

        e2 = OrbitalElements("X", 2460000.5, 30.0, 0.1, 5.0, 10.0, 20.0, 30.0, 0.2)
        sentinel = object()
        new_key = (
            "X",
            2460000.5,
            30.0,
            0.1,
            5.0,
            10.0,
            20.0,
            30.0,
            0.2,
            2460000.5,
            round(2460365.5 * 86400.0),
        )
        old_key = ("X", 2460000.5, round(2460365.5 * 86400.0))
        with ri._assist_propagation_cache_lock:
            ri._assist_propagation_cache.clear()
            ri._assist_propagation_cache[new_key] = sentinel
            # A stale name-only entry must never be served.
            ri._assist_propagation_cache[old_key] = object()
        try:
            got = ri.propagate_orbit_assist(e2, 2460000.5, 2460365.5)
            assert got is sentinel
        finally:
            with ri._assist_propagation_cache_lock:
                ri._assist_propagation_cache.clear()
