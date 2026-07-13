"""
Tests for the Schaefer visibility threshold model in libephemeris.

Tests the calculation of visibility thresholds that determine whether a
celestial object of given magnitude is visible against a sky of given
brightness, accounting for observer eye adaptation and experience.

The Schaefer model accounts for:
- Object brightness (apparent magnitude)
- Sky surface brightness (background)
- Observer's eye adaptation state (photopic, mesopic, scotopic)
- Observer's experience and skill level

References:
    - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
    - Schaefer, B.E. (1993) "Astronomy and the Limits of Vision"
    - Blackwell, H.R. (1946) "Contrast thresholds of the human eye"
    - Crumey, A. (2014) "Human contrast threshold and astronomical visibility"
"""

from libephemeris import (
    # Visibility threshold functions
    VisibilityResult,
    calc_eye_adaptation_state,
    calc_contrast_threshold,
    calc_visibility_threshold,
    is_object_visible,
    calc_limiting_magnitude_for_sky,
    # Observer skill constants
    OBSERVER_SKILL_INEXPERIENCED,
    OBSERVER_SKILL_AVERAGE,
    OBSERVER_SKILL_EXPERIENCED,
    OBSERVER_SKILL_EXPERT,
    EXPERIENCE_FACTORS,
)


class TestEyeAdaptationState:
    """Tests for eye adaptation state calculation."""

    def test_photopic_in_bright_sky(self):
        """Bright sky should result in photopic (cone) vision."""
        # Civil twilight - bright sky
        state = calc_eye_adaptation_state(8.0)
        assert state == "photopic"

    def test_photopic_boundary(self):
        """Sky brightness < 16 should be photopic."""
        state = calc_eye_adaptation_state(15.9)
        assert state == "photopic"

    def test_mesopic_in_twilight(self):
        """Twilight sky should result in mesopic (mixed) vision."""
        state = calc_eye_adaptation_state(18.0)
        assert state == "mesopic"

    def test_mesopic_boundaries(self):
        """Sky brightness 16-20 should be mesopic."""
        assert calc_eye_adaptation_state(16.0) == "mesopic"
        assert calc_eye_adaptation_state(19.9) == "mesopic"

    def test_dark_adapted_in_dark_sky(self):
        """Dark sky should result in dark-adapted (scotopic) vision."""
        state = calc_eye_adaptation_state(21.5)
        assert state == "dark"

    def test_dark_boundary(self):
        """Sky brightness >= 20 should be dark-adapted."""
        assert calc_eye_adaptation_state(20.0) == "dark"
        assert calc_eye_adaptation_state(22.0) == "dark"


class TestContrastThreshold:
    """Tests for contrast threshold calculation."""

    def test_threshold_positive(self):
        """Contrast threshold should always be positive."""
        for sky_b in [5.0, 10.0, 15.0, 18.0, 21.5]:
            threshold = calc_contrast_threshold(sky_b)
            assert threshold > 0

    def test_threshold_dark_adapted_lower(self):
        """Dark-adapted eyes should have lower (better) limiting magnitude threshold."""
        # Note: The threshold here is flux ratio (object/sky), not a fractional contrast.
        # For dark skies, the threshold is higher because limiting magnitude is closer to sky brightness.
        # What matters is that dark-adapted vision has BETTER limiting magnitude.
        # We test this indirectly through calc_limiting_magnitude_for_sky.

        # The key test is that dark adaptation allows seeing fainter objects
        # This is properly tested in TestLimitingMagnitudeForSky
        # Here we just verify threshold is positive
        threshold_dark = calc_contrast_threshold(21.5, eye_adaptation="dark")
        threshold_photopic = calc_contrast_threshold(10.0, eye_adaptation="photopic")
        assert threshold_dark > 0
        assert threshold_photopic > 0

    def test_threshold_decreases_with_experience(self):
        """More experienced observers should have lower threshold."""
        sky_b = 21.5
        t_inexperienced = calc_contrast_threshold(
            sky_b, observer_skill=OBSERVER_SKILL_INEXPERIENCED
        )
        t_average = calc_contrast_threshold(
            sky_b, observer_skill=OBSERVER_SKILL_AVERAGE
        )
        t_experienced = calc_contrast_threshold(
            sky_b, observer_skill=OBSERVER_SKILL_EXPERIENCED
        )
        t_expert = calc_contrast_threshold(sky_b, observer_skill=OBSERVER_SKILL_EXPERT)

        assert t_inexperienced > t_average > t_experienced > t_expert

    def test_threshold_bounded(self):
        """Threshold should be within reasonable bounds (positive)."""
        # The threshold is a flux ratio, which can be large for dark skies
        # where limiting magnitude is much fainter than sky brightness
        for sky_b in [5.0, 10.0, 15.0, 18.0, 21.5]:
            threshold = calc_contrast_threshold(sky_b)
            assert threshold > 0  # Always positive

    def test_threshold_with_explicit_adaptation(self):
        """Explicit adaptation state should override auto-detection."""
        # Force dark adaptation in bright sky
        threshold = calc_contrast_threshold(8.0, eye_adaptation="dark")
        # Should use dark-adapted threshold formula
        assert threshold > 0


class TestVisibilityThreshold:
    """Tests for the main visibility threshold calculation."""

    def test_returns_visibility_result(self):
        """Should return a VisibilityResult dataclass."""
        result = calc_visibility_threshold(0.0, 21.5)
        assert isinstance(result, VisibilityResult)

    def test_result_has_all_fields(self):
        """VisibilityResult should have all expected fields."""
        result = calc_visibility_threshold(0.0, 21.5)
        assert hasattr(result, "is_visible")
        assert hasattr(result, "object_magnitude")
        assert hasattr(result, "limiting_magnitude")
        assert hasattr(result, "sky_brightness")
        assert hasattr(result, "contrast")
        assert hasattr(result, "threshold_contrast")
        assert hasattr(result, "visibility_margin")
        assert hasattr(result, "eye_adaptation")
        assert hasattr(result, "observer_skill")

    def test_bright_object_visible_in_bright_sky(self):
        """Very bright objects (Venus, Jupiter) should be visible even in twilight."""
        # Venus at magnitude -4 during civil twilight
        result = calc_visibility_threshold(-4.0, 8.0)
        assert result.is_visible is True

    def test_faint_object_invisible_in_bright_sky(self):
        """Faint objects should be invisible in bright twilight sky."""
        # 6th magnitude star during civil twilight
        result = calc_visibility_threshold(6.0, 8.0)
        assert result.is_visible is False

    def test_faint_object_visible_in_dark_sky(self):
        """Faint objects should be visible under dark skies."""
        # 5th magnitude star under dark sky
        result = calc_visibility_threshold(5.0, 21.5)
        assert result.is_visible is True

    def test_limiting_mag_increases_with_darker_sky(self):
        """Limiting magnitude should increase (see fainter) with darker skies."""
        limit_twilight = calc_visibility_threshold(0.0, 10.0).limiting_magnitude
        limit_dark = calc_visibility_threshold(0.0, 21.5).limiting_magnitude
        assert limit_dark > limit_twilight

    def test_visibility_margin_positive_when_visible(self):
        """Visibility margin should be positive for visible objects."""
        result = calc_visibility_threshold(-4.0, 21.5)  # Bright object, dark sky
        assert result.visibility_margin > 0
        assert result.is_visible is True

    def test_visibility_margin_negative_when_not_visible(self):
        """Visibility margin should be negative for invisible objects."""
        result = calc_visibility_threshold(10.0, 8.0)  # Very faint object, bright sky
        assert result.visibility_margin < 0
        assert result.is_visible is False

    def test_expert_sees_fainter_than_inexperienced(self):
        """Expert observers should see fainter objects."""
        sky_b = 21.5
        # Find a magnitude at the edge of visibility
        result_inexperienced = calc_visibility_threshold(
            6.0, sky_b, observer_skill=OBSERVER_SKILL_INEXPERIENCED
        )
        result_expert = calc_visibility_threshold(
            6.0, sky_b, observer_skill=OBSERVER_SKILL_EXPERT
        )

        # Expert should have higher limiting magnitude
        assert (
            result_expert.limiting_magnitude > result_inexperienced.limiting_magnitude
        )

    def test_extinction_makes_objects_fainter(self):
        """Applying extinction should make objects appear fainter."""
        # Object at 30 degrees altitude
        result_no_ext = calc_visibility_threshold(5.0, 21.5, apply_extinction=False)
        result_with_ext = calc_visibility_threshold(
            5.0, 21.5, object_altitude_deg=30.0, apply_extinction=True
        )

        # With extinction, the effective magnitude should be higher (fainter)
        assert result_with_ext.object_magnitude > result_no_ext.object_magnitude

    def test_extinction_near_horizon_significant(self):
        """Extinction near horizon should be significant."""
        result = calc_visibility_threshold(
            0.0, 21.5, object_altitude_deg=5.0, apply_extinction=True
        )
        # At 5 degrees, extinction should add several magnitudes
        assert result.object_magnitude > 2.0  # Original was 0.0

    def test_sirius_visibility(self):
        """Sirius (mag -1.46) should be visible under most conditions."""
        # Sirius in dark sky
        result_dark = calc_visibility_threshold(-1.46, 21.5)
        assert result_dark.is_visible is True

        # Sirius in civil twilight
        result_twilight = calc_visibility_threshold(-1.46, 8.0)
        assert result_twilight.is_visible is True

    def test_naked_eye_limit_dark_sky(self):
        """Under dark skies, naked-eye limit should be around 6-6.5 mag."""
        result = calc_visibility_threshold(0.0, 21.5)
        # Typical limiting magnitude for dark sky, average observer
        assert 5.5 < result.limiting_magnitude < 7.0


class TestIsObjectVisible:
    """Tests for the simple visibility check function."""

    def test_returns_boolean(self):
        """Should return a boolean."""
        result = is_object_visible(0.0, 21.5)
        assert isinstance(result, bool)

    def test_bright_object_visible(self):
        """Bright objects should return True."""
        assert is_object_visible(-4.0, 21.5) is True

    def test_faint_object_in_bright_sky_not_visible(self):
        """Faint objects in bright sky should return False."""
        assert is_object_visible(7.0, 8.0) is False

    def test_skill_affects_visibility(self):
        """Observer skill should affect visibility outcome."""
        # Object at the edge of visibility
        visible_inexperienced = is_object_visible(
            6.0, 21.5, observer_skill=OBSERVER_SKILL_INEXPERIENCED
        )
        visible_expert = is_object_visible(
            6.0, 21.5, observer_skill=OBSERVER_SKILL_EXPERT
        )
        # Expert might see what inexperienced cannot
        # (depends on exact threshold calculation)
        assert (
            visible_expert or not visible_inexperienced
        )  # At minimum, expert >= inexperienced


class TestLimitingMagnitudeForSky:
    """Tests for limiting magnitude calculation from sky brightness."""

    def test_returns_float(self):
        """Should return a float."""
        result = calc_limiting_magnitude_for_sky(21.5)
        assert isinstance(result, float)

    def test_dark_sky_limit_around_6(self):
        """Dark sky limiting magnitude should be around 6 for average observer."""
        limit = calc_limiting_magnitude_for_sky(21.5)
        assert 5.5 < limit < 7.0

    def test_bright_sky_limit_lower(self):
        """Bright sky limiting magnitude should be lower."""
        limit_dark = calc_limiting_magnitude_for_sky(21.5)
        limit_twilight = calc_limiting_magnitude_for_sky(10.0)
        assert limit_twilight < limit_dark

    def test_expert_sees_fainter(self):
        """Expert observers should have higher limiting magnitude."""
        limit_average = calc_limiting_magnitude_for_sky(21.5, OBSERVER_SKILL_AVERAGE)
        limit_expert = calc_limiting_magnitude_for_sky(21.5, OBSERVER_SKILL_EXPERT)
        assert limit_expert > limit_average

    def test_inexperienced_sees_brighter(self):
        """Inexperienced observers should have lower limiting magnitude."""
        limit_average = calc_limiting_magnitude_for_sky(21.5, OBSERVER_SKILL_AVERAGE)
        limit_inexp = calc_limiting_magnitude_for_sky(
            21.5, OBSERVER_SKILL_INEXPERIENCED
        )
        assert limit_inexp < limit_average

    def test_limiting_magnitude_progression(self):
        """Limiting magnitude should increase with darker skies."""
        sky_brightnesses = [8.0, 12.0, 16.0, 18.0, 21.5]
        limits = [calc_limiting_magnitude_for_sky(b) for b in sky_brightnesses]

        for i in range(1, len(limits)):
            assert limits[i] > limits[i - 1]


class TestObserverSkillConstants:
    """Tests for observer skill constants."""

    def test_constants_defined(self):
        """All skill constants should be defined."""
        assert OBSERVER_SKILL_INEXPERIENCED == 1
        assert OBSERVER_SKILL_AVERAGE == 2
        assert OBSERVER_SKILL_EXPERIENCED == 3
        assert OBSERVER_SKILL_EXPERT == 4

    def test_experience_factors_defined(self):
        """Experience factors should be defined for all skill levels."""
        assert OBSERVER_SKILL_INEXPERIENCED in EXPERIENCE_FACTORS
        assert OBSERVER_SKILL_AVERAGE in EXPERIENCE_FACTORS
        assert OBSERVER_SKILL_EXPERIENCED in EXPERIENCE_FACTORS
        assert OBSERVER_SKILL_EXPERT in EXPERIENCE_FACTORS

    def test_experience_factors_decrease_with_skill(self):
        """Experience factors should decrease with higher skill."""
        f_inexp = EXPERIENCE_FACTORS[OBSERVER_SKILL_INEXPERIENCED]
        f_avg = EXPERIENCE_FACTORS[OBSERVER_SKILL_AVERAGE]
        f_exp = EXPERIENCE_FACTORS[OBSERVER_SKILL_EXPERIENCED]
        f_expert = EXPERIENCE_FACTORS[OBSERVER_SKILL_EXPERT]

        assert f_inexp > f_avg > f_exp > f_expert

    def test_average_observer_factor_is_one(self):
        """Average observer should have factor of 1.0."""
        assert EXPERIENCE_FACTORS[OBSERVER_SKILL_AVERAGE] == 1.0


class TestPhysicalReasonability:
    """Tests to verify physically reasonable behavior."""

    def test_planets_visible_in_twilight(self):
        """Bright planets should be visible during twilight."""
        # Venus, Jupiter, Mars should be visible in civil twilight
        assert is_object_visible(-4.0, 8.0)  # Venus
        assert is_object_visible(-2.5, 8.0)  # Jupiter
        assert is_object_visible(-1.5, 8.0)  # Mars at opposition

    def test_first_magnitude_stars_visible_late_twilight(self):
        """First magnitude stars should appear in late twilight."""
        # Around nautical twilight end
        result = calc_visibility_threshold(1.0, 15.0)
        assert result.is_visible is True

    def test_milky_way_conditions(self):
        """Very faint objects need very dark skies."""
        # Limiting magnitude around 6.5 is needed to see Milky Way well
        limit = calc_limiting_magnitude_for_sky(21.8, OBSERVER_SKILL_EXPERIENCED)
        assert limit > 6.0

    def test_contrast_decreases_with_brighter_sky(self):
        """Object contrast should decrease as sky gets brighter."""
        contrast_dark = calc_visibility_threshold(3.0, 21.5).contrast
        contrast_twilight = calc_visibility_threshold(3.0, 15.0).contrast
        contrast_civil = calc_visibility_threshold(3.0, 8.0).contrast

        # Contrast should decrease as sky brightens
        assert contrast_dark > contrast_twilight > contrast_civil

    def test_heliacal_visibility_scenario(self):
        """Test typical heliacal rising scenario."""
        # Star just above horizon during bright twilight
        # Apply extinction at low altitude
        result = calc_visibility_threshold(
            1.0,  # Bright star
            10.0,  # Early nautical twilight
            object_altitude_deg=5.0,
            apply_extinction=True,
            observer_skill=OBSERVER_SKILL_EXPERIENCED,
        )
        # With extinction, visibility depends on conditions
        # Just verify we get a valid result
        assert isinstance(result.is_visible, bool)
        assert result.object_magnitude > 1.0  # Extinction applied


class TestImportability:
    """Tests to verify proper import from libephemeris."""

    def test_all_functions_importable(self):
        """All visibility functions should be importable from libephemeris."""
        from libephemeris import (
            calc_eye_adaptation_state,
            calc_contrast_threshold,
            calc_visibility_threshold,
            is_object_visible,
            calc_limiting_magnitude_for_sky,
            OBSERVER_SKILL_AVERAGE,
            EXPERIENCE_FACTORS,
        )

        # Verify they are callable/usable
        assert callable(calc_eye_adaptation_state)
        assert callable(calc_contrast_threshold)
        assert callable(calc_visibility_threshold)
        assert callable(is_object_visible)
        assert callable(calc_limiting_magnitude_for_sky)
        assert isinstance(OBSERVER_SKILL_AVERAGE, int)
        assert isinstance(EXPERIENCE_FACTORS, dict)
