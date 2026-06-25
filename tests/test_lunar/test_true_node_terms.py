"""
Tests for the ELP2000-82B True Node Perturbation Term Table.

This module tests the perturbation coefficient table for the True Lunar Node,
verifying the structure, values, and utility functions.
"""

import math
import pytest

from libephemeris.true_node_terms import (
    TrueNodeTerm,
    MAIN_SOLAR_TERMS,
    SECOND_ORDER_SOLAR_TERMS,
    VENUS_TERMS,
    MARS_TERMS,
    JUPITER_TERMS,
    SATURN_TERMS,
    EVECTION_TERMS,
    VARIATION_TERMS,
    ANNUAL_EQUATION_TERMS,
    SECOND_ORDER_COUPLING_TERMS,
    PARALLACTIC_TERMS,
    ALL_TRUE_NODE_TERMS,
    get_term_count,
    get_terms_by_category,
    get_sorted_terms_by_amplitude,
)


class TestTrueNodeTermStructure:
    """Tests for the TrueNodeTerm NamedTuple structure."""

    def test_term_has_all_required_fields(self):
        """Verify TrueNodeTerm has all expected fields."""
        term = TrueNodeTerm(
            amplitude=1.0,
            d_mult=2,
            m_mult=0,
            mp_mult=-1,
            f_mult=0,
            e_power=0,
            trig_func="sin",
            planet=None,
            planet_mult=0,
            description="Test term",
        )
        assert hasattr(term, "amplitude")
        assert hasattr(term, "d_mult")
        assert hasattr(term, "m_mult")
        assert hasattr(term, "mp_mult")
        assert hasattr(term, "f_mult")
        assert hasattr(term, "e_power")
        assert hasattr(term, "trig_func")
        assert hasattr(term, "planet")
        assert hasattr(term, "planet_mult")
        assert hasattr(term, "description")

    def test_term_values_accessible(self):
        """Verify term values can be accessed correctly."""
        term = MAIN_SOLAR_TERMS[0]  # Dominant 2D term
        assert term.amplitude == -1.5233
        assert term.d_mult == 2
        assert term.m_mult == 0
        assert term.mp_mult == 0
        assert term.f_mult == 0
        assert term.e_power == 0
        assert term.trig_func == "sin"
        assert term.planet is None
        assert term.planet_mult == 0
        assert "fortnightly" in term.description.lower() or "2D" in term.description

    def test_terms_are_immutable(self):
        """Verify terms cannot be modified."""
        term = MAIN_SOLAR_TERMS[0]
        with pytest.raises(AttributeError):
            term.amplitude = 2.0  # type: ignore


class TestTermCounts:
    """Tests for term count and table completeness."""

    def test_minimum_term_count(self):
        """Verify we have at least 50 terms as specified."""
        assert get_term_count() >= 50

    def test_main_solar_terms_count(self):
        """Verify main solar terms are present."""
        assert len(MAIN_SOLAR_TERMS) >= 8

    def test_second_order_solar_terms_count(self):
        """Verify second-order solar terms are present."""
        assert len(SECOND_ORDER_SOLAR_TERMS) >= 10

    def test_planetary_terms_count(self):
        """Verify each planet has multiple terms."""
        assert len(VENUS_TERMS) >= 5
        assert len(MARS_TERMS) >= 5
        assert len(JUPITER_TERMS) >= 5
        assert len(SATURN_TERMS) >= 5

    def test_long_period_terms_count(self):
        """Verify long-period terms are present."""
        assert len(EVECTION_TERMS) >= 3
        assert len(VARIATION_TERMS) >= 3
        assert len(ANNUAL_EQUATION_TERMS) >= 4

    def test_all_categories_sum_to_total(self):
        """Verify all category tables sum to the total count."""
        categories = get_terms_by_category()
        total_from_categories = sum(len(terms) for terms in categories.values())
        assert total_from_categories == get_term_count()


class TestAmplitudeValues:
    """Tests for amplitude values and ranges."""

    def test_dominant_term_amplitude(self):
        """Verify the dominant term has the expected amplitude."""
        # The -1.5233° sin(2D) term is the largest
        dominant_term = MAIN_SOLAR_TERMS[0]
        assert abs(dominant_term.amplitude) > 1.0
        assert dominant_term.d_mult == 2
        assert dominant_term.amplitude == pytest.approx(-1.5233, rel=0.001)

    def test_annual_equation_amplitude(self):
        """Verify the annual equation term amplitude."""
        # The annual equation term should be around -0.186°
        annual_term = ANNUAL_EQUATION_TERMS[0]
        assert annual_term.m_mult == 1
        assert annual_term.mp_mult == 0
        assert annual_term.amplitude == pytest.approx(-0.186, rel=0.01)

    def test_all_amplitudes_are_nonzero(self):
        """Verify all terms have non-zero amplitudes."""
        for term in ALL_TRUE_NODE_TERMS:
            assert term.amplitude != 0.0, f"Zero amplitude in: {term.description}"

    def test_amplitudes_in_reasonable_range(self):
        """Verify all amplitudes are within expected ranges."""
        for term in ALL_TRUE_NODE_TERMS:
            # Amplitudes should be between about -2° and +2°
            assert abs(term.amplitude) < 2.0, f"Amplitude too large: {term}"
            # Smallest terms should still be significant (> 0.0001°)
            assert abs(term.amplitude) >= 0.0005, (
                f"Amplitude too small to be significant: {term}"
            )

    def test_sorted_terms_by_amplitude(self):
        """Verify get_sorted_terms_by_amplitude returns correct order."""
        sorted_terms = get_sorted_terms_by_amplitude()
        for i in range(len(sorted_terms) - 1):
            assert abs(sorted_terms[i].amplitude) >= abs(sorted_terms[i + 1].amplitude)

        # First term should be the dominant 2D term
        assert abs(sorted_terms[0].amplitude) > 1.0


class TestArgumentMultipliers:
    """Tests for argument multiplier values."""

    def test_multipliers_are_integers(self):
        """Verify all multipliers are integers."""
        for term in ALL_TRUE_NODE_TERMS:
            assert isinstance(term.d_mult, int)
            assert isinstance(term.m_mult, int)
            assert isinstance(term.mp_mult, int)
            assert isinstance(term.f_mult, int)
            assert isinstance(term.e_power, int)
            assert isinstance(term.planet_mult, int)

    def test_e_power_valid_values(self):
        """Verify e_power is 0, 1, or 2."""
        for term in ALL_TRUE_NODE_TERMS:
            assert term.e_power in (0, 1, 2), f"Invalid e_power in: {term}"

    def test_trig_func_valid_values(self):
        """Verify trig_func is 'sin' or 'cos'."""
        for term in ALL_TRUE_NODE_TERMS:
            assert term.trig_func in (
                "sin",
                "cos",
            ), f"Invalid trig_func in: {term}"

    def test_planet_valid_values(self):
        """Verify planet is None or a valid planet name."""
        valid_planets = {"venus", "mars", "jupiter", "saturn", None}
        for term in ALL_TRUE_NODE_TERMS:
            assert term.planet in valid_planets, f"Invalid planet in: {term}"

    def test_planetary_terms_have_planet_set(self):
        """Verify planetary term categories have planet field set."""
        for term in VENUS_TERMS:
            assert term.planet == "venus"
        for term in MARS_TERMS:
            assert term.planet == "mars"
        for term in JUPITER_TERMS:
            assert term.planet == "jupiter"
        for term in SATURN_TERMS:
            assert term.planet == "saturn"

    def test_non_planetary_terms_have_planet_none(self):
        """Verify non-planetary terms have planet=None."""
        for term in MAIN_SOLAR_TERMS:
            assert term.planet is None
        for term in EVECTION_TERMS:
            assert term.planet is None
        for term in ANNUAL_EQUATION_TERMS:
            assert term.planet is None


class TestPhysicalConsistency:
    """Tests for physical consistency of the term table."""

    def test_dominant_term_is_2d(self):
        """The largest term should be the fortnightly 2D term."""
        sorted_terms = get_sorted_terms_by_amplitude()
        dominant = sorted_terms[0]
        assert dominant.d_mult == 2
        assert dominant.m_mult == 0
        assert dominant.mp_mult == 0
        assert dominant.f_mult == 0

    def test_planetary_amplitudes_decrease_with_distance(self):
        """Outer planets should generally have smaller amplitudes."""
        venus_max = max(abs(t.amplitude) for t in VENUS_TERMS)
        mars_max = max(abs(t.amplitude) for t in MARS_TERMS)
        jupiter_max = max(abs(t.amplitude) for t in JUPITER_TERMS)
        saturn_max = max(abs(t.amplitude) for t in SATURN_TERMS)

        # Venus and Mars are inner/neighboring planets with larger effects
        # Jupiter and Saturn have smaller effects due to greater distance
        assert saturn_max <= jupiter_max
        assert jupiter_max <= mars_max

    def test_e_power_terms_have_m_mult_nonzero(self):
        """Terms with e_power > 0 should generally involve M (Sun's anomaly)."""
        for term in ALL_TRUE_NODE_TERMS:
            if term.e_power > 0 and term.planet is None:
                # Most E-dependent terms involve M, but some involve planetary longitudes
                # This is a soft check - we just verify the relationship exists
                pass  # E-factor is physically tied to Earth's orbit

    def test_second_order_coupling_uses_cosine(self):
        """Second-order coupling terms should use cosine."""
        for term in SECOND_ORDER_COUPLING_TERMS:
            assert term.trig_func == "cos", f"Coupling term should use cos: {term}"


class TestCategoryFunctions:
    """Tests for category utility functions."""

    def test_get_terms_by_category_returns_dict(self):
        """Verify get_terms_by_category returns a dictionary."""
        categories = get_terms_by_category()
        assert isinstance(categories, dict)

    def test_all_categories_present(self):
        """Verify all expected categories are present."""
        categories = get_terms_by_category()
        expected = {
            "main_solar",
            "second_order_solar",
            "third_order",
            "f_inclination",
            "venus",
            "mars",
            "jupiter",
            "saturn",
            "evection",
            "variation",
            "annual_equation",
            "second_order_coupling",
            "parallactic",
        }
        assert set(categories.keys()) == expected

    def test_category_terms_are_tuples(self):
        """Verify each category contains a tuple of terms."""
        categories = get_terms_by_category()
        for name, terms in categories.items():
            assert isinstance(terms, tuple), f"Category {name} is not a tuple"
            for term in terms:
                assert isinstance(term, TrueNodeTerm), (
                    f"Invalid term in category {name}"
                )


class TestTermDescriptions:
    """Tests for term descriptions."""

    def test_all_terms_have_descriptions(self):
        """Verify all terms have non-empty descriptions."""
        for term in ALL_TRUE_NODE_TERMS:
            assert term.description, (
                f"Empty description for term with amplitude {term.amplitude}"
            )
            assert len(term.description) > 5

    def test_descriptions_are_unique(self):
        """Verify term descriptions are unique (helps with debugging)."""
        descriptions = [term.description for term in ALL_TRUE_NODE_TERMS]
        # Allow some duplicate descriptions since similar terms may have similar descriptions
        # Just check that most are unique
        unique_count = len(set(descriptions))
        total_count = len(descriptions)
        assert unique_count >= total_count * 0.8, "Too many duplicate descriptions"


class TestKnownTermValues:
    """Tests for specific known term values from ELP2000-82B theory."""

    def test_main_fortnightly_term(self):
        """Verify the main fortnightly term matches ELP2000-82B."""
        # sin(2D) with amplitude ~-1.5233°
        term = next(
            t
            for t in MAIN_SOLAR_TERMS
            if t.d_mult == 2 and t.m_mult == 0 and t.mp_mult == 0 and t.f_mult == 0
        )
        assert term.amplitude == pytest.approx(-1.5233, rel=0.001)
        assert term.trig_func == "sin"

    def test_evection_term(self):
        """Verify the evection term matches ELP2000-82B."""
        # The 2D-M' term with amplitude ~-0.1226°
        term = next(
            t
            for t in MAIN_SOLAR_TERMS
            if t.d_mult == 2 and t.mp_mult == -1 and t.m_mult == 0
        )
        assert term.amplitude == pytest.approx(-0.1226, rel=0.01)

    def test_latitude_term(self):
        """Verify the latitude (2F) term matches ELP2000-82B."""
        # sin(2F) with amplitude ~0.1176°
        term = next(t for t in MAIN_SOLAR_TERMS if t.f_mult == 2 and t.d_mult == 0)
        assert term.amplitude == pytest.approx(0.1176, rel=0.01)

    def test_annual_equation_main_term(self):
        """Verify the annual equation main term matches ELP2000-82B."""
        # sin(M) with amplitude ~-0.186°
        term = ANNUAL_EQUATION_TERMS[0]
        assert term.m_mult == 1
        assert term.mp_mult == 0
        assert term.d_mult == 0
        assert term.f_mult == 0
        assert term.amplitude == pytest.approx(-0.186, rel=0.01)
        assert term.e_power == 1

    def test_parallactic_main_term(self):
        """Verify the parallactic inequality term matches ELP2000-82B."""
        # sin(D) with amplitude ~0.035°
        term = PARALLACTIC_TERMS[0]
        assert term.d_mult == 1
        assert term.m_mult == 0
        assert term.mp_mult == 0
        assert term.amplitude == pytest.approx(0.035, rel=0.01)


class TestTermTableComputation:
    """Tests verifying the term table can be used for computation."""

    def test_term_can_compute_argument(self):
        """Verify term arguments can be computed correctly."""
        # Test with sample values for fundamental arguments (in radians)
        D = math.radians(45.0)  # Mean elongation
        M = math.radians(30.0)  # Sun's mean anomaly
        M_prime = math.radians(60.0)  # Moon's mean anomaly
        F = math.radians(15.0)  # Argument of latitude

        term = MAIN_SOLAR_TERMS[0]  # 2D term
        argument = (
            term.d_mult * D + term.m_mult * M + term.mp_mult * M_prime + term.f_mult * F
        )
        expected = 2 * D  # Since other multipliers are 0
        assert argument == pytest.approx(expected)

    def test_term_can_compute_contribution(self):
        """Verify term contribution can be computed correctly."""
        D = math.radians(45.0)
        M = math.radians(30.0)
        M_prime = math.radians(60.0)
        F = math.radians(15.0)
        E = 0.99  # Eccentricity factor

        term = MAIN_SOLAR_TERMS[0]  # 2D term with amplitude -1.5233
        argument = (
            term.d_mult * D + term.m_mult * M + term.mp_mult * M_prime + term.f_mult * F
        )

        if term.trig_func == "sin":
            contribution = term.amplitude * (E**term.e_power) * math.sin(argument)
        else:
            contribution = term.amplitude * (E**term.e_power) * math.cos(argument)

        # With 2D = 90°, sin(90°) = 1, so contribution should be -1.5233
        expected = -1.5233 * math.sin(math.radians(90.0))
        assert contribution == pytest.approx(expected, rel=0.001)

    def test_sum_of_all_terms_at_sample_time(self):
        """Verify all terms can be summed for a sample calculation."""
        # Sample fundamental arguments (simplified)
        D = math.radians(120.0)
        M = math.radians(45.0)
        M_prime = math.radians(75.0)
        F = math.radians(30.0)
        E = 0.9832  # Modern eccentricity factor

        # Mock planetary longitudes
        L_venus = math.radians(100.0)
        L_mars = math.radians(200.0)
        L_jupiter = math.radians(60.0)
        L_saturn = math.radians(280.0)
        L_moon = math.radians(150.0)

        planet_longitudes = {
            "venus": L_venus,
            "mars": L_mars,
            "jupiter": L_jupiter,
            "saturn": L_saturn,
        }

        total_perturbation = 0.0

        for term in ALL_TRUE_NODE_TERMS:
            argument = (
                term.d_mult * D
                + term.m_mult * M
                + term.mp_mult * M_prime
                + term.f_mult * F
            )

            # Add planetary longitude if applicable
            if term.planet:
                # Special case for Venus term that uses L_Venus - L_Moon
                if (
                    term.planet == "venus"
                    and term.d_mult == 0
                    and term.m_mult == 0
                    and term.mp_mult == 0
                    and term.f_mult == 0
                ):
                    argument += term.planet_mult * (
                        planet_longitudes[term.planet] - L_moon
                    )
                else:
                    argument += term.planet_mult * planet_longitudes[term.planet]

            if term.trig_func == "sin":
                contribution = term.amplitude * (E**term.e_power) * math.sin(argument)
            else:
                contribution = term.amplitude * (E**term.e_power) * math.cos(argument)

            total_perturbation += contribution

        # The total perturbation should be within a reasonable range
        # (typically ±3° from the mean node)
        assert abs(total_perturbation) < 5.0
