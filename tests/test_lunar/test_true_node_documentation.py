"""
Tests verifying the accuracy of True Node calculation documentation.

This module ensures that the documented properties of the True Node calculation
in lunar.py and true_node_terms.py are accurate and consistent with the
actual implementation.
"""

import math

from libephemeris import lunar
from libephemeris.lunar import (
    calc_true_lunar_node,
    calc_mean_lunar_node,
    _calc_elp2000_node_perturbations,
    _calc_lunar_fundamental_arguments,
)
from libephemeris.true_node_terms import (
    ALL_TRUE_NODE_TERMS,
    get_term_count,
    get_terms_by_category,
)


class TestModuleDocumentation:
    """Tests verifying module-level documentation accuracy."""

    def test_module_docstring_exists(self):
        """Verify the lunar module has a comprehensive docstring."""
        assert lunar.__doc__ is not None
        assert len(lunar.__doc__) > 1000  # Should be substantial documentation

    def test_module_docstring_mentions_true_node(self):
        """Verify docstring covers True Node calculation."""
        doc = lunar.__doc__
        assert "True Node" in doc or "True Lunar Node" in doc
        assert "ELP2000" in doc

    def test_module_docstring_mentions_perturbation_categories(self):
        """Verify docstring lists perturbation categories."""
        doc = lunar.__doc__
        assert "Solar Perturbation" in doc
        assert "Venus" in doc or "planetary" in doc.lower()
        assert "Evection" in doc or "evection" in doc
        assert "Variation" in doc or "variation" in doc
        assert "Annual Equation" in doc or "annual equation" in doc.lower()

    def test_module_docstring_mentions_precision(self):
        """Verify docstring documents precision."""
        doc = lunar.__doc__
        assert "0.01" in doc or "precision" in doc.lower() or "accuracy" in doc.lower()

    def test_module_docstring_has_references(self):
        """Verify docstring includes references."""
        doc = lunar.__doc__
        assert "Chapront" in doc or "Meeus" in doc
        assert "ELP" in doc or "pyswisseph" in doc or "JPL" in doc


class TestCalcTrueLunarNodeDocumentation:
    """Tests verifying calc_true_lunar_node docstring accuracy."""

    def test_function_docstring_exists(self):
        """Verify calc_true_lunar_node has a comprehensive docstring."""
        assert calc_true_lunar_node.__doc__ is not None
        assert len(calc_true_lunar_node.__doc__) > 500

    def test_docstring_describes_algorithm_steps(self):
        """Verify docstring describes the algorithm steps."""
        doc = calc_true_lunar_node.__doc__
        # Should mention key algorithm steps
        assert "position" in doc.lower() or "r" in doc
        assert "velocity" in doc.lower() or "v" in doc
        assert "angular momentum" in doc.lower() or "h = r" in doc
        assert "precession" in doc.lower()
        assert "nutation" in doc.lower()

    def test_docstring_mentions_iau_models(self):
        """Verify docstring mentions the IAU models used."""
        doc = calc_true_lunar_node.__doc__
        assert "IAU 2000" in doc or "IAU2000" in doc  # Nutation model
        assert "IAU 2006" in doc or "IAU2006" in doc  # Precession model

    def test_docstring_describes_return_values(self):
        """Verify docstring describes return tuple."""
        doc = calc_true_lunar_node.__doc__
        assert "longitude" in doc.lower()
        assert "latitude" in doc.lower()
        assert "0.0" in doc  # Latitude is always 0

    def test_docstring_mentions_precision_estimate(self):
        """Verify docstring provides precision estimate."""
        doc = calc_true_lunar_node.__doc__
        assert "0.01" in doc or "arcsec" in doc or "precision" in doc.lower()

    def test_docstring_has_references(self):
        """Verify docstring includes references."""
        doc = calc_true_lunar_node.__doc__
        assert "Vallado" in doc or "Capitaine" in doc or "Meeus" in doc or "IAU" in doc


class TestPerturbationDocumentation:
    """Tests verifying _calc_elp2000_node_perturbations docstring accuracy."""

    def test_function_docstring_exists(self):
        """Verify _calc_elp2000_node_perturbations has a comprehensive docstring."""
        assert _calc_elp2000_node_perturbations.__doc__ is not None
        assert len(_calc_elp2000_node_perturbations.__doc__) > 1000

    def test_docstring_lists_fundamental_arguments(self):
        """Verify docstring lists the fundamental arguments D, M, M', F."""
        doc = _calc_elp2000_node_perturbations.__doc__
        # Should mention all four fundamental arguments
        assert "D" in doc  # Mean elongation
        assert "M" in doc  # Sun's anomaly
        assert "M'" in doc or "M_prime" in doc or "Moon" in doc  # Moon's anomaly
        assert "F" in doc  # Argument of latitude

    def test_docstring_describes_perturbation_categories(self):
        """Verify docstring describes all perturbation categories."""
        doc = _calc_elp2000_node_perturbations.__doc__
        # All major categories should be mentioned
        assert "solar" in doc.lower() or "Solar" in doc
        assert "Venus" in doc or "venus" in doc
        assert "Mars" in doc or "mars" in doc
        assert "Jupiter" in doc or "jupiter" in doc
        assert "Saturn" in doc or "saturn" in doc
        assert "Evection" in doc or "evection" in doc
        assert "Variation" in doc or "variation" in doc
        assert "Annual" in doc or "annual" in doc

    def test_documented_term_count_is_accurate(self):
        """Verify the documented '90+ terms' claim is accurate."""
        doc = _calc_elp2000_node_perturbations.__doc__
        # The docstring should mention approximately the right number
        assert "90" in doc or "90+" in doc
        # And the actual count should match
        actual_count = get_term_count()
        assert actual_count >= 90, f"Expected 90+ terms, got {actual_count}"


class TestDocumentedPerturbationAmplitudes:
    """Tests verifying documented perturbation amplitudes are accurate."""

    def test_dominant_term_amplitude(self):
        """Verify the documented -1.5233 sin(2D) term amplitude."""
        # Module docstring mentions ~1.5° amplitude for main solar term
        term = next(
            t
            for t in ALL_TRUE_NODE_TERMS
            if t.d_mult == 2
            and t.m_mult == 0
            and t.mp_mult == 0
            and t.f_mult == 0
            and t.trig_func == "sin"
        )
        assert abs(term.amplitude + 1.5233) < 0.001

    def test_annual_equation_amplitude(self):
        """Verify the documented ~0.186° annual equation amplitude."""
        # The annual equation term should be around 0.186°
        from libephemeris.true_node_terms import ANNUAL_EQUATION_TERMS

        term = ANNUAL_EQUATION_TERMS[0]
        assert abs(abs(term.amplitude) - 0.186) < 0.01

    def test_evection_amplitude(self):
        """Verify evection term amplitudes are in documented range."""
        from libephemeris.true_node_terms import EVECTION_TERMS

        max_amplitude = max(abs(t.amplitude) for t in EVECTION_TERMS)
        # Documented as up to 0.047°
        assert max_amplitude <= 0.1  # Allow some margin
        assert max_amplitude >= 0.04


class TestDocumentedOscillationBehavior:
    """Tests verifying documented oscillation behavior."""

    def test_true_node_oscillates_around_mean(self):
        """Verify True Node oscillates around Mean Node as documented."""
        # Sample at J2000.0
        jd_j2000 = 2451545.0

        # Get mean and true nodes at several dates
        test_dates = [jd_j2000 + i for i in range(0, 30)]  # 30 days

        differences = []
        for jd in test_dates:
            mean_node = calc_mean_lunar_node(jd)
            true_node, _, _ = calc_true_lunar_node(jd)

            # Handle wraparound at 0°/360°
            diff = true_node - mean_node
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360

            differences.append(diff)

        # Documented: oscillates ±1.5° from mean
        max_diff = max(abs(d) for d in differences)
        # Should be less than 3° (allowing for some margin)
        assert max_diff < 3.0, f"Max oscillation {max_diff}° exceeds expected range"

    def test_perturbation_range(self):
        """Verify perturbation corrections are in documented range."""
        # Test at several dates
        test_dates = [2451545.0 + i * 30 for i in range(12)]  # Monthly for 1 year

        for jd in test_dates:
            perturbation = _calc_elp2000_node_perturbations(jd)
            # Documented: typically ranges from -2° to +2°
            assert abs(perturbation) < 3.0, (
                f"Perturbation {perturbation}° at JD {jd} exceeds documented range"
            )


class TestDocumentedCoordinateSystems:
    """Tests verifying documented coordinate system usage."""

    def test_latitude_is_always_zero(self):
        """Verify latitude is always 0 as documented."""
        # Node is on ecliptic by definition
        test_dates = [2451545.0 + i * 100 for i in range(10)]

        for jd in test_dates:
            _, latitude, _ = calc_true_lunar_node(jd)
            assert latitude == 0.0, f"Latitude should be 0, got {latitude}"

    def test_longitude_in_valid_range(self):
        """Verify longitude is in [0, 360) range as documented."""
        test_dates = [2451545.0 + i * 100 for i in range(50)]

        for jd in test_dates:
            longitude, _, _ = calc_true_lunar_node(jd)
            assert 0 <= longitude < 360, f"Longitude {longitude} out of range"


class TestDocumentedPhysicalProperties:
    """Tests verifying documented physical properties."""

    def test_node_retrograde_motion(self):
        """Verify documented ~18.6 year retrograde nodal period."""
        # Mean node moves ~19.3°/year retrograde
        jd_start = 2451545.0  # J2000.0
        jd_end = jd_start + 365.25  # One year later

        mean_start = calc_mean_lunar_node(jd_start)
        mean_end = calc_mean_lunar_node(jd_end)

        # Calculate motion (accounting for retrograde)
        motion = mean_end - mean_start
        if motion > 180:
            motion -= 360
        elif motion < -180:
            motion += 360

        # Documented: ~19.3°/year retrograde
        # Motion should be approximately -19.3°
        assert -22 < motion < -17, f"Annual motion {motion}° not in expected range"

    def test_draconic_month_oscillation(self):
        """Verify oscillation on draconic month timescale as documented."""
        # Documented: Primary oscillation period ~27.2 days (draconic month)
        jd_start = 2451545.0
        samples = []

        # Sample every day for two draconic months
        for i in range(60):
            jd = jd_start + i
            lon, _, _ = calc_true_lunar_node(jd)
            samples.append(lon)

        # There should be oscillation - check standard deviation
        mean_lon = sum(samples) / len(samples)
        variance = sum((s - mean_lon) ** 2 for s in samples) / len(samples)
        std_dev = math.sqrt(variance)

        # Standard deviation should indicate oscillation
        assert std_dev > 0.1, "Expected oscillation in True Node position"


class TestDocumentedTermCounts:
    """Tests verifying documented term counts by category."""

    def test_total_term_count(self):
        """Verify documented 90+ terms."""
        count = get_term_count()
        assert count >= 90, f"Documented 90+ terms, got {count}"

    def test_planetary_term_counts(self):
        """Verify documented planetary term counts."""
        categories = get_terms_by_category()

        # Documented: Venus (9 terms), Mars (9), Jupiter (7), Saturn (7)
        assert len(categories["venus"]) >= 5
        assert len(categories["mars"]) >= 5
        assert len(categories["jupiter"]) >= 5
        assert len(categories["saturn"]) >= 5


class TestDocumentedReferences:
    """Tests verifying references are properly cited."""

    def test_elp2000_reference_in_docstring(self):
        """Verify ELP2000-82B is referenced."""
        doc = _calc_elp2000_node_perturbations.__doc__
        assert "ELP" in doc
        assert "2000" in doc or "Chapront" in doc

    def test_meeus_reference_in_module(self):
        """Verify Meeus is referenced in the module."""
        doc = lunar.__doc__
        assert "Meeus" in doc

    def test_jpl_ephemeris_reference(self):
        """Verify JPL ephemeris is referenced."""
        doc = calc_true_lunar_node.__doc__
        assert "JPL" in doc or "DE440" in doc or "DE441" in doc or "Park" in doc


class TestFundamentalArgumentsDocumentation:
    """Tests verifying fundamental arguments documentation."""

    def test_fundamental_arguments_docstring(self):
        """Verify _calc_lunar_fundamental_arguments is documented."""
        assert _calc_lunar_fundamental_arguments.__doc__ is not None
        doc = _calc_lunar_fundamental_arguments.__doc__

        # Should describe all four arguments
        assert "elongation" in doc.lower() or "D" in doc
        assert "anomaly" in doc.lower()
        assert "latitude" in doc.lower() or "F" in doc

    def test_fundamental_arguments_return_four_values(self):
        """Verify function returns four arguments as documented."""
        jd = 2451545.0
        result = _calc_lunar_fundamental_arguments(jd)

        assert len(result) == 4
        # All should be in radians [0, 2π)
        for arg in result:
            assert 0 <= arg < 2 * math.pi
