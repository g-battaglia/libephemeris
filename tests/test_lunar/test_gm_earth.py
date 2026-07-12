"""
Tests for gravitational parameters (GM) in lunar calculations.

The GM values used in orbital calculations must match the IAU 2015 Resolution B3
TDB-compatible values used by the DE440 ephemeris for consistency.

For the lunar apogee calculation, we need the combined Earth-Moon gravitational
parameter because the eccentricity vector formula uses μ = G(M_Earth + M_Moon)
for the two-body problem.

Reference:
    IAU 2015 Resolution B3: Recommended Nominal Values for Selected Solar and
    Planetary Quantities
    GM_Earth (TDB) = 398600.435436 km^3/s^2
    Earth/Moon mass ratio = 81.3005691
"""

import re
import math
from pathlib import Path


# The IAU 2015 Resolution B3 TDB-compatible GM value for Earth in km^3/s^2
# This is the value used by DE440/DE441 ephemeris files
GM_EARTH_IAU2015_KM3_S2 = 398600.435436

# Earth/Moon mass ratio (IAU 2015 Resolution B3)
EARTH_MOON_MASS_RATIO_IAU2015 = 81.3005691

# Derived: GM_Moon from Earth/Moon mass ratio
GM_MOON_IAU2015_KM3_S2 = GM_EARTH_IAU2015_KM3_S2 / EARTH_MOON_MASS_RATIO_IAU2015

# Combined Earth-Moon gravitational parameter for two-body problem
# μ = G(M_Earth + M_Moon) = GM_Earth + GM_Moon
GM_EARTH_MOON_IAU2015_KM3_S2 = GM_EARTH_IAU2015_KM3_S2 + GM_MOON_IAU2015_KM3_S2

# Conversion factors
AU_KM = 149597870.7  # km per AU (IAU 2012 definition)
SECONDS_PER_DAY = 86400.0


class TestGMEarthConstant:
    """Tests for the GM_Earth constant value used in lunar.py."""

    def test_gm_earth_value_in_source_code(self):
        """Verify that lunar.py builds its GM from the IAU 2015 GM_Earth value."""
        lunar_py = Path(__file__).parent.parent.parent / "libephemeris" / "lunar.py"
        content = lunar_py.read_text()

        # The shared constant embeds the IAU 2015 GM_Earth figure verbatim
        assert "398600.435436" in content, (
            "lunar.py should derive GM_EARTH_MOON_AU3_DAY2 from the IAU 2015 "
            "GM_Earth value 398600.435436 km^3/s^2"
        )

    def test_earth_moon_mass_ratio_in_source_code(self):
        """Verify that lunar.py uses the IAU 2015 Earth/Moon mass ratio."""
        lunar_py = Path(__file__).parent.parent.parent / "libephemeris" / "lunar.py"
        content = lunar_py.read_text()

        assert "81.3005691" in content, (
            "lunar.py should derive GM_EARTH_MOON_AU3_DAY2 from the IAU 2015 "
            "Earth/Moon mass ratio 81.3005691"
        )

    def test_uses_combined_earth_moon_mu(self):
        """The shared constant equals the combined Earth-Moon mu in AU^3/day^2."""
        from libephemeris.lunar import GM_EARTH_MOON_AU3_DAY2

        expected = GM_EARTH_MOON_IAU2015_KM3_S2 / (AU_KM**3) * (SECONDS_PER_DAY**2)
        assert math.isclose(GM_EARTH_MOON_AU3_DAY2, expected, rel_tol=1e-12), (
            f"GM_EARTH_MOON_AU3_DAY2 = {GM_EARTH_MOON_AU3_DAY2} does not equal "
            f"the combined IAU 2015 Earth-Moon mu {expected} AU^3/day^2"
        )
        # And the osculating-element code consumes the shared constant
        lunar_py = Path(__file__).parent.parent.parent / "libephemeris" / "lunar.py"
        content = lunar_py.read_text()
        pattern = r"mu\s*=\s*GM_EARTH_MOON_AU3_DAY2"
        assert re.search(pattern, content), (
            "osculating-element mu should come from the shared "
            "GM_EARTH_MOON_AU3_DAY2 constant"
        )

    def test_gm_value_documented_in_docstring(self):
        """Verify that the GM value is documented in the docstring."""
        lunar_py = Path(__file__).parent.parent.parent / "libephemeris" / "lunar.py"
        content = lunar_py.read_text()

        # Check that the IAU 2015 value is mentioned in docstring
        assert "398600.435436" in content, (
            "GM_Earth value 398600.435436 should be documented in lunar.py"
        )

        # Check that IAU 2015 is referenced
        assert "IAU 2015" in content, (
            "IAU 2015 Resolution B3 should be referenced in lunar.py"
        )

    def test_gm_conversion_to_au_day(self):
        """Test the conversion of GM from km^3/s^2 to AU^3/day^2."""
        # Convert from km^3/s^2 to AU^3/day^2
        # Formula: GM (AU^3/day^2) = GM (km^3/s^2) / (AU_km^3) * (sec_per_day^2)
        gm_au3_day2 = GM_EARTH_IAU2015_KM3_S2 / (AU_KM**3) * (SECONDS_PER_DAY**2)

        # Expected value approximately 8.887e-10 AU^3/day^2
        # This is used in orbital mechanics for Moon's orbit
        assert gm_au3_day2 > 0, "GM in AU^3/day^2 must be positive"
        assert 8.8e-10 < gm_au3_day2 < 9.0e-10, (
            f"GM conversion result {gm_au3_day2} is outside expected range"
        )

    def test_gm_relative_precision(self):
        """
        Test that the IAU 2015 value provides sufficient precision.

        The relative difference between IAU 2015 and older values (e.g., WGS-84)
        is small but significant for precise orbital calculations.
        """
        # Old WGS-84/EGM96 value (circa 1996)
        gm_old = 398600.4418

        # Relative difference
        rel_diff = abs(GM_EARTH_IAU2015_KM3_S2 - gm_old) / GM_EARTH_IAU2015_KM3_S2

        # The difference should be on the order of 1e-8 (about 16 ppb)
        assert rel_diff < 1e-6, "Relative difference is unexpectedly large"
        assert rel_diff > 1e-10, "Values should not be exactly equal"

        # The new value is slightly smaller (better satellite tracking data)
        assert GM_EARTH_IAU2015_KM3_S2 < gm_old, (
            "IAU 2015 value should be smaller than old WGS-84 value"
        )


class TestGMImpactOnOrbitalCalculations:
    """Tests for the impact of GM value on orbital calculations."""

    def test_orbital_velocity_formula(self):
        """
        Test that circular orbital velocity formula works with the GM value.

        v = sqrt(GM/r) where v is velocity, GM is gravitational parameter, r is radius

        For two-body problem, use combined Earth-Moon GM.
        """
        # Mean Moon distance in km
        moon_distance_km = 384400.0

        # Calculate circular orbital velocity using combined Earth-Moon GM
        v_km_s = math.sqrt(GM_EARTH_MOON_IAU2015_KM3_S2 / moon_distance_km)

        # Moon's actual mean orbital velocity is about 1.022-1.024 km/s
        assert 1.0 < v_km_s < 1.1, (
            f"Calculated velocity {v_km_s} km/s is outside expected range"
        )

    def test_orbital_period_formula(self):
        """
        Test that orbital period formula works with the GM value.

        T = 2*pi*sqrt(a^3/GM) where T is period, a is semi-major axis

        For two-body problem, use combined Earth-Moon GM.
        """
        # Moon's mean semi-major axis in km
        moon_sma_km = 384748.0

        # Calculate orbital period in seconds using combined Earth-Moon GM
        T_seconds = (
            2 * math.pi * math.sqrt(moon_sma_km**3 / GM_EARTH_MOON_IAU2015_KM3_S2)
        )

        # Convert to days
        T_days = T_seconds / SECONDS_PER_DAY

        # Moon's sidereal orbital period is about 27.32 days
        assert 27.0 < T_days < 28.0, (
            f"Calculated period {T_days} days is outside expected range"
        )


class TestEarthMoonCombinedGM:
    """Tests for the combined Earth-Moon gravitational parameter."""

    def test_combined_gm_is_larger_than_earth_only(self):
        """Verify that combined GM is larger than Earth-only GM."""
        assert GM_EARTH_MOON_IAU2015_KM3_S2 > GM_EARTH_IAU2015_KM3_S2, (
            "Combined Earth-Moon GM should be larger than Earth-only GM"
        )

    def test_moon_contribution_is_about_1_percent(self):
        """Verify the Moon's contribution to the combined GM is ~1.2%."""
        moon_fraction = GM_MOON_IAU2015_KM3_S2 / GM_EARTH_IAU2015_KM3_S2

        # Moon is about 1/81.3 of Earth's mass, so ~1.23% contribution
        assert 0.012 < moon_fraction < 0.013, (
            f"Moon's GM fraction {moon_fraction} is outside expected range (0.012-0.013)"
        )

    def test_combined_gm_value_matches_expected(self):
        """Verify the combined GM matches the expected value."""
        # Expected combined GM from IAU 2015 values
        # GM_Earth = 398600.435436 km³/s²
        # GM_Moon = 398600.435436 / 81.3005691 = 4902.80008 km³/s²
        # Total = 403503.2355 km³/s²
        expected_combined = 403503.2355

        assert abs(GM_EARTH_MOON_IAU2015_KM3_S2 - expected_combined) < 0.1, (
            f"Combined GM {GM_EARTH_MOON_IAU2015_KM3_S2} does not match expected "
            f"value {expected_combined} km^3/s^2"
        )

    def test_using_earth_only_gm_gives_incorrect_period(self):
        """
        Demonstrate that using Earth-only GM gives incorrect orbital period.

        This test shows why we need the combined Earth-Moon GM for the
        eccentricity vector calculation.
        """
        # Moon's mean semi-major axis in km
        moon_sma_km = 384748.0

        # Calculate period with Earth-only GM (incorrect)
        T_earth_only = 2 * math.pi * math.sqrt(moon_sma_km**3 / GM_EARTH_IAU2015_KM3_S2)
        T_earth_only_days = T_earth_only / SECONDS_PER_DAY

        # Calculate period with combined GM (correct)
        T_combined = (
            2 * math.pi * math.sqrt(moon_sma_km**3 / GM_EARTH_MOON_IAU2015_KM3_S2)
        )
        T_combined_days = T_combined / SECONDS_PER_DAY

        # The difference is about 0.2 days (significant!)
        difference_days = T_earth_only_days - T_combined_days

        # Using Earth-only GM gives a period that's too long by ~0.16 days
        assert 0.1 < difference_days < 0.3, (
            f"Period difference {difference_days} days is outside expected range"
        )

    def test_combined_gm_conversion_to_au_day(self):
        """Test the conversion of combined GM from km^3/s^2 to AU^3/day^2."""
        # Convert from km^3/s^2 to AU^3/day^2
        gm_au3_day2 = GM_EARTH_MOON_IAU2015_KM3_S2 / (AU_KM**3) * (SECONDS_PER_DAY**2)

        # Expected value approximately 8.997e-10 AU^3/day^2 (slightly larger than Earth-only)
        assert gm_au3_day2 > 0, "GM in AU^3/day^2 must be positive"
        assert 8.9e-10 < gm_au3_day2 < 9.2e-10, (
            f"Combined GM conversion result {gm_au3_day2} is outside expected range"
        )
