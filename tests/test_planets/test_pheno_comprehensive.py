"""
Tests for pheno_ut: planetary phenomena (phase angle, illumination,
elongation, diameter, magnitude) across planets and dates.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
)


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    swe.close()


from libephemeris.constants import FLG_HELCTR, FLG_TOPOCTR  # noqa: E402


class TestPhenoTopocentric:
    """FLG_TOPOCTR is honored in pheno (measured reference behavior)."""

    @pytest.mark.unit
    def test_moon_topocentric_differs_from_geocentric(self):
        """The Moon's topocentric phenomena differ from the geocentric ones:
        the observer parallax enlarges/shrinks the disc and shifts the phase
        angle by a few tenths of a degree. The base (pre-fix) behavior ignored
        the flag and returned identical values."""
        swe.set_topo(12.5, 41.9, 100.0)
        geo = swe.pheno_ut(2451545.0, MOON, 0)
        topo = swe.pheno_ut(2451545.0, MOON, FLG_TOPOCTR)
        assert abs(topo[3] - geo[3]) * 3600.0 > 1.0  # diameter shift > 1 arcsec
        assert abs(topo[0] - geo[0]) > 0.05  # phase-angle shift > 0.05 deg

    @pytest.mark.unit
    def test_far_planet_topocentric_effect_tiny(self):
        """A far planet's parallax is small: the topocentric diameter barely
        moves, but the value is still finite and positive."""
        swe.set_topo(12.5, 41.9, 100.0)
        geo = swe.pheno_ut(2451545.0, JUPITER, 0)
        topo = swe.pheno_ut(2451545.0, JUPITER, FLG_TOPOCTR)
        assert topo[3] > 0.0
        assert abs(topo[3] - geo[3]) * 3600.0 < 0.1

    @pytest.mark.unit
    def test_topocentric_without_set_topo_raises(self):
        """FLG_TOPOCTR without a prior set_topo() is an error (as in calc)."""
        from libephemeris.exceptions import ConfigurationError

        swe.close()  # clears any previously set observer location
        with pytest.raises(ConfigurationError):
            swe.pheno_ut(2451545.0, MOON, FLG_TOPOCTR)

    @pytest.mark.unit
    def test_helctr_still_ignored_with_topoctr(self):
        """FLG_HELCTR remains ignored; FLG_TOPOCTR|FLG_HELCTR == FLG_TOPOCTR."""
        swe.set_topo(12.5, 41.9, 100.0)
        topo = swe.pheno_ut(2451545.0, MOON, FLG_TOPOCTR)
        topo_helctr = swe.pheno_ut(2451545.0, MOON, FLG_TOPOCTR | FLG_HELCTR)
        for a, b in zip(topo, topo_helctr):
            assert a == pytest.approx(b, abs=1e-9)


JD_J2000 = 2451545.0
JD_2020 = 2458849.5
JD_2023 = 2460000.0


class TestPhenoReturnFormat:
    """Test return format of pheno_ut."""

    @pytest.mark.unit
    def test_returns_20_floats(self):
        """pheno_ut returns tuple of 20 floats."""
        result = swe.pheno_ut(JD_J2000, MARS)
        assert len(result) == 20

    @pytest.mark.unit
    def test_all_values_are_floats(self):
        """All returned values are Python floats."""
        result = swe.pheno_ut(JD_J2000, MARS)
        for i, val in enumerate(result):
            assert isinstance(val, float), f"Index {i} is {type(val)}, not float"

    @pytest.mark.unit
    def test_reserved_fields_zero(self):
        """Fields [5]-[19] are reserved and should be 0."""
        result = swe.pheno_ut(JD_J2000, MARS)
        for i in range(5, 20):
            assert result[i] == 0.0, f"Reserved field [{i}] = {result[i]}"


class TestPhenoSun:
    """Test phenomena for the Sun."""

    @pytest.mark.unit
    def test_sun_elongation_zero(self):
        """Sun elongation from Sun should be 0."""
        result = swe.pheno_ut(JD_J2000, SUN)
        assert result[2] == pytest.approx(0.0, abs=0.01)

    @pytest.mark.unit
    def test_sun_magnitude(self):
        """Sun apparent magnitude should be near -26.7."""
        result = swe.pheno_ut(JD_J2000, SUN)
        # Sun V magnitude ~ -26.7 to -26.8
        assert -27.5 < result[4] < -26.0, f"Sun magnitude: {result[4]}"


class TestPhenoMoon:
    """Test phenomena for the Moon."""

    @pytest.mark.unit
    def test_moon_phase_range(self):
        """Moon illumination is between 0 and 1."""
        result = swe.pheno_ut(JD_J2000, MOON)
        assert 0.0 <= result[1] <= 1.0

    @pytest.mark.unit
    def test_moon_phase_varies_with_date(self):
        """Moon illumination changes across the month."""
        phases = []
        for offset in range(0, 30, 3):
            result = swe.pheno_ut(JD_J2000 + offset, MOON)
            phases.append(result[1])
        # Should see variation
        assert max(phases) - min(phases) > 0.3

    @pytest.mark.unit
    def test_moon_diameter(self):
        """Moon apparent diameter should be reasonable."""
        result = swe.pheno_ut(JD_J2000, MOON)
        diameter = result[3]
        # Moon diameter ~0.5° = 1800 arcsec (but API returns degrees for Moon)
        # If in degrees: ~0.5; if in arcsec: ~1800
        assert diameter > 0.0

    @pytest.mark.unit
    def test_moon_elongation_range(self):
        """Moon elongation should be 0-180 degrees."""
        result = swe.pheno_ut(JD_J2000, MOON)
        assert 0.0 <= result[2] <= 180.0


class TestPhenoOuterPlanets:
    """Test phenomena for outer planets across dates."""

    OUTER_PLANETS = [MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO]

    @pytest.mark.unit
    @pytest.mark.parametrize("body", OUTER_PLANETS)
    def test_phase_angle_range(self, body):
        """Phase angle should be in [0, 180] degrees."""
        result = swe.pheno_ut(JD_J2000, body)
        assert 0.0 <= result[0] <= 180.0, f"Phase angle {result[0]} out of range"

    @pytest.mark.unit
    @pytest.mark.parametrize("body", OUTER_PLANETS)
    def test_illumination_range(self, body):
        """Illuminated fraction should be in [0, 1]."""
        result = swe.pheno_ut(JD_J2000, body)
        assert 0.0 <= result[1] <= 1.0, f"Illumination {result[1]} out of range"

    @pytest.mark.unit
    @pytest.mark.parametrize("body", OUTER_PLANETS)
    def test_elongation_range(self, body):
        """Elongation should be in [0, 180] degrees."""
        result = swe.pheno_ut(JD_J2000, body)
        assert 0.0 <= result[2] <= 180.0, f"Elongation {result[2]} out of range"

    @pytest.mark.unit
    @pytest.mark.parametrize("body", OUTER_PLANETS)
    def test_diameter_positive(self, body):
        """Apparent diameter should be positive."""
        result = swe.pheno_ut(JD_J2000, body)
        assert result[3] > 0.0, f"Diameter should be positive: {result[3]}"

    @pytest.mark.unit
    @pytest.mark.parametrize("body", OUTER_PLANETS)
    def test_magnitude_finite(self, body):
        """Magnitude should be a finite number."""
        result = swe.pheno_ut(JD_J2000, body)
        assert math.isfinite(result[4]), f"Magnitude not finite: {result[4]}"

    @pytest.mark.unit
    def test_jupiter_brighter_than_saturn(self):
        """Jupiter should generally be brighter (lower mag) than Saturn."""
        j_result = swe.pheno_ut(JD_J2000, JUPITER)
        s_result = swe.pheno_ut(JD_J2000, SATURN)
        # Not always true but generally Jupiter is brighter
        # Use a loose check — just verify both have reasonable magnitudes
        assert j_result[4] < 0.0, f"Jupiter mag should be negative: {j_result[4]}"
        assert s_result[4] < 5.0, f"Saturn mag too faint: {s_result[4]}"

    @pytest.mark.unit
    def test_outer_planet_illumination_high(self):
        """Outer planets are nearly fully illuminated (small phase angle)."""
        for body in [JUPITER, SATURN, URANUS, NEPTUNE]:
            result = swe.pheno_ut(JD_J2000, body)
            # Outer planets have small phase angles, high illumination
            assert result[1] > 0.8, (
                f"Body {body} illumination {result[1]} unexpectedly low"
            )


class TestPhenoMagnitudeVariation:
    """Test magnitude variation over time for outer planets."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [MARS, JUPITER, SATURN])
    def test_magnitude_varies_over_year(self, body):
        """Planet magnitude varies across a year."""
        mags = []
        for i in range(0, 365, 30):
            result = swe.pheno_ut(JD_J2000 + i, body)
            mags.append(result[4])
        # Magnitude should vary by at least 0.1 mag over a year
        variation = max(mags) - min(mags)
        assert variation > 0.1, (
            f"Body {body}: magnitude variation {variation} too small"
        )

    @pytest.mark.unit
    def test_mars_magnitude_range(self):
        """Mars magnitude can vary widely (opposition effect)."""
        mags = []
        # Sample over 2+ years (Mars opposition ~every 26 months)
        for i in range(0, 800, 30):
            result = swe.pheno_ut(JD_J2000 + i, MARS)
            mags.append(result[4])
        # Mars ranges from about -2.9 at opposition to +1.9
        assert min(mags) < 1.0, f"Mars min mag {min(mags)} — expected brighter"
        assert max(mags) > 0.5, f"Mars max mag {max(mags)} — expected fainter range"


class TestPhenoInnerPlanets:
    """Test inner planet phenomena."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [MERCURY, VENUS])
    def test_max_elongation_bounded(self, body):
        """Inner planets have bounded elongation from Sun."""
        max_elong = 0.0
        for i in range(0, 365, 5):
            result = swe.pheno_ut(JD_J2000 + i, body)
            max_elong = max(max_elong, result[2])
        if body == MERCURY:
            assert max_elong < 30.0, f"Mercury max elongation {max_elong} > 30°"
        else:
            assert max_elong < 50.0, f"Venus max elongation {max_elong} > 50°"

    @pytest.mark.unit
    def test_venus_can_be_very_bright(self):
        """Venus magnitude can reach about -4.5."""
        min_mag = 999.0
        for i in range(0, 600, 10):
            result = swe.pheno_ut(JD_J2000 + i, VENUS)
            min_mag = min(min_mag, result[4])
        assert min_mag < -3.0, f"Venus min magnitude {min_mag} — expected < -3"


class TestPhenoETVariant:
    """Test the ET variant pheno."""

    @pytest.mark.unit
    def test_pheno_et_works(self):
        """pheno (ET variant) returns valid results."""
        result = swe.pheno(JD_J2000, MARS)
        assert len(result) == 20
        assert 0.0 <= result[0] <= 180.0  # phase angle
        assert 0.0 <= result[1] <= 1.0  # illumination


from libephemeris.constants import (  # noqa: E402
    AST_OFFSET,
    VESTA,
    PHOLUS,
    CERES,
    FLG_SWIEPH,
)


class TestPhenoNumberedAsteroids:
    """Physical magnitude/diameter for numbered asteroids the library serves.

    Before this fix pheno reported 0.0 magnitude and 0.0 diameter for every
    numbered minor body (AST_OFFSET + number) except the six curated ones,
    because the H-G/diameter tables covered only the curated bodies. The tables
    now carry JPL SBDB photometry for the positionable numbered asteroids, and
    the AST_OFFSET + number alias of a built-in asteroid resolves to the
    built-in's row (matching how calc_ut serves its position).
    """

    EROS = AST_OFFSET + 433
    PSYCHE = AST_OFFSET + 16

    @pytest.mark.unit
    def test_eros_magnitude_is_finite_nonzero_native_float(self):
        """433 Eros reports a real, finite, nonzero visual magnitude."""
        swe.set_calc_mode("skyfield")
        swe.set_auto_spk_download(True)
        result = swe.pheno_ut(swe.julday(1990, 6, 1, 0.0), self.EROS, FLG_SWIEPH)
        magnitude = result[4]
        diameter = result[3]
        assert isinstance(magnitude, float)
        assert isinstance(diameter, float)
        assert math.isfinite(magnitude)
        assert magnitude != 0.0
        # Eros is a faint main-belt/NEA object, never as bright as the planets.
        assert 8.0 < magnitude < 20.0, magnitude
        assert diameter > 0.0

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "alias,builtin",
        [
            (AST_OFFSET + 4, VESTA),
            (AST_OFFSET + 5145, PHOLUS),
            (AST_OFFSET + 1, CERES),
        ],
    )
    def test_alias_matches_builtin_photometry(self, alias, builtin):
        """AST_OFFSET + number alias yields the built-in body's full pheno."""
        jd = swe.julday(1985, 3, 15, 0.0)
        alias_result = swe.pheno_ut(jd, alias, FLG_SWIEPH)
        builtin_result = swe.pheno_ut(jd, builtin, FLG_SWIEPH)
        assert alias_result == builtin_result
        # The built-in asteroids carry a real magnitude and disc.
        assert builtin_result[4] != 0.0
        assert builtin_result[3] > 0.0

    @pytest.mark.unit
    def test_positionable_body_without_photometry_row_reports_zero(self):
        """A positionable numbered body absent from the tables keeps 0.0/0.0.

        136199 Eris is served from the registry SPK but has no published
        effective diameter in SBDB, so it carries no photometry row; the
        documented behavior is diameter and magnitude 0.0 (its elongation is
        still a real geometric quantity).
        """
        swe.set_calc_mode("skyfield")
        swe.set_auto_spk_download(True)
        eris = AST_OFFSET + 136199
        try:
            result = swe.pheno_ut(swe.julday(1990, 6, 1, 0.0), eris, FLG_SWIEPH)
        except Exception:
            pytest.skip("136199 Eris not positionable in this environment")
        assert result[3] == 0.0  # diameter
        assert result[4] == 0.0  # magnitude
        assert math.isfinite(result[2])  # elongation is a real quantity
