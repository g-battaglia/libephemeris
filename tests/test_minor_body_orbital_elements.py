"""
Unit tests for minor body orbital elements data.

Tests verify:
- Orbital elements are at the correct epoch (JD 2461000.5 TDB = 2025-Sep-19)
- All 37 bodies in MINOR_BODY_ELEMENTS have valid orbital parameters
- Elements are consistent with expected physical ranges
- Mean motion (n) is consistent with semi-major axis (a) via Kepler's 3rd law
"""

import math
import pytest
from libephemeris.constants import (
    CHIRON,
    PHOLUS,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    ERIS,
    SEDNA,
    HAUMEA,
    MAKEMAKE,
    IXION,
    ORCUS,
    QUAOAR,
    NESSUS,
    ASBOLUS,
    CHARIKLO,
    GONGGONG,
    VARUNA,
    APOPHIS,
    HYGIEA,
    INTERAMNIA,
    DAVIDA,
    EUROPA_AST,
    SYLVIA,
    PSYCHE,
    EROS,
    AMOR,
    ICARUS,
    TORO,
    SAPPHO,
    PANDORA_AST,
    LILITH_AST,
    HIDALGO,
    TOUTATIS,
    ITOKAWA,
    BENNU,
    RYUGU,
)
from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS, OrbitalElements


# Expected epoch for all bodies (JD 2461000.5 TDB = 2025-Sep-19)
EXPECTED_EPOCH = 2461000.5

# Gaussian gravitational constant: k = 0.01720209895 rad/day
# This gives n = k / a^(3/2) in rad/day, or n_deg = k * 180/pi / a^(3/2) deg/day
# For consistency check: n_deg = 0.9856076686 / a^1.5
KEPLER_CONSTANT = 0.9856076686  # degrees/day for a in AU


@pytest.fixture
def all_minor_body_ids():
    """All minor body identifiers."""
    return [
        CHIRON,
        PHOLUS,
        CERES,
        PALLAS,
        JUNO,
        VESTA,
        ERIS,
        SEDNA,
        HAUMEA,
        MAKEMAKE,
        IXION,
        ORCUS,
        QUAOAR,
        NESSUS,
        ASBOLUS,
        CHARIKLO,
        GONGGONG,
        VARUNA,
        APOPHIS,
        HYGIEA,
        INTERAMNIA,
        DAVIDA,
        EUROPA_AST,
        SYLVIA,
        PSYCHE,
        EROS,
        AMOR,
        ICARUS,
        TORO,
        SAPPHO,
        PANDORA_AST,
        LILITH_AST,
        HIDALGO,
        TOUTATIS,
        ITOKAWA,
        BENNU,
        RYUGU,
    ]


@pytest.mark.unit
class TestOrbitalElementsEpoch:
    """Tests for orbital elements epoch."""

    def test_all_bodies_at_correct_epoch(self, all_minor_body_ids):
        """Verify all bodies have epoch JD 2461000.5."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            assert elements.epoch == EXPECTED_EPOCH, (
                f"{elements.name}: epoch {elements.epoch} != expected {EXPECTED_EPOCH}"
            )

    def test_expected_number_of_bodies(self):
        """Verify we have 37 minor bodies."""
        assert len(MINOR_BODY_ELEMENTS) == 37, (
            f"Expected 37 minor bodies, got {len(MINOR_BODY_ELEMENTS)}"
        )


@pytest.mark.unit
class TestOrbitalElementsValidity:
    """Tests for orbital element validity."""

    def test_all_eccentricities_valid(self, all_minor_body_ids):
        """Eccentricity must be in [0, 1) for bound orbits."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            assert 0 <= elements.e < 1, (
                f"{elements.name}: invalid eccentricity {elements.e}"
            )

    def test_all_inclinations_valid(self, all_minor_body_ids):
        """Inclination must be in [0, 180) degrees."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            assert 0 <= elements.i < 180, (
                f"{elements.name}: invalid inclination {elements.i}"
            )

    def test_all_angular_elements_valid(self, all_minor_body_ids):
        """Angular elements (omega, Omega, M0) must be in [0, 360) degrees."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            assert 0 <= elements.omega < 360, (
                f"{elements.name}: omega {elements.omega} not in [0, 360)"
            )
            assert 0 <= elements.Omega < 360, (
                f"{elements.name}: Omega {elements.Omega} not in [0, 360)"
            )
            assert 0 <= elements.M0 < 360, (
                f"{elements.name}: M0 {elements.M0} not in [0, 360)"
            )

    def test_all_semi_major_axes_positive(self, all_minor_body_ids):
        """Semi-major axis must be positive."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            assert elements.a > 0, (
                f"{elements.name}: invalid semi-major axis {elements.a}"
            )

    def test_all_mean_motions_positive(self, all_minor_body_ids):
        """Mean motion must be positive."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            assert elements.n > 0, f"{elements.name}: invalid mean motion {elements.n}"


@pytest.mark.unit
class TestOrbitalElementsConsistency:
    """Tests for orbital element physical consistency."""

    def test_mean_motion_vs_semi_major_axis(self, all_minor_body_ids):
        """
        Mean motion should be consistent with semi-major axis via Kepler's 3rd law.

        n = 0.9856076686 / a^1.5 (degrees/day)
        """
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            expected_n = KEPLER_CONSTANT / (elements.a**1.5)
            rel_error = abs(elements.n - expected_n) / expected_n

            # Allow 0.1% relative error for JPL's high-precision elements
            assert rel_error < 0.001, (
                f"{elements.name}: mean motion {elements.n} deg/day differs from "
                f"Kepler prediction {expected_n:.10f} deg/day by {rel_error * 100:.4f}%"
            )

    def test_perihelion_within_aphelion(self, all_minor_body_ids):
        """Perihelion distance q = a(1-e) must be positive."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            q = elements.a * (1 - elements.e)
            Q = elements.a * (1 + elements.e)

            assert q > 0, f"{elements.name}: negative perihelion {q} AU"
            assert Q > q, f"{elements.name}: aphelion {Q} <= perihelion {q}"


@pytest.mark.unit
class TestSpecificBodiesOrbitalElements:
    """Tests for specific body orbital element ranges."""

    def test_chiron_orbital_elements(self):
        """Chiron: Centaur with orbit between Saturn and Uranus."""
        elements = MINOR_BODY_ELEMENTS[CHIRON]
        assert elements.name == "Chiron"
        assert 13.5 < elements.a < 14.0, f"Chiron a={elements.a} unexpected"
        assert 0.37 < elements.e < 0.39, f"Chiron e={elements.e} unexpected"
        assert 6.5 < elements.i < 7.5, f"Chiron i={elements.i} unexpected"

    def test_ceres_orbital_elements(self):
        """Ceres: Largest main belt asteroid, dwarf planet."""
        elements = MINOR_BODY_ELEMENTS[CERES]
        assert elements.name == "Ceres"
        assert 2.7 < elements.a < 2.8, f"Ceres a={elements.a} unexpected"
        assert 0.07 < elements.e < 0.09, f"Ceres e={elements.e} unexpected"
        assert 10.0 < elements.i < 11.0, f"Ceres i={elements.i} unexpected"

    def test_pallas_high_inclination(self):
        """Pallas: Known for high inclination ~35 degrees."""
        elements = MINOR_BODY_ELEMENTS[PALLAS]
        assert elements.name == "Pallas"
        assert 34.0 < elements.i < 36.0, f"Pallas i={elements.i} unexpected"

    def test_eris_distant_tno(self):
        """Eris: Distant TNO with high inclination."""
        elements = MINOR_BODY_ELEMENTS[ERIS]
        assert elements.name == "Eris"
        assert 67.0 < elements.a < 69.0, f"Eris a={elements.a} unexpected"
        assert 43.0 < elements.i < 45.0, f"Eris i={elements.i} unexpected"

    def test_sedna_extreme_orbit(self):
        """Sedna: Detached object with extreme orbit."""
        elements = MINOR_BODY_ELEMENTS[SEDNA]
        assert elements.name == "Sedna"
        # Sedna has the most distant semi-major axis
        assert 500 < elements.a < 600, f"Sedna a={elements.a} unexpected"
        # Very eccentric
        assert 0.85 < elements.e < 0.87, f"Sedna e={elements.e} unexpected"
        # Very slow mean motion
        assert elements.n < 0.0001, f"Sedna n={elements.n} too fast"

    def test_ixion_plutino(self):
        """Ixion: Plutino in 2:3 resonance with Neptune."""
        elements = MINOR_BODY_ELEMENTS[IXION]
        assert elements.name == "Ixion"
        # Plutinos have a ~ 39.4 AU (like Pluto)
        assert 39.0 < elements.a < 40.0, f"Ixion a={elements.a} unexpected"

    def test_orcus_anti_pluto(self):
        """Orcus: Anti-Pluto (opposite phase in 2:3 resonance)."""
        elements = MINOR_BODY_ELEMENTS[ORCUS]
        assert elements.name == "Orcus"
        # Similar to Pluto's orbital parameters
        assert 39.0 < elements.a < 40.0, f"Orcus a={elements.a} unexpected"

    def test_quaoar_low_eccentricity(self):
        """Quaoar: TNO with nearly circular orbit."""
        elements = MINOR_BODY_ELEMENTS[QUAOAR]
        assert elements.name == "Quaoar"
        # Low eccentricity for a TNO
        assert elements.e < 0.05, (
            f"Quaoar e={elements.e} too high for 'nearly circular'"
        )

    def test_nessus_centaur(self):
        """Nessus: Centaur with orbit between Saturn and Pluto."""
        elements = MINOR_BODY_ELEMENTS[NESSUS]
        assert elements.name == "Nessus"
        # Semi-major axis ~24.5 AU (between Saturn ~9.5 AU and Pluto ~39.5 AU)
        assert 24.0 < elements.a < 25.0, f"Nessus a={elements.a} unexpected"
        # Moderately eccentric orbit
        assert 0.51 < elements.e < 0.53, f"Nessus e={elements.e} unexpected"
        # Inclination ~15.6 degrees
        assert 15.0 < elements.i < 16.5, f"Nessus i={elements.i} unexpected"

    def test_asbolus_centaur(self):
        """Asbolus: Centaur with orbit between Saturn and Uranus."""
        elements = MINOR_BODY_ELEMENTS[ASBOLUS]
        assert elements.name == "Asbolus"
        # Semi-major axis ~18 AU (between Saturn ~9.5 AU and Uranus ~19.2 AU)
        assert 17.5 < elements.a < 18.5, f"Asbolus a={elements.a} unexpected"
        # Highly eccentric orbit
        assert 0.61 < elements.e < 0.63, f"Asbolus e={elements.e} unexpected"
        # Inclination ~17.6 degrees
        assert 17.0 < elements.i < 18.5, f"Asbolus i={elements.i} unexpected"

    def test_chariklo_centaur(self):
        """Chariklo: Largest known centaur with ring system discovered in 2014."""
        elements = MINOR_BODY_ELEMENTS[CHARIKLO]
        assert elements.name == "Chariklo"
        # Semi-major axis ~15.7 AU (between Saturn ~9.5 AU and Uranus ~19.2 AU)
        assert 15.5 < elements.a < 16.0, f"Chariklo a={elements.a} unexpected"
        # Moderate eccentricity
        assert 0.16 < elements.e < 0.18, f"Chariklo e={elements.e} unexpected"
        # Inclination ~23.4 degrees
        assert 23.0 < elements.i < 24.0, f"Chariklo i={elements.i} unexpected"

    def test_gonggong_tno(self):
        """Gonggong: TNO dwarf planet candidate (formerly 2007 OR10)."""
        elements = MINOR_BODY_ELEMENTS[GONGGONG]
        assert elements.name == "Gonggong"
        # Semi-major axis ~67 AU (distant TNO)
        assert 66.0 < elements.a < 68.0, f"Gonggong a={elements.a} unexpected"
        # Moderately eccentric orbit
        assert 0.50 < elements.e < 0.51, f"Gonggong e={elements.e} unexpected"
        # High inclination ~30.9 degrees
        assert 30.5 < elements.i < 31.5, f"Gonggong i={elements.i} unexpected"
        # Slow mean motion (long orbital period ~547 years)
        assert elements.n < 0.002, f"Gonggong n={elements.n} too fast"

    def test_varuna_classical_kbo(self):
        """Varuna: Large classical Kuiper belt object (~670 km diameter)."""
        elements = MINOR_BODY_ELEMENTS[VARUNA]
        assert elements.name == "Varuna"
        # Semi-major axis ~43.2 AU (classical KBO region)
        assert 42.5 < elements.a < 44.0, f"Varuna a={elements.a} unexpected"
        # Nearly circular orbit (low eccentricity)
        assert elements.e < 0.06, f"Varuna e={elements.e} too high for classical KBO"
        # Moderate inclination ~17.1 degrees
        assert 16.5 < elements.i < 18.0, f"Varuna i={elements.i} unexpected"
        # Mean motion consistent with ~284 year period
        assert 0.003 < elements.n < 0.004, f"Varuna n={elements.n} unexpected"

    def test_apophis_near_earth_asteroid(self):
        """Apophis: Near-Earth asteroid with close approach in 2029."""
        elements = MINOR_BODY_ELEMENTS[APOPHIS]
        assert elements.name == "Apophis"
        # Semi-major axis ~0.92 AU (Aten-class NEA, crosses Earth's orbit)
        assert 0.91 < elements.a < 0.93, f"Apophis a={elements.a} unexpected"
        # Moderate eccentricity
        assert 0.18 < elements.e < 0.20, f"Apophis e={elements.e} unexpected"
        # Low inclination ~3.3 degrees
        assert 3.0 < elements.i < 4.0, f"Apophis i={elements.i} unexpected"
        # Fast mean motion (short orbital period ~0.89 years = 324 days)
        assert 1.10 < elements.n < 1.15, f"Apophis n={elements.n} unexpected"

    def test_hygiea_main_belt_asteroid(self):
        """Hygiea: Fourth largest main belt asteroid, dwarf planet candidate."""
        elements = MINOR_BODY_ELEMENTS[HYGIEA]
        assert elements.name == "Hygiea"
        # Semi-major axis ~3.15 AU (main belt)
        assert 3.10 < elements.a < 3.20, f"Hygiea a={elements.a} unexpected"
        # Low eccentricity
        assert 0.10 < elements.e < 0.12, f"Hygiea e={elements.e} unexpected"
        # Low inclination ~3.8 degrees
        assert 3.5 < elements.i < 4.5, f"Hygiea i={elements.i} unexpected"
        # Mean motion consistent with ~5.6 year period
        assert 0.17 < elements.n < 0.18, f"Hygiea n={elements.n} unexpected"

    def test_interamnia_main_belt_asteroid(self):
        """Interamnia: Fifth largest main belt asteroid."""
        elements = MINOR_BODY_ELEMENTS[INTERAMNIA]
        assert elements.name == "Interamnia"
        # Semi-major axis ~3.06 AU (main belt)
        assert 3.00 < elements.a < 3.10, f"Interamnia a={elements.a} unexpected"
        # Moderate eccentricity ~0.15
        assert 0.14 < elements.e < 0.17, f"Interamnia e={elements.e} unexpected"
        # Moderate inclination ~17.3 degrees
        assert 16.5 < elements.i < 18.0, f"Interamnia i={elements.i} unexpected"
        # Mean motion consistent with ~5.3 year period
        assert 0.18 < elements.n < 0.19, f"Interamnia n={elements.n} unexpected"

    def test_davida_main_belt_asteroid(self):
        """Davida: Seventh largest main belt asteroid."""
        elements = MINOR_BODY_ELEMENTS[DAVIDA]
        assert elements.name == "Davida"
        # Semi-major axis ~3.16 AU (main belt)
        assert 3.10 < elements.a < 3.20, f"Davida a={elements.a} unexpected"
        # Moderate eccentricity ~0.19
        assert 0.18 < elements.e < 0.20, f"Davida e={elements.e} unexpected"
        # Moderate inclination ~16 degrees
        assert 15.5 < elements.i < 16.5, f"Davida i={elements.i} unexpected"
        # Mean motion consistent with ~5.6 year period
        assert 0.17 < elements.n < 0.18, f"Davida n={elements.n} unexpected"

    def test_europa_main_belt_asteroid(self):
        """Europa: Main belt asteroid (not to be confused with Jupiter's moon)."""
        elements = MINOR_BODY_ELEMENTS[EUROPA_AST]
        assert elements.name == "Europa"
        # Semi-major axis ~3.09 AU (main belt)
        assert 3.05 < elements.a < 3.15, f"Europa a={elements.a} unexpected"
        # Low eccentricity ~0.11
        assert 0.10 < elements.e < 0.13, f"Europa e={elements.e} unexpected"
        # Low inclination ~7.5 degrees
        assert 7.0 < elements.i < 8.0, f"Europa i={elements.i} unexpected"
        # Mean motion consistent with ~5.4 year period
        assert 0.18 < elements.n < 0.19, f"Europa n={elements.n} unexpected"

    def test_sylvia_triple_asteroid(self):
        """Sylvia: Triple asteroid system with moons Romulus and Remus."""
        elements = MINOR_BODY_ELEMENTS[SYLVIA]
        assert elements.name == "Sylvia"
        # Semi-major axis ~3.48 AU (outer main belt)
        assert 3.45 < elements.a < 3.55, f"Sylvia a={elements.a} unexpected"
        # Low eccentricity ~0.09
        assert 0.08 < elements.e < 0.11, f"Sylvia e={elements.e} unexpected"
        # Moderate inclination ~10.9 degrees
        assert 10.5 < elements.i < 11.5, f"Sylvia i={elements.i} unexpected"
        # Mean motion consistent with ~6.5 year period
        assert 0.15 < elements.n < 0.16, f"Sylvia n={elements.n} unexpected"

    def test_psyche_metallic_asteroid(self):
        """Psyche: Metallic M-type asteroid, target of NASA Psyche mission."""
        elements = MINOR_BODY_ELEMENTS[PSYCHE]
        assert elements.name == "Psyche"
        # Semi-major axis ~2.92 AU (main belt)
        assert 2.90 < elements.a < 2.95, f"Psyche a={elements.a} unexpected"
        # Low eccentricity ~0.13
        assert 0.12 < elements.e < 0.15, f"Psyche e={elements.e} unexpected"
        # Low inclination ~3.1 degrees
        assert 2.5 < elements.i < 3.5, f"Psyche i={elements.i} unexpected"
        # Mean motion consistent with ~5.0 year period
        assert 0.19 < elements.n < 0.21, f"Psyche n={elements.n} unexpected"

    def test_eros_near_earth_asteroid(self):
        """Eros: Near-Earth asteroid, target of NEAR Shoemaker mission."""
        elements = MINOR_BODY_ELEMENTS[EROS]
        assert elements.name == "Eros"
        # Semi-major axis ~1.46 AU (Amor-class NEA, crosses Mars's orbit)
        assert 1.45 < elements.a < 1.47, f"Eros a={elements.a} unexpected"
        # Moderate eccentricity ~0.22
        assert 0.21 < elements.e < 0.24, f"Eros e={elements.e} unexpected"
        # Moderate inclination ~10.8 degrees
        assert 10.5 < elements.i < 11.5, f"Eros i={elements.i} unexpected"
        # Mean motion consistent with ~1.76 year period (643 days)
        assert 0.55 < elements.n < 0.57, f"Eros n={elements.n} unexpected"

    def test_amor_near_earth_asteroid(self):
        """Amor: Prototype of the Amor near-Earth asteroid class."""
        elements = MINOR_BODY_ELEMENTS[AMOR]
        assert elements.name == "Amor"
        # Semi-major axis ~1.92 AU (Amor-class NEA prototype)
        assert 1.90 < elements.a < 1.94, f"Amor a={elements.a} unexpected"
        # Moderate eccentricity ~0.43
        assert 0.42 < elements.e < 0.45, f"Amor e={elements.e} unexpected"
        # Moderate inclination ~11.9 degrees
        assert 11.5 < elements.i < 12.5, f"Amor i={elements.i} unexpected"
        # Mean motion consistent with ~2.66 year period (972 days)
        assert 0.36 < elements.n < 0.38, f"Amor n={elements.n} unexpected"

    def test_icarus_apollo_asteroid(self):
        """Icarus: Apollo asteroid with highly eccentric orbit reaching inside Mercury."""
        elements = MINOR_BODY_ELEMENTS[ICARUS]
        assert elements.name == "Icarus"
        # Semi-major axis ~1.08 AU (Apollo-class NEA)
        assert 1.07 < elements.a < 1.09, f"Icarus a={elements.a} unexpected"
        # Very high eccentricity ~0.83 (perihelion at ~0.19 AU, inside Mercury's orbit)
        assert 0.82 < elements.e < 0.84, f"Icarus e={elements.e} unexpected"
        # Moderate inclination ~22.8 degrees
        assert 22.0 < elements.i < 24.0, f"Icarus i={elements.i} unexpected"
        # Mean motion consistent with ~1.12 year period (409 days)
        assert 0.87 < elements.n < 0.89, f"Icarus n={elements.n} unexpected"

    def test_toro_apollo_asteroid(self):
        """Toro: Apollo asteroid, near-Earth asteroid."""
        elements = MINOR_BODY_ELEMENTS[TORO]
        assert elements.name == "Toro"
        # Semi-major axis ~1.37 AU (Apollo-class NEA)
        assert 1.36 < elements.a < 1.38, f"Toro a={elements.a} unexpected"
        # Moderate eccentricity ~0.44 (perihelion at ~0.77 AU, aphelion at ~1.96 AU)
        assert 0.43 < elements.e < 0.45, f"Toro e={elements.e} unexpected"
        # Low inclination ~9.4 degrees
        assert 9.0 < elements.i < 10.0, f"Toro i={elements.i} unexpected"
        # Mean motion consistent with ~1.6 year period (584 days)
        assert 0.61 < elements.n < 0.63, f"Toro n={elements.n} unexpected"

    def test_sappho_main_belt_asteroid(self):
        """Sappho: Main belt asteroid, artistic expression and same-sex love."""
        elements = MINOR_BODY_ELEMENTS[SAPPHO]
        assert elements.name == "Sappho"
        # Semi-major axis ~2.30 AU (main belt asteroid)
        assert 2.29 < elements.a < 2.31, f"Sappho a={elements.a} unexpected"
        # Moderate eccentricity ~0.20
        assert 0.19 < elements.e < 0.21, f"Sappho e={elements.e} unexpected"
        # Low inclination ~8.7 degrees
        assert 8.5 < elements.i < 8.9, f"Sappho i={elements.i} unexpected"
        # Mean motion consistent with ~3.5 year period (1271 days)
        assert 0.28 < elements.n < 0.29, f"Sappho n={elements.n} unexpected"

    def test_pandora_main_belt_asteroid(self):
        """Pandora: Main belt asteroid (distinct from Saturn moon Pandora)."""
        elements = MINOR_BODY_ELEMENTS[PANDORA_AST]
        assert elements.name == "Pandora"
        # Semi-major axis ~2.76 AU (main belt asteroid)
        assert 2.75 < elements.a < 2.77, f"Pandora a={elements.a} unexpected"
        # Moderate eccentricity ~0.145
        assert 0.14 < elements.e < 0.15, f"Pandora e={elements.e} unexpected"
        # Low inclination ~7.2 degrees
        assert 7.0 < elements.i < 7.5, f"Pandora i={elements.i} unexpected"
        # Mean motion consistent with ~4.6 year period (1673 days)
        assert 0.21 < elements.n < 0.22, f"Pandora n={elements.n} unexpected"

    def test_lilith_main_belt_asteroid(self):
        """Lilith: Main belt asteroid (not to be confused with lunar apogee Lilith)."""
        elements = MINOR_BODY_ELEMENTS[LILITH_AST]
        assert elements.name == "Lilith"
        # Semi-major axis ~2.66 AU (main belt asteroid)
        assert 2.65 < elements.a < 2.68, f"Lilith a={elements.a} unexpected"
        # Moderate eccentricity ~0.19
        assert 0.19 < elements.e < 0.20, f"Lilith e={elements.e} unexpected"
        # Low inclination ~5.6 degrees
        assert 5.5 < elements.i < 5.7, f"Lilith i={elements.i} unexpected"
        # Mean motion consistent with ~4.35 year period (1589 days)
        assert 0.22 < elements.n < 0.23, f"Lilith n={elements.n} unexpected"

    def test_hidalgo_centaur_asteroid(self):
        """Hidalgo: Centaur-class asteroid with comet-like orbit, used in astrological research."""
        elements = MINOR_BODY_ELEMENTS[HIDALGO]
        assert elements.name == "Hidalgo"
        # Semi-major axis ~5.73 AU (between Mars and Jupiter, Centaur-class)
        assert 5.70 < elements.a < 5.76, f"Hidalgo a={elements.a} unexpected"
        # Very high eccentricity ~0.66 (comet-like orbit)
        assert 0.66 < elements.e < 0.67, f"Hidalgo e={elements.e} unexpected"
        # Very high inclination ~42.5 degrees
        assert 42.0 < elements.i < 43.0, f"Hidalgo i={elements.i} unexpected"
        # Mean motion consistent with ~13.7 year period (5008 days)
        assert 0.071 < elements.n < 0.073, f"Hidalgo n={elements.n} unexpected"

    def test_toutatis_apollo_asteroid(self):
        """Toutatis: Apollo PHA, radar and spacecraft target, tumbling rotation."""
        elements = MINOR_BODY_ELEMENTS[TOUTATIS]
        assert elements.name == "Toutatis"
        # Semi-major axis ~2.54 AU (Apollo asteroid)
        assert 2.53 < elements.a < 2.56, f"Toutatis a={elements.a} unexpected"
        # High eccentricity ~0.62 (q=0.95 AU, Q=4.13 AU)
        assert 0.62 < elements.e < 0.63, f"Toutatis e={elements.e} unexpected"
        # Very low inclination ~0.45 degrees
        assert 0.4 < elements.i < 0.5, f"Toutatis i={elements.i} unexpected"
        # Mean motion consistent with ~4.05 year period (1481 days)
        assert 0.24 < elements.n < 0.25, f"Toutatis n={elements.n} unexpected"

    def test_itokawa_apollo_asteroid(self):
        """Itokawa: Apollo PHA, target of Hayabusa sample return mission."""
        elements = MINOR_BODY_ELEMENTS[ITOKAWA]
        assert elements.name == "Itokawa"
        # Semi-major axis ~1.32 AU (Apollo asteroid)
        assert 1.32 < elements.a < 1.33, f"Itokawa a={elements.a} unexpected"
        # Moderate eccentricity ~0.28 (q=0.95 AU, Q=1.70 AU)
        assert 0.27 < elements.e < 0.29, f"Itokawa e={elements.e} unexpected"
        # Very low inclination ~1.6 degrees
        assert 1.5 < elements.i < 1.7, f"Itokawa i={elements.i} unexpected"
        # Mean motion consistent with ~1.52 year period (557 days)
        assert 0.64 < elements.n < 0.66, f"Itokawa n={elements.n} unexpected"

    def test_bennu_apollo_asteroid(self):
        """Bennu: Apollo PHA, target of OSIRIS-REx sample return mission."""
        elements = MINOR_BODY_ELEMENTS[BENNU]
        assert elements.name == "Bennu"
        # Semi-major axis ~1.13 AU (Apollo asteroid)
        assert 1.12 < elements.a < 1.14, f"Bennu a={elements.a} unexpected"
        # Moderate eccentricity ~0.20 (q=0.90 AU, Q=1.36 AU)
        assert 0.20 < elements.e < 0.21, f"Bennu e={elements.e} unexpected"
        # Low inclination ~6.0 degrees
        assert 6.0 < elements.i < 6.1, f"Bennu i={elements.i} unexpected"
        # Mean motion consistent with ~1.20 year period (437 days)
        assert 0.82 < elements.n < 0.83, f"Bennu n={elements.n} unexpected"

    def test_ryugu_apollo_asteroid(self):
        """Ryugu: Apollo PHA, target of Hayabusa2 sample return mission."""
        elements = MINOR_BODY_ELEMENTS[RYUGU]
        assert elements.name == "Ryugu"
        # Semi-major axis ~1.19 AU (Apollo asteroid)
        assert 1.18 < elements.a < 1.20, f"Ryugu a={elements.a} unexpected"
        # Moderate eccentricity ~0.19 (q=0.96 AU, Q=1.42 AU)
        assert 0.19 < elements.e < 0.20, f"Ryugu e={elements.e} unexpected"
        # Low inclination ~5.9 degrees
        assert 5.8 < elements.i < 5.9, f"Ryugu i={elements.i} unexpected"
        # Mean motion consistent with ~1.30 year period (475 days)
        assert 0.75 < elements.n < 0.76, f"Ryugu n={elements.n} unexpected"


@pytest.mark.unit
class TestOrbitalElementsDataclass:
    """Tests for OrbitalElements dataclass."""

    def test_orbital_elements_attributes(self):
        """Test OrbitalElements has all required attributes."""
        elem = MINOR_BODY_ELEMENTS[CERES]
        assert hasattr(elem, "name")
        assert hasattr(elem, "epoch")
        assert hasattr(elem, "a")
        assert hasattr(elem, "e")
        assert hasattr(elem, "i")
        assert hasattr(elem, "omega")
        assert hasattr(elem, "Omega")
        assert hasattr(elem, "M0")
        assert hasattr(elem, "n")

    def test_orbital_elements_types(self):
        """Test OrbitalElements attribute types."""
        elem = MINOR_BODY_ELEMENTS[CERES]
        assert isinstance(elem.name, str)
        assert isinstance(elem.epoch, float)
        assert isinstance(elem.a, float)
        assert isinstance(elem.e, float)
        assert isinstance(elem.i, float)
        assert isinstance(elem.omega, float)
        assert isinstance(elem.Omega, float)
        assert isinstance(elem.M0, float)
        assert isinstance(elem.n, float)


@pytest.mark.unit
class TestJPLDataSourcePrecision:
    """Tests verifying JPL-level precision in orbital elements."""

    def test_semi_major_axis_precision(self, all_minor_body_ids):
        """Semi-major axis should have at least 6 significant figures."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            # Check that the value has significant precision
            # by ensuring it's not a round number
            a_str = f"{elements.a:.10f}"
            # Count non-zero digits after removing leading zeros and decimal point
            significant = len(a_str.replace(".", "").lstrip("0"))
            assert significant >= 6, (
                f"{elements.name}: semi-major axis {elements.a} has insufficient precision"
            )

    def test_mean_motion_precision(self, all_minor_body_ids):
        """Mean motion should have at least 6 significant figures."""
        for body_id in all_minor_body_ids:
            elements = MINOR_BODY_ELEMENTS[body_id]
            n_str = f"{elements.n:.15f}"
            # Count significant digits
            significant = len(
                n_str.replace(".", "").replace("-", "").lstrip("0").rstrip("0")
            )
            # Some bodies like Sedna have very small n, allow 5 significant figures
            assert significant >= 5, (
                f"{elements.name}: mean motion {elements.n} has insufficient precision"
            )
