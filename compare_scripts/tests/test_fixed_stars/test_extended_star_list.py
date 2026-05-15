"""
Tests for the extended FIXED_STARS list in comparison_utils.py.

Validates that all stars in the FIXED_STARS list are:
1. Recognized by libephemeris
2. Match pyswisseph positions within tolerance
3. Have correct astrological positions
"""

import pytest
import swisseph as swe
import libephemeris as ephem

from compare_scripts.comparison_utils import FIXED_STARS


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TOLERANCES
# ============================================================================

POSITION_TOL = 0.001  # degrees - 3.6 arcseconds


# ============================================================================
# TESTS
# ============================================================================


class TestExtendedStarListCoverage:
    """Test that the extended FIXED_STARS list has adequate coverage."""

    def test_minimum_star_count(self):
        """FIXED_STARS should contain at least 65 unique stars."""
        unique_stars = list(dict.fromkeys(FIXED_STARS))
        assert len(unique_stars) >= 65, (
            f"Expected at least 65 unique stars, got {len(unique_stars)}"
        )

    def test_royal_stars_included(self):
        """All four Royal Stars of Persia should be included."""
        royal_stars = {"Aldebaran", "Regulus", "Antares", "Fomalhaut"}
        included = royal_stars & set(FIXED_STARS)
        assert included == royal_stars, f"Missing Royal Stars: {royal_stars - included}"

    def test_behenian_stars_included(self):
        """Key Behenian stars should be included."""
        behenian_stars = {
            "Algol",
            "Alcyone",
            "Sirius",
            "Procyon",
            "Regulus",
            "Algorab",
            "Spica",
            "Arcturus",
            "Alphecca",
            "Antares",
            "Vega",
            "Deneb Algedi",
            "Altair",
            "Markab",
            "Scheat",
        }
        included = behenian_stars & set(FIXED_STARS)
        missing = behenian_stars - included
        assert len(missing) == 0, f"Missing Behenian stars: {missing}"

    def test_pleiades_cluster_included(self):
        """Pleiades cluster stars should be included."""
        pleiades = {
            "Alcyone",
            "Asterope",
            "Celaeno",
            "Electra",
            "Maia",
            "Merope",
            "Taygeta",
            "Atlas",
            "Pleione",
        }
        included = pleiades & set(FIXED_STARS)
        missing = pleiades - included
        assert len(missing) == 0, f"Missing Pleiades: {missing}"

    def test_hyades_cluster_included(self):
        """Hyades cluster stars should be included."""
        hyades = {"Prima Hyadum", "Secunda Hyadum", "Theta Tauri", "Ain"}
        included = hyades & set(FIXED_STARS)
        missing = hyades - included
        assert len(missing) == 0, f"Missing Hyades: {missing}"

    def test_zodiacal_stars_included(self):
        """Key zodiacal stars should be included."""
        zodiacal = {
            "Hamal",
            "Sheratan",  # Aries
            "Acubens",  # Cancer
            "Algieba",
            "Denebola",  # Leo
            "Vindemiatrix",  # Virgo
            "Zubenelgenubi",
            "Zubeneschamali",  # Libra
            "Kaus Australis",
            "Nunki",  # Sagittarius
            "Algedi",
            "Deneb Algedi",  # Capricorn
            "Sadalsuud",
            "Sadalmelik",  # Aquarius
            "Alrescha",  # Pisces
        }
        included = zodiacal & set(FIXED_STARS)
        missing = zodiacal - included
        assert len(missing) == 0, f"Missing zodiacal stars: {missing}"


class TestExtendedStarListValidity:
    """Test that all stars in FIXED_STARS are valid and calculable."""

    @pytest.mark.parametrize("star_name", list(dict.fromkeys(FIXED_STARS)))
    def test_star_recognized_by_libephemeris(self, star_name):
        """Each star should be recognized by libephemeris."""
        jd = 2451545.0  # J2000
        pos, name, retflag = ephem.fixstar(star_name, jd, 0)

        assert pos[0] != 0.0 or pos[1] != 0.0, (
            f"Star '{star_name}' returned zero position"
        )

    # Known catalog divergences: stars where libephemeris resolves a different
    # component than pyswisseph (e.g. Algedi alpha-1 vs alpha-2).
    KNOWN_DIVERGENCES = {"Algedi"}

    @pytest.mark.parametrize("star_name", list(dict.fromkeys(FIXED_STARS)))
    def test_star_matches_pyswisseph(self, star_name):
        """Each star should match pyswisseph position within tolerance."""
        if star_name in self.KNOWN_DIVERGENCES:
            pytest.xfail(
                f"Known catalog divergence: {star_name} wrong component"
            )

        jd = 2451545.0  # J2000

        try:
            result_swe = swe.fixstar_ut(star_name, jd, 0)
        except Exception as e:
            pytest.skip(f"Star {star_name} not in pyswisseph: {e}")
            return

        result_py = ephem.fixstar_ut(star_name, jd, 0)

        lon_swe = result_swe[0][0]
        lon_py = result_py[0][0]

        diff = angular_diff(lon_swe, lon_py)

        assert diff < POSITION_TOL, (
            f"Star '{star_name}' longitude diff {diff:.6f}° exceeds tolerance. "
            f"SWE: {lon_swe:.6f}°, libephemeris: {lon_py:.6f}°"
        )


class TestExtendedStarListAstrology:
    """Test astrologically significant star positions."""

    def test_alcyone_pleiades_position(self):
        """Alcyone (brightest Pleiad) should be in late Taurus."""
        jd = 2451545.0  # J2000
        pos, _, retflag = ephem.fixstar("Alcyone", jd, 0)
        # Alcyone should be around 29-30° Taurus (59-60° absolute)
        assert 59 < pos[0] < 61, (
            f"Alcyone position {pos[0]:.2f}° not in expected Taurus range"
        )

    def test_algol_demon_star_position(self):
        """Algol (Demon Star) should be in Taurus."""
        jd = 2451545.0  # J2000
        pos, _, retflag = ephem.fixstar("Algol", jd, 0)
        # Algol should be around 26° Taurus (56° absolute)
        assert 55 < pos[0] < 58, (
            f"Algol position {pos[0]:.2f}° not in expected Taurus range"
        )

    def test_prima_hyadum_hyades_position(self):
        """Prima Hyadum (first Hyad) should be in Gemini."""
        jd = 2451545.0  # J2000
        pos, _, retflag = ephem.fixstar("Prima Hyadum", jd, 0)
        # Prima Hyadum should be around 5° Gemini (65° absolute)
        assert 64 < pos[0] < 67, (
            f"Prima Hyadum position {pos[0]:.2f}° not in expected Gemini range"
        )

    def test_deneb_algedi_behenian_position(self):
        """Deneb Algedi (tail of the goat) should be in Aquarius."""
        jd = 2451545.0  # J2000
        pos, _, retflag = ephem.fixstar("Deneb Algedi", jd, 0)
        # Deneb Algedi should be around 23° Aquarius (323° absolute)
        assert 322 < pos[0] < 325, (
            f"Deneb Algedi position {pos[0]:.2f}° not in expected Aquarius range"
        )

    def test_alphecca_behenian_position(self):
        """Alphecca (gem in crown) should be in Scorpio."""
        jd = 2451545.0  # J2000
        pos, _, retflag = ephem.fixstar("Alphecca", jd, 0)
        # Alphecca should be around 12° Scorpio (222° absolute)
        assert 221 < pos[0] < 224, (
            f"Alphecca position {pos[0]:.2f}° not in expected Scorpio range"
        )

    def test_algorab_behenian_position(self):
        """Algorab (the crow) should be in Libra."""
        jd = 2451545.0  # J2000
        pos, _, retflag = ephem.fixstar("Algorab", jd, 0)
        # Algorab should be around 13° Libra (193° absolute)
        assert 192 < pos[0] < 195, (
            f"Algorab position {pos[0]:.2f}° not in expected Libra range"
        )
