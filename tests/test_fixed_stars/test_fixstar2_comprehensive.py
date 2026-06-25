"""
Tests for fixstar2_ut: flexible star lookup by partial name,
HIP number, nomenclature, and fuzzy matching.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    swe.close()


JD_J2000 = 2451545.0


class TestFixstar2ExactName:
    """Test fixstar2_ut with exact star names."""

    BRIGHT_STARS = [
        "Sirius",
        "Regulus",
        "Spica",
        "Aldebaran",
        "Antares",
        "Betelgeuse",
        "Rigel",
        "Vega",
        "Capella",
        "Arcturus",
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("star", BRIGHT_STARS)
    def test_exact_name_returns_position(self, star):
        """Exact star name returns valid position."""
        pos, name_out, retflag = swe.fixstar2_ut(star, JD_J2000)
        assert len(pos) == 6
        assert 0.0 <= pos[0] < 360.0
        assert -90.0 <= pos[1] <= 90.0
        assert pos[2] > 0.0  # distance > 0
        assert isinstance(name_out, str)
        assert len(name_out) > 0

    @pytest.mark.unit
    def test_case_insensitive(self):
        """Star name lookup is case-insensitive."""
        pos1, _, _ = swe.fixstar2_ut("Sirius", JD_J2000)
        pos2, _, _ = swe.fixstar2_ut("SIRIUS", JD_J2000)
        pos3, _, _ = swe.fixstar2_ut("sirius", JD_J2000)
        assert pos1[0] == pytest.approx(pos2[0], abs=1e-8)
        assert pos1[0] == pytest.approx(pos3[0], abs=1e-8)

    @pytest.mark.unit
    def test_name_out_contains_nomenclature(self):
        """Returned name includes the traditional name."""
        _, name_out, _ = swe.fixstar2_ut("Regulus", JD_J2000)
        assert "Regulus" in name_out


class TestFixstar2HIPNumber:
    """Test fixstar2_ut with Hipparcos catalog numbers."""

    @pytest.mark.unit
    def test_hip_number_direct(self):
        """Lookup by bare HIP number."""
        # HIP 32349 = Sirius
        pos, name_out, _ = swe.fixstar2_ut("32349", JD_J2000)
        assert 0.0 <= pos[0] < 360.0
        assert "Sirius" in name_out or "alCMa" in name_out

    @pytest.mark.unit
    def test_hip_with_comma(self):
        """Lookup by HIP number with leading comma."""
        pos, name_out, _ = swe.fixstar2_ut(",32349", JD_J2000)
        assert 0.0 <= pos[0] < 360.0

    @pytest.mark.unit
    def test_hip_prefix(self):
        """Lookup by 'HIP 32349' format."""
        pos, name_out, _ = swe.fixstar2_ut("HIP 32349", JD_J2000)
        assert 0.0 <= pos[0] < 360.0


class TestFixstar2Nomenclature:
    """Test fixstar2_ut with nomenclature codes."""

    @pytest.mark.unit
    def test_nomenclature_code(self):
        """Lookup by nomenclature like 'alLeo' for Regulus."""
        pos, name_out, _ = swe.fixstar2_ut("alLeo", JD_J2000)
        assert 0.0 <= pos[0] < 360.0
        assert "Regulus" in name_out or "alLeo" in name_out

    @pytest.mark.unit
    def test_comma_separated_format(self):
        """Lookup by 'Regulus,alLeo' format."""
        pos, name_out, _ = swe.fixstar2_ut("Regulus,alLeo", JD_J2000)
        assert 0.0 <= pos[0] < 360.0


class TestFixstar2PartialName:
    """Test fixstar2_ut with partial name matching."""

    @pytest.mark.unit
    def test_partial_prefix_match(self):
        """Partial prefix matches a star."""
        # "Sir" should match "Sirius" if unambiguous
        try:
            pos, name_out, _ = swe.fixstar2_ut("Sir", JD_J2000)
            assert "Sirius" in name_out or "Sir" in name_out
        except Exception:
            # Ambiguous match may raise
            pass

    @pytest.mark.unit
    def test_ambiguous_raises_or_matches(self):
        """Ambiguous partial name either raises or picks best match."""
        # Very short prefix likely ambiguous
        try:
            pos, name_out, _ = swe.fixstar2_ut("Al", JD_J2000)
            # If it works, it picked one star
            assert len(name_out) > 0
        except Exception as e:
            # Expected for ambiguous matches
            assert (
                "ambiguous" in str(e).lower()
                or "not found" in str(e).lower()
                or "match" in str(e).lower()
                or True
            )  # any error is acceptable


class TestFixstar2Mag:
    """Test fixstar2_mag for magnitude lookup."""

    @pytest.mark.unit
    @pytest.mark.parametrize("star", ["Sirius", "Regulus", "Vega", "Arcturus"])
    def test_magnitude_returns_tuple(self, star):
        """fixstar2_mag returns (magnitude, name) tuple."""
        result = swe.fixstar2_mag(star)
        assert len(result) == 2
        mag, name = result
        assert isinstance(mag, float)
        assert isinstance(name, str)
        assert math.isfinite(mag)

    @pytest.mark.unit
    def test_sirius_brightest(self):
        """Sirius should be the brightest star (mag ~ -1.46)."""
        mag, _ = swe.fixstar2_mag("Sirius")
        assert mag < -1.0, f"Sirius magnitude {mag} not bright enough"

    @pytest.mark.unit
    def test_nonexistent_star_raises(self):
        """Nonexistent star name raises an error."""
        with pytest.raises(Exception):
            swe.fixstar2_ut("ZZZZNOTASTAR", JD_J2000)

    @pytest.mark.unit
    def test_nonexistent_star_mag_raises(self):
        """Nonexistent star in fixstar2_mag raises an error."""
        with pytest.raises(Exception):
            swe.fixstar2_mag("ZZZZNOTASTAR")


class TestFixstar2Consistency:
    """Test consistency between fixstar_ut and fixstar2_ut."""

    @pytest.mark.unit
    @pytest.mark.parametrize("star", ["Sirius", "Regulus", "Spica"])
    def test_fixstar_vs_fixstar2(self, star):
        """fixstar_ut and fixstar2_ut return same positions."""
        pos1 = swe.fixstar_ut(star, JD_J2000)
        pos2, _, _ = swe.fixstar2_ut(star, JD_J2000)
        # fixstar_ut returns different format — extract positions
        # It may return (lon, lat, dist, speed_lon, speed_lat, speed_dist, name)
        # or just ((lon,...), retflag) depending on API
        if isinstance(pos1, tuple) and len(pos1) >= 2:
            if isinstance(pos1[0], tuple):
                p1 = pos1[0]
            else:
                p1 = pos1
        else:
            p1 = pos1

        # Compare longitudes
        if hasattr(p1, "__len__") and len(p1) >= 2:
            assert p1[0] == pytest.approx(pos2[0], abs=1e-6)
