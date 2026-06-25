"""
Tests for the Delta T model selector.

``set_delta_t_model()`` / ``get_delta_t_model()`` and the
``LIBEPHEMERIS_DELTAT_MODEL`` environment variable choose the Delta T model used
*after* the user-defined and IERS-observed priorities:

    - ``smh2016`` (default): Stephenson-Morrison-Hohenkerk 2016 via Skyfield.
    - ``espenak_meeus``: the classic NASA Espenak & Meeus (2006) polynomial set,
      a self-contained clean-room implementation.

The selector is clean-room: there is intentionally **no** ``swisseph`` mode
(libephemeris never imports pyswisseph). Exact Swiss parity, when needed for
validation, is injected externally via ``set_delta_t_userdef``.

References:
    Espenak, F. & Meeus, J. (2006). "Five Millennium Canon of Solar Eclipses:
    -1999 to +3000." NASA/TP-2006-214141.
"""

from __future__ import annotations

import os

import pytest

import libephemeris as ephem
from libephemeris import _config_toml as cfg


def _em_reference(year: float) -> float:
    """Espenak & Meeus (2006) Delta T in seconds (the published polynomials)."""
    y = year
    if y < -500:
        u = (y - 1820) / 100.0
        return -20 + 32 * u * u
    if y < 500:
        u = y / 100.0
        return (
            10583.6
            - 1014.41 * u
            + 33.78311 * u**2
            - 5.952053 * u**3
            - 0.1798452 * u**4
            + 0.022174192 * u**5
            + 0.0090316521 * u**6
        )
    if y < 1600:
        u = (y - 1000) / 100.0
        return (
            1574.2
            - 556.01 * u
            + 71.23472 * u**2
            + 0.319781 * u**3
            - 0.8503463 * u**4
            - 0.005050998 * u**5
            + 0.0083572073 * u**6
        )
    if y < 1700:
        t = y - 1600
        return 120 - 0.9808 * t - 0.01532 * t**2 + t**3 / 7129.0
    if y < 1800:
        t = y - 1700
        return (
            8.83 + 0.1603 * t - 0.0059285 * t**2 + 0.00013336 * t**3 - t**4 / 1174000.0
        )
    if y < 1860:
        t = y - 1800
        return (
            13.72
            - 0.332447 * t
            + 0.0068612 * t**2
            + 0.0041116 * t**3
            - 0.00037436 * t**4
            + 0.0000121272 * t**5
            - 0.0000001699 * t**6
            + 0.000000000875 * t**7
        )
    if y < 1900:
        t = y - 1860
        return (
            7.62
            + 0.5737 * t
            - 0.251754 * t**2
            + 0.01680668 * t**3
            - 0.0004473624 * t**4
            + t**5 / 233174.0
        )
    if y < 1920:
        t = y - 1900
        return (
            -2.79 + 1.494119 * t - 0.0598939 * t**2 + 0.0061966 * t**3 - 0.000197 * t**4
        )
    if y < 1941:
        t = y - 1920
        return 21.20 + 0.84493 * t - 0.076100 * t**2 + 0.0020936 * t**3
    if y < 1961:
        t = y - 1950
        return 29.07 + 0.407 * t - t**2 / 233.0 + t**3 / 2547.0
    if y < 1986:
        t = y - 1975
        return 45.45 + 1.067 * t - t**2 / 260.0 - t**3 / 718.0
    if y < 2005:
        t = y - 2000
        return (
            63.86
            + 0.3345 * t
            - 0.060374 * t**2
            + 0.0017275 * t**3
            + 0.000651814 * t**4
            + 0.00002373599 * t**5
        )
    if y < 2050:
        t = y - 2000
        return 62.92 + 0.32217 * t + 0.005589 * t**2
    if y < 2150:
        return -20 + 32 * ((y - 1820) / 100.0) ** 2 - 0.5628 * (2150 - y)
    u = (y - 1820) / 100.0
    return -20 + 32 * u * u


def _jan1_jd(year: int) -> float:
    return ephem.julday(year, 1, 1, 0.0)


@pytest.fixture(autouse=True)
def _restore_deltat_state():
    """Ensure global Delta T state is restored after every test."""
    yield
    ephem.set_delta_t_model(None)
    ephem.set_delta_t_userdef(None)
    os.environ.pop("LIBEPHEMERIS_DELTAT_MODEL", None)
    cfg.reset()


class TestModelSelectorBasics:
    @pytest.mark.unit
    def test_default_is_smh2016(self):
        ephem.set_delta_t_model(None)
        os.environ.pop("LIBEPHEMERIS_DELTAT_MODEL", None)
        assert ephem.get_delta_t_model() == "smh2016"

    @pytest.mark.unit
    def test_set_and_get(self):
        ephem.set_delta_t_model("espenak_meeus")
        assert ephem.get_delta_t_model() == "espenak_meeus"
        ephem.set_delta_t_model("smh2016")
        assert ephem.get_delta_t_model() == "smh2016"

    @pytest.mark.unit
    def test_case_insensitive(self):
        ephem.set_delta_t_model("Espenak_Meeus")
        assert ephem.get_delta_t_model() == "espenak_meeus"

    @pytest.mark.unit
    def test_none_defers_to_env(self):
        ephem.set_delta_t_model(None)
        os.environ["LIBEPHEMERIS_DELTAT_MODEL"] = "espenak_meeus"
        assert ephem.get_delta_t_model() == "espenak_meeus"

    @pytest.mark.unit
    def test_none_defers_to_toml(self, tmp_path):
        path = tmp_path / "libephemeris-config.toml"
        path.write_text(
            '[libephemeris]\ndeltat_model = "espenak_meeus"\n',
            encoding="utf-8",
        )
        cfg.reset()
        assert cfg.load_config(path) is True
        ephem.set_delta_t_model(None)
        os.environ.pop("LIBEPHEMERIS_DELTAT_MODEL", None)

        assert ephem.get_delta_t_model() == "espenak_meeus"

    @pytest.mark.unit
    def test_invalid_model_raises(self):
        with pytest.raises(ValueError):
            ephem.set_delta_t_model("bogus")

    @pytest.mark.unit
    def test_no_swisseph_mode(self):
        """Clean-room: there must be NO 'swisseph' Delta T mode in the core."""
        with pytest.raises(ValueError):
            ephem.set_delta_t_model("swisseph")


class TestEspenakMeeusValues:
    """The espenak_meeus model must reproduce the published NASA polynomials."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "year", [1700, 1800, 1900, 1950, 2000, 1000, 500, -500, 2100, 2200]
    )
    def test_matches_published_polynomials(self, year):
        ephem.set_delta_t_model("espenak_meeus")
        got = ephem.deltat(_jan1_jd(year)) * 86400.0
        expected = _em_reference(year + (1 - 0.5) / 12.0)  # decimal year, Jan
        assert got == pytest.approx(expected, abs=0.05)

    @pytest.mark.unit
    def test_year_2000_anchor(self):
        """ΔT(2000) is ~63.8-63.9 s for both models (well-observed era)."""
        ephem.set_delta_t_model("espenak_meeus")
        assert ephem.deltat(_jan1_jd(2000)) * 86400.0 == pytest.approx(63.87, abs=0.1)


class TestModelInteraction:
    @pytest.mark.unit
    def test_userdef_overrides_model(self):
        """User-defined Delta T has priority over the selected model."""
        ephem.set_delta_t_model("espenak_meeus")
        ephem.set_delta_t_userdef(0.001)  # days
        assert ephem.deltat(_jan1_jd(2000)) == pytest.approx(0.001, abs=1e-12)

    @pytest.mark.unit
    def test_models_agree_in_modern_era(self):
        """Both models are within a few seconds in the well-observed range."""
        jd = _jan1_jd(2000)
        ephem.set_delta_t_model("smh2016")
        smh = ephem.deltat(jd) * 86400.0
        ephem.set_delta_t_model("espenak_meeus")
        em = ephem.deltat(jd) * 86400.0
        assert abs(smh - em) < 1.0

    @pytest.mark.unit
    def test_models_diverge_in_far_future(self):
        """The models legitimately diverge where Delta T is extrapolated."""
        jd = _jan1_jd(2400)
        ephem.set_delta_t_model("smh2016")
        smh = ephem.deltat(jd) * 86400.0
        ephem.set_delta_t_model("espenak_meeus")
        em = ephem.deltat(jd) * 86400.0
        assert abs(smh - em) > 10.0

    @pytest.mark.unit
    def test_deltat_ex_uses_selected_model(self):
        """deltat_ex() must use the same selected model as deltat()."""
        jd = _jan1_jd(2400)
        ephem.set_delta_t_model("espenak_meeus")

        assert ephem.deltat_ex(jd, ephem.FLG_SWIEPH) == pytest.approx(
            ephem.deltat(jd),
            abs=1e-15,
        )
