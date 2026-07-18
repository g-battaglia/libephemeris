# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for libephemeris.minor_bodies.

This file targets the remaining uncovered lines/branches in the minor-body
calculation module: orbital-element validation, the exterior-planet secular
perturbation branches, near-resonance amplification, the (h,k)/(p,q) forced
element edge cases, Kepler-equation non-convergence paths, the short-period
perturbation guards, the SPK fast paths in the position helpers, and the
JPL SBDB / Horizons auto-download networking branches (all mocked, no
real network access).
"""

from __future__ import annotations

import json
import logging
import sys

import pytest
from unittest.mock import patch

pytestmark = pytest.mark.usefixtures("allow_mocked_network")

import libephemeris.minor_bodies as mb
from libephemeris import state
from libephemeris.constants import (
    CERES,
    ERIS,
    NAIF_CERES,
    NAIF_ERIS,
)
from libephemeris.minor_bodies import (
    OrbitalElements,
    JUPITER_A,
    JUPITER_N,
    SATURN_A,
    URANUS_A,
    NEPTUNE_A,
    NEPTUNE_N,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


class _FakeResponse:
    """Minimal context-manager stand-in for urllib's HTTP response."""

    def __init__(self, payload: object) -> None:
        self._bytes = json.dumps(payload).encode("utf-8")

    def read(self) -> bytes:
        return self._bytes

    def __enter__(self) -> "_FakeResponse":
        return self

    def __exit__(self, *exc: object) -> bool:
        return False


@pytest.fixture(autouse=True)
def _clear_spk_and_caches():
    """Isolate SPK registrations and SBDB caches around each test."""
    saved_map = dict(state._SPK_BODY_MAP)
    state._SPK_BODY_MAP.clear()
    state._SPK_BODY_MAP.update(saved_map)
    mb._ASTEROID_ELEMENTS_CACHE.clear()
    mb._ASTEROID_NAME_CACHE.clear()
    yield
    state._SPK_BODY_MAP.clear()
    state._SPK_BODY_MAP.update(saved_map)
    mb._ASTEROID_ELEMENTS_CACHE.clear()
    mb._ASTEROID_NAME_CACHE.clear()


# ---------------------------------------------------------------------------
# OrbitalElements.validate() and validate_elements()  (463-494, 506-510)
# ---------------------------------------------------------------------------


def test_validate_flags_all_invalid_branches():
    """Non-positive a, negative e, out-of-range i, non-positive n all warn."""
    bad = OrbitalElements(
        name="Bad",
        epoch=2451545.0,
        a=-1.0,  # a <= 0
        e=-0.5,  # e < 0  (and e < 1 with n <= 0)
        i=200.0,  # i outside [0, 180]
        omega=10.0,
        Omega=20.0,
        M0=30.0,
        n=-1.0,  # non-positive mean motion
    )
    issues = bad.validate()
    joined = " ".join(issues)
    assert "semi-major axis" in joined
    assert "eccentricity" in joined
    assert "inclination" in joined
    assert "mean motion" in joined


def test_validate_flags_hyperbolic_with_positive_a():
    """e >= 1 but a > 0 triggers the hyperbolic-orbit warning branch."""
    bad = OrbitalElements(
        name="Hyperbolic",
        epoch=2451545.0,
        a=1.0,  # a > 0
        e=1.5,  # e >= 1 with a > 0
        i=10.0,
        omega=10.0,
        Omega=20.0,
        M0=30.0,
        n=0.5,
    )
    issues = bad.validate()
    assert any("hyperbolic" in msg for msg in issues)


def test_validate_clean_elements_no_warnings():
    """Physically sane elements produce no warnings (empty list path)."""
    good = mb.MINOR_BODY_ELEMENTS[CERES]
    assert good.validate() == []


def test_validate_elements_logs_warnings(caplog):
    """validate_elements() forwards each issue to the module logger."""
    bad = OrbitalElements(
        name="LogMe",
        epoch=2451545.0,
        a=-2.0,
        e=0.1,
        i=10.0,
        omega=10.0,
        Omega=20.0,
        M0=30.0,
        n=0.2,
    )
    from libephemeris.logging_config import get_logger

    logger = get_logger()
    original_propagate = logger.propagate
    logger.propagate = True
    try:
        with caplog.at_level(logging.WARNING, logger="libephemeris"):
            mb.validate_elements(bad)
    finally:
        logger.propagate = original_propagate

    assert any("validation" in r.message for r in caplog.records)


def test_validate_elements_clean_no_log(caplog):
    """validate_elements() with valid elements logs nothing (issues empty)."""
    with caplog.at_level(logging.WARNING, logger="libephemeris"):
        mb.validate_elements(mb.MINOR_BODY_ELEMENTS[CERES])
    assert not any("validation" in r.message for r in caplog.records)


# ---------------------------------------------------------------------------
# _calc_laplace_coefficients near-singular denominator  (888->885)
# ---------------------------------------------------------------------------


def test_laplace_coefficient_near_unity_alpha_skips_singular_term():
    """alpha ~ 1 makes denom <= 1e-10 at psi=0, exercising the loop skip."""
    value = mb._calc_laplace_coefficients(0.999999, 1.5, 1)
    assert value > 0.0


# ---------------------------------------------------------------------------
# calc_secular_perturbation_rates exterior-planet branches
#   981->998, 1009->1025, 1036->1054, 1069->1105
# ---------------------------------------------------------------------------


def test_secular_rates_exterior_to_all_planets():
    """A body beyond Neptune hits every exterior-planet branch in one call."""
    d_omega, d_Omega, d_n = mb.calc_secular_perturbation_rates(
        mb.MINOR_BODY_ELEMENTS[ERIS]
    )
    assert isinstance(d_omega, float)
    assert isinstance(d_Omega, float)
    assert isinstance(d_n, float)


@pytest.mark.parametrize(
    "a_value",
    [JUPITER_A, SATURN_A, URANUS_A, NEPTUNE_A],
    ids=["jupiter", "saturn", "uranus", "neptune"],
)
def test_secular_rates_a_equals_planet_skips_both_branches(a_value):
    """a exactly equal to a planet's a takes neither the < nor > branch.

    This exercises the fall-through arcs 981->998, 1009->1025, 1036->1054
    and 1069->1105 where both the interior `if` and exterior `elif` are False.
    """
    el = OrbitalElements(
        name="OnPlanet",
        epoch=2451545.0,
        a=a_value,
        e=0.1,
        i=5.0,
        omega=10.0,
        Omega=20.0,
        M0=30.0,
        n=0.05,
    )
    d_omega, d_Omega, d_n = mb.calc_secular_perturbation_rates(el)
    assert isinstance(d_omega, float)
    assert isinstance(d_Omega, float)
    assert isinstance(d_n, float)


def test_secular_rates_dn_zero_for_near_parabolic(caplog):
    """e >= 0.99 forces the d_n = 0.0 else-branch (1111)."""
    el = OrbitalElements(
        name="HighE",
        epoch=2451545.0,
        a=39.0,
        e=0.995,
        i=10.0,
        omega=10.0,
        Omega=20.0,
        M0=30.0,
        n=0.004,
    )
    _, _, d_n = mb.calc_secular_perturbation_rates(el)
    assert d_n == 0.0


def test_secular_rates_tiny_mean_motion_skips_resonance_loop():
    """n <= 1e-20 makes both the d_n guard and the resonance guard False."""
    el = OrbitalElements(
        name="TinyN",
        epoch=2451545.0,
        a=39.0,
        e=0.2,
        i=10.0,
        omega=10.0,
        Omega=20.0,
        M0=30.0,
        n=1e-25,
    )
    _, _, d_n = mb.calc_secular_perturbation_rates(el)
    assert d_n == 0.0


def test_secular_rates_near_resonance_amplification():
    """A mean motion tuned near the 2:3 Neptune resonance amplifies rates."""
    # delta = |2n - 3*NEPTUNE_N| / n ~ 0.02 -> amplification factor 2.5
    n = 0.009123
    el = OrbitalElements(
        name="NearRes",
        epoch=2451545.0,
        a=39.0,
        e=0.1,
        i=5.0,
        omega=10.0,
        Omega=20.0,
        M0=30.0,
        n=n,
    )
    delta = abs(2 * n - 3 * NEPTUNE_N) / n
    assert 0.001 < delta < 0.05  # confirms the amplification window is selected
    d_omega, d_Omega, _ = mb.calc_secular_perturbation_rates(el)
    assert isinstance(d_omega, float)
    assert isinstance(d_Omega, float)


# ---------------------------------------------------------------------------
# _calc_forced_elements edge cases  (1234, 1284, 1293)
# ---------------------------------------------------------------------------


def test_forced_elements_static_planet_table():
    """jd_tt=None selects the static J2000.0 planet table (1234)."""
    el = OrbitalElements(
        name="Static",
        epoch=2451545.0,
        a=2.7,
        e=0.1,
        i=5.0,
        omega=10.0,
        Omega=20.0,
        M0=30.0,
        n=0.21,
    )
    result = mb._calc_forced_elements(el, jd_tt=None)
    assert len(result) == 6


def test_forced_elements_within_hill_sphere_skips_planet():
    """a equal to Jupiter's a lands inside Jupiter's Hill sphere (1284)."""
    el = OrbitalElements(
        name="CoOrbital",
        epoch=2451545.0,
        a=JUPITER_A,  # abs(a - a_p) = 0 < hill_radius
        e=0.1,
        i=5.0,
        omega=10.0,
        Omega=20.0,
        M0=30.0,
        n=0.08,
    )
    result = mb._calc_forced_elements(el)
    assert len(result) == 6


def test_forced_elements_negative_a_skips_via_alpha_guard():
    """Negative a yields alpha <= 0 for the inner case, hitting 1293."""
    el = OrbitalElements(
        name="NegA",
        epoch=2451545.0,
        a=-5.0,  # a < a_p for every planet, alpha = a/a_p < 0
        e=0.5,
        i=5.0,
        omega=10.0,
        Omega=20.0,
        M0=30.0,
        n=0.08,
    )
    g, h, k, s, p, q = mb._calc_forced_elements(el)
    # No planet contributes, so all accumulators stay zero.
    assert (g, h, k, s, p, q) == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


# ---------------------------------------------------------------------------
# apply_secular_perturbations free-vector zero branches  (1462, 1489)
# ---------------------------------------------------------------------------


def test_apply_secular_perturbations_zero_free_vectors():
    """Zero forced and zero osculating vectors collapse e_free/i_free to 0.

    A negative semi-major axis means no planet contributes (alpha <= 0), so
    g and s stay zero and the forced vectors are zero. With e=0 and i=0 the
    osculating vectors are also zero, so e_free and i_free fall below 1e-15,
    selecting beta=0.0 (1462) and gamma=0.0 (1489).
    """
    el = OrbitalElements(
        name="ZeroFree",
        epoch=2451545.0,
        a=-5.0,
        e=0.0,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        M0=10.0,
        n=0.05,
    )
    out = mb.apply_secular_perturbations(el, 2451545.0 + 4000.0)
    assert len(out) == 6


# ---------------------------------------------------------------------------
# Kepler-equation non-convergence paths  (2125-2133, 2168->2185, 2175, 2186-2188)
# ---------------------------------------------------------------------------


def test_elliptic_kepler_non_convergence_warns(caplog):
    """A negative tolerance prevents convergence, hitting the warning path."""
    with caplog.at_level(logging.WARNING, logger="libephemeris"):
        e_anom = mb.solve_kepler_equation_elliptic(1.0, 0.5, tol=-1.0)
    assert isinstance(e_anom, float)


def test_hyperbolic_kepler_non_convergence_loop_exhausts(caplog):
    """Negative tolerance exhausts the 50-iteration loop (2168->2185, 2186)."""
    with caplog.at_level(logging.WARNING, logger="libephemeris"):
        h_anom = mb.solve_kepler_equation_hyperbolic(1.0, 1.5, tol=-1.0)
    assert isinstance(h_anom, float)


def test_hyperbolic_kepler_zero_derivative_break():
    """e ~ 1 with M=0 drives fp below 1e-15, hitting the break (2175)."""
    h_anom = mb.solve_kepler_equation_hyperbolic(0.0, 1.0 + 1e-16, tol=1e-8)
    assert h_anom == 0.0


# ---------------------------------------------------------------------------
# _calc_short_period_correction guards  (2323, 2382)
# ---------------------------------------------------------------------------


def test_short_period_correction_returns_zero_for_hyperbolic():
    """a <= 0 (or e >= 1) returns the zero correction early (2323)."""
    el = OrbitalElements(
        name="Hyp",
        epoch=2451545.0,
        a=-1.0,
        e=1.5,
        i=10.0,
        omega=10.0,
        Omega=20.0,
        M0=30.0,
        n=0.5,
    )
    assert mb._calc_short_period_correction(el, 2451545.0 + 1000, 0.1, 0.2, 0.05) == (
        0.0,
        0.0,
        0.0,
    )


def test_short_period_correction_skips_near_synodic_resonance():
    """A mean motion equal to Jupiter's n zeroes dn, hitting the continue (2382)."""
    el = OrbitalElements(
        name="SynRes",
        epoch=2451545.0,
        a=6.0,
        e=0.1,
        i=5.0,
        omega=10.0,
        Omega=20.0,
        M0=30.0,
        n=JUPITER_N,  # dn = n_rad - n_p_rad = 0 for Jupiter
    )
    dx, dy, dz = mb._calc_short_period_correction(el, 2451545.0 + 1000, 0.1, 0.2, 0.05)
    assert all(isinstance(v, float) for v in (dx, dy, dz))


# ---------------------------------------------------------------------------
# Epoch-blend helpers  (2010-2011, 2035, 2050)
# ---------------------------------------------------------------------------


def test_closest_epoch_elements_returns_best():
    """_get_closest_epoch_elements returns the primary blend element (2010)."""
    el = mb._get_closest_epoch_elements(CERES, 2461000.5)
    assert el.name == "Ceres"


def test_epoch_blend_no_multi_entries_returns_base():
    """A body absent from the multi-epoch table returns base/None/0 (2035)."""
    # 111955 has single-epoch elements only (no multi-epoch table entry).
    base, secondary, weight = mb._get_epoch_elements_blend(111955, 2451545.0)
    assert secondary is None
    assert weight == 0.0


def test_epoch_blend_beyond_last_epoch_returns_nearest():
    """A date past the final multi-epoch entry returns nearest/None/0 (2050)."""
    multi = mb.MINOR_BODY_ELEMENTS_MULTI[CERES]
    jd = multi[-1].epoch + 100.0
    nearest, secondary, weight = mb._get_epoch_elements_blend(CERES, jd)
    assert secondary is None
    assert weight == 0.0
    assert nearest.epoch == multi[-1].epoch


# ---------------------------------------------------------------------------
# calc_minor_body_heliocentric SPK fast path  (2692-2700)
# ---------------------------------------------------------------------------


def test_heliocentric_uses_registered_spk(monkeypatch):
    """A registered SPK and a successful position dispatch the SPK branch."""
    monkeypatch.setitem(state._SPK_BODY_MAP, CERES, ("/fake.bsp", NAIF_CERES))
    monkeypatch.setattr(
        "libephemeris.spk.calc_spk_body_position",
        lambda *a, **k: (101.5, 3.25, 2.7, 0.0, 0.0, 0.0),
    )
    lon, lat, dist = mb.calc_minor_body_heliocentric(CERES, 2451545.0)
    assert (lon, lat, dist) == (101.5, 3.25, 2.7)


def test_heliocentric_spk_none_result_falls_back(monkeypatch):
    """A None SPK result falls through to the Keplerian path (2695->2708)."""
    monkeypatch.setitem(state._SPK_BODY_MAP, CERES, ("/fake.bsp", NAIF_CERES))
    monkeypatch.setattr(
        "libephemeris.spk.calc_spk_body_position",
        lambda *a, **k: None,
    )
    lon, lat, dist = mb.calc_minor_body_heliocentric(CERES, 2451545.0)
    assert 0.0 <= lon < 360.0
    assert dist > 0.0


def test_heliocentric_spk_exception_falls_back(monkeypatch):
    """An SPK exception is caught and the Keplerian path runs (2698-2700)."""
    monkeypatch.setitem(state._SPK_BODY_MAP, CERES, ("/fake.bsp", NAIF_CERES))

    def _raise(*a, **k):
        raise ValueError("spk boom")

    monkeypatch.setattr("libephemeris.spk.calc_spk_body_position", _raise)
    lon, lat, dist = mb.calc_minor_body_heliocentric(CERES, 2451545.0)
    assert 0.0 <= lon < 360.0
    assert dist > 0.0


# ---------------------------------------------------------------------------
# fetch_orbital_elements_from_sbdb error paths  (2824, 2836, 2856-2857)
# ---------------------------------------------------------------------------


def test_fetch_sbdb_zero_epoch_returns_none():
    """An epoch of 0 is treated as missing and returns None (2824)."""
    payload = {
        "orbit": {
            "epoch": 0,
            "elements": [
                {"name": "a", "value": "2.7"},
                {"name": "e", "value": "0.1"},
                {"name": "i", "value": "5"},
                {"name": "n", "value": "0.2"},
            ],
        },
        "object": {"fullname": "X"},
    }
    with patch("urllib.request.urlopen", return_value=_FakeResponse(payload)):
        assert mb.fetch_orbital_elements_from_sbdb(900001) is None


def test_fetch_sbdb_missing_core_element_returns_none():
    """A record missing a core element (n) returns None (2836)."""
    payload = {
        "orbit": {
            "epoch": 2461000.5,
            "elements": [
                {"name": "a", "value": "2.7"},
                {"name": "e", "value": "0.1"},
                {"name": "i", "value": "5"},
            ],
        },
        "object": {"fullname": "X"},
    }
    with patch("urllib.request.urlopen", return_value=_FakeResponse(payload)):
        assert mb.fetch_orbital_elements_from_sbdb(900002) is None


def test_fetch_sbdb_malformed_element_returns_none():
    """A KeyError while reading element values is caught (2856-2857)."""
    payload = {
        "orbit": {
            "epoch": 2461000.5,
            "elements": [{"name": "a"}],  # missing "value" -> KeyError
        },
        "object": {"fullname": "X"},
    }
    with patch("urllib.request.urlopen", return_value=_FakeResponse(payload)):
        assert mb.fetch_orbital_elements_from_sbdb(900003) is None


# ---------------------------------------------------------------------------
# get_asteroid_number SBDB parsing  (3018-3051, 3055-3056)
# ---------------------------------------------------------------------------


def test_get_asteroid_number_spkid_above_offset():
    """spkid >= 2_000_000 maps to spkid - 2_000_000 (3027-3028)."""
    with patch(
        "urllib.request.urlopen",
        return_value=_FakeResponse({"object": {"spkid": "2000433"}}),
    ):
        assert mb.get_asteroid_number("cov_name_a") == 433


def test_get_asteroid_number_spkid_direct():
    """spkid below the offset is used directly (3031)."""
    with patch(
        "urllib.request.urlopen",
        return_value=_FakeResponse({"object": {"spkid": "12345"}}),
    ):
        assert mb.get_asteroid_number("cov_name_b") == 12345


def test_get_asteroid_number_spkid_unparseable():
    """A non-numeric spkid raises ValueError and returns None (3032-3033)."""
    with patch(
        "urllib.request.urlopen",
        return_value=_FakeResponse({"object": {"spkid": "not-a-number"}}),
    ):
        assert mb.get_asteroid_number("cov_name_c") is None


def test_get_asteroid_number_from_designation():
    """A digit designation is parsed when spkid is absent (3036-3038)."""
    with patch(
        "urllib.request.urlopen",
        return_value=_FakeResponse({"object": {"des": "777"}}),
    ):
        assert mb.get_asteroid_number("cov_name_d") == 777


def test_get_asteroid_number_from_fullname():
    """A leading numeric token in fullname is parsed (3040-3044, 3049-3051)."""
    with patch(
        "urllib.request.urlopen",
        return_value=_FakeResponse({"object": {"fullname": "888 Foo"}}),
    ):
        assert mb.get_asteroid_number("cov_name_e") == 888


def test_get_asteroid_number_unparseable_returns_none():
    """No usable identifier yields None (3046)."""
    with patch(
        "urllib.request.urlopen",
        return_value=_FakeResponse({"object": {"fullname": "Foo Bar"}}),
    ):
        assert mb.get_asteroid_number("cov_name_f") is None


def test_get_asteroid_number_handles_type_error():
    """A non-dict payload raises TypeError and is caught (3055-3056)."""
    with patch("urllib.request.urlopen", return_value=_FakeResponse(5)):
        assert mb.get_asteroid_number("cov_name_g") is None


# ---------------------------------------------------------------------------
# is_spk_downloadable / is_in_jpl_major_body_index  (3182, 3200)
# ---------------------------------------------------------------------------


def test_is_spk_downloadable_returns_membership():
    """A non-blocked body in the maps is reported as downloadable (3182)."""
    assert mb.is_spk_downloadable(CERES) is True
    assert mb.is_spk_downloadable(424242) is False


def test_is_in_jpl_major_body_index_always_false():
    """Backward-compat shim always returns False (3200)."""
    assert mb.is_in_jpl_major_body_index(CERES) is False


# ---------------------------------------------------------------------------
# auto_download_asteroid_spk paths  (3282-3283, 3294-3354)
# ---------------------------------------------------------------------------


def test_auto_download_uses_spk_body_name_map(monkeypatch):
    """A body present only in SPK_BODY_NAME_MAP resolves its name path (3282-3283)."""
    # ERIS is in SPK_BODY_NAME_MAP but not MAJOR_ASTEROID_SPK_INFO.
    monkeypatch.setitem(state._SPK_BODY_MAP, ERIS, ("/eris.bsp", NAIF_ERIS))
    assert mb.auto_download_asteroid_spk(ERIS) == "/eris.bsp"


def test_auto_download_success(monkeypatch):
    """Full download path returns the downloaded SPK path (3294-3344)."""
    state._SPK_BODY_MAP.pop(CERES, None)
    monkeypatch.setattr(
        "libephemeris.spk.download_and_register_spk",
        lambda **k: "/downloaded.bsp",
    )
    assert mb.auto_download_asteroid_spk(CERES) == "/downloaded.bsp"


def test_auto_download_range_outside_horizons_limits(monkeypatch):
    """A requested range entirely outside Horizons limits returns None (3311-3320)."""
    state._SPK_BODY_MAP.pop(CERES, None)
    monkeypatch.setattr(
        "libephemeris.spk.download_and_register_spk",
        lambda **k: "/unused.bsp",
    )
    # Both bounds above HORIZONS_SPK_JD_MAX -> clamp collapses the window.
    assert (
        mb.auto_download_asteroid_spk(CERES, jd_start=1.0e7, jd_end=1.0e7 + 10.0)
        is None
    )


def test_auto_download_range_clamped_but_valid(monkeypatch):
    """A partially-out-of-range request is clamped yet still downloads (3314->3323)."""
    state._SPK_BODY_MAP.pop(CERES, None)
    monkeypatch.setattr(
        "libephemeris.spk.download_and_register_spk",
        lambda **k: "/clamped.bsp",
    )
    # jd_start below the Horizons minimum, jd_end well within range -> clamp
    # leaves a valid, non-empty window.
    result = mb.auto_download_asteroid_spk(CERES, jd_start=2200000.0, jd_end=2451545.0)
    assert result == "/clamped.bsp"


def test_auto_download_offline_padding_extends_calendar_clamp(monkeypatch):
    """Offline fitting can obtain and later discard one safe boundary day."""
    state._SPK_BODY_MAP.pop(CERES, None)
    received = {}

    def _download(**kwargs):
        received.update(kwargs)
        return "/padded.bsp"

    monkeypatch.setattr("libephemeris.spk.download_and_register_spk", _download)

    result = mb.auto_download_asteroid_spk(
        CERES,
        jd_start=2200000.0,
        jd_end=2800000.0,
        boundary_padding_days=1.0,
    )

    assert result == "/padded.bsp"
    assert received["start"] == "1599-12-31"
    assert received["end"] == "2500-01-02"


def test_auto_download_value_error_returns_none(monkeypatch):
    """A ValueError from the downloader is caught and returns None (3346-3348)."""
    state._SPK_BODY_MAP.pop(CERES, None)

    def _boom(**kwargs):
        raise ValueError("bad request")

    monkeypatch.setattr("libephemeris.spk.download_and_register_spk", _boom)
    assert mb.auto_download_asteroid_spk(CERES) is None


def test_auto_download_os_error_returns_none(monkeypatch):
    """An OSError/ConnectionError from the downloader returns None (3349-3354)."""
    state._SPK_BODY_MAP.pop(CERES, None)

    def _neterr(**kwargs):
        raise OSError("network down")

    monkeypatch.setattr("libephemeris.spk.download_and_register_spk", _neterr)
    assert mb.auto_download_asteroid_spk(CERES) is None


# ---------------------------------------------------------------------------
# is_spk_available_for_body ImportError guard  (3382-3383)
# ---------------------------------------------------------------------------


def test_is_spk_available_handles_import_error(monkeypatch):
    """A failed state import is treated as 'not available' (3382-3383).

    `from . import state` resolves via the package attribute first, so the
    attribute must be removed in addition to poisoning sys.modules to force
    the ImportError. monkeypatch restores both afterwards.
    """
    import libephemeris

    monkeypatch.setitem(sys.modules, "libephemeris.state", None)
    monkeypatch.delattr(libephemeris, "state")
    assert mb.is_spk_available_for_body(CERES) is False


# ---------------------------------------------------------------------------
# ensure_major_asteroid_spk name resolution  (3439-3441, 3454-3456)
# ---------------------------------------------------------------------------


def test_ensure_spk_name_falls_back_to_id_for_unknown_map_body(monkeypatch):
    """A SPK_BODY_NAME_MAP body absent from elements uses str(id) (3439)."""
    elems_no_eris = {k: v for k, v in mb.MINOR_BODY_ELEMENTS.items() if k != ERIS}
    monkeypatch.setattr(mb, "MINOR_BODY_ELEMENTS", elems_no_eris)
    monkeypatch.setitem(state._SPK_BODY_MAP, ERIS, ("/eris.bsp", NAIF_ERIS))
    # Already cached -> returns True without download.
    assert mb.ensure_major_asteroid_spk(ERIS) is True


def test_ensure_spk_name_falls_back_to_id_for_unmapped_body():
    """An entirely unmapped body uses str(id) and cannot download (3441)."""
    # 424242 is in no map; auto_download returns None -> overall False.
    assert mb.ensure_major_asteroid_spk(424242) is False


def test_ensure_spk_with_jd_passes_window(monkeypatch):
    """Providing jd computes a ±10-year window before downloading (3453-3456)."""
    state._SPK_BODY_MAP.pop(CERES, None)
    captured = {}

    def _fake_download(body_id, jd_start=None, jd_end=None):
        captured["jd_start"] = jd_start
        captured["jd_end"] = jd_end
        return "/cached.bsp"

    monkeypatch.setattr(mb, "auto_download_asteroid_spk", _fake_download)
    assert mb.ensure_major_asteroid_spk(CERES, jd=2451545.0) is True
    assert captured["jd_start"] == 2451545.0 - 3652.5
    assert captured["jd_end"] == 2451545.0 + 3652.5


# ---------------------------------------------------------------------------
# calc_asteroid_by_number SPK and special-mapping fall-through
#   3579->3583, 3592-3602
# ---------------------------------------------------------------------------


def test_calc_asteroid_by_number_uses_spk(monkeypatch):
    """A registered SPK for the offset body id dispatches the SPK branch (3592)."""
    # 70000 -> body_id 80000, which is in neither MINOR_BODY_ELEMENTS nor the
    # special low-number mapping, so the function reaches its own SPK branch.
    body_id = 70000 + 10000  # AST_OFFSET = 10000
    monkeypatch.setitem(state._SPK_BODY_MAP, body_id, ("/fake.bsp", 2070000))
    monkeypatch.setattr(
        "libephemeris.spk.calc_spk_body_position",
        lambda *a, **k: (44.4, 1.1, 1.5, 0.0, 0.0, 0.0),
    )
    lon, lat, dist = mb.calc_asteroid_by_number(70000, 2451545.0)
    assert (lon, lat, dist) == (44.4, 1.1, 1.5)


def test_calc_asteroid_by_number_spk_exception_falls_back(monkeypatch):
    """An SPK exception falls through to the SBDB fetch path (3600-3602)."""
    body_id = 70001 + 10000
    monkeypatch.setitem(state._SPK_BODY_MAP, body_id, ("/fake.bsp", 2070001))

    def _raise(*a, **k):
        raise ValueError("spk boom")

    monkeypatch.setattr("libephemeris.spk.calc_spk_body_position", _raise)
    monkeypatch.setattr(mb, "fetch_orbital_elements_from_sbdb", lambda *a, **k: None)
    with pytest.raises(ValueError, match="not found"):
        mb.calc_asteroid_by_number(70001, 2451545.0)


def test_calc_asteroid_by_number_special_mapping_fallthrough(monkeypatch):
    """A special-mapping body absent from elements falls through to SBDB (3579->3583)."""
    elems_no_ceres = {k: v for k, v in mb.MINOR_BODY_ELEMENTS.items() if k != CERES}
    monkeypatch.setattr(mb, "MINOR_BODY_ELEMENTS", elems_no_ceres)
    monkeypatch.setattr(mb, "fetch_orbital_elements_from_sbdb", lambda *a, **k: None)
    with pytest.raises(ValueError, match="not found"):
        mb.calc_asteroid_by_number(1, 2451545.0)
