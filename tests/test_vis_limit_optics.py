"""Observer-dependent corrections of ``vis_limit_mag`` (age + optical aid).

The reference API applies two observer corrections on top of the naked-eye
VISLIMIT limiting magnitude: the observer's *age* (via the dark-adapted eye
pupil) and, when ``HELFLAG_OPTICAL_PARAMS`` is set, an *optical instrument*
(binocular / telescope). Both are implemented as deltas that are exactly zero
for the naked-eye default ``(36, 1, 0, 0, 0, 0)`` so the calibrated naked-eye
behaviour is untouched.

The expected deltas below are taken from a black-box dump of the reference
``vis_limit_mag`` (192 cases at geopos ``(12.5, 41.9, 0)``, atmosphere
``(1013.25, 15, 40, 0)``), expressed relative to the *same-call naked
baseline* so this test needs no oracle at run time.

Precision notes (measured against the 192-row oracle):

* Snellen ratio: reproduced exactly (``5*log10(SN)``).
* Flag gating / naked baseline: exact.
* Age curve: max ~0.07 mag over the calibration ages, growing to ~0.15 mag at
  the extreme (age 90); the age effect is genuinely sky-brightness dependent,
  so a single fitted form cannot reach the oracle exactly.
* Optical aid: the reference's optical-aid model is a bespoke, sky-coupled
  formula (its magnification dependence is opposite to Schaefer 1990's
  published telescopic formula and its aperture dependence is stronger than
  the physical ``D**2`` light grasp). It is not reproducible clean-room below
  ~0.5 mag typical / ~1.3 mag worst case; the tolerances below reflect that
  measured floor, not a target.
"""

from __future__ import annotations

import math

import pytest

from libephemeris.constants import HELFLAG_OPTICAL_PARAMS
from libephemeris.heliacal import vis_limit_mag

GEOPOS = (12.5, 41.9, 0.0)
ATMO = (1013.25, 15.0, 40.0, 0.0)

# Live scenes from the oracle dump (object above the horizon).
JD_JUP = 2460390.2291666665  # Jupiter, bright twilight (photopic)
JD_SIR = 2460686.25  # Sirius, dark night

# Naked-eye baselines produced by this library (NOT the reference). These must
# stay put when the observer corrections are all no-ops. The values differ
# between the LEB and Skyfield backends by a few 1e-9 (ephemeris noise), so
# the guard uses a 1e-6 tolerance rather than bit identity.
NAKED_JUP = -0.45731009787402943
NAKED_SIR = 0.5869070076075849
NAKED_TOL = 1e-6

AGE_TOL = 0.16  # age-form floor (see module docstring)
OPTICS_TOL = 1.30  # optical-aid floor (see module docstring)
OPTICS_TOL_TIGHT = 0.30  # well-reproduced mid-size instruments


def _delta(jd, body, obs, flags):
    """Limiting-mag delta of ``obs`` vs the same-scene naked baseline."""
    base = vis_limit_mag(jd, GEOPOS, ATMO, (36, 1, 0, 0, 0, 0), body, 0)[1][0]
    got = vis_limit_mag(jd, GEOPOS, ATMO, obs, body, flags)[1][0]
    return got - base


class TestNakedBaselineUnchanged:
    """The naked-eye default must be unchanged by the observer corrections."""

    def test_jupiter_naked_bit_identical(self):
        got = vis_limit_mag(JD_JUP, GEOPOS, ATMO, (36, 1, 0, 0, 0, 0), "jupiter", 0)
        assert abs(got[1][0] - NAKED_JUP) < NAKED_TOL

    def test_sirius_naked_bit_identical(self):
        got = vis_limit_mag(JD_SIR, GEOPOS, ATMO, (36, 1, 0, 0, 0, 0), "sirius", 0)
        assert abs(got[1][0] - NAKED_SIR) < NAKED_TOL

    def test_naked_default_delta_is_zero(self):
        assert _delta(JD_JUP, "jupiter", (36, 1, 0, 0, 0, 0), 0) == 0.0


class TestSnellen:
    """Snellen ratio contributes exactly ``5*log10(SN)`` (sky-independent)."""

    @pytest.mark.parametrize("sn", [0.5, 1.0, 1.5, 2.0])
    def test_snellen_exact_jupiter(self, sn):
        got = _delta(JD_JUP, "jupiter", (36, sn, 0, 0, 0, 0), 0)
        assert got == pytest.approx(5.0 * math.log10(sn), abs=1e-9)

    @pytest.mark.parametrize("sn", [0.5, 1.5, 2.0])
    def test_snellen_exact_sirius(self, sn):
        # Sky-independent: identical term in a very different (dark) sky.
        got = _delta(JD_SIR, "sirius", (36, sn, 0, 0, 0, 0), 0)
        assert got == pytest.approx(5.0 * math.log10(sn), abs=1e-9)


class TestAgeCurve:
    """Age correction: sign, magnitude, clamp below ~23, monotonicity."""

    # Oracle deltas at the Jupiter (photopic) scene, relative to age 36.
    AGE_FACTS = {
        23: 0.196395,
        25: 0.166832,
        30: 0.091904,
        40: -0.062522,
        50: -0.223448,
        60: -0.391414,
        70: -0.567029,
        80: -0.750988,
        90: -0.944079,
    }

    @pytest.mark.parametrize("age,expected", sorted(AGE_FACTS.items()))
    def test_age_delta_matches_oracle(self, age, expected):
        got = _delta(JD_JUP, "jupiter", (age, 1, 0, 0, 0, 0), 0)
        assert got == pytest.approx(expected, abs=AGE_TOL)

    def test_age_applies_without_optical_flag(self):
        # Age is read from observer[0] even when HELFLAG_OPTICAL_PARAMS is off.
        got = _delta(JD_JUP, "jupiter", (60, 1, 0, 0, 0, 0), 0)
        assert got == pytest.approx(-0.391414, abs=AGE_TOL)

    def test_age_clamp_below_23(self):
        # Ages at/below ~23 share the age-23 delta (max dark-adapted pupil).
        d23 = _delta(JD_JUP, "jupiter", (23, 1, 0, 0, 0, 0), 0)
        d10 = _delta(JD_JUP, "jupiter", (10, 1, 0, 0, 0, 0), 0)
        d20 = _delta(JD_JUP, "jupiter", (20, 1, 0, 0, 0, 0), 0)
        assert d10 == d23
        assert d20 == d23
        assert d23 > 0.0  # younger than 36 sees fainter

    def test_age_monotonic_decreasing(self):
        deltas = [
            _delta(JD_JUP, "jupiter", (age, 1, 0, 0, 0, 0), 0)
            for age in range(23, 91, 1)
        ]
        assert all(b < a for a, b in zip(deltas, deltas[1:]))

    def test_age_zero_at_36(self):
        assert _delta(JD_JUP, "jupiter", (36, 1, 0, 0, 0, 0), 0) == 0.0

    def test_age_is_sky_dependent(self):
        # Same age, different sky brightness -> different delta (as the oracle).
        d_jup = _delta(JD_JUP, "jupiter", (60, 1, 0, 0, 0, 0), 0)
        d_sir = _delta(JD_SIR, "sirius", (60, 1, 0, 0, 0, 0), 0)
        assert d_jup != pytest.approx(d_sir, abs=0.05)


class TestOpticalAid:
    """Optical-instrument correction (only with HELFLAG_OPTICAL_PARAMS)."""

    # (observer, oracle delta) at the Jupiter scene. observer =
    # (age, snellen, is_binocular, magnification, aperture_mm, transmission).
    OPTICS_FACTS = [
        ((36, 1, 1, 50, 200, 0.8), 8.793116),
        ((36, 1, 1, 50, 200, 1.0), 9.241590),
        ((36, 1, 1, 100, 200, 0.8), 7.661777),
        ((36, 1, 1, 50, 100, 0.8), 6.156627),
        ((36, 1, 1, 10, 50, 0.8), 6.191654),
        ((36, 1, 0, 50, 200, 0.8), 8.096221),  # monocular
        ((36, 1, 1, 7, 50, 1.0), 7.239182),
    ]
    # Configs the fitted model reproduces tightly (binocular, realistic exit
    # pupils). Monocular and small-scope configs sit at the ~1 mag floor.
    OPTICS_TIGHT = {
        (36, 1, 1, 50, 200, 0.8),
        (36, 1, 1, 100, 200, 0.8),
        (36, 1, 1, 50, 100, 0.8),
    }

    @pytest.mark.parametrize("obs,expected", OPTICS_FACTS)
    def test_optics_delta_matches_oracle(self, obs, expected):
        got = _delta(JD_JUP, "jupiter", obs, HELFLAG_OPTICAL_PARAMS)
        tol = OPTICS_TOL_TIGHT if obs in self.OPTICS_TIGHT else OPTICS_TOL
        assert got == pytest.approx(expected, abs=tol)

    @pytest.mark.parametrize("obs,expected", OPTICS_FACTS)
    def test_optics_gives_large_positive_gain(self, obs, expected):
        got = _delta(JD_JUP, "jupiter", obs, HELFLAG_OPTICAL_PARAMS)
        assert got > 4.0  # every instrument here is a big gain

    def test_optics_ignored_without_flag(self):
        # Same optical observer, flag OFF -> optics are ignored (delta 0).
        got = _delta(JD_JUP, "jupiter", (36, 1, 1, 50, 200, 0.8), 0)
        assert got == 0.0

    def test_optics_ignored_when_magnification_not_meaningful(self):
        # A "1x, 7 mm" naked-eye-equivalent instrument is a no-op even with
        # the flag set (magnification <= 1).
        got = _delta(JD_JUP, "jupiter", (36, 1, 1, 1, 7, 1.0), HELFLAG_OPTICAL_PARAMS)
        assert got == 0.0

    def test_optics_transmission_monotonic(self):
        # Higher transmission -> fainter limit (larger delta).
        lo = _delta(JD_JUP, "jupiter", (36, 1, 1, 50, 200, 0.6), HELFLAG_OPTICAL_PARAMS)
        hi = _delta(JD_JUP, "jupiter", (36, 1, 1, 50, 200, 1.0), HELFLAG_OPTICAL_PARAMS)
        assert hi > lo

    def test_optics_aperture_monotonic(self):
        # Larger aperture (fixed magnification) -> fainter limit.
        small = _delta(
            JD_JUP, "jupiter", (36, 1, 1, 50, 100, 0.8), HELFLAG_OPTICAL_PARAMS
        )
        big = _delta(
            JD_JUP, "jupiter", (36, 1, 1, 50, 300, 0.8), HELFLAG_OPTICAL_PARAMS
        )
        assert big > small

    def test_age_not_combined_with_optics(self):
        # The reference ignores age once an instrument is used: age 80 + scope
        # equals age 36 + scope (Snellen still applies, here SN=1).
        old = _delta(
            JD_JUP, "jupiter", (80, 1, 1, 50, 200, 0.8), HELFLAG_OPTICAL_PARAMS
        )
        young = _delta(
            JD_JUP, "jupiter", (36, 1, 1, 50, 200, 0.8), HELFLAG_OPTICAL_PARAMS
        )
        assert old == pytest.approx(young, abs=1e-9)

    def test_snellen_additive_with_optics(self):
        # Snellen still applies on top of the optical gain.
        base = _delta(
            JD_JUP, "jupiter", (36, 1, 1, 50, 200, 0.8), HELFLAG_OPTICAL_PARAMS
        )
        with_sn = _delta(
            JD_JUP, "jupiter", (36, 1.2, 1, 50, 200, 0.8), HELFLAG_OPTICAL_PARAMS
        )
        assert with_sn - base == pytest.approx(5.0 * math.log10(1.2), abs=1e-9)


class TestBackendParity:
    """LEB and Skyfield twins must handle the observer identically."""

    def _both_backends(self, jd, body, obs, flags):
        import libephemeris.state as st

        leb = vis_limit_mag(jd, GEOPOS, ATMO, obs, body, flags)[1][0]
        orig = st.get_leb_reader
        st.get_leb_reader = lambda: None
        try:
            sky = vis_limit_mag(jd, GEOPOS, ATMO, obs, body, flags)[1][0]
        finally:
            st.get_leb_reader = orig
        return leb, sky

    def test_telescopic_call_backend_equal(self):
        leb, sky = self._both_backends(
            JD_JUP, "jupiter", (36, 1, 1, 50, 200, 0.8), HELFLAG_OPTICAL_PARAMS
        )
        assert leb == pytest.approx(sky, abs=1e-6)

    def test_age_call_backend_equal(self):
        leb, sky = self._both_backends(JD_JUP, "jupiter", (60, 1, 0, 0, 0, 0), 0)
        assert leb == pytest.approx(sky, abs=1e-6)
