"""Independent physical invariants for ``vis_limit_mag`` observer effects.

The observer model follows Schaefer (1990, PASP 102, 212): dark-adapted pupil
diameter controls naked-eye light grasp, while a telescope is governed by its
entrance aperture, magnification/exit pupil, transmission, and sky luminance.
No expected value in this module was captured from another ephemeris.
"""

from __future__ import annotations

import math

import pytest

from libephemeris import julday
from libephemeris.constants import HELFLAG_OPTICAL_PARAMS
from libephemeris.heliacal import vis_limit_mag

GEOPOS = (12.5, 41.9, 0.0)
ATMO = (1013.25, 15.0, 40.0, 0.0)

# Two ordinary round-number calendar instants chosen independently to exercise
# bright and dark sky paths.
JD_TWILIGHT = julday(2024, 3, 20, 18.0)
JD_NIGHT = julday(2025, 1, 10, 18.0)


def _delta(jd, body, obs, flags):
    """Return the observer correction against the age-36 naked baseline."""
    base = vis_limit_mag(jd, GEOPOS, ATMO, (36, 1, 0, 0, 0, 0), body, 0)[1][0]
    got = vis_limit_mag(jd, GEOPOS, ATMO, obs, body, flags)[1][0]
    return got - base


def _published_pupil(age: float) -> float:
    """Schaefer dark-adapted pupil diameter, with its age-23 cap."""
    capped_age = max(age, 23.0)
    return 7.0 * math.exp(-(capped_age**2) / 20000.0)


class TestNakedBaseline:
    """The default observer is a true no-op and produces finite output."""

    @pytest.mark.parametrize(
        "jd,body", [(JD_TWILIGHT, "jupiter"), (JD_NIGHT, "sirius")]
    )
    def test_naked_default_is_finite_and_repeatable(self, jd, body):
        first = vis_limit_mag(jd, GEOPOS, ATMO, (36, 1, 0, 0, 0, 0), body, 0)
        second = vis_limit_mag(jd, GEOPOS, ATMO, (36, 1, 0, 0, 0, 0), body, 0)
        assert math.isfinite(first[1][0])
        assert first == second
        assert _delta(jd, body, (36, 1, 0, 0, 0, 0), 0) == 0.0


class TestSnellen:
    """Snellen acuity contributes the published ``5*log10(SN)`` term."""

    @pytest.mark.parametrize("sn", [0.5, 1.0, 1.5, 2.0])
    @pytest.mark.parametrize(
        "jd,body", [(JD_TWILIGHT, "jupiter"), (JD_NIGHT, "sirius")]
    )
    def test_snellen_term(self, sn, jd, body):
        got = _delta(jd, body, (36, sn, 0, 0, 0, 0), 0)
        assert got == pytest.approx(5.0 * math.log10(sn), abs=1e-9)


class TestAgeCurve:
    """Age follows pupil collecting area, without an output-fitted curve."""

    @pytest.mark.parametrize("age", [10, 20, 23, 25, 30, 40, 50, 60, 70, 80, 90])
    def test_age_delta_is_pupil_area_magnitude(self, age):
        got = _delta(JD_TWILIGHT, "jupiter", (age, 1, 0, 0, 0, 0), 0)
        expected = 5.0 * math.log10(_published_pupil(age) / _published_pupil(36))
        assert got == pytest.approx(expected, abs=1e-12)

    def test_age_clamp_below_23(self):
        d23 = _delta(JD_TWILIGHT, "jupiter", (23, 1, 0, 0, 0, 0), 0)
        assert _delta(JD_TWILIGHT, "jupiter", (10, 1, 0, 0, 0, 0), 0) == d23
        assert _delta(JD_TWILIGHT, "jupiter", (20, 1, 0, 0, 0, 0), 0) == d23
        assert d23 > 0.0

    def test_age_monotonic_decreasing(self):
        deltas = [
            _delta(JD_TWILIGHT, "jupiter", (age, 1, 0, 0, 0, 0), 0)
            for age in range(23, 91)
        ]
        assert all(b < a for a, b in zip(deltas, deltas[1:]))

    def test_age_zero_at_reference_age(self):
        assert _delta(JD_TWILIGHT, "jupiter", (36, 1, 0, 0, 0, 0), 0) == 0.0

    def test_pupil_area_term_is_sky_independent(self):
        twilight = _delta(JD_TWILIGHT, "jupiter", (60, 1, 0, 0, 0, 0), 0)
        night = _delta(JD_NIGHT, "sirius", (60, 1, 0, 0, 0, 0), 0)
        assert twilight == pytest.approx(night, abs=1e-12)


class TestOpticalAid:
    """Telescope behavior follows light grasp and exit-pupil invariants."""

    CONFIGS = [
        (36, 1, 1, 50, 200, 0.8),
        (36, 1, 1, 100, 200, 0.8),
        (36, 1, 1, 50, 100, 0.8),
        (36, 1, 1, 10, 50, 0.8),
        (36, 1, 0, 50, 200, 0.8),
        (36, 1, 1, 7, 50, 1.0),
    ]

    @pytest.mark.parametrize("obs", CONFIGS)
    def test_real_instruments_improve_limit(self, obs):
        got = _delta(JD_TWILIGHT, "jupiter", obs, HELFLAG_OPTICAL_PARAMS)
        assert got > 0.0
        assert math.isfinite(got)

    def test_optics_ignored_without_flag(self):
        got = _delta(JD_TWILIGHT, "jupiter", (36, 1, 1, 50, 200, 0.8), 0)
        assert got == 0.0

    def test_optics_ignored_at_unit_magnification(self):
        got = _delta(
            JD_TWILIGHT,
            "jupiter",
            (36, 1, 1, 1, 7, 1.0),
            HELFLAG_OPTICAL_PARAMS,
        )
        assert got == 0.0

    def test_transmission_monotonic(self):
        lo = _delta(
            JD_TWILIGHT,
            "jupiter",
            (36, 1, 1, 50, 200, 0.6),
            HELFLAG_OPTICAL_PARAMS,
        )
        hi = _delta(
            JD_TWILIGHT,
            "jupiter",
            (36, 1, 1, 50, 200, 1.0),
            HELFLAG_OPTICAL_PARAMS,
        )
        assert hi > lo

    def test_aperture_monotonic_before_pupil_limit(self):
        small = _delta(
            JD_TWILIGHT,
            "jupiter",
            (36, 1, 0, 50, 100, 0.8),
            HELFLAG_OPTICAL_PARAMS,
        )
        big = _delta(
            JD_TWILIGHT,
            "jupiter",
            (36, 1, 0, 50, 200, 0.8),
            HELFLAG_OPTICAL_PARAMS,
        )
        assert big > small

    def test_binocular_summation_is_sqrt_two_snr(self):
        mono = _delta(
            JD_TWILIGHT,
            "jupiter",
            (36, 1, 0, 50, 200, 0.8),
            HELFLAG_OPTICAL_PARAMS,
        )
        bino = _delta(
            JD_TWILIGHT,
            "jupiter",
            (36, 1, 1, 50, 200, 0.8),
            HELFLAG_OPTICAL_PARAMS,
        )
        assert bino - mono == pytest.approx(1.25 * math.log10(2.0), abs=1e-12)

    def test_age_limits_effective_aperture(self):
        young = _delta(
            JD_TWILIGHT,
            "jupiter",
            (23, 1, 0, 10, 100, 0.8),
            HELFLAG_OPTICAL_PARAMS,
        )
        old = _delta(
            JD_TWILIGHT,
            "jupiter",
            (80, 1, 0, 10, 100, 0.8),
            HELFLAG_OPTICAL_PARAMS,
        )
        assert young > old

    def test_snellen_additive_with_optics(self):
        base = _delta(
            JD_TWILIGHT,
            "jupiter",
            (36, 1, 1, 50, 200, 0.8),
            HELFLAG_OPTICAL_PARAMS,
        )
        with_sn = _delta(
            JD_TWILIGHT,
            "jupiter",
            (36, 1.2, 1, 50, 200, 0.8),
            HELFLAG_OPTICAL_PARAMS,
        )
        assert with_sn - base == pytest.approx(5.0 * math.log10(1.2), abs=1e-9)


class TestBackendParity:
    """LEB and Skyfield twins must handle the observer identically."""

    @staticmethod
    def _both_backends(jd, body, obs, flags):
        import libephemeris.state as st

        leb = vis_limit_mag(jd, GEOPOS, ATMO, obs, body, flags)[1][0]
        orig = st.get_leb_reader
        st.get_leb_reader = lambda: None
        try:
            sky = vis_limit_mag(jd, GEOPOS, ATMO, obs, body, flags)[1][0]
        finally:
            st.get_leb_reader = orig
        return leb, sky

    @pytest.mark.parametrize(
        "obs,flags",
        [
            ((36, 1, 1, 50, 200, 0.8), HELFLAG_OPTICAL_PARAMS),
            ((60, 1, 0, 0, 0, 0), 0),
        ],
    )
    def test_backend_equal(self, obs, flags):
        leb, sky = self._both_backends(JD_TWILIGHT, "jupiter", obs, flags)
        assert leb == pytest.approx(sky, abs=1e-6)
