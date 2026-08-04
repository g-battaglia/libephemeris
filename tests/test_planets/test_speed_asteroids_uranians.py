"""Independent speed checks for asteroids and sourced fictitious bodies."""

from __future__ import annotations


import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    CHIRON,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_TRUEPOS,
    FLG_HELCTR,
    FLG_J2000,
    NODBIT_OSCU,
)
from libephemeris.exceptions import UnknownBodyError


@pytest.fixture(autouse=True)
def _reset_state():
    yield
    swe.close()


JD_J2000 = 2451545.0
DT = 0.01  # 0.01 days for numerical derivative

REVIEWED_FICTITIOUS_IDS = (*range(40, 49), *range(50, 59))
UNVERIFIED_FICTITIOUS_IDS = (49,)


def _numerical_speed(body, jd, flags, component=0):
    """Compute numerical derivative via central differences."""
    pos_m, _ = swe.calc_ut(jd - DT, body, flags)
    pos_p, _ = swe.calc_ut(jd + DT, body, flags)
    val_m = pos_m[component]
    val_p = pos_p[component]

    if component == 0:
        diff = val_p - val_m
        if diff > 180.0:
            diff -= 360.0
        elif diff < -180.0:
            diff += 360.0
        return diff / (2 * DT)

    return (val_p - val_m) / (2 * DT)


class TestAsteroidSpeedConsistency:
    """Test speed consistency for main asteroids."""

    ASTEROIDS = [CHIRON, CERES, PALLAS, JUNO, VESTA]

    @pytest.mark.unit
    @pytest.mark.parametrize("body", ASTEROIDS)
    def test_lon_speed(self, body):
        """Asteroid longitude speed matches numerical derivative."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        returned = pos[3]
        numerical = _numerical_speed(body, JD_J2000, flags, 0)
        assert returned == pytest.approx(numerical, abs=0.02), (
            f"Body {body}: returned {returned}, numerical {numerical}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("body", ASTEROIDS)
    def test_lat_speed(self, body):
        """Asteroid latitude speed matches numerical derivative."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        returned = pos[4]
        numerical = _numerical_speed(body, JD_J2000, flags, 1)
        assert returned == pytest.approx(numerical, abs=0.02), (
            f"Body {body}: lat speed returned {returned}, numerical {numerical}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("body", ASTEROIDS)
    def test_dist_speed(self, body):
        """Asteroid distance speed matches numerical derivative."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        returned = pos[5]
        numerical = _numerical_speed(body, JD_J2000, flags, 2)
        assert returned == pytest.approx(numerical, abs=1e-3), (
            f"Body {body}: dist speed returned {returned}, numerical {numerical}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("body", ASTEROIDS)
    def test_speed_magnitude_reasonable(self, body):
        """Asteroid speed magnitude is within physical limits."""
        flags = FLG_SWIEPH | FLG_SPEED
        pos, _ = swe.calc_ut(JD_J2000, body, flags)
        # Asteroids move slower than inner planets
        assert abs(pos[3]) < 1.0, f"Body {body}: lon speed {pos[3]} too fast"


class TestAsteroidTruePosSpeedSelfConsistency:
    """FLG_TRUEPOS asteroid speed is the true derivative of the served place.

    Regression guard for the SPK type-21 geometric-place speed on the Skyfield
    backend. The type-21 difference-table evaluation carries ~4e-7" of
    per-sample jitter; a one-second central difference amplified that into a
    ~0.007..0.02"/day longitude-speed bias, so the reported geometric speed
    stopped matching the body's own geometric positions (and drove the
    osculating nod_aps perihelion up to ~115" off the reference, versus the
    self-consistent LEB backend). The reported speed must equal a clean,
    wide-step central difference of the served TRUEPOS longitudes to well under
    0.001"/day, in every heliocentric/geocentric and of-date/J2000 frame.
    """

    ASTEROIDS = [CHIRON, CERES, PALLAS, JUNO, VESTA]
    # Whole-day UT epochs comfortably inside the 1920-2080 asteroid SPK safe
    # range (and away from any download-window edge).
    EPOCHS = [2435000.5, JD_J2000, 2455000.5]
    FRAMES = [0, FLG_HELCTR, FLG_J2000, FLG_HELCTR | FLG_J2000]

    # Ephemeris-grade sources the self-consistency bound applies to. The
    # multi-epoch Keplerian fallback (used when no SPK is provisioned, e.g.
    # offline CI) is explicitly NOT ephemeris-grade -- its piecewise anchors
    # are not smooth enough for a <0.001"/day derivative check -- so the guard
    # skips it rather than asserting a bound it never promised.
    _EPHEMERIS_GRADE = {"SPK", "LEB", "ASSIST", "Horizons"}

    @pytest.mark.unit
    @pytest.mark.parametrize("body", ASTEROIDS)
    @pytest.mark.parametrize("frame", FRAMES)
    @pytest.mark.parametrize("epoch", EPOCHS)
    def test_truepos_speed_matches_own_positions(self, body, frame, epoch):
        # The bug lived on the Skyfield type-21 path; LEB was always
        # self-consistent. Pin the backend so the guard exercises the fix, and
        # allow the SPK to be provisioned (the suite disables auto-download by
        # default). Where no SPK can be obtained (offline CI) the body falls to
        # the non-ephemeris-grade Keplerian model and the case is skipped.
        swe.set_calc_mode("skyfield")
        swe.set_auto_spk_download(True)
        base = FLG_SWIEPH | FLG_TRUEPOS | frame
        token = swe.start_tracing()
        try:
            reported = swe.calc_ut(epoch, body, base | FLG_SPEED)[0][3]
            source = swe.get_trace_results().get(body)
        finally:
            token.var.reset(token)
        if source not in self._EPHEMERIS_GRADE:
            pytest.skip(
                f"body {body} @ {epoch} served by {source}, not ephemeris-grade"
            )
        # Wide half-step: immune to the type-21 sample jitter that the buggy
        # one-second stencil amplified, so this reference derivative is clean.
        h = 0.05
        lon_p = swe.calc_ut(epoch + h, body, base)[0][0]
        lon_m = swe.calc_ut(epoch - h, body, base)[0][0]
        derivative = ((lon_p - lon_m + 180.0) % 360.0 - 180.0) / (2.0 * h)
        err_arcsec_per_day = abs(reported - derivative) * 3600.0
        assert err_arcsec_per_day < 0.001, (
            f"Body {body} frame {frame} @ {epoch}: reported {reported}, "
            f"own-position derivative {derivative} "
            f'({err_arcsec_per_day:.5f}"/day)'
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("body", ASTEROIDS)
    def test_osculating_perihelion_speed_is_own_derivative(self, body):
        """Osculating nod_aps perihelion speed differentiates its own longitude.

        The FLG_SPEED perihelion-longitude rate must be the central difference
        of the reported perihelion longitude. This is the invariant the
        geometric-speed bias broke (the osculating orbit is built from the
        body's geometric state, whose longitude speed was inconsistent).
        """
        swe.set_calc_mode("skyfield")
        swe.set_auto_spk_download(True)
        flags = FLG_SWIEPH | FLG_SPEED
        _, _, peri, _ = swe.nod_aps_ut(JD_J2000, body, NODBIT_OSCU, flags)
        reported = peri[3]
        h = 0.05
        _, _, peri_p, _ = swe.nod_aps_ut(JD_J2000 + h, body, NODBIT_OSCU, FLG_SWIEPH)
        _, _, peri_m, _ = swe.nod_aps_ut(JD_J2000 - h, body, NODBIT_OSCU, FLG_SWIEPH)
        derivative = ((peri_p[0] - peri_m[0] + 180.0) % 360.0 - 180.0) / (2.0 * h)
        assert reported == pytest.approx(derivative, abs=0.02), (
            f"Body {body}: perihelion-lon speed {reported} vs own derivative "
            f"{derivative}"
        )


class TestFictitiousBodies:
    """Only independently sourced built-ins enter the speed pipeline."""

    @pytest.mark.unit
    def test_harrington_speed_request_returns_finite_state(self) -> None:
        swe.set_calc_mode("skyfield")
        position, _ = swe.calc_ut(JD_J2000, swe.HARRINGTON, FLG_SWIEPH | FLG_SPEED)
        assert len(position) == 6
        assert all(math.isfinite(value) for value in position)

    @pytest.mark.unit
    @pytest.mark.parametrize("body", REVIEWED_FICTITIOUS_IDS)
    def test_reviewed_body_speed_request_returns_finite_state(self, body: int) -> None:
        swe.set_calc_mode("skyfield")
        position, _ = swe.calc_ut(JD_J2000, body, FLG_SWIEPH | FLG_SPEED)
        assert len(position) == 6
        assert all(math.isfinite(value) for value in position)

    @pytest.mark.unit
    @pytest.mark.parametrize("body", UNVERIFIED_FICTITIOUS_IDS)
    def test_unverified_body_speed_request_fails_closed(self, body: int) -> None:
        swe.set_calc_mode("skyfield")
        with pytest.raises(UnknownBodyError) as raised:
            swe.calc_ut(JD_J2000, body, FLG_SWIEPH | FLG_SPEED)
        assert raised.value.body_id == body
