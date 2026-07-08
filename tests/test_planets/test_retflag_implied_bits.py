"""Regression: retflag echoes the flag bits the reference API implies.

The reference's flag-plausibility step switches on (and echoes) flags the
request implies: J2000 and SIDEREAL output is referred to a mean equinox, so
FLG_NONUT is echoed; heliocentric/barycentric/true-position output skips
light deflection and annual aberration, so FLG_NOGDEFL | FLG_NOABERR are
echoed. Verified behaviourally against the reference oracle (Mars, JD
2451545.0): J2000 -> +NONUT, SIDEREAL -> +NONUT, HELCTR/TRUEPOS ->
+NOGDEFL|NOABERR, composing for combinations, including ECL_NUT.

Only the echo carries the implied bits; the computation flags are untouched
(adding NONUT to a sidereal request would change the ayanamsha realization).
"""

from __future__ import annotations

import pytest

import libephemeris as le
from libephemeris.constants import (
    ECL_NUT,
    FLG_EQUATORIAL,
    FLG_HELCTR,
    FLG_J2000,
    FLG_NOABERR,
    FLG_NOGDEFL,
    FLG_NONUT,
    FLG_SIDEREAL,
    FLG_TRUEPOS,
    MARS,
)

JD = 2451545.0

CASES = [
    (FLG_J2000, FLG_NONUT),
    (FLG_SIDEREAL, FLG_NONUT),
    (FLG_SIDEREAL | FLG_EQUATORIAL, FLG_NONUT),
    (FLG_HELCTR, FLG_NOGDEFL | FLG_NOABERR),
    (FLG_TRUEPOS, FLG_NOGDEFL | FLG_NOABERR),
    (FLG_HELCTR | FLG_J2000, FLG_NONUT | FLG_NOGDEFL | FLG_NOABERR),
]


class TestImpliedRetflagBits:
    @pytest.mark.unit
    @pytest.mark.parametrize("req,implied", CASES)
    def test_calc_ut_echoes_implied_bits(self, req: int, implied: int) -> None:
        _, retflag = le.calc_ut(JD, MARS, req)
        assert retflag & implied == implied
        # The requested bits themselves are still echoed too.
        assert retflag & req == req

    @pytest.mark.unit
    @pytest.mark.parametrize("req,implied", CASES)
    def test_calc_echoes_implied_bits(self, req: int, implied: int) -> None:
        _, retflag = le.calc(JD, MARS, req)
        assert retflag & implied == implied

    @pytest.mark.unit
    def test_no_spurious_bits_on_plain_requests(self) -> None:
        for req in (0, FLG_NONUT, FLG_EQUATORIAL):
            _, retflag = le.calc_ut(JD, MARS, req)
            assert not (retflag & FLG_NOGDEFL)
            assert not (retflag & FLG_NOABERR)
            if not req & FLG_NONUT:
                assert not (retflag & FLG_NONUT)

    @pytest.mark.unit
    def test_ecl_nut_echoes_implied_bits(self) -> None:
        _, retflag = le.calc_ut(JD, ECL_NUT, FLG_J2000)
        assert retflag & FLG_NONUT
        _, retflag = le.calc_ut(JD, ECL_NUT, FLG_HELCTR)
        assert retflag & (FLG_NOGDEFL | FLG_NOABERR) == FLG_NOGDEFL | FLG_NOABERR

    @pytest.mark.unit
    def test_implied_bits_are_computational_noops(self) -> None:
        """Why the reference can echo these bits: they are no-ops.

        Sidereal longitude is nutation-free by construction (the equinox
        nutation cancels between the position and the ayanamsha), so
        SIDEREAL == SIDEREAL|NONUT exactly -- which is precisely why the
        reference flags sidereal output as NONUT. The tropical longitudes,
        by contrast, must still differ by the nutation in longitude (sanity
        that the frame machinery is alive and the echo did not leak)."""
        le.set_sid_mode(le.SIDM_LAHIRI, 0.0, 0.0)
        sid_true, _ = le.calc_ut(JD, MARS, FLG_SIDEREAL)
        sid_mean, _ = le.calc_ut(JD, MARS, FLG_SIDEREAL | FLG_NONUT)
        assert sid_true[0] == pytest.approx(sid_mean[0], abs=1e-9)
        trop_true, _ = le.calc_ut(JD, MARS, 0)
        trop_mean, _ = le.calc_ut(JD, MARS, FLG_NONUT)
        delta = (trop_true[0] - trop_mean[0] + 180.0) % 360.0 - 180.0
        assert 1.0 < abs(delta) * 3600.0 < 60.0  # nutation-sized
