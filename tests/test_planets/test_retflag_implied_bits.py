"""Regression: retflag echoes the flag bits the reference API implies.

The reference's flag normalization switches on (and echoes) flags the
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
    MEAN_NODE,
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
    @pytest.mark.parametrize(
        "req,expect_speed,expect_speed3",
        [
            (le.FLG_SPEED3, False, True),
            (le.FLG_SPEED, True, False),
            (le.FLG_SPEED | le.FLG_SPEED3, True, False),  # SPEED wins
            (le.FLG_SIDEREAL | le.FLG_SPEED3, False, True),
        ],
    )
    def test_speed3_echo(self, req, expect_speed, expect_speed3) -> None:
        """The reference echoes SPEED3 when SPEED3 alone was requested and
        SPEED when SPEED is present (SPEED wins); the computed values are the
        same either way (identical 3-position differentiation)."""
        from libephemeris.constants import FLG_SPEED, FLG_SPEED3

        for fn in (le.calc_ut, le.calc):
            _, rf = fn(JD, MARS, req)
            assert bool(rf & FLG_SPEED) is expect_speed
            assert bool(rf & FLG_SPEED3) is expect_speed3
        # South node, context and the LEB backend echo the same.
        _, rf = le.calc_ut(JD, -MEAN_NODE, req)
        assert bool(rf & FLG_SPEED3) is expect_speed3
        _, rf = le.EphemerisContext().calc_ut(JD, MARS, req)
        assert bool(rf & FLG_SPEED3) is expect_speed3
        # Values identical to the plain-SPEED request.
        pos3, _ = le.calc_ut(JD, MARS, req)
        pos_s, _ = le.calc_ut(JD, MARS, (req & ~FLG_SPEED3) | FLG_SPEED)
        assert pos3 == pos_s

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


@pytest.mark.unit
class TestMosephEchoRegression:
    """MOSEPH-only requests echo FLG_MOSEPH in the retflag (measured
    black-box: MOSEPH -> 4, MOSEPH|SPEED -> 260, MOSEPH|SPEED3 -> 132,
    MOSEPH|SWIEPH -> 2 with SWIEPH winning). The computation still uses the
    DE440/DE441 path — only the echo changes."""

    def test_moseph_only_echoes_moseph(self):
        from libephemeris.constants import FLG_MOSEPH

        _, rf = le.calc_ut(JD, MARS, FLG_MOSEPH)
        assert rf == FLG_MOSEPH

    def test_moseph_speed_matrix(self):
        from libephemeris.constants import FLG_MOSEPH, FLG_SPEED, FLG_SPEED3

        assert le.calc_ut(JD, MARS, FLG_MOSEPH | FLG_SPEED)[1] == (
            FLG_MOSEPH | FLG_SPEED
        )
        assert le.calc_ut(JD, MARS, FLG_MOSEPH | FLG_SPEED3)[1] == (
            FLG_MOSEPH | FLG_SPEED3
        )

    def test_swieph_beats_moseph(self):
        from libephemeris.constants import FLG_MOSEPH, FLG_SWIEPH

        _, rf = le.calc_ut(JD, MARS, FLG_MOSEPH | FLG_SWIEPH)
        assert rf == FLG_SWIEPH

    def test_moseph_position_unchanged(self):
        """MOSEPH must not change the computed position (DE path always)."""
        from libephemeris.constants import FLG_MOSEPH

        p_mos, _ = le.calc_ut(JD, MARS, FLG_MOSEPH)
        p_def, _ = le.calc_ut(JD, MARS, 0)
        assert p_mos == p_def


@pytest.mark.unit
class TestFictitiousSpeedGate:
    """Fictitious/predicted bodies (49-58) return zero speed slots when
    FLG_SPEED is absent, like every other body class (measured black-box:
    the reference returns zeros)."""

    def test_no_speed_returns_zeros(self):
        for ipl in (49, 50, 51, 55, 56, 57, 58):
            r, _ = le.calc_ut(JD, ipl, 0)
            assert r[3:] == (0.0, 0.0, 0.0), (ipl, r[3:])

    def test_speed_flag_still_computes(self):
        from libephemeris.constants import FLG_SPEED

        r, _ = le.calc_ut(JD, 56, FLG_SPEED)
        assert r[3] != 0.0

    def test_helio_no_speed_zeros(self):
        from libephemeris.constants import FLG_HELCTR

        r, _ = le.calc_ut(JD, 55, FLG_HELCTR)
        assert r[3:] == (0.0, 0.0, 0.0)


class TestEclNutRetflagEphemerisBits:
    """ECL_NUT (-1) ephemeris-bit echo, measured on the reference API:
    calc() (TT) echoes only the ephemeris bits the caller passed (no forced
    default SWIEPH: calc(jd,-1,0) echoes 0), while calc_ut() adds the
    default (calc_ut(jd,-1,0) echoes 2)."""

    def test_calc_tt_no_forced_default(self):
        for f, want in [(0, 0), (256, 256), (128, 128), (64, 64), (32, 96)]:
            assert le.calc(JD, -1, f)[1] == want, f

    def test_calc_tt_explicit_bits_kept(self):
        for f, want in [(2, 2), (4, 4), (258, 258)]:
            assert le.calc(JD, -1, f)[1] == want, f

    def test_calc_ut_adds_default(self):
        for f, want in [(0, 2), (256, 258), (128, 130), (64, 66), (32, 98)]:
            assert le.calc_ut(JD, -1, f)[1] == want, f
