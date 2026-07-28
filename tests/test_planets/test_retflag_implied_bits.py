"""Public retflag normalization and implied-bit invariants.

Flag normalization switches on (and echoes) flags the request implies: J2000
and SIDEREAL output is referred to a mean equinox, so
FLG_NONUT is echoed; heliocentric/barycentric/true-position output skips
light deflection and annual aberration, so FLG_NOGDEFL | FLG_NOABERR are
echoed. The relations compose for combinations, including ECL_NUT.

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
        """The API echoes SPEED3 when SPEED3 alone was requested and
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
        SIDEREAL and SIDEREAL|NONUT agree to floating-point precision --
        which is precisely why the reference flags sidereal output as NONUT.
        The tropical longitudes, by contrast, must still differ by the
        nutation in longitude (sanity that the frame machinery is alive and
        the echo did not leak)."""
        le.set_sid_mode(le.SIDM_LAHIRI, 0.0, 0.0)
        sid_true, _ = le.calc_ut(JD, MARS, FLG_SIDEREAL)
        sid_mean, _ = le.calc_ut(JD, MARS, FLG_SIDEREAL | FLG_NONUT)
        # The Skyfield and LEB reduction paths reassociate the cancelling
        # rotations slightly differently (a few microarcseconds at J2000).
        assert sid_true[0] == pytest.approx(sid_mean[0], abs=1e-8)
        trop_true, _ = le.calc_ut(JD, MARS, 0)
        trop_mean, _ = le.calc_ut(JD, MARS, FLG_NONUT)
        delta = (trop_true[0] - trop_mean[0] + 180.0) % 360.0 - 180.0
        assert 1.0 < abs(delta) * 3600.0 < 60.0  # nutation-sized


@pytest.mark.unit
class TestMosephEchoRegression:
    """MOSEPH-only requests echo named flags; SWIEPH wins conflicts."""

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
    """The reviewed Harrington body gates speed slots behind FLG_SPEED."""

    def test_no_speed_returns_zeros(self):
        r, _ = le.calc_ut(JD, 50, 0)
        assert r[3:] == (0.0, 0.0, 0.0)

    def test_speed_flag_still_computes(self):
        from libephemeris.constants import FLG_SPEED

        r, _ = le.calc_ut(JD, 50, FLG_SPEED)
        assert r[3] != 0.0

    def test_helio_no_speed_zeros(self):
        from libephemeris.constants import FLG_HELCTR

        r, _ = le.calc_ut(JD, 50, FLG_HELCTR)
        assert r[3:] == (0.0, 0.0, 0.0)


class TestEclNutRetflagEphemerisBits:
    """ECL_NUT retflags: calc echoes the input flags, calc_ut adds a source.

    Measured reference behavior: the TT entry point returns the request
    unchanged (apart from the SPEED3->SPEED implication), while the UT entry
    point resolves a default ephemeris bit when none is requested.
    """

    def test_calc_tt_echoes_flags(self):
        for f, want in [(0, 0), (256, 256), (128, 128), (64, 64), (32, 96)]:
            assert le.calc(JD, -1, f)[1] == want, f

    def test_calc_tt_explicit_bits_kept(self):
        for f, want in [(2, 2), (4, 4), (258, 258)]:
            assert le.calc(JD, -1, f)[1] == want, f

    def test_calc_ut_adds_default(self):
        for f, want in [(0, 2), (256, 258), (128, 130), (64, 66), (32, 98)]:
            assert le.calc_ut(JD, -1, f)[1] == want, f


class TestCalcTTNormalBodyEphemerisBits:
    """TT echoes caller source bits; UT adds the default source bit."""

    def test_no_forced_default_matrix(self):
        le.set_topo(0.0, 0.0, 0.0)
        for f, want in [
            (0, 0),
            (256, 256),
            (8, 1544),
            (16384, 17920),
            (32768, 32768),
            (16, 1552),
            (64, 64),
            (32, 96),
            (65536, 65600),
        ]:
            assert le.calc(JD, 2, f)[1] == want, f

    def test_explicit_bits_kept(self):
        for f, want in [(2, 2), (4, 4), (258, 258)]:
            assert le.calc(JD, 2, f)[1] == want, f

    def test_calc_ut_keeps_default(self):
        assert le.calc_ut(JD, 2, 0)[1] == 2
        assert le.calc_ut(JD, 2, 256)[1] == 258

    def test_south_node_and_leb_paths(self):
        assert le.calc(JD, -10, 0)[1] == 0
        for mode in ("skyfield", "leb"):
            le.set_calc_mode(mode)
            assert le.calc(JD, 2, 0)[1] == 0, mode
        le.set_calc_mode("auto")


class TestCenterBodyRetflagEcho:
    """FLG_CENTER_BODY is consumed for Sun..Mars, echoed for other bodies.

    Measured reference behavior: for ipl 0-4 there is no satellite-system
    barycenter to resolve, and the retflag does NOT carry the bit back
    (calc and calc_ut alike); nodes, apogees and asteroids echo it
    unchanged. Backend-independent (LEB and Skyfield agree).
    """

    CB = 1048576  # FLG_CENTER_BODY

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [0, 1, 2, 3, 4])
    def test_stripped_for_sun_through_mars(self, body: int) -> None:
        assert le.calc_ut(JD, body, 2 | self.CB)[1] == 2
        assert le.calc(JD, body, 2 | self.CB)[1] == 2

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [10, 11, 12, 15, 17])
    def test_echoed_for_nodes_apogees_asteroids(self, body: int) -> None:
        assert le.calc_ut(JD, body, 2 | self.CB)[1] == 2 | self.CB
        assert le.calc(JD, body, 2 | self.CB)[1] == 2 | self.CB

    @pytest.mark.unit
    def test_backend_independent(self) -> None:
        for mode in ("skyfield", "leb"):
            le.set_calc_mode(mode)
            assert le.calc_ut(JD, 0, 2 | self.CB)[1] == 2, mode
            assert le.calc_ut(JD, 15, 2 | self.CB)[1] == 2 | self.CB, mode
        le.set_calc_mode("auto")


class TestJplhorRetflagConsumed:
    """FLG_JPLHOR / FLG_JPLHOR_APPROX are consumed, never echoed.

    This library performs no JPL-Horizons dpsi/deps Earth-orientation
    reduction (the flags are accepted for API compatibility only), so the
    measured retflag convention is mirrored: the bits are stripped and the
    position equals the plain request on both time arguments and backends.
    """

    @pytest.mark.unit
    @pytest.mark.parametrize("body", [0, 2, 4])
    def test_bits_stripped_and_position_plain(self, body: int) -> None:
        from libephemeris.constants import FLG_JPLHOR, FLG_JPLHOR_APPROX

        for bit in (FLG_JPLHOR, FLG_JPLHOR_APPROX):
            pos, rf = le.calc_ut(JD, body, 2 | bit)
            plain, rf_plain = le.calc_ut(JD, body, 2)
            assert rf == rf_plain == 2
            assert pos[0] == pytest.approx(plain[0], abs=1e-12)
            pos_tt, rf_tt = le.calc(JD, body, 2 | bit)
            assert rf_tt == 2
            assert pos_tt[0] == pytest.approx(le.calc(JD, body, 2)[0][0], abs=1e-12)
