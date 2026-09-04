"""SPK type 21: usable coverage, center chain and loud in-coverage failures.

Regression tests for the type-21 path (GitHub issue #49):

* an epoch inside the kernel's usable coverage is either served by the
  kernel or raises ``SPKEvaluationError`` -- it is never answered from the
  Keplerian approximation behind the caller's back, on any of the four
  dispatch paths (planet pipeline, direct path, heliocentric minor-body
  helpers, auto-download);
* the usable coverage is asymmetric: light-time band plus speed stencil at
  the start, the stencil alone at the end, and every epoch inside the span
  reported by ``get_spk_coverage()`` is served;
* the state chain honours the center the kernel declares (barycenter, Sun,
  or a body of the base ephemeris) and refuses a center it cannot resolve;
* ``EphemerisRangeError`` (outside the usable coverage) keeps the documented
  fallback at every caller;
* the gate and the evaluation both use TDB, the reader's scale.

The kernel is a synthetic reader with a hand-built linear trajectory, so
every expected position and velocity is known in closed form. The reader
double enforces the same bounds as the vendored reader: inclusive start,
exclusive end, one (target, center) pair.
"""

from __future__ import annotations

import math
import os

import numpy as np
import pytest

import libephemeris as eph
from libephemeris import minor_bodies, planets, spk, spk_auto, state, tracing
from libephemeris.constants import (
    AST_OFFSET,
    CHIRON,
    FLG_EQUATORIAL,
    FLG_ICRS,
    FLG_J2000,
    FLG_NOABERR,
    FLG_NOGDEFL,
    FLG_SPEED,
    FLG_TRUEPOS,
)
from libephemeris.exceptions import (
    CalculationError,
    EphemerisRangeError,
    SPKEvaluationError,
    SPKRequiredError,
)

NAIF = 20002060
J2000_JD = 2451545.0
STORED_START = 2458000.0  # TDB
STORED_END = 2459000.0  # TDB (exclusive, as in the reader)
AU_KM = 149597870.7
C_AU_DAY = 173.1446326846693
OBLIQUITY_J2000 = math.radians(23.4392911)

# Hand-built trajectory relative to the kernel center (ICRS): ~13 AU at J2000,
# ~24-25 AU over the stored span, ~5 km/s. Linear in time so the closed-form
# derivative is exact.
P0_AU = np.array([12.0, -5.0, 3.0])
V_AU_DAY = np.array([0.0012, 0.0026, 0.0006])


def relative_state(jd_tdb: float) -> tuple[np.ndarray, np.ndarray]:
    """Hand-built state relative to the center (AU, AU/day, ICRS)."""
    return P0_AU + V_AU_DAY * (jd_tdb - J2000_JD), V_AU_DAY.copy()


class _Segment:
    def __init__(self, target: int, center: int, start_jd: float, end_jd: float):
        self.target = target
        self.center = center
        self.data_type = 21
        self.frame = 1
        self.start_jd = start_jd
        self.end_jd = end_jd
        self.start_second = (start_jd - J2000_JD) * 86400.0
        self.end_second = (end_jd - J2000_JD) * 86400.0


class _SyntheticKernel:
    """Reader double: linear motion, the vendored reader's own bounds."""

    def __init__(self, center: int = 10, segments=None, fault_jd: float | None = None):
        self.segments = (
            list(segments)
            if segments is not None
            else [_Segment(NAIF, center, STORED_START, STORED_END)]
        )
        self.fault_jd = fault_jd
        self.evaluated: list[float] = []
        self.closed = False

    def compute_type21(self, center: int, target: int, jd: float, jd2: float = 0.0):
        # Two-part JD like the vendored reader: the fraction keeps the
        # sub-ULP spacing of the speed stencil.
        matched = [
            s for s in self.segments if s.target == target and s.center == center
        ]
        if not matched:
            raise ValueError("Invalid Target and/or Center")
        jd_total = jd + jd2
        if jd_total < min(s.start_jd for s in matched) or jd_total >= max(
            s.end_jd for s in matched
        ):
            raise ValueError("Invalid Time to evaluate")
        if self.fault_jd is not None and abs(jd_total - self.fault_jd) < 0.01:
            raise RuntimeError("Invalid data")
        self.evaluated.append(jd_total)
        elapsed = (jd - J2000_JD) + jd2
        pos_au, vel_au_d = P0_AU + V_AU_DAY * elapsed, V_AU_DAY.copy()
        return (pos_au * AU_KM).tolist(), (vel_au_d * AU_KM / 86400.0).tolist()

    def close(self) -> None:
        self.closed = True


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _isolate(monkeypatch):
    """Skyfield mode, no downloads, no LEB, clean SPK registries."""
    saved_map = dict(state._SPK_BODY_MAP)
    saved_kernels = dict(state._SPK_TYPE21_KERNELS)
    prev_mode = state.get_calc_mode()
    monkeypatch.delenv("LIBEPHEMERIS_LEB", raising=False)
    monkeypatch.delenv("LIBEPHEMERIS_MODE", raising=False)
    monkeypatch.setattr(state, "_LEB_FILE", None)
    monkeypatch.setattr(state, "_LEB_READER", None)
    monkeypatch.setattr(state, "_discover_leb_file", lambda: None)
    eph.set_calc_mode("skyfield")
    eph.set_strict_precision(None)
    eph.set_auto_spk_download(False)
    yield
    state._SPK_BODY_MAP.clear()
    state._SPK_BODY_MAP.update(saved_map)
    state._SPK_TYPE21_KERNELS.clear()
    state._SPK_TYPE21_KERNELS.update(saved_kernels)
    eph.set_strict_precision(None)
    eph.set_auto_spk_download(None)
    eph.set_calc_mode(prev_mode)


@pytest.fixture
def synthetic(tmp_path, monkeypatch):
    """Factory: install a synthetic type-21 kernel and register it.

    Returns ``(path, kernel)``. The file on disk is a placeholder (the reader
    is served from the cache), so type detection is pinned to 21 for it.
    """
    original_detect = spk._detect_spk_type
    pinned: dict[str, int] = {}
    monkeypatch.setattr(
        spk, "_detect_spk_type", lambda p: pinned.get(p) or original_detect(p)
    )

    def _install(
        center: int = 10,
        *,
        ipl: int = CHIRON,
        naif_id: int = NAIF,
        name: str = "synthetic.bsp",
        register: bool = True,
        **kernel_kwargs,
    ):
        path = os.path.join(str(tmp_path), name)
        with open(path, "wb") as fh:
            fh.write(b"\x00" * 16)
        kernel = _SyntheticKernel(center, **kernel_kwargs)
        pinned[path] = 21
        state._SPK_TYPE21_KERNELS[path] = kernel
        if register:
            spk.register_spk_body(ipl, path, naif_id)
        return path, kernel

    return _install


def _ts():
    return state.get_timescale()


def _calc_tt(jd_tt: float, flags: int):
    """calc() with the dispatch source recorded."""
    tok = tracing.start_tracing()
    try:
        pos, _ = eph.calc(jd_tt, CHIRON, flags)
        return pos, tracing.get_trace_results().get(CHIRON)
    finally:
        tracing._trace_data.reset(tok)


def _center_state(center: int, t):
    """Barycentric ICRS state of a kernel center from the base ephemeris."""
    if center == 0:
        return np.zeros(3), np.zeros(3)
    body = state.get_planets()[center].at(t)
    return np.asarray(body.position.au), np.asarray(body.velocity.au_per_d)


def _geocentric_state(center: int, t):
    """Hand-built geometric geocentric ICRS state: center + relative - Earth."""
    c_pos, c_vel = _center_state(center, t)
    r_pos, r_vel = relative_state(float(t.tdb))
    e_pos, e_vel = _center_state(399, t)
    return c_pos + r_pos - e_pos, c_vel + r_vel - e_vel


def _to_ecliptic_j2000(v: np.ndarray) -> np.ndarray:
    ce, se = math.cos(OBLIQUITY_J2000), math.sin(OBLIQUITY_J2000)
    return np.array([v[0], v[1] * ce + v[2] * se, -v[1] * se + v[2] * ce])


def _spherical(v: np.ndarray) -> tuple[float, float, float]:
    r = float(np.linalg.norm(v))
    lon = math.degrees(math.atan2(v[1], v[0])) % 360.0
    lat = math.degrees(math.asin(v[2] / r))
    return lon, lat, r


def _spherical_rates(v: np.ndarray, dv: np.ndarray) -> tuple[float, float, float]:
    """Closed-form d(lon, lat, dist)/dt in deg/day, deg/day, AU/day."""
    x, y, z = v
    r = float(np.linalg.norm(v))
    xy_sq = x * x + y * y
    dlon = math.degrees((x * dv[1] - y * dv[0]) / xy_sq)
    dr = float(np.dot(v, dv)) / r
    dlat = math.degrees((dv[2] * r - z * dr) / (r * r * math.sqrt(1.0 - (z / r) ** 2)))
    return dlon, dlat, dr


GEOMETRIC_ICRS = (
    FLG_EQUATORIAL | FLG_ICRS | FLG_J2000 | FLG_TRUEPOS | FLG_NOABERR | FLG_NOGDEFL
)


# ---------------------------------------------------------------------------
# Registration validation
# ---------------------------------------------------------------------------


class TestRegistration:
    def test_unknown_naif_id_rejected(self, synthetic):
        with pytest.raises(ValueError, match="no type 21 segment"):
            synthetic(10, naif_id=NAIF + 1)

    def test_unresolvable_center_rejected(self, synthetic):
        with pytest.raises(ValueError, match="cannot be resolved"):
            synthetic(2000001)
        assert CHIRON not in state._SPK_BODY_MAP

    def test_center_absent_from_base_ephemeris_rejected(self, synthetic):
        # 599 is a planetary-system id but DE440 carries the Jupiter
        # barycenter (5) only.
        with pytest.raises(ValueError, match="not available in the base ephemeris"):
            synthetic(599)

    def test_mixed_centers_rejected(self, synthetic):
        segments = [
            _Segment(NAIF, 10, STORED_START, STORED_START + 500.0),
            _Segment(NAIF, 0, STORED_START + 500.0, STORED_END),
        ]
        with pytest.raises(ValueError, match="different centers"):
            synthetic(10, segments=segments)

    @pytest.mark.parametrize("center", [0, 10, 3, 399])
    def test_resolvable_centers_register(self, synthetic, center):
        path, _ = synthetic(center)
        assert state._SPK_BODY_MAP[CHIRON] == (path, NAIF)


# ---------------------------------------------------------------------------
# Center chain: resulting state, position and velocity
# ---------------------------------------------------------------------------


class TestCenterChain:
    @pytest.mark.parametrize("center", [0, 10, 399])
    def test_barycentric_state_from_declared_center(self, synthetic, center):
        synthetic(center)
        t = _ts().tt_jd(2458500.0)
        target = spk.get_spk_type21_target(CHIRON, t.tt)
        assert target is not None and target.segment_center == center
        icrf = target.at(t)
        c_pos, c_vel = _center_state(center, t)
        r_pos, r_vel = relative_state(float(t.tdb))
        np.testing.assert_allclose(icrf.position.au, c_pos + r_pos, atol=1e-12)
        np.testing.assert_allclose(icrf.velocity.au_per_d, c_vel + r_vel, atol=1e-12)

    @pytest.mark.parametrize("center", [0, 10, 399])
    def test_direct_path_geometric_position_and_speed(self, synthetic, center):
        synthetic(center)
        t = _ts().tt_jd(2458500.0)
        lon, lat, dist, dlon, dlat, ddist = spk.calc_spk_body_position(
            t, CHIRON, FLG_J2000 | FLG_TRUEPOS | FLG_SPEED
        )
        geo, geo_vel = _geocentric_state(center, t)
        exp_lon, exp_lat, exp_dist = _spherical(_to_ecliptic_j2000(geo))
        exp_dlon, exp_dlat, exp_ddist = _spherical_rates(
            _to_ecliptic_j2000(geo), _to_ecliptic_j2000(geo_vel)
        )
        assert lon == pytest.approx(exp_lon, abs=1e-9)
        assert lat == pytest.approx(exp_lat, abs=1e-9)
        assert dist == pytest.approx(exp_dist, abs=1e-12)
        # One-second central difference with exact spacing: the closed-form
        # rate is recovered to the rounding of the sampled coordinates.
        assert dlon == pytest.approx(exp_dlon, abs=1e-8)
        assert dlat == pytest.approx(exp_dlat, abs=1e-8)
        assert ddist == pytest.approx(exp_ddist, abs=1e-9)

    @pytest.mark.parametrize("center", [0, 10, 399])
    def test_planet_pipeline_geometric_position_and_speed(self, synthetic, center):
        synthetic(center)
        jd_tt = 2458500.0
        t = _ts().tt_jd(jd_tt)
        (ra, dec, dist, dra, ddec, ddist), source = _calc_tt(
            jd_tt, GEOMETRIC_ICRS | FLG_SPEED
        )
        assert source == "SPK"
        geo, geo_vel = _geocentric_state(center, t)
        exp_ra, exp_dec, exp_dist = _spherical(geo)
        exp_dra, exp_ddec, exp_ddist = _spherical_rates(geo, geo_vel)
        assert ra == pytest.approx(exp_ra, abs=1e-9)
        assert dec == pytest.approx(exp_dec, abs=1e-9)
        assert dist == pytest.approx(exp_dist, abs=1e-12)
        assert dra == pytest.approx(exp_dra, abs=1e-7)
        assert ddec == pytest.approx(exp_ddec, abs=1e-7)
        assert ddist == pytest.approx(exp_ddist, abs=1e-9)

    def test_center_zero_differs_from_sun_by_the_sun_offset(self, synthetic):
        """The same relative state under center 0 and center 10 must differ by
        exactly the Sun's barycentric offset -- the error the previous
        heliocentric arithmetic would have hidden."""
        t = _ts().tt_jd(2458500.0)
        states = {}
        for center in (0, 10):
            synthetic(center, name=f"c{center}.bsp")
            states[center] = spk.get_spk_type21_target(CHIRON, t.tt).at(t).position.au
        sun = state.get_planets()[10].at(t).position.au
        np.testing.assert_allclose(states[10] - states[0], sun, atol=1e-12)
        assert np.linalg.norm(sun) > 1e-3  # the offset is not negligible

    def test_runtime_unresolvable_center_is_loud(self, synthetic):
        """A center the base ephemeris lacks fails explicitly on every path."""
        path, _ = synthetic(599, register=False)
        state._SPK_BODY_MAP[CHIRON] = (path, NAIF)
        t = _ts().tt_jd(2458500.0)
        with pytest.raises(SPKEvaluationError, match="not available"):
            spk.get_spk_type21_target(CHIRON, t.tt)
        with pytest.raises(SPKEvaluationError, match="not available"):
            spk.calc_spk_body_position(t, CHIRON, FLG_SPEED)
        with pytest.raises(SPKEvaluationError):
            eph.calc(t.tt, CHIRON, FLG_SPEED)


# ---------------------------------------------------------------------------
# Usable coverage: asymmetric edges, TDB, no exclusion inside the span
# ---------------------------------------------------------------------------


class TestUsableCoverage:
    def test_reported_span_is_asymmetric(self, synthetic):
        path, _ = synthetic(10)
        start, end = spk.get_spk_coverage(path)
        stencil = spk._type21_speed_stencil_days()
        distance_au = float(np.linalg.norm(relative_state(STORED_START)[0]))
        light_time = (distance_au + 1.1) / C_AU_DAY
        assert start == pytest.approx(STORED_START + light_time + stencil, abs=1e-9)
        assert end == pytest.approx(STORED_END - stencil, abs=1e-9)
        # Only the stencil (300 s) is lost at the end; the light-time band
        # (~0.145 d for ~24 AU) is subtracted at the start alone.
        assert STORED_END - end == pytest.approx(300.0 / 86400.0, abs=1e-9)
        assert 0.1 < start - STORED_START < 0.2

    def test_direct_path_edges(self, synthetic):
        path, kernel = synthetic(10)
        start, end = spk.get_spk_coverage(path)
        ts = _ts()
        for jd_tdb in (start + 1e-7, end - 1e-6):
            lon, _, _, dlon, _, _ = spk.calc_spk_body_position(
                ts.tdb_jd(jd_tdb), CHIRON, FLG_SPEED
            )
            assert 0.0 <= lon < 360.0 and dlon != 0.0
        # The reader was only ever asked for epochs inside its own bounds.
        assert all(STORED_START <= jd < STORED_END for jd in kernel.evaluated)
        for jd_tdb in (STORED_START, start - 1e-6, end, STORED_END, STORED_END + 1.0):
            with pytest.raises(EphemerisRangeError) as info:
                spk.calc_spk_body_position(ts.tdb_jd(jd_tdb), CHIRON, FLG_SPEED)
            assert "usable coverage" in str(info.value)
            assert info.value.start_jd == pytest.approx(start, abs=1e-9)
            assert info.value.end_jd == pytest.approx(end, abs=1e-9)
            assert info.value.requested_jd == pytest.approx(jd_tdb, abs=1e-9)

    @pytest.mark.parametrize("flags", [FLG_SPEED, FLG_SPEED | FLG_TRUEPOS])
    def test_planet_pipeline_serves_the_whole_reported_span(self, synthetic, flags):
        """Both stencils (1 s apparent, 300 s geometric) fit inside the span."""
        path, _ = synthetic(10)
        start, end = spk.get_spk_coverage(path)
        ts = _ts()
        for jd_tdb in (start + 1e-7, 0.5 * (start + end), end - 1e-6):
            pos, source = _calc_tt(ts.tdb_jd(jd_tdb).tt, flags)
            assert source == "SPK", jd_tdb
            assert math.isfinite(pos[3]) and pos[3] != 0.0

    def test_outside_the_reported_span_is_a_range_miss(self, synthetic, monkeypatch):
        path, _ = synthetic(10)
        start, end = spk.get_spk_coverage(path)
        ts = _ts()
        eph.set_strict_precision(False)
        # The out-of-coverage chain lands on ASSIST when its data files are
        # installed and on the Keplerian elements otherwise; pin the latter so
        # the assertions below do not depend on the machine.
        from libephemeris import rebound_integration

        monkeypatch.setattr(
            rebound_integration, "check_assist_data_available", lambda: False
        )
        one_second = 1.0 / 86400.0
        # Light-time band, exclusive end, beyond the end: never served.
        for jd_tdb, flags in (
            (STORED_START, 0),
            (STORED_START, FLG_SPEED),
            (STORED_END, 0),
            (STORED_END + 5.0, FLG_SPEED),
        ):
            pos, source = _calc_tt(ts.tdb_jd(jd_tdb).tt, flags)
            assert source == "Keplerian", (jd_tdb, flags)
            assert 0.0 <= pos[0] < 360.0
        # The gate is stencil-aware: half a second before the stored end a
        # position is served, but a speed request (samples at +/- 1 s) is a
        # range miss as a whole -- never an SPK position with a Keplerian
        # stencil sample.
        near_end = STORED_END - 0.5 * one_second
        assert _calc_tt(ts.tdb_jd(near_end).tt, 0)[1] == "SPK"
        assert _calc_tt(ts.tdb_jd(near_end).tt, FLG_SPEED)[1] == "Keplerian"
        # Likewise the 300 s geometric stencil is refused where the 1 s
        # apparent one still fits.
        jd_tdb = end + 1e-6
        assert _calc_tt(ts.tdb_jd(jd_tdb).tt, FLG_SPEED)[1] == "SPK"
        assert _calc_tt(ts.tdb_jd(jd_tdb).tt, FLG_SPEED | FLG_TRUEPOS)[1] == "Keplerian"
        # Strict mode reports the range miss instead of degrading.
        eph.set_strict_precision(True)
        with pytest.raises(SPKRequiredError, match="does not cover"):
            _calc_tt(ts.tdb_jd(STORED_START).tt, FLG_SPEED)

    def test_gate_uses_tdb(self, synthetic):
        """A TT epoch whose TDB lies just below a boundary is refused, even
        though the TT value itself is above it."""
        path, _ = synthetic(10)
        start, _ = spk.get_spk_coverage(path)
        ts = _ts()
        # Direct path: gated on the reported (stencil-trimmed) start.
        t = ts.tdb_jd(start - 5e-9)  # 0.43 ms below the boundary, in TDB
        # Precondition: at this date TDB runs behind TT by ~1.4 ms, so the TT
        # value sits above the (TDB) boundary.
        assert t.tt > start
        with pytest.raises(EphemerisRangeError) as info:
            spk.calc_spk_body_position(t, CHIRON, 0)
        assert info.value.requested_jd == pytest.approx(float(t.tdb), abs=1e-9)
        # Planet-pipeline gate for a position request: the stencil-free start.
        servable_start = start - spk._type21_speed_stencil_days()
        t = ts.tdb_jd(servable_start - 5e-9)
        assert t.tt > servable_start
        assert spk.get_spk_type21_target(CHIRON, t.tt, 0) is None
        # And 1 ms inside (TDB) is served, whatever the TT value.
        t_in = ts.tdb_jd(servable_start + 1e-8)
        assert spk.get_spk_type21_target(CHIRON, t_in.tt, 0) is not None

    def test_in_coverage_fault_raises_and_never_degrades(self, synthetic):
        mid = 0.5 * (STORED_START + STORED_END)
        synthetic(10, fault_jd=mid)
        t = _ts().tdb_jd(mid)
        with pytest.raises(SPKEvaluationError) as info:
            spk.calc_spk_body_position(t, CHIRON, FLG_SPEED)
        err = info.value
        assert "Invalid data" in str(err) and "not degraded" in str(err)
        assert err.body_id == CHIRON and err.naif_id == NAIF
        assert err.spk_file == state._SPK_BODY_MAP[CHIRON][0]
        assert err.requested_jd == pytest.approx(mid, abs=0.02)
        # Through calc(), in both strictness modes: no Keplerian answer.
        for strict in (True, False):
            eph.set_strict_precision(strict)
            with pytest.raises(SPKEvaluationError):
                eph.calc(t.tt, CHIRON, FLG_SPEED)
        # A neighbouring epoch the kernel serves is unaffected.
        pos, source = _calc_tt(_ts().tdb_jd(mid + 1.0).tt, FLG_SPEED)
        assert source == "SPK"

    def test_kernel_that_cannot_be_opened_is_loud(self, synthetic, monkeypatch):
        synthetic(10)
        monkeypatch.setattr(spk, "_load_type21_kernel", lambda _p: None)
        t = _ts().tt_jd(2458500.0)
        with pytest.raises(SPKEvaluationError, match="could not be opened"):
            spk.calc_spk_body_position(t, CHIRON, FLG_SPEED)
        with pytest.raises(SPKEvaluationError):
            eph.calc(t.tt, CHIRON, FLG_SPEED)


# ---------------------------------------------------------------------------
# The four callers: range miss falls back, evaluation failure propagates
# ---------------------------------------------------------------------------


class TestCallers:
    MID = 0.5 * (STORED_START + STORED_END)

    def test_planet_pipeline_direct_fallback_path(self, synthetic, monkeypatch):
        """planets._calc_body -> calc_spk_body_position (type-21 target None)."""
        synthetic(10, fault_jd=self.MID)
        monkeypatch.setattr(spk, "get_spk_type21_target", lambda *_a, **_k: None)
        with pytest.raises(SPKEvaluationError):
            eph.calc(self.MID, CHIRON, FLG_SPEED)

    def test_auto_download_path(self, synthetic, monkeypatch):
        """planets._try_auto_spk_download -> calc_spk_body_position."""
        installed = {}

        def fake_download(**kwargs):
            installed["path"], _ = synthetic(10, fault_jd=self.MID)
            return installed["path"]

        monkeypatch.setattr(spk, "download_and_register_spk", fake_download)
        t = _ts().tt_jd(self.MID)
        with pytest.raises(SPKEvaluationError):
            planets._try_auto_spk_download(t, CHIRON, FLG_SPEED)
        assert installed["path"] == state._SPK_BODY_MAP[CHIRON][0]

    def test_minor_body_heliocentric(self, synthetic):
        """minor_bodies.calc_minor_body_heliocentric."""
        synthetic(10, fault_jd=self.MID)
        with pytest.raises(SPKEvaluationError):
            minor_bodies.calc_minor_body_heliocentric(CHIRON, self.MID)
        # Outside the usable coverage: the documented Keplerian answer.
        outside = STORED_END + 5.0
        assert minor_bodies.calc_minor_body_heliocentric(
            CHIRON, outside
        ) == minor_bodies.calc_minor_body_heliocentric(CHIRON, outside, use_spk=False)

    def test_asteroid_by_number(self, synthetic, monkeypatch):
        """minor_bodies.calc_asteroid_by_number."""
        number = 70000
        synthetic(10, ipl=AST_OFFSET + number, fault_jd=self.MID)
        monkeypatch.setattr(
            minor_bodies, "fetch_orbital_elements_from_sbdb", lambda *a, **k: None
        )
        with pytest.raises(SPKEvaluationError):
            minor_bodies.calc_asteroid_by_number(number, self.MID)


# ---------------------------------------------------------------------------
# Exception contract and cached-kernel discovery
# ---------------------------------------------------------------------------


class TestContract:
    def test_exception_hierarchy(self):
        assert issubclass(SPKEvaluationError, CalculationError)
        assert not issubclass(SPKEvaluationError, EphemerisRangeError)
        assert not issubclass(SPKEvaluationError, ValueError)
        assert eph.SPKEvaluationError is SPKEvaluationError
        assert "SPKEvaluationError" in eph.__all__
        err = SPKEvaluationError(
            "m", body_id=1, naif_id=2, spk_file="f", requested_jd=3.0
        )
        assert str(err) == "m" and "naif_id=2" in repr(err)

    def test_discovery_tolerates_the_light_time_band(self, synthetic, monkeypatch):
        """A kernel downloaded for [start, end] is still found for that range
        although get_spk_coverage() reports the narrower usable span."""
        path, _ = synthetic(10, name="2060_synthetic.bsp", register=False)
        cache_dir = os.path.dirname(path)
        monkeypatch.setattr(spk_auto, "_is_valid_bsp", lambda _p: True)
        assert spk_auto.is_spk_cached("2060", STORED_START, STORED_END, cache_dir)
        assert (
            spk_auto._find_covering_spk("2060", STORED_START, STORED_END, cache_dir)
            == path
        )
        assert not spk_auto.is_spk_cached(
            "2060", STORED_START - 3.0, STORED_END, cache_dir
        )
        assert not spk_auto.is_spk_cached(
            "2060", STORED_START, STORED_END + 1.0, cache_dir
        )

    def test_probe_does_not_cache_an_unregistered_kernel(self, tmp_path, monkeypatch):
        """get_spk_coverage() on a scanned file opens and closes the reader."""
        path = os.path.join(str(tmp_path), "scan.bsp")
        with open(path, "wb") as fh:
            fh.write(b"\x00" * 16)
        kernel = _SyntheticKernel(10)
        monkeypatch.setattr(spk, "_detect_spk_type", lambda _p: 21)
        monkeypatch.setattr(
            spk,
            "_get_spktype21",
            lambda: type("R", (), {"open": staticmethod(lambda _p: kernel)}),
        )
        assert spk.get_spk_coverage(path) is not None
        assert kernel.closed
        assert path not in state._SPK_TYPE21_KERNELS
