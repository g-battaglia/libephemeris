# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for :mod:`libephemeris.planetary_moons`.

These tests drive the satellite-SPK registration, position, coverage and
cleanup code paths using lightweight fake kernel/target/observer objects so
that no real ``.bsp`` files or network access are required. They target the
specific lines/branches that the full suite leaves uncovered.
"""

from __future__ import annotations

import math
import os

import pytest

from libephemeris import planetary_moons as pm
from libephemeris import planets as planets_mod
from libephemeris import state
from libephemeris.constants import (
    FLG_BARYCTR,
    FLG_HELCTR,
    FLG_NOABERR,
    FLG_SIDEREAL,
    FLG_SPEED,
    FLG_TRUEPOS,
    MOON_IO,
    MOON_TITAN,
    NAIF_IO,
    PLMOON_OFFSET,
)
from libephemeris.exceptions import EphemerisRangeError


# ---------------------------------------------------------------------------
# Test doubles
# ---------------------------------------------------------------------------

# Ecliptic obliquity used to build vectors at a chosen ecliptic longitude.
_EPS = math.radians(23.4392911)


def _ecl_vector(longitude_deg: float, dist_au: float = 5.0) -> tuple:
    """Return an ICRF (x, y, z) AU vector at a given ecliptic longitude."""
    lam = math.radians(longitude_deg)
    xe, ye, ze = math.cos(lam), math.sin(lam), 0.0
    x = xe
    y = ye * math.cos(_EPS) - ze * math.sin(_EPS)
    z = ye * math.sin(_EPS) + ze * math.cos(_EPS)
    return (x * dist_au, y * dist_au, z * dist_au)


class _Vec:
    """Minimal stand-in for Skyfield's ``.position``/``.velocity`` wrapper."""

    def __init__(self, xyz):
        self.au = xyz
        self.au_per_d = xyz


class _State:
    """Minimal stand-in for a Skyfield astrometric state."""

    def __init__(self, position, velocity):
        self.position = _Vec(position)
        self.velocity = _Vec(velocity)


class _Target:
    """Fake SPK target returning a fixed (or per-time) position."""

    def __init__(self, position=(5.0, 0.0, 0.0), per_time=None, raise_value=False):
        self._position = position
        self._per_time = per_time
        self._raise_value = raise_value

    def at(self, t):
        if self._raise_value:
            raise ValueError("outside coverage")
        if self._per_time is not None:
            return _State(self._per_time(t), (0.0, 0.0, 0.0))
        return _State(self._position, (0.0, 0.0, 0.0))


class _Observer:
    """Fake observer with a fixed position and small velocity."""

    def __init__(self, position=(0.0, 0.0, 0.0), velocity=(0.001, 0.001, 0.001)):
        self._position = position
        self._velocity = velocity

    def at(self, t):
        return _State(self._position, self._velocity)


class _Segment:
    """Fake SPK segment for coverage / auto-detect scans."""

    def __init__(self, target, start_jd=2451000.0, end_jd=2452000.0):
        self.target = target
        self.start_jd = start_jd
        self.end_jd = end_jd


class _Spk:
    def __init__(self, segments):
        self.segments = segments


class _Kernel:
    """Fake Skyfield kernel: subscriptable by NAIF id, optional ``spk``."""

    def __init__(self, targets, segments=None, spk_raises=False, closed_raises=None):
        self._targets = dict(targets)
        self._segments = segments
        self._spk_raises = spk_raises
        self._closed_raises = closed_raises
        self.closed = False

    def __getitem__(self, naif_id):
        try:
            return self._targets[naif_id]
        except KeyError:
            raise KeyError(naif_id)

    @property
    def spk(self):
        if self._spk_raises:
            raise AttributeError("no spk")
        return _Spk(self._segments if self._segments is not None else [])

    def close(self):
        self.closed = True
        if self._closed_raises is not None:
            raise self._closed_raises


@pytest.fixture(autouse=True)
def _clean_moon_state():
    """Ensure the module's global registries are empty around each test."""
    pm.close_moon_kernels()
    yield
    pm.close_moon_kernels()


def _register_fake(monkeypatch, path, kernel):
    """Directly seed the module registries with a fake kernel by NAIF target."""
    pm._MOON_SPK_KERNELS[path] = kernel
    pm._MOON_SPK_BY_BODY[PLMOON_OFFSET + NAIF_IO] = path


# ---------------------------------------------------------------------------
# register_moon_spk: path resolution branches (248, 249->252, 253)
# ---------------------------------------------------------------------------


def test_register_resolves_file_in_data_dir(monkeypatch, tmp_path):
    """Relative file found in the data dir is promoted to its full path (248)."""
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    spk_path = data_dir / "jup.bsp"
    spk_path.write_bytes(b"FAKE")

    monkeypatch.setattr(state, "_get_data_dir", lambda: str(data_dir))

    kernel = _Kernel({NAIF_IO: _Target()})
    monkeypatch.setattr(pm, "get_loader", lambda: (lambda p: kernel))

    pm.register_moon_spk("jup.bsp", moons=[MOON_IO])

    # Registered under the resolved absolute path.
    assert pm._MOON_SPK_BY_BODY[MOON_IO] == str(spk_path)
    assert str(spk_path) in pm._MOON_SPK_KERNELS


def test_register_relative_existing_in_cwd(monkeypatch, tmp_path):
    """Relative file absent from data dir but present in CWD (arc 249->252)."""
    data_dir = tmp_path / "empty_data"
    data_dir.mkdir()
    monkeypatch.setattr(state, "_get_data_dir", lambda: str(data_dir))

    work = tmp_path / "work"
    work.mkdir()
    rel_name = "local.bsp"
    (work / rel_name).write_bytes(b"FAKE")
    monkeypatch.chdir(work)

    kernel = _Kernel({NAIF_IO: _Target()})
    monkeypatch.setattr(pm, "get_loader", lambda: (lambda p: kernel))

    pm.register_moon_spk(rel_name, moons=[MOON_IO])
    assert pm._MOON_SPK_BY_BODY[MOON_IO] == rel_name


def test_register_absolute_missing_raises(tmp_path):
    """Absolute path that does not exist raises at the final guard (253)."""
    missing = str(tmp_path / "nope.bsp")
    assert os.path.isabs(missing)
    with pytest.raises(FileNotFoundError):
        pm.register_moon_spk(missing)


# ---------------------------------------------------------------------------
# register_moon_spk: explicit moons list (263-273) and ValueError (297)
# ---------------------------------------------------------------------------


def test_register_with_moons_list(monkeypatch, tmp_path):
    """moons=[...] path handles None ids, missing targets, and success."""
    spk = tmp_path / "mix.bsp"
    spk.write_bytes(b"FAKE")

    # Kernel only has Io; Titan id will KeyError.
    kernel = _Kernel({NAIF_IO: _Target()})
    monkeypatch.setattr(pm, "get_loader", lambda: (lambda p: kernel))

    # 0 -> not a moon id (naif_id None -> continue, 265-266)
    # MOON_TITAN -> naif resolves but not in kernel (KeyError -> continue, 270-271)
    # MOON_IO -> registered (272-273)
    pm.register_moon_spk(str(spk), moons=[0, MOON_TITAN, MOON_IO])

    assert pm._MOON_SPK_BY_BODY[MOON_IO] == str(spk)
    assert pm._MOON_SPK_BY_BODY[PLMOON_OFFSET + NAIF_IO] == str(spk)
    assert MOON_TITAN not in pm._MOON_SPK_BY_BODY


def test_register_with_moons_list_none_found_raises(monkeypatch, tmp_path):
    """Explicit moons all missing from kernel -> ValueError (297)."""
    spk = tmp_path / "empty.bsp"
    spk.write_bytes(b"FAKE")

    kernel = _Kernel({})  # no targets at all
    monkeypatch.setattr(pm, "get_loader", lambda: (lambda p: kernel))

    with pytest.raises(ValueError):
        pm.register_moon_spk(str(spk), moons=[MOON_IO])


# ---------------------------------------------------------------------------
# register_moon_spk: auto-detect branches (284-285, 290-291)
# ---------------------------------------------------------------------------


def test_register_autodetect_spk_attribute_error(monkeypatch, tmp_path):
    """Auto-detect: kernel.spk access raises -> fall back to MOON_NAIF_MAP.

    Only Io chains successfully; every other NAIF id raises KeyError in the
    auto-detect loop, exercising the continue branch (290-291).
    """
    spk = tmp_path / "auto.bsp"
    spk.write_bytes(b"FAKE")

    kernel = _Kernel({NAIF_IO: _Target()}, spk_raises=True)
    monkeypatch.setattr(pm, "get_loader", lambda: (lambda p: kernel))

    pm.register_moon_spk(str(spk))  # moons=None -> auto-detect

    assert pm._MOON_SPK_BY_BODY[PLMOON_OFFSET + NAIF_IO] == str(spk)
    assert pm._MOON_SPK_BY_BODY[MOON_IO] == str(spk)


def test_register_autodetect_via_segments(monkeypatch, tmp_path):
    """Auto-detect over real segments registers the chained target."""
    spk = tmp_path / "seg.bsp"
    spk.write_bytes(b"FAKE")

    segments = [_Segment(NAIF_IO), _Segment(999)]  # 999 not in (400, 999)
    kernel = _Kernel({NAIF_IO: _Target()}, segments=segments)
    monkeypatch.setattr(pm, "get_loader", lambda: (lambda p: kernel))

    pm.register_moon_spk(str(spk))
    assert pm._MOON_SPK_BY_BODY[PLMOON_OFFSET + NAIF_IO] == str(spk)


# ---------------------------------------------------------------------------
# unregister_moon_spk (311-325)
# ---------------------------------------------------------------------------


def test_unregister_removes_bodies_and_kernel(monkeypatch, tmp_path):
    """unregister clears body keys (basename match) and closes the kernel."""
    spk = tmp_path / "u.bsp"
    spk.write_bytes(b"FAKE")
    path = str(spk)

    kernel = _Kernel({NAIF_IO: _Target()})
    monkeypatch.setattr(pm, "get_loader", lambda: (lambda p: kernel))
    pm.register_moon_spk(path, moons=[MOON_IO])

    assert path in pm._MOON_SPK_KERNELS
    pm.unregister_moon_spk(path)

    assert path not in pm._MOON_SPK_KERNELS
    assert MOON_IO not in pm._MOON_SPK_BY_BODY
    assert kernel.closed is True


def test_unregister_close_swallows_error(monkeypatch):
    """unregister tolerates a kernel whose close() raises (323-324)."""
    path = "/tmp/raisy.bsp"
    kernel = _Kernel({NAIF_IO: _Target()}, closed_raises=OSError("boom"))
    pm._MOON_SPK_KERNELS[path] = kernel
    pm._MOON_SPK_BY_BODY[MOON_IO] = path

    pm.unregister_moon_spk(path)  # should not raise
    assert path not in pm._MOON_SPK_KERNELS
    assert MOON_IO not in pm._MOON_SPK_BY_BODY


def test_unregister_without_cached_kernel(monkeypatch):
    """unregister of a body whose kernel is absent skips close (arc 320->exit)."""
    path = "phantom.bsp"
    pm._MOON_SPK_BY_BODY[MOON_IO] = path  # no kernel object cached

    pm.unregister_moon_spk(path)  # should not raise
    assert MOON_IO not in pm._MOON_SPK_BY_BODY
    assert path not in pm._MOON_SPK_KERNELS


# ---------------------------------------------------------------------------
# calc_moon_position: early None (450) and target resolution (462-470)
# ---------------------------------------------------------------------------


def test_calc_position_non_moon_id_returns_none():
    """A non-moon id yields None immediately (450)."""
    ts = state.get_timescale()
    t = ts.tt_jd(2451545.0)
    assert pm.calc_moon_position(t, 0, FLG_SPEED) is None


def test_calc_position_primary_keyerror_then_fallback(monkeypatch):
    """Primary kernel KeyErrors (462-463); fallback scan finds target (466-468)."""
    ts = state.get_timescale()
    t = ts.tt_jd(2451545.0)

    primary = _Kernel({})  # registered for Io but lacks the target -> KeyError
    fallback = _Kernel({NAIF_IO: _Target(position=_ecl_vector(0.0))})
    pm._MOON_SPK_KERNELS["primary.bsp"] = primary
    pm._MOON_SPK_KERNELS["fallback.bsp"] = fallback
    pm._MOON_SPK_BY_BODY[MOON_IO] = "primary.bsp"

    monkeypatch.setattr(pm, "get_planets", lambda: {"earth": _Observer()})
    monkeypatch.setattr(planets_mod, "get_ayanamsa_ut", lambda jd: 0.0)

    result = pm.calc_moon_position(t, MOON_IO, FLG_TRUEPOS | FLG_NOABERR)
    assert result is not None
    assert result[2] > 0  # positive distance


def test_calc_position_no_target_anywhere_returns_none(monkeypatch):
    """No kernel contains the target -> None after fallback scan (469-472)."""
    ts = state.get_timescale()
    t = ts.tt_jd(2451545.0)

    pm._MOON_SPK_KERNELS["a.bsp"] = _Kernel({})
    pm._MOON_SPK_KERNELS["b.bsp"] = _Kernel({})
    pm._MOON_SPK_BY_BODY[MOON_IO] = "a.bsp"

    assert pm.calc_moon_position(t, MOON_IO, 0) is None


# ---------------------------------------------------------------------------
# calc_moon_position: observer selection + aberration arc (482, 484, 508-509,
# 527->535)
# ---------------------------------------------------------------------------


def _setup_simple_io(monkeypatch, position=None):
    """Register a single kernel with Io and patch ephemeris helpers."""
    if position is None:
        position = _ecl_vector(123.0)
    kernel = _Kernel({NAIF_IO: _Target(position=position)})
    pm._MOON_SPK_KERNELS["k.bsp"] = kernel
    pm._MOON_SPK_BY_BODY[MOON_IO] = "k.bsp"
    monkeypatch.setattr(
        pm,
        "get_planets",
        lambda: {"earth": _Observer(), "sun": _Observer(position=(1.0, 0.0, 0.0))},
    )
    monkeypatch.setattr(planets_mod, "get_ayanamsa_ut", lambda jd: 30.0)


def test_calc_position_barycentric(monkeypatch):
    """Barycentric: observer None (482), obs_xyz zero (508-509), skip aberr."""
    ts = state.get_timescale()
    t = ts.tt_jd(2451545.0)
    _setup_simple_io(monkeypatch)

    result = pm.calc_moon_position(t, MOON_IO, FLG_BARYCTR)
    assert result is not None
    # No light-time / aberration applied; distance == |position| (~5 AU).
    assert result[2] == pytest.approx(5.0, abs=1e-6)


def test_calc_position_heliocentric(monkeypatch):
    """Heliocentric selects the Sun as observer (484)."""
    ts = state.get_timescale()
    t = ts.tt_jd(2451545.0)
    _setup_simple_io(monkeypatch)

    result = pm.calc_moon_position(t, MOON_IO, FLG_HELCTR)
    assert result is not None
    assert result[2] > 0


# ---------------------------------------------------------------------------
# calc_moon_position: speed wrap branches (566, 568) and sidereal (572-573)
# ---------------------------------------------------------------------------


def _wrap_target(prev_lon, cur_lon, next_lon):
    """Build a target that returns crafted longitudes at t-dt, t, t+dt."""
    dt = 1.0 / 86400.0
    base = 2451545.0

    def per_time(t):
        # Identify which sample this is by comparing tt to base +/- dt.
        delta = round((t.tt - base) / dt)
        if delta <= -1:
            return _ecl_vector(prev_lon)
        if delta >= 1:
            return _ecl_vector(next_lon)
        return _ecl_vector(cur_lon)

    return _Target(per_time=per_time)


def test_calc_position_speed_wrap_positive(monkeypatch):
    """Central-difference longitude wraps positive -> subtract (566)."""
    ts = state.get_timescale()
    t = ts.tt_jd(2451545.0)

    kernel = _Kernel({NAIF_IO: _wrap_target(1.0, 0.0, 359.0)})
    pm._MOON_SPK_KERNELS["k.bsp"] = kernel
    pm._MOON_SPK_BY_BODY[MOON_IO] = "k.bsp"
    monkeypatch.setattr(pm, "get_planets", lambda: {"earth": _Observer()})
    monkeypatch.setattr(planets_mod, "get_ayanamsa_ut", lambda jd: 0.0)

    result = pm.calc_moon_position(t, MOON_IO, FLG_SPEED | FLG_TRUEPOS | FLG_NOABERR)
    assert result is not None
    # Wrapped rate must be small in magnitude, not ~15M deg/day.
    assert abs(result[3]) < 1.0e6


def test_calc_position_speed_wrap_negative(monkeypatch):
    """Central-difference longitude wraps negative -> add (568)."""
    ts = state.get_timescale()
    t = ts.tt_jd(2451545.0)

    kernel = _Kernel({NAIF_IO: _wrap_target(359.0, 0.0, 1.0)})
    pm._MOON_SPK_KERNELS["k.bsp"] = kernel
    pm._MOON_SPK_BY_BODY[MOON_IO] = "k.bsp"
    monkeypatch.setattr(pm, "get_planets", lambda: {"earth": _Observer()})
    monkeypatch.setattr(planets_mod, "get_ayanamsa_ut", lambda jd: 0.0)

    result = pm.calc_moon_position(t, MOON_IO, FLG_SPEED | FLG_TRUEPOS | FLG_NOABERR)
    assert result is not None
    assert abs(result[3]) < 1.0e6


def test_calc_position_geocentric_with_aberration(monkeypatch):
    """Geocentric apparent path applies light-time + aberration (527-534)."""
    ts = state.get_timescale()
    t = ts.tt_jd(2451545.0)
    _setup_simple_io(monkeypatch, position=_ecl_vector(45.0))

    result = pm.calc_moon_position(t, MOON_IO, FLG_SPEED)
    assert result is not None
    assert result[2] > 0


def test_calc_position_sidereal(monkeypatch):
    """Sidereal flag subtracts the ayanamsa (572-573)."""
    ts = state.get_timescale()
    t = ts.tt_jd(2451545.0)
    _setup_simple_io(monkeypatch, position=_ecl_vector(123.0))
    # The sidereal path now routes through the canonical flag-aware helper
    # (_apply_sidereal_correction -> _get_ayanamsa_for_flags), matching the
    # planet path, instead of subtracting the mean get_ayanamsa_ut directly.
    monkeypatch.setattr(
        planets_mod,
        "_get_ayanamsa_for_flags",
        lambda jd, iflag, sid_mode=None: 30.0,
    )

    tropical = pm.calc_moon_position(t, MOON_IO, FLG_TRUEPOS | FLG_NOABERR)
    sidereal = pm.calc_moon_position(
        t, MOON_IO, FLG_TRUEPOS | FLG_NOABERR | FLG_SIDEREAL
    )
    assert tropical is not None and sidereal is not None
    expected = (tropical[0] - 30.0) % 360.0
    assert sidereal[0] == pytest.approx(expected, abs=1e-9)


def test_calc_position_outside_coverage_raises(monkeypatch):
    """A target whose .at() raises ValueError becomes EphemerisRangeError."""
    ts = state.get_timescale()
    t = ts.tt_jd(2451545.0)

    kernel = _Kernel({NAIF_IO: _Target(raise_value=True)})
    pm._MOON_SPK_KERNELS["k.bsp"] = kernel
    pm._MOON_SPK_BY_BODY[MOON_IO] = "k.bsp"
    monkeypatch.setattr(pm, "get_planets", lambda: {"earth": _Observer()})
    monkeypatch.setattr(planets_mod, "get_ayanamsa_ut", lambda jd: 0.0)

    with pytest.raises(EphemerisRangeError):
        pm.calc_moon_position(t, MOON_IO, 0)


# ---------------------------------------------------------------------------
# get_moon_coverage (596-613)
# ---------------------------------------------------------------------------


def test_get_moon_coverage_non_moon_returns_none():
    """Non-moon id -> None (596-598)."""
    assert pm.get_moon_coverage(0) is None


def test_get_moon_coverage_resolves_segments(monkeypatch):
    """Coverage resolves via canonical id and reads segment bounds (599-614)."""
    segments = [_Segment(NAIF_IO, 2450000.0, 2455000.0)]
    kernel = _Kernel({NAIF_IO: _Target()}, segments=segments)
    path = "cov.bsp"
    pm._MOON_SPK_KERNELS[path] = kernel
    # Key only under canonical id to exercise the second .get() lookup.
    pm._MOON_SPK_BY_BODY[PLMOON_OFFSET + NAIF_IO] = path

    coverage = pm.get_moon_coverage(MOON_IO)
    assert coverage == (2450000.0, 2455000.0)


def test_get_moon_coverage_no_kernel_returns_none():
    """Registered body but missing kernel object -> None (602-604)."""
    pm._MOON_SPK_BY_BODY[MOON_IO] = "ghost.bsp"  # no kernel cached
    assert pm.get_moon_coverage(MOON_IO) is None


# ---------------------------------------------------------------------------
# close_moon_kernels (633-634)
# ---------------------------------------------------------------------------


def test_close_moon_kernels_resets_globals():
    """close_moon_kernels closes kernels and clears both registries (633-642)."""
    kernel = _Kernel({NAIF_IO: _Target()})
    pm._MOON_SPK_KERNELS["x.bsp"] = kernel
    pm._MOON_SPK_BY_BODY[MOON_IO] = "x.bsp"

    pm.close_moon_kernels()

    assert kernel.closed is True
    assert pm._MOON_SPK_KERNELS == {}
    assert pm._MOON_SPK_BY_BODY == {}
