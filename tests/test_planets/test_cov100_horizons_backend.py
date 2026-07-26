# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Full-mock coverage tests for :mod:`libephemeris.horizons_backend`.

The Horizons backend does network IO; per project policy these tests use
FULL MOCKS to exercise every reachable line without touching the network.

Two seams are mocked:

1. The calc pipeline (:func:`horizons_calc_ut`, :func:`_to_ecliptic_output`,
   the deflector adapter) is driven with a :class:`FakeHorizonsClient` that
   returns canned :class:`StateVector` objects from ``fetch_state_vector`` and
   a dict from ``fetch_batch``.
2. The HTTP client (:class:`HorizonsClient`, ``_fetch_with_retry``,
   ``_parse_response``) is driven by patching the policy-gated ``net.open_url``
   choke point (imported into ``horizons_backend``) with a fake context manager
   whose ``.read()`` returns canned JSON bytes, and by patching ``time.sleep``
   so retry backoff is instantaneous.
"""

from __future__ import annotations

import json
import math
from unittest import mock

import pytest

pytestmark = pytest.mark.usefixtures("allow_mocked_network")

from libephemeris import horizons_backend as hb
from libephemeris.horizons_backend import (
    HorizonsClient,
    StateVector,
    _HorizonsDeflectorSource,
    _apply_deflection_horizons,
    _calc_analytical,
    _get_ssl_context,
    _to_ecliptic_output,
    horizons_calc_ut,
)
from libephemeris.exceptions import UnknownBodyError

JD = 2451545.0  # J2000.0

# Flag values (mirrored from libephemeris.constants for explicitness).
from libephemeris.constants import (  # noqa: E402
    FLG_BARYCTR,
    FLG_EQUATORIAL,
    FLG_HELCTR,
    FLG_ICRS,
    FLG_J2000,
    FLG_NOABERR,
    FLG_NOGDEFL,
    FLG_NONUT,
    FLG_RADIANS,
    FLG_SIDEREAL,
    FLG_SPEED,
    FLG_SWIEPH,
    FLG_TOPOCTR,
    FLG_TRUEPOS,
    FLG_XYZ,
    HARRINGTON,
)


# ---------------------------------------------------------------------------
# Fakes
# ---------------------------------------------------------------------------


def _sv(x, y, z, vx=0.0, vy=0.0, vz=0.0):
    """Build a StateVector quickly."""
    return StateVector(x, y, z, vx, vy, vz)


class FakeHorizonsClient:
    """Stand-in for HorizonsClient that never touches the network.

    Returns position-only canned vectors keyed by Horizons command. Positions
    are chosen so the deflection/aberration math stays well-conditioned (no
    division by zero, target far from the deflector line of sight).
    """

    # command -> barycentric position (AU); velocities small but nonzero.
    _POS = {
        "10": (0.001, 0.002, 0.0005, 1e-6, 2e-6, 3e-6),  # Sun
        "399": (0.5, 0.8, 0.3, 0.01, -0.02, 0.005),  # Earth
        "499": (1.2, -1.0, 0.4, 0.005, 0.006, 0.001),  # Mars (target)
        "5": (4.0, 2.0, 0.1, 0.001, 0.002, 0.0),  # Jupiter bary
        "6": (-8.0, 3.0, 0.5, 0.0005, 0.001, 0.0),  # Saturn bary
        "599": (4.0, 2.0, 0.1, 0.001, 0.002, 0.0),  # Jupiter planet (unused)
    }

    def __init__(self):
        self.calls = []

    def fetch_state_vector(self, command, jd, center="@0", time_type="TDB"):
        self.calls.append((command, jd, center, time_type))
        p = self._POS.get(command)
        if p is None:
            raise KeyError(command)
        return _sv(*p)

    def fetch_batch(self, requests, errors=None):
        out = {}
        for cmd, jd, center in requests:
            p = self._POS.get(cmd)
            if p is not None:
                out[(cmd, jd, center)] = _sv(*p)
            elif errors is not None:
                errors[(cmd, jd, center)] = KeyError(cmd)
        return out


class _FakeResp:
    """Context-manager double for the object returned by urlopen()."""

    def __init__(self, payload):
        self._payload = payload

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def read(self):
        return self._payload


def _good_json_bytes():
    """A well-formed Horizons VECTORS CSV response."""
    result = (
        "headers...\n"
        "$$SOE\n"
        "2451545.0, A.D. 2000-Jan-01, 1.2, -1.0, 0.4, 0.005, 0.006, 0.001, "
        "0.01, 1.6, 0.0,\n"
        "$$EOE\n"
        "footer...\n"
    )
    return json.dumps({"result": result}).encode("utf-8")


# ===========================================================================
# StateVector dataclass (lines 87-92, 96, 100)
# ===========================================================================


class TestStateVector:
    def test_attrs_and_properties(self):
        sv = StateVector(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
        assert sv.x == 1.0 and sv.y == 2.0 and sv.z == 3.0
        assert sv.vx == 4.0 and sv.vy == 5.0 and sv.vz == 6.0
        assert sv.pos == (1.0, 2.0, 3.0)
        assert sv.vel == (4.0, 5.0, 6.0)


# ===========================================================================
# _get_ssl_context (lines 120-124)
# ===========================================================================


class TestSSLContext:
    def test_builds_once_and_caches(self, monkeypatch):
        # Reset module-level cache so the build branch is exercised.
        monkeypatch.setattr(hb, "_SSL_CTX", None)
        sentinel = object()
        fake_ssl = mock.MagicMock()
        fake_ssl.create_default_context.return_value = sentinel
        fake_certifi = mock.MagicMock()
        fake_certifi.where.return_value = "/fake/ca.pem"
        with mock.patch.dict("sys.modules", {"ssl": fake_ssl, "certifi": fake_certifi}):
            ctx1 = _get_ssl_context()
            assert ctx1 is sentinel
            # Second call returns the cached object without rebuilding.
            ctx2 = _get_ssl_context()
            assert ctx2 is sentinel
        fake_ssl.create_default_context.assert_called_once()

    def test_double_checked_lock_already_built(self, monkeypatch):
        # Exercise the inner "already built" branch (122->128): the outer
        # check sees None, but by the time the lock is acquired another
        # thread has populated _SSL_CTX. Simulate with a fake lock whose
        # __enter__ sets the module global.
        monkeypatch.setattr(hb, "_SSL_CTX", None)
        sentinel = object()

        class _SettingLock:
            def __enter__(self):
                hb._SSL_CTX = sentinel
                return self

            def __exit__(self, *exc):
                return False

        with mock.patch.object(hb, "_SSL_CTX_LOCK", _SettingLock()):
            ctx = _get_ssl_context()
        assert ctx is sentinel


def test_ssl_context_unit_tests_leave_a_runtime_compatible_cache() -> None:
    """Mock sentinels must not leak into later live Horizons calculations."""
    assert hb._SSL_CTX is None or callable(getattr(hb._SSL_CTX, "wrap_socket", None))


# ===========================================================================
# HorizonsClient.fetch_state_vector + cache (lines 148-208)
# ===========================================================================


class TestFetchStateVector:
    def test_fetch_parses_and_caches(self):
        client = HorizonsClient(max_cache_size=4)
        with (
            mock.patch.object(
                hb, "open_url", return_value=_FakeResp(_good_json_bytes())
            ) as op,
            mock.patch.object(hb, "_get_ssl_context", return_value=None),
        ):
            sv1 = client.fetch_state_vector("499", JD, "@0", "TDB")
            assert sv1.x == 1.2 and sv1.vz == 0.001
            # Second identical call hits the cache (no new urlopen).
            sv2 = client.fetch_state_vector("499", JD, "@0", "TDB")
            assert sv2 is sv1
        assert op.call_count == 1

    def test_cache_lru_eviction(self):
        client = HorizonsClient(max_cache_size=1)
        with (
            mock.patch.object(
                hb, "open_url", return_value=_FakeResp(_good_json_bytes())
            ),
            mock.patch.object(hb, "_get_ssl_context", return_value=None),
        ):
            client.fetch_state_vector("499", JD, "@0", "TDB")
            client.fetch_state_vector("499", JD + 1.0, "@0", "TDB")
        # max_cache_size=1 forced eviction of the first entry.
        assert len(client._cache) == 1


# ===========================================================================
# HorizonsClient.fetch_batch (lines 222-250)
# ===========================================================================


class TestFetchBatch:
    def test_batch_cache_hit_and_parallel_fetch(self):
        client = HorizonsClient(max_workers=2)
        with (
            mock.patch.object(
                hb, "open_url", return_value=_FakeResp(_good_json_bytes())
            ),
            mock.patch.object(hb, "_get_ssl_context", return_value=None),
        ):
            # Prime the cache with one entry.
            client.fetch_state_vector("499", JD, "@0", "TDB")
            res = client.fetch_batch([("499", JD, "@0"), ("10", JD, "@0")])
        assert ("499", JD, "@0") in res
        assert ("10", JD, "@0") in res

    def test_batch_all_cached_short_circuit(self):
        client = HorizonsClient()
        with (
            mock.patch.object(
                hb, "open_url", return_value=_FakeResp(_good_json_bytes())
            ),
            mock.patch.object(hb, "_get_ssl_context", return_value=None),
        ):
            client.fetch_state_vector("499", JD, "@0", "TDB")
            res = client.fetch_batch([("499", JD, "@0")])
        assert res == {("499", JD, "@0"): res[("499", JD, "@0")]}

    def test_batch_skips_failed_fetch(self):
        client = HorizonsClient(max_workers=2)

        def _boom(*a, **k):
            raise OSError("network down")

        with mock.patch.object(client, "fetch_state_vector", side_effect=_boom):
            res = client.fetch_batch([("499", JD, "@0")])
        assert res == {}


# ===========================================================================
# _fetch_with_retry (lines 256-289) + _parse_response (291-336)
# ===========================================================================


class TestFetchWithRetry:
    def test_api_error_raises_keyerror_no_retry(self):
        client = HorizonsClient()
        payload = json.dumps({"error": "No matches found"}).encode("utf-8")
        with (
            mock.patch.object(hb, "open_url", return_value=_FakeResp(payload)) as op,
            mock.patch.object(hb, "_get_ssl_context", return_value=None),
        ):
            with pytest.raises(KeyError):
                client.fetch_state_vector("nope", JD)
        # KeyError is not retried.
        assert op.call_count == 1

    def test_retry_then_success(self):
        client = HorizonsClient()
        good = _FakeResp(_good_json_bytes())
        seq = [OSError("temp1"), OSError("temp2"), good]

        def _urlopen(*a, **k):
            item = seq.pop(0)
            if isinstance(item, Exception):
                raise item
            return item

        with (
            mock.patch.object(hb, "open_url", side_effect=_urlopen),
            mock.patch.object(hb, "_get_ssl_context", return_value=None),
            mock.patch("time.sleep", return_value=None),
        ):
            sv = client.fetch_state_vector("499", JD)
        assert sv.x == 1.2

    def test_all_attempts_fail_raises_connectionerror(self):
        client = HorizonsClient()
        with (
            mock.patch.object(hb, "open_url", side_effect=OSError("down")),
            mock.patch.object(hb, "_get_ssl_context", return_value=None),
            mock.patch("time.sleep", return_value=None),
        ):
            with pytest.raises(ConnectionError):
                client.fetch_state_vector("499", JD)


class TestParseResponse:
    def _parse(self, result_text):
        client = HorizonsClient()
        return client._parse_response({"result": result_text}, "499")

    def test_missing_markers_raises(self):
        with pytest.raises(ValueError):
            self._parse("no markers here at all")

    def test_empty_block_raises(self):
        with pytest.raises(ValueError):
            self._parse("$$SOE\n   \n$$EOE")

    def test_too_few_values_raises(self):
        # Only 3 numeric values -> < 7 required.
        with pytest.raises(ValueError):
            self._parse("$$SOE\n2451545.0, 1.0, 2.0\n$$EOE")

    def test_non_numeric_fields_skipped(self):
        # Calendar-date string field is skipped without error.
        txt = (
            "$$SOE\n"
            "2451545.0, A.D. 2000-Jan-01 00:00, 1.2, -1.0, 0.4, 0.005, "
            "0.006, 0.001\n"
            "$$EOE"
        )
        sv = self._parse(txt)
        assert sv.x == 1.2 and sv.z == 0.4


# ===========================================================================
# clear_cache / shutdown (lines 338-343)
# ===========================================================================


class TestClearCache:
    def test_clear_and_shutdown(self):
        client = HorizonsClient()
        with (
            mock.patch.object(
                hb, "open_url", return_value=_FakeResp(_good_json_bytes())
            ),
            mock.patch.object(hb, "_get_ssl_context", return_value=None),
        ):
            client.fetch_state_vector("499", JD)
        assert len(client._cache) == 1
        client.clear_cache()
        assert len(client._cache) == 0
        # shutdown() also clears (idempotent).
        client.shutdown()
        assert len(client._cache) == 0


# ===========================================================================
# horizons_calc_ut — flag/body dispatch (lines 387, 396, 402+)
# ===========================================================================


class TestCalcUtDispatch:
    def test_topoctr_raises(self):
        with pytest.raises(KeyError):
            horizons_calc_ut(None, JD, 4, FLG_TOPOCTR)

    def test_nonut_raises(self):
        with pytest.raises(KeyError):
            horizons_calc_ut(None, JD, 4, FLG_NONUT)

    def test_analytical_dispatch(self):
        (data, fl) = horizons_calc_ut(None, JD, 10, FLG_SWIEPH)
        assert len(data) == 6

    def test_hypothetical_geocentric_uses_deterministic_fallback(self):
        with pytest.raises(KeyError, match="Skyfield fallback"):
            horizons_calc_ut(None, JD, 40, FLG_SWIEPH)

    def test_reviewed_hypothetical_helio_uses_local_model(self):
        data, _flags = horizons_calc_ut(
            None, JD, 40, FLG_SWIEPH | FLG_HELCTR | FLG_SPEED
        )
        assert len(data) == 6
        assert all(math.isfinite(value) for value in data)

    @pytest.mark.parametrize(
        "body",
        [49],
    )
    def test_unverified_heliocentric_hypothetical_ids_fail_closed(self, body):
        with pytest.raises(UnknownBodyError) as raised:
            horizons_calc_ut(None, JD, body, FLG_SWIEPH | FLG_HELCTR)
        assert raised.value.body_id == body

    @pytest.mark.parametrize("body", [*range(40, 49), *range(50, 56), 57])
    def test_reviewed_heliocentric_hypothetical_ids_are_finite(self, body):
        data, _flags = horizons_calc_ut(
            None, JD, body, FLG_SWIEPH | FLG_HELCTR | FLG_SPEED
        )
        assert len(data) == 6
        assert all(math.isfinite(value) for value in data)

    def test_reviewed_white_moon_native_geocentric_is_finite(self):
        data, _flags = horizons_calc_ut(None, JD, 56, FLG_SWIEPH | FLG_SPEED)
        assert len(data) == 6
        assert all(math.isfinite(value) for value in data)

    def test_reviewed_waldemath_native_geocentric_is_finite(self):
        data, _flags = horizons_calc_ut(None, JD, 58, FLG_SWIEPH | FLG_SPEED)
        assert len(data) == 6
        assert all(math.isfinite(value) for value in data)

    def test_harrington_geocentric_raises(self):
        with pytest.raises(KeyError):
            horizons_calc_ut(None, JD, HARRINGTON, FLG_SWIEPH)

    def test_harrington_helio_dispatch(self):
        (data, fl) = horizons_calc_ut(None, JD, HARRINGTON, FLG_SWIEPH | FLG_HELCTR)
        assert len(data) == 6

    def test_unknown_body_raises(self):
        # Body 99 not in command map / analytical / uranian / moon-derived.
        with pytest.raises(KeyError):
            horizons_calc_ut(None, JD, 99, FLG_SWIEPH)


# ===========================================================================
# horizons_calc_ut — heliocentric path (lines 437-484)
# ===========================================================================


class TestCalcUtHeliocentric:
    def test_sun_heliocentric_zero(self):
        client = FakeHorizonsClient()
        (data, fl) = horizons_calc_ut(client, JD, 0, FLG_SWIEPH | FLG_HELCTR)
        assert data == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        # No HTTP for Sun-from-Sun.
        assert client.calls == []

    def test_planet_heliocentric_lighttime(self):
        client = FakeHorizonsClient()
        (data, fl) = horizons_calc_ut(client, JD, 4, FLG_SWIEPH | FLG_HELCTR)
        assert len(data) == 6
        # Sun + target + light-time retarded target fetches happened.
        assert any(c[0] == "10" for c in client.calls)
        target_calls = [call for call in client.calls if call[0] == "499"]
        assert len(target_calls) == 4

        # The first retarded epoch follows the heliocentric distance
        # |Mars(t)-Sun(t)|/c, not either body's distance from the SSB.
        sun = FakeHorizonsClient._POS["10"]
        mars = FakeHorizonsClient._POS["499"]
        helio_distance = math.sqrt(
            (mars[0] - sun[0]) ** 2 + (mars[1] - sun[1]) ** 2 + (mars[2] - sun[2]) ** 2
        )
        expected_retarded_jd = target_calls[0][1] - helio_distance / 173.14463267
        assert target_calls[1][1] == pytest.approx(expected_retarded_jd, abs=1e-12)

    def test_planet_heliocentric_truepos_skips_lighttime(self):
        client = FakeHorizonsClient()
        horizons_calc_ut(client, JD, 4, FLG_SWIEPH | FLG_HELCTR | FLG_TRUEPOS)
        # No retarded refetch: exactly one target fetch at jd_tt.
        target_calls = [c for c in client.calls if c[0] == "499"]
        assert len(target_calls) == 1

    def test_barycentric_lighttime(self):
        client = FakeHorizonsClient()
        (data, fl) = horizons_calc_ut(client, JD, 4, FLG_SWIEPH | FLG_BARYCTR)
        assert len(data) == 6

    def test_barycentric_truepos(self):
        client = FakeHorizonsClient()
        (data, fl) = horizons_calc_ut(
            client, JD, 4, FLG_SWIEPH | FLG_BARYCTR | FLG_TRUEPOS
        )
        assert len(data) == 6

    def test_heliocentric_zero_distance_skips_lighttime(self):
        # Target coincides with the Sun -> rel_pos magnitude == 0, so the
        # ``if dist > 0.0`` light-time block is skipped (arc 459->474).
        class _Coincident(FakeHorizonsClient):
            def fetch_state_vector(self, command, jd, center="@0", time_type="TDB"):
                self.calls.append((command, jd, center, time_type))
                # Both Sun ("10") and target ("499") share the same position.
                return _sv(0.001, 0.002, 0.0005)

        client = _Coincident()
        (data, fl) = horizons_calc_ut(client, JD, 4, FLG_SWIEPH | FLG_HELCTR)
        assert len(data) == 6

    def test_barycentric_zero_distance_skips_lighttime(self):
        # Target at SSB origin -> dist == 0, light-time block skipped
        # (arc 481->484).
        class _AtOrigin(FakeHorizonsClient):
            def fetch_state_vector(self, command, jd, center="@0", time_type="TDB"):
                self.calls.append((command, jd, center, time_type))
                return _sv(0.0, 0.0, 0.0)

        client = _AtOrigin()
        (data, fl) = horizons_calc_ut(client, JD, 4, FLG_SWIEPH | FLG_BARYCTR)
        assert len(data) == 6


# ===========================================================================
# horizons_calc_ut — geocentric full pipeline (lines 488-599)
# ===========================================================================


class TestCalcUtGeocentric:
    def test_full_pipeline_no_speed(self):
        client = FakeHorizonsClient()
        (data, fl) = horizons_calc_ut(client, JD, 4, FLG_SWIEPH)
        assert len(data) == 6
        # speed components default to zero-ish output (no SPEED flag).

    def test_full_pipeline_with_speed(self):
        client = FakeHorizonsClient()
        (data, fl) = horizons_calc_ut(client, JD, 4, FLG_SWIEPH | FLG_SPEED)
        assert len(data) == 6

    def test_truepos_geocentric(self):
        client = FakeHorizonsClient()
        (data, fl) = horizons_calc_ut(client, JD, 4, FLG_SWIEPH | FLG_TRUEPOS)
        assert len(data) == 6

    def test_nogdefl_noaberr(self):
        client = FakeHorizonsClient()
        (data, fl) = horizons_calc_ut(
            client, JD, 4, FLG_SWIEPH | FLG_NOGDEFL | FLG_NOABERR
        )
        assert len(data) == 6

    def test_speed_with_truepos_and_no_corrections(self):
        client = FakeHorizonsClient()
        (data, fl) = horizons_calc_ut(
            client,
            JD,
            4,
            FLG_SWIEPH | FLG_SPEED | FLG_TRUEPOS | FLG_NOGDEFL | FLG_NOABERR,
        )
        assert len(data) == 6

    def test_target_or_earth_missing_raises_connectionerror(self):
        class _PartialClient(FakeHorizonsClient):
            def fetch_batch(self, requests, errors=None):
                # Return Earth but not the target -> ConnectionError path.
                out = {}
                for cmd, jd, center in requests:
                    if cmd == "399":
                        out[(cmd, jd, center)] = _sv(*self._POS["399"])
                    elif errors is not None:
                        errors[(cmd, jd, center)] = ConnectionError(
                            "simulated rate limit"
                        )
                return out

        with pytest.raises(ConnectionError):
            horizons_calc_ut(_PartialClient(), JD, 4, FLG_SWIEPH)


# ===========================================================================
# _to_ecliptic_output — frame matrix (lines 659-737)
# ===========================================================================


class TestToEclipticOutput:
    POS = (1.2, -1.0, 0.4)
    VEL = (0.005, 0.006, 0.001)

    def _run(self, iflag):
        return _to_ecliptic_output(self.POS, self.VEL, JD, JD, iflag)

    def test_icrs(self):
        (out, fl) = self._run(FLG_ICRS)
        assert len(out) == 6

    def test_j2000_equatorial(self):
        (out, fl) = self._run(FLG_J2000 | FLG_EQUATORIAL)
        assert len(out) == 6

    def test_j2000_ecliptic(self):
        (out, fl) = self._run(FLG_J2000)
        assert len(out) == 6

    def test_equatorial_of_date(self):
        (out, fl) = self._run(FLG_EQUATORIAL)
        assert len(out) == 6

    def test_default_ecliptic_of_date(self):
        (out, fl) = self._run(0)
        assert len(out) == 6

    def test_sidereal(self):
        (out, fl) = self._run(FLG_SIDEREAL)
        lon = out[0]
        assert 0.0 <= lon < 360.0

    def test_xyz_output(self):
        (out, fl) = self._run(FLG_XYZ)
        # XYZ returns pos+vel concatenation (6 values).
        assert len(out) == 6

    def test_radians_output(self):
        (out_deg, _) = self._run(0)
        (out_rad, _) = self._run(FLG_RADIANS)
        # Radians value differs from degrees value.
        assert out_rad[0] != out_deg[0]


# ===========================================================================
# _to_ecliptic_output — frame-speed conventions (reference-free)
# ===========================================================================


class TestToEclipticOutputSpeedConventions:
    """The reported speed must be the exact derivative of the reported position.

    Reference-free invariant (no network, no ephemeris files): feed a synthetic
    ICRS state, re-run _to_ecliptic_output at jd ± h with the analytically
    extrapolated state, and central-difference the reported longitude. This
    pins the post-rework frame-speed conventions and would fail with any of:
      (a) the former spurious general-precession subtraction from J2000
          spherical speeds (error = exactly -0.1377"/day at this epoch);
      (b) the double subtraction under SIDEREAL|J2000;
      (c) the missing nutation-rate term dΔψ/dt for of-date sidereal
          (error ~0.05"/day).
    """

    POS = (1.2, -1.0, 0.4)
    VEL = (0.005, 0.006, 0.001)
    JD_T = 2465712.625
    H = 0.01
    SID_MODE = 1  # Lahiri

    def _out(self, jd, pos, vel, iflag):
        (out, _) = _to_ecliptic_output(pos, vel, jd, jd, iflag, self.SID_MODE)
        return out

    def _invariant_err_arcsec(self, iflag):
        jd, h = self.JD_T, self.H
        rep = self._out(jd, self.POS, self.VEL, iflag)[3]
        pm = tuple(self.POS[i] - self.VEL[i] * h for i in range(3))
        pp = tuple(self.POS[i] + self.VEL[i] * h for i in range(3))
        lon_m = self._out(jd - h, pm, self.VEL, iflag)[0]
        lon_p = self._out(jd + h, pp, self.VEL, iflag)[0]
        dl = (lon_p - lon_m) % 360.0
        if dl > 180.0:
            dl -= 360.0
        return (rep - dl / (2.0 * h)) * 3600.0

    @pytest.mark.parametrize(
        "label, iflag",
        [
            ("tropical_of_date", FLG_SPEED),
            ("j2000", FLG_SPEED | FLG_J2000),
            ("sidereal_of_date", FLG_SPEED | FLG_SIDEREAL),
            ("sidereal_j2000", FLG_SPEED | FLG_SIDEREAL | FLG_J2000),
            ("equatorial_of_date", FLG_SPEED | FLG_EQUATORIAL),
            ("sidereal_equatorial", FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL),
            ("equatorial_j2000", FLG_SPEED | FLG_EQUATORIAL | FLG_J2000),
        ],
    )
    def test_speed_is_derivative_of_reported_position(self, label, iflag):
        err = self._invariant_err_arcsec(iflag)
        assert abs(err) < 0.005, (
            f"{label}: reported dlon differs from the derivative of the "
            f'reported lon by {err:+.4f}"/day'
        )

    def test_sid_j2000_single_ayanamsa_drift(self):
        """SID|J2000 dlon = J2000 dlon − general precession rate, exactly once."""
        from libephemeris.fast_calc import _general_precession_rate_deg_day

        d_j2000 = self._out(self.JD_T, self.POS, self.VEL, FLG_SPEED | FLG_J2000)[3]
        d_sid = self._out(
            self.JD_T, self.POS, self.VEL, FLG_SPEED | FLG_SIDEREAL | FLG_J2000
        )[3]
        p = _general_precession_rate_deg_day(self.JD_T)
        assert abs((d_j2000 - d_sid) - p) * 3600.0 < 1e-6, (
            "SID|J2000 must subtract the ayanamsa drift exactly once "
            f"(got {(d_j2000 - d_sid) * 3600.0:.6f} arcsec/day vs p = "
            f"{p * 3600.0:.6f})"
        )


# ===========================================================================
# _HorizonsDeflectorSource (lines 616-627) + _apply_deflection_horizons
# (lines 647-656)
# ===========================================================================


class TestDeflectorSource:
    def test_eval_body_batch_hit(self):
        batch = {("10", JD, "@0"): _sv(0.001, 0.002, 0.0005)}
        src = _HorizonsDeflectorSource(batch, FakeHorizonsClient(), JD)
        pos, vel = src.eval_body(0, JD)  # body 0 -> command "10"
        assert pos == (0.001, 0.002, 0.0005)

    def test_eval_body_batch_miss_falls_back_to_client(self):
        # jd matches batch epoch but the command isn't in the batch dict.
        batch = {}
        client = FakeHorizonsClient()
        src = _HorizonsDeflectorSource(batch, client, JD)
        pos, vel = src.eval_body(5, JD)  # Jupiter -> command "5"
        assert pos == client._POS["5"][:3]

    def test_eval_body_different_jd_fetches(self):
        # jd != batch epoch -> always fetch from client.
        batch = {("10", JD, "@0"): _sv(9.0, 9.0, 9.0)}
        client = FakeHorizonsClient()
        src = _HorizonsDeflectorSource(batch, client, JD)
        pos, vel = src.eval_body(0, JD + 5.0)
        assert pos == client._POS["10"][:3]

    def test_eval_body_unknown_deflector_raises(self):
        src = _HorizonsDeflectorSource({}, FakeHorizonsClient(), JD)
        with pytest.raises(KeyError):
            src.eval_body(7, JD)  # 7 not a deflector


class TestApplyDeflectionHorizons:
    def test_normal_deflection(self):
        client = FakeHorizonsClient()
        batch = {
            ("10", JD, "@0"): _sv(*client._POS["10"]),
            ("5", JD, "@0"): _sv(*client._POS["5"]),
            ("6", JD, "@0"): _sv(*client._POS["6"]),
        }
        geo = (1.2 - 0.5, -1.0 - 0.8, 0.4 - 0.3)
        earth_bary = (0.5, 0.8, 0.3)
        out = _apply_deflection_horizons(geo, earth_bary, JD, 0.01, batch, client)
        assert len(out) == 3

    def test_connectionerror_returns_undeflected(self):
        class _Boom(FakeHorizonsClient):
            def fetch_state_vector(self, *a, **k):
                raise ConnectionError("deflector fetch failed")

        client = _Boom()
        geo = (0.7, -1.8, 0.1)
        earth_bary = (0.5, 0.8, 0.3)
        # Empty batch + different jd forces the source to call the client,
        # which raises ConnectionError -> warning + undeflected geo returned.
        out = _apply_deflection_horizons(geo, earth_bary, JD, 0.01, {}, client)
        assert out == geo


# ===========================================================================
# Independently derived analytical bodies
# ===========================================================================


class TestCalcAnalytical:
    def test_mean_node(self):
        (data, fl) = _calc_analytical(JD, 10, FLG_SWIEPH | FLG_SPEED)
        assert len(data) == 6
        assert data[3] < 0  # retrograde mean node

    def test_mean_apogee(self):
        (data, fl) = _calc_analytical(JD, 12, FLG_SWIEPH | FLG_SPEED)
        assert len(data) == 6

    def test_mean_node_sidereal(self):
        (data, fl) = _calc_analytical(JD, 10, FLG_SWIEPH | FLG_SIDEREAL)
        assert 0.0 <= data[0] < 360.0

    def test_non_analytical_body_raises(self):
        with pytest.raises(KeyError):
            _calc_analytical(JD, 99, FLG_SWIEPH)

    def test_mean_node_wrap_branch(self):
        # Force the speed wrap branch by patching the node function so the
        # ±0.5-day central-difference samples straddle the 0/360 boundary.
        # _calc_analytical now calls calc_mean_lunar_node three times:
        # position (jd_tt), then jd_tt-dt and jd_tt+dt for the speed.
        seq = iter([0.0, 359.9, 0.1])
        with mock.patch.object(hb, "_calc_analytical", wraps=hb._calc_analytical):
            import libephemeris.lunar as lunar

            with mock.patch.object(
                lunar, "calc_mean_lunar_node", side_effect=lambda jd: next(seq)
            ):
                (data, fl) = _calc_analytical(JD, 10, FLG_SWIEPH | FLG_SPEED)
        # Wrapped speed stays small, not a huge ~360 deg/day jump: the samples
        # 359.9 -> 0.1 must wrap to +0.2 deg over the day, so speed_lon ~ 0.2.
        # Without the wrap the raw diff (-359.8) would blow this bound.
        assert len(data) == 6
        assert abs(data[3]) < 1.0


class TestSiderealSpeedBranchResolvesGlobalMode:
    """The star-anchored speed branch must resolve sid_mode=None to the
    global mode. The module-level Horizons dispatch passes no sid_mode, so an
    unresolved None sent star-based modes down the generic precession branch.
    """

    POS = (1.2, -1.0, 0.4)
    VEL = (0.005, 0.006, 0.001)

    def test_none_matches_explicit_star_mode(self):
        import libephemeris as le
        from libephemeris.constants import SIDM_TRUE_CITRA

        le.set_sid_mode(SIDM_TRUE_CITRA)
        try:
            iflag = FLG_SWIEPH | FLG_SIDEREAL | FLG_SPEED
            (out_none, _) = _to_ecliptic_output(self.POS, self.VEL, JD, JD, iflag, None)
            (out_explicit, _) = _to_ecliptic_output(
                self.POS, self.VEL, JD, JD, iflag, SIDM_TRUE_CITRA
            )
            # Same branch, same ayanamsha: identical longitude AND speed.
            assert abs(out_none[0] - out_explicit[0]) < 1e-12
            assert abs(out_none[3] - out_explicit[3]) < 1e-12
        finally:
            le.set_sid_mode(0)

    def test_entry_snapshots_global_mode(self):
        """horizons_calc_ut resolves sid_mode at entry for sidereal requests
        (nodes go through _calc_analytical without HTTP)."""
        import libephemeris as le
        from libephemeris.constants import MEAN_NODE, SIDM_TRUE_CITRA

        le.set_sid_mode(SIDM_TRUE_CITRA)
        try:
            iflag = FLG_SWIEPH | FLG_SIDEREAL | FLG_SPEED
            (d_none, _) = hb.horizons_calc_ut(None, JD, MEAN_NODE, iflag)
            (d_expl, _) = hb.horizons_calc_ut(
                None, JD, MEAN_NODE, iflag, sid_mode=SIDM_TRUE_CITRA
            )
            assert abs(d_none[0] - d_expl[0]) < 1e-12
            assert abs(d_none[3] - d_expl[3]) < 1e-12
        finally:
            le.set_sid_mode(0)
