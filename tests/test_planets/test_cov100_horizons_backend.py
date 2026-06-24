# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
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
   ``_parse_response``) is driven by patching ``urllib.request.urlopen`` with
   a fake context manager whose ``.read()`` returns canned JSON bytes, and by
   patching ``time.sleep`` so retry backoff is instantaneous.
"""

from __future__ import annotations

import json
from unittest import mock

import pytest

from libephemeris import horizons_backend as hb
from libephemeris.horizons_backend import (
    HorizonsClient,
    StateVector,
    _HorizonsDeflectorSource,
    _apply_deflection_horizons,
    _calc_analytical,
    _calc_uranian,
    _get_ssl_context,
    _to_ecliptic_output,
    horizons_calc_ut,
)

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

    def fetch_batch(self, requests):
        out = {}
        for cmd, jd, center in requests:
            p = self._POS.get(cmd)
            if p is not None:
                out[(cmd, jd, center)] = _sv(*p)
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
    def test_builds_once_and_caches(self):
        # Reset module-level cache so the build branch is exercised.
        hb._SSL_CTX = None
        sentinel = object()
        fake_ssl = mock.MagicMock()
        fake_ssl.create_default_context.return_value = sentinel
        fake_certifi = mock.MagicMock()
        fake_certifi.where.return_value = "/fake/ca.pem"
        with mock.patch.dict(
            "sys.modules", {"ssl": fake_ssl, "certifi": fake_certifi}
        ):
            ctx1 = _get_ssl_context()
            assert ctx1 is sentinel
            # Second call returns the cached object without rebuilding.
            ctx2 = _get_ssl_context()
            assert ctx2 is sentinel
        fake_ssl.create_default_context.assert_called_once()

    def test_double_checked_lock_already_built(self):
        # Exercise the inner "already built" branch (122->128): the outer
        # check sees None, but by the time the lock is acquired another
        # thread has populated _SSL_CTX. Simulate with a fake lock whose
        # __enter__ sets the module global.
        hb._SSL_CTX = None
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


# ===========================================================================
# HorizonsClient.fetch_state_vector + cache (lines 148-208)
# ===========================================================================


class TestFetchStateVector:
    def test_fetch_parses_and_caches(self):
        client = HorizonsClient(max_cache_size=4)
        with mock.patch.object(
            hb.urllib.request, "urlopen", return_value=_FakeResp(_good_json_bytes())
        ) as op, mock.patch.object(hb, "_get_ssl_context", return_value=None):
            sv1 = client.fetch_state_vector("499", JD, "@0", "TDB")
            assert sv1.x == 1.2 and sv1.vz == 0.001
            # Second identical call hits the cache (no new urlopen).
            sv2 = client.fetch_state_vector("499", JD, "@0", "TDB")
            assert sv2 is sv1
        assert op.call_count == 1

    def test_cache_lru_eviction(self):
        client = HorizonsClient(max_cache_size=1)
        with mock.patch.object(
            hb.urllib.request, "urlopen", return_value=_FakeResp(_good_json_bytes())
        ), mock.patch.object(hb, "_get_ssl_context", return_value=None):
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
        with mock.patch.object(
            hb.urllib.request, "urlopen", return_value=_FakeResp(_good_json_bytes())
        ), mock.patch.object(hb, "_get_ssl_context", return_value=None):
            # Prime the cache with one entry.
            client.fetch_state_vector("499", JD, "@0", "TDB")
            res = client.fetch_batch(
                [("499", JD, "@0"), ("10", JD, "@0")]
            )
        assert ("499", JD, "@0") in res
        assert ("10", JD, "@0") in res

    def test_batch_all_cached_short_circuit(self):
        client = HorizonsClient()
        with mock.patch.object(
            hb.urllib.request, "urlopen", return_value=_FakeResp(_good_json_bytes())
        ), mock.patch.object(hb, "_get_ssl_context", return_value=None):
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
        with mock.patch.object(
            hb.urllib.request, "urlopen", return_value=_FakeResp(payload)
        ) as op, mock.patch.object(hb, "_get_ssl_context", return_value=None):
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

        with mock.patch.object(
            hb.urllib.request, "urlopen", side_effect=_urlopen
        ), mock.patch.object(hb, "_get_ssl_context", return_value=None), mock.patch(
            "time.sleep", return_value=None
        ):
            sv = client.fetch_state_vector("499", JD)
        assert sv.x == 1.2

    def test_all_attempts_fail_raises_connectionerror(self):
        client = HorizonsClient()
        with mock.patch.object(
            hb.urllib.request, "urlopen", side_effect=OSError("down")
        ), mock.patch.object(hb, "_get_ssl_context", return_value=None), mock.patch(
            "time.sleep", return_value=None
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
        with mock.patch.object(
            hb.urllib.request, "urlopen", return_value=_FakeResp(_good_json_bytes())
        ), mock.patch.object(hb, "_get_ssl_context", return_value=None):
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

    def test_uranian_geocentric_raises(self):
        # Uranian body without HELCTR -> KeyError.
        with pytest.raises(KeyError):
            horizons_calc_ut(None, JD, 40, FLG_SWIEPH)

    def test_uranian_helio_dispatch(self):
        (data, fl) = horizons_calc_ut(None, JD, 40, FLG_SWIEPH | FLG_HELCTR)
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
        assert any(c[0] == "499" for c in client.calls)

    def test_planet_heliocentric_truepos_skips_lighttime(self):
        client = FakeHorizonsClient()
        horizons_calc_ut(
            client, JD, 4, FLG_SWIEPH | FLG_HELCTR | FLG_TRUEPOS
        )
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
            def fetch_batch(self, requests):
                # Return Earth but not the target -> ConnectionError path.
                out = {}
                for cmd, jd, center in requests:
                    if cmd == "399":
                        out[(cmd, jd, center)] = _sv(*self._POS["399"])
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
# _calc_analytical (lines 740-789) + _calc_uranian (792-814)
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
        with mock.patch.object(
            hb, "_calc_analytical", wraps=hb._calc_analytical
        ):
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


class TestCalcUranian:
    def test_transpluto(self):
        (data, fl) = _calc_uranian(JD, 48, FLG_SWIEPH)
        assert len(data) == 6

    def test_uranian_planet(self):
        (data, fl) = _calc_uranian(JD, 40, FLG_SWIEPH)
        assert len(data) == 6

    def test_uranian_sidereal(self):
        (data, fl) = _calc_uranian(JD, 40, FLG_SWIEPH | FLG_SIDEREAL)
        assert 0.0 <= data[0] < 360.0
