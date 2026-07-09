# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
NASA JPL Horizons API backend for zero-install ephemeris.

Fetches state vectors from the Horizons REST API and computes apparent
positions using the same correction pipeline as the LEB/Skyfield paths
(light-time iteration, gravitational deflection, stellar aberration,
frame rotation).

This backend is used when:
- mode="horizons" — always use Horizons
- mode="auto" — no LEB file AND no DE440 locally available

Bodies not supported by Horizons (fixed stars, planetary moons, FLG_TOPOCTR)
raise KeyError to trigger Skyfield fallback.
"""

from __future__ import annotations

import json
import logging
import threading
import urllib.error
import urllib.request
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger("libephemeris")

# =============================================================================
# CONSTANTS
# =============================================================================

API_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"

# Map SE_* body IDs to Horizons COMMAND strings
_HORIZONS_COMMAND: Dict[int, str] = {
    0: "10",  # Sun
    1: "301",  # Moon
    2: "199",  # Mercury
    3: "299",  # Venus
    4: "499",  # Mars
    5: "599",  # Jupiter
    6: "699",  # Saturn
    7: "799",  # Uranus
    8: "899",  # Neptune
    9: "999",  # Pluto
    14: "399",  # Earth
    15: "2060;",  # Chiron (small body)
    # 16: "5145;", # Pholus
    17: "1;",  # Ceres
    18: "2;",  # Pallas
    19: "3;",  # Juno
    20: "4;",  # Vesta
}

# Bodies computed analytically (no HTTP needed)
_ANALYTICAL_BODIES = {10, 12}  # Mean Node, Mean Apogee

# Bodies computed from Moon state vectors (single HTTP for Moon)
_MOON_DERIVED_BODIES = {11, 13, 21, 22}  # True Node, Oscu Apogee, Interp Apogee/Perigee

# Heliocentric-only bodies (Uranians) — analytical, no HTTP
_URANIAN_BODIES = {40, 41, 42, 43, 44, 45, 46, 47, 48}

# Deflector bodies for gravitational light bending
_DEFLECTOR_NAIF = {
    0: "10",  # Sun
    5: "5",  # Jupiter barycenter
    6: "6",  # Saturn barycenter
}


# =============================================================================
# STATE VECTOR DATACLASS
# =============================================================================


class StateVector:
    """Barycentric ICRS state vector from Horizons."""

    __slots__ = ("x", "y", "z", "vx", "vy", "vz")

    def __init__(self, x: float, y: float, z: float, vx: float, vy: float, vz: float):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz

    @property
    def pos(self) -> Tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def vel(self) -> Tuple[float, float, float]:
        return (self.vx, self.vy, self.vz)


# =============================================================================
# HTTP CLIENT
# =============================================================================


_SSL_CTX = None
_SSL_CTX_LOCK = threading.Lock()


def _get_ssl_context():
    """Return a shared default SSL context, building it once on first use.

    Parsing the certifi CA bundle (disk read + parse) is non-trivial and the
    Horizons path issues many sequential fetches per body, so the context is
    cached at module level rather than rebuilt on every request.
    """
    global _SSL_CTX
    if _SSL_CTX is None:
        with _SSL_CTX_LOCK:
            if _SSL_CTX is None:
                import ssl

                import certifi

                _SSL_CTX = ssl.create_default_context(cafile=certifi.where())
    return _SSL_CTX


class HorizonsClient:
    """HTTP client for NASA JPL Horizons API with LRU cache.

    Thread-safe. Caches state vectors keyed by (jd, command, center).
    """

    def __init__(
        self,
        max_cache_size: int = 4096,
        timeout: int = 30,
        max_workers: int = 8,
    ):
        self._cache: OrderedDict = OrderedDict()
        self._cache_lock = threading.Lock()
        self._max_cache_size = max_cache_size
        self._timeout = timeout
        self._max_workers = max_workers

    def fetch_state_vector(
        self,
        command: str,
        jd: float,
        center: str = "@0",
        time_type: str = "TDB",
    ) -> StateVector:
        """Fetch a barycentric ICRS state vector from Horizons.

        Args:
            command: Horizons body command (e.g. "499" for Mars).
            jd: Julian Day.
            center: Horizons center (default "@0" = solar system barycenter).
            time_type: "TDB" or "UT".

        Returns:
            StateVector with position (AU) and velocity (AU/day).

        Raises:
            ConnectionError: Network/API error after retries.
            KeyError: Body not found on Horizons.
        """
        cache_key = (round(jd, 12), command, center, time_type)

        with self._cache_lock:
            if cache_key in self._cache:
                self._cache.move_to_end(cache_key)
                return self._cache[cache_key]

        # Build URL
        params = {
            "format": "json",
            "COMMAND": f"'{command}'",
            "EPHEM_TYPE": "'VECTORS'",
            "CENTER": f"'{center}'",
            "TLIST": str(jd),
            "TLIST_TYPE": "'JD'",
            "VEC_TABLE": "'2'",
            "OUT_UNITS": "'AU-D'",
            "VEC_CORR": "'NONE'",
            "CSV_FORMAT": "'YES'",
            "REF_SYSTEM": "'ICRF'",
            # x-y plane = Earth mean equator (ICRF); without this Horizons
            # defaults to REF_PLANE=ECLIPTIC and the whole pipeline (which
            # assumes equatorial ICRS vectors) would be misframed by ~23.4°.
            "REF_PLANE": "'FRAME'",
            "TIME_TYPE": f"'{time_type}'",
        }
        query = "&".join(f"{k}={v}" for k, v in params.items())
        url = f"{API_URL}?{query}"

        # Fetch with retry
        sv = self._fetch_with_retry(url, command)

        with self._cache_lock:
            self._cache[cache_key] = sv
            if len(self._cache) > self._max_cache_size:
                self._cache.popitem(last=False)

        return sv

    def fetch_batch(
        self,
        requests: List[Tuple[str, float, str]],
        errors: Optional[Dict[Tuple[str, float, str], Exception]] = None,
    ) -> Dict[Tuple[str, float, str], StateVector]:
        """Fetch multiple state vectors in parallel.

        Args:
            requests: List of (command, jd, center) tuples.
            errors: Optional dict populated with the exception of each
                failed fetch, keyed like the results — lets callers chain
                the real cause instead of raising a generic failure.

        Returns:
            Dict mapping (command, jd, center) -> StateVector.
        """
        # Filter out cached results.
        # The batch interface is TDB-only (all callers pass jd_tt): the
        # cache key hardcodes the time type on purpose. A future UT-based
        # caller must extend the request tuples, NOT reuse this method,
        # or it would read TDB-cached vectors as if they were UT.
        to_fetch = []
        results = {}
        for cmd, jd, center in requests:
            cache_key = (round(jd, 12), cmd, center, "TDB")
            with self._cache_lock:
                if cache_key in self._cache:
                    self._cache.move_to_end(cache_key)
                    results[(cmd, jd, center)] = self._cache[cache_key]
                else:
                    to_fetch.append((cmd, jd, center))

        if not to_fetch:
            return results

        with ThreadPoolExecutor(max_workers=self._max_workers) as pool:
            futures = {}
            for cmd, jd, center in to_fetch:
                fut = pool.submit(self.fetch_state_vector, cmd, jd, center)
                futures[fut] = (cmd, jd, center)

            for fut in as_completed(futures):
                key = futures[fut]
                try:
                    results[key] = fut.result()
                except (OSError, ValueError, KeyError) as exc:
                    # Skip the failed fetch and record the post-retry cause so
                    # the caller can chain it. Log at DEBUG only: a missing
                    # deflector (Sun/Jupiter/Saturn) is recovered gracefully by
                    # _apply_deflection_horizons, so a WARNING here would be
                    # spurious on the normal path — horizons_calc_ut raises the
                    # aggregated WARNING/error at the decision point when
                    # target/Earth are actually missing.
                    if errors is not None:
                        errors[key] = exc
                    logger.debug(
                        "Horizons batch fetch failed for %s @ JD %.6f (%s): %s",
                        key[0],
                        key[1],
                        key[2],
                        exc,
                    )

        return results

    def _fetch_with_retry(
        self, url: str, command: str, max_retries: int = 2
    ) -> StateVector:
        """Fetch URL with exponential backoff retry."""
        import time

        ssl_ctx = _get_ssl_context()

        last_err = None
        for attempt in range(max_retries + 1):
            try:
                req = urllib.request.Request(
                    url,
                    headers={"User-Agent": "libephemeris/1.0"},
                )
                with urllib.request.urlopen(
                    req, timeout=self._timeout, context=ssl_ctx
                ) as resp:
                    data = json.loads(resp.read().decode("utf-8"))

                if "error" in data:
                    raise KeyError(f"Horizons API error for {command}: {data['error']}")

                return self._parse_response(data, command)

            except KeyError:
                raise  # don't retry API errors (body not found etc.)
            except (OSError, ValueError, KeyError) as e:
                last_err = e
                if attempt < max_retries:
                    time.sleep(0.5 * (2**attempt))

        raise ConnectionError(
            f"Horizons API request failed after {max_retries + 1} attempts "
            f"for {command}: {last_err}. "
            f"Check internet connection, or download local ephemeris: "
            f"python -m libephemeris download"
        )

    def _parse_response(self, data: dict, command: str) -> StateVector:
        """Parse Horizons JSON response to extract state vector."""
        result_text = data.get("result", "")

        # Find data between $$SOE and $$EOE markers
        soe_idx = result_text.find("$$SOE")
        eoe_idx = result_text.find("$$EOE")
        if soe_idx < 0 or eoe_idx < 0:
            raise ValueError(
                f"Cannot find $$SOE/$$EOE markers in Horizons response for {command}"
            )

        data_block = result_text[soe_idx + 5 : eoe_idx].strip()
        # VEC_TABLE='2' + CSV_FORMAT='YES' fixes the record shape:
        # JDTDB, Calendar Date, X, Y, Z, VX, VY, VZ,
        lines = [line.strip() for line in data_block.split("\n") if line.strip()]
        if not lines:
            raise ValueError(f"Empty data block in Horizons response for {command}")

        # Parse the first (and only, single-epoch TLIST) data record
        # positionally by CSV field: accumulating every parseable float
        # across the block would let a stray numeric token shift the
        # coordinate indices silently.
        fields = [f.strip() for f in lines[0].split(",")]
        if len(fields) < 8:
            raise ValueError(
                f"Malformed data record in Horizons response for {command}: "
                f"got {len(fields)} CSV fields, need at least 8"
            )
        try:
            x, y, z, vx, vy, vz = (float(fields[i]) for i in range(2, 8))
        except ValueError as exc:
            raise ValueError(
                f"Malformed data record in Horizons response for {command}: "
                f"non-numeric state-vector field ({exc})"
            ) from exc

        return StateVector(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)

    def clear_cache(self) -> None:
        with self._cache_lock:
            self._cache.clear()

    def shutdown(self) -> None:
        self.clear_cache()


# =============================================================================
# CALCULATION PIPELINE
# =============================================================================


def _get_body_command(body_id: int) -> str:
    """Get Horizons COMMAND for a body, or raise KeyError."""
    if body_id in _HORIZONS_COMMAND:
        return _HORIZONS_COMMAND[body_id]
    raise KeyError(f"Body {body_id} not supported by Horizons backend")


def _is_horizons_body(body_id: int) -> bool:
    """Check if a body can be computed via Horizons (directly or analytically)."""
    return (
        body_id in _HORIZONS_COMMAND
        or body_id in _ANALYTICAL_BODIES
        or body_id in _MOON_DERIVED_BODIES
        or body_id in _URANIAN_BODIES
    )


def horizons_calc_ut(
    client: HorizonsClient,
    jd_ut: float,
    body_id: int,
    iflag: int,
    sid_mode: int | None = None,
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """Calculate body position via Horizons API.

    Same return format as fast_calc.fast_calc_ut().

    Args:
        client: HorizonsClient instance.
        jd_ut: Julian Day UT.
        body_id: SE_* body constant.
        iflag: The reference-ephemeris calculation flags.
        sid_mode: Sidereal mode override; reads the global state when None.
            An explicit value lets EphemerisContext dispatch through
            Horizons without swapping the global sidereal mode across the
            network call (thread-safe, like fast_calc's sid_mode).

    Returns:
        ((lon, lat, dist, dlon, dlat, ddist), iflag)

    Raises:
        KeyError: Body not supported by Horizons.
    """
    from .constants import (
        FLG_HELCTR,
        FLG_BARYCTR,
        FLG_NOABERR,
        FLG_NOGDEFL,
        FLG_SPEED,
        FLG_TOPOCTR,
        FLG_TRUEPOS,
        FLG_NONUT,
    )

    # Unsupported flags → fallback to Skyfield
    if iflag & FLG_TOPOCTR:
        raise KeyError("FLG_TOPOCTR not supported by Horizons backend")
    # Analytical bodies — no HTTP needed. _calc_analytical is FLG_NONUT-aware
    # (it conditionally omits Δψ), so it must be dispatched BEFORE the generic
    # FLG_NONUT rejection below, otherwise that handling is unreachable.
    if body_id in _ANALYTICAL_BODIES:
        return _calc_analytical(jd_ut, body_id, iflag, sid_mode)

    if iflag & FLG_NONUT:
        # The state-vector frames below always include nutation; defer to
        # Skyfield (same policy as the LEB backend).
        raise KeyError("FLG_NONUT not supported by Horizons backend")

    # Uranian hypotheticals — analytical, heliocentric only
    if body_id in _URANIAN_BODIES:
        if not (iflag & FLG_HELCTR):
            raise KeyError(
                f"Uranian body {body_id} geocentric not supported via Horizons"
            )
        return _calc_uranian(jd_ut, body_id, iflag, sid_mode)

    # Bodies not in Horizons command map
    if body_id not in _HORIZONS_COMMAND:
        raise KeyError(f"Body {body_id} not in Horizons command map")

    # Convert UT to TT (approximate, good enough for Horizons queries)
    from .time_utils import deltat

    delta_t = deltat(jd_ut)
    jd_tt = jd_ut + delta_t

    command = _HORIZONS_COMMAND[body_id]

    # Heliocentric or barycentric — simpler pipeline (no aberration or
    # deflection, but light-time retardation still applies unless TRUEPOS:
    # the target is evaluated at the retarded epoch while the observer
    # stays at observation time, matching the LEB/Skyfield pipeline).
    if iflag & FLG_HELCTR or iflag & FLG_BARYCTR:
        import math as _math

        c_au_day = 173.14463267
        if iflag & FLG_HELCTR:
            if body_id == 0:
                # Sun heliocentric = Sun from Sun = (0, 0, 0)
                return ((0.0, 0.0, 0.0, 0.0, 0.0, 0.0), iflag)
            sun_sv = client.fetch_state_vector("10", jd_tt, "@0", "TDB")
            tgt_sv = client.fetch_state_vector(command, jd_tt, "@0", "TDB")
            rel_pos = (
                tgt_sv.x - sun_sv.x,
                tgt_sv.y - sun_sv.y,
                tgt_sv.z - sun_sv.z,
            )
            rel_vel = (
                tgt_sv.vx - sun_sv.vx,
                tgt_sv.vy - sun_sv.vy,
                tgt_sv.vz - sun_sv.vz,
            )
            if not (iflag & FLG_TRUEPOS):
                dist = _math.sqrt(rel_pos[0] ** 2 + rel_pos[1] ** 2 + rel_pos[2] ** 2)
                if dist > 0.0:
                    lt = dist / c_au_day
                    tgt_lt = client.fetch_state_vector(command, jd_tt - lt, "@0", "TDB")
                    rel_pos = (
                        tgt_lt.x - sun_sv.x,
                        tgt_lt.y - sun_sv.y,
                        tgt_lt.z - sun_sv.z,
                    )
                    rel_vel = (
                        tgt_lt.vx - sun_sv.vx,
                        tgt_lt.vy - sun_sv.vy,
                        tgt_lt.vz - sun_sv.vz,
                    )
            # Speed slots are 0.0 without FLG_SPEED (matching the reference
            # ephemeris and the geocentric branch below, which zeroes velocity).
            if not (iflag & FLG_SPEED):
                rel_vel = (0.0, 0.0, 0.0)
            return _to_ecliptic_output(rel_pos, rel_vel, jd_tt, jd_ut, iflag, sid_mode)

        # Barycentric: the SSB is inertial, so retardation is just an
        # earlier evaluation epoch.
        sv = client.fetch_state_vector(command, jd_tt, "@0", "TDB")
        if not (iflag & FLG_TRUEPOS):
            dist = _math.sqrt(sv.x**2 + sv.y**2 + sv.z**2)
            if dist > 0.0:
                lt = dist / c_au_day
                sv = client.fetch_state_vector(command, jd_tt - lt, "@0", "TDB")
        # Speed slots are 0.0 without FLG_SPEED (matching the reference
        # ephemeris and the geocentric branch below, which passes a zero
        # velocity in that case).
        bary_vel = sv.vel if (iflag & FLG_SPEED) else (0.0, 0.0, 0.0)
        return _to_ecliptic_output(sv.pos, bary_vel, jd_tt, jd_ut, iflag, sid_mode)

    # Deflection is skipped for true geometric positions and when explicitly
    # disabled; only prefetch the deflectors when it will actually run.
    apply_deflection = not (iflag & FLG_NOGDEFL) and not (iflag & FLG_TRUEPOS)

    # Geocentric apparent — full pipeline
    # Prefetch: target + Earth (always); deflectors (Sun, Jupiter, Saturn) only
    # when deflection will run.
    prefetch_cmds = [
        (command, jd_tt, "@0"),  # target barycentric
        ("399", jd_tt, "@0"),  # Earth barycentric
    ]
    if apply_deflection:
        prefetch_cmds += [
            ("10", jd_tt, "@0"),  # Sun barycentric (deflector)
            ("5", jd_tt, "@0"),  # Jupiter barycenter (deflector)
            ("6", jd_tt, "@0"),  # Saturn barycenter (deflector)
        ]
    batch_errors: Dict[Tuple[str, float, str], Exception] = {}
    batch = client.fetch_batch(prefetch_cmds, errors=batch_errors)

    target_sv = batch.get((command, jd_tt, "@0"))
    earth_sv = batch.get(("399", jd_tt, "@0"))

    if target_sv is None or earth_sv is None:
        # Chain the real post-retry cause (rate limit, DNS, parse error)
        # instead of hiding it behind a generic failure message.
        cause = batch_errors.get((command, jd_tt, "@0")) or batch_errors.get(
            ("399", jd_tt, "@0")
        )
        msg = f"Failed to fetch target/Earth for body {body_id}"
        if cause is not None:
            msg = f"{msg}: {cause}"
        raise ConnectionError(msg) from cause

    # Geometric geocentric
    geo = (
        target_sv.x - earth_sv.x,
        target_sv.y - earth_sv.y,
        target_sv.z - earth_sv.z,
    )

    # Light-time correction (single iteration, sufficient for arcsecond precision)
    if not (iflag & FLG_TRUEPOS):
        import math

        c_au_day = 173.14463267  # speed of light in AU/day
        dist = math.sqrt(geo[0] ** 2 + geo[1] ** 2 + geo[2] ** 2)
        lt = dist / c_au_day

        # Re-fetch target at retarded time
        target_lt = client.fetch_state_vector(command, jd_tt - lt, "@0", "TDB")
        geo = (
            target_lt.x - earth_sv.x,
            target_lt.y - earth_sv.y,
            target_lt.z - earth_sv.z,
        )
        dist = math.sqrt(geo[0] ** 2 + geo[1] ** 2 + geo[2] ** 2)
        lt = dist / c_au_day
    else:
        lt = 0.0

    # Gravitational deflection (suppressed for true geometric positions, matching
    # fast_calc and the reference ephemeris, which gate deflection on both
    # NOGDEFL and TRUEPOS — like the light-time and aberration steps).
    if apply_deflection:
        geo = _apply_deflection_horizons(geo, earth_sv.pos, jd_tt, lt, batch, client)

    # Aberration (relativistic when light-time is known, matching fast_calc)
    if not (iflag & FLG_NOABERR) and not (iflag & FLG_TRUEPOS):
        from .fast_calc import _apply_aberration

        geo = _apply_aberration(geo, earth_sv.vel, lt)

    # Speeds only on request — they cost several extra HTTP round-trips
    if not (iflag & FLG_SPEED):
        return _to_ecliptic_output(geo, (0.0, 0.0, 0.0), jd_tt, jd_ut, iflag, sid_mode)

    # Velocity via numerical derivative of the apparent position
    # Compute position at jd + dt to get d(apparent_pos)/dt
    dt = 1.0 / 86400.0  # 1 second in days
    jd_tt2 = jd_tt + dt

    target_sv2 = client.fetch_state_vector(command, jd_tt2, "@0", "TDB")
    earth_sv2 = client.fetch_state_vector("399", jd_tt2, "@0", "TDB")

    geo2 = (
        target_sv2.x - earth_sv2.x,
        target_sv2.y - earth_sv2.y,
        target_sv2.z - earth_sv2.z,
    )

    # Apply same corrections to geo2
    lt2 = 0.0
    if not (iflag & FLG_TRUEPOS):
        dist2 = math.sqrt(geo2[0] ** 2 + geo2[1] ** 2 + geo2[2] ** 2)
        lt2 = dist2 / c_au_day
        target_lt2 = client.fetch_state_vector(command, jd_tt2 - lt2, "@0", "TDB")
        geo2 = (
            target_lt2.x - earth_sv2.x,
            target_lt2.y - earth_sv2.y,
            target_lt2.z - earth_sv2.z,
        )
        # Update lt2 from the retarded position exactly like the first
        # epoch: the aberration formula scales with light-time, so an
        # asymmetric (geometric vs retarded) lt between the two epochs
        # would leak a spurious velocity into the finite difference.
        dist2 = math.sqrt(geo2[0] ** 2 + geo2[1] ** 2 + geo2[2] ** 2)
        lt2 = dist2 / c_au_day

    if apply_deflection:
        # Use same deflector positions (good enough for dt=1s)
        geo2 = _apply_deflection_horizons(
            geo2,
            earth_sv2.pos,
            jd_tt2,
            lt2,
            batch,
            client,
        )

    if not (iflag & FLG_NOABERR) and not (iflag & FLG_TRUEPOS):
        from .fast_calc import _apply_aberration

        geo2 = _apply_aberration(geo2, earth_sv2.vel, lt2)

    # Velocity = (pos2 - pos1) / dt in ICRS
    geo_vel = (
        (geo2[0] - geo[0]) / dt,
        (geo2[1] - geo[1]) / dt,
        (geo2[2] - geo[2]) / dt,
    )

    return _to_ecliptic_output(geo, geo_vel, jd_tt, jd_ut, iflag, sid_mode)


class _HorizonsDeflectorSource:
    """Adapter exposing Horizons state vectors via LEBReader's eval_body API.

    Lets the Horizons pipeline reuse fast_calc._apply_gravitational_deflection
    (the verified NOVAS/Skyfield formula) instead of maintaining a second
    deflection implementation.  Deflector positions come from the prefetch
    batch when available, otherwise from the (cached) client.
    """

    def __init__(self, batch: dict, client: "HorizonsClient", jd_batch: float):
        self._batch = batch
        self._client = client
        self._jd_batch = jd_batch

    def eval_body(
        self, body_id: int, jd_tt: float
    ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        command = _DEFLECTOR_NAIF.get(body_id)
        if command is None:
            raise KeyError(f"body {body_id} is not a Horizons deflector")
        if jd_tt == self._jd_batch:
            sv = self._batch.get((command, jd_tt, "@0"))
            if sv is not None:
                return sv.pos, sv.vel
        sv = self._client.fetch_state_vector(command, jd_tt, "@0", "TDB")
        return sv.pos, sv.vel


def _apply_deflection_horizons(
    geo: Tuple[float, float, float],
    earth_bary: Tuple[float, float, float],
    jd_tt: float,
    light_time: float,
    batch: dict,
    client: "HorizonsClient",
) -> Tuple[float, float, float]:
    """Apply gravitational deflection using Horizons-fetched deflector data.

    Delegates the geometry to fast_calc._apply_gravitational_deflection so
    all three backends share one verified PPN formula.  Network failures for
    a deflector degrade gracefully (that deflector is skipped inside the
    shared routine via its KeyError/ValueError handling).
    """
    from .fast_calc import _apply_gravitational_deflection

    source = _HorizonsDeflectorSource(batch, client, jd_tt)
    try:
        return _apply_gravitational_deflection(
            geo,
            earth_bary,
            jd_tt,
            light_time,
            source,  # type: ignore[arg-type]
        )
    except ConnectionError:
        logger.warning(
            "Horizons deflector fetch failed; returning undeflected position"
        )
        return geo


def _to_ecliptic_output(
    pos_icrs: Tuple[float, float, float],
    vel_icrs: Tuple[float, float, float],
    jd_tt: float,
    jd_ut: float,
    iflag: int,
    sid_mode: int | None = None,
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """Convert ICRS Cartesian to ecliptic spherical output."""
    from .constants import (
        FLG_EQUATORIAL,
        FLG_J2000,
        FLG_SIDEREAL,
        FLG_SPEED,
        FLG_XYZ,
        FLG_RADIANS,
        FLG_ICRS,
    )
    from .fast_calc import (
        _VEL_H,
        _cartesian_to_spherical,
        _cartesian_velocity_to_spherical,
        _rotate_equatorial_to_ecliptic,
        _rotate_icrs_to_ecliptic_j2000,
        _get_skyfield_frame_data,
        _mat3_vec3,
        _prec_matrix,
        _spherical_to_cartesian_with_velocity,
        _general_precession_rate_deg_day,
    )
    import math

    def _rot(p: Tuple[float, float, float], jd: float) -> Tuple[float, float, float]:
        """Rotate an ICRS vector into the requested output frame at epoch jd."""
        if iflag & FLG_ICRS:
            # Output in ICRS — no rotation
            return p
        if iflag & FLG_J2000:
            if iflag & FLG_EQUATORIAL:
                # Equatorial J2000: ICRS vectors as-is (frame bias ~0.02",
                # same convention as the LEB pipeline)
                return p
            # J2000 ecliptic
            return _rotate_icrs_to_ecliptic_j2000(*p)
        if (iflag & FLG_EQUATORIAL) and (iflag & FLG_SIDEREAL):
            # Sidereal + equatorial: mean equator of date (precession only, no
            # nutation) and NO ayanamsa, mirroring the Skyfield/LEB path in
            # fast_calc (the ayanamsa is applied to ecliptic longitude only).
            return _mat3_vec3(_prec_matrix(jd), p)
        if iflag & FLG_EQUATORIAL:
            # True equatorial of date — apply precession-nutation matrix
            pn_mat, _, _, _ = _get_skyfield_frame_data(jd)
            return _mat3_vec3(pn_mat, p)
        # Default: ecliptic of date
        pn_mat, _, _, eps_true = _get_skyfield_frame_data(jd)
        q = _mat3_vec3(pn_mat, p)
        # Equatorial -> ecliptic (eps_true is in radians)
        return _rotate_equatorial_to_ecliptic(*q, eps_true)

    # ICRS / J2000 output frames are fixed in time; the of-date frames
    # (ecliptic/equator of date, mean equator for SID+EQ) rotate with the
    # equinox, so their reported speed must include the frame-rotation rate.
    _fixed_frame = bool(iflag & (FLG_ICRS | FLG_J2000))

    pos = _rot(pos_icrs, jd_tt)
    lon, lat, dist = _cartesian_to_spherical(*pos)

    if _fixed_frame or not (iflag & FLG_SPEED):
        # Fixed frame: rotating the ICRS velocity with the (constant) matrix
        # gives the exact derivative of the reported position — no correction.
        # (Same convention as fast_calc post-rework: NO general-precession
        # subtraction for J2000 spherical output; the former subtraction here
        # made the J2000 dlon differ from the derivative of the J2000
        # longitude by exactly the precession rate.)
        # Without FLG_SPEED the callers pass a zero velocity that must stay 0.
        vel = _rot(vel_icrs, jd_tt)
        dlon, dlat, ddist = _cartesian_velocity_to_spherical(*pos, *vel)
    else:
        # Of-date frame with FLG_SPEED: central-difference the reported place
        # through the frame at jd_tt ± _VEL_H, extrapolating the ICRS state
        # analytically — mirroring fast_calc._pipeline_icrs, so the reported
        # speed is the exact derivative of the reported position (captures the
        # precession/nutation equinox motion, ~0.09"/day in longitude).
        vel = _rot(vel_icrs, jd_tt)  # for the XYZ (Cartesian) output below
        _wm = _rot(
            (
                pos_icrs[0] - vel_icrs[0] * _VEL_H,
                pos_icrs[1] - vel_icrs[1] * _VEL_H,
                pos_icrs[2] - vel_icrs[2] * _VEL_H,
            ),
            jd_tt - _VEL_H,
        )
        _wp = _rot(
            (
                pos_icrs[0] + vel_icrs[0] * _VEL_H,
                pos_icrs[1] + vel_icrs[1] * _VEL_H,
                pos_icrs[2] + vel_icrs[2] * _VEL_H,
            ),
            jd_tt + _VEL_H,
        )
        _sm = _cartesian_to_spherical(*_wm)
        _sp = _cartesian_to_spherical(*_wp)
        _inv2h = 1.0 / (2.0 * _VEL_H)
        _dl = (_sp[0] - _sm[0]) % 360.0
        if _dl > 180.0:
            _dl -= 360.0
        dlon = _dl * _inv2h
        dlat = (_sp[1] - _sm[1]) * _inv2h
        ddist = (_sp[2] - _sm[2]) * _inv2h
        if iflag & FLG_XYZ:
            # Cartesian velocity in the of-date frame: same central difference
            # on the Cartesian components (fast_calc's XYZ convention).
            vel = tuple((_wp[i] - _wm[i]) * _inv2h for i in range(3))  # type: ignore[assignment]

    # Sidereal correction (ecliptic longitude only; equatorial sidereal output
    # uses the mean-equator-of-date frame handled above, with no ayanamsa — the
    # same rule as planets._calc_body_pctr and fast_calc).
    #
    # Only Pipeline-A-like bodies (planets/asteroids/helio/bary state vectors)
    # reach _to_ecliptic_output; the deferred ecliptic-direct bodies
    # (nodes/apogees) are handled by _calc_analytical and never arrive here, so
    # the simple J2000 handling below — without fast_calc's _deferred_sid_j2k
    # rebuild — is correct for every body that does.
    if (iflag & FLG_SIDEREAL) and not (iflag & FLG_EQUATORIAL):
        from .planets import _get_ayanamsa_for_flags

        # The longitude is on the TRUE ecliptic of date (the default branch
        # above applies the precession-nutation matrix), so it already carries
        # the nutation-in-longitude term Δψ. Subtract the TRUE ayanamsa
        # (mean + Δψ), exactly as fast_calc (lines 1996-2003) and
        # planets._apply_sidereal_correction do. Using the mean ayanamsa here
        # (the old get_ayanamsa_ut call) left a spurious ~9-17" offset versus
        # the LEB/Skyfield backend for the same request. _get_ayanamsa_for_flags
        # correctly returns the mean ayanamsa for FLG_NONUT/FLG_J2000 (where the
        # position carries no nutation), so the J2000 ecliptic branch above
        # stays correct.
        ayan = _get_ayanamsa_for_flags(jd_ut, iflag, sid_mode)
        lon = (lon - ayan) % 360.0

        # Sidereal speed correction (ayanamsa drift): the sidereal longitude
        # subtracts the ayanamsa, so its rate subtracts the ayanamsa drift
        # (mirroring fast_calc post-rework). Done here so BOTH the spherical
        # return below and the XYZ rebuild use the sidereal rate. Fires for
        # every sidereal SPEED request:
        #   - of-date: the position subtracted the TRUE ayanamsa (mean + Δψ),
        #     so the rate drops the general-precession rate PLUS the
        #     nutation-in-longitude rate dΔψ/dt (~0.05"/day);
        #   - J2000: the mean ayanamsa still drifts in the fixed frame, so the
        #     rate drops the general-precession rate only (no nutation, and no
        #     second frame term — the former "J2000 frame conversion" block
        #     below double-counted it).
        # Skipped only when FLG_SPEED is absent (dlon is 0.0 and must stay 0.0).
        if iflag & FLG_SPEED:
            from .fast_calc import _STAR_BASED_MODES

            if sid_mode in _STAR_BASED_MODES:
                # Star-anchored modes: the ayanamsa drift is the anchor star's
                # apparent-longitude drift (annual-aberration-dominated), not the
                # IAU 2006 precession polynomial. Central-difference the actual
                # ayanamsa (_get_ayanamsa_for_flags returns mean for FLG_J2000/
                # FLG_NONUT and true otherwise), matching the Skyfield and LEB
                # backends (see fast_calc's sidereal speed correction).
                # Unwrap the delta into (-180, 180]: the ayanamsha is mod 360,
                # and the two samples can straddle the 0/360 wrap (e.g.
                # SIDM_GALCENT_COCHRANE near JD 2533810), which would turn the
                # raw -360 difference into a ~1.6e7 deg/day speed spike.
                _dt_aya = 1.0 / 86400.0
                _aya_p = _get_ayanamsa_for_flags(jd_ut + _dt_aya, iflag, sid_mode)
                _aya_m = _get_ayanamsa_for_flags(jd_ut - _dt_aya, iflag, sid_mode)
                _d_aya = (_aya_p - _aya_m + 180.0) % 360.0 - 180.0
                dlon -= _d_aya / (2.0 * _dt_aya)
            else:
                dlon -= _general_precession_rate_deg_day(jd_tt)
                if not (iflag & FLG_J2000):
                    _, _dpsi_m, _, _ = _get_skyfield_frame_data(jd_tt - _VEL_H)
                    _, _dpsi_p, _, _ = _get_skyfield_frame_data(jd_tt + _VEL_H)
                    dlon -= math.degrees(_dpsi_p - _dpsi_m) / (2.0 * _VEL_H)

    # NOTE: J2000 spherical output needs NO precession-rate subtraction here.
    # The J2000 frame is fixed, so the fixed-matrix rotation of the ICRS
    # velocity already yields the exact derivative of the reported J2000
    # longitude (same convention as fast_calc post-rework). The former
    # subtraction here made the J2000 dlon differ from that derivative by
    # exactly the general-precession rate, and double-counted the term for
    # SID|J2000 on top of the ayanamsa drift above.

    # XYZ output
    if iflag & FLG_XYZ:
        if (iflag & FLG_SIDEREAL) and not (iflag & FLG_EQUATORIAL):
            # Sidereal ecliptic Cartesian: rebuild the vectors from the
            # ayanamsa-shifted longitude, mirroring fast_calc's XYZ+SIDEREAL
            # post-processing (shared helper). The bare pos/vel were rotated to
            # ecliptic-of-date but NOT by the ayanamsa, so returning them as-is
            # would yield tropical Cartesian output while the LEB/Skyfield
            # backends return sidereal (ayanamsa-rotated) Cartesian.
            return (
                _spherical_to_cartesian_with_velocity(
                    lon, lat, dist, dlon, dlat, ddist
                ),
                iflag,
            )
        return (pos + vel, iflag)  # type: ignore

    # Radians
    if iflag & FLG_RADIANS:
        lon = math.radians(lon)
        lat = math.radians(lat)
        dlon = math.radians(dlon)
        dlat = math.radians(dlat)

    return ((lon, lat, dist, dlon, dlat, ddist), iflag)


def _calc_analytical(
    jd_ut: float, body_id: int, iflag: int, sid_mode: int | None = None
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """Calculate analytical body (Mean Node, Mean Apogee/Lilith).

    Mirrors the canonical Mean Node / Mean Apogee branches of
    ``planets._calc_body`` so the Horizons backend produces identical results
    to the LEB/Skyfield path for every flag combination. The previous version
    returned a bare mean-ecliptic-of-date spherical 6-tuple and ignored
    FLG_J2000 / FLG_EQUATORIAL / FLG_XYZ / FLG_RADIANS, omitted the nutation-in-
    longitude term Δψ that the of-date output carries (a ~17" tropical error),
    and left the sidereal longitude speed uncorrected for the ayanamsha drift.
    All of this is pure analytical math (no DE440 / HTTP), so it stays valid in
    Horizons mode.
    """
    import math

    from .time_utils import deltat
    from .constants import FLG_SIDEREAL, FLG_NONUT, FLG_EQUATORIAL, FLG_SPEED
    from .constants import _MOON_MEAN_DIST_AU, _MOON_MEAN_APOG_DIST_AU
    from .planets import (
        _apply_output_flags,
        _apply_sidereal_correction,
        _maybe_equatorial_convert,
        _nutation_rate_deg_per_day,
    )
    from .cache import get_cached_nutation

    jd_tt = jd_ut + deltat(jd_ut)
    is_sidereal = bool(iflag & FLG_SIDEREAL)
    # SIDEREAL+EQUATORIAL output is on the mean ecliptic (no nutation), so Δψ is
    # neither added to the position nor its rate to the speed.
    _sid_eq = is_sidereal and bool(iflag & FLG_EQUATORIAL)
    add_nutation = not (iflag & FLG_NONUT) and not _sid_eq

    if body_id == 10:  # Mean Node
        from .lunar import calc_mean_lunar_node

        lon = calc_mean_lunar_node(jd_tt)
        lat = 0.0
        dist = _MOON_MEAN_DIST_AU
    elif body_id == 12:  # Mean Apogee (Lilith)
        from .lunar import calc_mean_lilith_with_latitude

        lon, lat = calc_mean_lilith_with_latitude(jd_tt)
        dist = _MOON_MEAN_APOG_DIST_AU
    else:
        raise KeyError(f"Body {body_id} not analytical")

    # The analytical formulae return the MEAN ecliptic of date; the reference
    # outputs the TRUE ecliptic, so add Δψ to the longitude unless suppressed.
    if add_nutation:
        dpsi_rad, _ = get_cached_nutation(jd_tt)
        lon = (lon + math.degrees(dpsi_rad)) % 360.0

    # Speed via ±0.5-day central difference (matching the canonical path).
    dlon, dlat = 0.0, 0.0
    if iflag & FLG_SPEED:
        dt = 0.5
        if body_id == 10:
            lon_prev = calc_mean_lunar_node(jd_tt - dt)
            lon_next = calc_mean_lunar_node(jd_tt + dt)
        else:
            lon_prev, lat_prev = calc_mean_lilith_with_latitude(jd_tt - dt)
            lon_next, lat_next = calc_mean_lilith_with_latitude(jd_tt + dt)
            dlat = (lat_next - lat_prev) / (2.0 * dt)
        lon_diff = lon_next - lon_prev
        if lon_diff > 180:
            lon_diff -= 360.0
        elif lon_diff < -180:
            lon_diff += 360.0
        dlon = lon_diff / (2.0 * dt)
        # True-ecliptic output: the speed carries the nutation rate, exactly
        # as the position carries Δψ.
        if add_nutation:
            dlon += _nutation_rate_deg_per_day(jd_tt, dt)

    # Sidereal correction (flag-aware ayanamsa + ayanamsha-rate speed term),
    # ecliptic longitude only — equatorial sidereal output uses the mean-equator
    # frame in _maybe_equatorial_convert with no ayanamsa.
    if is_sidereal and not (iflag & FLG_EQUATORIAL):
        lon, dlon = _apply_sidereal_correction(lon, dlon, jd_ut, iflag, sid_mode)

    result = (lon, lat, dist, dlon, dlat, 0.0)
    # FLG_J2000 (precession) and FLG_EQUATORIAL conversion, then the output
    # representation flags (FLG_XYZ / FLG_RADIANS).
    result = _maybe_equatorial_convert(result, jd_tt, iflag)
    result = _apply_output_flags(result, iflag)
    return (result, iflag)


def _calc_uranian(
    jd_ut: float, body_id: int, iflag: int, sid_mode: int | None = None
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """Calculate Uranian hypothetical body (heliocentric only)."""
    from .time_utils import deltat
    from .constants import FLG_SIDEREAL, FLG_J2000, FLG_EQUATORIAL, FLG_SPEED

    jd_tt = jd_ut + deltat(jd_ut)

    from .hypothetical import calc_uranian_planet, calc_transpluto

    if body_id == 48:
        lon, lat, dist, dlon, dlat, ddist = calc_transpluto(jd_tt)
    else:
        lon, lat, dist, dlon, dlat, ddist = calc_uranian_planet(body_id, jd_tt)

    # Zero the speed slots when FLG_SPEED is absent so callers never receive
    # unrequested velocity data — mirroring the state-vector and analytical
    # paths. (The sidereal speed correction below is itself gated on FLG_SPEED,
    # so it leaves these zeros untouched.)
    if not (iflag & FLG_SPEED):
        dlon = dlat = ddist = 0.0

    # calc_uranian_planet / calc_transpluto return a J2000-frame ecliptic
    # longitude. Precess it to the ecliptic OF DATE unless FLG_J2000 is set
    # (carrying the velocity when FLG_SPEED so the of-date speed is the true
    # derivative of the of-date position), add the of-date nutation Δψ, then
    # apply the sidereal correction — mirroring the canonical heliocentric
    # Uranian path in planets._calc_body so both backends agree.
    if not (iflag & FLG_J2000):
        if iflag & FLG_SPEED:
            from .planets import _precess_ecliptic_state

            lon, lat, dlon, dlat = _precess_ecliptic_state(
                lon, lat, dlon, dlat, jd_tt, to_j2000=False
            )
        else:
            from .astrometry import _precess_ecliptic

            lon, lat = _precess_ecliptic(lon, lat, 2451545.0, jd_tt)

    from .planets import _add_of_date_nutation

    lon, dlon = _add_of_date_nutation(lon, dlon, jd_tt, iflag)

    if (iflag & FLG_SIDEREAL) and not (iflag & FLG_EQUATORIAL):
        from .planets import _apply_sidereal_correction

        lon, dlon = _apply_sidereal_correction(lon, dlon, jd_ut, iflag, sid_mode)

    # Equatorial conversion and output-representation flags, mirroring the
    # canonical heliocentric Uranian path in planets._calc_body. The bare
    # return above produced ecliptic spherical DEGREES for every request, so
    # FLG_EQUATORIAL / FLG_XYZ / FLG_RADIANS were silently ignored — a
    # tens-of-degrees / wrong-representation divergence from the LEB/Skyfield
    # backend. already_j2000=True because the J2000 framing was handled by the
    # precession gate above: the converter must skip its internal date→J2000
    # precession but still pick the J2000 obliquity for EQ+J2000 output (the
    # canonical path passes the same combination). Both helpers are pure
    # analytical math (obliquity/nutation via erfa), so no DE440/HTTP is needed.
    from .planets import _apply_output_flags, _maybe_equatorial_convert

    result = (lon, lat, dist, dlon, dlat, ddist)
    result = _maybe_equatorial_convert(result, jd_tt, iflag, already_j2000=True)
    result = _apply_output_flags(result, iflag)
    return (result, iflag)
