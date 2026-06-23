# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
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
from typing import Dict, List, Tuple

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
    ) -> Dict[Tuple[str, float, str], StateVector]:
        """Fetch multiple state vectors in parallel.

        Args:
            requests: List of (command, jd, center) tuples.

        Returns:
            Dict mapping (command, jd, center) -> StateVector.
        """
        # Filter out cached results
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
                except (OSError, ValueError, KeyError):
                    pass  # skip failed fetches

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
        # CSV format: JDTDB, Calendar Date, X, Y, Z, VX, VY, VZ, LT, RG, RR,
        lines = [line.strip() for line in data_block.split("\n") if line.strip()]
        if not lines:
            raise ValueError(f"Empty data block in Horizons response for {command}")

        # Parse the first (and only) data line
        # CSV fields are comma-separated
        values = []
        for line in lines:
            for field in line.split(","):
                field = field.strip()
                if field:
                    try:
                        values.append(float(field))
                    except ValueError:
                        pass  # skip non-numeric (calendar date string)

        # We need at least: JDTDB, X, Y, Z, VX, VY, VZ (7 values)
        if len(values) < 7:
            raise ValueError(
                f"Insufficient data in Horizons response for {command}: "
                f"got {len(values)} values, need at least 7"
            )

        # Values: [JDTDB, X, Y, Z, VX, VY, VZ, ...]
        return StateVector(
            x=values[1],
            y=values[2],
            z=values[3],
            vx=values[4],
            vy=values[5],
            vz=values[6],
        )

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
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """Calculate body position via Horizons API.

    Same return format as fast_calc.fast_calc_ut().

    Args:
        client: HorizonsClient instance.
        jd_ut: Julian Day UT.
        body_id: SE_* body constant.
        iflag: Swiss Ephemeris flags.

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
    if iflag & FLG_NONUT:
        # The of-date frames below always include nutation; defer to Skyfield
        # (same policy as the LEB backend).
        raise KeyError("FLG_NONUT not supported by Horizons backend")

    # Analytical bodies — no HTTP needed
    if body_id in _ANALYTICAL_BODIES:
        return _calc_analytical(jd_ut, body_id, iflag)

    # Uranian hypotheticals — analytical, heliocentric only
    if body_id in _URANIAN_BODIES:
        if not (iflag & FLG_HELCTR):
            raise KeyError(
                f"Uranian body {body_id} geocentric not supported via Horizons"
            )
        return _calc_uranian(jd_ut, body_id, iflag)

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
                    tgt_lt = client.fetch_state_vector(
                        command, jd_tt - lt, "@0", "TDB"
                    )
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
            return _to_ecliptic_output(rel_pos, rel_vel, jd_tt, jd_ut, iflag)

        # Barycentric: the SSB is inertial, so retardation is just an
        # earlier evaluation epoch.
        sv = client.fetch_state_vector(command, jd_tt, "@0", "TDB")
        if not (iflag & FLG_TRUEPOS):
            dist = _math.sqrt(sv.x**2 + sv.y**2 + sv.z**2)
            if dist > 0.0:
                lt = dist / c_au_day
                sv = client.fetch_state_vector(command, jd_tt - lt, "@0", "TDB")
        return _to_ecliptic_output(sv.pos, sv.vel, jd_tt, jd_ut, iflag)

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
    batch = client.fetch_batch(prefetch_cmds)

    target_sv = batch.get((command, jd_tt, "@0"))
    earth_sv = batch.get(("399", jd_tt, "@0"))

    if target_sv is None or earth_sv is None:
        raise ConnectionError(f"Failed to fetch target/Earth for body {body_id}")

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
    # fast_calc and Swiss Ephemeris app_pos_etc_plan, which gate deflection on
    # both NOGDEFL and TRUEPOS — like the light-time and aberration steps).
    if apply_deflection:
        geo = _apply_deflection_horizons(geo, earth_sv.pos, jd_tt, lt, batch, client)

    # Aberration (relativistic when light-time is known, matching fast_calc)
    if not (iflag & FLG_NOABERR) and not (iflag & FLG_TRUEPOS):
        from .fast_calc import _apply_aberration

        geo = _apply_aberration(geo, earth_sv.vel, lt)

    # Speeds only on request — they cost several extra HTTP round-trips
    if not (iflag & FLG_SPEED):
        return _to_ecliptic_output(geo, (0.0, 0.0, 0.0), jd_tt, jd_ut, iflag)

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

    return _to_ecliptic_output(geo, geo_vel, jd_tt, jd_ut, iflag)


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
            geo, earth_bary, jd_tt, light_time, source  # type: ignore[arg-type]
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
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """Convert ICRS Cartesian to ecliptic spherical output."""
    from .constants import (
        FLG_EQUATORIAL,
        FLG_J2000,
        FLG_SIDEREAL,
        FLG_XYZ,
        FLG_RADIANS,
        FLG_ICRS,
    )
    from .fast_calc import (
        _cartesian_to_spherical,
        _cartesian_velocity_to_spherical,
        _rotate_equatorial_to_ecliptic,
        _rotate_icrs_to_ecliptic_j2000,
        _get_skyfield_frame_data,
        _mat3_vec3,
        _prec_matrix,
    )
    import math

    pos = pos_icrs
    vel = vel_icrs

    if iflag & FLG_ICRS:
        # Output in ICRS — no rotation
        pass
    elif iflag & FLG_J2000:
        if iflag & FLG_EQUATORIAL:
            # Equatorial J2000: ICRS vectors as-is (frame bias ~0.02",
            # same convention as the LEB pipeline)
            pass
        else:
            # J2000 ecliptic
            pos = _rotate_icrs_to_ecliptic_j2000(*pos)
            vel = _rotate_icrs_to_ecliptic_j2000(*vel)
    elif (iflag & FLG_EQUATORIAL) and (iflag & FLG_SIDEREAL):
        # Sidereal + equatorial: mean equator of date (precession only, no
        # nutation) and NO ayanamsa, mirroring the Skyfield/LEB path in
        # fast_calc (the ayanamsa is applied to ecliptic longitude only).
        p_mat = _prec_matrix(jd_tt)
        pos = _mat3_vec3(p_mat, pos)
        vel = _mat3_vec3(p_mat, vel)
    elif iflag & FLG_EQUATORIAL:
        # True equatorial of date — apply precession-nutation matrix
        pn_mat, dpsi, deps, eps_true = _get_skyfield_frame_data(jd_tt)
        pos = _mat3_vec3(pn_mat, pos)
        vel = _mat3_vec3(pn_mat, vel)
    else:
        # Default: ecliptic of date
        pn_mat, dpsi, deps, eps_true = _get_skyfield_frame_data(jd_tt)
        # ICRS -> equatorial of date
        pos = _mat3_vec3(pn_mat, pos)
        vel = _mat3_vec3(pn_mat, vel)
        # Equatorial -> ecliptic (eps_true is in radians)
        pos = _rotate_equatorial_to_ecliptic(*pos, eps_true)
        vel = _rotate_equatorial_to_ecliptic(*vel, eps_true)

    # Convert to spherical
    lon, lat, dist = _cartesian_to_spherical(*pos)
    dlon, dlat, ddist = _cartesian_velocity_to_spherical(*pos, *vel)

    # Sidereal correction (ecliptic longitude only; equatorial sidereal output
    # uses the mean-equator-of-date frame handled above, with no ayanamsa — the
    # same rule as planets._calc_body_pctr and fast_calc).
    if (iflag & FLG_SIDEREAL) and not (iflag & FLG_EQUATORIAL):
        from .planets import get_ayanamsa_ut

        ayan = get_ayanamsa_ut(jd_ut)
        lon = (lon - ayan) % 360.0

    # XYZ output
    if iflag & FLG_XYZ:
        if (iflag & FLG_SIDEREAL) and not (iflag & FLG_EQUATORIAL):
            # Sidereal ecliptic Cartesian: rebuild the vectors from the
            # ayanamsa-shifted longitude, mirroring fast_calc's XYZ+SIDEREAL
            # post-processing. The bare pos/vel were rotated to ecliptic-of-date
            # but NOT by the ayanamsa, so returning them as-is would yield
            # tropical Cartesian output while the LEB/Skyfield backends return
            # sidereal (ayanamsa-rotated) Cartesian for the same request.
            lon_r = math.radians(lon)
            lat_r = math.radians(lat)
            cos_lat = math.cos(lat_r)
            sin_lat = math.sin(lat_r)
            cos_lon = math.cos(lon_r)
            sin_lon = math.sin(lon_r)
            x = dist * cos_lat * cos_lon
            y = dist * cos_lat * sin_lon
            z = dist * sin_lat
            if dlon != 0.0 or dlat != 0.0 or ddist != 0.0:
                dlon_r = math.radians(dlon)
                dlat_r = math.radians(dlat)
                vx = (
                    -dist * cos_lat * sin_lon * dlon_r
                    - dist * sin_lat * cos_lon * dlat_r
                    + cos_lat * cos_lon * ddist
                )
                vy = (
                    dist * cos_lat * cos_lon * dlon_r
                    - dist * sin_lat * sin_lon * dlat_r
                    + cos_lat * sin_lon * ddist
                )
                vz = dist * cos_lat * dlat_r + sin_lat * ddist
            else:
                vx = vy = vz = 0.0
            return ((x, y, z, vx, vy, vz), iflag)
        return (pos + vel, iflag)  # type: ignore

    # Radians
    if iflag & FLG_RADIANS:
        lon = math.radians(lon)
        lat = math.radians(lat)
        dlon = math.radians(dlon)
        dlat = math.radians(dlat)

    return ((lon, lat, dist, dlon, dlat, ddist), iflag)


def _calc_analytical(
    jd_ut: float, body_id: int, iflag: int
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """Calculate analytical body (Mean Node, Mean Apogee)."""
    from .time_utils import deltat
    from .constants import FLG_SIDEREAL

    jd_tt = jd_ut + deltat(jd_ut)

    if body_id == 10:  # Mean Node
        from .lunar import calc_mean_lunar_node

        lon = calc_mean_lunar_node(jd_tt)
        lat = 0.0
        dist = 0.002569  # mean lunar distance in AU (approximate)
    elif body_id == 12:  # Mean Apogee (Lilith)
        from .lunar import calc_mean_lilith_with_latitude

        lon, lat = calc_mean_lilith_with_latitude(jd_tt)
        dist = 0.002710  # mean apogee distance
    else:
        raise KeyError(f"Body {body_id} not analytical")

    # Speed via finite difference
    dt = 1.0 / 86400.0  # 1 second
    dlat = 0.0
    if body_id == 10:
        from .lunar import calc_mean_lunar_node

        lon2 = calc_mean_lunar_node(jd_tt + dt)
        dlon = (lon2 - lon) / dt
        if abs(dlon) > 180 / dt:
            dlon = ((lon2 - lon + 180) % 360 - 180) / dt
    elif body_id == 12:
        from .lunar import calc_mean_lilith_with_latitude

        lon2, lat2 = calc_mean_lilith_with_latitude(jd_tt + dt)
        dlon = ((lon2 - lon + 180) % 360 - 180) / dt
        dlat = (lat2 - lat) / dt
    else:  # pragma: no cover - only bodies 10/12 reach here (others raise KeyError above)
        dlon = 0.0

    # Sidereal
    if iflag & FLG_SIDEREAL:
        from .planets import get_ayanamsa_ut

        ayan = get_ayanamsa_ut(jd_ut)
        lon = (lon - ayan) % 360.0

    return ((lon, lat, dist, dlon, dlat, 0.0), iflag)


def _calc_uranian(
    jd_ut: float, body_id: int, iflag: int
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """Calculate Uranian hypothetical body (heliocentric only)."""
    from .time_utils import deltat
    from .constants import FLG_SIDEREAL

    jd_tt = jd_ut + deltat(jd_ut)

    from .hypothetical import calc_uranian_planet, calc_transpluto

    if body_id == 48:
        lon, lat, dist, dlon, dlat, ddist = calc_transpluto(jd_tt)
    else:
        lon, lat, dist, dlon, dlat, ddist = calc_uranian_planet(body_id, jd_tt)

    if iflag & FLG_SIDEREAL:
        from .planets import get_ayanamsa_ut

        ayan = get_ayanamsa_ut(jd_ut)
        lon = (lon - ayan) % 360.0

    return ((lon, lat, dist, dlon, dlat, ddist), iflag)
