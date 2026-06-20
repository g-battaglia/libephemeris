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
        cache_key = (round(jd, 12), command, center)

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
            cache_key = (round(jd, 12), cmd, center)
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
                    # Do not raise: a failed deflector (Sun/Jupiter/Saturn) must
                    # not abort the whole calculation. But do not fail silently
                    # either — a dropped fetch degrades precision (e.g. missing
                    # gravitational deflection), so make it visible in the log.
                    cmd, jd, center = key
                    logger.warning(
                        "Horizons fetch failed for command=%s jd=%s center=%s: %s",
                        cmd,
                        jd,
                        center,
                        exc,
                    )

        return results

    def _fetch_with_retry(
        self, url: str, command: str, max_retries: int = 2
    ) -> StateVector:
        """Fetch URL with exponential backoff retry."""
        import time

        last_err = None
        for attempt in range(max_retries + 1):
            try:
                req = urllib.request.Request(
                    url,
                    headers={"User-Agent": "libephemeris/1.0"},
                )
                with urllib.request.urlopen(req, timeout=self._timeout) as resp:
                    data = json.loads(resp.read().decode("utf-8"))

                if "error" in data:
                    # Permanent API error (e.g. body not found). Signalled as
                    # KeyError per this method's documented contract; not retried.
                    raise KeyError(f"Horizons API error for {command}: {data['error']}")

                return self._parse_response(data, command)

            except KeyError:
                raise  # permanent (body not found / API error): do not retry
            except (OSError, ValueError) as e:
                # Transient: network errors and malformed/partial responses.
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
        FLG_EQUATORIAL,
        FLG_HELCTR,
        FLG_BARYCTR,
        FLG_J2000,
        FLG_NOABERR,
        FLG_NOGDEFL,
        FLG_SIDEREAL,
        FLG_SPEED,
        FLG_TOPOCTR,
        FLG_TRUEPOS,
        FLG_XYZ,
        FLG_RADIANS,
        FLG_ICRS,
        FLG_NONUT,
    )

    # Unsupported flags → fallback to Skyfield
    if iflag & FLG_TOPOCTR:
        raise KeyError("FLG_TOPOCTR not supported by Horizons backend")

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

    # Heliocentric or barycentric — simpler pipeline
    if iflag & FLG_HELCTR:
        if body_id == 0:
            # Sun heliocentric = Sun from Sun = (0, 0, 0)
            return ((0.0, 0.0, 0.0, 0.0, 0.0, 0.0), iflag)
        sv = client.fetch_state_vector(command, jd_tt, center="@10", time_type="TDB")
        return _to_ecliptic_output(sv.pos, sv.vel, jd_tt, iflag)

    if iflag & FLG_BARYCTR:
        sv = client.fetch_state_vector(command, jd_tt, center="@0", time_type="TDB")
        return _to_ecliptic_output(sv.pos, sv.vel, jd_tt, iflag)

    # Geocentric apparent — full pipeline
    # Prefetch: target + Earth + deflectors (Sun, Jupiter, Saturn)
    prefetch_cmds = [
        (command, jd_tt, "@0"),  # target barycentric
        ("399", jd_tt, "@0"),  # Earth barycentric
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

    # Gravitational deflection
    if not (iflag & FLG_NOGDEFL):
        geo = _apply_deflection_horizons(geo, earth_sv.pos, jd_tt, lt, batch)

    # Aberration
    if not (iflag & FLG_NOABERR) and not (iflag & FLG_TRUEPOS):
        from .fast_calc import _apply_aberration

        geo = _apply_aberration(geo, earth_sv.vel)

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
    if not (iflag & FLG_TRUEPOS):
        dist2 = math.sqrt(geo2[0] ** 2 + geo2[1] ** 2 + geo2[2] ** 2)
        lt2 = dist2 / c_au_day
        target_lt2 = client.fetch_state_vector(command, jd_tt2 - lt2, "@0", "TDB")
        geo2 = (
            target_lt2.x - earth_sv2.x,
            target_lt2.y - earth_sv2.y,
            target_lt2.z - earth_sv2.z,
        )

    if not (iflag & FLG_NOGDEFL):
        # Use same deflector positions (good enough for dt=1s)
        geo2 = _apply_deflection_horizons(
            geo2,
            earth_sv2.pos,
            jd_tt2,
            lt if not (iflag & FLG_TRUEPOS) else 0.0,
            batch,
        )

    if not (iflag & FLG_NOABERR) and not (iflag & FLG_TRUEPOS):
        from .fast_calc import _apply_aberration

        geo2 = _apply_aberration(geo2, earth_sv2.vel)

    # Velocity = (pos2 - pos1) / dt in ICRS
    geo_vel = (
        (geo2[0] - geo[0]) / dt,
        (geo2[1] - geo[1]) / dt,
        (geo2[2] - geo[2]) / dt,
    )

    return _to_ecliptic_output(geo, geo_vel, jd_tt, iflag)


def _apply_deflection_horizons(
    geo: Tuple[float, float, float],
    earth_bary: Tuple[float, float, float],
    jd_tt: float,
    light_time: float,
    batch: dict,
) -> Tuple[float, float, float]:
    """Apply gravitational deflection using Horizons-fetched deflector positions."""
    from .fast_calc import _vec3_sub, _vec3_dist, _mat3_vec3
    import math

    c_au_day = 173.14463267
    deflectors = [
        ("10", 1.32712440041279419e11),  # Sun GM
        ("5", 1.26712764945480000e8),  # Jupiter GM
        ("6", 3.79406260288322009e7),  # Saturn GM
    ]

    result = list(geo)

    for defl_cmd, gm in deflectors:
        key = (defl_cmd, jd_tt, "@0")
        if key not in batch:
            # The Sun ("10") dominates gravitational deflection (~1.8" near the
            # limb); skipping it silently would degrade precision invisibly.
            if defl_cmd == "10":
                logger.warning(
                    "Horizons: Sun deflector unavailable at jd=%s; "
                    "gravitational deflection is degraded for this position.",
                    jd_tt,
                )
            continue

        defl_sv = batch[key]
        defl_pos = defl_sv.pos

        # Deflector relative to Earth
        e = (
            defl_pos[0] - earth_bary[0],
            defl_pos[1] - earth_bary[1],
            defl_pos[2] - earth_bary[2],
        )
        e_dist = math.sqrt(e[0] ** 2 + e[1] ** 2 + e[2] ** 2)

        # Body relative to deflector
        q = (
            result[0] - e[0],
            result[1] - e[1],
            result[2] - e[2],
        )
        q_dist = math.sqrt(q[0] ** 2 + q[1] ** 2 + q[2] ** 2)

        geo_dist = math.sqrt(result[0] ** 2 + result[1] ** 2 + result[2] ** 2)

        if e_dist < 1e-20 or q_dist < 1e-20 or geo_dist < 1e-20:
            continue

        # PPN deflection angle
        # δθ ≈ (1+γ) GM / (c² e_dist) * (unit_geo + unit_e) / (1 + cos(angle))
        two_gm_c2 = 2.0 * gm / (c_au_day * c_au_day * 1.495978707e8)  # km -> AU

        dot_eq = (result[0] * e[0] + result[1] * e[1] + result[2] * e[2]) / (
            geo_dist * e_dist
        )

        if dot_eq > 0.9999:
            continue  # body behind deflector, skip

        factor = two_gm_c2 / (e_dist * (1.0 + dot_eq + 1e-30))

        for i in range(3):
            unit_geo_i = result[i] / geo_dist
            unit_e_i = e[i] / e_dist
            result[i] += factor * (unit_geo_i - dot_eq * unit_e_i)

    return tuple(result)  # type: ignore[return-value]


def _to_ecliptic_output(
    pos_icrs: Tuple[float, float, float],
    vel_icrs: Tuple[float, float, float],
    jd_tt: float,
    iflag: int,
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
        FLG_NONUT,
    )
    from .fast_calc import (
        _cartesian_to_spherical,
        _cartesian_velocity_to_spherical,
        _rotate_equatorial_to_ecliptic,
        _rotate_icrs_to_ecliptic_j2000,
        _get_skyfield_frame_data,
        _mat3_vec3,
        _mean_obliquity_iau2006,
    )
    import math

    pos = pos_icrs
    vel = vel_icrs

    if iflag & FLG_ICRS:
        # Output in ICRS — no rotation
        pass
    elif iflag & FLG_J2000:
        # J2000 ecliptic
        pos = _rotate_icrs_to_ecliptic_j2000(pos)
        vel = _rotate_icrs_to_ecliptic_j2000(vel)
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
        # Equatorial -> ecliptic
        pos = _rotate_equatorial_to_ecliptic(pos, eps_true)
        vel = _rotate_equatorial_to_ecliptic(vel, eps_true)

    # Convert to spherical
    lon, lat, dist = _cartesian_to_spherical(pos)
    dlon, dlat, ddist = _cartesian_velocity_to_spherical(pos, vel)

    # Sidereal correction
    if iflag & FLG_SIDEREAL:
        from .ayanamsha import get_ayanamsha_ut

        ayan = get_ayanamsha_ut(jd_tt)
        lon = (lon - ayan) % 360.0

    # XYZ output
    if iflag & FLG_XYZ:
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
    from .constants import FLG_SPEED, FLG_SIDEREAL

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
    else:
        dlon = 0.0

    # Sidereal
    if iflag & FLG_SIDEREAL:
        from .ayanamsha import get_ayanamsha_ut

        ayan = get_ayanamsha_ut(jd_tt)
        lon = (lon - ayan) % 360.0

    return ((lon, lat, dist, dlon, 0.0, 0.0), iflag)


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
        from .ayanamsha import get_ayanamsha_ut

        ayan = get_ayanamsha_ut(jd_tt)
        lon = (lon - ayan) % 360.0

    return ((lon, lat, dist, dlon, dlat, ddist), iflag)
