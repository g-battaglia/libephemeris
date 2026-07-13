# Proposal: Horizons Live Backend (Zero-Install Ephemeris)

**Status:** Draft — pending approval  
**Author:** AI pair-programming session  
**Date:** 2026-03-26  
**Priority:** Medium  
**Estimated effort:** 8-10 days  

---

## 1. Motivation

When a casual user does:

```bash
pip install libephemeris
python -c "import libephemeris as swe; print(swe.swe_calc_ut(2451545.0, swe.SE_MARS, swe.SEFLG_SPEED))"
```

…the library currently needs to download DE440 (~114 MB) before the first calculation can run. This is a friction point for users who want to quickly evaluate the library, run a one-off chart, or use it in lightweight/serverless environments.

**Goal:** Make `swe_calc_ut()` work immediately after `pip install` with zero local files, by transparently fetching ephemeris data from the NASA JPL Horizons REST API over HTTP.

---

## 2. Context: Current Architecture

### 2.1 Existing Backends

LibEphemeris currently has two computation backends, selected via `set_calc_mode()` / `LIBEPHEMERIS_MODE`:

| Backend | Data Source | Latency | Precision | File Size |
|---------|-----------|---------|-----------|-----------|
| **Skyfield** (default) | JPL DE440/441 `.bsp` files | ~0.5 ms/call | ~0.001" | 31 MB – 3.1 GB |
| **LEB** (fast path) | Precomputed `.leb` Chebyshev files | ~5 μs/call | ~0.001" | 5–50 MB |

Both require local files. Mode `"auto"` uses LEB if configured, otherwise Skyfield.

### 2.2 Calculation Pipeline (Skyfield path)

The call chain for `swe_calc_ut(jd, body, flags)` is:

```
swe_calc_ut()                          # planets.py:760
  ├── LEB fast path                    # planets.py:829-848
  │   └── fast_calc.fast_calc_ut()     # fast_calc.py:1270
  ├── [PROPOSED: Horizons path]        # <-- insertion point
  └── Skyfield path                    # planets.py:850+
      ├── get_cached_time_ut1(jd)      # cache.py
      └── _calc_body(t, body, flags)   # planets.py:1534
          ├── Dispatcher by body type
          │   ├── Standard planets → Skyfield .at(t).observe().apparent()
          │   ├── Lunar nodes → lunar.py (analytical or Skyfield Moon vectors)
          │   ├── Lilith variants → lunar.py (analytical or Skyfield Moon vectors)
          │   ├── Minor bodies → SPK / ASSIST / Keplerian
          │   ├── Fixed stars → fixed_stars.py + Skyfield Earth position
          │   ├── Angles → houses.py (pure GAST + obliquity, no ephemeris)
          │   ├── Arabic parts → arithmetic on cached values
          │   └── Harrington → independently sourced Keplerian propagation
          └── Post-processing: equatorial convert, sidereal, XYZ, radians
```

### 2.3 Key State Management (`state.py`)

Global state includes:
- `_CALC_MODE`: `"auto"` | `"skyfield"` | `"leb"` (proposed: add `"horizons"`)
- `_LEB_READER`: singleton LEBReader instance
- `_PLANETS`: loaded Skyfield SpiceKernel (DE440/441)
- `get_planets()`: lazy-loads DE440, **triggers download if missing**

The critical insight: `get_planets()` is called by `_calc_body()` unconditionally at line 1583. The Horizons path must intercept **before** `_calc_body()` to avoid triggering the DE440 download.

### 2.4 Coordinate Transform Functions (Reusable)

`fast_calc.py` contains standalone math functions decoupled from LEBReader:

| Function | Line | Description | Reusable? |
|----------|------|-------------|-----------|
| `_mean_obliquity_iau2006` | 84 | IAU 2006 obliquity polynomial | Yes |
| `_cartesian_to_spherical` | 112 | (x,y,z) → (lon°, lat°, dist AU) | Yes |
| `_cartesian_velocity_to_spherical` | 126 | (pos, vel) → (dlon°/d, dlat°/d, ddist AU/d) | Yes |
| `_rotate_equatorial_to_ecliptic` | 167 | R_x(ε) rotation | Yes |
| `_rotate_icrs_to_ecliptic_j2000` | 188 | ICRS → J2000 ecliptic (fixed ε) | Yes |
| `_apply_aberration` | 195 | Special-relativistic stellar aberration | Yes |
| `_apply_gravitational_deflection` | 359 | PPN light bending (Sun+Jup+Sat) | **Refactor needed** |
| `_get_skyfield_frame_data` | 460 | Precession-nutation matrix + dpsi/deps | Yes (needs Skyfield timescale) |
| `_get_precession_matrix` | 504 | ICRS → mean equatorial of date | Yes (needs Skyfield timescale) |
| `_mat3_vec3` | 540 | 3x3 matrix × vector | Yes |
| `_cotrans` | 552 | Spherical coordinate transform | Yes |
| `_precess_ecliptic` | 589 | Lieske ecliptic precession | Yes |
| `_calc_ayanamsa_from_leb` | 654 | General precession + J2000 offset | Yes (rename, pass mode) |
| `_vec3_sub` / `_vec3_dist` | 102/107 | Vector arithmetic | Yes |

**13 of 15 functions** are fully reusable or only depend on Skyfield's timescale (which is lightweight — no BSP file needed). Only `_apply_gravitational_deflection` is coupled to LEBReader's `eval_body()` method and needs refactoring.

### 2.5 Body Dependency Classification

Each body type in `_calc_body()` has been classified by its data dependency:

| Body | IDs | Classification | What It Needs |
|------|-----|----------------|---------------|
| Standard planets | 0-9, 14 | **ephemeris-dependent** | DE440 planet positions via Skyfield |
| SE_MEAN_NODE | 10 | **analytical** | ERFA/IERS lunar arguments and independent geometry |
| SE_TRUE_NODE | 11 | **ephemeris-dependent** | Moon state vectors (r, v) for angular momentum |
| SE_MEAN_APOG (Lilith) | 12 | **analytical** | ERFA/IERS lunar arguments and independent geometry |
| SE_OSCU_APOG (True Lilith) | 13 | **ephemeris-dependent** | Moon state vectors for eccentricity vector |
| SE_INTP_APOG | 21 | **analytical** | Versioned Delaunay compatibility series and pinned refinement |
| SE_INTP_PERG | 22 | **analytical** | Separate Delaunay compatibility series and pinned refinement |
| SE_CHIRON, SE_PHOLUS | 15-16 | **ephemeris-dependent** | Minor body positions |
| SE_CERES..SE_VESTA | 17-20 | **ephemeris-dependent** | Minor body positions |
| Asteroids | SE_AST_OFFSET+N | **ephemeris-dependent** | Minor body positions |
| Historical hypothetical bodies | 40–58 | **analytical** | Runtime element/convention registry with per-body provenance status |
| Fixed stars | FIXSTAR_OFFSET+N | **ephemeris-dependent** | Earth position for parallax/aberration |
| Angles (ASC, MC, Vertex) | ANGLE_OFFSET+N | **analytical** | GAST + obliquity from pyerfa (zero file deps) |
| Arabic parts | ARABIC_OFFSET+N | **analytical** | Arithmetic on cached longitudes |
| Planetary moons | SE_MOON_* | **ephemeris-dependent** | Satellite SPK + DE440 |

---

## 3. JPL Horizons API Overview

### 3.1 Endpoint

```
GET https://ssd.jpl.nasa.gov/api/horizons.api
```

Version 1.3 (2025 June). Free, no API key required. No documented rate limits (public best-effort service).

### 3.2 Relevant Configuration for This Backend

We use `EPHEM_TYPE='VECTORS'` to get Cartesian state vectors in ICRS:

```
?format=json
&COMMAND='499'              # Mars (planet center)
&EPHEM_TYPE='VECTORS'
&CENTER='500@399'           # Geocentric
&TLIST='2451545.0'          # Single JD
&TLIST_TYPE='JD'
&VEC_TABLE='2'              # x,y,z,vx,vy,vz
&OUT_UNITS='AU-D'           # AU and days
&VEC_CORR='NONE'            # Geometric (we apply corrections locally)
&CSV_FORMAT='YES'
&REF_SYSTEM='ICRF'
&TIME_TYPE='UT'             # or 'TDB' for swe_calc()
```

Response contains state vector data between `$$SOE` and `$$EOE` markers in the `result` field (JSON string).

### 3.3 Body Command Mapping

| libephemeris Body | Horizons COMMAND | Notes |
|-------------------|-----------------|-------|
| SE_SUN (0) | `'10'` | Sun center |
| SE_MOON (1) | `'301'` | Moon center |
| SE_MERCURY (2) | `'199'` | Mercury center |
| SE_VENUS (3) | `'299'` | Venus center |
| SE_EARTH (14) | `'399'` | Earth center |
| SE_MARS (4) | `'499'` | Mars center |
| SE_JUPITER (5) | `'599'` | Jupiter center (not barycenter 5) |
| SE_SATURN (6) | `'699'` | Saturn center |
| SE_URANUS (7) | `'799'` | Uranus center |
| SE_NEPTUNE (8) | `'899'` | Neptune center |
| SE_PLUTO (9) | `'999'` | Pluto center |
| SE_CHIRON (15) | `'2060;'` | Small-body syntax |
| SE_PHOLUS (16) | `'5145;'` | Small-body syntax |
| SE_CERES (17) | `'1;'` | Asteroid #1 |
| SE_PALLAS (18) | `'2;'` | Asteroid #2 |
| SE_JUNO (19) | `'3;'` | Asteroid #3 |
| SE_VESTA (20) | `'4;'` | Asteroid #4 |
| SE_AST_OFFSET+N | `'N;'` | Generic asteroid by number |

Note: Horizons planet center IDs (199, 299, 599, etc.) return the **physical center** of the planet directly, so no COB (center-of-body) correction is needed — unlike the Skyfield path which uses system barycenters and must apply COB corrections from `planet_centers_*.bsp`.

### 3.4 Bodies NOT Supported by Horizons

These must be handled analytically or trigger fallback:

- **Analytical (no HTTP needed):** Mean Node, Mean Apogee/Lilith, Harrington, Angles, Arabic Parts
- **Computed from Moon vectors (1 HTTP fetch):** True Node, Osculating Apogee, Interpolated Apogee/Perigee
- **Fallback to Skyfield (triggers DE440 download):** Fixed stars, Planetary moons, SEFLG_TOPOCTR

---

## 4. Technical Design

### 4.1 Activation Logic

```
User calls swe_calc_ut(jd, body, flags)
  │
  ├── LEB path (if .leb configured and mode != "skyfield")
  │
  ├── Horizons path (NEW):
  │   ├── mode == "horizons" → always use Horizons
  │   ├── mode == "auto" AND no DE440 locally → use Horizons
  │   ├── mode == "auto" AND DE440 exists → skip, use Skyfield
  │   └── mode == "skyfield" → skip
  │
  └── Skyfield path (existing, triggers DE440 download if missing)
```

On first Horizons activation, a one-time info message is logged:

```
INFO: No local ephemeris file found. Using NASA JPL Horizons API
      (requires internet, slower). For offline/faster calculations:
        python -m libephemeris download    # downloads DE440 (~114 MB)
```

### 4.2 Computation Pipeline

For ephemeris-dependent bodies, the Horizons backend must replicate the full Skyfield pipeline to achieve identical precision (~0.001"):

```
Step 1: Fetch geometric state vectors via HTTP (parallel)
        ├── Target body: ICRS barycentric geometric (CENTER='@0', VEC_CORR='NONE')
        ├── Earth: ICRS barycentric geometric
        ├── Sun: ICRS barycentric geometric          ┐
        ├── Jupiter bary: ICRS barycentric geometric  ├── Deflectors
        └── Saturn bary: ICRS barycentric geometric   ┘

Step 2: Light-time iteration (local)
        ├── Compute geometric geocentric: target_bary - earth_bary
        ├── Initial light-time estimate: |geo| / c
        ├── Iterate: re-fetch target at (jd - lt) → new geo → new lt
        │   (Horizons TLIST supports fractional JD, so each iteration
        │    is a new HTTP request. Typically converges in 2-3 iterations.
        │    Cache prevents redundant fetches.)
        └── Final: light-time-corrected geocentric position

Step 3: Gravitational deflection (local, PPN model)
        ├── Deflectors: Sun (GM=1.327e11), Jupiter (GM=1.267e8), Saturn (GM=3.794e7)
        ├── Deflector positions at observation time (from Step 1 cache)
        ├── Deflector positions at closest-approach time (additional fetch if needed)
        └── Apply PPN formula: δr = (2GM/c²) × [...] 
        Uses refactored _apply_gravitational_deflection() from fast_calc.py

Step 4: Stellar aberration (local, special-relativistic)
        ├── Earth velocity from barycentric state vector (Step 1)
        └── Apply _apply_aberration() from fast_calc.py (reuse as-is)

Step 5: Frame rotation (local)
        ├── ICRS → true equatorial of date: precession-nutation matrix
        │   via _get_skyfield_frame_data() (needs only Skyfield timescale, no BSP)
        ├── Equatorial → ecliptic: _rotate_equatorial_to_ecliptic(ε_true)
        └── Alternative paths for SEFLG_J2000, SEFLG_EQUATORIAL, SEFLG_ICRS

Step 6: Spherical conversion (local)
        ├── _cartesian_to_spherical() → (lon°, lat°, dist AU)
        └── _cartesian_velocity_to_spherical() → (dlon°/d, dlat°/d, ddist AU/d)

Step 7: Post-processing (local)
        ├── Sidereal correction (if SEFLG_SIDEREAL): subtract ayanamsha
        └── Output format (SEFLG_XYZ, SEFLG_RADIANS)
```

### 4.3 Light-Time Iteration Strategy

The naive approach (re-fetch target at retarded time for each iteration) would require 2-3 extra HTTP requests per body. This is expensive.

**Optimization: Single-fetch with interpolation.** Instead of iterating:

1. Fetch target at 3 closely-spaced times: `[jd - 0.02, jd, jd + 0.02]` (using `TLIST`)
2. Compute geometric light-time from the `jd` position
3. Interpolate (quadratic) the target position at `jd - lt`
4. One iteration is sufficient for sub-milliarcsecond precision

This reduces the total HTTP requests per body from 3-4 to **1** (with 3 time steps in a single TLIST request).

### 4.4 Batch Prefetch Strategy

When the first `swe_calc_ut()` call arrives for a given JD:

```python
# Pseudo-code for the batch prefetch optimization
def _prefetch_for_jd(jd, iflag):
    """Pre-fetch all likely-needed bodies for this JD in parallel."""
    bodies_to_fetch = [
        # Always: Earth + 3 deflectors (barycentric geometric)
        ('399', '@0', 'NONE'),  # Earth
        ('10', '@0', 'NONE'),   # Sun (deflector)
        ('5', '@0', 'NONE'),    # Jupiter bary (deflector)
        ('6', '@0', 'NONE'),    # Saturn bary (deflector)
        # Standard chart bodies (barycentric geometric for LT iteration)
        ('10', '@0', 'NONE'),   # Sun
        ('301', '@0', 'NONE'),  # Moon
        ('199', '@0', 'NONE'),  # Mercury
        ('299', '@0', 'NONE'),  # Venus
        ('499', '@0', 'NONE'),  # Mars
        ('599', '@0', 'NONE'),  # Jupiter
        ('699', '@0', 'NONE'),  # Saturn
        ('799', '@0', 'NONE'),  # Uranus
        ('899', '@0', 'NONE'),  # Neptune
        ('999', '@0', 'NONE'),  # Pluto
    ]
    # De-duplicate and fire all in parallel via ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = {pool.submit(fetch_state_vector, *req): req for req in unique(bodies_to_fetch)}
        for future in as_completed(futures):
            cache.store(futures[future], future.result())
```

This way a full chart (Sun through Pluto + nodes + Lilith + angles) completes in:
- **Wave 1:** ~300-600ms (14 unique HTTP requests, 8 parallel)
- **Subsequent calls same JD:** 0ms (cache hit)

### 4.5 Flag Dispatch

| SEFLG Flag | Effect on Horizons Request | Effect on Local Pipeline |
|------------|--------------------------|------------------------|
| *(default)* | CENTER='@0', VEC_CORR='NONE' | Full pipeline (LT + deflect + aberr + frame) |
| SEFLG_HELCTR | CENTER='@10' (Sun) | Skip LT iteration, skip deflection, skip aberration |
| SEFLG_BARYCTR | CENTER='@0' (SSB) | Skip LT iteration, skip deflection, skip aberration |
| SEFLG_TRUEPOS | Same as default | Skip aberration (keep LT + deflection) |
| SEFLG_NOABERR | Same as default | Skip aberration step |
| SEFLG_NOGDEFL | Same as default | Skip deflection step |
| SEFLG_SIDEREAL | Same as default | Subtract ayanamsha after ecliptic conversion |
| SEFLG_EQUATORIAL | Same as default | Skip ecliptic rotation, output RA/DEC |
| SEFLG_J2000 | Same as default | Use J2000 obliquity instead of date |
| SEFLG_SPEED | Same as default | Compute velocity from state vector (always available) |
| SEFLG_TOPOCTR | **raise KeyError** | Fallback to Skyfield (Earth orientation params needed) |
| SEFLG_XYZ | Same as default | Output Cartesian instead of spherical |
| SEFLG_RADIANS | Same as default | Convert degrees to radians |
| SEFLG_NONUT | Same as default | Use mean obliquity, skip nutation in frame rotation |
| SEFLG_ICRS | Same as default | Skip frame rotation entirely |
| SEFLG_MOSEPH | Stripped (ignored) | — |

### 4.6 Body Routing Table

| Body | Route | HTTP Requests | Local Computation |
|------|-------|--------------|-------------------|
| SE_SUN..SE_PLUTO, SE_EARTH | Horizons VECTORS | 1 target + Earth + 3 deflectors (cached) | Full pipeline |
| SE_CHIRON, SE_PHOLUS | Horizons VECTORS | 1 target (+ cached deps) | Full pipeline |
| SE_CERES..SE_VESTA | Horizons VECTORS | 1 target (+ cached deps) | Full pipeline |
| SE_AST_OFFSET+N | Horizons VECTORS | 1 target (+ cached deps) | Full pipeline |
| SE_MEAN_NODE | **Analytical** | 0 | `lunar.calc_mean_lunar_node()` + nutation |
| SE_MEAN_APOG | **Analytical** | 0 | `lunar.calc_mean_lilith_with_latitude()` + nutation |
| SE_TRUE_NODE | Horizons (Moon) | 1 Moon geometric (+ cached Earth) | Angular momentum r×v → node longitude |
| SE_OSCU_APOG | Horizons (Moon) | 1 Moon geometric (+ cached Earth) | Eccentricity vector → apogee |
| SE_INTP_APOG | Mixed | Moon-state samples (+ cached Earth) | Symmetric JPL-state smoothing |
| SE_INTP_PERG | Mixed | Moon-state samples (+ cached Earth) | Symmetric JPL-state smoothing |
| Harrington | **Analytical** | 0 | Published Harrington (1988) elements |
| Unsupported compatibility IDs | **Unsupported** | 0 | `UnknownBodyError` |
| Angles (ASC, MC, etc.) | **Analytical** | 0 | `houses.py` (GAST + obliquity via pyerfa) |
| Arabic parts | **Analytical** | 0 | Arithmetic on cached longitudes |
| Fixed stars | **raise KeyError** | — | Fallback → Skyfield (triggers DE440 download) |
| Planetary moons | **raise KeyError** | — | Fallback → Skyfield (triggers DE440 download) |

### 4.7 Caching Strategy

**In-memory LRU cache** (`collections.OrderedDict`, max 2048 entries):

- **Key:** `(jd_rounded_to_12_decimals, naif_command, center, vec_corr)`
- **Value:** `StateVector(x, y, z, vx, vy, vz, lt, range, rdot)`
- **Eviction:** LRU when max size exceeded
- **Lifetime:** Process lifetime (cleared on `close()`)
- **Thread safety:** `threading.Lock` on cache access

The cache is effective because:
- A typical chart calls `swe_calc_ut()` 10-15 times with the **same JD**
- Earth + deflector positions are shared across all body calculations for a given JD
- The Moon state vector is shared between True Node, Osculating Apogee, and Interpolated Apogee/Perigee

### 4.8 HTTP Client Design

```python
class HorizonsClient:
    API_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
    
    def __init__(self, max_cache_size=2048, max_workers=8, timeout=30):
        self._cache: OrderedDict           # LRU cache
        self._cache_lock: threading.Lock   # Thread safety
        self._timeout: int                 # HTTP timeout in seconds
        self._max_workers: int             # ThreadPoolExecutor parallelism
    
    def fetch_state_vector(self, command, jd, center='@0', vec_corr='NONE',
                           time_type='UT') -> StateVector:
        """Fetch a single state vector, with caching and retry."""
    
    def fetch_batch(self, requests: list[FetchRequest]) -> dict[tuple, StateVector]:
        """Fetch multiple state vectors in parallel via ThreadPoolExecutor."""
    
    def shutdown(self):
        """Clean up resources."""
    
    def clear_cache(self):
        """Clear the in-memory cache."""
```

**HTTP specifics:**
- Uses `urllib.request` (stdlib, no extra dependencies)
- SSL via `certifi` (already a transitive dependency of Skyfield)
- User-Agent: `libephemeris/{version}`
- Retry: max 2 retries with backoff (0.5s, 1.0s)
- Timeout: 30s per request (configurable)
- Response parsing: JSON → extract `result` string → find `$$SOE`/`$$EOE` → parse CSV line

### 4.9 Error Handling

| Error Condition | Handling | User-Facing Message |
|----------------|---------|-------------------|
| Network timeout | Retry 2× with exponential backoff (0.5s, 1s), then raise | `ConnectionError: "Horizons API request timed out after 3 attempts. Check your internet connection, or download local ephemeris: python -m libephemeris download"` |
| DNS resolution failure | Raise immediately (no retry) | `ConnectionError: "Cannot reach ssd.jpl.nasa.gov. No internet connection? Download local ephemeris: python -m libephemeris download"` |
| HTTP 429 (rate limit) | Backoff 2s, retry up to 3× | `ConnectionError: "Horizons API rate limit exceeded. Please wait or download local ephemeris."` |
| Horizons API error (JSON `error` field) | Raise with Horizons message | `ValueError: "Horizons API error: {message}"` |
| Body not found on Horizons | raise `KeyError` → triggers Skyfield fallback | (silent fallback, debug log only) |
| Unsupported flag (TOPOCTR) | raise `KeyError` → triggers Skyfield fallback | (silent fallback, debug log only) |
| Malformed response (parse failure) | Retry once, then raise | `ValueError: "Failed to parse Horizons API response"` |
| SSL certificate error | Use certifi CA bundle | Should not occur with certifi |

---

## 5. Files to Create

### 5.1 `libephemeris/horizons_backend.py` (~600-800 lines)

New module containing:

**a) `StateVector` dataclass** — holds `(x, y, z, vx, vy, vz)` in AU and AU/day, ICRS frame.

**b) `HorizonsClient` class** — HTTP client with LRU cache, connection pooling, retry logic, parallel fetch. See section 4.8.

**c) Body command mapping** — `_HORIZONS_COMMAND_MAP` dict mapping `SE_*` constants to Horizons `COMMAND` strings. See section 3.3.

**d) URL builder** — `_build_vectors_url()` constructs the full API URL with proper parameter encoding (semicolons as `%3B`, spaces as `%20`).

**e) Response parser** — `_parse_vectors_response()` extracts CSV data between `$$SOE`/`$$EOE` markers from the JSON `result` string.

**f) Calculation pipeline** — `horizons_calc_ut()` and `horizons_calc_tt()` implement the full pipeline from section 4.2, reusing math functions from `fast_calc.py`.

**g) Body routing** — dispatches to independent analytical functions for Mean
Node, Mean Apogee, Harrington, Angles, and Arabic Parts. It fetches Moon vectors
for state-derived lunar points and rejects unsupported hypothetical IDs.

**h) Flag dispatch** — `_resolve_horizons_params()` maps SEFLG flags to Horizons API parameters and local pipeline options. See section 4.5.

---

## 6. Files to Modify

### 6.1 `libephemeris/state.py`

**a) Extend valid calc modes:**

```python
# Line 171 (current)
_VALID_CALC_MODES = ("auto", "skyfield", "leb")

# Proposed
_VALID_CALC_MODES = ("auto", "skyfield", "leb", "horizons")
```

**b) Add ephemeris-local-availability check:**

```python
def _is_ephemeris_locally_available() -> bool:
    """Check if DE440/441 BSP exists locally WITHOUT triggering download."""
    data_dir = _get_data_dir()
    eph_file = _get_effective_ephemeris_file()
    if os.path.exists(os.path.join(data_dir, eph_file)):
        return True
    if _EPHEMERIS_PATH and os.path.exists(os.path.join(_EPHEMERIS_PATH, eph_file)):
        return True
    return False
```

**c) Add Horizons client singleton management:**

```python
_HORIZONS_CLIENT: Optional["HorizonsClient"] = None
_HORIZONS_WARNED: bool = False

def get_horizons_client() -> Optional["HorizonsClient"]:
    """Get Horizons client if mode requires it."""
    mode = get_calc_mode()
    if mode == "skyfield":
        return None
    if mode == "horizons":
        return _get_or_create_horizons_client()
    if mode == "auto":
        # Use Horizons only when DE440 is not available locally
        if not _is_ephemeris_locally_available():
            return _get_or_create_horizons_client()
    return None

def _get_or_create_horizons_client():
    global _HORIZONS_CLIENT, _HORIZONS_WARNED
    if _HORIZONS_CLIENT is None:
        if not _HORIZONS_WARNED:
            logger.info(
                "No local ephemeris file found. Using NASA JPL Horizons API "
                "(requires internet, slower). For offline/faster calculations:\n"
                "  python -m libephemeris download"
            )
            _HORIZONS_WARNED = True
        from .horizons_backend import HorizonsClient
        _HORIZONS_CLIENT = HorizonsClient()
    return _HORIZONS_CLIENT
```

**d) Extend `close()` to clean up Horizons client** (around line 1264):

```python
# Add to close():
global _HORIZONS_CLIENT, _HORIZONS_WARNED
if _HORIZONS_CLIENT is not None:
    _HORIZONS_CLIENT.shutdown()
_HORIZONS_CLIENT = None
_HORIZONS_WARNED = False
```

### 6.2 `libephemeris/planets.py`

**Insert Horizons path in `swe_calc_ut()`** (after LEB block at line 848, before Skyfield block at line 850):

```python
    # --- END LEB fast path --- (existing line 848)

    # --- Horizons path: use API when no local ephemeris ---
    from .state import get_horizons_client

    h_client = get_horizons_client()
    if h_client is not None:
        try:
            from . import horizons_backend

            result = horizons_backend.horizons_calc_ut(
                h_client, tjdut, planet, flags
            )
            get_logger().debug("body=%d jd=%.1f source=Horizons", planet, tjdut)
            return result
        except KeyError as _hz_err:
            get_logger().debug(
                "body=%d jd=%.1f source=Horizons->fallback reason=%s",
                planet, tjdut, _hz_err,
            )
    # --- END Horizons path ---

    # Validate JD range for bodies that use the JPL ephemeris (existing line 850)
```

**Same pattern in `swe_calc()`** (after LEB block at line 938).

### 6.3 `libephemeris/fast_calc.py`

**Refactor `_apply_gravitational_deflection()`** (line 359) to decouple from LEBReader:

```python
# CURRENT signature:
def _apply_gravitational_deflection(
    geo, earth_bary, jd_tt, light_time, reader
):
    # ... uses reader.eval_body(naif_id, jd) to get deflector positions ...

# PROPOSED signature:
def _apply_gravitational_deflection(
    geo: Tuple[float, float, float],
    earth_bary: Tuple[float, float, float],
    jd_tt: float,
    light_time: float,
    get_deflector_pos: Callable[[int, float], Tuple[float, float, float]],
):
    """Apply PPN gravitational light deflection.
    
    Args:
        geo: Light-time-corrected geocentric ICRS position (AU).
        earth_bary: Earth barycentric ICRS position (AU).
        jd_tt: Julian Day TT of observation.
        light_time: Light travel time in days.
        get_deflector_pos: Callable(naif_id, jd_tt) -> (x, y, z) in AU,
                           returns barycentric ICRS position of a deflector body.
    """
```

The LEB caller wraps `reader.eval_body()` into the callable. The Horizons caller wraps its HTTP fetch + cache. The PPN math (~60 lines) remains **identical**.

The existing call site in `_pipeline_icrs()` (line 837) is updated to pass a lambda:

```python
# Current (line 837):
geo = _apply_gravitational_deflection(geo, observer, jd_tt, lt, reader)

# Proposed:
geo = _apply_gravitational_deflection(
    geo, observer, jd_tt, lt,
    lambda naif_id, jd: reader.eval_body(naif_id, jd)[:3],
)
```

### 6.4 `libephemeris/lunar.py`

**Add alternative entry points** for True Node and Osculating Apogee that accept pre-computed Moon state vectors instead of calling `get_planets()` internally:

```python
def calc_true_lunar_node_from_vectors(
    jd_tt: float, moon_r: tuple, moon_v: tuple
) -> float:
    """Compute True Node longitude from pre-supplied Moon state vectors.
    
    Args:
        jd_tt: Julian Day TT.
        moon_r: Moon geocentric ecliptic position (x, y, z) in AU.
        moon_v: Moon geocentric ecliptic velocity (vx, vy, vz) in AU/day.
    
    Returns:
        True Node longitude in degrees [0, 360).
    """
    # Same angular momentum h = r × v logic, just using provided vectors
    # instead of calling get_planets() + Skyfield
```

Similarly for `calc_true_lilith_from_vectors()` and `calc_osculating_perigee_from_vectors()`.

This keeps the existing functions unchanged (no risk of regression) and adds parallel entry points for the Horizons backend.

### 6.5 `libephemeris/__init__.py`

- Update module docstring to mention `"horizons"` mode (lines 12-19)
- `"horizons"` is already valid via the `set_calc_mode()` / `get_calc_mode()` API — no new exports needed

---

## 7. Performance Characteristics

### 7.1 Latency Comparison

| Scenario | LEB | Skyfield | Horizons |
|----------|-----|----------|----------|
| Single `swe_calc_ut()` call | ~5 μs | ~0.5 ms | ~300 ms (first), ~0 ms (cached) |
| Full chart (15 bodies, same JD) | ~75 μs | ~7.5 ms | ~600 ms (first), ~0 ms (cached) |
| 1000 charts (different JDs) | ~75 ms | ~7.5 s | ~600 s (serial) / ~75 s (8 parallel) |

### 7.2 HTTP Request Budget (Full Chart, Single JD)

```
First call to swe_calc_ut(jd, SE_SUN, ...):
  → Triggers batch prefetch for jd
  → 14 unique HTTP requests (8 parallel, 2 waves):
      Wave 1 (8 parallel): Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn, Uranus
      Wave 2 (6 parallel): Neptune, Pluto, Earth(bary), Sun(bary/defl), Jup(bary/defl), Sat(bary/defl)
  → Total: ~300-600 ms

Subsequent calls (SE_MOON, SE_MARS, ... SE_PLUTO):
  → All cache hits → ~0 ms each

Analytical bodies (Mean Node, Mean Apogee, Angles):
  → Zero HTTP requests → <1 ms each

True Node / Osculating Apogee:
  → Moon already cached from batch → ~0 ms
```

### 7.3 Memory Footprint

- LRU cache at max capacity (2048 entries × ~200 bytes per StateVector): ~400 KB
- Thread pool: 8 threads × ~1 MB stack: ~8 MB
- Negligible compared to DE440 loaded in memory (~114 MB)

---

## 8. Test Plan

### 8.1 Unit Tests (Mocked HTTP)

**`tests/test_horizons_client.py`:**

| Test | Description |
|------|-------------|
| `test_build_vectors_url` | Verify URL construction with proper encoding (semicolons, spaces) |
| `test_parse_vectors_response` | Parse real Horizons JSON response fixtures |
| `test_parse_vectors_csv` | Parse CSV data between $$SOE/$$EOE markers |
| `test_cache_lru_eviction` | Verify LRU eviction at max_cache_size |
| `test_cache_hit` | Verify cache returns stored values |
| `test_retry_on_timeout` | Verify retry with backoff on timeout |
| `test_retry_on_429` | Verify retry with backoff on rate limit |
| `test_error_on_network_failure` | Verify ConnectionError with helpful message |
| `test_body_command_mapping` | Verify all SE_* → COMMAND mappings |
| `test_flag_to_params` | Verify SEFLG_* → Horizons parameter mapping |
| `test_unsupported_body_raises_keyerror` | Fixed stars, planetary moons raise KeyError |
| `test_unsupported_flag_raises_keyerror` | SEFLG_TOPOCTR raises KeyError |
| `test_analytical_bodies_no_http` | Mean Node, Mean Apogee, Angles don't trigger HTTP |
| `test_batch_deduplication` | Verify duplicate requests are collapsed |
| `test_parallel_fetch` | Verify ThreadPoolExecutor is used |

### 8.2 Integration Tests (Live HTTP, marked `@pytest.mark.slow`)

**`tests/test_horizons_precision.py`:**

| Test | Description | Tolerance |
|------|-------------|-----------|
| `test_planets_vs_skyfield` | All 10 planets (Sun-Pluto) at 5 dates, compare vs Skyfield | 0.001" lon/lat |
| `test_chiron_vs_skyfield` | Chiron at 3 dates | 0.01" |
| `test_asteroids_vs_skyfield` | Ceres, Pallas, Juno, Vesta at 3 dates | 0.01" |
| `test_true_node_vs_skyfield` | True Node at 5 dates | 0.001" |
| `test_oscu_apogee_vs_skyfield` | Osculating Apogee at 5 dates | 0.01" |
| `test_helctr_flag` | Heliocentric positions match Skyfield | 0.001" |
| `test_baryctr_flag` | Barycentric positions match Skyfield | 0.001" |
| `test_sidereal_flag` | Sidereal positions match Skyfield | 0.001" |
| `test_equatorial_flag` | RA/DEC match Skyfield | 0.001" |
| `test_j2000_flag` | J2000 ecliptic match Skyfield | 0.001" |
| `test_speed_values` | Velocity components match Skyfield | 0.001"/day |

### 8.3 Fallback Tests

**`tests/test_horizons_fallback.py`:**

| Test | Description |
|------|-------------|
| `test_auto_mode_with_local_de440` | With DE440 present, Horizons is NOT used |
| `test_auto_mode_without_de440` | Without DE440, Horizons IS used transparently |
| `test_horizons_mode_explicit` | `set_calc_mode("horizons")` always uses Horizons |
| `test_skyfield_mode_explicit` | `set_calc_mode("skyfield")` never uses Horizons |
| `test_fixed_star_fallback` | Fixed star triggers Skyfield fallback |
| `test_topoctr_fallback` | SEFLG_TOPOCTR triggers Skyfield fallback |
| `test_info_message_once` | Info message logged only on first use |
| `test_close_resets_client` | `close()` clears Horizons client |

---

## 9. Dependencies

**No new dependencies.** Everything is implementable with:

| Module | Source | Already Used? |
|--------|--------|--------------|
| `urllib.request` | stdlib | Yes (in `spk.py`) |
| `json` | stdlib | Yes |
| `concurrent.futures` | stdlib | No (new usage) |
| `collections.OrderedDict` | stdlib | Yes |
| `ssl` | stdlib | Yes (in `spk.py`) |
| `certifi` | transitive dep of Skyfield | Yes |
| `threading` | stdlib | Yes (in `state.py`) |

---

## 10. Limitations and Known Constraints

1. **SEFLG_TOPOCTR** — Not implementable via Horizons without replicating Earth orientation parameters. Fallback to Skyfield (triggers DE440 download on first use of this flag).

2. **Fixed stars** — Horizons does not have a generic fixed-star catalog query equivalent to the Hipparcos-based `fixed_stars.py`. Fallback to Skyfield.

3. **Planetary moons** (Galilean, Titan, etc.) — Require satellite-specific SPK kernels not available via the Horizons VECTORS API. Fallback to Skyfield.

4. **Offline usage** — The backend requires internet connectivity. When offline with no local DE440, the error message guides the user to `python -m libephemeris download`.

5. **Horizons API availability** — NASA's Horizons is a best-effort public service. Extended outages, maintenance windows, or rate limiting could affect availability. The info message at first use encourages downloading local files.

6. **Precision for very distant dates** — Horizons may use different Delta T models or Earth orientation parameters than Skyfield for extreme dates (before 1550 or after 2650). For the casual user this is irrelevant.

7. **Asteroids not in Horizons** — Very recently discovered or rare objects might not be in Horizons database. Fallback to Skyfield's Keplerian propagation.

---

## 11. Future Extensions

- **Disk cache** — Optional SQLite persistence layer for frequently-queried JDs across sessions.
- **SEFLG_TOPOCTR support** — Use Horizons' `SITE_COORD` parameter with geodetic coordinates. Requires validation of precision parity with Skyfield's topocentric path.
- **WebSocket/streaming API** — If Horizons adds a persistent connection API, could dramatically reduce latency.
- **Partial DE440** — Download only the segments needed for a specific date range (~10 MB instead of 114 MB) as a middle ground between zero-install and full DE440.

---

## 12. Implementation Order

| Phase | Description | Effort | Depends On |
|-------|-------------|--------|------------|
| 1 | Refactor `_apply_gravitational_deflection` in `fast_calc.py` | 0.5 day | — |
| 2 | Add `_from_vectors()` entry points in `lunar.py` | 1 day | — |
| 3 | Create `horizons_backend.py`: HTTP client, parsing, caching | 2-3 days | — |
| 4 | Create `horizons_backend.py`: calculation pipeline + body routing | 2 days | Phase 1, 2, 3 |
| 5 | Integrate in `state.py` + `planets.py` (mode, auto-detection) | 1 day | Phase 4 |
| 6 | Unit tests (mocked HTTP) | 1-2 days | Phase 5 |
| 7 | Integration tests (live HTTP, cross-validation vs Skyfield) | 1-2 days | Phase 5 |
| **Total** | | **8-10 days** | |

Phases 1, 2, 3 can be worked on in parallel as they have no mutual dependencies.
