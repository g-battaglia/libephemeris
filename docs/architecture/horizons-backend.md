# Horizons API Backend

> **Requires:** Internet connection to ssd.jpl.nasa.gov

## Overview

The Horizons backend enables zero-install ephemeris computation by fetching
state vectors from the [NASA JPL Horizons REST API](https://ssd.jpl.nasa.gov/horizons/).
When no local ephemeris files (DE440 or LEB) are available, the library
transparently fetches data from Horizons and computes apparent positions
using the same correction pipeline as the Skyfield/LEB paths.

## Calculation Modes

LibEphemeris supports 4 calculation modes, configured via `set_calc_mode()`
or the `LIBEPHEMERIS_MODE` environment variable:

| Mode | Behavior | Fails when |
|------|----------|-----------|
| `"auto"` (default) | LEB → Horizons (if no DE440) → Skyfield | never |
| `"leb"` | Sealed LEB-only persistent data; only traced local analytical/Keplerian models are allowed where LEB has no meaningful channel | no LEB resolvable, corrupt/missing required data, or core range miss |
| `"horizons"` | Prefer Horizons; unsupported bodies/flags fall back to Skyfield | no internet |
| `"skyfield"` | Always Skyfield/DE440 | DE440 not available |

### Auto Mode Flow

```
calc_ut(jd, body, flags)
    |
    +-> LEB fast path (if .leb file configured)
    |   |-> success: return result
    |   |-> KeyError/ValueError: fall through
    |
    +-> Horizons path (if mode="horizons" OR (mode="auto" AND no DE440))
    |   |-> success: return result
    |   |-> KeyError: fall through (unsupported body/flag)
    |
    +-> Skyfield path (default fallback)
        |-> success: return result
        |-> exception: propagate to caller
```

This flow describes `auto`. Forced `leb` mode stops at the LEB boundary: it
does not create a Horizons client, open a DE kernel or planet-center BSP, use a
registered/automatic SPK, or invoke ASSIST. A curated body that is intentionally
served by an in-code model is returned with explicit `Keplerian` or
`Analytical` provenance.

### Configuration

```python
import libephemeris as swe

# Explicit mode
swe.set_calc_mode("horizons")

# Or via environment variable
# export LIBEPHEMERIS_MODE=horizons
```

## Supported Bodies

| Category | Bodies | Source |
|----------|--------|--------|
| Standard planets | Sun, Moon, Mercury-Pluto, Earth | Horizons VECTORS API |
| Asteroids | Chiron, Ceres, Pallas, Juno, Vesta | Horizons small-body syntax |
| Mean Node / Mean Apogee | MEAN_NODE (10), MEAN_APOG (12) | ERFA/IERS fundamental arguments (no HTTP) |
| Reviewed hypothetical bodies | IDs 40–47, 50–53, and 56 | Cited local runtime models (no HTTP) |

### Not Supported (fallback to Skyfield)

- **True Node, Osculating Apogee, Interpolated Apogee/Perigee** (11, 13, 21, 22) — require Moon state vectors
- **Fixed stars** — no Horizons equivalent
- **Planetary moons** — require satellite-specific SPK
- **Unverified historical hypothetical IDs** — IDs 48, 49, 54, 55, 57, and 58 are
  recognised but raise `UnknownBodyError`; they do not fall through to HTTP
- **FLG_TOPOCTR** — requires Earth orientation parameters
- **Hypothetical center conversions** — non-native heliocentric/geocentric or
  barycentric requests fall through to Skyfield for the Earth/Sun vector
  subtraction; they are never classified as unknown bodies

## Pipeline Architecture

For geocentric apparent positions, the Horizons backend replicates the
full Skyfield/LEB pipeline:

```
1. Fetch barycentric ICRS state vectors (parallel HTTP)
   - Target body
   - Earth (for geocentric conversion)
   - Sun, Jupiter, Saturn (gravitational deflectors)

2. Compute geometric geocentric vector
   target_bary - earth_bary

3. Light-time correction (single iteration)
   Re-fetch target at jd - light_time

4. Gravitational deflection (PPN model)
   Sun (GM=1.327e11) + Jupiter + Saturn
   Uses reusable _apply_deflection_horizons()

5. Stellar aberration (special-relativistic)
   Uses fast_calc._apply_aberration()

6. Frame rotation
   ICRS -> equatorial of date (precession-nutation via Skyfield timescale)
   equatorial -> ecliptic (true obliquity)

7. Spherical conversion
   Cartesian -> (lon, lat, dist) in degrees/AU
```

## HTTP Client

`HorizonsClient` in `libephemeris/horizons_backend.py`:

- **LRU cache**: 4096 entries, keyed by (JD, command, center)
- **Parallel fetch**: ThreadPoolExecutor with 8 workers
- **Retry**: 2 retries with exponential backoff (0.5s, 1.0s)
- **Timeout**: 30 seconds per request
- **Thread-safe**: Lock on cache access

A typical astrological chart (15 bodies, same JD) requires:
- **First call**: ~5 parallel HTTP requests (target + Earth + 3 deflectors)
- **Subsequent calls same JD**: 0 HTTP (all cached)
- **Total latency**: ~300-600ms for first chart, ~0ms for cached

## Flag Support

| Flag | Horizons behavior |
|------|-------------------|
| FLG_SPEED | Velocity from state vector differences |
| FLG_HELCTR | center='@10' (Sun center) |
| FLG_BARYCTR | center='@0' (SSB) |
| FLG_SIDEREAL | Subtract ayanamsha after ecliptic conversion |
| FLG_EQUATORIAL | Skip ecliptic rotation, output RA/Dec |
| FLG_J2000 | Use J2000 ecliptic frame |
| FLG_NOABERR | Skip aberration step |
| FLG_NOGDEFL | Skip deflection step |
| FLG_TRUEPOS | Skip light-time + aberration |
| FLG_TOPOCTR | **raises KeyError** -> fallback to Skyfield |

## Precision

Measured against Skyfield/DE440 reference (15,400 tests, 200 dates x 13 bodies x 6 flags):

| Mode | Max error | Notes |
|------|-----------|-------|
| Geocentric (default/sidereal/equatorial/J2000/no_aberr) | **0.0003"** | Excellent |
| Heliocentric | **0.027"** | Systematic offset: Horizons center=@10 vs Skyfield SSB-Sun |

### Velocity Precision

Velocities are computed via numerical differentiation of the apparent position
(dt=1 second). This is slightly less precise than the analytical Chebyshev
derivative used by the LEB/Skyfield paths.

| Body | Speed diff vs Skyfield/DE440 | Notes |
|------|--------------------------|-------|
| Moon | ~0.001°/day | Largest because Moon moves ~12°/day |
| Inner planets | < 0.0001°/day | Negligible |
| Outer planets | < 0.00001°/day | Negligible |

The Moon velocity difference of ~0.001°/day (~3.6"/day) is an inherent
limitation of the numerical derivative approach and cannot be improved
without analytical access to the ephemeris polynomials (which defeats
the purpose of the Horizons backend).

For astrological applications this is irrelevant: a 0.001°/day velocity
error translates to ~0.04 seconds of arc over a 1-minute time step.

## Testing

```bash
# Horizons vs Skyfield precision (needs internet)
leph test horizons precision          # 200 dates, ~45s
leph test horizons precision-quick    # 50 dates, ~15s

# Horizons vs LEB2 cross-validation
leph test horizons vs-leb             # 100 dates, ~30s

# Equivalent poe task shortcuts
poe test:horizons:core                # precision-quick (50 dates, ~15s)
poe test:horizons:fast                # precision (200 dates, ~45s)
poe test:horizons:full                # vs-leb (LEB2 cross-validation)
```

## Error Handling

| Error | Behavior |
|-------|----------|
| Network timeout | Retry 2x with backoff, then raise `ConnectionError` |
| DNS failure | Raise immediately with helpful message |
| Body not found on Horizons | Raise `KeyError` -> triggers Skyfield fallback |
| Unsupported flag (TOPOCTR) | Raise `KeyError` -> triggers Skyfield fallback |
| Horizons API error | Raise `KeyError` with Horizons error message |

## Limitations

1. **Requires internet** — offline use requires LEB or Skyfield/DE440
2. **Latency** — ~300-600ms for first chart (vs ~75us for LEB, ~7ms for Skyfield)
3. **Horizons API availability** — NASA best-effort public service, no SLA
4. **Heliocentric offset** — ~0.01-0.03" systematic difference vs Skyfield
   due to different Sun center definitions (physical center vs SSB offset)
5. **Date range** — Horizons covers most dates but may have gaps for extreme dates
6. **Moon velocity** — ~0.001°/day difference vs analytical Chebyshev derivative (Skyfield/DE440).
   Numerical differentiation of apparent position is inherently less precise than
   analytical Chebyshev polynomial derivatives. Irrelevant for astrological use.
7. **Unsupported bodies** — True Node (11), Osculating Apogee (13), Interpolated
   Apogee/Perigee (21, 22) fall through to Skyfield because they require Moon
   state vectors computed in a specific frame. Fixed stars and planetary moons
   also fall through.
