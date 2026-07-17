# Migration Guide: pyswisseph to LibEphemeris

This guide helps users migrate from `pyswisseph` (Python bindings to Swiss Ephemeris) to `libephemeris`.

## Table of Contents

- [Overview](#overview)
- [Quick Migration](#quick-migration)
- [API Differences](#api-differences)
- [Numerical model differences](#numerical-model-differences)
- [Known Compatibility Gaps](#known-compatibility-gaps)
- [Thread Safety with EphemerisContext](#thread-safety-with-ephemeriscontext)
- [FLG_MOSEPH (Moshier Ephemeris Flag)](#flg_moseph-moshier-ephemeris-flag)
- [Calculation Backend](#calculation-backend)
- [Migration Checklist](#migration-checklist)
- [Reporting Issues](#reporting-issues)

## Overview

`libephemeris` is designed as a **drop-in replacement** for `pyswisseph` in most common use cases. It provides:

- Same function names — call them on the `libephemeris` module exactly as you called them on `swisseph` (e.g. `swe.calc_ut`)
- Same flag and body constants (`FLG_*`, `SIDM_*`, `SUN`, `MOON`, etc.)
- Same return value structures
- Same behavior for global state management

However, there are important differences to be aware of when migrating.

---

## Quick Migration

### Minimal Change

In most cases, you can simply replace your import:

```python
# Before (pyswisseph)
import swisseph as swe

# After (libephemeris)
import libephemeris as swe
```

Your existing code should continue to work for common planetary calculations, house systems, and ayanamshas.

### Constants Import

```python
# Before (pyswisseph)
import swisseph as swe
planet = swe.SUN
flag = swe.FLG_SPEED

# After (libephemeris) - Option 1: Same style (constants have identical names)
import libephemeris as swe
planet = swe.SUN
flag = swe.FLG_SPEED

# After (libephemeris) - Option 2: Explicit constants import (recommended)
import libephemeris as swe
from libephemeris.constants import SUN, FLG_SPEED
```

---

## API Differences

### Function Names

Function names are identical to pyswisseph. After `import libephemeris as swe`,
existing call sites keep working unchanged — `swe.calc_ut(...)`, `swe.houses(...)`,
`swe.julday(...)`, `swe.set_sid_mode(...)`, `swe.get_ayanamsa_ut(...)`, and so on.
Functions are exposed **unprefixed**; there is no `swe_`-prefixed alias form.

### Constant Names

LibEphemeris uses the **same constant names** as pyswisseph — no prefix changes needed:

| pyswisseph | LibEphemeris |
|------------|--------------|
| `swe.SUN` | `SUN` (0) |
| `swe.MOON` | `MOON` (1) |
| `swe.FLG_SPEED` | `FLG_SPEED` |
| `swe.FLG_SIDEREAL` | `FLG_SIDEREAL` |
| `swe.FLG_TOPOCTR` | `FLG_TOPOCTR` |
| `swe.SIDM_LAHIRI` | `SIDM_LAHIRI` (computed predefined mode) |

### House Cusp Array Indexing

Both libraries return the same shape: a 12-element cusp tuple, 0-indexed
(cusp 1 at index 0), plus an 8-element `ascmc` tuple. No index shifting is
needed when migrating:

```python
# pyswisseph: 12 elements, 0-indexed
cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, b'P')
cusp1_swe = cusps_swe[0]

# libephemeris: identical
cusps, ascmc = ephem.houses(jd, lat, lon, b'P')
cusp1 = cusps[0]  # First house cusp
```

---

## Numerical model differences

LibEphemeris uses NASA JPL DE ephemerides, IAU/ERFA reductions, and independently
published models. Results are therefore not guaranteed to be bit-identical to
another ephemeris. Compatibility checks are performed ephemerally; per-date
per-date comparison values and delta tables are not stored in this project.

Planetary coordinates are validated against JPL states, house geometry against
the published spherical definitions, and frame transforms against ERFA. Speeds
are direct JPL state derivatives where available and documented numerical
derivatives otherwise.

### Ayanamshas

Only fixed-epoch modes, True Citra, and Galactic Centre 0 Sagittarius are
computed natively. Every predefined base mode 0–46 is operational without a
warn-and-degrade fallback. Review the [sidereal mode
classification](../reference/ayanamsha.md) for each mode's defining anchor;
use `SIDM_USER` for a custom epoch and offset.

### Lunar Nodes

| Component | Model | Notes |
|-----------|-------|-------|
| Mean Node | ERFA/IERS Delaunay node argument | Standards-derived mean curve |
| True Node | JPL angular-momentum geometry | Instantaneous orbital plane |

The True Node is derived from the Moon's JPL position and velocity and is
validated through independent orbital-plane identities.

### Lilith (Lunar Apogee)

| Component | Model | Notes |
|-----------|-------|-------|
| Mean Apogee (Mean Lilith) | ERFA/IERS Delaunay identities | Conventional inclined mean orbit |
| Osculating Apogee (True Lilith) | JPL eccentricity-vector geometry | Instantaneous orbit |
| Interpolated Apogee/Perigee | Separate Delaunay series and hash-pinned refinement | rc7 compatibility contract |

**Note:** True Lilith (`OSCU_APOG`, body 13) comes directly from the
eccentricity vector of the JPL DE440/DE441 geocentric lunar state; no fitted
perturbation series is added. See [Osculating Lunar Apogee](../methodology/true-lilith.md).

---

## Known Compatibility Gaps

The following features are present in pyswisseph but have limited or different implementation in LibEphemeris:

### Eclipses and other measured differences

Solar/lunar eclipse and occultation entry points are implemented. Remaining
numeric or convention differences are maintained in the evidence-backed
[Known Differences](../comparison/known-differences.md) catalog rather than a
blanket "partial" classification.
- `lun_eclipse_when()` / `lun_eclipse_when()`
- `lun_occult_when_glob()` / `lun_occult_when_glob()`

### Fixed Star Velocities

Fixed star velocities are computed from proper motion, precession, and frame
rates. Note the 3-tuple return shape:

```python
xx, star_name, retflag = ephem.fixstar_ut("Aldebaran", jd, FLG_SPEED)
# xx[3], xx[4], xx[5] carry the longitude/latitude/distance rates
```

### Date Range Limitations

LibEphemeris uses JPL DE ephemerides with specific date ranges:

| Ephemeris | Date Range | Size | Notes |
|-----------|-----------|------|-------|
| DE440s | 1849-2150 | ~31 MB | Lightweight subset of DE440 |
| DE440 (default) | 1550-2650 | ~114 MB | ICRF 3.0, recommended |
| DE441 | -13200 to 17191 | ~3.4 GB | Extended version of DE440 |

### Precision Tiers

LibEphemeris organizes the current-generation files into three precision tiers:

| Tier | File | Use Case |
|------|------|----------|
| `base` | de440s.bsp | Lightweight, modern-era usage |
| `medium` | de440.bsp | General purpose **(DEFAULT)** |
| `extended` | de441.bsp | Historical/far-future research |

```python
from libephemeris import set_precision_tier

set_precision_tier("extended")   # uses de441.bsp
```

```bash
export LIBEPHEMERIS_PRECISION=extended
```

### Selecting an Ephemeris File Directly

You can also bypass the tier system and select a specific file:

```python
from libephemeris import set_ephemeris_file, set_ephe_path

# Set custom ephemeris file
set_ephemeris_file("de441.bsp")  # Extended range version of DE440

# Or specify a custom directory
set_ephe_path("/path/to/ephemeris/files")
```

You can also select the ephemeris file via the `LIBEPHEMERIS_EPHEMERIS` environment variable:

```bash
export LIBEPHEMERIS_EPHEMERIS=de441.bsp
```

### Resolution Priority

Resolution priority (highest to lowest):

1. `LIBEPHEMERIS_EPHEMERIS` environment variable
2. `set_ephemeris_file()` / `set_jpl_file()` programmatic call
3. Precision tier (`LIBEPHEMERIS_PRECISION` env var or `set_precision_tier()`)
4. Default: `de440.bsp` (medium tier)

---

## Thread Safety with EphemerisContext

**This is a major difference from pyswisseph.**

The compatibility module-level API exposes process-global configuration state.
LibEphemeris preserves that public behavior and also offers an explicitly
thread-safe alternative through `EphemerisContext`.

### Global API (pyswisseph-compatible, NOT thread-safe)

```python
import libephemeris as swe

# This works like pyswisseph - uses global state
swe.set_topo(12.5, 41.9, 0)  # Rome
swe.set_sid_mode(swe.SIDM_TRUE_CITRA)  # independently sourced stellar mode

pos, _ = swe.calc_ut(2451545.0, swe.SUN, swe.FLG_SIDEREAL)
```

**Warning**: Do NOT use the global API in multi-threaded applications. Different threads modifying global state (topo, sidereal mode) will cause race conditions.

### Context API (thread-safe)

For multi-threaded applications, use `EphemerisContext`:

```python
from libephemeris import EphemerisContext
from libephemeris.constants import SUN, MOON, FLG_SIDEREAL

# Each thread creates its own context with isolated state
ctx = EphemerisContext()
ctx.set_topo(12.5, 41.9, 0)    # Rome - isolated to this context
ctx.set_sid_mode(swe.SIDM_TRUE_CITRA)  # isolated to this context

# Thread-safe calculations
sun, _ = ctx.calc_ut(2451545.0, SUN, FLG_SIDEREAL)
moon, _ = ctx.calc_ut(2451545.0, MOON, FLG_SIDEREAL)
cusps, ascmc = ctx.houses(2451545.0, 41.9, 12.5, ord('P'))
```

### Multi-Threading Example

```python
from libephemeris import EphemerisContext
from libephemeris.constants import SUN, MOON, FLG_SIDEREAL
from concurrent.futures import ThreadPoolExecutor

def calculate_chart(location: dict, jd: float) -> dict:
    """Thread-safe chart calculation function."""
    # Each thread creates its own context
    ctx = EphemerisContext()
    ctx.set_topo(location['lon'], location['lat'], 0)
    ctx.set_sid_mode(swe.SIDM_TRUE_CITRA)

    sun, _ = ctx.calc_ut(jd, SUN, FLG_SIDEREAL)
    moon, _ = ctx.calc_ut(jd, MOON, FLG_SIDEREAL)
    cusps, ascmc = ctx.houses(jd, location['lat'], location['lon'], ord('P'))

    return {
        'location': location['name'],
        'sun': sun[0],
        'moon': moon[0],
        'asc': ascmc[0]
    }

# Different locations
locations = [
    {'name': 'Rome', 'lon': 12.5, 'lat': 41.9},
    {'name': 'London', 'lon': -0.1, 'lat': 51.5},
    {'name': 'Tokyo', 'lon': 139.7, 'lat': 35.7},
    {'name': 'New York', 'lon': -74.0, 'lat': 40.7},
    {'name': 'Sydney', 'lon': 151.2, 'lat': -33.9},
]

# Calculate concurrently - each thread has its own isolated context
jd = 2451545.0  # J2000.0
with ThreadPoolExecutor(max_workers=5) as executor:
    results = list(executor.map(
        lambda loc: calculate_chart(loc, jd),
        locations
    ))

for r in results:
    print(f"{r['location']}: Sun={r['sun']:.2f}, Moon={r['moon']:.2f}, Asc={r['asc']:.2f}")
```

### EphemerisContext Methods

| Method | Description |
|--------|-------------|
| `__init__(ephe_path, ephe_file)` | Create context with optional ephemeris path/file |
| `set_topo(lon, lat, alt)` | Set observer location (isolated to this context) |
| `get_topo()` | Get current observer location |
| `set_sid_mode(mode, t0, ayan_t0)` | Set sidereal mode (isolated to this context) |
| `get_sid_mode(full=False)` | Get current sidereal mode |
| `calc_ut(tjd_ut, ipl, iflag)` | Calculate planetary position (UT) |
| `calc(tjd, ipl, iflag)` | Calculate planetary position (TT/ET) |
| `calc_pctr(tjd_ut, ipl, iplctr, iflag)` | Planet-centric position |
| `houses(tjd_ut, lat, lon, hsys)` | Calculate house cusps and angles |
| `close()` | Class method - close all shared ephemeris resources |

### Resource Sharing

`EphemerisContext` is designed for both thread safety and memory efficiency:

- **Isolated State**: Each context has its own observer location, sidereal mode, and angles cache
- **Shared Resources**: Expensive resources (ephemeris files ~16MB+, timescale data) are shared across all contexts
- **Thread-Safe Loading**: Uses double-checked locking pattern for lazy initialization

```python
# Multiple contexts share the same ephemeris file (memory efficient)
ctx1 = EphemerisContext()
ctx2 = EphemerisContext()
ctx3 = EphemerisContext()

# Each context has isolated state
ctx1.set_sid_mode(swe.SIDM_TRUE_CITRA)
ctx2.set_sid_mode(swe.SIDM_J2000)
ctx3.set_sid_mode(27) # True Citra

# But they all share the same ephemeris data in memory
```

---

## FLG_MOSEPH (Moshier Ephemeris Flag)

The `FLG_MOSEPH` flag is accepted for API compatibility but does not select a
Moshier semi-analytical ephemeris. Backend selection remains controlled by the
active calculation mode: LEB, Horizons, or the local Skyfield/JPL pipeline.
Code that used `FLG_MOSEPH` continues to run, but the flag does not replace the
active source with Moshier.

```python
# This still works, but FLG_MOSEPH does not select another backend:
pos, _ = swe.calc_ut(jd, swe.SUN, swe.FLG_MOSEPH | swe.FLG_SPEED)
# Equivalent to:
pos, _ = swe.calc_ut(jd, swe.SUN, swe.FLG_SPEED)
```

---

## Calculation Backend

LibEphemeris uses **JPL DE440/DE441** as its planetary data source. The default
`"auto"` mode resolves positions via LEB (if configured), Horizons API (if no
local DE kernel), or Skyfield. Compatibility flags are accepted according to
measured public behavior; no claim about the reference API's internal backend
is needed.

For performance-critical workloads, LibEphemeris also supports an optional **LEB** (LibEphemeris Binary) backend that provides ~14x faster evaluation using precomputed Chebyshev approximations. LEB is entirely opt-in and not needed for correctness.

The calculation mode controls which backend is used:

| Mode | Behavior |
|------|----------|
| `"auto"` **(default)** | Try LEB first, then Horizons API (if no local DE440), then Skyfield |
| `"skyfield"` | Always Skyfield/DE440 |
| `"leb"` | Require LEB and keep it as the only persistent ephemeris source; curated non-LEB bodies may use a traced local analytical/Keplerian model |
| `"horizons"` | Prefer Horizons API (requires internet); unsupported bodies/flags fall back to Skyfield |

```python
from libephemeris import set_calc_mode

set_calc_mode("skyfield")   # Force pure JPL/Skyfield
set_calc_mode("leb")        # Sealed LEB (never JPL/Horizons/SPK)
set_calc_mode("horizons")   # Force Horizons API
set_calc_mode("auto")       # Default: LEB -> Horizons -> Skyfield
```

The default `"auto"` mode tries the bundled/configured/local LEB first, then
Horizons API, then Skyfield. The wheel contains the reviewed base core; reviewed
medium and extended cores are available through the SHA-256-pinned downloader.
Additional modular LEB1/LEB2 files can be generated locally from the selected
JPL kernel.

See the [LEB Technical Guide](../leb/guide.md) for details.

---

## Migration Checklist

- [ ] Replace `import swisseph as swe` with `import libephemeris as swe`
- [ ] Constants use the same names as pyswisseph (`SUN`, `MOON`, `FLG_SPEED`, etc.) — no changes needed
- [ ] Check house cusp array indexing (0-based in LibEphemeris)
- [ ] Verify date range is within ephemeris coverage (1550-2650 for DE440)
- [ ] For multi-threaded apps: migrate to `EphemerisContext` API
- [ ] Update tests for relaxed tolerances on star-based ayanamshas (< 0.06 degrees)
- [ ] Review any application-specific tolerance against the current
      [precision report](../comparison/precision.md)

---

## Reporting Issues

If you encounter compatibility issues not covered in this guide, please report them at:

https://github.com/g-battaglia/libephemeris/issues

Include:
1. A minimal public-API call that demonstrates the semantic issue
2. The documented behavior you expected
3. The LibEphemeris error or behavior category
4. Python version and LibEphemeris version

Do not attach reference-distribution files or persist numerical comparison
output in the issue. Maintainers reproduce public-API comparisons ephemerally.
