# LEB (LibEphemeris Binary) — Complete Technical Guide

> **Version:** 2.2 — March 2026
> **Status:** Production-ready (LEB1 and LEB2 formats), **31 core bodies <0.001" precision (+ 31 exotic minor bodies)**
> **Source of truth:** This document. See also [Algorithms & Theory](algorithms.md) for detailed mathematical foundations.
> **Quick reference:** [Generation Quickstart](quickstart.md) — step-by-step commands for generating LEB1 and LEB2 files.
> LEB2 compressed format details: see `proposals/leb2-implementation-plan.md` and `release-notes/v1.0.0.md`.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Architecture](#2-architecture)
3. [Binary File Format](#3-binary-file-format)
4. [Reader (`leb_reader.py`)](#4-reader)
5. [Calculation Pipelines (`fast_calc.py`)](#5-calculation-pipelines)
6. [Generator (`scripts/generate_leb.py`)](#6-generator)
   - [6.15 Group Generation & Merge](#615-group-generation--merge)
7. [Integration with LibEphemeris](#7-integration-with-libephemeris)
8. [Thread Safety and `EphemerisContext`](#8-thread-safety-and-ephemeriscontext)
9. [Body Catalog](#9-body-catalog)
10. [Precision and Validation](#10-precision-and-validation)
11. [Performance](#11-performance)
12. [Commands Reference](#12-commands-reference)
13. [LEB2 Compressed Format](#13-leb2-compressed-format)
14. [Troubleshooting](#14-troubleshooting)
15. [Internals Deep-Dive](#15-internals-deep-dive)

---

## 1. Overview

### What is LEB?

LEB is a precomputed binary ephemeris format that stores Chebyshev polynomial
approximations of celestial body positions. When activated, it replaces the
default Skyfield/JPL pipeline as the primary data source for `calc_ut()`
and `calc()`, providing approximately **14x speedup** for typical
astrological chart calculations.

### Design Principles

- **Zero new runtime dependencies.** The reader uses only `mmap`, `struct`,
  `math`, and `dataclasses` (all stdlib). `numpy`, `erfa`, and `spktype21`
  are only needed at generation time.
- **Silent transparent fallback.** If a body is not in the `.leb` file, or
  the Julian Day is out of range, or an unsupported flag combination is
  requested, the system silently falls back to the full Skyfield pipeline.
  Callers never need to know whether LEB is active.
- **API-transparent.** The public API (`calc_ut`, `calc`, `calc_ut`,
  `calc`, `EphemerisContext.calc_ut`, etc.) remains unchanged. LEB is
  activated by setting a file path; no code changes are needed.
- **Immutable after init.** Once a `LEBReader` is constructed, all its data
  structures are read-only. This makes it inherently thread-safe.

### Two Data Paths

| Path | Data Source | Activation | Speed |
|------|------------|------------|-------|
| **Skyfield path** | JPL DE440/DE441 via Skyfield | Always available as fallback | ~120 µs/eval |
| **LEB path** | Precomputed `.leb` / `.leb2` file | `set_leb_file()`, auto-discovery, or auto-download | ~5-8 µs/eval |

The default `auto` mode uses the LEB path when available (including auto-download), falling back to Skyfield otherwise.

---

## 2. Architecture

### Module Map

```
libephemeris/
  leb_format.py    Format constants, struct layouts, dataclasses, serialization helpers
  leb_reader.py    LEBReader class: mmap, Clenshaw evaluation, delta-T, star catalog
  fast_calc.py     Four calculation pipelines (A/A'/B/C), flag dispatch, sidereal/ayanamsa,
                   gravitational deflection, COB corrections

scripts/
  generate_leb.py  CLI generator: Chebyshev fitting, vectorized evaluation, binary assembly

data/leb/
  ephemeris_base.leb      Base tier (de440s, 1850-2150, ~53 MB)
  ephemeris_medium.leb    Medium tier (de440, 1550-2650, ~175 MB)
  ephemeris_extended.leb  Extended tier (de441, -5000 to 5000, ~1.6 GB)

tests/test_leb/
  test_leb_format.py      Format constants and serialization
  test_leb_reader.py      Reader, Clenshaw, segment lookup
  test_fast_calc.py       Pipeline A/A'/B/C, flags, sidereal
  test_generate_leb.py    Generator functions, fitting, verification
  test_context_leb.py     EphemerisContext LEB integration
  conftest.py             Shared fixtures

tests/test_leb/compare/
  conftest.py                          Shared infrastructure (tolerances, helpers, fixtures)
  test_compare_leb_planets.py          Planets (19 test files for medium tier)
  base/                                Base tier comparison tests (8 test files)
  extended/                            Extended tier tests
  crosstier/                           Cross-tier consistency tests
```

### Data Flow

```
User calls calc_ut(jd, MARS, FLG_SPEED)
  |
  v
planets.py: calc_ut()
  |-- state.get_leb_reader() -> LEBReader or None
  |
  |-- [LEB active] fast_calc.fast_calc_ut(reader, jd, ipl, iflag)
  |     |-- Delta-T: deltat(jd) -> jd_tt  (NOT reader.delta_t)
  |     |-- Body lookup: reader._bodies[ipl] -> BodyEntry
  |     |-- Pipeline dispatch by coord_type:
  |     |     COORD_ICRS_BARY         -> _pipeline_icrs()                  [Pipeline A]
  |     |     COORD_ICRS_BARY_SYSTEM  -> _pipeline_icrs(is_system_bary=True) [Pipeline A']
  |     |     COORD_ECLIPTIC          -> _pipeline_ecliptic()              [Pipeline B]
  |     |     COORD_HELIO_ECL         -> _pipeline_helio()                 [Pipeline C]
  |     |-- Gravitational deflection (Sun, Jupiter, Saturn) [Pipeline A/A']
  |     |-- Sidereal correction (if FLG_SIDEREAL)
  |     |-- Return (lon, lat, dist, dlon, dlat, ddist), iflag
  |     |
  |     |-- [KeyError/ValueError] -> fall through to Skyfield
  |
  |-- [Skyfield fallback] _calc_body(t, ipl, iflag)
```

### Activation

```bash
# Method 0: Download pre-generated LEB file (easiest)
libephemeris download leb-base       # ~53 MB, 1850-2150
libephemeris download leb-medium     # ~175 MB, 1550-2650
# Auto-discovered from ~/.libephemeris/leb/ — no further configuration needed
```

```python
# Method 1: Programmatic
from libephemeris import set_leb_file
set_leb_file("/path/to/ephemeris.leb")   # enable
set_leb_file(None)                        # disable

# Method 2: Environment variable
export LIBEPHEMERIS_LEB=/path/to/ephemeris.leb

# Method 3: Per-context (thread-safe)
ctx = EphemerisContext()
ctx.set_leb_file("/path/to/ephemeris.leb")
```

**Resolution priority** (highest to lowest):
1. `EphemerisContext._leb_file` (per-context)
2. Global `set_leb_file()` call
3. `LIBEPHEMERIS_LEB` environment variable
4. Auto-discovery: `~/.libephemeris/leb/{tier}_core.leb2` (LEB2) or `~/.libephemeris/leb/ephemeris_{tier}.leb` (LEB1)

### Calculation Mode (`LIBEPHEMERIS_MODE`)

The calculation mode controls how `calc_ut()` and `calc()` resolve
positions. Set via `set_calc_mode()` or the `LIBEPHEMERIS_MODE` environment
variable.

| Mode | Behavior |
|------|----------|
| `auto` (default) | LEB if configured, then Horizons (if no local DE440), then Skyfield |
| `skyfield` | Always use Skyfield, even if a `.leb` file is configured |
| `leb` | Require LEB (auto-discovered or auto-downloaded if needed); unsupported bodies/flags fall back to Skyfield |
| `horizons` | Prefer Horizons API; unsupported bodies/flags fall back to Skyfield |

```python
from libephemeris import set_calc_mode, get_calc_mode

set_calc_mode("skyfield")  # Force Skyfield for benchmarking/validation
set_calc_mode("leb")       # Require LEB (error if unavailable)
set_calc_mode("horizons")  # Prefer Horizons API
set_calc_mode("auto")      # Default behavior
set_calc_mode(None)        # Reset to env var / default
```

```bash
export LIBEPHEMERIS_MODE=skyfield   # Force Skyfield
export LIBEPHEMERIS_MODE=leb        # Require LEB
export LIBEPHEMERIS_MODE=horizons   # Prefer Horizons API
export LIBEPHEMERIS_MODE=auto       # Default (same as unset)
```

In `leb` mode, bodies not present in the `.leb` file still fall through to
Skyfield (the mode only requires that a valid `.leb` file is loaded, not
that every body is in it).

---

## 3. Binary File Format

**Source file:** `libephemeris/leb_format.py`

### Magic and Version

```
Magic:   b"LEB1" (4 bytes)
Version: 1 (uint32)
```

### Overall Layout

```
Offset      Content                          Size
────────────────────────────────────────────────────
0x0000      File Header                      64 bytes
0x0040      Section Directory                N × 24 bytes (N=5 currently)
variable    Section 0: Body Index            body_count × 52 bytes
variable    Section 1: Chebyshev Data        variable (bulk of the file)
variable    Section 2: Nutation Data         header(40B) + segments
variable    Section 3: Delta-T Table         header(8B) + entries(16B each)
variable    Section 4: Star Catalog          n_stars × 64 bytes
(reserved)  Section 5: Orbital Elements      not yet used
```

### 3.1 File Header (64 bytes)

```c
// Struct format: "<4sIIIdddI20s" (little-endian)
struct FileHeader {
    char     magic[4];           // b"LEB1"
    uint32_t version;            // 1
    uint32_t section_count;      // Number of sections (currently 5)
    uint32_t body_count;         // Number of bodies in the file
    float64  jd_start;           // Start of date coverage (JD TT)
    float64  jd_end;             // End of date coverage (JD TT)
    float64  generation_epoch;   // JD when file was generated
    uint32_t flags;              // Reserved (0)
    char     reserved[20];       // Padding to 64 bytes
};
```

**Python dataclass:** `leb_format.FileHeader`
**Constant:** `HEADER_SIZE = 64`

### 3.2 Section Directory Entry (24 bytes)

```c
// Struct format: "<IIQQ"
struct SectionEntry {
    uint32_t section_id;    // SECTION_BODY_INDEX=0, ..., SECTION_STARS=4
    uint32_t reserved;      // Padding
    uint64_t offset;        // Absolute byte offset from file start
    uint64_t size;          // Section size in bytes
};
```

**Python dataclass:** `leb_format.SectionEntry`
**Constant:** `SECTION_DIR_SIZE = 24`

### 3.3 Body Index Entry (52 bytes)

```c
// Struct format: "<iIIdddIIQ"
struct BodyEntry {
    int32_t  body_id;        // SE_* constant (0=Sun, 1=Moon, ...)
    uint32_t coord_type;     // 0=ICRS_BARY, 1=ECLIPTIC, 2=HELIO_ECL
    uint32_t segment_count;  // Number of Chebyshev segments
    float64  jd_start;       // Body coverage start (JD TT)
    float64  jd_end;         // Body coverage end (JD TT)
    float64  interval_days;  // Segment width in days
    uint32_t degree;         // Chebyshev polynomial degree
    uint32_t components;     // Always 3 (x/y/z or lon/lat/dist)
    uint64_t data_offset;    // Absolute byte offset to first coefficient
};
```

**Python dataclass:** `leb_format.BodyEntry`
**Constant:** `BODY_ENTRY_SIZE = 52`

### 3.4 Chebyshev Coefficient Storage

Coefficients are stored as contiguous `float64` arrays in **component-major**
order. For a body with `degree=13` and `components=3`:

```
Segment layout (3 × 14 × 8 = 336 bytes per segment):
  [c0_x, c1_x, ..., c13_x,     ← 14 coefficients for component 0 (x or lon)
   c0_y, c1_y, ..., c13_y,     ← 14 coefficients for component 1 (y or lat)
   c0_z, c1_z, ..., c13_z]     ← 14 coefficients for component 2 (z or dist)
```

The byte size of one segment is computed by:
```python
segment_byte_size(degree, components) = components * (degree + 1) * 8
```

### 3.5 Nutation Section

**Header (40 bytes):**
```c
// Struct format: "<dddIIII"
struct NutationHeader {
    float64  jd_start;
    float64  jd_end;
    float64  interval_days;   // 16.0 days
    uint32_t degree;          // 16
    uint32_t components;      // 2 (dpsi, deps in radians)
    uint32_t segment_count;
    uint32_t reserved;
};
```

Followed by `segment_count` segments, each containing `2 × 17 × 8 = 272 bytes`
of Chebyshev coefficients for dpsi and deps (IAU 2006/2000A nutation).

### 3.6 Delta-T Section

**Header (8 bytes):**
```c
// Struct format: "<II"
struct DeltaTHeader {
    uint32_t n_entries;
    uint32_t reserved;
};
```

**Entries (16 bytes each):**
```c
// Struct format: "<dd"
struct DeltaTEntry {
    float64 jd;          // Julian Day
    float64 delta_t;     // TT - UT1 in days
};
```

Sampled every 30 days. The reader uses linear interpolation between entries.

### 3.7 Star Catalog Entry (64 bytes)

```c
// Struct format: "<iddddddd4s"
struct StarEntry {
    int32_t  star_id;
    float64  ra_j2000;     // Right ascension (degrees)
    float64  dec_j2000;    // Declination (degrees)
    float64  pm_ra;        // Proper motion RA (deg/yr, includes cos(dec))
    float64  pm_dec;       // Proper motion Dec (deg/yr)
    float64  parallax;     // Parallax (arcsec)
    float64  rv;           // Radial velocity (km/s)
    float64  magnitude;    // Visual magnitude
    char     reserved[4];
};
```

### 3.8 Coordinate Types

| Value | Constant | Meaning | Bodies |
|-------|----------|---------|--------|
| 0 | `COORD_ICRS_BARY` | ICRS barycentric planet center (x, y, z) in AU | Sun, Moon, Mercury, Venus, Mars, Earth, Chiron, Ceres-Vesta |
| 1 | `COORD_ECLIPTIC` | Ecliptic of date (lon, lat, dist) in deg/deg/AU | Mean/True Node, Mean/Oscu Apogee, Interp Apogee/Perigee |
| 2 | `COORD_HELIO_ECL` | Heliocentric ecliptic (lon, lat, dist) | Cupido-Poseidon, Transpluto |
| 3 | `COORD_GEO_ECLIPTIC` | Geocentric ecliptic of date — **reserved, not used** | None |
| 4 | `COORD_ICRS_BARY_SYSTEM` | ICRS system barycenter (x, y, z) in AU, COB at runtime | Jupiter, Saturn, Uranus, Neptune, Pluto |

**Why `COORD_ICRS_BARY_SYSTEM`?** Outer planets (Jupiter-Pluto) have moons whose
gravitational influence creates high-frequency oscillations in the planet center
position relative to the system barycenter (the Center-of-Body or COB correction).
These oscillations are difficult to fit with Chebyshev polynomials — they require
very short intervals and high degree, producing large files and residual fitting
errors. The solution is to store the smooth system barycenter in the `.leb` file
and apply the COB correction at runtime. The COB correction uses either
`planet_centers.bsp` (SPK segments, <0.001" precision) or analytical moon theory
corrections as fallback (<0.01"). See [Algorithms & Theory](algorithms.md) for details.

---

## 4. Reader

**Source file:** `libephemeris/leb_reader.py`

### 4.1 LEBReader Class

```python
class LEBReader:
    def __init__(self, path: str) -> None
    def __enter__(self) -> "LEBReader"
    def __exit__(self, *args) -> None
    def has_body(self, body_id: int) -> bool
    def eval_body(self, body_id: int, jd: float) -> ((pos), (vel))
    def eval_nutation(self, jd_tt: float) -> (dpsi, deps)
    def delta_t(self, jd: float) -> float
    def get_star(self, star_id: int) -> StarEntry
    def warm(self, jd_start: float, jd_end: float) -> None
    def close(self) -> None

    # Properties
    path -> str
    jd_range -> (jd_start, jd_end)
```

### 4.2 Memory Mapping

The file is opened read-only via `mmap.mmap(fd, 0, access=mmap.ACCESS_READ)`.
All struct reads use `struct.unpack_from()` which operates directly on the
mmap'd buffer — **zero-copy** for coefficient access.

Pages are loaded **on demand** by the kernel (demand paging).  No
`madvise(MADV_WILLNEED)` is issued at open time, so the mmap does not
pre-fault pages into physical RAM.  This keeps idle RSS near zero
regardless of file size.

**Selective preloading:** Both `LEBReader` and `LEB2Reader` expose a
`warm(jd_start, jd_end)` method that calls `madvise(MADV_WILLNEED)` on
the page-aligned byte ranges of segments/chunks overlapping the given
Julian Day range.  This allows pre-faulting only the date ranges that
will be needed (e.g. 1800-2200 CE ≈ 11 MB for extended tier) without
loading the entire file (~1.6 GB for the extended tier).  `CompositeLEBReader.warm()` delegates
to all constituent readers.

Automatic preloading can be enabled via TOML configuration:

```toml
[libephemeris]
mmap_preload = true          # default: false
mmap_preload_start = 1800    # start year
mmap_preload_end = 2200      # end year
```

**Page cache release (containerised environments):** In cgroup v2
environments (Docker, Railway, Kubernetes), file-backed page cache is
counted in `memory.current`.  All readers expose a `cool()` method
(complement of `warm()`) that calls `madvise(MADV_DONTNEED)` to advise
the kernel that mmap pages can be reclaimed from the cgroup counter.
The `release_data_cache()` function does the same for all files in the
data directory via `posix_fadvise(FADV_DONTNEED)`.

```python
import libephemeris

reader = libephemeris.get_leb_reader()
if reader:
    reader.cool()                    # release mmap pages
libephemeris.release_data_cache()    # release data file pages
```

Both are advisory hints — the kernel ignores them on desktop systems
with available RAM.  `cool()` is idempotent and safe on closed readers.
It does not clear Python-level caches (`_eval_cache`, `_chunk_cache`).

**Resource safety:** The `__init__` wraps `_parse()` in try/except. If parsing
fails, `close()` is called to release the mmap and file handle before
re-raising the exception (`leb_reader.py:181-186`).

### 4.3 Parsing (at construction time)

1. **Header** — Validates magic (`b"LEB1"`) and version (1).
2. **Section directory** — Builds `_sections: Dict[int, SectionEntry]`.
3. **Body index** — Builds `_bodies: Dict[int, BodyEntry]`.
4. **Nutation header** — Stores `_nutation: NutationHeader` and data offset.
5. **Delta-T table** — Loads into two Python lists: `_delta_t_jds` and
   `_delta_t_vals`.
6. **Star catalog** — Builds `_stars: Dict[int, StarEntry]`.

### 4.4 Body Evaluation (`eval_body`)

The core evaluation path (the hot path, ~1.5 us per call):

```python
def eval_body(self, body_id: int, jd: float):
    body = self._bodies[body_id]

    # 1. O(1) segment lookup (no binary search needed)
    seg_idx = int((jd - body.jd_start) / body.interval_days)
    seg_idx = clamp(seg_idx, 0, body.segment_count - 1)

    # 2. Map JD to normalized tau in [-1, 1]
    seg_start = body.jd_start + seg_idx * body.interval_days
    seg_mid = seg_start + 0.5 * body.interval_days
    tau = 2.0 * (jd - seg_mid) / body.interval_days
    tau = clamp(tau, -1.0, 1.0)

    # 3. Zero-copy coefficient read from mmap
    byte_offset = body.data_offset + seg_idx * n_coeffs * 8
    coeffs = struct.unpack_from(f"<{n_coeffs}d", self._mm, byte_offset)

    # 4. Clenshaw evaluation per component (value + derivative)
    for c in range(3):
        comp_coeffs = coeffs[c * deg1 : (c + 1) * deg1]
        val, deriv = _clenshaw_with_derivative(comp_coeffs, tau)
        pos.append(val)
        vel.append(deriv * 2.0 / body.interval_days)

    # 5. Longitude wrapping for ecliptic bodies
    if coord_type in (COORD_ECLIPTIC, COORD_HELIO_ECL):
        pos[0] = pos[0] % 360.0

    return tuple(pos), tuple(vel)
```

**Key design choices:**

- **O(1) segment lookup** via integer division (not binary search). This is
  possible because all segments for a given body have equal width
  (`interval_days`).

- **Clenshaw algorithm** for Chebyshev evaluation (not Horner). Clenshaw is
  numerically stable for Chebyshev polynomials and requires only 2 temporary
  variables. Implementation at `leb_reader.py:58-79`.

- **Analytical derivative** via the Chebyshev derivative recurrence relation
  (`_deriv_coeffs` at `leb_reader.py:82-112`). The recurrence computes
  derivative coefficients d_k from the original coefficients c_k:
  ```
  d_{n-1} = 2n * c_n
  d_k     = d_{k+2} + 2(k+1) * c_{k+1}    for k = n-2, ..., 1
  d_0     = d_2/2 + c_1
  ```
  Then evaluates the derivative polynomial via Clenshaw. The derivative is
  in d/d(tau) units; it is scaled by `2 / interval_days` to get d/d(jd).

### 4.5 Delta-T Interpolation

```python
def delta_t(self, jd: float) -> float:
    # Binary search (bisect_right) for the interval
    idx = bisect_right(jds, jd) - 1
    # Linear interpolation (sufficient for 30-day spacing)
    t = (jd - jds[idx]) / (jds[idx+1] - jds[idx])
    return vals[idx] + t * (vals[idx+1] - vals[idx])
```

Values are clamped at boundaries (returns first/last value for out-of-range JDs).

### 4.6 Nutation Evaluation

Same pattern as body evaluation: O(1) segment lookup, tau mapping, Clenshaw
on 2 components (dpsi, deps). Returns values in **radians**.

---

## 5. Calculation Pipelines

**Source file:** `libephemeris/fast_calc.py`

### 5.1 Entry Points

```python
def fast_calc_ut(reader, tjd_ut, ipl, iflag, *,
                 sid_mode=None, sid_t0=None, sid_ayan_t0=None)
    -> ((lon, lat, dist, dlon, dlat, ddist), iflag)

def fast_calc_tt(reader, tjd_tt, ipl, iflag, *,
                 sid_mode=None, sid_t0=None, sid_ayan_t0=None)
    -> ((lon, lat, dist, dlon, dlat, ddist), iflag)
```

Both functions:
1. Check for unsupported flags and raise `KeyError` to trigger fallback.
2. Snapshot sidereal state at entry (thread-safe).
3. Convert UT to TT via `deltat()` from `time_utils` (high-precision
   Skyfield model, **not** `reader.delta_t()` which uses linear interpolation
   on a sparse table and can introduce up to ~0.004s error near 1985).
   The `fast_calc_tt` path uses `reader.delta_t()` only for the reverse
   TT→UT approximation needed by sidereal ayanamsa.
4. Delegate to `_fast_calc_core()`.

### 5.2 Flag Handling

**Flags that trigger immediate Skyfield fallback (raise KeyError):**

| Flag | Reason |
|------|--------|
| `FLG_ICRS` | ICRS output frame not implemented on the LEB path |

**Flags handled natively by LEB:**

| Flag | Effect |
|------|--------|
| `FLG_SPEED` | Compute velocities (analytical Chebyshev derivative for all pipelines) |
| `FLG_HELCTR` | Use Sun as observer instead of Earth; skip aberration |
| `FLG_BARYCTR` | Use SSB as observer; skip aberration |
| `FLG_TRUEPOS` | Skip light-time correction and aberration |
| `FLG_NOABERR` | Skip aberration only |
| `FLG_EQUATORIAL` | Output in equatorial coordinates instead of ecliptic |
| `FLG_J2000` | Output in J2000 frame instead of of-date |
| `FLG_SIDEREAL` | Apply ayanamsa correction |
| `FLG_TOPOCTR` | Topocentric offset from `set_topo()` / context topo |
| `FLG_XYZ` | Cartesian output |
| `FLG_RADIANS` | Radian output |
| `FLG_NONUT` | No-nutation frame (mean ecliptic/equator of date) |
| `FLG_MOSEPH` | Silently stripped (always uses JPL data) |

Individual body classes may still fall back per body when the loaded
file lacks the data they need.

### 5.3 Pipeline A: ICRS Barycentric Bodies (Planet Centers)

**Bodies:** Sun, Moon, Mercury, Venus, Mars, Earth, Chiron, Ceres, Pallas,
Juno, Vesta (11 bodies)

**Stored as:** (x, y, z) in AU, ICRS barycentric frame (`COORD_ICRS_BARY`)

**Algorithm (`_pipeline_icrs`, line 681):**

1. **Body position:** `reader.eval_body(ipl, jd_tt)` -> (x, y, z) in AU
2. **Observer selection:**
   - Geocentric (default): Earth position from LEB
   - Heliocentric (`FLG_HELCTR`): Sun position from LEB
   - Barycentric (`FLG_BARYCTR`): origin (0, 0, 0)
3. **Geometric vector:** target - observer
4. **Light-time correction** (unless `FLG_TRUEPOS`):
   - 3 fixed-point iterations
   - `lt = dist / C_LIGHT_AU_DAY` (173.14 AU/day)
   - Re-evaluate body at `jd_tt - lt`
5. **Gravitational deflection** (unless Moon, helio, bary, truepos, or noaberr):
   - PPN formula matching Skyfield's `apparent(deflectors=(10, 599, 699))`
   - Three deflectors: Sun (mass ratio 1.0), Jupiter (1047.3), Saturn (3497.9)
   - Evaluates deflector positions at both observation time and closest-approach time
   - See [Algorithms & Theory](algorithms.md#9-gravitational-deflection) for details
6. **Aberration** (unless disabled or helio/bary/truepos):
   - Classical first-order formula using Earth velocity
   - `u' = u + v/c - u*(u.v/c)`, renormalize
7. **Coordinate transform** (one of four paths):
   - **True ecliptic of date** (default, most common):
     - Precess ICRS -> true equatorial of date via the Vondrák 2011
       long-term precession + LEB-stored nutation (`vondrak_pn_matrix`,
       see `precession_vondrak.py`)
     - Rotate equatorial -> ecliptic using true obliquity (mean + nutation deps)
   - **J2000 ecliptic** (`FLG_J2000`):
     - Rotate ICRS -> the reference J2000 ecliptic frame: frame bias, then
       the IAU 2006 J2000 obliquity 84381.406"
       (`_rotate_icrs_to_ecliptic_j2000`)
   - **True equatorial of date** (`FLG_EQUATORIAL`):
     - Precess ICRS -> equatorial of date via precession-nutation matrix
   - **J2000 equatorial** (`FLG_EQUATORIAL | FLG_J2000`):
     - Rotate ICRS -> mean equator/equinox of J2000 (frame bias only),
       then convert to spherical

**Velocity** is the exact time-derivative of the *reported apparent* position
in the requested frame, so every speed slot differentiates its own position
slot. It is built as a central difference of the analytically-extrapolated
apparent place, which keeps the geometry exact while capturing all the rate
terms a raw Chebyshev-derivative rotation would drop:

1. Geometric velocity from the exact Chebyshev derivative (`eval_body()`),
   minus the observer velocity, with the light-time *rate* factor
   `1 + d(lt)/dt` applied (not just the velocity at the retarded time).
2. The barycenter-to-centre offset rate `d(COB)/dt` for the outer planets
   (evaluated at the light-time retarded epoch for heliocentric/barycentric
   output, matching the position).
3. The apparent corrections — annual aberration, gravitational deflection
   and the precession-nutation frame rotation — carried through by a central
   difference of the apparent vector at `jd_tt +/- h`, with the geometric
   state extrapolated analytically (so no Chebyshev-fit error leaks into the
   difference) and the deflectors evaluated at each sample epoch. A deflector
   is skipped when the observer sits at that body (the far-field formula is
   singular there).
4. Frame handling matches the position: J2000/sidereal speeds pick up the
   equinox-drift and nutation-rate terms; sidereal-of-date adds the true
   ayanamsha rate.

An earlier revision computed velocity as a plain analytical Chebyshev
derivative rotated like the position, dropping the aberration, deflection,
light-time and frame-rotation *rates*; that biased the geocentric Moon speed
by ~4"/day and heliocentric/planetocentric speeds by up to a few tenths of an
arcsec/day. The current scheme keeps the analytic geometry (no extra full
pipeline runs for the geometric part) while making the reported speed the
true derivative of the reported position (residual <= 0.001"/day geocentric).

### 5.3.1 Pipeline A': ICRS System Barycenter Bodies (Runtime COB)

**Bodies:** Jupiter, Saturn, Uranus, Neptune, Pluto (5 bodies)

**Stored as:** (x, y, z) in AU, system barycenter in ICRS (`COORD_ICRS_BARY_SYSTEM`)

**Why a separate pipeline?** Outer planets have moons whose gravitational
influence creates high-frequency oscillations in the planet center position.
These oscillations cannot be accurately captured by Chebyshev polynomials
without impractically short intervals. The solution stores the smooth system
barycenter and applies the Center-of-Body (COB) correction at runtime.

**Algorithm (`_pipeline_icrs` with `is_system_bary=True`):**

Same as Pipeline A, with one critical addition:

- **COB correction** is applied via `_apply_cob_correction()` at **observer time**
  (`jd_tt`), not at retarded time (`jd_tt - lt`). This matches Skyfield's
  `_SpkCenterTarget._observe_from_bcrs()` behavior.
- The light-time iteration operates on the raw barycenter position (smooth),
  then COB is applied after convergence.
- COB uses `planet_centers.bsp` (SPK segments for NAIF IDs 599, 699, 799, 899, 999)
  when available (<0.001" precision), with analytical moon theory corrections
  as fallback (<0.01").

**Architectural limitation:** For nearby asteroids (Ceres, Pallas, Juno,
Vesta), the ICRS->ecliptic pipeline amplifies errors by `1/geocentric_distance`.
This produces latitude velocity errors of 0.19-0.71 deg/day -- an inherent
property of the coordinate transformation, not the Chebyshev fitting.

### 5.4 Pipeline B: Ecliptic Direct Bodies

**Bodies:** Mean Node, True Node, Mean Apogee, Osculating Apogee,
Interpolated Apogee, Interpolated Perigee (6 bodies)

**Stored as:** (lon, lat, dist) in degrees/degrees/AU, ecliptic of date

**Algorithm (`_pipeline_ecliptic`, line 508):**

1. **Direct read:** `reader.eval_body(ipl, jd_tt)` -> (lon, lat, dist) and
   (dlon, dlat, ddist) from Chebyshev analytical derivative
2. **Coordinate transforms** (if requested):
   - **True equatorial of date** (`FLG_EQUATORIAL`):
     - Rotate ecliptic -> equatorial using true obliquity
     - Velocity via finite difference on the ecliptic coords (dt=0.001 days)
   - **J2000 equatorial** (`FLG_EQUATORIAL | FLG_J2000`):
     - Precess ecliptic-of-date -> J2000 ecliptic, then rotate to J2000 equatorial
     - Velocity via finite difference
   - **J2000 ecliptic** (`FLG_J2000`):
     - Precess ecliptic-of-date -> J2000 ecliptic
     - Velocity via finite difference

**No light-time or aberration** is applied to these bodies (they are
Earth-relative analytical quantities).

### 5.5 Pipeline C: Heliocentric Bodies

**Bodies:** Cupido, Hades, Zeus, Kronos, Apollon, Admetos, Vulkanus,
Poseidon, Transpluto (9 bodies)

**Stored as:** (lon, lat, dist) in degrees/degrees/AU, heliocentric J2000
ecliptic

**Algorithm:** `_pipeline_helio()`, a distinct pipeline. For `FLG_HELCTR`
requests the stored heliocentric coordinates are evaluated directly; for the
default geocentric output the pipeline performs a real heliocentric ->
geocentric conversion — it subtracts the Earth's position (with the
light-time correction, matching the Skyfield path) and then applies the
same of-date frame transforms as Pipeline B.

### 5.6 Sidereal Correction

Applied after the pipeline, for any body, when `FLG_SIDEREAL` is set and
the output is ecliptic (not equatorial):

```python
aya = _calc_ayanamsa_from_leb(reader, jd_tt, sid_mode, sid_t0, sid_ayan_t0)
lon = (lon - aya) % 360.0

# Speed correction: subtract IAU 2006 general precession rate
T = (jd_tt - J2000) / 36525.0
prec_rate = (5028.796195 + 2 * 1.1054348 * T) / (3600.0 * 36525.0)  # deg/day
dlon -= prec_rate
```

**Ayanamsa computation** (`_calc_ayanamsa_from_leb`, line 335):

- Supports **29 formula-based sidereal modes** (stored in `_AYANAMSHA_J2000` dict)
- Supports `SIDM_USER` (mode 255) with custom t0/ayan_t0
- **Star-based modes** (17, 27-36, 39-40, 42) raise `KeyError` to trigger
  Skyfield fallback
- True ayanamsa = mean ayanamsa + nutation in longitude (from LEB nutation data)

### 5.7 Utility Functions

| Function | Purpose | Location |
|----------|---------|----------|
| `_mean_obliquity_iau2006(jd_tt)` | IAU 2006 mean obliquity polynomial | line 68 |
| `_vec3_sub(a, b)` | 3-vector subtraction | line 86 |
| `_vec3_dist(v)` | Euclidean distance | line 91 |
| `_cartesian_to_spherical(x, y, z)` | Cartesian -> (lon, lat, dist) deg | line 96 |
| `_rotate_equatorial_to_ecliptic(x, y, z, eps)` | Frame rotation | line 110 |
| `_rotate_icrs_to_ecliptic_j2000(x, y, z)` | ICRS -> ecliptic J2000 | line 131 |
| `_apply_aberration(geo, earth_vel)` | Classical aberration formula | line 138 |
| `_get_leb_frame_data(reader, jd_tt)` | Vondrák 2011 PN matrix + LEB nutation | fast_calc.py |
| `_mat3_vec3(mat, vec)` | 3x3 matrix * 3-vector | line 221 |
| `_cotrans(lon, lat, eps)` | Ecliptic <-> equatorial spherical | line 233 |
| `_precess_ecliptic(lon, lat, from_jd, to_jd)` | Ecliptic precession | line 270 |
| `_apply_cob_correction(pos, ipl, jd_tt)` | Center-of-Body correction for outer planets | line 286 |
| `_apply_gravitational_deflection(...)` | PPN gravitational deflection (Sun/Jupiter/Saturn) | line 348 |
| `_cartesian_velocity_to_spherical(...)` | Cartesian velocity -> spherical | line 448 |

---

## 6. Generator

**Source file:** `scripts/generate_leb.py`

### 6.1 Overview

The generator creates `.leb` files by:
1. Evaluating celestial body positions at Chebyshev nodes
2. Fitting Chebyshev polynomials to the sampled values
3. Verifying the fit against intermediate test points
4. Assembling all data into the binary format

### 6.2 CLI Usage

```bash
# Tier-based (recommended)
python scripts/generate_leb.py --tier base --verify
python scripts/generate_leb.py --tier medium --verify
python scripts/generate_leb.py --tier extended --verify

# Custom range
python scripts/generate_leb.py --output custom.leb --start 1900 --end 2100

# Options
  --tier {base,medium,extended}   Preset configuration
  --output PATH                   Output file path
  --start YEAR                    Start year
  --end YEAR                      End year
  --workers N                     Parallel workers (default: CPU count)
  --verify                        Run post-generation validation
  --verify-samples N              Samples per body for verification (default: 500)
  --bodies 0,1,2                  Comma-separated body IDs (default: all)
  --quiet                         Suppress progress output
```

### 6.3 Tier Configurations

| Tier | Ephemeris | Years | Output | Approx Size |
|------|-----------|-------|--------|-------------|
| `base` | de440s.bsp | 1850-2150 | `ephemeris_base.leb` | ~53 MB |
| `medium` | de440.bsp | 1550-2650 | `ephemeris_medium.leb` | ~175 MB |
| `extended` | de441.bsp | -5000 to 5000 | `ephemeris_extended.leb` | ~1.6 GB |

### 6.4 Chebyshev Fitting

**Node computation** (`chebyshev_nodes`, line 257):
```python
# Type I Chebyshev nodes on [-1, 1]
nodes = cos(pi * (arange(n) + 0.5) / n)
```

**Segment fitting** (`fit_segment`, line 262):
1. Map Chebyshev nodes from [-1, 1] to [jd_start, jd_end]
2. Evaluate body function at each node
3. Fit using `numpy.polynomial.chebyshev.chebfit(nodes, values, degree)`

**Verification** (`verify_segment`, line 297):
- 10 uniform test points per segment (not on Chebyshev nodes)
- Compare fitted polynomial vs reference function
- Track maximum error across all components

### 6.5 Vectorized Evaluation (Key Optimization)

The major performance optimization is **batching all JDs across all segments
into a single Skyfield evaluation call**:

**Step 1 — Precompute all JDs** (`_compute_all_segment_jds`, line 333):
```
For each segment i:
  - (degree+1) Chebyshev fit nodes
  - N_VERIFY (10) uniform verification points
Total: n_segments * (degree + 1 + 10) JDs
```

**Step 2 — Single vectorized Skyfield call** (`_eval_body_icrs_vectorized`, line 598):
```python
# For inner planets (Sun, Moon, Mercury, Venus, Earth):
target = planets[target_name]
t_arr = ts.tt_jd(all_jds)                    # vectorized Time
positions = target.at(t_arr).position.au.T    # (N, 3) in one call

# For outer planets (Mars, Jupiter, Saturn, Uranus, Neptune, Pluto):
# Barycenter + SPK center offset or COB correction
bary_vals = barycenter.at(t_arr).position.au.T
offset_vals = center_segment.at(t_arr).position.au.T
positions = bary_vals + offset_vals
```

**Step 3 — Batch fit and verify** (`_fit_and_verify_from_values`, line 374):
```
For each segment:
  - Extract pre-evaluated values at Chebyshev nodes
  - chebfit() per component
  - Compare pre-evaluated verification values against fitted polynomial
```

This eliminates the per-JD overhead of Skyfield's time conversion, SPK
evaluation, and Python function call overhead. Speedup: **~150x for planets**.

### 6.6 SPK Boundary Overshoot

**Problem:** The last Chebyshev segment may extend its fit nodes beyond
`jd_end`. For example, Uranus with 128-day intervals: if the last segment
starts at JD X, its nodes extend to X+128. If de440s.bsp ends at JD 2506352.5
(~2150-01-22), nodes can overshoot by up to 128 days.

**Solution: Linear extrapolation** (`_eval_target_vectorized`, line 517):
- Identify which JDs are outside the SPK valid range (with 1-day safety margin)
- For in-range JDs: vectorized Skyfield evaluation
- For out-of-range JDs: evaluate position + velocity at the boundary, then
  extrapolate: `pos(jd) = pos(boundary) + vel(boundary) * (jd - boundary)`

**Why not clamp?** Clamping (replacing out-of-range JDs with the boundary
value) would move the Chebyshev nodes, corrupting the polynomial fit for
the entire segment. Linear extrapolation preserves the node positions and
produces smooth, reasonable values for the short extrapolation distance
(up to ~128 days).

### 6.7 Asteroid Generation

**Single path** (`generate_body_icrs_asteroid`, line ~1329):

Uses `spktype21` exclusively (~36x faster than old scalar path):
1. Open the SPK type 21 file via `spktype21.SPKType21.open()`
2. Verify it covers the requested date range (raises `RuntimeError` if not)
3. Get Sun barycentric positions via vectorized Skyfield
4. For each JD (scalar loop, ~65 us/eval):
   ```python
   pos_km, _ = kernel.compute_type21(center_id, target_id, jd)
   pos_au = pos_km / 149597870.7  # km -> AU
   bary_pos = pos_au + sun_bary   # helio -> SSB barycentric
   ```
5. Fit and verify from the collected values

**No Keplerian fallback.** If the SPK file is unavailable or doesn't cover
the requested range, the function raises `RuntimeError`. The generator
excludes asteroids without SPK rather than producing inaccurate data.

### 6.8 Per-Body Date Ranges for Asteroids

**Problem:** JPL Horizons SPK files for minor bodies cover approximately
1600-2500 CE. For the base tier (1850-2150), this is sufficient. For the
medium tier (1550-2650), SPK coverage may be slightly narrower at the edges.
For the extended tier (-5000 to 5000), asteroids cannot cover the full range.

Previously, asteroids were excluded entirely if their SPK didn't cover the
full tier range. Now, the generator uses **per-body date ranges**: each
asteroid covers only its actual SPK range, while planets and analytical
bodies cover the full tier range.

**How it works** (`assemble_leb`, step 0):

1. For each asteroid, attempt to download SPK for the full tier range
2. If the SPK covers the full range → use global `jd_start`/`jd_end`
3. If the SPK is narrower → discover actual range via `_get_asteroid_spk_range()`
4. Intersect SPK range with tier range (with 1-day safety margin)
5. If overlap is less than 20 years → exclude the asteroid
6. Otherwise → generate with the per-body range and write it to `BodyEntry`

```python
# Per-body range discovery:
spk_range = _get_asteroid_spk_range(spk_file, bid)
eff_start = max(spk_jd_start, jd_start) + 1.0  # safety margin
eff_end = min(spk_jd_end, jd_end) - 1.0
body_jd_ranges[bid] = (eff_start, eff_end)
```

**At runtime:** when a body's JD is outside its per-body range, `eval_body()`
raises `ValueError`, which `calc_ut()`/`calc()` catches to fall
through to Skyfield. This is transparent to the caller.

**Coverage by tier:**

| Tier | Planets | Asteroids | Analytical |
|------|---------|-----------|------------|
| Base (1850-2150) | Full | Full (SPK covers) | Full |
| Medium (1550-2650) | Full | ~1600-2500 (per-body) | Full |
| Extended (-5000 to 5000) | Full | ~1600-2500 (per-body) | Full |

### 6.9 Ecliptic Body Generation

**Longitude unwrapping** (`_generate_segments_unwrap`, line 1041):

Ecliptic bodies (lunar nodes, Lilith) have longitude that can wrap around
0/360 degrees. Fitting a polynomial across a 359->1 discontinuity would
produce wildly wrong results.

Solution:
1. Evaluate at Chebyshev nodes
2. `numpy.unwrap(radians(lon))` -> removes 2pi jumps
3. Convert back to degrees and fit
4. Verification re-wraps with `% 360`

### 6.10 Analytical Body Functions

| Body ID | Function | Source |
|---------|----------|--------|
| 10 (Mean Node) | `calc_mean_lunar_node(jd)` | `lunar.py` |
| 11 (True Node) | `calc_true_lunar_node(jd)` | `lunar.py` |
| 12 (Mean Apogee) | `calc_mean_lilith_with_latitude(jd)` | `lunar.py` |
| 13 (Oscu Apogee) | `calc_true_lilith(jd)` | `lunar.py` |
| 21 (Interp Apogee) | `calc_interpolated_apogee(jd)` | `lunar.py` |
| 22 (Interp Perigee) | `calc_interpolated_perigee(jd)` | `lunar.py` |
| 40-47 (Uranians) | `calc_uranian_planet(body_id, jd)` | `hypothetical.py` |
| 48 (Transpluto) | `calc_transpluto(jd)` | `hypothetical.py` |

These are pure-Python analytical functions (~315 us/eval). They are evaluated
**sequentially** — one body at a time with per-body progress bars.

### 6.11 Nutation Generation

Uses vectorized `erfa.nut06a()` (IAU 2006/2000A):
```python
jd1 = np.full_like(all_jds, 2451545.0)  # J2000 epoch
jd2 = all_jds - 2451545.0               # offset
dpsi, deps = erfa.nut06a(jd1, jd2)      # radians, vectorized
```
Parameters: interval=16 days, degree=16, 2 components.

### 6.12 Progress Bars

Lightweight `ProgressBar` class (stdlib only, line 38):
- Throttled to 100ms redraws to avoid terminal flicker
- Terminal width detection via `shutil.get_terminal_size()`
- Per-body bars for sequential generation
- Output goes to `sys.stdout`

### 6.13 Execution Strategy

```
Phase 1: ICRS planets (11 bodies)
  Sequential, vectorized Skyfield -> very fast (~seconds for 300yr)

Phase 2: ICRS asteroids + exotics (5 classic + 31 exotic)
  Sequential, spktype21 scalar loop -> tens of seconds (classic) + minutes (exotics)

Phase 3: Analytical bodies (15 bodies)
  Sequential, per-body progress bars -> ~2-3 min for 300yr
```

All three phases run sequentially in the main process.
Parallelization via `ProcessPoolExecutor` was removed because it caused
deadlocks on macOS due to numpy/BLAS/Accelerate initialization in spawned
processes. The group workflow (Section 6.15) provides the primary speedup
mechanism: regenerate only the group that changed.

### 6.14 File Assembly

The `assemble_leb()` function (line 1328):
1. Pre-download asteroid SPKs with full-range coverage
2. Generate body coefficients (phases 1-3 above)
3. Generate nutation coefficients
4. Generate Delta-T sparse table
5. Generate star catalog (from `STAR_CATALOG` in `fixed_stars.py`)
6. Calculate section sizes and offsets
7. Allocate single `bytearray(total_size)`
8. Write header, section directory, body index, coefficients, nutation,
   Delta-T, and star catalog
9. Write buffer to file in one `f.write(buf)` call

### 6.15 Group Generation & Merge

Generating all 62 bodies in a single process can be slow. The **group workflow**
splits generation into four independent runs — one per body group (planets,
asteroids, exotics, analytical) — then merges the partial files into a single
`.leb`. This allows regenerating only the group that changed (e.g. after
updating asteroid SPK files).

#### Body Groups

The `BODY_GROUPS` dict in `generate_leb.py` maps group names to body ID lists:

```python
BODY_GROUPS: dict[str, List[int]] = {
    "planets": sorted(_PLANET_MAP.keys()),  # 11 planets (ICRS_BARY, vectorized)
    # Classic asteroids = SPK bodies minus the exotic registry (stays correct if
    # a 6th classic is ever added to _ASTEROID_NAIF without touching EXOTIC_IDS).
    "asteroids": sorted(set(_ASTEROID_NAIF) - set(EXOTIC_IDS)),
    "exotics": list(EXOTIC_IDS),  # centaurs/TNOs/NEAs (ICRS_BARY, spktype21)
    "analytical": sorted(
        bid
        for bid in BODY_PARAMS
        if bid not in _PLANET_MAP and bid not in _ASTEROID_NAIF
    ),  # ecliptic/helio analytical bodies
}
```

| Group | Bodies | Method | Typical time (base) |
|---|---|---|---|
| `planets` | Sun, Moon, Mercury–Pluto, Earth | Vectorized Skyfield | ~1 s |
| `asteroids` | Chiron, Ceres, Pallas, Juno, Vesta | spktype21 (scalar) | ~15–60 s |
| `exotics` | 31 centaurs / TNOs / main-belt / NEAs (Pholus, Eris, Chariklo, Apophis, …) | spktype21 (scalar) | ~several min |
| `analytical` | Mean/true nodes, mean/true Lilith, 8 Uranians, mean apogee/perigee, osc. apogee | Sequential (scalar) | ~2–3 min |

#### CLI: `--group`

```bash
# Generate only the planets group (partial file)
python scripts/generate_leb.py --tier base --group planets

# Output: data/leb/ephemeris_base_planets.leb  (auto-suffixed)
```

When `--group` is used without an explicit `--output`, the output path is
auto-suffixed by `_group_output_path()`:

```
data/leb/ephemeris_base.leb  +  planets  →  data/leb/ephemeris_base_planets.leb
```

Each partial file is a valid `.leb` with the standard header, section directory,
and body index — it just contains fewer bodies. Nutation, Delta-T, and star
catalog are only generated when their respective bodies are present (nutation
is generated in every group run; stars likewise).

#### CLI: `--merge`

```bash
# Merge the four partial files into one complete file
python scripts/generate_leb.py --tier base --merge \
  data/leb/ephemeris_base_planets.leb \
  data/leb/ephemeris_base_asteroids.leb \
  data/leb/ephemeris_base_exotics.leb \
  data/leb/ephemeris_base_analytical.leb \
  --verify
```

The `merge_leb_files()` function:

1. **Validates JD range consistency** — all inputs must cover the same
   `[jd_start, jd_end]` range (tolerance: 0.5 days).
2. **Checks for duplicate bodies** — raises `ValueError` if any body ID
   appears in more than one input file.
3. **Copies raw coefficient blobs** — zero re-computation; each body's
   Chebyshev coefficients are copied verbatim from the source file.
4. **Takes auxiliary sections from first provider** — nutation, Delta-T,
   and star catalog are taken from the first input file that contains them.
5. **Writes merged file** — single `bytearray` allocation, one `f.write()`.

#### Complete Group Workflow

The recommended workflow for regenerating a tier:

```bash
# Step by step
python scripts/generate_leb.py --tier base --group planets
python scripts/generate_leb.py --tier base --group asteroids
python scripts/generate_leb.py --tier base --group exotics
python scripts/generate_leb.py --tier base --group analytical
python scripts/generate_leb.py --tier base --merge \
  data/leb/ephemeris_base_planets.leb \
  data/leb/ephemeris_base_asteroids.leb \
  data/leb/ephemeris_base_exotics.leb \
  data/leb/ephemeris_base_analytical.leb \
  --verify

# Or all at once via leph
leph leb generate base groups
```

#### Regenerating a Single Group

If only one group needs to be regenerated (e.g. after updating asteroid SPK
files), regenerate just that group and re-merge:

```bash
python scripts/generate_leb.py --tier base --group asteroids
python scripts/generate_leb.py --tier base --merge \
  data/leb/ephemeris_base_planets.leb \
  data/leb/ephemeris_base_asteroids.leb \
  data/leb/ephemeris_base_exotics.leb \
  data/leb/ephemeris_base_analytical.leb \
  --verify
```

#### macOS Deadlock History

Early versions used `ProcessPoolExecutor` to parallelize analytical body
generation across multiple workers. This caused persistent deadlocks on
macOS due to numpy/BLAS/Accelerate initialization in spawned child
processes. Switching to the `spawn` start method (via
`multiprocessing.get_context("spawn")`) partially helped but still hung
intermittently. The parallelization was ultimately **removed entirely** —
analytical bodies now run sequentially in the main process. The group
workflow (generate each group independently, then merge) provides the
primary mechanism for selective regeneration without redoing everything.

---

## 7. Integration with LibEphemeris

### 7.1 Global State (`state.py`)

```python
# state.py globals (line 157-163)
_LEB_FILE: Optional[str] = None      # Path to .leb file
_LEB_READER: Optional["LEBReader"] = None  # Cached reader instance
_CALC_MODE: Optional[str] = None     # None = check env var
_CALC_MODE_ENV_VAR = "LIBEPHEMERIS_MODE"
_VALID_CALC_MODES = ("auto", "skyfield", "leb", "horizons")
```

**`set_calc_mode(mode)`** (line 170):
- Sets the calculation mode (`"auto"`, `"skyfield"`, `"leb"`, `"horizons"`, or `None`)
- `None` resets to environment variable / default

**`get_calc_mode()`** (line 209):
- Returns the effective mode: programmatic override > env var > `"auto"`

**`set_leb_file(filepath)`** (line 231):
- Closes existing reader if any
- Sets `_LEB_FILE` and clears `_LEB_READER` (lazy re-creation)

**`get_leb_reader()`** (line 256):
- Respects `get_calc_mode()`:
  - `"skyfield"` → always returns `None`
  - `"leb"` → returns reader or raises `RuntimeError`
  - `"auto"` → returns reader if configured, else `None`
- Returns cached `_LEB_READER` if available
- Otherwise checks `_LEB_FILE` and `LIBEPHEMERIS_LEB` env var
- Creates `LEBReader(path)` on first access
- Logs warning and returns `None` if file is invalid/corrupt (in `auto` mode)

### 7.2 Dispatch in `calc_ut()` (`planets.py:772-782`)

```python
def calc_ut(tjd_ut, ipl, iflag):
    # ... ECL_NUT handling, MOSEPH stripping ...

    # --- LEB fast path ---
    reader = get_leb_reader()
    if reader is not None:
        try:
            return fast_calc.fast_calc_ut(reader, tjd_ut, ipl, iflag)
        except (KeyError, ValueError):
            pass  # fall through to Skyfield
    # --- END LEB fast path ---

    # Full Skyfield pipeline
    return _calc_body(t, ipl, iflag)
```

The same pattern is used in `calc()` (line 835-845), dispatching to
`fast_calc.fast_calc_tt()`.

### 7.3 Ayanamsa in LEB Mode

When `get_ayanamsa_ut()` or `get_ayanamsa_ex_ut()` is called:
```python
# planets.py line 2250
reader = get_leb_reader()
if reader is not None:
    try:
        from .fast_calc import _calc_ayanamsa_from_leb
        return _calc_ayanamsa_from_leb(reader, jd_tt, sid_mode)
    except KeyError:
        pass  # star-based mode, fall back to Skyfield
```

### 7.4 Exported API

```python
# __init__.py
from .state import set_leb_file, get_leb_reader, set_calc_mode, get_calc_mode
```

All four are accessible as `libephemeris.set_leb_file()`,
`libephemeris.get_leb_reader()`, `libephemeris.set_calc_mode()`, and
`libephemeris.get_calc_mode()`.

---

## 8. Thread Safety and EphemerisContext

### 8.1 LEBReader Thread Safety

`LEBReader` is **inherently thread-safe** because:
- All data structures are populated at construction time and never mutated
- `mmap` read access is safe from multiple threads
- `struct.unpack_from()` is a pure function
- Clenshaw evaluation uses only local variables
- `warm()` only calls `madvise()` (a read-only advisory hint) and is safe
  to call from any thread

### 8.2 EphemerisContext LEB Integration (`context.py`)

`EphemerisContext` provides per-context LEB support:

```python
class EphemerisContext:
    def __init__(self):
        self._leb_file: Optional[str] = None
        self._leb_reader: Optional["LEBReader"] = None
        self.sidereal_mode: int = 1   # Lahiri
        self.sidereal_t0: float = 2451545.0
        self.sidereal_ayan_t0: float = 0.0

    def set_leb_file(self, filepath)   # Per-context LEB
    def get_leb_reader(self)           # Per-context reader (falls back to global)

    def calc_ut(self, tjd_ut, ipl, iflag):
        reader = self.get_leb_reader()
        if reader is None:
            reader = state.get_leb_reader()  # global fallback
        if reader is not None:
            return fast_calc.fast_calc_ut(
                reader, tjd_ut, ipl, iflag,
                sid_mode=self.sidereal_mode,        # thread-safe
                sid_t0=self.sidereal_t0,             # thread-safe
                sid_ayan_t0=self.sidereal_ayan_t0,   # thread-safe
            )
        # ... Skyfield fallback ...
```

**Critical:** Sidereal parameters are passed as **keyword arguments** to
`fast_calc_ut()` / `fast_calc_tt()`, not read from global state. This
prevents race conditions when multiple threads use different sidereal modes.

The sidereal state snapshot happens once at the entry point (`fast_calc_ut`
line 662-667 for global API, or passed explicitly by `EphemerisContext`).

---

## 9. Body Catalog

### 9.1 BODY_PARAMS Table

**Source:** `leb_format.py:149-185`

This is the **single source of truth** for all Chebyshev parameters:

```python
BODY_PARAMS: dict[int, tuple[float, int, int, int]] = {
    # body_id: (interval_days, degree, coord_type, components)
    ...
}
```

| Body ID | Name | Interval | Degree | Coord Type | Components |
|---------|------|----------|--------|------------|------------|
| 0 | Sun | 32 | 13 | ICRS_BARY | 3 |
| 1 | Moon | 4 | 13 | ICRS_BARY | 3 |
| 2 | Mercury | 16 | 15 | ICRS_BARY | 3 |
| 3 | Venus | 16 | 13 | ICRS_BARY | 3 |
| 4 | Mars | 16 | 13 | ICRS_BARY | 3 |
| 5 | Jupiter | 32 | 13 | **ICRS_BARY_SYSTEM** | 3 |
| 6 | Saturn | 32 | 13 | **ICRS_BARY_SYSTEM** | 3 |
| 7 | Uranus | 64 | 13 | **ICRS_BARY_SYSTEM** | 3 |
| 8 | Neptune | 64 | 13 | **ICRS_BARY_SYSTEM** | 3 |
| 9 | Pluto | 64 | 11 | **ICRS_BARY_SYSTEM** | 3 |
| 10 | Mean Node | 8 | 13 | ECLIPTIC | 3 |
| 11 | True Node | 8 | 13 | ECLIPTIC | 3 |
| 12 | Mean Apogee | 8 | 13 | ECLIPTIC | 3 |
| 13 | Oscu Apogee | **4** | **15** | ECLIPTIC | 3 |
| 14 | Earth | 4 | 13 | ICRS_BARY | 3 |
| 15 | Chiron | 8 | 13 | ICRS_BARY | 3 |
| 17 | Ceres | 8 | 13 | ICRS_BARY | 3 |
| 18 | Pallas | 8 | 13 | ICRS_BARY | 3 |
| 19 | Juno | 8 | 13 | ICRS_BARY | 3 |
| 20 | Vesta | 8 | 13 | ICRS_BARY | 3 |
| 21 | Interp Apogee | **1** | **17** | ECLIPTIC | 3 |
| 22 | Interp Perigee | **1** | **17** | ECLIPTIC | 3 |
| 40 | Cupido | 256 | 7 | HELIO_ECL | 3 |
| 41 | Hades | 256 | 7 | HELIO_ECL | 3 |
| 42 | Zeus | 256 | 7 | HELIO_ECL | 3 |
| 43 | Kronos | 256 | 7 | HELIO_ECL | 3 |
| 44 | Apollon | 256 | 7 | HELIO_ECL | 3 |
| 45 | Admetos | 256 | 7 | HELIO_ECL | 3 |
| 46 | Vulkanus | 256 | 7 | HELIO_ECL | 3 |
| 47 | Poseidon | 256 | 7 | HELIO_ECL | 3 |
| 48 | Transpluto | 256 | 7 | HELIO_ECL | 3 |

**Total: 31 core bodies** (listed above). `BODY_PARAMS` additionally
registers **31 exotic minor bodies** (`EXOTIC_IDS`, see §9.3) — Pholus
(id 16), the TNOs, and extra centaur / main-belt / NEA asteroids (ids
10010–235088) — for **62 bodies** in all.

Bodies marked **bold** differ from the original design document:
- **Jupiter-Pluto (5-9):** Use `COORD_ICRS_BARY_SYSTEM` (4) instead of
  `COORD_ICRS_BARY` (0). Stores pure system barycenters; COB correction
  applied at runtime. This was the key fix for outer planet precision.
- **OscuApogee (13), InterpApogee (21), InterpPerigee (22):** Tightened
  to interval=4, degree=15 (from interval=8, degree=13) for sub-arcsecond
  precision on these high-frequency analytical bodies.

### 9.2 Parameter Design Rationale

Parameters were aggressively tuned in the precision improvement work to
minimize fitting error for all bodies. The key insight is that even
slow-moving outer planets are stored in ICRS barycentric coordinates, so
their *apparent* position (as seen from Earth) varies faster than the raw
orbital motion due to Earth's own motion and parallax effects.

- **Moon, Earth (interval=4, degree=13):** Earth is the observer for all
  geocentric calculations — any error in Earth's position propagates to
  every other body. Moon is the fastest-moving body (~13 deg/day). Both
  use 4-day intervals for maximum precision.
- **Mercury (interval=16, degree=15):** Most eccentric orbit, needs highest
  degree to capture orbital variations.
- **Venus, Mars (interval=16, degree=13):** Halved from 32 days to reduce
  latitude error caused by ICRS→ecliptic pipeline amplification at close
  approach (Venus at ~0.26 AU, Mars at ~0.37 AU).
- **Jupiter, Saturn (interval=32, degree=13):** Use `COORD_ICRS_BARY_SYSTEM`
  (pure system barycenter, COB at runtime). This eliminated the high-frequency
  moon oscillation fitting problem that previously caused 3.95" errors.
- **Uranus, Neptune (interval=64, degree=13):** Also use `COORD_ICRS_BARY_SYSTEM`.
  Halved from 128 days and degree increased from 9 to 13.
- **Pluto (interval=64, degree=11):** Uses `COORD_ICRS_BARY_SYSTEM`. The
  long, slow orbit fits well at 64-day intervals (0.0005" fit error).
- **Asteroids (interval=8, degree=13):** Reduced from 32 days. Eccentric
  and perturbed orbits need short intervals for sub-arcsecond accuracy.
- **Hypotheticals (interval=256, degree=7):** Pure analytical Keplerian
  orbits (no perturbations), so very long intervals with a low degree fit
  to ~0" error while keeping the segment count small.
- **OscuApogee (interval=4, degree=15):** Tightened from interval=8,
  degree=13 to achieve <0.001" precision on this high-frequency
  analytical body.
- **InterpApogee, InterpPerigee (interval=1, degree=17):** 1-day
  intervals — the interpolated apsides oscillate too rapidly for longer
  segments.
- **Other ecliptic bodies (interval=8, degree=13):** Unchanged -- already tight.
- **Nutation (interval=16, degree=16):** Halved from 32 days to reduce
  obliquity error, which affects latitude of all bodies.

### 9.3 Bodies NOT in LEB (Skyfield Fallback)

When LEB is active and a body is **not** present in the loaded `.leb`
file, the library raises `KeyError`/`ValueError` internally, which is
caught by `calc_ut()` / `calc()`, and the request falls through to the
full Skyfield pipeline. This is completely transparent to the caller.

The same fallback is triggered by:

- **Unsupported flags:** `FLG_ICRS` (the only flag that always falls
  back; `FLG_TOPOCTR`, `FLG_XYZ`, `FLG_RADIANS` and `FLG_NONUT` are
  handled on the LEB path — see §5.2)
- **JD out of range:** Julian Day outside the body's coverage in the file
- **Star-based sidereal modes:** e.g., `SIDM_TRUE_REVATI` (requires
  fixed star position not available in LEB fast path)

#### Exotic minor bodies (served from LEB when present)

The LEB catalog is **not** limited to the 31 core bodies. `BODY_PARAMS`
also registers **31 exotic minor bodies** (`EXOTIC_IDS` in
`libephemeris/exotic_bodies.py`) — centaurs, trans-Neptunian objects,
large main-belt asteroids, and near-Earth asteroids — generated into the
LEB `exotics` group (`ephemeris_{tier}_exotics.leb` for LEB1,
`{tier}_exotics.leb2` for LEB2):

| Class | Bodies | Count |
|-------|--------|-------|
| TNOs | Eris, Sedna, Haumea, Makemake, Quaoar, Orcus, Ixion, Gonggong, Varuna | 9 |
| Centaurs | Pholus (id 16), Nessus, Asbolus, Chariklo, Hidalgo | 5 |
| Main-belt | Hygiea, Interamnia, Davida, Europa, Sylvia, Psyche, Sappho, Pandora, Lilith-ast | 9 |
| NEAs | Eros, Amor, Apophis, Itokawa, Ryugu, Toro, Toutatis, Icarus | 8 |

When the loaded LEB file contains the `exotics` group these bodies are
served **directly from LEB** (SPK-derived Chebyshev coefficients over
each object's JPL SPK coverage window). They fall back to the
SPK → Skyfield pipeline **only** when the loaded file lacks the body —
e.g. the exotics group was not generated/downloaded, the requested JD is
outside the object's SPK window, or on the extended tier, where the 8
chaotic NEAs are intentionally excluded from generation
(`EXOTIC_EXTENDED_IDS` keeps the 23 regular exotics only).

#### Bodies that always fall back to Skyfield

These are never stored as LEB Chebyshev data:

| Category | Bodies | IDs | Count | How computed |
|----------|--------|-----|-------|--------------|
| Additional hypotheticals | Nibiru, Harrington, Leverrier, Adams, Lowell, Pickering, Vulcan, Selena, Proserpina, Waldemath | 49–58 | 10 | Keplerian from `hypothetical.py` |
| Fixed stars | full star catalog (Regulus, Spica, Aldebaran, …) | FIXSTAR_OFFSET + n | 1447 | `fixed_stars.py` (see §9.4) |
| Planetary moons | Io, Europa, Ganymede, Callisto, Titan, Triton, Charon, etc. | MOON_OFFSET + n | 21 | SPK via `planetary_moons.py` |
| Astrological angles | Ascendant, MC, Descendant, IC, Vertex, Antivertex | 9000–9005 | 6 | `angles.py` (house-based) |
| Arabic parts | Pars Fortunae, Pars Spiritus, Pars Amoris, Pars Fidei | 9100–9103 | 4 | `arabic_parts.py` (derived) |
| Nutation/obliquity | ECL_NUT | -1 | 1 | LEB nutation section (§9.4) |

Any minor body outside the `exotics` registry (e.g. Bennu) likewise falls
back to the SPK → Skyfield pipeline.

**Total bodies NOT in LEB Chebyshev data:** ~1489 (10 hypotheticals + 1447
stars + 21 moons + 10 angles/parts + 1 nutation).

**Why the core/exotics split:** the 31 **core** bodies are used in nearly
every chart; their compressed LEB2 `core` group is small enough (~10.6 MB)
to bundle with the wheel. The 31 **exotic** bodies are SPK-derived and much
larger on disk (the base `exotics` group is ~59 MB compressed), so they ship
as a separate optional group served from LEB once downloaded. Fixed stars,
planetary moons, and chart angles/parts stay on the Skyfield pipeline — they
are cheap to compute analytically or vary with observer/house settings, so
precomputing them would add size with little benefit. The Skyfield fallback
provides identical accuracy at the cost of ~120 µs per evaluation (vs ~8 µs
for LEB bodies).

### 9.4 Auxiliary LEB Sections (Non-Chebyshev Data)

In addition to the Chebyshev polynomial coefficients for its bodies, each
`.leb` file contains three auxiliary data sections that accelerate other
parts of the calculation pipeline:

#### Nutation (Section 2)

Stores Chebyshev polynomial approximations for **dpsi** (nutation in
longitude) and **deps** (nutation in obliquity) using the IAU 2006/2000A
model. These are used by Pipelines A/A' for the precession-nutation matrix,
by ecliptic body pipelines for true obliquity, and by sidereal ayanamsha
calculations.

- **Parameters:** interval=16 days, degree=16
- **Precision:** sub-milliarcsecond (matches `erfa.nut06a()`)
- **Access:** `reader.eval_nutation(jd_tt)` → `(dpsi, deps)` in radians
- **Evaluation time:** ~0.8 µs

#### Delta-T (Section 3)

Stores a sparse table of historical Delta-T values (TT − UT1) at regular
intervals. Used by `fast_calc_tt()` for reverse UT↔TT lookups, but
**not** used by `fast_calc_ut()` for the forward UT→TT conversion (which
uses `deltat()` for higher precision — see §5).

- **Format:** array of `(jd, delta_t_days)` pairs
- **Interpolation:** linear between adjacent entries
- **Access:** `reader.delta_t(jd)` → Delta-T in days
- **Evaluation time:** ~0.3 µs
- **Caveat:** Linear interpolation introduces up to ~0.004s error near 1985.
  This is why `fast_calc_ut()` uses `deltat()` instead.

#### Star Catalog (Section 4)

Stores a snapshot of the fixed star catalog (J2000 positions, proper motions,
parallax, radial velocity) for the full star catalog defined in
`fixed_stars.py` (1447 entries).
This is **read-only reference data** — star positions are not stored as
Chebyshev polynomials because proper motion is a simple linear correction
that doesn't benefit from polynomial approximation.

- **Format:** per-star entries with RA, Dec, pmRA, pmDec, parallax, radial
  velocity, visual magnitude
- **Access:** `reader.get_star(hip_number)` → star data dict
- **Note:** Fixed star calculations still go through the full Skyfield
  pipeline when called via `calc_ut()`. The star catalog in LEB is
  used internally for sidereal ayanamsha calculations involving reference
  stars.

---

## 10. Precision and Validation

### 10.1 Generation-Time Verification

Every body is verified during generation. For each Chebyshev segment,
10 uniformly-spaced test points are evaluated and compared against the
reference function. The maximum error across all segments is reported.

Typical generation-time errors:

| Body | Max Error | Unit |
|------|-----------|------|
| Sun | <1e-12 | AU (~0.00002") |
| Moon | <5e-11 | AU (~0.001") |
| Mercury | <1e-11 | AU |
| Venus | <1e-12 | AU |
| Mars | <1e-12 | AU |
| Jupiter | <1e-12 | AU |
| Mean Node | <1e-12 | degrees |
| True Node | <1e-9 | degrees (~4e-6") |
| Interp Apogee | <1e-7 | degrees (~4e-4") |

### 10.2 End-to-End Precision (vs Skyfield Reference)

The compare test suite (`tests/test_leb/compare/`) validates LEB output
against Skyfield for all 31 core bodies across hundreds of dates per tier.

**All 31 core bodies achieve <0.001 arcsecond geocentric position precision**
on all three tiers (base, medium, and extended). (The 31 exotic minor bodies
of §9.3 are SPK-derived and best-effort; NEAs in particular use looser
targets.) This was accomplished through:

1. `COORD_ICRS_BARY_SYSTEM` storage for outer planets (eliminates COB
   oscillation fitting errors)
2. PPN gravitational deflection (Sun, Jupiter, Saturn)
3. Runtime COB correction at observer time (not retarded time)
4. `deltat()` for UT->TT conversion (not reader's sparse table)
5. Asteroid pipeline via `_SpkType21Target` VectorFunction wrapper
6. Tightened Chebyshev parameters for bodies 13, 21, 22

#### Base Tier (1850-2150, verified)

| Group | Bodies | Worst Case Body | Max Error |
|-------|--------|-----------------|-----------|
| Planets (11) | Sun-Pluto, Earth | Moon | 0.000332" |
| Asteroids (5) | Chiron, Ceres-Vesta | Juno | 0.000045" |
| Ecliptic (6) | Nodes, Lilith | OscuApog | 0.000049" |
| Hypothetical (9) | Uranians, Transpluto | all | ~0.000000" |

**Tests passed:** 404 comparison tests (planets, asteroids, hypothetical,
lunar, velocities, distances, flags, sidereal).

#### Medium Tier (1550-2650, verified)

| Group | Bodies | Worst Case Body | Max Error |
|-------|--------|-----------------|-----------|
| Planets (11) | Sun-Pluto, Earth | Moon | 0.000325" |
| Asteroids (5) | Chiron, Ceres-Vesta | Vesta | 0.000036" |
| Ecliptic (6) | Nodes, Lilith | OscuApog | 0.000075" |
| Hypothetical (9) | Uranians, Transpluto | all | ~0.000000" |

**Tests passed:** 904 comparison tests (all planet/asteroid/hypothetical/lunar
tests plus eclipses, crossings, stations, rise/transit, elongation, ayanamsha,
nutation, houses, gauquelin).

#### Extended Tier (-5000 to 5000 CE, verified)

| Group | Bodies | Worst Case Body | Max Error |
|-------|--------|-----------------|-----------|
| Planets (11) | Sun-Pluto, Earth | Mars | 0.000010" |
| Asteroids (5) | Chiron, Ceres-Vesta | Pallas | 0.000018" |
| Ecliptic (6) | Nodes, Lilith | OscuApog | 0.054" * |
| Hypothetical (9) | Uranians, Transpluto | all | ~0.000000" |

\* Ecliptic body precision is limited by Meeus polynomial degradation
beyond ±20 centuries from J2000.0. OscuApogee uses analytical lunar
formulas that lose accuracy at extreme dates. Within ±1000 CE, ecliptic
body errors are <0.001".

**Tests passed:** 261 comparison tests (planets, hypothetical, velocities,
lunar, flags, ancient/future sub-ranges, boundary dates).

**File size:** ~1.6 GB (de441.bsp, -5000 to 5000 CE, 10,000 years).

#### Test Tolerances (as configured in `conftest.py`)

| Field | Base | Medium | Extended | Notes |
|-------|------|--------|----------|-------|
| `POSITION_ARCSEC` | 0.001 | 0.001 | 0.001 | All planets including outer |
| `ASTEROID_ARCSEC` | 0.001 | 0.001 | 0.001 | |
| `ECLIPTIC_ARCSEC` | 0.001 | 0.001 | 0.1 | Meeus limit at extreme dates |
| `EQUATORIAL_ARCSEC` | 0.02 | 0.02 | 0.02 | Heliocentric amplification |
| `J2000_ARCSEC` | 0.001 | 0.001 | 0.001 | |
| `SIDEREAL_ARCSEC` | 0.001 | 0.001 | 0.001 | |
| `HYPOTHETICAL_ARCSEC` | 0.001 | 0.001 | 0.001 | |
| `DISTANCE_AU` | 5e-6 | 5e-6 | 5e-6 | |

### 10.3 Architectural Limitations

**Heliocentric/equatorial Moon amplification (~0.01"):** When computing
Moon's heliocentric position, the geocentric error is amplified by the
ratio of heliocentric to geocentric distance. This produces ~0.01" errors
in heliocentric coordinates. The `EQUATORIAL_ARCSEC` tolerance of 0.02"
accommodates this.

**Asteroid latitude velocity (0.19-1.7 deg/day):** The ICRS->ecliptic
pipeline amplifies velocity errors by `1/geocentric_distance` for nearby
asteroids (Ceres, Pallas, Juno, Vesta). This is handled in tests via a
separate `ASTEROID_SPEED_LAT_DEG_DAY` tolerance of 1.7 deg/day. Chiron
is less affected due to greater distance.

These are inherent to the coordinate transformation and cannot be fixed
without changing the storage format.

### 10.4 Asteroid Precision Caveat

Asteroid precision depends entirely on how the LEB file was generated:

- **With spktype21 SPK:** Position errors <1" — excellent
- **With scalar calc() fallback:** Position errors can reach ~1500" — unacceptable
- **With Keplerian fallback (no SPK):** Even worse

Always ensure asteroid SPK files cover the full tier date range before
generating. Use `_spk_covers_range()` to verify.

---

## 11. Performance

### 11.1 Per-Evaluation Costs

| Operation | Time |
|-----------|------|
| Clenshaw evaluation (1 component, degree 13) | ~0.5 us |
| `eval_body()` (3 components + derivatives) | ~1.5 us |
| `eval_nutation()` | ~0.8 us |
| `delta_t()` | ~0.3 us |
| Pipeline A full (ICRS -> ecliptic of date, with speed) | ~8 us |
| Pipeline B full (ecliptic direct, with speed) | ~2 us |
| Skyfield calc_ut() for comparison | ~120 us |
| **Speedup (Pipeline A)** | **~14x** |

### 11.2 Generation Performance

| Configuration | Time |
|---------------|------|
| Old scalar (1 worker, 300yr) | ~18 min |
| Vectorized (10yr, 1 worker) | 18.0s |
| Vectorized (10yr, 4 workers) | 6.4s |
| Vectorized (5yr, 4 workers) | 3.8s |
| Vectorized (300yr, no asteroids, 4 workers) | 170s (~2.8 min) |
| Vectorized (300yr, with spktype21, 4 workers) | ~3-4 min (estimated) |

### 11.3 Vectorization Speedups

| Method | Scalar Cost | Vectorized Cost | Speedup |
|--------|------------|-----------------|---------|
| Skyfield `target.at(t)` (Sun) | 61 us | 0.4 us | **150x** |
| Skyfield `target.at(t)` (Moon) | 121 us | 0.7 us | **170x** |
| `erfa.nut06a()` (nutation) | 32 us | ~0.5 us | **~50x** |
| `spktype21.compute_type21()` | 65 us | N/A (scalar) | 36x vs calc |
| True lunar node (pure Python) | 315 us | N/A | not vectorizable |
| Interp perigee (pure Python) | 329 us | N/A | not vectorizable |

---

## 12. Commands Reference

### 12.1 leph Tasks (Developer CLI)

```bash
# Group generation (recommended — avoids fork-deadlock, allows partial regen)
leph leb generate base groups      # All four groups + merge in one command
leph leb generate medium groups    # Medium tier
leph leb generate extended groups  # Extended tier

# Verify existing LEB files
leph leb verify base               # Verify base tier
leph leb verify medium             # Verify medium tier

# Download (convenience wrappers for libephemeris CLI)
libephemeris download leb-base       # Download base tier LEB (~53 MB)
libephemeris download leb-medium     # Download medium tier LEB (~175 MB)
libephemeris download leb-extended   # Download extended tier LEB (~1.6 GB)

# Release — upload to GitHub Releases (requires: gh auth login)
leph release leb 1.0.0             # Upload all tiers + update hashes in download.py
leph release leb-dry-run 1.0.0     # Show what would be uploaded (no changes)

# Testing
leph test leb-format all           # All LEB format tests (excludes @slow)
leph test leb-format precision     # Full precision suite (slow, all tiers)
```

### 12.2 Direct pytest

```bash
# By file
pytest tests/test_leb/test_leb_format.py -v
pytest tests/test_leb/test_leb_reader.py -v
pytest tests/test_leb/test_fast_calc.py -v
pytest tests/test_leb/test_generate_leb.py -v
pytest tests/test_leb/test_leb_precision.py -v
pytest tests/test_leb/test_context_leb.py -v

# By marker
pytest tests/test_leb/ -v -m "not slow"
pytest tests/test_leb/ -v -m "precision"
```

### 12.3 Manual Generation

```bash
# Generate with custom options
python scripts/generate_leb.py \
  --tier base \
  --workers 8 \
  --verify \
  --verify-samples 1000

# Generate specific bodies only
python scripts/generate_leb.py \
  --output test.leb \
  --start 2000 --end 2030 \
  --bodies 0,1,2,3,4 \
  --workers 4

# Group generation (recommended for base tier)
python scripts/generate_leb.py --tier base --group planets
python scripts/generate_leb.py --tier base --group asteroids
python scripts/generate_leb.py --tier base --group exotics
python scripts/generate_leb.py --tier base --group analytical
python scripts/generate_leb.py --tier base --merge \
  data/leb/ephemeris_base_planets.leb \
  data/leb/ephemeris_base_asteroids.leb \
  data/leb/ephemeris_base_exotics.leb \
  data/leb/ephemeris_base_analytical.leb \
  --verify

# Regenerate only the asteroids group, then re-merge
python scripts/generate_leb.py --tier base --group asteroids
python scripts/generate_leb.py --tier base --merge \
  data/leb/ephemeris_base_planets.leb \
  data/leb/ephemeris_base_asteroids.leb \
  data/leb/ephemeris_base_exotics.leb \
  data/leb/ephemeris_base_analytical.leb \
  --verify

# Quick validation after regeneration
python3 -c "
from libephemeris.leb_reader import LEBReader
reader = LEBReader('data/leb/ephemeris_base.leb')
pos, vel = reader.eval_body(0, 2451545.0)  # Sun at J2000
print(f'Sun ICRS: x={pos[0]:.10f} y={pos[1]:.10f} z={pos[2]:.10f} AU')
reader.close()
"
```

### 12.4 Release Workflow

LEB files are hosted as assets on the GitHub release tagged `data-v1`.
The release script (`scripts/release_leb.py`) handles upload and hash updates.

**Prerequisites:**

- `gh` CLI installed and authenticated (`gh auth login`)
- LEB files generated locally (in `data/leb/` or `~/.libephemeris/leb/`)

**Full workflow (generate + release):**

```bash
# 1. Generate LEB file(s)
leph leb generate medium groups    # recommended group workflow

# 2. Dry run — verify what would be uploaded
leph release leb-dry-run 1.0.0

# 3. Upload to GitHub Releases + auto-update hashes in download.py
leph release leb 1.0.0

# 4. Commit the updated download.py
git add libephemeris/download.py
git commit -m "update LEB medium tier hash after regeneration"
```

**What the release script does:**

1. Searches for `.leb` files in `data/leb/` (repo) and `~/.libephemeris/leb/`
2. Computes SHA256 hash and file size for each file
3. Uploads to GitHub release `data-v1` via `gh release upload --clobber`
4. With `--update-hashes`: updates `DATA_FILES` in `libephemeris/download.py`
   with the new SHA256 and size_mb values

**Direct script usage (without leph):**

```bash
python scripts/release_leb.py --version 1.0.0 --tier medium --update-hashes
python scripts/release_leb.py --version 1.0.0 --dry-run
python scripts/release_leb.py --version 1.0.0 --tier base --tag data-v1
```

### 12.5 User Download

End users download pre-generated LEB files via the library CLI (no `leph` needed):

```bash
# Install the library
pip install libephemeris

# Download LEB for the desired tier
libephemeris download leb-base       # ~53 MB, 1850-2150
libephemeris download leb-medium     # ~175 MB, 1550-2650
libephemeris download leb-extended   # ~1.6 GB, -5000 to +5000
```

Files are saved to `~/.libephemeris/leb/` and auto-discovered at runtime
based on the active precision tier — no further configuration needed.

Programmatic download is also available:

```python
from libephemeris import download_leb_for_tier

download_leb_for_tier("medium")  # downloads + activates
```

---

## 13. LEB2 Compressed Format

### 13.1 Overview

LEB2 is a compressed variant of the LEB format that uses error-bounded lossy
compression to reduce file sizes by 4-10x while maintaining <0.001" precision.
The compression is transparent: `open_leb()` auto-detects the format via magic
bytes (`LEB1` vs `LEB2`), and the runtime API is identical.

**Motivation:** LEB1 files are too large to bundle in a wheel: the base-tier
core set totals ~102 MB across the separate group files (over PyPI's 100 MB
limit), and even the merged base-tier file is ~53 MB. LEB2 compresses the core
body set to ~10.6 MB, enabling `pip install libephemeris` to include precomputed
ephemeris with zero additional downloads.

**Dependency:** `zstandard` (required, ~200 KB wheel).

### 13.2 Binary File Format

**Magic and version:**

```
Magic:   b"LEB2" (4 bytes)
Version: 2 (uint32)   — current chunked format (LEB2_VERSION)
         1            — legacy monolithic format, still readable (LEB2_VERSION_V1)
```

**Format versions.** The current conversion pipeline writes **v2
(chunked)**: each body's coefficient stream is split into ~10-year
temporal chunks, each compressed independently, with a per-body chunk
index; the reader decompresses only the chunk containing the requested
JD (lazy, per-chunk cache), instead of the whole body. **v1
(monolithic)** stored each body as a single compressed blob;
`LEB2Reader` keys off the header version and still reads v1 files
through the monolithic decode path. Sections 13.3–13.8 below describe
the shared per-body compression primitives (mantissa truncation,
coeff-major reorder, byte shuffle, zstd), which apply to each v2 chunk
exactly as they applied to a whole v1 body.

**Overall layout:**

```
Offset      Content                          Size
────────────────────────────────────────────────────
0x0000      File Header                      64 bytes (same struct as LEB1)
0x0040      Section Directory                N × 24 bytes
variable    Section 0: Body Index            body_count × 68 bytes (CompressedBodyEntry)
variable    Section 6: Compressed Chebyshev  per-body compressed blobs (contiguous)
variable    Section 2: Nutation Data         uncompressed (same as LEB1)
variable    Section 3: Delta-T Table         uncompressed (same as LEB1)
variable    Section 4: Star Catalog          uncompressed (same as LEB1)
```

**Key differences from LEB1:**

- Magic is `b"LEB2"` instead of `b"LEB1"`
- Body index uses `CompressedBodyEntry` (68 bytes) instead of `BodyEntry` (52 bytes)
- Chebyshev data is in section 6 (`SECTION_COMPRESSED_CHEBYSHEV`) instead of section 1
- Header `flags` field encodes `COMPRESSION_ZSTD_TRUNC_SHUFFLE = 1`
- Auxiliary sections (nutation, delta-T, stars) remain uncompressed

### 13.3 CompressedBodyEntry (68 bytes)

```c
// Struct format: "<iIIdddIIQQQ"
struct CompressedBodyEntry {
    // First 52 bytes: identical to LEB1 BodyEntry
    int32_t  body_id;
    uint32_t coord_type;
    uint32_t segment_count;
    float64  jd_start;
    float64  jd_end;
    float64  interval_days;
    uint32_t degree;
    uint32_t components;
    uint64_t data_offset;        // offset to COMPRESSED blob

    // Additional 16 bytes (LEB2-specific)
    uint64_t compressed_size;    // size of compressed blob in bytes
    uint64_t uncompressed_size;  // size of raw coefficients in bytes
};
```

**Python dataclass:** `leb_format.CompressedBodyEntry`
**Constant:** `COMPRESSED_BODY_ENTRY_SIZE = 68`

Each body's compressed blob is independently decompressible. Only the bodies
actually queried at runtime get decompressed.

### 13.4 Compression Pipeline

```
Raw float64 coefficients
    ↓
[1] Mantissa truncation  — zero unneeded mantissa bits per coefficient order
    ↓                      (high-order coefficients need very few bits)
[2] Coefficient-major     — reorder (segments, coeffs) → (coeffs, segments)
    reorder                 so same-order coefficients are contiguous
    ↓
[3] Byte shuffle          — transpose byte lanes (HDF5/Blosc-style)
    ↓
[4] zstd level 19         — the regularized data compresses well
    ↓
Compressed blob
```

**Step 1 — Mantissa truncation** (lossy, error-bounded):

For each coefficient order `k`, the minimum mantissa bits are computed:

```
bits_needed(k) = ceil(-log2(target / max|c_k|))
```

If `max|c_k| < target`, the coefficient is zeroed entirely (0 bits).
The excess bits in the IEEE 754 mantissa are zeroed via bitmask:

```python
mask = 0xFFFFFFFFFFFFFFFF << (52 - keep_bits)
uint64_repr &= mask
```

Example for Moon (degree 13, base tier):

| Coefficient | Magnitude | Bits needed | Wasted bits zeroed |
|-------------|-----------|-------------|-------------------|
| c0 | ~1.0 AU | 28 | 24 |
| c1 | ~0.01 AU | 23 | 29 |
| c5 | ~1e-8 AU | 4 | 48 |
| c6-c13 | < 5e-9 AU | 0 | 52 (zeroed entirely) |

**Step 2 — Coefficient-major reorder:**

Transposes from segment-major (how LEB1 stores data) to coefficient-major,
grouping all c0 values together, all c1 together, etc. This creates long
runs of values with similar exponents and zeroed mantissa bits.

**Step 3 — Byte shuffle** (HDF5/Blosc technique):

Transposes byte lanes so all byte-0 positions cluster, all byte-1 positions
cluster, etc. The zeroed mantissa bytes compress to nearly nothing.

**Step 4 — zstd level 19:**

Maximum practical compression level. Generation is slow but one-time;
decompression runs at ~5 GB/s.

**Decompression** is the exact reverse: zstd → unshuffle → segment-major
reorder. The truncated floats are valid float64 values — no special
handling at evaluation time.

### 13.5 Per-Body Precision Targets

The default target is `5e-9 AU` (≈0.001"). Bodies with small geocentric
distances use tighter targets because positional errors are amplified into
angular errors by `1/distance`, and because Moon/Earth positions feed into
light-time, deflection, and aberration corrections for **all** other bodies.

**Source:** `BODY_TARGET_AU` in `leb_compression.py`

| Body | Target (AU) | Reason |
|------|-------------|--------|
| Moon (1) | 1e-12 | d_geo ~0.002 AU, amplification ~500x |
| Earth (14) | 1e-12 | Used in corrections for every body |
| Sun (0) | 1e-10 | Deflector for all bodies |
| Mercury (2) | 1e-10 | d_geo ~0.55 AU, fast orbit |
| Venus (3) | 1e-10 | d_geo ~0.27 AU, closest to Earth |
| Mars (4) | 1e-10 | d_geo ~0.37 AU |
| All others | 5e-9 | Default (0.001" at 1 AU) |

### 13.6 Modular File Architecture

LEB2 files are organized into **body groups** instead of one monolithic file:

| Group | Bodies | Base size | Description |
|-------|--------|-----------|-------------|
| `core` | 14 | 10.6 MB | Sun-Pluto, Earth, Mean/True Node, Mean Apogee |
| `asteroids` | 5 | 8.7 MB | Chiron, Ceres, Pallas, Juno, Vesta |
| `exotics` | 31 | 59.0 MB | Centaurs, TNOs, main-belt & NEA minor bodies (Pholus, Eris, …) |
| `apogee` | 3 | 11.4 MB | Oscu Apogee, Interp Apogee/Perigee |
| `uranians` | 9 | 2.1 MB | Cupido-Transpluto |

The `exotics` group (`LEB2_GROUPS` in `libephemeris/leb_groups.py`) is by far
the largest and ships as a separate optional download; the pip wheel bundles
only the `core` group.

**File naming convention:** `{tier}_{group}.leb2` (e.g. `base_core.leb2`,
`medium_asteroids.leb2`).

**Output directory:** `data/leb2/`

### 13.7 CompositeLEBReader

The `CompositeLEBReader` (`leb_composite.py`) wraps multiple LEB readers
and dispatches `eval_body()` to the reader containing the requested body.

**Auto-discovery:** when `set_leb_file()` is called with a group file
(e.g. `base_core.leb2`), companion files in the same directory with the
same tier prefix are automatically loaded.

```python
import libephemeris as swe

swe.set_leb_file("data/leb2/base_core.leb2")  # companions auto-discovered
swe.set_calc_mode("leb")

# Bodies from different files work transparently
swe.calc_ut(2451545.0, swe.SUN, swe.FLG_SPEED)     # from core
swe.calc_ut(2451545.0, swe.CHIRON, swe.FLG_SPEED)   # from asteroids
swe.calc_ut(2451545.0, 40, swe.FLG_SPEED)              # from uranians
```

**How it works:**

1. `from_file_with_companions(path)` extracts the tier prefix from the filename
2. Discovers all `{prefix}_*.leb` files in the same directory
3. Opens each with `open_leb()` (auto-detects LEB1/LEB2)
4. Builds a `body_id → reader` dispatch dict
5. Auxiliary data (nutation, delta-T, stars) comes from the first reader that has it

**Programmatic usage:**

```python
from libephemeris.leb_composite import CompositeLEBReader

# From a directory (opens all .leb files)
reader = CompositeLEBReader.from_directory("data/leb2/")

# From a single file + companions
reader = CompositeLEBReader.from_file_with_companions("data/leb2/base_core.leb2")
```

### 13.8 LEB2Reader

**Source:** `libephemeris/leb2_reader.py`

`LEB2Reader` provides the same interface as `LEBReader` with one key
difference: coefficients are decompressed **lazily** on first access
per body, then cached in memory.

```
LEB1 hot path:  jd → O(1) index → struct.unpack_from(mmap, offset) → Clenshaw
LEB2 hot path:  jd → O(1) index → struct.unpack_from(cache[body_id], offset) → Clenshaw
                                    └── bytes buffer, populated once per body
```

After first access, performance is identical to LEB1 (~1.5 µs/eval).
First-access cost per body: ~0.2-1 ms (zstd decompression at ~5 GB/s).

The Clenshaw evaluation functions (`_clenshaw`, `_clenshaw_with_derivative`,
`_deriv_coeffs`) are imported from `leb_reader.py` — zero code duplication.

**Thread safety:** same as LEBReader. The `_cache` dict assignment is atomic
under CPython's GIL; once populated, the `bytes` buffer is immutable.

### 13.9 Key Modules

| Module | Purpose |
|--------|---------|
| `leb_compression.py` | `compress_body()`, `decompress_body()`, `compute_mantissa_bits()`, shuffle/unshuffle |
| `leb2_reader.py` | `LEB2Reader` — lazy per-body decompression, same interface as `LEBReader` |
| `leb_composite.py` | `CompositeLEBReader` — wraps multiple readers, dispatches by body_id |
| `leb_reader.py` | `open_leb()` factory — auto-detects LEB1/LEB2 via magic bytes |
| `scripts/generate_leb2.py` | CLI: `convert`, `convert-all`, `generate`, `verify` |
| `scripts/test_leb2_precision.py` | Fast precision test: all bodies × 6 flags × N dates per tier |

### 13.10 Generation Workflow

LEB2 files are produced by **converting** existing LEB1 files. The conversion
applies the compression pipeline (§13.4) to each body's raw coefficients.

**Prerequisites:** a valid LEB1 file for the target tier must exist in `data/leb/`.

**Recommended workflow (base tier):**

```bash
# Step 1: Generate LEB1 (if not already present)
leph leb generate base groups

# Step 2: Convert LEB1 → LEB2 (all 5 groups)
leph leb2 convert base

# Step 3: Verify LEB2 against LEB1 reference
leph leb2 verify base

# Step 4: Run precision tests
leph test leb2-format precision-base
```

**From scratch (no LEB1):** the `generate` subcommand creates a temporary
LEB1 file, then converts it:

```bash
python scripts/generate_leb2.py generate --tier base --group core -o data/leb2/base_core.leb2
```

### 13.11 Commands Reference

**leph tasks (developer CLI):**

```bash
# Convert LEB1 → LEB2 (all groups for a tier)
leph leb2 convert base              # Base tier → data/leb2/base_{core,asteroids,exotics,apogee,uranians}.leb2
leph leb2 convert medium            # Medium tier
leph leb2 convert extended          # Extended tier

# Verify against LEB1 reference
leph leb2 verify base               # 500 samples per body, compare vs LEB1

# Unit tests
leph test leb2-format all           # Compression round-trip + reader tests

# Precision tests (end-to-end via calc_ut)
leph test leb2-format precision-base       # Base tier (~15s)
leph test leb2-format precision-all        # All tiers (~45s)
```

**Direct CLI (`scripts/generate_leb2.py`):**

```bash
# Convert a single group
python scripts/generate_leb2.py convert data/leb/ephemeris_base.leb \
  -o data/leb2/base_core.leb2 --group core

# Convert all groups at once
python scripts/generate_leb2.py convert-all data/leb/ephemeris_base.leb \
  -o data/leb2/ --tier-name base

# Generate from scratch (Skyfield → Chebyshev → compress)
python scripts/generate_leb2.py generate --tier base --group core \
  -o data/leb2/base_core.leb2

# Verify against LEB1
python scripts/generate_leb2.py verify data/leb2/base_core.leb2 \
  --reference data/leb/ephemeris_base.leb --samples 500
```

### 13.12 Measured Compression Results (Base Tier)

Measured from the shipped `data/leb2/base_*.leb2` files (chunk-index sums).
"Payload" is the Chebyshev coefficient data alone; the compression ratio is
payload-to-payload. Each on-disk group file additionally carries ~1.9 MB of
*uncompressed* auxiliary sections (nutation series, star catalog, ΔT table),
which is why a tiny group like Uranians produces a 2.1 MB file from a
0.7 MB payload: the auxiliary sections dominate.

| Group | Bodies | Payload (raw) | Payload (compressed) | Ratio | File on disk |
|-------|--------|---------------|----------------------|-------|--------------|
| Core | 14 | 44.5 MB | 8.7 MB | 5.1x | 10.6 MB |
| Asteroids | 5 | 23.0 MB | 6.7 MB | 3.4x | 8.7 MB |
| Apogee | 3 | 31.6 MB | 9.5 MB | 3.3x | 11.4 MB |
| Uranians | 9 | 0.7 MB | 0.1 MB | 6.2x | 2.1 MB |
| **Core-set subtotal** | **31** | **99.8 MB** | **25.0 MB** | **4.0x** | **32.8 MB** |
| Exotics | 31 | 199.7 MB | 56.9 MB | 3.5x | 59.0 MB |
| **All groups** | **62** | **299.5 MB** | **81.9 MB** | **3.7x** | **91.8 MB** |

---

## 14. Troubleshooting

### 14.1 "Body X not in LEB file"

The body is not one of the 62 bodies in `BODY_PARAMS`. LEB silently falls
back to Skyfield. This is expected behavior, not an error.

### 14.2 "JD outside range"

The requested date is outside the body's coverage range in the LEB file.
Falls back to Skyfield silently. Note that different bodies may have
different date ranges — asteroids typically cover ~1600-2500 CE (limited
by JPL Horizons SPK availability), while planets cover the full tier range.
Check `reader._bodies[body_id].jd_start` and `.jd_end` for per-body ranges.

### 14.3 Large Asteroid Errors (~1500")

The LEB file was generated without proper SPK coverage for asteroids.
Regenerate with:
```bash
export LIBEPHEMERIS_AUTO_SPK=1
leph leb generate base groups
```

### 14.4 "Failed to open LEB file"

- Check the file path exists and is readable
- Check the file is a valid LEB file (magic bytes = `b"LEB1"`)
- Check the file was not truncated during generation

### 14.5 Performance Not Improved

- Ensure `set_leb_file()` is called before any `calc_ut()` calls
- Check that `get_leb_reader()` returns a non-None value
- Verify the body you're computing is in the LEB file
- If using `FLG_ICRS`, LEB always falls back to Skyfield
  (`FLG_TOPOCTR`, `FLG_XYZ`, `FLG_RADIANS` and `FLG_NONUT` are handled
  on the LEB path)

---

## 15. Internals Deep-Dive

### 15.1 Outer Planet Position Computation

LibEphemeris supports three strategies for outer planet centers:

1. **Direct target** — Inner planets (Sun, Moon, Mercury, Venus, Earth) have
   direct segments in DE440. `planets["sun"]` works directly.

2. **SpkCenterTarget** — Outer planets use a separate `planet_centers.bsp`
   file with NAIF IDs 599 (Jupiter), 699 (Saturn), 799 (Uranus), 899
   (Neptune), 999 (Pluto). Position = barycenter + center offset from SPK.
   Precision: <0.001".

3. **CobCorrectedTarget** — Analytical Center-of-Body corrections from
   moon theory. Position = barycenter + analytical offset. Precision: <0.01".
   Used as fallback when `planet_centers.bsp` is unavailable.

The generator uses strategy 2 or 3 transparently via
`_eval_body_icrs_vectorized()`. For vectorized evaluation, the COB path
requires a scalar loop for `get_cob_offset()` (the offset function is
pure Python and not vectorizable), but the expensive barycenter evaluation
is still vectorized.

### 15.2 _PLANET_FALLBACK Map

```python
# planets.py:131-138
_PLANET_FALLBACK = {
    "mars": "mars barycenter",
    "jupiter": "jupiter barycenter",
    "saturn": "saturn barycenter",
    "uranus": "uranus barycenter",
    "neptune": "neptune barycenter",
    "pluto": "pluto barycenter",
}
```

Used when `planets[target_name]` raises `KeyError` (outer planets don't
have direct segments in DE440).

### 15.3 Asteroid NAIF IDs

```python
# generate_leb.py:243-249
_ASTEROID_NAIF = {
    15: 2060,      # Chiron
    17: 2000001,   # Ceres
    18: 2000002,   # Pallas
    19: 2000003,   # Juno
    20: 2000004,   # Vesta
}
```

SPK files from JPL Horizons use the `20000000 + N` convention for small
bodies, where N is the asteroid number. The `_spk_covers_range()` function
checks both conventions when searching for segments.

### 15.4 _SPK_BODY_MAP

```python
# state.py:135-137
_SPK_BODY_MAP: dict[int, tuple[str, int]] = {}
# Maps: body_id -> (spk_file_path, naif_id)
```

Populated by `auto_download_asteroid_spk()` and `download_and_register_spk()`.
The generator reads this to find the SPK file for each asteroid.

### 15.5 Precession-Nutation Matrix

`fast_calc.py` builds the ICRS -> true-equator-of-date rotation from the
**Vondrák 2011 long-term precession** combined with the nutation angles via
`vondrak_pn_matrix()` (`precession_vondrak.py`; valid ±200,000 years, unlike
the few-century IAU 2006 polynomials it replaced). The nutation angles come
from the LEB-stored Chebyshev coefficients on the pure-LEB path
(`_get_leb_frame_data`) and from Skyfield's IAU 2000A model on the Skyfield
path (`_get_skyfield_frame_data`) — both feed the same single Vondrák
matrix source, so the two backends cannot diverge at remote epochs.

```python
# fast_calc.py (LEB path)
dpsi, deps = <nutation from LEB Chebyshev coefficients>
pn_mat, eps_true_rad = vondrak_pn_matrix(jd_tt, dpsi, deps)
```

### 15.6 IAU 2006 General Precession (for Sidereal)

Used in both `_calc_ayanamsa_from_leb()` and the speed correction:

```python
# fast_calc.py
_PREC_COEFFS = (5028.796195, 1.1054348, 0.00007964, -0.000023857, -0.0000000383)
# arcsec/century polynomial: P(T) = sum(c_i * T^(i+1))
```

The sidereal speed correction subtracts the instantaneous precession rate:
```
dP/dT = c0 + 2*c1*T + 3*c2*T^2 + ...  (arcsec/century)
```
converted to deg/day and subtracted from dlon.

### 15.7 Aberration Formula

Classical first-order stellar aberration (`fast_calc.py:138-183`):

```
u = geo / |geo|                  # unit vector to body
v = earth_vel / c                # Earth velocity in units of c
u' = u + v - u*(u.v)            # aberrated unit vector
result = normalize(u') * |geo|   # scale back to original distance
```

This matches the pyswisseph implementation and provides ~0.01" accuracy
(sufficient for astrological purposes; the rigorous formula differs by
<1 milliarcsecond).

### 15.8 Full Segment Width Invariant

**Critical implementation detail:** All segments have the same width
(`interval_days`), including the last segment. The last segment may extend
beyond `jd_end` by up to `interval_days` worth. This is necessary because
the reader's O(1) segment lookup assumes uniform width:

```python
seg_idx = int((jd - body.jd_start) / body.interval_days)
```

If the last segment were truncated, the tau mapping would be wrong for all
dates in that segment. The generator handles this by:
- Using full-width segments for fitting (nodes extend beyond jd_end)
- Verifying only within `[seg_start, min(seg_end, jd_end)]`
- Using linear extrapolation for SPK overshoot

---

## Appendix A: File Size Estimation

```
For body with interval I, degree D, components C, over range R days:

segments = ceil(R / I)
bytes_per_segment = C * (D + 1) * 8
body_total = segments * bytes_per_segment + 52 (index entry)

Example: Moon, base tier (300yr = 109,573 days):
  segments = ceil(109573 / 4) = 27,394
  bytes/seg = 3 * 14 * 8 = 336
  total = 27,394 * 336 + 52 = 9,204,416 bytes (~8.8 MB)

Example: Uranus, base tier:
  segments = ceil(109573 / 64) = 1,713
  bytes/seg = 3 * 14 * 8 = 336
  total = 1,713 * 336 + 52 = 575,620 bytes (~0.6 MB)
```

## Appendix B: Adding a New Body to LEB

1. Add entry to `BODY_PARAMS` in `leb_format.py`:
   ```python
   NEW_BODY_ID: (interval_days, degree, coord_type, components),
   ```

2. Add name to `BODY_NAMES` in `generate_leb.py`:
   ```python
   NEW_BODY_ID: "Body Name",
   ```

3. Add evaluation function:
   - For ICRS: add to `_PLANET_MAP` or `_ASTEROID_NAIF`
   - For ecliptic: add to `eval_funcs` in `generate_body_ecliptic()`
   - For heliocentric: add to `generate_body_helio()`

4. Regenerate LEB files:
   ```bash
   leph leb generate base groups
   ```

5. Run precision tests:
   ```bash
   leph test leb-format precision
   ```

## Appendix C: Key Constants

```python
# leb_format.py
MAGIC = b"LEB1"
VERSION = 1
HEADER_SIZE = 64
SECTION_DIR_SIZE = 24
BODY_ENTRY_SIZE = 52
STAR_ENTRY_SIZE = 64
NUTATION_HEADER_SIZE = 40

# fast_calc.py
C_LIGHT_AU_DAY = 173.1446326846693
J2000 = 2451545.0
OBLIQUITY_J2000_DEG = 23.4392911

# generate_leb.py
NUTATION_INTERVAL = 16.0   # days
NUTATION_DEGREE = 16
DELTA_T_INTERVAL = 30.0    # days
N_VERIFY = 10              # verification points per segment
```
