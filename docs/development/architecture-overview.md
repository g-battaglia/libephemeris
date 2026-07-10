# Architecture Overview

Codebase metrics, performance bottleneck analysis, and future development vision
for LibEphemeris. For the implemented LEB binary ephemeris system, see the
[LEB Technical Guide](../leb/guide.md).

---

## Table of Contents

1. [Codebase Overview](#1-codebase-overview)
2. [Performance Bottleneck Analysis](#2-performance-bottleneck-analysis)
3. [Future Rust Port via PyO3](#3-future-rust-port-via-pyo3)
4. [Appendix: Full Module Inventory](#appendix-full-module-inventory)

---

## 1. Codebase Overview

### Project Metrics

| Metric                      | Value                                  |
|-----------------------------|----------------------------------------|
| Package version             | 3.0.0rc4                               |
| Python support              | 3.12 - 3.14                            |
| Library source files        | 70+ `.py` files in `libephemeris/`     |
| Total library LOC           | ~100,000 lines                         |
| Largest module              | `eclipse.py` (15,700+ lines)           |
| Unit test files             | ~150+ in `tests/`                      |
| Comparison test files       | ~120 in `compare_scripts/tests/`       |
| Public API entries (`__all__`)| 500+                                 |
| Constants defined           | ~300+ (`SUN`, `MOON`, `FLG_*`, `SIDM_*`, etc.) |
| Exception classes           | 14 (1 base + 4 categories + 9 specific)|
| Runtime dependencies        | 7                                      |

### Architecture Summary

LibEphemeris does NOT implement its own orbital integrator. It delegates
all positional astronomy to **Skyfield** (which reads JPL DE440/DE441 SPK
kernels via **jplephem**). The library provides:

- A Swiss Ephemeris-compatible Python API (1:1 with pyswisseph)
- 25 astrological house systems (26 codes)
- 47 sidereal/ayanamsha modes
- Solar/lunar eclipse search and contact point calculation
- Besselian element computation
- Lunar node and Lilith calculations (4 methods)
- Fixed star catalog (1,447 stars)
- Hypothetical/Uranian body support (8 Hamburg School + Transpluto)
- Minor body/asteroid support with auto-SPK download
- Atmospheric extinction and heliacal visibility models

### Runtime Dependencies

| Dependency     | Purpose                                              |
|----------------|------------------------------------------------------|
| `skyfield`     | Core ephemeris engine (JPL DE440 via Skyfield)       |
| `skyfield-data`| Skyfield data files                                  |
| `astroquery`   | Auto SPK kernel downloads from JPL Horizons          |
| `certifi`      | SSL certificates for HTTPS to JPL                    |
| `zstandard`    | Zstandard compression for LEB2 format                |
| `pyerfa`       | IAU ERFA nutation, obliquity, precession, aberration |
| `click`        | CLI framework (`libephemeris` and `leph` commands)   |

> Note: `spktype21` is vendored at `libephemeris/vendor/spktype21.py` (upstream unmaintained since 2018, NumPy 2.x fix applied) — it is not a separate package dependency.

### Module Breakdown by LOC

| Module              | LOC    | Domain                                  |
|---------------------|--------|-----------------------------------------|
| `eclipse.py`        | 15,700+ | Solar/lunar eclipses, occultations     |
| `planets.py`        | 6,160  | Core planetary calculations             |
| `fixed_stars.py`    | 6,127  | Star catalog, proper motion             |
| `houses.py`         | 5,374  | 25 house systems (26 codes)             |
| `hypothetical.py`   | 3,651  | Uranian planets, Transpluto             |
| `lunar.py`          | 3,033  | Lunar nodes, Lilith (4 methods)         |
| `minor_bodies.py`   | 2,649  | Asteroids, TNO, Kepler equation         |
| `heliacal.py`       | 2,414  | Heliacal rising/setting                 |
| `crossing.py`       | 1,914  | Sun/Moon/planet crossings               |
| `spk_auto.py`       | 1,879  | Automatic SPK download/cache            |
| `utils.py`          | 1,722  | Angular math, coordinate transforms     |
| `extinction.py`     | 1,623  | Atmospheric extinction model            |
| `state.py`          | 1,492  | Global state management                 |
| `true_node_terms.py`| 1,412  | Lunar node perturbation coefficients    |
| `time_utils.py`     | 1,323  | Julian day, Delta-T, calendar           |
| `spk.py`            | 1,304  | SPK kernel support                      |
| `iers_data.py`      | 1,286  | IERS Delta-T data                       |
| `__init__.py`       | 1,142  | Public API surface (~690 exports)       |
| `constants.py`      | 1,071  | All body IDs, FLG_*, SIDM_*, etc. constants     |
| `exceptions.py`     | 962    | Exception hierarchy                     |
| `context.py`        | 547    | Thread-safe EphemerisContext             |
| `profiling.py`      | 518    | Profiling infrastructure                |
| `astrometry.py`     | 511    | IAU precession/nutation/aberration      |
| `cache.py`          | 155    | LRU cache for nutation/obliquity        |

---

## 2. Performance Bottleneck Analysis

### 2.1 Current Performance Profile

The library's own code is **pure Python** — no Cython, mypyc, or compiled
extensions in the package source, so its hot paths run as interpreted Python.
(Its scientific dependencies — NumPy, pyerfa — are compiled, but the ephemeris
algorithms themselves are not, which is the relevant factor for the per-op cost
below.)

#### Cost per Operation

| Operation                       | Cost     | Why                                           |
|---------------------------------|----------|-----------------------------------------------|
| `calc_ut()` with `FLG_SPEED` | ~1ms  | 3 full Skyfield pipelines (9-12 Chebyshev evals + Python object overhead) |
| `calc_ut()` without speed   | ~350us   | 1 Skyfield pipeline                           |
| `houses()` (Placidus)       | ~500us   | Up to 50 trig iterations + GAST from Skyfield |
| `sol_eclipse_when_glob()`       | ~500ms-2s| 600-1050 Skyfield pipelines (iterative search)|
| Full chart (10 planets + houses) | ~12ms   | 30x `calc_ut` + 1x `houses`          |

### 2.2 Where CPU Time Is Spent (Priority Order)

#### 1. Skyfield `.observe()` calls (dominant cost, ~90%+ of total time)

Each `.observe()` involves light-time iteration (2-3 iterations of JPL
Chebyshev polynomial evaluation via jplephem). This is the single largest
cost factor in every operation.

#### 2. Central difference velocity (3x multiplier)

The velocity computation in `_calc_body()` (planets.py:1909-1910) calls
itself recursively twice for central difference:

```python
result_prev, _ = _calc_body(t_prev, ipl, flags_no_speed)
result_next, _ = _calc_body(t_next, ipl, flags_no_speed)
velocity = (result_next - result_prev) / (2 * dt)
```

This means every `calc_ut()` with `FLG_SPEED` costs **3x** the
position-only cost. Each of the 3 calls runs the full Skyfield pipeline:
`.at().observe().apparent()`.

**Total Skyfield calls per `calc_ut()` with FLG_SPEED:**

| Step                    | `.at()` | `.observe()` | `.apparent()` |
|-------------------------|---------|--------------|---------------|
| Main position           | 1       | 1            | 1             |
| Velocity t-dt (recursive)| 1     | 1            | 1             |
| Velocity t+dt (recursive)| 1     | 1            | 1             |
| **Total**               | **3**   | **3**        | **3**         |

Each `.observe()` internally iterates light-time (~2-3 iterations), so
the true count of JPL ephemeris segment evaluations is approximately
**9-12 Chebyshev polynomial evaluations** per single `calc_ut()` call.

#### 3. Eclipse search (combinatorial explosion)

`sol_eclipse_when_glob()` call breakdown:

| Component           | `calc_ut()` calls | Direct Skyfield pipelines |
|---------------------|-----------------------|---------------------------|
| Find New Moon       | ~42                   | ~126                      |
| Check eclipse       | ~3                    | ~18                       |
| Refine maximum      | 0                     | ~90-180                   |
| Contact times       | 0                     | ~369-729                  |
| **Total**           | **~45**               | **~600-1050**             |

#### 4. No position caching

When calculating a full chart (10 planets at the same JD), the Earth
observer position `earth.at(t)` is recomputed **30 times** (10 planets x
3 pipelines each) despite being identical every time.

#### 5. `erfa.nut06a()` (mitigated by cache)

Would be expensive (~0.02ms per call with 1365 nutation terms), but the
LRU cache in `cache.py` (128 entries) makes this negligible for same-JD
calculations. Cache hit: ~0.0001ms (200x speedup).

### 2.3 Existing Caching Infrastructure

| Cached function              | Cache size | Underlying call    | Speedup |
|------------------------------|------------|--------------------|---------| 
| `get_cached_nutation(jd_tt)` | 128 entries| `erfa.nut06a()`    | 200x    |
| `get_cached_obliquity(jd_tt)`| 128 entries| `erfa.obl06()` + nutation | 200x |

**What is NOT cached:**
- Planet positions (every `calc_ut()` hits Skyfield every time)
- Skyfield Time objects (created fresh each call)
- Earth position for observer (recomputed per planet)

### 2.4 numpy Usage in Hot Paths

All numpy usage in the hot paths is **scalar or small 3-element vector
operations** (e.g., `np.array(offset)` for tuple→array conversion). No
vectorized batch computation (`np.dot`, `np.matmul`, etc.) is used anywhere
in the critical path.

---

## 3. Future Rust Port via PyO3

> **Note:** The LEB-reader-only Rust port strategy is also documented in
> the [LEB Technical Guide](../leb/guide.md).

This section describes a vision for porting the full library (~7,800 lines)
to Rust, beyond just the LEB reader.

### 3.1 Strategy

Once the Python implementation is validated, port the **reader + pipeline**
to Rust. The same `.leb` file format works for both Python and Rust. The
generator stays in Python (it runs once, not performance-critical).

### 3.2 Rust Crate Structure

```
engine/
├── Cargo.toml
├── pyproject.toml          # maturin build config
├── src/
│   ├── lib.rs              # PyO3 module entry point
│   ├── leb.rs              # .leb reader, mmap, index        (~200 lines)
│   ├── chebyshev.rs        # Clenshaw eval + derivative      (~80 lines)
│   ├── calc.rs             # Main calc pipeline               (~400 lines)
│   ├── coords.rs           # ICRS/ecliptic/equatorial transforms (~300 lines)
│   ├── precession.rs       # IAU 2006 polynomial              (~150 lines)
│   ├── light_time.rs       # 3-iteration light-time           (~60 lines)
│   ├── aberration.rs       # Annual aberration                (~100 lines)
│   ├── ayanamsha.rs        # 47 sidereal systems              (~800 lines)
│   ├── houses/             # 25 house systems (26 codes)                 (~2000 lines)
│   │   ├── mod.rs
│   │   ├── placidus.rs
│   │   ├── koch.rs
│   │   └── ...
│   ├── eclipse.rs          # Eclipse search + contacts         (~1500 lines)
│   ├── crossing.rs         # Crossings, stations               (~500 lines)
│   ├── lunar.rs            # Lunar nodes, Lilith               (~400 lines)
│   ├── hypothetical.rs     # Kepler equation, hypotheticals    (~300 lines)
│   ├── stars.rs            # Fixed stars + proper motion       (~200 lines)
│   ├── constants.rs        # All body IDs, FLG_*                 (~400 lines)
│   ├── error.rs            # Error types                       (~100 lines)
│   └── pymodule.rs         # PyO3 bindings                     (~800 lines)
│
│   Total: ~7,800 lines Rust
```

### 3.3 Rust Dependencies

```toml
[dependencies]
pyo3 = { version = "0.23", features = ["extension-module"] }
memmap2 = "0.9"              # Memory-mapped .leb files
byteorder = "1.5"            # Binary parsing
thiserror = "2"              # Typed errors
parking_lot = "0.12"         # Fast Mutex/RwLock

[dev-dependencies]
approx = "0.5"               # Float comparison in tests
criterion = "0.5"            # Benchmarks
```

Zero math dependencies: all trigonometry uses `f64::sin()`, `f64::cos()`, etc.

### 3.4 Expected Speedup After Rust Port

| Operation                     | Python .leb | Rust .leb  | Total vs Skyfield |
|-------------------------------|-------------|------------|-------------------|
| `calc_ut()` + speed       | ~75us       | ~1-5us     | **200-1000x**     |
| `calc_ut()` mean node     | ~5us        | ~0.1us     | **1000x**         |
| `houses()` Placidus       | ~400us      | ~1-5us     | **100-500x**      |
| Eclipse search                | ~100ms      | ~0.5-2ms   | **500-2000x**     |
| Full chart (10 planets+houses)| ~1ms        | ~20-50us   | **240-600x**      |

### 3.5 Python Integration (Transparent Fallback)

```python
# libephemeris/planets.py
try:
    from libephemeris_engine import Ephemeris as _NativeEphemeris
    _engine = _NativeEphemeris("path/to/ephemeris.leb")
    _HAS_ENGINE = True
except ImportError:
    _HAS_ENGINE = False

def calc_ut(tjd_ut, ipl, iflag):
    if _HAS_ENGINE:
        return _engine.calc_ut(tjd_ut, ipl, iflag)
    reader = state.get_leb_reader()
    if reader is not None:
        return fast_calc.fast_calc_ut(reader, tjd_ut, ipl, iflag)
    return _calc_body_skyfield(...)  # original code, unchanged
```

Three tiers of performance:
1. **Rust engine** (~1-5us) - if `libephemeris_engine` is installed
2. **Python .leb reader** (~25-75us) - if .leb file is configured
3. **Skyfield** (~350-1050us) - fallback, always available

### 3.6 Rust Port Timeline (After Python Validation)

| Phase | Weeks | Deliverable                                        |
|-------|-------|----------------------------------------------------|
| 1     | 1-3   | `leb.rs` + `chebyshev.rs` + `calc.rs` (core pipeline) |
| 2     | 3-4   | `coords.rs` + `precession.rs` + `aberration.rs`   |
| 3     | 4-6   | `houses/` (all 25 systems) + `ayanamsha.rs`        |
| 4     | 6-9   | `eclipse.rs` + `crossing.rs`                       |
| 5     | 9-11  | `lunar.rs` + `stars.rs` + `hypothetical.rs`        |
| 6     | 11-13 | `pymodule.rs` (PyO3 bindings) + integration         |
| 7     | 13-15 | Test suite + validation vs pyswisseph + benchmarks  |

---

## Appendix: Full Module Inventory

### A.1 Coordinate Frames by Body Type

This table documents the coordinate frame produced by each body's
calculation function, as found in the source code analysis:

| Body                | Source Function                | Uses Skyfield? | Frame                  | Geocentric? |
|---------------------|--------------------------------|----------------|------------------------|-------------|
| Sun-Pluto           | Skyfield `.at().observe().apparent()` | Yes     | True ecliptic of date  | Yes         |
| `MEAN_NODE`      | `calc_mean_lunar_node()`       | No             | Mean ecliptic of date  | N/A         |
| `TRUE_NODE`      | `calc_true_lunar_node()`       | Yes            | True ecliptic of date  | Yes         |
| `MEAN_APOG`      | `calc_mean_lilith_with_latitude()` | No         | Mean ecliptic of date  | N/A         |
| `OSCU_APOG`      | `calc_true_lilith()`           | Yes            | True ecliptic of date  | Yes         |
| `INTP_APOG`      | `calc_interpolated_apogee()`   | No             | Mean ecliptic of date  | N/A         |
| `INTP_PERG`      | `calc_interpolated_perigee()`  | No             | Mean ecliptic of date  | N/A         |
| Uranians (40-47)     | `calc_uranian_planet()`        | No             | Heliocentric ecliptic  | Helio       |
| Transpluto (48)     | `calc_transpluto()`            | No             | Heliocentric ecliptic  | Helio       |
| Fixed stars         | `calc_fixed_star_position()`   | Yes            | True ecliptic of date  | Yes         |
| Chiron, Ceres-Vesta | SPK kernel via Skyfield        | Yes            | True ecliptic of date  | Yes         |

### A.2 Numerical Algorithms Used

| Algorithm                  | Where Used                            | Parameters                |
|----------------------------|---------------------------------------|---------------------------|
| Clenshaw (Chebyshev eval)  | .leb reader, JPL kernel evaluation    | degree 9-16               |
| Newton-Raphson             | `crossing.py`, `eclipse.py` (New Moon)| tol: 0.001-0.1 arcsec    |
| Brent's method             | `crossing.py` (root finding)          | tol: 1e-6 deg/day        |
| Bisection                  | `eclipse.py` (contact times)          | 60 iterations             |
| Golden section search      | `eclipse.py` (min separation)         | 60 iterations             |
| Fixed-point iteration      | `planets.py` (light-time, 3 iter)     | Convergence in 3 iter     |
| Placidus iteration         | `houses.py`                           | max 50 iter, tol 1e-7 deg|
| Gauquelin iteration        | `houses.py`                           | max 100 iter, tol 1e-8   |
| Pullen SR Newton-Raphson   | `houses.py`                           | max 20 iter, tol 1e-10   |
| Kepler equation            | `minor_bodies.py`, `hypothetical.py`  | Elliptic, hyperbolic, parabolic |
| Vondrák 2011 precession    | `precession_vondrak.py` (`erfa.ltp`/`ltpb`) | valid ±200,000 yr (IAU 2006 polynomial kept as no-erfa fallback) |
| IAU 2006/2000A nutation    | `cache.py` -> `erfa.nut06a()`         | 678 lunisolar + 687 planetary terms (1365 total) |

### A.3 House Systems Implemented (25)

| Code | Name                    | Method        | Iterative? | Polar safe? |
|------|-------------------------|---------------|------------|-------------|
| P    | Placidus                | Time-based    | Yes (50)   | No          |
| K    | Koch (GOH)              | Birthplace    | No         | No          |
| O    | Porphyry                | Arc trisection| No         | Yes         |
| R    | Regiomontanus           | Equator proj. | No         | Unstable    |
| C    | Campanus                | Prime vertical| No         | Unstable    |
| E/A  | Equal (from Asc)        | Geometric     | No         | Yes         |
| W    | Whole Sign              | Geometric     | No         | Yes         |
| X    | Meridian (Zariel)       | RA-based      | No         | Yes         |
| H    | Horizontal              | Azimuthal     | No         | Unstable    |
| T    | Polich-Page (Topocentric)| Modified Regio| No        | Unstable    |
| B    | Alcabitius              | Arabic        | No         | Unstable    |
| M    | Morinus                 | RA projection | No         | Yes         |
| U    | Krusinski-Pisa          | Multi-step    | No         | Unstable    |
| G    | Gauquelin (36 sectors)  | Iterative     | Yes (100)  | No          |
| V    | Vehlow Equal            | Geometric     | No         | Yes         |
| Y    | APC                     | Parallel circle| No        | Unstable    |
| F    | Carter Poli-Equatorial  | RA-based      | No         | Unstable    |
| N    | Natural Gradient        | Fixed (0=H1)  | No         | Yes         |
| S    | Sripati (Indian)        | Midpoints     | No         | Yes         |
| L    | Pullen SD               | Sinusoidal    | No         | Yes         |
| Q    | Pullen SR               | Newton-Raphson| Yes (20)   | Yes         |
| D    | Equal from MC           | Geometric     | No         | Yes         |
| I    | Sunshine (Treindl)      | Sun arc       | No         | Unstable    |
| i    | Sunshine (Makransky)    | Sun arc       | No         | Unstable    |
| J    | Savard-A                | Equator proj. | No         | Unstable    |

### A.4 Ayanamsha Systems (47)

All 47 systems are implemented in `planets.py::_calc_ayanamsa()` (~375 lines).
They use the IAU 2006 general precession polynomial (5th degree) as the
precession backbone, with system-specific reference epochs and offsets.

Categories:
- **Formula-based** (Lahiri, Fagan-Bradley, Krishnamurti, etc.): `aya = offset + precession_polynomial(T)`
- **Star-based** ("True" modes): compute apparent ecliptic longitude of reference star via Skyfield, subtract target position
- **Galactic-based**: calibrated formulas or galactic pole computation
- **User-defined** (`SIDM_USER`): `aya = ayan_t0 + [p(T_now) - p(T0)]`
