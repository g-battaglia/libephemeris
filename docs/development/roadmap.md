# Historical Roadmap and Status Record

This document preserves the development plan and completion record assembled in
March 2026.

> Last updated: March 2026

> **Historical snapshot — not the current release checklist.** Later releases
> changed test commands, counts, download manifests, module locations, and some
> fallback behavior. Completed items below describe the state when they were
> closed; open items must be revalidated against current code before being
> treated as planned work. Use the [Architecture Overview](architecture-overview.md),
> [Testing Guide](testing.md), and current issue tracker for present behavior.

## Table of Contents

1. [LEB: Per-Body Asteroid Coverage](#1-leb-per-body-asteroid-coverage)
2. [LEB: Generation Speed](#2-leb-generation-speed)
3. [Architecture: SPK+Kepler vs LEB Mode Decoupling](#3-architecture-spkkepler-vs-leb-mode-decoupling)
4. [Precision: Keplerian Fallback Improvement](#4-precision-keplerian-fallback-improvement)
5. [Open Items](#open-items)

---

## 1. LEB: Per-Body Asteroid Coverage

**Status: COMPLETED**

### Context

Asteroid SPK files (downloaded from JPL Horizons) cover a limited range (~1600-2500 CE), while planetary ephemerides (DE440/DE441) cover much wider ranges (DE441: -13200 to +17191). The LEB generator previously attempted to download SPK data for the entire tier range: for the "extended" tier (-5000 to +5000), it would request 10,000 years of coverage for Ceres — Horizons would refuse, and the asteroid would be silently excluded from the LEB file or fall back to Keplerian propagation (errors on the order of degrees).

### Solution

Each body in the LEB file has its own independent date range. Planets cover the full tier range; asteroids cover only the range actually supported by their SPK data (~1600-2500). At runtime, if the requested date falls outside an asteroid's range in the LEB, the library falls back to on-the-fly computation (SPK if available, otherwise Keplerian). The LEB binary format already supports per-body ranges — each body entry has independent `jd_start` and `jd_end` fields (`leb_format.py:BodyEntry`).

### Completed Subtasks

- **1.1 Discover maximum SPK range from Horizons** — Implemented `_get_asteroid_spk_range()` in `scripts/generate_leb.py` that opens SPK type 21 files, finds segments for the target NAIF ID (trying both 2M and 20M conventions), and returns the effective `(jd_start, jd_end)`.
- **1.2 Generate asteroids with independent range** — Modified `assemble_leb()` to discover effective SPK range, intersect with tier range (1-day margin), exclude asteroids with less than 20 years of useful coverage, and write per-body ranges into `BodyEntry`.
- **1.3 Update reader for per-body range** — Already implemented. The reader (`leb_reader.py:eval_body()`) uses per-body `jd_start`/`jd_end` fields and raises `ValueError` when out of range, caught by `calc_ut()`/`calc()` for Skyfield fallback.
- **1.4 Transparent fallback for out-of-range bodies** — Already implemented. The dispatch in `calc_ut()` and `calc()` catches both `KeyError` (body not in LEB) and `ValueError` (JD out of per-body range) for transparent per-body fallback.
- **1.5 Abort if SPK missing (strict mode for generator)** — Removed scalar Keplerian fallback from `generate_body_icrs_asteroid()`. The generator now raises `RuntimeError` if SPK is unavailable; `assemble_leb()` excludes asteroids without SPK with an `(EXCLUDED)` message.
- **1.6 Improve verify_leb() for asteroids and ecliptic bodies** — Rewrote
  `verify_leb()` for independently sourced body types: planets, JPL asteroids,
  and ecliptic bodies. The historical unsupported-hypothetical verification
  branch was later retired.
- **1.7 Poe commands for medium and extended with adaptive range** — Commands for group generation already exist and work automatically with per-body range since the logic is entirely in the generator.
- **1.8 Update documentation** — Documented per-body range concept, SPK coverage limits, fallback behavior, and effective range tables in `docs/leb/guide.md` and `README.md`.

---

## 2. LEB: Generation Speed

**Status: COMPLETED**

Vectorization completed for all 6 ecliptic bodies. A single Skyfield call `(moon - earth).at(t_array)` serves all bodies, followed by vectorized numpy computation for each.

**Batch functions implemented in `scripts/generate_leb.py`:**

- `_calc_mean_lilith_batch()` — vectorized `_calc_mean_apse_analytical()`
- `_calc_lunar_fundamental_arguments_batch()` — vectorized Delaunay arguments
- `_calc_elp2000_apogee_perturbations_batch()` — vectorized 40+ term series
- `_calc_elp2000_perigee_perturbations_batch()` — vectorized 61-term series
- `_eval_ecliptic_bodies_batch()` — core: single Skyfield call, 6 bodies
- `generate_ecliptic_bodies_vectorized()` — orchestrator
- `_fit_and_verify_from_values_unwrap()` — fitting with longitude unwrapping

**Measured times (base tier, 300 years):**

| Group | Before | After | Speedup |
|-------|--------|-------|---------|
| Planets (11 bodies) | ~15s | ~15s | — |
| Asteroids (5 bodies) | ~3-5 min | ~3-5 min | — (spktype21 limited) |
| Ecliptics (6 bodies) | ~5-8 min | **~38s** | **~10x** |
| **Supported groups total (historical measurement)** | **~10-15 min** | **~5 min** | **~2-3x** |

These historical timings predate the clean-room removal of unsupported
hypothetical and fitted lunar models. Current verification is source-based and
uses only independently derived data.

---

## 3. Architecture: SPK+Kepler vs LEB Mode Decoupling

**Status: COMPLETED**

### Context

The library has four calculation modes (default: `auto`):

1. **Auto mode** (default): tries a bundled, auto-discovered, or reviewed LEB2 installation, then Horizons API (if no local DE440), then Skyfield.
2. **Skyfield mode**: queries Skyfield in real time (SPK + frame rotation + nutation). Precise but slower (~120µs/call).
3. **LEB mode**: requires a valid locally generated LEB1 or reviewed LEB2 file (configured or auto-discovered). Calls served from precomputed Chebyshev coefficients (~2µs/call, ~14x speedup). Falls back to Skyfield for unsupported bodies/flags.
4. **Horizons mode**: prefers NASA JPL Horizons REST API. Falls back to Skyfield for unsupported bodies/flags.

### Goal

Make the library fully transparent without LEB. A user can use `libephemeris` by installing only the package and never touching `.leb` files, with the guarantee that all planets work (via Skyfield/DE440), all asteroids work (via auto-download SPK + Keplerian fallback), and there are no LEB-related errors or warnings.

### Completed Subtasks

- **3.1 Verify LEB-free mode is complete** — Already implemented. `_LEB_FILE` and `_LEB_READER` initialize to `None`; `get_leb_reader()` returns `None` when no LEB is configured; the `if reader is not None:` guard in `calc_ut()`/`calc()` skips the entire LEB block. All non-LEB tests exercise this path.
- **3.2 Document both modes** — Documented in `README.md` (new "Binary Ephemeris Mode (LEB)" section) and `docs/leb/guide.md` (Calculation Mode section with mode table, examples, env var documentation).
- **3.3 Environment variable for explicit mode** — Implemented `LIBEPHEMERIS_MODE` with four values: `auto` (default: LEB → Horizons → Skyfield), `skyfield` (force Skyfield), `leb` (require LEB with auto-discovery), `horizons` (prefer Horizons API). Added `set_calc_mode()`/`get_calc_mode()` in `state.py`, exported in `__init__.py`, reset in `close()`.
- **3.4 Graceful handling of missing LEB** — Already implemented. `get_leb_reader()` handles all edge cases: missing file (warning + `None`), corrupt file (warning + `None`), range too narrow (`ValueError` caught for per-body fallback), no file specified (`None` without warning).
- **3.5 Distribution of pre-generated LEB files** — The wheel bundles the
  reviewed, hash-pinned base core. `download_leb_for_tier()` and
  `download_leb2_for_tier()` install pinned base, medium, and extended core
  assets; locally generated LEB1 and LEB2 files remain available through
  explicit configuration.

---

## 4. Precision: Keplerian Fallback Improvement

**Status: MOSTLY COMPLETED**

### Context

The current Keplerian fallback (`minor_bodies.py`) uses osculating orbital elements at a fixed epoch (JD 2461000.5 = Sep 2025) with first-order secular perturbations from Jupiter, Saturn, Uranus, and Neptune. Only `omega`, `Omega`, and `n` are perturbed; `a`, `e`, `i` are constant. Includes resonant libration correction for plutinos (Ixion, Orcus).

Typical errors: main belt asteroids 10-30" over months, degrees over decades; TNOs 1-3' over months, degrees over years. Root cause: constant `a`, `e`, `i` + absence of short-period perturbations.

A hook for REBOUND/ASSIST already exists in `rebound_integration.py`, providing the most precise available method (sub-arcsecond), but requires optional dependencies (`pip install rebound assist`) and data files (~1 GB).

The current fallback cascade in `planets.py:1643-1728`:

```
1. SPK kernel (registered)            -> sub-arcsecond
2. Auto SPK download (from Horizons)  -> sub-arcsecond
3. Strict precision check             -> error if active
4. REBOUND/ASSIST (if installed)      -> sub-arcsecond
5. Keplerian with perturbations       -> degrees over decades (LAST RESORT)
```

### 4.1 Second-Order Secular Perturbations — COMPLETED

Implemented the complete Laplace-Lagrange secular theory (Murray & Dermott Ch. 7) with vectorial (h,k)/(p,q) evolution for eccentricity and inclination.

Functions added in `minor_bodies.py`:
- `_calc_forced_elements()` (~120 lines) — computes forced eccentricity vectors (h_f, k_f) and inclination vectors (p_f, q_f) plus proper frequencies (g, s) from Jupiter, Saturn, Uranus, and Neptune contributions using Laplace coefficients b_{3/2}^{(1)} and b_{3/2}^{(2)}.
- `apply_secular_perturbations()` now returns 6 values: `(omega_pert, Omega_pert, M_pert, n_pert, e_pert, i_pert)` with eccentricity oscillating around the forced value (proper period ~23,000 years for main belt asteroids).

Impact: minimal on 100-year benchmarks due to the ~23,000-year proper period, but significant on millennial timescales. The theory is correct and necessary for completeness.

### 4.2 Short-Period Perturbations — INVESTIGATED, DEPRIORITIZED

**Priority:** Low — the empirical approach does not work, and the analytical theory is too complex for the achievable gain. The SPK + LEB path is superior.

Investigation completed with two empirical approaches (Fourier fit of SPK - Keplerian residuals):
1. **Direct fit over 800 years** — RMS reduction only 2-7% for main belt asteroids. Residuals are dominated by secular drift (~degrees), not periodic oscillations (~arcminutes).
2. **Fit with polynomial detrending (degree 5)** — Removes secular drift before fitting. Fourier amplitudes are correct (~300" for Ceres, ~780" for Pallas) but pointwise validation shows catastrophic results: corrections worsen positions at short term (5" -> 19' at 1 month for Ceres).

The fundamental problem: Keplerian propagation with secular perturbations produces errors that grow as a power of time, not as stationary periodic oscillations. A static fit captures the mean amplitude, which is wrong for any specific date. For precision beyond the current level (arcminutes over decades), the only viable path is numerical integration (REBOUND/ASSIST). Generation script `scripts/generate_short_period_corrections.py` kept as reference, not integrated.

### 4.3 Multi-Epoch Orbital Elements — COMPLETED

Generated orbital elements from SPK type 21 state vectors at 50-year intervals (1650-2450 CE) for 6 bodies: Chiron, Pholus, Ceres, Pallas, Juno, Vesta (17 epochs each, ~600 lines in `MINOR_BODY_ELEMENTS_MULTI`). `_get_closest_epoch_elements()` selects the nearest epoch considering both the original single-epoch element and all multi-epoch entries.

**Benchmark results:**

| Offset | Before | After | Improvement |
|--------|--------|-------|-------------|
| epoch | 0.002" | 0.002" | — |
| 1 month | 7.4" | 7.5" | — |
| 25 years | 3.3 deg | **2.6'** | **~75x** |
| 50 years | 5.3 deg | 3.6 deg | ~1.5x |
| 100 years | 10.7 deg | 3.5 deg | ~3x |

The dramatic improvement at 25 years (~75x) is due to the multi-epoch table with entries every 50 years, placing the test point ~0 years from the nearest epoch.

### 4.4 REBOUND/ASSIST Integration as Intermediate Fallback — COMPLETED

The module `rebound_integration.py` (~890 lines) is complete and functional:
- `propagate_orbit_assist()` — ephemeris-quality integration with ASSIST
- `propagate_orbit_rebound()` — 2-body integration with REBOUND
- `propagate_trajectory()` — multi-point with automatic fallback
- `compare_with_keplerian()` — precision comparison

The hook in `planets.py:1677-1710` calls `check_assist_data_available()` and uses ASSIST before Keplerian when installed with data files. Completed: cached availability check, reset function, `AssistEphemConfig` search paths (`~/.libephemeris/assist/`, `ASSIST_DIR`, `./data/`), `download_assist_data()` helper, conditional end-to-end tests (47 REBOUND tests pass, 9 skipped without data files), macOS rpath fix.

Dependencies: `rebound>=4.0.0`, `assist>=1.1.0` (optional, in `pyproject.toml` as `nbody` extra).

### 4.5 Caching of REBOUND Results for Asteroids — DEPRIORITIZED

**Priority:** Low — performance optimization. The bundled reviewed base core
and locally generated LEB files cover the primary use case. REBOUND/ASSIST
fallback is adequate for sporadic queries.

### 4.6 Systematic Keplerian Precision Validation — COMPLETED

Created `tests/test_keplerian_precision_benchmark.py` with SPK comparison in ecliptic J2000 frame for 5 main asteroids (Ceres, Pallas, Juno, Vesta, Chiron) across 11 time intervals (epoch to 100 years). Regression tests: epoch < 1", 1 month < 60", 1 year < 5'.

---

## Open Items

Two items remain before this roadmap is fully closed:

### 1. ASSIST Data Files Download + End-to-End Verification

**Status:** Awaiting user approval

ASSIST data files (~714 MB total: `linux_p1550p2650.440` ~98 MB + `sb441-n16.bsp` ~616 MB) need to be downloaded and the full end-to-end path verified with real planetary perturbations. Nine conditional tests are currently skipped due to missing data files. Also pending: documentation for `pip install libephemeris[nbody]`.

### 2. Extended-tier LEB core — completed

The reviewed, SHA-256-pinned extended core is installable from the data
release. Optional extended companion groups still require local generation.

### 3. ASSIST Performance Verification for Single Evaluations

**Status:** Pending

ASSIST integration is slow for single evaluations (~ms per call). Performance characteristics need to be measured and documented to set user expectations. The REBOUND results caching approach (task 4.5) was deprioritized since the LEB system covers the primary use case.

### 4. Automatic Page Cache Release in `reset_session()`

**Status:** Proposed

**Problem:** libephemeris uses mmap for LEB file access (zero-copy,
~80x faster than Skyfield).  In cgroup v2 environments (Docker,
Railway, Kubernetes), file-backed page cache is counted in
`memory.current` and billed. The libephemeris's mmap design introduces
a container-specific regression: idle memory grows to ~1 GB+ even 
though the process heap is only ~220 MB.

**Current workaround:** Applications call `reader.cool()` +
`libephemeris.release_data_cache()` explicitly (added in v1.4.0).
The astrologer-api does this via middleware after each HTTP response.

**Proposed fix:** Call `cool()` automatically inside `reset_session()`
(`state.py`).  This is the natural cleanup point between calculation
sessions — kerykeion already calls it between charts.  No config,
no user action required.

**Considerations before implementing:**
- Interaction with `warm()`: warm pre-faults pages, then reset_session
  would mark them reclaimable.  Semantically OK because the kernel
  doesn't reclaim immediately — pages stay cached unless there's
  memory pressure.
- Batch processing: thousands of consecutive `calc_ut()` calls
  without reset_session would accumulate page cache.  This is the
  expected fast path — cool() should NOT run per-calculation.
- `release_data_cache()` walks the entire data directory — too heavy
  for reset_session.  Only `cool()` (one madvise syscall) belongs
  there.  Data directory cleanup remains application-side.

### 5. Port Eclipse/Fixed Stars/Heliacal Modules to LEB Path

**Status:** Proposed — analysis complete

**Problem:** `eclipse.py`, `fixed_stars.py`, and `heliacal.py` bypass
the LEB path and call Skyfield `get_planets()` directly.  When
invoked, Skyfield loads the DE kernel (up to 3.1 GB for DE441) into
memory.  The LEB files already contain all the position data these
modules need.

**Analysis document:** `docs/analysis/skyfield-to-leb-porting.md`

**Estimated effort:** fixed_stars (1-2 days), heliacal (2-3 days),
eclipse (2-4 weeks).  Requires implementing FLG_XYZ and
FLG_TOPOCTR in the LEB fast_calc path.
