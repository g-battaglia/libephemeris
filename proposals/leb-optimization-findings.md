# LEB Optimization: Findings & Known Issues

## Parameter Optimization Results

### Applied Changes (base tier)
| Body | Old (int/deg) | New (int/deg) | Savings |
|------|---------------|---------------|---------|
| Pluto | 32/13 | 64/11 | 0.63 MB |
| 9 Uranians (40-48) | 32/13 | 256/7 | 9.5 MB |
| **Total** | | | **~10.1 MB** |

Base file: 107 MB → 97 MB (-9.3%)

### Parameters That Do NOT Work (tested and failed)
| Body | Proposed | Fitting Error | Reason |
|------|----------|---------------|--------|
| Mercury | 32/13 | 1.4" | Eccentric orbit (e=0.206), high frequency |
| Venus | 32/11 | 15.9" | Needs short intervals for perturbations |
| Mars | 32/11 | 2.8" | Same reasons as Venus |
| Jupiter | 64/11 | 0.26" | Perturbations from Saturn are significant |
| Saturn | 64/11 | 0.044" | Mutual perturbations with Jupiter |
| Uranus | 128/9 | 0.026" | Degree 9 insufficient for perturbations |
| Neptune | 128/9 | 0.007" | Borderline, but exceeds worst-case target |
| Asteroids (16d) | 16/13 | 15-321" | Irregular orbits, strong perturbations |
| Mean Node | 128/7 | 0.82" | Precession adds complexity |
| True Node | 32/11 | 24.6" | High-frequency lunar perturbations |

### Key Insights
1. The precision margins reported by the accuracy sweep (e.g. "3800× margin for Jupiter") are misleading because the sweep does not test at moments of minimum geocentric distance (worst case).
2. The fitting errors of the ORIGINAL parameters already exceed 0.001" for some bodies (Jupiter 0.0015", Uranus 0.0025") — the documented target is not met at worst case even now.
3. Barycentric ICRS trajectories are NOT as smooth as hypothesized — planetary perturbations require short intervals and high degrees.
4. Ecliptic bodies (lunar nodes, apogees) have complex oscillations from lunar perturbations that require conservative parameters.

## Pre-existing Bugs Discovered

### 1. Geocentric Uranians not supported (bodies 40-47)
**File**: `libephemeris/planets.py:1854-1879`
**Problem**: `swe_calc(jd, 40, SEFLG_SPEED)` without `SEFLG_HELCTR` falls through to `raise UnknownBodyError`. The `if SE_CUPIDO <= ipl <= SE_POSEIDON:` block only has the `if is_helio:` path and does not handle the geocentric case (unlike Transpluto/body 48 which has geocentric conversion).
**Impact**: 16 test failures in the LEB test suite + all Uranians SKIP'd in the accuracy sweep.
**Severity**: Medium — Uranians are used almost exclusively in heliocentric mode in astrology, but the lack of geocentric support is an incompatibility with pyswisseph.

### 2. Heliocentric Sun (SEFLG_HELCTR) returns incorrect data
**File**: `tests/test_leb/compare/base/test_base_flags.py:85`
**Problem**: `swe_calc(jd, SE_SUN, SEFLG_SPEED | SEFLG_HELCTR)` produces an error of 646800" (~180 degrees).
**Impact**: 1 test failure in the flags suite.

### 3. TrueNode distance failure
**File**: `tests/test_leb/compare/base/test_base_lunar.py`
**Problem**: True Node distance exceeds the DISTANCE_AU tolerance.
**Impact**: 1 test failure.

---

## Read Speed Investigation / LEB Evaluation (2026-05-29)

Investigation of the question: "is a Rust module worth it to supercharge LEB file reading/opening?". Conclusion: **no for the real-world use case**; the actual margin is limited and elsewhere. All prototyped code was **discarded** after measurement — this document preserves the findings.

### What is already optimal (nothing to gain)
- **LEB file I/O**: already zero-copy `mmap` (`leb_reader.py`, `leb2_reader.py`), **singleton** reader opened once (`state.get_leb_reader`), coefficients read with `struct.unpack_from` (C). Opening/reading is NOT the bottleneck.
- **LEB2 decompression**: zstd (C) + numpy reorder (C), amortized and cached per body/chunk — does not happen per-query.

### Where the real bottleneck is
- The only pure-Python per-query loop is Chebyshev evaluation (Clenshaw), measured at **~0.65 µs** (docstrings claimed ~1.5 µs, optimistic by 2×).
- The full apparent pipeline (`fast_calc._pipeline_icrs`) is dominated by **per-date trigonometric transformations** (precession-nutation matrix, aberration, deflection), NOT by Chebyshev. So the single-date case (birth chart, ~15 bodies) has no removable Python hot loop.

### Rust: discarded
Even porting `eval_body` to Rust/PyO3, the estimated gain on full `calc_ut` was only ~1.2-1.4× (Amdahl: Chebyshev is half the cost, and ~8-11 FFI traversals/call remain), at the cost of breaking zero-toolchain `pip install` and introducing a wheel maintenance matrix. Not worth it.

### Measured prototypes, discarded
| Intervention | Measured real gain | Outcome |
|-------------|-------------------|---------|
| "Fused" Clenshaw (value + derivative in one recurrence) | ~2.26× on evaluation, **~2% on `calc_ut`** | Discarded: touches pyswiss-critical code for invisible gain. NB: it is also more *accurate* — vs exact truth (`fractions.Fraction`) mean error 1.47e-13 vs 2.06e-13 of the two-pass form. |
| `eval_body_many` (vectorized numpy Chebyshev, multi-date) | **7-8×** (Sun 7.5×, Moon 8.4×, Mars 7.2×), bit-identical | Discarded: no batch consumer in the project. Valid IF multi-date raw position scans are needed in the future. |
| `calc_ut_many` (pre-warm cache + real pipeline) | **0.59× (slower)** | Discarded: light-time and deflection look for *retarded* dates that are not pre-computable → cache hits don't cover the overhead. The apparent pipeline is not vectorizable without rewriting pyswiss-critical math. |

### Note on numpy in the scalar path
numpy on a **single point** is SLOWER than pure Python (call overhead > degree-13 loop): `chebval` 1 pt ~3 µs vs Python Clenshaw ~0.65 µs. numpy pays off **only** in batches (~0.03 µs/pt on 1000 points). The current choice of "pure Python, no numpy" in `leb_reader.py` is therefore already correct for single-call.

### Confirmed precision rule ("NASA grade")
Evaluate a numerical change against **exact mathematical truth** (`fractions.Fraction` / `mpmath`), not against the previous float implementation: the old code is also an approximation, not the reference.

### If batch speed is needed in the future
Build on a vectorized `eval_body_many` (raw Chebyshev positions, 7-8×, bit-identical) and apply only the necessary apparent corrections on top. Do NOT expect a generic `calc_ut_many` to be fast.
