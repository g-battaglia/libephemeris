# Test Suite Performance Analysis

> **Historical snapshot — not current test status.** This document records one
> diagnostic run from 2026-03-30 and preserves the bottleneck investigation
> that followed. Its command names, counts, failures, timings, and suite layout
> are not current release gates. Use the
> [Testing Guide](../development/testing.md) and `leph test --help` for today's
> supported commands, and rely on fresh command output for current status.

**Date**: 2026-03-30
**Run**: `poe test:full` (sequential, no `-n auto`)
**Duration**: 8h 48m 45s (31,725s)
**Results**: 17,463 passed, 43 failed, 107 skipped, 9 xfailed / 17,621 collected

---

## 1. Test Results: Are Things OK?

**Almost entirely OK.** The 43 failures break down as:

| Failure Group                    | Count | Cause                                                              | Status                                         |
| -------------------------------- | ----- | ------------------------------------------------------------------ | ---------------------------------------------- |
| `test_changelog.py`              | 1     | Invalid changelog section type                                     | Minor bug, non-critical                        |
| `test_context_thread_safety.py`  | 1     | Concurrent topocentric contexts                                    | Real bug, likely race condition                |
| `test_extended_asteroids.py`     | 25    | Asteroids outside SPK range (1920-2080) tested over -5000 to +5000 | **Known limitation** (documented in CLAUDE.md) |
| `test_compare_leb_asteroids.py`  | 15    | Same asteroid SPK range issue                                      | **Known limitation**                           |
| `test_heliocentric_positions.py` | 1     | Earth helio vs Sun geocentric mismatch                             | Precision bug                                  |

The **40 asteroid failures are expected** — JPL SPK21 coverage is 1920-2080 CE, but extended tier tests span -5000 to +5000. Outside SPK range, Keplerian fallback produces catastrophically wrong results. These tests should be filtered or skipped.

**4,664,832 warnings** were emitted (mostly Skyfield RuntimeWarning from `relativity.py` at extreme dates).

---

## 2. Performance Bottlenecks

### Bottleneck #1: `CompareHelper` Resets State on Every Call in Inner Loops

The single largest time sink. The LEB comparison tests iterate over hundreds of dates per test, calling both `compare.skyfield()` and `compare.leb()` at each step:

```python
# test_extended_lunar.py — runs for each of 6 ecliptic bodies × 6 test classes
for jd in ext_dates_500:  # 500 dates!
    ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
    leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)
```

Each `compare.skyfield()` call:

1. Sets `_LEB_FILE = None`, `_LEB_READER = None`
2. Calls `set_precision_tier(tier)` — may trigger BSP reload
3. Calls `set_calc_mode("skyfield")`
4. Executes full Skyfield computation (DE441 for extended tier)
5. Calls `set_calc_mode(None)` in finally block

Each `compare.leb()` call:

1. Sets `_LEB_FILE = path`, **`_LEB_READER = None`** — **forces LEB file re-open**
2. Calls `set_calc_mode("auto")`
3. Executes LEB computation (reader re-created from scratch)
4. Cleanup: `_LEB_FILE = None`, `_LEB_READER = None`, `set_calc_mode(None)`

**Impact on `test_extended_lunar.py` alone:**

- 6 bodies × 6 test classes × 500 dates × 2 calls = **36,000 `calc_ut` invocations**
- Each with full mode switching + LEB reader re-creation
- Estimated: **3-5 hours** for this single test file

The `_LEB_READER = None` on every call is the critical issue — it forces `open_leb()` (file open + header parse + index load) on every iteration of the inner loop.

### Bottleneck #2: `reset_ephemeris_state` Autouse Fixture Destroys Caches

This fixture runs **before and after every single test** (17,621 times):

```python
@pytest.fixture(autouse=True)
def reset_ephemeris_state():
    ephem.set_sid_mode(SIDM_FAGAN_BRADLEY)
    clear_caches()  # <-- NUKES ALL LRU CACHES
    yield
    # ... restore state
```

`clear_caches()` empties:

- Nutation LRU cache
- Obliquity LRU cache
- Time UT1 conversion LRU cache
- Time TT conversion LRU cache
- Observer-at cache

**Impact**: Sequential date computations within a test benefit from cache hits (same nutation/obliquity for nearby dates). But the cache is wiped between every test, destroying locality. For the comparison tests that iterate over 200-500 dates, the first few calls after a cache clear are significantly slower.

### Bottleneck #3: Skyfield DE441 Is Intrinsically Slow

Extended tier tests use DE441, which covers 10,000 years (-13200 to +17191). Each `calc_ut` via Skyfield involves:

1. JD → TT conversion (with Delta-T lookup)
2. SPK segment interpolation (DE441 has large segments)
3. Nutation computation (Meeus polynomial series)
4. Aberration, gravitational deflection corrections
5. Full ecliptic/equatorial coordinate pipeline

Estimated ~10-50ms per Skyfield call. With ~100,000+ total Skyfield calls across the suite, this accounts for hours.

### Bottleneck #4: No Parallelization

`poe test:full` runs bare `pytest` without `-n auto`. All 17,621 tests execute **sequentially on a single CPU core**. The machine has multiple cores available (98% CPU utilization on one core).

### Bottleneck #5: Redundant Computation in Test Structure

The comparison tests split each body into 6 separate test classes (longitude, latitude, distance, speed_lon, speed_lat, speed_dist). Each class iterates over the **same dates** calling the **same function** with the **same arguments**, but only checks one component of the 6-tuple result.

For `test_extended_lunar.py`:

- 6 classes × 6 bodies = 36 tests
- Each test calls `calc_ut()` 500×2 = 1,000 times
- But all 6 classes compute identical results — they just assert on different indices

A single test verifying all 6 components would reduce Skyfield calls by **6×**.

---

## 3. Estimated Time Distribution

| Test File/Area                             | # Tests    | Dates/Test | calc_ut Calls | Est. Time   |
| ------------------------------------------ | ---------- | ---------- | ----------------- | ----------- |
| `test_extended_lunar.py`                   | 36         | 500        | 36,000            | 3-5h        |
| `test_leb_precision.py`                    | ~2,400     | 150-200    | ~480,000          | 1-2h        |
| `test_extended_velocities.py`              | ~60        | 200        | 24,000            | 1-2h        |
| `test_extended_ancient/future/sidereal.py` | ~100       | 100-200    | ~40,000           | 30m-1h      |
| `test_extended_distances.py`               | ~35        | 200        | 14,000            | 15-30m      |
| `test_flag_pair_robustness.py`             | 940        | 5          | 9,400             | ~10m        |
| `test_flag_combinations_stress.py`         | 892        | varies     | ~5,000            | ~10m        |
| All other tests (~13,000)                  | 13,000     | 1-9        | ~50,000           | ~30m        |
| **Total**                                  | **17,621** |            | **~650,000+**     | **~8h 49m** |

The **LEB compare/extended tests** dominate: ~80% of total time in ~300 tests.

---

## 4. Recommended Fixes (Priority Order)

### P0: Fix `CompareHelper` to Keep LEB Reader Open

**Impact: 2-5× speedup on compare tests (saves ~4h)**

Stop resetting `_LEB_READER = None` on every call. The reader should stay open for the duration of the test:

```python
def leb(self, fn, *args, **kwargs):
    ephem.state._LEB_FILE = self.leb_path
    # DON'T set _LEB_READER = None — let it persist
    ephem.set_calc_mode("auto")
    try:
        return fn(*args, **kwargs)
    finally:
        ephem.set_calc_mode(None)
```

The reader cleanup should happen once in `teardown()`, not in every call.

### P1: Consolidate Test Classes to Avoid Redundant Computation

**Impact: ~6× fewer Skyfield calls for compare tests (saves ~3-5h)**

Instead of 6 separate test classes each iterating 500 dates:

```python
# BEFORE: 6 tests × 500 dates × 2 calls = 6,000 calc_ut per body
class TestLongitude: ...   # iterates 500 dates
class TestLatitude: ...    # iterates same 500 dates (redundant!)
class TestDistance: ...    # iterates same 500 dates (redundant!)
class TestSpeedLon: ...   # iterates same 500 dates (redundant!)
class TestSpeedLat: ...   # iterates same 500 dates (redundant!)
class TestSpeedDist: ...  # iterates same 500 dates (redundant!)
```

Consolidate into one test per body that checks all 6 components:

```python
# AFTER: 1 test × 500 dates × 2 calls = 1,000 calc_ut per body
class TestExtLunarPrecision:
    def test_all_components(self, compare, ext_dates_500, body_id, body_name):
        for jd in ext_dates_500:
            ref, _ = compare.skyfield(...)
            leb, _ = compare.leb(...)
            # Check all 6 components from same call
            assert lon_err < tol_lon
            assert lat_err < tol_lat
            assert dist_err < tol_dist
            assert speed_lon_err < tol_speed_lon
            assert speed_lat_err < tol_speed_lat
            assert speed_dist_err < tol_speed_dist
```

### P2: Run Tests in Parallel

**Impact: 4-8× speedup with `-n auto`**

Change `poe test:full` to use pytest-xdist:

```toml
"test:full" = "pytest -n auto"
```

This is already available (`pytest-xdist>=3.5.0` is installed) and used by `test:full:fast`.

### P3: Don't `clear_caches()` in `reset_ephemeris_state` by Default

**Impact: 1.5-2× speedup globally**

The `clear_caches()` call before every test destroys LRU cache locality. Options:

- Remove it from the autouse fixture; only tests that actually change ephemeris files need it.
- Or make the full reset a non-autouse fixture, applying it only to tests that modify global state.

### P4: Skip Asteroid Tests Outside SPK Range

**Impact: eliminates 40 known failures + wasted computation**

Add explicit skips or filter dates to 1920-2080 for asteroid bodies in extended tier tests. The current failures are expected and waste time computing meaningless results.

### P5: Reduce Date Counts for Extended Tier

**Impact: proportional time reduction**

500 dates per test for extended lunar tests is excessive for regression detection. 100-200 dates would catch the same precision issues with 2.5-5× fewer Skyfield calls.

---

## 5. Expected Outcome

| Fix                                | Time Saved        | Cumulative     |
| ---------------------------------- | ----------------- | -------------- |
| P0: Keep LEB reader open           | ~2-3h             | ~6h            |
| P1: Consolidate test classes       | ~3-5h             | ~2-3h          |
| P2: Parallel execution (`-n auto`) | 4-8× on remaining | **~20-40 min** |
| P3: Keep caches across tests       | ~15-30% further   | **~15-30 min** |

**Realistic target with P0+P1+P2**: from **8h 49m → ~20-40 minutes**.
