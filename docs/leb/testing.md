# LEB vs Skyfield — Comparison Test Guide

> **Updated:** March 2026 — reflects <0.001" precision achieved for all 31 bodies.

This guide explains how to run and interpret the tests that compare LEB (LibEphemeris Binary) calculation results with Skyfield (direct calculation from NASA JPL ephemerides).

## Prerequisites

### LEB Files

Tests require pre-generated `.leb` files. There are two possible locations:

| Location | Path | Notes |
|----------|------|-------|
| **Local** (default) | `data/leb/ephemeris_{tier}.leb` | Used automatically by tests |
| **External** | `/Volumes/Data/libephemeris/leb/ephemeris_{tier}.leb` | Must be set via env var |

To use an external LEB file (e.g. on a separate disk):

```bash
export LIBEPHEMERIS_LEB=/Volumes/Data/libephemeris/leb/ephemeris_medium.leb
```

If the LEB file is not found, tests are **skipped** (they don't fail).

### Available Tiers

| Tier | Ephemeris | Range | File |
|------|-----------|-------|------|
| `base` | de440s.bsp | 1850–2150 | `ephemeris_base.leb` (~53 MB) |
| `medium` | de440.bsp | 1550–2650 | `ephemeris_medium.leb` (~175 MB) |
| `extended` | de441.bsp | -5000–+5000 | `ephemeris_extended.leb` (~1.6 GB) |

### Dependencies

```bash
uv pip install -e ".[dev]"
```

Asteroid tests require network access to automatically download SPK21 files from JPL Horizons (handled by `set_auto_spk_download(True)` in the test setup).

---

## Quick Commands

### Medium tier (default)

```bash
# All medium comparison tests (full, ~90 seconds with -n 4)
leph test leb-format vs-skyfield medium

# With parallelism (direct pytest)
LIBEPHEMERIS_LEB=/path/to/ephemeris_medium.leb \
  pytest tests/test_leb/compare/ -m "leb_compare" -v --tb=short -n 4
```

### Base tier

```bash
# All base tests (~90 seconds with -n 4)
leph test leb-format vs-skyfield base

# With parallelism (direct pytest)
pytest tests/test_leb/compare/base/ -m "leb_compare_base" -v --tb=short -n 4
```

### Extended tier

```bash
leph test leb-format vs-skyfield extended
```

### Single test

```bash
# A specific file
pytest tests/test_leb/compare/test_compare_leb_planets.py -v -n 4

# A specific test
pytest tests/test_leb/compare/test_compare_leb_planets.py::TestPlanetLongitude::test_longitude[6-Saturn] -v

# By keyword
pytest tests/test_leb/compare/ -m "leb_compare" -k "asteroid" -v -n 4
```

---

## Test Structure

### Directory

```
tests/test_leb/compare/
├── conftest.py                          # Shared infrastructure (tolerances, helpers, fixtures)
├── test_compare_leb_planets.py          # Lon, lat, dist, speed for ICRS planets
├── test_compare_leb_asteroids.py        # Position, speed, distance for asteroids
├── test_compare_leb_hypothetical.py     # Uranian bodies (Cupido, Hades, Zeus, ...)
├── test_compare_leb_velocities.py       # Speed lon/lat/dist for all 31 bodies
├── test_compare_leb_distances.py        # Geocentric and heliocentric distance
├── test_compare_leb_crossings.py        # cross_ut, solcross_ut, ...
├── test_compare_leb_eclipses_solar.py   # Solar eclipses
├── test_compare_leb_eclipses_lunar.py   # Lunar eclipses
├── test_compare_leb_nutation.py         # Nutation
├── test_compare_leb_deltat.py           # Delta-T
├── test_compare_leb_ayanamsha.py        # Ayanamsha (27 sidereal modes)
├── test_compare_leb_sidereal.py         # Sidereal positions
├── test_compare_leb_observations.py     # Equatorial coordinates, J2000
├── test_compare_leb_houses.py           # Houses (Placidus, Koch, ...)
├── test_compare_leb_rise_transit.py     # Rise, set, transit
├── test_compare_leb_stations.py         # Stations (retrogradation)
├── test_compare_leb_elongation.py       # Elongation
├── test_compare_leb_gauquelin.py        # Gauquelin sectors
├── test_compare_leb_lunar.py            # Lunar-specific functions
├── base/                                # Base tier tests (de440s)
│   ├── conftest.py                      # Base tier tolerances and fixtures
│   ├── test_base_planets.py
│   ├── test_base_asteroids.py
│   ├── test_base_velocities.py
│   ├── test_base_distances.py
│   ├── test_base_hypothetical.py
│   ├── test_base_sidereal.py
│   ├── test_base_flags.py
│   └── test_base_lunar.py
├── extended/                            # Extended tier tests (de441)
└── crosstier/                           # Cross-tier consistency tests
```

### Pytest Markers

| Marker | Meaning |
|--------|---------|
| `leb_compare` | Medium tier test |
| `leb_compare_base` | Base tier test |
| `leb_compare_extended` | Extended tier test |
| `leb_compare_crosstier` | Cross-tier test |
| `slow` | Intensive tests (100-200 dates per body) |

---

## How Comparison Works

### CompareHelper

The core of the infrastructure is the `CompareHelper` class in `conftest.py`. It runs the same function in two modes:

```python
# Skyfield mode (reference): direct calculation from NASA ephemerides
ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)

# LEB mode: calculation via precomputed Chebyshev polynomials
leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)
```

Internally, `CompareHelper`:
1. **Saves** the global libephemeris state (mode, tier, LEB file)
2. **skyfield()**: forces `set_calc_mode("skyfield")` and the correct precision tier
3. **leb()**: forces the specified LEB file and `set_calc_mode("auto")`
4. **teardown()**: restores the original state

### Typical Test Pattern

```python
@pytest.mark.leb_compare
@pytest.mark.slow
@pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
def test_longitude(self, compare, test_dates_200, body_id, body_name):
    max_err = 0.0
    worst_jd = 0.0

    for jd in test_dates_200:
        ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
        leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

        err = lon_error_arcsec(ref[0], leb[0])
        if err > max_err:
            max_err = err
            worst_jd = jd

    assert max_err < TOLS.POSITION_ARCSEC, (
        f'{body_name}: max lon error = {max_err:.4f}" at JD {worst_jd:.1f}'
    )
```

The test:
1. Iterates over N dates uniformly distributed across the tier range
2. For each date, calculates with Skyfield and with LEB
3. Measures the maximum error (worst case)
4. Compares against the tier tolerance

### Test Dates

Dates are generated uniformly across the tier range with a 30-day margin at the edges:

| Fixture | N dates | Usage |
|---------|---------|-------|
| `test_dates_200` | 200 | Planet and asteroid positions |
| `test_dates_100` | 100 | Velocities, distances |
| `test_dates_50` | 50 | Equatorial, sidereal |
| `test_dates_20` | 20 | Quick tests |

For the base tier: `base_dates_300`, `base_dates_150`, `base_dates_100`, `base_dates_50`.

### Asteroid Date Filtering

Asteroids have valid SPK data only for ~1900-2100 CE. Dates are filtered automatically:

```python
dates = filter_asteroid_dates(test_dates_200, body_id)
```

For non-asteroid bodies, all dates are returned unchanged.

---

## Tolerances

### Structure

Tolerances are defined in the `TierTolerances` dataclass and configured per tier in `TIER_DEFAULTS`:

```python
TOLS = TierTolerances.for_tier("medium")  # Loads medium tier defaults
```

### Override via Environment Variables

```bash
# Override for a specific tier
LEB_TOL_BASE_POSITION_ARCSEC=10.0 pytest ...

# Global override (fallback for all tiers)
LEB_TOL_POSITION_ARCSEC=10.0 pytest ...
```

Priority order (highest to lowest):
1. Explicit override via kwargs
2. `LEB_TOL_{TIER}_{FIELD}` (per-tier env var)
3. `LEB_TOL_{FIELD}` (global env var)
4. `TIER_DEFAULTS[tier]` (code defaults)
5. Dataclass defaults

### Current Tolerances

All 31 bodies achieve <0.001 arcsecond geocentric position precision on both
base and medium tiers.

#### Position

| Field | Base | Medium | Unit | Notes |
|-------|------|--------|------|-------|
| `POSITION_ARCSEC` | 0.001 | 0.001 | arcsec | All planets including outer |
| `ASTEROID_ARCSEC` | 0.001 | 0.001 | arcsec | |
| `ECLIPTIC_ARCSEC` | 0.001 | 0.001 | arcsec | Nodes, Lilith |
| `HYPOTHETICAL_ARCSEC` | 0.001 | 0.001 | arcsec | Uranians (error ~0) |
| `EQUATORIAL_ARCSEC` | 0.02 | 0.02 | arcsec | Heliocentric amplification |
| `J2000_ARCSEC` | 0.001 | 0.001 | arcsec | |
| `SIDEREAL_ARCSEC` | 0.001 | 0.001 | arcsec | |
| `DISTANCE_AU` | 5e-6 | 5e-6 | AU | |

#### Velocity

| Field | Base | Medium | Unit | Notes |
|-------|------|--------|------|-------|
| `SPEED_LON_DEG_DAY` | 0.045 | 0.045 | deg/day | OscuApogee dominates |
| `SPEED_LAT_DEG_DAY` | 0.004 | 0.004 | deg/day | |
| `SPEED_DIST_AU_DAY` | 1.2e-4 | 1e-4 | AU/day | |
| `ASTEROID_SPEED_LON_DEG_DAY` | 0.15 | 0.15 | deg/day | |
| `ASTEROID_SPEED_LAT_DEG_DAY` | 1.7 | 1.7 | deg/day | Architectural limit |
| `ASTEROID_SPEED_DIST_AU_DAY` | 5e-3 | 5e-3 | AU/day | |

#### Ecliptic Per-Body Overrides

| Body | Lon (arcsec) | Speed (deg/day) |
|------|-------------|-----------------|
| Mean Node | 0.001 | 0.0001 |
| True Node | 0.001 | 0.01 |
| Mean Apogee | 0.001 | 0.0001 |
| Oscu Apogee | 0.001 | 0.05 |
| Interp Apogee | 0.001 | 0.01 |
| Interp Perigee | 0.001 | 0.01 |

#### Timing (indirect functions)

| Field | Value | Unit | Tested function |
|-------|-------|------|-----------------|
| `CROSSING_SUN_SEC` | 1.0 | sec | `solcross_ut` |
| `CROSSING_MOON_SEC` | 5.0 | sec | `mooncross_ut` |
| `CROSSING_PLANET_SEC` | 30.0 | sec | `cross_ut` |
| `ECLIPSE_TIMING_SEC` | 1.0 | sec | `sol_eclipse_*` |
| `STATION_TIMING_SEC` | 1.0 | sec | Retrograde stations |
| `RISE_TRANSIT_SEC` | 1.0 | sec | `rise_trans` |

---

## Tested Bodies (31 total)

### Pipeline A — ICRS (11 bodies)

| ID | Name | Notes |
|----|------|-------|
| 0 | Sun | |
| 1 | Moon | |
| 2 | Mercury | |
| 3 | Venus | |
| 4 | Mars | |
| 5 | Jupiter | |
| 6 | Saturn | |
| 7 | Uranus | |
| 8 | Neptune | |
| 9 | Pluto | |
| 14 | Earth | |

### Pipeline A — Asteroids (5 bodies)

| ID | Name | Notes |
|----|------|-------|
| 15 | Chiron | SPK valid only 1900–2100 |
| 17 | Ceres | SPK valid only 1900–2100 |
| 18 | Pallas | SPK valid only 1900–2100, orbital inclination 34.8° |
| 19 | Juno | SPK valid only 1900–2100 |
| 20 | Vesta | SPK valid only 1900–2100 |

### Pipeline B — Ecliptic (6 bodies)

| ID | Name | Notes |
|----|------|-------|
| 10 | MeanNode | Error ~0 (analytical formula) |
| 11 | TrueNode | |
| 12 | MeanApogee | Error ~0 (analytical formula) |
| 13 | OscuApogee | Highest velocity error |
| 21 | InterpApogee | |
| 22 | InterpPerigee | |

### Pipeline B — Hypothetical/Uranian (9 bodies)

| ID | Name |
|----|------|
| 40 | Cupido |
| 41 | Hades |
| 42 | Zeus |
| 43 | Kronos |
| 44 | Apollon |
| 45 | Admetos |
| 46 | Vulkanus |
| 47 | Poseidon |
| 48 | Transpluto |

---

## LEB File Regeneration

If you modify Chebyshev parameters or the calculation pipeline, you must regenerate the LEB files.

### Group generation (recommended)

On macOS always use group generation (avoids multiprocessing deadlocks):

```bash
# Base tier
leph leb generate base groups

# Medium tier
leph leb generate medium groups

# Extended tier
leph leb generate extended groups
```

Each command runs in sequence:
1. `planets` — Sun-Pluto, Earth (11 bodies)
2. `asteroids` — Chiron, Ceres, Pallas, Juno, Vesta (5 bodies)
3. `analytical` — Nodes, Lilith, Uranians (15 bodies)
4. `merge` — Merges the 3 partial files + verification

### Single-body generation (lowest memory)

If group generation still uses too much memory, use single-body mode. Each of the 31 bodies is generated in its own subprocess (one at a time), then all partial files are merged:

```bash
# Base tier (direct CLI)
python scripts/generate_leb.py --tier base --single-body

# Medium tier
python scripts/generate_leb.py --tier medium --single-body

# Extended tier
python scripts/generate_leb.py --tier extended --single-body
```

This is slower than group mode but uses minimal memory (~1 body in memory at a time).

### Direct generation (Linux)

```bash
python scripts/generate_leb.py --tier base --verify
python scripts/generate_leb.py --tier medium --verify
python scripts/generate_leb.py --tier extended --verify
```

### Copy to external disk

After generation, the file is in `data/leb/`. To use it from tests with env var:

```bash
cp data/leb/ephemeris_medium.leb /Volumes/Data/libephemeris/leb/
```

---

## xfailed and Skipped Tests

### xfail (expected failures)

| Test | Reason |
|------|--------|
| Jupiter/Saturn geocentric crossing | Pre-existing bug in `crossing.py` solver (not LEB) |
| Mars 180° geocentric crossing | `RuntimeError: Maximum iterations reached` (not LEB) |
| Saturn 180°/270° heliocentric crossing | `RuntimeError: Heliocentric crossing search diverged` (not LEB) |

### skip

Tests are skipped if the LEB file for the corresponding tier is not found.

---

## Interpreting Failures

### Position error (arcsec)

```
AssertionError: Moon: max lon error = 0.0015" at JD 2451544.5
assert 0.0015 < 0.001
```

The error is in **arcseconds**. To convert:
- To degrees: divide by 3600 (0.0015" = 0.0000004°)
- To arcminutes: divide by 60 (0.0015" = 0.000025')

### Velocity error (deg/day)

```
AssertionError: Pallas: max lat speed error = 0.450000 deg/day at JD 2451544.5
assert 0.45 < 0.40
```

### Distance error (AU)

```
AssertionError: Pluto: max dist error = 4.50e-05 AU at JD 2451544.5
assert 4.5e-05 < 3e-05
```

1 AU = ~150 million km. An error of 3e-5 AU = ~4500 km.

### What to do if a test fails

1. **Check the date** — large errors at extreme dates (start/end of tier) suggest SPK contamination or edge Chebyshev segments
2. **Check the body** — Saturn, Uranus, asteroids have known architectural limits
3. **Widen the tolerance?** — only if the error is close to the limit and the safety margin (2x) is confirmed across many dates
4. **Regenerate the LEB?** — if you changed parameters in `leb_format.py` or `fast_calc.py`
