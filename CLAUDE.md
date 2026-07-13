## Project Overview

LibEphemeris is a high-precision astronomical ephemeris library for Python, providing Swiss Ephemeris-compatible API using NASA JPL DE440 ephemeris via Skyfield. 

VERY IMPORTANT:
*Keep always 1:1 compatibility with PySwissEphemeris.*

## Commands

Uses `uv` for dependencies and `poe` (poethepoet) for task running.

The `leph` dev CLI is repo tooling and is not shipped in wheels: from a
source checkout use `./leph`, `python -m libephemeris.dev_cli`, or
`poe leph -- <args>`. No `leph` console entry point is installed.

```bash
uv pip install -e ".[dev]"        # Install with dev dependencies
poe lint                          # Ruff linter with auto-fix
poe format                        # Ruff formatter
poe typecheck                     # mypy
```

### Important: Never Run the Full Test Suite

**Do NOT run `poe test` or `poe test:full` during development.** Always run targeted tests:

```bash
pytest tests/test_file.py -v                          # Single file
pytest tests/test_file.py::TestClass::test_method -v  # Single test
pytest tests/test_file.py -k "pattern" -v             # By keyword within one file

# Skyfield backend (main unit test suite)
poe test:skyfield:core       # Essential ~900 tests, ~20s (parallel)
poe test:skyfield:fast       # All unit tests ~16,000, ~1 min (parallel, no @slow)

# LEB backend (precomputed Chebyshev ephemeris, ~14x faster)
poe test:leb:core            # LEB essential ~900 tests, ~20s (parallel)
poe test:leb:fast            # LEB unit tests ~16,000, ~1 min (parallel, no @slow) [RECOMMENDED]

# Feature-specific suites (via leph CLI)
./leph test lunar all          # All lunar tests (nodes, Lilith, perigee, apogee), no @slow
./leph test lunar perigee      # Perigee tests only
./leph test lunar apogee       # Apogee tests only
./leph test lunar lilith       # Lilith tests only (7 files)
./leph test leb2-format all    # LEB2 format unit tests (compression + reader)
```

## Code Style

- `from __future__ import annotations` at top of every module
- Imports grouped: stdlib, third-party, local (relative imports)
- Line length 88, Python 3.12+, double quotes, Ruff formatter
- Google-style docstrings with Args/Returns/Raises
- Naming: `snake_case` functions, `PascalCase` classes, `SCREAMING_SNAKE_CASE` constants, `_underscore` private
- Public API uses canonical bare names matching the upstream reference API (no `swe_`/`SE_`/`SEFLG_` prefixes). Single allowed exception: `SE_FNAME_DE431` (mirrored from upstream).
- Always return native Python floats (not numpy types)
- Exceptions in `libephemeris/exceptions.py`: `Error`, `CoordinateError`, `UnknownBodyError`, `EphemerisRangeError`, `PolarCircleError`

## Ephemeris File Selection

Three precision tiers: `base` (de440s.bsp, 1849-2150), `medium` (de440.bsp, 1550-2650, **default**), `extended` (de441.bsp, -13200 to +17191).

`FLG_MOSEPH` is accepted for API compatibility: it is echoed in the retflag like the reference (MOSEPH-only echoes MOSEPH), but the computation always uses JPL DE440/DE441 via Skyfield.

## LEB vs Skyfield Precision (Measured)

LEB Chebyshev approximation error vs Skyfield reference, per body group and tier:

| Body Group | Base (1860-2140) | Medium (1560-2640) | Extended (-5000 to +5000) |
|---|---|---|---|
| **Planets** (Sun-Pluto, Earth) | <0.001" | <0.001" | <0.001" |
| **Moon** | <0.001" (0.000332") | <0.001" (0.000325") | <0.001" |
| **Asteroids** (Chiron-Vesta) | <0.001" (0.000045") | <0.001" (0.000036") | <0.001" |
| **Ecliptic** (Nodes, Lilith) | <0.001" (0.000049") | <0.001" (0.000075") | <0.001" modern, up to 0.054" extreme dates |
| **IntpApog/IntpPerig** | Covered by the ecliptic tolerance | Covered by the ecliptic tolerance | Covered by the 0.1" extreme-date tolerance |
| **Uranians** | ~0.000000" | ~0.000000" | ~0.000000" |

**Precession at extreme dates**: the apparent-place reduction uses the **Vondrák 2011** long-term precession (valid ±200,000 years), via pyerfa — see `libephemeris/precession_vondrak.py`. This replaced the IAU 2006 precession (only valid a few centuries from J2000, ~36" off for the Sun at year -3000) and improves independently measured black-box agreement at extreme dates. Modern results are unchanged (Vondrák ≡ IAU 2006 to <1 mas near J2000).

**Known limitations at extreme dates**: (1) the nutation series still adds ~0.003" of degradation beyond ±2000 years (a physical limit, not a Chebyshev fit error; test floor 0.005" for extended extreme-date tests). (2) Versus pyswisseph specifically, the black-box reference solution differs from libephemeris DE441, creating an irreducible floor of tens-to-hundreds of arcsec for planets at deep-BCE dates that Vondrák does NOT remove (e.g. Mars ~600" at -3000). For the Sun this floor is small (~6" at -3000 once precession and ΔT are accounted for); no internal reference data source is inferred.

**Asteroid SPK coverage**: Safe range 1920-2080 CE. Outside this range, SPK data is unavailable and calculations use Keplerian fallback (catastrophically wrong for LEB compare tests). Test dates are filtered to this range.

## Binary Ephemeris Mode (LEB)

Precomputed `.leb` files with Chebyshev polynomial approximations (~14x speedup). Automatic fallback to Skyfield for unsupported bodies/flags — the only flag that always falls back is `FLG_ICRS`; `FLG_TOPOCTR`, `FLG_XYZ`, `FLG_RADIANS` and `FLG_NONUT` are handled on the LEB path (some body classes may still fall back per body).

```python
from libephemeris import set_leb_file
set_leb_file("/path/to/ephemeris.leb")  # or env LIBEPHEMERIS_LEB=...
```

Group workflow is **recommended** for generation (avoids macOS multiprocessing deadlocks):

```bash
poe leb:generate:base      # All 4 groups + merge for base tier
poe leb:generate:medium    # All 4 groups + merge for medium tier
poe leb:generate:extended  # All 4 groups + merge for extended tier
```

Runtime always uses a **single merged file** for LEB1. See `docs/leb/guide.md` for details.

## LEB2 Compressed Format

LEB2 uses error-bounded lossy compression (mantissa truncation + coeff-major
reorder + byte shuffle + zstd) to achieve 4-10x compression. The core
end-to-end gate requires <0.001" angular agreement with LEB1; the direct
verifier reports normalized per-component errors in native stored units
(AU for Cartesian ICRS; degrees/degrees/AU for ecliptic data), after checking
the declared inventory and per-body metadata against LEB1.

### Architecture

| Module | Purpose |
|--------|---------|
| `libephemeris/leb_compression.py` | Compression primitives: `compress_body()`, `decompress_body()`, `compute_mantissa_bits()` |
| `libephemeris/leb2_reader.py` | `LEB2Reader` — lazy per-body decompression, same interface as `LEBReader` |
| `libephemeris/leb_composite.py` | `CompositeLEBReader` — wraps multiple LEB files, dispatches by body_id |
| `libephemeris/leb_reader.py` | `open_leb()` factory auto-detects LEB1/LEB2 via magic bytes |
| `scripts/generate_leb2.py` | LEB2 conversion engine (invoked via `leph leb2`) |
| `scripts/test_leb2_precision.py` | Fast precision test: 14 core bodies x 6 flags x N dates per tier |

### Body Groups

| Group | Bodies | Base size |
|-------|--------|-----------|
| `core` | Sun-Pluto, Earth, Mean/True Node, Mean Apogee (14) | ~10.6 MB |
| `asteroids` | Chiron, Ceres, Pallas, Juno, Vesta (5) | ~8.7 MB |
| `apogee` | OscuApog, IntpApog, IntpPerig (3) | ~11.4 MB |
| `uranians` | Cupido-Transpluto (9) | ~2.1 MB |
| `exotics` | Centaurs, TNOs, NEAs (31) | ~59.0 MB |

### Per-body Precision Targets (`BODY_TARGET_AU` in `leb_compression.py`)

Moon/Earth use 1e-12 AU (not default 5e-9) because small geocentric distance amplifies errors through the pipeline (light-time, deflection, aberration). Inner planets use 1e-10 AU.

### Key Commands

```bash
poe leb2:convert:base              # Convert LEB1 -> LEB2 (all 5 groups)
./leph leb2 convert base-core      # Core group only (~10.6 MB)
./leph leb2 convert base-exotics   # Exotic registry group only
poe leb2:verify:base               # Verify base_core.leb2 against LEB1
./leph test leb2-format all             # Unit tests (compression + reader)
./leph test leb2-format precision-base  # Fast precision test
./leph test leb2-format precision-all   # Core companions, all tiers
```

### Full documentation
- `docs/leb/guide.md` — Complete LEB technical guide (section 13 for LEB2)
- `proposals/leb2-implementation-plan.md` — Implementation plan with benchmarks
- `release-notes/v1.0.0.md` — Release notes with measured results

## Horizons API Backend

Zero-install ephemeris via NASA JPL Horizons REST API. Used automatically in `"auto"` mode when no local DE440/LEB files are available, or explicitly via `set_calc_mode("horizons")`.

### Calculation Modes

| Mode | Flow | Fails when |
|------|------|-----------|
| `"auto"` (default) | LEB → Horizons (if no DE440) → Skyfield | never (always has fallback) |
| `"leb"` | Require LEB (auto-discovered or auto-downloaded if needed); unsupported bodies/flags fall back to Skyfield | no LEB resolvable |
| `"horizons"` | Prefer Horizons; unsupported bodies/flags fall back to Skyfield | no internet |
| `"skyfield"` | Always Skyfield/DE440 | DE440 not downloaded |

Set via `set_calc_mode()` or env var `LIBEPHEMERIS_MODE`.

### Architecture

| Module | Purpose |
|--------|---------|
| `libephemeris/horizons_backend.py` | `HorizonsClient` (LRU cache, parallel fetch, retry) + `horizons_calc_ut()` pipeline |
| `libephemeris/state.py` | `get_horizons_client()` — singleton with auto-detection |
| `libephemeris/planets.py` | Horizons dispatch between LEB and Skyfield paths |

### Supported Bodies
- Planets (Sun-Pluto, Earth): via Horizons VECTORS API
- Asteroids (Chiron, Ceres, Pallas, Juno, Vesta): via Horizons small-body syntax
- Mean Node, Mean Apogee: analytical (no HTTP)
- Uranians: analytical heliocentric (no HTTP)
- NOT supported: fixed stars, planetary moons, `FLG_TOPOCTR` -> fallback to Skyfield

### Full documentation
- `proposals/horizons-live-backend.md` — Original proposal with detailed design

## Lunar Calibration Workflow

The calibration/generation tooling (which fits coefficients against the reference
ephemeris as a black-box oracle) now lives in the **separate `validation/` repo**,
not in this library. Workflow:

1. From `validation/`: run the perigee calibration (`validation/calibrate/calibrate_perigee_perturbations.py`)
2. Paste coefficients into `_calc_elp2000_perigee_perturbations()` in `lunar.py`
3. From `validation/`: regenerate the residual table (`validation/calibrate/generate_lunar_apse_corrections.py --write`),
   which rewrites `lunar_apse_corrections.py` (the live residual table) in this repo
4. `./leph test lunar perigee` (the reference-free perigee tests in this repo)

See `docs/methodology/interpolated-perigee.md` for the full methodology.
