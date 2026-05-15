## Project Overview

LibEphemeris is a high-precision astronomical ephemeris library for Python, providing Swiss Ephemeris-compatible API using NASA JPL DE440 ephemeris via Skyfield. 

VERY IMPORTANT:
*Keep always 1:1 compatibility with PySwissEphemeris.*

## Commands

Uses `uv` for dependencies and `poe` (poethepoet) for task running.

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
pytest -k "pattern" -v                                # By keyword

# Unit test suites
poe test:unit               # All unit tests (no slow), sequential
poe test:unit:fast           # All unit tests (no slow), parallel (-n auto)
poe test:unit:full           # All unit tests including slow

# Unit tests in LEB mode (precomputed Chebyshev ephemeris, ~14x faster)
poe test:unit:leb            # LEB mode, no slow, sequential
poe test:unit:leb:fast       # LEB mode, no slow, parallel (-n auto)
poe test:unit:leb:full       # LEB mode, including slow

# Feature-specific suites
poe test:lunar              # All lunar tests (no slow)
poe test:lunar:perigee      # Perigee tests only
poe test:lunar:apogee       # Apogee tests only
poe test:lunar:lilith       # Lilith tests only
poe test:leb                # LEB-specific unit tests (no slow)
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

`SEFLG_MOSEPH` is accepted for API compatibility but silently ignored. All calculations always use JPL DE440/DE441 via Skyfield.

## LEB vs Skyfield Precision (Measured)

LEB Chebyshev approximation error vs Skyfield reference, per body group and tier:

| Body Group | Base (1860-2140) | Medium (1560-2640) | Extended (-5000 to +5000) |
|---|---|---|---|
| **Planets** (Sun-Pluto, Earth) | <0.001" | <0.001" | <0.001" |
| **Moon** | <0.001" (0.000332") | <0.001" (0.000325") | <0.001" |
| **Asteroids** (Chiron-Vesta) | <0.001" (0.000045") | <0.001" (0.000036") | <0.001" |
| **Ecliptic** (Nodes, Lilith) | <0.001" (0.000049") | <0.001" (0.000075") | <0.001" modern, <0.005" extreme dates |
| **IntpApog/IntpPerig** | ~1-2° (pre-regen) | ~1-2° (pre-regen) | ~1-2° (pre-regen) |
| **Uranians** | ~0.000000" | ~0.000000" | ~0.000000" |

**Known limitation**: At extreme dates (beyond ±2000 years from J2000), Meeus nutation polynomial degradation adds ~0.003" to ecliptic body errors. This is a physical limit of the nutation series, not a Chebyshev fit error. Test tolerance floor: 0.005" for extended extreme-date tests only.

**Asteroid SPK coverage**: Safe range 1920-2080 CE. Outside this range, SPK data is unavailable and calculations use Keplerian fallback (catastrophically wrong for LEB compare tests). Test dates are filtered to this range.

## Binary Ephemeris Mode (LEB)

Precomputed `.leb` files with Chebyshev polynomial approximations (~14x speedup). Automatic fallback to Skyfield for unsupported bodies/flags (`SEFLG_TOPOCTR`, `SEFLG_XYZ`, `SEFLG_RADIANS`, `SEFLG_NONUT`).

```python
from libephemeris import set_leb_file
set_leb_file("/path/to/ephemeris.leb")  # or env LIBEPHEMERIS_LEB=...
```

Group workflow is **recommended** for generation (avoids macOS multiprocessing deadlocks):

```bash
poe leb:generate:base:groups      # All 3 groups + merge for base tier
poe leb:generate:medium:groups    # All 3 groups + merge for medium tier
poe leb:generate:extended:groups  # All 3 groups + merge for extended tier
```

Runtime always uses a **single merged file** for LEB1. See `docs/leb/guide.md` for details.

## LEB2 Compressed Format

LEB2 uses error-bounded lossy compression (mantissa truncation + coeff-major reorder + byte shuffle + zstd) to achieve 4-10x compression while maintaining <0.001" precision vs LEB1.

### Architecture

| Module | Purpose |
|--------|---------|
| `libephemeris/leb_compression.py` | Compression primitives: `compress_body()`, `decompress_body()`, `compute_mantissa_bits()` |
| `libephemeris/leb2_reader.py` | `LEB2Reader` — lazy per-body decompression, same interface as `LEBReader` |
| `libephemeris/leb_composite.py` | `CompositeLEBReader` — wraps multiple LEB files, dispatches by body_id |
| `libephemeris/leb_reader.py` | `open_leb()` factory auto-detects LEB1/LEB2 via magic bytes |
| `scripts/generate_leb2.py` | CLI: `convert`, `convert-all`, `generate`, `verify` |
| `scripts/test_leb2_precision.py` | Fast precision test: all bodies x 6 flags x N dates per tier |

### Body Groups

| Group | Bodies | Base size |
|-------|--------|-----------|
| `core` | Sun-Pluto, Earth, Mean/True Node, Mean Apogee (14) | ~10.6 MB |
| `asteroids` | Chiron, Ceres, Pallas, Juno, Vesta (5) | ~8.7 MB |
| `apogee` | OscuApog, IntpApog, IntpPerig (3) | ~11.4 MB |
| `uranians` | Cupido-Transpluto (9) | ~2.1 MB |

### Per-body Precision Targets (`BODY_TARGET_AU` in `leb_compression.py`)

Moon/Earth use 1e-12 AU (not default 5e-9) because small geocentric distance amplifies errors through the pipeline (light-time, deflection, aberration). Inner planets use 1e-10 AU.

### Key Commands

```bash
poe leb2:convert:base              # Convert LEB1 -> LEB2 (all groups)
poe leb2:convert:base:core         # Core group only
poe leb2:verify:base               # Verify against LEB1
poe test:leb2                      # Unit tests (27)
poe test:leb2:precision:base       # Fast precision test (~15s)
poe test:leb2:precision:all        # All tiers (~45s)
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
- NOT supported: fixed stars, planetary moons, SEFLG_TOPOCTR -> fallback to Skyfield

### Full documentation
- `proposals/horizons-live-backend.md` — Original proposal with detailed design

## Lunar Calibration Workflow

1. `poe calibrate-perigee` (or `poe calibrate-perigee:quick`)
2. Paste coefficients into `_calc_elp2000_perigee_perturbations()` in `lunar.py`
3. `poe generate-lunar-corrections` (regenerates `lunar_corrections.py`)
4. `poe test:lunar:perigee`

See `docs/methodology/interpolated-perigee.md` for the full methodology.
