## Project Overview

LibEphemeris is a high-precision astronomical ephemeris library for Python, providing Swiss Ephemeris-compatible API using NASA JPL DE440 ephemeris via Skyfield. 

VERY IMPORTANT:
*Keep always 1:1 compatibility with PySwissEphemeris.*

## Development Policy

The project is licensed AGPL-3.0-only.

Development rules for new and modified code:

- Never inspect, retrieve, possess, read, translate, adapt, or copy Swiss
  Ephemeris source code, source comments, documentation prose, algorithms, or
  data files. Reference-distribution files must not enter this repository.
- The reference API may be used only for ephemeral behavioral comparison:
  call its public API and compare outputs. Never persist, fit, or commit
  those outputs.
- Run `uv run python scripts/check_provenance.py` for integrity checks.

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

| Body Group | Base (1860-2140) | Medium (1560-2640) | Extended (-13200 to +17191) |
|---|---|---|---|
| **Planets** (Sun-Pluto, Earth) | <0.001" | <0.001" | <0.001" |
| **Moon** | <0.001" (0.000332") | <0.001" (0.000325") | <0.001" |
| **Asteroids** (Chiron-Vesta) | <0.001" (0.000045") | <0.001" (0.000036") | <0.001" |
| **Ecliptic** (Nodes, Lilith) | <0.001" (0.000049") | <0.001" (0.000075") | <0.001" modern, up to 0.054" extreme dates |
| **IntpApog/IntpPerig** | Covered by the ecliptic tolerance | Covered by the ecliptic tolerance | Covered by the 0.1" extreme-date tolerance |
| **Uranians** | ~0.000000" | ~0.000000" | ~0.000000" |

The table compares like-for-like states (system barycenters for the outer
planets). In non-sealed `skyfield` mode the runtime additionally applies the
`planet_centers` offset for Saturn and Pluto (planet center vs barycenter,
up to ~0.05"/0.08" of geocentric longitude), so LEB-vs-Skyfield apparent
places differ by that physical offset for those two bodies.

**Precession at extreme dates**: the apparent-place reduction uses the
**Vondrák 2011** long-term precession through ERFA; see
`libephemeris/precession_vondrak.py`. Near J2000 it converges with the IAU 2006
model, while its published validity interval makes it suitable for the extended
DE441 tier.

**Known limitations at extreme dates**: nutation and Earth-orientation models
have narrower published validity intervals than DE441. Treat deep-time apparent
coordinates as model-dependent even when the underlying barycentric state is
within kernel coverage. Compatibility checks may use the public reference API
only ephemerally; no per-date outputs or fitted corrections belong in this
repository.

**Asteroid SPK coverage**: Safe range 1920-2080 CE. Outside this range, SPK data is unavailable and calculations use Keplerian fallback (catastrophically wrong for LEB compare tests). Test dates are filtered to this range.

## Binary Ephemeris Mode (LEB)

Precomputed `.leb` files with Chebyshev polynomial approximations (~14x speedup). In non-sealed `auto` mode unsupported bodies/flags can fall back to Skyfield, and `FLG_ICRS` is completed via the in-memory `LEBVectorEphemeris` adapter; `FLG_TOPOCTR`, `FLG_XYZ`, `FLG_RADIANS` and `FLG_NONUT` are handled on the LEB path. In sealed `leb` mode there is **no** fallback: a core body outside coverage raises `EphemerisRangeError`, an absent core body raises `UnknownBodyError`, and unsupported bodies (moons, fixed stars) fail explicitly rather than opening JPL/Skyfield.

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
| `core` | Sun-Pluto, Earth, Mean/True Node, Mean Apogee (14) | ~10.2 MB |
| `asteroids` | Chiron, Ceres, Pallas, Juno, Vesta (5) | ~2.15 MB |
| `apogee` | OscuApog, IntpApog, IntpPerig (3) | ~9.8 MB |
| `exotics` | Centaurs, TNOs, NEAs (31) | ~29.4 MB |
| `uranians` | Hamburg bodies Cupido-Poseidon (40-47, 8) | ~46 KB |

`uranians` is companion-only: fitted from the runtime Neely (1980) propagation
in `libephemeris.hypothetical` (never from legacy coefficients), generated from
the standalone `ephemeris_{tier}_uranians.leb` partial (merged mains carry no
fictitious IDs). Since 3.0.0rc15 the runtime trusts it on presence under a
`DATA_FILES` name; its SHA-256 is verified once, by the installer that writes
it. With a trusted reader active, the Uranian branch in `planets.py`
sources body and Earth positions from the LEB (same transform chain, ~4x
faster geocentric speed); `fast_calc` still rejects IDs 40-58 from persisted
channels. Other hypothetical bodies (48-58) remain runtime-only; each registry
row carries a source annotation checked by the integrity gate.

### Per-body Precision Targets (`BODY_TARGET_AU` in `leb_compression.py`)

Moon/Earth use 1e-12 AU (not default 5e-9) because small geocentric distance amplifies errors through the pipeline (light-time, deflection, aberration). Inner planets use 1e-10 AU. Uranians (40-47) use 1e-12 because their channels store degrees natively, and the default AU-calibrated target would allow ~5e-9 deg of angular error.

### Key Commands

```bash
poe leb2:convert:base              # Convert LEB1 -> LEB2 (all 5 groups)
./leph leb2 convert base-core      # Core group only (~10.2 MB)
./leph leb2 convert base-exotics   # Exotic registry group only
./leph leb2 convert base-uranians  # Hamburg companion (from the standalone partial)
poe leb2:verify:base               # Verify base_core.leb2 against LEB1
./leph leb2 verify base-uranians   # Verify the Hamburg companion (--group/--tier)
./leph test leb2-format all             # Unit tests (compression + reader)
./leph test leb2-format precision-base  # Fast precision test
./leph test leb2-format precision-all   # Core companions, all tiers

# Regenerate the uranians partial (companion-only, never merged)
python scripts/generate_leb.py --tier base --group uranians
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
| `"leb"` | **Sealed** — LEB only, no fallback to Skyfield/JPL/BSP/network. `FLG_ICRS` and the event paths are completed in-memory via `LEBVectorEphemeris`; a core range miss raises `EphemerisRangeError`, an absent core body raises `UnknownBodyError`; unsupported bodies (moons, fixed stars) fail explicitly | a canonical group is missing/corrupted, or the request is outside every LEB interval |
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

## Independent Lunar Model Workflow

Mean lunar points use ERFA's implementation of the IERS 2003 Delaunay
arguments and conventional inclined-orbit geometry. Interpolated apsides are
anchored to the actual JPL DE440 apsis passages: a Delaunay series fitted at
the passages plus residual tables, regenerated only by
`scripts/generate_lunar_apse_model.py` and pinned by SHA-256.

For lunar changes, document the JPL/IAU or primary-literature basis, run
the targeted lunar tests, and run
`uv run python scripts/check_provenance.py` before release.

See `docs/methodology/independence.md` and
`docs/methodology/interpolated-perigee.md` for the methodology.
