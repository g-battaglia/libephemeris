# Architecture Overview

LibEphemeris is a pure-Python compatibility layer and astronomical calculation
engine. It combines local NASA JPL ephemerides, a precomputed LEB fast path, a
remote Horizons backend, analytical models, and optional minor-body data behind
one canonical bare-name API.

## Public entry points and state

`libephemeris/__init__.py` defines the supported public surface in `__all__`.
Most applications use module-level functions such as `calc_ut()`, `houses()`,
and `set_sid_mode()`. Their mutable configuration lives in `state.py`.

`EphemerisContext` provides an isolated state container for concurrent or
independent calculations. Contexts share expensive immutable resources where
possible, but retain their own topocentric location and sidereal settings.

Two cleanup operations serve different purposes:

- `reset_session()` resets topocentric/sidereal session state without closing
  kernels, the LEB reader, the timescale, or reusable caches.
- `close()` performs full resource teardown and is appropriate when changing
  ephemeris resources or releasing file handles.

## Calculation dispatch

`set_calc_mode()` and `LIBEPHEMERIS_MODE` select one of four modes:

| Mode | Contract |
|---|---|
| `auto` | Prefer LEB, then use Horizons when the local DE kernel is absent, then Skyfield. |
| `leb` | Require a resolved LEB file; unsupported bodies or flags may use the general local pipeline. |
| `horizons` | Prefer NASA JPL Horizons; unsupported operations may use the local pipeline. |
| `skyfield` | Use the local DE440/DE441 pipeline through Skyfield. |

The selected mode is not necessarily the source of every derived quantity.
Nodes, lunar apsides, hypothetical bodies, houses, fixed stars, and coordinate
reductions have dedicated analytical or catalog-backed implementations. Use the
[tracing API](../guides/tracing.md) when an application needs to record which
backend supplied a body.

## Ephemeris tiers

Precision tiers coordinate the primary JPL kernel and supported date span:

| Tier | Primary kernel | Approximate range | Role |
|---|---|---|---|
| `base` | DE440s | 1849–2150 | Small, bundled LEB2 core and modern-date use. |
| `medium` | DE440 | 1550–2650 | Default general-purpose tier. |
| `extended` | DE441 | −13200–+17191 | Deep-history and far-future calculations. |

The wheel contains the reviewed `base_core.leb2`, covering 14 core bodies.
Reviewed medium and extended cores are distributed through exact SHA-256 pins;
other modular LEB1/LEB2 files must be generated and verified locally from JPL
kernels. JPL DE and planet-center kernels remain optional downloads from their
original source. Missing bodies continue through the normal fallback chain.

## Core data flow

For a typical planetary call, the engine:

1. converts UT to TT when needed using the selected Delta-T model;
2. resolves the requested body and calculation source;
3. obtains barycentric/geocentric state data from LEB, Horizons, or the local
   DE kernel;
4. applies light-time, aberration, deflection, body-center, observer, frame,
   nutation, and sidereal transformations required by the flags;
5. returns the reference-compatible six-value tuple and retflag.

LEB stores precomputed Chebyshev coefficients. Its reader supports the legacy
monolithic layout and the LEB2 v2 ten-year chunked layout. Chunks are
decompressed lazily; `warm()` can preload a selected date interval and
`cool()`/`release_data_cache()` can advise the operating system to reclaim
cached pages.

## Major modules

| Module | Responsibility |
|---|---|
| `planets.py` | Public planetary calculation pipeline and flag handling. |
| `state.py` | Global configuration, resource discovery, modes, and tiers. |
| `context.py` | Thread-safe context facade. |
| `fast_calc.py` | LEB calculation and reduction pipeline. |
| `leb_reader.py`, `leb2_reader.py` | LEB1/LEB2 readers and composite dispatch. |
| `horizons_backend.py` | Direct NASA JPL Horizons API backend. |
| `time_utils.py`, `astrometry.py` | Time scales and ERFA/Vondrák reductions. |
| `houses.py` | House cusps, angles, positions, and polar handling. |
| `eclipse.py`, `heliacal.py`, `crossing.py` | Event searches and visibility calculations. |
| `fixed_stars.py` | Catalog lookup, proper motion, and fixed-star positions. |
| `spk.py`, `spk_auto.py` | Minor-body SPK registration, direct download, and cache management. |
| `rebound_integration.py` | Optional REBOUND/ASSIST propagation. |

## Minor-body fallback

Minor bodies use the highest-quality available source:

1. a matching body in the active LEB reader;
2. a registered local SPK kernel;
3. direct SPK download from JPL Horizons when auto-download is enabled;
4. optional ASSIST N-body propagation when dependencies and data are present;
5. Keplerian propagation when strict precision permits it.

Strict precision prevents silent degradation for bodies expected to have SPK
coverage. REBOUND and ASSIST are an explicit GPL-licensed `[nbody]` opt-in and <!-- provenance-ok -->
are never part of the core or `[all]` installation.

## Shared astronomical models

- Planetary/lunar state vectors originate from NASA JPL DE440/DE441 or
  Horizons products.
- Nutation and frame operations use PyERFA, with the long-term Vondrák 2011
  precession/ecliptic geometry where remote epochs require it.
- Delta T is selected centrally and shared by positions and house angles.
- Outer-planet barycenter-to-body-center corrections use available JPL
  planet-center kernels and documented analytical satellite models.
- Fixed stars use the project's Hipparcos-derived catalog and proper motions.

See [Computational Methodology](../methodology/overview.md) for model choices and
[Independence and Methodology](../methodology/independence.md) for the
relationship to the reference implementation.

## Performance invariants

- `set_ephe_path()` is idempotent when the path is unchanged.
- `reset_session()` preserves expensive resources between independent jobs.
- LEB2 reads and decompresses only the body/date chunks requested.
- LRU caches hold reusable time, nutation, obliquity, and calculation data.
- Public results are normalized to native Python floats even when numerical
  dependencies use NumPy internally.

Performance numbers depend on hardware, file format, date locality, and flags;
benchmark them in the target environment rather than treating historical
microbenchmarks as architecture guarantees.
