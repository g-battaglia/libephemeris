## Project Overview

LibEphemeris is a pure-Python astronomical ephemeris library with a Swiss Ephemeris-compatible API, powered by NASA JPL DE440/DE441.

**Always maintain 1:1 compatibility with PySwissEphemeris.**

## Commands

```bash
uv pip install -e ".[dev]"              # Install with dev deps

leph code lint                          # Ruff linter with auto-fix
leph code format                        # Ruff formatter
leph code typecheck                     # mypy

leph test skyfield essential            # ~490 tests, ~20s (fast sanity check)
leph test leb-backend unit-fast         # ~5890 tests, ~1 min [RECOMMENDED]
pytest tests/test_file.py -v            # Single file
pytest tests/test_file.py::test_name -v # Single test
```

**Never run the full test suite.** Always use `leph test <subgroup>` or `pytest` on specific files. Full CLI reference (dev, prod, poe shortcuts, TAB completion setup): see `CLI.md`.

## Code Style

- `from __future__ import annotations` at top of every module
- Line length 88, Python 3.12+, double quotes, Ruff formatter
- Google-style docstrings, `snake_case` functions, `PascalCase` classes, `_underscore` private
- `swe_` prefix for Swiss Ephemeris-compatible functions
- Always return native Python floats (not numpy types)

## Architecture

Three tiers: `base` (de440s, 1849-2150), `medium` (de440, 1550-2650, **default**), `extended` (de441, -13200 to +17191).

Four calc modes (`set_calc_mode()` / `LIBEPHEMERIS_MODE`): `auto` (default: LEB->Horizons->Skyfield), `leb`, `horizons`, `skyfield`.

LEB = precomputed Chebyshev polynomials (~14x speedup). Set via `set_leb_file()` or `LIBEPHEMERIS_LEB` env var.

## Performance API

These are libephemeris-only features, not part of the SwissEph-compatible API surface.

- **`reset_session()`** — lightweight state reset preserving file handles, LEB reader, Skyfield timescale, and LRU caches. Only resets topo, sidereal mode, angles cache. Use between independent calculations instead of `close()`. 1750x faster than `close()` for consecutive calcs.
- **`set_ephe_path()` idempotent** — no-op when path unchanged (avoids teardown).
- **LEB2 v2 chunked format** — 10-year temporal chunks per body. Reader auto-detects v1 (monolithic) vs v2 (chunked). Cold-start 33x faster (47ms vs 1568ms). Generate with `scripts/generate_leb2.py convert`.
- **Lazy mmap (no global `MADV_WILLNEED`)** — pages are loaded on demand by the kernel. Use `reader.warm(jd_start, jd_end)` for selective preloading of specific date ranges. Configurable via TOML keys `mmap_preload`, `mmap_preload_start`, `mmap_preload_end`.

Key files for performance work:
- `libephemeris/state.py` — `reset_session()`, `set_ephe_path()`
- `libephemeris/fast_calc.py` — LEB calculation pipeline
- `libephemeris/leb2_reader.py` — v1/v2 reader with chunked decompression
- `libephemeris/leb_format.py` — chunk index structs (ChunkEntry, CHUNK_*)
- `libephemeris/leb_compression.py` — `compress_body_chunked()`
- `scripts/generate_leb2.py` — LEB2 v2 writer
