# libephemeris-postgres

SPDX-License-Identifier: AGPL-3.0-only

`libephemeris-postgres` provides a small PostgreSQL-backed
`libephemeris.segment_source.SegmentSource`. It stores immutable LEB2 artifacts
as metadata and pages of little-endian float64 coefficients. The core library
continues to perform all numerical evaluation.

## Install

```bash
uv pip install \
  'libephemeris-postgres @ git+https://github.com/g-battaglia/libephemeris@TAG#subdirectory=packages/libephemeris-postgres'
```

For local development, install the package from this directory. Its
`[tool.uv.sources]` entry uses the repository root as an editable
`libephemeris` dependency; this is intentionally not a uv workspace.

## Database setup

Run the schema with an administrator URL:

```bash
export LIBEPHEMERIS_PG_ADMIN_URL='postgresql://...'
python -m libephemeris_postgres schema
```

Use separate roles. The importer role owns the schema. The runtime role needs
only `USAGE` on `libephemeris` and `SELECT` on its tables, and should have
`ALTER ROLE libephemeris_runtime SET default_transaction_read_only = on`.
No triggers are used: immutable datasets are published under a new UUID and
runtime permissions provide the immutability boundary.

## Import and verify

Each artifact commits independently, so an interrupted import can be resumed:

```bash
python -m libephemeris_postgres import --tier extended \
  extended_core.leb2 extended_asteroids.leb2 extended_exotics.leb2 extended_apogee.leb2
python -m libephemeris_postgres import --tier extended --resume DATASET_UUID FILES...
python -m libephemeris_postgres verify --dataset DATASET_UUID FILES...
```

The importer reads only public reader accessors, writes pages of 32 segments,
and runs `ANALYZE` after a complete dataset. A failed import may leave dead
space; run normal PostgreSQL vacuum maintenance before retrying a large load.
Verification compares checksums, float metadata, stars, and every page byte
for byte using a streaming reader.

## Runtime

Set `LIBEPHEMERIS_PG_URL` to a read-only runtime DSN and configure the core:

```bash
export LIBEPHEMERIS_PG_URL='postgresql://...'
export LIBEPHEMERIS_TIER_SOURCE_MEDIUM='libephemeris_postgres:open_tier'
```

Optional settings are `LIBEPHEMERIS_PG_DATASET_BASE`,
`LIBEPHEMERIS_PG_DATASET_MEDIUM`, and `LIBEPHEMERIS_PG_DATASET_EXTENDED`, plus
`LIBEPHEMERIS_PG_POOL_MAX` (default 2), `LIBEPHEMERIS_PG_TIMEOUT_SECONDS`
(default 5), and `LIBEPHEMERIS_PG_CACHE_PAGES` (default 1024). Pools are lazy
and process-local. The page cache is bounded; a miss prefetches the matching
page for every series in the artifact. Cache updates use a short lock;
concurrent cold misses may issue duplicate reads rather than serialize I/O.

Plan a connection budget as
`replicas * WEB_CONCURRENCY * (1 + active spawned subprocesses) * POOL_MAX`.
Keep that value below PostgreSQL `max_connections`, leaving room for admin and
maintenance sessions.
