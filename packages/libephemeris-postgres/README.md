# libephemeris-postgres

SPDX-License-Identifier: AGPL-3.0-only

Store original, pinned LEB2 files in PostgreSQL as 64 KiB byte blocks. The
native LibEphemeris reader parses and evaluates the same bytes as local files;
no astronomical metadata or decompressed coefficient tables are duplicated.

## Install and provision

This directory builds a separate distribution; PostgreSQL dependencies are
not part of the core wheel. For local development its uv source points to the
repository root, without a workspace.

```bash
uv pip install \
  'libephemeris-postgres @ git+https://github.com/g-battaglia/libephemeris@TAG#subdirectory=packages/libephemeris-postgres'
export LIBEPHEMERIS_PG_ADMIN_URL='postgresql://...'
python -m libephemeris_postgres upload medium_core.leb2 medium_asteroids.leb2 \
  medium_exotics.leb2 medium_apogee.leb2
```

Upload creates the schema automatically. Only names and SHA-256 hashes in the
installed LibEphemeris manifest are accepted. Use the matching feature revision
for core and provider until the byte-backed reader is released. Upload commits every 256 blocks (16 MiB); rerunning the same command
resumes missing blocks. A session advisory lock serializes uploads of the same
hash. Stored bytes are hashed through a server-side streaming cursor before
publication. Completed files are never modified by the uploader. An incomplete
upload whose existing blocks are corrupt is rejected, not published; an
administrator must repair that incomplete upload before retrying.

The schema uses `files` and `blocks`. Payloads retain LEB2 compression, with
TOAST recompression disabled. It is not a migration of the earlier experimental
coefficient-page schema: provision a fresh database for this format.

Use separate credentials: the uploader owns the schema; the runtime role needs
only `USAGE` on `libephemeris` and `SELECT` on both tables. Set that role's
`default_transaction_read_only = on`. Administrative access can still change
stored data; it must not modify published hashes.

## Runtime

```bash
export LIBEPHEMERIS_PG_URL='postgresql://...'
export LIBEPHEMERIS_LEB_SOURCE='libephemeris_postgres:open_reader'
export LIBEPHEMERIS_PG_TIERS='medium,extended'
```

The provider composes remote selected tiers and reviewed local unsourced tiers
through the current precision tier. Base typically remains local. Runtime
selection uses installed manifest hashes, not dataset UUIDs or mutable aliases.

Each process uses one lazy read-only connection, with a 5-second connection and
statement timeout. Forked children create their own connection. Reads serialize
on the connection; one slice batches its missing blocks. Each file caches at
most 256 blocks (16 MiB), clearing the cache when full. Closing a reader clears
its bytes, not the shared connection. Missing/truncated blocks, damaged chunks
and transport errors raise the provider's `CoefficientSourceError`, without
scientific fallback or credential-bearing diagnostics. `ping()` probes the
transport independently of cached bytes.

Budget one connection per application worker and active spawned subprocess,
leaving capacity for maintenance and provisioning.
