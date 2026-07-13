# LEB2: Compressed Binary Ephemeris Format

**Status:** Implemented; historical proposal sanitized after the clean-room
review.

This document records the architectural intent of LEB2 without preserving
metadata for retired generated artifacts. Historical file sizes, hashes, body
counts, residual tables, release links, and download instructions were removed
because some early files included bodies or models that no longer meet the
project's provenance standard.

## Current clean-room boundary

- The wheel may contain only the reviewed `base_core.leb2` asset.
- Every other LEB1 or LEB2 file must be generated locally from NASA JPL data
  and independently sourced analytical models.
- Published LEB artifacts from earlier releases are retired and are not valid
  conversion or verification inputs.
- Compatibility-oracle output is never stored in LEB files, fixtures,
  coefficients, thresholds, or reports.
- Unsupported hypothetical-body IDs must not be materialized into a generated
  file. Harrington (ID 50) is the sole built-in hypothetical orbit and is
  traced to Harrington (1988).

## Motivation

LEB1 stores Chebyshev coefficients in a memory-mappable binary layout. LEB2
keeps the same evaluation model while compressing coefficient blocks so a
reviewed core can be shipped without requiring a separate data acquisition
step.

## Format design

LEB2 provides:

- a versioned header and section directory;
- a per-body index with explicit coverage and coordinate type;
- independently compressed coefficient blocks;
- bounded, lazy decompression into a thread-safe cache;
- unchanged Clenshaw evaluation after decompression;
- optional chunking so only the requested date interval is decompressed; and
- validation of magic, version, bounds, offsets, lengths, and checksums before
  coefficient data is evaluated.

Nutation, time-scale support data, and catalog sections follow their own
documented provenance and packaging rules. A generated file is acceptable only
when every included section can be reproduced from cited independent sources.

## Reader architecture

`open_leb()` detects LEB1 or LEB2 and returns a reader with the shared runtime
interface. The LEB2 reader memory-maps the container, validates its index, and
decompresses a body or temporal chunk on first access. Subsequent calls use the
same read-only coefficient buffer and evaluation path.

Mutable reader caches are protected by locks. Closing a reader invalidates its
mapping and cache; malformed or truncated input fails closed rather than
falling through as a routine out-of-range condition.

## Generation and conversion

Generation starts from the tier's NASA JPL kernel plus independently sourced
models. The local workflow fits coefficients, verifies mathematical and
source-based invariants, writes LEB1 or LEB2, and records provenance without
embedding public compatibility-oracle observations.

Conversion from LEB1 to LEB2 is allowed only for a locally generated LEB1 file
whose complete provenance is known. It changes serialization and compression,
not the scientific source of the coefficients.

## Verification requirements

- Round-trip and reader-equivalence tests use synthetic fixtures or locally
  generated, provenance-approved inputs.
- Compression error bounds derive from declared numerical targets, not from
  compatibility-oracle residual fitting.
- Packaging tests allow only the reviewed bundled core.
- The repository provenance gate must pass before any generated artifact is
  considered releasable.

See `docs/leb/guide.md` for the active local-generation workflow and
`NOTICE.md` for the controlling clean-room policy.
