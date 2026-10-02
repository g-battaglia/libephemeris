# PostgreSQL coefficient backend

`db` is an explicit calculation mode, not a file cache or an extension of
`auto`. A published PostgreSQL dataset supplies the same native coefficients
that a LEB reader would evaluate. Mathematical reductions remain in the engine.
For explicit base/local and wider/DB selection, see
[tier routing and calculation sessions](tier-routing.md).

## Installation and configuration

Install the optional driver and pool:

```sh
pip install 'libephemeris[postgres]'
```

Select a published immutable version:

```python
import os
import libephemeris as ephemeris

ephemeris.set_db_config(
    os.environ["LIBEPHEMERIS_DB_URL"],
    "01234567-89ab-cdef-0123-456789abcdef",
)
ephemeris.set_calc_mode("db")
position, flags = ephemeris.calc_ut(2451545.0, ephemeris.MARS)
```

Alternatively configure `LIBEPHEMERIS_MODE=db`, `LIBEPHEMERIS_DB_URL` and
`LIBEPHEMERIS_DB_DATASET`. TOML keys under `[libephemeris]` are `mode`, `db_url`
and `db_dataset`. Explicit configuration overrides environment settings, which
in turn override TOML. `set_db_config(None)` clears the explicit override.
Do not commit passwords in TOML or shell history; prefer an injected secret.
`libephemeris config` redacts the entire `db_url` value before printing it,
including URI, query-parameter and libpq keyword connection-string forms.

Configuring a DSN does not change the calculation mode. `auto` never consults
PostgreSQL implicitly. `calc()`, `calc_ut()` and `EphemerisContext` preserve their
existing arguments, result shapes, flag normalization and model conventions.

## Explicit provisioning

Schema installation and import are offline administrative actions, never
implicit runtime migrations. Use a separate provisioning role:

```sh
# LIBEPHEMERIS_DB_URL points to the administrative connection for these steps.
python -m libephemeris.db schema
python -m libephemeris.db import \
    --dataset 01234567-89ab-cdef-0123-456789abcdef --tier base \
    /provisioning/base_core.leb2 /provisioning/base_asteroids.leb2
```

`scripts/import_leb_db.py` delegates to the same command. Commands also accept
`--dsn` before the subcommand, but injected environment secrets are preferable.
The schema is packaged at `libephemeris/db/schema_v1.sql`.

The importer supports LEB1 and both LEB2 versions. It preserves the decoded
binary64 values, without new truncation, fits or interpolation. Nutation,
Delta-T, stars, source headers/directories and non-coefficient sections are
retained as well. Reserved sections remain opaque. Original compressed body
blobs are replaced by decoded segments, so this is a scientific-data import,
not a byte-for-byte backup of the original container.

Files are explicit and belong to one tier. Duplicate bodies, mixed tiers,
conflicting auxiliary realizations and invalid coefficient inventories fail
closed. Companion auxiliary sections may repeat if they agree. All source
hashes and metadata are recorded in the manifest. Repeating an import with the
same UUID and identical manifest is idempotent; changing its artifacts is not.

Publication is atomic: failed imports roll back, and unpublished datasets are
not readable at runtime. After successful bulk publication, provisioning runs
`ANALYZE` on dataset/series/segment tables so a fresh large import is not served
with stale planner estimates. This is never a runtime write. Failed large COPY
transactions can leave reclaimable table/WAL space; rollback is not disk compaction. PostgreSQL guards published scientific rows against
ordinary INSERT/UPDATE/DELETE and validates complete segment inventories before
publication. A privileged database owner can of course disable triggers; this
is not protection against a malicious administrator.

Give the runtime role only USAGE on the `libephemeris` schema and SELECT on its
tables. The runtime pool additionally sets read-only sessions. Keep published
versions used by deployed instances; there is deliberately no mutable `latest`
alias or automatic garbage collection.

## Runtime ownership and efficiency

The only persistent backend resources are configuration and a lazy connection
pool (maximum four connections per process, ten-second acquisition/connect and
statement timeouts). Create workers before opening DB pools. An inherited pool
is rejected after fork rather than reused unsafely. Close/reconfigure while
workers are idle; an in-flight operation may fail if its pool is closed.

An outer operation owns its reader and inputs. Nested engine calls borrow that
owner. Inputs are discarded on both success and failure; repeating a request
performs fresh DB reads. No reader, coefficient or evaluated DB state enters
an application-global LRU, frame cache, vector-adapter cache, observer cache or
Besselian eclipse cache. Besselian geometry bypasses global cache reads/writes
in DB/routed mode; ordinary file caches are invalidated when sources change.
`close()` releases the pool while preserving explicit DB configuration;
`reset_session()` does not turn operation inputs into reusable state. The public
synchronous `calculation_session()` can group unchanged calls in one owner;
operation-local frame reductions are bounded and discarded with it. A caught
fatal storage error still prevents successful outer session exit.

Headers are read together. The reader prefetches segments for the target,
observer, gravitational deflectors and nutation using a paired-key batch query.
Neighboring segments cover many retarded/speed epochs; any additional epoch is
resolved by an exact indexed DB read. Prefetch is never a reason to extrapolate
or silently omit a required target state. Stored input buffers are bounded to
256 segments within an operation, so long searches do not accumulate a dataset.

The initial implementation uses one row per segment and uncompressed `bytea`
payloads. This trades disk space for small indexed reads and avoids decompressing
a ten-year chunk on every request. Transfer uses PostgreSQL binary results, not
hexadecimal JSON/text. No astronomical computation runs in SQL.

There is no guarantee of one query for every possible calculation, nor does a
loop of independent `calc_ut()` calls automatically become a public batch. The
storage contract already supports batched keys; explicit cross-body public
batching, temporal windows and microblock compression are future optimizations
to evaluate with actual workload measurements.

On Railway, place DB and replicas in the same region and use the private
connection address. Sum pool limits across all worker processes and replicas.
PostgreSQL's own shared buffers are centralized DB internals, not replica-local
application caches. Service replicas need no ephemeris volume.

## Source and network boundaries

`db` never discovers or opens `.leb`, `.leb2` or BSP files at runtime and never
falls back to Horizons/JPL when a DB operation fails. Local analytical models
already allowed by the LEB contract remain local models, not DB approximation
fallbacks. The integrated temporal models and fixed-star catalog keep their
existing rules; importing Delta-T/star sections does not silently replace them.
Normal package/code/model assets still exist: this mode eliminates ephemeris
file dependencies, not every filesystem read in Python and its dependencies.

The automatic network policy seals HTTP/downloads in DB mode while permitting
the explicitly configured PostgreSQL connection through a separate capability.
Explicit `network_policy=sealed` denies DB connections too. Explicit provisioning
commands select the allow policy. The connection address is not obtained from
an arbitrary fallback service.

Tracing reports `DB` for persisted DB states. Common body coverage reports DB
ranges with no invented filename. The file-oriented `get_leb_reader()` returns
None; `get_leb_inventory()` does not claim DB readiness. `get_current_file_data()`
returns its empty file record instead of reporting an old local kernel.

Transport failures use `DBError`; corrupt/missing published segments and invalid
metadata use `DBDataError`. Neither is a ValueError/KeyError fallback signal.
Ordinary missing-channel/range signals follow the existing reader protocol;
public source guards classify them into the existing public error hierarchy.
An absent optional nutation/Delta-T section is also the native reader's ordinary
missing-input signal, allowing only the already-declared local model. A missing
or corrupt segment in a declared series remains `DBDataError`, never that signal.
Connection credentials are not included in backend error messages.

## Explicit storage and numerical contracts

PostgreSQL-specific code lives in `libephemeris/db/`. Generic public scopes
and mixed-tier dispatch live in `operations.py` and `routing.py` respectively:

| File | Responsibility |
|---|---|
| `contract.py` | Data-only Series record, explicit validation/layout functions, store interface |
| `kernel.py` | Pure evaluation, epoch normalization, interpolation and read planning |
| `schema_v1.sql` | Versioned SQL schema and publication/immutability guards |
| `store.py` | Optional synchronous PostgreSQL transport, no astronomy |
| `reader.py` | Operation-owned adapter to the existing mathematical reader API |
| `backend.py` | Configuration, explicit operation lifetime and thin Python API shim |
| `policy.py` | Narrow DB transport authorization |
| `importer.py`, `__main__.py` | Explicit offline conversion/provisioning |

Storage methods are `metadata(dataset_id)`,
`fetch_segments(dataset_id, [(body_id, segment_index), ...])`,
`delta_t_points(dataset_id, jd)` and `star(dataset_id, star_id)`.
Records contain UUID strings, integers, float64 fields and bytes. There is no
ORM model, Python object serialization, callback or Skyfield object in the
storage contract. The small Python decorator is only a lifetime shim at existing
entry points; the backend has a data-only `DatabaseOperation` record passed to
`begin_operation(operation)` / `finish_operation(operation)`. Nested calls
borrow the outer owner, and finishing is idempotent. An owner must release
its inputs on every exit path, including failures. The synchronous compatibility
routing slot is thread-local, not part of the numerical kernel or an
asynchronous task abstraction.

### Data-only numerical boundary

`Series` is only a record: no methods, inheritance or reflection. Its free
functions receive all inputs explicitly. `kernel.py` has no classes,
decorators, driver objects, callbacks, I/O, caches or ambient state.

`DBReader` adapts those functions to the existing engine interface.
`PostgresStore` and `Artifact` bind the optional driver and native file parser.
Those adapters are integration machinery, separate from the numerical model.

The interfaces are explicit:

| Interface | Inputs and outputs |
|---|---|
| `Series` | Signed body ID, unsigned channel/count/degree fields, binary64 epochs |
| `segment_index(series, jd)` | Metadata and epoch → unsigned index, or typed error |
| `decode_segment(series, bytes)` | Metadata and little-endian bytes → finite coefficient sequence |
| `segment_coordinates(series, jd)` | Metadata and epoch → index and normalized epoch |
| `evaluate_body(series, coefficients, tau)` | Explicit inputs → three positions and three velocities |
| `evaluate_nutation(series, coefficients, tau)` | Explicit inputs → two radian angles |
| `interpolate_delta_t(points, jd)` | At most two ordered samples and epoch → time difference |
| `required_segment_keys(series_by_body, jd, body_id, flags)` | Metadata map and request → sorted body/index pairs |
| Operation-owned byte buffers | Segment-keyed byte map discarded on outer scope exit |
| `CoefficientStore` | Metadata, segment and auxiliary-data transport interface |

Evaluation functions require validated metadata/coefficients; every evaluator
must validate at these boundaries. The advertised positive span must reach the last
stored segment without extending past `segment_count * interval_days`; a
partial final segment is allowed. The upper comparison permits only the sum
of binary64 ULPs of the two endpoints and grid capacity, not a relative date
tolerance. Non-finite grid products or unreachable extra segments are invalid.
Lookup, decode and evaluation remain separate steps, so replacing a driver
never changes astronomy. Equivalent evaluators must preserve arithmetic order,
byte decoding, explicit error categories and operation-owned input lifetimes.

Series fields, in explicit SQL order:

```text
body_id: i32, coord_type: u32, segment_count: u32,
jd_start: f64, jd_end: f64, interval_days: f64,
degree: u32, components: u32
```

`body_id=-1` is reserved for nutation. Body coefficient channels retain their
native LEB coordinate types and units; nutation has two components in radians.
One payload contains `components * (degree + 1)` IEEE-754 binary64 values,
little-endian, ordered by component and then increasing polynomial degree.
Position and velocity are evaluated with the unchanged LEB Clenshaw recurrences.
Decode each eight-byte word explicitly as a little-endian binary64 value.

Index selection uses the actual stored interval:

```text
index = min(trunc((jd - jd_start) / interval_days), segment_count - 1)
midpoint = jd_start + index * interval_days + 0.5 * interval_days
tau = clamp(2 * (jd - midpoint) / interval_days, -1, 1)
```

Validate finite dates and inclusive coverage first. The final endpoint belongs
to the last segment. Preserve arithmetic order at rounded segment boundaries;
do not enable floating-point reassociation/FMA merely as a storage optimization.
Velocity has the derivative scale `2 / interval_days`. Angular longitude uses
Euclidean modulo, not signed remainder, for negative angles. Normalize into
`[0, 360)` and canonicalize zero to positive zero when comparing binary
representations; the direct kernel tests cover the negative-angle case.

`tests/test_db/fixtures/contract_v1.json` contains synthetic coefficients and
mathematical expected values that exercise the contract without an engine or
transport dependency. Any conforming reader can reuse the schema and existing
published datasets; no coefficient regeneration or fitting is required.

## Validation

```sh
uv run --no-sync pytest tests/test_db -m 'not db_integration' -q
# Only against an explicitly configured disposable database:
LIBEPHEMERIS_TEST_DB_URL='...' \
    uv run --no-sync pytest tests/test_db -m db_integration -q
```

The integration test installs the schema and publishes small synthetic versions
in that disposable database; it does not delete published versions afterward.
It skips when the explicit test DSN is absent. Unit tests compare native file
and DB-adapter behavior ephemerally, never by storing external reference outputs.
Do not run the entire project test suite for this change.
