# PostgreSQL vs native LEB: Apple Container manual validation

Manual run on 2026-09-30, working tree on `feat/postgres-db-backend`.
These measurements are observations of this implementation, not a Railway SLA
or results for a different implementation. No implementation changes were needed
during this run, and no Swiss Ephemeris reference outputs were used or stored.

## Environment and isolation

Client: Apple M5 Pro, 18 logical CPUs, 48 GiB RAM, macOS 26.5.2,
Python 3.12.13, LibEphemeris 3.2.1, Psycopg 3.3.6.
Apple Container CLI 1.0.0. PostgreSQL 17.11, official `postgres:17` arm64 image,
**2 virtual CPUs / 2 GiB RAM**, `shared_buffers=256MB`, `max_connections=80`,
`track_io_timing=on`, and `pg_stat_statements` enabled. Normal durability settings
were retained; no `fsync=off` or similar benchmark shortcuts.

A newly named container, database and administrative/runtime roles were used.
Only a random loopback host port was published. No host ephemeris/data volume
was mounted, and the pre-existing `vector-global` container was not modified.
The host Python client reads native files only for import and the comparison
oracle; actual DB API calls do not need those files.

The runtime role received only schema USAGE and table SELECT. Runtime sessions
were additionally verified read-only, with the library's four-connection pool
and ten-second acquisition/statement limits.

## Real artifacts and import

| Artifact | Format | Segment rows incl. nutation | Decoded coefficient bytes | Import time |
|---|---|---:|---:|---:|
| `ephemeris_base.leb` | LEB1, full 53-body dataset | 1,032,414 | 374,311,728 | 7.11 s |
| `medium_core.leb2` | Chunked LEB2 | 508,493 | 170,150,528 | 5.57 s |
| `extended_core.leb2` | Chunked LEB2 | 14,048,478 | 4,700,863,552 | 167.65 s |

The bundled `base_core.leb2` was also imported as a separate dataset for its
own paired performance comparison. Different artifacts were **not** assumed
to contain identical coefficients. Each DB result was compared to its exact
source artifact. All datasets were explicitly published under distinct UUIDs.

The initial full-base database occupied 504,668,687 bytes. With all tiers,
indexes, retained sections, synthetic integration fixtures and auxiliary data,
`pg_database_size()` was 7,834,203,827 bytes (about 7.30 GiB; excludes cluster WAL).
The extended core source was 1,153,900,064 bytes: decoded coefficients alone are
about 4.07 times its compressed file size. Uncompressed per-segment storage has
an observable disk/index cost; this is not a compressed-container backup.

## Correctness and manual fault coverage

Streaming hashes compared **every payload byte**, in segment-index order per
series, for all three large datasets: **15,589,385 segments**, all matching.
This is not merely a sample of imported rows. Full-base hashing took 1.73 s,
medium 0.90 s, extended 21.23 s. Coefficients were decoded by the native artifact
parser and independently retrieved with binary PostgreSQL cursors.

Reader evaluation covered 1,350 full-base and 495 each medium/extended
body/nutation cases, including exact coverage endpoints, interval boundaries,
adjacent representable epochs and deterministic random epochs. All matched
native file evaluations exactly. Also compared all **1,447 full-base stars**
and five Delta-T interpolation/clamping epochs.

The public API matrix checked **6,454 cases**: 6,448 matching successes,
four matching `EphemerisRangeError` and two matching `UnknownBodyError` cases.
It included planets, lunar points, asteroids and fictitious bodies, `calc` and
`calc_ut`, SPEED, J2000, EQUATORIAL, XYZ, RADIANS, TRUEPOS, HELCTR, BARYCTR,
TOPOCTR, SIDEREAL, NONUT, NOABERR/NOGDEFL, invalid dates and missing bodies.
Ephemeris suffixes `.leb`, `.leb2`, `.bsp`, `.spk` were blocked at both
`builtins.open` and `io.open` during the DB matrix. DB routing and frame-owner
cleanup were asserted after every case.

An additional 16 API cases covered fixed stars/batching, orbital elements,
orbital distances, nodes/apsides, phenomena, planet-centric calculations,
eclipse maximum/circumstances, rise times and stellar ayanamshas. Another
400 requests used four differently configured topocentric/sidereal contexts
on eight threads, all exactly matching independently configured file contexts.

Medium added 99 public cases in years 1600, 2000 and 2600. Extended added
198 cases in years -12000, -5000, 0, 2000, 10000 and 17000. All succeeded with
zero numerical difference. These test years are samples, not a claim to test
every possible ancient/future API or scientific model.

Actual database fault injection verified:

- Idempotent reimport; duplicate-artifact import fails and rolls back.
- Published-row deletion is rejected even for the administrative role.
- The runtime role cannot write; explicit sealed policy blocks DB access.
- Truncated, NaN-containing and missing published segments raise `DBDataError`.
  Only inside this disposable database, the owner temporarily disabled user
  triggers to inject each fault, then restored the exact original bytes and
  re-enabled triggers. Every recovery result matched the original calculation;
  failed operation readers were closed and their maps cleared.
- Leasing all four pool connections makes an extra request fail in **10.003 s**.
- An administrative exclusive table lock makes a query fail in **10.002 s**.
  Both recover normally once the resource is released.
- Absent and unpublished UUIDs are rejected.
- A 1,000-segment scan never retains more than 256 input segments and clears
  inputs when closed.
- Stopping the actual container with SIGINT produces `DBError` on the existing
  connection, without fallback; restarting it recovers the same published UUID.

The existing `tests/test_db` selection also passed against this container:
**99 passed**. This was a targeted selection, not the full test suite.

## Latency benchmark

Final latency run was separate from pytest/fault injection. Client-side
`perf_counter_ns()` measures complete library calls, not just SQL execution.
Flags: `FLG_SWIEPH | FLG_SPEED`. Mars is the single target; charts are ten
independent `calc_ut` calls for Sun through Pluto. Date samples are deterministic
(`Random(9384)`, JD uniformly between 2450000 and 2460000).

Each scenario starts with fresh library state, then 25 warm-up calls. There are
300 single-target samples, 100 chart samples and ten first-request samples.
`reset_session()` runs outside the timed warm requests. The table shows median
and p95 in milliseconds from the final run; raw aggregates are in the sibling
JSON report.

| Workload | Full `.leb` | DB, same full dataset | Base `.leb2` | DB, same core dataset |
|---|---:|---:|---:|---:|
| Mars, same epoch: median / p95 | 0.022 / 0.024 | 0.921 / 1.116 | 0.021 / 0.035 | 1.034 / 1.249 |
| Mars, varying epoch | 0.167 / 0.353 | 1.048 / 1.345 | 0.145 / 0.171 | 0.951 / 1.543 |
| Ten-body chart, varying epoch | 0.905 / 1.441 | 9.385 / 11.941 | 0.689 / 0.948 | 9.560 / 11.923 |
| First request after `close()` | 1.664 / 1.876 | 10.045 / 10.857 | 2.439 / 2.678 | 9.099 / 9.349 |

Default-tier medium comparison: Mars varying epoch **0.147 ms file vs
0.741 ms DB** median (p95 0.168 vs 1.396 ms); chart **0.693 vs 8.851 ms**
(p95 0.847 vs 11.302 ms). Median overhead is about **5.0x** for Mars and
**12.8x** for a ten-body chart. Full-base medians are about **6.3x / 10.4x**;
base-core about **6.6x / 13.9x**.

File mode legitimately retains native scientific caches; DB mode intentionally
does not retain DB inputs/evaluations across operations. The particularly large
same-epoch gap is therefore a cache-policy comparison, not slow Clenshaw math.
Varying epochs are the more useful ordinary-request baseline.

An experimental explicit `begin_operation`/`finish_operation` around the chart
reuses transient inputs **only within that chart**: full-base median 4.713 ms,
base-core 3.894 ms, medium 4.231 ms. This is an internal ownership experiment,
**not a new public batch API**, and no such wrapping was added to production code.

Global eclipse searches from JD 2460310.5 returned identical results. Five-run
medians: solar **19.55 ms file / 42.24 ms DB**, lunar **10.58 / 27.39 ms**.
Their larger scoped operation amortizes reads; the standalone-planet overhead
ratio does not apply directly to every API.

## Query profile and concurrency

A measured ordinary Mars request made **two queries**: one metadata read and
one segment batch, 18 segment rows / 5,856 coefficient bytes. Ten independent
bodies made **20 queries**, 171 rows / 55,536 bytes. One explicitly scoped chart
made **9 queries**, 36 rows / 11,904 bytes. Read counts are concrete measurements
for these flags/epochs, not guarantees for every calculation.

An `EXPLAIN (ANALYZE, BUFFERS)` probe used `Nested Loop` plus composite
`segments_pkey` index scans, not a coefficient-table sequential scan. The local
`SELECT 1` median round trip was 0.125 ms. Probe timings are illustrative,
not averages from the latency benchmark.

Longer 10,000-request trials on one process and its shared four-connection pool:

| Threads | Requests/s | Median ms | p95 ms | p99 ms |
|---:|---:|---:|---:|---:|
| 1 | 1,156 | 0.815 | 1.066 | 1.176 |
| 4 | 2,545 | 1.533 | 2.171 | 2.484 |
| 16 | 2,414 | 6.492 | 7.215 | 8.470 |

No failures or leaked operation routes. Beyond four threads, queues increased
latency without increasing throughput. Four independently spawned processes
handled 8,000 requests at **3,904 requests/s including startup**; each created
its own pool and used one active connection for its serial workload. These are
short local load trials (seconds), not an hours-long production capacity test.

## Interpretation, reproduction and cleanup

The tested transport is numerically equivalent to its source artifacts and
fails closed under the injected faults. It currently buys centralized data and
volume-free replicas at a measurable latency/storage cost. It is **not faster
than local LEB**. Two sequential SQL round trips per ordinary planet, repeated
headers and repeated common segment transfer are important optimization targets;
scoped cross-body batching is promising without violating cache ownership.

Railway private-network latency, TLS, CPU throttling, concurrent applications,
server sizing and disk/cache pressure were not reproduced. Increasing RTT by
1 ms adds roughly 2 ms to an ordinary two-query request, or 20 ms to a chart's
20 independent queries, before queueing and computation effects. This is a
round-trip budget estimate, not a measured Railway benchmark.

First-request tests reset file handles/pools, **not OS page cache or PostgreSQL
shared buffers**. Python imports are already loaded. No system-wide cache purge
was used because that would affect unrelated services. The host was not a
dedicated benchmark machine, so small differences and tail percentiles should
not be treated as stable platform constants or extrapolated to other engines.

Manual harnesses and complete aggregate results were retained in the temporary
working directory reported by the session. To reproduce, start a fresh
`postgres:17` container with the above settings and a loopback port; provision
schema/import each explicitly selected artifact using `python -m libephemeris.db`,
create the SELECT-only role, then run the harness phases `manual.py validate`,
`faults.py`, `tiers.py`, `manual.py benchmark`, and `diagnostics.py` against that
**disposable** database. The fault script must never run on a deployed database.
Harness configuration supplies `admin_dsn`, `runtime_dsn`, `dataset`, `source`,
container `name`, and subsequently imported tier UUIDs. Passwords must be
injected anew; saved configuration was sanitized after completion.

After measurement, the dedicated container was stopped/deleted and its isolated
PostgreSQL storage discarded. No existing container/database was stopped or
removed. No ephemeris artifact was modified and no commit was made.
