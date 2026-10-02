# Mixed tier routing — Apple Container trial

The DB baseline was saved separately in `d7e77ee4`. This report measures the
subsequent uncommitted routing/session integration, not a production deployment.
Machine-readable results: [tier-routing-apple-container.json](tier-routing-apple-container.json).

## Setup and scope

Apple M5 Pro, 48 GiB, macOS 26.5.2; Python 3.12.13 for Kerykeion timings and
3.12.12 for API fault checks; Psycopg 3.3.6. PostgreSQL 17 ran in a dedicated
Apple Container with 2 CPUs/2 GiB and a loopback connection. Runtime access used
a SELECT-only role. No scientific volume was mounted into the DB container.
Existing unrelated containers and source artifacts were not modified.

Base used four reviewed local LEB2 groups. Medium imported core, asteroids,
exotics and apogee (~25.1 seconds). Extended imported **core, asteroids and
exotics only** (370.624 seconds): the four-group trial exceeded a 600-second
budget and rolled back. Extended apogee is therefore **not covered by the real
DB trial**. The paired file-tier reader deliberately used precisely the same
artifact sets, not a more complete reference distribution.

The disposable database occupied 38,257,907,379 bytes after the attempts. This
includes dead space from the rolled-back large import and is not a steady-state
storage estimate. Rollback prevents publication; it does not compact storage.

An initial post-import query hit the runtime statement timeout. Explicit
`ANALYZE` refreshed planner estimates and the complete rerun succeeded. The
provisioning importer now refreshes dataset/series/segment statistics after
successful bulk publication; calculations never perform maintenance writes.

## Correctness and source purity

**1,170 public calculation cases were bit-identical** to the paired native
file reader: major planets, Earth, nodes/apogee, Ceres/Chiron; epochs from
−12000 through +17000, real tier boundaries, speed, equatorial/Cartesian,
J2000/true-position and heliocentric/barycentric flags. Owner cleanup was
checked after each case. Additional local trials passed 138 NONUT/local-model
cases and four fixed-star flag combinations.

Kerykeion ten-point charts matched in 2000, 1600 and 1000, with actual `LEB`
for the modern chart and `DB` for the historical charts. Numerical chart fields
were compared independently of source/coverage metadata. Forty-day Mercury/Mars
station searches matched in modern, historical and base-boundary ranges; the
historical/boundary trials contained one/two stations respectively.

Modern calculations, modern charts, modern searches and base inventory produced
**zero DB queries**. One public scope used all three tiers, retained at most 256
remote segments in aggregate and cleared its owner on exit. Explicit `sealed`
policy allowed local base but rejected remote access.

Real API readiness passed base and complete medium, and correctly rejected
required extended because its apogee group was absent. A complete base catalog
cannot conceal an incomplete or unreviewed required remote source.

In the disposable DB only, segment triggers were temporarily disabled to inject
truncated and missing Mars input. Each fault produced `DBDataError`, prevented
a partial chart even after internal catches, cleaned ownership, left local base
usable and recovered exactly after byte-for-byte restoration. Actual container
shutdown left base readiness/charts usable while historical work returned a
sanitized 503 without fallback. Restart recovered exact results. The observed
warm-connection outage surfaced in ~1.6 ms; this is not a general timeout bound.

## Warm client-call timings

Charts use Kerykeion's existing lock/session and ten active points, with ecliptic
and equatorial computations. Each chart row has 100 samples over twenty repeated
dates; Mars has 300 varying dates. Historical station searches have ten samples
with slightly shifted forty-day windows. Pools/file caches were warmed; these
are complete client-call timings, not SQL-only or evaluator-only numbers.

| Workload | File median / p95 (ms) | Routed median / p95 (ms) | Routed queries | Coefficient bytes |
|---|---:|---:|---:|---:|
| Chart, 2000 | 1.751 / 1.927 | 2.758 / 2.980 | 0 | 0 |
| Chart, 1600 | 1.825 / 2.076 | 9.302 / 10.189 | 10.05 | 11,920.8 |
| Chart, 1000 | 1.944 / 2.072 | 10.444 / 11.219 | 12.05 | 11,984.8 |
| Mars, modern | 0.197 / 0.226 | 0.241 / 0.261 | 0 | 0 |
| Station search, 1600 | 8.170 / 8.260 | 18.578 / 18.887 | 12 | 13,360 |

Queries/bytes are average per chart/search, with bytes counting returned
coefficient payloads only, not metadata, PostgreSQL protocol or TLS overhead.
The extra queries near some dates are exact support-epoch reads.

Pure `db` against the same medium import measured **7.625 / 8.316 ms** for a
chart, 9.05 queries and the same average 11,920.8 coefficient bytes. Routed mode
also verifies tier/publication identity and considers higher-priority local
coverage. Operation-local DB frames are reused inside the chart/search and
cleared afterward; no scientific DB cache survives into the next request.

Local routing is not free: modern Mars median increased ~23%, and chart median
~58%, despite no remote contact. Mixed proxies remain non-cacheable, while
concrete file-only frame caches retain their existing behavior. These measured
costs should not be described as equal to pure LEB performance.

## Limits and reproducibility

No Railway, TLS, other-engine, full ASGI load/capacity trial or extended-apogee DB parity
is claimed. This is a short, nondedicated-host trial using modified local builds,
not the API's released dependency pins. The API source-contract tests also ran
against locally built wheels in an isolated interpreter, without synchronizing
shared environments or changing deployment pins.

The retained temporary harness contains native parity, chart/search timing and
disposable fault scripts. Credentials and the dedicated container are removed
after verification; reproduction needs a fresh disposable database. Never run
fault injection against a deployed dataset. Provision all required groups,
refresh statistics, and benchmark the actual network/workload before rollout.
