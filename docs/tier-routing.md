# Explicit tier routing and calculation sessions

`routed` is an opt-in coefficient mode. It preserves the existing per-body,
per-epoch **base → medium → extended** selection, but each declared tier can
use local LEB artifacts or an immutable PostgreSQL dataset. `precision` is the
maximum eligible tier, not a date threshold or a request-wide tier pin.

The original `auto`, `leb`, `db`, `skyfield` and `horizons` modes retain their
selection rules. `auto` never consults routing or PostgreSQL implicitly.
Astronomical reductions, public calculation signatures and native coefficient
values are unchanged. See [the DB contract](db-backend.md) for provisioning.

## Configuration

A typical mixed deployment declares every local file explicitly:

```toml
[libephemeris]
mode = "routed"
precision = "extended"
network_policy = "auto"

[libephemeris.tier_routes.base]
backend = "leb"
files = [
  "/app/ephemeris/base_core.leb2",
  "/app/ephemeris/base_asteroids.leb2",
  "/app/ephemeris/base_exotics.leb2",
  "/app/ephemeris/base_apogee.leb2",
]

[libephemeris.tier_routes.medium]
backend = "db"
dataset_id = "11111111-1111-4111-8111-111111111111"

[libephemeris.tier_routes.extended]
backend = "db"
dataset_id = "22222222-2222-4222-8222-222222222222"
```

Replace the illustrative UUIDs with published versions. Inject the common DSN
through `LIBEPHEMERIS_DB_URL`, not a committed password. Tier UUID overrides
are `LIBEPHEMERIS_DB_DATASET_BASE`, `LIBEPHEMERIS_DB_DATASET_MEDIUM` and
`LIBEPHEMERIS_DB_DATASET_EXTENDED`. Each DB tier requires a distinct UUID whose
published `tier` matches its declaration. Undeclared tiers are disabled.

Programmatic configuration uses immutable data-only records:

```python
import os
import libephemeris as ephe

ephe.set_tier_routes(
    {
        "base": ephe.TierRoute("leb", files=("/data/base_core.leb2",)),
        "medium": ephe.TierRoute(
            "db", dataset_id=os.environ["LIBEPHEMERIS_DB_DATASET_MEDIUM"]
        ),
    },
    db_url=os.environ["LIBEPHEMERIS_DB_URL"],
)
ephe.set_calc_mode("routed")
```

This smaller example deliberately declares only core local coverage; it does
not discover asteroid, exotic or apogee companions. Production users exposing
those products must provision the corresponding groups.

Precedence is setter → environment → TOML. `set_tier_routes(None)` restores
configuration lookup without changing mode. Validation copies records before
replacing the previous policy and does no file/DB I/O. Changing configuration,
precision or closing resources is an idle-worker action, not hot reload.

The automatic network policy seals HTTP and downloads while permitting the
explicit PostgreSQL transport. Explicit `sealed` also blocks PostgreSQL.
Missing files are never repaired during calculations; build/provision them
separately using the reviewed artifact manifest.

## Lazy selection and failures

The portable `select_tier()` function consumes explicit closed intervals and
returns a selected tier, a missing metadata requirement, or no coverage. The
Python adapter fetches only the next required source. A fully local calculation
has no DB connection, driver invocation, metadata discovery or coefficient
query, including its supporting bodies and auxiliary epochs.

Light-time, speed, searches and series evaluate their actual dates separately.
They may cross tiers within one calculation/session. There is no smoothing,
extrapolation, hardcoded calendar cutoff or forced single tier per request.

Only verified absence or insufficient coverage advances to another tier. A
configured source that cannot be read, has the wrong tier/schema/publication,
or has corrupt declared segments is fatal. It does not trigger file, JPL,
Horizons, lower-precision or partial-result fallback. The already-approved
local analytical and temporal models retain their existing precedence.

`DBError` represents transport availability, `DBDataError` invalid DB data,
and `RoutingDataError` failed declared file data. These are not ordinary
`ValueError`/`KeyError` fallback signals. An absent optional channel keeps the
existing native missing-channel behavior; a damaged declared channel does not.

## Synchronous ownership

Independent calculations manage their own inputs automatically. To share
bounded inputs across a logical chart or search, use the optional public scope:

```python
with ephe.calculation_session():
    mars = ephe.calc_ut(2451545.0, ephe.MARS)
    jupiter = ephe.calc_ut(2451545.0, ephe.JUPITER)
```

Nested scopes borrow one thread-local owner. All DB readers, metadata,
coefficient buffers and operation-local frame results are cleared on outer
exit, including failures. Repeated independent sessions reread remote inputs.
The aggregate retained remote coefficient budget is **256 segments across all
tiers**, not 256 per dataset. There is no long-lived SQL transaction or checked
out connection for the duration of the public scope.

Mixed proxies are always non-cacheable in persistent scientific LRUs. Concrete
file-only frame readers may use the existing local cache; DB frame reductions
are reused only inside the operation and cleared with its reader. File handles
and connection pools can survive; evaluated DB states cannot. `reset_session()`
preserves configuration/resources but never extends input ownership.

The first fatal source error poisons the owner. Catching it inside a factory
cannot make the outer scope exit successfully; further coefficient calls fail
with the retained category. Ordinary caught input/coverage errors do not poison
it. Cancellation and interruption are not replaced by a previously caught
storage error.

Open and close scopes **inside synchronous calculation threads**, never across
an asynchronous `await`. Integration locks/state-reset conventions remain the
caller's responsibility. For example, a chart factory should enter this scope
inside its existing ephemeris lock, then exit it before resetting global state.

## Diagnostics and source provenance

Tracing reports actual `LEB`, `DB` or `Mixed` reductions. Body coverage records
include serving tier/dataset where known; a DB source never invents a filename.
For an all-tier record, `jd_start`/`jd_end` are the outer envelope and
`intervals` retain the actual closed windows. `BodyCoverage.contains()` rejects
gaps inside that envelope; numerical/vector evaluation raises the public
`EphemerisRangeError` rather than evaluating an unsupported edge. Ordinary
coverage misses still permit only the existing curated local models.
Scientific precision is independent of transport. Publication alone does not
prove reviewed status: the source-artifact manifest must match reviewed hashes.

```python
# Explicit health probe: does not inspect optional medium/extended sources.
status = ephe.get_runtime_inventory("base")
```

The default probe ceiling is configured precision and **may contact DB tiers**.
Diagnostics omit DSNs and driver error text. They return readiness, declared
sources and body coverage, not reusable calculation inputs. `get_leb_reader()`
and `get_leb_inventory()` remain file-only; in routed mode they expose only
declared local routes. `get_current_file_data()` returns the empty file record
in DB/routed mode, never a retained JPL kernel from an earlier mode. Analytic
lunar range policy checks the active mode before any file metadata, so such a
kernel cannot narrow an otherwise independent analytic model.

HTTP services can gate readiness at base while optional DB tiers are down, or
require wider datasets explicitly. A bounded health TTL is not a coefficient
cache. Availability errors should become sanitized 503 responses; data and
configuration faults are server errors, not invalid-date 400 responses. Preserve
these categories through causal wrappers and subprocess messages.

## Validation and explicit interfaces

The routing module introduces storage policy, not a new astronomical kernel.
`TierRoute`/`RouteDecision` are data-only records; the closed-interval selector
has no I/O, callbacks, driver or ambient state. The existing DB binary64 contract
and native arithmetic remain the pure numerical contract. File/DB adapters,
lifetime guards and thread-local slots are integration machinery.

Targeted routing tests cover configuration, source purity, exact boundaries,
additional epochs, mixed sources, nested/thread ownership, aggregate budgets,
swallowed fatal errors and persistent-cache exclusion. See the
[Apple Container routing trial](benchmarks/tier-routing-apple-container.md)
for numerical parity, chart timings and explicit limitations. Routing has
measurable local dispatch overhead; zero DB queries does not mean zero cost.
