# LibEphemeris performance plan

Revision 3. Incorporates external review and the owner's explicit RAM-reduction
requirement.

Goal: reduce the time and memory cost of LibEphemeris 4.x on realistic
workloads **without changing numerical output**. The work is incremental:
every phase must justify itself with the results it actually obtains, not
with a desired final figure.

Status: planning only; no implementation changes. Reported baseline measurements
must be reproduced in phase 0 before optimisation starts.

---

## 1. Ground rules

### 1.1 The precision contract

Every change belongs to exactly one class.

| Class | Meaning | Gate |
| --- | --- | --- |
| **A — bit-identical** | Same IEEE-754 operations on the same operands in the same order, or pure caching, storage and layout changes. Every output word equals the baseline. | Bit-identity differential gate (section 3.2), zero differing words, on every platform of the wheel matrix. |
| **A-candidate** | Believed to be bit-identical but the argument is not trivial. | Treated as not approved until a written proof exists. Tests alone are a search for counterexamples, not a proof. |
| **B — rounding-level** | Model unchanged, but individual results may differ in low-order bits. | Not scheduled. Each item needs owner approval, its own measured comparison, a documented error bound and regenerated golden data. |

Consequences:

- No truncated series, no shortened models, no reduced-precision arithmetic,
  no fast-math, no `target-cpu=native`, no change in FMA usage.
- A cache key must contain every input of the cached computation as exact
  bit patterns.
- Reusing a result is allowed only when the operands that produced it are
  bit-identical to the ones that would be used.

### 1.2 What belongs in the plan

- **Performance priorities** are chosen from profiles. An optimisation with no
  profile evidence is listed as an experiment with no promised gain.
- **Integrity, correctness, error handling and memory bounds** are requirements
  in their own right. They stay in the plan whether or not they show up as a
  cost in a profile.
- **Figures** in this document are either *measured* (section 2) or
  *estimates*. Estimates are hypotheses to be replaced by measurements after
  each step; none of them is a commitment.

### 1.3 Engineering rules

- Minimal code, explicit data, few abstractions.
- One step per PR, each with benchmark results before and after and a
  bit-identity report.
- Comments, tests, commits and docs in English.
- Downstream consumer suites are run unchanged after every merged step.

### 1.4 RAM reduction is an acceptance requirement

- Reducing actual RAM use is a project objective alongside speed, not merely
  keeping caches bounded. A faster workload that uses more RAM is a trade-off,
  not a success on both objectives.
- Measure each mode separately against phase 0 on identical data and workloads.
  Report steady-state and peak process RSS, anonymous memory, resident file-backed
  pages, and total physical memory attributable to 1/2/4/8 worker processes.
  Account for shared pages without counting them once per worker where the OS
  exposes that information; otherwise report the measurement limitation.
- The phase 1 default configuration must demonstrate a reduction in RAM on the
  representative workloads before the memory objective is declared achieved.
  Record peak regressions and modes/workloads without a reduction explicitly;
  gains in `auto` do not establish gains in `leb`. Lower anonymous memory alone
  does not establish lower total resident memory.
- Select the cache default using both latency and memory measurements. Larger
  caches may remain explicit throughput options. No final default or numerical
  memory target is fixed before the baseline and budget sweep are available.

---

## 2. Measured baseline

Published 4.0.0 wheel, Apple M5 Pro, CPython 3.12. The reviewed repository at
commit `1abe74d1d5a7693b69a16392795b0521be9892d2`
differs from tag `python-v4.0.0` only in `README.md` and `docs/README.md`, and
the locally installed extension has the same SHA-256 as the PyPI wheel, so
these reported measurements describe that code baseline. Phase 0 records the
wheel, data-file and build identities and independently reproduces the results.

The comparison column is the reference API's Python binding without its data
files (its built-in analytic route), on the same machine.

### 2.1 Per call, consecutive dates (step 0.37 day), µs

| Case | `leb` | `auto` (default) | `skyfield` | reference API |
| --- | --- | --- | --- | --- |
| Position, no speed | 5 | — | — | 13 |
| Position with `FLG_SPEED` | 33–37 | 105–150 | 400–450 | 10–16 |
| `FLG_SPEED \| FLG_EQUATORIAL` | 18–22 | — | — | 15 |
| `FLG_SPEED \| FLG_J2000` | 11–19 | — | — | — |
| `FLG_SPEED \| FLG_HELCTR` | 41–54 | — | — | — |
| `FLG_TOPOCTR` (no speed) | 24 | — | — | — |
| `FLG_TOPOCTR \| FLG_SPEED` | 290–340 | — | — | 82 |
| True node with speed | 5.5 | 5.5 | 96 | 41 |
| Placidus houses | 22 | — | — | 5 |
| `sidtime` | 19.5 | — | — | 2.5 |
| `get_ayanamsa_ut` | 13 | — | — | 0.6 |
| Same date repeated (memo hit) | 1.35 | 88 | 150 | 4.6 |

### 2.2 Anatomy of one chart at consecutive dates, `leb`, µs

Eleven bodies with speed, then Placidus houses.

| Part | Cost |
| --- | --- |
| First body of a new date | 37 |
| Each of the other ten bodies | 11.8 (118 in total) |
| Houses after the bodies | 25 |
| **Chart** | **180** |

Two facts follow. The frame memo already works: bodies after the first do not
rebuild the frame. Houses get no benefit from the preceding position calls
(25 µs after them, 22.6 µs alone).

### 2.3 Full chart at random dates in 1900–2100, µs

| | `leb` | `auto` | `skyfield` | reference API |
| --- | --- | --- | --- | --- |
| Full chart | 862 | 1808 | 2103 | 182 |
| Bodies only, with speed | 831 | 1764 | 2067 | 175 |
| Bodies only, no speed | 720 | 769 | 367 | 110 |
| Mars only, with speed | 271 | 352 | 409 | 20 |

The gap between 862 µs (random) and 180 µs (consecutive) is the cost of
reaching coefficients for a new date.

### 2.4 Memory and first call

| | `leb` | `auto` | `skyfield` | reference API |
| --- | --- | --- | --- | --- |
| RSS after import | 25 MB | 25 MB | 25 MB | 18 MB |
| RSS after first call | 72 MB | 670 MB | 772 MB | 18 MB |
| RSS after mixed use | 73–250 MB | 770 MB | 880 MB | 18 MB |
| First call | 9 ms | 269 ms | 165 ms | 0.05 ms |

In `auto`, `vmmap` shows 550 MB of live large heap allocations and 190 MB of
freed-but-dirty ones.

### 2.5 `leb` and `auto` are not interchangeable

At the same date Mars has the same position word in both modes, while Jupiter
differs by about 0.01″ and Pluto differs as well: `auto` applies
barycentre-to-centre corrections from the planet-centre kernels. Choosing a
mode is a scientific choice as well as a performance one.

### 2.6 README

The README states ~9 µs per position with speed and ~0.11 ms per chart. The
published wheel measures 33–37 µs and ~180 µs at consecutive dates, and
862 µs per chart at random dates.

---

## 3. Phase 0 — measurement and verification infrastructure

Nothing else starts before this is in place.

### 3.1 Benchmark harness — `tools/bench/`

A dependency-free Python driver and a Rust example for the core, same JSON
output. Median and p95 over at least five runs.

Workloads:

1. `seq_position` — one body, consecutive dates, per flag combination.
2. `random_chart` — section 2.3, on three date ranges.
3. `chart_anatomy` — section 2.2: first body, other bodies, houses.
4. `same_date_chart` — memo path.
5. `time_series` — one body, 100 000 consecutive steps.
6. `houses`, `sidtime`, `ayanamsa`, `fixstar`, `topo`.
7. `cold_start`, in three separately reported conditions:
   - first run ever, any one-time preparation still to do;
   - new process, data prepared, operating-system page cache cold;
   - new process, page cache warm.
8. `memory` — after import, first call, first chart and 100 000 random-date
   calls, reporting steady-state and peak RSS, **anonymous memory** and
   **resident file-backed pages** separately. Include peaks during file opening
   and decoding, and 1/2/4/8 worker processes with shared-page accounting as
   described in section 1.4. The global cache budget is per process, not shared
   between workers.
9. `threads` — `random_chart` on 1, 2, 4, 8 threads.

All workloads in `leb`, `auto` and `skyfield`, on `base` and `medium` tiers,
against both the published wheel and the working tree, on the same platform
and the same data files.

### 3.2 Bit-identity differential gate

- Baseline answers from the published 4.0.0 wheel.
- Deterministic corpus of at least two million calls: every public position
  entry point, houses, sidereal time, ayanamsha, stars, nodes and apsides,
  phenomena, rise and set, crossings, a reduced eclipse grid; every body
  class; every accepted flag combination; dates across the tier span including
  segment and chunk boundaries; all three modes.
- Every returned `f64` compared by bit pattern (signed zeros and NaN payloads
  included); returned flags, warnings and error types compared exactly.
- Run per platform. A baseline is valid only for the platform that produced
  it, because parts of the engine call the platform math library.

### 3.3 Existing gates

`cargo test` for the core, `pytest python/tests`, the golden harness and
validation suite, the wheel workflow, and the unchanged test suites of the
downstream consumers.

### 3.4 Profiling recipe

Document the procedure used here (tight loop, `sample` on macOS or `perf` on
Linux, top-of-stack table) and add a profiling profile that inherits `release`
with line tables only.

### 3.5 Exit criterion

The harness reproduces section 2 within noise, and the gate reports zero
differences when the working tree is compared with itself and with the wheel.

---

## 4. Phase 1 — data access without touching arithmetic

All Class A. This is where the eliminable work is largest and best understood.

### 4.1 Step 1a — configurable decoded-coefficient cache

**Evidence.** At random dates 85 % of the samples are in zstd decompression
(`decompress_body_window_words_limited`, `core/src/providers/leb/reader.rs`).
LEB2 chunks span a nominal ten years per body. The decoded cache holds at most
16 entries (`MAX_CACHED_CHUNKS`, `core/src/providers/leb/index.rs:17`) in
1 MiB (`DEFAULT_CHUNK_CACHE_BYTES`, `index.rs:48`), with a 64 KiB window per
entry (`index.rs:20`). One chart needs more entries than that, so the cache
holds roughly one date. The budget cannot be set from Python.

**Change.**

- One **global** byte budget shared by all images of a process, instead of a
  budget per image. Eviction by bytes, least recently used; no fixed entry
  count.
- Configuration through an environment variable, the TOML file and a Python
  setter, with the value reported by `libephemeris status -vv`.
- The budget may be large enough to hold a whole tier decoded. The decoded
  size of each tier is unknown today and is the first thing to measure.
- Reuse one decompression context per thread; decode straight into the
  retained allocation.
- Make cache hits free of a process-wide lock if the `threads` workload shows
  contention.
- The integrity guarantee is unchanged: compressed bytes are verified exactly
  as today, and decoded words exist only in process memory.
- Initial experiment: **64 MiB, global per process**, not a selected release
  default. Compare the existing default with 16, 32, 64 and 128 MiB using both
  latency and memory workloads, including multiple workers. Select the smallest
  budget that meets the measured latency/memory trade-off and the RAM acceptance
  requirement in section 1.4, rather than choosing only by throughput. Holding
  a whole tier decoded is an explicit option after its size has been measured.
- **The budget is not the total memory.** It bounds retained decoded
  coefficients only. Decompression buffers, coefficients still borrowed by
  calls in progress and cache metadata are accounted and reported separately
  by the `memory` workload and by `status -vv`.
- **The shared budget must not become a global bottleneck.** Accounting uses
  atomic counters; lookups are sharded or read-mostly; decompression happens
  outside any shared lock. The `threads` workload is the check.

**Estimate.** With a budget that holds the dates in use, a random-date chart
should approach the consecutive-date cost (862 µs towards ~180 µs). Warm-up
cost and the budget needed for each tier are to be measured.

**Decision this step informs.** How much of the random-date problem a cache
alone solves, and therefore whether phase 3 is needed for steady state at all
or only for start-up and multi-process sharing.

### 4.2 Step 1b — stop copying whole files to the heap

**Evidence.** SPK kernels are read with `read_to_end`
(`core/src/providers/spk/reader.rs:322`) and LEB files with `std::fs::read`
(`core/src/providers/leb/reader.rs:2095`). In `auto` this is `de440.bsp`
(114 MB), `planet_centers_medium.bsp` (200 MB) and more. A further 190 MB of
freed large allocations suggests a transient second copy during open.

**Change.**

1. **Open only what a request needs.** In `auto`, a kernel is opened when a
   request first routes to it. Bodies served from LEB must not cause `de440`
   to be loaded. Planet-centre kernels are still needed in `auto` for the
   outer planets; that is by design (section 2.5).
2. **Find and remove the transient second copy** on open, using allocation
   stack logging.
3. **Avoid the whole-file heap image for large kernels.** Two mechanisms
   exist in the code and must be compared on measurements:
   - the bounded/positional readers (`spk/bounded.rs`, `spk/positional.rs`,
     `leb/access.rs`), which read records on demand and give no exposure to
     in-place modification beyond what a normal file read has;
   - memory mapping (`LebFile::open_mapped`, `core/src/providers/leb/reader.rs:2120`),
     which is faster and shares pages between processes but depends on the
     immutability contract described in section 6.2.
   **Decision: start with positional reads plus a bounded record cache.**
   Measure memory, start-up and random access. Introduce mapping only if it
   brings a concrete benefit that justifies the immutability work, and then
   only together with the operational guarantees of section 6.2.

   Positional reads remove the truncated-mapping hazard, but they do **not**
   give a coherent view of a file that is modified while it is being read:
   records fetched at different times could come from different contents.
   The data-update rules of section 6.2 (replace by rename, never rewrite in
   place, identity recorded and rechecked) therefore apply to positional
   reads as well, and a retained record cache is bound to the file identity
   it was read from.
4. **Audit the other whole-file reads** (star catalogue, photometry tables,
   IERS data, the classic DE reader) and leave small files alone.
5. **Bound every cache by bytes** and report totals in `status -vv`.

Whatever verification of file content exists today is kept as it is. If
profiling shows that whole-file hashing dominates start-up, changing *when*
verification happens is an integrity decision to be made explicitly, not a
side effect of this step.

**Estimate.** Anonymous memory in `auto` should fall from hundreds of MB to
tens; first call from hundreds of ms to a few ms plus whatever verification
costs. Both to be measured in the three cold-start conditions.

### 4.3 Step 1c — reuse of identical results and removal of copies

Only items that are bit-identical by construction.

1. Constants evaluated on every call with constant arguments are computed
   once (for example `frame_sin_cos_words(VONDRAK_EPS0_RAD)`,
   `core/src/positions/frame_velocity.rs:955`).
2. Within one call, the matrix and rate entry points of the frame both
   evaluate the periodic pole terms at the same epoch
   (`analytic_output_frame_words_with_nutation` and
   `analytic_output_frame_matrix_rate_with_nutation`). First confirm with a
   profile that the terms really are evaluated twice; if so, evaluate once
   and pass them on.
3. Within one call, states of Earth, the Sun and the deflectors read at
   bit-identical epochs by both the position path and the rate path are read
   once.
4. In the topocentric rate stencil (`local_model_velocity`,
   `core/src/positions/serving.rs:460`), nutation requested more than once at
   the same epoch is evaluated once. The set of distinct epochs is unchanged.
5. In `sidtime` and houses, confirm whether the nutation series is entered
   more than once per public call (`core/src/runtime/api/houses.rs:415`,
   `:548`; `public_houses.rs:582`) and, if so, evaluate it once.
6. Hot callers use the fixed-size `eval_raw_inline` variants instead of the
   `Vec`-returning `eval_raw` (`core/src/positions/mod.rs:235`).
7. Internal maps keyed by integers stop using SipHash
   (`core/src/providers/leb/index.rs:87`, `:112`); hashing is about a quarter
   of the no-speed path in the profile.

**Estimate, decomposed against section 2.2.** None of this is established.

| Part | Today | What could reduce it | Status |
| --- | --- | --- | --- |
| First body, 37 µs | Half of it is two-word sine/cosine | Items 1–2 here; phase 4 for the rest | Unquantified until item 2 is confirmed |
| Other ten bodies, 11.8 µs each | Two-word polynomial rates, deflection, hashing, allocation | Items 3, 6, 7 | No breakdown yet; needs its own profile |
| Houses, 25 µs | Nutation series is about 80 % | Item 5, if a duplicate exists | Residual with the current implementation estimated at 15–17 µs; this is not a proven lower bound, and further bit-identical optimisation of the series is simply not quantified today |

A chart at 180 µs becomes about 165 µs if the first body reaches 22 µs and
about 157 µs if houses reach 17 µs. Anything below that requires savings on
the ten other bodies or on fixed overhead that have **not yet been
identified**. The first task of this step is therefore a profile of "body
after the first" to see whether those savings exist.

### 4.4 Exit criterion for phase 1

Re-run every workload, publish the numbers, re-profile, and rewrite the
remainder of this plan from what the new profiles show.

Publish the latency/memory budget sweep and the chosen default with its
justification. Check the RAM acceptance requirement in section 1.4, including
opening peaks and multi-process totals. A cache-only intermediate step may
increase memory, but must report that trade-off; it does not by itself complete
the memory objective. If a mode fails that objective, retain it as outstanding
work rather than reporting overall success from another mode's improvement.

---

## 5. Phase 2 — batch API (Class A, additive)

1. `calc_ut_many(jd, bodies, flags)` and `calc_ut_series(jds, body, flags)`
   with TT variants, returning arrays (NumPy is already a dependency of the
   binding).
2. One resolved session and one GIL release per call.
3. **Per-element contract, to be specified before implementation:**
   - returned flags per element, since the scalar call returns them per body;
   - warnings: collected per element and emitted in element order, or returned
     as data — decide which, and keep parity with the scalar calls;
   - errors: **in the first version one failing element fails the whole
     call**, and the error identifies the index, body and epoch of that
     element. A separate explicit mode with structured per-element errors is
     added only if partial results turn out to be needed;
   - memo behaviour identical to the scalar path.
4. Every element equals the corresponding scalar call bit for bit; the test
   runs the bit-identity corpus through both.

The saving is the per-call overhead that the scalar path pays eleven times.
It must not be counted again on top of the frame reuse that already exists.
No figure is promised before a prototype is measured.

---

## 6. Phase 3 — persistent decoded store (prototype, Class A)

Only if phase 1 shows a remaining need. Its justification is specific:
start-up without decode work and sharing of pages between worker processes.
Steady-state speed is expected to come from step 1a already.

### 6.1 Design outline

A decoded image with the same index sections as LEB2 and the coefficient
payload stored uncompressed and aligned, produced from the verified compressed
file by an explicit CLI command, at the end of `download`, or optionally on
first use. `extended` is opt-in because of its size. The bundled `base` data
stays compressed in the wheel. Evaluation code receives the same coefficient
words as today.

### 6.2 Integrity and immutability — to be decided before any code

The first draft was incomplete here. The plan must state which guarantee is
offered and how it is enforced.

1. **Content.** Proposed: a checksum per block, checked the first time each
   block is used in a process. To be defined as part of the same decision:
   - how the **block-checksum table and the index themselves** are
     authenticated (for example one digest over header, index and checksum
     table, checked at open and bound to the pinned hash of the compressed
     source);
   - what happens on a mismatch (reject the decoded file, fall back to the
     compressed one, rebuild).
   Alternatives to weigh: full verification at every open (strongest,
   slowest), or trusting the local cache with a weaker stated guarantee.
2. **Immutability.** Opening read-only and publishing by atomic rename do not
   prevent later writes to the same inode, as
   `core/src/providers/leb/storage.rs` already notes. A sentence in the
   documentation is not a mechanism. Options to evaluate and combine:
   - decoded files live in a directory owned by the library, are created with
     read-only permissions, and are only ever replaced by rename, never
     rewritten in place;
   - library code paths that update data (download, regeneration) are audited
     to guarantee that rule;
   - identity of the mapped file (device, inode, size, modification time) is
     recorded at verification and rechecked at open;
   - positional reads instead of mapping wherever the library does not control
     the file, which removes the truncation hazard altogether.
   What remains undetected (a block changed in place after it was verified in
   a running process) must be written down as the residual risk.
3. **Fallback.** Read-only installations and sealed environments keep the
   compressed path with the cache of step 1a.

### 6.3 Exit criterion

A prototype measured in the three cold-start conditions and in the
multi-process memory workload, with the integrity design reviewed.

---

## 7. Phase 4 — deeper CPU work, only if the new profiles justify it

| Item | Class | Notes |
| --- | --- | --- |
| Evaluate the 22 independent pole-term recurrences in parallel lanes (`vondrak_pair_words`, `frame_velocity.rs:1000`) | A | Each lane performs the same operations in the same order. Verified word for word against the current scalar code and with the bit-identity gate on every platform. |
| Parallel lanes for the three components in the two-word polynomial evaluators | A | Same conditions. |
| Restructure the nutation series loop for memory layout, keeping summation order and the same sine/cosine calls | A | Gain unknown; measure before keeping. |
| Stop the sine/cosine recurrence at its fixed point | A-candidate | Terms alternate in sign and round-to-nearest is not symmetric at binade boundaries, so "later terms are smaller" does not by itself prove they are no-ops. Requires a written proof for this two-word addition; otherwise not done. |
| Profile-guided optimisation of the wheels | A | Changes layout and inlining only; proven per wheel by the gate. |
| Thread-local or sharded memo, lazy construction of optional helpers | A | Only if the `threads` workload or a profile shows a cost. Experiments, no promised gain. |

---

## 8. Class B items — not scheduled, owner decision required

Each needs its own measured comparison and error bound before it can even be
proposed. "Only the last bits change" and "the error is equal or smaller" are
claims to be demonstrated per item, not assumptions.

| Id | Change | Motivation | First evidence |
| --- | --- | --- | --- |
| B1 | Analytic rates in `auto`. `BoundProvider::WithCenters` reports `has_analytic_kinematics() == false` (`core/src/positions/context_serving.rs:559`), so every body takes the sampled route and re-runs the pipeline at each stencil sample. | `auto` is the default and costs 105–150 µs per position with speed against 33–37 µs in `leb`. | Mars longitude rate differs between the analytic (`leb`) and sampled (`auto`) routes by about 1e-13 °/day at one date. One data point, not a bound. |
| B2 | Analytic rates on the direct SPK route. | 400–450 µs per position with speed. | None yet. |
| B3 | Houses and sidereal time take nutation from the stored LEB channel. | Would remove the series evaluation from houses and sidereal time. | None yet; the two sources are numerically different by construction. |
| B4 | Plain binary64 rates. | — | Not recommended. |

Class A work still helps `auto` and `skyfield`, because every stencil sample
gets cheaper through phase 1.

---

## 9. Documentation

1. Replace the README performance table with harness output, stating workload
   (consecutive or random dates), mode, tier, flags and hardware. Any
   published comparison with the reference API states which route each side
   used: it is a comparison between configurations, not necessarily at equal
   model and accuracy (in section 2 the reference side ran its built-in
   analytic route without data files).
2. Document the modes honestly: `leb` is the fastest route; `auto` applies
   planet-centre corrections and gives different results for the outer
   planets. Neither is presented as a drop-in replacement for the other.
3. Document the cache budget and the memory figures, distinguishing anonymous
   memory from resident file-backed pages.

---

## 10. Order of work

| Phase | Content | Decision at the end |
| --- | --- | --- |
| 0 | Harness, bit-identity gate, profiling recipe | Baseline reproduced |
| 1a | Configurable global decoded cache | Is a cache enough for steady state? |
| 1b | Lazy open, no whole-file heap copies | Is mapping needed, or do positional reads suffice? |
| 1c | Reuse of identical results, copies, hashing | Do savings exist on bodies after the first? |
| — | Re-measure, re-profile, rewrite the rest of the plan | — |
| 2 | Batch API | Per-element contract agreed first |
| 3 | Decoded store prototype | Integrity design agreed first |
| 4 | Lanes, series layout, PGO, early exit | Only items a profile supports |
| — | Class B | Only on approval |

No final target is set. The expected direction, as hypotheses:

- random-date chart in `leb` moving from 862 µs towards the consecutive-date
  cost of ~180 µs through phase 1a;
- a further reduction below ~157 µs only if phase 1c finds savings on the
  bodies after the first;
- anonymous memory in `auto` falling by roughly an order of magnitude through
  phase 1b.

---

## 11. Definition of done for every PR

The implementation checks below apply when the relevant code or harness exists.
This initial plan-only draft is validated by reviewing the document and checking
the diff; it does not claim benchmark, bit-identity or test-suite results.

- [ ] Harness results before and after, same machine and data.
- [ ] Bit-identity gate: zero differing words, all modes, each platform.
- [ ] `cargo test`, `pytest python/tests`, golden harness green.
- [ ] Downstream consumer suites green and untouched.
- [ ] Error behaviour unchanged, or the change is described and tested.
- [ ] Memory workload shows bounded use over 100 000 random-date calls, with
      anonymous and file-backed figures, steady-state/peak RSS and multi-process
      totals reported separately; any memory increase is identified as a trade-off.
- [ ] Phase 1 completion demonstrates RAM reduction against phase 0 as specified
      in section 1.4, with unmet modes/workloads explicitly left outstanding.
- [ ] No change to model constants, series lengths or arithmetic precision.
- [ ] Estimates in this plan replaced by the measured values.

---

## 12. Decisions

Recommended by the external review and adopted here as working decisions.
The owner can overrule any of them.

| Question | Decision |
| --- | --- |
| Cache budget | 64 MiB global per process is an experimental value, not the release default. Compare the current default with 16/32/64/128 MiB; select the smallest budget justified by latency and actual RAM use, including multiple workers. Whole-tier retention is an explicit option after decoded sizes are measured. |
| RAM acceptance | Demonstrate lower RAM on representative workloads against phase 0; report steady-state, peak and multi-process memory per mode. Faster execution with higher RAM is a trade-off, not completion of the memory objective. |
| Large kernels | Positional reads with a bounded cache first. Mapping only if measurements show a benefit worth the immutability work. |
| Persistent decoded store | Deferred to the results of phase 1. Worth pursuing mainly if start-up and total memory across processes remain a problem. |
| Batch errors | First version fails the whole call, naming index, body and epoch. Per-element structured errors only as a later explicit mode. |
| Class B (B1–B3) | Deferred until Class A is measured. B1 is evaluated first afterwards, because it concerns rates in the default mode. |

## 13. Starting scope

Phase 0, then step 1a (global configurable cache) and step 1b (lazy open and
positional reads). Everything else stays in this plan as options to be
selected from the measured results; none of it needs to be decided now.

## 14. Phase 0 implementation checklist

- [x] Isolated published and candidate wheel environments, with matching Python
  and NumPy and recorded extension and wheel SHA-256.
- [x] Pinned input inventory, provider preflight, timing and memory drivers,
  and a native timing example.
- [x] Streaming bit comparison, compressed captures, coverage summary and
  tests that detect changes to bits, result shape, flags, warnings and errors.
- [x] Add the 10,000-call PR check and manual full mode to the wheel workflow.
- [x] Document commands and profiling in [tools/bench/README.md](tools/bench/README.md).
- [ ] Complete numerical qualification across the wheel matrix and record
  discrepancies, incomplete checks and consumer test results.
- [ ] Complete the baseline timing and RAM matrix. Deferring the cache budget
  selection until both can be compared against these measurements.

Local reports and logs are under `target/bench/`. The first five-repeat
candidate run on macOS arm64 (Python 3.12, base/leb, 500 charts per repeat)
measured **857 µs** per random-date chart and **171 µs** per consecutive-date
chart, compared with the reported 862 and 180 µs. On the published wheel, a
single base/leb worker reached approximately **56 MB RSS** and **34 MB physical
footprint** after 100,000 random-date charts. It used an explicitly staged
base core. These numbers qualify that specific setup and workload. They do not
reproduce the original unspecified mixed-use memory workload.

The paired 1900–2100 timing run covered all six mode/tier pairs with five
repeats of 300 charts each. Base/leb random-chart medians were 876 µs for
the published wheel and 874 µs for the rebuilt wheel. Hosted timings are
intentionally advisory. The data-v4 body indexes report 44,544,768 bytes
(42.5 MiB) of decoded coefficients in the base core and 163,322,016 bytes
(155.8 MiB) in the medium core, before metadata, other tables, scratch
space or any per-process duplication.

Validation on macOS arm64:

- The 10,000-call comparison passed. Eight comparator regression tests passed.
- The Rust release suite passed with `bundled-data,nbody,python-internals`
  after provisioning the SHA-256-pinned data-v3 input through
  `LIBEPHEMERIS_TEST_DATA_V3_LEB2`. It includes 3,238 passing library tests
  plus integration and documentation tests. The installed-wheel Python
  suite passed 1,482 tests.
- The external scientific checkout's selected historical calendar lane
  passed 30 tests. The full external qualification campaign remains open.
- Kerykeion offline passed 13,436 tests, kerykeion-parity's source suite
  passed 51, and the Astrologer Studio backend passed 2,139. The Astrologer
  API suite had missing qualification fixtures for its dependency pair and
  three tracing failures in the full suite. The same focused tracing suite
  passed with both wheels, and both wheels lacked the qualification fixtures.
- The extended macOS corpus found a coverage discrepancy on planet-centric
  calls in `auto/medium`. For example, at TT JD 2642370.864629844
  (`0x414428e16eac30d4`), `calc_pctr(jd, 5, 14, 130082)` returned flags
  114274 with a position on the published wheel and raised
  `EphemerisRangeError` in the rebuilt wheel during the run. Both return
  identical results when that call is run alone. Keep this qualification
  open and preserve the sequence captures for reproduction.

The first local full run found 70 inconsistent `auto/medium` planet-centric
results. The staging helper used hard links, so creating and removing links
from a concurrent benchmark changed inode ctime during the comparison.
That could cause the reader to admit a new allocation for identical SPK
bytes. A focused reproduction on the published wheel induced that metadata
change and alternated between a value and `EphemerisRangeError`. The runtime
bug has a separate fix in [PR #88](https://github.com/g-battaglia/libephemeris/pull/88).

The harness now copies each input to an independent inode and validates the
copy. The repeated local full run completed **2,100,000 calls with zero
differences**, using the same published and rebuilt wheels as the first run.
Its captures and comparison JSON are in `target/bench/full-private/`.

The hosted 10,000-call check passed on Linux x86-64, Linux ARM64 and Windows
x86-64. Both Linux jobs passed the [full hosted comparison](https://github.com/g-battaglia/libephemeris/actions/runs/43853954157).
The Windows job was canceled after the 90-minute job limit elapsed. Its logs
show all LEB assets provisioned, followed by the cancellation. The
[diagnostic run](https://github.com/g-battaglia/libephemeris/actions/runs/43867334410)
adds progress during the corpus capture. A sample of the random-date base/leb workload on the
published wheel put Zstandard decompression and LEB coefficient unpacking at the top of
the stacks, supporting the cache experiment once qualification is resolved.

The hosted qualification now reports pinned data provisioning and prints
each capture's progress by mode and tier every 50,000 calls. The diagnostic
change passed all nine comparator tests and the 10,000-call macOS smoke
comparison with zero differences. It leaves the corpus and capture schema
unchanged.

The complementary timing matrix passed on macOS for both isolated wheels.
Each wheel completed all 42 mode/tier combinations of repeated-date charts,
time series, houses, sidereal time, ayanamsha, stars and topocentric calls.
Five fresh processes ran each case, with 100 operations per timing repeat
and 1,000 calls per series. Reports are
`target/bench/family-matrix-{baseline,candidate}.json`. These runs shared
the host with memory qualification, so they establish workload coverage.
Dedicated performance comparisons use the raw timing samples under
controlled conditions.

The base-tier macOS memory comparison used one worker per wheel and 100,000
random-date charts per mode. Published/rebuilt steady RSS was 109.6/108.8 MB
in `auto` and 190.1/190.2 MB in `skyfield`. Peak RSS was 136.1/135.2 MB
in `auto` and 190.1/190.2 MB in `skyfield`. Reports and their raw `vmmap`
snapshots are in `target/bench/base-native-memory-{baseline,candidate}*`.

The [Windows 32 MiB cache run](https://github.com/g-battaglia/libephemeris/actions/runs/43863141515)
finished all 2,100,000 calls and recorded 537 differences, all in `auto`.
525 changed between a result and `EphemerisRangeError` in both directions.
The other 12 were angular output or angular rate changes at Mars dates
outside the LEB coverage. These retained the same returned flags, distance,
and radial rate. The focused `auto-edge` profile runs the first 30,000
requests of each `auto` tier with the full input set. The wheel workflow
also runs the published wheel twice as a control. The focused profile
omits the preceding LEB runs in the full sequence. The unchanged-build
Windows full comparison and diagnosis remain open.

The exact Windows wheel binaries embed Rust standard-library source
revisions `58cf8f9c35b355199727627d07ae34a4e2846d67` (published, 1.99.1)
and `0ef1bc386059596dada8fb65b50425aab2ea637a` (candidate, 1.100.0).
These identify the compiler and standard library difference independently
of the runner's current `rustc -Vv`. The optional wheel workflow toolchain
matches the published build for a controlled rerun.

The separate [SPK storage Windows run](https://github.com/g-battaglia/libephemeris/actions/runs/43867334424)
completed 2,100,000 calls with 562 differences. Comparing the two published
4.0.0 captures from the cache and SPK runs against each other also completed
2,100,000 calls. They used the same extension SHA-256, reported platform
and Python, input hashes and request sequence. The published wheel
disagreed with itself on 546 `auto` requests: 536 value/range-error
changes and ten cases with values on both sides. This cross-run control
establishes that both types of difference can occur before any optimization.
The focused paired controls and the effect of the center-retention fix
remain under investigation. The strict bit comparator still fails on any
observed difference.

The unchanged rebuild's first 1,400,000 requests also exhibited 544
differences, all in `auto`. Counts were 534 typed range-error differences
and ten value differences. That run reached Skyfield medium before the
90-minute job limit and retained partial files. These counts cover only
the completed LEB and Auto prefixes.

A macOS diagnostic reopened the same center supplement after changing
only the ctime of a private copy. Across 24 repetitions, the published
wheel's context calls returned three bit-distinct apparent Mars vectors
at each of the sampled outside-LEB dates, JD 2396756.5, 2506333.5 and
2688954.5. With `FLG_NOGDEFL`, each request produced one vector. The
context's planet-centric call returned both values and range errors.
With the [center retention correction](https://github.com/g-battaglia/libephemeris/pull/88),
each request produced identical bits in every repetition and every
planet-centric call returned a value.

The [Windows fixed-wheel repeat](https://github.com/g-battaglia/libephemeris/actions/runs/43882915978)
compares two independent captures of the wheel with the center correction.
The [optimization check](https://github.com/g-battaglia/libephemeris/actions/runs/43883394649)
applies the same isolated correction to the reference and optimized wheels
and uses Rust 1.99.1 for both builds. The comparison with the released
4.0.0 wheel remains a separate strict run.

On macOS, the composed optimization stack passed all 2,100,000 calls
against the corrected reference with both wheels using indexed LEB and
the candidate using positional centers. Raw captures are in
`target/bench/corrected-reference-indexed-full-mac/`. The companion
[Windows owned-reader run](https://github.com/g-battaglia/libephemeris/actions/runs/43884616052)
uses the same correction and the candidate's 32 MiB process cache.
The corresponding macOS owned-reader comparison also passed 2,100,000
calls with zero differences. Its report records the 32 MiB cache and
positional centers as enabled on the optimized wheel. Raw captures are
in `target/bench/corrected-reference-owned-full-mac/`.

The full external Swiss comparison ran the fast and slow pytest lanes
against both wheel overlays with network policy sealed. Each run evaluated
21,359 cases. Both had 147 skips and the same 62 failing case IDs. The
rebuilt wheel had one additional failure in a Chiron SPK input-stability
check. The existing check passed against both wheels when rerun alone.
The source of the transient inventory change remains unresolved. The
external standalone rounds and verification commands were outside these
runs. The full logs are `target/bench/validation-full-leb{,-baseline}.log`.

On the hosted Linux ARM64 runner, the paired base/LEB comparison completed
100,000 charts per worker at 1/2/4/8 processes. At eight processes, published
and rebuilt summed RSS was 397.3/399.9 MB, and their summed PSS was
254.1/254.0 MB. Each worker retained about 29.6 MB of anonymous RSS.
The archived reports include each process's anonymous, file-backed, peak
and proportional memory counters at all four checkpoints.
Linux x86-64 likewise recorded 432.5/432.7 MB summed RSS and
256.4/256.4 MB summed PSS for the two wheels at eight processes.
[Linux memory run](https://github.com/g-battaglia/libephemeris/actions/runs/43875470097)

A paired macOS base/LEB run compared the Python wrapper with the native
release example using the same dates and pinned inputs. Across five
fresh-process repeats of 5,000 charts, Python/native consecutive medians
were 109.5/107.9 µs and random medians were 539.5/545.2 µs. These
measurements put the consecutive path's observed difference at about
1.6 µs per chart and show no benefit in the random-date comparison.
They do not justify a new batch API on the basis of wrapper overhead.
See `target/bench/binding-paired/` for the samples and per-call p95.

The fresh-process startup report now includes import, first-call and first-
chart medians and p95s. The warm-page condition explicitly primes the staged
inputs. First use has only the staging reads. Both controls passed a
five-repeat base/leb run on macOS, where the cold-page condition reports a
missing prerequisite. On Linux the driver evicts only the staged files
and verifies zero resident bytes with `fincore` before each cold sample.
Three tool cases refuse a nonzero count or incomplete probe output.

The remaining baseline PSS and anonymous/file splits run as separate
Linux cases for [LEB medium](https://github.com/g-battaglia/libephemeris/actions/runs/43886523936),
[Auto base](https://github.com/g-battaglia/libephemeris/actions/runs/43886529057),
[Auto medium](https://github.com/g-battaglia/libephemeris/actions/runs/43886534322),
[Skyfield base](https://github.com/g-battaglia/libephemeris/actions/runs/43886540807)
and [Skyfield medium](https://github.com/g-battaglia/libephemeris/actions/runs/43886546316).
The SPK candidate's Linux counterparts cover
[Auto medium](https://github.com/g-battaglia/libephemeris/actions/runs/43886630970)
and [Skyfield medium](https://github.com/g-battaglia/libephemeris/actions/runs/43886636910).
The positional candidate is also measured in
[Auto medium](https://github.com/g-battaglia/libephemeris/actions/runs/43886693099).
Each case uses 1/2/4/8 processes and 100,000 charts per process.

For the step 1c investigation, the native profiling example separately
attributes the first body and the ten subsequent bodies. A base/leb
50,000-chart consecutive run produced medians of 25.5, 75.3 and 16.0 µs
for the first body, the ten subsequent bodies and houses. This was a
profiling build with other qualification running on the host.
Of 692 first-body samples, 327 were in the two-word sine/cosine function.
Of 2,034 subsequent-body samples, 646 were in coefficient evaluation or
the two-word polynomial rate path and 258 were in deflection arithmetic.
The observations guide the reuse work. The native trace and raw timings
are under `target/bench/chart-anatomy-profile/`.

## 15. Step 1a experiment

An opt-in process pool retains whole decoded chunks across runtime-owned LEB
images. Each shard evicts its least recently used entries by byte capacity.
Explicit indexed readers retain their existing resource policy. The legacy
default stays in place pending the combined latency and RAM measurements.

Mac arm64, 1900–2100 random charts, five repeats of 1,000; RSS after a separate
5,000-chart run. All rows use the same rebuilt wheel and pinned input bytes.

| Process budget | Base median | Base RSS | Medium median | Medium RSS |
| --- | ---: | ---: | ---: | ---: |
| Legacy | 536 µs | 46.7 MB | 552 µs | 84.0 MB |
| 8 MiB | 731 µs | 62.9 MB | 743 µs | 94.3 MB |
| 16 MiB | 234 µs | 66.8 MB | 226 µs | 104.3 MB |
| 32 MiB | 126 µs | 72.5 MB | 132 µs | 109.8 MB |
| 64 MiB | 127 µs | 72.6 MB | 132 µs | 110.2 MB |
| 128 MiB | 129 µs | 72.6 MB | 132 µs | 109.9 MB |

The 32 MiB pool retained 24,720,192 coefficient bytes for this date range.
It is the smallest tested budget that reached the latency plateau for this
workload. Every opt-in budget increased RSS here. The separate positional
reader experiment saved roughly 8 MB in base on a short run, with higher
latency. These results justify proceeding with lower-copy data access and
measuring the two changes together before selecting a default.

- [x] Configurable shared budget, whole chunks, per-shard LRU.
- [x] Per-process retention diagnostics and isolated benchmark sweep.
- [x] Ten local Python checks and the Rust pool lifetime/resize test passed.
- [x] Candidate with 64 MiB matched 10,000 published-wheel calls bit for bit.
- [x] The macOS 32 MiB full comparison passed 2,100,000 calls with zero differences.
- [ ] Finish [cross-platform qualification](https://github.com/g-battaglia/libephemeris/actions/runs/43863141515).
- [x] Measure 100,000-chart process memory and 1/2/4/8 thread throughput.
- [ ] Re-profile cold misses and quantify decoder context reuse and workspace.
- [ ] Select a default after the data access changes and RAM acceptance check.

The 32 MiB base/leb run retained roughly 24.7 MB of coefficient words per
worker after 100,000 charts. Summed RSS was 76.4, 152.4, 305.4 and 626.9 MB
for 1, 2, 4 and 8 independent workers. These sums include shared pages.
The paired thread experiment returned median batch averages of 143, 108, 94
and 130 µs per chart at 1, 2, 4 and 8 threads (legacy: 557, 348, 264 and
387 µs). Throughput peaked at four threads in this experiment.

The installed-wheel suite passed 1,485 tests in each of the legacy and
32 MiB configurations. The standalone CLI regression also verified that
`status -vv` reports the budget before any coefficient acquisition.

### Decoder-workspace experiment

The isolated `codex/leb-decoder-workspace` branch reused one direct zstd
context per thread for owned, bounded, single-frame decoding. Indexed
resource admissions and the streaming route used the existing decoder
paths. All 243 LEB tests passed, including recovery from a corrupt frame
after a successful decode.

Compared with `31f55eee`, the experimental `2f6d6e61` wheel moved random
chart medians from 555.5 to 553.5 µs in base and 558.2 to 555.1 µs in medium.
Consecutive medians moved from 110.5 to 110.7 µs and from 115.8 to 116.5 µs.
These are five independent 2,000-chart repeats per mode and tier on macOS
arm64 with the legacy coefficient cache. Zstandard reported 95,984 retained
bytes per decoder. An eight-thread control retained 671,888 bytes across
seven threads that decoded and returned to zero after those threads joined.
The timing changes are below 1%. The release candidate keeps its existing
decoder lifetime. Reports are under `target/bench/leb-workspace-qualification/`.
