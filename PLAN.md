# LibEphemeris performance plan

Revision 4. The qualification summary at the front records the latest checks
and deployment-shaped measurements. Later sections preserve the experiment log.

Goal: reduce the time and memory cost of LibEphemeris 4.x on realistic
workloads **without changing numerical output**. The work is incremental:
every phase must justify itself with the results it actually obtains, not
with a desired final figure.

Draft [PR #87](https://github.com/g-battaglia/libephemeris/pull/87)
holds the baseline harness, the SPK correction and the Phase 1 changes.
Keep one PR for this effort and review individual changes by commit.
Section 14 tracks the remaining baseline checks and links the measured
results. The Windows Auto reference and its released-wheel limitation are
recorded in the qualification summary.

---

## Current qualification

The draft combines the measurement harness and the measured Phase 1 changes.
The release comparison is direct against the published 4.0.0 wheel on
macOS and Linux. The Windows full Auto gate uses `d86cea0f`, which adds the
SPK role-retention correction to 4.0.0. A Windows control using the released
build configuration disagreed across independent processes in Auto, so
bit identity to the published Windows wheel remains unqualified for Auto.
The Windows extended indexed LEB gate compares directly to the published wheel.

The [extended indexed wheel matrix](https://github.com/g-battaglia/libephemeris/actions/runs/44271681814)
matched 350,000 scalar calls per platform on Linux ARM64, Linux x86-64
and Windows x86-64. The profile uses the cumulative base, medium and
extended inventory and defaults to IERS disabled on hosted runners.
Local runs use the application image's observed IERS snapshot.

| Mode and case | Measured memory result | Timing context |
| --- | --- | --- |
| LEB base, Linux ARM64, eight workers | Summed PSS 259 → 186 MB | Default bundled reader |
| LEB medium, macOS, eight workers | Physical footprint 526 → 198 MB | Random chart median 597 → 626 µs, warm first call 3.66 → 16.19 ms |
| LEB medium, Linux ARM64, eight workers | Default reader flat. Opt-in compressed sharing PSS 443 → 183 MB and cgroup anon + file 516 → 255 MB | Opt-in warm first call 13.6 → 30.9 ms on the ordinary owned reader |
| Auto medium, Linux x86-64, eight workers | Default storage PSS 4.76 → 3.16 GB. Experimental positional comparison 4.76 → 1.55 GB | Separate opt-in positional center policy |
| Skyfield medium, Linux x86-64, eight workers | Summed PSS 5.29 → 2.73 GB | Paired macOS chart timing separately showed roughly 40% lower median |
| LEB extended/indexed, application, one web worker | Container usage 279.8 → 277.9 MB at idle and 306.6 → 307.8 MB at nominal 10 requests/s | CPU 39.9 → 39.6 ms/success, HTTP p95 169.6 → 167.2 ms at that load |

The application comparison uses the same pinned image contents and data,
apart from the wheel, on a local native Linux ARM64 container limited to
one vCPU and approximately 2 GB. Both wheels run with sealed networking,
the indexed reader and prewarm. Each variant ran three times with the same
seeded open-loop arrivals and shifted synthetic dates. At four concurrent
requests the closed-loop medians were 26.50 and 27.25 successful requests/s,
with p95 of 285 and 272 ms. The 930 open-loop responses per wheel had
zero errors. Across paired responses, 882 digests matched and the remaining
48 came from the endpoint which evaluates the current instant.

The one-worker application results place both wheels within a few MB of each
other, with similar end-to-end latency. With two web workers
the candidate used about 553 MB at idle and 610 MB at nominal 10 requests/s,
with comparable or lower goodput. The shared Linux medium reader and the
owned-cache tuning remain specific to their documented configurations.

The reference API was also exercised through the same application and its
supported backend switch, with two runs on compressed ephemeris files and
one analytical planetary fallback run with ancillary data. On the file
route, usage was 272.3 MB at idle and 293.2 MB at nominal 10 requests/s.
At that offered rate the candidate recorded 39.6 CPU-ms/success and 167 ms
HTTP p95, and the reference recorded 39.0 CPU-ms/success and 205 ms p95.
At four concurrent clients the reference served 29.75 successful requests/s
at 230 ms p95 and the candidate served 27.25 at 272 ms p95. Both routes
responded successfully across the measured mix. The backend switch compares
application configurations using different ephemeris models.

The separate 1900–2100 random-chart microbenchmark used the cumulative
inventory and the image's IERS snapshot. The application wheel pair moved
from 830 to 526 µs at the median and from 943 to 787 µs at p95, with
15,000 charts per wheel. The reference API measured 37 µs with compressed
files and 101 µs in its fileless analytical mode on the same dates and
call sequence. The latter figures compare different ephemeris models.
After 100,000 charts per worker the released and candidate wheels used
34.2 and 41.8 MB summed PSS at one worker, and 132.0 and 169.5 MB at five.
The higher retained-coefficient capacity is a speed and memory trade-off.

In a diagnostic container, the measured native entry points accounted
for 1.67 of 15.54 container CPU-seconds over 370 successful requests.
The timing wraps the public native functions in the web process. Native
class methods and subprocesses are outside the attribution. `rise_trans`
used 0.68 CPU-seconds and `calc_ut` used 0.47 in that window.

The application cache screen tested 2, 6, 8, 16 and 32 MiB in the indexed
reader. At nominal 5 requests/s, the 2 MiB setting reduced mean container
usage from 303.2 to 299.9 MB, with p95 moving from 132.0 to 138.1 ms.
The other settings were within approximately 1 MB of the default. In the
separate chart benchmark, the 2 MiB setting achieved a 659 µs median and
724 µs p95, with steady PSS of 34.3 MB at one worker and 132.3 MB at five.
Keep the 8 MiB default for its chart median and document the smaller
setting for memory-constrained indexed workloads. The decoded
coefficient lengths for the four extended images alone sum to 18.0 GB,
compared with 2.87 GB on disk. A decoded archive, batch API or deeper CPU
work needs separate workload evidence. The RAM target for this application
configuration remains open.

### Local follow-up on indexed decoding

Two separate prototypes used the same twelve pinned images. Each indexed
timing case used five fresh processes per arm, 3,000 charts per process and
alternating A/B order. All measurements in this follow-up ran locally.

- A 16 KiB indexed window improved the macOS 2 MiB cache's modern
  random-chart median from 592 to 493 µs. At the default 8 MiB, the
  modern median moved from 402 to 475 µs and p95 for dates within a
  20-year interval moved from 132 to 266 µs. Keep the 64 KiB window.
- Reconstructing adjacent rows within each shuffled byte lane uses the
  existing output allocation. On macOS ARM64 at 8 MiB, the modern median
  moved from 418 to 382 µs and the p95 from 653 to 561 µs. Peak RSS was
  21.87 and 21.86 MB. On local Linux ARM64 with the image's IERS data,
  the modern median moved from 515 to 498 µs and the broad-coverage
  median from 769 to 695 µs. The consecutive-date median was effectively
  flat on both systems.
- The macOS owned-reader comparison used five repeats of 2,000 charts.
  It measured random-chart medians of
  598 → 505 µs for base and 588 → 504 µs for medium. Consecutive
  medians stayed within 1% of the earlier build.
- With the same two locally built Linux wheels in the API image at
  nominal 10 requests/s, idle usage was 277.9 MB in both arms and
  mean load usage was 306.4 → 306.6 MB. CPU was 41.0 → 41.2 ms
  per success. Each arm completed 650 requests. The p95 ranged
  from 160 to 183 ms before the change and from 164 to 167 ms
  after it. The isolated engine measurements justify retaining the rewrite.
- In a separate application run, limiting glibc to two allocation
  arenas kept idle usage near 278 MB and lowered mean usage under load
  by about 2.4 MB. Latency varied across the two runs per arm, so the
  deployment keeps its ordinary allocator configuration.

The cumulative extended indexed comparison and the ordinary LEB smoke
case matched 350,000 and 10,000 calls on both local operating systems.
Raw reports and local reproduction scripts are under
`target/bench/local-indexed-screen/`.

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
- Keep implementation steps reviewable by commit in PR #87. Record before
  and after measurements and bit-identity evidence for each step.
- Comments, tests, commits and docs in English.
- Run the downstream consumer suites unchanged on the composed candidate.

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

**Observed.** For random dates in 1900–2100 on macOS, 32 MiB reduced
the paired chart median from 589.7 to 138.2 µs in base and from 591.9
to 149.4 µs in medium. Larger budgets reached the same timing plateau
in the sweep. Every tested opt-in budget raised memory in the paired
macOS run, and the medium full-coverage workload lost throughput at
32 MiB. The default remains the per-image policy. Section 15 records
the budget sweep and the warm-up results.

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

**Observed.** In Linux x86-64 Auto medium at eight workers, summed PSS
moved from 4,760 to 3,158 MB with default SPK storage and to 1,553 MB in
the separate positional-center pair. In the default storage pair the
median first call was 470 → 346 ms on first use, 458 → 350 ms with
warm pages and 560 → 502 ms with verified cold pages. Section 14 records
the other startup cases and their input conditions.

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

**Observed.** In the native macOS base/LEB chart profile, medians were
23.38 µs for the first body, 69.58 µs for the ten subsequent bodies
combined and 14.75 µs for houses. In a separate paired macOS benchmark,
reusing reference-obliquity constants changed the first-body median
from 24.33 to 23.88 µs and the full chart median from 110.92 to 110.42 µs.
The median for the ten subsequent bodies stayed at 71.71 µs in that pair.
The measured paths and remaining CPU candidates are recorded in the
individual profiles and section 14.

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
| 2 | Batch API | Revisit if a measured workload shows significant binding overhead or a separate API need |
| 3 | Decoded store prototype | Integrity design agreed first |
| 4 | Lanes, series layout, PGO, early exit | Only items a profile supports |
| — | Class B | Only on approval |

The phase 1 comparisons give these results:

- In the macOS LEB 1900–2100 workload, a 32 MiB pool changed random-chart
  medians from 589.7 to 138.2 µs in base and 591.9 to 149.4 µs in medium
  against the corrected reference. Broader medium coverage and measured
  multi-process RAM favor retaining the existing default.
- In the native LEB base chart profile, medians were 23.38 µs for the
  first body, 69.58 µs for the other ten and 14.75 µs for houses.
  Decompression still dominates indexed random charts. The scalar
  Python/native comparison does not currently justify a batch API for speed.
- In the paired macOS runs, the SPK lookup change reduced Skyfield chart
  latency by about 40%. The [default wheel matrix](https://github.com/g-battaglia/libephemeris/actions/runs/43949734696)
  qualified the improved lookup against the corrected reference.

The following comparisons used eight processes and 100,000 charts per
process. Each row is a separate paired run. Values are decimal MB.

| Mode and path | Host | Steady memory | Scope |
| --- | --- | ---: | --- |
| LEB base, bundled snapshot | macOS ARM64 | Footprint 276.0 → 206.8 | Bundled base comparison |
| LEB base, bundled snapshot | Linux x86-64 | PSS 261.0 → 187.7 | Bundled base comparison |
| LEB medium, default | macOS ARM64 | Footprint 526.4 → 198.3 | Full coverage against PyPI 4.0.0 |
| LEB medium, default | Linux x86-64 | PSS 477.6 → 477.5 | 1900–2100 against the corrected reference |
| LEB medium, shared opt-in | Linux ARM64 | PSS 443.4 → 182.9 | Full coverage, same wheel with sharing off and on |
| Auto medium, default SPK storage | Linux x86-64 | PSS 4,760.3 → 3,157.7 | Published wheel against the SPK storage candidate |
| Auto medium, positional-center opt-in | Linux x86-64 | PSS 4,760.2 → 1,553.1 | Separate positional-center comparison |
| Skyfield medium, SPK storage | Linux x86-64 | PSS 5,290.8 → 2,728.5 | Separate Skyfield comparison |

Section 14 records paths, inputs and verification for each pair. In the
matched macOS medium comparison the default snapshot changed median chart
time from 597 to 626 µs and the warm first call from 3.7 to 16.2 ms. In
the separate eight-worker Linux ARM64 opt-in run, the per-worker batch
averages were about 1.2 ms per chart and the median first call changed
from 13.6 to 30.9 ms. The Linux measurement at one worker showed no
material RAM saving.

---

## 11. Definition of done

Apply these checks to the relevant changes in the consolidated draft.
Section 14 and the step reports record completed runs and open requirements.

- [x] Harness results before and after, same machine and data. Rechecked
      60 timing pairs and the 24 macOS baseline memory pairs for matched
      data, dates, flags, settings, and workload size.
- [ ] Bit-identity gate: zero differing words, all modes, each platform.
- [ ] `cargo test`, `pytest python/tests`, golden harness green.
- [ ] Downstream consumer suites green and untouched.
- [ ] Error behaviour unchanged, or the change is described and tested.
- [ ] Memory workload shows bounded use over 100 000 random-date calls, with
      anonymous and file-backed figures, steady-state/peak RSS and multi-process
      totals reported separately; any memory increase is identified as a trade-off.
- [ ] Phase 1 completion demonstrates RAM reduction against phase 0 as specified
      in section 1.4, with unmet modes/workloads explicitly left outstanding.
- [x] No change to model constants, series lengths or arithmetic precision.
      The source diff leaves the data and series kernels intact. LEB cache
      entries and SPK records retain the decoded binary64 words. The
      reference-obliquity memo calls the same two-word routine with the
      same constant. The Auto center-retention correction affects source
      selection and is qualified separately.
- [x] Replace Phase 1 estimates with the paired measurements in sections 4
      and 10. Individual reports retain workloads, inputs and raw samples.

---

## 12. Decisions

Recommended by the external review and adopted here as working decisions.
The owner can overrule any of them.

| Question | Decision |
| --- | --- |
| Cache budget | Keep legacy retention as the default. All four process-pool budgets raised measured RAM at every macOS process count. Offer 32 MiB explicitly for the measured narrow-window throughput gain. Medium full-coverage latency and the increase in RAM remain trade-offs. |
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

### Reproducible environments and measurements

- [x] Install the published wheel and the rebuilt wheel in separate
  environments. Record the revision, wheel and extension digests, platform,
  Python, compiler and build options.
- [x] Stage pinned inputs into private directories. Record their coverage,
  use sealed calls with explicit settings, and report unavailable inputs.
- [x] Share workload requests and JSON formats across a Python driver and
  a native Rust example. Document timing, memory and native profiling in
  [tools/bench/README.md](tools/bench/README.md).
- [x] Collect the macOS timing pairs across all six mode and tier cases,
  including chart components, 100,000-position series, random date ranges
  and 1/2/4/8 threads. Record five independent samples and the p95 unit.
- [x] Archive the 1/2/4/8-process macOS RAM pairs across all six cases.
  Retain per-process RSS, peak RSS, footprint and raw `vmmap` snapshots.
- [x] Archive the Linux LEB base RAM pair with its released companion
  files. Keep the available RSS, anonymous, file, peak and proportional
  counters separate.
- [x] Archive the Linux startup matrix for all six mode and tier cases.
  Compare first use, explicitly warm pages and verified cold pages.

### Numerical and consumer validation

- [x] Generate a deterministic 2,100,000-call corpus and capture it in
  separate processes. Archive compressed captures, request coverage and
  raw bitwise results. Probe interval boundaries, state changes, warnings
  and typed errors. Keep the named satellite comparison as a separate run.
- [x] Verify that the comparator detects modified binary64 words, signed
  zeros, NaN payloads, flags, structure, warnings, errors and truncation.
  The smoke and full request generators now also require a real warning
  from the sidereal-mode transition and report the count on each wheel.
  The local 10,000-call release/candidate smoke matched, including one
  warning-emitting call on each side.
  An additional 25,518-call prefix across all six mode/tier cases against
  the published wheel matched exactly, with six real warning calls per
  wheel. Captures are in `target/bench/spk-registration-lookup/warning-prefix/`.
- [x] Run a 10,000-call PR check and support full manual runs on both
  Linux architectures and Windows. Keep the macOS verification local.
- [ ] Pass the full comparison against the published reference on all
  required platforms. Linux and macOS passed. The Windows published wheel
  disagrees with independent runs of itself in Auto; the corrected
  reference and optimized checks are reported separately.
- [x] Run Rust release tests with `bundled-data,nbody,python-internals`,
  installed-wheel Python tests and harness tests. Preserve their logs.
- [x] Execute the external validation and downstream comparisons on the
  unchanged suites. Record missing prerequisites and compare failure
  inventories against the published run.

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
  passed 30 tests. The complete campaign results are recorded later in this section.
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
bug was originally investigated in [PR #88](https://github.com/g-battaglia/libephemeris/pull/88).
The consolidated draft includes its correction alongside the optimization
stack. The first consolidated runtime matched the qualified build at
`b565e146`. Later code changes have their own measurement and
qualification records below.

The harness now copies each input to an independent inode and validates the
copy. The repeated local full run completed **2,100,000 calls with zero
differences**, using the same published and rebuilt wheels as the first run.
Its captures and comparison JSON are in `target/bench/full-private/`.

The saved full comparison was replayed to inventory actual position
requests. It identifies 36 mode/tier/entry-point cases, ten target or
diagnostic categories and 16,406 distinct requested flag masks.
The planet-centric cases use Earth as their center (ID 14). The inventory
preserves exact masks and includes requests that produced errors.
The replay still reports zero differences on 2,100,000 requests.
See `target/bench/full-private/comparison-with-inventory.json`.

ID 9599 is in the planetary-moon namespace, but has no named satellite
in the library's public catalog. A supplemental macOS Skyfield medium
capture registered Charon's alias and three named Pluto satellite IDs
from the local `plu017.bsp` fixture. Its manifest records the SHA-256 and
descriptor coverage. Module positions and the direct moon API used five
epochs, including segment edges and immediate neighbors, with four flag
sets. All 320 requests matched bit for bit on the unchanged and composed
wheels against the published wheel. The published capture contains 88
values. The initial command and raw files are under `target/bench/moon-bits/`.
The same comparison is now available as `bitdiff.py run --profile satellites`.
The isolated local runs reproduced the result for both candidate wheels.
The verified NAIF fixture is provisioned by the
[supplemental wheel check](https://github.com/g-battaglia/libephemeris/actions/runs/43897311620).
All 320 calls matched on Windows x86-64 and both Linux architectures,
with 88 values per capture on each job.

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
The independent archive replay verifies identical executable and wheel
digests, and compares every request after removing only the ephemeral
input staging paths from the comparison headers. Its results for the
published control, strict rebuild and corrected control are under
`target/bench/windows-reference-audit.json`. The strict comparator still
fails on any observed difference.

Replaying the SPK Windows comparison with the complete difference counts
confirmed 556 changes between values and `EphemerisRangeError`, all on
planet-centric calls in Auto. Six Auto context position results differed
with values on both sides. No warnings changed. A replay of the Windows
owned-reader comparison against the corrected reference again found zero
differences. Each replay is saved as `replay-with-difference-counts.json`
in its corresponding local artifact directory.

The [focused Windows center-retention run](https://github.com/g-battaglia/libephemeris/actions/runs/43878253315)
completed both 60,000-request comparisons. The same published binary
produced 594 differences in its self-comparison, including 592 changes
between values and range errors in both directions and two Mars vectors.
Between the published wheel and the wheel with the center fix, 1,537
requests changed from a range error to a value. Of these, 1,509 were
planet-centric. The other 12 differences were Mars positions and rates
at the three outside-LEB epochs below. The returned flags, warnings,
distances and radial rates agreed. This compares the fix with the
unstable release. The fixed wheel's own repeat is tracked separately.
The raw captures are in `target/bench/windows-fixed-edge-43878253315/`.

The [Rust 1.99.1 repeat](https://github.com/g-battaglia/libephemeris/actions/runs/43879212102)
found 562 changes in the published control, all between values and
range errors. The fixed build produced values for 1,543 old range errors
plus the same 12 angular output differences. The two corrected
captures, built with Rust 1.99.1 and 1.100.0, matched bit for bit on all
60,000 requests. Request sequences and pinned input hashes were checked
before replay. These runs had different CPU fingerprints. The complete,
same-runner corrected comparison remains separate. The second capture
is archived locally in `target/bench/windows-fixed-edge-43879212102/`.

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
passed the full 2,100,000-call comparison of two independent captures of
the wheel with the center correction.
The [optimization check](https://github.com/g-battaglia/libephemeris/actions/runs/43883394649)
applies the same isolated correction to the reference and optimized wheels
and uses Rust 1.99.1 for both builds. The indexed and positional
configuration of the optimized wheel matched all 2,100,000 requests.
The same corrected-reference corpus passed on Linux ARM64 and x86-64
for [indexed readers](https://github.com/g-battaglia/libephemeris/actions/runs/43909331459)
and [the 32 MiB owned reader](https://github.com/g-battaglia/libephemeris/actions/runs/43909331865).
Each configuration matched all 2,100,000 calls on each architecture.
A separate
[strict run](https://github.com/g-battaglia/libephemeris/actions/runs/43882398843)
with the unmodified runtime and the published Windows Rust toolchain
still found 572 differences over 2,100,000 calls on Windows. The
archived captures replayed to 566 value/range-error transitions and six
value differences, all in Auto, with zero warning changes. Both Linux jobs
passed. The strict Windows gate against the released wheel remains open.
The wheel workflow on PR #87 supports reproducing the separate corrected
comparison via `bitdiff_profile=full` and `corrected_spk_reference=true`.
The default selection still compares against the published wheel.
The [default-configuration corrected comparison](https://github.com/g-battaglia/libephemeris/actions/runs/43918145944)
at `18ce19be` used legacy LEB retention and the default center reader.
The archived ARM64 and x86-64 Linux captures each matched all 2,100,000
requests, with 350,000 per mode and tier. The Windows job also passed.
The later SPK lookup change has a separate run in the same PR.

On macOS, the composed optimization stack passed all 2,100,000 calls
against the corrected reference with both wheels using indexed LEB and
the candidate using positional centers. Raw captures are in
`target/bench/corrected-reference-indexed-full-mac/`. The companion
[Windows owned-reader run](https://github.com/g-battaglia/libephemeris/actions/runs/43884616052)
used the same correction and the candidate's 32 MiB process cache.
It passed all 2,100,000 calls with zero differences.
The corresponding macOS owned-reader comparison also passed 2,100,000
calls with zero differences. Its report records the 32 MiB cache and
positional centers as enabled on the optimized wheel. Raw captures are
in `target/bench/corrected-reference-owned-full-mac/`.

The composed stack also passed the Rust release suite with
`bundled-data,nbody,python-internals` and the pinned v3 test fixture.
It reports 3,246 passing library tests, 18 ignored library tests, and
passing integration and documentation tests. The installed-wheel Python
suite passed all 1,486 tests with the bundled test core pinned as in CI.
The first Python run selected the base core in the user's larger data
directory and failed the inventory test on the uncovered-date review
status. The unchanged and optimized wheels both reproduce that focused
failure with the larger directory selection. Both pass when the explicit
test core is pinned outside that directory. The raw qualification logs
are `target/bench/optimized-corrected-{rust,python-pinned}-suite.log`.
Using a wheel-only overlay verified by the candidate extension hash,
Kerykeion's unchanged offline suite passed 13,436 tests with 1,392 skips
and three deselections. Kerykeion-parity's unchanged source suite passed
51 tests. Logs are `target/bench/optimized-corrected-kerykeion*.log`.

The unchanged native API consumer suite ran in two independent copies of
its environment with the published and composed wheels installed. The
extension SHA-256 was checked before each run. Both wheels had identical
outcomes on all 106 cases, including the same one failure comparing direct
wheel results with the engine embedded in the consumer. A direct replay
captured identical binary64 words and flags from the two standalone
wheels and the same difference from the embedded engine. Local artifacts
are `target/bench/astrologer-v7-*-installed.{log,xml}` and
`target/bench/astrologer-v7-*-direct-probe.json`.

The paired upstream API v6 run also produced identical outcomes by
test ID. Each wheel had 2,439 passes, four skips, three tracing failures
and 95 setup errors due to missing qualified references for the existing
dependency pair. The initial launcher used stdin as Python's main module,
which made spawned workers unable to import their entry point. Both
qualified runs use the ordinary `python -m pytest` entry point after a
separate digest check. The unchanged suite and its reference fixtures were
left intact. Reports are `target/bench/astrologer-v6-*-qualified.{log,xml}`.

The Windows positional and indexed comparisons started before the
timeout increase. Their Linux jobs passed. The Windows jobs reached the
90-minute limit, including the
[positional run](https://github.com/g-battaglia/libephemeris/actions/runs/43870941467)
where the candidate finished and the published capture was still in the
last configuration. The corrected-reference Windows runs retain the
150-minute job allowance.

The full external reference API comparison ran the fast and slow pytest lanes
against both wheel overlays with network policy sealed. Each run evaluated
21,359 cases. Both had 147 skips and the same 62 failing case IDs. The
rebuilt wheel had one additional failure in a Chiron SPK input-stability
check. The existing check passed against both wheels when rerun alone.
The source of the transient inventory change remains unresolved. The
external standalone rounds and verification commands were outside these
runs. The full logs are `target/bench/validation-full-leb{,-baseline}.log`.

The external `leb all` campaigns finished against the unchanged validation
checkout at `5a2051f2`. The published wheel and composed candidate have
identical outcomes on all 21,359 pytest IDs, including 62 failures and
147 skips in each run. All 273 standalone rounds and all 19 verification
scripts have matching exit statuses. The 269 round wrapper records and
19 verification records match in assertions, failed predicates, call
counts, qualifiers and reported errors. Four rounds reached the same
timeout on both wheels. The complete parity JSON and execution evidence
match, including 3,330 assertions and 253 failures in each run. The result
records the existing scientific failures with their original thresholds.
Logs and the paired audit are under `target/bench/validation-complete-*.log`
and `target/bench/validation-campaign-comparison.json`.

On the hosted Linux ARM64 runner, the paired base/LEB comparison completed
100,000 charts per worker at 1/2/4/8 processes. At eight processes, published
and rebuilt summed RSS was 397.3/399.9 MB, and their summed PSS was
254.1/254.0 MB. Each worker retained about 29.6 MB of anonymous RSS.
The archived reports include each process's anonymous, file-backed, peak
and proportional memory counters at all four checkpoints.
Linux x86-64 likewise recorded 432.5/432.7 MB summed RSS and
256.4/256.4 MB summed PSS for the two wheels at eight processes.
[Linux memory run](https://github.com/g-battaglia/libephemeris/actions/runs/43875470097)

The [Linux medium/LEB baseline case](https://github.com/g-battaglia/libephemeris/actions/runs/43886523936)
completed on ARM64 and x86-64 with its tier-specific release files.
At eight processes, published/rebuilt summed PSS was 470.8/470.4 MB on
ARM64 and 472.9/473.1 MB on x86-64. Summed RSS was 614.0/615.5 and
649.4/649.2 MB respectively. The wheel artifacts contain all four
checkpoints at 1/2/4/8 processes.

The [Linux base/Auto baseline case](https://github.com/g-battaglia/libephemeris/actions/runs/43886529057)
also completed with both wheels and matching inputs. At eight processes,
published/rebuilt summed PSS was 1002.2/1002.2 MB on ARM64 and
1004.5/1004.7 MB on x86-64. Summed RSS was 1144.2/1146.4 and
1179.9/1180.3 MB respectively.

The macOS paired baseline has completed all LEB checkpoints at 1/2/4/8
processes, with 100,000 random charts per process and the companion files
for the selected tier. At eight processes, published/rebuilt summed RSS
was 450.2/454.3 MB for base and 668.1/665.5 MB for medium. Sampled
per-process physical footprints summed to 275.6/274.9 and 514.5/521.7 MB
respectively. Each pair used identical input hashes, epochs and flags.
These per-tier installations are distinct from the mixed-installation
tests used to select an indexed cache budget. Raw Mac snapshots and
reports are in `target/bench/phase0-macos-memory/`.

The paired macOS Auto base run also completed 1/2/4/8 processes with
100,000 charts per process. At eight, published/rebuilt summed RSS was
1199.4/1199.5 MB and sampled physical footprint was 1054.3/1057.2 MB.
The four-process baseline showed 377.8 MB RSS against 603.1 MB rebuilt,
with similar footprints of 526.6/523.6 MB. The raw `vmmap` summaries
record about 191 MiB of writable pages as swapped out for the published
capture and zero in the rebuilt capture. Both reports retain the original
measurements and this case is flagged for residency pressure when
interpreting its RSS.

For Auto medium, both macOS wheels completed the same companion-inclusive
case at 1/2/4/8 processes. At eight, published/rebuilt summed RSS was
4957.6/4959.5 MB and sampled physical footprint was 4814.3/4819.6 MB.
This pair uses the same input hashes. Its figures are separate from the
earlier storage screens with different staging layouts.

The two Skyfield pairs complete the macOS RAM matrix. At eight processes,
published/rebuilt steady summed RSS was 1,538.4/1,530.0 MB in base and
8,012.9/8,020.1 MB in medium. Summed sampled per-process physical
footprint was 1,416.5/1,408.2 and 7,897.8/7,901.9 MB respectively.
All six pairs used matching input manifests and 100,000 charts per
process. Each kept the four phase checkpoints at 1/2/4/8 workers.
Group PSS is unavailable on macOS.

The ARM64 portion of the [Linux base/Skyfield run](https://github.com/g-battaglia/libephemeris/actions/runs/43886540807)
completed with 1249.6/1249.0 MB summed RSS and 1112.2/1111.2 MB summed
PSS at eight processes. The x86-64 job is a separate measurement.

The added macOS timing matrix completed for both wheels with five
independent samples per case. It covers the two wider random-date
windows, 36 cases for 1/2/4/8 threads and the consecutive workloads,
and all six first-use and warm-page startup cases. Each time-series
sample contains 100,000 calls. The pairing was checked for the input
inventory, seed, effective date range, flags, case order and sample
counts. These coverage runs shared the machine with RAM and numerical
qualification. Raw data is under `target/bench/phase0-macos-timing/`.

The README chart table now uses the reproducible base-tier run at
`target/bench/readme-baseline.json`. Its published LEB medians were
118 µs for consecutive dates and 606 µs for random dates. The earlier
timing run measured 171 and 876 µs under different host load. File
digests, workload parameters and route probes agree across both runs.

The separate Mars sweep completed both builds at all six mode and tier
pairs with masks 2, 258, 266, 322, 4354 and 33026. Each case retained
five fresh-process samples of 1,200 consecutive dates. The initial
heliocentric Skyfield probe misclassified the workload by tracing a
whole chart. It now traces only Mars for single-position workloads,
using the requested flags, and keeps the whole-chart trace for chart
workloads. The paired traces match in
`target/bench/phase0-seq-position-flags-qualified/`.

The unchanged external `leb all` campaign is running sequentially against
the published and composed wheel overlays with network policy sealed.
This adds the standalone rounds, verification scripts and parity
matrix to the external pytest cases. Raw console logs are in
`target/bench/validation-complete-{published,optimized}.log`.
The paired Step 1 timing sweep runs after the baseline memory and
external-validation processes. It samples each workload with the corrected
reference and the default optimized wheel, and tests the indexed,
process-cache and positional settings on their applicable routes. Its
commands and reports are under
`target/bench/phase1-composed-timing/`.

The LEB portion completed all 30 report pairs. Both builds use matching
data, routes, epoch sequences and five fresh processes per case. Random
chart medians in µs are shown below, with 300 charts per process.

| Tier and window | Corrected reference | Optimized default | Indexed | 32 MiB owned pool |
| --- | ---: | ---: | ---: | ---: |
| Base, 1900–2100 | 589.7 | 606.0 | 593.4 | 138.2 |
| Base, full coverage | 592.9 | 586.1 | 612.1 | 137.4 |
| Medium, 1900–2100 | 591.9 | 580.6 | 510.6 | 149.4 |
| Medium, full coverage | 603.5 | 605.3 | 674.7 | 925.5 |

The indexed ceiling is 2 MiB in base and 6 MiB in medium. Across the
1900–2100 samples, the 32 MiB pool's p95 for the first 100 charts per
process was 939.2 µs in base and 1,039.6 µs in medium. The last 100
charts measured 147.7 and 233.6 µs. Medium full coverage still had
a p95 of 1,074.6 µs in the last 100 charts. Report the date-window and
warm-up dependence alongside the median speedup.

The same sweep also measured 300 random charts per batch at each thread
count, in five fresh processes per case. The table gives the median of
the batch averages in µs per chart.

| LEB tier and reader | 1 thread | 2 threads | 4 threads | 8 threads |
| --- | ---: | ---: | ---: | ---: |
| Base, optimized default | 561.9 | 353.7 | 265.5 | 389.9 |
| Base, 32 MiB owned pool | 188.1 | 130.7 | 105.8 | 153.8 |
| Base, 2 MiB indexed | 573.4 | 527.4 | 547.9 | 654.4 |
| Medium, optimized default | 600.0 | 363.6 | 311.9 | 406.0 |
| Medium, 32 MiB owned pool | 226.7 | 189.8 | 160.9 | 186.0 |
| Medium, 6 MiB indexed | 522.4 | 495.5 | 495.0 | 473.4 |

The 32 MiB pool improved the batch average at all four thread counts in
this window. The default and 32 MiB owned pool reached their best result
at four threads on this machine. The indexed reader's shared decode gate
serializes cold scratch reservations. The timings alone cannot attribute
the scaling to that gate. A separate worker-sampling run compared the
indexed medium path at one and eight threads and the owned alternatives
at eight. The provenance and traces are under
`target/bench/leb-thread-profile/`.
The four contention samples completed over 64,000 random medium-tier
charts per batch. At eight threads, mutex waits accounted for
1,961 of 15,968 sampled worker frames with legacy retention,
8,445 of 16,384 with the 32 MiB owned pool, and 13,978 of 16,904
with the 6 MiB indexed reader. The one-thread indexed control had none
of these mutex waits. The dominant indexed caller was the coefficient
path containing the cold decode gate. With the 32 MiB pool, the most
frequent caller was chunk directory lookup (4,816 wait samples).
An isolated shared-read directory experiment passed 244 bundled LEB
tests, including single publication under concurrent first access.
Its paired timing used LEB medium, five repeats of 6,000 charts each.
At 32 MiB, median batch time per chart at 1 / 2 / 4 / 8 threads was
155.1 / 113.3 / 97.9 / 169.7 µs with the current lock and
154.0 / 101.9 / 79.0 / 95.8 µs with concurrent cache-hit reads.
The legacy-retention medians moved by 0.4–2.2% across thread counts.
The local 10,000-call gate with the candidate pool active matched,
including the expected sidereal warning.
The shared-read change is in PR #87 for combined qualification.
Raw comparison artifacts are in `target/bench/leb-shared-directory-read/`.

Auto base also completed with matching input and workload settings.
The random-chart median was 1,183 µs on the corrected reference,
1,066 µs on the default candidate and 943 µs with positional centers.
Consecutive chart medians were 863, 766 and 606 µs respectively.
The corresponding Auto medium random-chart medians were 1,229,
1,254 and 1,110 µs. Consecutive medians were 735, 753 and 636 µs.
On full Auto medium coverage the random chart medians were
1,224, 1,222 and 1,063 µs. The default storage changes reduce RAM
in the paired memory measurements. The timing change depends on
the workload and the center reader setting.
The full sweep passed its input and workload checks for all 60 report
pairs. Skyfield's default storage changes primarily improve startup.
The table gives the median in ms, corrected reference → composed build.

| Tier | First call | Chart after first call | First call, warm pages | Chart after first call, warm pages | Random chart | Consecutive chart |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Base | 20.69 → 13.89 | 4.98 → 1.65 | 20.62 → 13.43 | 4.93 → 1.61 | 1.535 → 1.549 | 1.565 → 1.612 |
| Medium | 105.16 → 80.64 | 10.80 → 1.44 | 104.89 → 84.60 | 10.57 → 1.50 | 1.543 → 1.603 | 1.572 → 1.569 |

At full coverage the Skyfield random-chart median was 1.486 → 1.506 ms
in base and 1.549 → 1.546 ms in medium. The 3.9% increase in the
1900–2100 medium random-chart median did not repeat at full coverage.
The first-use condition includes the mandatory staging reads. Only the
warm-page trials explicitly prime the data.
The local comparison audit is `target/bench/phase1-composed-timing/comparison.json`.
All four follow-up native profiles completed, with the same optimized
source and pinned input manifests. Each trace samples three seconds.
The symbol counts below are exclusive leaf samples. The macOS summary
omits leaves observed fewer than five times and the counts are not
independent timing measurements.

| Workload | Calls or charts | Main-thread samples | Observed hot leaves |
| --- | ---: | ---: | --- |
| LEB base, consecutive chart components | 80,000 charts | 2,277 | Two-word sine/cosine (225), polynomial position and rate (201), cached relative velocity source closure (194) |
| LEB medium, indexed random charts with 6 MiB | 35,000 charts | 2,500 | Zstd sequences (719), LEB window decompression (481), Huffman decoding (268), positional read (57) |
| Auto medium, random charts with positional centers | 35,000 charts | 2,487 | Zstd sequences (409), LEB window decompression (349), LEB segment evaluation (171), SPK Chebyshev evaluation (138) |
| Skyfield medium, random charts | 70,000 charts | 2,513 | Allocator free (340) and tiny allocation (210), SPK registration sort (126), SPK Chebyshev evaluation (76) |

The base LEB native chart-component medians were 23.38 µs for the first
body, 69.58 µs for the other ten combined and 14.75 µs for houses.
The LEB and Auto random-date profiles support measuring the decode
trade-off before prioritizing arithmetic changes. The Skyfield trace
attributes registration-list allocation and sorting to single-target
lookups in routing and SPK center selection. The lookup change now in
draft PR #87 scans live mappings for the smallest matching signed body
key, preserving the sorted list's choice even when several library
aliases share a target.

The macOS comparison interleaved reference, candidate, candidate and
reference. Each run used five fresh processes and 500 charts per case
in the 1900–2100 window. The two trials per build are shown in µs.

| Workload | Reference medians | Candidate medians |
| --- | ---: | ---: |
| Skyfield base, random charts | 1,534.7 / 1,512.8 | 916.5 / 901.8 |
| Skyfield base, consecutive charts | 1,577.8 / 1,550.9 | 942.3 / 928.6 |
| Skyfield medium, random charts | 1,618.3 / 1,595.3 | 982.6 / 967.6 |
| Skyfield medium, consecutive charts | 1,586.7 / 1,573.4 | 943.2 / 934.6 |

The Auto differences were within the observed variation of the repeated
reference runs. The 150-test SPK suite passed, including duplicate,
replaced and removed aliases. The named-satellite bit gate matched all
320 requests with the pinned fixture, including 88 successful
calculations. The isolated macOS lookup comparison against the preceding
composed build matched all 2,100,000 calls, with 350,000 per mode and tier.
The [full corrected-reference matrix](https://github.com/g-battaglia/libephemeris/actions/runs/43921283876)
on the revision with the lookup change matched all 2,100,000 calls on
Windows x86-64 and both Linux architectures. The archived comparison
reports and capture trailers were checked locally.
Raw timings and both local gates are under
`target/bench/spk-registration-lookup/`.
The latest runtime's Rust release suite passed 3,247 library tests with
the three requested features and the pinned v3 test input. Integration
and documentation tests passed. The isolated wheel's Python suite
passed 1,486 tests. The unchanged consumer and external replays are
under `target/bench/spk-registration-lookup/current-suites/`.
The SPK revision's scientific replay matched the published comparison
on all 21,359 pytest IDs, 273 standalone round statuses, 19 verification
statuses and the complete parity reports. The 269 round wrapper records
and 19 verification wrapper records also matched; four round timeouts
recurred in both runs.

The combined default wheel at `39692f97` matched all 2,100,000 macOS
calls against corrected reference `d86cea0f`. The comparator counted six
warning-emitting calls on each wheel, with no differences in their
results. The complete comparison had 350,000 calls per mode and tier.
Rust release tests on the corresponding runtime passed 3,248 library
tests after the policy-switch test was updated to check contiguous
bytes instead of requiring a heap copy. Integration and documentation
tests passed. The wheel passed all 1,486 Python tests. Raw output is in
`target/bench/combined-final-candidate/`.
The archived Linux captures for the combined
[default](https://github.com/g-battaglia/libephemeris/actions/runs/43928262541)
and [32 MiB plus positional](https://github.com/g-battaglia/libephemeris/actions/runs/43928798173)
cases each contain 2,100,000 identical requests per architecture, with
six matching warning calls. The archived default Windows captures
were also checked locally and contain 2,100,000 identical requests,
again including all six warning calls. The 32 MiB plus positional
Windows archives also contain 2,100,000 identical requests and six
matching warning calls. The archived pytest reports contain 1,454
passes and 32 skips per Python version, on both 3.12 and 3.14.
The unchanged combined-wheel consumer replays passed on Kerykeion
(13,436 passed, 1,392 skipped, three deselected), Kerykeion Parity
(51 passed), and Astrologer Studio (2,139 passed, 43 skipped).
The API v6 replay matched every outcome in the installed published
baseline across all 2,541 test IDs, including 95 setup errors and
three failures. The API v7 replay matched all 106 installed-baseline
outcomes, including one failure. The first-line exception messages
matched for all 99 non-passing cases across the two suites. The
combined scientific replay matched the published wheel on all
21,359 pytest IDs, including the same 62 failures and 147 skips.
The 273 standalone round status codes and 19 verification status
codes matched. The 269 available round wrapper records and all 19
verification wrapper records also matched. Four round timeouts
recurred in both campaigns. The complete parity result and execution
record matched. The paired audit is under
`target/bench/combined-final-candidate/validation-campaign-comparison.json`.
Raw samples, traces, build hash and pinned inputs are under
`target/bench/phase1-composed-profiles/`.

A paired macOS base/LEB run compared the Python wrapper with the native
release example using the same dates and pinned inputs. Across five
fresh-process repeats of 5,000 charts, Python/native consecutive medians
were 109.5/107.9 µs and random medians were 539.5/545.2 µs. These
measurements put the consecutive path's observed difference at about
1.6 µs per chart and show no benefit in the random-date comparison.
They do not justify a new batch API on the basis of wrapper overhead.
See `target/bench/binding-paired/` for the samples and per-call p95.

The fresh-process startup report includes import, first-call and first-chart
medians and p95s. The warm-page condition explicitly primes the staged inputs.
The early first-use reports had empty process caches and reused the staged
files between the route probe and repetitions. The revised driver takes a
fresh private copy for every first-use sample, so only staging can populate
its input pages before launch. The [first Linux startup run](https://github.com/g-battaglia/libephemeris/actions/runs/43887347685)
completed first-use and warm-page measurements on the published x86-64
wheel. The runner lacked `fincore`. The revised driver evicts only the
private staged inputs and verifies zero resident pages directly through
`mincore`. Tests cover nonzero residency, the final partial page and the
per-page presence bit. A failed cold prerequisite retains an explicit
report and the driver continues the other startup measurements.
The [revised startup workflow](https://github.com/g-battaglia/libephemeris/actions/runs/43895329805)
completed the published/rebuilt base/LEB measurements on Linux x86-64 and
ARM64. All five cold samples per build and architecture recorded zero
resident input bytes immediately before spawning the measured process.

All six default mode/tier pairs were repeated in the ARM64 Linux guest
with the revised `first_use` staging, private ext4 input files, released
companions and five fresh processes per arm. The installed candidate
and PyPI wheel used matching inputs and calculation routes. Median
first-call times in milliseconds follow.

| Mode | Base, PyPI → candidate | Medium, PyPI → candidate |
| --- | ---: | ---: |
| LEB | 2.12 → 1.97 | 5.53 → 5.68 |
| Auto | 38.23 → 30.16 | 183.68 → 150.81 |
| Skyfield | 29.41 → 17.76 | 154.06 → 97.18 |

The raw samples and p95s are in
`target/bench/linux-medium-lab/phase1-corrected-first-use-ext4.jsonl`.

The paired baseline PSS and anonymous/file splits cover
[LEB medium](https://github.com/g-battaglia/libephemeris/actions/runs/43886523936),
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

The [LEB base RAM follow-up with companions](https://github.com/g-battaglia/libephemeris/actions/runs/43903192938)
completed on both architectures with four matched LEB input hashes.
The additional Linux startup pairs cover
[LEB medium](https://github.com/g-battaglia/libephemeris/actions/runs/43903679912),
[Auto base](https://github.com/g-battaglia/libephemeris/actions/runs/43903679854),
[Auto medium](https://github.com/g-battaglia/libephemeris/actions/runs/43903679339),
[Skyfield base](https://github.com/g-battaglia/libephemeris/actions/runs/43903679363)
and [Skyfield medium](https://github.com/g-battaglia/libephemeris/actions/runs/43910066861).
The initial Skyfield medium dispatch stayed queued with zero jobs while
the cancellation endpoint reported that it had completed. The replacement
dispatch reached both Linux runners.
All six cases completed on both Linux architectures. The archived pairs
match their input manifests. Each of the 120 cold launches across both
builds reported zero input bytes resident before spawning. The newer
reports count one date, one first position and one first chart per fresh
process. The verified summary is
`target/bench/linux-startup-matrix-summary.json`.

Additional Linux startup measurements compare the published wheel against
[the SPK storage changes](https://github.com/g-battaglia/libephemeris/actions/runs/43905660043)
and [the opted-in positional reader](https://github.com/g-battaglia/libephemeris/actions/runs/43905659756).
Both use Auto medium and keep the first call and first chart separate.
All four architecture and candidate pairs completed with five fresh
processes per condition and zero resident input bytes in each cold sample.
The table shows median first call and first chart in ms, as measured
within each published/candidate pair.

| Candidate and runner | First use, call | Warm, call | Cold, call | Cold, chart |
| --- | ---: | ---: | ---: | ---: |
| SPK storage, ARM64 | 407.76 → 315.76 | 406.53 → 318.96 | 848.53 → 783.91 | 18.95 → 2.34 |
| SPK storage, x86-64 | 470.23 → 345.74 | 458.40 → 349.74 | 560.35 → 502.27 | 23.58 → 1.81 |
| Positional reader, ARM64 | 401.49 → 275.13 | 409.49 → 277.40 | 826.50 → 798.04 | 18.54 → 2.20 |
| Positional reader, x86-64 | 461.28 → 287.73 | 460.82 → 290.88 | 564.27 → 417.50 | 23.52 → 1.74 |

The paired Auto medium SPK comparison completed on both Linux
architectures with identical pinned input manifests. Steady summed PSS
for the published wheel and the SPK candidate, in decimal MB.

| Processes | ARM64 published → candidate | x86-64 published → candidate |
| ---: | ---: | ---: |
| 1 | 601.0 → 400.0 | 602.8 → 403.0 |
| 2 | 1,198.2 → 796.2 | 1,198.5 → 798.2 |
| 4 | 2,385.5 → 1,583.3 | 2,386.9 → 1,585.8 |
| 8 | 4,757.9 → 3,154.9 | 4,760.3 → 3,157.7 |

At eight processes summed anonymous RSS fell from 4,740.8 to 3,138.1 MB
on ARM64 and from 4,741.8 to 3,139.1 MB on x86-64. Summed resident file
pages were 162.2 → 163.7 MB and 193.2 → 193.8 MB respectively.

The Auto medium published/rebuilt control held steady at eight processes,
4,757.1 → 4,757.0 MB PSS on ARM64 and 4,759.6 → 4,759.6 MB on x86-64.
With the additional positional center reader enabled, the paired
comparison reduced PSS from 4,757.9 to 1,550.8 MB on ARM64 and from
4,760.2 to 1,553.1 MB on x86-64.
The paired Skyfield medium runs reduced PSS from 5,288.6 to 2,726.9 MB on
ARM64 and from 5,290.8 to 2,728.5 MB on x86-64 at eight processes.
Summed anonymous RSS fell by about 2,562 MB on both architectures.
These are separate pairs with the inputs and routing recorded in each
report.

The latest baseline PSS and RSS pairs use identical published/rebuilt
manifests within each runner. These are the steady summed PSS values in
decimal MB after 100,000 random charts per process at eight processes.
The earlier LEB base pair staged just the core. Other LEB and Auto pairs
include the released companion files for the tested tier.

| Mode and tier | ARM64 published → rebuilt | x86-64 published → rebuilt |
| --- | ---: | ---: |
| LEB base, core only | 254.1 → 254.0 | 256.4 → 256.4 |
| LEB base with companions | 257.2 → 257.2 | 259.5 → 259.5 |
| LEB medium | 470.8 → 470.4 | 472.9 → 473.1 |
| Auto base | 1,002.2 → 1,002.2 | 1,004.5 → 1,004.7 |
| Auto medium | 4,757.1 → 4,757.0 | 4,759.6 → 4,759.6 |
| Skyfield base | 1,112.2 → 1,111.2 | 1,114.2 → 1,113.6 |
| Skyfield medium | 5,287.6 → 5,287.2 | 5,289.8 → 5,290.2 |

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

The step 1b audit checked additional whole-file reads in the native
runtime. The fixed-star catalog and photometry are compiled tables.
The native JPL loader selects BSP kernels. Classic DE filename constants
remain part of the compatibility surface. The optional IERS cache does
read its three text inputs before parsing. In the local cache their sizes
were 2,460,544, 1,352 and 13,419 bytes. The sealed benchmark configurations
disable IERS loading. These paths do not explain the measured whole-image
SPK or LEB allocations.

The ordinary sealed LEB position path already requests a planetary or
lunar subset where the body and flags allow it. Reader acquisition
inspects partition metadata before deferring irrelevant body-only
partitions, with the explicit primary and required auxiliary readers
still attached. The owned LEB reader moves the file buffer into shared
storage, and bundled base uses static bytes. Removing a second full-file
copy on those paths therefore offers no equivalent to the SPK change.
External LEB images still retain their owned buffer. The indexed reader
and the authenticated positional experiment target that allocation.

The small-budget owned-reader screen uses the composed wheel with the
released LEB companions. In the completed zero-byte case, steady
single-worker RSS after 10,000 charts moved from 49.4 to 48.1 MB in base
and from 87.2 to 85.5 MB in medium. Random-chart medians moved from
580 to 17,660 µs and from 601 to 19,501 µs respectively. Releasing all
decoded coefficients between calls is excluded from default selection.
This screen shared the host with other work. The positive-budget results
follow below. Raw timings and memory snapshots are under
`target/bench/leb-owned-small-budget-screening/`.

At 1 MiB, random-chart medians were 13.7 ms in base and 15.3 ms in
medium. At 10,000 charts, steady RSS was 52.0 and 88.0 MB. The process
pool divides its budget among eight shards and cannot retain complete
chunks larger than a shard's allowance. The follow-up experiment
reuses the existing bounded coefficient windows for those chunks.
The prototype is isolated in `codex/leb-process-cache-windows`.
Its first build passed 245 LEB Rust tests and matched all 10,000 smoke
requests with the candidate's 1 MiB cache active.

The unmodified process cache at 2 MiB measured 13.4/13.6 ms median
per random chart in base/medium and 54.2/94.0 MB steady RSS. At 4 MiB,
these were 0.84/0.84 ms and 62.5/96.7 MB. Both RSS readings follow
10,000 charts on the same companion-inclusive workload.

In the separate 1 MiB prototype screen, bounded windows improved random
chart medians from 13.73 to 0.77 ms in base and from 13.60 to 0.77 ms in
medium. At 1,000 charts, steady RSS moved 49.9 → 47.7 MB and
87.3 → 85.2 MB. Consecutive chart medians were still about 2.10 ms.
This screen compares the same small configured budget on both builds.
The subsequent same-wheel sweep paired 2 and 3 MiB with legacy
retention. Medians reflect five fresh processes of 1,500 charts each.
Memory is steady single-worker RSS after 10,000 charts.

| Setting | Base random median | Medium random median | Base RSS | Medium RSS |
| --- | ---: | ---: | ---: | ---: |
| Legacy | 599 µs | 599 µs | 49.5 MB | 87.3 MB |
| 2 MiB windows | 659 µs | 663 µs | 53.4 MB | 90.8 MB |
| 3 MiB windows | 1,062 µs | 1,073 µs | 59.6 MB | 94.7 MB |

The small-budget variants fail the latency and RAM criteria for a default.
Keep this prototype outside the draft optimization chain. Local
cache-screen artifacts are under `target/bench/leb-process-cache-windows/`.

The companion-inclusive medium-only follow-up compared authenticated
positional acquisition with the ordinary owned path, using the same
experimental wheel on both sides. In a separate check the positional
path retained a descriptor for the pinned core after serving the Mars
probe. At eight processes and 100,000 charts per process, summed RSS
fell from 671.3 to 615.1 MB and summed sampled physical footprint fell
from 520.4 to 456.5 MB. At one worker RSS fell from 84.4 to 70.4 MB.
The same-wheel random-chart median rose from 609.1 to 637.7 µs and
p95 fell from 711.8 to 683.6 µs. The consecutive medians were 115.8
and 121.2 µs. These separate timing samples used five fresh processes
of 3,000 charts each over the covered part of 1550–2650. The digest
check also increased first-call time in the earlier startup trial.
Keep automatic positional acquisition experimental.
Raw results are in `target/bench/managed-leb-medium-follow-up/`.

The full owned-reader RAM sweep measured legacy retention and 16/32/64/128
MiB on both LEB tiers after the local timing and profiling commands.
All 40 cases passed the local build, data, workload and checkpoint
audit. Each used 100,000 charts per worker at 1/2/4/8 processes.
Commands and artifacts are under `target/bench/owned-cache-ram-matrix/`.
The sweep reports macOS process measurements separately from Linux PSS.
A [hosted LEB medium run](https://github.com/g-battaglia/libephemeris/actions/runs/43918300315)
compared the published release with the 32 MiB candidate using Linux PSS
at the same four process counts. Both architectures passed the covered
10,000-call smoke and matching-data checks. At eight processes summed
PSS increased from 475.2 to 671.2 MB on ARM64 and from 477.4 to
673.8 MB on x86-64. The increase largely appears in anonymous RSS.

RSS and sampled physical footprint increased at every process count in
both tiers for all four budgets. Values below represent eight processes,
in decimal MB. The sum of each process's lifetime peak RSS equaled the
steady summed RSS in these eight-process cases. Peaks need not coincide.

| Pool budget | Base summed RSS | Base summed footprint | Medium summed RSS | Medium summed footprint | Retained coefficients per worker |
| --- | ---: | ---: | ---: | ---: | ---: |
| Legacy | 458.9 | 277.5 | 758.5 | 608.6 | — |
| 16 MiB | 843.6 | 485.0 | 1,138.1 | 816.1 | 16.07 |
| 32 MiB | 627.8 | 502.0 | 921.6 | 796.1 | 24.72 |
| 64 MiB | 626.0 | 500.3 | 920.9 | 795.4 | 24.72 |
| 128 MiB | 625.9 | 500.2 | 924.2 | 798.7 | 24.72 |

The 32 MiB budget was the smallest that reached the narrow-window latency
plateau. Increasing it retained no extra decoded data in this workload.
The existing default meets the lower-RAM preference for the cache choice.
Lower RAM than Phase 0 for default LEB remains a separate requirement.
The 16 MiB single-process base `vmmap` sample shows regions classified as
`MALLOC_SMALL (empty)` with 2.8 → 24.7 MiB resident. This is allocator
residency outside the decoded retention count. It does not identify
which allocations caused it. An isolated owned-reader prototype reuses
the bounded temporary output buffer on the existing direct Zstd path.
It keeps the indexed working-budget path and streaming fallback. Its
243 LEB Rust tests and the 10,000-call smoke comparison passed. The
paired macOS screen used five processes of 1,500 charts for timing and
one process of 100,000 charts for memory. With legacy retention, random
medians rose from 548.3 to 568.0 µs in base and from 552.3 to 566.5 µs
in medium. RSS moved from 57.7 to 57.0 MB and from 95.1 to 94.6 MB,
respectively. Physical footprint grew slightly in both cases.
At 16 MiB, random medians rose from 216.3 to 225.6 µs in base and from
217.4 to 232.6 µs in medium. Base RSS increased from 74.8 to 98.1 MB,
while medium sampled footprint increased from 93.1 to 106.3 MB.
The experiment is excluded from the draft. Raw timing, memory and
correctness reports are under `target/bench/leb-decode-buffer/`.

A second isolated candidate keeps ordinary ownership semantics for
discovered base-tier data. It compares the opened `base_core.leb2`
stream with the static bytes already bundled in the executable, through
EOF. An exact match uses those immutable bytes with the discovered
pathname and a fresh cache. Other inputs keep the ordinary acquisition
path, including explicit files and indexed readers. This could reduce
physical memory across processes without relaxing the immutable snapshot
contract. The 245 bundled LEB tests passed, including an exact match,
same-size replacement and explicit-file control. The macOS smoke gate
matched 10,000 calls. After 100,000 charts per process, the paired
measurements gave the following steady sums in decimal MB.

| Processes | RSS, reference → snapshot | Sampled footprint, reference → snapshot |
| ---: | ---: | ---: |
| 1 | 57.6 → 57.7 | 38.0 → 24.9 |
| 2 | 113.2 → 115.2 | 70.6 → 59.6 |
| 4 | 230.5 → 231.1 | 139.1 → 102.7 |
| 8 | 459.3 → 459.3 | 276.0 → 206.8 |

At eight workers the sampled footprint after the first chart was
308.0 → 225.5 MB. Summed RSS includes the executable's resident pages
in each process. In the [Linux follow-up](https://github.com/g-battaglia/libephemeris/actions/runs/43924104899),
ARM64 steady summed PSS changed from 38.5 / 72.3 / 133.4 / 258.8 MB
to 37.2 / 61.2 / 104.0 / 186.0 MB at 1 / 2 / 4 / 8 workers. At
eight workers anonymous RSS fell from 241.4 to 159.5 MB. The x86-64
steady summed PSS moved from 40.4 / 73.7 / 137.3 / 261.0 MB to
39.1 / 62.1 / 105.1 / 187.7 MB at the same process counts. At eight
workers anonymous RSS fell from 242.4 to 160.6 MB. Both Linux runs
matched all 10,000 smoke cases. The base snapshot change is in PR #87
for [combined default qualification](https://github.com/g-battaglia/libephemeris/actions/runs/43928262541)
and [combined 32 MiB and positional qualification](https://github.com/g-battaglia/libephemeris/actions/runs/43928798173).
Mac random-chart medians were 541.7 → 576.6 µs and consecutive medians
were 106.8 → 111.4 µs across five independent repeats. The median first
call moved from 1.49 to 1.92 ms on first use and from 1.45 to 1.95 ms
with warm pages.

The follow-up alternated reference, candidate, candidate and reference
with five repeats of 2,400 charts per run. Random-chart medians were
583.4 / 585.3 µs for the reference and 584.0 / 586.1 µs for the
candidate. Consecutive medians were 112.1 / 113.2 and 111.8 / 111.2 µs.
The initial steady-state slowdown did not recur. Warm-page first-call
medians were 1.54 / 1.56 ms on the reference and 1.90 / 1.85 ms on the
candidate. The startup cost remains visible during qualification.
Raw results are under `target/bench/leb-static-base-snapshot/`.

The medium-tier investigation adds a local lab probe for a file-backed
snapshot. It opens a separate APFS clone of a verified staged medium
core in a private directory. It maps the clone read-only and removes
its path before publishing the reader to the runtime. It retains no
writable alias. The source and clone are distinct inodes with
copy-on-write data. [Apple clonefile documentation](https://github.com/apple-oss-distributions/xnu/blob/main/bsd/man/man2/clonefile.2)

A native core-only preflight matched position words, returned flags
and house output over 1,000 random charts. At one process, steady RSS
was 48.2 → 19.0 MB and sampled footprint was 41.5 → 5.4 MB. The existing
mapped reader used a different decoded-window policy in that preflight.

The matched-window experiment ran 100,000 random charts per worker
on the explicit medium core and matched the bitwise output digests at
every process count. The figures are steady sums in decimal MB.

| Processes | RSS, owned → clone | Sampled footprint, owned → clone |
| ---: | ---: | ---: |
| 1 | 55.0 → 24.8 | 45.7 → 7.2 |
| 2 | 109.6 → 49.6 | 88.7 → 13.9 |
| 4 | 220.2 → 99.5 | 180.8 → 27.6 |
| 8 | 441.5 → 199.9 | 358.4 → 54.4 |

Truncating a separately staged source after acquisition preserved all
1,200 native charts for both the owned and clone readers. All output
digests matched the unmodified-source control.

The isolated wheel exercises the same snapshot through ordinary
discovery with the companion files installed. It matched 10,000
requests including the sidereal warning. After 3,000 charts per worker,
eight-process summed RSS was 618.6 → 375.2 MB and sampled footprint
was 493.8 → 194.8 MB. In this screen, random-chart medians were
581.3 → 582.3 µs and consecutive medians were 113.8 → 113.5 µs.
The warm-page first-call median was 3.85 → 2.35 ms and the subsequent
first-chart median was 0.60 → 2.15 ms.

The discovery comparison completed 100,000 charts per worker with
the same route probes on both arms. Steady sums are in decimal MB.

| Processes | RSS, owned → clone | Sampled footprint, owned → clone |
| ---: | ---: | ---: |
| 1 | 84.7 → 54.5 | 67.2 → 24.7 |
| 2 | 169.3 → 108.9 | 130.5 → 51.7 |
| 4 | 338.5 → 217.9 | 264.7 → 102.0 |
| 8 | 667.6 → 430.4 | 522.7 → 197.2 |

Random medians were 582.0 µs for both arms and consecutive medians
were 114.0 → 113.4 µs. The warm-page first-call median moved from
3.81 to 2.47 ms and the subsequent first-chart median from 0.59 to
2.19 ms. The experiment assumes exclusive access to its scratch
directory from clone creation through unlink. Whether the runtime
can establish that invariant for automatic discovery remains open.
Default medium RAM remains open. Raw files are under
`target/bench/leb-medium-private-snapshot/`.

A subsequent local wheel (`3a1bb6f5`) tried a fresh scratch directory
automatically and used the original owned reader when the scratch
directory was unavailable. The two paths matched all 10,000 requests
and the warning. With 3,000 charts per worker, eight-process summed
RSS was 618.7 → 376.1 MB and sampled footprint was 492.9 → 194.7 MB.
The test found no scratch entries left behind after successful acquisition.
The wheel remains a lab artifact pending an immutability decision.

A second macOS prototype (`cb02fa89`) validates the digest of the
private clone against the pinned medium release, keeps an unlinked
read-only descriptor, and uses bounded positional reads with a 6 MiB
coefficient cache. Explicit contiguous byte borrows materialize and
verify the image on demand. Acquisition falls back to the owned reader
when the clone or its pinned digest is unavailable. The two paths
matched all 10,000 requests and the warning. At eight processes after
3,000 charts per worker, summed RSS was 618.3 → 446.7 MB and sampled
footprint was 492.7 → 320.6 MB. Over the 1900–2100 workload, random
median/p95 were 581.8/614.6 → 519.5/648.4 µs and consecutive
median/p95 were 113.9/118.9 → 118.6/124.0 µs. Warm-page first-call
median was 3.67 → 15.95 ms. Native tests verified that contiguous
borrows allocate once across reader clones, that the admitted clone
survives rewriting and truncating its original pathname, and that
another same-length file with a wrong digest falls back to the
ordinary reader. All 350,000 medium/LEB calls against the ordinary
combined runtime at `39692f97` matched, including the existing warning.

On the full common 1550–2650 coverage, the same wheel ran 100,000
charts per worker. The table uses steady sums in decimal MB.

| Processes | RSS, owned → positional snapshot | Sampled footprint, owned → positional snapshot |
| ---: | ---: | ---: |
| 1 | 84.9 → 65.4 | 65.4 → 31.4 |
| 2 | 169.8 → 127.6 | 130.9 → 72.1 |
| 4 | 339.7 → 157.0 | 261.8 → 198.4 |
| 8 | 668.9 → 523.3 | 525.5 → 238.9 |

At eight workers, the sums of individual process peaks were
668.9 → 523.3 MB RSS and 543.3 → 396.5 MB footprint. Peaks across
different processes need not be simultaneous. macOS does not report PSS.
On full coverage, random median/p95 were 585.4/613.0 → 620.9/654.6 µs
and consecutive median/p95 were 114.1/119.8 → 118.4/124.0 µs.
Warm-page first-call median was 3.83 → 16.15 ms. Full-coverage latency
differs from the shorter 1900–2100 workload reported above. The
prototype remains local while a matched-cache comparison isolates
the storage and cache-policy effects.

The matched-cache wheel (`e526c136`) uses the ordinary 1 MiB
coefficient cache setting. Its 3,000-chart screen sampled full
coverage and matched 10,000 calls bit for bit, including the warning.
At eight processes, summed RSS was 618.6 → 321.6 MB and sampled
footprint was 493.1 → 195.7 MB. Random-chart median/p95 were
590.4/619.4 → 619.5/651.5 µs and consecutive median/p95 were
114.3/119.1 → 119.9/124.8 µs. The warm-page first-call median was
3.70 → 16.10 ms. It also matched 350,000 calls against the combined
runtime with the same warning. Over 100,000 charts per worker on full
coverage, the steady eight-process RSS was 672.8 → 379.8 MB and
sampled footprint was 523.2 → 197.1 MB. At one, two and four workers
the footprint was 65.3 → 24.9, 131.0 → 49.2 and 261.9 → 100.3 MB.
The sum of individual eight-process peak footprints was 547.0 →
253.7 MB. Full-coverage random median/p95 changed from
608.3/685.9 to 637.2/708.0 µs and consecutive median/p95 from
114.6/120.2 to 119.8/125.3 µs. Warm-page first-call median was
3.83 → 16.32 ms.

Qualification found that an explicit 32 MiB process cache could
panic while attaching the positional snapshot's shared cache.
The revised local runtime (`a4ca7bd4`) only attempts the snapshot
through automatic discovery with the legacy cache. Direct opens and
explicit process-cache settings follow their previous acquisition path.
The explicit 32 MiB probe and 10,000-call comparison against the
combined runtime passed, including the warning. The corrected route
matched all 350,000 default-configuration calls. At first call, the
resident malloc regions in the automatic acquisition probe measured
8.2 MB versus 45.5 MB with a forced scratch fallback.

The macOS qualification candidate retains the 1 MiB policy and
preserves the ordinary reader diagnostics. A verified clone is admitted
only for the pinned medium core discovered in LEB mode. A failed scratch
or verification attempt leaves acquisition and errors to the ordinary
reader. A contiguous byte borrow verifies and materializes the image
on demand. The reduced steady and peak footprint is the reason to
qualify it for the default despite the measured first-use and chart
latency cost. Default medium memory on the other platforms remains
open.

The local wheel at `17c32b93` matched 350,000 LEB medium calls in
the default configuration and 10,000 with the explicit 32 MiB pool,
including the warning in both comparisons. The automatic discovery
probe began with 8.45 MB of resident malloc regions after the first
call, versus 45.60 MB on the scratch fallback and 46.10 MB with an
explicit path. The Rust suite passed 3,248 library tests plus
integration and documentation tests. The three data-dependent
snapshot tests also passed separately. The wheel passed 1,486 Python
tests with the repository's bundled base fixture, as used in the
wheel workflow, and all 18 harness tests. The out-of-range inventory
test has the same existing reviewed-label failure on the preceding
wheel and this wheel when pointed at the separately installed v4
base core. Neither the installed data nor the test was modified.

The interleaved five-repeat screen comparing the matched-cache
prototype directly with PyPI 4.0.0 measured random and consecutive
medians of 602.9 → 633.1 and 114.9 → 120.0 µs over full coverage.
The same wheel's scratch fallback measured 600.7 and 114.5 µs.
Warm-page first-call medians were 3.79 ms for PyPI and 3.76 →
20.07 ms in the fallback/snapshot comparison.

The exact macOS wheel at `17c32b93` matched all 2,100,000 calls
against `d86cea0f` in the full corrected-reference corpus. Each
mode/tier pair had 350,000 calls and the six warning-emitting calls
matched. The capture headers and trailers were verified along with
the comparison report. The [three-platform PR smoke](https://github.com/g-battaglia/libephemeris/actions/runs/43949666788)
also matched 10,000 calls and the warning against the published wheel
on Windows and both Linux architectures. The PR merge commit and the
branch head for that run had identical Git trees. The PR smoke skips
the installed-wheel pytest suites; the full manual matrix runs them.

The final macOS wheel repeated the same-wheel fallback/snapshot RAM
comparison on full coverage at 100,000 charts per worker.

| Processes | Steady RSS, fallback → automatic | Sampled footprint, fallback → automatic |
| ---: | ---: | ---: |
| 1 | 84.9 → 47.8 MB | 65.4 → 24.9 MB |
| 2 | 170.0 → 92.6 MB | 131.2 → 51.2 MB |
| 4 | 335.7 → 191.7 MB | 260.4 → 97.6 MB |
| 8 | 670.4 → 375.9 MB | 530.4 → 194.2 MB |

At eight workers, the sum of individual peak footprints moved from
544.9 to 250.0 MB. Random median/p95 moved from 597.0/626.8 to
626.2/657.2 µs and consecutive median/p95 from 114.4/119.2 to
119.2/124.2 µs. Warm-page first-call median was 3.66 → 16.19 ms.

The direct comparison staged the same inputs and dates for the
published PyPI 4.0.0 wheel and the exact `17c32b93` wheel. At eight
processes and 100,000 full-coverage charts per worker, steady summed
RSS fell from 674.1 to 380.5 MB and sampled footprint from 526.4
to 198.3 MB. First-call sampled footprint was 520.5 → 223.3 MB.
The sum of individual peak footprints fell from 550.2 to 254.5 MB.
The wheel hashes, routes and four checkpoints are recorded in
`target/bench/leb-medium-private-snapshot/positional-production/published-memory.json`.
The final wheel also repeated the unchanged downstream suites. Outcomes
matched the preceding combined wheel on all 14,828 Kerykeion test IDs,
51 Kerykeion Parity IDs and 2,182 Astrologer Studio IDs. The Studio
comparison matches the five randomly generated UUID cases by format
and version. The API v6 and v7 replays matched all 2,647 outcomes on both
the preceding wheel and the published wheel, including the messages of
the 99 existing non-passing cases. The exact wheel's external scientific
replay has matched all 21,359 published pytest outcomes, with identical
fixture digests and harness settings. The same comparison matched all
273 standalone round statuses, all 19 verification statuses and both
parity reports. The wrapper evidence matched for the 269 rounds that
produced it and all 19 verification cases. Four shared round timeouts
had no wrapper evidence. The paired report is
`target/bench/leb-medium-private-snapshot/positional-production/validation-campaign-comparison.json`.

A local warm-page acquisition profile (`15e92dbb`) measured one
private core acquisition per process across 12 repeats for each digest
buffer size. With the existing 64 KiB buffer, median clone acquisition
was 0.28 ms, full digest verification 14.95 ms, index parsing 0.02 ms,
and the complete first call 15.94 ms. At 256 KiB, verification
was 14.33 ms and the first call 15.32 ms. At 1 MiB these were
14.26 and 15.25 ms. The draft retains the 64 KiB buffer. The
remaining first-call cost comes primarily from reading and checking
the complete pinned digest. Raw samples are in
`target/bench/leb-medium-private-snapshot/open-profile/`.

The [combined indexed qualification](https://github.com/g-battaglia/libephemeris/actions/runs/43929882577)
also measured 100,000 charts per process on Linux ARM64 with the
final candidate and a 6 MiB indexed budget. At 1 / 2 / 4 / 8 workers
the steady summed PSS for `.leb` medium fell from 77.3 / 138.9 /
257.5 / 492.7 MB to 34.4 / 64.2 / 122.0 / 232.5 MB. For base it
moved from 39.1 / 60.8 / 104.6 / 186.3 MB to 34.4 / 65.1 / 122.4 /
232.4 MB. The existing bundled base reader uses less PSS at multiple
processes in this comparison. The x86-64 memory run confirmed the result
at eight workers with 493.6 → 235.5 MB in medium and 187.5 → 236.3 MB
in base. The archived ARM64 and x86-64 captures each matched 2,100,000
corrected-reference calls in the indexed and positional configuration,
with six warning-emitting calls on each side. The verified Windows
archives also contain 2,100,000 identical requests and six matching
warning calls. The archived Windows pytest reports contain 1,454
passes and 32 skips on both Python versions.

The [final default matrix](https://github.com/g-battaglia/libephemeris/actions/runs/43949734696)
has qualified Linux ARM64 at `0319bacb`: all 2,100,000 calls matched the
corrected reference and all six warning-emitting calls matched. Capture
headers, bodies and trailers were checked. Installed-wheel tests on
Python 3.12 and 3.14 each had 1,454 passes and 32 skips.
On the same data and 100,000 random medium charts per worker over
1900–2100, steady summed PSS for 1 / 2 / 4 / 8 workers was
67.23 / 126.74 / 242.40 / 475.03 MB on the reference and
64.67 / 124.77 / 243.13 / 474.47 MB on the candidate. The
eight-process default footprint was effectively unchanged.
The x86-64 job's comparison report also records 2,100,000 identical
calls and six matching warnings. Its 1 / 2 / 4 / 8-worker reference
PSS was 67.42 / 127.80 / 245.48 / 477.60 MB versus 67.47 / 127.78 /
245.52 / 477.51 MB on the candidate. Both reports were read from the
job archive. All x86-64 capture records and trailers were checked.
The x86-64 installed-wheel tests also had 1,454 passes and 32 skips
on each Python version.
The final Windows job completed successfully. Its decoded archive at
`0319bacb` also matched all 2,100,000 corrected-reference calls and
six warning-emitting calls. Both Python versions had 1,454 passes and
32 skips. The strict published-wheel Windows gate remains open.

The Linux lab uses an isolated ARM64 guest and the same installed
candidate wheel. Its guest filesystem declined reflink creation.
A copied unnamed snapshot retained the exact pinned digest after
the isolated source was truncated. The owned, direct indexed and
private indexed paths matched position, house and flag bits over
1,000 full-coverage charts. The private path moved process PSS from
68.3 to 37.8 MB, while `mincore` found an additional 32.1 MB of
resident cache pages in the private file. These unmapped pages are
outside process PSS. This screen does not establish a physical-memory
gain on Linux. Raw evidence is under `target/bench/linux-medium-lab/`.

With 10,000 charts per path, the private indexed copy kept 7.7 MB
of file pages on 1900–2100 and 35.5 MB on full coverage. Dropping
pages every 64 charts retained bitwise parity and reduced those
steady samples to 3.0 and 4.7 MB. The full-coverage batch average
rose from 610 µs on owned storage to 947 µs with periodic eviction.
The lab moved on to sharing an immutable image across processes.

A separate local screen used a Linux memfd with write, shrink, grow
and seal changes disabled by the kernel. Eight freshly executed
workers opened and independently verified the same inode. All
position, flag and house bits matched the ordinary reader over
3,000 full-coverage charts per worker. Summed worker PSS was
426.7 → 165.8 MB. The container's sampled `anon + file` charges
were 499.6 → 247.9 MB, including the experiment's owner process.
The 37.3 MB shared segment is included in `file`. The sealed arm
uses the stock indexed evaluator with a 1 MiB cache and separately
holds a read-only mapping to account for its shared pages.
[Linux file seals](https://man7.org/linux/man-pages/man2/F_ADD_SEALS.2const.html)
and [cgroup memory counters](https://docs.kernel.org/admin-guide/cgroup-v2.html#memory-interface-files).

At eight workers and 100,000 charts per worker over full coverage,
the owned and sealed chart digests matched. Worker PSS changed from
430.1 to 166.3 MB and the container's `anon + file` charges from
503.1 to 248.5 MB. The median of per-worker batch averages rose
from 1,295 to 1,387 µs per chart. The mapped inode retained its
pinned digest and all seals after the initial owner exited and
another worker continued to hold it.

The controller supplied the owner descriptor for the 100,000-chart
repeats. All worker digests matched the owned path at each process
count. The table reports decimal MB and includes the owner process
in the sealed cgroup charges. `file` includes shared memory.

| Workers | Owned → sealed summed worker PSS | Owned → sealed cgroup `anon + file` |
| ---: | ---: | ---: |
| 1 | 65.1 → 62.8 | 142.3 → 149.9 |
| 2 | 119.3 → 78.8 | 194.7 → 164.3 |
| 4 | 225.2 → 109.2 | 295.4 → 160.9 |
| 8 | 430.1 → 166.3 | 503.1 → 248.5 |

The four-worker repeat overlapped with the prototype build. Its
timings are excluded from the latency evaluation.

A local native candidate (`fe894aaf`) built on Linux ARM64 and used
the ordinary contiguous cache for both storage paths. Admission
requires the kernel seals and an exact match to the original pinned
source. With an explicit descriptor supplied by the lab and 3,000
charts per worker, eight workers used the same sealed inode and
retained identical position, flag and house bits. Worker PSS was
431.6 → 163.4 MB and cgroup `anon + file` was 504.6 → 245.9 MB. An
unsealed file with the exact pinned bytes fell back to owned storage
and produced the same digest.

A separate 2,000-chart, eight-worker run used a Python lab registry
to share the descriptor among independently started workers. Exactly
one worker created the image. After it exited, a replacement reused
the inode from a live peer and returned the same chart digest. With
the native evaluator, summed worker PSS was 433.1 → 168.9 MB and
cgroup `anon + file` was 507.5 → 243.2 MB. This lab has no dedicated
owner process. Running the same check as UID/GID 65534 in the guest
also confirmed one creator, reuse of the same inode after that worker
exited, and identical chart digests. PSS was 444.6 → 191.0 MB and
cgroup `anon + file` was 479.8 → 226.3 MB in that run.

The native 100,000-chart repeat used eight processes per arm and an
additional descriptor owner for the sealed arm. All worker digests
matched. Summed worker PSS was 428.6 → 163.3 MB, including 37.3 MB
of shared memory in the sealed arm. The cgroup sample was 501.6 →
245.7 MB for `anon + file`, including the descriptor owner. The
median of per-worker batch averages was 1,275 → 1,269 µs per chart
in this single run. The first-call median was 13.4 → 34.9 ms. The
fallback also matched when the descriptor was absent or a correctly
sealed image of the expected length differed near the end from the
pinned original. The subsequent automatic-discovery qualification is
documented later in this section.

An unprivileged repeat with 10,000 full-coverage charts per worker
used the Python registry and native evaluator, with no separate
owner process. All chart digests matched. Sampled cgroup charges
(`anon + file`, MB) for 1 / 2 / 4 / 8 workers were
143.2 / 196.7 / 303.5 / 517.3 → 143.8 / 160.6 / 194.2 / 261.4.
At one worker, the median of per-worker batch averages was about
601 µs for both paths. First-call timing inside the sealed worker
rose from 6.7 to 17.1 ms. The Python registry added a separate
22.5 ms including redundant hashing before the Rust verification.
This supports targeting multiprocess deployments and retaining the
ordinary path for single-process defaults until the acquisition
cost can be reduced. Raw files are in `target/bench/linux-medium-lab/`.

Five further unprivileged repeats per count (one and eight workers,
2,000 charts each) removed the duplicate Python digest while retaining
the Rust admission check. All chart digests matched. At one worker,
the median registry stage fell to 8.4 ms across these repeats. The
ordinary and sealed first-call medians, measured after the registry
stage, were 7.1 and 17.6 ms. At eight workers, the median first-call
time across the per-run medians was 14.6 → 35.3 ms, again excluding
the lab registry stage. A native opt-in registry prototype is in
qualification. It removes the external coordination step by publishing
live descriptors from the library, and preserves ordinary acquisition
if the private directory or sealed descriptor cannot be used.
The first installed-wheel run (`1190ed38`, ARM64, UID 65534) matched
all eight worker digests over 3,000 full-coverage charts each. One
worker created the image. A replacement reused it after that worker
exited. After every holder exited, a fresh worker recreated it. A
missing, public or symlinked directory, a busy directory lock and an
incorrect sealed image all fell back to the owned path with identical
results. Worker PSS was 443.0 → 182.6 MB, and cgroup `anon + file`
was 514.9 → 254.4 MB. First-call median was 11.1 → 31.0 ms. The
100,000-chart repeat matched all worker digests and the lifecycle
and fallback checks. PSS was 443.4 → 182.9 MB, and cgroup
`anon + file` was 515.5 → 254.9 MB. The median of per-worker
batch averages was 1,230 → 1,214 µs per chart in this single run.
First-call median was 13.6 → 30.9 ms.

PR #87 now includes the opt-in `LIBEPHEMERIS_LEB_SHARED_DIR` path.
It targets sealed LEB medium with default owned coefficients, based
on the multiprocess memory gain and the startup tradeoff. Workers
must share a private runtime directory and an accessible `/proc`
namespace. The library holds a bounded startup lock while discovering
or creating the immutable image. Every worker checks its source
against the pinned bytes. Host-local qualification observed correct
fallback on missing or unsafe settings and a bounded busy-lock wait.
The new `leb-medium` numerical profile matched 10,000 requests and
one warning-emitting call against the preceding wheel. The standard
benchmark verified the same mapped inode across 1 / 2 / 4 / 8 worker
counts. The [full hosted opt-in run](https://github.com/g-battaglia/libephemeris/actions/runs/43965840168)
includes both Linux memory comparisons. A repeatable installed-wheel
lifecycle check uses independent processes and exercises creator exit,
recreation after the last holder, an occupied lock and incorrect sealed
bytes. It passed on the unprivileged ARM64 candidate and rejected the
preceding wheel with sharing unavailable. The
[published-wheel medium check](https://github.com/g-battaglia/libephemeris/actions/runs/43966964284)
passed the lifecycle scenarios on Linux ARM64 and x86-64. The archived
captures and local replay matched 10,000 requests and one warning call
against PyPI 4.0.0 on each architecture. Python 3.12 and 3.14 each had
1,454 installed-wheel passes and 32 skips on both Linux jobs. The same
run passed the Windows medium check against PyPI without the Linux
opt-in, with 10,000 matching calls, one warning call and the same Python
test counts. Windows Rust test targets type-checked successfully.

Startup on the same ARM64 guest wheel was also measured without an
existing shared image, using five fresh processes per condition. The
shared reader created an image on each first call. Medians in ms were
6.55 → 20.23 at first use, 7.64 → 20.87 with warm pages, and
10.06 → 24.68 with verified cold pages. Both cold arms recorded zero
resident input bytes before each process. In an independent seven-pair
warm-page run with a live owner, the first call was 7.23 → 18.62 ms.
Every new opt-in worker mapped the owner's inode. The explicit-file,
indexed-reader and configured-process-cache controls each kept their
ordinary path and returned the same position bits. Raw files are
`native-registry-startup-pairs.json` and `shared-live-startup-pair.json`
under `target/bench/linux-medium-lab/`.

The macOS wheel rebuilt after the Linux integration passed 1,486 Python
tests. Rust on the same source passed 3,248 library tests plus integration
and 24 documentation tests when configured with the pinned data-v3
fixture. The rebuilt wheel at `42acb11a` also matched the corrected
reference on 2,100,000 requests across all six mode/tier pairs, including
all six warning-emitting calls. Reports are under
`target/bench/macos-after-linux-sharing/` and
`target/bench/linux-medium-lab/macos-post-integration-rust-pinned.log`.

In the local ARM64 guest, the first parallel release test build was
killed by SIGKILL. With one build job, LTO disabled and four codegen
units, the Rust 1.100.0 library suite ran to completion on both sides of
the Linux sharing change. Both builds had 3,246 passes, 18 ignored
tests and identical failures in the obliquity/nutation, heliacal and
topocentric pinned-value tests. The observed and expected values match
between `6cd28816` and `54c3547f`. The failures are retained as a guest
qualification limitation. Raw logs and the word-for-word comparison are
under `target/bench/linux-medium-lab/`.

## 15. Step 1a experiment

An opt-in process pool retains whole decoded chunks across runtime-owned LEB
images. Each shard evicts its least recently used entries by byte capacity.
Explicit indexed readers retain their existing resource policy. Legacy
retention remains the default based on the combined latency and RAM
measurements recorded in Section 14.

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
- [x] Qualify the 32 MiB setup against the corrected reference on macOS,
  both Linux architectures and Windows. Section 14 links the full runs.
- [ ] Pass strict published-release qualification on Windows. The
  [original cache trial](https://github.com/g-battaglia/libephemeris/actions/runs/43863141515)
  includes the released wheel's Auto instability.
- [x] Measure 100,000-chart process memory and 1/2/4/8 thread throughput.
- [x] Re-profile cold misses and measure decoder context and output-buffer
  reuse. Both experiments remain outside the draft.
- [x] Select legacy retention as the default after the paired latency and
  memory measurements. Keep 32 MiB available as an explicit option.

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

## 16. Step 1b opening-copy experiment

The first storage change moves the Vec received by the SPK parser into a
shared owner. Acquiring a kernel no longer copies the whole Vec into an Arc
slice. This preserves the public immutable byte view and record evaluator.
The regression test checks that the parser and its clones retain the original
buffer address. The earlier digest and weak-owner test still passes.

Mac arm64, medium tier, controlled 1,000-chart samples, default LEB retention:

| Mode | After first call before → after | After 1,000 charts before → after |
| --- | ---: | ---: |
| auto | 535.5 → 335.1 MB RSS | 676.6 → 474.8 MB RSS |
| skyfield | 794.6 → 474.1 MB RSS | 996.0 → 675.5 MB RSS |

This removes the transient full-file copy and the allocator pages it leaves
resident. The owned kernel image and its decoded segments still consume RAM.
Positional SPK records and lazy acquisition remain part of step 1b.

The opening-copy wheel (commit `abd6e262`, wheel SHA-256
`c06461ae417708b565c1eb0121c142cfefe466ba76fd6feb77778abcd06d8348`)
matched the 4.0.0 wheel on all 2,100,000 macOS arm64 corpus calls with 32 MiB
of LEB coefficient capacity. The capture is under
`target/bench/spk-copy-qualification/full/`.

### Decode only the selected SPK record

Registry evaluation now retains type 2 and 3 directory results and decodes
the selected record from the shared immutable image. Each installed image
holds at most 128 cached records and 1 MiB of decoded coefficient words.
Active references, allocation overhead, the raw image and type 21 segments
are accounted for separately. The 153 SPK tests passed, including a regression
that compares eager and on-demand evaluation bits across record joins and at
the closed endpoints.

Mac arm64, medium tier, 1900–2100, legacy LEB retention. One process per
memory sample with 5,000 random charts. Latency is the pooled median over
five independent processes with 5,000 charts each.

| Mode | RSS after first call | RSS after 5,000 charts | Physical footprint after 5,000 charts | Median per chart |
| --- | ---: | ---: | ---: | ---: |
| auto | 335.3 → 272.3 MB | 475.7 → 274.6 MB | 460.5 → 259.2 MB | 1,172 → 1,189 µs |
| skyfield | 474.4 → 354.6 MB | 675.8 → 355.4 MB | 660.9 → 340.3 MB | 1,555 → 1,591 µs |

The comparison uses the opening-copy wheel and the on-demand record wheel
(commit `7009b012`, SHA-256
`4c80d771e503c043603365298476fd424f31d5c504ee5778d64fc33336db8385`).
Raw reports and vmmap output are under
`target/bench/spk-copy-qualification/`. The new wheel's 2,100,000-call
comparison passed with zero differences with the configured LEB cache of
32 MiB. The combined wheel (commit `e76dd9f7`) passed all 1,486 Python tests
with that setting. Its 10,000-call smoke comparison passed locally and the
PR wheel checks passed on Linux x86-64, Linux ARM64 and Windows x86-64. The full Rust release
suite passed 3,241 library tests plus integration and documentation tests
with `bundled-data,nbody,python-internals` and the pinned data-v3 fixture.

The sealed external planet-center and planet-centric selections ran in
`auto` and `skyfield` with both wheels. Each run had 34 passes and two
matching import failures in the legacy planet-map introspection tests.
The complete scientific comparison is recorded in section 14.

### 100,000-chart memory comparison

Mac arm64, medium tier, one process, random dates in 1900–2100 with seed
20261116 and legacy LEB coefficient retention. The published 4.0.0 wheel and
the combined SPK wheel at `e76dd9f7` used the same pinned files and requests.
Values are decimal MB with one run per configuration.

| Mode | RSS after first call | RSS after 100,000 charts | Peak RSS | Physical footprint after 100,000 charts |
| --- | ---: | ---: | ---: | ---: |
| leb | 79.2 → 79.5 | 84.0 → 84.2 | 84.0 → 84.2 | 66.8 → 64.9 |
| auto | 543.3 → 279.7 | 485.2 → 284.7 | 685.4 → 284.7 | 466.2 → 269.4 |
| skyfield | 802.1 → 362.1 | 1,002.0 → 362.9 | 1,002.0 → 362.9 | 987.3 → 347.8 |

In a separate same-wheel comparison, the existing indexed LEB reader lowered
base RSS from 57.0 to 51.5 MB and medium RSS from 94.3 to 52.4 MB after
100,000 charts. Physical footprint moved from 35.5 to 24.4 MB and from 74.6
to 26.4 MB respectively. The raw reports and vmmap captures are under
`target/bench/spk-records-qualification/`. The indexed reader's latency
and the process and thread matrix are still open before selecting a default.

The next storage experiment is positional I/O for managed data under the
immutable-data contract. The public owned-image path supplies callers with
a stable byte slice for the life of the kernel.

## 17. Step 1b positional automatic-center prototype

`LIBEPHEMERIS_EXPERIMENTAL_POSITIONAL_CENTERS=1` selects the prototype for
scalar `auto` calls. It defers acquisition of a published planet-center
supplement until the first offset request. The open descriptor is scanned
against the published SHA-256 pin before its summary index can be used.
Record offsets and split-time Chebyshev evaluation share the existing code.
Explicit image and complete session APIs retain their owned byte views.
Native DE reads and PCTR continue through their ordinary acquisition paths.

The first prototype used an identity check before and after every offset
evaluation. On macOS arm64, medium tier and legacy LEB retention, on the
same wheel (`503ae7a5`) and data:

| Metric | Owned | Positional, checks on every sample |
| --- | ---: | ---: |
| RSS after 5,000 random charts | 274.2 MB | 73.8 MB |
| Physical footprint after 5,000 random charts | 258.9 MB | 58.4 MB |
| Median random chart | 1,189 µs | 1,947 µs |
| Median first call with OS pages warm | 74.7 ms | 69.3 ms |

The focused scalar comparison matched all 20,000 calls across base and
medium bit for bit. A four-second sample during the `auto/medium` workload
recorded 1,677 samples in `fstat`. The revision checks the descriptor before
and after each cache miss. Cache hits use the previously admitted immutable
words, so a cached result can retain the original value after an external
in-place write. New reads refuse ordinary changes detected by the descriptor
size and modification time. This remains dependent on replace-by-rename and
does not protect against a writer that restores the observed metadata.
The pinned descriptor survives atomic replacement of its path.

Run the paired benchmark in `tools/bench/benchmark.py` with and without
`--positional-centers` on the same installed wheel, using
`--mode auto --tier medium --workload random_chart --workload memory
--workload cold_start --workers 1 --calls 5000 --memory-calls 5000`.
For a scalar-only replay, use `bitdiff.py run --profile auto-centers
--candidate-positional-centers` with the same wheel as both inputs. Raw reports
and the sampling trace are in `target/bench/spk-positional-qualification/`.
The first prototype also matched the published wheel on all 2,100,000
corpus calls on macOS with the 32 MiB LEB cache enabled. The revised
wheel (`dd21b595`, SHA-256
`54628e7c9c84b6cbb069b441905a3ad5bf29ab560b697adb0ad90688c845c5d6`)
again matched 20,000 focused calls. All 155 SPK tests and 1,486 Python tests
passed. The revised wheel also matched all 2,100,000 macOS calls with the
32 MiB LEB pool, zero differences. Section 14 records the composed
hosted comparison.

| Metric, same revised wheel and setup | Owned | Positional, checks around uncached reads |
| --- | ---: | ---: |
| RSS after 5,000 random charts | 274.3 MB | 73.5 MB |
| Physical footprint after 5,000 random charts | 259.1 MB | 58.1 MB |
| Median random chart | 1,193 µs | 1,040 µs |
| p95 random chart | 1,247 µs | 1,093 µs |
| Median first call with OS pages warm | 75.2 ms | 68.5 ms |

### Multiprocess memory at 100,000 charts per worker

Mac arm64, medium tier, the revised wheel and pins as above, legacy LEB
retention. Each worker ran the same 100,000 random dates independently.

| Processes | Owned sum of RSS | Positional sum of RSS | Median physical footprint per process |
| --- | ---: | ---: | ---: |
| 1 | 284.5 MB | 84.1 MB | 269.3 → 65.0 MB |
| 2 | 569.1 MB | 168.5 MB | 269.3 → 65.1 MB |
| 4 | 1,137.0 MB | 336.3 MB | 269.4 → 65.0 MB |
| 8 | 2,269.4 MB | 668.1 MB | 268.1 → 66.0 MB |

These are per-process RSS totals, which count shared mappings once per
process. macOS did not expose PSS for these runs. The footprint column is a
per-process median. Raw data and vmmap outputs are the two
`*-memory-workers` reports in `target/bench/spk-positional-qualification/`.

### Existing indexed LEB reader

The separate base and medium `leb` trial used the same revised wheel and data
for both storage policies, five fresh processes and 2,000 charts per
timing repeat. The 100,000-chart memory trial above supplies the RSS figures.

| Tier | RSS after 100,000 charts | Median random chart | Median consecutive chart |
| --- | ---: | ---: | ---: |
| base | 57.0 → 51.5 MB | 556.5 → 812.0 µs | 110.7 → 114.1 µs |
| medium | 94.3 → 52.4 MB | 556.6 → 816.2 µs | 116.4 → 121.0 µs |

The indexed policy reduces retained heap, with a roughly 46% cost on this
random-date workload. The medium memory saving warrants checking a larger
coefficient cache under the positional reader's live allocation limits.
No release default has been selected.

## 18. Step 1b coefficient windows for indexed LEB

A four-second profile of `leb/medium` put full-frame
coefficient reconstruction at the top of the indexed reader profile
(1,400 samples). `pread` and `fstat` had 86 and 34 samples. The change
validates and decompresses the full frame but reconstructs and retains
only the selected window of at most 64 KiB. Positional admission still
uses the existing frame, compressed-range, window and working limits.
The retained working lease covers the actual selected word allocation.

All 243 LEB tests passed, including cross-window bit equality and releasing
the indexed retention and working leases on cooling. The two indexed wheels
matched all 10,000 smoke calls. The standard installed-wheel Python suite
passed 1,486 tests. With indexed mode explicit, both Python wheels
also passed all 1,486 tests.

Initial macOS arm64 screening, matching indexed policies, same pinned data,
five fresh processes with 1,000 charts each. Auto memory jobs were running
during the screening, so the timings are preliminary.

| Tier | Full decoded chunks, random | Selected windows, random | Full decoded chunks, consecutive | Selected windows, consecutive |
| --- | ---: | ---: | ---: | ---: |
| base | 842.4 µs | 613.1 µs | 117.2 µs | 117.7 µs |
| medium | 851.1 µs | 623.0 µs | 124.4 µs | 124.1 µs |

Use `benchmark.py --mode leb --leb-reader indexed` with the same pinned set.
Use `bitdiff.py run --leb-reader indexed` to set the same reader on both
sides. Raw profiles are under `target/bench/indexed-sampling/` and the
comparison and timing reports are under
`target/bench/indexed-windows-qualification/`. The quiet timing repeat,
100,000-chart memory confirmation and full corpus are still open.

A repeat over 3,000 charts per process confirmed the reduction for indexed
random charts, 845.7 → 609.5 µs in base and 848.8 → 621.8 µs in medium.
Consecutive medians were 116.2 → 116.8 µs in base and 122.0 → 125.6 µs
in medium. Other qualification was active on the host in this repeat.

At 100,000 random charts, one indexed worker reached 46.9 MB RSS and
23.9 MB footprint in base, and 46.7 MB RSS and 24.3 MB footprint in medium.
The prior indexed reader reached 51.5 and 52.4 MB RSS respectively in the
same workload. The [hosted run](https://github.com/g-battaglia/libephemeris/actions/runs/43873289950)
records the platform checks for this revision.

## 19. Step 1b indexed cache capacity

The indexed reader's 16-entry limit could retain at most 1 MiB of 64 KiB
coefficient windows, even with an 8 MiB byte budget. The follow-up raises
only the indexed entry limit to 256. The shared indexed retention ceiling
and working reservations are still charged by actual coefficient capacity.
The entry limit bounds map overhead independently of those reservations.
The 244 LEB tests passed. A targeted regression fills the entry limit and
then lowers the byte budget, checking both accounting paths.

Screening on macOS arm64, pinned 1900–2100 coverage, same revised wheel
(`ef47ed5e`) for both readers. Timing uses five 2,000-chart repeats. RSS
is measured separately after 5,000 charts. Long memory jobs also ran on
this host, so these are screening measurements.

| Tier | Owned → indexed random median | Owned → indexed consecutive median | Owned → indexed RSS |
| --- | ---: | ---: | ---: |
| base | 728.4 → 421.0 µs | 141.0 → 116.9 µs | 49.7 → 46.9 MB |
| medium | 744.9 → 419.0 µs | 151.9 → 122.9 µs | 87.1 → 46.5 MB |

The candidate matched the 10,000-call indexed smoke corpus. The same-policy
full comparison, the 100,000-chart comparisons and the quiet timing
repeat follow below. Reports are under `target/bench/indexed-capacity-qualification/`.

The same revised wheel completed 100,000 random charts with the following
process memory. This longer run supersedes the 5,000-chart memory snapshot
for selecting a budget.

| Tier | Owned → indexed RSS, default 8 MiB | Owned → indexed footprint |
| --- | ---: | ---: |
| base | 57.1 → 65.4 MB | 34.8 → 37.2 MB |
| medium | 94.8 → 70.4 MB | 75.5 → 43.8 MB |

The base tier requires a smaller budget to meet the RAM target. Use
`benchmark.py --leb-reader indexed --indexed-cache-mib 2` (or 4, 6, 8)
to vary the shared indexed ceiling on the same staged data. This writes a
temporary `leb_limits.cache_bytes` setting for the worker. The timing
and memory results for lower budgets follow below.

The full macOS indexed comparison with the published release passed all
2,100,000 calls with zero differences. The host had finished the other
benchmark runs before the following same-wheel timing sweep. Each timing
cell represents five fresh processes with 3,000 charts each.

| Reader bound | Base random median / p95 | Medium random median / p95 | Base RSS after 100,000 charts |
| --- | ---: | ---: | ---: |
| owned | 556.8 / 587.1 µs | 566.1 / 596.9 µs | 57.1 MB |
| indexed, 2 MiB | 581.8 / 619.4 µs | 591.9 / 629.7 µs | 45.3 MB |
| indexed, 4 MiB | 567.1 / 621.5 µs | 577.1 / 630.0 µs | 56.0 MB |
| indexed, 6 MiB | 492.7 / 619.9 µs | 502.2 / 632.0 µs | 48.6 MB |
| indexed, 8 MiB | 396.4 / 620.0 µs | 401.5 / 630.0 µs | 65.4 MB |

The 6 MiB base configuration reduced peak RSS from 57.1 to 48.6 MB and
steady physical footprint from 34.8 to 32.8 MB in the single-worker runs.
Consecutive chart medians rose from 110.6 to 115.0 µs in base and from
116.4 to 121.6 µs in medium with the 8 MiB budget. These are candidate
budgets for automatically discovered pinned data, with the random p95 and
consecutive latency regressions recorded. Section 14 records the later
comparisons with the complete companion groups.

Two follow-up experiments left the performance candidate unchanged. The
paired scalar Vondrak rate loop measured 105.8 → 106.1 µs per consecutive
chart on `f5530185`, with all 698 relevant Rust tests and 10,000 smoke calls
passing. Reconstructing coefficient words from pre-sliced byte planes
(`86f88a06`) measured 478.3 → 511.7 µs per indexed random chart at 6 MiB.
All 244 LEB tests and the 10,000-call indexed comparison passed. The
reports are in `target/bench/frame-rate-trig-qualification/` and
`target/bench/leb-word-unpack-qualification/`.

With both tiers and the released companion groups installed, the same
wheel completed 100,000 charts per process at an indexed ceiling of 6 MiB.
This table reports the sum of steady RSS in MB. The sum can count shared
resident pages more than once.

| Tier | Processes | Owned → indexed RSS |
| --- | ---: | ---: |
| base | 1 | 56.3 → 48.6 MB |
| base | 2 | 114.6 → 116.2 MB |
| base | 4 | 229.1 → 254.6 MB |
| base | 8 | 455.1 → 515.7 MB |
| medium | 1 | 94.7 → 61.8 MB |
| medium | 2 | 189.4 → 129.3 MB |
| medium | 4 | 378.2 → 257.6 MB |
| medium | 8 | 753.6 → 377.1 MB |

The eight-process repeats at 2 and 4 MiB cover the same mixed installation
and 100,000 random charts per process:

| Indexed ceiling | Base summed RSS | Medium summed RSS |
| --- | ---: | ---: |
| 2 MiB | 419.1 MB | 421.1 MB |
| 4 MiB | 505.3 MB | 509.4 MB |
| 6 MiB | 515.7 MB | 377.1 MB |

For base, 2 MiB lowered eight-process RSS by about 8% relative to owned.
Its separate quiet timing run had about 4.5% higher median chart latency.
The two additional reports are `mixed-indexed-{2,4}-8x100k.json`.

For eight medium workers the sum of each process's peak RSS was
753.6 → 513.4 MB. Peaks can occur at different times. The full provisioned
base layout does not support a default switch based on RSS. The indexed
reader and its byte ceiling remain explicit choices. Reports are
`mixed-{owned,indexed}-100k.json` in the indexed capacity directory.

The [paired Linux memory comparison](https://github.com/g-battaglia/libephemeris/actions/runs/43879751526)
completed on ARM64 and x86-64 using one identical candidate wheel per
architecture for both readers. Each run used all eight pinned LEB files
and 100,000 random charts per process. At eight processes:

| Linux architecture | Tier | Owned → indexed summed PSS | Owned → indexed anonymous RSS |
| --- | --- | ---: | ---: |
| ARM64 | base | 254.3 → 229.8 MB | 237.4 → 212.6 MB |
| ARM64 | medium | 560.4 → 229.1 MB | 543.8 → 212.0 MB |
| x86-64 | base | 257.4 → 232.2 MB | 238.8 → 213.5 MB |
| x86-64 | medium | 563.5 → 231.8 MB | 544.9 → 213.2 MB |

The indexed runs used 6 MiB. The artifacts contain per-process peaks,
file residency and each checkpoint for 1/2/4/8 workers. File-backed RSS
was similar across readers. At eight processes the medium PSS reduction
was about 59% on both architectures and the base reduction was about 10%.
These are memory results for the pinned release layout and sampled epochs.
The Mac base latency and RSS tradeoffs still apply to default selection.

A follow-up prototype for automatically discovered release files hashes the
opened descriptor before admitting positional reads (`2c798dd9`).
Comparing five fresh-process runs, this raised median first-call latency
from 1.85 to 5.34 ms for base and from 5.73 to 22.62 ms for medium.
In its separate one-process 100,000-chart comparison, steady RSS moved
55.4 → 64.4 MB in base and 96.6 → 53.8 MB in medium. The prototype's 10,000
base LEB calls matched the published wheel. The startup and base RSS
tradeoffs keep automatic acquisition out of this change. Raw reports are
in `target/bench/managed-leb-measurement/`.

Source inspection found existing exact-epoch memoization in
`nut00a_cached_centuries` and existing memoization for the frame matrices.
The native sidereal wrappers can call nutation twice, but a repeated epoch
uses the existing series memo. No additional cache is proposed for those
values.

## Step 1c reference obliquity

The native chart profile showed the two-word sine/cosine accounting for
327 of 692 first-body samples. The fixed Vondrak reference rotation is
now acquired once, using the same trigonometric function and arguments.
All 23 frame-velocity tests passed. The isolated wheel matched all 10,000
published-protocol smoke inputs against its parent.

Mac arm64, base/leb, legacy coefficient retention, five fresh processes
with 5,000 consecutive dates. The same pinned bytes were used in both
builds. The first-body median moved from 24.33 to 23.88 µs and complete
chart median from 110.92 to 110.42 µs. The median for the ten subsequent
bodies stayed at 71.71 µs. The raw results in
`target/bench/frame-constant-qualification/` qualify this small reuse.

The full macOS wheel comparison with the published 4.0.0 passed all
2,100,000 calls with zero differences. Capture files and the report are
in `target/bench/frame-constant-qualification/full/`.
