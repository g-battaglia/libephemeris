# Performance qualification

Run from a source checkout. The timing driver and the wheel comparison use the
Python standard library. The native timing example uses the existing Rust
workspace dependencies. Artifacts belong under `target/bench/` locally.

## Isolated environments

The following example uses CPython 3.12 and uv. On Windows use
`Scripts/python.exe` in each environment. Keep the selected data directory
stable throughout a run. The drivers stage only named, hash-checked source
files in a temporary data directory, use sealed network policy and disable IERS
by default, with automatic SPK downloads disabled. Use the release inputs rather than locally
regenerated files with matching names.

```sh
uv venv --python 3.12 target/bench/baseline
uv venv --python 3.12 target/bench/candidate
uvx --from pip pip download --only-binary=:all: --no-deps \
  --dest target/bench/published libephemeris==4.0.0
uvx --from maturin==1.9.4 maturin build --release --locked \
  -i target/bench/candidate/bin/python --out target/bench/wheels
uv pip install --python target/bench/baseline/bin/python \
  target/bench/published/*.whl numpy==2.4.6
uv pip install --python target/bench/candidate/bin/python \
  target/bench/wheels/*.whl numpy==2.4.6 pytest
```

`DATA_DIR` below is the directory with pinned inputs, normally
`$HOME/.libephemeris`. Smoke needs only `leb/base_core.leb2`. Full needs the
base and medium LEB2 groups, `de440s.bsp`, `de440.bsp` and the base and medium
planet-centre kernels. The wheel workflow stages those inputs when manually
dispatched with `bitdiff_profile=full`. Use the installed wheel for measurement
and a dedicated profiling build for sampling.

For sealed extended LEB, use `bitdiff.py run --profile leb-extended-indexed`
with the same interpreter, wheel, data and output options as the full gate.
This profile adds 350,000 scalar calls and requires the released LEB2 groups
for base, medium and extended to exercise the cumulative inventory. It
selects indexed readers and needs no SPK kernel.
Pass `--iers-data /path/to/captured/iers_cache` to stage and hash the three
observed-clock files from a deployment snapshot. Both processes use the same
snapshot, and downloads stay disabled. The wheel workflow can run the
profile with `bitdiff_profile=leb-extended-indexed` and IERS disabled.
For timing and memory, use `--mode leb --tier extended --leb-reader indexed
--companions` with the same optional IERS input.

## Time and memory

```sh
DATA_DIR="$HOME/.libephemeris"
target/bench/baseline/bin/python tools/bench/benchmark.py \
  --data "$DATA_DIR" --mode leb --tier base --window 1900-2100 \
  --wheel target/bench/published/*.whl --commit python-v4.0.0 \
  --output target/bench/baseline-time.json
target/bench/candidate/bin/python tools/bench/benchmark.py \
  --data "$DATA_DIR" --mode leb --tier base --window 1900-2100 \
  --wheel target/bench/wheels/*.whl --commit "$(git rev-parse HEAD)" \
  --build-options 'maturin build --release --locked' \
  --output target/bench/candidate-time.json
cargo build --release --locked -p libephemeris-core \
  --features bundled-data,nbody,python-internals --example performance_bench
target/bench/candidate/bin/python tools/bench/benchmark.py \
  --data "$DATA_DIR" --mode leb --tier base --window 1900-2100 \
  --native target/release/examples/performance_bench --commit "$(git rev-parse HEAD)" \
  --build-options 'cargo build --release --locked --features bundled-data,nbody,python-internals' \
  --output target/bench/native-time.json
target/bench/candidate/bin/python tools/bench/benchmark.py \
  --data "$DATA_DIR" --mode leb --tier base --workload memory \
  --output target/bench/candidate-memory.json
```

Omit `--mode` and `--tier` for the base/medium and leb/auto/skyfield matrix.
The defaults time random charts, consecutive charts and chart components.
Repeat `--workload` to select other cases. `--window 1550-2650` and
`--window coverage` probe broader windows. Dates are restricted to the
intersection of the relevant raw coverage and an interior margin. The report
records the date interval and a Python provider trace. The single-position
workloads trace the measured Mars call with its actual flags. Chart workloads
trace all chart bodies. A failed calculation or unexpected provider fails
the run. Use `--companions` to include published companion groups. Compare
reports using the same options and inputs.

Each timing repeat is a fresh process. The default is five repeats, with
median and nearest-rank p95 over pooled individual timings. Thread throughput
uses batch averages per chart. Memory collects checkpoints for 1/2/4/8
independent workers, each running 100,000 charts by default. Use `--workers`
or `--memory-calls` to select a smaller diagnostic run. Threading and process
memory are measured through the Python binding. Measure performance with
competing jobs stopped, including validation jobs, and retain raw samples.

For a cache sweep on the candidate, repeat the same workloads with
`--leb-cache-mib 16`, `32`, `64` and `128`. Each run verifies that the
worker activated the specified budget and includes the worker's retained
coefficient bytes at each memory checkpoint. For the numerical comparison,
`--candidate-cache-mib 64` enables that budget only on the candidate side.
Omit these flags for the legacy baseline. `--leb-reader indexed` selects
the independently budgeted positional LEB reader for comparison.

On Linux, use `bitdiff.py run --profile leb-medium --candidate-shared-leb`
with the same interpreter, wheel, data and output arguments as the other
comparisons. The comparison gives the candidate a private registry directory
and requires its first medium LEB calculation to map the sealed image. The
full profile also accepts this switch. The wheel workflow dispatch input is
`candidate_shared_leb` and applies on the two Linux runners. For memory
and timing, pass `--mode leb --tier medium --shared-leb-medium` to
`benchmark.py`. The memory report verifies that all measured workers hold the
same sealed inode.
The same dispatch checks owner exit, recreation, busy-lock fallback and bad
sealed bytes on each Linux wheel. To rerun the lifecycle check locally, use
the installed candidate interpreter with `tools/bench/shared_leb_check.py
--data "$DATA_DIR" --output target/bench/shared-leb-lifecycle.json`.

Use `--series-calls` to change the length of the consecutive single-body
time series, which defaults to 100,000 calls.

To reproduce the README's published-wheel chart cases, use the prepared
baseline interpreter and run the base-tier chart workloads in all three
modes. The report retains medians, p95s and raw timing samples.

```sh
target/bench/baseline/bin/python tools/bench/benchmark.py \
  --data "$DATA_DIR" --tier base --window 1900-2100 \
  --workload seq_chart --workload random_chart --flags 258 --calls 300 --repeats 5 \
  --wheel target/bench/published/*.whl --commit python-v4.0.0 \
  --output target/bench/readme-baseline.json
```

The cold-start workload records import time, first calculation and first
full chart in each fresh process. With `page_cache_warm`, the driver reads
all staged inputs immediately before each timed process. With `first_use`,
each process receives its own staged input files. The route probe and
previous repetitions use different files. Both conditions start
with empty process-local decoder and coefficient caches. A previously
decoded persistent store is outside these configurations. On Linux,
`page_cache_cold` syncs and evicts the private staged inputs before each
timed process. It calls `mincore` on private file mappings and requires
zero resident bytes immediately before spawn. Each sample records that
observation.
Other hosts report this condition as an unmet prerequisite. [Linux file
cache advice](https://man7.org/linux/man-pages/man2/posix_fadvise.2.html)
and [residency probe](https://man7.org/linux/man-pages/man2/mincore.2.html).

Each startup sample uses one date for one position call followed by one
chart. Memory samples perform the same two startup steps before the
separately counted random-chart workload.

Linux uses `/proc/PID/status` and `smaps_rollup`. Sum of PSS apportions shared
pages and is the available estimate of the group's attributable resident
memory. Sum of RSS double-counts sharing. macOS records `ps` RSS, the OS peak
RSS, and raw `vmmap -summary` output with physical footprint and allocator
resident bytes. Precise anonymous/file RSS and group PSS remain null on
macOS. Allocator resident bytes give a narrower view than all anonymous
memory. On Windows unavailable memory fields are null. Use a platform with
the required counters when qualifying a RAM claim.

For hosted checks, enable `linux_memory` in the manual wheel workflow
and choose a `linux_memory_case` such as `auto/medium`. The selected
mode and tier run against the published and rebuilt wheels with 1, 2, 4
and 8 workers and 100,000 charts per worker. Pinned companion inputs
are staged for LEB and Auto. Raw Linux reports are in the wheel artifacts.
Enable `linux_startup` to measure all three startup conditions on the
selected case.

## Bit identity

```sh
target/bench/candidate/bin/python tools/bench/bitdiff.py run --profile smoke \
  --data "$DATA_DIR" --output target/bench/smoke \
  --baseline target/bench/baseline/bin/python \
  --baseline-wheel target/bench/published/*.whl --baseline-commit python-v4.0.0 \
  --candidate target/bench/candidate/bin/python \
  --candidate-wheel target/bench/wheels/*.whl --candidate-commit "$(git rev-parse HEAD)" \
  --candidate-build-options 'maturin build --release --locked'
target/bench/candidate/bin/python tools/bench/bitdiff.py compare \
  target/bench/smoke/baseline.jsonl.gz target/bench/smoke/candidate.jsonl.gz
```

Smoke is 10,000 calls in base/leb. Full is 2,100,000 calls over six mode/tier
pairs. It combines the module and context position calls, boundary probes,
representative event and derived APIs, state changes, single public flag
bits, a sweep of 14 independent position bits and dense date samples. This
covers the enumerated masks, not their full Cartesian product with every
entry point, body, provider and date. Mode-dependent refusals are compared
as typed errors. Compare successful coverage per API when interpreting a
result. The dump preserves float words, including NaN payloads and signed
zero, and checks result structure, returned flags and Python warnings.
The sidereal-mode transition in smoke and full includes an unsupported
projection bit that emits a real warning. The capture requires that
warning and the comparison reports how many calls emitted warnings on
each wheel.
Capture files and the coverage/error summary can be retained as CI artifacts.
Exit status 1 means a difference, and 2 means no completed comparison.
Missing or incompatible required inputs are reported as
`missing_prerequisite` in `comparison.json`, including when the
wheel workflow cannot stage a pinned file or a capture rejects its wheel.

The report lists exact requested flag masks, target IDs grouped by body
class, and observed center IDs for each mode, tier and position entry
point. This inventory includes requests that returned errors. Re-run
`bitdiff.py compare` on earlier captures to generate it for those runs.
On a difference, the report counts all affected requests by mode, tier
and API, including changes between values and typed errors. The warning
count can overlap result changes. The report also retains the first ten
differences with their full requests and outcomes.

The manual `auto-edge` profile runs the first 30,000 calls of each
`auto` tier with pinned LEB companions and SPKs. It omits the LEB and
Skyfield runs. The wheel workflow also compares the published wheel
against itself in independent processes and retains the control captures
in `reports/bitdiff/published-control/`.

The manual `satellites` profile compares 320 calls for named Pluto
satellites and the Charon alias in Skyfield medium. It checks module
positions, planet-centric positions and the direct moon API at J2000
and around the segment edges. Local runs use `--profile satellites` and
expect the pinned `plu017.bsp` next to `de440.bsp` and the medium planet
centers file in the data directory. The manual workflow fetches and
verifies the fixture. [NAIF satellite archive](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/).

Set the manual workflow's `rust_toolchain` to `1.99.1` to build with the
toolchain identified in the published Windows 4.0.0 wheel.
Each wheel report records any Rust source revisions embedded in the
extension. The host's `rustc -Vv` is recorded independently.
For a separate full comparison against the pinned SPK correction, select
`bitdiff_profile=full` and `corrected_spk_reference=true`. The workflow
builds the corrected 4.0.0 reference and the candidate with the selected
toolchain on every runner and retains both wheels and complete captures.
Leave the option disabled for the strict comparison against the published
4.0.0 wheel. The PR smoke check uses that published wheel.
The full manual matrix allows 150 minutes per job. Smoke and focused
comparisons retain the 90-minute cap.

Run validation and downstream suites using the candidate wheel in isolated
consumer environments. Where a consumer lock selects a different build, use
an isolated import overlay containing only the candidate package and its
distribution metadata, and check the imported extension SHA-256 before the
suite. Use the consumer's existing tests and record omitted or unavailable
suites explicitly. The external scientific validation checkout can select the
candidate via that overlay and `LIBEPH_VAL_PYTHON`.

## Profiling

Build the Rust example with `--profile profiling`, which inherits release
settings and adds line tables. For the binding, build a separate wheel with
`maturin build --profile profiling --locked`. Use a fixed workload with enough
calls to sample the worker during calculation. For a long-lived worker
`PID`, use `sample PID 3` on macOS or `perf record -g -p PID -- sleep 3` on
Linux. Inspect the hot stacks, distinguishing initial acquisition from the
steady-state calculation. Keep timing reports from the release build.
