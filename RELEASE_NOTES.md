# LibEphemeris 4.0.0 — native Rust rewrite

**Stable release preparation:** Python `4.0.0`; Rust core/data `1.0.0`.
Source tag: `python-v4.0.0`. Earlier tags and published packages remain immutable.
Registry publication and platform availability must be checked separately from
these source notes. Preparing Python packages does not publish Rust crates.

## Promotion from Python 4.0.0b4

The stable preparation retains the b4 numerical implementation, coefficient
bytes and data-v4 identities. Changes are package identities, exact internal
version pins, stable metadata, release instructions and identity assertions;
no numerical tolerance or reference value is changed. The optional offline
LEB tools retain their existing independent version.

Python still requires 3.12+. The wheel matrix covers macOS ARM64, Linux
x86_64/aarch64 and Windows x86_64; other platforms need the source toolchain.
A stable version label does not remove the documented scientific limitations.

## Empty-path compatibility

An empty explicit path selects the current directory, matching the historic
Python `abspath("")` behavior. `None` restores default discovery; `"."` remains
an explicit current-directory selection. Normalization happens in the Rust
runtime, so Python and Rust callers cannot retain an invalid empty path.
Repeated empty/dot selections preserve readers. This supersedes the prepared
b3 interpretation without moving its source tag. Do not upload superseded b3
artifacts. Numerical formulas and coefficient bytes are unchanged.

## Changes since Python 4.0.0b2

Completed checksum ranges are released from the Linux page cache. Explicit
reader cooling preserves warm decoded coefficients, and owned images are
identified by shared allocation rather than full byte comparison. Consumers
can query an opaque reader-generation token to invalidate coverage caches after
provisioning or configuration changes. Historical regression inputs are bound
to their original data-v3 image; expected numbers and tolerances are unchanged.
Numerical models and data-v4 coefficients are unchanged. Reader selection remains
explicit: this release introduces no process-memory cap or concurrency limit.

## Changes since Python 4.0.0b1

The public IERS cache-age constant now has a concrete integer declaration in
its Python stub, with a regression against duplicate import bindings. UTC
live-clock test references retain integer ticks through the seconds conversion;
exact-rational checks prevent a false one-ULP interval failure without loosening
the bounds or changing the runtime clock. Release preparation retries unsettled
GitHub job metadata without accepting incomplete verification. The Python
package identity advances independently of Rust core/data. Numerical models,
coefficient bytes and the previously published releases are unchanged.

## A complete engine rewrite

LibEphemeris v4 replaces the Python calculation engine with an independently
authored **Rust engine**. It is not a wrapper around the old Python engine.
Python applications keep `import libephemeris`, using thin PyO3 bindings,
compatibility modules and shipped type stubs. Rust applications can use
`libephemeris-core` directly without Python.

The repository is reorganized around explicit responsibilities:

- `core/`: calculations, configuration, native data readers and CLI;
- `python/`: bindings, public imports, typing and Python integration tests;
- `core/data/`: immutable bundled coefficients, separate from calculation code;
- `tools/`: table generators, optional offline LEB tools and release checks.

There is one numerical engine rather than parallel Python and Rust engines.
The old implementation remains available in the v3 Git history and releases;
it is not installed alongside the new engine.

## Optimized native execution and resource management

The rewrite moves calculation kernels, JPL/SPK evaluation, LEB reading and
N-body integration into Rust. Immutable resources and decoded caches are shared
rather than recreated for each observing context. Indexed LEB reading can bound
decoded coefficient storage without loading an entire image into owned memory.
Native calls release the Python GIL, supporting bounded worker-based execution.

These changes provide a cleaner architecture and more efficient execution and
resource reuse. Actual gains depend on the workload, data route, cache state
and platform: no universal speedup or process-memory ceiling is claimed.
Calls remain synchronous, and cancelling an asynchronous wrapper does not
interrupt a calculation already running in Rust.

One local macOS ARM64 / CPython 3.12.13 check compared the published v3.2.2
with this beta using identical data-v4, sealed LEB mode, IERS off and a fixed
64-second Delta-T. After 200 warm-up calls, five samples of 1,000 Mars
`calc_ut` calls with `FLG_SPEED` at JD `2451545 + i/100` had median times of
0.1459 s (v3) and 0.02387 s (v4), approximately **6.1× faster for that workload**.
This is a small warmed calculation comparison, not a cold-start, peak-memory
or cross-platform benchmark.

## Independent of Swiss Ephemeris

The v4 calculation engine is **independent of Swiss Ephemeris**: it does not
require its engine, source code, FFI library or ephemeris data files.
Compatibility terminology describes interfaces, not implementation provenance,
affiliation or endorsement. This does not imply that earlier LibEphemeris
versions depended on or copied Swiss Ephemeris.

Skyfield and PyERFA are no longer Python runtime calculation backends.
REBOUND/ASSIST implementations are not bundled or linked; native Rust supplies
the supported integration capabilities. Historical names such as `skyfield`,
`rebound_integration` and `AssistEphemConfig` remain compatibility labels.
NASA/JPL inputs, ERFA-derived components, numerical tables and other included
third-party material retain their [documented provenance and terms](THIRD_PARTY_NOTICES.md).

## Python API compatibility

The public Python calculation surface is preserved: familiar names, argument
order, flags and ordinary result layouts remain available. Existing applications
continue to import `libephemeris`; using familiar APIs does not require Swiss
Ephemeris to be installed.

Compatibility does **not** promise identical dependency internals, exception
graphs, arbitrary Python numeric-object behavior or bit-for-bit results across
independent numerical implementations. Native records, resource lifecycle,
configuration capture and independently justified model corrections have
explicit migration notes. See [MIGRATION.md](MIGRATION.md) before upgrading.

The existing detailed compatibility baseline is v3.2.1. The v3.2.2 update
changed the data inventory without changing the public API; that fact is not
a claim that a new exhaustive v3.2.2 behavioral audit has been completed.
Precision and physical accuracy remain distinct; consult the
[limitations](docs/limitations.md) for model and data domains.

## The published data-v4 inventory

The beta uses the same immutable [data-v4 release](https://github.com/g-battaglia/libephemeris/releases/tag/data-v4)
introduced by [LibEphemeris v3.2.2](https://github.com/g-battaglia/libephemeris/releases/tag/v3.2.2).
No ephemerides were regenerated for this transfer.

The bundled `base_core.leb2` is 10,232,271 bytes, SHA-256
`b891de02f82b68346046da67389fd3713aa81ab299287c632cb20cc396397964`.
Its stored core interval is **JD 2396757.5–2506332.5**, one additional day at
each end compared with data-v3. Per-body intervals and the surrounding data
needed by apparent-place or event calculations remain authoritative.

All twelve LEB2 asset identities follow the unchanged
[upstream manifest](core/data/leb2/data-v4-manifest.json). Only the base core is
bundled; other groups and tiers remain explicit external inputs. Planet-center
SPK assets retain their separate data-v2rc1 identities. Library version v4 and
the data-v4 release tag are separate versioning schemes.

## Installing v4

Python requires **3.12+**. Use an available wheel for your OS/architecture, or
[build a local wheel](docs/installation.md#build-from-source) with the Rust and
C toolchains. Do not assume that wheels exist for every platform or that a
source preparation has already been published to PyPI or crates.io.

Once the corresponding registry artifacts are published, select Python 4.0.0
explicitly to avoid installing an older stable release. For offline
operation, begin with the bundled base core, choose `leb` mode and seal the
network. Follow the [installation guide](docs/installation.md) and
[migration guide](MIGRATION.md) for configuration and first calculations.

Public license declarations are unchanged. For inquiries about private
licensing, use the [project contact](README.md#license); requests are evaluated
case by case without a promise of an alternative license.
