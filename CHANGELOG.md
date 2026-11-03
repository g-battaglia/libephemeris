# Changelog

## Unreleased — cross-platform fixes

- Discover the native user profile when `HOME` is absent, including Windows
  `USERPROFILE`, while preserving an explicit nonempty `HOME` override.
- Runtime coverage queries no longer label an out-of-range date with a serving
  file; explicit reader inspection still identifies its input image.
- Keep the Unix-only optional fingerprint capture out of Windows test builds,
  and use portable timeout cleanup for the threaded minor-body subprocess test.
- Check Windows Rust test-target compilation in CI without executing Rust tests.


## 4.0.0b1 (Python) · 1.0.0-beta.1 (Rust core) — Unreleased

See the [beta release notes](RELEASE_NOTES.md) for the full rewrite and migration overview.

### Native engine

- Rust owns calculations, configuration, data acquisition and the CLI.
  Python contains native bindings, re-exports and type stubs.
- Local LEB1/LEB2, native JPL/SPK and explicit Horizons calculation routes retain
  tiered data selection and actual per-body coverage checks.
- Planetary/lunar/stellar and minor-body calculations, houses, sidereal frames,
  events, photometry and visibility use the native engine.
- Independent Rust Solar and prescribed-ephemeris N-body integration preserves
  historical compatibility API names without REBOUND/ASSIST implementation
  dependencies. Force datasets remain external inputs.
- Canonical sidereal time separates UT1 rotation from the continuous TT mean
  celestial origin. Primary ayanamsha getters retain the public observation
  clock, dynamic-anchor, flag and error semantics.

### Interfaces and resources

- `libephemeris_core::api` provides concrete Sessions, typed errors, closed enums
  and complete named result records. Native data/provider integration remains
  in `context`, `config` and `native`.
- Python retains the pinned v3.2.1 public calculation surface, with documented
  native input, error and handle differences. Numerical precision, not cross-
  engine binary identity, is required; see [migration](MIGRATION.md).
- Shared immutable resources, bounded coefficient caches and optional indexed
  LEB access support retained contexts. Global close/reset semantics remain
  application lifecycle concerns.
- Native work releases the Python GIL but remains synchronous; bounded whole-job
  offload is documented for asynchronous services.

### Tools and documentation

- The bundled base and all twelve LEB2 download descriptors use the published
  data-v4 inventory from v3.2.2. Base stored coverage is JD
  2396757.5–2506332.5, extending the previous core interval by one day per side.

- One native CLI supports configuration, inspection, explicit provisioning,
  lossless conversion and maintenance products.
- Optional standalone LEB fitting/generation tools use native sources and do not
  add Python astronomy to the runtime package.
- Scientific comparisons and reference environments live in a separate
  qualification repository; product builds and tests do not require it.
- Documentation is organized into guides, reference, data, sources and
  development; release tools and generators have separate directories.
- Clarified the distinction between project-owned code, third-party material and
  historical N-body implementations; corrected the ERFA three-clause attribution.
  Private licensing inquiries are evaluated case by case; public license terms
  remain unchanged.
- Project license metadata consistently uses `AGPL-3.0`; upstream grants remain
  unchanged. Rust and Python distribution versions are independent.
- The Python package carries its README and upstream license texts. The
  binding Cargo package is not intended for separate crates.io publication.
- Stack-specific source, tests and examples live under `core/` and `python/`;
  common project metadata, sources, documentation and tools remain at the root.
- An optional data-only package preserves static, zero-copy embedding of the
  unchanged base image while keeping both Rust archives within registry limits.
- Published-manifest/golden-value regressions use their exact identified data;
  current-data tests verify real digest-based trust instead of trusting filenames.
- The beginner guide explains sky coordinates, clocks and runnable first examples;
  documentation navigation separates learning, API use, models and development.

### Release boundaries

The [known limitations](docs/limitations.md) remain applicable.
Model domains, optional data coverage and unqualified platforms are not waived
by preparation work or a future beta label.
