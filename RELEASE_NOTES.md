# Release Notes

The authoritative, machine-checkable change history lives in
[`CHANGELOG.md`](CHANGELOG.md). Long-form notes for each release live under
[`release-notes/`](release-notes/). This page points to the candidate currently
prepared in the source tree.

---

## 3.0.0rc8 — prepared review candidate (2026-07-13)

`3.0.0rc8` is the candidate currently prepared in this source tree; it has not
yet been published to PyPI. It is a clean-room review and release-hardening pass
over rc7, with focused fixes for context entry points, forced backend selection,
nakshatra boundaries, pre-UTC Julian dates, packaging, and public API
documentation.

After publication it will be installable with
`pip install --pre libephemeris==3.0.0rc8`. Until then, install this source tree
to evaluate rc8.

Highlights:

- `EphemerisContext.calc()` and `calc_ut()` now preserve backend, flag, and
  tracing semantics across LEB, Skyfield, Horizons, SPK, ASSIST, and Keplerian
  paths.
- Explicit Skyfield or Horizons mode bypasses a cached LEB reader.
- `split_deg()` handles nakshatra boundary-adjacent values consistently.
- Pre-1972 Julian-calendar input is mapped to the correct civil instant.
- The documented top-level export surface now includes all intentionally public
  performance helpers, eclipse extensions, moon identifiers, and constants.
- Clean-room and packaging checks cover ignored/untracked worktree artifacts,
  data namespaces, optional dependencies, wheels, and source distributions.
- API documentation is generated from the live public exports and validated by
  strict Sphinx and focused contract tests.

For the complete rc8 fixes and upgrade details, see
[`release-notes/v3.0.0rc8.md`](release-notes/v3.0.0rc8.md),
[`CHANGELOG.md`](CHANGELOG.md), and the
[`migration guide`](docs/guides/migration-guide.md).

---

## Earlier releases

Per-version narrative notes are under [`release-notes/`](release-notes/):

- **3.0.0rc7 and earlier v3 prereleases** — correctness, licensing,
  provenance, precision, and compatibility history.
- **2.x** — canonical bare-name API and removal of legacy prefixed aliases.
- **1.x** — the LEB binary-ephemeris and performance line.
- **0.x** — early development.

For the complete structured history, see [`CHANGELOG.md`](CHANGELOG.md).
