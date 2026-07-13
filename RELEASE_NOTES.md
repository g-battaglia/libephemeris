# Release Notes

The authoritative, machine-checkable change history lives in
[`CHANGELOG.md`](CHANGELOG.md). Long-form notes for each release live under
[`release-notes/`](release-notes/). This page points to the candidate currently
prepared in the source tree.

---

## 3.0.0rc8 — prepared review candidate (2026-07-13)

`3.0.0rc8` is the candidate currently prepared in this source tree; it has not
yet been published to PyPI. It restores the complete rc7 compatibility surface
while retaining the provenance, package-content, and hash-pin hardening from
the review. It also fixes context/backend edge cases, LEB tier installation,
and outer-planet center coverage.

After publication it will be installable with
`pip install --pre libephemeris==3.0.0rc8`. Until then, install this source tree
to evaluate rc8.

Highlights:

- All predefined ayanamsha modes 0–46 calculate without a J2000 warning or
  degradation; direct, Horizons, and LEB paths share the same definitions.
- Historical hypothetical bodies 40–58 and their rc7 module-level classes and
  constants are available again in every calculation mode. Old LEB1 channels
  fall through to the runtime models.
- `INTP_APOG` and `INTP_PERG` retain the rc7 Delaunay-series behavior through
  an immutable, SHA-256-pinned compatibility refinement.
- Default medium installations discover a pinned bundled LEB fast path;
  base/medium/extended downloads work through both LEB entry points, and local
  legacy `.leb` files remain discoverable.
- The regenerated medium planet-center asset retains all ten JPL descriptors,
  extends Saturn/Uranus/Neptune through 1550–2650, and removes speed spikes at
  center-kernel coverage edges.
- `EphemerisContext.calc()` and `calc_ut()` preserve backend, flag, and tracing
  semantics across LEB, Skyfield, Horizons, SPK, ASSIST, and Keplerian paths.
- Horizons range checks for analytical lunar points no longer trigger a local
  DE440 download; `get_tid_acc()` starts at the reference-compatible `-25.80`.
- The documented export surface, generated API reference, packaging audits,
  provenance scans, and source-status tables now match the runtime behavior.
- Compatibility items are individually classified with per-mode source/audit
  status. All implementations are derived from published standards.

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
