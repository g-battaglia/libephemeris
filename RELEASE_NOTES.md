# Release Notes

The authoritative, machine-checkable change history lives in
[`CHANGELOG.md`](CHANGELOG.md). Long-form notes for each release live under
[`release-notes/`](release-notes/). This page points to the release currently
prepared in the source tree.

---

## 3.1.0 — current source-tree release (2026-08-10)

`3.1.0` retires the `uranians` LEB companion (fictitious bodies are always
served by their runtime analytical models) and stops sealed `leb` mode from
warning on the by-design fictitious-body dispatch. Details in
[`CHANGELOG.md`](CHANGELOG.md) and
[`release-notes/v3.1.0.md`](release-notes/v3.1.0.md). This page only
records which version the tree currently carries; it does not duplicate
the changelog.

After publication it is installable with `pip install libephemeris==3.1.0`.

For the complete, structured change history across the `3.0.0rc*` series and
upgrade details, see [`CHANGELOG.md`](CHANGELOG.md), the per-release notes under
[`release-notes/`](release-notes/), and the
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
