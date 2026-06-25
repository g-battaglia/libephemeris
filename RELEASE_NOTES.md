# Release Notes

The authoritative, machine-checkable change history lives in
[`CHANGELOG.md`](CHANGELOG.md). Long-form, narrative release notes for each
version live under [`release-notes/`](release-notes/). This file front-pages the
current release.

---

## 3.0.0rc1 — the v3.0.0 release candidate (2026-06-25)

**v3.0.0 is the dual-licensed (AGPL-3.0-only OR commercial), clean-room
provenance release** that re-grounds the whole library on long-term-valid models
(Vondrák 2011 precession & obliquity, long-term sidereal-time house cusps, a
multi-era Delta T with a selectable model) and ships a full review-driven
correctness sweep across eclipses, houses, fixed stars, minor bodies, and the
LEB / Horizons backends — while keeping the v2 canonical bare-name public API and
1:1 reference parity except for documented intentional divergences.

This is the first **release candidate** of the v3 line, published on the PyPI
pre-release channel and intended to be promoted to `3.0.0` final unchanged if it
proves clean.

```bash
pip install --pre libephemeris==3.0.0rc1
```

Highlights:

- **Dual licensing** — `AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial`, on
  a clean-room provenance footing with no copyleft code (see `LICENSING.md`).
- **Vondrák 2011 long-term precession & obliquity** across the whole pipeline
  (valid ±200,000 years); modern results unchanged to sub-milliarcsecond.
- **Multi-era Delta T** with a selectable model (`set_delta_t_model`).
- **Long-term house cusps** and **true time-derivative cusp speeds**
  (`houses_ex2`).
- **Eclipse / occultation / fixed-star / minor-body correctness sweep** and a
  functional live **Horizons** backend.
- A hardened **test & validation** regime (100% line-coverage campaign,
  oracle-free invariants, independent cross-checks).

**Upgrading from v2?** Several observable behavior changes are deliberate
(default sidereal ayanamsha, eclipse retflags/obscuration, error policy, fixed
stars, remote-epoch positions, house speeds). Re-check any pinned values against
the migration table in
[`release-notes/v3.0.0rc1.md`](release-notes/v3.0.0rc1.md),
[`docs/guides/migration-guide.md`](docs/guides/migration-guide.md), and
[`docs/comparison/intentional-divergences.md`](docs/comparison/intentional-divergences.md).

Full detail:
[release-notes/v3.0.0rc1.md](release-notes/v3.0.0rc1.md) ·
[CHANGELOG.md](CHANGELOG.md).

---

## Earlier releases

Per-version narrative notes are under [`release-notes/`](release-notes/):

- **2.0.0** (2026-05-14) — removal of the legacy `swe_`/`SE_`/`SEFLG_` prefixed
  aliases in favour of the canonical bare-name public API; see `CHANGELOG.md`.
- **1.x** — the LEB binary-ephemeris line (LEB1/LEB2 compressed format, page-cache
  management, performance work). See `release-notes/v1.*.md`.
- **0.x** — early development. See `release-notes/v0.*.md`.

For the complete, structured history of every release, see
[`CHANGELOG.md`](CHANGELOG.md).
