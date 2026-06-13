# LibEphemeris Commercial License

LibEphemeris is dual-licensed. The default public grant is **AGPL-3.0-only**
(see [LICENSE](LICENSE) and [LICENSING.md](LICENSING.md)). Organizations that
cannot or do not want to comply with the AGPL — for example, embedding
LibEphemeris in a closed-source product, or running a SaaS without source
disclosure — can instead obtain a commercial license, identified by the SPDX
expression `LicenseRef-LibEphemeris-Commercial`.

## How to obtain a commercial license

Commercial-licensing terms are arranged **directly with the copyright
holder**. There is no online purchase or self-serve form — please get in
touch and we will work out the terms for your use case:

- **Licensor:** Giacomo Battaglia
- **Contact:** <kerykeion.astrology@gmail.com>

A commercial license typically covers the right to use, modify and
redistribute LibEphemeris in object or source form as part of your products
without the AGPL's source-disclosure obligations, scoped to your
distribution model (per product / per organization, distributed software vs.
SaaS).

## What the grant covers

The grant covers only code owned by the Licensor. The vendored components
listed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) keep their own
licenses; they are all permissive (MIT) and their notices must accompany
distributions. **The package contains no copyleft code.** No rights are
granted to the "Swiss Ephemeris" name (Astrodienst AG) or to the
LibEphemeris name beyond accurate attribution.

## Independence & provenance

The independence/provenance record that backs the commercial offering:

- [NOTICE.md](NOTICE.md) — independent-implementation statement, calibration
  disclosure, AI-authorship note.
- `docs/methodology/provenance-sweep-2026-06.md` — Swiss Ephemeris
  fingerprint remediation; `docs/methodology/galilean-clean-room-2026-06.md`
  — the Galilean clean-room rewrite that retired the last copyleft module.
- `docs/methodology/hypothetical-bodies.md` — published-source provenance
  for the hypothetical-body elements.
- Verification baselines under `tasks/results/`.
