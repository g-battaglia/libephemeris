# Licensing

LibEphemeris is **dual-licensed**: every release of the same codebase is
offered under two alternative grants, and you choose the one you use it
under.

## 1. Open source — AGPL-3.0-only

The default license is the [GNU Affero General Public License v3](LICENSE)
(AGPL-3.0-only). It is free for any use — including commercial use — as
long as you comply with its terms, most notably: if you distribute the
software or make a modified version available to users over a network,
you must make the corresponding source code available under the AGPL.

The packages published on PyPI are AGPL-3.0-only.

## 2. Commercial license

Organizations that cannot or do not want to comply with the AGPL (for
example: embedding LibEphemeris in a closed-source product, or running a
SaaS without source disclosure) can obtain a commercial license from the
copyright holder instead.

Commercial-licensing terms are arranged directly with the copyright holder
on request — please get in touch:

- Contact: Giacomo Battaglia — <kerykeion.astrology@gmail.com>

## How the dual license works

- Giacomo Battaglia is the copyright holder of LibEphemeris (provenance:
  [NOTICE.md](NOTICE.md)).
- Each owned source file carries the SPDX expression
  `AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial`: the recipient may
  elect either grant.
- **Exceptions**: a few vendored/adapted files keep their upstream
  licenses (see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)); they are
  all permissive (MIT). There is no copyleft code in the package: the
  Galilean satellite module was rewritten clean-room in June 2026 from the
  published Lieske 1998 / Meeus theory and is now dual-licensed like the
  rest of the project (see
  [docs/methodology/galilean-clean-room-2026-06.md](docs/methodology/galilean-clean-room-2026-06.md)).
- Versions up to and including the 2.x series were published under
  AGPL-3.0-only; those grants remain valid for those published versions.

## Contributions

External contributions are accepted only with a contributor agreement that
licenses the contribution to the maintainer for distribution under **both**
grants (AGPL-3.0-only and the commercial license). Opening a pull request
constitutes agreement unless stated otherwise; a formal CLA flow will be
added before external contributions are merged.

## Relationship to Swiss Ephemeris

LibEphemeris is an **independent implementation** with an API that is
signature-compatible with pyswisseph; it does not contain Swiss Ephemeris
code and has no license relationship with Astrodienst AG (details and
provenance record: [NOTICE.md](NOTICE.md)). "Swiss Ephemeris" is a product
of Astrodienst AG; the name is used here nominatively only.
pyswisseph appears solely in the `dev` extra as a black-box test oracle
and is never a runtime dependency.

---

*This document describes the project's licensing model; it is not legal
advice. Commercial-licensing terms are arranged directly with the
copyright holder — contact <kerykeion.astrology@gmail.com>.*
