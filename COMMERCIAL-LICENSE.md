# LibEphemeris Commercial License — DRAFT

> **STATUS: DRAFT TEMPLATE — NOT A LICENSE OFFER.**
> This outline is an engineering placeholder for the commercial grant
> referenced by the SPDX expression `LicenseRef-LibEphemeris-Commercial`.
> It must be completed, reviewed and approved by legal counsel before any
> commercial license is sold. Until then, the only effective public grant
> is AGPL-3.0-only (see [LICENSE](LICENSE) and [LICENSING.md](LICENSING.md)).

- **Licensor:** Giacomo Battaglia
- **Contact:** <giacomo@libephemeris.dev>
- **Identifier:** `LicenseRef-LibEphemeris-Commercial`

## Outline to be completed with counsel

1. **Grant** — non-exclusive, non-transferable right to use, modify and
   redistribute LibEphemeris in object or source form as part of the
   licensee's products, without AGPL source-disclosure obligations.
2. **Scope tiers** — per-product / per-organization; SaaS vs distributed
   software; redistribution to end users.
3. **Third-party components** — the grant covers only code owned by the
   Licensor. Vendored components keep their own licenses (see
   [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)); in particular
   `libephemeris/moon_theories/galilean.py` remains LGPL-3.0 and its
   notices/source must accompany distributions until replaced.
4. **Trademarks** — no rights to the "Swiss Ephemeris" name (Astrodienst
   AG) or to the LibEphemeris name beyond accurate attribution.
5. **Warranty disclaimer and limitation of liability.**
6. **Support and updates** — what the fee includes (versions, duration).
7. **Term and termination** — effect of termination on shipped products.
8. **Governing law and venue.**

## Counsel review checklist (commercial edition)

- Final license/EULA text for the items above.
- Independence/provenance dossier: [NOTICE.md](NOTICE.md),
  `docs/methodology/provenance-sweep-2026-06.md`, verification logs under
  `tasks/results/`.
- Calibration-data position (tables fitted to pyswisseph output as a
  black-box oracle — see "Calibration Data Disclosure" in NOTICE.md).
- API-compatibility position (*Google LLC v. Oracle America*, 593 U.S. 1
  (2021); EU Directive 2009/24/EC, interoperability provisions).
- ESA Hipparcos data attribution terms; EU database right considerations
  for compiled catalogs.
- Contributor rights audit (sole-author status; AI-assisted output
  ownership statement in NOTICE.md).
