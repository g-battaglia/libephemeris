# Licensing

LibEphemeris is licensed under the **Apache License, Version 2.0**
(`Apache-2.0`) — a permissive open-source license. See the [LICENSE](LICENSE)
file for the full text.

You may use, modify, and redistribute LibEphemeris — including in
closed-source and commercial products — subject to the terms of the Apache
License 2.0, most notably preservation of copyright, license, and attribution
notices and the [NOTICE](NOTICE.md) file.

The packages published on PyPI are Apache-2.0.

## Copyright

Copyright (c) 2025-2026 Giacomo Battaglia.

Every owned source file carries the SPDX expression
`Apache-2.0`. The independent-implementation and provenance record is in
[NOTICE.md](NOTICE.md).

## Third-party components

A few vendored/adapted files keep their own upstream licenses — they are all
permissive (MIT) and Apache-2.0 compatible. The full inventory, license
texts, and the runtime/dev dependency licenses are in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

**Optional `nbody` extra:** the opt-in `libephemeris[nbody]` feature pulls in
`rebound` and `assist`, which are licensed **GPL-3.0-or-later**. They are not
part of the core install and are never bundled in any LibEphemeris artifact;
installing the extra forms a combined work that is, for that user, governed by
the GPL. The core library and all required runtime dependencies are permissive.
See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for details.

## Contributions

Contributions are accepted under the terms of the Apache License 2.0: unless
you explicitly state otherwise, any contribution you intentionally submit for
inclusion in the work shall be licensed as Apache-2.0, without any additional
terms or conditions (Apache-2.0, Section 5). No separate contributor license
agreement is required.

## Relationship to Swiss Ephemeris

LibEphemeris is an **independent implementation** with an API that is
signature-compatible with pyswisseph; independence from Swiss Ephemeris
code is the project's working standard, enforced by provenance sweeps and
remediated whenever a finding surfaces (details, disclosures and the
current remediation record: [NOTICE.md](NOTICE.md) and
docs/methodology/independence-remediation-2026-07.md). While that record
lists open items (repository history, output-calibrated data sets), this
document does not assert a categorical absence of obligations toward
Astrodienst AG. "Swiss Ephemeris" is a product
of Astrodienst AG; the name is used here nominatively only.
pyswisseph is used only as a black-box test oracle in the separate
validation tooling; it is not declared among this package's dependencies
(not even in the `dev` extra) and is never a runtime dependency.

---

*This document describes the project's licensing; it is not legal advice.*
