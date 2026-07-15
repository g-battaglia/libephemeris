# Licensing

LibEphemeris is licensed under the **Apache License, Version 2.0**
(`Apache-2.0`) — a permissive open-source license. See the [LICENSE](LICENSE)
file for the full text.

You may use, modify, and redistribute LibEphemeris — including in
closed-source and commercial products — subject to the terms of the Apache
License 2.0, most notably preservation of copyright, license, and attribution
notices and the [NOTICE](NOTICE.md) file.

Project-authored files in the current release tree are offered under
Apache-2.0; vendored files and external data retain the terms identified in
`THIRD_PARTY_NOTICES.md`. Published distributions retain the terms and content
under which they were released. PyPI distributions are immutable, so consult
the metadata and bundled notices for the exact version you use. The v3
remediation does not retroactively relicense removed material or earlier
artifacts.

## Copyright

Copyright (c) 2025-2026 Giacomo Battaglia.

Every project-owned Python source file in the shipped `libephemeris` package
carries the SPDX expression `Apache-2.0`.
Attribution and data-source details are in [NOTICE.md](NOTICE.md).

## Third-party components

A few vendored/adapted files keep their own upstream licenses — they are all
permissive (MIT) and Apache-2.0 compatible. The full inventory, license
texts, and the runtime/dev dependency licenses are in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

**Optional `nbody` extra:** the opt-in `libephemeris[nbody]` feature pulls in
`rebound` and `assist`, which are licensed **GPL-3.0-or-later**. They are not
part of the core install and are never bundled in any LibEphemeris artifact;
users and redistributors must comply with those packages' GPL terms. Whether a
particular installation or distribution is a combined work is fact-specific
and is not determined by this notice. The core library declares no
strong-copyleft (GPL/LGPL/AGPL) runtime dependency; `certifi` is a separately
installed MPL-2.0 dependency.
See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for details.

## Contributions

Contributions are accepted under the terms of the Apache License 2.0: unless
you explicitly state otherwise, any contribution you intentionally submit for
inclusion in the work shall be licensed as Apache-2.0, without any additional
terms or conditions (Apache-2.0, Section 5). No separate contributor license
agreement is required.

## Relationship to Swiss Ephemeris

The current LibEphemeris release line is maintained as an independently
implemented, signature-compatible API (see [NOTICE.md](NOTICE.md)). "Swiss
Ephemeris" is a product of Astrodienst AG; the name is used here nominatively
only. The reference API's Python binding may be used only for ephemeral
behavioral comparisons in a separate validation environment and is not a
dependency of this package. Earlier development history included material now
removed by the v3 remediation; Apache-2.0 applies only where the project has
the necessary rights.

---

*This document describes the project's licensing; it is not legal advice.*
