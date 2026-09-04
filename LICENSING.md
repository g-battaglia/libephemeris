# Licensing

LibEphemeris is licensed under the **GNU Affero General Public License v3.0**
(`AGPL-3.0-only`). See the [LICENSE](LICENSE) file for the full text.

You may use, modify, and redistribute LibEphemeris under the terms of the
AGPL v3. If you modify the program and make it available to users
interacting with it remotely over a network, you must make the complete
corresponding source available to those users under the AGPL.

## License history

Versions through `v2.0.2` and v3 pre-releases through `v3.0.0rc1` were
released under AGPL-3.0-only.

Pre-releases `v3.0.0rc2` through `v3.0.0rc8` incorrectly declared
`Apache-2.0` in their PyPI metadata. An internal audit determined that
the relicensing was premature: unresolved components prevent an Apache-2.0
grant. Those releases are to be treated as **AGPL-3.0-only** and have been
**yanked** from PyPI.

Starting with `v3.0.0rc10`, the project metadata correctly declares
`AGPL-3.0-only`.

## Copyright

Copyright (c) 2025-2026 Giacomo Battaglia.

Every project-owned Python source file in the shipped `libephemeris` package
carries the SPDX expression `AGPL-3.0-only`.
Attribution and data-source details are in [NOTICE.md](NOTICE.md).

## Third-party components

A few vendored/adapted files keep their own upstream licenses (all MIT).
The full inventory, license texts, and the runtime/dev dependency licenses
are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

**Optional `nbody` extra:** the opt-in `libephemeris[nbody]` feature pulls in
`rebound` and `assist`, which are licensed **GPL-3.0-or-later**. GPL-3.0 and
AGPL-3.0 are compatible per their respective Section 13 provisions; the
combined work is governed by the AGPL-3.0. These packages are not part of
the core install and are never bundled in any LibEphemeris artifact;
installing the extra forms a combined work that is, for that user, governed
by the AGPL. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for
details.

## Contributions

Contributions are accepted only under the Contributor License Agreement in
[CONTRIBUTING.md](CONTRIBUTING.md). By submitting a pull request or any other
contribution you irrevocably assign its copyright to the maintainer, Giacomo
Battaglia (with an exclusive licence as a fallback where such an assignment is
not effective under local law), and you keep a perpetual licence to reuse your
own work.

Centralizing copyright is what allows the maintainer to re-license the project
— in whole or in part, including your contribution — under any other terms,
open-source or proprietary, at their sole discretion. The project itself stays
publicly available under AGPL-3.0-only, and published AGPL releases remain
AGPL: those grants are irrevocable. Your authorship is credited in the Git
history, in [AUTHORS](AUTHORS), and, where appropriate, in release notes.

On GitHub, a CLA bot checks each pull request against the recorded list of
contributors who have agreed (`.clabot`) and comments with instructions when an
author is not yet on it.

## Relationship to Swiss Ephemeris

LibEphemeris provides an API that is signature-compatible with the reference
ephemeris API (see [NOTICE.md](NOTICE.md)). "Swiss Ephemeris" is a product
of Astrodienst AG; the name is used here nominatively only. The reference
API's Python binding may be used for ephemeral compatibility comparisons in
a separate validation environment, and is never a dependency of this
package.

---

*This document describes the project's licensing; it is not legal advice.*
