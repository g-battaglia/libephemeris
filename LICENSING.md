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

Contributions are accepted under the terms of the AGPL-3.0-only: unless
you explicitly state otherwise, any contribution you intentionally submit for
inclusion in the work shall be licensed as AGPL-3.0-only, without any
additional terms or conditions. No separate contributor license agreement
is required.

## Relationship to Swiss Ephemeris

LibEphemeris provides an API that is signature-compatible with the reference
ephemeris API (see [NOTICE.md](NOTICE.md)). "Swiss Ephemeris" is a product
of Astrodienst AG; the name is used here nominatively only. The reference
API's Python binding may be used for ephemeral compatibility comparisons in
a separate validation environment, and is never a dependency of this
package.

---

*This document describes the project's licensing; it is not legal advice.*
