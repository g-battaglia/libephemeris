# Third-Party Notices

This file inventories third-party software and data that LibEphemeris
**contains** (vendored or adapted code shipped inside the package) or
**depends on** (PyPI dependencies installed alongside it), together with
their licenses. It is maintained as part of the Apache-2.0 compliance
process — see [LICENSING.md](LICENSING.md) and [NOTICE.md](NOTICE.md).

## Code shipped inside the package (vendored / adapted)

| Component | Path | Upstream | License |
|---|---|---|---|
| spktype21 0.1.0 | `libephemeris/vendor/spktype21.py` | Shushi Uetsuki (whiskie14142) | MIT |
| TASS 1.7 port | `libephemeris/moon_theories/tass17.py`, `libephemeris/moon_theories/tass17_data.py` (generated periodic-term tables) | Johannes Gajdosik's Stellarium implementation of Vienne & Duriez TASS 1.7 | MIT |

Each of these files carries its own `SPDX-License-Identifier` header and
documents its upstream provenance in the module docstring. They are **not**
covered by the project's Apache-2.0 license, but they are all permissively
(MIT) licensed; the package contains no copyleft code.

The Galilean satellite module (`libephemeris/moon_theories/galilean.py`)
was, through v2.1.0, adapted from PyMeeus and licensed LGPL-3.0. It was
**rewritten clean-room in June 2026** from the published theory (Lieske,
J.H. 1998, "Galilean satellite ephemerides E5", A&AS 129, 205; Meeus 1998,
*Astronomical Algorithms*, ch. 44) and is now owned by the project under
the Apache-2.0 license. See
[docs/methodology/galilean-clean-room-2026-06.md](docs/methodology/galilean-clean-room-2026-06.md)
for the process and independence record.

### MIT license text (spktype21, TASS 1.7 port)

> Permission is hereby granted, free of charge, to any person obtaining a
> copy of this software and associated documentation files (the
> "Software"), to deal in the Software without restriction, including
> without limitation the rights to use, copy, modify, merge, publish,
> distribute, sublicense, and/or sell copies of the Software, and to
> permit persons to whom the Software is furnished to do so, subject to
> the following conditions:
>
> The above copyright notice and this permission notice shall be included
> in all copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
> OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
> MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
> IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
> CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
> TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
> SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Runtime dependencies (installed from PyPI, not redistributed here)

| Package | License | Notes |
|---|---|---|
| skyfield | MIT | Ephemeris access, time scales |
| skyfield-data | MIT | Bundled astronomy data files |
| jplephem | MIT | Transitive via skyfield; SPK kernel reader |
| sgp4 | MIT | Transitive via skyfield |
| numpy | BSD-3-Clause | Transitive |
| astroquery | BSD-3-Clause | JPL Horizons / SPK downloads |
| pyerfa | BSD-3-Clause | Bundles the ERFA C library, distributed under the ERFA license (permissive, SOFA-derived, with naming restrictions) |
| click | BSD-3-Clause | CLI framework |
| certifi | MPL-2.0 | Used unmodified; MPL file-level copyleft imposes no obligation on LibEphemeris code |
| zstandard | BSD-3-Clause | LEB2 compression |

Dependencies are installed from PyPI by the user or installer under their
own licenses; they are not part of LibEphemeris's license grant.

## Optional `nbody` extra — GPL-3.0 (not bundled)

The opt-in `libephemeris[nbody]` (and `libephemeris[all]`) extra installs two
**GPL-3.0-or-later** packages used by the shipped `rebound_integration.py`
module for ephemeris-quality n-body propagation of minor bodies:

| Package | License | Notes |
|---|---|---|
| rebound | GPL-3.0-or-later | N-body integrator; imported lazily, only when the extra is installed |
| assist | GPL-3.0-or-later | Ephemeris-quality REBOUND extension (JPL small-body integrator) |

These packages are **not** part of the core install, are **never bundled or
redistributed** in any LibEphemeris artifact (they are installed separately
from PyPI under their own GPL license), and are imported only on demand. The
LibEphemeris source itself is Apache-2.0. Because Apache-2.0 is one-way
compatible with the GPL, a user who chooses to install `libephemeris[nbody]`
forms a combined runtime work that is, **for that user**, governed by the
GPLv3. The default/core library and all of its required runtime dependencies
remain fully permissive, and the shipped wheel contains no copyleft code.

## Data sources

- **JPL DE440 / DE441** planetary and lunar ephemerides (NASA/JPL/Caltech):
  U.S. government work, freely usable; downloaded at runtime, not bundled.
- **ESA Hipparcos Catalogue** (ESA SP-1200, 1997) and **van Leeuwen 2007**
  reduction: star positions and proper motions (facts of nature) compiled
  into `libephemeris/fixed_stars.py`, with attribution in the module.
- **IAU Working Group on Star Names** (WGSN) designations.
- **IERS** Earth-orientation / Delta-T / leap-second tables: public data,
  fetched at runtime.
- `libephemeris/data/fictitious_orbits.csv`: this project's own compilation
  of hypothetical-body orbital elements. Most rows cite primary historical
  publications (Witte & Sieggrün 1928, Neely 1988, Strubell 1952, Hoyt
  1980, Weston, peer-reviewed papers); a few rows without a known
  publication are disclosed interoperability values recovered by black-box
  fits against reference-API output (NOTICE.md, Calibration Data
  Disclosure).
- `libephemeris/data/leb2/base_core.leb2`: generated by this project from
  JPL DE440 (Chebyshev coefficient tables).

## Development-only tools (not shipped, never linked)

- **pyswisseph** (AGPL-3.0 / Swiss Ephemeris dual license): used
  exclusively as a black-box comparison oracle in the separate validation
  tooling. It is not declared among this package's dependencies (not even in
  the `dev` extra), is never imported by the shipped package, and is not a
  runtime dependency. See [NOTICE.md](NOTICE.md) for the independence
  statement.
- **Swiss Ephemeris reference data files** (`sefstars.txt`, `seorbel.txt`,
  `*.se1`): used only as local oracle inputs for comparison tests; they
  are gitignored (`data/reference/`), downloaded on demand by
  `scripts/download_reference_data.py`, and never bundled in any artifact.
