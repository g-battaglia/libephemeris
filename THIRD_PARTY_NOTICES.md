# Third-Party Notices

Project-owned Rust and binding code is licensed AGPL-3.0; the complete
text is in [COPYING](COPYING). Components below retain their upstream terms.
The Python package includes base LEB coefficients and generated numerical
catalog/clock tables. Engine `src/` paths below are relative to `core/`;
`tools/` paths are workspace-relative. The optional `core/data/` package exports
the unchanged coefficient bytes without calculation code. Full DE/SPK and ASSIST datasets remain external inputs.

Scientific definitions, upstream software and external datasets are distinct
attribution and licensing subjects. Citations do not grant rights to upstream
software or datasets. External NASA/JPL, catalog and observation data retain
their providers' terms. Conditions on included third-party material apply to
that material; they do not determine ownership of independently authored
project code.

## Code and numerical data included in the project

| Component | Path | Origin | Terms |
| --- | --- | --- | --- |
| TASS 1.7 periodic terms | `src/tables/tass17_tables.rs`, `tools/generators/data/tass17.json` | Johannes Gajdosik's Stellarium implementation of the Vienne–Duriez theory | MIT, retained below |
| Earth-orientation/precession routines and coefficient tables | `src/time/erfa.rs`, `src/tables/nutation_tables.rs` | ERFA 2.0.1 / IAU SOFA-derived material | BSD-style ERFA conditions, retained below |
| Nutation and historical/IERS Delta-T numerical tables | `src/tables/iau2000_tables.rs`, `src/tables/smh2016_tables.rs`, `tools/generators/data/iau2000_nutation.npz` | Skyfield bundled data | MIT, retained below |
| Resolved Cargo dependencies | `licenses/` | Packages recorded in `Cargo.lock` | Individual original texts and notices included |

The TASS tables implement the published Vienne & Duriez (1995, A&A 297,
588–605) theory with the later IMCCE TASS 1.7 revision. They are generated
from the attributed numerical snapshot in `tools/generators/data/tass17.json`
by `tools/generators/gen_tass17_tables.py`. That snapshot retains the exact
historical v3.2.1 source identity; regeneration requires no v3 checkout or
execution. Their original MIT attribution is retained. Fixed-star catalog measurements and nomenclature retain their
Hipparcos, IAU and WGSN provenance in the generated table headers.

Retained upstream copyright notice:

- Stellarium TASS 1.7 implementation: Copyright (c) 2005 Johannes
  Gajdosik.

### MIT license text (TASS 1.7 tables)

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

## ERFA long-term precession

The native Earth-orientation routines and tables include material derived
from ERFA `ltpecl`, `ltpequ`, `ltp` and `ltpb` v2.0.1 (peeled commit
`9915ba38c9365f8b0738269b8c2ac1fdd5f8dee3`), whose source is at
<https://github.com/liberfa/erfa/tree/9915ba38c9365f8b0738269b8c2ac1fdd5f8dee3/src>.
The published Vondrák coefficients and IAU series are retained in native
numerical evaluation. Rust validation and arithmetic are adaptations; ERFA
is not linked at runtime.
These Earth-orientation components are derived from ERFA/SOFA; they are not
SOFA itself or an endorsed standards implementation. This attribution does not
extend to independently authored project components, which retain the crate
license.

### Retained ERFA copyright, conditions and disclaimer

Copyright (C) 2013-2023, NumFOCUS Foundation. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of the Standards Of Fundamental Astronomy Board, the
   International Astronomical Union nor the names of its contributors may
   be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

## Historical N-body references

Earlier project revisions bundled REBOUND 4.6.0 and ASSIST 1.2.3 under
GPL-3.0-or-later, with their complete license texts and copyright notices.
REBOUND credits Hanno Rein, Shangfei Liu and the contributors named in its
source files. ASSIST credits Matthew Holman, Arya Akmal, Hanno Rein and the
contributors named in its source files. The historical bundled revisions were
REBOUND `e853676f202355b7c3e62f96501fb2c2d3c8c862` and ASSIST
`9493cc36d07e5cc13107fe713c67e1139111edd5`. Their notices remain applicable to
the earlier distributions that contain that code.

The current N-body implementation is project-owned Rust under AGPL-3.0.
The REBOUND/ASSIST implementations, FFI adapter and native build requirements
have been removed from the current sources and distributions. Public API names
such as `rebound_integration`, `AssistEphemConfig` and the five integrator
selectors remain for compatibility; their version queries identify the Rust
implementations. The `nbody` feature enables these capabilities, and
`native-nbody` remains its compatibility alias. Selected classic-DE and
asteroid SPK datasets remain external inputs. This change does not alter the
project license or the terms of other included dependencies.

## Skyfield numerical tables

Only numerical data are retained; Skyfield is not a runtime computation backend.
The original MIT text below is copied from the installed pinned reference's
Skyfield distribution license. Table input identities are recorded in their
generated Rust headers.
The bundled `tools/generators/data/iau2000_nutation.npz` has the same identity
as the Skyfield 1.54 dataset recorded in the Rust nutation-table header. The
available Skyfield 1.55 distribution supplies byte-identical numerical data;
including it changes neither coefficients nor the runtime model. Offline
regeneration reads this local file without importing Skyfield or NumPy.

### Retained Skyfield MIT license

Copyright © 2013–2018 Brandon Rhodes and available under the MIT license:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Cargo dependency licenses and notices

`licenses/index.json` inventories all non-workspace
packages resolved by the current Cargo lockfile, including optional and
target-specific dependencies. Each entry records the package/version, upstream
source, declared license expression and paths to its actual copied license,
copyright, authorship and notice files. Original filenames and text bytes are
preserved; the index does not choose between alternative licenses.

This includes PyO3 0.29.3, pyo3-log 0.13.4, rust-numpy 0.29.0, serde_json,
transport/TLS dependencies and the
bundled Zstandard sources. Their terms remain those in the corresponding
upstream files. A package's presence in the all-features lockfile does not mean
that every target or artifact links it.

Regenerate the bundle with `python3 tools/release/package_licenses.py` after a dependency
change, then run `python3 tools/release/package_licenses.py --check`. The check reads
Cargo metadata offline and verifies exact upstream text copies for the current
resolution. `pyproject.toml` includes these files, this notice and the project
AGPL text in its distribution license-file selection.
Final wheel/sdist content checks must also inspect the resulting artifacts.

The isolated native allocation regression uses an original project-owned
system-allocator wrapper under AGPL-3.0. It requires no third-party
allocation-instrumentation package and is not linked into the runtime library
or Python extension.

## Project-owned SPK reader

The native DAF/SPK reader (`src/providers/spk/`) is project-owned code written
against the public NAIF DAF/SPK Required Reading. No third-party SPK reader
implementation is copied into the crate.
