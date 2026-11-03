# Offline LEB tools

These tools use the installed LibEphemeris v4 native extension and local data.
They do not require LibEphemeris v3, Skyfield or a separate Python ephemeris
implementation. From the v4 checkout, install the native package and generation
dependencies with:

```sh
uv sync --group generate
```

The native package and the optional `libephemeris-leb-tools` distribution can
also be installed from their respective local wheels. The tools depend on
LibEphemeris 4.0.0b1, NumPy and zstandard; the runtime does not depend on the
tools. Installed commands and spawned workers use the installed modules and
need no checkout or v3 installation. Build the separate tools artifacts with:

```sh
uv build --out-dir target/leb-tools tools/leb_generation
```

For private licensing inquiries, see the [project contact](https://github.com/g-battaglia/libephemeris/blob/v4.0.0b1/README.md#license);
requests are evaluated case by case.

## Installed generation command

Both the native executable and Python console support:

```sh
libephemeris leb generate medium full --data-dir /path/to/local/data \
  --auto-sources --start 2000 --end 2001 --output-dir /path/to/output
libephemeris leb generate medium body SUN,MOON,EARTH --kernel /path/to/de440.bsp \
  --start-jd 2451545 --end-jd 2451911 -o /path/to/three-bodies.leb
```

The native executable finds `libephemeris-generate` on PATH. The Python
console selects that executable beside its own interpreter, including when
invoked by an absolute path without activating the environment. The helper is
a separate process with literal arguments and inherited standard streams;
its exit status is preserved. There is no installation or checkout fallback.
Sources/output options must be explicit: a child cannot inherit an embedding
Rust Runtime's in-memory configuration.

| Mode | Result |
| --- | --- |
| `full` | Every eligible body in one verified file. |
| `single` | The same inventory, fitting each body in a fresh sequential process. |
| `groups` | Planets, asteroids, exotics and analytical parts, then verified merge; retain all parts. |
| `body NAME[,NAME...]` | Exact names/IDs, case-insensitive; no silently omitted channel. |
| `planets`, `asteroids`, `exotics`, `analytical` | One complete named part. |
| `merge` | Verify and merge the four exact named parts without refitting. |
| `all` in place of tier/mode | Base, medium and extended full generation, stopping on the first failure. |

Without `--output`, files go to `--output-dir` or `DATA_DIR/leb` as
`ephemeris_TIER.leb`; group parts append `_GROUP`. All destinations are checked
before a run. Each publication is atomic; a later failure may leave earlier
verified parts or tiers. `--force` permits replacement after verification.
Merge checks exact inventories, global domains and shared auxiliary bytes.
Grouped knot-aligned bodies retain their actual domains beneath one common
global header; no body is extrapolated. Auxiliary sections are stored in the
planets part and retained in the merged file.

`all` requires `--data-dir`, optionally `--output-dir`, and rejects a single
`--kernel` or `--output`. Each tier selects its own existing main kernel.
The default is one fitting worker; `--workers N` is explicit, and `single`
requires one worker. `--group`, `--bodies`, `--single` and `--tier` cannot
override positional selections. Body-only outputs may need observer/core
channels from other parts before complete apparent-coordinate calculations.

The same local-source, source-domain and verification rules below apply.
Missing requested channels are errors. Calendar domains use exact January 1
endpoints, without v3's historical one-day margins. Complete nominal tier
regeneration still requires sufficient local source/force coverage; the mode
name is not a claim of astronomical qualification over that whole domain.

The direct installed interfaces remain available as
`libephemeris-generate-leb` (the options of `tools/generate_leb.py`) and
`libephemeris-generate-leb2` (the options of `tools/generate_leb2.py`).

### Regenerating one group

Keep the four part files produced by `groups`. To refresh one group, generate
that named part with the same tier, date range, source selections and clock
configuration, then rebuild the main file from all four parts:

```sh
libephemeris leb generate medium planets --kernel /path/to/de440.bsp \
  --start 2000 --end 2001 --output-dir /path/to/output --force
libephemeris leb generate medium merge --output-dir /path/to/output --force
```

Only the selected part is refitted. Merge streams the other parts' original
coefficient words and body domains, checks the complete named inventories and
shared auxiliary data, then verifies the combined container through Rust.
Each part and the main file publish atomically; a failed merge leaves the
previous main file intact and retains the newly verified part for inspection.

The main file and a replacement part contain overlapping body IDs, so they are
not disjoint merge inputs. The planets part owns nutation, Delta-T and the star
snapshot in the grouped workflow; preserve that clock selection when maintaining
the parts. There is no implicit winner for duplicate bodies or conflicting
auxiliary sections in a disjoint merge.

To update an existing main without the other three part files, select one named
group and add `--update-main MAIN --force`:

```sh
libephemeris leb generate medium planets --kernel /path/to/de440.bsp \
  --start 2000 --end 2001 --output-dir /path/to/parts \
  --update-main /path/to/ephemeris_medium.leb --force
```

This fits and verifies the selected part, then replaces exactly its complete
body inventory in the existing LEB1. It retains every unselected body's domain
and coefficient words, adds any selected body absent from the main, and sorts
the final body index by ID. The main and replacement must have identical global
date ranges, including the common knot-aligned start. An incomplete selected
group is an error. The group output must differ from the main and source files.

All auxiliary sections come from the original main, including opaque sections.
Clock options and auxiliary generation for the new part do not change the
main's nutation, Delta-T or catalogue. Use `groups` or `full` when those shared
tables also need regeneration. `--update-main` accepts one of `planets`,
`asteroids`, `exotics` or `analytical`; it requires an existing LEB1 and explicit
`--force`. Convert the updated LEB1 separately if LEB2 products are needed.

The combined file passes native verification in a sibling temporary directory
before atomic replacement. The main and partial hashes are checked again before
publication; inputs must stay immutable during the operation. A failed update
preserves the main and retains a successfully verified new part for inspection.

## Native lossless conversion

Convert an existing LEB1 with the native executable or installed Python console
entry. This command needs no offline Python generation dependencies:

```sh
libephemeris leb2 convert base --data-dir /path/to/data --output-dir /path/to/output --json
libephemeris leb2 convert medium-core --input /path/to/input.leb --output /path/to/core.leb2
```

Without explicit paths, the source is `DATA_DIR/leb/ephemeris_TIER.leb` and
outputs are `DATA_DIR/leb2/TIER_GROUP.leb2`. Whole-tier conversion creates the
complete groups present in the input and refuses partial groups or unknown
bodies before writing. Core owns the auxiliary sections; if core is absent,
the first whole-tier product retains them. An explicit companion-group target
selects only its body channels. The extended exotic inventory is checked
against its declared subset.

Every coefficient bit is retained, including signed zero and subnormal words.
The default decoded payload ceiling is 1 MiB per chunk; use
`--max-chunk-bytes N` to change it. Segments are never split. Reorder, shuffle,
compression and verification need additional chunk buffers, decoder workspace
and bounded index metadata; the ceiling is not a process RSS limit. Auxiliary
sections are copied and compared in blocks of at most 64 KiB. The input inode
must remain immutable throughout conversion.

The native reader decompresses every output chunk and compares all original
coefficient words and retained auxiliary bytes before publication. Existing
files are refused unless `--force` is supplied. Publication is atomic for each
file; a failure in a later group can leave earlier verified files published.
Native conversion always uses lossless Zstandard level 3. The offline Python
converter below has additional fitting/compression options and may produce
different compressed sizes from the same coefficient words.

## Offline Python conversion

Convert a local LEB1 file into a chunked LEB2 file:

```sh
uv run --group generate libephemeris-generate-leb2 convert input.leb -o output.leb2
```

Conversion preserves every coefficient bit by default. Reordering, byte
shuffling and Zstandard compression are reversible. The native indexed Rust
reader checks the complete output, including all compressed chunks, before an
atomic publication. Failed conversion preserves any existing output; use
`--force` to replace it after successful verification.

`--group core|asteroids|exotics|apogee` selects an exact named inventory;
`--tier base|medium|extended` authenticates its tier, including the extended
exotic subset. Shared nutation, Delta-T and stellar sections are retained in
unfiltered/core files; named companions contain their body channels.
Coefficients and auxiliary sections are streamed, including from large LEB1
files. No data file is downloaded or modified.

Explicit `--lossy` enables mantissa truncation. `--precision` then limits the
sum of coefficient errors in each stored component's units; it cannot relax a
body's default budget. That bound does not establish apparent-coordinate or
rate precision. The base-file diagnostic has 60 failing public-call pairs under
the fixed rate criteria, despite small position differences, so lossless
conversion is the default. See the diagnostic
and the lossless verification.

Stored structure and source-precision checks are separate. The standalone
console tools verify generated output with the native reader; cross-engine and
apparent-coordinate comparisons are maintained in the separate qualification
repository.

Generate a core LEB1 from an exact existing local JPL kernel:

```sh
uv run --group generate libephemeris-generate-leb --tier medium \
  --kernel /path/to/de440.bsp --start 2000 --end 2001 \
  --group core -o target/generated/core-2000.leb
```

Generate LEB2 directly from the same native sources, or verify an existing file:

```sh
uv run --group generate libephemeris-generate-leb2 generate --tier medium \
  --kernel /path/to/de440.bsp --start 2000 --end 2001 \
  --group core -o target/generated/medium_core.leb2
uv run --group generate libephemeris-generate-leb2 verify target/generated/medium_core.leb2 \
  --kernel /path/to/de440.bsp --tier medium
uv run --group generate libephemeris-generate-leb2 verify target/generated/medium_core.leb2 --structural-only
```

Direct generation fits a private temporary LEB1, streams its conversion and
checks the final compressed reader against the selected native source before
atomic publication. It accepts the same explicit SPK/N-body sources, domain,
body inventory and auxiliary choices as LEB1 generation. Lossless remains the
default. Explicit lossy generation must also pass the final source criteria;
coefficient budgets alone cannot admit it. A failed final check preserves an
existing output and removes this operation's temporary files.

Both source generators accept `--workers N` to fit independent bodies in
parallel processes. Each body opens its own exact native sources, checks the
input hashes and exits after returning its temporary coefficient-file metadata.
`--single` uses one fresh process per body sequentially, containing the lifetime
of native allocations; combine it only with `--workers 1`. Without either
option, fitting remains in the current process. Worker counts are explicit;
the default is one to avoid multiplying memory use by the machine's CPU count.
All modes preserve requested body order and use the same fitting/verification
criteria. Auxiliary sections are generated once, and the final assembled reader
is checked against the parent's retained source before atomic publication.
A worker refusal or changed input preserves an existing destination. These
options apply to source generation, not merge or conversion.

`--auto-sources` assembles missing minor-body inputs from existing local SPKs.
Select the main explicitly with `--kernel`, or use `--data-dir` and the tier's
standard main filename. SPK search defaults to `DATA_DIR/spk`; repeat `--spk-dir`
to select other directories. Directories are not searched recursively. Native
descriptors identify actual targets, the two supported asteroid NAIF offsets,
J2000 frames and complete center-chain coverage. Names and file sizes never
establish target identity or coverage. Every selected image is hash-checked
and printed; `--spk BODY PATH NAIF` always overrides discovery for that body.
The DAF/SPK header identifies a candidate; non-SPK cache files are reported and
left unused. Unreadable files or invalid metadata in an identified SPK are
errors. ASSIST force binaries are separate inputs, even when their filenames
end in `.bsp`.

For ordinary SPK generation, selection requires one continuous requested
domain (or the explicit `--available-domains` overlap). One candidate must
contain all other eligible continuous source windows. Multiple copies of the
same bytes/target count once; distinct solutions with equal or incomparable
coverage require an explicit `--spk`. Width is a declared selection policy,
not an accuracy ranking. No failed numerical read triggers a different image,
download or Keplerian substitute.

On `extended`, automatic assembly routes the declared regular exotic inventory
through SPK-seeded N-body integration, except the six declared asteroid
perturbers that retain their SPK route. Distinct eligible J2000 seed solutions
require an explicit `--spk`; seed width does not choose an initial condition.
The default existing force files are `DATA_DIR/assist/linux_m13000p17000.441`
and `DATA_DIR/assist/sb441-n16.bsp`. Other tiers use the existing `.440` force
file when `--nbody` explicitly requests integration. Explicit `--assist-planets`
and `--assist-asteroids` override their individual selections. Missing forces
refuse generation, and actual guarded common coverage still limits each body.
Neither a wider main nor a wider planet force file supplies absent asteroid
forces. All paths are native, including force loading and integration.

`--group all` requests every supported core, asteroid, exotic and apogee channel.
Its extended-tier inventory retains the declared regular exotics; an explicit
`--bodies` list is never silently shortened. The default group remains `core`.
Automatic source verification also uses every actual stored body domain:

```sh
uv run --group generate libephemeris-generate-leb2 generate --tier medium \
  --data-dir /path/to/local/data --auto-sources --available-domains \
  --group all --workers 2 -o target/generated/medium.leb2
uv run --group generate libephemeris-generate-leb2 verify target/generated/medium.leb2 \
  --tier medium --data-dir /path/to/local/data --auto-sources
```

Both commands accept the same explicit source overrides. Verification must use
the same physical selections as generation; a different cached orbital solution
does not establish a fitting error. Nothing is downloaded by these commands.

By default every requested body must cover the complete requested interval.
Explicit `--available-domains` instead stores each body's complete continuous
overlap with that interval. This allows a full core and a shorter minor-body
channel in the same file. The global header retains the requested interval;
each body index records its actual closed endpoints, also printed by the CLI.
The native reader refuses dates outside that body's domain. This option never
omits a requested body, bridges a source gap, selects only the largest island
or extrapolates a source. A missing overlap or multiple disjoint overlaps
refuse generation and preserve an existing destination. N-body coverage already
includes its force-stage guard. Serial and isolated fits use identical domains.

For example, when a selected Chiron image begins after the requested core:

```sh
uv run --group generate libephemeris-generate-leb2 generate --tier medium \
  --kernel /path/to/de440.bsp --start 1850 --end 1950 --bodies SUN CHIRON \
  --spk 15 /path/to/2060_190001_210001.bsp 20002060 --available-domains \
  --workers 2 -o target/generated/partial-chiron.leb2
```

That two-body example checks raw fitting/storage only; complete apparent
minor-body calls also require their generated observer core. The public source
checker samples each actual body domain with its unchanged 30-day margin and
refuses channels too short for that declared comparison.

Actual process tests preserve serial coefficient values for Cartesian and
interpolated-apse channels, explicit SPK and SPK-seeded N-body sources, and both
supported Python interpreters. A 2000 core generated directly as LEB2 with
`--single` passes 4,005 off-grid position / 2,763 Cartesian-rate checks and
all 1,344 complete aligned-clock native-source comparisons. These cover the
declared grid; full-tier and auxiliary qualification remain separate. See the
process-generation check.

`verify` checks every compressed chunk through actual native evaluation, then
checks raw positions and Cartesian rates against the explicit local sources.
`--structural-only` omits the source comparison and says so in its result. It
does not accept source image selections or establish astronomical precision.
Prefetching pages with `warm` alone does not decompress or validate a payload.
Neither verification mode changes a tolerance; auxiliary-model precision and
complete apparent-coordinate precision require their separate comparisons.

Every LEB channel is fitted on a TT grid. Rust converts Cartesian samples to
split TDB epochs before evaluating their geometric ICRS sources, and scales
analytic velocities by d(TDB)/d(TT). Angular lunar models already use TT.
Source TDB coverage is mapped inward to representable TT endpoints; a rounded
endpoint never claims an out-of-source clock instant. Raw TDB state APIs remain
available separately. Python performs only offline polynomial fitting and
serialization. Fitting uses actual representable sample epochs,
with fixed conservative intervals/degrees, continuous longitude unwrapping,
bounded batches and temporary coefficient files. Source coverage includes the
entire target/center chain; gaps, invalid frames and unavailable data are
refusals. No value is clamped or extrapolated and no file is downloaded.

Interpolated apsides have native piecewise-linear residual tables. Their
polynomial grid aligns with every residual and outer taper knot, using the
exact grid metadata returned by Rust. This avoids fitting across a derivative
break. An arbitrary requested start can require an earlier body boundary; the
file header describes the actual union and body indices retain their explicit
domains. The source samples and tolerances are unchanged.

Nine independent points per segment check positions; Cartesian fits also check
analytic native velocities. Before publication the native LEB reader checks 64
samples per body, including endpoints. `--verify-samples` changes the sampling
count, never a tolerance. `--start-jd`/`--end-jd` select exact dates instead of
calendar years; `--bodies` selects explicit public IDs or constant names.
`--skip-aux` omits shared nutation, Delta-T and catalogue snapshots. Otherwise
these use the native IAU 2006/2000A model, the selected native Delta-T model
(default SMH2016, IERS off) and the native catalogue. The established LEB star
snapshot has zero parallax/radial velocity; runtime star calculation uses the
complete native catalogue.

### Explicit local clock inputs

Generation and source verification accept `--delta-t-model smh2016` or
`--delta-t-model espenak_meeus`. To use existing IERS observations, supply both
`--iers-delta-t` and `--iers-cache-dir DIR`. The directory is explicit; no file
is downloaded. Native IERS selection uses the available local rows and retains
the selected model outside their coverage. An empty directory therefore does
not establish observational coverage.

```sh
libephemeris leb generate base full --kernel /path/to/de440s.bsp \
  --start 2000 --end 2001 --iers-delta-t --iers-cache-dir /path/to/iers_cache \
  --output-dir /path/to/output
```

`--delta-t-userdef-days VALUE` selects a finite constant **TT minus UT1 in
days**. It takes priority over IERS and the model. Omit it for automatic
selection; the native `DELTAT_AUTOMATIC` sentinel is not an explicit value.
Delta-T auxiliary entries also store days. `--skip-aux` omits that section;
it does not reset the source clock configuration.

Each fitting process reopens the same selected inputs. The tool records hashes
for `finals2000A.data`, `leap_seconds.dat` and `deltat.data`, including their
absence, and rejects changed, removed or newly added inputs before publication.
Merge preserves existing auxiliary bytes and rejects clock-selection options.
Structural-only verification also rejects those options because it does not
compare astronomical sources.

Nutation fitting also checks rotation rates. A degree-32 interpolant on a
denser source grid checks the derivative of the stored degree-16 polynomial;
nine off-grid points check actual native values and relative vector-rate
convergence at 1e-6. If a segment fails, the uniform section interval is halved
and fitting restarts. Undefined samples or failure after eight refinements
refuse generation. The denser derivative is a convergence estimate, not a
rigorous error bound; complete public comparisons remain necessary.

The complete medium 2000 generation includes all **53 supported bodies** and
passes 50,751 off-grid positions / 22,329 Cartesian rates, 3,392 / 3,008 final
native-reader source samples and all 5,088 complete TT/UT position/rate pairs.
Its former value-only nutation fit had four small speed failures: 23 intervals
of 15.913 days became 46 intervals of 7.957 days, retaining degree 16 and every
source/model/PLAN criterion. The earlier failures
and corrected complete comparison
remain separate. Automatic assembly also generates the medium core/five-asteroid
selection (1,536 complete pairs) and a real DE441 core/Pholus selection with
local DE441/N16 forces (1,152 complete pairs), all passing their declared grids.
The latter is a short SPK-orbit agreement check, not full deep-time N-body
accuracy. Full tier/auxiliary coverage remains pending.

Merge disjoint LEB1 body inventories, split named LEB2 groups, or recheck a LEB1:

```sh
uv run --group generate libephemeris-generate-leb --merge planets.leb analytical.leb -o merged.leb
uv run --group generate libephemeris-generate-leb2 convert-all merged.leb -o target/generated --tier-name medium
uv run --group generate libephemeris-generate-leb --verify-only merged.leb --kernel /path/to/de440.bsp --tier medium
```

Merge requires identical global date ranges and refuses duplicate bodies or
conflicting shared auxiliary data. `convert-all` produces standard
`medium_core.leb2` etc. for the complete groups present in the source; partial
inventories are errors. Every destination is verified before atomic publication.

The 2000 core example passes all raw fits and native-reader source checks.
Lossless LEB2 conversion exactly repeats all 1,344 complete TT/UT/frame pairs.
The direct source comparison also passes **1,344/1,344 pairs** when both sides
use the same physical TT epoch and shared planet-center corrections. Native
JPL's UT entry uses its independent bundled clock; LEB's UT entry uses the
configured model Delta-T. Pairing native JPL UT with LEB TT at that clock's
explicit epoch explains all 176 original differences. No position residual
selects a time, and no model or tolerance was changed. The
aligned source check
and original comparison
remain separate. Forced LEB's system-barycenter model adds further differences.
Source generation is **not yet fully qualified** across all date/body ranges;
angular raw rates and full-tier/auxiliary precision still need their full scope.

The current source implementation supports planets and angular lunar channels.
The 2000 apogee group passes 27,180 off-grid position checks and 192 native-read
source samples. All 3,600 complete position/rate pairs pass, including exact
native residual knots and their immediately adjacent representable epochs.
Lossless conversion exactly repeats its sampled raw and public outputs. See
the boundary check.
The former one-day grid crossed source knots and correctly refused publication;
the regression now covers an arbitrary start and both sides of a native-model
break. Wider-date and auxiliary qualification still remain.

Minor bodies can use explicit existing SPK images. Supply the public body ID,
exact path and the NAIF target actually stored in that file; repeat `--spk` for
each body. For example, the local ASSIST asteroid kernel uses target 2000001
for Ceres:

```sh
uv run --group generate libephemeris-generate-leb --tier medium \
  --kernel /path/to/de440.bsp --start 2000 --end 2001 --bodies CERES \
  --spk 17 /path/to/sb441-n16.bsp 2000001 -o target/generated/ceres-2000.leb
```

The native `GeometricSpkSource` retains the exact primary and supporting main
images, composes SSB states in Rust, and intersects all target/center coverage
windows. A wider alternate image never replaces the primary target. Source
images are hashed and reported; a missing center, a gap or an uncovered range
refuses generation. No download or N-body fallback is implicit in `--spk`.

Use a local SPK only as a J2000 seed for N-body generation:

```sh
uv run --group generate libephemeris-generate-leb --tier medium \
  --kernel /path/to/de440.bsp --start 1850 --end 1852 --bodies CHIRON \
  --spk 15 /path/to/2060_190001_210001.bsp 20002060 --nbody CHIRON \
  --assist-planets /path/to/linux_p1550p2650.440 \
  --assist-asteroids /path/to/sb441-n16.bsp -o target/generated/chiron-1850.leb
```

The SPK must contain the explicit target in heliocentric J2000 coordinates at
J2000. Rust reads its position and analytic velocity, adds the Sun state from
the selected force ephemeris, and integrates SSB ICRS states with IAS15. Forward
and backward sweeps start independently from that same seed. Every asteroid
perturber remains active except a coincident target's own force. The installed
wheel includes the native REBOUND/ASSIST engines; their Python packages and v3
are not generation dependencies.

Both force files are explicit existing paths, hashed with the seed. The
requested domain must fit their actual common continuous coverage with the
established 32-day force-stage guard. Deep-time work must explicitly select an
available DE441 force file; the generator never silently chooses a shorter
dataset or clips the requested body interval. No SPK value beyond the anchor
initializes a later fitting sample. Output cannot overwrite any selected source
or force file, including with `--force`.

The native source matches all 22 independently read/propagated reference states
for Chiron and Eris, including nonmonotone requests, duplicates and both guarded
force boundaries. The check uses the pinned v3 generator's SPK/ASSIST model;
it does not certify complete-force JPL accuracy over centuries. Three newly
generated body files pass 1,863 off-grid position/rate checks and 192 fresh
native-reader source samples per observable. Chiron 1850–1852 and Eris
2200–2202 lie outside their seed SPK spans. With the matching generated cores,
their lossless conversions each preserve 1,440 complete TT/UT/frame pairs;
every public target is traced to LEB. See the
source comparison,
Chiron conversion
and Eris conversion.
These establish source preservation and fitting/conversion fidelity. The
inherited reduced-force absolute precision discrepancy remains open.

For complete minor-body apparent calls, merge the part with the generated core
over the same global domain. The public source checker requires that core and
confirms the actual LEB target trace. Chiron/Ceres/Pallas/Juno/Vesta generation
over 2000 passes 2,070 off-grid position and 2,070 Cartesian-rate checks; the
combined 19-body file passes all 1,824 complete TT/UT/frame comparisons with
explicit clock alignment. The independent type-2 JPL reader also passes 24
primary/support states including actual common coverage endpoints. See the
source check
and independent geometric check.

Standard named groups are consumed directly by the native runtime: select
`medium_core.leb2` and its sibling `medium_asteroids.leb2`, `medium_apogee.leb2`
etc. are discovered in the same directory. The indexed forced-LEB comparison
passes 2,112 complete pairs for all 22 generated core/asteroid/apogee bodies,
under unchanged criteria with no target substitute. Floating differences are
retained in the group comparison.

Automatic local assembly is supported under the explicit policies above.
Full tier/auxiliary precision, complete regeneration and the remaining
maintenance generators still need completion before claiming all v3 generation
capabilities. Unsupported sources fail explicitly.

The TT-grid correction found during per-body coverage validation replaces an
earlier generator mistake: Cartesian TDB samples had been stored under TT
labels. Historical 2000/core and N-body fitting reports above retain their
executed results, but require regeneration with the corrected clock before
they qualify current generation. The 1900 core plus shorter Chiron channel
passes 3,429 off-grid position/rate checks, 768/768 final-reader source samples
and all 1,152 complete aligned-clock position/rate pairs in each format. LEB1
uses parallel workers; direct LEB2 uses a fresh sequential process per body.
The original four lunar failures are retained separately; thresholds and epochs
are unchanged. See LEB1,
LEB2
and the earlier clock failure.

## Planet-center SPK generation

The native Rust CLI also assembles companion planet-center SPKs from existing
local NASA/JPL satellite kernels. It requires no Python generation package.

Prepare the source cache with an explicit native download command:

```sh
libephemeris download planet-centers-sources --tier base --data-dir /path/to/data
```

Omit `--tier` to prepare all tiers. Sources go to
`DATA_DIR/tmp/planet_centers`, the generator's default cache. Existing files
are reused after their available digest pins and SPK directory are checked.
Transfers stream to temporary files and publish atomically after validation;
`--force` replaces the BSP inputs. The pinned `naif0012.tls` companion is
prepared for external SPICE tooling; native TDB assembly does not consume it.
Base sources have no fixed digest in the catalog: their actual hashes identify
their bytes, while the generator checks their consumed segments and coverage.

Assemble one tier from that cache or select another existing cache explicitly:

```sh
libephemeris generate planet-centers-base --data-dir /path/to/data --json
libephemeris generate planet-centers-all --cache-dir /path/to/satellite-kernels \
  --output-dir /path/to/generated --json
```

The default cache is `DATA/tmp/planet_centers`; `--data-dir DATA` selects the
data root and the default output directory. The single-profile commands are
`planet-centers-base`, `planet-centers-medium`, `planet-centers-extended`;
`planet-centers-spk` aliases medium. `--force` explicitly permits replacement.
Every product is published atomically; failures preserve its previous file.

Base requires `jup204.bsp`, `sat319.bsp`, `ura083.bsp`, `nep050.bsp` and
`plu017.bsp`. Medium/extended require `jup365.bsp`, `sat441xl_part-1.bsp`,
`sat441xl_part-2.bsp`, `ura111xl-799.bsp`, `nep097xl-899.bsp` and `plu060.bsp`.
Wide-source SHA-256 identities are fixed; every consumed file is hashed before
and after copying. Missing sources cause an explicit error. This command does
not download inputs.

The writer copies all J2000 type-2/3 segments for planet centers 599, 699,
799, 899 and 999 relative to their system barycenters. It clips base descriptors
to 1850–2150 and medium to 1550–2650 in TDB, preserving complete coefficient
records and their edge neighbours. Extended retains the actual source domains.
There is no refitting, extrapolation or invented coverage.

The JSON report lists actual disjoint coverage intervals and `covers_tier` for
each body. For example, the pinned Jupiter source starts on 1600-01-10, and
Pluto covers 1800-01-02 through 2199-12-30: neither spans the whole medium tier.
`--require-full-tier` makes base/medium generation refuse incomplete coverage;
it does not apply to extended or all, whose source domains are body-specific.

An independent JPL-only comparison
checks every retained coefficient word and samples both native and public
jplephem readers at identical split TDB epochs. That scientific comparison is
maintained in the separate qualification repository; the CLI's JSON product
reports describe the generated output here.
The writer is self-contained. The separate native `generate keplerian-elements`
command builds multi-epoch element tables, and `leb generate` dispatches to the
optional standalone fitting application described above.
