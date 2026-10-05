# data-v4 validation

The 2026-10-05 rebuild completed all 344 phases: 151 per-body generation
checkpoints, their source comparisons, merged LEB1 exports and LEB2 conversion
checks. The release distributes the twelve LEB2 files consumed by the downloader.
Fifteen additional LEB1 exports remain available in the local build inventory.

## Integrity and reproducibility

All 27 exports passed a second structural inspection and full SHA-256 check
against the completed build manifest and checksum list. The generator and
library source hashes still matched their recorded generation inputs before
release metadata was changed. The twelve old local LEB2 files matched the
published data-v3 hashes, making them a reproducible comparison baseline.

The build used Python 3.12.13, NumPy 2.3.5, Skyfield 1.54, skyfield-data 7.0.0,
pyerfa 2.0.1.5, jplephem 2.24, zstandard 0.25.0, REBOUND 4.6.0 and ASSIST 1.2.3.
Sources are registered JPL kernels and the project's independently sourced
models. No reference-distribution outputs enter this build or comparison.

## What changed

- The canonical body inventory is unchanged: 53 bodies for base and medium,
  45 for extended. The retired uranians group is not reintroduced.
- Base coverage extends from JD 2396758.5–2506331.5 to
  2396757.5–2506332.5; medium extends from 2287185.5–2688952.5 to
  2287184.5–2688953.5. Extended and the narrower source-bounded companion
  intervals retain their previous coverage.
- Every file has a new digest; timestamps and auxiliary sections also affect
  bytes, so a changed digest alone does not imply changed body positions.
- At 40,113 sampled stored states across 151 tier/body combinations,
  medium/extended asteroid and exotic channels were numerically identical.
  Sampling combined a uniform grid, a seeded random grid, 1900–2100 dates,
  J2000 and both coverage edges. This is not an exhaustive accuracy proof.
- The public `calc()` comparison used explicit old/new readers and identical
  code, flags and dates: 4,928 calls, with no newly introduced errors. The
  largest sampled Sun/Moon/planet longitude change was 0.000492 arcseconds.
  With automatic cumulative routing, changed base channels also serve dates
  inside base coverage when medium or extended is configured.

## Differences requiring interpretation

Gonggong's base channel changed most: up to 8.85 arcseconds in the stored
barycentric direction and 8.45 arcseconds in the sampled public longitude.
At the worst stored-direction sample, comparison against the selected JPL
kernel gave a maximum-component angular error estimate of 6.72 arcseconds
for the old file and 0.0000192 arcseconds for the new compressed file.
Targeted checks at that epoch, J2000 and a modern epoch also confirmed closer
agreement for Asbolus, Varuna, Toutatis and Icarus. These are different metrics
from the public apparent geocentric longitude; they must not be conflated.

Stored mean-apogee channels differ by up to 1.43 arcseconds and interpolated
apsides by about 0.0642 arcseconds. The current public API serves these points
from its declared runtime models. The sampled public interpolated-apsides
results were unchanged; sampled mean-node/apogee changes stayed below
0.0000017 arcseconds, including the changed auxiliary frame interpolation.

The extended Asbolus trajectory is unchanged in the sampled old/new states.
Its generation verifier reports 3.5613 arcseconds of angular error estimate
against the finite Horizons window, near 1622–1624 CE. The generator explicitly
uses a 4-arcsecond verification budget for that numerical model. This does not
claim a precision improvement or extend that exception to direct-source bodies.

## Installation and rollback

The data-v4 tag supersedes data-v3 for libephemeris 3.2.2 downloads. The older
release remains immutable so previous library versions keep their pinned URLs
and rollback remains possible. Installing a library update does not rewrite an
existing user cache: run `libephemeris download leb2-extended` (or the desired
tier) to verify and replace its files through the atomic downloader.
