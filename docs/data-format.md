# LEB binary format

LEB stores polynomial channels and auxiliary tables. All fields are
little-endian, without struct padding. Python exposes retained native reader
handles; Rust exposes `leb::LebFile`. See [data](data.md#large-leb-files) for
storage policy and immutable-file requirements.

## Header and sections

| Container | Accepted versions | Coefficients |
| --- | --- | --- |
| `LEB1` | 1 | Raw binary64 |
| `LEB2` | 1 | One compressed blob per body |
| `LEB2` | 2 | Independently compressed chunks |

LEB2 accepts compression flag **1** (`ZSTD_TRUNC_SHUFFLE`); zero and unknown
schemes are rejected. Compression is not inferred from the filename.

The 64-byte header (`<4sIIIdddI20s`) contains magic, version, section count,
body count, start JD, end JD, generation epoch, flags and 20 reserved bytes.
Each section-directory entry is 24 bytes (`<IIQQ`): section ID, reserved,
offset and size. IDs are 0 body index, 1 coefficients, 2 nutation, 3 Delta-T,
4 stars; orbital elements are reserved.

The published data-v4 bundled image has SHA-256
`b891de02f82b68346046da67389fd3713aa81ab299287c632cb20cc396397964`,
LEB2 v2/flag 1, five sections and fourteen bodies, JD 2396757.5–2506332.5.
Compared with data-v3, stored coverage includes one additional day at each end.
Always inspect a selected file's actual metadata; its filename does not
establish that identity.

## Body index and units

LEB1 entries are 52 bytes (`<iIIdddIIQ`): body ID, coordinate type, segment
count, start/end JD, interval in days, degree, components, absolute data offset.
LEB2 entries add compressed and uncompressed sizes (68 bytes,
`<iIIdddIIQQQ`). Counts must fit the enclosing section.

| Coordinate type | Stored channels |
| --- | --- |
| 0 `ICRS_BARY` | Planet-centre ICRS barycentric xyz in AU |
| 1 `ECLIPTIC` | Longitude/latitude in degrees, distance in AU |
| 2 retired | Historical channel; never reuse this ID; direct longitude reads still normalize |
| 3 reserved | Unused |
| 4 `ICRS_BARY_SYSTEM` | System-barycentre xyz in AU; planet-centre correction belongs to the position pipeline |

Body evaluation epochs are TT. Raw polynomial derivatives use TT days.
These stored quantities are not automatically apparent positions or velocities.
Use stored intervals verbatim; generators can fit non-integer day intervals
to an integer number of segments.

## Coefficient storage

Each segment contains component-major coefficients, degree + 1 binary64 words
per component. LEB1 segment bytes are `components * (degree + 1) * 8`.

LEB2 v1 decompresses a whole body. Its expected size is
`segment_count * components * (degree + 1) * 8`, checked before decoding.

LEB2 v2 starts each body with a 16-byte chunk header (`<IIfI`): chunk count,
reserved, nominal chunk interval in days, reserved. It is followed by 48-byte
entries (`<ddQQQII`): start/end JD, absolute blob offset, compressed size,
uncompressed size, first segment, segment count. Validate both directory bounds
and decoded sizes. The nominal chunk interval is normally 3652.5 days.

Chunk JD bounds guide lookup, but segment ownership is authoritative. After
binary-searching chunk dates, the reader checks the segment range and adjacent
chunks before reporting corruption. Rounded chunk boundaries must not change
which polynomial supplies a query.

## LEB2 compression

Generator order:

1. Explicit, lossy mantissa truncation using native-coordinate error budgets.
2. Lossless coefficient-major transpose.
3. Lossless eight-lane byte shuffle.
4. Zstandard compression.

Reader order: bounded zstd decode → unshuffle → segment-major reorder.
Decoded values are ordinary binary64; the reader does not repeat truncation.
The decoded length must equal the declared, validated size. Decode failures,
short frames and invalid offsets are corruption errors, not range misses.
A generator's component budget does not establish an apparent-place error bound.

## Auxiliary sections

| Section | Layout and meaning |
| --- | --- |
| Nutation | 40-byte `<dddIIII` header: start/end JD, interval, degree, components (2), segment count, reserved; then raw binary64 dpsi/deps in radians at TT |
| Delta-T | `<II` count/reserved header, then binary64 `(JD, Delta-T days)` pairs |
| Stars | 64-byte `<iddddddd4s` records: ID, J2000 RA/Dec (degrees), proper motions (degrees/year), parallax (arcseconds), radial velocity (km/s), magnitude, padding |

Zero nutation segments mean absent data and must not be evaluated. Delta-T
lookup clamps to endpoints and linearly interpolates within the table using
the right insertion point; a zero-width pair returns its left value. Star count
is section size divided by 64; a short trailing remainder is ignored.

## Evaluation and boundaries

Each body's stored interval is closed. Reject extrapolation and nonfinite
queries. The initial record estimate uses `(jd - start) / interval`; actual
ownership is the greatest rounded nominal record start not exceeding the
represented query. The first endpoint belongs to record zero and the last to
the final record. Correct quotient disagreement at rounded joins rather than
changing the stored interval.

```text
segment_start = jd_start + index * interval_days
midpoint = segment_start + 0.5 * interval_days
tau = clamp(2 * (jd - midpoint) / interval_days, -1, 1)
```

A split epoch selects the same record as its collapsed word, retaining the
ordered `(hi - midpoint) + lo` displacement in tau. Re-split retarded epochs
after stepping; derivative support must stay in the selected record.

Clenshaw evaluation, for each component:

```text
b1 = b2 = 0
for k from degree down to 1:
    b = c[k] + (2*tau)*b1 - b2
    b2 = b1
    b1 = b
value = c[0] + tau*b1 - b2
```

Scale normalized-time derivatives by `2 / interval_days`. Ecliptic longitudes
wrap into `[0, 360)`; ICRS channels do not. Reject zero record counts, invalid
intervals and resulting NaNs. Preserve the native operation association when
reproducing binary64 results.

## Errors

Missing files, unknown magic, unsupported versions/compression, corrupt
metadata/payloads, missing bodies and out-of-range epochs have distinct native
errors (`LebFileNotFound`, `LebUnknownFormat`, `LebUnsupportedVersion`,
`LebUnsupportedCompression`, `LebCorrupt`, `LebBodyNotFound`, `LebOutOfRange`).
Storage refusals use `ResourceLimit`, `InvalidLebReader` or
`LebBytesUnavailable`; closed Python handles report reader errors even for
presence probes. None of these conditions authorizes silently changing a
forced provider.

For writing/conversion, use the [offline tools](../tools/leb_generation/README.md).
Their tests and the native reader define executable checks for these layouts.
