# Command line

The Python console command and Rust binary run the same native application.
Exit codes: **0** success, **1** operational failure, **2** usage error.

## Inspect and configure

```sh
libephemeris --version
libephemeris init
libephemeris config --json
libephemeris status --tier base --data-dir /path/to/data --json
libephemeris status -vv
libephemeris leb info /path/to/base_core.leb2 --json
```

`-v` adds paths; `-vv` adds detailed inventory. Missing assets appear in status
rather than being provisioned. Configured IERS access can still prepare its
cache or acquire observations: use `LIBEPHEMERIS_NETWORK_POLICY=sealed` for
strictly offline inspection. `leb info` reports actual body/channel coverage.

For scripted setup, use `init --non-interactive -o libephemeris-config.toml`.
Existing output needs `--force`. See [configuration](data.md#configuration).

## Download only what you need

```sh
libephemeris download auto
libephemeris download base --no-progress
libephemeris download medium --data-dir /path/to/data
libephemeris download leb2-base
libephemeris download iers
libephemeris download assist --no-asteroids
```

- `auto` follows the configured tier; `base`, `medium`, `extended` are explicit.
- `leb-{tier}` / `leb2-{tier}` select local coefficients.
- `assist` means external force data, not the retired executable. It accepts
  `--no-planets`, `--no-asteroids` and `--target-dir`.
- `planet-centers-sources --tier base` acquires planet-centre sources.
- `all --yes` requests the full inventory and can use substantial disk space.
- `--data-dir`, `--quiet`, `--no-progress` and `--force` control destination,
  output and replacement. Existing valid files are reused; verified downloads
  are published atomically. A later failure can leave earlier completed files.

Explicit download commands authorize network access. Historical colon aliases
remain supported, but use the space-separated commands for new scripts.

## Verify and convert

```sh
libephemeris leb verify base --data-dir /path/to/data --json
libephemeris leb2 verify /path/to/base_core.leb2 --json
libephemeris leb2 verify /path/to/base_core.leb2 --reference /path/to/base.leb
libephemeris leb2 convert base-core --input /path/to/base.leb \
  --output /path/to/output/base_core.leb2 --max-chunk-bytes 1048576
```

Conversion preserves coefficient bits and verifies each output before
publication. Existing destinations need `--force`. Keep source files immutable
while indexed readers use them. Container verification is not a measurement
of apparent-position accuracy.

## Generate new data (optional)

Ordinary installation does not require generating or replacing datasets.
The separate [offline LEB tools](../tools/leb_generation/README.md) document
polynomial generation from explicit local sources.

```sh
libephemeris generate keplerian-dry-run --spk-dir /path/to/spk --body Ceres
libephemeris generate keplerian-elements --spk-dir /path/to/spk \
  --body Ceres Juno --start 1650 --end 2450 --spacing 10 \
  --format json --output elements.json --require-complete
libephemeris generate planet-centers-base --data-dir /path/to/data --json
```

Element generation reads local Sun-centred J2000 type-2/3/21 SPKs without
acquisition. Years label `2451545 + (year - 2000) * 365.25` JD TDB, not calendar
January 1. `--kernel BODY=PATH` resolves ambiguity and `--require-complete`
refuses omissions. Planet-centre generation copies existing JPL segments
without refitting; reported per-body coverage governs use.
