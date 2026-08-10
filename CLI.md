# CLI Reference

libephemeris has three CLI interfaces:

| CLI | Purpose | Entry point |
|-----|---------|-------------|
| `leph` | Developer tools: tests, generation, lint, release | repo launcher `./leph` |
| `libephemeris` | End-user: download data, check status | `libephemeris.cli:main` |
| `poe` | Curated shortcuts that delegate to `leph` | poethepoet (pyproject.toml) |

All three support `-h` / `--help` on every command and subcommand.

## Setup

### Install in dev mode

```bash
uv pip install -e ".[dev]"
```

This installs the production `libephemeris` console script. The developer CLI
is intentionally repo-only; invoke its checked-in launcher from the source
checkout.

### Running commands

If the virtualenv is activated, use directly:

```bash
./leph test skyfield essential
libephemeris download medium
poe lint
```

If the virtualenv is **not** activated, prefix with `uv run`:

```bash
uv run ./leph test skyfield essential
uv run libephemeris download medium
uv run poe lint
```

### TAB completion + `leph` shell function

Click's built-in completion doesn't work with `uv run` because it requires
`leph` to be in PATH. The `leph completion` command generates custom scripts
that solve this: they create a `leph` shell function (so you can type `leph`
directly without activating the venv) and wire up TAB completion through
`uv run` transparently.

#### Quick test (try it now, no permanent changes)

```bash
eval "$(uv run ./leph completion zsh)"
```

After running this, `leph` and TAB completion work in that terminal session.

#### Permanent setup

**zsh** (most macOS / modern Linux):

```bash
# Option A: append directly to .zshrc
uv run ./leph completion zsh >> ~/.zshrc

# Option B: separate file (cleaner)
uv run ./leph completion zsh > ~/.leph-completion.zsh
echo 'source ~/.leph-completion.zsh' >> ~/.zshrc
```

**bash**:

```bash
uv run ./leph completion bash >> ~/.bashrc
```

**fish**:

```bash
uv run ./leph completion fish > ~/.config/fish/conf.d/leph.fish
```

Then reload your shell (`source ~/.zshrc`, `source ~/.bashrc`, or restart fish).

#### What it does

The generated script does two things:

1. **Shell function**: creates a repo-anchored wrapper around
   `uv run --project /path/to/libephemeris python -m libephemeris.dev_cli`, so
   `leph test skyfield essential` works from any directory without activating
   the virtualenv. Regenerate the completion file after moving the checkout.

2. **TAB completion**: wires up completions that work at every nesting level,
   with descriptions:

```
leph <TAB>                         # test, code, leb, download, ...
leph test <TAB>                    # skyfield, leb-backend, leb-format, ...
leph test skyfield <TAB>           # essential, smoke, unit, unit-fast, ...
leph test leb-format vs-skyfield <TAB>  # base, medium, extended, crosstier, ...
```

#### libephemeris (production CLI) completion

The production CLI uses the same Click mechanism. If you want completion for
`libephemeris` too:

```bash
# zsh (requires venv activated or libephemeris in PATH)
eval "$(_LIBEPHEMERIS_COMPLETE=zsh_source libephemeris)"

# bash
eval "$(_LIBEPHEMERIS_COMPLETE=bash_source libephemeris)"

# fish
_LIBEPHEMERIS_COMPLETE=fish_source libephemeris | source
```

---

## Dev CLI (`leph`)

The dev CLI has over 100 commands organized in 11 subgroups. Every command name
makes it clear WHAT is being tested/generated and WHICH backend is used.

### `leph code` — Code quality

```bash
leph code lint                # Ruff linter with auto-fix
leph code format              # Ruff formatter (line-length 88)
leph code format-black        # Black formatter (legacy)
leph code typecheck           # mypy static type checker
```

### `leph test` — Run tests

**Never run the full test suite.** Always pick a targeted subcommand.

#### `leph test skyfield` — Unit tests via Skyfield/DE440

Positions are computed in real time from DE440 binary kernels via Skyfield.

```bash
leph test skyfield essential      # ~900 tests, ~20s — fast sanity check
leph test skyfield smoke          # ~1,700 tests, ~30s — broader sanity check
leph test skyfield unit           # ~16,000 tests, ~2 min — all unit tests, sequential
leph test skyfield unit-fast      # ~16,000 tests, ~1 min — all unit tests, parallel
leph test skyfield unit-full      # ALL tests including @slow (~10+ min)
leph test skyfield all            # unit + compare, sequential
leph test skyfield all-fast       # unit + compare, parallel
leph test skyfield all-full       # everything including @slow, sequential
leph test skyfield all-full-fast  # everything including @slow, parallel
leph test skyfield progress       # dots-only output (CI-friendly)
```

#### `leph test leb-backend` — Sealed runtime via LEB precomputed data

This curated suite verifies the source-sealed LEB contract: no-JPL routing,
vector adapters, product minor bodies, fixed stars, events, inventory, network
policy, and comparison precision. The generic tree deliberately remains under
`skyfield`; it contains JPL-only APIs and dates outside the bundled base tier.

```bash
leph test leb-backend essential   # Minimal sealed contract (parallel)
leph test leb-backend unit        # Curated contract, sequential and verbose
leph test leb-backend unit-fast   # Curated contract, parallel [RECOMMENDED]
leph test leb-backend unit-full   # Curated contract including @slow
```

#### `leph test lunar` — Lunar module

```bash
leph test lunar all       # All lunar tests (nodes, Lilith, perigee, apogee)
leph test lunar perigee   # Interpolated-model + osculating perigee contracts
leph test lunar apogee    # Mean/interpolated/osculating apogee contracts
leph test lunar lilith    # Mean + osculating Lilith geometry and precision
```

#### `leph test leb-format` — LEB binary format internals

Tests the file format itself (not using LEB as a backend).

```bash
leph test leb-format all                         # Reader/writer/fast_calc tests
leph test leb-format precision                   # Chebyshev error sweep, all tiers
leph test leb-format precision-quick             # Medium tier only

# LEB vs Skyfield accuracy per tier
leph test leb-format vs-skyfield base            # Base tier
leph test leb-format vs-skyfield base-quick      # Base tier, no @slow
leph test leb-format vs-skyfield medium          # Medium tier
leph test leb-format vs-skyfield medium-quick    # Medium tier, no @slow
leph test leb-format vs-skyfield extended        # Extended tier
leph test leb-format vs-skyfield extended-quick  # Extended tier, no @slow
leph test leb-format vs-skyfield crosstier       # Cross-tier consistency
leph test leb-format vs-skyfield all             # All tiers combined
leph test leb-format vs-skyfield legacy          # Legacy flat tests
leph test leb-format vs-skyfield legacy-quick    # Legacy, no @slow
```

#### `leph test leb2-format` — LEB2 compressed format

```bash
leph test leb2-format all                # Compression + reader unit tests
leph test leb2-format precision-base     # Core LEB2 vs LEB1, base tier (~15s)
leph test leb2-format precision-medium   # Core LEB2 vs LEB1, medium tier (~15s)
leph test leb2-format precision-extended # Core LEB2 vs LEB1, extended tier (~15s)
leph test leb2-format precision-all      # Core companions, all tiers (~45s)
```

#### `leph test horizons` — Horizons API precision (requires internet)

```bash
leph test horizons precision        # 13 bodies, 200 dates, 6 flags (~45s)
leph test horizons precision-quick  # 50 dates (~15s)
leph test horizons vs-leb           # Cross-validate against LEB2
```

#### `leph test coverage` — Coverage reports

```bash
leph test coverage run    # Standard run (term + XML + HTML), no @slow
leph test coverage full   # Including @slow
```

### `leph leb` — LEB1 binary ephemeris

#### Generation

Recommended workflow (avoids macOS multiprocessing deadlocks):

```bash
leph leb generate medium groups    # planets + asteroids + exotics + analytical + merge
```

All generation modes per tier (`base`, `medium`, `extended`):

```bash
leph leb generate <tier> groups      # Recommended: group-by-group then merge
leph leb generate <tier> full        # All eligible bodies at once (may deadlock on macOS)
leph leb generate <tier> single      # One eligible body at a time (lowest memory)
leph leb generate <tier> planets     # Planets group only
leph leb generate <tier> asteroids   # Asteroids group only
leph leb generate <tier> exotics     # Exotic registry group only
leph leb generate <tier> analytical  # Analytical group only
leph leb generate <tier> merge       # Merge partial files
leph leb generate <tier> body <name> # Specific body (e.g. 'moon', '1,2,3')
leph leb generate all                # All three tiers sequentially
```

#### Verification

```bash
leph leb verify base
leph leb verify medium
leph leb verify extended
```

### `leph leb2` — LEB2 compressed ephemeris

```bash
# Convert LEB1 -> LEB2
leph leb2 convert base              # All active groups (includes exotics)
leph leb2 convert medium
leph leb2 convert extended
leph leb2 convert base-core         # One group only (~10.2 MB, for PyPI)
leph leb2 convert base-exotics      # The exotic registry group only
leph leb2 convert medium-apogee     # Any tier/group combination is available

# Verify the exact 14-body core in native stored component units
leph leb2 verify base               # base_core.leb2 against LEB1
leph leb2 verify medium
leph leb2 verify extended
```

### `leph download` — Data files

```bash
# Full bootstrap for developers (downloadable prerequisites only)
leph download all             # DE/SPK + IERS + planet-center source kernels + ASSIST

# SPK kernels (required for Skyfield backend and LEB generation)
leph download spk-base        # DE440s + asteroid SPKs (1850-2150)
leph download spk-medium      # DE440 + asteroid SPKs (1900-2100)
leph download spk-extended    # de441 + max-range SPKs (1600-2500)

# Time / Earth-orientation data
leph download iers            # finals2000A.data + leap_seconds.dat + deltat.data

# Source kernels used by `leph generate planet-centers-*`
leph download planet-centers-sources   # All tiers
leph download planet-centers-sources --tier medium

# ASSIST n-body data
leph download assist           # ~120 MB (requires libephemeris[nbody])
```

### `leph diag` — Diagnostics

```bash
leph diag positions-base       # Print all body positions for base tier
leph diag positions-medium     # Print all body positions for medium tier
leph diag positions-extended   # Print all body positions for extended tier
leph diag download-data        # Download all required files for current tier
```

### `leph generate` — Derived data

```bash
# Planet center-of-body SPKs (sub-arcsecond gas giants)
leph generate planet-centers-base
leph generate planet-centers-medium
leph generate planet-centers-extended
leph generate planet-centers-all       # All three tiers
leph generate planet-centers-spk       # Legacy alias for medium

# Keplerian orbital elements
leph generate keplerian-elements       # Multi-epoch elements for 37 bodies
leph generate keplerian-dry-run        # Preview available SPKs
```

### Lunar model reproducibility

Mean lunar points are evaluated directly from ERFA/IERS fundamental arguments.
Interpolated apogee and perigee are evaluated at runtime from the active NASA
JPL state with a deterministic symmetric smoothing kernel. There are no lunar
compatibility-model artifacts or oracle-output generators to rebuild. See
`docs/methodology/lunar-apsides.md`.

### `leph release` — Release management

The historical LEB upload tasks are retired. The cumulative `data-v3` GitHub
release contains all four SHA-256-pinned LEB2 groups for base, medium, and extended;
the wheel bundles a byte-identical copy of `base_core.leb2`.
Maintainers rebuild the complete matrix with
`./regenerate-leb.sh all`; the legacy release workflow must not be used.

### `leph manual` — Documentation builds

```bash
# Pandoc workflow (requires: brew install pandoc tectonic)
leph manual build              # EPUB + PDF, Italian + English
leph manual build-epub         # EPUB only
leph manual build-pdf          # PDF only
leph manual build-it           # Italian only
leph manual build-en           # English only

# Ebooklib workflow (no pandoc required, Kobo-compatible)
leph manual generate-epub      # All EPUBs
leph manual generate-epub-it   # Italian EPUB
leph manual generate-epub-en   # English EPUB
```

---

## Production CLI (`libephemeris`)

For end-users: download data files and check library status.

```bash
# Download data
libephemeris download auto              # Exactly what the configured mode/tier needs
libephemeris download base              # DE440s + planet centers + SPKs (1850-2150)
libephemeris download medium            # DE440 + planet centers + SPKs (default)
libephemeris download extended          # DE441 + planet centers + SPKs
libephemeris download assist            # ASSIST n-body data (~714 MB)

# Status (comprehensive: version, config, all data files)
libephemeris status                     # Formatted text output
libephemeris status --json              # Machine-readable JSON

# Configuration guide (env vars, file locations, Python API)
libephemeris config                     # Show all settings and how to change them

# Generate a TOML config file (interactive wizard)
libephemeris init                       # Interactive wizard -> libephemeris-config.toml
libephemeris init -o config.toml        # Custom output path
libephemeris init --non-interactive     # All defaults, no prompts
libephemeris init --force               # Overwrite existing file

# Options available on all download commands
libephemeris download medium --force           # Re-download even if present
libephemeris download medium --no-progress     # Suppress progress bars
libephemeris download medium --quiet           # Suppress all output
```

The historical monolithic LEB download commands remain retired. In `auto` or
`leb` mode, `libephemeris download auto` installs every SHA-256-pinned LEB2
group cumulatively through the configured tier: 4 files for base, 8 through
medium, or 12 through extended. The wheel supplies the bundled
`base_core.leb2` and the remaining files come from the immutable `data-v3` release.

### `libephemeris status` output

The `status` command shows a comprehensive overview:

- **Configuration**: version, calc mode, precision tier, LEB file, data directory
- **Ephemeris Kernels**: de440s.bsp, de440.bsp, de441.bsp with [OK]/[--] and sizes
- **Planet Center Corrections**: per-tier BSP files, active tier marked with `*`
- **Locally generated LEB1 Binary Ephemeris**: per-tier .leb files, active one marked with `*`
- **LEB2 Compressed Ephemeris**: per-tier group counts
  (core/asteroids/exotics/apogee)
- **SPK Asteroid Cache**: directory path, file count, total size
- **ASSIST N-body Data**: planet ephemeris + asteroid perturbers
- **IERS Earth Orientation Data**: finals2000A, leap seconds, delta T with age in days
- **Commands**: setup hints with correct command syntax

### Backward compatibility

The reader can decode legacy LEB files that contain the retired hypothetical
group: pre-3.1.0 `{tier}_uranians.leb2` companions left on disk are recognized
by the tier guard and ignored harmlessly — fictitious-body channels are never
consulted, and the `<tier>-uranians` selectors are gone. Since 3.0.0rc15 a
file named after a manifest artifact is trusted on presence; its SHA-256 is
verified once, when the installer writes it.

Old colon-separated syntax still works:

```bash
libephemeris download:medium      # Same as: libephemeris download medium
```

### `libephemeris config`

The `config` command is a reference guide showing every configurable setting:

- **Data directory**: current path, default, env var
- **Precision tier**: base/medium/extended with kernel sizes and date ranges
- **Ephemeris file**: DE kernel override
- **Calculation mode**: auto/skyfield/leb/horizons with descriptions
- **LEB binary ephemeris**: bundled base core and locally generated LEB1/LEB2 listings
- **SPK asteroid cache**: directory, auto-download toggle
- **Planet center corrections**: per-tier BSP files
- **IERS Earth orientation data**: file names, env vars
- **Delta T model**: current model, env var, Python API, and TOML key
- **LEB mmap preloading**: current toggle and TOML year range
- **ASSIST n-body data**: file names, sizes, install requirements
- **.env file**: location and override env var
- **TOML config file**: path, loaded values, generate command
- **Logging**: log level env var
- **Other settings**: strict precision mode

Each section shows the current value, relevant env var, and Python API equivalent.

### `libephemeris init`

Interactive wizard that generates a `libephemeris-config.toml` file.
The file can be committed to version control so every developer on the
project uses the same library settings.

```bash
libephemeris init                       # Interactive wizard
libephemeris init -o config.toml        # Custom output path
libephemeris init --non-interactive     # All defaults, no prompts
libephemeris init --force               # Overwrite existing file
```

Generated file example:

```toml
[libephemeris]
precision = "extended"       # base | medium | extended
mode = "leb"                 # auto | skyfield | leb | horizons
network_policy = "sealed"    # auto | allow | sealed (LEB defaults sealed)
auto_spk = true              # auto-download SPK for minor bodies
strict_precision = true      # require SPK for major asteroids
iers_auto_download = false   # auto-download IERS data
iers_delta_t = false         # use observed Delta T from IERS
deltat_model = "smh2016"     # smh2016 | espenak_meeus
log_level = "WARNING"        # DEBUG | INFO | WARNING | ERROR | CRITICAL
mmap_preload = false         # pre-fault mmap pages for a date range
mmap_preload_start = 1800    # start year (when mmap_preload = true)
mmap_preload_end = 2200      # end year (when mmap_preload = true)
```

`network_policy = "sealed"` blocks every calculation-time HTTP path before
DNS or socket work. Explicit `libephemeris download ...` commands are the
provisioning boundary and remain authorized. Use `libephemeris leb info FILE`
to inspect the exact bodies and per-body Julian Day ranges stored in an
artifact.

### Configuration resolution order

When multiple sources set the same value, highest priority wins:

1. **`set_*()`** function calls (programmatic, in code)
2. **Environment variables** (including `.env` file)
3. **TOML config file** (`libephemeris-config.toml`)
4. **Built-in defaults**

TOML file search order:

1. `LIBEPHEMERIS_CONFIG` env var (explicit path)
2. `./libephemeris-config.toml` (project directory)
3. `~/.libephemeris/config.toml` (global user config)

---

## Poe shortcuts

Curated aliases for the most common workflows. Every shortcut delegates to `leph`.

Backend shortcuts expose the documented **core / fast / full** scopes below.

```bash
# Code quality
poe lint                                # -> leph code lint
poe format                              # -> leph code format
poe typecheck                           # -> leph code typecheck

# Skyfield backend: core / fast / full
poe test:skyfield:core                  # -> leph test skyfield essential (~900, ~20s)
poe test:skyfield:fast                  # -> leph test skyfield unit-fast (~16,000, ~1 min)
poe test:skyfield:full                  # -> leph test skyfield all-full-fast (+@slow)

# LEB backend: core / fast / full
poe test:leb:core                       # -> minimal sealed LEB contract
poe test:leb:fast                       # -> curated sealed contract [RECOMMENDED]
poe test:leb:full                       # -> curated contract (+@slow)

# Horizons backend: core / fast / full
poe test:horizons:core                  # -> leph test horizons precision-quick (~15s)
poe test:horizons:fast                  # -> leph test horizons precision (~45s)
poe test:horizons:full                  # -> leph test horizons vs-leb

# Coverage
poe coverage                            # -> leph test coverage run

# LEB1 generation + verification (resumable group workflow)
poe leb:generate:base                   # -> leph leb generate base groups
poe leb:generate:medium                 # -> leph leb generate medium groups
poe leb:generate:extended               # -> leph leb generate extended groups
poe leb:verify:base                     # -> leph leb verify base
poe leb:verify:medium                   # -> leph leb verify medium
poe leb:verify:extended                 # -> leph leb verify extended

# LEB2 conversion + verification
poe leb2:convert:base                   # -> leph leb2 convert base
poe leb2:convert:medium                 # -> leph leb2 convert medium
poe leb2:convert:extended               # -> leph leb2 convert extended
poe leb2:verify:base                    # -> leph leb2 verify base
poe leb2:verify:medium                  # -> leph leb2 verify medium
poe leb2:verify:extended                # -> leph leb2 verify extended
```

---

## Direct pytest usage

For single-file or single-test runs, use pytest directly:

```bash
pytest tests/test_file.py -v                          # One file
pytest tests/test_file.py::TestClass::test_method -v  # One test
pytest tests/test_file.py -k "pattern" -v             # By keyword
pytest tests/test_file.py -v --calc-mode leb           # One file with LEB backend
```

---

## Architecture

### Dev CLI source: `libephemeris/dev_cli/`

| Module | Commands | Purpose |
|--------|----------|---------|
| `__init__.py` | 1 | Root group, registers all subgroups |
| `cmd_code.py` | 4 | lint, format, format-black, typecheck |
| `cmd_test.py` | 41 | All test suites with backend/target naming |
| `cmd_leb.py` | 31 | LEB1 generate (per-tier, per-group, body) + verify |
| `cmd_leb2.py` | 21 | LEB2 convert + verify (per tier) |
| `cmd_download.py` | 7 | SPK, LEB, ASSIST downloads |
| `cmd_diag.py` | 4 | Tier diagnostics, data download |
| `cmd_generate.py` | 7 | Planet-center SPKs and Keplerian elements |
| `cmd_release.py` | 5 | Retired LEB publication commands (fail closed) |
| `cmd_manual.py` | 8 | Manual build (EPUB/PDF, pandoc/ebooklib) |
| `cmd_completion.py` | 3 | Shell completion scripts (zsh, bash, fish) |

### Production CLI source: `libephemeris/cli/__init__.py`

Click-based CLI with `download` subgroup (JPL tier data and ASSIST),
comprehensive `status` command (with `--json`), and `config` reference guide.
Backward compatible with old colon-separated syntax via alias rewriting.

### Shared: `libephemeris/cli/shared.py`

Tier metadata (`TIER_INFO`), download help text generators — single source of
truth reused by both CLIs.
