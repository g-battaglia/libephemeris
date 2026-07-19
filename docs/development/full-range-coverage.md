# Full Range Coverage

> **Status: implemented; retained as a historical design record.** The
> future-tense sections below describe the work as originally planned, not the
> current implementation. Local SPK discovery now prefers the widest kernel,
> the extended tier requests the supported 1600-2500 Horizons SPK interval,
> `scripts/download_max_range_spk.py` and `leph download spk-extended` are
> available, and ASSIST configuration is documented. The major-planet DE441
> range remains about -13200 to +17191. All former published LEB downloads are
> retired; the bundled base core and locally generated LEB1/LEB2 files are the
> supported accelerated paths. A later data-v3 audit corrected one important
> assumption in this historical plan: the complete ASSIST force model is
> bounded by the intersection with `sb441-n16.bsp`, not by planet-only DE441.
> The pinned guarded interval is JD -1200493.5 through 5008210.5, and generated
> rows are classified as a numerical model rather than a direct minor-body
> ephemeris.

Historical plan for extending minor-body coverage across the DE441 range.

## Table of Contents

- [Problem Statement](#problem-statement)
- [Current Situation](#current-situation)
- [Target Coverage](#target-coverage)
- [Bugs Found During Investigation](#bugs-found-during-investigation)
- [Implementation Steps](#implementation-steps)
- [Verification](#verification)
- [JPL Horizons SPK Limits](#jpl-horizons-spk-limits)

## Problem Statement

The extended precision tier uses DE441 which covers -13200 to +17191 for major
planets, but minor bodies (asteroids, TNOs, centaurs) fall back to low-precision
Keplerian propagation (~10-30 arcsec error) outside a narrow SPK window.

## Current Situation

| Source | Range | Precision | Status |
|--------|-------|-----------|--------|
| SPK (on disk, `*_190001_210001`) | 1900-2100 | Sub-arcsecond | Registered (by chance) |
| SPK (on disk, `*_185001_215001`) | 1850-2150 | Sub-arcsecond | Present but **ignored** |
| SPK (max from JPL Horizons) | 1600-2500 | Sub-arcsecond | **Not downloaded** |
| ASSIST (DE441 + sb441-n16 common interval) | JD -1200493.5 to 5008210.5 after guard | Numerical model; body/date dependent | Configured for data-v3 generation |
| Keplerian fallback | Model-dependent | Accuracy is body/date dependent; no universal bound | Traced outside stored/source-backed intervals where a curated model exists |

## Target Coverage

```
DE441 range:   -13200 ==============================================> +17191
SPK (Horizons):            |1600 ============ 2500|
ASSIST complete model:       |-1200493.5 ============== 5008210.5|
Keplerian:        outside each body's stored/source-backed interval
```

## Bugs Found During Investigation

### 1. `discover_local_spks` registers random SPK file (first-wins)

**File:** `libephemeris/spk_auto.py`, `discover_local_spks()` (~line 1634)

`os.listdir()` returns files in filesystem-dependent order. When multiple SPK
files exist for the same body (e.g. `2060_190001_210001.bsp` and
`2060_185001_215001.bsp`), whichever file is listed first gets registered.
The wider-range file may be silently ignored.

**Status:** Pending fix

### 2. `spk_date_range` for extended tier is wrong

**File:** `libephemeris/state.py`, line 64

Currently set to `("1550-01-01", "2650-01-01")`, but JPL Horizons only accepts
SPK requests in the range **1600-01-01 to 2500-01-01** (verified empirically).
This causes `_try_auto_spk_download()` to request impossible ranges, which fail
silently and fall through to Keplerian.

**Status:** Pending fix

### 3. Diagnostic `_get_source()` was lying about data source

**File:** `scripts/_tier_diagnostic.py`

Was reporting "SPK" for all minor bodies regardless of whether the date was
actually covered by the SPK file. Now checks actual SPK file coverage via
`get_spk_coverage()`.

**Status:** Fixed

## Implementation Steps

### Step 1: Fix `discover_local_spks` to prefer widest SPK

**Status:** Pending

**File:** `libephemeris/spk_auto.py`

When a body is already registered and a new SPK file containing the same body
is found, compare coverage ranges. Re-register if the new file has wider
coverage.

```python
# Pseudocode for the fix:
if ipl in state._SPK_BODY_MAP:
    existing_file, _ = state._SPK_BODY_MAP[ipl]
    existing_coverage = get_spk_coverage(existing_file)
    new_coverage = get_spk_coverage(filepath)
    if new_coverage and existing_coverage:
        new_span = new_coverage[1] - new_coverage[0]
        existing_span = existing_coverage[1] - existing_coverage[0]
        if new_span > existing_span:
            # Re-register with wider file
            spk.register_spk_body(ipl, filepath, target_naif)
```

### Step 2: Update `spk_date_range` for tiers

**Status:** Pending

**File:** `libephemeris/state.py`

Update the extended tier to reflect the actual JPL Horizons maximum:

```python
"extended": PrecisionTier(
    name="extended",
    ephemeris_file="de441.bsp",
    spk_date_range=("1600-01-01", "2500-01-01"),  # Was 1550-2650
    description="Extended range (-13200 to +17191), ~3.1 GB",
),
```

Also update base and medium tiers if their ranges exceed what Horizons accepts:
- **base:** `("1850-01-01", "2150-01-01")` — OK, within 1600-2500
- **medium:** `("1900-01-01", "2100-01-01")` — OK, within 1600-2500

### Step 3: Create max-range SPK download script

**Status:** Pending

**File:** new `scripts/download_max_range_spk.py`

Downloads SPK files with the maximum range (1600-2500) for all 21 bodies in
`SPK_BODY_NAME_MAP`. Replaces existing narrower-range files.

Bodies (from `libephemeris/constants.py`):
- Chiron (2060), Pholus (5145), Ceres (Ceres;), Pallas (Pallas;),
  Juno (Juno;), Vesta (Vesta;), Eris (136199), Sedna (90377),
  Haumea (136108), Makemake (136472), Ixion (28978), Orcus (90482),
  Quaoar (50000), Varuna (20000), Nessus (7066), Asbolus (8405),
  Chariklo (10199), Gonggong (225088), Apophis (99942),
  Hygiea (Hygiea;), Eros (433)

### Step 4: Add poe task for max-range download

**Status:** Pending

**File:** `pyproject.toml`

```toml
"spk:download:maxrange" = "python scripts/download_max_range_spk.py"
```

### Step 5: Document ASSIST configuration (historical proposal)

**Status:** Pending

The original proposal assumed that adding the long DE441 planet file extended
the complete force model to -13000..+17000. That is false: ASSIST also consumes
the asteroid perturber SPK, and its narrower coverage is binding. Current LEB
generation reads both files, reserves an IAS15 boundary guard, and serializes
only their common interval.

**Requirements:**

```bash
pip install rebound assist

# Download ASSIST data files (~1 GB total):
# Option A: Use existing de441.bsp (already ~3.1 GB)
# Option B: Download Linux-format DE441:
curl https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de441/linux_m13000p17000.441 \
    -o data/linux_m13000p17000.441

# Required: asteroid perturbers for n-body integration:
curl https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp \
    -o data/sb441-n16.bsp

# Set environment variable:
export ASSIST_DIR=/path/to/data/
```

The `AssistEphemConfig` in `rebound_integration.py` searches for ephemeris files
in this order: `de441.bsp`, `de440.bsp`, `linux_m13000p17000.441`,
`linux_p1550p2650.440`.

## Verification

After implementation, run the tier diagnostics to verify:

```bash
leph diag positions-extended
```

Expected results:
- Dates 1600-2500: all minor bodies show **SPK**
- Dates inside the guarded common ASSIST source interval may show **ASSIST** or
  a stored `numerical-model` LEB trajectory, depending on mode and body.
- Dates outside the actual per-body source interval show **Keplerian** when the
  curated deterministic model exists; they are never inferred from the
  planet-only DE441 range.

## JPL Horizons SPK Limits

Tested empirically against the Horizons API on 2026-02-20:

| Parameter | Value |
|-----------|-------|
| Earliest START_TIME | `1600-01-01` |
| Latest STOP_TIME | `2500-01-01` |
| Max span | 900 years |
| Applies to | All 21 bodies in SPK_BODY_NAME_MAP |
| Bodies using name syntax | Ceres;, Pallas;, Juno;, Vesta;, Hygiea; |

Requesting dates outside this range returns:
`"START time outside set SPK limits"` or `"STOP time outside set SPK limits"`.
