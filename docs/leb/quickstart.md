# LEB Generation Quickstart

> Quick reference for generating and regenerating LEB1 and LEB2 files.
> For full technical details, see [guide.md](guide.md).

---

## Prerequisites

```bash
uv pip install -e ".[dev]"
```

For asteroids, JPL Horizons SPK files are needed (downloaded automatically
during generation if `LIBEPHEMERIS_AUTO_SPK=1`).

---

## 1. Generate LEB1 Files

### Via leph (recommended)

```bash
# Group workflow (recommended — avoids macOS deadlocks, allows partial regen)
leph leb generate base groups       # 4 groups + merge → ephemeris_base.leb
leph leb generate medium groups     # same for medium tier
leph leb generate extended groups   # same for extended tier
```

### Step-by-step group workflow

```bash
# 1. Generate each group separately (via direct CLI)
python scripts/generate_leb.py --tier base --group planets      # → ephemeris_base_planets.leb    (~1s)
python scripts/generate_leb.py --tier base --group asteroids    # → ephemeris_base_asteroids.leb   (~15-60s)
python scripts/generate_leb.py --tier base --group exotics      # → ephemeris_base_exotics.leb     (~2-5min)
python scripts/generate_leb.py --tier base --group analytical   # → ephemeris_base_analytical.leb (~2-3min)

# 2. Merge + verify
python scripts/generate_leb.py --tier base --merge \
  data/leb/ephemeris_base_planets.leb \
  data/leb/ephemeris_base_asteroids.leb \
  data/leb/ephemeris_base_exotics.leb \
  data/leb/ephemeris_base_analytical.leb \
  --verify
```

### Regenerate a single group

If only asteroids changed (e.g. SPK update):

```bash
python scripts/generate_leb.py --tier base --group asteroids    # regenerate asteroids only
python scripts/generate_leb.py --tier base --merge \
  data/leb/ephemeris_base_planets.leb \
  data/leb/ephemeris_base_asteroids.leb \
  data/leb/ephemeris_base_exotics.leb \
  data/leb/ephemeris_base_analytical.leb \
  --verify                                                      # re-merge all 4 partial files
```

### Direct CLI

```bash
# Full tier
python scripts/generate_leb.py --tier base --verify

# Custom range
python scripts/generate_leb.py --output custom.leb --start 1900 --end 2100 --verify

# Specific bodies only
python scripts/generate_leb.py --output test.leb --start 2000 --end 2030 --bodies 0,1,2,3,4

# Single group + manual merge
python scripts/generate_leb.py --tier base --group planets
python scripts/generate_leb.py --tier base --group asteroids
python scripts/generate_leb.py --tier base --group exotics
python scripts/generate_leb.py --tier base --group analytical
python scripts/generate_leb.py --tier base --merge \
  data/leb/ephemeris_base_planets.leb \
  data/leb/ephemeris_base_asteroids.leb \
  data/leb/ephemeris_base_exotics.leb \
  data/leb/ephemeris_base_analytical.leb \
  --verify
```

### Output files

| Tier | File | Size | Coverage |
|------|------|------|----------|
| Base | `data/leb/ephemeris_base.leb` | ~53 MB | 1850-2150 |
| Medium | `data/leb/ephemeris_medium.leb` | ~175 MB | 1550-2650 |
| Extended | `data/leb/ephemeris_extended.leb` | ~1.6 GB | -5000 / +5000 |

### LEB1 body groups

| Group | Bodies | Method | Time (base) |
|-------|--------|--------|-------------|
| `planets` | Sun-Pluto, Earth (11) | Vectorized Skyfield | ~1s |
| `asteroids` | Chiron, Ceres, Pallas, Juno, Vesta (5) | spktype21 scalar | ~15-60s |
| `exotics` | Centaurs, TNOs, NEAs (31) | spktype21 scalar | ~2-5min |
| `analytical` | Nodes, Lilith, Uranians, Apogees (15) | Scalar Python | ~2-3min |

The four groups total 62 bodies.

---

## 2. Convert LEB1 → LEB2

> LEB2 requires an existing LEB1 file as input. Generate one first
> (section 1) if you don't have one.

### Via leph

```bash
# All 5 groups for a tier
leph leb2 convert base              # → data/leb2/base_{core,asteroids,exotics,apogee,uranians}.leb2

# Medium / Extended
leph leb2 convert medium
leph leb2 convert extended
```

### Direct CLI

```bash
# Single group
python scripts/generate_leb2.py convert data/leb/ephemeris_base.leb \
  -o data/leb2/base_core.leb2 --group core

# All groups at once
python scripts/generate_leb2.py convert-all data/leb/ephemeris_base.leb \
  -o data/leb2/ --tier-name base
```

### Generate LEB2 from scratch (no pre-existing LEB1)

Creates a temporary LEB1 and converts it automatically:

```bash
python scripts/generate_leb2.py generate --tier base --group core \
  -o data/leb2/base_core.leb2
```

### Output files (base tier)

| Group | File | Size | Bodies |
|-------|------|------|--------|
| core | `data/leb2/base_core.leb2` | 10.6 MB | Sun-Pluto, Earth, Nodes, Mean Apogee (14) |
| asteroids | `data/leb2/base_asteroids.leb2` | 8.7 MB | Chiron, Ceres, Pallas, Juno, Vesta (5) |
| exotics | `data/leb2/base_exotics.leb2` | 59.0 MB | Centaurs, TNOs, NEAs (31) |
| apogee | `data/leb2/base_apogee.leb2` | 11.4 MB | OscuApog, IntpApog, IntpPerig (3) |
| uranians | `data/leb2/base_uranians.leb2` | 2.1 MB | Cupido-Transpluto (9) |
| **Total** | | **91.7 MB** | **62 bodies** |

### LEB2 groups vs LEB1 groups

LEB2 groups differ from LEB1 groups:

| LEB1 | LEB2 | Difference |
|------|------|------------|
| `planets` (11) | `core` (14) | Core also includes Mean/True Node and Mean Apogee |
| `asteroids` (5) | `asteroids` (5) | Identical |
| `exotics` (31) | `exotics` (31) | Identical |
| `analytical` (15) | `apogee` (3) + `uranians` (9) | Analytical split into 2 groups |

---

## 3. Verification and Testing

### Verify LEB1

```bash
# Verification built into generation (add --verify)
python scripts/generate_leb.py --tier base --verify --verify-samples 1000

# Quick manual smoke test
python3 -c "
from libephemeris.leb_reader import LEBReader
reader = LEBReader('data/leb/ephemeris_base.leb')
pos, vel = reader.eval_body(0, 2451545.0)  # Sun at J2000
print(f'Sun ICRS: x={pos[0]:.10f} y={pos[1]:.10f} z={pos[2]:.10f} AU')
reader.close()
"
```

### Verify LEB2

```bash
# Against LEB1 reference
leph leb2 verify base

# Direct CLI
python scripts/generate_leb2.py verify data/leb2/base_core.leb2 \
  --reference data/leb/ephemeris_base.leb --samples 500
```

### Test suites

```bash
# LEB1
leph test leb-format all           # Unit tests (no @slow)
leph test leb-format precision     # Full precision suite

# LEB2
leph test leb2-format all          # Compression round-trip + reader tests
leph test leb2-format precision-base  # End-to-end precision, base (~15s)
leph test leb2-format precision-all   # All tiers (~45s)

# With LEB active (any format)
leph test leb-backend unit         # Unit tests in LEB mode
leph test leb-backend unit-fast    # Same, parallel
```

### Individual test files

```bash
pytest tests/test_leb/test_leb_format.py -v
pytest tests/test_leb/test_leb_reader.py -v
pytest tests/test_leb/test_fast_calc.py -v
pytest tests/test_leb/test_leb_compression.py -v
pytest tests/test_leb/test_leb2_reader.py -v
```

---

## 4. Release

```bash
# 1. Generate the LEB files
leph leb generate medium groups

# 2. Dry run
leph release leb-dry-run 1.0.0

# 3. Upload + update hashes
leph release leb 1.0.0

# 4. Commit
git add libephemeris/download.py
git commit -m "update LEB medium tier hash after regeneration"
```

---

## 5. Runtime Usage

```python
from libephemeris import set_leb_file, set_calc_mode

# LEB1 (single file)
set_leb_file("data/leb/ephemeris_base.leb")

# LEB2 (auto-discovers companion files)
set_leb_file("data/leb2/base_core.leb2")

# Or via environment variable
# export LIBEPHEMERIS_LEB=/path/to/file.leb

# Force LEB mode (error if no file available)
set_calc_mode("leb")
```

Auto-discovery: files in `~/.libephemeris/leb/` are found automatically.

```bash
# End-user download
libephemeris download leb-base       # ~53 MB
libephemeris download leb-medium     # ~175 MB
libephemeris download leb-extended   # ~1.6 GB
```

---

## Command Summary

| Action | Command |
|--------|---------|
| **LEB1 base (recommended)** | `leph leb generate base groups` |
| **LEB1 medium** | `leph leb generate medium groups` |
| **LEB1 extended** | `leph leb generate extended groups` |
| **LEB2 base (all groups)** | `leph leb2 convert base` |
| **LEB2 medium** | `leph leb2 convert medium` |
| **Verify LEB2** | `leph leb2 verify base` |
| **Test LEB1** | `leph test leb-format all` |
| **Test LEB2** | `leph test leb2-format all` |
| **Test LEB2 precision** | `leph test leb2-format precision-base` |
| **Release** | `leph release leb 1.0.0` |
