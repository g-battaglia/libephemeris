# Swiss Ephemeris Bug: SIDEREAL + J2000 for Lunar Nodes and Apsides

## Summary

pyswisseph silently ignores `FLG_J2000` for four lunar bodies when
`FLG_SIDEREAL` is also set:

| Body | `SIDEREAL + J2000` in pyswisseph |
|------|----------------------------------|
| `SE_MEAN_NODE` (10) | J2000 applied correctly |
| `SE_MEAN_APOG` (12) | J2000 applied correctly |
| `SE_TRUE_NODE` (11) | **J2000 silently ignored** |
| `SE_OSCU_APOG` (13) | **J2000 silently ignored** |
| `SE_INTP_APOG` (21) | **J2000 silently ignored** |
| `SE_INTP_PERG` (22) | **J2000 silently ignored** |

LibEphemeris intentionally corrects this: `FLG_J2000` is honored for
**all** bodies uniformly.

**Impact for users:** If you use `FLG_SIDEREAL` without `FLG_J2000`
(the vast majority of use cases), there is zero difference. The divergence
only affects the specific combination `FLG_SIDEREAL | FLG_J2000` on
these four bodies.

---

## How we detected it

During systematic validation of all sidereal flag combinations against
pyswisseph, we observed an internal inconsistency:

- `SE_MEAN_NODE` with `SIDEREAL | J2000` returns a **different** longitude
  than with `SIDEREAL` alone (J2000 precession is applied).
- `SE_TRUE_NODE` with `SIDEREAL | J2000` returns the **same** longitude
  as with `SIDEREAL` alone (J2000 precession is not applied).

Bodies from the same physical family respond differently to the same flag
combination. This inconsistency was confirmed across multiple dates, all
43 ayanamsha modes, and both longitude and latitude components.

No error or warning is emitted when the user requests `SIDEREAL | J2000`
for the affected bodies. The `J2000` flag is accepted but silently
discarded.

---

## Why this is incorrect

### 1. Ayanamsha and J2000 precession are distinct operations

- **Ayanamsha** shifts the longitude zero point from the vernal equinox to
  a sidereal reference. It is a 1D rotation along the ecliptic longitude
  coordinate.

- **J2000 precession** changes the reference plane from the ecliptic of
  date to the ecliptic of J2000.0. It is a 3D rotation that accounts for
  the ~47 arcseconds per century drift of the ecliptic plane.

These two operations are geometrically independent and composable. Applying
ayanamsha does not substitute for J2000 precession, and vice versa.

### 2. Internal inconsistency

If the design intent were that `SIDEREAL` and `J2000` are incompatible for
lunar nodes/apsides, the behavior should be consistent across all bodies
of the same type. Instead:

- Mean bodies (`SE_MEAN_NODE`, `SE_MEAN_APOG`) apply J2000 correctly.
- True/osculating/interpolated bodies do not.

This pattern indicates a code-path issue, not a deliberate API decision.

### 3. The error grows with distance from J2000

The delta between the correct and incorrect results grows proportionally
with time distance from J2000.0, consistent with missing ecliptic plane
precession:

| Epoch | TrueNode delta | Physical meaning |
|-------|---------------|------------------|
| J2000.0 | ~0.004° | Frame bias only |
| 2024 CE | ~0.34° | 24 years of ecliptic precession |
| J1900 | ~1.40° | 100 years of ecliptic precession |
| 3000 CE | ~14.0° | 1000 years of ecliptic precession |

If the difference were merely a convention choice, it would not produce
a systematic, time-dependent error of this form.

### 4. Physical sanity check fails

The true lunar node oscillates around the mean node with an amplitude of
approximately ±1.5°. At any epoch, the two should be within ~2° of each
other.

With the pyswisseph behavior at 0 CE:
- `|TrueNode - MeanNode|` with `SIDEREAL | J2000` = **~29°**
- This is physically impossible.

With the LibEphemeris fix:
- `|TrueNode - MeanNode|` with `SIDEREAL | J2000` = **~1.04°**
- This is physically correct.

---

## Numerical evidence

All measurements use Lahiri ayanamsha (`SE_SIDM_LAHIRI`).

### A. SID+J2K vs SID-only (LibEphemeris, after fix)

| Epoch | Body | SID+J2K lon | SID lon | Delta |
|-------|------|-------------|---------|-------|
| 2024 | TrueNode | 15.280° | 15.618° | -0.339° |
| 2024 | OscuApog | 190.790° | 191.128° | -0.339° |
| 2024 | MeanNode | 15.774° | 16.113° | -0.339° |
| 2024 | MeanApog | 169.639° | 169.978° | -0.339° |

All four bodies show a consistent ~0.339° J2000 precession shift — as
expected for a uniform coordinate transformation.

### B. TrueNode vs MeanNode physical sanity

| Epoch | With fix | Without fix (SE behavior) |
|-------|----------|--------------------------|
| 2024 CE | 0.49° | 0.49° |
| 0 CE | 1.04° | 28.85° |
| 3000 CE | 1.37° | 15.37° |

---

## The LibEphemeris fix

### Algorithm

For all Pipeline B bodies (MeanNode, MeanApog, TrueNode, OscuApog,
IntpApog, IntpPerg), when both `FLG_SIDEREAL` and `FLG_J2000` are
set:

1. Compute tropical ecliptic-of-date position.
2. Subtract **mean ayanamsha** (not true, because the J2000 ecliptic
   frame has no nutation component).
3. Precess from ecliptic of date to J2000 ecliptic.
4. If equatorial output is also requested, rotate to J2000 equatorial
   using J2000 obliquity.

This is the same pipeline already used correctly for `SE_MEAN_NODE` and
`SE_MEAN_APOG`. The fix simply extends it to the remaining four bodies.

### Code locations

The fix touches two files:

- **`libephemeris/fast_calc.py`** (LEB binary ephemeris path):
  - Removed `_SID_J2K_SKIP_BODIES` — J2000 precession is no longer
    suppressed for any body.
  - Extended `_deferred_sid_j2k` pattern to all ecliptic-direct bodies
    (was previously limited to mean bodies only).
  - Removed `_J2K_SKIP` from ayanamsha selection — mean ayanamsha is
    used for all bodies when J2000 is requested.

- **`libephemeris/planets.py`** (Skyfield computation path):
  - Removed `_eff_flags = iflag & ~FLG_J2000` logic from TrueNode,
    OscuApog, and IntpApog/IntpPerg handlers. These bodies now use
    `iflag` directly, matching the MeanNode/MeanApog pattern.

### Test coverage

- **`tests/test_sidereal/test_se_bug_j2k_nodes.py`**: Dedicated tests
  verifying J2000 is applied, physical sanity checks, LEB vs Skyfield
  consistency, and documented SE divergence magnitude.

- **`compare_scripts/tests/test_compare_sidereal_regression.py`**:
  Updated from "verify J2K suppression" to "verify intentional
  divergence" for the four affected bodies.

- **`tests/test_leb/compare/extended/test_extended_sidereal.py`**:
  Updated to verify LEB vs Skyfield agreement for SID+J2K on true
  bodies (was previously verifying J2K suppression).

---

## Note on methodology

This analysis was performed entirely through **black-box behavioral
observation** of pyswisseph. We compared the outputs of different flag
combinations across multiple dates, bodies, and ayanamsha modes. We did
not inspect the Swiss Ephemeris source code. The conclusions about the
nature of the bug are derived from the observed behavior, the mathematical
properties of the coordinate transformations involved, and the internal
inconsistency between mean and true body handling.
