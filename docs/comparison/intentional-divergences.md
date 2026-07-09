# Intentional divergences from Swiss Ephemeris

libephemeris aims for 1:1 behavioural compatibility with pyswisseph. In a small
number of cases where Swiss Ephemeris has a demonstrable behavioural bug or returns
a physically unbounded quantity, libephemeris deliberately departs from it to
return the physically correct value. Each case is documented here so the difference
is never a surprise.

---

## 1. `SIDEREAL | J2000` for lunar nodes and apsides

**Status:** intentional divergence since the leb/precision branch.

### Summary

pyswisseph silently ignores `FLG_J2000` for four lunar bodies when `FLG_SIDEREAL`
is also set:

| Body | `SIDEREAL + J2000` in pyswisseph |
|------|----------------------------------|
| `MEAN_NODE` (10) | J2000 applied correctly |
| `MEAN_APOG` (12) | J2000 applied correctly |
| `TRUE_NODE` (11) | **J2000 silently ignored** |
| `OSCU_APOG` (13) | **J2000 silently ignored** |
| `INTP_APOG` (21) | **J2000 silently ignored** |
| `INTP_PERG` (22) | **J2000 silently ignored** |

libephemeris honours `FLG_J2000` for **all** bodies uniformly.

**Impact for users:** if you use `FLG_SIDEREAL` without `FLG_J2000` (the vast
majority of use cases), there is zero difference. The divergence only affects the
specific combination `FLG_SIDEREAL | FLG_J2000` on these four bodies.

### How it was detected

During systematic validation of all sidereal flag combinations, an internal
inconsistency appeared: `MEAN_NODE` with `SIDEREAL | J2000` returns a *different*
longitude than with `SIDEREAL` alone (J2000 precession applied), while `TRUE_NODE`
with `SIDEREAL | J2000` returns the *same* longitude as with `SIDEREAL` alone
(J2000 precession not applied). Bodies of the same physical family respond
differently to the same flag combination. This was confirmed across multiple dates,
all ayanamsha modes, and both longitude and latitude. No error or warning is
emitted; the `J2000` flag is accepted but silently discarded.

### Why this is incorrect

1. **Ayanamsha and J2000 precession are distinct operations.** Ayanamsha shifts the
   longitude zero point along the ecliptic (a 1D rotation). J2000 precession changes
   the reference plane from the ecliptic of date to the ecliptic of J2000.0 (a 3D
   rotation accounting for the ~47"/century drift of the ecliptic plane). They are
   geometrically independent and composable; applying one does not substitute for
   the other.
2. **Internal inconsistency.** If `SIDEREAL` and `J2000` were intended to be
   incompatible for these bodies, the behaviour should be consistent across the
   family. Instead mean bodies apply J2000 and true/osculating/interpolated bodies
   do not — a code-path issue, not an API decision.
3. **The error grows with distance from J2000**, consistent with missing
   ecliptic-plane precession:

   | Epoch | TrueNode delta | Physical meaning |
   |-------|---------------|------------------|
   | J2000.0 | ~0.004° | Frame bias only |
   | 2024 CE | ~0.34° | 24 years of ecliptic precession |
   | J1900 | ~1.40° | 100 years of ecliptic precession |
   | 3000 CE | ~14.0° | 1000 years of ecliptic precession |

4. **Physical sanity check fails.** The true node oscillates around the mean node
   with amplitude ~±1.5°; the two should be within ~2° at any epoch. With the
   pyswisseph behaviour at 0 CE, `|TrueNode − MeanNode|` under `SIDEREAL | J2000`
   reaches ~29° — physically impossible. With the libephemeris fix it returns to
   ~1.04°.

### Numerical evidence

All measurements use Lahiri ayanamsha (`SIDM_LAHIRI`).

**SID+J2K vs SID-only (libephemeris, after fix):**

| Epoch | Body | SID+J2K lon | SID lon | Delta |
|-------|------|-------------|---------|-------|
| 2024 | TrueNode | 15.280° | 15.618° | −0.339° |
| 2024 | OscuApog | 190.790° | 191.128° | −0.339° |
| 2024 | MeanNode | 15.774° | 16.113° | −0.339° |
| 2024 | MeanApog | 169.639° | 169.978° | −0.339° |

All four bodies show a consistent ~0.339° J2000 precession shift — as expected for a
uniform coordinate transformation.

**TrueNode vs MeanNode physical sanity:**

| Epoch | With fix | Without fix (SE behaviour) |
|-------|----------|---------------------------|
| 2024 CE | 0.49° | 0.49° |
| 0 CE | 1.04° | 28.85° |
| 3000 CE | 1.37° | 15.37° |

### The libephemeris fix

For all Pipeline-B bodies (MeanNode, MeanApog, TrueNode, OscuApog, IntpApog,
IntpPerg), when both `FLG_SIDEREAL` and `FLG_J2000` are set:

1. Compute the tropical ecliptic-of-date position.
2. Subtract **mean ayanamsha** (not true — the J2000 ecliptic frame has no nutation
   component).
3. Precess from ecliptic of date to J2000 ecliptic.
4. If equatorial output is also requested, rotate to J2000 equatorial using J2000
   obliquity.

This is the same pipeline already used correctly for `MEAN_NODE` and `MEAN_APOG`;
the fix extends it to the remaining four bodies.

**Code locations.** `libephemeris/fast_calc.py` (LEB path): removed
`_SID_J2K_SKIP_BODIES` and `_J2K_SKIP`, extended `_deferred_sid_j2k` to all
ecliptic-direct bodies. `libephemeris/planets.py` (Skyfield path): removed the
`_eff_flags = iflag & ~FLG_J2000` suppression from the TrueNode, OscuApog and
IntpApog/IntpPerg handlers so they use `iflag` directly.

**Test coverage.** `validation/compare_scripts/tests/test_sidereal/test_se_bug_j2k_nodes.py`
(J2000 applied, physical sanity, LEB-vs-Skyfield consistency, documented SE divergence
magnitude) and `validation/compare_scripts/tests/test_compare_sidereal_regression.py`
— both in the separate, non-shipped `validation/` repository (not part of a clean
checkout of this library) — and `tests/test_leb/compare/extended/test_extended_sidereal.py`
(updated to verify the intentional divergence / LEB-vs-Skyfield agreement).

### Methodology note

This analysis was performed entirely through **black-box behavioural observation**
of pyswisseph — comparing the outputs of different flag combinations across dates,
bodies and ayanamsha modes. The Swiss Ephemeris source code was not inspected. The
conclusions follow from the observed behaviour, the mathematical properties of the
coordinate transformations involved, and the internal inconsistency between mean and
true body handling.

---

## 2. Total-eclipse obscuration

**Status: divergence removed in the v3.0.0 review (round 2).** The
reference-compatible entry points (`sol_eclipse_how` `attr[2]`,
`sol_eclipse_when_loc` `attr[2]`, `sol_eclipse_where` `attr[2]`,
`lun_occult_where` `attr[2]`) now report exactly what the reference API
reports for total eclipses: the lunar/solar disc *area ratio*
`(R_moon/R_sun)² ≈ 1.05–1.12 > 1` (e.g. ~1.1167 at the Dallas maximum of
2024-04-08). An earlier release clamped this to 1.0 (the physically-bounded
covered fraction); that clamp broke 1:1 parity and has been reverted.

The convenience extensions that are **not** part of the reference API —
`sol_eclipse_obscuration_at_loc` and `sol_eclipse_how_details`
`max_obscuration` — still return the physically-bounded covered fraction
(1.0 when the Sun is fully covered), as their own documentation states.

Annular eclipses report `(R_moon/R_sun)² < 1` (the ring-residual area
fraction, identical in both conventions) and partial eclipses use the
standard two-disc lens overlap; both agree with the reference to ~1e-3. At a
no-eclipse instant the obscuration is **0.0** from every entry point
(`sol_eclipse_how`, `sol_eclipse_where`, `lun_occult_where`).

---

## 3. Of-date mean obliquity at deep-BCE epochs

For the of-date mean obliquity libephemeris uses the **angle between the Vondrák
2011 of-date ecliptic pole and equator pole** — the *same two poles* that build
the precession matrix (`erfa.ltp`/`ltpb`; `sidereal_longterm.precession_matrix`).
Deriving the obliquity from those poles makes the precession and the
equator-of-date → ecliptic-of-date rotation **one self-consistent frame**: a
direction lying in the mean ecliptic of date (the Sun, by definition) reduces to
~0 ecliptic latitude.

The reference engine instead reports Vondrák's **direct `ε_A` obliquity series**
(the published `p_A`/`ε_A` polynomial+periodic terms). That series is a genuine,
high-quality Vondrák fit — it tracks the IAU 2006 obliquity to < 1 mas near J2000,
better than the pole angle — but it is a *separate* fit from the pole series, so
it is not the angle between the pole vectors that define the of-date ecliptic.
Pairing the direct series with the pole-based precession (which both this library
and the reference use) tilts the of-date ecliptic away from the pole-defined
ecliptic by up to ~6.5″ at −3000, which surfaces as a **spurious ecliptic latitude
on the Sun** of the same size. libephemeris avoids that by taking the obliquity
from the poles; the reference does not, so its Sun sits ~6″ off its own of-date
ecliptic at −3000. (The direct series is still available in libephemeris, for
reference parity and provenance, as `sidereal_longterm.mean_obliquity_series_rad`;
it is simply not used for any equator↔ecliptic transform.)

The trade is therefore *physical correctness for bit-parity of a non-physical
value*, and the divergence is bounded, latitude-only, and identically zero in the
modern era:

| Year | Reference obliquity (`ε_A` series) | pole-angle − `ε_A` | IAU 2006 − `ε_A` | Sun latitude (lib − reference) |
|------|-----------------|-----------------|------------------|-----------------------------|
| −3000 | 24.021270° | −6.475″ | +5.718″ | −5.999″ |
| −1000 | 23.814592° | −1.040″ | +0.285″ | −1.048″ |
| 0 | 23.695022° | −0.206″ | −0.010″ | −0.208″ |
| 2000 | 23.439279° | 0.000″ | 0.000″ | 0.000″ |
| 3000 | 23.309726° | +0.048″ | +0.010″ | −0.008″ |

The pole-angle and the direct series agree to < 0.001″ across 1900–2100
(identically 0 at J2000), so modern obliquity, positions and house cusps are
unchanged. Ecliptic **longitude** is essentially unaffected (< 0.03 mas modern; the
of-date-sidereal-time / ARMC readout shifts ≤ 0.03″ at −3000), and the whole
deviation is confined to ecliptic latitude. Independently verified with
pyerfa + astropy: reducing the DE441 Sun through `erfa.ltpb` + the pole angle
gives |lat| < 1.4″ at every epoch 2000 → −4000 in both backends (vs a ramp to
~+9″ at −4000 with the direct series); astropy's DE-kernel Sun latitude at J2000
matches to 0.0001″. At −3000 the ≤ ~6″ latitude divergence is already far below
the ephemeris-generation floor on the planets (e.g. Mars ≈ 600″ at −3000).

Implementation: `sidereal_longterm.mean_obliquity_rad` returns the angle between
its own Vondrák `_ecliptic_pole`/`_equator_pole` vectors (identical to
`erfa.ltpecl`/`erfa.ltpequ` to < 1e-4 mas); `precession_vondrak` and every
equator↔ecliptic-of-date path (both backends, plus the house cusps) share this one
realization.

---

## 4. `nod_aps` MEAN speed slots

For the MEAN node/apsis method the reference engine fills the latitude-speed
slot (`xperi[4]` / `xaphe[4]`) with the **latitude evaluated one day later**,
not with a rate: probing shows `lat_spd == lat(t + 1 day)` to full printed
precision at every epoch, while the longitude and distance slots contain a
genuine forward difference over one day, `(x(t+1d) − x(t)) / 1d` (both
verified black-box, e.g. Moon MEAN perigee at J2000: reference
`lat = −3.419715`, `lat_spd = −3.408674 = lat(t+1d)`).

libephemeris reproduces the one-day forward-difference convention for the
longitude and distance speeds (exact match for the Moon), but returns the
true forward-difference **rate** for the latitude speed — a value like
−3.4°/day for a point whose latitude changes by 0.011°/day is a defect, not
a convention. Consumers of the reference's MEAN `lat_spd` are reading a
position, so no astrological result can depend on it being reproduced.

For the MEAN apsides of the *planets*, libephemeris also keeps the true
derivative semantics for the longitude and distance speed channels. This was
checked independently by deriving perihelion/aphelion positions with
`FLG_SPEED` disabled and applying central step halving plus Richardson
extrapolation (`validation/verify/verify_nodaps_mean_speed_derivative.py`):
Mercury, Venus and Mars at 1950, J2000 and 2025 agree with the numerical
derivative to <3e-7°/day for longitude speed and <4e-9 AU/day for distance
speed. The reference engine differs by up to ~0.15°/day on those channels.

The certified divergence therefore applies only to MEAN speed channels. The
node/apside positions themselves are not covered by this certification.

## 5. `deltat_ex` has no side effect on `get_tid_acc`

In the reference API, calling `deltat_ex(jd, FLG_MOSEPH)` in automatic mode
*mutates* the stored tidal acceleration (a later `get_tid_acc()` returns
−25.58). libephemeris computes the same ΔT value (the tid-acc adjustment is
applied functionally) but keeps `deltat_ex` pure: `get_tid_acc()` only
changes when `set_tid_acc()` is called. Only introspection via
`get_tid_acc()` after `deltat_ex` can observe the difference.

## 6. `nod_aps` node positions: transient light-time spikes

With default flags the reference engine applies its apparent-position
correction to each node/apsis point with a per-point numerical scheme that
occasionally glitches: on isolated days (~15 days out of 3653 sampled over
1995-2005 for Jupiter) one node — and only one — jumps by up to ~3.5" while
its reported speed spikes an order of magnitude, and with `FLG_TRUEPOS` the
same call returns exactly antipodal nodes again (0.000000"). The median
default-flag deviation is 0.015" (p99 ≈ 0.43").

libephemeris applies its aberration correction analytically, so its nodes
are exactly antipodal at all times and agree with the reference to ≤0.03"
outside those spike days. The spikes are not reproduced. Validation
comparisons of node positions should either use `FLG_TRUEPOS` or tolerate
the reference's p99 (~0.5").

## 7. Fixed-star distance and distance-speed channels

For fixed stars, libephemeris returns `dist` as the apparent geocentric distance
of the same star state used for the angular coordinates, and `dist_spd` as the
central finite-difference derivative of that distance with a half-day step. This
keeps the fixed-star channel semantics consistent with the planetary API: a
speed slot is the derivative of the corresponding position slot.

The reference engine uses a different radial-distance convention for some
catalog entries. Reproducing that convention would make the same `dist_spd`
channel mean one thing for planets and another for fixed stars, so it is not
copied.

Independent verification (`validation/verify/verify_star_distance_astropy.py`)
uses Astropy/SkyCoord space motion plus the geocentric Earth vector as the
oracle. Across Spica, Sirius, Regulus, Aldebaran, Vega, Arcturus, Rigil Kent,
Altair, Antares and one zero-parallax catalog entry at 1950, J2000 and 2025,
libephemeris agrees with the Astropy distance to <1.3e-7 relative and with the
distance-speed derivative to <1.7e-5 AU/day.

## 8. `EQUATORIAL | J2000` for hypothetical bodies (RESOLVED — parity, not a divergence)

**Status:** resolved during the v3 speed/frame round; verified parity.

Before v3, libephemeris rotated the hypothetical orbital bodies
(Cupido…Poseidon, Transpluto/Isis) from J2000 ecliptic to the equator with
the **true obliquity of date** under `FLG_EQUATORIAL | FLG_J2000` — a
frame-mixed output ~3" away from a consistent J2000 equatorial frame in
2023 (growing with distance from J2000), and inconsistent with every other
body class. The v3 round fixed this to use the **J2000 obliquity**, like
all other bodies, in both backends.

An earlier revision of this section certified the fix as a divergence,
attributing the frame-mixed rotation to the reference API. Direct
measurement against the reference oracle later disproved that attribution:
the reference's reported `EQ|J2000` equals rotating its own J2000-ecliptic
output with the **J2000 obliquity** for all bodies including the
hypotheticals (Cupido 0.006", Isis 0.008" from the eps_J2000 rotation at
2005, vs 0.6–0.9" from the of-date one; the gap grows to >100" by 2300).
The pre-v3 mixed frame was therefore a plain libephemeris bug, and the fix
**restored** parity. `EQ|J2000` output now matches the reference for every
body class; nothing is intentionally divergent here.

## 9. `SIDEREAL | EQUATORIAL` speed frame (RESOLVED — parity, not a divergence)

**Status:** resolved during the v3 speed/frame round; verified parity.

For `FLG_SIDEREAL | FLG_EQUATORIAL`, both engines report the **position** on
the *mean* equator of date (precession only, no nutation), with no ayanamsha
subtraction from RA. libephemeris computes the accompanying RA/Dec **speed**
in the mean-equator frame of the reported position, in both backends, so
every speed slot is the exact time-derivative of the corresponding position
slot (residual ≤0.034"/day vs a central derivative for the Moon).

An earlier revision of this section certified this as a divergence, citing
an old in-code note that the reference returned the *true*-equator rate
(~0.8"/day RA / ~1.8"/day Dec away from the derivative for the Moon).
Direct measurement against the reference oracle later disproved that note:
the reference's reported SID|EQ Moon speed (dRA 45197.42, dDec 19297.82
"/day at JD 2460000.5, Lahiri) sits on the derivative-of-reported-position
side and matches libephemeris to **0.05"/day** — not the true-equator rate.
The two engines agree; nothing is intentionally divergent here.

## 10. Sripati `house_pos` in the house-12 wrap region

**Status:** intentional divergence (documented after the v3 review round).

`house_pos` with the Sripati system (`hsys = 'S'`) places each body between
1.0 and 13.0 by shifting the Porphyry result forward by half a house, so a body
deep in house 12 lands in `[12.0, 13.0)` before the wrap back to `[1.0, 1.5)`.
libephemeris wraps only the true overflow (`if hpos >= 13.0: hpos -= 12.0`),
**preserving the fraction** so the house-12 region reports values like 12.03,
12.27, … (internally consistent with the engine's own Sripati cusps, which put
those longitudes in house 12).

The reference API instead collapses the whole house-12 wrap region onto house 1:
for ~43 of 360 body longitudes (the entire Sripati house-12 region plus the
first half of house 1) it returns exactly `1.0`, dropping the fraction. This is
internally inconsistent with the reference's *own* Sripati cusps — the same
longitudes fall between its cusp 12 and cusp 1, yet its `house_pos` reports
house 1. Verified by direct differential measurement at armc 100°, lat 51.5°
(`lib house_pos = 12.03` vs `reference = 1.00` at bodylon 144°; the pattern
reproduces at every lat/armc tested).

libephemeris keeps the geometrically-consistent, fraction-preserving value
because the reference's collapse loses information (a body in house 12 is
reported as house 1) and contradicts the cusps it returns for the same chart.
Consumers that need bit-parity with the reference's collapse for the house-12
region should apply `1.0 if hpos > 12.0 else hpos` to the result. Every other
house region agrees exactly; only the Sripati house-12 wrap region differs.
