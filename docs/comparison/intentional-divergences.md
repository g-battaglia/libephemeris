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

**Test coverage.** `tests/test_sidereal/test_se_bug_j2k_nodes.py` (J2000 applied,
physical sanity, LEB-vs-Skyfield consistency, documented SE divergence magnitude);
`compare_scripts/tests/test_compare_sidereal_regression.py` and
`tests/test_leb/compare/extended/test_extended_sidereal.py` (updated to verify the
intentional divergence / LEB-vs-Skyfield agreement).

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

For the of-date mean obliquity, libephemeris uses the true angle between the
Vondrák-2011 of-date equator and ecliptic poles (valid ±200,000 years), rather than
the IAU 2006 obliquity polynomial (valid only a few centuries from J2000). At
deep-BCE epochs Swiss Ephemeris reports an obliquity that matches neither model
exactly — it sits between the rigorous Vondrák pole angle and the IAU 2006
extrapolation. libephemeris keeps the physically-consistent Vondrák value rather
than reproducing Swiss's value bit-for-bit. The deviation is bounded and benign,
and — because obliquity does not affect ecliptic longitude — confined to ecliptic
latitude:

| Year | Swiss obliquity | Vondrák − Swiss | IAU 2006 − Swiss | Sun latitude (lib − Swiss) |
|------|-----------------|-----------------|------------------|-----------------------------|
| −3000 | 24.021270° | −6.475″ | +5.718″ | −5.999″ |
| −1000 | 23.814592° | −1.040″ | +0.285″ | −1.048″ |
| 0 | 23.695022° | −0.206″ | −0.010″ | −0.208″ |
| 2000 | 23.439279° | 0.000″ | 0.000″ | 0.000″ |
| 3000 | 23.309726° | +0.048″ | +0.010″ | −0.008″ |

The deviation is identically zero in the modern era, sub-arcsecond within recorded
history, and at −3000 (≤ ~6″ in latitude) already well below the
ephemeris-generation floor on the planets (e.g. Mars ≈ 600″ at −3000). The of-date
mean obliquity is taken from ERFA's Vondrák long-term routines
(`erfa.ltpequ` / `erfa.ltpecl`).

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

## 8. `EQUATORIAL | J2000` for hypothetical bodies (Uranians, Transpluto)

**Status:** intentional divergence since the v3 speed/frame round.

For the hypothetical orbital bodies (Cupido…Poseidon, Transpluto/Isis), the
reference API computes `FLG_EQUATORIAL | FLG_J2000` by rotating the body's
J2000 *ecliptic* coordinates to the equator with the **true obliquity of
date** — a frame-mixed output (J2000 ecliptic longitude zero-point, of-date
equator tilt) that is ~3" away from a consistent J2000 equatorial frame in
2023 (Cupido 3.1", Transpluto 1.4"; the offset grows with distance from
J2000). Every other body class (planets, nodes, apsides, Lilith, White Moon,
fixed stars) already gets a consistent J2000 frame: ecliptic → equator with
the J2000 obliquity.

libephemeris rotates these bodies with the **J2000 obliquity** as well, in
both backends (Skyfield path in `planets.py`, LEB fast path in
`fast_calc.py` Pipeline C), so `EQ|J2000` means the same frame for every
body. RA/Dec speeds are carried through the same fixed rotation, so the
reported speed remains the derivative of the reported position.

Verification: rotating the body's own `FLG_J2000` ecliptic output to the
equator with eps_J2000 = 23.4392911° reproduces the returned `EQ|J2000`
RA/Dec to <0.001" for Cupido, Isis, Mars, White Moon and True Node alike
(scratch check `uranian_eqj2000.py`, JD 2460000.5, both backends). Only the
combination `FLG_EQUATORIAL | FLG_J2000` on hypothetical bodies is affected;
of-date equatorial, ecliptic and J2000-ecliptic outputs are unchanged.

## 9. `SIDEREAL | EQUATORIAL` speed frame (mean equator)

**Status:** intentional divergence since the v3 speed/frame round.

For `FLG_SIDEREAL | FLG_EQUATORIAL`, both the reference API and libephemeris
report the **position** on the *mean* equator of date (precession only, no
nutation), with no ayanamsha subtraction from RA. The reference API, however,
computes the accompanying RA/Dec **speed** in the *true*-equator (plain
equatorial) frame — the same rate it returns without `FLG_SIDEREAL`. The
result is a speed that does not differentiate the reported position: for the
Moon the reference's SID|EQ rate differs from the derivative of its own
SID|EQ position by ~0.8"/day in RA and ~1.8"/day in Dec (the mean-vs-true
equator rate gap; smaller for slower bodies).

libephemeris computes the SID|EQ speed in the **mean-equator frame of the
reported position**, in both backends: the LEB fast path differentiates the
apparent place through the same mean-equator (P-matrix) transform used for
the position, and the Skyfield path keeps `FLG_SIDEREAL` set on its
central-difference samples so they are taken in the reported frame. This is
the project's certified true-rate convention — every speed slot is the exact
time-derivative of the corresponding position slot.

Numerical evidence (Moon, JD 2460000.5, Lahiri): reported dRA/dDec vs the
central derivative of the reported RA/Dec agree to 0.016"/day (LEB) and
0.034"/day (Skyfield backend), while re-computing the speed in the
true-equator frame — the reference convention — would be off by −0.760"/day
in RA and +1.843"/day in Dec from that derivative. Only the speed slots of
the `FLG_SIDEREAL | FLG_EQUATORIAL` combination are affected; positions are
identical to the reference convention, and plain-equatorial output (true
equator, no sidereal flag) is unchanged.
