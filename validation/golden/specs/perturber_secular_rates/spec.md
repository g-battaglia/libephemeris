# Specification — perturbing-planet secular rates (B-11)

**Status: revised design contract; implementation not included.** This document
specifies the project-owned perturbing-planet element and frame contract. It does
not change runtime code, golden files, checkers, or allowances. The algebra and
frame portions may be reviewed independently, but **B-11 closure remains
BLOCKED** until the separately reviewed numerical oracle in §9 is supplied.

## 1. Scope and scientific decision

The local minor-body model uses Jupiter, Saturn, Uranus, and Neptune as
perturbers. This ticket replaces the undocumented first-order planet-element
rates used by the forced eccentricity and inclination vectors. It does not
replace the local test-particle secular model, its project-defined resonance
heuristics, or any SPK, Horizons, or direct-vector backend.

The one coefficient system is the project-owned transcription in
`libephemeris/planetary_mean_elements.py`: Simon et al. (1994) mean planetary
elements as represented by Meeus (1998), Table 31.A. The source frame is the
**mean ecliptic and mean equinox of date**. The source time argument is

```text
T = (jd_tt - 2451545.0) / 36525.0
```

where `jd_tt` is Julian Date in TT and `T` is Julian centuries from J2000.0 TT.
The target minor-body orbital-element vectors are in the **mean ecliptic of
J2000**. Those two frames must not be mixed as scalar angles.

The implementation shall evaluate the **full registered polynomial** for every
consumed field `a`, `e`, `i`, `Omega`, and `varpi` at the requested `T`. It shall
not replace the polynomial by its J2000 intercept plus first derivative. In
particular, all nonzero polynomial terms in `a` are applied. The argument of
perihelion is derived from the two evaluated longitudes:

```text
omega_raw(T) = varpi(T) - Omega(T)
omega(T)     = norm360(omega_raw(T))
norm360(x)  = x - 360 * floor(x / 360)       # [0, 360)
```

When a derivative is required for a diagnostic, it is the derivative of the
full polynomial, not a constant rate:

```text
dx/dT = c1 + 2*c2*T + 3*c3*T**2
domega/dT = dvarpi/dT - dOmega/dT
```

No mean-longitude coefficient table is part of this contract. Mean longitude
and the existing Neptune resonance/libration path are explicitly out of scope;
that path must not consume a newly invented `L` field or rate.

## 2. Registered full coefficient contract

Each row below is an ascending-power coefficient tuple
`(c0, c1, c2, c3)` for `x(T) = sum(c[k] * T**k)`. A shorter tuple has exactly
zero coefficients for omitted higher powers. Units are:

* `a`: AU;
* `e`: dimensionless;
* `i`, `Omega`, and `varpi`: degrees.

The corresponding non-constant coefficients per `T` are in AU/cy for `a`,
dimensionless/cy for `e`, and deg/cy for angles. These are the four
project-registered rows, copied here to make the implementation contract
executable without a hidden linearization.

| planet | `a(T)` | `e(T)` | `i(T)` | `Omega(T)` | `varpi(T)` |
| --- | --- | --- | --- | --- | --- |
| Jupiter | `(5.202603191, 0.0000001913)` | `(0.04849485, 0.000163244, -0.0000004719, -0.00000000197)` | `(1.303270, -0.0054966, 0.00000465, -0.000000004)` | `(100.464441, 1.0209550, 0.00040117, 0.000000569)` | `(14.331309, 1.6126668, 0.00103127, -0.000004569)` |
| Saturn | `(9.554909596, -0.0000021389)` | `(0.05550862, -0.000346818, -0.0000006456, 0.00000000338)` | `(2.488878, -0.0037363, -0.00001516, 0.000000089)` | `(113.665524, 0.8770979, -0.00012067, -0.000002380)` | `(93.056787, 1.9637694, 0.00083757, 0.000004899)` |
| Uranus | `(19.218446062, -0.0000000372, 0.00000000098)` | `(0.04629590, -0.000027337, 0.0000000790, 0.00000000025)` | `(0.773196, 0.0007744, 0.00003749, -0.000000092)` | `(74.005947, 0.5211258, 0.00133982, 0.000018516)` | `(173.005159, 1.4863784, 0.00021450, 0.000000433)` |
| Neptune | `(30.110386869, -0.0000001663, 0.00000000069)` | `(0.00898809, 0.000006408, -0.0000000008, -0.00000000005)` | `(1.769952, -0.0093082, -0.00000708, 0.000000028)` | `(131.784057, 1.1022057, 0.00026006, -0.000000636)` | `(48.123691, 1.4262677, 0.00037918, -0.000000003)` |

At `T = 0`, the derived argument values and first derivatives are:

| planet | `omega(0)` deg | `domega/dT(0)` deg/cy |
| --- | ---: | ---: |
| Jupiter | `273.8668680` | `+0.5917118` |
| Saturn | `339.3912630` | `+1.0866715` |
| Uranus | `98.9992120` | `+0.9652526` |
| Neptune | `276.3396340` | `+0.3240620` |

The Neptune value is the exact subtraction of the registered linear terms,
`1.4262677 - 1.1022057 = +0.3240620 deg/cy`. It is not `+0.3240613`.
Normalization applies only to a returned angle; it never modulo-reduces a
signed derivative.

The table is intentionally limited to `a`, `e`, `i`, `Omega`, and `varpi`.
There are no `L` rows, aliases, claims of `L` registration, or resonance
coefficients in this specification. The existing
`NEPTUNE_MEAN_LONGITUDE_J2000`/`NEPTUNE_N` resonance implementation remains a
separate, unchanged contract and is not a consumer of this table.

## 3. Evaluator and tuple semantics

`mean_orbital_elements()` is the project-owned evaluator and table authority.
The implementation shall use that same table/evaluation rule for the
perturbers, rather than maintaining a second rounded static list. In particular:

* `jd_tt=None` in `_calc_forced_elements()` means exactly
  `jd_tt = 2451545.0`, hence exactly `T = 0`;
* the `None` and explicit-J2000 branches call/use the same evaluator and table;
* no rounded J2000 compatibility row survives as a separate implementation
  branch;
* `a(T)` is passed through with its nonzero terms, including the outer-planet
  terms shown in §2;
* `mu` remains the existing separately sourced mass ratio, and the Neptune
  `a_threshold` remains `20.0` AU. Those values are not inferred from the
  polynomial table.

The existing seven-tuple shape is retained so callers do not acquire a new
positional API:

```text
(a_planet_au, mu, e, omega_deg_j2000, node_deg_j2000,
 i_deg_j2000, a_threshold_au)
```

Its angle fields have the following exact meanings:

* `omega_deg_j2000` is the argument of perihelion in the mean ecliptic J2000;
* `node_deg_j2000` is the longitude of the ascending node in that same frame;
* `i_deg_j2000` is inclination to that same ecliptic.

They are not source-frame scalars copied from the date polynomial. They are
obtained from the single vector transformation in §4. If an implementation
uses an internal typed record instead of this tuple, it must carry the same
fields and the two transformed vectors; it must never carry a date-frame
`omega`/`Omega` pair alongside J2000 vectors.

`_calc_forced_elements()` then retains its existing vector semantics:

```text
varpi_j2000 = omega_deg_j2000 + node_deg_j2000
h_j, k_j    = e_j * sin(varpi_j2000), e_j * cos(varpi_j2000)
p_j, q_j    = sin(i_j/2) * sin(node_j), sin(i_j/2) * cos(node_j)
```

The synthetic-angle tests must make `omega + node` materially different from
`omega`, so using argument of perihelion as longitude of perihelion cannot pass
accidentally. The `OrbitalElements.omega` field and
`rebound_integration.elements_to_rebound_orbit()` remain argument-of-perihelion
consumers and are regression guards.

## 4. Independent date-ecliptic to J2000-ecliptic frame contract

### 4.1 Vector construction in the source frame

For each perturber, first evaluate the five source elements at `T` and form
vectors in a right-handed date-ecliptic Cartesian frame whose axes are:

* `x`: the mean equinox of date;
* `y`: +90 degrees of mean ecliptic longitude from that equinox;
* `z`: the north pole of the mean ecliptic of date.

Use the evaluated `varpi` and `Omega`, not independently propagated scalar
rates:

```text
evec_date = e * (cos(varpi), sin(varpi), 0)
pole_date = (sin(i)*sin(Omega),
             -sin(i)*cos(Omega),
              cos(i))
```

`evec_date` is the eccentricity vector. `pole_date` is the unit orbital pole
for the standard perifocal-to-ecliptic rotation used by the project. The
inclination vector `(p, q)` is derived only after the pole has reached the
J2000 frame; it is not rotated as two unrelated scalar angles.

### 4.2 Matrix direction, order, and convention

All vectors are **column vectors** and all displayed matrices act on the left.
The project passive x-axis rotation is

```text
Rx(eps) = [[1, 0,       0      ],
           [0, cos(eps), sin(eps)],
           [0, -sin(eps), cos(eps)]]
```

Under this convention, ecliptic-to-equatorial is `Rx(-eps)` and
equatorial-to-ecliptic is `Rx(+eps)`.

Let:

* `eps_date = vondrak_mean_obliquity_rad(jd_tt)`, the project's registered
  Vondrak 2011 mean obliquity for the date;
* `eps_j2000` be the registered J2000 mean obliquity used by the existing
  ecliptic-J2000 state paths (`23.4392911` degrees, with no nutation);
* `P_date = vondrak_precession_matrix(jd_tt, frame_bias=False)`, whose rows map
  a mean-equatorial J2000 vector to the mean equator/equinox of date, as
  documented in `precession_vondrak.py`.

The required date-ecliptic to J2000-ecliptic matrix is

```text
C_date_to_j2000 = Rx(+eps_j2000) @ transpose(P_date) @ Rx(-eps_date)
```

and it is applied as `v_j2000 = C_date_to_j2000 @ v_date`. The transpose is
required because `P_date` maps J2000 to date and is an orthogonal rotation.
The reverse witness is

```text
C_j2000_to_date = Rx(+eps_date) @ P_date @ Rx(-eps_j2000)
```

so `C_j2000_to_date @ C_date_to_j2000` must be identity to binary64 rounding.
There is no nutation matrix, aberration, light-time correction, or longitude
subtraction in this transformation.

The J2000 target obliquity is stated explicitly because the minor-body J2000
Cartesian paths in `spk.py`, `rebound_integration.py`, and `fast_calc.py` use
the registered IAU J2000 value above. A future migration of all those paths to
the Vondrak pole-angle value would be a separate frame contract, not an
implicit part of B-11.

### 4.3 Recovering the tuple from transformed vectors

Apply `C_date_to_j2000` to both `evec_date` and `pole_date`. Then derive the
canonical tuple fields from those vectors:

```text
e_j2000       = norm(evec_j2000)
varpi_j2000   = norm360(deg(atan2(evec_y, evec_x)))
i_j2000       = deg(acos(clamp(pole_z, -1, 1)))
node_j2000    = norm360(deg(atan2(pole_x, -pole_y)))
omega_j2000  = norm360(varpi_j2000 - node_j2000)
```

The normal vector is unit-normalized before the `acos`/`atan2` operations if
rounding has changed its norm. For a nonzero source inclination, its node is
well-defined. A zero-inclination pole has no defined node; that case is not
present among these four perturbers and must be handled by the project's
existing zero-inclination convention rather than by inventing a date/J2000
angle mix.

This construction makes all scalar angles consumed by `_calc_forced_elements`
J2000 angles, while preserving the source-frame polynomial as the sole source
of the perturbing vectors. Merely subtracting a precession longitude is not an
implementation of this contract.

## 5. Source domain and runtime policy

`planetary_mean_elements.py` and the current provenance registry identify the
Simon et al. (1994)/Meeus Table 31.A coefficient source, TT epoch, Julian-century
unit, angle units, and date-ecliptic frame. They do **not** register validated
start and end dates for the full mean-element expressions. This specification
therefore does not invent endpoints or an accuracy claim. The exact source
validity interval must be recorded from the reviewed publication before any
implementation review claims source-domain accuracy.

For every finite input accepted by the existing public path, runtime
compatibility remains defined and the evaluator returns the full polynomial
value. Outside the source's eventually recorded validity interval, the result
shall be labelled **model extrapolation**. Full-polynomial evaluation removes
the additional algebraic error caused by a J2000 linear truncation; it does not
validate the polynomials physically in deep time. No hidden date refusal, new
exception, or silent backend switch is authorized by B-11.

The existing SPK, Horizons, and optional ASSIST/REBOUND dispatch policy remains
in force. When one of those paths supplies a state, B-11's local perturbing
coefficient evaluator is not part of that state calculation.

## 6. Affected and unaffected surfaces

### Affected local-model surfaces

The implementation review shall inventory and test:

* `libephemeris/planetary_mean_elements.py`, including its full table and
  evaluator contract;
* `minor_bodies._get_planet_elements_at_time()` and
  `minor_bodies._calc_forced_elements()`;
* `minor_bodies.apply_secular_perturbations()` and the local Keplerian position
  path `calc_minor_body_position()`;
* `minor_bodies.calc_minor_body_heliocentric()` when it falls through to the
  local model;
* `planets._keplerian_position_at()` and `_calc_keplerian_fallback()`, plus
  public `calc_ut()` callers that select that fallback;
* `rebound_integration.elements_to_rebound_orbit()` and its argument/node
  semantics; and
* the local-model portions of element, phenomenon, node/apsis, and
  distance-extrema paths when they reach the minor-body fallback.

The existing test-particle secular-rate formulas (`d_omega`, `d_Omega`, and the
project-defined `d_n`) are not silently re-sourced by this specification. A
change to those formulas requires its own reviewed contract. The four
perturber `a(T)` values are used where the forced-vector coefficients require a
perturber semi-major axis.

### Unaffected backend guards

The following must not start calling the new evaluator merely because B-11 is
implemented:

* successful SPK fast paths;
* successful Horizons/direct-vector/Skyfield state paths;
* propagation with `include_perturbations=False`;
* fixed-star, lunar, and major-planet state paths; and
* the public meaning of `OrbitalElements.omega`, `OrbitalElements.Omega`, and
  the REBOUND tuple conversion.

`calc_minor_body_position()` exposes position `(x, y, z)` in AU in ecliptic
J2000. `calc_minor_body_heliocentric()` exposes `(longitude, latitude,
distance)` and no public velocity. This contract makes no speed or velocity
accuracy claim. Although `calc_ut()` may expose speed slots through its own
finite-difference machinery and `PropagationResult` has AU/day velocity fields,
those are separate surfaces and are not B-11 acceptance components.

## 7. Deterministic implementation guards and semantic tests

The implementation test set shall include all of the following deterministic
guards. A monkeypatch sentinel must raise if the full polynomial evaluator is
called unexpectedly.

1. **Disabled perturbations.** `apply_secular_perturbations(...,
   include_perturbations=False)` returns before the evaluator, forced-vector
   calculation, and local secular-rate calculation. The corresponding
   `calc_minor_body_position(..., include_perturbations=False)` path also makes
   no evaluator call, including its short-period branch.
2. **`None` identity.** `_calc_forced_elements(elements, jd_tt=None)` and
   `_calc_forced_elements(elements, jd_tt=2451545.0)` use the same evaluator
   and produce identical six-tuples. A call sentinel records the exact J2000
   JD used by both branches.
3. **SPK no-call.** With a registered successful SPK fixture, call
   `calc_minor_body_heliocentric(..., use_spk=True)` while the evaluator is a
   raising sentinel. The SPK result must be returned without touching the
   evaluator. A deliberately unsuccessful SPK fixture is a different test and
   may fall through to the local model.
4. **Horizons/direct-vector no-call.** Force a successful Horizons/direct
   vector state dispatch in the existing planet/backend harness and repeat the
   raising-sentinel test. The state path must not import or evaluate the local
   mean-element table. The test covers both the Horizons route and the direct
   vector route where those fixtures are available; it must not use a failed
   backend as evidence of no-call behavior.
5. **Tuple semantics.** Synthetic nonzero `omega`, `Omega`, `i`, and `e` values
   verify that the forced eccentricity vector uses `omega + Omega`, while the
   inclination vector uses the node. Separate synthetic date-frame vectors
   verify that the tuple angles are recovered only after the matrix conversion.
6. **Frame round-trip.** For arbitrary nonzero vectors and every engineering
   date stratum in §9, independently compute date-ecliptic -> J2000-ecliptic ->
   date-ecliptic. Check matrix orthogonality, matrix direction, preservation of
   vector norms, and component recovery to a binary64 conditioning bound.
7. **Full-polynomial identity.** At positive and negative nonzero `T`, assert
   that every returned field equals the complete registered coefficient
   polynomial, including the nonzero `a` terms. Check an angular wrap at a
   360-degree boundary and the exact Neptune `+0.3240620 deg/cy` J2000
   derivative.

## 8. Implementation acceptance versus scientific closure

### 8.1 Algebra/frame implementation acceptance

A role-C implementation may be approved for algebra and frame correctness when
it demonstrates, without compatibility outputs:

* coefficient identity for all five fields and all four perturbers;
* full-polynomial evaluation and full-polynomial derivative diagnostics;
* exact `jd_tt=None`/J2000 identity and application of nonzero `da/dT` terms;
* the derived argument-of-perihelion rule and corrected Neptune derivative;
* the tuple/vector semantics in §§3–4;
* date-ecliptic/J2000 matrix direction, order, no-nutation policy, and round-trip;
* all no-call and disabled-perturbation guards in §7; and
* provenance references that identify the table locator, TT conversion, units,
  source frame, Vondrak/ERFA matrix convention, and the unregistered source
  domain rather than claiming one.

This is an internal mathematical contract. Passing it does not establish that
the local secular model predicts physical minor-body positions.

### 8.2 Why B-11 closure is still BLOCKED

The current local secular model is a reduced, project-authored approximation.
Its Laplace-coefficient quadrature, Hill-sphere skip, high-eccentricity and
high-inclination corrections, resonance heuristics, and mean-motion correction
are not a complete published dynamical solution with a registered output
validity/accuracy envelope. The project therefore cannot accept changed public
golden positions, longitudes, latitudes, or distances from an observed maximum,
compatibility output, or a permissive blanket Level-2 tolerance.

B-11 closure remains **BLOCKED** until a separately reviewed numerical oracle
certifies at least the internal perturbing vectors and explicitly decides
whether a scientifically sourced public-state gate exists. If that review
finds no physical source/validity contract for the current secular positions,
acceptance is limited to internal coefficient, frame, tuple, and forced-vector
correctness. Public golden differences remain blocked; no Level-2 allowance is
proposed here.

## 9. Required independent oracle protocol

The following protocol defines the pending oracle closely enough for a separate
review. It is an independent calculation, not a persisted output comparison.

### 9.1 Independent evaluator and frame calculation

The oracle shall:

1. implement the five coefficient polynomials independently (for example with
   `math.fsum(c[k] * T**k for k in range(len(c)))`, not by calling the project
   evaluator), using the coefficient table in §2;
2. derive `omega` from independently evaluated `varpi` and `Omega`;
3. construct `evec_date` and `pole_date` independently;
4. obtain the registered Vondrak precession matrix directly from ERFA's
   `ltp()` semantics (or an independently coded equivalent of the published
   matrix), use the independently coded `Rx` matrices, and apply exactly the
   column-vector order in §4; and
5. independently recover J2000 `e`, `omega`, `Omega`, and `i` from the transformed
   vectors before comparing the candidate tuple.

The oracle must not call the candidate helper, fit a correction, use a public
compatibility result, or persist per-date outputs.

### 9.2 Exact finite grid

The grid is deliberately an engineering grid, not a claim about source
validity. Every date is TT:

```text
JD(T) = 2451545.0 + 36525.0*T

S0: T = {0}
S1: T = {-1, -0.25, 0.25, 1}
S2: T = {-10, -5, 5, 10}
S3: T = {-100, 100}       # model-extrapolation conditioning only
```

The local-body set is every named entry currently in
`MINOR_BODY_ELEMENTS`, exactly (the symbolic names below are the names used by
the constants module):

```text
CHIRON, PHOLUS, CERES, PALLAS, JUNO, VESTA,
ERIS, SEDNA, HAUMEA, MAKEMAKE, IXION, ORCUS, QUAOAR,
NESSUS, ASBOLUS, CHARIKLO, GONGGONG, VARUNA,
APOPHIS, HYGIEA, INTERAMNIA, DAVIDA, EUROPA_AST, SYLVIA,
PSYCHE, EROS, AMOR, ICARUS, TORO, SAPPHO, PANDORA_AST,
LILITH_AST, HIDALGO, TOUTATIS, ITOKAWA, BENNU, RYUGU
```

This is 37 entries. The oracle must enumerate the actual current table rather
than silently omitting a newly registered entry; the list above is the explicit
review snapshot for this specification.

The oracle runs the direct perturbing-table and forced-vector checks for all
listed bodies and all strata. It additionally includes one synthetic body with
`a` just outside each of the four perturber Hill-radius guards, one synthetic
body with zero eccentricity, and one with zero inclination, to exercise the
condition and convention branches. Resonance/libration coefficients and mean
longitude are not part of this oracle.

For each body/date, the exact request variants are:

* `_calc_forced_elements(jd_tt=None)` and explicit `jd_tt=JD(0)`;
* explicit `jd_tt=JD(T)` for every `T` in `S0`–`S3`;
* `apply_secular_perturbations` with
  `include_perturbations=True` and `False` for `S0`, `S1`, and `S2`;
* `calc_minor_body_position` with perturbations enabled and disabled for the
  same body/date subset; and
* successful SPK and successful Horizons/direct-vector dispatch fixtures for
  the no-call guards.

The last three variants are control-flow tests only. They do not create a
public physical-accuracy gate. The exact compared components are:

* source scalars: `a`, `e`, `i`, `Omega`, `varpi`, derived `omega`;
* transformed vectors: all three components of `evec_j2000` and `pole_j2000`;
* forced output: `g`, `h_forced`, `k_forced`, `s`, `p_forced`, `q_forced`;
* local position controls: `x`, `y`, `z` in AU, and the separately exposed
  longitude, latitude, and distance; and
* no velocity component. A speed/velocity result is not asserted because the
  public minor-body helper has no velocity API in this contract.

### 9.3 Conditioning-derived comparison bounds

The oracle shall derive, record, and apply a bound for each value from binary64
rounding and the declared algorithm, never from an observed maximum. Let
`u = 2**-53` and

```text
gamma(n) = (n*u) / (1 - n*u)
```

for an operation path with `n` rounded operations. For a degree-`d` polynomial,
the scalar bound is computed from the coefficient magnitude and condition
number, using a two-evaluation bound such as

```text
B_poly = 2 * gamma(2*d) * sum(abs(c[k] * T**k) for k in range(d + 1))
          + 2*u*abs(value)
```

with the actual counted operation path recorded by the review. The independent
power-sum evaluator may use a larger declared `n`; it must use that larger
bound rather than an observed residual.

For a three-term matrix dot product, use the corresponding `gamma(3)` dot
product bound, propagate the input scalar bounds through the absolute matrix,
and add the bound for the three output rows. The matrices are orthogonal, so
there is no physical amplification beyond the explicitly calculated component
conditioning. For vector normalization, `acos`, `atan2`, and division, propagate
the Jacobian/condition number at the actual vector; if a norm or a forced
frequency is below the documented zero-branch threshold (`1e-20` for `g`/`s`),
assert the exact zero-branch behavior and do not claim an angular bound for an
undefined direction.

For the independent forced-vector formula, use the same declared finite
quadrature (200 panels, or 500 when `alpha > 0.7`), explicit Hill/threshold
guards, and explicit high-`e`/high-`i` correction branches as the current
project equation, but implement the summation independently. The review shall
count quadrature, numerator, denominator, and division operations and propagate
their `gamma(n)` and condition-number bounds. A discrepancy beyond that
conditioning-derived bound is an implementation failure; no empirical
allowance may replace it.

These bounds certify agreement between two implementations of the stated
algebra. They do not certify an N-body truth model. A future physical oracle,
if separately sourced and reviewed, must declare its initial-state convention,
frame, integration tolerances, covered dates, and component-specific bounds
before it can unblock public golden differences.

## 10. Golden families, targeted commands, and provenance work

Potentially changed public families are `pheno_nodaps` (especially
`elements_minor`, `nodaps_*`, `pheno_*`, and extrema paths) and local-model
`positions`. No changed cell or allowance is authorized by this specification;
public golden review stays blocked as stated in §8.2.

Run only targeted commands during implementation:

```text
pytest tests/test_minor_bodies/test_secular_perturbations.py -q
pytest tests/test_minor_bodies/test_cov100_minor_bodies.py -q
pytest tests/test_secular_epoch_continuity.py -q
pytest tests/test_precession_vondrak.py -q
pytest validation/compare_scripts/tests/test_compare_minor_bodies.py -q
(cd validation && ./run.sh golden check --families pheno_nodaps,positions --sha v3.2.0 --jobs 6)
uv run python scripts/check_provenance.py
```

Do not run the full project test suite for this ticket. The implementation
commit, if later authorized, must update the existing `minor-body-runtime`
provenance entry and the relevant methodology text with the exact coefficient
table locator, TT conversion, units, date-ecliptic source frame, the matrix
source/order, and the explicit unregistered-domain/extrapolation statement.
This specification commit itself adds no runtime, checker, allowance, or
provenance-registry change.

## References

* Simon, J. L. et al. (1994), “Numerical expressions for precession formulae
  and mean elements for the Moon and the planets,” *Astronomy & Astrophysics*
  282, 663–683, DOI/bibliographic record:
  <https://ui.adsabs.harvard.edu/abs/1994A%26A...282..663S/abstract>.
* Meeus, J. (1998), *Astronomical Algorithms*, 2nd ed., Chapter 31,
  Table 31.A.
* Vondrak, J., Capitaine, N. & Wallace, P. (2011), “New precession expressions,
  valid for long time intervals,” *Astronomy & Astrophysics* 534, A22, with
  corrigendum, DOI: <https://doi.org/10.1051/0004-6361/201117274>.
* ERFA/pyerfa `ltp` and the project's registered Vondrak mean-obliquity
  realization, as wrapped by
  `/Users/giacomo/dev/libephemeris/libephemeris/precession_vondrak.py`.
* Project-owned coefficient table and evaluator:
  `/Users/giacomo/dev/libephemeris/libephemeris/planetary_mean_elements.py`.
* Project-owned frame consumers and semantic paths:
  `/Users/giacomo/dev/libephemeris/libephemeris/minor_bodies.py`,
  `/Users/giacomo/dev/libephemeris/libephemeris/astrometry.py`,
  `/Users/giacomo/dev/libephemeris/libephemeris/rebound_integration.py`, and
  `/Users/giacomo/dev/libephemeris/libephemeris/fast_calc.py`.
