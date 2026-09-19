# Exact finite-source strict shadow certificates

This document records the project-authored mathematics used by the private
`_finite_source_certificates` checker. The fixed, independently reviewed
implementation contract is
`validation/golden/specs/occultation_contact_residuals/finite-source-strict-verifier.md`
at validation commit `03292c3390b9a843e0b5547c69c54e5e933b69ec`
(SHA-256 `e96d262c53b0b5af0edf70ca3f68f2a98f81987160ec4a8c84d69ee401faff46`).
The underlying project mathematical proofs are recorded in the audit files
`B16_FULL_CORE_DOMAIN_CANDIDATE.md` and
`B16_TRUNCATED_CONE_CERTIFICATES_CANDIDATE.md`, each with an independent
high-reasoning review. Those documents and this page derive the predicates
from the declared geometric model; no reference implementation or reference
output supplies a formula or fitted parameter.

## Exact model and domain

Let `B` and `p` be source and occulter sphere centres, with radii
`r>m>0`, all in the same declared Euclidean coordinate frame and units.
The receiver is the closed solid ellipsoid

```text
K = {x in R³ : x^T A x <= 1},
```

where `A` is an exact symmetric positive-definite matrix. The caller also
supplies `R>0`, and the checker proves `A-I/R²` positive semidefinite before
using `R` as a Euclidean enclosure of `K`. Set

```text
u=p-B,  U=u·u,  s=r+m,  d=r-m,
v=x-p,  w=v·u,  z=u×v,  T=z·z.
```

`U>s²` proves the source and occulter closed balls disjoint. The strict
inequalities `||B||²>(R+r)²` and `||p||²>(R+m)²` then suffice to prove the
whole receiver disjoint from both balls. A failed enclosure or receiver
guard is an **uncertified domain**, since these sufficient checks do not
decide intersection of an arbitrary ellipsoid and sphere. The `U` guard
failure is a separate **excluded theorem domain**. Valid `r<=m` or zero
radii need a different theorem and remain **unsupported slopes**. None of
these classifications is a physical shadow miss.

## Oriented targets

Under the proved domain, the finite-ray penumbral target and both core
nappes are the following exact closed sets. Write `H_p=U-s²>0` and
`H_c=U-d²>0`:

```text
penumbra:  g_p=w+ms >= 0,  a_p=Um+sw >= 0,  a_p² >= H_p T;
umbra:     g_c=w-md >= 0,  a_u=Um-dw >= 0,  a_u² >= H_c T;
antumbra:                    a_a=dw-Um >= 0,  a_a² >= H_c T.
```

The antumbral gate is redundant in this positive-slope domain. The physical
core is the union of the two core nappes, so one nappe's miss cannot prove
a core miss. Keeping the affine side and gate checks prevents a positive
squared residual from admitting the opposite, physically unlit nappe.

## Why strict payload checks are sound

For each target, a rational point inside `K` with strict side, gate where
present, and cone inequality is a strict interior overlap. The reviewed
whole-solid disjointness and convex segment argument then gives strict
shadow reach at the receiver surface. A boundary point or a zero residual
does not establish strict reach.

For a strict miss, take exact rational multipliers `tau>=0`, `beta`, and
`lambda>=0` for a gated target, with `beta·beta<=tau²H`; antumbra requires
`lambda=0`. Form the affine function

```text
F(x)=tau*a(x)+beta·z(x)+lambda*g(x),
```

omitting the gate term for antumbra. On its target, `a>=sqrt(H)||z||`;
Cauchy–Schwarz therefore gives
`tau*a+beta·z >= (tau*sqrt(H)-||beta||)||z|| >= 0`, and the gate term is
nonnegative. Write `F(x)=c+n·x`, obtaining `c=F(0)` and
`n_i=F(e_i)-c` from the supplied exact payload. The ellipsoid support
identity is

```text
max_{x in K} n·x = sqrt(n^T A^-1 n).
```

Thus `c<0` and `c²>n^T A^-1 n` prove `F<0` everywhere on `K`, separating
the entire target from the receiver. Every comparison is exact rational
arithmetic; the code never rounds an angle or accepts an independently
supplied axis or affine separator. Equality is not a strict miss.

The checker validates caller-supplied payloads only. It does not produce
certificates, decide exact contact, authenticate astronomical states,
bound ephemeris error, select data assets, search in time, or change a
public result. Those require separate specifications and verification.
