# Exact finite-source algebraic contact certificates

This page documents the private, disconnected
`_finite_source_algebraic_contact` checker. Its independently reviewed
contract is
`validation/golden/specs/occultation_contact_residuals/finite-source-algebraic-contact-verifier.md`
at validation commit `658875fd084ee339873d708687ff99d365b7bd9c`
(SHA-256 `b2aca7a287a9fb41711e2b6fe74c2cdc1f92bc396c4fa05934e722bb4a5508cd`).
The mathematical soundness argument is the project-authored and
independently reviewed audit
`B16_ALGEBRAIC_CONTACT_CERTIFICATE_CANDIDATE.md`. The exact root encoding,
sign method and singular discriminator are recorded in the separately
reviewed audit `B16_ALGEBRAIC_VERIFIER_ENCODING_DRAFT.md`. No reference
ephemeris algorithm, data or persisted output supplies this checker.

## Proof domain and one-root representation

The caller provides the exact rational finite-source geometry from the
strict checker. This checker reconstructs it from its base fields, so a
caller-supplied inverse, axis or domain flag cannot establish contact.
The source and Moon balls must be disjoint from the entire solid
receiving ellipsoid, and the reviewed positive-slope three-target model
must apply. Each decision is for one penumbral, umbral or antumbral target;
the module does not combine nappes or classify a public eclipse.

One exact integer polynomial `P` and a rational open interval with exactly
one root select a real algebraic number `theta`. All dual coefficients are
rational polynomials in that same `theta`. `P` is primitive and squarefree
with positive leading coefficient; irreducibility is unnecessary because
the interval selects one real root. A degree-one polynomial is enough for
rational coefficients. Equivalent values can have different encodings;
this disconnected verifier checks value truth rather than wire identity.

The verifier decides a polynomial sign at `theta` without a tolerance.
First, the gcd with `P` and an exact Sturm root count decide whether the
expression vanishes at the selected root. Otherwise it bisects the
isolating interval with exact root counts until rational interval Horner
evaluation has a strict sign. For a nonzero value, this terminates because
the interval shrinks to `theta`. This is an unbounded exact algorithm; it
is not a production time or memory guarantee.

## Contact certificate

For one target, write the oriented cone and optional physical gate as
exact rational affine forms `a(x)=a0+v^T x`, `z(x)=z0+Zx` and
`g(x)=g0+t^T x`. With exact dual coefficients `tau,beta,gamma`, the
checker requires `tau>=0`, `||beta||²<=H tau²` and `gamma>=0` when the
target has a gate. It constructs

```text
F(x)=tau a(x)+beta·z(x)+gamma g(x)=n·x+c,
b=A^-1 n,
```

and checks `c<0` and `c²=n·b`. This makes `F` tangent to the ellipsoid
at the unique point `x*=-b/c`. The target's exact side, squared cone and
gate checks are performed through the polynomial numerators

```text
A_num=a0 c-v·b,
Z_num=z0 c-Zb,
G_num=g0 c-t·b.
```

Since `c<0`, membership requires `A_num<=0`,
`A_num²-H||Z_num||²>=0`, and `G_num<=0` when gated. Thus the verifier
never divides by an algebraic `c` and never accepts a separately supplied
support point. The dual condition proves `F>=0` throughout the target;
the support condition proves `F<=0` throughout the ellipsoid; membership
at the tangent point proves exact per-target `CONTACT`. A failed payload
is neither a miss nor a reach proof.

The focused synthetic cases include a positive singular multiplier with
`A=I`, `r=11/4`, `m=1/4`, `B=(-55/12,5/4,0)`,
`p=(5/12,5/4,0)` and penumbral contact at `(3/5,4/5,0)`. Here the
multiplier is `1/225`, `a=9`, the physical gate is `5/3>0`, and the
exact dual coefficients are `tau=2/25`, `beta=(0,0,8/25)`. The code's
independent checker does not use this generator multiplier; the case
ensures that a true singular contact is not discarded.

This module neither obtains astronomical states nor generates a contact
candidate. It has no runtime work cap and cannot be called by a public
path until a separately reviewed boundary limits the **entire** invocation
and returns typed `UNRESOLVED` on exhaustion. Authenticated source
receipts, complete replay and time-dependent contacts remain separate.
