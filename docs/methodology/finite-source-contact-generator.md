# Exact finite-source contact candidate generation

This page records the private, disconnected contact generator in
`_finite_source_contact_generator.py` and its exact-arithmetic helpers. Its
reviewed implementation contract is
`validation/golden/specs/occultation_contact_residuals/finite-source-contact-generator.md`
at validation commit `c9f609891aea4c041a642198450e3f6c60ace43a`, SHA-256
`9ea40fb03e873f4149985a09010cebdd5c8427d759c830bd4998b37fa9270b7c`.
The project-authored mathematical derivation and reviewed design are in audit
`B16_CANONICAL_CONTACT_GENERATOR_DRAFT.md`, SHA-256
`6dca6e6964cebe89a52275f27c1cd71226911bc8a09053693345a651471964bf`.
The audit rereview approved only a non-normative mathematical design. The
validation specification and its independent rereview supplied the
disconnected implementation contract; neither approved a runtime route.

## Exact fixed-state problem

The input is one rational, positive-definite receiving ellipsoid and the
source/occulter balls of the previously reviewed three-target finite-source
geometry. Construction repeats the full domain guard from the six base fields;
cached inverses and derived axes supplied by a caller carry no authority.
Each penumbral, umbral or antumbral target is handled separately. The output is
either one independently verified algebraic contact payload with a complete
branch trace, or a complete exhaustion trace. Exhaustion is an internal
absence-of-contact result, not a physical `MISS` certificate.

For the target's signed affine cone forms `a(x)=a0+v·x` and `z(x)=z0+Zx`,
the generator expands the cone boundary as

```text
f(x)=a(x)^2-H ||z(x)||^2=x^T C x+2d^T x+e,
C=v v^T-H Z^T Z,  d=a0 v-H Z^T z0,  e=a0^2-H ||z0||^2.
```

The receiver boundary is `q(x)=x^T A x=1`. At a smooth tangency the
positive multiplier `lambda` obeys `A x=lambda(C x+d)`. Setting
`M=A-lambda C`, `D=det(M)` and `N=lambda adj(M)d` gives two rational
polynomials, `Pq=N^T A N-D^2` and
`Pf=N^T C N+2D d^T N+eD^2`. Exact factorization of their gcd, real-root
isolation and sign decisions enumerate every positive nonsingular candidate.
The zero-determinant multiplier is transferred to the separate singular
branch, without dividing by `D`.

At the unique positive singular multiplier, the exact rank-two linear system
is reduced to `x=x0+t v0`. Substitution into `q=1` produces a quadratic over
`Q[lambda]`. Its discriminant decides zero, one or two real points. If a
point's free coordinate lies outside `Q[lambda]`, the producer constructs an
exact quadratic extension, selects the first positive integer `c` for which
`theta=lambda+c t` is primitive, and expresses the payload in `Q[theta]`.
Exact rational multiplication matrices and characteristic polynomials certify
field degree and coordinate conversion. The core apex is checked separately
as a rational point with its own dual support and physical target test.

The generator checks stationarity, receiver equality, cone equality, signed
side and physical gate before constructing a smooth payload. Accepted
payloads are then checked by the separately implemented
`_finite_source_algebraic_contact.py` verifier. A discrepancy in a defining
identity or in the verifier's decision is an invariant failure; ordinary
side, gate or apex dual failure is recorded as a rejected candidate.

## Representation and source boundary

The producer uses the project's existing `python-flint` dependency for exact
rational polynomial factorization, matrices and characteristic polynomials.
Its real-root counts and all decisive sign choices use independently written
rational Sturm and interval arithmetic. No floating-point midpoint decides
an equality, sign or branch. Each accepted payload uses a primitive integer
minimal polynomial with positive leading coefficient, the first isolating
centered dyadic cell, and reduced rational coordinate polynomials in one
selected real embedding. These are mathematical in-memory records; no wire
serialization or authenticated replay is specified here.

The generator has no astronomical state acquisition, public routing, temporal
search, strict reach/miss producer, source receipt or resource cap. It cannot
serve a runtime classifier until a separately reviewed boundary caps the
entire generator and verifier invocation and converts exhaustion of resources
to a typed unresolved result. It neither consults reference implementation
materials nor stores reference API outputs.
