# Disconnected finite-source strict certificate proposer

The private `_finite_source_strict_proposer` implements the project-authored
finite-policy proposal sequence in the reviewed validation contract
`validation/golden/specs/occultation_contact_residuals/finite-source-strict-proposer.md`
(SHA-256 `a0300389313401ef69789a39e3ead7b490edbecef1b216ebdbbc731dbb9b516b`,
validation commit `8754dba`). It consumes one exact rational
`_ValidatedGeometry`, one target, and explicit iteration and Decimal-precision
limits. It reads no astronomical state or reference data.

The proposal sequence extracts the target's exact affine side, transverse,
and optional gate forms at the origin and coordinate basis. An exact rational
bound on these forms fixes the Decimal step size. A fresh half-even Decimal
context computes ellipsoid best responses and projected dual iterates. The
iterates are candidate construction only: at each selected checkpoint, the
current and averaged primal and dual candidates are rationalized, repaired
inside their respective feasibility sets, and sent to the separate exact
`strict_reach` or `strict_miss` verifier. All four attempts at a checkpoint
are committed together. Only an exact verifier acceptance yields `REACH` or
`MISS`; exhausted iterations and Decimal faults remain `UNRESOLVED`.

The policy limits iteration count and Decimal precision. It does not bound
input size, exact rational growth, verification cost, or whole-call resources.
This module has no source receipt, public route, time search, or contact
equality decision. The two core nappes remain separate targets.
