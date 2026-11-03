//! Numerical foundation kernels: angles, Chebyshev evaluation, exact arithmetic,
//! and small linear algebra.
//!
//! These are source-independent mathematical kernels with documented definitions,
//! not release-algorithm ports: exact identities plus stated rounding and overflow
//! behavior are the contract. Release-convention confirmation (for example the
//! Chebyshev coefficient ordering against the LEB reader) stays with the owning
//! reader tasks and is not claimed here.
//!
//! Ordered layers:
//! - [`angles`]: radian/degree normalization with the half-turn and snap
//!   conventions, plus spherical/Cartesian transforms with the pole `lon = 0`
//!   rule over [`crate::math::linalg::Vec3`].
//! - [`chebyshev`]: Clenshaw value evaluation, exact basis-recurrence
//!   `x`-derivatives, and [`crate::types::Segment`] scaling of derivative
//!   rows to the time domain.
//! - [`exact`]: reduced `i128` rational arithmetic and outward-rounded
//!   floating-point intervals over `+`, `-`, `*`, `/`, `sqrt`.
//! - [`linalg`]: three-vectors and 3x3 matrices (dot/cross/norm, active
//!   right-handed rotations, transpose, adjugate inverse, orthogonality
//!   checks) in the written association order.
//!
//! Frames, units, epochs: angles are radians unless the name says `deg`;
//! the Chebyshev argument is the normalized segment coordinate in `[-1, 1]`;
//! matrices are row-major `m[row][col]`; rationals and intervals are
//! dimensionless. There are no epochs in this module: time scaling enters
//! only through the caller-supplied segment at the [`chebyshev`] boundary.
//!
//! Refusals: typed errors, never silent substitution — empty coefficient
//! slices, non-finite arguments, invalid derivative orders, overflowing
//! segment-derivative scales, zero-vector normalize/spherical conversions,
//! negative spherical radii, zero rational denominators, integer overflow,
//! inverted or non-finite interval bounds, zero-spanning interval divisors,
//! and negative interval square roots.
//!
//! Domains and limits: every entry accepts all finite `f64` inputs and
//! rejects the rest; rationals are exact inside the checked `i128` range
//! and report overflow outside it; intervals enclose the five stated
//! operations only and offer no transcendental extensions.
//!

pub mod angles;
pub mod chebyshev;
pub mod exact;
pub mod linalg;
