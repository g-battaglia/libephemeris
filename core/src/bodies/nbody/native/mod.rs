// SPDX-License-Identifier: AGPL-3.0
//! Independently derived Rust numerical foundations and N-body engine.
//!
//! The optional `nbody` feature selects this implementation beneath the public
//! runtime, binding and offline-source interfaces. `native-nbody` is a compatible
//! feature alias. No GPL implementation is an implementation input. Qualification
//! of actual built artifacts remains separate from successful compilation.
//!
//! The standalone test target also compiles this module without the parent crate:
//! `rustc --edition=2024 --test src/bodies/nbody/native/mod.rs -o /tmp/native-nbody`.
//! This checks that the foundation itself has no FFI or dependency on the C engine.

#[cfg(feature = "nbody")]
pub(super) mod bridge;
mod catalogue;
mod classic_de;
mod conic;
mod encounter;
// Keep the standalone mathematical target dependency-free. The combined
// provider uses the crate's existing SPK implementation under N-body builds.
#[cfg(feature = "nbody")]
mod ephemeris;
#[cfg(feature = "nbody")]
mod evaluator;
mod forces;
mod kepler;
mod newtonian;
#[cfg(feature = "nbody")]
mod poles;
mod radau;
mod radau_tables;
mod solar;
mod state;
#[cfg(feature = "nbody")]
mod trajectory;

#[cfg(test)]
mod tests;
