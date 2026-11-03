//! Native SDK for provider integration and the engine's complete Rust capabilities.
//!
//! This namespace retains native session configuration, resource/acquisition
//! interfaces, model-specific operations and unconverted native error/result
//! types. Use [`crate::api`] for the deliberately smaller cross-engine primary
//! contract. The native SDK is independently owned by LibEphemeris; its extra
//! operations are not part of that primary contract.
//!
//! `SessionBuilder::isolated` builds an explicitly configured local calculation
//! owner. `crate::config::Runtime` provides mutable application configuration
//! and lazy acquisition; its snapshots are these same native sessions. The
//! feature-gated Python implementation bridge uses the same engine seam, so
//! no calculation or Python exception is reconstructed through the primary API.

pub use crate::engine_api::*;
