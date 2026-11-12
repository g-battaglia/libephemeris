// SPDX-License-Identifier: AGPL-3.0
//! Distribution metadata shared by the native CLI and the Python binding.
//!
//! The Rust executable identifies the core package. A binding can install its
//! own immutable distribution label before exposing the shared CLI/status
//! surface. This metadata never selects a model, provider or scientific path.

#[cfg(feature = "python-internals")]
use std::sync::OnceLock;

#[cfg(feature = "python-internals")]
static BINDING_VERSION: OnceLock<&'static str> = OnceLock::new();

/// Return the embedding distribution's version, or the standalone core version.
pub fn distribution_version() -> &'static str {
    #[cfg(feature = "python-internals")]
    if let Some(version) = BINDING_VERSION.get() {
        return version;
    }
    env!("CARGO_PKG_VERSION")
}

/// Install the binding's immutable version label before exposing its CLI.
///
/// Reinitializing the same binding is harmless. A conflicting label is refused
/// rather than changing process-wide metadata underneath an existing caller.
#[cfg(feature = "python-internals")]
#[doc(hidden)]
pub fn set_distribution_version(version: &'static str) -> Result<(), &'static str> {
    install_version(&BINDING_VERSION, version)
}

#[cfg(feature = "python-internals")]
fn install_version(
    slot: &OnceLock<&'static str>,
    version: &'static str,
) -> Result<(), &'static str> {
    if slot.set(version).is_ok() || slot.get() == Some(&version) {
        Ok(())
    } else {
        Err("distribution version is already initialized with a different label")
    }
}

#[cfg(all(test, feature = "python-internals"))]
mod tests {
    use super::*;

    #[test]
    fn embedding_label_is_immutable_and_reinitialization_is_idempotent() {
        // A local slot avoids mutating process metadata for other unit tests.
        let slot = OnceLock::new();
        assert_eq!(install_version(&slot, "4.0.0b1"), Ok(()));
        assert_eq!(install_version(&slot, "4.0.0b1"), Ok(()));
        assert!(install_version(&slot, "different").is_err());
        assert_eq!(slot.get(), Some(&"4.0.0b1"));
    }
}
