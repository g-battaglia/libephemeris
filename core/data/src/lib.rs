// SPDX-License-Identifier: AGPL-3.0
//! Immutable source bytes for LibEphemeris's optional bundled base image.
//!
//! This crate performs no calculations, I/O, conversion or decompression. The
//! engine validates and interprets the exact LEB2 container, including each
//! channel's actual domain. A dataset's coverage is not an accuracy guarantee.
//! Scientific attribution and redistribution notices accompany the archive.

#![deny(missing_docs)]

/// The unchanged published base-core LEB2 image, backed by static read-only data.
///
/// Its data-v4 SHA-256 is `b891de02f82b68346046da67389fd3713aa81ab299287c632cb20cc396397964`.
/// Use the engine's `bundled-data` feature rather than interpreting these opaque
/// bytes yourself. No allocation or runtime asset download is required.
pub static BASE_CORE_LEB2: &[u8] = include_bytes!("../leb2/base_core.leb2");
