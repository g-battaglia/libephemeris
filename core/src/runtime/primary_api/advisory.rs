//! Owned, backend-independent advisories for position results.
//!
//! No advisory is inferred from diagnostic prose. The producer supplies the
//! structured detail and its actual display text. Absence is an empty vector.
//! The contract enum is closed: adding a variant requires a coordinated new
//! major version of each implementation of this primary API.

use std::fmt;

/// One non-fatal producer advisory, separate from the native diagnostic buffer.
#[derive(Debug, Clone, PartialEq)]
pub struct Advisory {
    /// Original structured advisory data; units are documented per variant.
    pub detail: AdvisoryDetail,
    /// Original producer text, retained without parsing or reconstruction.
    pub message: String,
}

/// Structured non-fatal producer information, independent of engine error types.
///
/// This is a closed contract enum, intentionally not `#[non_exhaustive]`.
/// Adding a variant is a breaking change requiring a new major in both API
/// implementations. A producer must not invent a variant to label a diagnostic.
#[derive(Debug, Clone, PartialEq)]
pub enum AdvisoryDetail {
    /// A lunar polynomial is outside its preferred or declared validity range.
    LunarPolynomial {
        /// Epoch offset from J2000, in Julian centuries of 36,525 days.
        centuries_from_j2000: f64,
        /// Whether the producer reports the declared validity range exceeded.
        outside_valid_range: bool,
    },
    /// The producer replaced an unknown sidereal selection with its default.
    UnknownSiderealMode {
        /// Actual unrecognized integer selection.
        mode: i32,
    },
    /// The producer omitted unsupported sidereal projection bits.
    UnsupportedSiderealProjection {
        /// Actual omitted bit subset, dimensionless.
        bits: i32,
    },
    /// A configuration source was ignored after a parsing failure.
    MalformedConfiguration {
        /// Actual configuration source path.
        path: String,
        /// Original reason supplied by its parser.
        detail: String,
    },
    /// An integrator excluded a target's identical perturbing force source.
    SelfInteractionExcluded {
        /// Actual excluded source name, not a reconstructed body identifier.
        perturber: String,
    },
}

impl fmt::Display for Advisory {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.message)
    }
}
