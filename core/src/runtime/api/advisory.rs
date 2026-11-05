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

impl Advisory {
    /// Convert actual engine metadata without formatting or class inference.
    /// Numeric payloads are moved unchanged, including signed zero/NaN bits.
    /// This adapter is crate-private, not part of the primary contract.
    pub(crate) fn from_engine(warning: crate::errors::Warning) -> Self {
        use crate::errors::Warning;

        // Capture the existing display route before moving owned payloads.
        let message = warning.to_string();
        let detail = match warning {
            Warning::MeeusPolynomial {
                centuries_from_j2000,
                outside_valid_range,
            } => AdvisoryDetail::LunarPolynomial {
                centuries_from_j2000,
                outside_valid_range,
            },
            Warning::UnknownSiderealMode { mode } => AdvisoryDetail::UnknownSiderealMode { mode },
            Warning::UnsupportedSiderealProjection { bits } => {
                AdvisoryDetail::UnsupportedSiderealProjection { bits }
            }
            Warning::MalformedConfig { path, detail } => {
                AdvisoryDetail::MalformedConfiguration { path, detail }
            }
            Warning::ReboundSelfInteraction { perturber } => {
                AdvisoryDetail::SelfInteractionExcluded { perturber }
            }
        };
        Self { detail, message }
    }
}

impl fmt::Display for Advisory {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.message)
    }
}

#[cfg(test)]
mod tests {
    use super::{Advisory, AdvisoryDetail};
    use crate::errors::Warning;

    #[test]
    fn conversion_retains_every_payload_and_original_message() {
        let inputs = [
            Warning::UnknownSiderealMode { mode: -193 },
            Warning::UnsupportedSiderealProjection { bits: 0x4000 },
            Warning::MalformedConfig {
                path: "actual-config.toml".to_owned(),
                detail: "actual parser detail".to_owned(),
            },
            Warning::ReboundSelfInteraction {
                perturber: "actual target".to_owned(),
            },
        ];
        for input in inputs {
            let original = input.to_string();
            let output = Advisory::from_engine(input.clone());
            assert_eq!(output.message, original);
            assert_eq!(output.to_string(), original);
            match (input, output.detail) {
                (
                    Warning::UnknownSiderealMode { mode: a },
                    AdvisoryDetail::UnknownSiderealMode { mode: b },
                ) => assert_eq!(a, b),
                (
                    Warning::UnsupportedSiderealProjection { bits: a },
                    AdvisoryDetail::UnsupportedSiderealProjection { bits: b },
                ) => assert_eq!(a, b),
                (
                    Warning::MalformedConfig { path: a, detail: c },
                    AdvisoryDetail::MalformedConfiguration { path: b, detail: d },
                ) => {
                    assert_eq!(a, b);
                    assert_eq!(c, d);
                }
                (
                    Warning::ReboundSelfInteraction { perturber: a },
                    AdvisoryDetail::SelfInteractionExcluded { perturber: b },
                ) => assert_eq!(a, b),
                _ => panic!("advisory payload changed during conversion"),
            }
        }
    }

    #[test]
    fn conversion_preserves_float_payload_bits() {
        for bits in [
            0_u64,
            0x8000_0000_0000_0000,
            0x4026_0000_0000_0000,
            0x7ff8_0000_0000_0042,
        ] {
            for outside_valid_range in [false, true] {
                let warning = Warning::MeeusPolynomial {
                    centuries_from_j2000: f64::from_bits(bits),
                    outside_valid_range,
                };
                let message = warning.to_string();
                let output = Advisory::from_engine(warning);
                assert_eq!(output.message, message);
                let AdvisoryDetail::LunarPolynomial {
                    centuries_from_j2000,
                    outside_valid_range: actual,
                } = output.detail
                else {
                    panic!("wrong advisory kind");
                };
                assert_eq!(centuries_from_j2000.to_bits(), bits);
                assert_eq!(actual, outside_valid_range);
            }
        }
    }
}
