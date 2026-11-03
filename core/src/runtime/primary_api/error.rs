//! Primary errors retain the genuine engine cause without changing Python errors.
use super::ErrorKind;
use std::{fmt, sync::Arc};

/// Owned primary failure and its optional real backend metadata.
///
/// `source()` retains the original engine error, including its structured fields
/// and causal chain. Primary classification never parses diagnostic wording.
#[derive(Clone, Debug)]
pub struct Error {
    kind: ErrorKind,
    operation: Option<String>,
    native_status: Option<i32>,
    native_scalar: Option<f64>,
    diagnostic: Option<String>,
    source: Option<Arc<crate::errors::Error>>,
}
impl Error {
    /// Construct an error from explicitly supplied backend metadata.
    ///
    /// # Arguments
    ///
    /// * `kind` - Closed broad failure class.
    /// * `operation` - Public operation name, when supplied.
    /// * `native_status` - Actual backend integer status, when supplied.
    /// * `native_scalar` - Actual backend scalar, retaining its binary64 bits.
    /// * `diagnostic` - Actual diagnostic text, when supplied.
    ///
    /// # Returns
    ///
    /// An owned error record. This constructor has no native cause; failures
    /// produced by calculation methods retain their original engine source.
    pub fn new(
        kind: ErrorKind,
        operation: Option<String>,
        native_status: Option<i32>,
        native_scalar: Option<f64>,
        diagnostic: Option<String>,
    ) -> Self {
        Self {
            kind,
            operation,
            native_status,
            native_scalar,
            diagnostic,
            source: None,
        }
    }
    /// Closed broad failure class; native structured details remain in `source()`.
    pub fn kind(&self) -> ErrorKind {
        self.kind
    }
    /// Actual public operation where the error crossed this boundary.
    pub fn operation(&self) -> Option<&str> {
        self.operation.as_deref()
    }
    /// Actual native integer status, only when the engine error supplies one.
    pub fn native_status(&self) -> Option<i32> {
        self.native_status
    }
    /// Actual rejected native scalar, preserving its binary64 bits.
    pub fn native_scalar(&self) -> Option<f64> {
        self.native_scalar
    }
    /// Actual engine diagnostic, without rewritten exception formatting.
    pub fn diagnostic(&self) -> Option<&str> {
        self.diagnostic.as_deref()
    }
    pub(crate) fn from_engine(operation: &str, error: crate::errors::Error) -> Self {
        let kind = classify(&error);
        let native_scalar = scalar(&error);
        let native_status = status(&error);
        let diagnostic = Some(error.to_string());
        Self {
            kind,
            operation: Some(operation.to_owned()),
            native_status,
            native_scalar,
            diagnostic,
            source: Some(Arc::new(error)),
        }
    }
}
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(source) = &self.source {
            return fmt::Display::fmt(source.as_ref(), f);
        }
        if let Some(message) = &self.diagnostic {
            return f.write_str(message);
        }
        write!(f, "{:?}", self.kind)
    }
}
impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.source
            .as_deref()
            .map(|value| value as &(dyn std::error::Error + 'static))
    }
}

// Classification follows structured variants, never diagnostic wording.
fn classify(error: &crate::errors::Error) -> ErrorKind {
    use crate::errors::{Error as E, ErrorCategory};
    match error {
        E::NormalizedTimeOutOfRange { .. }
        | E::EpochOutOfSegment { .. }
        | E::LebOutOfRange { .. }
        | E::SpkOutOfRange { .. }
        | E::SpkRegisteredOutOfCoverage { .. }
        | E::EphemerisRange { .. } => ErrorKind::OutOfRange,
        E::UnsupportedFlagCombination { .. }
        | E::SealedChannelNotServed { .. }
        | E::HouseSystemUnserved { .. }
        | E::MoonTheoryUnserved { .. }
        | E::LebUnsupportedVersion { .. }
        | E::LebUnsupportedCompression { .. }
        | E::SpkUnsupportedType { .. }
        | E::SpkUnsupportedFrame { .. }
        | E::SpkUnsupportedNdNi { .. } => ErrorKind::Unsupported,
        E::NBodyBackendUnavailable { .. }
        | E::DeltaTModelUnavailable { .. }
        | E::StarMagnitudeUnavailable { .. }
        | E::LebBytesUnavailable
        | E::LebReaderClosed
        | E::AyanamshaAnchorUnserved { .. }
        | E::GauquelinMethodUnserved { .. } => ErrorKind::Unavailable,
        E::BackendRange { error, .. } | E::Caused { error, .. } => classify(error),
        E::PublicCrossing {
            cause: Some(cause), ..
        } => classify(cause),
        _ => match error.category() {
            ErrorCategory::InputValidation => ErrorKind::InvalidInput,
            ErrorCategory::DataNotFound => ErrorKind::MissingData,
            ErrorCategory::Calculation => ErrorKind::Calculation,
            ErrorCategory::Configuration => ErrorKind::Configuration,
        },
    }
}

fn scalar(error: &crate::errors::Error) -> Option<f64> {
    use crate::errors::Error as E;
    match error {
        E::NonFiniteArgument { value } | E::InvalidObserver { value, .. } => Some(*value),
        E::BackendRange { error, .. } => scalar(error),
        E::Caused { error, cause } => scalar(error).or_else(|| scalar(cause)),
        E::PublicCrossing {
            cause: Some(cause), ..
        } => scalar(cause),
        _ => None,
    }
}

fn status(error: &crate::errors::Error) -> Option<i32> {
    use crate::errors::Error as E;
    match error {
        E::NBodyNative { code, .. } => Some(*code),
        E::HttpStatus { status, .. } => Some(i32::from(*status)),
        E::BackendRange { error, .. } => status(error),
        E::Caused { error, cause } => status(error).or_else(|| status(cause)),
        E::PublicCrossing {
            cause: Some(cause), ..
        } => status(cause),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::errors::Error as E;

    #[test]
    fn structured_range_support_and_availability_remain_distinct() {
        let cases = [
            (
                E::LebOutOfRange {
                    body_id: 0,
                    jd: 1.0,
                    start: 2.0,
                    end: 3.0,
                },
                ErrorKind::OutOfRange,
            ),
            (
                E::UnsupportedFlagCombination {
                    entry: "actual native entry",
                    flags: 8,
                },
                ErrorKind::Unsupported,
            ),
            (
                E::NBodyBackendUnavailable {
                    backend: "actual backend",
                },
                ErrorKind::Unavailable,
            ),
        ];
        for (source, kind) in cases {
            let error = Error::from_engine("actual operation", source.clone());
            assert_eq!(error.kind(), kind);
            assert_eq!(error.to_string(), source.to_string());
        }
    }

    #[test]
    fn real_native_status_survives_a_causal_wrapper() {
        let error = Error::from_engine(
            "actual operation",
            E::PublicCrossing {
                message: "actual outer diagnostic".into(),
                cause: Some(Box::new(E::NBodyNative {
                    stage: 2,
                    code: -17,
                    diagnostic: "actual native diagnostic".into(),
                })),
            },
        );
        assert_eq!(error.native_status(), Some(-17));
        assert_eq!(error.to_string(), "actual outer diagnostic");
        let http = Error::from_engine(
            "actual request",
            E::HttpStatus {
                status: 503,
                purpose: "actual request".into(),
            },
        );
        assert_eq!(http.native_status(), Some(503));
    }
}
