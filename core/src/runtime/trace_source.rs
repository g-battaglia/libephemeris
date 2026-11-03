//! Owned provider names for successful trace publications.

/// The actual provider label consumed by a calculation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TraceSource {
    /// An owned UTF-8 name, independent of the caller's string lifetime.
    Named(String),
}

impl TraceSource {
    /// Retain a known provider name without inferring another selection.
    pub fn named(source: impl Into<String>) -> Self {
        Self::Named(source.into())
    }

    /// Borrow the retained provider name.
    pub fn as_named(&self) -> Option<&str> {
        match self {
            Self::Named(source) => Some(source),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_own_their_storage_and_compare_by_value() {
        let mut input = String::from("LEB");
        let source = TraceSource::named(&input);
        input.clear();
        let retained = source.clone();
        drop(source);
        assert_eq!(retained.as_named(), Some("LEB"));
        assert_eq!(retained, TraceSource::Named("LEB".into()));
        assert_ne!(retained, TraceSource::named("JPL"));
    }
}
