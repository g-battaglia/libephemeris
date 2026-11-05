//! Complete owned transport for the primary node/apse operation.

/// Four six-component node/apse vectors in the request's frame and units.
///
/// Default spherical components are degrees/degrees/AU, followed by matching
/// daily rates. Cartesian/radian options retain their existing units. FOPOINT
/// places the second focus in aphelion. Producer NaNs and signed zeros retain
/// their original words; this record performs no numerical transformation.
#[derive(Debug, Clone, PartialEq)]
pub struct NodesApsides {
    /// Ascending node, including all three rates.
    pub ascending: [f64; 6],
    /// Descending node, including all three rates.
    pub descending: [f64; 6],
    /// Perihelion, including all three rates.
    pub perihelion: [f64; 6],
    /// Aphelion or second focus, including all three rates.
    pub aphelion: [f64; 6],
    /// Actual producer integer status; absent when none was returned.
    pub status: Option<i32>,
    /// Actual diagnostic channel; None is distinct from an empty diagnostic.
    pub diagnostic: Option<String>,
}

impl NodesApsides {
    /// Move the four native arrays without evaluating or normalizing any word.
    pub(crate) fn from_vectors(vectors: [[f64; 6]; 4]) -> Self {
        let [ascending, descending, perihelion, aphelion] = vectors;
        Self {
            ascending,
            descending,
            perihelion,
            aphelion,
            status: None,
            diagnostic: None,
        }
    }
}
