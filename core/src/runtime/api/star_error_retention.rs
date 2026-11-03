//! Ordinary owned results for prepared catalog and fixed-star calculations.
use crate::errors::Error;

/// A completed fixed-star result independent of context/provider lifetimes.
/// Errors are ordinary native values; no interpreter presentation is retained.
#[derive(Debug)]
pub struct RetainedStarResult<T> {
    /// Native result, including its actual typed data/domain failure.
    pub result: Result<T, Error>,
}
impl<T> RetainedStarResult<T> {
    /// Return the native catalog/calculation result unchanged.
    pub fn into_catalog_result(self) -> Result<T, Error> {
        self.result
    }
    pub(crate) fn success(value: T) -> Self {
        Self { result: Ok(value) }
    }
    pub(crate) fn failure(error: Error) -> Self {
        Self { result: Err(error) }
    }
    pub(crate) fn from_serving(
        value: crate::positions::star_serving::CapturedStarResult<T>,
    ) -> Self {
        Self {
            result: value.result,
        }
    }
}
