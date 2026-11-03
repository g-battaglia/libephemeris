//! Consumed atmospheric roles for selected visibility geometry.

use crate::errors::Result;
pub use crate::vislim::VisibilityOpticalModel;
use crate::vislim::VisibilityPlace;

/// Atmospheric operands remain owned and unconverted until their Source role.
/// Stored scalar geometry consumes pressure/temperature in each public
/// horizontal call. Native planet vectors do not; the Native stellar public
/// horizon call does. Model construction occurs only after the below-horizon
/// decision in both routes.
pub trait VisibilityAtmosphericAcquisition: Send + Sync {
    /// Public ECL2HOR with the actual retained raw pressure/temperature.
    fn horizontal(
        &self,
        jd: f64,
        site: [f64; 3],
        place: VisibilityPlace,
    ) -> Result<VisibilityPlace>;
    /// Acquire the actual defining factory and create one model after the
    /// above-horizon guard. No phase/separation read precedes this operation.
    fn model(&self) -> Result<Box<dyn VisibilityOpticalModel>>;
}
