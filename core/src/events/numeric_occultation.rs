//! Fixed-epoch numeric lunar-occultation composition over owned acquisitions.
//! The geometry callback precedes the center-point circumstance callback;
//! caller UT1 and the exact numeric body identity are retained throughout.

use super::{OccultWhereOutcome, SolHowCircumstances, WhereGeometry, compose_occult_where_outcome};
use crate::errors::Error;

/// Call-owned numeric occultation acquisition. Geometry consumes the public
/// Moon/target scalar positions; local circumstances consume an independent
/// bundled UT1 WGS84 apparent-place channel at the actual returned center.
/// Neither method may retain a provider borrow beyond the owning call.
pub(crate) trait NumericOccultWhereSource {
    /// Global true-equatorial Moon/target shadow geometry at UT1 Julian Day.
    /// Longitude east and latitude north are geodetic degrees; the core-shadow
    /// diameter is kilometres. Errors propagate without a no-event substitute.
    fn global(&mut self, jd_ut1: f64, body: i32, flags: i32) -> Result<WhereGeometry, Error>;
    /// Local apparent ecliptic target/Moon circumstances at the actual center,
    /// with uncapped finite-disc WHERE obscuration and independent horizon TT.
    /// The three site words are east longitude/north latitude degrees/metres.
    fn local(
        &mut self,
        jd_ut1: f64,
        body: i32,
        flags: i32,
        site: [f64; 3],
    ) -> Result<SolHowCircumstances, Error>;
}

/// Compose fixed-instant lunar-occultation character and center circumstances.
/// Global geometry is acquired once, followed by local circumstances even for
/// character zero. The result contains ten geoposition and twenty attribute
/// words; `attr[3]` is replaced by the actual global core-shadow diameter km.
/// This leaf neither validates epochs early nor applies a local altitude gate.
/// Public Moon-self rejection/negative-id folding precede this internal leaf.
pub(crate) fn numeric_occult_where_via(
    source: &mut impl NumericOccultWhereSource,
    jd_ut1: f64,
    body: i32,
    flags: i32,
) -> Result<OccultWhereOutcome, Error> {
    let geometry = source.global(jd_ut1, body, flags)?;
    let circumstances = source.local(
        jd_ut1,
        body,
        flags,
        [geometry.lon_deg(), geometry.lat_deg(), 0.0],
    )?;
    Ok(compose_occult_where_outcome(
        &geometry,
        circumstances.attr(),
    ))
}

/// Local finite-disc circumstances from ecliptic positions in degrees and AU.
/// Target and lunar angular radii must be real. The horizon and Saros
/// callbacks supply local angles and cycle indices after domain validation.
/// Non-finite coordinates, zero distance and impossible radii return errors.
pub(crate) fn numeric_occult_how_via(
    target: (f64, f64, f64),
    moon: (f64, f64, f64),
    radius_au: f64,
    is_solar: bool,
    horizontal: &mut dyn FnMut() -> Result<(f64, f64, f64), Error>,
    height_m: f64,
    saros: &mut dyn FnMut() -> Result<(f64, f64), Error>,
) -> Result<SolHowCircumstances, Error> {
    if radius_au > 0.0 {
        occult_radius_domain(radius_au, target.2)?;
    }
    occult_radius_domain(super::ECL_RMOON_AU, moon.2)?;
    super::sol_how_core_lazy(
        target, moon, radius_au, is_solar, true, horizontal, height_m, saros,
    )
}

/// Check that a finite physical radius has a real angular image.
fn occult_radius_domain(radius: f64, distance: f64) -> Result<(), Error> {
    if distance == 0.0 {
        return Err(Error::RationalDivisionByZero);
    }
    let ratio = radius / distance;
    if !ratio.is_finite() || !(-1.0..=1.0).contains(&ratio) {
        return Err(Error::MathDomain { value: ratio });
    }
    Ok(())
}

#[cfg(test)]
#[path = "numeric_occultation/tests.rs"]
mod tests;
