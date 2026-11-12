//! Closed application-facing Rust API over the independent native engine.
//!
//! Choose this namespace with `use libephemeris_core::api as eph`. Pure
//! operations are free functions; configured calculations borrow an immutable
//! concrete Session. Native providers, data setup and the feature-gated Python
//! implementation seam remain outside this namespace. No backend trait or
//! shared astronomy is involved. All results own their complete native words.
mod constants;
pub use constants::*;
mod selectors;
pub use selectors::*;
mod types;
pub use types::*;
mod error;
pub use error::Error;
mod advisory;
pub use advisory::{Advisory, AdvisoryDetail};
mod pure;
pub use pure::*;
mod operations;

/// Immutable configured calculation owner. Cloning retains its real data owners.
#[derive(Clone, Debug)]
pub struct Session {
    pub(crate) engine: crate::engine_api::Session,
}
/// Explicit deterministic builder; data setup is supplied by native SDK context.
#[derive(Clone, Debug)]
pub struct SessionBuilder {
    engine: crate::engine_api::SessionBuilder,
}
impl Default for SessionBuilder {
    fn default() -> Self {
        Self::new()
    }
}
impl SessionBuilder {
    /// Create explicit local defaults without environment resolution or file I/O.
    ///
    /// # Returns
    ///
    /// A builder selecting sealed local LEB mode and networking, SMH2016 Delta-T
    /// with IERS disabled, a zero-degree/zero-metre observer, Fagan/Bradley
    /// sidereal mode, and automatic tidal acceleration and Delta-T. Data owners
    /// can be supplied through the native context/Runtime adapters outside api.
    pub fn new() -> Self {
        Self {
            engine: crate::engine_api::SessionBuilder::new(),
        }
    }
    /// Set the observer coordinates, checked when building the session.
    ///
    /// # Arguments
    ///
    /// * `longitude` - East longitude in degrees.
    /// * `latitude` - North latitude in degrees.
    /// * `height` - Height above the reference surface in metres.
    ///
    /// # Returns
    ///
    /// This builder for chained configuration; invalid domains fail in `build`.
    pub fn topo(&mut self, longitude: f64, latitude: f64, height: f64) -> &mut Self {
        self.engine.topo(longitude, latitude, height);
        self
    }
    /// Select the native sidereal mode and its defining epoch/offset.
    ///
    /// # Arguments
    ///
    /// * `mode` - Sidereal mode and projection bits.
    /// * `epoch` - Defining Julian day interpreted by the native selected mode.
    /// * `offset` - Defining ayanamsha in degrees.
    ///
    /// # Returns
    ///
    /// This builder; the native mode-selection diagnostics are retained at build.
    pub fn sid_mode(&mut self, mode: i32, epoch: f64, offset: f64) -> &mut Self {
        self.engine.sid_mode(mode, epoch, offset);
        self
    }
    /// Override tidal acceleration for the native Delta-T model.
    ///
    /// # Arguments
    ///
    /// * `value` - Acceleration in arcseconds per Julian century squared.
    ///
    /// # Returns
    ///
    /// This builder; a non-finite value fails in `build`.
    pub fn tid_acc(&mut self, value: f64) -> &mut Self {
        self.engine.tid_acc(value);
        self
    }
    /// Override TT minus UT1 without replacing independently defined frame clocks.
    ///
    /// # Arguments
    ///
    /// * `days` - Explicit TT minus UT1 in days.
    ///
    /// # Returns
    ///
    /// This builder; a non-finite value fails in `build`.
    pub fn delta_t_override(&mut self, days: f64) -> &mut Self {
        self.engine.delta_t_override(days);
        self
    }
    /// Restore the configured Delta-T model by clearing its scalar override.
    ///
    /// # Returns
    ///
    /// This builder with the existing model and frame-clock policy retained.
    pub fn clear_delta_t_override(&mut self) -> &mut Self {
        self.engine.clear_delta_t_override();
        self
    }
    /// Set the atmospheric temperature gradient.
    ///
    /// # Arguments
    ///
    /// * `value` - Lapse rate in K/m.
    ///
    /// # Returns
    ///
    /// This builder; a non-finite value fails in `build`.
    pub fn lapse_rate(&mut self, value: f64) -> &mut Self {
        self.engine.lapse_rate(value);
        self
    }
    /// Validate the configured settings and retain their native data owners.
    ///
    /// # Returns
    ///
    /// An immutable session sharing the configured images. Building performs no
    /// file I/O or scientific calculation.
    ///
    /// # Errors
    ///
    /// Returns the native invalid-observer or non-finite Delta-T, tidal,
    /// lapse-rate or sidereal-setting error. Context settings not consumed by
    /// the builder retain their own checks at the consuming operation.
    pub fn build(&self) -> Result<Session, Error> {
        self.engine
            .build()
            .map(|engine| Session { engine })
            .map_err(|e| Error::from_engine("SessionBuilder::build", e))
    }
    pub(crate) fn from_context(context: crate::context::EphemerisContext) -> Self {
        Self {
            engine: crate::engine_api::SessionBuilder::from_context(context),
        }
    }
}
impl Session {
    pub(crate) fn from_engine(engine: crate::engine_api::Session) -> Self {
        Self { engine }
    }
}
impl Position {
    /// Component zero, in the actual requested frame/units; no conversion.
    pub fn longitude(&self) -> f64 {
        self.values[0]
    }
    /// Component one, in the actual requested frame/units; no conversion.
    pub fn latitude(&self) -> f64 {
        self.values[1]
    }
    /// Component two, in the actual requested frame/units; no conversion.
    pub fn distance(&self) -> f64 {
        self.values[2]
    }
}
impl Houses {
    /// Actual Ascendant word, in the producing method's angle units.
    pub fn ascendant(&self) -> f64 {
        self.angles[0]
    }
    /// Actual Midheaven word, in the producing method's angle units.
    pub fn mc(&self) -> f64 {
        self.angles[1]
    }
}
impl HousesWithSpeeds {
    /// Actual Ascendant word, in the producing method's angle units.
    pub fn ascendant(&self) -> f64 {
        self.angles[0]
    }
    /// Actual Midheaven word, in the producing method's angle units.
    pub fn mc(&self) -> f64 {
        self.angles[1]
    }
}
