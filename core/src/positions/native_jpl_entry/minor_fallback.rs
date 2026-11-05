//! Call-scoped minor-body policy after registered/refreshed SPK coverage declines.
//!
//! SPK/main-image acquisition and the SBDB cache remain in the outer entry.
//! These inputs only borrow its resources and callbacks; they neither retain
//! owners nor retry SPK. The ASSIST Session lives only inside its lazy attempt.

use super::super::{
    local_jpl_keplerian_fallback, local_jpl_lunar_frame, local_jpl_numbered_asteroid_fallback,
    local_jpl_sidereal_word, native_source_result,
};
use super::optional_assist_result;
use crate::constants::{FLG_EQUATORIAL, FLG_J2000, FLG_SIDEREAL, FLG_TRUEPOS};
use crate::context::EphemerisContext;
use crate::engine_api::{MinorAssistAcquisition, NativeMinorAssist, Session};
use crate::errors::{Error, Result};
use crate::minor_bodies::SbdbElementsProvider;
use crate::nbody::OwnedTrajectoryError;
use crate::positions::{JplLeg, PositionResult};

/// Borrowed inputs for one fallback evaluation, never a new resource owner.
/// The concrete native supplier remains distinct from the public callback:
/// only it can evaluate an owned trajectory at physical emission epochs.
pub(super) struct MinorBodyFallback<'a, 'registry> {
    pub(super) leg: &'a JplLeg<'registry>,
    pub(super) context: &'a EphemerisContext,
    pub(super) minor_assist: Option<&'a dyn MinorAssistAcquisition>,
    pub(super) owned_minor: Option<&'a NativeMinorAssist<'a>>,
    pub(super) original_ut: Option<f64>,
    pub(super) ayanamsha_at: &'a dyn Fn(f64) -> Result<f64>,
    pub(super) sbdb_elements: SbdbElementsProvider<'a>,
}

impl MinorBodyFallback<'_, '_> {
    /// Enforce the configured precision gate before trying ASSIST, then use
    /// curated or SBDB Keplerian elements only if optional acquisition declines.
    /// TOPOCTR is still echoed for Keplerian output without inventing parallax.
    pub(super) fn position(&self, tt: f64, body: i32, flags: i32) -> Result<(PositionResult, i32)> {
        self.require_precision(tt, body)?;
        if let Some(supplier) = self.minor_assist {
            let attempt = self.assist_position(supplier, tt, body, flags);
            if let Some(output) = optional_assist_result(attempt)? {
                return native_source_result(Ok((output, flags)), "ASSIST");
            }
        }
        let mut output = if crate::minor_bodies::minor_body_elements(body).is_some() {
            local_jpl_keplerian_fallback(self.leg, tt, body, flags)?
        } else {
            local_jpl_numbered_asteroid_fallback(self.leg, tt, body, flags, self.sbdb_elements)?
        };
        local_jpl_sidereal_word(&mut output, tt, self.original_ut, flags, self.ayanamsha_at)?;
        local_jpl_lunar_frame(&mut output, tt, flags)?;
        native_source_result(Ok((output, flags)), "Keplerian")
    }

    /// Availability is deliberately not memoized here: strict-policy probing
    /// and the subsequent ASSIST attempt are separate public callback events.
    /// A registered coverage miss must retain its more specific strict error.
    fn require_precision(&self, tt: f64, body: i32) -> Result<()> {
        if self.context.effective_strict_precision()
            && !self.context.effective_auto_spk_download()
            && crate::constants::get_spk_body_info_from_map(body).is_some()
            && !crate::constants::is_spk_auto_download_blocked(body)
        {
            let available = if let Some(supplier) = self.minor_assist {
                if supplier.contains(body)? {
                    match supplier.available() {
                        Ok(available) => available,
                        Err(error) if crate::nbody::assist_fallback(&error) => false,
                        Err(error) => return Err(error),
                    }
                } else {
                    false
                }
            } else {
                false
            };
            if !available {
                if self.leg.registry().body_info(body).is_some() {
                    return Err(Error::SpkRegisteredOutOfCoverage {
                        body,
                        requested_jd: tt,
                    });
                }
                return Err(Error::SpkRequired { body });
            }
        }
        Ok(())
    }

    /// Keep the Session inside the actual ASSIST attempt. A concrete acquired
    /// orbit's emission or frame failure cannot become a Keplerian fallback;
    /// opaque suppliers retain their initial-acquisition error classification.
    fn assist_position(
        &self,
        supplier: &dyn MinorAssistAcquisition,
        tt: f64,
        body: i32,
        flags: i32,
    ) -> std::result::Result<Option<PositionResult>, OwnedTrajectoryError> {
        use OwnedTrajectoryError::Initial;
        let curved = self.owned_minor.filter(|_| flags & FLG_TRUEPOS == 0);
        if !supplier.contains(body).map_err(Initial)? || !supplier.available().map_err(Initial)? {
            return Ok(None);
        }
        let session = Session::from_context(self.context.clone());
        let equatorial = flags & FLG_EQUATORIAL != 0;
        // Held-PV suppliers include the inertial frame in their analytic
        // derivative. The concrete owned orbit keeps its existing curved-ray
        // projection and final-frame ordering.
        let inertial_ecliptic =
            curved.is_none() && !equatorial && flags & FLG_J2000 != 0 && flags & FLG_SIDEREAL == 0;
        let mut output = if let Some(owned) = curved {
            session.assist_minor_owned_position(owned, body, tt, flags)?
        } else if equatorial {
            session
                .assist_minor_equatorial(supplier, body, tt, flags)
                .map_err(Initial)?
        } else if inertial_ecliptic {
            session
                .assist_minor_inertial_ecliptic(supplier, body, tt, flags)
                .map_err(Initial)?
        } else {
            session
                .assist_minor_ecliptic(supplier, body, tt, flags)
                .map_err(Initial)?
        };
        let frame_error = |error| {
            if curved.is_some() {
                OwnedTrajectoryError::after_acquisition(error)
            } else {
                Initial(error)
            }
        };
        local_jpl_sidereal_word(&mut output, tt, self.original_ut, flags, self.ayanamsha_at)
            .map_err(frame_error)?;
        if !equatorial && !inertial_ecliptic {
            local_jpl_lunar_frame(&mut output, tt, flags).map_err(frame_error)?;
        }
        Ok(Some(output))
    }
}
