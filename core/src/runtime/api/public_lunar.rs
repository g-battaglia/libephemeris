//! Scalar-TT osculating lunar positions over the selected vector ephemeris.

use super::{Error, Session};
use crate::constants::{EARTH, MOON};
use crate::context::ProviderMode;
use crate::positions::StateProvider;

/// Geocentric ecliptic-of-date position in AU and velocity in AU/TT-day.
pub(super) type LunarOsculatingState = ((f64, f64, f64), (f64, f64, f64));

impl Session {
    /// Osculating lunar perigee from the actual eccentricity vector.
    /// Input is TT Julian day; return true ecliptic-of-date longitude and
    /// latitude in degrees and perigee range in AU, p/(1+e). The state and
    /// provider domain are shared with [`Self::true_lunar_point`]. Singular
    /// conic geometry, invalid time and resource errors propagate.
    pub fn calc_osculating_perigee(&self, jd_tt: f64) -> Result<(f64, f64, f64), Error> {
        let whole = jd_tt.floor();
        let (position, velocity) = self.osculating_lunar_state_parts(whole, jd_tt - whole)?;
        crate::lunar::osculating_perigee_from_state(position, velocity)
    }

    /// Compare the two equivalent osculating-apogee formulas at TT JD.
    /// Positions are degrees/degrees/AU; differences retain those units.
    /// Both formulas use a single acquired geometric lunar state.
    pub fn compare_true_lilith_methods(
        &self,
        jd_tt: f64,
    ) -> Result<crate::lunar::TrueLilithComparison, Error> {
        crate::lunar::compare_true_lilith_methods_via(
            &|tt| {
                let whole = tt.floor();
                self.osculating_lunar_state_parts(whole, tt - whole)
            },
            jd_tt,
        )
    }

    /// Mean lunar ascending-node position and daily rates at TT JD.
    /// Coordinates are degrees/degrees/AU. SPEED3 differentiates the same
    /// analytic model; active main-kernel range and invalid time are checked.
    pub fn calc_mean_lunar_node_state(
        &self,
        jd_tt: f64,
        speed3: bool,
    ) -> Result<(f64, f64, f64, f64, f64, f64), Error> {
        crate::lunar::calc_mean_lunar_node_state(
            jd_tt,
            self.mean_lunar_coverage(jd_tt, crate::constants::MEAN_NODE)?,
            speed3,
        )
    }

    /// Mean lunar apogee position and daily rates at TT JD.
    /// Units and range policy follow [`Self::calc_mean_lunar_node_state`].
    pub fn calc_mean_lilith_state(
        &self,
        jd_tt: f64,
        speed3: bool,
    ) -> Result<(f64, f64, f64, f64, f64, f64), Error> {
        crate::lunar::calc_mean_lilith_state(
            jd_tt,
            self.mean_lunar_coverage(jd_tt, crate::constants::MEAN_APOG)?,
            speed3,
        )
    }

    /// Osculating ascending node or apogee at a scalar TT Julian day.
    ///
    /// `apogee == false` selects the ascending node (zero latitude), otherwise
    /// the eccentricity-vector apogee. Return longitude/latitude in degrees
    /// on the true ecliptic of date and conic distance in AU. The caller first
    /// acquires the computation ephemeris; this call retains its immutable
    /// LEB/SPK owners. Earth is selected before Moon. Native SPK keeps scalar
    /// TT's independent split TDB sampling and vector frame; LEB evaluates
    /// both geometric states at scalar TT and rotates their difference by the
    /// same vector frame. No configured Delta-T, public position flag, speed
    /// stencil, alternate provider or stored nutation table is substituted.
    ///
    /// Missing aliases, selected-reader coverage, raw SPK ranges, non-finite
    /// epochs, frame faults and osculating-state failures propagate. This
    /// method has no public scalar-range envelope or network fallback. Python
    /// epoch protocols and backend object/error reconstruction belong to the
    /// binding at their actual acquisition stages.
    pub fn true_lunar_point(&self, jd_tt: f64, apogee: bool) -> Result<(f64, f64, f64), Error> {
        let whole = jd_tt.floor();
        self.true_lunar_point_parts(whole, jd_tt - whole, apogee)
    }

    /// Osculating lunar position retaining a scalar Time's original TT split.
    ///
    /// Output, selected providers and errors follow [`Self::true_lunar_point`].
    /// Native SPK reads retain `whole` and `tt_fraction`; frames and LEB readers
    /// consume their scalar sum. No resplitting, model clock or public flag
    /// conversion changes the acquired fraction.
    pub fn true_lunar_point_parts(
        &self,
        whole: f64,
        tt_fraction: f64,
        apogee: bool,
    ) -> Result<(f64, f64, f64), Error> {
        let state = self.osculating_lunar_state_parts(whole, tt_fraction)?;
        if apogee {
            crate::lunar::true_lilith_from_state(state.0, state.1)
        } else {
            crate::lunar::true_lunar_node_from_state(state.0, state.1)
        }
    }
    /// Acquire the owned geocentric lunar state once for osculating models.
    /// TT parts are Julian days. Native SPK uses the Session's shared immutable
    /// main selection, preserving explicit registrations and active snapshots.
    /// LEB uses the same configured Earth/Moon coefficients and date frame.
    /// Missing data, coverage and frame errors propagate without fallback.
    pub(super) fn osculating_lunar_state_parts(
        &self,
        whole: f64,
        tt_fraction: f64,
    ) -> Result<LunarOsculatingState, Error> {
        let jd_tt = whole + tt_fraction;
        let state = if self.ctx.mode() == ProviderMode::NativeJpl {
            let selected = self.star_main.snapshot(&self.ctx, None)?;
            crate::positions::LocalSpkStateProvider::new(selected.spk_registry(), selected.mode())?
                .public_true_moon_state_parts(whole, tt_fraction)?
        } else {
            crate::positions::context_serving::with_sealed_serving_stack(
                &self.ctx,
                self.delta_t_model(),
                self.ctx.effective_iers_delta_t(),
                false,
                |stack| {
                    let provider = stack.serving.provider();
                    for (body, name) in [(EARTH, "earth"), (MOON, "moon")] {
                        if !provider.has_body(body) {
                            return Err(missing_vector_alias(body, name));
                        }
                    }
                    let relative = provider
                        .relative_geometric_state(MOON, EARTH, jd_tt)
                        .map_err(|error| provider.vector_error(error))?;
                    let frame = crate::planetary_moons::moon_spk_frame(jd_tt)?;
                    let position =
                        crate::positions::reference_frames::matvec3(frame, relative.pos())?;
                    let velocity =
                        crate::positions::reference_frames::matvec3(frame, relative.vel())?;
                    Ok((
                        (position[0], position[1], position[2]),
                        (velocity[0], velocity[1], velocity[2]),
                    ))
                },
            )?
        };
        Ok(state)
    }
}

/// A recognized missing reader alias, after the actual owned membership check.
/// The public `from None` suppresses the caught KeyError while retaining that
/// same key/error as context. No provider, date or coefficient is substituted.
fn missing_vector_alias(body: i32, _name: &str) -> Error {
    Error::UnknownBody { body_id: body }
}
