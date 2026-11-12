//! Hypothetical-family reduction over the entry's borrowed model and clock feeds.
//!
//! No support state or observer offset is prepared up front. Each branch reads
//! only its original dependencies, so declined topocentric models can continue
//! through ordinary JPL dispatch and support failures retain their precedence.

use super::super::{
    PublicModelDirectory, hypothetical_topo_offset_j2000, hypothetical_topo_site,
    local_jpl_symbolic_center_state_via, model_serving, native_model_ut_at,
};
use super::{
    local_jpl_predicted_center_state_with_models, local_jpl_uranian_bary_state_with_models,
};
use crate::constants::{
    CUPIDO, FLG_BARYCTR, FLG_HELCTR, FLG_TOPOCTR, ISIS, NIBIRU, POSEIDON, WALDEMATH, WHITE_MOON,
};
use crate::errors::{Error, Result};
use crate::positions::{AnalyticalModelCallback, PositionResult, strip_output_flags};
use crate::sidereal::{SealedAyanamshaServing, SiderealSelection};

/// Only the borrowed geometric model and frame-clock dependencies of this
/// family. Resource acquisition, range-error promotion and source publication
/// stay with the entry; this view owns no ephemeris, callback or prepared state.
pub(super) struct HypotheticalModels<'a, 'registry> {
    pub(super) models: Option<&'a AnalyticalModelCallback<'a>>,
    pub(super) directory: &'a PublicModelDirectory<'a, 'registry>,
    pub(super) original_ut: Option<f64>,
    pub(super) model_deltat: &'a dyn Fn(f64) -> Result<f64>,
    pub(super) ayanamsha: &'a SealedAyanamshaServing<'a>,
    pub(super) selection: SiderealSelection,
}

impl HypotheticalModels<'_, '_> {
    /// Return a reduced model only when its original branch serves the request.
    /// `None` preserves fallthrough, including topocentric cases without a site;
    /// it is not a fallback for a failed support read or numerical calculation.
    pub(super) fn position(
        &self,
        tt: f64,
        body: i32,
        flags: i32,
        site: Option<(f64, f64, f64)>,
    ) -> Result<Option<PositionResult>> {
        let Self {
            models,
            directory,
            original_ut,
            model_deltat,
            ayanamsha,
            selection,
        } = *self;
        if (CUPIDO..=POSEIDON).contains(&body)
            && flags & FLG_BARYCTR != 0
            && flags & FLG_TOPOCTR == 0
        {
            let state = local_jpl_uranian_bary_state_with_models(models, directory, tt, body)?;
            let aya = |u, f| crate::planetary_moons::moon_spk_ayanamsha(u, f, ayanamsha);
            let output = model_serving::finish_uranian_bary_state(
                state,
                &|sample_tt| native_model_ut_at(sample_tt, original_ut),
                tt,
                strip_output_flags(flags),
                &aya,
            )?;
            return Ok(Some(output));
        }
        if (CUPIDO..=ISIS).contains(&body) {
            let aya = |u, f| crate::planetary_moons::moon_spk_ayanamsha(u, f, ayanamsha);
            let state_at = |id, epoch: f64| {
                // Keep separate TT day and fraction for each support sample.
                let whole = epoch.floor();
                directory.state(id, whole, epoch - whole)
            };
            let output = model_serving::classic_model_state_source_time_with_models(
                models,
                &state_at,
                &|epoch| directory.earth(epoch),
                &|sample_tt| native_model_ut_at(sample_tt, original_ut),
                tt,
                body,
                strip_output_flags(flags),
                site,
                &aya,
            )?;
            return Ok(Some(output));
        }
        if body == WHITE_MOON || body == WALDEMATH {
            // Recenter the symbolic geocentric model, then apply nutation and output
            // frame shaping. Observer parallax is subtracted in the frame of date before
            // any sidereal rotation or requested J2000 precession.
            let topo_site = hypothetical_topo_site(flags, site)?;
            let is_geocentric_topo =
                flags & FLG_TOPOCTR != 0 && flags & (FLG_HELCTR | FLG_BARYCTR) == 0;
            if !is_geocentric_topo || topo_site.is_some() {
                let state = local_jpl_symbolic_center_state_via(
                    tt,
                    body,
                    flags,
                    model_deltat,
                    &|id, day, fraction| directory.state(id, day, fraction),
                    &|epoch| directory.earth(epoch),
                )?;
                let aya = |u, f| crate::planetary_moons::moon_spk_ayanamsha(u, f, ayanamsha);
                let numerics = crate::positions::horizons_output::HorizonsNumerics::new(
                    model_deltat,
                    &aya,
                    selection.full().0,
                    crate::lunar::LunarCoverage::Unenforced,
                );
                let offset = topo_site
                    .map(|site| hypothetical_topo_offset_j2000(tt, site))
                    .transpose()?;
                let (output, _) = numerics.hypothetical_source_time_state_using_ut(
                    &|| native_model_ut_at(tt, original_ut),
                    tt,
                    body,
                    strip_output_flags(flags),
                    state,
                    offset,
                )?;
                return Ok(Some(output));
            }
        }
        if body == NIBIRU {
            // Acquire Earth-from-Sun support for the geocentric model before its center reduction.
            if flags & (FLG_HELCTR | FLG_BARYCTR) == 0 {
                directory.earth(tt)?;
            }
            return Err(Error::UnknownBody { body_id: body });
        }
        if crate::horizons::is_predicted_hypothetical_body(body) {
            // Recenter the retarded heliocentric model with native PPN and annual
            // aberration, then apply the shared hypothetical-body frame reduction.
            // Observer parallax precedes sidereal rotation and J2000 output precession.
            let topo_site = hypothetical_topo_site(flags, site)?;
            let is_geocentric_topo =
                flags & FLG_TOPOCTR != 0 && flags & (FLG_HELCTR | FLG_BARYCTR) == 0;
            if !is_geocentric_topo || topo_site.is_some() {
                let (state, rate_origin) = local_jpl_predicted_center_state_with_models(
                    models,
                    directory.leg,
                    directory,
                    tt,
                    body,
                    flags,
                )?;
                let aya = |u, f| crate::planetary_moons::moon_spk_ayanamsha(u, f, ayanamsha);
                let numerics = crate::positions::horizons_output::HorizonsNumerics::new(
                    model_deltat,
                    &aya,
                    selection.full().0,
                    crate::lunar::LunarCoverage::Unenforced,
                );
                let offset = topo_site
                    .map(|site| hypothetical_topo_offset_j2000(tt, site))
                    .transpose()?;
                let (output, _) = numerics.hypothetical_source_time_state_with_rate_origin(
                    &|| native_model_ut_at(tt, original_ut),
                    tt,
                    body,
                    strip_output_flags(flags),
                    state,
                    offset,
                    rate_origin,
                )?;
                return Ok(Some(output));
            }
        }
        Ok(None)
    }
}
