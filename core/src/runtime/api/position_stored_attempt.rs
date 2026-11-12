//! One ordinary stored-data attempt, without retry policy or result packaging.

use super::{
    FastPositionAnchors, PositionAcquisitions, PositionRequest, Session, earth_origin_trace,
    opaque_live_body_chain, stored_state_with_clock_at, tropical,
};
use crate::constants::{EARTH, ECL_NUT, FLG_SIDEREAL, FLG_SPEED};
use crate::context::ProviderMode;
use crate::errors::Error;
use crate::houses::sidereal_selection_for_context;
use crate::positions::{self, context_serving as serving};
use crate::sidereal::{DynamicAnchorProvider, resolve_sidereal_route, strip_to_base_mode};

#[cfg(test)]
use super::super::owned_apparent_session_tests::record_owned_apparent_admission;

impl Session {
    /// Attempt stored body serving after raw flag admission and the diagnostic
    /// snapshot. The coordinator owns early routes, retries and final packaging.
    ///
    /// The live anchor owner stays borrowed from the coordinator through any
    /// later fallback. This attempt's fast clock adapter contains only borrowed
    /// feeds and a Copy clock cell; retiring it neither releases a source nor
    /// invokes a callback. Derivative serving views remain inside the serving
    /// stack's closure, including all observer, clock and projection operations.
    ///
    /// `main_image_supplied` is observed only by the admission test probe. The
    /// actual main-image capability is not consumed by the stored-data route.
    pub(super) fn ordinary_stored_position_attempt(
        &self,
        request: PositionRequest,
        acquisitions: &PositionAcquisitions<'_>,
        live: Option<&dyn DynamicAnchorProvider>,
        main_image_supplied: bool,
    ) -> Result<([f64; 6], i32), Error> {
        let PositionRequest {
            epoch,
            body,
            input_flags,
        } = request;
        let (jd, is_ut) = epoch.into_entry();
        #[cfg(not(test))]
        let _ = main_image_supplied;

        let (selection, _) = sidereal_selection_for_context(&self.ctx);
        let route = resolve_sidereal_route(
            input_flags,
            strip_to_base_mode(self.ctx.sidereal_mode()),
            selection.bits(),
        );
        let delta = |ut| self.deltat(ut);
        let delta_value_rate = |tt| {
            self.ctx.delta_t_value_and_rate_days(
                tt,
                self.delta_t_model(),
                self.ctx.effective_iers_delta_t(),
            )
        };
        let fast_owner = live.map(|live| FastPositionAnchors {
            live,
            delta: &delta,
            delta_value_rate: &delta_value_rate,
            want_rate: input_flags & (FLG_SPEED | crate::constants::FLG_SPEED3) != 0,
            clock: std::cell::Cell::new(None),
        });
        let fast_live = fast_owner
            .as_ref()
            .map(|owner| owner as &dyn crate::sidereal::DynamicAnchorProvider);
        let owned_apparent_id = (self.ctx.mode() == ProviderMode::SealedLeb
            && !self.ctx.leb_files().is_empty()
            && self
                .ctx
                .leb_files()
                .iter()
                .all(crate::leb::LebFile::is_owned_image)
            && acquisitions.anchor_rows.is_none()
            && acquisitions.models.is_none())
        .then_some(self.owned_apparent_id)
        .flatten();
        self.with_serving_stack(|stack| {
            // The existing constructor has already admitted the actual native
            // provider and strict handle state. Opaque models/anchors never opt in.
            // Main-image suppliers are consumed only by the separate native
            // continuation, never by this strict stored physical core.
            #[cfg(test)]
            record_owned_apparent_admission(main_image_supplied, owned_apparent_id);
            let owned_serving = stack.serving.with_owned_apparent_token(owned_apparent_id);
            let owned_stack = stack.with_derivative_serving(&owned_serving);
            let stack = &owned_stack;
            let sampled_serving = stack.serving.with_sampled_derivatives();
            let sampled_stack = stack.with_derivative_serving(&sampled_serving);
            let stack = if opaque_live_body_chain(&self.ctx, input_flags, live) {
                &sampled_stack
            } else {
                stack
            };
            // Nutation returns before body serving; topocentric Earth is the
            // exact zero origin. Neither early domain consumes an observer.
            if !matches!(body, ECL_NUT | EARTH) {
                let site = self
                    .ctx
                    .observer()
                    .map(|o| (o.lon_deg(), o.lat_deg(), o.height_m()));
                positions::resolve_topo_request(
                    positions::normalize_calc_flags(input_flags, false),
                    positions::remap_ast_offset(body),
                    site,
                    None,
                )?;
            }
            let site_clock = |tt| {
                positions::configured_site_ut1_partner_with_table(
                    tt,
                    is_ut.then_some(jd),
                    stack.delta_cfg,
                    &|epoch| stack.serving.provider().delta_t(epoch),
                )
            };
            let raw = |tt, id, calc_flags| {
                tropical(stack, &self.ctx, tt, id, calc_flags, &site_clock, jd)
            };
            if route.is_some() {
                return earth_origin_trace(body, input_flags, || {
                    serving::with_fixed_epoch_engine_with_native_motion(
                        stack,
                        &self.ctx,
                        &raw,
                        fast_live,
                        live,
                        |engine| {
                            if is_ut {
                                positions::calc_ut(engine, jd, body, input_flags)
                            } else {
                                positions::calc(engine, jd, body, input_flags)
                            }
                        },
                    )
                });
            }
            // Preserve the normal entry's time/flag normalization, with lunar
            // dispatch supplied as its tropical TT core. XYZ receives direct
            // stored Cartesian words for ordinary bodies, not a roundtrip.
            let engine = positions::CalcEngine {
                ut1_to_tt_jd: &|ut| Ok(ut + self.deltat(ut)?),
                compute_tt: &raw,
                nutation_obliquity_ut: &positions::ecl_nut_result_ut,
                nutation_obliquity_tt: &|tt, nut_flags| {
                    positions::ecl_nut_result(tt, nut_flags, &crate::erfa::nut06a)
                },
                jpl_fallback_active: false,
            };
            let tropical_compute_flags = input_flags & !FLG_SIDEREAL;
            let cartesian = |tt, id, calc_flags| {
                if crate::lunar::is_lunar_point(id) {
                    return Ok(None);
                }
                match positions::serving_channel_gate(stack.serving.provider(), id) {
                    Err(Error::SealedChannelNotServed { .. }) => Ok(None),
                    Err(error) => Err(positions::escalate_sealed_provider_miss(
                        self.ctx.mode(),
                        stack.serving.provider(),
                        id,
                        tt,
                        error,
                    )),
                    Ok(_) => stored_state_with_clock_at(
                        stack,
                        &self.ctx,
                        tt,
                        id,
                        calc_flags,
                        true,
                        &site_clock,
                        Some(jd),
                    )
                    .map(Some),
                }
            };
            let (values, tropical_returned_flags) = earth_origin_trace(body, input_flags, || {
                if is_ut {
                    positions::calc_ut_with_cartesian(
                        &engine,
                        &cartesian,
                        jd,
                        body,
                        tropical_compute_flags,
                    )
                } else {
                    positions::calc_with_cartesian(
                        &engine,
                        &cartesian,
                        jd,
                        body,
                        tropical_compute_flags,
                    )
                }
            })?;
            // SIDEREAL was masked for the tropical calculation only; restore
            // its implied NONUT echo without changing the computation flags.
            Ok((
                values,
                tropical_returned_flags
                    | (input_flags & FLG_SIDEREAL)
                    | positions::implied_retflag_bits(input_flags),
            ))
        })
    }
}
