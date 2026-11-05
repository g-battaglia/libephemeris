//! Native provider fallback after a stored-data lookup declines.

use super::*;
use crate::context::ProviderMode;
use crate::positions::MainImageAcquisition;

impl Session {
    /// Evaluate the physical provider chain after a stored-data miss.
    /// Native images remain owned by this context; a remote cache hit acquires no local image.
    pub(super) fn position_after_fast(
        &self,
        request: PositionRequest,
        acquisitions: &PositionAcquisitions<'_>,
        live: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
        main_image: Option<&dyn MainImageAcquisition>,
        owned_minor: Option<&super::NativeMinorAssist<'_>>,
    ) -> Result<([f64; 6], i32), Error> {
        let PositionRequest {
            epoch,
            body,
            input_flags,
        } = request;
        // Keep the public scalar untouched; the selected provider owns clocks.
        let (jd, is_ut) = epoch.into_entry();
        if self.ctx.mode() == ProviderMode::SealedLeb {
            return self.with_serving_stack(|stack| {
                if crate::lunar::is_lunar_point(body) {
                    return lunar_point_route_with_anchors_prepared(
                        stack,
                        &self.ctx,
                        jd,
                        body,
                        input_flags,
                        is_ut,
                        false,
                        live,
                        None,
                    );
                }
                if crate::minor_bodies::minor_body_elements(positions::remap_ast_offset(body))
                    .is_some()
                    || (CUPIDO..=WALDEMATH).contains(&body)
                {
                    return model_entry_with_models_and_anchors(
                        acquisitions.models,
                        stack,
                        &self.ctx,
                        jd,
                        body,
                        input_flags,
                        is_ut,
                        live,
                    );
                }
                positions::calc_stored_vector_retry(
                    stack,
                    &self.ctx,
                    jd,
                    body,
                    input_flags,
                    is_ut,
                    live,
                )
            });
        }
        let native_ctx = self.ctx.clone().with_mode(ProviderMode::NativeJpl);
        let native_request = positions::ServingRequest {
            ctx: &native_ctx,
            delta_t_model: self.delta_t_model(),
            iers_enabled: self.ctx.effective_iers_delta_t(),
            observer: self.ctx.observer(),
            tier_policy: positions::TierPolicy::AttachOrder,
            skip_closed: false,
        };
        let native = |epoch: f64, target: i32, request_flags: i32, ut_entry: bool| {
            if let Some(owned) = owned_minor {
                return owned.calc_native_position_with_main_image(
                    acquisitions.models,
                    live,
                    native_request,
                    epoch,
                    target,
                    request_flags,
                    ut_entry,
                    main_image,
                );
            }
            positions::calc_native_jpl_entry_with_main_image(
                acquisitions.models,
                live,
                acquisitions.minor_assist,
                native_request,
                epoch,
                target,
                request_flags,
                ut_entry,
                main_image,
            )
        };
        let remote = self.ctx.mode() == ProviderMode::HorizonsRemote
            || (self.ctx.mode() == ProviderMode::Auto
                && !selected_local_main(&self.ctx, main_image));
        if !remote {
            return native(jd, body, input_flags, is_ut);
        }
        let remote_ctx = if self.ctx.horizons_client().is_some() {
            self.ctx.clone()
        } else {
            self.ctx.clone().with_horizons_config(Default::default())?
        };
        let remote_request = positions::ServingRequest {
            ctx: &remote_ctx,
            delta_t_model: self.delta_t_model(),
            iers_enabled: self.ctx.effective_iers_delta_t(),
            observer: self.ctx.observer(),
            tier_policy: positions::TierPolicy::AttachOrder,
            skip_closed: false,
        };
        positions::calc_horizons_with_fallback(
            remote_request,
            jd,
            body,
            input_flags,
            is_ut,
            &native,
        )
    }
}

/// Whether the selected local main is present or already retained.
/// A partial minor-body registry does not establish main-planet coverage.
fn selected_local_main(
    ctx: &EphemerisContext,
    main_image: Option<&dyn MainImageAcquisition>,
) -> bool {
    let filename =
        positions::select_de_kernel(ctx.ephemeris_file_override(), None, None, ctx.tier());
    let file = std::path::Path::new(filename);
    let present = if file.is_absolute() {
        file.is_file()
    } else {
        [ctx.jpl_data_dir(), ctx.paths().jpl_dir.as_ref()]
            .into_iter()
            .flatten()
            .any(|directory| directory.join(file).is_file())
    };
    present || main_image.is_some_and(|acquire| acquire.has_current(ctx))
}
