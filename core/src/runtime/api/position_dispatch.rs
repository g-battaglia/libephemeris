//! Specialized stored/model routes selected before the ordinary body attempt.

use super::{
    PositionAcquisitions, PositionRequest, Session, curated_entry_with_anchors,
    model_entry_with_models_and_anchors, model_point_entry_with_anchors,
    true_point_entry_with_anchors,
};
use crate::constants::{
    CUPIDO, INTP_APOG, INTP_PERG, MEAN_APOG, MEAN_NODE, OSCU_APOG, TRUE_NODE, WALDEMATH,
};
use crate::context::ProviderMode;
use crate::errors::Error;
use crate::positions::{self, context_serving as serving};
use crate::sidereal::DynamicAnchorProvider;

impl Session {
    /// Select only the body families that own their complete serving route.
    ///
    /// Selection and serving stay ordered: curated sealed minors, hypothetical
    /// models, true lunar points, then mean/interpolated lunar points. Returning
    /// `None` does not acquire a serving stack or evaluate a model/anchor feed.
    pub(super) fn specialized_position(
        &self,
        request: PositionRequest,
        acquisitions: &PositionAcquisitions<'_>,
        live: Option<&dyn DynamicAnchorProvider>,
    ) -> Result<Option<([f64; 6], i32)>, Error> {
        let PositionRequest {
            epoch,
            body,
            input_flags,
        } = request;
        let (jd, is_ut) = epoch.into_entry();
        if self.ctx.mode() == ProviderMode::SealedLeb
            && crate::minor_bodies::minor_body_elements(positions::remap_ast_offset(body)).is_some()
        {
            return serving::with_sealed_serving_stack(
                &self.ctx,
                self.delta_t_model(),
                self.ctx.effective_iers_delta_t(),
                true,
                |stack| {
                    curated_entry_with_anchors(stack, &self.ctx, jd, body, input_flags, is_ut, live)
                },
            )
            .map(Some);
        }
        if (CUPIDO..=WALDEMATH).contains(&body) {
            return serving::with_sealed_serving_stack(
                &self.ctx,
                self.delta_t_model(),
                self.ctx.effective_iers_delta_t(),
                true,
                |stack| {
                    model_entry_with_models_and_anchors(
                        acquisitions.models,
                        stack,
                        &self.ctx,
                        jd,
                        body,
                        input_flags,
                        is_ut,
                        live,
                    )
                },
            )
            .map(Some);
        }
        if matches!(body, TRUE_NODE | OSCU_APOG) || body == -TRUE_NODE {
            return serving::with_sealed_serving_stack(
                &self.ctx,
                self.delta_t_model(),
                self.ctx.effective_iers_delta_t(),
                true,
                |stack| {
                    true_point_entry_with_anchors(
                        stack,
                        &self.ctx,
                        jd,
                        body,
                        input_flags,
                        is_ut,
                        live,
                    )
                },
            )
            .map(Some);
        }
        let source_body = crate::lunar::south_node_north_source(body).unwrap_or(body);
        if matches!(source_body, MEAN_NODE | MEAN_APOG | INTP_APOG | INTP_PERG) {
            return self
                .with_serving_stack(|stack| {
                    model_point_entry_with_anchors(
                        stack,
                        &self.ctx,
                        jd,
                        body,
                        input_flags,
                        is_ut,
                        live,
                    )
                })
                .map(Some);
        }
        Ok(None)
    }
}
