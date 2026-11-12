//! Call-scoped assembly of the retained provider and borrowed serving callbacks.
//!
//! Provider selection remains lazy at the call boundary. The selected provider,
//! nutation closure, serving view and time policy live through the callback;
//! none is newly retained by the caller or moved into the numerical kernels.

use super::ServingStack;
use crate::context::{EphemerisContext, ProviderMode};
use crate::errors::{Error, Result};
use crate::positions::{self, SealedServing};
use crate::time::DeltaTModel;

/// Assemble the actual bound provider/serving/nutation/time policy and run one
/// synchronous borrowed callback. No state can escape. Non-LEB modes still
/// fail explicitly, rather than using geometric JPL words as apparent places.
pub fn with_sealed_serving_stack<'c, R>(
    ctx: &'c EphemerisContext,
    model: DeltaTModel,
    iers: bool,
    skip_closed: bool,
    f: impl FnOnce(&ServingStack<'_, 'c>) -> Result<R>,
) -> Result<R> {
    let mode = ctx.mode();
    if !matches!(mode, ProviderMode::SealedLeb | ProviderMode::Auto) {
        return Err(Error::ConfigurationMissing {
            missing_config: "non_leb_leg_wiring",
            suggestion: "the full non-LEB serving pipeline is not bound; no geometric-state substitution",
        });
    }
    if !skip_closed {
        positions::require_leb(!ctx.leb_files().is_empty())?;
    }
    let provider = ctx.leb_provider(skip_closed)?;
    let nutation_at = |jd| provider.nutation_at(jd);
    let serving = SealedServing::new_stored(provider, &nutation_at);
    let delta_cfg = ctx.delta_t_config(model, iers);
    f(&ServingStack {
        mode,
        serving: &serving,
        nutation_at: &nutation_at,
        stored_nutation_feed: true,
        delta_cfg: &delta_cfg,
    })
}
