//! Apparent satellite places over independently retained satellite/main images.

use super::{Error, Session};
use crate::planetary_moons::MoonSpkRegistry;

impl Session {
    /// Evaluate a registered satellite at scalar TT Julian day.
    ///
    /// The satellite registry and this Session's main observer registry remain
    /// separately owned throughout the calculation. No file acquisition or
    /// registration occurs here. `moon_id` is a supported legacy or canonical
    /// satellite id; `flags` control observer, light time, aberration, nutation,
    /// sidereal reduction and one-second central-stencil speed calculation.
    ///
    /// Returns ecliptic-of-date longitude/latitude in degrees, radius in AU and
    /// their daily rates, or None for an unknown/unregistered target. Explicit
    /// sealed mode refuses JPL. Non-finite TT, missing main/site settings, frame
    /// and SPK coverage/evaluation errors propagate from the satellite kernel.
    /// The Session's configured Delta-T affects model-defined sidereal epochs;
    /// target TT, split TDB and frame/rotation clocks remain independent.
    pub fn planetary_moon_native_tt(
        &self,
        moons: &MoonSpkRegistry,
        moon_id: i32,
        jd_tt: f64,
        flags: i32,
    ) -> Result<Option<[f64; 6]>, Error> {
        let delta = self
            .ctx
            .delta_t_config(self.delta_t_model(), self.ctx.effective_iers_delta_t());
        crate::planetary_moons::moon_position_native_tt(
            &self.ctx, moons, moon_id, jd_tt, flags, &delta,
        )
    }
}
