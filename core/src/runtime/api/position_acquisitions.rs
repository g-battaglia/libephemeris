//! Borrowed native model, catalog-anchor and numerical minor-body providers.
use super::Session;
use crate::errors::Error;
use crate::positions::AnalyticalModelCallback;
use crate::sidereal::{DynamicAnchor, DynamicAnchorProvider};
use crate::stars::AnchorRowAcquisition;

/// Optional synchronous native providers for a position calculation.
/// The caller retains ownership until the calculation returns. Assembly performs
/// no I/O; a provider is consulted only when the selected body or frame needs it.
#[derive(Default)]
pub struct PositionAcquisitions<'a> {
    /// Analytical orbital elements at the geometric model stage.
    pub models: Option<&'a AnalyticalModelCallback<'a>>,
    /// Catalog rows used by dynamic sidereal anchors.
    pub anchor_rows: Option<&'a dyn AnchorRowAcquisition>,
    /// Numerical integration after registered and automatic SPK sources decline.
    pub minor_assist: Option<&'a dyn super::MinorAssistAcquisition>,
}

/// Lazy independent Source computation provider for one live position call.
/// Its caller has already acquired observation TT. Forced LEB keeps scalar
/// vector queries; Auto uses the separately selected Native main on a miss.
/// Completed cached rows do not acquire an image. No static row follows error.
pub(super) struct PositionAnchors<'a> {
    pub session: &'a Session,
    pub rows: Option<&'a dyn AnchorRowAcquisition>,
}

impl DynamicAnchorProvider for PositionAnchors<'_> {
    fn has_analytic_anchor_rates(&self) -> bool {
        self.rows.is_none() && self.session.ctx.mode() == crate::context::ProviderMode::SealedLeb
    }

    fn anchor_rate_deg_day(
        &self,
        anchor: DynamicAnchor,
        tt: f64,
        noaberr: bool,
        nogdefl: bool,
    ) -> Result<Option<f64>, Error> {
        if self.rows.is_some() {
            return Ok(None);
        }
        if self.session.ctx.mode() != crate::context::ProviderMode::SealedLeb {
            return Ok(None);
        }
        crate::positions::star_serving::with_context_star_serving_captured_with_main(
            &self.session.ctx,
            &self.session.star_main,
            None,
            self.session.delta_t_model(),
            self.session.ctx.effective_iers_delta_t(),
            0,
            None,
            |serving, _| {
                serving
                    .anchor_provider()
                    .anchor_rate_deg_day(anchor, tt, noaberr, nogdefl)
            },
        )
        .result
    }

    fn anchor_longitude_deg(
        &self,
        anchor: DynamicAnchor,
        tt: f64,
        noaberr: bool,
        nogdefl: bool,
    ) -> Result<f64, Error> {
        if let Some(rows) = self.rows {
            return self
                .session
                .star_anchor_provider_captured(rows, |provider| {
                    provider.anchor_longitude_deg(anchor, tt, noaberr, nogdefl)
                })
                .result;
        }
        let native_context =
            (self.session.ctx.mode() == crate::context::ProviderMode::Auto).then(|| {
                self.session
                    .ctx
                    .clone()
                    .with_mode(crate::context::ProviderMode::NativeJpl)
            });
        crate::positions::star_serving::with_context_star_serving_captured_with_main(
            native_context.as_ref().unwrap_or(&self.session.ctx),
            &self.session.star_main,
            None,
            self.session.delta_t_model(),
            self.session.ctx.effective_iers_delta_t(),
            0,
            None,
            |serving, _| {
                serving
                    .anchor_provider()
                    .anchor_longitude_deg(anchor, tt, noaberr, nogdefl)
            },
        )
        .result
    }
}

/// Dynamic fast-LEB helper's independent observation clock. The stored body
/// and nutation keep their model TT; only the delegated live anchor receives
/// global TT from `fast_TT - configured_delta(fast_TT)`.
#[derive(Clone, Copy)]
pub(super) struct AcquiredAnchorClock {
    pub tt: f64,
    pub ut1: f64,
    pub global_tt: f64,
    pub configured_delta_rate: f64,
}

pub(super) struct FastPositionAnchors<'a> {
    pub live: &'a dyn DynamicAnchorProvider,
    pub delta: &'a dyn Fn(f64) -> Result<f64, Error>,
    pub delta_value_rate: &'a dyn Fn(f64) -> Result<(f64, f64), Error>,
    pub want_rate: bool,
    pub clock: std::cell::Cell<Option<AcquiredAnchorClock>>,
}

impl DynamicAnchorProvider for FastPositionAnchors<'_> {
    fn has_analytic_anchor_rates(&self) -> bool {
        self.live.has_analytic_anchor_rates()
    }

    fn anchor_rate_deg_day(
        &self,
        anchor: DynamicAnchor,
        tt: f64,
        noaberr: bool,
        nogdefl: bool,
    ) -> Result<Option<f64>, Error> {
        if !self.live.has_analytic_anchor_rates() {
            return Ok(None);
        }
        // The preceding central scalar anchor acquired these exact clock words
        // and the derivative of that same selected source. No second lookup.
        let Some(clock) = self.clock.get().filter(|v| v.tt.to_bits() == tt.to_bits()) else {
            return Err(Error::ConfigurationMissing {
                missing_config: "anchor_clock_rate",
                suggestion: "analytic anchor rates consume their central native clock acquisition",
            });
        };
        let Some(rate) =
            self.live
                .anchor_rate_deg_day(anchor, clock.global_tt, noaberr, nogdefl)?
        else {
            return Ok(None);
        };
        let clock_rate = crate::smh2016::timescale_tt_rate_from_ut1(clock.ut1)?
            * (1.0 - clock.configured_delta_rate);
        Ok(Some(rate * clock_rate))
    }
    fn anchor_longitude_deg(
        &self,
        anchor: DynamicAnchor,
        tt: f64,
        noaberr: bool,
        nogdefl: bool,
    ) -> Result<f64, Error> {
        let (delta, delta_rate) = if self.want_rate && self.live.has_analytic_anchor_rates() {
            (self.delta_value_rate)(tt)?
        } else {
            ((self.delta)(tt)?, 0.0)
        };
        let ut = tt - delta;
        let observation_tt = crate::time::default_timescale_tt_from_ut1(ut)?;
        self.clock.set(Some(AcquiredAnchorClock {
            tt,
            ut1: ut,
            global_tt: observation_tt,
            configured_delta_rate: delta_rate,
        }));
        self.live
            .anchor_longitude_deg(anchor, observation_tt, noaberr, nogdefl)
    }
}
