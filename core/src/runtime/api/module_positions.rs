//! Typed direct provider helpers for useful numeric module interfaces.

use super::{Error, Session};
use crate::horizons::{HorizonsClient, HorizonsTimeType, StateVector};
use std::num::NonZeroUsize;

/// Construct a bounded Horizons vector client without opening a connection.
///
/// `timeout` is seconds; cache and worker capacities must be positive. The
/// cache owns plain AU/AU-day states and each fetch enforces the caller policy.
/// Invalid capacities or timeout are ordinary native configuration errors.
pub fn new_horizons_client(
    max_cache_size: usize,
    timeout: f64,
    max_workers: usize,
) -> Result<HorizonsClient, Error> {
    let positive = |value, name| {
        NonZeroUsize::new(value).ok_or_else(|| Error::EphemerisLoaderValue {
            message: format!("{name} must be positive"),
        })
    };
    HorizonsClient::with_config(
        timeout,
        positive(max_cache_size, "max_cache_size")?,
        positive(max_workers, "max_workers")?,
    )
}

/// Fetch one geometric ICRF vector using its explicit epoch scale.
///
/// Position is AU and velocity AU/day. `time_type` is `TDB` or `UT`; invalid
/// scale names are refused before transport. Cache, HTTP, parsing and sealed
/// policy errors propagate without a local model substitution.
pub fn fetch_horizons_state(
    client: &HorizonsClient,
    command: &str,
    jd: f64,
    center: &str,
    time_type: &str,
    policy: crate::context::EffectiveNetworkPolicy,
) -> Result<StateVector, Error> {
    let time_type = match time_type {
        "TDB" => HorizonsTimeType::Tdb,
        "UT" => HorizonsTimeType::Ut,
        _ => {
            return Err(Error::EphemerisLoaderValue {
                message: format!("unsupported Horizons time type: {time_type}"),
            });
        }
    };
    client.fetch_vectors_cached(command, jd, center, time_type, policy)
}

impl Session {
    /// Direct remote-provider position at UT1 JD using this Session's settings.
    ///
    /// The supplied client owns transport/cache state independently of the
    /// Session. Positions are degrees/degrees/AU and daily rates unless flags
    /// select radians or Cartesian output. Unsupported provider combinations
    /// remain an explicit refusal, rather than falling back to another model.
    pub fn horizons_backend_calc_ut(
        &self,
        client: &HorizonsClient,
        jd_ut: f64,
        body: i32,
        flags: i32,
        sid_mode: Option<i32>,
    ) -> Result<crate::positions::horizons_serving::HorizonsOutcome, Error> {
        let coverage = if crate::horizons::is_analytical_body(body)
            && flags
                & (crate::constants::FLG_HELCTR
                    | crate::constants::FLG_BARYCTR
                    | crate::constants::FLG_TOPOCTR
                    | crate::constants::FLG_ICRS)
                == 0
        {
            Some(self.mean_lunar_coverage(jd_ut + self.deltat(jd_ut)?, body)?)
        } else {
            None
        };
        let policy = match self.ctx.effective_network_policy() {
            crate::context::EffectiveNetworkPolicy::Sealed => crate::context::NetworkPolicy::Sealed,
            crate::context::EffectiveNetworkPolicy::Allow => crate::context::NetworkPolicy::Allow,
        };
        let mut context = self
            .ctx
            .clone()
            .with_mode(crate::context::ProviderMode::HorizonsRemote)
            .with_network_policy(Some(policy));
        if let Some(mode) = sid_mode {
            context.set_sidereal_mode(mode, 0.0, 0.0);
        }
        crate::positions::calc_horizons_backend_ut(
            crate::positions::ServingRequest {
                ctx: &context,
                delta_t_model: self.delta_t_model(),
                iers_enabled: context.effective_iers_delta_t(),
                observer: context.observer(),
                tier_policy: crate::positions::TierPolicy::AttachOrder,
                skip_closed: false,
            },
            client,
            jd_ut,
            body,
            flags,
            coverage,
        )
    }
}

/// Convert a civil calendar label between Julian and Gregorian calendars.
/// The input calendar changes at 1582-10-15; `calendar` selects the output
/// (`j`/`g`, ASCII case-insensitive). Hour is decimal hours. Invalid codes or
/// a non-finite hour are errors; all numerical conversion uses the time core.
pub fn convert_calendar_date(
    year: i32,
    month: i32,
    day: i32,
    hour: f64,
    calendar: &str,
) -> Result<crate::time::CalendarDate, Error> {
    let mut codes = calendar.chars();
    let code = codes.next().unwrap_or('\0');
    if codes.next().is_some() {
        return Err(Error::InvalidCalendarCode { code });
    }
    crate::time::date_conversion(
        year,
        month,
        day,
        hour,
        crate::time::CalendarTarget::from_code(code)?,
    )
}
