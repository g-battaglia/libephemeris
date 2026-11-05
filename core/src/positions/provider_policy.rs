//! Provider selection and coverage-error policy.
//!
//! Planning uses resolved configuration and availability only. Coverage errors
//! preserve the requested body and distinguish it from failed support reads.
//! The served-place adapter delegates calculation to the existing serving path.

use super::{
    BackendLeg, JplFileDecision, LebMissOutcome, LegAvailability, PositionResult, RoutePlan,
    SealedServing, StateProvider, serving_compute_tt,
};
use crate::constants::{
    EARTH, FICT_OFFSET, INTP_APOG, INTP_PERG, MEAN_APOG, MEAN_NODE, PLUTO, SUN, WALDEMATH,
};
use crate::context::{EffectiveNetworkPolicy, ProviderMode, Tier};
use crate::errors::{Error, Result};

/// JPL DE kernel filename selected by the configured data tier.
pub fn tier_de_kernel(tier: Tier) -> &'static str {
    match tier {
        Tier::Base => "de440s.bsp",
        Tier::Medium => "de440.bsp",
        Tier::Extended => "de441.bsp",
    }
}

/// Select the effective JPL filename from already-acquired settings.
///
/// A nonblank environment value takes precedence, then the explicit setting,
/// then nonblank TOML, then the active tier's kernel. Selection is deterministic
/// and performs no filesystem or network access.
pub fn select_de_kernel<'a>(
    explicit: Option<&'a str>,
    env: Option<&'a str>,
    toml: Option<&'a str>,
    tier: Tier,
) -> &'a str {
    if let Some(value) = env {
        if !value.trim().is_empty() {
            return value;
        }
    }
    if let Some(value) = explicit {
        return value;
    }
    if let Some(value) = toml {
        if !value.trim().is_empty() {
            return value;
        }
    }
    tier_de_kernel(tier)
}

/// Whether a basename matches `de` plus digits, optional `s`, and `.bsp`,
/// case-insensitively. Recognized absent names can be downloaded on demand;
/// other absent names select the tier default. The basename uses `/` separators.
pub fn is_recognized_de_kernel(filename: &str) -> bool {
    let basename = filename.rsplit('/').next().unwrap_or(filename);
    let lower = basename.to_ascii_lowercase();
    let stem = lower.strip_suffix(".bsp");
    let Some(stem) = stem else {
        return false;
    };
    let Some(number) = stem.strip_prefix("de") else {
        return false;
    };
    let number = number.strip_suffix('s').unwrap_or(number);
    !number.is_empty() && number.bytes().all(|b| b.is_ascii_digit())
}

/// Use a recognized DE name or existing local file as requested.
/// Other absent names select the tier default, without downloading that name.
pub fn decide_jpl_file(requested_recognized: bool, requested_present: bool) -> JplFileDecision {
    if requested_recognized || requested_present {
        JplFileDecision::UseRequested
    } else {
        JplFileDecision::TierDefault
    }
}

/// Whether a core body's coverage miss must fail in sealed LEB mode.
/// Includes standard planets, Earth and geocentrically defined lunar points;
/// other families can continue to their declared analytical or propagated model.
pub fn is_leb_core_body(body_id: i32) -> bool {
    (SUN..=MEAN_APOG).contains(&body_id) || body_id == EARTH || is_leb_model_window_point(body_id)
}

/// Whether an interpolated lunar apsis uses the fitted model's time window.
/// The model owns its numerical limits; provider policy records membership only.
fn is_leb_model_window_point(body_id: i32) -> bool {
    body_id == INTP_APOG || body_id == INTP_PERG
}

/// Whether a core body requires its own stored LEB channel.
///
/// A missing channel without declared coverage is [`Error::UnknownBody`].
/// Analytical lunar points instead use their model's explicit validity window.
pub fn is_leb_stored_state_body(body_id: i32) -> bool {
    (SUN..=PLUTO).contains(&body_id) || body_id == EARTH
}

/// Whether a lunar point uses a runtime analytical/interpolated model rather
/// than its own stored channel. These points receive an explicit validity-window
/// check before evaluation.
pub fn is_leb_model_served_point(body_id: i32) -> bool {
    body_id == MEAN_NODE || body_id == MEAN_APOG || is_leb_model_window_point(body_id)
}

impl BackendLeg {
    /// Trace source tag for the backend that actually completed the calculation.
    pub fn as_source_tag(self) -> &'static str {
        match self {
            Self::Leb => "LEB",
            Self::NativeJpl => "NativeJpl",
            Self::Horizons => "Horizons",
        }
    }
}

impl LegAvailability {
    /// Snapshot both legs explicitly.
    pub fn new(leb_ready: bool, local_jpl_ready: bool) -> Self {
        Self {
            leb_ready,
            local_jpl_ready,
        }
    }

    /// Whether a sealed LEB reader serves core states.
    pub fn leb_ready(self) -> bool {
        self.leb_ready
    }

    /// Whether a local JPL source is provisioned.
    pub fn local_jpl_ready(self) -> bool {
        self.local_jpl_ready
    }
}

impl RoutePlan {
    /// Legs in attempt order; the first serving leg wins, per-leg miss
    /// policy decides fail-closed vs fall-through.
    pub fn legs(&self) -> &[BackendLeg] {
        &self.legs[..self.len]
    }

    /// Effective network gate for the planned attempts.
    pub fn network_gate(self) -> EffectiveNetworkPolicy {
        self.net
    }
}

/// Plan the backend chain from resolved mode, availability and network policy.
///
/// Sealed LEB offers only stored data; native JPL offers only local SPK.
/// Horizons offers its remote leg then native JPL on an eligible miss.
/// Auto offers available LEB first, remote vectors when no local JPL image
/// is available, then native JPL. The effective network policy is enforced
/// by acquisition; planning itself reads no files, sockets or mutable globals.
pub fn plan_route(
    mode: ProviderMode,
    avail: LegAvailability,
    net: EffectiveNetworkPolicy,
) -> RoutePlan {
    let mut legs = [BackendLeg::NativeJpl; 3];
    let mut len = 0;
    let mut push = |leg: BackendLeg| {
        legs[len] = leg;
        len += 1;
    };
    match mode {
        ProviderMode::SealedLeb => push(BackendLeg::Leb),
        ProviderMode::NativeJpl => push(BackendLeg::NativeJpl),
        ProviderMode::HorizonsRemote => {
            push(BackendLeg::Horizons);
            push(BackendLeg::NativeJpl);
        }
        ProviderMode::Auto => {
            if avail.leb_ready {
                push(BackendLeg::Leb);
            }
            if !avail.local_jpl_ready {
                push(BackendLeg::Horizons);
            }
            push(BackendLeg::NativeJpl);
        }
    }
    RoutePlan { legs, len, net }
}

/// Whether a LEB miss may continue to the next backend.
/// Explicit sealed LEB mode applies [`leb_miss_outcome`] instead.
pub fn leb_miss_fallthrough(mode: ProviderMode) -> bool {
    mode != ProviderMode::SealedLeb
}

/// Classify one body's LEB miss under the configured mode.
/// Sealed core-body misses return [`Error::EphemerisRange`], or
/// [`Error::UnknownBody`] when a required stored channel has no coverage.
/// Other modes can continue through their backend plan.
pub fn leb_miss_outcome(
    mode: ProviderMode,
    body_id: i32,
    jd: f64,
    coverage: Option<(f64, f64)>,
) -> LebMissOutcome {
    if let Some((start, end)) = coverage {
        if start <= jd && jd <= end {
            return LebMissOutcome::Served;
        }
    }
    if leb_miss_fallthrough(mode) || !is_leb_core_body(body_id) {
        return LebMissOutcome::FallThrough;
    }
    match coverage {
        None if !is_leb_stored_state_body(body_id) => LebMissOutcome::FallThrough,
        None => LebMissOutcome::Fail(Error::UnknownBody { body_id }),
        Some((start, end)) => LebMissOutcome::Fail(Error::EphemerisRange {
            requested_jd: jd,
            start_jd: start,
            end_jd: end,
            body: Some(body_id),
        }),
    }
}

/// Promote an eligible sealed provider miss to a public coverage error.
///
/// A core body with known requested coverage produces [`Error::EphemerisRange`];
/// a required stored channel without coverage produces [`Error::UnknownBody`].
/// For a non-core body, simultaneous absence of requested and Sun support can
/// name the requested body rather than the missing support channel.
/// Non-miss, non-sealed, in-range and non-finite errors propagate unchanged.
pub fn escalate_sealed_provider_miss<P: StateProvider>(
    mode: ProviderMode,
    provider: &P,
    body_id: i32,
    jd_tt: f64,
    err: Error,
) -> Error {
    if mode != ProviderMode::SealedLeb || !jd_tt.is_finite() {
        return err;
    }
    if (FICT_OFFSET..=WALDEMATH).contains(&body_id) {
        return err;
    }
    match err {
        Error::LebOutOfRange {
            body_id: err_body,
            jd,
            start,
            end,
        } => {
            if is_leb_core_body(body_id) {
                let (start, end) = provider.body_coverage(body_id).unwrap_or((start, end));
                Error::EphemerisRange {
                    requested_jd: jd,
                    start_jd: start,
                    end_jd: end,
                    body: Some(body_id),
                }
            } else if err_body != body_id {
                // A support-read miss inside the pipeline (e.g. the Earth
                // re-read): it names the failing support read, never the
                // requested minor body — propagate it unchanged.
                Error::LebOutOfRange {
                    body_id: err_body,
                    jd,
                    start,
                    end,
                }
            } else {
                match provider.body_coverage(SUN) {
                    Some((s, e)) if s <= jd_tt && jd_tt <= e => Error::LebOutOfRange {
                        body_id,
                        jd,
                        start,
                        end,
                    },
                    _ => Error::EphemerisRange {
                        requested_jd: jd_tt,
                        start_jd: start,
                        end_jd: end,
                        body: Some(body_id),
                    },
                }
            }
        }
        Error::LebBodyNotFound { .. } => {
            if is_leb_core_body(body_id) {
                match leb_miss_outcome(mode, body_id, jd_tt, provider.body_coverage(body_id)) {
                    LebMissOutcome::Fail(escalated) => escalated,
                    LebMissOutcome::Served | LebMissOutcome::FallThrough => {
                        Error::LebBodyNotFound { body_id }
                    }
                }
            } else if provider.body_coverage(body_id).is_some() {
                Error::LebBodyNotFound { body_id }
            } else {
                match provider.body_coverage(SUN) {
                    Some((s, e)) if s <= jd_tt && jd_tt <= e => Error::LebBodyNotFound { body_id },
                    Some((s, e)) => Error::EphemerisRange {
                        requested_jd: jd_tt,
                        start_jd: s,
                        end_jd: e,
                        body: Some(body_id),
                    },
                    None => Error::LebBodyNotFound { body_id },
                }
            }
        }
        miss @ Error::SealedChannelNotServed { .. } if !is_leb_core_body(body_id) => {
            // A retired or unsupported stored frame can refuse before its
            // raw state read. Coverage still decides whether the complete
            // local-model retry has the support it needs.
            match provider.body_coverage(body_id) {
                Some((start, end)) if jd_tt < start || jd_tt > end => {
                    match provider.body_coverage(SUN) {
                        Some((s, e)) if s <= jd_tt && jd_tt <= e => miss,
                        _ => Error::EphemerisRange {
                            requested_jd: jd_tt,
                            start_jd: start,
                            end_jd: end,
                            body: Some(body_id),
                        },
                    }
                }
                None => match provider.body_coverage(SUN) {
                    Some((start, end)) if jd_tt < start || jd_tt > end => Error::EphemerisRange {
                        requested_jd: jd_tt,
                        start_jd: start,
                        end_jd: end,
                        body: Some(body_id),
                    },
                    _ => miss,
                },
                _ => miss,
            }
        }
        other => other,
    }
}

/// Compute a stored apparent place and apply public sealed-coverage errors.
///
/// [`serving_compute_tt`] supplies the value. Eligible misses pass through
/// [`escalate_sealed_provider_miss`]; other errors are unchanged so the caller
/// can apply its backend plan.
pub fn serving_compute_tt_escalated<P: StateProvider>(
    mode: ProviderMode,
    serving: &SealedServing<'_, P>,
    jd_tt: f64,
    body_id: i32,
    calc_flags: i32,
) -> Result<(PositionResult, i32)> {
    match serving_compute_tt(serving, jd_tt, body_id, calc_flags) {
        Ok(ok) => Ok(ok),
        Err(err) => Err(escalate_sealed_provider_miss(
            mode,
            serving.provider,
            body_id,
            jd_tt,
            err,
        )),
    }
}

/// Require an available reader for sealed LEB mode; absence is a typed error.
pub fn require_leb(leb_ready: bool) -> Result<()> {
    if leb_ready {
        return Ok(());
    }
    Err(Error::ConfigurationMissing {
        missing_config: "leb_reader",
        suggestion: "bind a reviewed LEB core (set_leb_dir/attach_leb_file) or select another calculation mode",
    })
}

/// Reject native JPL/SPK acquisition in sealed LEB mode.
/// Only the active stored reader and declared analytical models are eligible.
pub fn forbid_jpl_in_sealed_mode(mode: ProviderMode) -> Result<()> {
    if mode != ProviderMode::SealedLeb {
        return Ok(());
    }
    Err(Error::ConfigurationMissing {
        missing_config: "jpl_source_access",
        suggestion: "sealed LEB mode never opens JPL kernels; use the bound LEB reader or switch calculation mode",
    })
}
