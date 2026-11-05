//! Apparition search orchestration over civil-day provider samples.
//!
//! A complete batch is acquired before any decision consumes it. Lookback,
//! streak and conjunction gates are deliberately separate from the provider's
//! twilight stencil and from the final visibility-window searches.

use super::{
    ELONG_GAP_DEG, ELONG_START_GAP_DEG, HELIACAL_MAX_DAYS, HELIACAL_SEARCH_BATCH_DAYS,
    HeliacalEphemeris, HeliacalEvent, HeliacalSearchRequest, MORNING_LAST_STREAK_DAYS,
    RISING_SETTING_STREAK_DAYS, refine_heliacal_time,
};
use crate::errors::Result;

/// Acquire every twilight sample in list order before the caller reads
/// elongations or makes streak decisions. Keep the provider's seasonal state
/// and callback order independent of early search success.
fn batch_twilight_samples<P: HeliacalEphemeris>(
    provider: &P,
    jd_days: &[f64],
    morning: bool,
) -> Vec<(bool, f64)> {
    jd_days
        .iter()
        .map(|&jd_day| provider.twilight_sample(jd_day, morning))
        .collect()
}

/// Look backward in descending day order, consuming the complete batch first.
/// The 999-degree sentinel must not satisfy a conjunction-gap test.
fn lookback_streak<P: HeliacalEphemeris>(
    provider: &P,
    jd_start: f64,
    morning: bool,
    days: usize,
) -> (usize, f64, bool) {
    let jd_days: Vec<f64> = (1..=days).map(|i| jd_start - i as f64).collect();
    let samples = batch_twilight_samples(provider, &jd_days, morning);
    let mut consecutive_invisible = 0_usize;
    let mut min_elong = 999.0_f64;
    let mut seen_visible = false;
    for (&jd_lb, &(vis, _)) in jd_days.iter().zip(samples.iter()) {
        if !vis {
            consecutive_invisible += 1;
            min_elong = min_elong.min(provider.elongation(jd_lb));
        } else {
            seen_visible = true;
            break;
        }
    }
    (consecutive_invisible, min_elong, seen_visible)
}

/// Rising-shaped scan shared only by morning-first and evening-first events.
/// The high-precision flag disables the pre-window gap rule, not refinement.
#[allow(clippy::too_many_arguments)]
fn search_rising_shape<P: HeliacalEphemeris>(
    provider: &P,
    jd_start: f64,
    morning: bool,
    is_star: bool,
    high_precision: bool,
    consecutive_invisible: usize,
    min_elong: f64,
    seen_visible: bool,
) -> Result<Option<f64>> {
    let mut consecutive_invisible = consecutive_invisible;
    let mut min_elong = min_elong;
    let mut seen_visible = seen_visible;
    // Every 100-day batch finishes before the first streak decision.
    let mut batch_start = 0_usize;
    while batch_start < HELIACAL_MAX_DAYS {
        let batch_end = (batch_start + HELIACAL_SEARCH_BATCH_DAYS).min(HELIACAL_MAX_DAYS);
        let jd_days: Vec<f64> = (batch_start..batch_end)
            .map(|d| jd_start + d as f64)
            .collect();
        let samples = batch_twilight_samples(provider, &jd_days, morning);
        for (&jd_day, &(vis, jd_vis)) in jd_days.iter().zip(samples.iter()) {
            if !vis {
                consecutive_invisible += 1;
                min_elong = min_elong.min(provider.elongation(jd_day));
            } else {
                let current_elong = provider.elongation(jd_day);
                let pre_window_gap = !seen_visible
                    && !high_precision
                    && min_elong <= ELONG_START_GAP_DEG
                    && current_elong >= min_elong + 2.0;
                if consecutive_invisible >= RISING_SETTING_STREAK_DAYS
                    && (is_star || min_elong <= ELONG_GAP_DEG || pre_window_gap)
                {
                    let jd_event = refine_heliacal_time(
                        |jd| provider.refinement_visible(jd),
                        jd_vis,
                        morning,
                    )?;
                    return Ok(Some(jd_event));
                }
                consecutive_invisible = 0;
                min_elong = 999.0;
                seen_visible = true;
            }
        }
        batch_start = batch_end;
    }
    Ok(None)
}

/// Search morning-first visibility from a UT Julian day, using the six-day
/// lookback and conditional days 7..30 extension before the forward scan.
/// Returns a refined UT epoch or `None`; refinement rejects nonfinite anchors.
pub fn search_heliacal_rising<P: HeliacalEphemeris>(
    jd_start: f64,
    provider: &P,
    is_star: bool,
    high_precision: bool,
) -> Result<Option<f64>> {
    let (consecutive_invisible, mut min_elong, mut seen_visible) =
        lookback_streak(provider, jd_start, true, 6);
    if !seen_visible && min_elong > ELONG_GAP_DEG {
        // The extended lookback also acquires all days before deciding.
        let jd_days: Vec<f64> = (7..=30).map(|i| jd_start - i as f64).collect();
        let samples = batch_twilight_samples(provider, &jd_days, true);
        for (&jd_lb, &(vis, _)) in jd_days.iter().zip(samples.iter()) {
            if vis {
                seen_visible = true;
                break;
            }
            min_elong = min_elong.min(provider.elongation(jd_lb));
        }
    }
    search_rising_shape(
        provider,
        jd_start,
        true,
        is_star,
        high_precision,
        consecutive_invisible,
        min_elong,
        seen_visible,
    )
}

/// Search evening-last visibility from a UT Julian day. Five invisible
/// evenings confirm the last visible day; a trailing day is still refined.
/// Returns a UT epoch or `None`; refinement rejects nonfinite anchors.
pub fn search_heliacal_setting<P: HeliacalEphemeris>(
    jd_start: f64,
    provider: &P,
) -> Result<Option<f64>> {
    let mut last_visible_jd = 0.0_f64;
    let mut consecutive_invisible = 0_usize;
    let mut batch_start = 0_usize;
    while batch_start < HELIACAL_MAX_DAYS {
        let batch_end = (batch_start + HELIACAL_SEARCH_BATCH_DAYS).min(HELIACAL_MAX_DAYS);
        let jd_days: Vec<f64> = (batch_start..batch_end)
            .map(|d| jd_start + d as f64)
            .collect();
        let samples = batch_twilight_samples(provider, &jd_days, false);
        for &(vis, jd_vis) in samples.iter() {
            if vis {
                last_visible_jd = jd_vis;
                consecutive_invisible = 0;
            } else {
                consecutive_invisible += 1;
                if consecutive_invisible >= RISING_SETTING_STREAK_DAYS && last_visible_jd > 0.0 {
                    let jd_event = refine_heliacal_time(
                        |jd| provider.refinement_visible(jd),
                        last_visible_jd,
                        false,
                    )?;
                    return Ok(Some(jd_event));
                }
            }
        }
        batch_start = batch_end;
    }
    // Setting refines its trailing visible day; morning-last does not.
    if last_visible_jd > 0.0 {
        let jd_event =
            refine_heliacal_time(|jd| provider.refinement_visible(jd), last_visible_jd, false)?;
        return Ok(Some(jd_event));
    }
    Ok(None)
}

/// Search first evening visibility from a UT Julian day. A clear seven-day
/// conjunction gate starts the scan at day 1; otherwise use the rising shape.
/// Returns a refined UT epoch or `None`; refinement rejects nonfinite anchors.
pub fn search_evening_first<P: HeliacalEphemeris>(
    jd_start: f64,
    provider: &P,
    is_star: bool,
    high_precision: bool,
) -> Result<Option<f64>> {
    if !is_star && provider.elongation(jd_start) <= ELONG_GAP_DEG {
        // All seven gate evenings precede the conjunction decision.
        let gate_days: Vec<f64> = (0..7).map(|i| jd_start - i as f64).collect();
        let gate_results = batch_twilight_samples(provider, &gate_days, false);
        let gate_clear = !gate_results.iter().any(|&(vis, _)| vis);
        if gate_clear {
            // A genuine conjunction starts the forward search at day 1.
            let mut batch_start = 1_usize;
            while batch_start < HELIACAL_MAX_DAYS {
                let batch_end = (batch_start + HELIACAL_SEARCH_BATCH_DAYS).min(HELIACAL_MAX_DAYS);
                let jd_days: Vec<f64> = (batch_start..batch_end)
                    .map(|d| jd_start + d as f64)
                    .collect();
                let samples = batch_twilight_samples(provider, &jd_days, false);
                for &(vis, jd_vis) in samples.iter() {
                    if vis {
                        let jd_event = refine_heliacal_time(
                            |jd| provider.refinement_visible(jd),
                            jd_vis,
                            false,
                        )?;
                        return Ok(Some(jd_event));
                    }
                }
                batch_start = batch_end;
            }
            return Ok(None);
        }
    }
    let (consecutive_invisible, min_elong, seen_visible) =
        lookback_streak(provider, jd_start, false, 6);
    search_rising_shape(
        provider,
        jd_start,
        false,
        is_star,
        high_precision,
        consecutive_invisible,
        min_elong,
        seen_visible,
    )
}

/// Search last morning visibility from a UT Julian day, requiring three
/// invisible mornings through the conjunction gap for confirmation. An
/// unconfirmed trailing day returns raw, without refinement. Returns `None`
/// when never visible; refinement rejects nonfinite anchors.
pub fn search_morning_last<P: HeliacalEphemeris>(
    jd_start: f64,
    provider: &P,
) -> Result<Option<f64>> {
    let mut last_visible_jd = 0.0_f64;
    let mut found_visible = false;
    let mut consecutive_invisible = 0_usize;
    let mut streak_min_elong = 999.0_f64;
    let mut batch_start = 0_usize;
    while batch_start < HELIACAL_MAX_DAYS {
        let batch_end = (batch_start + HELIACAL_SEARCH_BATCH_DAYS).min(HELIACAL_MAX_DAYS);
        let jd_days: Vec<f64> = (batch_start..batch_end)
            .map(|d| jd_start + d as f64)
            .collect();
        let samples = batch_twilight_samples(provider, &jd_days, true);
        for (&(vis, jd_vis), &jd_day) in samples.iter().zip(jd_days.iter()) {
            if vis {
                last_visible_jd = jd_vis;
                found_visible = true;
                consecutive_invisible = 0;
                streak_min_elong = 999.0;
            } else if found_visible {
                consecutive_invisible += 1;
                streak_min_elong = streak_min_elong.min(provider.elongation(jd_day));
                if consecutive_invisible >= MORNING_LAST_STREAK_DAYS
                    && streak_min_elong <= ELONG_GAP_DEG
                    && last_visible_jd > 0.0
                {
                    let jd_event = refine_heliacal_time(
                        |jd| provider.refinement_visible(jd),
                        last_visible_jd,
                        true,
                    )?;
                    return Ok(Some(jd_event));
                }
            }
        }
        batch_start = batch_end;
    }
    if last_visible_jd > 0.0 {
        return Ok(Some(last_visible_jd));
    }
    Ok(None)
}

/// Dispatch a validated event request over the provider. Returns its UT epoch
/// or `None`; an optimum already before the start restarts one day later.
/// Provider verdicts are infallible; refinement can reject a nonfinite anchor.
pub fn search_heliacal_event<P: HeliacalEphemeris>(
    request: &HeliacalSearchRequest,
    provider: &P,
) -> Result<Option<f64>> {
    let run = |jd_start: f64| -> Result<Option<f64>> {
        match request.event {
            HeliacalEvent::Rising => search_heliacal_rising(
                jd_start,
                provider,
                request.is_star,
                request.high_precision(),
            ),
            HeliacalEvent::Setting => search_heliacal_setting(jd_start, provider),
            HeliacalEvent::EveningFirst => search_evening_first(
                jd_start,
                provider,
                request.is_star,
                request.high_precision(),
            ),
            HeliacalEvent::MorningLast => search_morning_last(jd_start, provider),
        }
    };
    let jd_event = run(request.jd_start)?;
    if let Some(t) = jd_event {
        if t > 0.0 && t < request.jd_start - 1e-6 {
            return run(request.jd_start + 1.0);
        }
    }
    Ok(jd_event)
}
