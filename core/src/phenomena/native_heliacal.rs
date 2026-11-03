//! Fallible heliacal event and visibility-window walks over retained acquisitions.
use crate::errors::{Error, Result};
use crate::heliacal::{HELIACAL_MAX_DAYS, HELIACAL_SEARCH_BATCH_DAYS, HeliacalEvent};

// Test-only observations are made after the ordinary acquisition, never by
// asking the provider to reproduce a predicate or an elongation.
#[cfg(test)]
mod observation {
    use std::cell::RefCell;

    thread_local! {
        static EVENTS: RefCell<Option<Vec<serde_json::Value>>> = const { RefCell::new(None) };
    }

    pub(super) fn record(value: impl FnOnce() -> serde_json::Value) {
        let _ = EVENTS.try_with(|events| {
            if let Some(events) = events.borrow_mut().as_mut() {
                events.push(value());
            }
        });
    }

    pub(super) fn begin() {
        EVENTS.with(|events| assert!(events.replace(Some(Vec::new())).is_none()));
    }

    pub(super) fn finish() -> Vec<serde_json::Value> {
        EVENTS.with(|events| events.take().expect("active heliacal observation"))
    }

    pub(super) fn word(value: f64) -> serde_json::Value {
        serde_json::json!([value, value.to_bits()])
    }
}

/// Shared epoch resolution for search refinement and public window boundaries.
const TIME_TOLERANCE_DAYS: f64 = 1e-6;

/// One owned search's acquisition interface. A batch is evaluated completely
/// before its streak decisions; a failed acquisition remains an error, never
/// an invisible-day or exhausted-search substitute. All epochs are UT days.
pub trait HeliacalAcquisition {
    /// Ordered civil-day twilight detections: first visibility epoch or zero.
    fn twilight_batch(&mut self, days: &[f64], morning: bool) -> Result<Vec<(bool, f64)>>;
    /// Geocentric solar elongation in degrees at the unrounded requested day.
    fn elongation(&mut self, jd: f64) -> Result<f64>;
    /// Full Moon-inclusive visibility predicate at a refinement epoch.
    fn visible(&mut self, jd: f64) -> Result<bool>;
}

/// Resolve a narrow positive twilight interval between coarse grid points.
/// The seed is the strongest eligible sample; its two neighbouring grid
/// intervals bound the local maximum. Callback epochs/spacing are UT days and
/// values are V-magnitude margins. Excluded geometry is `None`; actual errors
/// propagate. Resolution is the search's existing 1e-6-day refinement scale.
pub(crate) fn refine_twilight_peak(
    seed: f64,
    spacing: f64,
    mut margin: impl FnMut(f64) -> Result<Option<f64>>,
) -> Result<Option<f64>> {
    let mut best = None;
    let mut objective = |jd| {
        let Some(value) = margin(jd)? else {
            return Ok(f64::INFINITY);
        };
        if !value.is_finite() {
            return Err(Error::InvalidEclipseSearch {
                entry: "heliacal_ut",
                detail: "twilight refinement margin is nonfinite".into(),
            });
        }
        if best.is_none_or(|(_, previous)| value > previous) {
            best = Some((jd, value));
        }
        Ok(-value)
    };
    objective(seed)?;
    let peak = crate::events::golden_section_min(
        &mut objective,
        seed - spacing,
        seed + spacing,
        TIME_TOLERANCE_DAYS,
        30,
    )?;
    objective(peak)?;
    Ok(best.filter(|(_, value)| *value > 0.0).map(|(jd, _)| jd))
}

/// Refine a visibility transition in a +/-0.1-day bracket. Morning visibility
/// moves the upper endpoint; evening visibility moves the lower endpoint.
/// Stop after 30 samples or a bracket narrower than 1e-6 day. No finite gate
/// precedes the first provider read; provider failures propagate unchanged.
fn refine(feed: &mut impl HeliacalAcquisition, epoch: f64, morning: bool) -> Result<f64> {
    let (mut lo, mut hi) = (epoch - 0.1, epoch + 0.1);
    for _ in 0..30 {
        let mid = (lo + hi) / 2.0;
        let visible = feed.visible(mid)?;
        #[cfg(test)]
        observation::record(
            || serde_json::json!({"phase":"refine", "mid":observation::word(mid), "visible":visible, "morning":morning}),
        );
        if visible == morning {
            hi = mid;
        } else {
            lo = mid;
        }
        if hi - lo < TIME_TOLERANCE_DAYS {
            break;
        }
    }
    Ok((lo + hi) / 2.0)
}

/// Require one verdict for every supplied civil day without silently trimming
/// or duplicating a batch. The provider owns all samples before this returns.
fn batch(
    feed: &mut impl HeliacalAcquisition,
    days: &[f64],
    morning: bool,
) -> Result<Vec<(bool, f64)>> {
    let values = feed.twilight_batch(days, morning)?;
    if values.len() != days.len() {
        return Err(Error::InvalidEclipseSearch {
            entry: "heliacal_ut",
            detail: "twilight batch returned a different number of days".into(),
        });
    }
    #[cfg(test)]
    observation::record(|| {
        serde_json::json!({"phase":"batch", "morning":morning,
        "days":days.iter().map(|&jd| observation::word(jd)).collect::<Vec<_>>(),
        "verdicts":values.iter().map(|&(visible,hit)| (visible,observation::word(hit))).collect::<Vec<_>>()})
    });
    Ok(values)
}

/// Six-day lookback in descending epoch order. All twilight samples precede
/// the elongation reads; the streak ends at the first visible lookback day.
fn lookback(
    feed: &mut impl HeliacalAcquisition,
    start: f64,
    morning: bool,
) -> Result<(usize, f64, bool)> {
    let days: Vec<_> = (1..=6).map(|i| start - i as f64).collect();
    let values = batch(feed, &days, morning)?;
    let (mut count, mut minimum, mut seen) = (0, 999.0, false);
    for (jd, (visible, _)) in days.into_iter().zip(values) {
        if visible {
            seen = true;
            #[cfg(test)]
            observation::record(
                || serde_json::json!({"phase":"lookback-visible", "jd":observation::word(jd), "count":count, "minimum":observation::word(minimum), "seen":seen}),
            );
            break;
        }
        count += 1;
        let elongation = feed.elongation(jd)?;
        if elongation < minimum {
            minimum = elongation;
        }
        #[cfg(test)]
        observation::record(
            || serde_json::json!({"phase":"lookback-invisible", "jd":observation::word(jd), "elongation":observation::word(elongation), "count":count, "minimum":observation::word(minimum), "seen":seen}),
        );
    }
    Ok((count, minimum, seen))
}

/// First visibility following a conjunction gap. The six-day lookback may
/// extend to days7..30 for a morning search with no previous visibility.
/// A pre-window gap admits the next apparition only when the minimum is at
/// most30 degrees, elongation has grown by2 degrees, and HIGH_PRECISION is off.
fn first(
    feed: &mut impl HeliacalAcquisition,
    start: f64,
    morning: bool,
    star: bool,
    high: bool,
) -> Result<f64> {
    let (mut count, mut minimum, mut seen) = lookback(feed, start, morning)?;
    #[cfg(test)]
    observation::record(
        || serde_json::json!({"phase":"first-initial", "start":observation::word(start), "morning":morning, "star":star, "high":high, "count":count, "minimum":observation::word(minimum), "seen":seen}),
    );
    if morning && !seen && minimum > 10.0 {
        let days: Vec<_> = (7..31).map(|i| start - i as f64).collect();
        let values = batch(feed, &days, true)?;
        for (jd, (visible, _)) in days.into_iter().zip(values) {
            if visible {
                seen = true;
                #[cfg(test)]
                observation::record(
                    || serde_json::json!({"phase":"extended-lookback-visible", "jd":observation::word(jd), "count":count, "minimum":observation::word(minimum), "seen":seen}),
                );
                break;
            }
            let elongation = feed.elongation(jd)?;
            if elongation < minimum {
                minimum = elongation;
            }
            #[cfg(test)]
            observation::record(
                || serde_json::json!({"phase":"extended-lookback-invisible", "jd":observation::word(jd), "elongation":observation::word(elongation), "count":count, "minimum":observation::word(minimum), "seen":seen}),
            );
        }
    }
    for begin in (0..HELIACAL_MAX_DAYS).step_by(HELIACAL_SEARCH_BATCH_DAYS) {
        let days: Vec<_> = (begin..(begin + HELIACAL_SEARCH_BATCH_DAYS).min(HELIACAL_MAX_DAYS))
            .map(|i| start + i as f64)
            .collect();
        let values = batch(feed, &days, morning)?;
        for (jd, (visible, hit)) in days.into_iter().zip(values) {
            if !visible {
                count += 1;
                let elongation = feed.elongation(jd)?;
                if elongation < minimum {
                    minimum = elongation;
                }
            } else {
                let elongation = feed.elongation(jd)?;
                let pre_window = !seen && !high && minimum <= 30.0 && elongation >= minimum + 2.0;
                #[cfg(test)]
                observation::record(
                    || serde_json::json!({"phase":"first-candidate", "jd":observation::word(jd), "hit":observation::word(hit), "elongation":observation::word(elongation), "count":count, "minimum":observation::word(minimum), "seen":seen, "pre_window":pre_window, "admitted":count>=5 && (star || minimum<=10.0 || pre_window)}),
                );
                if count >= 5 && (star || minimum <= 10.0 || pre_window) {
                    return refine(feed, hit, morning);
                }
                count = 0;
                minimum = 999.0;
                seen = true;
            }
        }
    }
    Ok(0.0)
}

/// Forward event search with bounded800-day,100-day batch chronology.
/// Returns zero for the actual completed no-event walk. Evening first tests
/// the current plus six previous evenings before its conjunction fast path.
/// Morning last preserves the unrefined trailing hit when no confirming gap
/// occurs; evening last refines that trailing hit. A classified event before
/// the caller start retries exactly once from the following day.
pub fn search(
    feed: &mut impl HeliacalAcquisition,
    start: f64,
    event: HeliacalEvent,
    star: bool,
    flags: i32,
) -> Result<f64> {
    #[cfg(test)]
    observation::record(
        || serde_json::json!({"phase":"search", "start":observation::word(start), "event":event as i32, "star":star, "flags":flags}),
    );
    let hit = walk(feed, start, event, star, flags)?;
    #[cfg(test)]
    observation::record(
        || serde_json::json!({"phase":"walk-result", "hit":observation::word(hit), "retry":0.0<hit && hit<start-1e-6}),
    );
    if 0.0 < hit && hit < start - 1e-6 {
        walk(feed, start + 1.0, event, star, flags)
    } else {
        Ok(hit)
    }
}

/// One bounded search without the outer previous-day retry.
fn walk(
    feed: &mut impl HeliacalAcquisition,
    start: f64,
    event: HeliacalEvent,
    star: bool,
    flags: i32,
) -> Result<f64> {
    let high = flags & crate::constants::HELFLAG_HIGH_PRECISION != 0;
    match event {
        HeliacalEvent::Rising => first(feed, start, true, star, high),
        HeliacalEvent::EveningFirst => {
            if !star
                && {
                    let elongation = feed.elongation(start)?;
                    #[cfg(test)]
                    observation::record(
                        || serde_json::json!({"phase":"evening-conjunction", "jd":observation::word(start), "elongation":observation::word(elongation), "within_ten_degrees":elongation<=10.0}),
                    );
                    elongation <= 10.0
                }
            {
                let days: Vec<_> = (0..7).map(|i| start - i as f64).collect();
                if !batch(feed, &days, false)?.iter().any(|(v, _)| *v) {
                    for begin in (1..HELIACAL_MAX_DAYS).step_by(HELIACAL_SEARCH_BATCH_DAYS) {
                        let days: Vec<_> = (begin
                            ..(begin + HELIACAL_SEARCH_BATCH_DAYS).min(HELIACAL_MAX_DAYS))
                            .map(|i| start + i as f64)
                            .collect();
                        for (visible, hit) in batch(feed, &days, false)? {
                            if visible {
                                #[cfg(test)]
                                observation::record(
                                    || serde_json::json!({"phase":"evening-fast-candidate", "hit":observation::word(hit)}),
                                );
                                return refine(feed, hit, false);
                            }
                        }
                    }
                    return Ok(0.0);
                }
            }
            first(feed, start, false, star, high)
        }
        HeliacalEvent::Setting | HeliacalEvent::MorningLast => {
            let morning = event == HeliacalEvent::MorningLast;
            let (mut last, mut count, mut seen, mut minimum) = (0.0, 0, false, 999.0);
            for begin in (0..HELIACAL_MAX_DAYS).step_by(HELIACAL_SEARCH_BATCH_DAYS) {
                let days: Vec<_> = (begin
                    ..(begin + HELIACAL_SEARCH_BATCH_DAYS).min(HELIACAL_MAX_DAYS))
                    .map(|i| start + i as f64)
                    .collect();
                let values = batch(feed, &days, morning)?;
                for (jd, (visible, hit)) in days.into_iter().zip(values) {
                    if visible {
                        last = hit;
                        seen = true;
                        count = 0;
                        minimum = 999.0;
                        #[cfg(test)]
                        observation::record(
                            || serde_json::json!({"phase":"last-visible", "jd":observation::word(jd), "hit":observation::word(hit), "count":count, "minimum":observation::word(minimum), "seen":seen, "morning":morning}),
                        );
                    } else if !morning || seen {
                        count += 1;
                        if morning {
                            let elongation = feed.elongation(jd)?;
                            if elongation < minimum {
                                minimum = elongation;
                            }
                        }
                        if last > 0.0
                            && ((!morning && count >= 5)
                                || (morning && count >= 3 && minimum <= 10.0))
                        {
                            #[cfg(test)]
                            observation::record(
                                || serde_json::json!({"phase":"last-admitted", "jd":observation::word(jd), "hit":observation::word(last), "count":count, "minimum":observation::word(minimum), "seen":seen, "morning":morning}),
                            );
                            return refine(feed, last, morning);
                        }
                    }
                }
            }
            if last > 0.0 && !morning {
                refine(feed, last, false)
            } else {
                Ok(last)
            }
        }
    }
}

/// Visibility-window start/optimum/end in UT days. Scan145 points across
/// +/-0.12 day, refine the strongest grid neighbourhood for20 golden-section
/// steps, then bracket each edge with25 bisections to1e-6 day. Every margin
/// failure propagates at its actual sample; no sentinel hides a failed read.
/// A nonfinite margin or a candidate with no positive visibility interval is
/// an [`Error::InvalidEclipseSearch`], rather than a fabricated event epoch.
/// Epochs are UT Julian days; a nonfinite epoch returns
/// [`Error::NonFiniteArgument`] after the first fallible margin acquisition.
pub fn visibility_window(
    epoch: f64,
    mut evaluate: impl FnMut(f64) -> Result<f64>,
) -> Result<[f64; 3]> {
    let mut margin = |jd| {
        let value = evaluate(jd)?;
        if !value.is_finite() {
            return Err(Error::InvalidEclipseSearch {
                entry: "heliacal_ut",
                detail: "visibility margin is nonfinite".into(),
            });
        }
        Ok(value)
    };
    let (mut best_jd, mut best) = (epoch, margin(epoch)?);
    if !epoch.is_finite() {
        return Err(Error::NonFiniteArgument { value: epoch });
    }
    for i in 0..=144 {
        let t = epoch - 0.12 + (2.0 * 0.12) * i as f64 / 144.0;
        let m = margin(t)?;
        if m > best {
            best = m;
            best_jd = t;
        }
    }
    let step = (2.0 * 0.12) / 144.0;
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let (mut lo, mut hi) = (best_jd - step, best_jd + step);
    let (mut a, mut b) = (hi - (hi - lo) / phi, lo + (hi - lo) / phi);
    let (mut fa, mut fb) = (margin(a)?, margin(b)?);
    for _ in 0..20 {
        if fa > fb {
            hi = b;
            b = a;
            fb = fa;
            a = hi - (hi - lo) / phi;
            fa = margin(a)?;
        } else {
            lo = a;
            a = b;
            fa = fb;
            b = lo + (hi - lo) / phi;
            fb = margin(b)?;
        }
        if hi - lo < TIME_TOLERANCE_DAYS {
            break;
        }
    }
    let refined_jd = 0.5 * (lo + hi);
    let refined_margin = margin(refined_jd)?;
    if refined_margin > best {
        best_jd = refined_jd;
        best = refined_margin;
    }
    if best <= 0.0 {
        return Err(Error::InvalidEclipseSearch {
            entry: "heliacal_ut",
            detail: "candidate has no positive visibility interval".into(),
        });
    }
    let start = crossing(best_jd, best_jd - 0.12, &mut margin)?;
    let end = crossing(best_jd, best_jd + 0.12, &mut margin)?;
    Ok([
        if best_jd < start { best_jd } else { start },
        best_jd,
        if best_jd > end { best_jd } else { end },
    ])
}

/// Fallible zero crossing from the positive interior towards one bracket edge.
fn crossing(
    mut inner: f64,
    mut outer: f64,
    margin: &mut impl FnMut(f64) -> Result<f64>,
) -> Result<f64> {
    let (inside, outside) = (margin(inner)?, margin(outer)?);
    if inside <= 0.0 {
        return Ok(inner);
    }
    if outside > 0.0 {
        return Ok(outer);
    }
    for _ in 0..25 {
        let mid = 0.5 * (inner + outer);
        if margin(mid)? > 0.0 {
            inner = mid;
        } else {
            outer = mid;
        }
        if (outer - inner).abs() < TIME_TOLERANCE_DAYS {
            break;
        }
    }
    Ok(0.5 * (inner + outer))
}

#[cfg(test)]
mod tests {
    use super::*;
    struct Failed;
    impl HeliacalAcquisition for Failed {
        fn twilight_batch(&mut self, _: &[f64], _: bool) -> Result<Vec<(bool, f64)>> {
            Err(Error::UnknownBody { body_id: 991 })
        }
        fn elongation(&mut self, _: f64) -> Result<f64> {
            panic!("after failed batch")
        }
        fn visible(&mut self, _: f64) -> Result<bool> {
            panic!("after failed batch")
        }
    }
    #[test]
    fn failed_acquisition_is_not_invisible_or_no_event() {
        assert!(matches!(
            search(&mut Failed, f64::NAN, HeliacalEvent::Rising, false, 0),
            Err(Error::UnknownBody { body_id: 991 })
        ));
    }
    #[test]
    fn fallible_window_keeps_first_failed_sample() {
        let mut calls = 0;
        let value = visibility_window(f64::INFINITY, |_| {
            calls += 1;
            Err(Error::UnknownBody { body_id: 992 })
        });
        assert_eq!(calls, 1);
        assert!(matches!(value, Err(Error::UnknownBody { body_id: 992 })));
    }
    #[test]
    fn nonpositive_visibility_is_not_a_finite_event() {
        for value in [-99.0, 0.0] {
            assert!(matches!(
                visibility_window(2451545.0, |_| Ok(value)),
                Err(Error::InvalidEclipseSearch { .. })
            ));
        }
    }
    #[test]
    fn undefined_margin_refuses_at_the_consumed_sample() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut calls = 0;
            assert!(matches!(
                visibility_window(2451545.0, |_| {
                    calls += 1;
                    Ok(value)
                }),
                Err(Error::InvalidEclipseSearch { .. })
            ));
            assert_eq!(calls, 1);
        }
    }
    #[test]
    fn finite_margin_cannot_admit_an_undefined_epoch() {
        for epoch in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut calls = 0;
            assert!(matches!(
                visibility_window(epoch, |_| {
                    calls += 1;
                    Ok(1.0)
                }),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert_eq!(calls, 1);
        }
    }
    #[test]
    fn narrow_positive_interval_is_refined_before_admission() {
        let peak = 2451545.0 + 0.0006;
        let curve = |jd: f64| Ok(1.0 - ((jd - peak) / 0.0001).powi(2));
        let times = visibility_window(2451545.0, curve).unwrap();
        assert!(times[0] < times[1] && times[1] < times[2]);
        assert!(curve(times[1]).unwrap() > 0.0);
        assert!((times[1] - peak).abs() * 86400.0 < 1.0);
    }
    #[test]
    fn twilight_peak_detects_a_visible_interval_between_negative_grid_points() {
        let seed = 2451545.0;
        let spacing = 15.0 / 1440.0;
        let peak = seed + 3.0 / 1440.0;
        let curve = |jd: f64| 1.0 - ((jd - peak) / (30.0 / 86400.0)).powi(2);
        assert!(
            [-spacing, 0.0, spacing]
                .iter()
                .all(|dt| curve(seed + dt) < 0.0)
        );
        let found = refine_twilight_peak(seed, spacing, |jd| Ok(Some(curve(jd))))
            .unwrap()
            .unwrap();
        assert!(curve(found) > 0.0);
        assert!((found - peak).abs() * 86400.0 < 1.0);
    }
    #[test]
    fn twilight_peak_keeps_exclusion_failure_and_invisibility_distinct() {
        assert!(
            refine_twilight_peak(2451545.0, 0.01, |_| Ok(None))
                .unwrap()
                .is_none()
        );
        assert!(
            refine_twilight_peak(2451545.0, 0.01, |_| Ok(Some(-0.01)))
                .unwrap()
                .is_none()
        );
        assert!(matches!(
            refine_twilight_peak(2451545.0, 0.01, |_| Err(Error::UnknownBody {
                body_id: 991
            })),
            Err(Error::UnknownBody { .. })
        ));
        assert!(matches!(
            refine_twilight_peak(2451545.0, 0.01, |_| Ok(Some(f64::NAN))),
            Err(Error::InvalidEclipseSearch { .. })
        ));
    }
    struct Marked {
        batches: Vec<Vec<f64>>,
        reads: Vec<f64>,
    }
    impl HeliacalAcquisition for Marked {
        fn twilight_batch(&mut self, days: &[f64], _: bool) -> Result<Vec<(bool, f64)>> {
            self.batches.push(days.to_vec());
            Ok(days.iter().map(|_| (false, 0.0)).collect())
        }
        fn elongation(&mut self, jd: f64) -> Result<f64> {
            self.reads.push(jd);
            Ok(5.0)
        }
        fn visible(&mut self, _: f64) -> Result<bool> {
            panic!("invisible chronology must not refine")
        }
    }
    #[test]
    fn evening_conjunction_checks_seven_days_then_starts_day_one() {
        let mut feed = Marked {
            batches: vec![],
            reads: vec![],
        };
        assert_eq!(
            search(&mut feed, 100.0, HeliacalEvent::EveningFirst, false, 0).unwrap(),
            0.0
        );
        assert_eq!(feed.batches[0], [100.0, 99.0, 98.0, 97.0, 96.0, 95.0, 94.0]);
        assert_eq!(feed.batches[1][0], 101.0);
        assert_eq!(*feed.batches.last().unwrap().last().unwrap(), 899.0);
        assert_eq!(feed.reads, [100.0]);
    }

    #[test]
    fn observation_preserves_reads_results_and_failed_short_circuit() {
        let mut ordinary = Marked {
            batches: vec![],
            reads: vec![],
        };
        let expected = search(&mut ordinary, 100.0, HeliacalEvent::EveningFirst, false, 0).unwrap();
        let mut observed = Marked {
            batches: vec![],
            reads: vec![],
        };
        observation::begin();
        let actual = search(&mut observed, 100.0, HeliacalEvent::EveningFirst, false, 0).unwrap();
        let events = observation::finish();
        assert_eq!(actual.to_bits(), expected.to_bits());
        assert_eq!(observed.batches, ordinary.batches);
        assert_eq!(observed.reads, ordinary.reads);
        assert!(
            events
                .iter()
                .any(|row| row["phase"] == "evening-conjunction")
        );
        observation::begin();
        let failed = search(&mut Failed, 100.0, HeliacalEvent::Rising, false, 0);
        let events = observation::finish();
        assert!(matches!(failed, Err(Error::UnknownBody { body_id: 991 })));
        assert_eq!(events.len(), 1);
        assert_eq!(events[0]["phase"], "search");
    }

    /// Replay only the twelve original differing public requests. The collector
    /// observes reads already made by the real finder; it adds no visibility
    /// evaluation. Historical result differences and actual failures are saved
    /// before the terminal equality assertion, never replaced with a sentinel.
    #[test]
    #[ignore = "explicit original-request observation; needs existing read-only astronomical files"]
    fn observe_original_heliacal_apparitions() {
        use serde_json::{Value, json};
        use sha2::{Digest, Sha256};
        use std::collections::{BTreeMap, btree_map::Entry};
        use std::io::{BufRead, BufReader, Read, Write};
        use std::path::{Path, PathBuf};

        fn fingerprint(path: &Path) -> Value {
            let mut file = std::fs::File::open(path).unwrap();
            let metadata = file.metadata().unwrap();
            let mut hash = Sha256::new();
            let mut buffer = [0_u8; 65536];
            loop {
                let count = file.read(&mut buffer).unwrap();
                if count == 0 {
                    break;
                }
                hash.update(&buffer[..count]);
            }
            assert_eq!(metadata.len(), file.metadata().unwrap().len());
            json!({"path":path,"bytes":metadata.len(),"sha256":hash.finalize().iter().map(|byte| format!("{byte:02x}")).collect::<String>(),
                "modified_ns":metadata.modified().unwrap().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos().to_string()})
        }
        fn words(value: &Value) -> Vec<f64> {
            value
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_f64().unwrap())
                .collect()
        }
        fn runtime(row: &Value) -> crate::config::Runtime {
            let config = &row["config"];
            let data = config["data_dir"]
                .as_str()
                .map(PathBuf::from)
                .unwrap_or_else(|| {
                    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                        .parent()
                        .expect("source workspace")
                        .join("target/v4/comparison-data")
                });
            let mut context = crate::context::EphemerisContext::new().with_jpl_data_dir(Some(data));
            context.set_iers_auto_download(Some(false));
            let runtime = crate::config::Runtime::from_context(context);
            runtime.set_network_policy(Some("sealed")).unwrap();
            runtime.set_auto_spk_download(Some(false));
            runtime
                .set_precision_tier(config["tier"].as_str().unwrap())
                .unwrap();
            runtime.set_calc_mode(config["mode"].as_str()).unwrap();
            runtime.set_leb_file(config["leb_file"].as_str().map(PathBuf::from));
            runtime.set_iers_delta_t_enabled(config["iers_delta_t_enabled"].as_bool());
            runtime.set_delta_t_userdef(config["delta_t_userdef"].as_f64());
            if let Some(path) = config["ephe_path"].as_str() {
                runtime.set_ephe_path(Some(PathBuf::from(path)));
            }
            if let Some(name) = config["ephemeris_file"].as_str() {
                runtime.set_ephemeris_file(name);
            }
            if let Some(sidereal) = config["sidereal"].as_array() {
                runtime.set_sid_mode(
                    sidereal[0].as_i64().unwrap() as i32,
                    sidereal[1].as_f64().unwrap(),
                    sidereal[2].as_f64().unwrap(),
                );
            }
            if let Some(observer) = config["observer"].as_array() {
                runtime
                    .set_topo(
                        observer[0].as_f64().unwrap(),
                        observer[1].as_f64().unwrap(),
                        observer[2].as_f64().unwrap(),
                    )
                    .unwrap();
            }
            for setup in row["setup"].as_array().into_iter().flatten() {
                let args = setup["args"].as_array();
                match setup["function"].as_str().unwrap() {
                    "set_topo" => {
                        let args = args.unwrap();
                        runtime
                            .set_topo(
                                args[0].as_f64().unwrap(),
                                args[1].as_f64().unwrap(),
                                args[2].as_f64().unwrap(),
                            )
                            .unwrap();
                    }
                    "set_sid_mode" => {
                        let args = args.unwrap();
                        runtime.set_sid_mode(args[0].as_i64().unwrap() as i32, 0.0, 0.0);
                    }
                    "set_iers_cache_dir" => runtime.set_iers_cache_dir(Some(PathBuf::from(
                        args.unwrap()[0].as_str().unwrap(),
                    ))),
                    "set_iers_auto_download" => {
                        runtime.set_iers_auto_download(Some(args.unwrap()[0].as_bool().unwrap()));
                    }
                    "clear_iers_cache" => runtime.clear_iers_cache(),
                    "load_iers_data" => assert!(runtime.load_iers_data(false).unwrap()),
                    operation => panic!("unsupported explicit heliacal capture setup: {operation}"),
                }
            }
            runtime
        }

        let root = PathBuf::from(
            std::env::var_os("HELIACAL_ORIGINAL_CAPTURE_DIR")
                .expect("HELIACAL_ORIGINAL_CAPTURE_DIR"),
        );
        let output = PathBuf::from(
            std::env::var_os("HELIACAL_OBSERVATION_OUTPUT").expect("HELIACAL_OBSERVATION_OUTPUT"),
        );
        let mut stream = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)
            .unwrap();
        let captures = [
            (
                "native-public.jsonl",
                "61e0dc716625b9de4b342a0424b3d3c95a555ff3c31e95b1a91bccbaef34fee4",
                [
                    "search/native/Venus/1",
                    "search/native/Mercury/3",
                    "search/native/Sirius/2",
                    "search/leb/Venus/1",
                    "search/leb/Mercury/3",
                    "search/leb/Sirius/2",
                ],
            ),
            (
                "native-events.jsonl",
                "8c6e21686df63e586927c71a541bf596153c9e56add70ea3dfbab21a472fd192",
                [
                    "rising/mars-30N/flags-default",
                    "rising/mars-30N/flags-0",
                    "setting/mars-30N/flags-default",
                    "setting/mars-30N/flags-0",
                    "setting/mercury-30N/flags-default",
                    "setting/mercury-30N/flags-0",
                ],
            ),
        ];
        let mut rows = Vec::new();
        let mut before = BTreeMap::new();
        let current_input = std::env::var_os("HELIACAL_CURRENT_INPUT").map(PathBuf::from);
        if let Some(path) = &current_input {
            let input_identity = fingerprint(path);
            let expected = std::env::var("HELIACAL_CURRENT_INPUT_SHA256")
                .expect("HELIACAL_CURRENT_INPUT_SHA256");
            assert_eq!(input_identity["sha256"], expected);
            before.insert(path.clone(), input_identity);
            let mut identifiers = std::collections::BTreeSet::new();
            for line in BufReader::new(std::fs::File::open(path).unwrap()).lines() {
                let row: Value = serde_json::from_str(&line.unwrap()).unwrap();
                assert!(identifiers.insert(row["id"].as_str().unwrap().to_owned()));
                rows.push(row);
            }
            assert_eq!(
                rows.len(),
                33,
                "complete current provider/legacy event scope"
            );
        } else {
            for (filename, expected_sha, ids) in captures {
                let path = root.join(filename);
                let identity = fingerprint(&path);
                assert_eq!(identity["sha256"], expected_sha);
                before.insert(path.clone(), identity);
                let mut selected = BTreeMap::new();
                for line in BufReader::new(std::fs::File::open(path).unwrap()).lines() {
                    let row: Value = serde_json::from_str(&line.unwrap()).unwrap();
                    if ids.contains(&row["id"].as_str().unwrap()) {
                        assert!(
                            selected
                                .insert(row["id"].as_str().unwrap().to_owned(), row)
                                .is_none()
                        );
                    }
                }
                for id in ids {
                    rows.push(selected.remove(id).expect("original request"));
                }
            }
        }
        let requests = rows.len();
        for row in &rows {
            assert_eq!(row["engine"], "v4");
            assert_eq!(row["function"], "heliacal_ut");
            assert_eq!(row["version"], "4.0.0");
            for input in row["inputs"].as_object().unwrap().values() {
                let path = PathBuf::from(input["path"].as_str().unwrap());
                match before.entry(path) {
                    Entry::Vacant(entry) => {
                        let actual = fingerprint(entry.key());
                        assert_eq!(actual["sha256"], input["sha256"]);
                        assert_eq!(actual["bytes"], input["bytes"]);
                        entry.insert(actual);
                    }
                    Entry::Occupied(entry) => {
                        assert_eq!(entry.get()["sha256"], input["sha256"]);
                    }
                }
            }
        }
        let source_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        for name in [
            "src/phenomena/native_heliacal.rs",
            "src/phenomena/heliacal.rs",
            "src/phenomena/heliacal/atmosphere.rs",
            "src/phenomena/heliacal/composition.rs",
            "src/phenomena/heliacal/event_window.rs",
            "src/phenomena/heliacal/optical.rs",
            "src/phenomena/heliacal/pheno_window.rs",
            "src/phenomena/heliacal/provider.rs",
            "src/phenomena/heliacal/search.rs",
            "src/runtime/api/native_heliacal.rs",
            "src/runtime/api/numeric_heliacal.rs",
            "src/runtime/context.rs",
            "src/runtime/config.rs",
            "src/tables/star_catalog_tables.rs",
        ] {
            let path = source_root.join(name);
            before.insert(path.clone(), fingerprint(&path));
        }
        let executable = std::env::current_exe().unwrap();
        before.insert(executable.clone(), fingerprint(&executable));
        writeln!(stream, "{}", json!({"kind":"header", "scope":if current_input.is_some(){"current-provider-legacy-heliacal-apparition-gates"}else{"twelve-original-heliacal-apparition-gates"}, "public_qualified":false,"inputs_before":before,"environment":{
            "network":"sealed", "auto_spk_download":false, "iers_auto_download":false,
            "data_root_per_request":"explicit historical config, otherwise target/v4/comparison-data"},
            "historical_binary_claim":false,"returned_flags":"not part of heliacal_ut result"})).unwrap();
        stream.flush().unwrap();
        let mut exact = 0;
        for row in rows {
            let runtime = runtime(&row);
            let args = row["args"].as_array().unwrap();
            let flags = row["kwargs"]["flags"].as_i64().unwrap_or(2) as i32;
            observation::begin();
            let trace = crate::trace::start_session();
            let result = runtime.shared_session().and_then(|session| {
                crate::binding_heliacal::heliacal_ut_words(
                    &session,
                    args[0].as_f64().unwrap(),
                    &words(&args[1]),
                    &words(&args[2]),
                    &words(&args[3]),
                    args[4].as_str().unwrap(),
                    args[5].as_i64().unwrap() as i32,
                    flags,
                )
            });
            let events = observation::finish();
            let providers = crate::trace::get_trace_results();
            assert!(trace.reset());
            let expected = words(&row["result"]);
            assert_eq!(expected.len(), 3);
            let (actual, error, matches) = match result {
                Ok(actual) => {
                    let matches = actual
                        .iter()
                        .zip(&expected)
                        .all(|(a, b)| a.to_bits() == b.to_bits());
                    (json!(actual.map(observation::word)), Value::Null, matches)
                }
                Err(error) => (
                    Value::Null,
                    json!({"display":error.to_string(),"debug":format!("{error:?}")}),
                    false,
                ),
            };
            exact += usize::from(matches);
            writeln!(stream, "{}", json!({"kind":"request", "request":row,"flags":flags,"actual_words":actual,"error":error,"historical_result_word_exact":matches,"events":events,"actual_providers":providers,"clock_scope":"actual UT search epochs and explicit configured scalar clock; internal observer TT is not observed by this hook"})).unwrap();
            stream.flush().unwrap();
        }
        let after: BTreeMap<_, _> = before
            .keys()
            .map(|path| (path.clone(), fingerprint(path)))
            .collect();
        let unchanged = before == after;
        writeln!(stream, "{}", json!({"kind":"footer", "requests":requests, "historical_word_exact":exact,"inputs_unchanged":unchanged,"inputs_after":after,"public_qualified":false})).unwrap();
        stream.flush().unwrap();
        assert!(unchanged, "consumed inputs changed");
        assert_eq!(
            exact, requests,
            "actual errors/differences preserved in observation output"
        );
    }
}
