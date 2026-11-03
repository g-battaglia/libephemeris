//! Contact refinement and per-contact optional-geometry composition.
//!
//! These solvers deliberately retain distinct bracket, endpoint and degradation
//! contracts. Sharing their numerical policies would change callback order and
//! the meaning of an unresolved contact (`0.0` or `None`).

use super::{
    BESSELIAN_CONTACT_MAX_WINDOW_DAYS, BESSELIAN_CONTACT_STEP_DAYS, CONTACT_BISECT_TOLERANCE_DAYS,
    DETAILS_CONTACT_ITERATIONS, DetailsContactBlock, EclipseHowDetailsEphemeris, Error,
    GOLDEN_SECTION_ITERATIONS, LOCAL_CONTACT_BISECT_ITERATIONS, LOCAL_CONTACT_SEP_TOLERANCE_DEG,
    LOCAL_CONTACT_VALID_TOLERANCE_DEG, MemoTriple, Result, eclipse_leaf_reraises,
};

/// Validate a contact time: returns `0.0` for non-positive times or times
/// farther than `max_offset_days` from maximum.
pub fn validate_contact_time(jd_contact: f64, jd_max: f64, max_offset_days: f64) -> f64 {
    if jd_contact <= 0.0 {
        return 0.0;
    }
    if (jd_contact - jd_max).abs() > max_offset_days {
        return 0.0;
    }
    jd_contact
}

/// Bisect a sign change of `residual` on `[lo, hi]`: exact-zero endpoints return the endpoint,
/// an unbracketed interval returns the `0.0` sentinel, otherwise 60
/// halvings to the `1e-8`-day bracket. Provider errors propagate (each
/// caller maps them per the release `try` placement).
pub fn bisect_sign_change(
    residual: &mut impl FnMut(f64) -> Result<f64>,
    lo: f64,
    hi: f64,
) -> Result<f64> {
    let f_lo = residual(lo)?;
    let f_hi = residual(hi)?;
    if f_lo == 0.0 {
        return Ok(lo);
    }
    if f_hi == 0.0 {
        return Ok(hi);
    }
    if (f_lo > 0.0) == (f_hi > 0.0) {
        return Ok(0.0);
    }
    let (mut lo, mut hi, mut f_lo) = (lo, hi, f_lo);
    for _ in 0..GOLDEN_SECTION_ITERATIONS {
        let mid = 0.5 * (lo + hi);
        let f_mid = residual(mid)?;
        if (f_mid > 0.0) == (f_lo > 0.0) {
            lo = mid;
            f_lo = f_mid;
        } else {
            hi = mid;
        }
        if hi - lo < CONTACT_BISECT_TOLERANCE_DAYS {
            break;
        }
    }
    Ok(0.5 * (lo + hi))
}

/// Read a Besselian residual, retaining the protected key, value and
/// arithmetic exception hierarchies. A caught failure returns `None`;
/// coverage misses, corrupt ephemerides and other failures propagate.
pub(in crate::events) fn contact_provider_value(
    residual: &mut impl FnMut(f64) -> Result<f64>,
    jd: f64,
) -> Result<Option<f64>> {
    match residual(jd) {
        Ok(value) => Ok(Some(value)),
        Err(err) => {
            if besselian_contact_caught(&err) {
                Ok(None)
            } else {
                Err(err)
            }
        }
    }
}

/// Contact-solve core over optional residual reads (shared structure
/// with [`super::solve_besselian_contact`]): `None` is the degraded-sentinel
/// path. A `None` at the initial read or within bisection returns `0.0`;
/// bracket expansion skips that sample. Uncaught errors propagate.
pub(in crate::events) fn solve_contact_opt(
    read: &mut impl FnMut(f64) -> Result<Option<f64>>,
    jd_max: f64,
    search_before: bool,
    max_window_days: f64,
) -> Result<f64> {
    let f_max = match read(jd_max)? {
        Some(value) => value,
        None => return Ok(0.0),
    };
    if f_max >= 0.0 {
        return Ok(0.0);
    }
    let sign = if search_before { -1.0 } else { 1.0 };
    let mut far = 0.0;
    let mut step = BESSELIAN_CONTACT_STEP_DAYS;
    while step <= max_window_days {
        let t = jd_max + sign * step;
        if let Some(ft) = read(t)? {
            if ft > 0.0 {
                far = t;
                break;
            }
        }
        step *= 1.6;
    }
    if far == 0.0 {
        return Ok(0.0);
    }
    let (mut lo, mut hi) = if search_before {
        (far, jd_max)
    } else {
        (jd_max, far)
    };
    // Inline `_root_bisect` over optional reads (same endpoints,
    // sentinel and iteration structure as [`bisect_sign_change`]); a
    // degraded read inside the bracket degrades to `0.0`.
    let f_lo = match read(lo)? {
        Some(value) => value,
        None => return Ok(0.0),
    };
    let f_hi = match read(hi)? {
        Some(value) => value,
        None => return Ok(0.0),
    };
    if f_lo == 0.0 {
        return Ok(lo);
    }
    if f_hi == 0.0 {
        return Ok(hi);
    }
    if (f_lo > 0.0) == (f_hi > 0.0) {
        return Ok(0.0);
    }
    let mut f_lo = f_lo;
    for _ in 0..GOLDEN_SECTION_ITERATIONS {
        let mid = 0.5 * (lo + hi);
        let f_mid = match read(mid)? {
            Some(value) => value,
            None => return Ok(0.0),
        };
        if (f_mid > 0.0) == (f_lo > 0.0) {
            lo = mid;
            f_lo = f_mid;
        } else {
            hi = mid;
        }
        if hi - lo < CONTACT_BISECT_TOLERANCE_DAYS {
            break;
        }
    }
    Ok(0.5 * (lo + hi))
}

impl<F> MemoTriple<'_, F>
where
    F: FnMut(f64) -> Result<(f64, f64, f64)>,
{
    pub(in crate::events) fn slot(&mut self, slot: usize, jd: f64) -> Result<Option<f64>> {
        for (key, triple) in &self.memo {
            if *key == jd.to_bits() {
                return Ok(Some(match slot {
                    0 => triple.0,
                    1 => triple.1,
                    _ => triple.2,
                }));
            }
        }
        match (self.inner)(jd) {
            Ok(triple) => {
                let value = match slot {
                    0 => triple.0,
                    1 => triple.1,
                    _ => triple.2,
                };
                self.memo.push((jd.to_bits(), triple));
                Ok(Some(value))
            }
            Err(err) => {
                if besselian_contact_caught(&err) {
                    Ok(None)
                } else {
                    Err(err)
                }
            }
        }
    }
}

/// Global first penumbral tangency before the supplied approximate maximum.
///
/// `jd_max` and the result are UT1 Julian days. `residuals(jd_ut1)` supplies
/// `(penumbral_edge, umbral_edge, shadow_axis)` signed residuals in kilometres,
/// negative during the corresponding intersection. This entry selects the
/// first component, with a bounded expanding bracket and bisection search.
/// Successful triples are memoized by exact epoch bits for this call only.
///
/// Returns zero if the feature is absent at maximum, no bracket is found, or
/// a selected acquisition/domain failure degrades the solve. Range/corruption
/// and uncaught errors propagate. Zero therefore means absent **or unresolved**,
/// not proof of physical nonexistence. The callback owns finite, consistent
/// residuals and its epoch domain; this reducer adds no eager epoch validation.
pub fn solar_contact_c1_via(
    jd_max: f64,
    residuals: &mut impl FnMut(f64) -> Result<(f64, f64, f64)>,
) -> Result<f64> {
    let mut memo = MemoTriple {
        inner: residuals,
        memo: Vec::new(),
    };
    let mut read = |jd: f64| memo.slot(0, jd);
    solve_contact_opt(&mut read, jd_max, true, BESSELIAN_CONTACT_MAX_WINDOW_DAYS)
}

/// Global first umbral/antumbral tangency before maximum, in UT1 Julian days.
///
/// Uses component 1 of the residual tuple documented on [`solar_contact_c1_via`].
/// Its callback units, bounded search, call-local memoization, propagated errors
/// and absent-or-unresolved zero sentinel apply unchanged. This is a global
/// Earth tangency, not the observer's local second contact.
pub fn solar_contact_c2_via(
    jd_max: f64,
    residuals: &mut impl FnMut(f64) -> Result<(f64, f64, f64)>,
) -> Result<f64> {
    let mut memo = MemoTriple {
        inner: residuals,
        memo: Vec::new(),
    };
    let mut read = |jd: f64| memo.slot(1, jd);
    solve_contact_opt(&mut read, jd_max, true, BESSELIAN_CONTACT_MAX_WINDOW_DAYS)
}

/// Global last umbral/antumbral tangency after maximum, in UT1 Julian days.
///
/// Uses component 1 of the residual tuple documented on [`solar_contact_c1_via`].
/// Its callback/domain/error contract and zero sentinel apply unchanged. This
/// is not the observer's local third contact.
pub fn solar_contact_c3_via(
    jd_max: f64,
    residuals: &mut impl FnMut(f64) -> Result<(f64, f64, f64)>,
) -> Result<f64> {
    let mut memo = MemoTriple {
        inner: residuals,
        memo: Vec::new(),
    };
    let mut read = |jd: f64| memo.slot(1, jd);
    solve_contact_opt(&mut read, jd_max, false, BESSELIAN_CONTACT_MAX_WINDOW_DAYS)
}

/// Global last penumbral tangency after maximum, in UT1 Julian days.
///
/// Uses component 0 of the residual tuple documented on [`solar_contact_c1_via`].
/// Its callback/domain/error contract, call-local memoization and zero sentinel
/// apply unchanged. This is last contact anywhere on Earth, not at one site.
pub fn solar_contact_c4_via(
    jd_max: f64,
    residuals: &mut impl FnMut(f64) -> Result<(f64, f64, f64)>,
) -> Result<f64> {
    let mut memo = MemoTriple {
        inner: residuals,
        memo: Vec::new(),
    };
    let mut read = |jd: f64| memo.slot(0, jd);
    solve_contact_opt(&mut read, jd_max, false, BESSELIAN_CONTACT_MAX_WINDOW_DAYS)
}

/// One legacy local contact solve:
/// bisection of `separation - target_sep` with the
/// achievability pre-check (an unachievable target reports `0.0`
/// without iterating), `50` iterations, the `1e-6` deg match stop,
/// the `1e-8` day width stop, and the `0.001` deg post-solve
/// validation (`0.0` when unreliable). `increasing` selects the
/// exit-side (`false`) or entry-side (`true`) branch update.
pub fn find_local_contact_time(
    separation_deg: &mut impl FnMut(f64) -> Result<f64>,
    jd_start: f64,
    jd_end: f64,
    target_sep_deg: f64,
    increasing: bool,
) -> Result<f64> {
    for value in [jd_start, jd_end, target_sep_deg] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    if jd_end <= jd_start {
        return Err(Error::InvalidSearchInterval {
            a: jd_start,
            b: jd_end,
        });
    }
    let sep_start = separation_deg(jd_start)?;
    let sep_end = separation_deg(jd_end)?;
    if !sep_start.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: sep_start });
    }
    if !sep_end.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: sep_end });
    }
    // Achievability pre-check: a midpoint probe
    // distinguishes "contact never reached" (`0.0`) from a bracket the
    // plain endpoint test would misread.
    if increasing {
        if sep_start > target_sep_deg || sep_end < target_sep_deg {
            let sep_mid = separation_deg(0.5 * (jd_start + jd_end))?;
            if !sep_mid.is_finite() {
                return Err(Error::NonFiniteSearchValue { x: sep_mid });
            }
            if sep_start.min(sep_end).min(sep_mid) > target_sep_deg {
                return Ok(0.0);
            }
        }
    } else if sep_start < target_sep_deg || sep_end > target_sep_deg {
        let sep_mid = separation_deg(0.5 * (jd_start + jd_end))?;
        if !sep_mid.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: sep_mid });
        }
        if sep_start.min(sep_end).min(sep_mid) > target_sep_deg {
            return Ok(0.0);
        }
    }
    let (mut jd_low, mut jd_high) = (jd_start, jd_end);
    for _ in 0..LOCAL_CONTACT_BISECT_ITERATIONS {
        let jd_mid = 0.5 * (jd_low + jd_high);
        let sep = separation_deg(jd_mid)?;
        if !sep.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: sep });
        }
        if (sep - target_sep_deg).abs() < LOCAL_CONTACT_SEP_TOLERANCE_DEG {
            return Ok(jd_mid);
        }
        if increasing {
            if sep < target_sep_deg {
                jd_low = jd_mid;
            } else {
                jd_high = jd_mid;
            }
        } else if sep > target_sep_deg {
            jd_low = jd_mid;
        } else {
            jd_high = jd_mid;
        }
        if jd_high - jd_low < CONTACT_BISECT_TOLERANCE_DAYS {
            break;
        }
    }
    let result = 0.5 * (jd_low + jd_high);
    let final_sep = separation_deg(result)?;
    if !final_sep.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: final_sep });
    }
    if (final_sep - target_sep_deg).abs() > LOCAL_CONTACT_VALID_TOLERANCE_DEG {
        return Ok(0.0);
    }
    Ok(result)
}

/// Refine a monotone local disc-contact equation to representable time.
///
/// The caller supplies a finite separation in degrees and a known entry or
/// exit bracket in UT1 Julian days. Bisect until the separation equals its
/// target or no interior binary64 epoch remains; return the evaluated epoch
/// with the smallest absolute angular residual. A fixed angular match limit
/// is insufficient near grazing contacts, where it corresponds to a large
/// timing error. At most sixty samples are taken. Provider failures and
/// nonfinite sampled separations propagate without substituting an event.
pub fn details_contact_time(
    separation_deg: &mut impl FnMut(f64) -> Result<f64>,
    jd_start: f64,
    jd_end: f64,
    target_sep_deg: f64,
    is_increasing: bool,
) -> Result<f64> {
    let (mut jd_low, mut jd_high) = (jd_start, jd_end);
    let mut best = None;
    for _ in 0..DETAILS_CONTACT_ITERATIONS {
        let jd_mid = jd_low + (jd_high - jd_low) / 2.0;
        if jd_mid == jd_low || jd_mid == jd_high {
            break;
        }
        let sep = separation_deg(jd_mid)?;
        if !sep.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd_mid });
        }
        let residual = sep - target_sep_deg;
        if best.is_none_or(|(_, error): (f64, f64)| residual.abs() < error) {
            best = Some((jd_mid, residual.abs()));
        }
        if residual == 0.0 {
            return Ok(jd_mid);
        }
        if (residual < 0.0) == is_increasing {
            jd_low = jd_mid;
        } else {
            jd_high = jd_mid;
        }
    }
    if let Some((epoch, _)) = best {
        Ok(epoch)
    } else {
        // An already adjacent bracket has no interior representable epoch.
        // Evaluate its actual endpoints rather than synthesize an epoch.
        let low = (separation_deg(jd_low)? - target_sep_deg).abs();
        let high = (separation_deg(jd_high)? - target_sep_deg).abs();
        if !low.is_finite() || !high.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd_low });
        }
        Ok(if low <= high { jd_low } else { jd_high })
    }
}

/// Solve a disc contact and evaluate the solar horizon and position angle.
/// Coverage failures propagate; unavailable optional geometry is omitted.
pub(in crate::events) fn details_contact_block(
    eph: &impl EclipseHowDetailsEphemeris,
    separation: &mut impl FnMut(f64) -> Result<f64>,
    lo_jd: f64,
    hi_jd: f64,
    threshold_deg: f64,
    increasing: bool,
) -> Result<Option<DetailsContactBlock>> {
    let jd = match details_contact_time(separation, lo_jd, hi_jd, threshold_deg, increasing) {
        Ok(jd) => jd,
        Err(e) if eclipse_leaf_reraises(&e) => return Err(e),
        Err(e) if details_contact_caught(&e) => return Ok(None),
        Err(e) => return Err(e),
    };
    let alt_az = match eph.details_sun_alt_az_deg(jd) {
        Ok(pair) => pair,
        Err(e) if eclipse_leaf_reraises(&e) => return Err(e),
        Err(e) if details_contact_caught(&e) => {
            return Ok(Some(DetailsContactBlock {
                jd,
                alt_az: None,
                position_angle_deg: None,
            }));
        }
        Err(e) => return Err(e),
    };
    let position_angle_deg = match eph.details_position_angle_deg(jd) {
        Ok(angle) => Some(angle),
        Err(e) if eclipse_leaf_reraises(&e) => return Err(e),
        Err(e) if details_contact_caught(&e) => None,
        Err(e) => return Err(e),
    };
    Ok(Some(DetailsContactBlock {
        jd,
        alt_az: Some(alt_az),
        position_angle_deg,
    }))
}

/// Geometric-domain, non-finite and arithmetic failures which leave a
/// contact unresolved. Data-range and corruption errors always propagate.
fn details_contact_caught(error: &Error) -> bool {
    matches!(
        error.primary_error(),
        Error::MathDomain { .. }
            | Error::NonFiniteArgument { .. }
            | Error::IntegerOverflow
            | Error::RationalDivisionByZero
            | Error::NumericalOverflow { .. }
            | Error::LebBodyNotFound { .. }
    )
}

/// Besselian phase failures may leave an individual contact unresolved;
/// ephemeris range and corrupt-data errors remain fatal.
pub(in crate::events) fn besselian_contact_caught(error: &Error) -> bool {
    !eclipse_leaf_reraises(error.primary_error()) && details_contact_caught(error)
}
