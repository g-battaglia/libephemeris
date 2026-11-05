//! Evaluation implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

impl ChebState {
    /// Position in astronomical units.
    pub fn position_au(&self) -> [f64; 3] {
        self.position_km.map(|v| v / AUNIT_TO_KM)
    }

    /// Velocity in AU/day.
    pub fn velocity_au_per_day(&self) -> [f64; 3] {
        self.velocity_km_s
            .map(|v| v * SECONDS_PER_DAY / AUNIT_TO_KM)
    }
}

/// Evaluate a decoded type2/3 segment at TDB seconds past J2000.
/// Directory coverage is closed. The final endpoint belongs to the last
/// record; other points outside the directory are errors, never extrapolated.
/// Each record's actual midpoint/radius defines its polynomial argument.
/// Type2 velocity analytically differentiates the polynomial in km/s; type3
/// uses its stored velocity coefficients. Inputs and coefficients are finite.
pub fn eval_chebyshev(segment: &ChebSegment, et: f64) -> Result<ChebState> {
    if !et.is_finite() {
        return Err(Error::NonFiniteArgument { value: et });
    }
    let end = chebyshev_directory_end(segment)?;
    if et < segment.init_et || et > end {
        return Err(Error::SpkOutOfRange {
            target: 0,
            et,
            start: segment.init_et,
            end,
        });
    }
    let relative = et - segment.init_et;
    let index = ((relative / segment.int_len).floor() as usize).min(segment.n_records - 1);
    eval_chebyshev_record(segment, index, relative - index as f64 * segment.int_len)
}

/// Validate the finite closed directory span before selecting a record.
fn chebyshev_directory_end(segment: &ChebSegment) -> Result<f64> {
    let end = segment.init_et + segment.n_records as f64 * segment.int_len;
    let components = if segment.seg_type == TYPE_CHEBYSHEV_POS {
        3
    } else {
        6
    };
    if segment.record_words < 2 + components
        || (segment.record_words - 2) % components != 0
        || segment.n_records.checked_mul(segment.record_words) != Some(segment.words.len())
        || segment.n_records == 0
        || !segment.init_et.is_finite()
        || !segment.int_len.is_finite()
        || segment.int_len <= 0.0
        || !end.is_finite()
    {
        return Err(Error::SpkCorrupt {
            what: "SPK Chebyshev directory".to_owned(),
            offset: segment.data_offset,
        });
    }
    Ok(end)
}

/// Evaluate one checked record using its stored midpoint and radius.
/// The directory offset remains separate from its large origin until the
/// small distance to the actual record midpoint has been formed. Positions
/// are km; type2 differentiated velocities are km/s, while type3 stores its
/// velocity coefficients directly. Nonfinite data and results are errors.
fn eval_chebyshev_record(segment: &ChebSegment, index: usize, offset: f64) -> Result<ChebState> {
    eval_chebyshev_record_with_domain(segment, index, offset, None).map(|(state, _, _)| state)
}

/// Evaluate once, retaining the open-domain flag and a separate finite dP/dTDB.
/// Type3's original velocity still comes from its independent coefficient sets.
fn eval_chebyshev_record_with_domain(
    segment: &ChebSegment,
    index: usize,
    offset: f64,
    seconds: Option<(f64, f64)>,
) -> Result<(ChebState, bool, Option<[f64; 3]>)> {
    if !offset.is_finite() {
        return Err(Error::NonFiniteArgument { value: offset });
    }
    let base = index
        .checked_mul(segment.record_words)
        .ok_or_else(|| Error::SpkCorrupt {
            what: "SPK Chebyshev record".to_owned(),
            offset: segment.data_offset,
        })?;
    let rec = segment
        .words
        .get(base..base + segment.record_words)
        .ok_or_else(|| Error::SpkCorrupt {
            what: "SPK Chebyshev record".to_owned(),
            offset: segment.data_offset,
        })?;
    let directory = super::decode::ChebyshevDirectory {
        init_et: segment.init_et,
        int_len: segment.int_len,
        record_words: segment.record_words,
        n_records: segment.n_records,
    };
    eval_selected_chebyshev_record(
        &directory,
        segment.seg_type,
        rec,
        index,
        offset,
        segment.data_offset + base * 8,
        seconds,
    )
}

/// Evaluate one already selected type-2/3 record with the shared arithmetic.
/// The validated directory and record retain their original origin/index so
/// positional reads never collapse a small time offset into a large epoch.
pub(super) fn eval_selected_chebyshev_record(
    directory: &super::decode::ChebyshevDirectory,
    seg_type: i32,
    rec: &[f64],
    index: usize,
    offset: f64,
    data_offset: usize,
    seconds: Option<(f64, f64)>,
) -> Result<(ChebState, bool, Option<[f64; 3]>)> {
    let mid = rec[0];
    let radius = rec[1];
    if !mid.is_finite() || !radius.is_finite() || radius <= 0.0 {
        return Err(Error::SpkCorrupt {
            what: "SPK Chebyshev record".to_owned(),
            offset: data_offset,
        });
    }
    let delta = crate::types::TwoPartEpoch {
        hi: directory.init_et - mid,
        lo: offset,
    }
    .add_offset(index as f64 * directory.int_len)?
    .value();
    let tau = delta / radius;
    if !tau.is_finite() || tau.abs() > 1.0 {
        return Err(Error::SpkOutOfRange {
            target: 0,
            et: mid + delta,
            start: mid - radius,
            end: mid + radius,
        });
    }
    let (state, position_tangent) = eval_chebyshev_values_and_tangent(seg_type, rec, tau, radius)?;
    let stored_interior = seconds.is_none_or(|(hi, lo)| {
        let (_, start_upper) = outward_sum(mid, -radius);
        let (end_lower, _) = outward_sum(mid, radius);
        split_seconds_interior(hi, lo, start_upper, end_lower)
    });
    Ok((
        state,
        offset > 0.0 && offset < directory.int_len && tau.abs() < 1.0 && stored_interior,
        position_tangent,
    ))
}

/// Evaluate a checked type-2/3 record for a bounded source reader.
///
/// The reader owns its directory and epoch selection. Positions are km,
/// velocities km/s, `radius` seconds and `tau` dimensionless. Invalid layouts,
/// nonfinite arguments and extrapolation are rejected before evaluation.
#[cfg(any(test, feature = "nbody"))]
pub(super) fn eval_chebyshev_record_values(
    seg_type: i32,
    rec: &[f64],
    tau: f64,
    radius: f64,
) -> Result<ChebState> {
    let components = match seg_type {
        TYPE_CHEBYSHEV_POS => 3,
        TYPE_CHEBYSHEV_POS_VEL => 6,
        _ => return Err(Error::SpkUnsupportedType { seg_type }),
    };
    if rec.len() < 2 + components || (rec.len() - 2) % components != 0 || radius <= 0.0 {
        return Err(Error::SpkCorrupt {
            what: "SPK checked coefficient record".into(),
            offset: 0,
        });
    }
    for value in [tau, radius] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    if tau.abs() > 1.0 {
        return Err(Error::SpkOutOfRange {
            target: 0,
            et: rec[0] + tau * radius,
            start: rec[0] - radius,
            end: rec[0] + radius,
        });
    }
    eval_chebyshev_values_and_tangent(seg_type, rec, tau, radius).map(|(state, _)| state)
}

/// Shared coefficient recurrence after the owning reader's admission checks.
/// Type-3 velocity remains distinct from the analytic position tangent.
fn eval_chebyshev_values_and_tangent(
    seg_type: i32,
    rec: &[f64],
    tau: f64,
    radius: f64,
) -> Result<(ChebState, Option<[f64; 3]>)> {
    let components = if seg_type == TYPE_CHEBYSHEV_POS { 3 } else { 6 };
    let sets = (rec.len() - 2) / components;
    let set = |k: usize| {
        let off = 2 + k * sets;
        &rec[off..off + sets]
    };
    let mut position_km = [0.0; 3];
    let mut velocity_km_s = [0.0; 3];
    let mut position_tangent_km_s = [0.0; 3];
    for c in 0..3 {
        let (value, derivative) = spk_chebyshev_value_and_derivative(set(c), tau)?;
        position_km[c] = value;
        position_tangent_km_s[c] = derivative / radius;
        velocity_km_s[c] = if seg_type == TYPE_CHEBYSHEV_POS {
            position_tangent_km_s[c]
        } else {
            spk_chebyshev_value_and_derivative(set(c + 3), tau)?.0
        };
    }
    for value in position_km.into_iter().chain(velocity_km_s) {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    Ok((
        ChebState {
            position_km,
            velocity_km_s,
        },
        position_tangent_km_s
            .iter()
            .all(|value| value.is_finite())
            .then_some(position_tangent_km_s),
    ))
}

/// Evaluate a type 2/3 segment at the uncollapsed TDB JD `jd + fraction`.
///
/// The two JD parts are separately converted to seconds. Three quotient/
/// remainder operations combine their interval indices and offsets without
/// first adding a tiny fraction to a large epoch. Directory coverage is closed;
/// its final endpoint uses the last record. No extrapolation is permitted.
/// Outputs are km and km/s; coefficient and record errors propagate.
pub fn eval_chebyshev_two_part_jd(
    segment: &ChebSegment,
    jd: f64,
    fraction: f64,
) -> Result<ChebState> {
    eval_chebyshev_two_part_jd_with_domain(segment, jd, fraction).map(|(state, _, _)| state)
}

/// The same selected state plus an open directory/stored-polynomial domain.
/// Joins and the final closed endpoint are uncertified; no C1 match is assumed.
pub(super) fn eval_chebyshev_two_part_jd_with_domain(
    segment: &ChebSegment,
    jd: f64,
    fraction: f64,
) -> Result<(ChebState, bool, Option<[f64; 3]>)> {
    let (index, offset, directory_interior, hi, lo) =
        chebyshev_record_at_two_part_jd(segment, jd, fraction)?;
    let (state, polynomial_interior, tangent) =
        eval_chebyshev_record_with_domain(segment, index, offset, Some((hi, lo)))?;
    Ok((state, directory_interior && polynomial_interior, tangent))
}

/// The exact directory-selected record and its small epoch offset.
/// Shared by evaluation and continuous-record admission.
fn chebyshev_record_at_two_part_jd(
    segment: &ChebSegment,
    jd: f64,
    fraction: f64,
) -> Result<(usize, f64, bool, f64, f64)> {
    let seconds = jd_parts_seconds(jd, fraction)?;
    chebyshev_directory_end(segment)?;
    chebyshev_record_coordinates_prepared(
        segment.init_et,
        segment.int_len,
        segment.n_records,
        seconds,
    )
}

/// Select a record with the same quotient/remainder order as owned reads.
#[cfg(feature = "nbody")]
pub(super) fn chebyshev_record_coordinates(
    init_et: f64,
    int_len: f64,
    n_records: usize,
    jd: f64,
    fraction: f64,
) -> Result<(usize, f64, bool, f64, f64)> {
    let seconds = jd_parts_seconds(jd, fraction)?;
    chebyshev_record_coordinates_prepared(init_et, int_len, n_records, seconds)
}

fn chebyshev_record_coordinates_prepared(
    init_et: f64,
    int_len: f64,
    n_records: usize,
    (a, b, hi, lo): (f64, f64, f64, f64),
) -> Result<(usize, f64, bool, f64, f64)> {
    let end = init_et + n_records as f64 * int_len;
    if !split_seconds_covers(hi, lo, init_et, end) {
        return Err(Error::SpkOutOfRange {
            target: 0,
            et: hi + lo,
            start: init_et,
            end,
        });
    }
    // Reduce the two epoch components separately, retaining a sub-ULP
    // fractional offset even at epochs many billions of seconds from J2000.
    let interval = int_len;
    let first = a - init_et;
    let q1 = first.div_euclid(interval);
    let q2 = b.div_euclid(interval);
    let remainder = first.rem_euclid(interval) + b.rem_euclid(interval);
    let q3 = remainder.div_euclid(interval);
    let quotient = q1 + q2 + q3;
    let mut offset = remainder.rem_euclid(interval);
    if !quotient.is_finite() || quotient < 0.0 || quotient > n_records as f64 {
        return Err(Error::SpkOutOfRange {
            target: 0,
            et: hi + lo,
            start: init_et,
            end,
        });
    }
    let mut index = quotient as usize;
    if index == n_records {
        index -= 1;
        offset += interval;
    }
    let interior = directory_record_interior(init_et, int_len, index, hi, lo);
    Ok((index, offset, interior, hi, lo))
}

/// Outward enclosure of a real sum of finite binary64 words.
/// Nonfinite bounds prevent certification without introducing errors.
fn outward_sum(a: f64, b: f64) -> (f64, f64) {
    let sum = a + b;
    (sum.next_down(), sum.next_up())
}

/// Enclose init + k*width for the validated positive directory width,
/// including conversion of an arbitrary usize k.
fn directory_boundary_enclosure(init: f64, width: f64, index: usize) -> (f64, f64) {
    let rounded = index as f64;
    let (lower, upper) = if (index as u128) <= (1_u128 << 53) {
        // Every nonnegative integer through 2^53 is exactly representable.
        (rounded, rounded)
    } else {
        (rounded.next_down(), rounded.next_up())
    };
    let low_product = (lower * width).next_down();
    let high_product = (upper * width).next_up();
    (
        outward_sum(init, low_product).0,
        outward_sum(init, high_product).1,
    )
}

fn split_seconds_interior(hi: f64, lo: f64, start: f64, end: f64) -> bool {
    (hi > start || (hi == start && lo > 0.0)) && (hi < end || (hi == end && lo < 0.0))
}

fn directory_record_interior(init: f64, width: f64, index: usize, hi: f64, lo: f64) -> bool {
    let Some(next) = index.checked_add(1) else {
        return false;
    };
    let (_, start_upper) = directory_boundary_enclosure(init, width, index);
    let (end_lower, _) = directory_boundary_enclosure(init, width, next);
    // Adjacent words inside either outward endpoint enclosure remain
    // uncertified too; no coefficient or source-state value is changed.
    split_seconds_interior(hi, lo, start_upper, end_lower)
}

#[cfg(test)]
mod directory_domain_tests {
    use super::{directory_boundary_enclosure, directory_record_interior};

    #[test]
    fn invalid_epoch_precedes_an_invalid_owned_directory() {
        let segment = super::ChebSegment {
            seg_type: 2,
            init_et: 0.0,
            int_len: 0.0,
            record_words: 0,
            n_records: 0,
            words: Vec::new(),
            data_offset: 0,
        };
        for (jd, fraction) in [(f64::NAN, 0.0), (2451545.0, f64::INFINITY)] {
            assert!(matches!(
                super::eval_chebyshev_two_part_jd(&segment, jd, fraction),
                Err(super::Error::NonFiniteArgument { .. })
            ));
        }
        assert!(matches!(
            super::eval_chebyshev_two_part_jd(&segment, 2451545.0, 0.0),
            Err(super::Error::SpkCorrupt { .. })
        ));
    }

    #[test]
    fn non_dyadic_metadata_seam_and_its_uncertain_neighbors_are_not_certified() {
        // The actual 0.1 binary64 word is 3602879701896397 / 2^55.
        // init + 4*width is therefore exactly 0.5 + 2^-55, not the
        // rounded scalar boundary 0.5. Its normalized two-word epoch would
        // pass a naive strict comparison with that rounded lower bound.
        let init = 0.1;
        let width = 0.1;
        let hi = 0.5;
        let lo = 2.0_f64.powi(-55);
        assert_eq!(init + 4.0 * width, hi);
        assert!(lo > 0.0);
        for index in [3, 4] {
            assert!(!directory_record_interior(init, width, index, hi, lo));
        }
        assert!(!directory_record_interior(
            init,
            width,
            3,
            hi.next_down(),
            0.0
        ));
        assert!(!directory_record_interior(
            init,
            width,
            4,
            hi.next_up(),
            0.0
        ));
        assert!(directory_record_interior(init, width, 3, 0.45, 0.0));
        assert!(directory_record_interior(init, width, 4, 0.55, 0.0));

        if usize::BITS > 53 {
            let index = usize::try_from((1_u128 << 53) + 1).unwrap();
            let (lower, upper) = directory_boundary_enclosure(0.0, 1.0, index);
            let rounded = index as f64;
            assert!(lower < rounded && rounded < upper);
        }
    }
}

/// Offsets in TDB days that stay on the selected type2/3 polynomial.
/// Directory joins belong to the following record; the last endpoint remains
/// closed. The stored polynomial domain can be narrower than its directory.
pub(super) fn chebyshev_record_offsets(
    segment: &ChebSegment,
    jd: f64,
    fraction: f64,
) -> Result<(f64, f64)> {
    let (index, _, _, _, _) = chebyshev_record_at_two_part_jd(segment, jd, fraction)?;
    let (_, _, hi, lo) = jd_parts_seconds(jd, fraction)?;
    let base = index * segment.record_words;
    let midpoint = segment.words[base];
    let radius = segment.words[base + 1];
    if !midpoint.is_finite() || !radius.is_finite() || radius <= 0.0 {
        return Err(Error::SpkCorrupt {
            what: "SPK Chebyshev record".to_owned(),
            offset: segment.data_offset + base * 8,
        });
    }
    let start = segment.init_et + index as f64 * segment.int_len;
    let end = segment.init_et + (index + 1) as f64 * segment.int_len;
    let relative = |boundary: f64| ((boundary - hi) - lo) / SECONDS_PER_DAY;
    let lower = relative(start.max(midpoint - radius));
    let mut upper = relative(end.min(midpoint + radius));
    if index + 1 != segment.n_records && end <= midpoint + radius {
        upper = upper.next_down();
    }
    Ok((lower, upper))
}

/// Offsets in TDB days owned by the selected type21 difference line.
/// The descriptor clips the first/last lines; an interior final epoch belongs
/// to its following line, independently of the polynomial reference epoch.
pub(super) fn emda_record_offsets(
    segment: &EmdaSegment,
    descriptor: Descriptor,
    jd: f64,
    fraction: f64,
) -> Result<(f64, f64)> {
    let et = emda_epoch_seconds(jd, fraction)?;
    let index = emda_record_index(segment, et)?;
    let (_, _, hi, lo) = jd_parts_seconds(jd, fraction)?;
    let relative = |boundary: f64| ((boundary - hi) - lo) / SECONDS_PER_DAY;
    let lower = if index == 0 {
        relative(descriptor.start_et)
    } else {
        emda_boundary_offsets(segment.epochs[index - 1], jd, fraction)?
            .1
            .max(relative(descriptor.start_et))
    };
    let upper = if index + 1 == segment.n_records {
        relative(descriptor.end_et)
    } else {
        emda_boundary_offsets(segment.epochs[index], jd, fraction)?
            .0
            .min(relative(descriptor.end_et))
    };
    Ok((lower, upper))
}

/// Adjacent offset words on either side of the *actual scalar* type21
/// selector's boundary. Normalizing JD parts before the seconds product can
/// round to a final epoch from either side of its real-valued TDB boundary.
/// Locate that selector cell rather than mixing it with the split descriptor
/// test. The bracket follows the two binary64 additions/products' forward
/// error scale; no state, derivative or comparison residual enters it.
fn emda_boundary_offsets(boundary: f64, jd: f64, fraction: f64) -> Result<(f64, f64)> {
    let day = jd - J2000_JD_TDB;
    let nominal = (boundary / SECONDS_PER_DAY - day) - fraction;
    let mut radius = (16.0
        * f64::EPSILON
        * day
            .abs()
            .max(fraction.abs())
            .max((boundary / SECONDS_PER_DAY).abs()))
    .max(f64::MIN_POSITIVE);
    let (mut lower, mut upper) = (nominal - radius, nominal + radius);
    while emda_epoch_seconds(jd, fraction + lower)? >= boundary
        || emda_epoch_seconds(jd, fraction + upper)? < boundary
    {
        radius *= 2.0;
        lower = nominal - radius;
        upper = nominal + radius;
        if !lower.is_finite() || !upper.is_finite() {
            return Err(Error::MathDomain { value: boundary });
        }
    }
    loop {
        let middle = lower + (upper - lower) / 2.0;
        if middle == lower || middle == upper {
            break;
        }
        if emda_epoch_seconds(jd, fraction + middle)? < boundary {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    Ok((lower, upper))
}

/// Convert separate JD parts to separate seconds and their error-free sum.
/// Both products are finite-checked. The error-free trailing sum word retains
/// fractions near a closed descriptor endpoint.
pub(super) fn jd_parts_seconds(jd: f64, fraction: f64) -> Result<(f64, f64, f64, f64)> {
    for value in [jd, fraction] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let a = (jd - J2000_JD_TDB) * SECONDS_PER_DAY;
    let b = fraction * SECONDS_PER_DAY;
    let hi = a + b;
    if !a.is_finite() || !b.is_finite() || !hi.is_finite() {
        return Err(Error::NonFiniteArgument { value: hi });
    }
    let recovered_b = hi - a;
    let lo = (a - (hi - recovered_b)) + (b - recovered_b);
    Ok((a, b, hi, lo))
}

/// Type-21 epoch in ephemeris seconds past J2000 from split TDB JD parts.
/// Type21's recurrence consumes scalar seconds. Normalize in days before
/// scaling, while coverage selection retains both original input parts.
pub(super) fn emda_epoch_seconds(jd: f64, fraction: f64) -> Result<f64> {
    for value in [jd, fraction] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let et = ((jd - J2000_JD_TDB) + fraction) * SECONDS_PER_DAY;
    if !et.is_finite() {
        return Err(Error::NonFiniteArgument { value: et });
    }
    Ok(et)
}

/// Compare a normalized two-word seconds sum against closed scalar endpoints.
/// A nonzero residual distinguishes fractions neighboring an endpoint even
/// when their rounded leading sum word equals that endpoint.
pub(super) fn split_seconds_covers(hi: f64, lo: f64, start: f64, end: f64) -> bool {
    (hi > start || (hi == start && lo >= 0.0)) && (hi < end || (hi == end && lo <= 0.0))
}

/// Select the latest covering descriptor at an uncollapsed two-part JD.
/// Coverage errors retain the target and stored span; split fractions cannot
/// leak through a rounded ET endpoint into descriptor extrapolation.
pub(super) fn select_two_part_jd(
    index: &SpkIndex,
    target: i32,
    jd: f64,
    fraction: f64,
) -> Result<&Segment> {
    select_two_part_jd_with_domain(index, target, jd, fraction).map(|(segment, _)| segment)
}

/// Select once and retain the open descriptor-priority region at this epoch.
/// Boundaries of higher-priority descriptors cannot certify a local derivative.
pub(super) fn select_two_part_jd_with_domain(
    index: &SpkIndex,
    target: i32,
    jd: f64,
    fraction: f64,
) -> Result<(&Segment, bool)> {
    let (_, _, hi, lo) = jd_parts_seconds(jd, fraction)?;
    let mut priority_boundary = false;
    for segment in index
        .segments
        .iter()
        .rev()
        .filter(|segment| segment.descriptor.target == target)
    {
        let d = segment.descriptor;
        priority_boundary |= lo == 0.0 && (hi == d.start_et || hi == d.end_et);
        if split_seconds_covers(hi, lo, d.start_et, d.end_et) {
            let open = (hi > d.start_et || (hi == d.start_et && lo > 0.0))
                && (hi < d.end_et || (hi == d.end_et && lo < 0.0));
            return Ok((segment, open && !priority_boundary));
        }
    }
    let spans = index.coverage(target);
    match (spans.first(), spans.last()) {
        (Some(&(start, _)), Some(&(_, end))) => Err(Error::SpkOutOfRange {
            target,
            et: hi + lo,
            start,
            end,
        }),
        _ => Err(Error::SpkBodyNotFound { target }),
    }
}

/// Evaluate `f(s) = sum c[k] T_k(s)` and its analytic derivative `df/ds`.
///
/// Coefficients are in ascending degree order; `s` is dimensionless and may
/// be any finite number. Backward Clenshaw steps group the subtraction before
/// adding the coefficient: `w0 = c[k] + (2*s*w1 - w2)`. The final value is
/// `c[0] + (s*w0 - w1)`. Differentiating those states gives
/// `dw0 = 2*w1 + dw1*(2*s) - dw2`, and `df/ds = w0 + s*dw0 - dw1`.
/// The derivative states use the previous value state at each iteration, so
/// both recurrences can be streamed together without saving intermediate rows.
///
/// Empty coefficients and non-finite inputs are errors; finite arithmetic may
/// overflow. Constants return their value and an exactly zero derivative without
/// forming `2*s`. Units are the coefficient units and coefficient units per s;
/// the caller supplies segment-time scaling. No finite differences are used.
fn spk_chebyshev_value_and_derivative(coeffs: &[f64], s: f64) -> Result<(f64, f64)> {
    if coeffs.is_empty() {
        return Err(Error::EmptyCoefficients);
    }
    if !s.is_finite() {
        return Err(Error::NonFiniteArgument { value: s });
    }
    if let Some(&value) = coeffs.iter().find(|c| !c.is_finite()) {
        return Err(Error::NonFiniteArgument { value });
    }
    if coeffs.len() == 1 {
        return Ok((coeffs[0], 0.0));
    }
    let s2 = 2.0 * s;
    let (mut w0, mut w1) = (0.0, 0.0);
    let (mut dw0, mut dw1) = (0.0, 0.0);
    for &coefficient in coeffs[1..].iter().rev() {
        let w2 = w1;
        w1 = w0;
        w0 = s2 * w1 - w2 + coefficient;
        let dw2 = dw1;
        dw1 = dw0;
        dw0 = 2.0 * w1 + s2 * dw1 - dw2;
    }
    Ok((s * w0 - w1 + coeffs[0], w0 + s * dw0 - dw1))
}

impl EmdaSegment {
    /// Difference-table dimension of this segment.
    pub fn maxdim(&self) -> usize {
        self.maxdim
    }

    /// Number of difference-line records.
    pub fn n_records(&self) -> usize {
        self.n_records
    }
}

/// Index of the record serving `et`: the first record whose final epoch
/// is strictly greater than `et`, or the last record when none is.
///
/// Interior record boundaries belong to the following difference line.
/// The final descriptor endpoint uses the last record. Callers establish
/// closed descriptor coverage before the recurrence.
pub fn emda_record_index(segment: &EmdaSegment, et: f64) -> Result<usize> {
    if !et.is_finite() {
        return Err(Error::NonFiniteArgument { value: et });
    }
    Ok(segment
        .epochs
        .iter()
        .position(|&e| e > et)
        .unwrap_or(segment.n_records - 1))
}

/// Reference state and epoch of record `index`: `(position_km,
/// velocity_km_s, tl_seconds)`. Used by tests and diagnostics;
/// evaluation at `tl` reproduces this state exactly.
pub fn emda_record_reference(
    segment: &EmdaSegment,
    index: usize,
) -> Result<([f64; 3], [f64; 3], f64)> {
    let rec = segment
        .records
        .get(index * segment.dlsize..(index + 1) * segment.dlsize);
    let Some(rec) = rec else {
        return Err(Error::SpkCorrupt {
            what: "SPK type 21 record".to_string(),
            offset: segment.data_offset,
        });
    };
    let m = segment.maxdim;
    Ok((
        [rec[m + 1], rec[m + 3], rec[m + 5]],
        [rec[m + 2], rec[m + 4], rec[m + 6]],
        rec[0],
    ))
}

/// Evaluate a decoded type 21 segment at ephemeris seconds past J2000.
///
/// Krogh's difference-line recurrence over the record selected by
/// [`emda_record_index`], transcribed operation-for-operation from the
/// NAIF-specified record layout: stepsize ratios from `G`, reciprocal
/// weights `W`, position interpolation from the reference state plus
/// the difference-table correction, then velocity interpolation. Units
/// are km and km/s, like [`eval_chebyshev`]. A zero stepsize entry in
/// the used prefix of `G` makes the record unevaluable
/// ([`Error::SpkCorrupt`]); unused tail entries may be zero.
pub fn eval_emda(segment: &EmdaSegment, et: f64) -> Result<ChebState> {
    if !et.is_finite() {
        return Err(Error::NonFiniteArgument { value: et });
    }
    let index = emda_record_index(segment, et)?;
    eval_emda_reader_record(segment, et, index)
}

/// Evaluate an already selected immutable EMDA record, in km and km/s.
///
/// `index` is zero based and must identify one decoded record. Invalid
/// indices and nonfinite epochs are refused before recurrence arithmetic.
/// Coefficient/stepsize validation is shared by every entry point.
pub(super) fn eval_emda_reader_record(
    segment: &EmdaSegment,
    et: f64,
    index: usize,
) -> Result<ChebState> {
    if !et.is_finite() {
        return Err(Error::NonFiniteArgument { value: et });
    }
    let m = segment.maxdim;
    if index >= segment.n_records {
        return Err(Error::SpkCorrupt {
            what: "SPK type 21 record".to_owned(),
            offset: segment.data_offset,
        });
    }
    let base = index * segment.dlsize;
    let rec = &segment.records[base..base + segment.dlsize];
    let tl = rec[0];
    let refpos = [rec[m + 1], rec[m + 3], rec[m + 5]];
    let refvel = [rec[m + 2], rec[m + 4], rec[m + 6]];
    // Difference table in Fortran column-major order: DT[row][col].
    let dt = |row: usize, col: usize| -> f64 { rec[m + 7 + row + col * m] };
    let kqmax1 = rec[4 * m + 7] as usize;
    let kq = [
        rec[4 * m + 8] as usize,
        rec[4 * m + 9] as usize,
        rec[4 * m + 10] as usize,
    ];

    let delta = et - tl;
    let mut tp = delta;
    let mq2 = kqmax1 - 2;
    let mut fc = vec![0.0; kqmax1];
    let mut wc = vec![0.0; kqmax1 - 1];
    fc[0] = 1.0;
    for j in 1..=mq2 {
        let g = rec[j];
        if g == 0.0 {
            return Err(Error::SpkCorrupt {
                what: "SPK type 21 stepsize".to_string(),
                offset: segment.data_offset + (base + j) * 8,
            });
        }
        fc[j] = tp / g;
        wc[j - 1] = delta / g;
        tp = delta + g;
    }
    let mut w: Vec<f64> = (1..=kqmax1).map(|i| 1.0 / i as f64).collect();
    let mut ks = kqmax1 - 1;
    let mut ks1 = ks - 1;
    let mut jx = 0usize;
    while ks >= 2 {
        jx += 1;
        for j in 1..=jx {
            w[j + ks - 1] = fc[j] * w[j + ks1 - 1] - wc[j - 1] * w[j + ks - 1];
        }
        ks = ks1;
        ks1 -= 1;
    }
    let mut position_km = [0.0; 3];
    for (i, slot) in position_km.iter_mut().enumerate() {
        let mut sum = 0.0;
        for j in (1..=kq[i]).rev() {
            sum += dt(j - 1, i) * w[j + ks - 1];
        }
        *slot = refpos[i] + delta * (refvel[i] + delta * sum);
    }
    for j in 1..=jx {
        w[j + ks - 1] = fc[j] * w[j + ks1 - 1] - wc[j - 1] * w[j + ks - 1];
    }
    ks -= 1;
    let mut velocity_km_s = [0.0; 3];
    for (i, slot) in velocity_km_s.iter_mut().enumerate() {
        let mut sum = 0.0;
        for j in (1..=kq[i]).rev() {
            sum += dt(j - 1, i) * w[j + ks - 1];
        }
        *slot = refvel[i] + delta * sum;
    }
    Ok(ChebState {
        position_km,
        velocity_km_s,
    })
}

/// State of `target` relative to its segment-declared center, with the
/// serving descriptor.
///
/// Selects the latest-wins covering segment (closed coverage), gates
/// its frame, and dispatches on its data type to the type-2/3 and
/// type-21 evaluators. A coverage miss is [`Error::SpkOutOfRange`] bounded by
/// the merged stored span (never a silent extrapolation); an unknown
/// target is [`Error::SpkBodyNotFound`].
pub fn eval_relative(
    bytes: &[u8],
    endian: DafEndian,
    index: &SpkIndex,
    target: i32,
    et: f64,
) -> Result<(ChebState, Descriptor)> {
    let Some(segment) = index.select(target, et)? else {
        let spans = index.coverage(target);
        let (Some(&(first, _)), Some(&(_, last))) = (spans.first(), spans.last()) else {
            return Err(Error::SpkBodyNotFound { target });
        };
        return Err(Error::SpkOutOfRange {
            target,
            et,
            start: first,
            end: last,
        });
    };
    if segment.descriptor.frame != FRAME_J2000 {
        return Err(Error::SpkUnsupportedFrame {
            frame: segment.descriptor.frame,
        });
    }
    let state = match segment.descriptor.seg_type {
        TYPE_CHEBYSHEV_POS | TYPE_CHEBYSHEV_POS_VEL => {
            let decoded = decode_chebyshev(bytes, segment, endian)?;
            eval_chebyshev(&decoded, et)?
        }
        TYPE_EMDA => {
            let decoded = decode_emda(bytes, segment, endian)?;
            eval_emda(&decoded, et)?
        }
        other => {
            return Err(Error::SpkUnsupportedType { seg_type: other });
        }
    };
    Ok((state, segment.descriptor))
}

/// SSB-centered ICRS/J2000 state of `target` in km and km/s.
///
/// The declared center chain consists of
/// the kernel-relative state from [`eval_relative`] plus the
/// barycentric state of the declared center. Center 0 contributes the
/// identically-zero state without consulting `provider`; any other
/// center resolves through `provider` (an unavailable center surfaces
/// the provider's typed error, never invented data).
///
/// Relative state is added to the physically resolved center state.
pub fn eval_ssb(
    bytes: &[u8],
    endian: DafEndian,
    index: &SpkIndex,
    target: i32,
    provider: &dyn CenterProvider,
    et: f64,
) -> Result<ChebState> {
    // Full chain gate first:
    // a mixed-center or non-J2000 kernel fails here, before any
    // segment is decoded or served. Registration  runs
    // `resolve_chain` once for an early failure; evaluation re-checks
    // standalone so the composed state is always single-center J2000.
    let chain = resolve_chain(index, target)?;
    let (relative, descriptor) = eval_relative(bytes, endian, index, target, et)?;
    debug_assert_eq!(descriptor.center, chain.center);
    debug_assert_eq!(descriptor.frame, chain.frame);
    if descriptor.center == CENTER_SSB {
        return Ok(relative);
    }
    let center = provider.center_state(descriptor.center, et)?;
    let mut position_km = [0.0; 3];
    let mut velocity_km_s = [0.0; 3];
    for i in 0..3 {
        position_km[i] = center.position_km[i] + relative.position_km[i];
        velocity_km_s[i] = center.velocity_km_s[i] + relative.velocity_km_s[i];
    }
    Ok(ChebState {
        position_km,
        velocity_km_s,
    })
}

/// Center-relative state at the uncollapsed TDB JD `jd + fraction`.
/// Type 2/3 records preserve both epoch parts; type 21 retains its scalar-ET
/// recurrence after strict split coverage selection. Output is km and km/s.
/// Byte-slice calls decode on each evaluation and do not retain caller memory.
pub fn eval_relative_two_part_jd(
    bytes: &[u8],
    endian: DafEndian,
    index: &SpkIndex,
    target: i32,
    jd: f64,
    fraction: f64,
) -> Result<(ChebState, Descriptor)> {
    let segment = select_two_part_jd(index, target, jd, fraction)?;
    if segment.descriptor.frame != FRAME_J2000 {
        return Err(Error::SpkUnsupportedFrame {
            frame: segment.descriptor.frame,
        });
    }
    let state = match segment.descriptor.seg_type {
        TYPE_CHEBYSHEV_POS | TYPE_CHEBYSHEV_POS_VEL => {
            eval_chebyshev_two_part_jd(&decode_chebyshev(bytes, segment, endian)?, jd, fraction)?
        }
        TYPE_EMDA => eval_emda(
            &decode_emda(bytes, segment, endian)?,
            emda_epoch_seconds(jd, fraction)?,
        )?,
        seg_type => return Err(Error::SpkUnsupportedType { seg_type }),
    };
    Ok((state, segment.descriptor))
}

/// SSB-relative state at an uncollapsed two-part TDB JD, in km and km/s.
/// Relative type 2/3 states preserve the fraction. External centers use
/// [`CenterProvider::center_state_two_part_jd`]; the native registry-backed
/// provider retains both parts recursively. Legacy scalar-only providers
/// retain their documented rounded-seconds default.
pub fn eval_ssb_two_part_jd(
    bytes: &[u8],
    endian: DafEndian,
    index: &SpkIndex,
    target: i32,
    provider: &dyn CenterProvider,
    jd: f64,
    fraction: f64,
) -> Result<ChebState> {
    let chain = resolve_chain(index, target)?;
    let (relative, _) = eval_relative_two_part_jd(bytes, endian, index, target, jd, fraction)?;
    if chain.center == CENTER_SSB {
        return Ok(relative);
    }
    let center = provider.center_state_two_part_jd(chain.center, jd, fraction)?;
    Ok(ChebState {
        position_km: std::array::from_fn(|i| center.position_km[i] + relative.position_km[i]),
        velocity_km_s: std::array::from_fn(|i| center.velocity_km_s[i] + relative.velocity_km_s[i]),
    })
}

#[cfg(test)]
mod primitive_tests {
    use super::*;

    fn synthetic_position_segment(init: f64, interval: f64, coeffs: &[f64]) -> ChebSegment {
        let mut words = vec![init + interval / 2.0, interval / 2.0];
        for _ in 0..3 {
            words.extend_from_slice(coeffs);
        }
        ChebSegment {
            seg_type: TYPE_CHEBYSHEV_POS,
            init_et: init,
            int_len: interval,
            record_words: words.len(),
            n_records: 1,
            words,
            data_offset: 0,
        }
    }

    #[test]
    fn chebyshev_record_support_preserves_split_join_ownership_for_both_types() {
        for seg_type in [TYPE_CHEBYSHEV_POS, TYPE_CHEBYSHEV_POS_VEL] {
            let components = if seg_type == TYPE_CHEBYSHEV_POS { 3 } else { 6 };
            let mut words = Vec::new();
            for midpoint in [4.0, 12.0] {
                words.extend([midpoint, 4.0]);
                words.extend(std::iter::repeat_n(1.0, components));
            }
            let segment = ChebSegment {
                seg_type,
                init_et: 0.0,
                int_len: 8.0,
                record_words: 2 + components,
                n_records: 2,
                words,
                data_offset: 0,
            };
            let at_join = 8.0 / SECONDS_PER_DAY;
            for fraction in [at_join.next_down(), at_join, at_join.next_up()] {
                let bounds = chebyshev_record_offsets(&segment, J2000_JD_TDB, fraction).unwrap();
                assert!(bounds.0 <= 0.0 && bounds.1 >= 0.0);
                let before = fraction < at_join;
                assert_eq!(bounds.1 < 1.0 / SECONDS_PER_DAY, before);
                assert_eq!(bounds.0 > -1.0 / SECONDS_PER_DAY, !before);
                assert!(eval_chebyshev_two_part_jd(&segment, J2000_JD_TDB, fraction).is_ok());
            }
        }
    }

    #[test]
    fn type21_record_support_matches_scalar_selector_at_sub_ulp_neighbors() {
        let boundary = 10_000_000.0;
        let segment = EmdaSegment {
            maxdim: 1,
            n_records: 2,
            dlsize: 15,
            records: vec![0.0; 30],
            epochs: vec![boundary, boundary + 16.0],
            data_offset: 0,
        };
        let descriptor = Descriptor {
            start_et: boundary - 16.0,
            end_et: boundary + 16.0,
            target: 2000001,
            center: 0,
            frame: FRAME_J2000,
            seg_type: TYPE_EMDA,
            initial_address: 1,
            final_address: 32,
        };
        let fraction = boundary / SECONDS_PER_DAY;
        for central in [fraction.next_down(), fraction, fraction.next_up()] {
            let bounds = emda_record_offsets(&segment, descriptor, J2000_JD_TDB, central).unwrap();
            assert!(
                bounds.0 <= 0.0 && bounds.1 >= 0.0,
                "{central:?}: {bounds:?}"
            );
            let index =
                emda_record_index(&segment, emda_epoch_seconds(J2000_JD_TDB, central).unwrap())
                    .unwrap();
            for offset in [bounds.0, bounds.1] {
                let sampled = emda_epoch_seconds(J2000_JD_TDB, central + offset).unwrap();
                assert_eq!(emda_record_index(&segment, sampled).unwrap(), index);
            }
        }
    }

    #[test]
    fn record_midpoint_defines_the_polynomial_origin() {
        let init = 10_000_000_000_000_000.0;
        let segment = synthetic_position_segment(init, 3.0, &[7.0, 2.0]);
        // The stored midpoint is init+2 after binary64 encoding. Its actual
        // polynomial coordinate is exactly zero, independently of directory.
        let state = eval_chebyshev(&segment, segment.words[0]).unwrap();
        assert_eq!(state.position_km, [7.0; 3]);
        assert_eq!(state.velocity_km_s, [2.0 / 1.5; 3]);
    }

    #[test]
    fn scalar_and_split_directories_refuse_extrapolation_and_nonfinite_epochs() {
        let model = synthetic_position_segment(0.0, 8.0, &[1.0, 2.0]);
        for et in [-1.0, 8.001, 16.0] {
            assert!(matches!(
                eval_chebyshev(&model, et),
                Err(Error::SpkOutOfRange { .. })
            ));
            assert!(matches!(
                eval_chebyshev_two_part_jd(&model, J2000_JD_TDB, et / SECONDS_PER_DAY),
                Err(Error::SpkOutOfRange { .. })
            ));
        }
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, f64::MAX] {
            assert!(eval_chebyshev_two_part_jd(&model, bad, 0.0).is_err());
            assert!(eval_chebyshev_two_part_jd(&model, J2000_JD_TDB, bad).is_err());
            assert!(emda_epoch_seconds(bad, 0.0).is_err());
        }
        assert_eq!(eval_chebyshev(&model, 0.0).unwrap().position_km, [-1.0; 3]);
        assert_eq!(eval_chebyshev(&model, 8.0).unwrap().position_km, [3.0; 3]);
    }

    #[test]
    fn spk_differentiated_clenshaw_matches_exact_quadratic() {
        for (s, value, derivative) in [
            (-1.0, 9.0, -17.0),
            (0.0, 2.0, 3.0),
            (0.25, 3.375, 8.0),
            (1.0, 15.0, 23.0),
        ] {
            assert_eq!(
                spk_chebyshev_value_and_derivative(&[7.0, 3.0, 5.0], s).unwrap(),
                (value, derivative)
            );
        }
        assert_eq!(
            spk_chebyshev_value_and_derivative(&[3.0], f64::MAX).unwrap(),
            (3.0, 0.0)
        );
        assert!(spk_chebyshev_value_and_derivative(&[], 0.0).is_err());
        assert!(spk_chebyshev_value_and_derivative(&[f64::NAN], 0.0).is_err());
        assert!(spk_chebyshev_value_and_derivative(&[1.0], f64::INFINITY).is_err());
    }
}
