//! Evaluation implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

/// Evaluate one segment's coefficient words at normalized epoch `tau`.
///
/// Shared by the LEB1 raw path and the LEB2 decoded-blob path (the body-evaluation layer):
/// `seg_coeffs` holds exactly one segment (`components` blocks of
/// `degree + 1` f64). Clenshaw evaluates values and first derivatives
/// with constant auxiliary space; multiplication by `2 / interval_days`
/// converts normalized-coordinate derivatives to rates per day.
/// Ecliptic longitude wraps via `rem_euclid(360.0)`.
///
/// Units follow the stored frame: AU xyz and AU/day for ICRS bodies,
/// deg/deg/AU and (deg/day, deg/day, AU/day) for ecliptic bodies.
pub(super) fn eval_one_segment(
    entry: &BodyEntry,
    seg_coeffs: &[f64],
    tau: f64,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let width = segment_component_width(entry, seg_coeffs)?;
    let mut pos = vec![0.0; entry.components() as usize];
    let mut vel = vec![0.0; entry.components() as usize];
    eval_segment_into(entry, seg_coeffs, tau, width, &mut pos, &mut vel)?;
    Ok((pos, vel))
}

/// Three-component form of the ordinary stored-channel evaluator. Selection,
/// coefficient validation, Clenshaw arithmetic and angular normalization are
/// shared with the public vector route; only the output storage is inline.
pub(super) fn eval_one_segment_inline(
    entry: &BodyEntry,
    seg_coeffs: &[f64],
    tau: f64,
) -> Result<CartesianState> {
    debug_assert_eq!(entry.components(), 3);
    let width = segment_component_width(entry, seg_coeffs)?;
    let mut pos = [0.0; 3];
    let mut vel = [0.0; 3];
    eval_segment_into(entry, seg_coeffs, tau, width, &mut pos, &mut vel)?;
    Ok((pos, vel))
}

/// Validate before allocating output buffers, including malformed shape/counts.
fn segment_component_width(entry: &BodyEntry, coefficients: &[f64]) -> Result<usize> {
    let invalid = || Error::LebCorrupt {
        what: format!("LEB coefficient blob (body {})", entry.body_id()),
        offset: usize::try_from(entry.data_offset()).unwrap_or(usize::MAX),
    };
    let width = (entry.degree() as usize)
        .checked_add(1)
        .ok_or_else(invalid)?;
    let count = (entry.components() as usize)
        .checked_mul(width)
        .ok_or_else(invalid)?;
    if coefficients.len() != count {
        return Err(invalid());
    }
    Ok(width)
}

/// Fill caller-owned output slices with the existing per-component operation order.
fn eval_segment_into(
    entry: &BodyEntry,
    coefficients: &[f64],
    tau: f64,
    width: usize,
    pos: &mut [f64],
    vel: &mut [f64],
) -> Result<()> {
    debug_assert_eq!(pos.len(), entry.components() as usize);
    debug_assert_eq!(vel.len(), pos.len());
    let scale = 2.0 / entry.interval_days();
    for component in 0..pos.len() {
        let row = &coefficients[component * width..(component + 1) * width];
        let (value, derivative) = crate::math::chebyshev::chebyshev_value_and_derivative(row, tau)?;
        pos[component] = value;
        vel[component] = derivative * scale;
    }
    if is_ecliptic_coord(entry.coord_type()) {
        pos[0] = pos[0].rem_euclid(360.0);
    }
    Ok(())
}

/// Shared multi-epoch batch tail: one scalar evaluation per epoch, in
/// input order, failing fast with the first element's error (the
/// release `leb_vector._at` loop contract). No segment, chunk, or
/// coefficient path of its own — every element runs the caller's
/// scalar evaluator.
pub(super) fn batch_eval(
    jds: &[f64],
    eval: impl Fn(f64) -> Result<(Vec<f64>, Vec<f64>)>,
) -> Result<Vec<(Vec<f64>, Vec<f64>)>> {
    let mut out = Vec::with_capacity(jds.len());
    for &jd in jds {
        out.push(eval(jd)?);
    }
    Ok(out)
}

/// Split-epoch form of [`batch_eval`].
pub(super) fn batch_eval_epochs(
    epochs: &[TwoPartEpoch],
    eval: impl Fn(TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)>,
) -> Result<Vec<(Vec<f64>, Vec<f64>)>> {
    let mut out = Vec::with_capacity(epochs.len());
    for &epoch in epochs {
        out.push(eval(epoch)?);
    }
    Ok(out)
}

/// Value-only Clenshaw evaluation of one nutation segment
/// (`(dpsi, deps)` channels, radians, no velocity), the release
/// `eval_nutation` tail shared by both formats.
pub(super) fn eval_nutation_words(
    header: &NutationHeader,
    seg_coeffs: &[f64],
    tau: f64,
) -> Result<(f64, f64)> {
    const WHAT: &str = "LEB nutation table";
    let deg1 = (header.degree() as usize)
        .checked_add(1)
        .ok_or_else(|| Error::LebCorrupt {
            what: WHAT.to_string(),
            offset: 0,
        })?;
    if header.components() != 2 || seg_coeffs.len() != 2 * deg1 {
        return Err(Error::LebCorrupt {
            what: WHAT.to_string(),
            offset: 0,
        });
    }
    let dpsi = chebyshev_value(&seg_coeffs[0..deg1], tau)?;
    let deps = chebyshev_value(&seg_coeffs[deg1..2 * deg1], tau)?;
    Ok((dpsi, deps))
}

/// Look up Delta-T (TT minus UT1), in days, at Julian Day `jd`.
///
/// The equally sized `jds` and `vals` slices contain a sparse table of
/// `(Julian Day, Delta-T in days)` pairs, with epochs in ascending order.
/// `jd` must use the same epoch-coordinate convention as `jds`; this function
/// has no time-scale metadata and performs no TT/UT1/UTC conversion. Sparse
/// spacing alone does not guarantee interchangeable clocks or physical accuracy.
/// Endpoint values clamp requests outside the table; interior requests use
/// right-biased binary search followed by linear interpolation. An adjacent
/// zero-width interval returns its left value. A finite request against a
/// singleton table returns its sole value.
///
/// An empty table, unequal slice lengths, or a singleton lacking an applicable
/// endpoint is [`Error::LebCorrupt`]. A non-finite requested epoch is
/// [`Error::NonFiniteArgument`]. These guards do not validate or alter the
/// stored table words; in particular, an interpolated non-finite value retains
/// its ordinary floating-point result. Empty-table refusal takes priority over
/// the other guards.
pub fn eval_delta_t(jds: &[f64], vals: &[f64], jd: f64) -> Result<f64> {
    if jds.is_empty() {
        return Err(Error::LebCorrupt {
            what: "LEB Delta-T table (absent)".to_string(),
            offset: 0,
        });
    }
    if jds.len() != vals.len() {
        return Err(Error::LebCorrupt {
            what: "LEB Delta-T table (unequal lengths)".to_string(),
            offset: 0,
        });
    }
    eval_delta_t_rows(jds.len(), jd, |index| Ok((jds[index], vals[index])))
}

/// One interpolation rule for contiguous tables and bounded positional rows.
/// Input records have the same sorted-epoch contract as `eval_delta_t`.
pub(super) fn eval_delta_t_rows(
    count: usize,
    jd: f64,
    row: impl Fn(usize) -> Result<(f64, f64)>,
) -> Result<f64> {
    if count == 0 {
        return Err(Error::LebCorrupt {
            what: "LEB Delta-T table (absent)".into(),
            offset: 0,
        });
    }
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    let first = row(0)?;
    if jd <= first.0 {
        return Ok(first.1);
    }
    let last = row(count - 1)?;
    if jd >= last.0 {
        return Ok(last.1);
    }
    if count < 2 {
        return Err(Error::LebCorrupt {
            what: "LEB Delta-T table (no interpolation interval)".into(),
            offset: 0,
        });
    }
    let (mut lo, mut hi) = (0, count);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if row(mid)?.0 <= jd {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    let index = lo.saturating_sub(1).min(count - 2);
    let left = row(index)?;
    let right = row(index + 1)?;
    let span = right.0 - left.0;
    if span == 0.0 {
        return Ok(left.1);
    }
    // Normalize the offset before blending, retaining the established order.
    let t = (jd - left.0) / span;
    Ok(left.1 + t * (right.1 - left.1))
}

/// Evaluate the requested body center and its barycentric state from an
/// already-evaluated stored LEB state.
///
/// Non-system-barycenter bodies
/// return the stored state untouched. For system barycenters,
/// `offset_at(jd)` supplies the planet-center-minus-barycenter offset in
/// AU (`None` when no segment covers the epoch): the position becomes
/// stored + offset, and the velocity gains the centered-difference rate
/// of the same offset with step
/// [`CENTER_OFFSET_VELOCITY_STEP_DAYS`], in the exact release operation
/// order (`stored_vel + (after - before) * (1 / (2h))`). Side epochs
/// without an offset contribute zero — mirroring the release, whose
/// side evaluations fall back to the passed-in zero vector on segment
/// failure. When the center epoch has no offset, or when the offset
/// position is bitwise identical to the stored one, the
/// stored state is returned unchanged: retaining the system barycenter
/// is the explicit data-availability fallback, never a silent
/// substitution. A non-finite `jd_tt` skips the correction.
pub fn eval_body_center(
    stored_pos: [f64; 3],
    stored_vel: [f64; 3],
    is_system_bary: bool,
    jd_tt: f64,
    offset_at: impl Fn(f64) -> Option<[f64; 3]>,
) -> ([f64; 3], [f64; 3]) {
    if !is_system_bary || !jd_tt.is_finite() {
        return (stored_pos, stored_vel);
    }
    let Some(offset_now) = offset_at(jd_tt).and_then(finite_offset) else {
        return (stored_pos, stored_vel);
    };
    let center_pos = apply_center_offset(stored_pos, offset_now);
    if center_pos == stored_pos {
        return (stored_pos, stored_vel);
    }
    let h = CENTER_OFFSET_VELOCITY_STEP_DAYS;
    let before = offset_at(jd_tt - h)
        .and_then(finite_offset)
        .unwrap_or([0.0; 3]);
    let after = offset_at(jd_tt + h)
        .and_then(finite_offset)
        .unwrap_or([0.0; 3]);
    let scale = 1.0 / (2.0 * h);
    let center_vel = [
        stored_vel[0] + (after[0] - before[0]) * scale,
        stored_vel[1] + (after[1] - before[1]) * scale,
        stored_vel[2] + (after[2] - before[2]) * scale,
    ];
    (center_pos, center_vel)
}

/// Exact inclusive TT window of an interpolated-apse model, derived from
/// the generated residual-grid metadata (first grid JD through the JD of
/// its last sample). `None` for any other body or when the grid metadata
/// is unavailable. Source: release `lunar._interpolated_apse_model_window`.
pub fn interpolated_apse_model_window(body_id: i32) -> Option<(f64, f64)> {
    match body_id {
        21 => Some((APOGEE_MODEL_JD_START, APOGEE_MODEL_JD_END)),
        22 => Some((PERIGEE_MODEL_JD_START, PERIGEE_MODEL_JD_END)),
        _ => None,
    }
}

/// Evaluate one probe epoch through the container's own flavor path:
/// monolithic coefficients for LEB1, decoded blobs for LEB2 v1/v2.
pub(super) fn eval_verify_probe(
    index: &FileIndex,
    bytes: &[u8],
    body_id: i32,
    jd: f64,
) -> Result<(Vec<f64>, Vec<f64>)> {
    if index.header().format() == LebFormat::Leb1 {
        index.eval_body(bytes, body_id, jd)
    } else {
        index.eval_leb2_body(bytes, body_id, jd)
    }
}
