// SPDX-License-Identifier: AGPL-3.0
//! Smooth polynomial support inside actual binary64 record ownership.

use super::{Error, Result, select_segment};
use std::cmp::Ordering;

/// Find the maximal closed TT interval owned by the selected polynomial.
///
/// The caller has established `start <= jd <= end`. Coordinate evaluation and
/// normalized polynomial time use the same selector as these bounds. Nominal
/// starts are search hints; the predicate also accounts for clipped coverage,
/// coincident rounded starts and the selector's finite-quotient domain.
pub(super) fn smooth_segment_interval(
    start: f64,
    end: f64,
    interval: f64,
    count: u32,
    jd: f64,
) -> Result<(f64, f64)> {
    let selected = |epoch| select_segment(start, interval, count, epoch).map(|value| value.0);
    let record = selected(jd)?;
    let relative_record = |epoch| match selected(epoch) {
        Ok(index) => Ok(index.cmp(&record)),
        // A distant finite coverage endpoint can overflow the subtraction or
        // quotient even though the central epoch is valid. It lies beyond all
        // selectable records, so bracket the last valid word before it. Actual
        // central errors and all other metadata errors still propagate.
        Err(Error::NonFiniteArgument { value }) if value == f64::INFINITY && epoch > jd => {
            Ok(Ordering::Greater)
        }
        Err(error) => Err(error),
    };
    let lower = if record == 0 {
        start
    } else {
        first_selected_epoch(start, jd, start + record as f64 * interval, &|epoch| {
            relative_record(epoch).map(|order| order != Ordering::Less)
        })?
    };
    let upper = if relative_record(end)? == Ordering::Equal {
        end
    } else {
        let first_next =
            first_selected_epoch(jd, end, start + (record + 1) as f64 * interval, &|epoch| {
                relative_record(epoch).map(|order| order == Ordering::Greater)
            })?;
        first_next.next_down()
    };
    Ok((lower, upper))
}

/// Locate the first true word of a monotone selector predicate.
///
/// `upper` is known true. Nearby nominal words resolve ordinary astronomical
/// grids cheaply. Ordered-bit bisection also handles cancellation near zero
/// and intervals shorter than an epoch ulp, without an unbounded ulp walk.
fn first_selected_epoch(
    lower: f64,
    upper: f64,
    hint: f64,
    selected: &dyn Fn(f64) -> Result<bool>,
) -> Result<f64> {
    let mut probe = hint.clamp(lower, upper);
    for _ in 0..4 {
        if selected(probe)? {
            let previous = probe.next_down();
            if previous < lower || !selected(previous)? {
                return Ok(probe);
            }
            probe = previous;
        } else {
            probe = probe.next_up();
        }
    }
    if selected(lower)? {
        return Ok(lower);
    }
    let mut left = ordered_key(lower);
    let mut right = ordered_key(upper);
    // At most 64 halvings: these are finite binary64 words, not a tolerance.
    while right - left > 1 {
        let middle = left + (right - left) / 2;
        if selected(from_ordered_key(middle))? {
            right = middle;
        } else {
            left = middle;
        }
    }
    Ok(from_ordered_key(right))
}

fn ordered_key(value: f64) -> u64 {
    let bits = value.to_bits();
    if value.is_sign_negative() {
        !bits
    } else {
        bits ^ (1 << 63)
    }
}

fn from_ordered_key(key: u64) -> f64 {
    f64::from_bits(if key & (1 << 63) == 0 {
        !key
    } else {
        key ^ (1 << 63)
    })
}
