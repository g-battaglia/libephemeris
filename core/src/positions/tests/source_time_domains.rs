// SPDX-License-Identifier: AGPL-3.0
//! Focused tests for the parent position-routing module.

use super::*;
use std::cell::Cell;

#[test]
fn source_observation_keeps_ten_nan_iterations_and_direct_at_distinction() {
    let reads = Cell::new(0);
    let splits = Cell::new(0);
    let sample = |day: f64, fraction: f64| {
        assert!(day.is_nan() && fraction.is_nan());
        reads.set(reads.get() + 1);
        Ok([f64::NAN; 3])
    };
    let split = |_: f64, _: f64, _: f64| {
        splits.set(splits.get() + 1);
        Ok((f64::NAN, f64::NAN))
    };
    let result = retarded_earth_minus_sun_samples_with_split(
        f64::NAN,
        f64::NAN,
        [f64::NAN; 3],
        [f64::NAN; 3],
        &sample,
        &split,
    );
    assert!(matches!(
        result,
        Err(Error::ConvergenceFailed {
            algorithm: "spk_light_time",
            iterations: 10,
            ..
        })
    ));
    assert_eq!(reads.get(), 10);
    assert_eq!(splits.get(), 10);
    let first_error =
        |_: f64, _: f64| -> Result<[f64; 3]> { Err(Error::SpkBodyNotFound { target: 399 }) };
    assert!(matches!(
        retarded_earth_minus_sun_samples_with_split(
            f64::NAN,
            f64::NAN,
            [f64::NAN; 3],
            [f64::NAN; 3],
            &first_error,
            &split
        ),
        Err(Error::SpkBodyNotFound { target: 399 })
    ));
}
