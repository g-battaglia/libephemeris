//! UT1-to-TT clock sanity bounds for modern epochs.
//!
//! The input is a UT1 Julian day and Delta-T is TT minus UT1 in seconds.
//! These broad physical bounds check the default clock near J2000 and 2020;
//! precise matched-data clock comparisons are covered separately.

use libephemeris_core::time::default_timescale_from_ut1;

#[test]
fn body08_original_j2000_skyfield_delta_t_assertions() {
    let (_tt, delta_t_seconds) = default_timescale_from_ut1(2451545.0).unwrap();
    assert!(60.0 < delta_t_seconds && delta_t_seconds < 70.0);
}

#[test]
fn body08_original_modern_skyfield_delta_t_assertions() {
    let (_tt, delta_t_seconds) = default_timescale_from_ut1(2458849.5).unwrap();
    assert!(65.0 < delta_t_seconds && delta_t_seconds < 75.0);
}
