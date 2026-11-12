//! Public compatibility contract for explicit legacy boundary calibration.

use libephemeris_core::errors::Error;
#[allow(deprecated)]
use libephemeris_core::time::legacy_sidereal;
use libephemeris_core::time::{
    SIDEREAL_J2000_JD, SIDEREAL_MODERN_END_JD, SIDEREAL_MODERN_START_JD, SiderealDeltaT,
    mean_sidereal_time_deg, sidereal_continuity_offsets,
};

fn angular_difference(a: f64, b: f64) -> f64 {
    (a - b + 180.0).rem_euclid(360.0) - 180.0
}

#[test]
#[allow(deprecated)] // This test deliberately exercises the explicit legacy contract.
fn supplied_offsets_have_their_documented_branch_effects() {
    let base = SiderealDeltaT::try_new(0.0008, 0.125, -0.25).unwrap();
    let lower = SiderealDeltaT::try_new(0.0008, 0.375, -0.25).unwrap();
    let upper = SiderealDeltaT::try_new(0.0008, 0.125, 0.25).unwrap();
    assert_eq!(base.offset_lo_deg(), 0.125);
    assert_eq!(base.offset_hi_deg(), -0.25);
    assert_eq!(lower.offset_lo_deg(), 0.375);
    assert_eq!(upper.offset_hi_deg(), 0.25);
    for (jd, expected_lower, expected_upper) in [
        (SIDEREAL_MODERN_START_JD.next_down(), -0.25, 0.0),
        (SIDEREAL_MODERN_START_JD, -0.25, 0.0),
        (SIDEREAL_MODERN_START_JD.next_up(), 0.0, 0.0),
        (SIDEREAL_J2000_JD, 0.0, 0.0),
        (SIDEREAL_MODERN_END_JD.next_down(), 0.0, 0.0),
        (SIDEREAL_MODERN_END_JD, 0.0, -0.5),
        (SIDEREAL_MODERN_END_JD.next_up(), 0.0, -0.5),
    ] {
        let original = legacy_sidereal::mean_sidereal_time_deg(jd, &base).unwrap();
        for (changed, expected) in [(&lower, expected_lower), (&upper, expected_upper)] {
            let value = legacy_sidereal::mean_sidereal_time_deg(jd, changed).unwrap();
            // These dyadic degree corrections are exact at the returned scale.
            assert_eq!(angular_difference(value, original), expected, "JD={jd}");
        }
    }
}

#[test]
#[allow(deprecated)] // Compatibility and explicit legacy entry points must agree.
fn legacy_offsets_retain_bits_and_canonical_time_uses_only_its_clock() {
    let offsets = legacy_sidereal::sidereal_continuity_offsets(0.001, 0.002).unwrap();
    let compatibility_offsets = sidereal_continuity_offsets(0.001, 0.002).unwrap();
    assert_ne!(offsets.0, 0.0);
    assert_ne!(offsets.1, 0.0);
    assert_eq!(offsets.0.to_bits(), compatibility_offsets.0.to_bits());
    assert_eq!(offsets.1.to_bits(), compatibility_offsets.1.to_bits());
    let dt = SiderealDeltaT::try_new(0.0008, offsets.0, offsets.1).unwrap();
    for jd in [
        SIDEREAL_MODERN_START_JD.next_down(),
        SIDEREAL_MODERN_START_JD,
        SIDEREAL_MODERN_START_JD.next_up(),
        SIDEREAL_J2000_JD,
        SIDEREAL_MODERN_END_JD.next_down(),
        SIDEREAL_MODERN_END_JD,
        SIDEREAL_MODERN_END_JD.next_up(),
        SIDEREAL_J2000_JD - 2001.0 * 36525.0,
        SIDEREAL_J2000_JD + 2001.0 * 36525.0,
        SIDEREAL_J2000_JD + 10000.0 * 36525.0,
    ] {
        let explicit = legacy_sidereal::mean_sidereal_time_deg(jd, &dt).unwrap();
        let compatibility = mean_sidereal_time_deg(jd, &dt).unwrap();
        assert!(explicit.is_finite(), "JD={jd}");
        assert!(compatibility.is_finite(), "JD={jd}");
        let different_offsets = SiderealDeltaT::try_new(0.0008, 123.0, -456.0).unwrap();
        assert_eq!(
            compatibility.to_bits(),
            mean_sidereal_time_deg(jd, &different_offsets)
                .unwrap()
                .to_bits(),
            "canonical clock must not consume legacy calibration at JD={jd}"
        );
    }
}

#[test]
#[allow(deprecated)] // Legacy typed errors remain part of the supported utility.
fn explicit_legacy_entries_preserve_nonfinite_error_order_and_payloads() {
    let dt = SiderealDeltaT::try_new(0.001, -0.0, 0.0).unwrap();
    assert_eq!(dt.offset_lo_deg().to_bits(), (-0.0_f64).to_bits());
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for result in [
            legacy_sidereal::mean_sidereal_time_deg(bad, &dt),
            legacy_sidereal::sidereal_continuity_offset_deg(bad, 0.001),
            legacy_sidereal::sidereal_continuity_offset_deg(SIDEREAL_J2000_JD, bad),
        ] {
            match result {
                Err(Error::NonFiniteArgument { value }) => {
                    assert_eq!(value.to_bits(), bad.to_bits())
                }
                result => panic!("unexpected error class: {result:?}"),
            }
        }
    }
    match legacy_sidereal::sidereal_continuity_offsets(f64::INFINITY, f64::NEG_INFINITY) {
        Err(Error::NonFiniteArgument { value }) => assert_eq!(value, f64::INFINITY),
        result => panic!("lower boundary must fail first: {result:?}"),
    }
}
