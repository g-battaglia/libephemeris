//! Exact polynomial-input admission, independently of coefficient evaluation.
use super::*;

fn body_entry(start: f64, interval: f64, count: u32) -> BodyEntry {
    BodyEntry {
        body_id: 0,
        coord_type: COORD_ICRS_BARY,
        segment_count: count,
        jd_start: start,
        jd_end: start + interval * f64::from(count),
        interval_days: interval,
        degree: 3,
        components: 3,
        data_offset: 0,
        compressed_size: None,
        uncompressed_size: None,
    }
}

#[test]
fn scalar_split_input_keys_admit_exact_records_and_midpoint_remainders() {
    for start in [0.0, -1_000_000.5, 2_451_544.5, 2_451_545.123_456_7] {
        for interval in [0.0625, 1.0, 2.7] {
            let entry = body_entry(start, interval, 4);
            for fraction in [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.5, 4.0] {
                let jd = start + fraction * interval;
                let epoch = TwoPartEpoch::split(jd).unwrap();
                assert!(
                    scalar_and_split_share_argument(&entry, jd, epoch),
                    "start={start}, interval={interval}, jd={jd}"
                );
            }
        }
    }
    let entry = body_entry(2_451_545.123_456_7, 1.0, 4);
    let midpoint = entry.jd_start + 0.5 * entry.interval_days;
    assert_ne!(TwoPartEpoch::split(midpoint).unwrap().lo, 0.0);
    let jd = entry.jd_start + 0.25;
    assert!(scalar_and_split_share_argument(
        &entry,
        jd,
        TwoPartEpoch::split(jd).unwrap()
    ));
}

#[test]
fn scalar_split_input_keys_keep_original_parts_and_zero_signs() {
    let entry = body_entry(2_451_544.5, 1.0, 4);
    let jd = 2_451_545.125;
    for epoch in [
        TwoPartEpoch { hi: jd, lo: -0.0 },
        TwoPartEpoch {
            hi: jd - 0.125,
            lo: 0.125,
        },
        TwoPartEpoch {
            hi: jd,
            lo: 2.0_f64.powi(-50),
        },
    ] {
        assert_eq!(epoch.value().to_bits(), jd.to_bits());
        assert!(!scalar_and_split_share_argument(&entry, jd, epoch));
    }
    let entry = body_entry(-1.0, 2.0, 1);
    assert!(scalar_and_split_share_argument(
        &entry,
        0.0,
        TwoPartEpoch::split(0.0).unwrap()
    ));
    assert!(!scalar_and_split_share_argument(
        &entry,
        -0.0,
        TwoPartEpoch::split(-0.0).unwrap()
    ));
}

#[test]
fn scalar_split_input_keys_reject_noncartesian_and_invalid_queries() {
    let valid = body_entry(0.0, 1.0, 4);
    for coord in [COORD_ECLIPTIC, COORD_GEO_ECLIPTIC, COORD_HELIO_ECL_RETIRED] {
        let mut entry = valid.clone();
        entry.coord_type = coord;
        assert!(!scalar_and_split_share_argument(
            &entry,
            1.0,
            TwoPartEpoch::split(1.0).unwrap()
        ));
    }
    for components in [0, 2, 4] {
        let mut entry = valid.clone();
        entry.components = components;
        assert!(!scalar_and_split_share_argument(
            &entry,
            1.0,
            TwoPartEpoch::split(1.0).unwrap()
        ));
    }
    for jd in [-1.0, 5.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(!scalar_and_split_share_argument(
            &valid,
            jd,
            TwoPartEpoch { hi: jd, lo: 0.0 }
        ));
    }
}

#[test]
fn scalar_split_input_keys_preserve_split_only_midpoint_failure() {
    let mut entry = body_entry(1e308, 1.7e308, 1);
    entry.jd_end = 1.6e308;
    let jd = 1.2e308;
    let epoch = TwoPartEpoch::split(jd).unwrap();
    // Ordered midpoint overflows. Scalar tau clamps the infinity to -1;
    // the split route must still refuse to construct its infinite midpoint.
    let (_, tau) = select_segment(entry.jd_start, entry.interval_days, 1, jd).unwrap();
    assert_eq!(tau.value(), -1.0);
    assert!(matches!(
        select_segment_two_part(entry.jd_start, entry.interval_days, 1, epoch),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(!scalar_and_split_share_argument(&entry, jd, epoch));
}
