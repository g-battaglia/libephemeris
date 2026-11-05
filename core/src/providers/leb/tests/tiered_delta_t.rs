//! Synthetic container controls for tiered Delta-T selection.
//! No ephemeris payload or source astronomical expectation is stored.

use super::*;

/// Make an uncompressed Delta-T image using only header and pair words.
fn table_image(pairs: &[(f64, f64)]) -> Vec<u8> {
    aux_image(vec![(SECTION_DELTA_T, delta_t_payload(pairs))])
}

/// Build concrete parsed tier composites; input order cannot override priority.
fn tiered_tables(images: &[(LebTier, Vec<u8>)]) -> TieredIndex<'_> {
    TieredIndex::new(
        images
            .iter()
            .map(|(tier, bytes)| {
                (
                    *tier,
                    CompositeIndex::new(vec![(
                        FileIndex::decode(bytes).unwrap(),
                        bytes.as_slice(),
                    )])
                    .unwrap(),
                )
            })
            .collect(),
    )
    .unwrap()
}

#[test]
fn tiered_delta_t_absent_and_empty_carriers_continue_in_priority_order() {
    let images = [
        (LebTier::Extended, table_image(&[(0.0, 3.0)])),
        (LebTier::Base, aux_image(vec![])),
        (LebTier::Medium, table_image(&[])),
    ];
    let tiered = tiered_tables(&images);
    assert_eq!(tiered.delta_t(8.0).unwrap().to_bits(), 3.0_f64.to_bits());
}

#[test]
fn tiered_delta_t_first_selected_table_preserves_clamps_knots_and_zero() {
    let images = [
        (LebTier::Medium, table_image(&[(0.0, 9.0)])),
        (LebTier::Base, table_image(&[(4.0, -0.0), (12.0, 2.0)])),
    ];
    let tiered = tiered_tables(&images);
    for (jd, expected) in [
        (0.0, -0.0_f64),
        (4.0, -0.0),
        (8.0, 1.0),
        (12.0, 2.0),
        (16.0, 2.0),
    ] {
        assert_eq!(tiered.delta_t(jd).unwrap().to_bits(), expected.to_bits());
    }
}

#[test]
fn tiered_delta_t_single_nan_node_does_not_hide_selected_interval_error() {
    let images = [
        (LebTier::Base, table_image(&[(f64::NAN, 1.0)])),
        (LebTier::Medium, table_image(&[(0.0, 2.0)])),
    ];
    let direct = FileIndex::decode(&images[0].1)
        .unwrap()
        .delta_t(&images[0].1, 8.0)
        .unwrap_err();
    assert!(
        matches!(&direct, Error::LebCorrupt { what, offset: 0 } if what == "LEB Delta-T table (no interpolation interval)")
    );
    assert_eq!(tiered_tables(&images).delta_t(8.0).unwrap_err(), direct);
}

#[test]
fn tiered_delta_t_selected_lazy_capacity_error_does_not_select_lower_tier() {
    let mut payload = delta_t_payload(&[(0.0, 1.0)]);
    payload[0..4].copy_from_slice(&2_u32.to_le_bytes());
    let images = [
        (LebTier::Base, aux_image(vec![(SECTION_DELTA_T, payload)])),
        (LebTier::Medium, table_image(&[(0.0, 2.0)])),
    ];
    // Native structural decode is lazy for this auxiliary table. The source
    // constructor refuses this inflated count before a reader is installed.
    let direct = FileIndex::decode(&images[0].1)
        .unwrap()
        .delta_t(&images[0].1, 8.0)
        .unwrap_err();
    assert!(matches!(&direct, Error::LebCorrupt { .. }));
    assert_eq!(tiered_tables(&images).delta_t(8.0).unwrap_err(), direct);
}

#[test]
fn tiered_delta_t_range_selection_precedes_selected_table_errors() {
    let images = [
        (LebTier::Base, table_image(&[(f64::NAN, 1.0)])),
        (LebTier::Medium, table_image(&[(0.0, 2.0)])),
    ];
    let mut images = images;
    images[0].1[24..32].copy_from_slice(&4.0_f64.to_le_bytes());
    assert_eq!(
        tiered_tables(&images).delta_t(8.0).unwrap().to_bits(),
        2.0_f64.to_bits()
    );
}

#[test]
fn tiered_delta_t_selected_nonfinite_input_refusal_is_not_rewritten() {
    let mut image = table_image(&[(0.0, 1.0)]);
    image[16..24].copy_from_slice(&f64::NEG_INFINITY.to_le_bytes());
    image[24..32].copy_from_slice(&f64::INFINITY.to_le_bytes());
    let images = [(LebTier::Base, image)];
    for value in [f64::NEG_INFINITY, f64::INFINITY] {
        let error = tiered_tables(&images).delta_t(value).unwrap_err();
        assert!(
            matches!(error, Error::NonFiniteArgument { value: actual } if actual.to_bits() == value.to_bits())
        );
    }
}

#[test]
fn tiered_delta_t_stored_nan_value_is_success_and_not_fallback() {
    let nan = f64::from_bits(0x7ff8_0000_0000_0042);
    let images = [
        (LebTier::Base, table_image(&[(0.0, nan)])),
        (LebTier::Medium, table_image(&[(0.0, 2.0)])),
    ];
    assert_eq!(
        tiered_tables(&images).delta_t(8.0).unwrap().to_bits(),
        nan.to_bits()
    );
}

#[test]
fn tiered_delta_t_nan_request_and_exhaustion_keep_no_coverage_error() {
    let images = [(LebTier::Base, table_image(&[(0.0, 1.0)]))];
    for jd in [f64::NAN, 17.0] {
        assert!(
            matches!(tiered_tables(&images).delta_t(jd), Err(Error::LebCorrupt { what, offset: 0 }) if what.starts_with("LEB tiered Delta-T (no coverage at JD "))
        );
    }
}

#[test]
fn tiered_delta_t_short_header_is_terminal_even_when_probe_reports_absence() {
    let images = [
        (
            LebTier::Base,
            aux_image(vec![(SECTION_DELTA_T, vec![0; 4])]),
        ),
        (LebTier::Medium, table_image(&[(0.0, 2.0)])),
    ];
    assert!(
        matches!(tiered_tables(&images).delta_t(8.0), Err(Error::LebCorrupt { what, .. }) if what == "LEB Delta-T header")
    );
}

#[test]
fn tiered_delta_t_empty_header_must_fit_its_declared_section() {
    let mut image = table_image(&[]);
    // Physical header bytes remain present, but the directory claims only
    // four bytes. Zero entries do not make a negative payload capacity valid.
    image[104..112].copy_from_slice(&4_u64.to_le_bytes());
    let images = [
        (LebTier::Base, image),
        (LebTier::Medium, table_image(&[(0.0, 2.0)])),
    ];
    assert!(
        matches!(tiered_tables(&images).delta_t(8.0), Err(Error::LebCorrupt { what, .. }) if what == "LEB Delta-T table")
    );
}
