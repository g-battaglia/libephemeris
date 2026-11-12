use super::*;

fn polynomial_image(specifications: &[(i32, i32, f64, f64)]) -> Vec<u8> {
    let specs: Vec<_> = specifications
        .iter()
        .map(|&(target, center, start, end)| SegSpec {
            target,
            center,
            start,
            end,
            frame: FRAME_J2000,
            seg_type: TYPE_CHEBYSHEV_POS,
            name: "SYNTHETIC RECORD DOMAIN",
        })
        .collect();
    let mut image = build_image(DafEndian::Little, &[specs]);
    for (ordinal, &(_, _, start, end)) in specifications.iter().enumerate() {
        let initial = image.len() / 8 + 1;
        // One constant polynomial; its exact interval is the descriptor.
        let words = [
            (start + end) / 2.0,
            (end - start) / 2.0,
            1.0,
            2.0,
            3.0,
            start,
            end - start,
            5.0,
            1.0,
        ];
        for word in words {
            image.extend_from_slice(&word.to_le_bytes());
        }
        let at = RECORD_BYTES + 3 * 8 + ordinal * SPK_SUMMARY_WORDS * 8;
        put_i32(&mut image, at + 32, initial as i32, DafEndian::Little);
        let final_word = image.len() / 8;
        put_i32(&mut image, at + 36, final_word as i32, DafEndian::Little);
    }
    image
}

#[test]
fn ssb_smooth_support_intersects_actual_center_records_only() {
    let mut registry = SpkRegistry::new();
    registry
        .register_image(
            "chain",
            polynomial_image(&[
                (501, 5, 0.0, 8.0 * SECONDS_PER_DAY),
                (5, 0, 0.0, 2.0 * SECONDS_PER_DAY),
                (606, 0, 0.0, SECONDS_PER_DAY / 16.0),
            ]),
        )
        .unwrap();
    registry.register_body(9001, "chain", 501).unwrap();
    let bounds = registry
        .geometric_ssb_record_offsets(9001, J2000_JD_TDB, 1.0)
        .unwrap();
    assert_eq!(bounds, (-1.0, 1.0));
    assert!(
        registry
            .eval_geometric_ssb_au_two_part_jd(9001, J2000_JD_TDB, 1.0)
            .is_ok()
    );
}

#[test]
fn later_descriptor_endpoints_bound_the_selected_record_without_envelope_widening() {
    let mut registry = SpkRegistry::new();
    registry
        .register_image(
            "priority",
            polynomial_image(&[
                (499, 0, 0.0, 8.0 * SECONDS_PER_DAY),
                (499, 0, 3.0 * SECONDS_PER_DAY, 5.0 * SECONDS_PER_DAY),
            ]),
        )
        .unwrap();
    registry.register_body(4, "priority", 499).unwrap();
    let left = registry
        .geometric_ssb_record_offsets(4, J2000_JD_TDB, 1.0)
        .unwrap();
    assert_eq!(left, (-1.0, 2.0_f64.next_down()));
    let middle = registry
        .geometric_ssb_record_offsets(4, J2000_JD_TDB, 4.0)
        .unwrap();
    assert_eq!(middle, (-1.0, 1.0));
    let right = registry
        .geometric_ssb_record_offsets(4, J2000_JD_TDB, 6.0)
        .unwrap();
    assert_eq!(right, ((-1.0_f64).next_up(), 2.0));
}

#[test]
fn absent_automatic_offset_is_not_a_smooth_parent_domain_across_its_start() {
    let mut registry = SpkRegistry::new();
    registry
        .register_image(
            "main",
            polynomial_image(&[(5, 0, 0.0, 8.0 * SECONDS_PER_DAY)]),
        )
        .unwrap();
    registry.register_body(5, "main", 5).unwrap();
    registry
        .register_image(
            "centers",
            polynomial_image(&[(599, 5, 3.0 * SECONDS_PER_DAY, 5.0 * SECONDS_PER_DAY)]),
        )
        .unwrap();
    registry
        .register_automatic_planet_centers("centers")
        .unwrap();
    let before = registry
        .geometric_ssb_record_offsets(599, J2000_JD_TDB, 1.0)
        .unwrap();
    assert_eq!(before, (-1.0, 2.0_f64.next_down()));
    let inside = registry
        .geometric_ssb_record_offsets(599, J2000_JD_TDB, 4.0)
        .unwrap();
    assert_eq!(inside, (-1.0, 1.0));
}
