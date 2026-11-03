//! Prepared provider metadata must retain the original scalar dispatch rules.

use super::*;

/// Reconstruct selection through the original composite metadata interface.
fn original_selection(index: &TieredIndex<'_>, body: i32, jd: f64) -> Option<usize> {
    index.tiers.iter().enumerate().find_map(|(pos, (_, tier))| {
        tier.body_coverage(body)
            .filter(|(start, end)| *start <= jd && jd <= *end)
            .map(|_| pos)
    })
}

#[test]
fn prepared_tier_bounds_preserve_priority_endpoints_gaps_and_missing_bodies() {
    let images = [
        (LebTier::Extended, leb1_const_image(0, 3.0, 6.0, 8.0)),
        (LebTier::Medium, leb1_const_image(0, 2.0, 1.0, 4.0)),
        (LebTier::Base, leb1_const_image(0, 1.0, 0.0, 2.0)),
    ];
    let index = TieredIndex::new(
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
    .unwrap();
    for jd in [
        f64::NAN,
        f64::NEG_INFINITY,
        f64::INFINITY,
        -0.0,
        0.0,
        1.0,
        2.0,
        4.0,
        6.0,
        8.0,
    ]
    .into_iter()
    .flat_map(|day| [day.next_down(), day, day.next_up()])
    {
        for body in [0, 1, 999] {
            assert_eq!(
                index.selected(body, jd),
                original_selection(&index, body, jd)
            );
        }
    }
    // Out-of-range evaluation still enters the highest-priority original reader.
    assert_eq!(
        index.eval_body_any(0, 5.0).unwrap_err(),
        index.tiers[0].1.eval_body_any(0, 5.0).unwrap_err()
    );
}

#[test]
fn prepared_tier_bounds_preserve_unusual_metadata_words() {
    let bytes = leb1_const_image(0, 1.0, 0.0, 2.0);
    for (start, end) in [
        (-0.0, 0.0),
        (0.0, -0.0),
        (f64::NEG_INFINITY, f64::INFINITY),
        (f64::from_bits(0x7ff8_0000_0000_0011), 2.0),
        (0.0, f64::NAN),
    ] {
        // Isolate dispatch metadata from parsing policy and numerical evaluation.
        let mut file = FileIndex::decode(&bytes).unwrap();
        file.bodies[0].jd_start = start;
        file.bodies[0].jd_end = end;
        let index = TieredIndex::new(vec![(
            LebTier::Base,
            CompositeIndex::new(vec![(file, bytes.as_slice())]).unwrap(),
        )])
        .unwrap();
        let bounds = index.candidates[&0][0].bounds.unwrap();
        assert_eq!(bounds.0.to_bits(), start.to_bits());
        assert_eq!(bounds.1.to_bits(), end.to_bits());
        for jd in [
            f64::NAN,
            f64::NEG_INFINITY,
            f64::INFINITY,
            -0.0,
            0.0,
            1.0,
            2.0,
        ] {
            assert_eq!(index.selected(0, jd), original_selection(&index, 0, jd));
        }
    }
}

#[test]
fn prepared_body_locations_match_original_rows_and_first_source_priority() {
    let bytes = leb1_const_image(0, 1.0, 0.0, 2.0);
    let mut first = FileIndex::decode(&bytes).unwrap();
    first.bodies.clear();
    // Non-first rows, both ends of the fixed vocabulary and arbitrary IDs.
    for (row, body) in [-1, 14, 0, 22, i32::MAX, -999].into_iter().enumerate() {
        let mut entry = FileIndex::decode(&bytes).unwrap().bodies.remove(0);
        entry.body_id = body;
        entry.coord_type = if row % 2 == 0 {
            COORD_ICRS_BARY
        } else {
            COORD_ICRS_BARY_SYSTEM
        };
        entry.jd_start = row as f64;
        entry.jd_end = row as f64 + 1.0;
        first.bodies.push(entry);
    }
    let second = FileIndex::decode(&bytes).unwrap();
    let index = CompositeIndex::new(vec![
        (first.clone(), bytes.as_slice()),
        (second, bytes.as_slice()),
    ])
    .unwrap();
    for body in [-999, -1, 0, 1, 14, 21, 22, 23, i32::MAX] {
        let original = first.find_body(body).ok();
        assert_eq!(index.serving_entry(body), original);
        assert_eq!(index.has_body(body), original.is_some());
        assert_eq!(index.serving_index(body), original.map(|_| 0));
        assert_eq!(
            index.body_coverage(body),
            original.map(|entry| (entry.jd_start(), entry.jd_end()))
        );
    }
    let mut inventory = index.body_ids();
    inventory.sort_unstable();
    assert_eq!(inventory, [-999, -1, 0, 14, 22, i32::MAX]);
}
