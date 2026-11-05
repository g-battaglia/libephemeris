// SPDX-License-Identifier: AGPL-3.0
//! Record support must agree with evaluation, including non-dyadic joins.

use super::*;

fn entry(start: f64, end: f64, interval: f64, count: u32) -> BodyEntry {
    BodyEntry {
        body_id: 7,
        coord_type: COORD_ICRS_BARY,
        segment_count: count,
        jd_start: start,
        jd_end: end,
        interval_days: interval,
        degree: 1,
        components: 3,
        data_offset: 0,
        compressed_size: None,
        uncompressed_size: None,
    }
}

fn assert_owned_interval(entry: &BodyEntry, jd: f64) -> (f64, f64) {
    let selected = |epoch| {
        select_segment(
            entry.jd_start,
            entry.interval_days,
            entry.segment_count,
            epoch,
        )
        .unwrap()
        .0
    };
    let owner = selected(jd);
    let (start, end) = entry.smooth_interval(jd).unwrap().unwrap();
    assert!(entry.jd_start <= start && start <= jd);
    assert!(jd <= end && end <= entry.jd_end);
    for epoch in [start, jd, end] {
        assert_eq!(selected(epoch), owner);
        assert_eq!(
            select_segment_two_part(
                entry.jd_start,
                entry.interval_days,
                entry.segment_count,
                TwoPartEpoch::split(epoch).unwrap(),
            )
            .unwrap()
            .0,
            owner,
        );
    }
    if start > entry.jd_start {
        assert_ne!(selected(start.next_down()), owner);
    }
    if end < entry.jd_end {
        let same_next_record = matches!(
            select_segment(
                entry.jd_start,
                entry.interval_days,
                entry.segment_count,
                end.next_up(),
            ),
            Ok((index, _)) if index == owner
        );
        assert!(!same_next_record);
    }
    (start, end)
}

#[test]
fn rounded_non_dyadic_join_remains_in_its_selected_record() {
    // Synthetic metadata reproduces multiplication/division disagreement;
    // no astronomical coefficients or captured output values are required.
    let grid = entry(2_396_757.5, 2_506_330.5, 63.966725043782844, 1713);
    let join = 2_433_282.5_f64;
    assert_eq!(
        select_segment(grid.jd_start, grid.interval_days, grid.segment_count, join)
            .unwrap()
            .0,
        571,
    );
    assert_eq!(assert_owned_interval(&grid, join).0, join);
    for jd in [join.next_down(), join.next_up(), grid.jd_start, grid.jd_end] {
        assert_owned_interval(&grid, jd);
    }
}

#[test]
fn record_support_matches_selector_across_uniform_grids_and_clipped_coverage() {
    for grid in [
        entry(0.0, 6.0, 2.0, 3),
        entry(2_451_545.0, 2_451_546.13, 0.1, 20),
        entry(-1e16, 1e16, 1e16, 2),
        entry(1.0, 1.0_f64.next_up(), f64::MIN_POSITIVE, 10),
    ] {
        for n in 0..=grid.segment_count {
            let join = grid.jd_start + f64::from(n) * grid.interval_days;
            for jd in [join.next_down(), join, join.next_up()] {
                if grid.jd_start <= jd && jd <= grid.jd_end {
                    assert_owned_interval(&grid, jd);
                }
            }
        }
        assert_owned_interval(&grid, grid.jd_start);
        assert_owned_interval(&grid, grid.jd_end);
        assert_eq!(
            grid.smooth_interval(grid.jd_start.next_down()).unwrap(),
            None
        );
        assert_eq!(grid.smooth_interval(grid.jd_end.next_up()).unwrap(), None);
    }
}

#[test]
fn nutation_support_uses_the_same_selected_record_bounds() {
    let grid = entry(2_396_757.5, 2_506_330.5, 63.966725043782844, 1713);
    let mut payload = Vec::new();
    for value in [grid.jd_start, grid.jd_end, grid.interval_days] {
        push_f64(&mut payload, value);
    }
    for value in [1, 2, grid.segment_count, 0] {
        push_u32(&mut payload, value);
    }
    payload.resize(
        NUTATION_HEADER_SIZE + grid.segment_count as usize * 4 * 8,
        0,
    );
    let file = aux_image(vec![(SECTION_NUTATION, payload)]);
    let index = FileIndex::decode(&file).unwrap();
    let join = 2_433_282.5_f64;
    for jd in [
        join.next_down(),
        join,
        join.next_up(),
        grid.jd_start,
        grid.jd_end,
    ] {
        assert_eq!(
            index
                .nutation_smooth_interval_source(LebSource::Slice(&file), jd)
                .unwrap(),
            Some(assert_owned_interval(&grid, jd)),
        );
    }
}

#[test]
fn non_dyadic_seams_across_epoch_scales_have_maximal_owned_support() {
    for start in [-1e12, -0.3, 0.0, 1.0, 2_451_545.0, 1e12] {
        for interval in [0.1, 0.3, 1.7, 32.0, 63.966725043782844] {
            let grid = entry(start, start + 12.5 * interval, interval, 20);
            for n in 0..=13 {
                let seam = start + f64::from(n) * interval;
                for jd in [seam.next_down(), seam, seam.next_up()] {
                    if start <= jd && jd <= grid.jd_end {
                        assert_owned_interval(&grid, jd);
                    }
                }
            }
        }
    }
}

#[test]
fn distant_overflowing_coverage_does_not_exclude_a_valid_central_epoch() {
    for count in [2, 3] {
        let grid = entry(-f64::MAX, f64::MAX, f64::MAX, count);
        let (_, end) = assert_owned_interval(&grid, 0.0);
        assert!(end < grid.jd_end);
        assert!(matches!(
            select_segment(grid.jd_start, grid.interval_days, count, end.next_up()),
            Err(Error::NonFiniteArgument { value }) if value == f64::INFINITY
        ));
    }
}

#[test]
fn retained_split_tail_does_not_change_collapsed_word_record_ownership() {
    let grid = entry(2_396_757.5, 2_506_330.5, 63.966725043782844, 1713);
    let jd = 2_433_282.5_f64;
    let tail = (jd.next_up() - jd) / 4.0;
    for low in [-tail, tail] {
        let epoch = TwoPartEpoch { hi: jd, lo: low };
        assert_eq!(epoch.value(), jd);
        let owner =
            select_segment_two_part(grid.jd_start, grid.interval_days, grid.segment_count, epoch)
                .unwrap()
                .0;
        assert_eq!(owner, 571);
        assert_owned_interval(&grid, epoch.value());
    }
}

#[test]
fn nominal_polynomial_join_is_not_promoted_to_smooth_interior_by_roundoff() {
    let grid = entry(-3_100_015.5, 4_899_984.5, 4.0, 2_000_000);
    let join = 2_415_020.5_f64;
    let selected = |epoch| {
        select_segment(grid.jd_start, grid.interval_days, grid.segment_count, epoch)
            .unwrap()
            .0
    };
    assert_ne!(selected(join.next_down()), selected(join));
    assert_eq!(
        assert_owned_interval(&grid, join.next_down()).1,
        join.next_down()
    );
    let (start, end) = assert_owned_interval(&grid, join);
    assert_eq!(start, join);
    assert!(join < end);
    // The next representable epoch is a genuine interior point, preserving
    // the distinction between a polynomial join and its neighboring samples.
    let (start, end) = assert_owned_interval(&grid, join.next_up());
    assert!(start < join.next_up() && join.next_up() < end);
}

#[test]
fn selector_preserves_outer_clamps_and_handles_coincident_starts() {
    for count in [1, 10, u32::MAX] {
        let start = 1.0_f64;
        let interval = f64::MIN_POSITIVE;
        // All nominal starts collapse to the first epoch on this tiny grid.
        for jd in [start.next_down(), start] {
            assert_eq!(select_segment(start, interval, count, jd).unwrap().0, 0);
        }
        assert_eq!(
            select_segment(start, interval, count, start.next_up())
                .unwrap()
                .0,
            (count - 1) as usize,
        );
    }
    for (jd, expected) in [(-1.0, 0), (0.0, 0), (2.0, 1), (6.0, 2), (7.0, 2)] {
        assert_eq!(select_segment(0.0, 2.0, 3, jd).unwrap().0, expected);
    }
}
