// SPDX-License-Identifier: AGPL-3.0
//! Extraction checks use independent synthetic DAF/record construction.

use super::*;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "leph-spk-subset-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn source(
    root: &Path,
    name: &str,
    endian: DafEndian,
    kind: i32,
    target: i32,
    value: f64,
) -> SpkSubsetInput {
    let sets = if kind == 2 { 3 } else { 6 };
    let records: Vec<_> = (0..8)
        .map(|i| ChebRec {
            mid: 16.0 * i as f64 + 8.0,
            radius: 8.0,
            sets: (0..sets)
                .map(|j| vec![value + j as f64 + i as f64, 0.25])
                .collect(),
        })
        .collect();
    let bytes = cheb_image(
        endian, kind, target, 5, 0.0, 128.0, 0.0, 16.0, &records, None, name,
    );
    let path = root.join(name);
    std::fs::write(&path, bytes).unwrap();
    SpkSubsetInput {
        path,
        target,
        center: 5,
        expected_sha256: None,
    }
}

#[test]
fn extracted_records_keep_states_for_both_types_and_endian_formats() {
    let root = Scratch::new();
    for endian in [DafEndian::Little, DafEndian::Big] {
        for kind in [2, 3] {
            let input = source(&root.0, "source.bsp", endian, kind, 599, 10.0);
            let output = root.0.join("subset.bsp");
            let report = write_spk_subset(
                std::slice::from_ref(&input),
                &output,
                Some((33.25, 70.5)),
                &[(599, 5, 33.25, 70.5)],
                true,
            )
            .unwrap();
            assert_eq!(report.segments.len(), 1);
            let old = std::fs::read(&input.path).unwrap();
            let new = std::fs::read(&output).unwrap();
            let old_index = SpkIndex::parse(&old).unwrap();
            let new_index = SpkIndex::parse(&new).unwrap();
            for epoch in [
                33.25,
                40.0,
                48.0_f64.next_down(),
                48.0,
                48.0_f64.next_up(),
                64.0,
                70.5,
            ] {
                let a = eval_relative(&old, endian, &old_index, 599, epoch)
                    .unwrap()
                    .0;
                let b = eval_relative(&new, DafEndian::Little, &new_index, 599, epoch)
                    .unwrap()
                    .0;
                assert_eq!(a, b);
            }
            assert_eq!(report.bytes as usize, new.len());
            assert!(
                report.segments[0].final_address
                    < old_index
                        .select(599, 64.0)
                        .unwrap()
                        .unwrap()
                        .descriptor
                        .final_address
            );
            assert!(eval_relative(&new, DafEndian::Little, &new_index, 599, 33.0).is_err());
        }
    }
}

#[test]
fn multiple_summary_records_and_later_segment_priority_survive() {
    let root = Scratch::new();
    let mut inputs: Vec<_> = (0..30)
        .map(|i| {
            source(
                &root.0,
                &format!("{i}.bsp"),
                DafEndian::Little,
                2,
                599 + i,
                i as f64,
            )
        })
        .collect();
    inputs.push(source(&root.0, "later.bsp", DafEndian::Big, 3, 599, 900.0));
    let output = root.0.join("many.bsp");
    let report = write_spk_subset(&inputs, &output, None, &[], false).unwrap();
    assert_eq!(report.segments.len(), 31);
    let bytes = std::fs::read(&output).unwrap();
    let index = SpkIndex::parse(&bytes).unwrap();
    for target in 599..629 {
        assert!(index.select(target, 64.0).unwrap().is_some());
    }
    let result = eval_relative(&bytes, DafEndian::Little, &index, 599, 64.0).unwrap();
    assert_eq!(result.1.seg_type, 3);
    assert_eq!(result.0.position_km[0], 903.75);
}

#[test]
fn required_union_gaps_and_wrong_hashes_preserve_existing_destination() {
    let root = Scratch::new();
    let mut input = source(&root.0, "source.bsp", DafEndian::Little, 2, 599, 1.0);
    let output = root.0.join("owner.bsp");
    std::fs::write(&output, b"owner content").unwrap();
    assert!(
        write_spk_subset(
            std::slice::from_ref(&input),
            &output,
            None,
            &[(599, 5, 0.0, 130.0)],
            true
        )
        .is_err()
    );
    let second = source(&root.0, "second.bsp", DafEndian::Little, 2, 599, 2.0);
    let mut early = std::fs::read(&input.path).unwrap();
    let mut late = std::fs::read(&second.path).unwrap();
    put_f64(&mut early, RECORD_BYTES + 24 + 8, 32.0, DafEndian::Little);
    put_f64(&mut late, RECORD_BYTES + 24, 64.0, DafEndian::Little);
    std::fs::write(&input.path, early).unwrap();
    std::fs::write(&second.path, late).unwrap();
    assert!(
        write_spk_subset(
            &[input.clone(), second],
            &output,
            None,
            &[(599, 5, 0.0, 128.0)],
            true
        )
        .is_err()
    );
    input.expected_sha256 = Some("00".repeat(32));
    assert!(write_spk_subset(&[input], &output, None, &[], true).is_err());
    assert_eq!(std::fs::read(&output).unwrap(), b"owner content");
    assert!(std::fs::read_dir(&root.0).unwrap().all(|entry| {
        !entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .ends_with(".spk-tmp")
    }));
}

#[test]
fn native_planet_center_command_builds_then_preserves_outputs_on_refusal() {
    let root = Scratch::new();
    for (name, target, center) in [
        ("jup204.bsp", 599, 5),
        ("sat319.bsp", 699, 6),
        ("ura083.bsp", 799, 7),
        ("nep050.bsp", 899, 8),
        ("plu017.bsp", 999, 9),
    ] {
        let input = source(&root.0, name, DafEndian::Little, 2, target, center as f64);
        let mut bytes = std::fs::read(&input.path).unwrap();
        put_i32(
            &mut bytes,
            RECORD_BYTES + 24 + 20,
            center,
            DafEndian::Little,
        );
        std::fs::write(input.path, bytes).unwrap();
    }
    let output = root.0.join("output");
    let cache = root.0.to_str().unwrap();
    let out = output.to_str().unwrap();
    let command = [
        "generate",
        "planet-centers-base",
        "--cache-dir",
        cache,
        "--output-dir",
        out,
        "--json",
    ];
    assert_eq!(crate::cli::run(command), 0);
    let product = output.join("planet_centers_base.bsp");
    let original = std::fs::read(&product).unwrap();
    let index = SpkIndex::parse(&original).unwrap();
    for target in [599, 699, 799, 899, 999] {
        assert!(index.select(target, 64.0).unwrap().is_some());
    }
    assert_eq!(crate::cli::run(command), 1);
    assert_eq!(std::fs::read(&product).unwrap(), original);
    let mut strict = command.to_vec();
    strict.extend(["--force", "--require-full-tier"]);
    assert_eq!(crate::cli::run(strict), 1);
    assert_eq!(std::fs::read(&product).unwrap(), original);
    std::fs::write(root.0.join("sat319.bsp"), b"corrupt").unwrap();
    let mut forced = command.to_vec();
    forced.push("--force");
    assert_eq!(crate::cli::run(forced), 1);
    assert_eq!(std::fs::read(&product).unwrap(), original);
}

#[test]
fn extraction_cannot_replace_or_alias_an_input() {
    let root = Scratch::new();
    let input = source(&root.0, "source.bsp", DafEndian::Little, 2, 599, 1.0);
    let original = std::fs::read(&input.path).unwrap();
    assert!(write_spk_subset(std::slice::from_ref(&input), &input.path, None, &[], true).is_err());
    let alias = root.0.join("alias.bsp");
    std::fs::hard_link(&input.path, &alias).unwrap();
    assert!(write_spk_subset(std::slice::from_ref(&input), &alias, None, &[], true).is_err());
    assert_eq!(
        write_spk_subset(&[input], &alias, None, &[], false)
            .unwrap_err()
            .kind(),
        std::io::ErrorKind::AlreadyExists
    );
    assert_eq!(std::fs::read(&alias).unwrap(), original);
}

#[test]
fn a_nonfinite_coefficient_refuses_publication_and_removes_its_temporary() {
    let root = Scratch::new();
    let input = source(&root.0, "source.bsp", DafEndian::Little, 2, 599, f64::NAN);
    let output = root.0.join("new.bsp");
    assert!(write_spk_subset(&[input], &output, None, &[], false).is_err());
    assert!(!output.exists());
    assert_eq!(std::fs::read_dir(&root.0).unwrap().count(), 1);
}

#[test]
fn invalid_record_radius_preserves_the_destination() {
    let root = Scratch::new();
    let input = source(&root.0, "source.bsp", DafEndian::Little, 2, 599, 1.0);
    let mut bytes = std::fs::read(&input.path).unwrap();
    let index = SpkIndex::parse(&bytes).unwrap();
    let segment = index.select(599, 64.0).unwrap().unwrap();
    let radius_offset = (segment.descriptor.initial_address as usize - 1) * 8 + 8;
    put_f64(&mut bytes, radius_offset, -8.0, DafEndian::Little);
    std::fs::write(&input.path, bytes).unwrap();
    let output = root.0.join("owner.bsp");
    std::fs::write(&output, b"owner content").unwrap();
    assert!(write_spk_subset(&[input], &output, None, &[], true).is_err());
    assert_eq!(std::fs::read(output).unwrap(), b"owner content");
    assert_eq!(std::fs::read_dir(&root.0).unwrap().count(), 2);
}
