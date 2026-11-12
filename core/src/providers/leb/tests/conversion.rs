//! Independent hand-packed container and filesystem publication controls.
use super::*;

struct Sandbox(PathBuf);
impl Sandbox {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "leb-convert-tests-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn input(&self, bytes: &[u8]) -> PathBuf {
        let path = self.0.join("source.leb");
        fs::write(&path, bytes).unwrap();
        path
    }
    fn output(&self) -> PathBuf {
        self.0.join("output.leb2")
    }
    fn assert_no_temporary(&self) {
        assert!(fs::read_dir(&self.0).unwrap().all(|entry| {
            !entry
                .unwrap()
                .file_name()
                .to_string_lossy()
                .starts_with(".leb2-convert-")
        }));
    }
}
impl Drop for Sandbox {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn u32_at(bytes: &mut [u8], at: usize, value: u32) {
    bytes[at..at + 4].copy_from_slice(&value.to_le_bytes());
}
fn u64_at(bytes: &mut [u8], at: usize, value: u64) {
    bytes[at..at + 8].copy_from_slice(&value.to_le_bytes());
}
fn f64_at(bytes: &mut [u8], at: usize, value: f64) {
    bytes[at..at + 8].copy_from_slice(&value.to_le_bytes());
}

/// This encoder writes LEB1 offsets/fields directly; it calls no converter,
/// reorder, shuffle or LEB2 encoder. The asymmetric finite IEEE words make
/// component/segment permutation, mantissa truncation and zero-sign loss visible.
fn fixture() -> Vec<u8> {
    let index_at = 64 + 4 * 24;
    let coefficients_at = index_at + 2 * 52;
    let raw_size = (5 + 3) * 6 * 8;
    let delta_at = coefficients_at + raw_size;
    let opaque_at = delta_at + 40;
    let mut bytes = vec![0; opaque_at + COPY_BLOCK + 3];
    bytes[..4].copy_from_slice(b"LEB1");
    for (at, value) in [(4, 1), (8, 4), (12, 2)] {
        u32_at(&mut bytes, at, value);
    }
    f64_at(&mut bytes, 16, 10.0);
    f64_at(&mut bytes, 24, 112.0);
    f64_at(&mut bytes, 32, 2450000.25);
    bytes[44..64].fill(0x5a);
    for (ordinal, id, offset, size) in [
        (0, 0, index_at, 104),
        (1, 1, coefficients_at, raw_size),
        (2, 3, delta_at, 40),
        (3, 99, opaque_at, COPY_BLOCK + 3),
    ] {
        let at = 64 + ordinal * 24;
        u32_at(&mut bytes, at, id);
        u32_at(&mut bytes, at + 4, 0x1234_5678);
        u64_at(&mut bytes, at + 8, offset as u64);
        u64_at(&mut bytes, at + 16, size as u64);
    }
    // Deliberately unequal component widths/domains, ordered ids7 then1.
    for (ordinal, id, coord, segments, start, end, interval, degree, components, offset) in [
        (0, 7i32, 0, 5, 10.0, 20.0, 2.0, 2, 2, coefficients_at),
        (
            1,
            1i32,
            1,
            3,
            100.0,
            112.0,
            4.0,
            1,
            3,
            coefficients_at + 240,
        ),
    ] {
        let at = index_at + ordinal * 52;
        bytes[at..at + 4].copy_from_slice(&id.to_le_bytes());
        u32_at(&mut bytes, at + 4, coord);
        u32_at(&mut bytes, at + 8, segments);
        f64_at(&mut bytes, at + 12, start);
        f64_at(&mut bytes, at + 20, end);
        f64_at(&mut bytes, at + 28, interval);
        u32_at(&mut bytes, at + 36, degree);
        u32_at(&mut bytes, at + 40, components);
        u64_at(&mut bytes, at + 44, offset as u64);
    }
    let patterns = [
        0u64,
        0x8000_0000_0000_0000,
        1,
        0x8000_0000_0000_0001,
        0x3ff0_0000_0000_0001,
        0xbfd0_1234_5678_9abc,
    ];
    for ordinal in 0..raw_size / 8 {
        u64_at(
            &mut bytes,
            coefficients_at + ordinal * 8,
            patterns[ordinal % patterns.len()] ^ ((ordinal / 6) as u64),
        );
    }
    u32_at(&mut bytes, delta_at, 2);
    for (ordinal, jd, delta) in [(0, 10.0, 0.001), (1, 112.0, 0.002)] {
        f64_at(&mut bytes, delta_at + 8 + ordinal * 16, jd);
        f64_at(&mut bytes, delta_at + 16 + ordinal * 16, delta);
    }
    for (ordinal, byte) in bytes[opaque_at..].iter_mut().enumerate() {
        *byte = (ordinal % 251) as u8;
    }
    bytes
}

fn options() -> Leb2ConvertOptions {
    Leb2ConvertOptions {
        chunk_interval_days: 8.0,
        max_chunk_bytes: 96,
        ..Leb2ConvertOptions::default()
    }
}

fn check_words_and_values(input: &LebFile, output: &LebFile) -> u64 {
    let mut total = 0;
    for entry in output.index().unwrap().bodies() {
        let original = input
            .index()
            .unwrap()
            .bodies()
            .iter()
            .find(|body| body.body_id == entry.body_id)
            .unwrap();
        assert_eq!(entry.coord_type, original.coord_type);
        assert_eq!(entry.components, original.components);
        assert_eq!(entry.degree, original.degree);
        assert_eq!(entry.segment_count, original.segment_count);
        assert_eq!(entry.jd_start.to_bits(), original.jd_start.to_bits());
        assert_eq!(entry.jd_end.to_bits(), original.jd_end.to_bits());
        assert_eq!(
            entry.interval_days.to_bits(),
            original.interval_days.to_bits()
        );
        let (_, chunks) = decode::decode_chunk_index_source(
            output.source().unwrap(),
            entry.data_offset,
            entry.body_id,
        )
        .unwrap();
        for (ordinal, chunk) in chunks.values.iter().enumerate() {
            let frame = output
                .read_range(chunk.blob_offset, chunk.compressed_size)
                .unwrap();
            let decoded = decompress_chunk(
                frame,
                chunk,
                entry.degree,
                entry.components,
                entry.body_id,
                ordinal,
            )
            .unwrap();
            let segment_bytes = (entry.degree as u64 + 1) * entry.components as u64 * 8;
            let raw = input
                .read_range(
                    original.data_offset + chunk.segment_start as u64 * segment_bytes,
                    chunk.uncompressed_size,
                )
                .unwrap();
            assert_eq!(decoded, raw);
            assert!(chunk.uncompressed_size <= 96);
            total += 1;
        }
        for segment in 0..entry.segment_count {
            for fraction in [0.0, 0.25, 0.5, 0.75, 1.0] {
                let jd = entry.jd_start + (f64::from(segment) + fraction) * entry.interval_days;
                let a = input.eval_body_any(entry.body_id, jd).unwrap();
                let b = output.eval_body_any(entry.body_id, jd).unwrap();
                let bits = |values: (Vec<f64>, Vec<f64>)| {
                    values
                        .0
                        .into_iter()
                        .chain(values.1)
                        .map(f64::to_bits)
                        .collect::<Vec<_>>()
                };
                assert_eq!(bits(a), bits(b));
            }
        }
    }
    total
}

#[test]
fn lossless_conversion_preserves_all_chunk_words_metadata_auxiliary_and_public_rates() {
    let sandbox = Sandbox::new();
    let bytes = fixture();
    let source = sandbox.input(&bytes);
    let destination = sandbox.output();
    let report = convert_leb1_to_leb2(&source, &destination, &options()).unwrap();
    assert_eq!(report.body_ids, [7, 1]);
    assert_eq!(report.chunks, 5);
    assert_eq!(
        report.bytes_written,
        fs::metadata(&destination).unwrap().len()
    );
    assert_eq!(fs::read(&source).unwrap(), bytes);
    let input = LebFile::open(&source).unwrap();
    let output = LebFile::open(&destination).unwrap();
    assert_eq!(output.format().unwrap(), LebFormat::Leb2V2);
    assert_eq!(
        output
            .index()
            .unwrap()
            .header()
            .generation_epoch()
            .to_bits(),
        input.index().unwrap().header().generation_epoch().to_bits()
    );
    assert_eq!(&output.bytes().unwrap()[44..64], &bytes[44..64]);
    assert_eq!(check_words_and_values(&input, &output), 5);
    for id in [SECTION_DELTA_T, 99] {
        let a = input.index().unwrap().section(id).unwrap();
        let b = output.index().unwrap().section(id).unwrap();
        assert_eq!(
            input.read_range(a.offset, a.size).unwrap(),
            output.read_range(b.offset, b.size).unwrap()
        );
        let ordinal = output
            .index()
            .unwrap()
            .sections
            .iter()
            .position(|s| s.section_id == id)
            .unwrap();
        assert_eq!(
            &output.bytes().unwrap()[64 + ordinal * 24 + 4..64 + ordinal * 24 + 8],
            &0x1234_5678_u32.to_le_bytes()
        );
    }
    sandbox.assert_no_temporary();
}

#[test]
fn conversion_selection_and_auxiliary_omission_are_explicit() {
    let sandbox = Sandbox::new();
    let source = sandbox.input(&fixture());
    let mut options = options();
    options.body_ids = Some(vec![1]);
    options.include_auxiliary = false;
    let report = convert_leb1_to_leb2(&source, &sandbox.output(), &options).unwrap();
    assert_eq!(report.body_ids, [1]);
    assert_eq!(report.chunks, 2);
    let output = LebFile::open(&sandbox.output()).unwrap();
    assert_eq!(output.index().unwrap().sections().len(), 2);
    assert!(output.index().unwrap().section(99).is_none());
    assert!(output.index().unwrap().section(SECTION_DELTA_T).is_none());
    check_words_and_values(&LebFile::open(&source).unwrap(), &output);
}

#[test]
fn conversion_rejects_invalid_options_budgets_ambiguity_shapes_and_nonfinite_words() {
    let sandbox = Sandbox::new();
    let original = fixture();
    let source = sandbox.input(&original);
    let mut invalid_options = Vec::new();
    for interval in [
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        f64::MAX,
        f64::MIN_POSITIVE,
    ] {
        invalid_options.push(Leb2ConvertOptions {
            chunk_interval_days: interval,
            ..options()
        });
    }
    invalid_options.push(Leb2ConvertOptions {
        max_chunk_bytes: 47,
        ..options()
    });
    invalid_options.push(Leb2ConvertOptions {
        body_ids: Some(vec![1, 1]),
        ..options()
    });
    invalid_options.push(Leb2ConvertOptions {
        body_ids: Some(vec![12345]),
        ..options()
    });
    for options in invalid_options {
        assert!(convert_leb1_to_leb2(&source, &sandbox.output(), &options).is_err());
        assert!(!sandbox.output().exists());
        sandbox.assert_no_temporary();
    }
    let index_at = 64 + 4 * 24;
    let data_at = index_at + 2 * 52;
    let mut corruptions = Vec::new();
    let mut nonfinite = original.clone();
    u64_at(&mut nonfinite, data_at, f64::NAN.to_bits());
    corruptions.push(nonfinite);
    let mut duplicate_body = original.clone();
    u32_at(&mut duplicate_body, index_at + 52, 7);
    corruptions.push(duplicate_body);
    let mut duplicate_section = original.clone();
    u32_at(&mut duplicate_section, 64 + 3 * 24, SECTION_DELTA_T);
    corruptions.push(duplicate_section);
    let mut conflict = original.clone();
    u32_at(&mut conflict, 64 + 3 * 24, SECTION_COMPRESSED_CHEBYSHEV);
    corruptions.push(conflict);
    let mut bad_shape = original.clone();
    u32_at(&mut bad_shape, index_at + 40, 0);
    corruptions.push(bad_shape);
    let mut short = original.clone();
    u64_at(&mut short, index_at + 44, original.len() as u64 - 8);
    corruptions.push(short);
    let mut overlap = original.clone();
    u64_at(&mut overlap, index_at + 44, 64);
    corruptions.push(overlap);
    for bytes in corruptions {
        fs::write(&source, &bytes).unwrap();
        fs::write(sandbox.output(), b"preserve existing output").unwrap();
        let options = Leb2ConvertOptions {
            overwrite: true,
            ..options()
        };
        assert!(convert_leb1_to_leb2(&source, &sandbox.output(), &options).is_err());
        assert_eq!(
            fs::read(sandbox.output()).unwrap(),
            b"preserve existing output"
        );
        assert_eq!(fs::read(&source).unwrap(), bytes);
        sandbox.assert_no_temporary();
    }
}

#[test]
fn conversion_preflight_and_exclusive_publication_refuse_existing_or_raced_targets() {
    let sandbox = Sandbox::new();
    let source = sandbox.input(&fixture());
    let destination = sandbox.output();
    validate_leb2_conversion_target(&source, &destination, false).unwrap();
    fs::write(&destination, b"raced target").unwrap();
    assert!(convert_leb1_to_leb2(&source, &destination, &options()).is_err());
    let (temporary, mut file) = temporary(&sandbox.0).unwrap();
    file.write_all(b"verified temporary").unwrap();
    file.sync_all().unwrap();
    assert!(publish(&temporary.0, &destination, false).is_err());
    assert_eq!(fs::read(&destination).unwrap(), b"raced target");
    drop(file);
    drop(temporary);
    sandbox.assert_no_temporary();
    let nested = sandbox.0.join("missing/deeper/out.leb2");
    validate_leb2_conversion_target(&source, &nested, false).unwrap();
    assert!(!sandbox.0.join("missing").exists());
    convert_leb1_to_leb2(&source, &nested, &options()).unwrap();
    assert_eq!(
        LebFile::open(&nested).unwrap().format().unwrap(),
        LebFormat::Leb2V2
    );
    assert_eq!(fs::read(&destination).unwrap(), b"raced target");
    let options = Leb2ConvertOptions {
        overwrite: true,
        ..options()
    };
    convert_leb1_to_leb2(&source, &destination, &options).unwrap();
    assert_eq!(
        LebFile::open(&destination).unwrap().format().unwrap(),
        LebFormat::Leb2V2
    );
    sandbox.assert_no_temporary();
}

#[cfg(unix)]
#[test]
fn conversion_refuses_canonical_hardlink_symlink_and_dangling_output_aliases() {
    use std::os::unix::fs::symlink;
    let sandbox = Sandbox::new();
    let original = fixture();
    let source = sandbox.input(&original);
    let options = Leb2ConvertOptions {
        overwrite: true,
        ..options()
    };
    assert!(convert_leb1_to_leb2(&source, &source, &options).is_err());
    let hardlink = sandbox.0.join("hardlink.leb");
    fs::hard_link(&source, &hardlink).unwrap();
    assert!(convert_leb1_to_leb2(&source, &hardlink, &options).is_err());
    let link = sandbox.0.join("symlink.leb");
    symlink(&source, &link).unwrap();
    assert!(convert_leb1_to_leb2(&source, &link, &options).is_err());
    let dangling = sandbox.0.join("dangling.leb");
    symlink(sandbox.0.join("does-not-exist"), &dangling).unwrap();
    assert!(validate_leb2_conversion_target(&source, &dangling, false).is_err());
    assert!(convert_leb1_to_leb2(&source, &dangling, &Leb2ConvertOptions::default()).is_err());
    assert!(
        fs::symlink_metadata(&dangling)
            .unwrap()
            .file_type()
            .is_symlink()
    );
    assert_eq!(fs::read(&source).unwrap(), original);
    sandbox.assert_no_temporary();
}

#[cfg(unix)]
#[test]
fn conversion_guard_tracks_acquired_inode_across_atomic_source_path_replacement() {
    let sandbox = Sandbox::new();
    let original = fixture();
    let source = sandbox.input(&original);
    let input = LebFile::open_indexed(&source, IndexedLimits::default()).unwrap();
    let LebSource::Positional(acquired) = input.source().unwrap() else {
        unreachable!()
    };
    let identity = acquired.acquired_metadata().clone();
    let old_alias = sandbox.0.join("old-inode.leb");
    fs::hard_link(&source, &old_alias).unwrap();
    let replacement = sandbox.0.join("replacement.leb");
    let mut changed = original.clone();
    changed[44] ^= 1;
    fs::write(&replacement, changed).unwrap();
    fs::rename(replacement, &source).unwrap();
    check_source(&input).unwrap();
    assert!(protect_input(&source, &identity, &old_alias).is_err());
    assert_eq!(
        &*input
            .source()
            .unwrap()
            .range(0, original.len(), "original snapshot")
            .unwrap(),
        original.as_slice()
    );
    assert!(!same_inode(&identity, &fs::metadata(source).unwrap()));
}

#[test]
fn conversion_nominal_span_keeps_single_wider_segments_losslessly() {
    let sandbox = Sandbox::new();
    let mut bytes = fixture();
    let index_at = 64 + 4 * 24;
    f64_at(&mut bytes, index_at + 20, 25010.0);
    f64_at(&mut bytes, index_at + 28, 5000.0);
    f64_at(&mut bytes, 24, 25010.0);
    let source = sandbox.input(&bytes);
    let options = Leb2ConvertOptions {
        max_chunk_bytes: 96,
        ..Leb2ConvertOptions::default()
    };
    let report = convert_leb1_to_leb2(&source, &sandbox.output(), &options).unwrap();
    assert_eq!(report.chunks, 7); // Five wide segments, two byte-bounded chunks.
    assert_eq!(
        check_words_and_values(
            &LebFile::open(&source).unwrap(),
            &LebFile::open(&sandbox.output()).unwrap()
        ),
        7
    );
    assert_eq!(fs::read(source).unwrap(), bytes);
}
