use super::*;

fn core_bytes() -> Vec<u8> {
    // These tests exercise the published manifest, not a same-named local image.
    std::fs::read(crate::test_data::published_core_path()).unwrap()
}

#[test]
fn manifest_length_rejection_leaves_immutable_digest_uncomputed() {
    let mut bytes = core_bytes();
    bytes.push(0);
    let file = LebFile::from_bytes("base_core.leb2".into(), bytes).unwrap();
    let asset = crate::download::asset_by_name("base_core.leb2").unwrap();
    for _ in 0..2 {
        assert!(!asset.matches_image(&file).unwrap());
        assert!(file.image.as_ref().unwrap().sha256.get().is_none());
    }
    // Explicit SHA remains available and identifies the actual full image.
    assert_ne!(file.sha256().unwrap(), asset.sha256());
    let mut closed = file;
    closed.close();
    assert!(matches!(
        asset.matches_image(&closed),
        Err(Error::LebReaderClosed)
    ));
}

#[test]
fn matching_length_requires_digest_and_detects_changed_content() {
    let asset = crate::download::asset_by_name("base_core.leb2").unwrap();
    let bytes = core_bytes();
    let file = LebFile::from_bytes("base_core.leb2".into(), bytes.clone()).unwrap();
    assert!(asset.matches_image(&file).unwrap());
    assert_eq!(
        file.image.as_ref().unwrap().sha256.get().unwrap(),
        asset.sha256()
    );

    let mut changed = bytes;
    *changed.last_mut().unwrap() ^= 1;
    let other = LebFile::from_bytes("base_core.leb2".into(), changed).unwrap();
    assert_eq!(other.len(), file.len());
    assert!(!asset.matches_image(&other).unwrap());
    assert!(other.image.as_ref().unwrap().sha256.get().is_some());
}

#[test]
fn indexed_images_keep_first_hash_and_cached_identity_behavior() {
    use std::io::Write;

    let path = std::env::temp_dir().join(format!(
        "leb-manifest-{}-{}.leb2",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let mut output = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&path)
        .unwrap();
    output.write_all(&core_bytes()).unwrap();
    output.write_all(&[0]).unwrap();
    output.flush().unwrap();
    let file = LebFile::open_indexed(&path, IndexedLimits::default()).unwrap();
    let asset = crate::download::asset_by_name("base_core.leb2").unwrap();
    assert!(!asset.matches_image(&file).unwrap());
    assert!(file.image.as_ref().unwrap().sha256.get().is_some());

    // Deliberately violate the immutable-inode contract in a controlled test:
    // cached SHA keeps its existing semantics; explicit probing detects change.
    output.write_all(&[0]).unwrap();
    output.flush().unwrap();
    assert!(!asset.matches_image(&file).unwrap());
    assert!(file.validate_identity().is_err());
    drop(output);
    drop(file);
    std::fs::remove_file(path).unwrap();
}
