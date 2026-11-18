//! Optional immutable base coefficients for local calculations without setup.

use super::LebFile;
use crate::errors::Result;
use std::io::Read;
use std::path::Path;
use std::sync::OnceLock;

// The data-only package retains the exact static bytes; the engine still owns
// parsing, validation, digest publication and the shared numerical reader.
use libephemeris_core_data::BASE_CORE_LEB2 as BASE;

/// Share the bundled 1850–2150 base coefficient image and decoded chunk cache.
///
/// Available with feature `bundled-data`. The complete container is validated
/// once and retains its own body-specific coverage limits; evaluation never
/// extrapolates outside those limits. Corrupt-container errors propagate.
/// Returned clones retain the immutable data independently of runtime resets.
pub fn bundled_base_core() -> Result<LebFile> {
    static READER: OnceLock<Result<LebFile>> = OnceLock::new();
    bundled_reader(&READER, || {
        LebFile::from_static_bytes("bundled/base_core.leb2".into(), BASE)
    })
}

pub(crate) fn bundled_base_core_with_process_cache(
    pool: std::sync::Arc<super::process_cache::ProcessCache>,
) -> Result<LebFile> {
    static READER: OnceLock<Result<LebFile>> = OnceLock::new();
    bundled_reader(&READER, || {
        LebFile::from_static_bytes("bundled/base_core.leb2".into(), BASE)
            .map(|file| file.with_process_cache(pool))
    })
}

/// Reuse the executable's immutable payload for an exact installed copy.
/// Retain the discovered pathname and a fresh reader/cache. An inconclusive
/// probe leaves ordinary owned acquisition and its errors to the caller.
pub(crate) fn exact_installed_base(path: &Path) -> Result<Option<LebFile>> {
    if path.file_name().and_then(|name| name.to_str()) != Some("base_core.leb2") {
        return Ok(None);
    }
    let Ok(mut file) = super::storage::open_regular_file(path) else {
        return Ok(None);
    };
    if file.metadata().map(|info| info.len()).ok() != Some(BASE.len() as u64) {
        return Ok(None);
    }
    let mut buffer = [0; 64 * 1024];
    for expected in BASE.chunks(buffer.len()) {
        let actual = &mut buffer[..expected.len()];
        if file.read_exact(actual).is_err() || actual != expected {
            return Ok(None);
        }
    }
    if !matches!(file.read(&mut buffer[..1]), Ok(0)) {
        return Ok(None);
    }
    let image = LebFile::from_static_bytes(path.to_string_lossy().into_owned(), BASE)?;
    pin_digest(&image);
    Ok(Some(image))
}

fn bundled_reader(
    cell: &'static OnceLock<Result<LebFile>>,
    open: impl FnOnce() -> Result<LebFile>,
) -> Result<LebFile> {
    cell.get_or_init(|| {
        let file = open()?;
        pin_digest(&file);
        Ok(file)
    })
    .clone()
}

fn pin_digest(file: &LebFile) {
    // The test below checks the executable's static input against this pin.
    file.image
        .as_ref()
        .expect("new static image")
        .sha256
        .set(crate::download::BASE_CORE_SHA256.to_owned())
        .expect("new image has no cached digest");
}

#[cfg(test)]
mod tests {
    use super::*;
    use sha2::{Digest, Sha256};

    #[test]
    fn discovered_base_reuses_only_an_exact_snapshot() {
        let dir = std::env::temp_dir().join(format!(
            "libephemeris-discovered-base-{}",
            std::process::id()
        ));
        std::fs::create_dir(&dir).unwrap();
        let path = dir.join("base_core.leb2");
        std::fs::write(&path, BASE).unwrap();

        let mut context = crate::context::EphemerisContext::new();
        let discovered = context.open_leb_file(&path).unwrap();
        assert_eq!(discovered.path(), path.to_string_lossy());
        assert_eq!(discovered.allocated_image_bytes(), 0);
        assert_eq!(discovered.bytes().unwrap().as_ptr(), BASE.as_ptr());

        context.set_leb_file(Some(path.clone()));
        let explicit = context.open_leb_file(&path).unwrap();
        assert_eq!(explicit.allocated_image_bytes(), BASE.len());
        assert_eq!(explicit.bytes().unwrap(), BASE);

        use std::io::Write;
        let mut changed = std::fs::OpenOptions::new().write(true).open(&path).unwrap();
        changed.write_all(b"BAD!").unwrap();
        drop(changed);
        assert!(exact_installed_base(&path).unwrap().is_none());
        context.set_leb_file(None);
        assert!(context.open_leb_file(&path).is_err());
        assert_eq!(discovered.bytes().unwrap(), BASE);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn bundled_coefficients_keep_their_published_identity_and_parse() {
        assert_eq!(BASE.len() as u64, crate::download::BASE_CORE_BYTES);
        assert_eq!(
            Sha256::digest(BASE)
                .iter()
                .map(|byte| format!("{byte:02x}"))
                .collect::<String>(),
            crate::download::BASE_CORE_SHA256
        );
        let first = bundled_base_core().unwrap();
        let second = bundled_base_core().unwrap();
        assert_eq!(
            first
                .image
                .as_ref()
                .unwrap()
                .sha256
                .get()
                .map(String::as_str),
            Some(crate::download::BASE_CORE_SHA256)
        );
        assert_eq!(first.sha256().unwrap(), crate::download::BASE_CORE_SHA256);
        assert_eq!(first.path, second.path);
        assert_eq!(first.allocated_image_bytes(), 0);
        assert_eq!(first.bytes().unwrap().as_ptr(), BASE.as_ptr());
        assert_eq!(second.bytes().unwrap().as_ptr(), BASE.as_ptr());
        let owned = LebFile::from_bytes("control".into(), BASE.to_vec()).unwrap();
        for jd in [2451545.0, 2460000.25, 2500000.5] {
            for body in [0, 1, 4, 9] {
                assert_eq!(
                    first.eval_body_any(body, jd).unwrap(),
                    owned.eval_body_any(body, jd).unwrap()
                );
            }
        }
    }
}
