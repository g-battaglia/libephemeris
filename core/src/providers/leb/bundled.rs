//! Optional immutable base coefficients for local calculations without setup.

use super::LebFile;
use crate::errors::Result;
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
    READER
        .get_or_init(|| {
            let file = LebFile::from_static_bytes("bundled/base_core.leb2".into(), BASE)?;
            // These are the executable's exact static bytes, not a pathname or
            // filesystem-metadata claim. The test below verifies the build input
            // against this pin; no caller can seed a digest for another image.
            file.image
                .as_ref()
                .expect("new static image")
                .sha256
                .set(crate::download::BASE_CORE_SHA256.to_owned())
                .expect("new image has no cached digest");
            Ok(file)
        })
        .clone()
}

#[cfg(test)]
mod tests {
    use super::*;
    use sha2::{Digest, Sha256};

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
