//! Network-free lifecycle controls for a valid LEB container lacking the Sun.
//! The default fixture has no coefficient records: these tests exercise missing
//! data, ownership and error propagation, without astronomical expected vectors.

use libephemeris_core::constants::*;
use libephemeris_core::context::{EphemerisContext, NetworkPolicy, ProviderMode};
use libephemeris_core::errors::Error;
use libephemeris_core::leb::LebFile;
use libephemeris_core::native::Session;
use sha2::{Digest, Sha256};
use std::sync::Arc;

const MEDIUM_SHA256: &str = "1c63b9cc2d53421d2ad0cc931c45b721810e11f51fe277b1a179efe0b6e80b20";
const FLAGS: i32 = FLG_SWIEPH | FLG_SPEED | FLG_TRUEPOS | FLG_BARYCTR | FLG_XYZ;

fn context(mode: ProviderMode) -> EphemerisContext {
    let (path, bytes) = if let Some(path) = std::env::var_os("LIBEPHEMERIS_TEST_MEDIUM_ASTEROIDS") {
        let path = std::path::PathBuf::from(path);
        let bytes = std::fs::read(&path).expect("read the supplied medium asteroid fixture");
        assert_eq!(bytes.len(), 6_476_624);
        assert_eq!(
            Sha256::digest(&bytes)
                .iter()
                .map(|byte| format!("{byte:02x}"))
                .collect::<String>(),
            MEDIUM_SHA256
        );
        (path.display().to_string(), bytes)
    } else {
        // LEB2 v1: one empty body-index section, covering the test epochs.
        // No coefficients are generated, substituted or evaluated.
        let mut bytes = vec![0u8; 88];
        bytes[..4].copy_from_slice(b"LEB2");
        bytes[4..8].copy_from_slice(&1u32.to_le_bytes());
        bytes[8..12].copy_from_slice(&1u32.to_le_bytes());
        bytes[16..24].copy_from_slice(&2_400_000.0f64.to_le_bytes());
        bytes[24..32].copy_from_slice(&2_500_000.0f64.to_le_bytes());
        bytes[40..44].copy_from_slice(&1u32.to_le_bytes());
        bytes[72..80].copy_from_slice(&88u64.to_le_bytes());
        ("medium_asteroids.leb2".to_string(), bytes)
    };
    let mut ctx = EphemerisContext::new()
        .with_mode(mode)
        .with_network_policy(Some(NetworkPolicy::Sealed))
        .with_delta_t_model(Some(libephemeris_core::time::DeltaTModel::Smh2016))
        .with_iers_delta_t(Some(false));
    ctx.attach_leb_file(LebFile::from_bytes(path, bytes).unwrap());
    ctx
}

#[test]
fn missing_sealed_core_body_is_a_native_error() {
    let session = Session::from_context(context(ProviderMode::SealedLeb));
    assert!(matches!(
        session
            .calc(2451545.0, SUN, FLAGS)
            .unwrap_err()
            .primary_error(),
        Error::UnknownBody { body_id: SUN }
    ));
}

#[test]
fn automatic_provider_cannot_download_under_sealed_network_policy() {
    let session = Session::from_context(context(ProviderMode::Auto));
    let error = session.calc(2451545.0, SUN, FLAGS).unwrap_err();
    assert!(
        matches!(error.primary_error(), Error::NetworkSealed { .. }),
        "{error:?}"
    );
    assert!(session.context().spk_registry().is_empty());
}

#[test]
fn closing_one_context_does_not_retire_another_reader_owner() {
    let mut original = context(ProviderMode::SealedLeb);
    let retained = original.clone();
    original.close();
    assert!(original.leb_files().is_empty());
    assert_eq!(retained.leb_files().len(), 1);
    let error = Session::from_context(retained)
        .calc(2451545.0, SUN, FLAGS)
        .unwrap_err();
    assert!(matches!(
        error.primary_error(),
        Error::UnknownBody { body_id: SUN }
    ));
}

#[test]
fn simultaneous_native_errors_remain_independent() {
    let session = Arc::new(Session::from_context(context(ProviderMode::SealedLeb)));
    let threads: Vec<_> = (0..8)
        .map(|i| {
            let session = Arc::clone(&session);
            std::thread::spawn(move || {
                session
                    .calc(2451545.0 + f64::from(i), SUN, FLAGS)
                    .unwrap_err()
            })
        })
        .collect();
    for thread in threads {
        assert!(matches!(
            thread.join().unwrap().primary_error(),
            Error::UnknownBody { body_id: SUN }
        ));
    }
    assert_eq!(session.context().leb_files().len(), 1);
}
