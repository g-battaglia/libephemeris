//! Full public stored-TT fallback and observer-priority regressions.
use super::*;
use sha2::{Digest, Sha256};

/// Acquire the exact historical image with explicit clocks and optional site.
fn canonical_context(with_site: bool) -> EphemerisContext {
    let path = std::env::var_os("LIBEPHEMERIS_TEST_LEB2")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(crate::test_data::core_leb_path);
    let bytes = std::fs::read(&path).expect("provision the pinned historical fixture");
    let digest: String = Sha256::digest(&bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect();
    assert!(crate::test_data::admitted_core_hash(&(digest)));
    let mut ctx = EphemerisContext::new()
        .with_mode(ProviderMode::SealedLeb)
        .with_delta_t_model(Some(crate::time::DeltaTModel::EspenakMeeus))
        .with_iers_delta_t(Some(false))
        .with_tidal_acceleration(Some(f64::NAN))
        .with_network_policy(Some(crate::context::NetworkPolicy::Sealed));
    ctx.attach_leb_file(
        crate::leb::LebFile::from_bytes(path.display().to_string(), bytes).unwrap(),
    );
    if with_site {
        ctx.set_observer(Observer::new(11.0, 45.0, 100.0).unwrap());
    }
    ctx
}

#[test]
fn canonical_tt_fallback_public_session_matches_pinned_precision() {
    let result = Session::from_context(canonical_context(true))
        .calc(2396760.125, MARS, 33026)
        .unwrap();
    assert_eq!(result.returned_flags, 33026);
    let words = [
        0x40543d7be029daf2,
        0x400a9602cca064d9,
        0x3fe41d7c4baaeb04,
        0xbfd1c073bee7a000,
        0x3f8cade74114d000,
        0x3f714ee408e53c00,
    ];
    for (slot, (value, word)) in result.values.into_iter().zip(words).enumerate() {
        let reference = f64::from_bits(word);
        let limit = match slot {
            0 | 1 => 0.001 / 3600.0,
            2 => reference.abs() * 1e-8,
            3 | 4 if reference.abs() <= 1e-6 => 1e-6,
            _ => reference.abs() * 1e-6,
        };
        assert!(value.is_finite());
        assert!(
            (value - reference).abs() <= limit,
            "component {slot}: {value:?} vs pinned {reference:?}; physical limit {limit}"
        );
    }
}

#[test]
fn canonical_tt_missing_observer_precedes_calendar_and_table_alternative() {
    let error = Session::from_context(canonical_context(false))
        .calc(2396760.125, MARS, 33026)
        .unwrap_err();
    assert!(matches!(
        error,
        Error::ConfigurationMissing {
            missing_config: "observer_location",
            ..
        }
    ));
}
