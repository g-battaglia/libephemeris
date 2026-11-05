//! Resolution precedence and construction effects, without provider acquisition.

use super::*;

fn sources(
    environment: &[(&str, &str)],
    dotenv: &[(&str, &str)],
    table: &str,
) -> ConfigurationSources {
    ConfigurationSources {
        environment: environment
            .iter()
            .map(|(key, value)| ((*key).into(), (*value).into()))
            .collect(),
        dotenv: dotenv
            .iter()
            .map(|(key, value)| ((*key).into(), (*value).into()))
            .collect(),
        table: table.parse().unwrap(),
        ..ConfigurationSources::default()
    }
}

#[test]
fn mode_and_tier_keep_explicit_captured_and_default_precedence() {
    let environment = [
        ("LIBEPHEMERIS_MODE", "leb"),
        ("LIBEPHEMERIS_PRECISION", "base"),
    ];
    let dotenv = [
        ("LIBEPHEMERIS_MODE", "skyfield"),
        ("LIBEPHEMERIS_PRECISION", "extended"),
    ];
    let table = "mode='horizons'\nprecision='medium'\n";
    for (env, dot, toml, mode, tier, expected_mode, expected_tier) in [
        (
            &environment[..],
            &dotenv[..],
            table,
            Some(ProviderMode::Auto),
            Some(Tier::Extended),
            ProviderMode::Auto,
            Tier::Extended,
        ),
        (
            &environment[..],
            &dotenv[..],
            table,
            None,
            None,
            ProviderMode::SealedLeb,
            Tier::Base,
        ),
        (
            &[][..],
            &dotenv[..],
            table,
            None,
            None,
            ProviderMode::NativeJpl,
            Tier::Extended,
        ),
        (
            &[][..],
            &[][..],
            table,
            None,
            None,
            ProviderMode::HorizonsRemote,
            Tier::Medium,
        ),
        (
            &[][..],
            &[][..],
            "",
            None,
            None,
            ProviderMode::Auto,
            Tier::Medium,
        ),
    ] {
        let resolved = sources(env, dot, toml)
            .resolve(&EphemerisContext::new(), mode, tier)
            .unwrap();
        assert_eq!(resolved.mode(), expected_mode);
        assert_eq!(resolved.tier(), expected_tier);
    }
}

#[test]
fn tier_application_restores_selected_filename_without_clearing_leb_override() {
    let sources = sources(&[("LIBEPHEMERIS_EPHEMERIS", "de441.bsp")], &[], "");
    let mut explicit = EphemerisContext::new().with_tier(Tier::Base);
    explicit.set_leb_file(Some(PathBuf::from("/not-opened/custom.leb2")));
    explicit.set_ephemeris_file(Some("custom.bsp"));
    let resolved = sources.resolve(&explicit, None, Some(Tier::Base)).unwrap();
    assert_eq!(resolved.ephemeris_file_override(), Some("custom.bsp"));
    assert_eq!(resolved.leb_file(), explicit.leb_file());
    assert_eq!(explicit.ephemeris_file_override(), Some("custom.bsp"));
    explicit.set_ephemeris_file(None);
    let resolved = sources.resolve(&explicit, None, Some(Tier::Base)).unwrap();
    assert_eq!(resolved.ephemeris_file_override(), Some("de441.bsp"));
    assert_eq!(resolved.leb_file(), explicit.leb_file());
}

#[test]
fn iers_directory_derivation_configured_override_and_explicit_choice_stay_distinct() {
    let mut explicit = EphemerisContext::new();
    explicit.set_jpl_data_dir(Some(PathBuf::from("/not-created/data")));
    let plain = ConfigurationSources::default()
        .resolve(&explicit, None, None)
        .unwrap();
    assert_eq!(
        plain.iers_cache_dir(),
        Some(&PathBuf::from("/not-created/data/iers_cache"))
    );
    let configured = sources(
        &[("LIBEPHEMERIS_IERS_CACHE_DIR", "/not-created/configured")],
        &[],
        "",
    );
    let selected = configured.resolve(&explicit, None, None).unwrap();
    assert_eq!(
        selected.iers_cache_dir(),
        Some(&PathBuf::from("/not-created/configured"))
    );
    explicit.set_iers_cache_dir(Some(PathBuf::from("/not-created/explicit")));
    let selected = configured.resolve(&explicit, None, None).unwrap();
    assert_eq!(selected.iers_cache_dir(), explicit.iers_cache_dir());
    assert_eq!(selected.jpl_data_dir(), explicit.jpl_data_dir());
}

#[test]
fn independent_boolean_overrides_keep_each_captured_layer() {
    let environment = [
        ("LIBEPHEMERIS_IERS_DELTA_T", "1"),
        ("LIBEPHEMERIS_AUTO_SPK", "1"),
        ("LIBEPHEMERIS_STRICT_PRECISION", "1"),
        ("LIBEPHEMERIS_IERS_AUTO_DOWNLOAD", "1"),
    ];
    let dotenv = [
        ("LIBEPHEMERIS_IERS_DELTA_T", "0"),
        ("LIBEPHEMERIS_AUTO_SPK", "0"),
        ("LIBEPHEMERIS_STRICT_PRECISION", "0"),
        ("LIBEPHEMERIS_IERS_AUTO_DOWNLOAD", "0"),
    ];
    let table =
        "iers_delta_t=true\nauto_spk=true\nstrict_precision=true\niers_auto_download=true\n";
    for (level, expected) in [(0, false), (1, true), (2, false), (3, true)] {
        let captured = sources(
            if level <= 1 { &environment } else { &[] },
            if level <= 2 { &dotenv } else { &[] },
            table,
        );
        let value = (level == 0).then_some(false);
        let explicit = EphemerisContext::new()
            .with_iers_delta_t(value)
            .with_auto_spk_download(value)
            .with_strict_precision(value)
            .with_iers_auto_download(value);
        let resolved = captured.resolve(&explicit, None, None).unwrap();
        assert_eq!(resolved.iers_delta_t_override(), Some(expected));
        assert_eq!(resolved.auto_spk_download_override(), Some(expected));
        assert_eq!(resolved.strict_precision_override(), Some(expected));
        assert_eq!(resolved.iers_auto_download_override(), Some(expected));
    }
}

#[test]
fn resolution_keeps_the_seed_observer_and_attached_image_owner() {
    let image = crate::leb::LebFile::from_bytes(
        "attached-base.leb2".into(),
        crate::test_data::core_bytes().to_vec(),
    )
    .unwrap();
    let mut explicit = EphemerisContext::new();
    explicit.set_leb_file(Some(PathBuf::from("attached-base.leb2")));
    explicit.attach_leb_file(image);
    explicit.set_observer(crate::context::Observer::new(12.5, 41.9, 30.0).unwrap());
    let pointer = explicit.leb_files()[0].bytes().unwrap().as_ptr();
    let resolved = ConfigurationSources::default()
        .resolve(&explicit, None, None)
        .unwrap();
    assert_eq!(resolved.leb_files()[0].bytes().unwrap().as_ptr(), pointer);
    let observer = resolved.observer().unwrap();
    assert_eq!(
        (observer.lon_deg(), observer.lat_deg(), observer.height_m()),
        (12.5, 41.9, 30.0)
    );
    assert_eq!(explicit.leb_files()[0].bytes().unwrap().as_ptr(), pointer);
}

#[test]
fn invalid_reader_mode_and_network_keep_consumed_error_precedence() {
    let mut sources = sources(
        &[
            ("LIBEPHEMERIS_LEB_READER", "invalid-reader"),
            ("LIBEPHEMERIS_MODE", "invalid-mode"),
            ("LIBEPHEMERIS_NETWORK_POLICY", "invalid-policy"),
        ],
        &[],
        "",
    );
    let explicit = EphemerisContext::new();
    assert!(matches!(
        sources.resolve(&explicit, None, None),
        Err(crate::errors::Error::InvalidLebReader { .. })
    ));
    sources.environment.remove("LIBEPHEMERIS_LEB_READER");
    assert!(matches!(
        sources.resolve(&explicit, None, None),
        Err(crate::errors::Error::InvalidCalcMode { .. })
    ));
    sources.environment.remove("LIBEPHEMERIS_MODE");
    assert!(matches!(
        sources.resolve(&explicit, None, None),
        Err(crate::errors::Error::InvalidNetworkPolicy { .. })
    ));
}
