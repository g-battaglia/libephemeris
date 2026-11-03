//! Primary ayanamsha entries retain the public native observation clock.
use libephemeris_core::{
    api,
    context::{EphemerisContext, NetworkPolicy, ProviderMode},
    errors, native,
};
use std::error::Error as _;

#[test]
fn primary_mean_and_extended_entries_preserve_public_native_clock_bits() {
    for mode in [ProviderMode::SealedLeb, ProviderMode::NativeJpl] {
        for (selection, epoch, offset) in [
            (api::SIDM_FAGAN_BRADLEY, 0.0, 0.0),
            (api::SIDM_LAHIRI, 0.0, 0.0),
            (api::SIDM_USER, 2451545.0, 24.2),
            (
                api::SIDM_USER | libephemeris_core::constants::SIDBIT_USER_UT,
                2451545.0,
                24.2,
            ),
        ] {
            let mut context = EphemerisContext::new()
                .with_mode(mode)
                .with_user_delta_t(Some(0.125))
                .with_iers_delta_t(Some(false));
            let _ = context.set_sidereal_mode(selection, epoch, offset);
            let primary = context.primary_session_builder().build().unwrap();
            let original = native::Session::from_context(context);
            for day in [2415020.5, 2451545.0, 2460000.5] {
                for (tt, ut) in [(true, false), (false, true)] {
                    let expected = original.public_ayanamsha(day, 0, ut, false).unwrap();
                    let actual = if tt {
                        primary.get_ayanamsa(day)
                    } else {
                        primary.get_ayanamsa_ut(day)
                    }
                    .unwrap();
                    assert_eq!(
                        actual.value.to_bits(),
                        expected.value.to_bits(),
                        "mean clock differs: {mode:?}, {selection}, {day}, UT={ut}"
                    );
                    for flags in [0, api::FLG_SWIEPH | api::FLG_SPEED | api::FLG_NONUT] {
                        let expected = original.public_ayanamsha(day, flags, ut, true).unwrap();
                        let actual = if tt {
                            primary.get_ayanamsa_ex(day, flags)
                        } else {
                            primary.get_ayanamsa_ex_ut(day, flags)
                        }
                        .unwrap();
                        assert_eq!(
                            actual.value.to_bits(),
                            expected.value.to_bits(),
                            "extended clock differs: {mode:?}, {selection}, {day}, UT={ut}, {flags}"
                        );
                        assert_eq!(actual.returned_flags, expected.returned_flags);
                        assert_eq!(actual.diagnostic, expected.diagnostic);
                    }
                }
            }
        }
    }
}

/// Compare the full exposed result, or the genuine structured native failure.
/// Mean records deliberately expose no flags/diagnostic channel.
fn check_entry(
    primary: &api::Session,
    original: &native::Session,
    day: f64,
    is_ut: bool,
    flags: Option<i32>,
) -> bool {
    let (operation, actual) = match (is_ut, flags) {
        (false, None) => (
            "get_ayanamsa",
            primary.get_ayanamsa(day).map(|value| (value.value, None)),
        ),
        (true, None) => (
            "get_ayanamsa_ut",
            primary
                .get_ayanamsa_ut(day)
                .map(|value| (value.value, None)),
        ),
        (false, Some(flags)) => (
            "get_ayanamsa_ex",
            primary
                .get_ayanamsa_ex(day, flags)
                .map(|value| (value.value, Some((value.returned_flags, value.diagnostic)))),
        ),
        (true, Some(flags)) => (
            "get_ayanamsa_ex_ut",
            primary
                .get_ayanamsa_ex_ut(day, flags)
                .map(|value| (value.value, Some((value.returned_flags, value.diagnostic)))),
        ),
    };
    let expected = original.public_ayanamsha(day, flags.unwrap_or(0), is_ut, flags.is_some());
    match (actual, expected) {
        (Ok((value, channels)), Ok(expected)) => {
            assert_eq!(
                value.to_bits(),
                expected.value.to_bits(),
                "{operation}, {day}"
            );
            assert_eq!(
                channels,
                flags.map(|_| (expected.returned_flags, expected.diagnostic)),
                "{operation}, {day}"
            );
            true
        }
        (Err(actual), Err(expected)) => {
            assert_eq!(actual.operation(), Some(operation));
            let source = actual
                .source()
                .unwrap()
                .downcast_ref::<errors::Error>()
                .unwrap();
            // Debug retains structured fields; unlike PartialEq it handles NaN inputs.
            assert_eq!(format!("{source:?}"), format!("{expected:?}"));
            assert_eq!(actual.diagnostic(), Some(expected.to_string().as_str()));
            false
        }
        (actual, expected) => panic!("{operation}: {actual:?} differs from {expected:?}"),
    }
}

#[test]
fn primary_entries_keep_unknown_mode_diagnostics_and_terminal_errors() {
    for (mode, epoch, should_succeed) in [
        (-1, 0.0, true),
        (
            api::SIDM_USER | libephemeris_core::constants::SIDBIT_USER_UT,
            f64::NAN,
            false,
        ),
        (libephemeris_core::constants::SIDM_TRUE_CITRA, 0.0, false),
    ] {
        // Sealed, empty resources make a missing dynamic anchor an actual error,
        // never a download, static substitute or optional test skip.
        let mut context = EphemerisContext::new()
            .with_mode(ProviderMode::SealedLeb)
            .with_network_policy(Some(NetworkPolicy::Sealed))
            .with_iers_delta_t(Some(false));
        let _ = context.set_sidereal_mode(mode, epoch, 24.2);
        let primary = context.primary_session_builder().build().unwrap();
        let original = native::Session::from_context(context);
        for is_ut in [false, true] {
            for flags in [None, Some(0), Some(api::FLG_NONUT)] {
                assert_eq!(
                    check_entry(&primary, &original, 2451545.0, is_ut, flags),
                    should_succeed
                );
            }
        }
    }
}

#[test]
#[ignore = "requires explicit LIBEPHEMERIS_AYANAMSHA_DATA_ROOT with existing de440.bsp and leb/base_core.leb2"]
fn primary_entries_retain_real_spk_anchors_and_reader_clocks() {
    let root = std::path::PathBuf::from(
        std::env::var_os("LIBEPHEMERIS_AYANAMSHA_DATA_ROOT").expect("select existing local data"),
    );
    for mode in [ProviderMode::SealedLeb, ProviderMode::NativeJpl] {
        let mut context = EphemerisContext::new()
            .with_mode(mode)
            .with_network_policy(Some(NetworkPolicy::Sealed))
            .with_user_delta_t(Some(0.125))
            .with_iers_delta_t(Some(false));
        if mode == ProviderMode::SealedLeb {
            context.attach_leb_file(
                libephemeris_core::leb::LebFile::open(&root.join("leb/base_core.leb2")).unwrap(),
            );
        } else {
            let key = context
                .spk_registry_mut()
                .register_file(&root.join("de440.bsp"))
                .unwrap();
            libephemeris_core::spk::register_all_naif_targets(context.spk_registry_mut(), &key)
                .unwrap();
        }
        for selection in [
            api::SIDM_LAHIRI,
            api::SIDM_USER | libephemeris_core::constants::SIDBIT_ECL_DATE,
            libephemeris_core::constants::SIDM_TRUE_CITRA,
            libephemeris_core::constants::SIDM_ALDEBARAN_15TAU,
        ] {
            let _ = context.set_sidereal_mode(selection, 2451545.0, 24.2);
            let primary = context.primary_session_builder().build().unwrap();
            let original = native::Session::from_context(context.clone());
            for day in [2415020.5, 2451545.0, 2460000.5] {
                for is_ut in [false, true] {
                    for flags in [None, Some(0), Some(api::FLG_NONUT | api::FLG_TRUEPOS)] {
                        assert!(
                            check_entry(&primary, &original, day, is_ut, flags),
                            "{mode:?}, {selection}, {day}"
                        );
                    }
                }
            }
        }
    }
}
