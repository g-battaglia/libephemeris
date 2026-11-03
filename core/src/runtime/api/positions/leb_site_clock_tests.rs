//! Stored-state site-clock regression and lazy acquisition contracts.
use super::*;
use sha2::{Digest, Sha256};

/// An owned immutable historical core with an explicit configured model.
fn context(with_site: bool) -> EphemerisContext {
    let path = std::env::var_os("LIBEPHEMERIS_TEST_LEB2")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(crate::test_data::core_leb_path);
    let bytes = std::fs::read(&path).expect("provision the pinned historical fixture");
    let digest: String = Sha256::digest(&bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    assert!(crate::test_data::admitted_core_hash(&(digest)));
    let mut ctx = EphemerisContext::new()
        .with_mode(ProviderMode::SealedLeb)
        .with_delta_t_model(Some(crate::time::DeltaTModel::EspenakMeeus))
        .with_iers_delta_t(Some(false))
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
fn december_tt_and_original_ut_match_source_within_precision_contract() {
    let session = Session::from_context(context(true));
    let flags = FLG_SWIEPH | FLG_TOPOCTR | FLG_SPEED;
    for (is_ut, expected) in [
        (
            true,
            [
                0x406e7a0481c466eb,
                0xbfc89a10d1b1e93e,
                0x400416b0ba5e568a,
                0x3fe6ecd8c2c50000,
                0xbf83b1f37654ba80,
                0xbf62c0187b044000,
            ],
        ),
        (
            false,
            [
                0x406e79ee33ace23e,
                0xbfc898dedbcc8d0e,
                0x400416b548d7c3ac,
                0x3fe6ecf566771000,
                0xbf839cc2b8ca5a00,
                0xbf62b75a0174c000,
            ],
        ),
    ] {
        let result = if is_ut {
            session.calc_ut(2506300.5, MARS, flags)
        } else {
            session.calc(2506300.5, MARS, flags)
        }
        .unwrap();
        assert_eq!(result.returned_flags, flags);
        for (slot, (actual, reference)) in result.values.into_iter().zip(expected).enumerate() {
            let reference = f64::from_bits(reference);
            let tolerance = match slot {
                0 | 1 => 0.001 / 3600.0,
                2 => reference.abs() * 1e-8,
                _ => reference.abs() * 1e-6,
            };
            assert!(
                (actual - reference).abs() <= tolerance,
                "UT entry {is_ut}, slot {slot}: {actual} vs {reference} (limit {tolerance})"
            );
        }
    }
}

#[test]
fn december_inverse_consumes_the_previous_month_delta_t() {
    let session = Session::from_context(context(true));
    let tt = 2506300.5;
    let delta = |ut| session.deltat(ut);
    let partner = fast_ut1_partner(tt, tt, false, &delta).unwrap();
    assert_eq!(partner.to_bits(), 0x41431f1e3f838c36);
    assert_ne!(partner.to_bits(), (tt - delta(tt).unwrap()).to_bits());
}

#[test]
fn missing_observer_precedes_site_clock_and_provider_acquisition() {
    let ctx = context(false);
    serving::with_sealed_serving_stack(
        &ctx,
        crate::time::DeltaTModel::EspenakMeeus,
        false,
        true,
        |stack| {
            let unused = |_: f64| -> Result<f64, Error> { panic!("unused site clock") };
            assert!(matches!(
                stored_state_with_clock(
                    stack,
                    &ctx,
                    2451545.125,
                    MARS,
                    FLG_TOPOCTR,
                    false,
                    &unused
                ),
                Err(Error::ConfigurationMissing {
                    missing_config: "observer_location",
                    ..
                })
            ));
            Ok(())
        },
    )
    .unwrap();
}

#[test]
fn geocentric_and_noop_requests_acquire_no_site_clock() {
    let ctx = context(false);
    serving::with_sealed_serving_stack(
        &ctx,
        crate::time::DeltaTModel::EspenakMeeus,
        false,
        true,
        |stack| {
            let unused = |_: f64| -> Result<f64, Error> { panic!("unused site clock") };
            stored_state_with_clock(stack, &ctx, 2451545.125, MARS, FLG_SWIEPH, false, &unused)?;
            tropical(
                stack,
                &ctx,
                2451545.125,
                MEAN_NODE,
                FLG_TOPOCTR,
                &unused,
                2451545.125,
            )?;
            Ok(())
        },
    )
    .unwrap();
}

#[test]
fn site_clock_failure_is_terminal_before_provider_reads() {
    let ctx = context(true);
    let marker = Error::ConfigurationMissing {
        missing_config: "test_site_clock",
        suggestion: "retain the actual clock failure",
    };
    serving::with_sealed_serving_stack(
        &ctx,
        crate::time::DeltaTModel::EspenakMeeus,
        false,
        true,
        |stack| {
            let clock = |_: f64| Err(marker.clone());
            assert_eq!(
                stored_state_with_clock(stack, &ctx, f64::MAX, MARS, FLG_TOPOCTR, false, &clock),
                Err(marker.clone())
            );
            Ok(())
        },
    )
    .unwrap();
}

#[test]
fn public_missing_observer_precedes_nonfinite_configured_clock() {
    let session = Session::from_context(context(false).with_user_delta_t(Some(f64::NAN)));
    let flags = FLG_SWIEPH | FLG_TOPOCTR | FLG_SPEED;
    for result in [
        session.calc_ut(2506300.5, MARS, flags),
        session.calc(2506300.5, MARS, flags),
    ] {
        assert!(matches!(
            result,
            Err(Error::ConfigurationMissing {
                missing_config: "observer_location",
                ..
            })
        ));
    }
}

#[test]
fn public_raw_flag_conflict_precedes_nonfinite_configured_clock() {
    let session = Session::from_context(context(false).with_user_delta_t(Some(f64::NAN)));
    let flags = FLG_BARYCTR | FLG_MOSEPH;
    for (entry, result) in [
        ("calc_ut", session.calc_ut(2506300.5, MARS, flags)),
        ("calc", session.calc(2506300.5, MARS, flags)),
    ] {
        assert!(matches!(
            result,
            Err(Error::UnsupportedFlagCombination { entry: actual, flags: actual_flags })
                if actual == entry && actual_flags == flags
        ));
    }
}

#[test]
fn public_early_domains_do_not_require_an_observer() {
    let session = Session::from_context(context(false));
    for is_ut in [false, true] {
        let nutation = if is_ut {
            session.calc_ut(2506300.5, ECL_NUT, FLG_TOPOCTR)
        } else {
            session.calc(2506300.5, ECL_NUT, FLG_TOPOCTR)
        }
        .expect("nutation has no observer acquisition");
        assert!(nutation.values.iter().all(|value| value.is_finite()));
        let flags = FLG_SWIEPH | FLG_TOPOCTR | FLG_SPEED;
        let earth = if is_ut {
            session.calc_ut(2506300.5, EARTH, flags)
        } else {
            session.calc(2506300.5, EARTH, flags)
        }
        .expect("topocentric Earth is the exact zero origin");
        assert_eq!(earth.values.map(f64::to_bits), [0; 6]);
        assert_eq!(earth.returned_flags, flags);
    }
}

#[test]
fn public_reader_acquisition_precedes_ordinary_observer_and_clock() {
    let ctx = EphemerisContext::new()
        .with_mode(ProviderMode::SealedLeb)
        .with_delta_t_model(Some(crate::time::DeltaTModel::EspenakMeeus))
        .with_iers_delta_t(Some(false))
        .with_user_delta_t(Some(f64::NAN));
    let session = Session::from_context(ctx);
    let flags = FLG_SWIEPH | FLG_TOPOCTR | FLG_SPEED;
    for result in [
        session.calc_ut(2506300.5, MARS, flags),
        session.calc(2506300.5, MARS, flags),
    ] {
        assert!(matches!(
            result,
            Err(Error::ConfigurationMissing {
                missing_config: "leb_reader",
                ..
            })
        ));
    }
}
