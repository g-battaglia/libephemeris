use super::*;

/// Global-solar fixture over the sealed asset (no site: gamma is
/// geocentric; Espenak-Meeus Delta-T). `SKIP`s (never fails) when
/// the asset is absent.
macro_rules! bind_solglob23 {
    ($provider:ident, $nutation_at:ident, $serving:ident, $deltat_at:ident, $bound:ident) => {
        let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
        let bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(_) => {
                eprintln!("evt14_slice23_maxtime_tests: SKIP (missing {path})");
                return;
            }
        };
        let index = crate::leb::FileIndex::decode(&bytes).expect("core LEB decodes");
        let cfg = crate::time::DeltaTConfig {
            model: crate::time::DeltaTModel::EspenakMeeus,
            user_override_days: None,
            tidal_acceleration: crate::constants::TIDAL_DEFAULT,
            iers_enabled: false,
            iers_lookup_seconds: None,
            smh2016_lookup_seconds: None,
        };
        let $provider = crate::positions::SingleLeb::new(&index, &bytes);
        let $nutation_at = |jd: f64| index.eval_nutation(&bytes, jd);
        let $serving = crate::positions::SealedServing::new(&$provider, &$nutation_at);
        let $deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &cfg);
        let $bound =
            SealedSolarGlobalServing::try_new(&$serving, &$deltat_at, crate::constants::FLG_SWIEPH)
                .expect("global-solar serving binds");
    };
}

/// Site-bound solar fixture over the sealed asset (Dallas: the
/// banked 2024 local-maximum site). `SKIP`s when absent.
macro_rules! bind_solloc23 {
    ($provider:ident, $nutation_at:ident, $serving:ident, $deltat_at:ident, $bound:ident) => {
        let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
        let bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(_) => {
                eprintln!("evt14_slice23_maxtime_tests: SKIP (missing {path})");
                return;
            }
        };
        let index = crate::leb::FileIndex::decode(&bytes).expect("core LEB decodes");
        let cfg = crate::time::DeltaTConfig {
            model: crate::time::DeltaTModel::EspenakMeeus,
            user_override_days: None,
            tidal_acceleration: crate::constants::TIDAL_DEFAULT,
            iers_enabled: false,
            iers_lookup_seconds: None,
            smh2016_lookup_seconds: None,
        };
        let $provider = crate::positions::SingleLeb::new(&index, &bytes);
        let $nutation_at = |jd: f64| index.eval_nutation(&bytes, jd);
        let $serving = crate::positions::SealedServing::new(&$provider, &$nutation_at);
        let $deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &cfg);
        let $bound = SealedSolarLocalServing::try_new(
            &$serving,
            &$deltat_at,
            crate::constants::FLG_SWIEPH,
            -96.7970,
            32.7767,
            131.0,
        )
        .expect("Dallas serving binds");
    };
}

#[test]
fn saros_inex_numbers_pin_banked_2024_maximum() {
    // Release documents solar Saros 139 / Inex 50 for 2024-04-08.
    assert_eq!(saros_number(2460409.2619832945, EclipseKind::Solar), 139);
    assert_eq!(inex_number(2460409.2619832945, EclipseKind::Solar), 50);
    assert!(parse_eclipse_kind("SOLAR").is_ok());
    assert!(parse_eclipse_kind("comet").is_err());
}

#[test]
fn sol_global_max_time_recovers_banked_maximum() {
    bind_solglob23!(provider, nutation_at, serving, deltat_at, bound);
    // Banked 2024-04-08 global maximum (slice-22 release word),
    // approached from three hours early with the release ±0.125 d
    // window: the golden search must land back on the maximum with
    // a sub-lunar-disc gamma.
    let (jd_max, gamma) = sol_global_max_time_sealed(&bound, 2460409.2619832945 - 0.125, 0.125)
        .expect("global max solves");
    assert!(
        (jd_max - 2460409.2619832945).abs() < 2e-6,
        "global max returns to the banked instant, got {jd_max}"
    );
    assert!(
        gamma < 1.0,
        "central-eclipse gamma pierces inside the disc, got {gamma}"
    );
}

#[test]
fn sol_local_max_time_recovers_banked_dallas_maximum() {
    bind_solloc23!(provider, nutation_at, serving, deltat_at, bound);
    // Banked 2024-04-08 Dallas local maximum (solhow release word),
    // approached from three hours early: the separation minimum
    // must land back on the local maximum inside totality.
    let (jd_max, sep) = sol_local_max_time_sealed(&bound, 2460409.279538682 - 0.125, 0.125)
        .expect("local max solves");
    assert!(
        (jd_max - 2460409.279538682).abs() < 2e-6,
        "local max returns to the banked instant, got {jd_max}"
    );
    assert!(
        sep < 0.1,
        "totality separation is a small fraction of a degree, got {sep}"
    );
}

#[test]
fn sol_max_time_refuses_bad_inputs() {
    bind_solglob23!(provider, nutation_at, serving, deltat_at, bound);
    assert!(sol_global_max_time_sealed(&bound, f64::NAN, 0.125).is_err());
    assert!(sol_global_max_time_sealed(&bound, 2460409.26, 0.0).is_err());
    assert!(sol_global_max_time_sealed(&bound, 2460409.26, -0.125).is_err());
}
