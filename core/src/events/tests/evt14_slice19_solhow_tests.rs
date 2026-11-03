use super::*;

/// Local-solar fixture over the sealed asset at one observer triple
/// (sealed core, Espenak-Meeus Delta-T).
macro_rules! bind_solloc19 {
    ($bytes:ident, $index:ident, $cfg:ident, $provider:ident, $nutation_at:ident,
         $serving:ident, $deltat_at:ident, $bound:ident, $lon:expr, $lat:expr, $alt:expr) => {
        let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
        let $bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(_) => {
                eprintln!("evt14_slice19_tests: SKIP (missing {path})");
                return;
            }
        };
        let $index = crate::leb::FileIndex::decode(&$bytes).expect("core LEB decodes");
        let $cfg = crate::time::DeltaTConfig {
            model: crate::time::DeltaTModel::EspenakMeeus,
            user_override_days: None,
            tidal_acceleration: crate::constants::TIDAL_DEFAULT,
            iers_enabled: false,
            iers_lookup_seconds: None,
            smh2016_lookup_seconds: None,
        };
        let $provider = crate::positions::SingleLeb::new(&$index, &$bytes);
        let $nutation_at = |jd: f64| $index.eval_nutation(&$bytes, jd);
        let $serving = crate::positions::SealedServing::new(&$provider, &$nutation_at);
        let $deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &$cfg);
        let $bound = SealedSolarLocalServing::try_new(
            &$serving,
            &$deltat_at,
            crate::constants::FLG_SWIEPH,
            $lon,
            $lat,
            $alt,
        )
        .expect("local-solar serving binds");
    };
}

#[test]
fn non_finite_instants_refuse() {
    bind_solloc19!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        bound,
        -96.7970,
        32.7767,
        131.0
    );
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            sol_how_sealed(&bound, bad),
            Err(Error::NonFiniteArgument { .. })
        ));
    }
}

#[test]
fn total_maximum_carries_type_and_observer_independent_core_width() {
    bind_solloc19!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        dallas,
        -96.7970,
        32.7767,
        131.0
    );
    bind_solloc19!(
        bytes2,
        index2,
        cfg2,
        provider2,
        nutation_at2,
        serving2,
        deltat_at2,
        denver,
        -104.9903,
        39.7392,
        1609.0
    );
    let _ = (
        &bytes2,
        &index2,
        &cfg2,
        &provider2,
        &nutation_at2,
        &serving2,
        &deltat_at2,
    );
    // Banked 2024-04-08 Dallas local maximum: total character with
    // the visible bit (midday Sun), magnitude at or above unity,
    // finite azalt triple.
    let jd_total = 2460409.279538682;
    let (retflag, attr) = sol_how_sealed(&dallas, jd_total).unwrap();
    assert_ne!(retflag & (ECL_TOTAL | ECL_PARTIAL | ECL_ANNULAR), 0);
    assert_ne!(retflag & ECL_VISIBLE, 0);
    assert!(attr[0] >= 1.0, "total magnitude {}", attr[0]);
    assert!(attr[4].is_finite() && attr[5].is_finite() && attr[6].is_finite());
    assert!(attr[6] > 0.0);
    // The core-shadow diameter at the global shadow-center point is
    // observer-independent: Denver sees the same word bitwise.
    let (_, denver_attr) = sol_how_sealed(&denver, jd_total).unwrap();
    assert_eq!(attr[3], denver_attr[3]);
}

#[test]
fn night_side_and_quiet_instants_zero_the_eclipse_slots() {
    bind_solloc19!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        tokyo,
        139.6917,
        35.6895,
        40.0
    );
    // Banked 2024-04-08 total maximum is night in Tokyo: retflag 0
    // with the eclipse slots zeroed while azalt/separation stay.
    let (retflag, attr) = sol_how_sealed(&tokyo, 2460409.279538682).unwrap();
    assert_eq!(retflag, 0);
    for i in [0, 1, 2, 3, 8, 9, 10] {
        assert_eq!(attr[i], 0.0, "slot {i}");
    }
    assert!(attr[4].is_finite() && attr[7].is_finite());
    // Quiet JD between seasons: no eclipse from Dallas either.
    bind_solloc19!(
        bytes2,
        index2,
        cfg2,
        provider2,
        nutation_at2,
        serving2,
        deltat_at2,
        dallas,
        -96.7970,
        32.7767,
        131.0
    );
    let _ = (
        &bytes2,
        &index2,
        &cfg2,
        &provider2,
        &nutation_at2,
        &serving2,
        &deltat_at2,
    );
    let (quiet_flag, quiet_attr) = sol_how_sealed(&dallas, 2460140.0).unwrap();
    assert_eq!(quiet_flag, 0);
    for i in [0, 1, 2, 3, 8, 9, 10] {
        assert_eq!(quiet_attr[i], 0.0, "quiet slot {i}");
    }
}
