use super::*;

/// Global-lunar fixture over the sealed asset (same shape as the
/// slice-5b fixture: sealed core, Espenak-Meeus Delta-T).
macro_rules! bind_lunglob17 {
    ($bytes:ident, $index:ident, $cfg:ident, $provider:ident, $nutation_at:ident,
         $serving:ident, $deltat_at:ident, $bound:ident) => {
        let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
        let $bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(_) => {
                eprintln!("evt14_slice17_tests: SKIP (missing {path})");
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
        let $bound =
            SealedLunarGlobalServing::try_new(&$serving, &$deltat_at, crate::constants::FLG_SWIEPH)
                .expect("global-lunar serving binds");
    };
}

#[test]
fn non_finite_maxima_refuse_on_all_contact_arms() {
    bind_lunglob17!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        bound
    );
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            lun_contact_p1_sealed(&bound, bad),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            lun_contact_p4_sealed(&bound, bad),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            lun_contact_u1_sealed(&bound, bad),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            lun_contact_u2_sealed(&bound, bad),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            lun_contact_u3_sealed(&bound, bad),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            lun_contact_u4_sealed(&bound, bad),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            lun_total_duration_sealed(&bound, bad),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            lun_umbral_duration_sealed(&bound, bad),
            Err(Error::NonFiniteArgument { .. })
        ));
    }
}

#[test]
fn off_eclipse_maximum_returns_zero_words() {
    bind_lunglob17!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        bound
    );
    // Off-eclipse JD: every contact/duration word is 0.0 (the
    // 32-case sealed sweep agrees word-for-word with the release
    // here; this pins the no-eclipse shape, not a golden value).
    let jd_off = 2459959.5;
    assert_eq!(lun_contact_p1_sealed(&bound, jd_off).unwrap(), 0.0);
    assert_eq!(lun_contact_u1_sealed(&bound, jd_off).unwrap(), 0.0);
    assert_eq!(lun_total_duration_sealed(&bound, jd_off).unwrap(), 0.0);
    assert_eq!(lun_umbral_duration_sealed(&bound, jd_off).unwrap(), 0.0);
}

#[test]
fn total_eclipse_maximum_resolves_totality() {
    bind_lunglob17!(
        bytes,
        index,
        cfg,
        provider,
        nutation_at,
        serving,
        deltat_at,
        bound
    );
    // Banked 2022-11-08 total maximum: totality resolves with a
    // positive duration in the expected ~85-minute neighborhood
    // .
    let jd_total = 2459891.957725105;
    let total = lun_total_duration_sealed(&bound, jd_total).unwrap();
    assert!(total > 80.0 && total < 90.0, "total {total}");
    let u2 = lun_contact_u2_sealed(&bound, jd_total).unwrap();
    let u3 = lun_contact_u3_sealed(&bound, jd_total).unwrap();
    assert!(u2 > 0.0 && u3 > u2, "u2 {u2} u3 {u3}");
}
