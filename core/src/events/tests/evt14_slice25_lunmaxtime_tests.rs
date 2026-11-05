use super::*;

/// Global-lunar fixture over the sealed asset (same shape as the
/// slice-17 fixture: sealed core, Espenak-Meeus Delta-T).
macro_rules! bind_lunglob25 {
    ($provider:ident, $nutation_at:ident, $serving:ident, $deltat_at:ident, $bound:ident) => {
        let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
        let bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(_) => {
                eprintln!("evt14_slice25_lunmaxtime_tests: SKIP (missing {path})");
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
            SealedLunarGlobalServing::try_new(&$serving, &$deltat_at, crate::constants::FLG_SWIEPH)
                .expect("global-lunar serving binds");
    };
}

#[test]
fn lun_max_time_recovers_banked_total_maximum() {
    bind_lunglob25!(provider, nutation_at, serving, deltat_at, bound);
    // Banked 2022-11-08 total maximum (luncontact release word),
    // approached from three hours early: the depth minimum must
    // land back on the maximum inside the contact-family ceiling.
    let jd_max =
        lun_eclipse_max_time_sealed(&bound, 2459891.957725105 - 0.125).expect("lunar max solves");
    assert!(
        (jd_max - 2459891.957725105).abs() < 2e-6,
        "lunar max returns to the banked instant, got {jd_max}"
    );
}

#[test]
fn lun_max_time_off_eclipse_stays_finite() {
    bind_lunglob25!(provider, nutation_at, serving, deltat_at, bound);
    // Quiet JD: no eclipse anywhere, but the fixed-window search
    // still minimizes to a finite instant (never a fabricated
    // event word).
    let jd_max = lun_eclipse_max_time_sealed(&bound, 2459959.5).expect("quiet max solves");
    assert!(jd_max.is_finite(), "quiet max is finite, got {jd_max}");
    assert!(
        (jd_max - 2459959.5).abs() <= 0.3,
        "quiet max stays inside the release window, got {jd_max}"
    );
}

#[test]
fn lun_max_time_refuses_non_finite() {
    bind_lunglob25!(provider, nutation_at, serving, deltat_at, bound);
    assert!(lun_eclipse_max_time_sealed(&bound, f64::NAN).is_err());
    assert!(lun_eclipse_max_time_sealed(&bound, f64::INFINITY).is_err());
}
