use super::*;

/// Local-solar fixture over the sealed asset at one observer triple
/// (sealed core, Espenak-Meeus Delta-T).
/// `SKIP`s (never fails) when the asset is absent.
macro_rules! bind_solloc21 {
    ($provider:ident, $nutation_at:ident, $serving:ident, $deltat_at:ident, $bound:ident,
         $lon:expr, $lat:expr, $alt:expr) => {
        let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
        let bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(_) => {
                eprintln!("evt14_slice21_solmagloc_tests: SKIP (missing {path})");
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
            $lon,
            $lat,
            $alt,
        )
        .expect("local-solar serving binds");
    };
}

#[test]
fn sol_magloc_non_finite_instants_refuse() {
    bind_solloc21!(
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
        assert!(
            matches!(
                sol_eclipse_magnitude_at_loc_sealed(&bound, bad),
                Err(Error::NonFiniteArgument { .. })
            ),
            "magnitude refuses {bad:?}"
        );
        assert!(
            matches!(
                sol_eclipse_obscuration_at_loc_sealed(&bound, bad),
                Err(Error::NonFiniteArgument { .. })
            ),
            "obscuration refuses {bad:?}"
        );
    }
}

#[test]
fn sol_magloc_totality_covers_disc() {
    bind_solloc21!(
        provider,
        nutation_at,
        serving,
        deltat_at,
        dallas,
        -96.7970,
        32.7767,
        131.0
    );
    // Banked 2024-04-08 total local maximum over Dallas.
    let jd_total = 2460409.279538682;
    let magnitude =
        sol_eclipse_magnitude_at_loc_sealed(&dallas, jd_total).expect("magnitude solves");
    let obscuration =
        sol_eclipse_obscuration_at_loc_sealed(&dallas, jd_total).expect("obscuration solves");
    assert!(
        magnitude >= 1.0,
        "totality covers the diameter, got {magnitude}"
    );
    assert_eq!(obscuration, 1.0, "totality covers the disc");
}

#[test]
fn sol_magloc_partial_ranks_magnitude_above_obscuration() {
    bind_solloc21!(
        provider,
        nutation_at,
        serving,
        deltat_at,
        denver,
        -104.9903,
        39.7392,
        1609.0
    );
    // Banked 2024-04-08 total local maximum seen partial from Denver.
    let jd_total = 2460409.279538682;
    let magnitude =
        sol_eclipse_magnitude_at_loc_sealed(&denver, jd_total).expect("magnitude solves");
    let obscuration =
        sol_eclipse_obscuration_at_loc_sealed(&denver, jd_total).expect("obscuration solves");
    assert!(
        0.0 < magnitude && magnitude < 1.0,
        "Denver sees a partial, got {magnitude}"
    );
    assert!(
        0.0 < obscuration && obscuration < magnitude,
        "area fraction trails the diameter fraction, got {obscuration} vs {magnitude}"
    );
}

#[test]
fn sol_magloc_quiet_and_night_instants_are_zero() {
    bind_solloc21!(
        provider,
        nutation_at,
        serving,
        deltat_at,
        dallas,
        -96.7970,
        32.7767,
        131.0
    );
    bind_solloc21!(
        provider2,
        nutation_at2,
        serving2,
        deltat_at2,
        tokyo,
        139.6917,
        35.6895,
        40.0
    );
    for (site, jd) in [(&dallas, 2460140.0), (&tokyo, 2460409.279538682)] {
        assert_eq!(
            sol_eclipse_magnitude_at_loc_sealed(site, jd).expect("magnitude degrades"),
            0.0
        );
        assert_eq!(
            sol_eclipse_obscuration_at_loc_sealed(site, jd).expect("obscuration degrades"),
            0.0
        );
    }
}
