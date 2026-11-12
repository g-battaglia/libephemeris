use super::*;

/// Occult-global fixture over the sealed asset for one numeric
/// body id (Espenak-Meeus Delta-T). `SKIP`s (never fails) when
/// the asset is absent.
macro_rules! bind_occglob24 {
    ($provider:ident, $nutation_at:ident, $serving:ident, $deltat_at:ident, $star:ident, $bound:ident, $body:expr) => {
        let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
        let bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(_) => {
                eprintln!("evt14_slice24_occwhere_tests: SKIP (missing {path})");
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
        let $star = crate::stars::SealedStarServing::new(
            &$provider,
            &$nutation_at,
            &$deltat_at,
            crate::sidereal::SiderealSelection::unset(),
        );
        let $bound = SealedOccultGlobalServing::try_new_planet(
            &$serving,
            &$star,
            &$deltat_at,
            crate::constants::FLG_SWIEPH,
            $body,
        )
        .expect("occult-global serving binds");
    };
}

#[test]
fn occult_where_wiring_matches_where_core_words() {
    bind_occglob24!(
        provider,
        nutation_at,
        serving,
        deltat_at,
        star,
        bound,
        crate::constants::VENUS
    );
    // Banked Venus occultation maximum (slice-14 release tret[0]
    // word): the fixed-instant composition must report an event
    // in progress with the where-core words wired through
    //  — character bits, center lon/lat and dcore.
    let jd_max = 2460027.938786889;
    let outcome = bound
        .lun_occult_where_at_ut(&nutation_at, jd_max)
        .expect("where resolves at the banked maximum");
    assert_ne!(outcome.retflag(), 0, "occultation in progress at maximum");
    let geom = bound
        .where_geometry_at_ut(jd_max)
        .expect("where core reads");
    assert_eq!(outcome.retflag(), geom.retc());
    assert_eq!(outcome.geopos()[0], geom.lon_deg());
    assert_eq!(outcome.geopos()[1], geom.lat_deg());
    assert_eq!(outcome.attr()[3], geom.dcore()[0]);
    assert!(
        outcome.geopos()[2..].iter().all(|w| *w == 0.0),
        "reserved geopos slots stay zero"
    );
}

#[test]
fn occult_where_reports_quiet_instant_as_no_event() {
    bind_occglob24!(
        provider,
        nutation_at,
        serving,
        deltat_at,
        star,
        bound,
        crate::constants::VENUS
    );
    // Quiet JD with no occultation anywhere (slice-23 quiet arm):
    // retflag 0 with the closest-approach center point carried.
    let outcome = bound
        .lun_occult_where_at_ut(&nutation_at, 2460140.0)
        .expect("where resolves off-event");
    assert_eq!(outcome.retflag(), 0);
    assert!(
        outcome.geopos()[0].is_finite() && outcome.geopos()[1].is_finite(),
        "closest-approach point carried, got {:?}",
        outcome.geopos()
    );
}

#[test]
fn occult_where_refuses_moon_and_nonfinite() {
    bind_occglob24!(
        provider,
        nutation_at,
        serving,
        deltat_at,
        star,
        bound,
        crate::constants::VENUS
    );
    assert!(
        bound
            .lun_occult_where_at_ut(&nutation_at, f64::NAN)
            .is_err()
    );
    bind_occglob24!(
        provider2,
        nutation_at2,
        serving2,
        deltat_at2,
        star2,
        moon,
        crate::constants::MOON
    );
    assert!(matches!(
        moon.lun_occult_where_at_ut(&nutation_at2, 2460027.938786889),
        Err(Error::InvalidEclipseSearch { .. })
    ));
}
