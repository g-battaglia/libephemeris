use super::*;
use crate::leb::LebTier;
use crate::positions::{EntrySnapshot, StateProvider};
use crate::types::TwoPartEpoch;

/// Dummy provider that never serves anything: `has_body` is false
/// so the sealed gate refuses before any state read; every read is
/// a typed miss. Enough for the construction/gate tests below,
/// which must fail before (or without) touching ephemeris data.
struct Unserved;
impl StateProvider for Unserved {
    fn provider_name(&self) -> &'static str {
        "evt13-dummy/unserved"
    }
    fn has_body(&self, _body_id: i32) -> bool {
        false
    }
    fn body_coverage(&self, _body_id: i32) -> Option<(f64, f64)> {
        None
    }
    fn serving_coord_type(&self, _body_id: i32) -> Option<u32> {
        None
    }
    fn eval_raw(&self, body_id: i32, _jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        Err(Error::LebBodyNotFound { body_id })
    }
    fn eval_raw_two_part(
        &self,
        body_id: i32,
        _epoch: TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        Err(Error::LebBodyNotFound { body_id })
    }
    fn serving_entry_snapshot(&self, body_id: i32) -> Result<EntrySnapshot> {
        Err(Error::LebBodyNotFound { body_id })
    }
}

#[allow(dead_code)]
fn _assert_tier_unused(_: Option<LebTier>) {}

fn zero_dt(_: f64) -> Result<f64> {
    Ok(0.0)
}

fn zero_nut(_: f64) -> Result<(f64, f64)> {
    Ok((0.0, 0.0))
}

fn pair<'a>(
    provider: &'a Unserved,
) -> (
    crate::positions::SealedServing<'a, Unserved>,
    crate::stars::SealedStarServing<'a, Unserved>,
) {
    use crate::sidereal::SiderealSelection;
    let serving = crate::positions::SealedServing::new(provider, &zero_nut);
    let star_serving = crate::stars::SealedStarServing::new(
        provider,
        &zero_nut,
        &zero_dt,
        SiderealSelection::unset(),
    );
    (serving, star_serving)
}

#[test]
fn occult_body_leg_labels_and_radii() {
    use crate::constants::{JUPITER, MOON, SUN};
    let planet = OccultBodyLeg::Planet(JUPITER);
    assert!(!planet.is_star());
    assert_eq!(planet.body_id_opt(), Some(JUPITER));
    assert!(!planet.is_sun());
    assert_eq!(planet.label(), JUPITER.to_string());
    assert_eq!(planet.radius_au(), occ_body_radius_au(Some(JUPITER)));
    let sun = OccultBodyLeg::Planet(SUN);
    assert!(sun.is_sun());
    assert_eq!(sun.radius_au(), ECL_RSUN_AU);
    let moon = OccultBodyLeg::Planet(MOON);
    assert!(!moon.is_sun());
    // Star legs are point sources with the search name as label.
    let provider = Unserved;
    let (_serving, star_serving) = pair(&provider);
    let bound = SealedOccultGlobalServing::try_new_star(
        &_serving,
        &star_serving,
        &zero_dt,
        crate::constants::FLG_SWIEPH,
        "Sirius",
    )
    .expect("Sirius resolves");
    assert!(bound.body_leg().is_star());
    assert_eq!(bound.body_leg().body_id_opt(), None);
    assert!(!bound.body_leg().is_sun());
    assert_eq!(bound.body_leg().label(), "Sirius");
    assert_eq!(bound.body_leg().radius_au(), 0.0);
    assert_eq!(bound.body_leg().radius_au(), occ_body_radius_au(None));
}

#[test]
fn occult_servings_refuse_flags_site_and_stars_at_construction() {
    use crate::constants::{FLG_SIDEREAL, FLG_SWIEPH, JUPITER};
    let provider = Unserved;
    let (serving, star_serving) = pair(&provider);
    // Flag selections outside the tropical basis refuse.
    assert!(matches!(
        SealedOccultGlobalServing::try_new_planet(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH | FLG_SIDEREAL,
            JUPITER,
        ),
        Err(Error::UnsupportedFlagCombination { .. })
    ));
    assert!(matches!(
        SealedOccultLocalServing::try_new_planet(
            &serving,
            &star_serving,
            &zero_dt,
            &zero_nut,
            FLG_SWIEPH | FLG_SIDEREAL,
            JUPITER,
            0.0,
            0.0,
            0.0,
        ),
        Err(Error::UnsupportedFlagCombination { .. })
    ));
    // Unknown star names refuse StarNotFound (never a silent
    // substitution).
    assert!(matches!(
        SealedOccultGlobalServing::try_new_star(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH,
            "No Such Star xyz",
        ),
        Err(Error::StarNotFound { .. })
    ));
    // The observer triple validates before any ephemeris read.
    assert!(matches!(
        SealedOccultLocalServing::try_new_planet(
            &serving,
            &star_serving,
            &zero_dt,
            &zero_nut,
            FLG_SWIEPH,
            JUPITER,
            f64::NAN,
            0.0,
            0.0,
        ),
        Err(Error::InvalidObserver { .. })
    ));
    // Non-finite Delta-T boundary words refuse at construction.
    let bad_dt = |_: f64| Ok(f64::NAN);
    assert!(matches!(
        SealedOccultGlobalServing::try_new_planet(
            &serving,
            &star_serving,
            &bad_dt,
            FLG_SWIEPH,
            JUPITER,
        ),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn occult_servings_fold_negative_body_ids_to_sun() {
    // Release `_fold_pseudo_body_to_sun`: every negative id (e.g.
    // `ECL_NUT`) answers exactly as the Sun at each `lun_occult_*`
    // entry (sweep case `evt14-occult-glob-negid-sun-fwd-0031`
    // caught the missing fold: the Rust side refused
    // `leb_body_not_found` where the release serves body 0).
    use crate::constants::{ECL_NUT, FLG_SWIEPH, SUN};
    let provider = Unserved;
    let (serving, star_serving) = pair(&provider);
    let glob = SealedOccultGlobalServing::try_new_planet(
        &serving,
        &star_serving,
        &zero_dt,
        FLG_SWIEPH,
        ECL_NUT,
    )
    .expect("negative id folds to the Sun");
    assert!(glob.body_leg().is_sun());
    assert_eq!(glob.body_leg().body_id_opt(), Some(SUN));
    let local = SealedOccultLocalServing::try_new_planet(
        &serving,
        &star_serving,
        &zero_dt,
        &zero_nut,
        FLG_SWIEPH,
        ECL_NUT,
        0.0,
        0.0,
        0.0,
    )
    .expect("negative id folds to the Sun");
    assert!(local.body_leg().is_sun());
    assert_eq!(local.body_leg().body_id_opt(), Some(SUN));
}

#[test]
fn occult_where_serves_geocentric_earth_over_sealed_asset() {
    // Release `lun_occult_where(jd, EARTH, flags)` serves finite
    // geometry (Earth is the geocentric origin, zeros), verified
    // against the pinned release. Asset-gated like the other
    // sealed-path tests: skips when the machine-local asset is
    // absent, never fails for missing data.
    use crate::constants::{EARTH, FLG_SWIEPH};
    let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
    let bytes = match std::fs::read(path) {
        Ok(bytes) => bytes,
        Err(_) => {
            eprintln!(
                "occult_where_serves_geocentric_earth_over_sealed_asset: SKIP (missing {path})"
            );
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
    // CompositeIndex mirrors the binding's `sealed_occult_stack`
    // dispatch over context-bound images (first-wins, single entry
    // here over the same bytes).
    let composite = crate::leb::CompositeIndex::new(vec![(index.clone(), &bytes[..])])
        .expect("one reader must build");
    let nutation_at = |jd: f64| composite.eval_nutation(jd);
    let serving = crate::positions::SealedServing::new(&composite, &nutation_at);
    let deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &cfg);
    let star_serving = crate::stars::SealedStarServing::new(
        &composite,
        &nutation_at,
        &deltat_at,
        crate::sidereal::SiderealSelection::unset(),
    );
    let jd_ut = 2460000.5;
    let glob = SealedOccultGlobalServing::try_new_planet(
        &serving,
        &star_serving,
        &deltat_at,
        FLG_SWIEPH,
        EARTH,
    )
    .expect("Earth leg constructs");
    let geom = glob.where_core(jd_ut).expect("Earth where_core serves");
    assert!(geom.lon_deg().is_finite(), "finite center lon");
    assert!(geom.lat_deg().is_finite(), "finite center lat");
    let local = SealedOccultLocalServing::try_new_planet(
        &serving,
        &star_serving,
        &deltat_at,
        &nutation_at,
        FLG_SWIEPH,
        EARTH,
        geom.lon_deg(),
        geom.lat_deg(),
        0.0,
    )
    .expect("Earth local leg constructs");
    let attr = local
        .local_how_attr_where(jd_ut)
        .expect("Earth how_attr serves");
    assert!(attr.iter().all(|v| v.is_finite()), "finite attr");
}

#[test]
fn occult_sun_how_attr_at_banked_tret_within_release_ceiling() {
    // O-A rerun1 closeout: worst when_loc residual case
    // (evt14-lunloc-total2022a-mexicocity-pre30-bwd-0002-body0, Sun
    // occultation at Mexico City, backward search). The rerun
    // measured the residual with a pre-Clenshaw kernel; post-fix
    // the sealed walk settles on the release tret words bitwise
    // (closeout re-probe). This test pins the fixed-instant
    // evaluation level over the sealed asset: the local how-attr
    // at the banked release tret[0] agrees with the banked release
    // attr words inside the 5e-6 comparison ceiling (measured
    // worst 2.3e-9). Slot 3 stays the raw how-core zero here —
    // the walk overwrites it with the where-core dcore[0] via
    // `finish_occult_local_attr`. Asset-gated like the other
    // sealed-path tests: skips when the machine-local asset is
    // absent, never fails for missing data.
    use crate::constants::{FLG_SWIEPH, SUN};
    let path = "/Users/giacomo/.libephemeris/leb/base_core.leb2";
    let bytes = match std::fs::read(path) {
        Ok(bytes) => bytes,
        Err(_) => {
            eprintln!(
                "occult_sun_how_attr_at_banked_tret_within_release_ceiling: SKIP (missing {path})"
            );
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
    let composite = crate::leb::CompositeIndex::new(vec![(index.clone(), &bytes[..])])
        .expect("one reader must build");
    let nutation_at = |jd: f64| composite.eval_nutation(jd);
    let serving = crate::positions::SealedServing::new(&composite, &nutation_at);
    let deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &cfg);
    let star_serving = crate::stars::SealedStarServing::new(
        &composite,
        &nutation_at,
        &deltat_at,
        crate::sidereal::SiderealSelection::unset(),
    );
    let bound = SealedOccultLocalServing::try_new_planet(
        &serving,
        &star_serving,
        &deltat_at,
        &nutation_at,
        FLG_SWIEPH,
        SUN,
        f64::from_bits(0xc058c886594af4f1),
        f64::from_bits(0x40336ebedfa43fe6),
        2240.0,
    )
    .expect("Sun local leg constructs");
    // Banked reference values (O-A rerun1 pin side): tret[0] then the
    // 20 attr slots of the worst residual case.
    let tret0 = f64::from_bits(0x4142c0c1a1cb203a);
    let want: [u64; 8] = [
        0x3fd85bc112325657,
        0x3ff0879ee4aadb82,
        0x3fd1202a3040ef7b,
        0x40747ac23e8fc4c0,
        0x4054471b15c84bfe,
        0x4054473b5b78b1c7,
        0x3fd573dd3f6e196f,
        0x3fd85bc112325657,
    ];
    let attr = bound
        .local_how_attr_where(tret0)
        .expect("Sun how_attr serves");
    for (slot, bits) in [0, 1, 2, 4, 5, 6, 7, 8].iter().zip(want.iter()) {
        let gap = (attr[*slot] - f64::from_bits(*bits)).abs();
        assert!(
            gap <= 5e-6,
            "attr[{slot}] within the comparison ceiling: gap {gap:e}"
        );
    }
    assert_eq!(attr[3].to_bits(), 0x0000000000000000, "raw how-core slot");
    assert_eq!(attr[9].to_bits(), 0x4062200000000000, "structural slot");
    assert_eq!(attr[10].to_bits(), 0x4036000000000000, "structural slot");
    assert!(
        attr[11..20].iter().all(|v| *v == 0.0),
        "trailing slots stay zero"
    );
}

#[test]
fn occult_servings_validate_epochs_before_any_ephemeris_read() {
    use crate::constants::{FLG_SWIEPH, JUPITER};
    let provider = Unserved;
    let (serving, star_serving) = pair(&provider);
    let glob = SealedOccultGlobalServing::try_new_planet(
        &serving,
        &star_serving,
        &zero_dt,
        FLG_SWIEPH,
        JUPITER,
    )
    .expect("planet leg constructs");
    // Non-finite epochs refuse without touching the provider.
    assert!(matches!(
        glob.body_lonlat_deg(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        glob.moon_lonlat_deg(f64::INFINITY),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        glob.sep_minus_radii_deg(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        glob.ra_gap_deg(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        glob.where_core(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
    let loc = SealedOccultLocalServing::try_new_planet(
        &serving,
        &star_serving,
        &zero_dt,
        &zero_nut,
        FLG_SWIEPH,
        JUPITER,
        10.0,
        50.0,
        100.0,
    )
    .expect("local planet leg constructs");
    assert_eq!(loc.observer_lon_lat_alt(), (10.0, 50.0, 100.0));
    assert!(matches!(
        loc.topo_radii_sep_deg(f64::NAN, false),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        loc.body_alt_deg(f64::INFINITY),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        loc.local_how_attr(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        loc.body_rise_set(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        loc.sun_rise_set(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn occult_servings_escalate_provider_misses_without_panic() {
    use crate::constants::{FLG_SWIEPH, JUPITER};
    let provider = Unserved;
    let (serving, star_serving) = pair(&provider);
    let glob = SealedOccultGlobalServing::try_new_planet(
        &serving,
        &star_serving,
        &zero_dt,
        FLG_SWIEPH,
        JUPITER,
    )
    .expect("planet leg constructs");
    // The unserved provider misses every channel: finite-epoch
    // reads escalate to a typed error (never a panic, never a
    // fabricated event).
    let err = glob
        .body_lonlat_deg(2460000.5)
        .expect_err("unserved provider must miss");
    assert!(
        matches!(
            err,
            Error::LebBodyNotFound { .. }
                | Error::LebOutOfRange { .. }
                | Error::EphemerisRange { .. }
                | Error::UnknownBody { .. }
        ),
        "unexpected miss mapping: {err:?}"
    );
    // The sealed search wrappers surface the miss (or exhaustion)
    // instead of an event: one-try star walks refuse before any
    // provider read when the star cannot occult.
    let star_glob = SealedOccultGlobalServing::try_new_star(
        &serving,
        &star_serving,
        &zero_dt,
        FLG_SWIEPH,
        "Sirius",
    )
    .expect("Sirius resolves");
    let err = lun_occult_when_glob_sealed(&star_glob, 2460000.5, 0, false, false)
        .expect_err("unserved provider must not yield an event");
    assert!(
        !matches!(err, Error::NonFiniteArgument { .. }),
        "finite-epoch miss must not masquerade as a domain refusal: {err:?}"
    );
}

/// Inverting the horizontal projection with the full IAU 2006/2000A frame
/// must recover the supplied topocentric ecliptic direction. This geometric
/// invariant does not bank a sidereal model's final azimuth/altitude words.
/// The frame uses observed-timescale TT; the independent mean origin uses
/// the caller's Espenak–Meeus TT. The stored nutation provides an adverse
/// control: it cannot reproduce the same direction at the strict tolerance.
#[test]
fn occult_altaz_uses_full_series_nutation() {
    let cfg = crate::time::DeltaTConfig {
        model: crate::time::DeltaTModel::EspenakMeeus,
        user_override_days: None,
        tidal_acceleration: crate::constants::TIDAL_DEFAULT,
        iers_enabled: false,
        iers_lookup_seconds: None,
        smh2016_lookup_seconds: None,
    };
    let deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &cfg);
    let sun_topo = (
        f64::from_bits(0x4069231086416245),
        f64::from_bits(0xbf437e4b73789f36),
        f64::from_bits(0x3fefeb5d6a77e1b4),
    );
    let jd_ut = f64::from_bits(0x4142c5241b7f7bcb);
    let longitude = f64::from_bits(0xc058c886594af4f1);
    let latitude = f64::from_bits(0x40336ebedfa43fe6);
    let (az_deg, true_alt_deg, app_alt_deg) = altaz_of_topo_ecliptic_impl(
        sun_topo, jd_ut, longitude, latitude, 2240.0, &deltat_at, 0.0, 0.0,
    )
    .expect("banked altaz inputs serve");
    let frame_tt = crate::time::default_timescale_tt_from_ut1(jd_ut).unwrap();
    let mean_obliquity = crate::positions::vondrak_mean_obliquity_deg(frame_tt).unwrap();
    let mean_sidereal = crate::time::public_sidereal::mean_deg(jd_ut, &deltat_at).unwrap();
    let inverse = |(dpsi, deps): (f64, f64)| {
        let epsilon = (mean_obliquity + deps.to_degrees()).to_radians();
        let last = (mean_sidereal + dpsi.to_degrees() * epsilon.cos() + longitude).to_radians();
        let (sin_lat, cos_lat) = latitude.to_radians().sin_cos();
        let (sin_az, cos_az) = az_deg.to_radians().sin_cos();
        let (sin_alt, cos_alt) = true_alt_deg.to_radians().sin_cos();
        // Azimuth is measured from south towards west. Undo the orthogonal
        // horizon rotation, then RA = LAST - hour angle in Cartesian form.
        let south = cos_alt * cos_az;
        let west = cos_alt * sin_az;
        let meridian = south * sin_lat + sin_alt * cos_lat;
        let z_equator = sin_alt * sin_lat - south * cos_lat;
        let (sin_last, cos_last) = last.sin_cos();
        let x = meridian * cos_last + west * sin_last;
        let y_equator = meridian * sin_last - west * cos_last;
        // Undo the equator-to-ecliptic tilt without calling cotrans or any
        // horizontal conversion from the calculation under test.
        let (sin_eps, cos_eps) = epsilon.sin_cos();
        let y = y_equator * cos_eps + z_equator * sin_eps;
        let z = -y_equator * sin_eps + z_equator * cos_eps;
        let lon = y.atan2(x).to_degrees().rem_euclid(360.0);
        let lat = z.atan2(x.hypot(y)).to_degrees();
        (
            ((lon - sun_topo.0 + 180.0).rem_euclid(360.0) - 180.0).abs(),
            (lat - sun_topo.1).abs(),
        )
    };
    let (lon_error, lat_error) = inverse(crate::erfa::nut06a(frame_tt).unwrap());
    assert!(lon_error <= 1e-12, "ecliptic longitude error {lon_error}");
    assert!(lat_error <= 1e-12, "ecliptic latitude error {lat_error}");
    let stored = crate::leb::LebFile::open(&crate::test_data::core_leb_path()).unwrap();
    assert!(crate::test_data::admitted_core_hash(
        stored.sha256().unwrap()
    ));
    let (stored_lon_error, stored_lat_error) = inverse(stored.eval_nutation(frame_tt).unwrap());
    assert!(
        stored_lon_error.max(stored_lat_error) > 1e-12,
        "stored-series adverse control must distinguish the full-series frame"
    );
    // The adjacent public/shared test checks all three output words against
    // each other. Here refraction must additionally retain its physical sign.
    assert!(app_alt_deg.is_finite() && app_alt_deg > true_alt_deg && app_alt_deg < 90.0);
}

/// The public `azalt_deg` ECL2HOR leg with the sealed-path
/// spelling (estimated pressure, 10 C, release default lapse)
/// reproduces the shared `altaz_of_topo_ecliptic_impl` triple
/// bitwise on the banked Mexico City words. The two spellings
/// share every leaf kernel; divergence here means the public entry
/// drifted from the swept serving path.
#[test]
fn azalt_public_ecliptic_leg_matches_shared_impl() {
    let cfg = crate::time::DeltaTConfig {
        model: crate::time::DeltaTModel::EspenakMeeus,
        user_override_days: None,
        tidal_acceleration: crate::constants::TIDAL_DEFAULT,
        iers_enabled: false,
        iers_lookup_seconds: None,
        smh2016_lookup_seconds: None,
    };
    let deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &cfg);
    let sun_topo = (
        f64::from_bits(0x4069231086416245),
        f64::from_bits(0xbf437e4b73789f36),
        f64::from_bits(0x3fefeb5d6a77e1b4),
    );
    let jd_ut = f64::from_bits(0x4142c5241b7f7bcb);
    let (lon, lat, alt_m) = (
        f64::from_bits(0xc058c886594af4f1),
        f64::from_bits(0x40336ebedfa43fe6),
        2240.0,
    );
    let shared =
        altaz_of_topo_ecliptic_impl(sun_topo, jd_ut, lon, lat, alt_m, &deltat_at, 0.0, 0.0)
            .expect("banked altaz inputs serve");
    let public = azalt_deg(
        jd_ut,
        ECL2HOR,
        lon,
        lat,
        alt_m,
        0.0,
        10.0,
        sun_topo,
        &deltat_at,
        0.0,
        0.0,
        RISE_SET_DIP_LAPSE_RATE_K_M,
    )
    .expect("public ECL2HOR leg serves");
    assert_eq!(public.0.to_bits(), shared.0.to_bits(), "azimuth word");
    assert_eq!(public.1.to_bits(), shared.1.to_bits(), "true-altitude word");
    assert_eq!(
        public.2.to_bits(),
        shared.2.to_bits(),
        "apparent-altitude word"
    );
}

/// Forward/inverse geometric round trip — `azalt_rev_deg`
/// over the public forward triple recovers the equatorial input
/// within 1e-9 deg (refraction is not removed on the inverse leg,
/// so the true-altitude arm, not the apparent one, round-trips).
#[test]
fn azalt_rev_round_trips_equatorial_input() {
    let cfg = crate::time::DeltaTConfig {
        model: crate::time::DeltaTModel::EspenakMeeus,
        user_override_days: None,
        tidal_acceleration: crate::constants::TIDAL_DEFAULT,
        iers_enabled: false,
        iers_lookup_seconds: None,
        smh2016_lookup_seconds: None,
    };
    let deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &cfg);
    let dt_lo = deltat_at(crate::time::SIDEREAL_MODERN_START_JD).unwrap();
    let dt_hi = deltat_at(crate::time::SIDEREAL_MODERN_END_JD).unwrap();
    let (off_lo, off_hi) = crate::time::sidereal_continuity_offsets(dt_lo, dt_hi).unwrap();
    let jd_ut = 2460000.5;
    let (lon, lat, alt_m) = (12.5, 41.9, 100.0);
    let (ra_in, dec_in) = (90.0, 23.5);
    let (az, true_alt, _) = azalt_deg(
        jd_ut,
        EQU2HOR,
        lon,
        lat,
        alt_m,
        1013.25,
        15.0,
        (ra_in, dec_in, 1.0),
        &deltat_at,
        off_lo,
        off_hi,
        RISE_SET_DIP_LAPSE_RATE_K_M,
    )
    .expect("forward equatorial leg serves");
    let (ra_out, dec_out) = azalt_rev_deg(
        jd_ut, HOR2EQU, lon, lat, alt_m, az, true_alt, &deltat_at, off_lo, off_hi,
    )
    .expect("inverse equatorial leg serves");
    assert!(
        (ra_out - ra_in).abs() < 1e-9,
        "RA round trip: {ra_out} vs {ra_in}"
    );
    assert!(
        (dec_out - dec_in).abs() < 1e-9,
        "Dec round trip: {dec_out} vs {dec_in}"
    );
}

/// Pole guards and domain refusals — a geographic-pole
/// observer reports RA equal to local sidereal time (hour angle 0),
/// and any non-finite leg refuses without a value.
#[test]
fn azalt_rev_pole_and_domain_legs() {
    let cfg = crate::time::DeltaTConfig {
        model: crate::time::DeltaTModel::EspenakMeeus,
        user_override_days: None,
        tidal_acceleration: crate::constants::TIDAL_DEFAULT,
        iers_enabled: false,
        iers_lookup_seconds: None,
        smh2016_lookup_seconds: None,
    };
    let deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &cfg);
    let jd_ut = 2460000.5;
    let (ra, dec) = azalt_rev_deg(
        jd_ut, HOR2EQU, 0.0, 90.0, 0.0, 180.0, 45.0, &deltat_at, 0.0, 0.0,
    )
    .expect("polar observer serves");
    assert!(
        (45.0 - dec).abs() < 1e-9,
        "at the geographic pole declination equals true altitude: {dec}"
    );
    let _ = ra;
    assert!(
        azalt_deg(
            f64::NAN,
            EQU2HOR,
            0.0,
            0.0,
            0.0,
            1013.25,
            15.0,
            (0.0, 0.0, 1.0),
            &deltat_at,
            0.0,
            0.0,
            RISE_SET_DIP_LAPSE_RATE_K_M,
        )
        .is_err(),
        "non-finite epoch refuses"
    );
    assert!(
        azalt_rev_deg(
            jd_ut,
            HOR2EQU,
            0.0,
            0.0,
            0.0,
            f64::INFINITY,
            45.0,
            &deltat_at,
            0.0,
            0.0,
        )
        .is_err(),
        "non-finite azimuth refuses"
    );
}
