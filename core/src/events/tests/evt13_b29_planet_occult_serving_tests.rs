use super::*;
use crate::leb::{
    COORD_ICRS_BARY, CompositeIndex, FileIndex, HEADER_SIZE, LEB1_BODY_ENTRY_SIZE, MAGIC_LEB1,
    SECTION_DIR_SIZE,
};

fn push_u32(v: &mut Vec<u8>, x: u32) {
    v.extend_from_slice(&x.to_le_bytes());
}

fn push_i32(v: &mut Vec<u8>, x: i32) {
    v.extend_from_slice(&x.to_le_bytes());
}

fn push_u64(v: &mut Vec<u8>, x: u64) {
    v.extend_from_slice(&x.to_le_bytes());
}

fn push_f64(v: &mut Vec<u8>, x: f64) {
    v.extend_from_slice(&x.to_le_bytes());
}

/// Hand-packed LEB1 image over JD `[0, 2]`: one section, degree-1
/// segments with constant positions (zero velocity), ICRS
/// barycentric AU. Same layout as the positions-fixture precedent;
/// no release coefficient is stored anywhere here.
fn fixture_icrs_const(states: &[(i32, [f64; 3])]) -> Vec<u8> {
    let n = states.len();
    let mut file = Vec::with_capacity(HEADER_SIZE);
    file.extend_from_slice(&MAGIC_LEB1);
    push_u32(&mut file, 1);
    push_u32(&mut file, 1);
    push_u32(&mut file, n as u32);
    push_f64(&mut file, 0.0);
    push_f64(&mut file, 2.0);
    push_f64(&mut file, 0.0);
    push_u32(&mut file, 0);
    file.extend_from_slice(&[0u8; 20]);
    debug_assert_eq!(file.len(), HEADER_SIZE);
    let body_at = HEADER_SIZE + SECTION_DIR_SIZE;
    let blob_at = body_at + n * LEB1_BODY_ENTRY_SIZE;
    let mut section = Vec::with_capacity(SECTION_DIR_SIZE);
    push_u32(&mut section, 0);
    push_u32(&mut section, 0);
    push_u64(&mut section, body_at as u64);
    push_u64(&mut section, (n * LEB1_BODY_ENTRY_SIZE) as u64);
    file.extend_from_slice(&section);
    for (i, (body_id, _)) in states.iter().enumerate() {
        // 2 segments over [0, 2], 1-day interval, degree 1, 3
        // components.
        let mut entry = Vec::with_capacity(LEB1_BODY_ENTRY_SIZE);
        push_i32(&mut entry, *body_id);
        push_u32(&mut entry, COORD_ICRS_BARY);
        push_u32(&mut entry, 2);
        push_f64(&mut entry, 0.0);
        push_f64(&mut entry, 2.0);
        push_f64(&mut entry, 1.0);
        push_u32(&mut entry, 1);
        push_u32(&mut entry, 3);
        push_u64(&mut entry, (blob_at + i * 96) as u64);
        file.extend_from_slice(&entry);
    }
    for (_, pos) in states {
        for _ in 0..2 {
            for component in pos {
                push_f64(&mut file, *component);
                push_f64(&mut file, 0.0);
            }
        }
    }
    debug_assert_eq!(file.len(), blob_at + n * 96);
    file
}

/// Venus/Jupiter conjunction fixture: Venus at 0.2 AU and Jupiter
/// at 1.0 AU from Earth along nearly the same geocentric direction
/// (discs overlap by a factor of three), so the sealed global walk
/// opens on a hit at the first sample.
fn occult_fixture() -> Vec<u8> {
    use crate::constants::{EARTH, JUPITER, SUN, VENUS};
    fixture_icrs_const(&[
        (SUN, [0.0, 0.0, 0.0]),
        (EARTH, [1.0, 0.0, 0.0]),
        (VENUS, [1.2, 0.0001, 0.0]),
        (JUPITER, [2.0, 0.0003, 0.0]),
    ])
}

fn zero_dt(_: f64) -> Result<f64> {
    Ok(0.0)
}

fn zero_nut(_: f64) -> Result<(f64, f64)> {
    Ok((0.0, 0.0))
}

/// Run `f` over the bound sealed serving pair on the conjunction
/// fixture (borrowing stays inside the closure).
fn with_servings<R>(
    f: impl FnOnce(
        crate::positions::SealedServing<'_, CompositeIndex<'_>>,
        crate::stars::SealedStarServing<'_, CompositeIndex<'_>>,
    ) -> R,
) -> R {
    use crate::sidereal::SiderealSelection;
    let file = occult_fixture();
    let owned = FileIndex::decode(&file).expect("fixture must decode");
    let composite = CompositeIndex::new(vec![(owned, &file[..])]).expect("one reader must build");
    let serving = crate::positions::SealedServing::new(&composite, &zero_nut);
    let star_serving = crate::stars::SealedStarServing::new(
        &composite,
        &zero_nut,
        &zero_dt,
        SiderealSelection::unset(),
    );
    f(serving, star_serving)
}

#[test]
fn planet_occult_servings_refuse_requests_flags_site_and_stars_at_construction() {
    use crate::constants::{FLG_SIDEREAL, FLG_SWIEPH, JUPITER, MOON, SUN, VENUS};
    with_servings(|serving, star_serving| {
        // Sun/Moon can never occult: the release constructors refuse.
        // (`Result::err` needs no `Debug` on the serving, unlike
        // `expect_err`.)
        let err = SealedPlanetOccultGlobalServing::try_new_planet(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH,
            SUN,
            JUPITER,
        )
        .err()
        .expect("construction must refuse");
        assert!(format!("{err:?}").contains("Sun cannot be an occulting body"));
        let err = SealedPlanetOccultGlobalServing::try_new_planet(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH,
            MOON,
            JUPITER,
        )
        .err()
        .expect("construction must refuse");
        assert!(format!("{err:?}").contains("Moon cannot be an occulting body"));
        // Classical-set membership gates both legs.
        let err = SealedPlanetOccultGlobalServing::try_new_planet(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH,
            12345,
            JUPITER,
        )
        .err()
        .expect("construction must refuse");
        assert!(format!("{err:?}").contains("Invalid occulting planet ID: 12345"));
        let err = SealedPlanetOccultGlobalServing::try_new_planet(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH,
            VENUS,
            12345,
        )
        .err()
        .expect("construction must refuse");
        assert!(format!("{err:?}").contains("Invalid occulted planet ID: 12345"));
        // Same-body and missing-target requests refuse.
        let err = SealedPlanetOccultGlobalServing::try_new_planet(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH,
            VENUS,
            VENUS,
        )
        .err()
        .expect("construction must refuse");
        assert!(format!("{err:?}").contains("cannot be the same"));
        let err = SealedPlanetOccultGlobalServing::try_new_star(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH,
            VENUS,
            "",
        )
        .err()
        .expect("construction must refuse");
        assert!(format!("{err:?}").contains("Either occulted_planet ID or starname"));
        // Unknown star names refuse StarNotFound (never a silent
        // substitution).
        assert!(matches!(
            SealedPlanetOccultGlobalServing::try_new_star(
                &serving,
                &star_serving,
                &zero_dt,
                FLG_SWIEPH,
                VENUS,
                "No Such Star xyz",
            ),
            Err(Error::StarNotFound { .. })
        ));
        // Flag selections outside the tropical basis refuse on both
        // servings.
        assert!(matches!(
            SealedPlanetOccultGlobalServing::try_new_planet(
                &serving,
                &star_serving,
                &zero_dt,
                FLG_SWIEPH | FLG_SIDEREAL,
                VENUS,
                JUPITER,
            ),
            Err(Error::UnsupportedFlagCombination { .. })
        ));
        assert!(matches!(
            SealedPlanetOccultLocalServing::try_new_planet(
                &serving,
                &star_serving,
                &zero_dt,
                FLG_SWIEPH | FLG_SIDEREAL,
                VENUS,
                JUPITER,
                10.0,
                50.0,
                100.0,
            ),
            Err(Error::UnsupportedFlagCombination { .. })
        ));
        // The observer triple validates before any ephemeris read.
        assert!(matches!(
            SealedPlanetOccultLocalServing::try_new_planet(
                &serving,
                &star_serving,
                &zero_dt,
                FLG_SWIEPH,
                VENUS,
                JUPITER,
                f64::NAN,
                50.0,
                100.0,
            ),
            Err(Error::InvalidObserver { .. })
        ));
        // Non-finite Delta-T boundary words refuse at construction.
        let bad_dt = |_: f64| Ok(f64::NAN);
        assert!(matches!(
            SealedPlanetOccultGlobalServing::try_new_planet(
                &serving,
                &star_serving,
                &bad_dt,
                FLG_SWIEPH,
                VENUS,
                JUPITER,
            ),
            Err(Error::NonFiniteArgument { .. })
        ));
    });
}

#[test]
fn planet_occult_acquisition_serves_finite_words_over_fixture() {
    use crate::constants::{FLG_SWIEPH, JUPITER, VENUS};
    with_servings(|serving, star_serving| {
        let glob = SealedPlanetOccultGlobalServing::try_new_planet(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH,
            VENUS,
            JUPITER,
        )
        .expect("planet pair constructs");
        assert_eq!(glob.occulting_id(), VENUS);
        assert_eq!(glob.occulted_id_opt(), Some(JUPITER));
        assert_eq!(glob.target_starname(), "");
        // Geocentric equatorial words are finite with the fixture
        // distances (Venus 0.2 AU, Jupiter 1.0 AU from Earth).
        let (occ_ra, occ_dec) = glob.occulting_ra_dec_deg(1.0).unwrap();
        assert!(occ_ra.is_finite() && occ_dec.is_finite());
        let (tgt_ra, tgt_dec) = glob.target_ra_dec_deg(1.0).unwrap();
        assert!(tgt_ra.is_finite() && tgt_dec.is_finite());
        let occ_r = glob.occulting_radius_deg(1.0).unwrap();
        let tgt_r = glob.target_radius_deg(1.0).unwrap();
        assert!(occ_r > 0.0 && tgt_r > 0.0);
        // Radii stay near the naive geometric values (Venus 0.2 AU
        // distant, Jupiter 1.0 AU): the sealed pipeline only shapes
        // the words, it does not rescale them.
        assert!((occ_r - 0.0116).abs() < 0.01);
        assert!((tgt_r - 0.0268).abs() < 0.01);
        let sep = planet_occult_separation_via(&glob, 1.0).expect("separation serves");
        assert!(sep.is_finite());
        assert!(
            planet_occult_threshold_met(sep, occ_r, tgt_r),
            "fixture conjunction must overlap: sep={sep}"
        );
        // Star targets resolve to finite equatorial words with the
        // point-source radius.
        let star_glob = SealedPlanetOccultGlobalServing::try_new_star(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH,
            VENUS,
            "Sirius",
        )
        .expect("Sirius resolves");
        assert_eq!(star_glob.occulted_id_opt(), None);
        assert_eq!(star_glob.target_starname(), "Sirius");
        let (sra, sdec) = star_glob.target_ra_dec_deg(1.0).unwrap();
        assert!(sra.is_finite() && sdec.is_finite());
        assert_eq!(
            star_glob.target_radius_deg(1.0).unwrap(),
            OCC_POINT_SOURCE_RADIUS_DEG
        );
    });
}

#[test]
fn planet_occult_global_walk_finds_fixture_conjunction() {
    use crate::constants::{ECL_PARTIAL, ECL_TOTAL, FLG_SWIEPH, JUPITER, VENUS};
    with_servings(|serving, star_serving| {
        let glob = SealedPlanetOccultGlobalServing::try_new_planet(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH,
            VENUS,
            JUPITER,
        )
        .expect("planet pair constructs");
        let (retflag, tret) = planet_occult_when_glob_sealed(&glob, 1.0, 0)
            .expect("fixture conjunction must be found");
        assert!(
            retflag == ECL_TOTAL || retflag == ECL_PARTIAL,
            "unexpected retflag {retflag}"
        );
        assert!(tret[0] > 0.0 && tret[0] <= 2.0);
        assert!(tret[2] > 0.0 && tret[3] > tret[2]);
        assert!((tret[2] <= tret[0]) && (tret[0] <= tret[3]));
    });
}

#[test]
fn planet_occult_servings_validate_epochs_before_any_ephemeris_read() {
    use crate::constants::{FLG_SWIEPH, JUPITER, VENUS};
    with_servings(|serving, star_serving| {
        let glob = SealedPlanetOccultGlobalServing::try_new_planet(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH,
            VENUS,
            JUPITER,
        )
        .expect("planet pair constructs");
        // Non-finite epochs refuse without touching the provider.
        assert!(matches!(
            glob.occulting_ra_dec_deg(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            glob.target_ra_dec_deg(f64::INFINITY),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            glob.occulting_radius_deg(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            glob.target_radius_deg(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        let loc = SealedPlanetOccultLocalServing::try_new_planet(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH,
            VENUS,
            JUPITER,
            10.0,
            50.0,
            100.0,
        )
        .expect("local planet pair constructs");
        assert_eq!(loc.observer_lon_lat_alt(), (10.0, 50.0, 100.0));
        assert!(matches!(
            loc.occulting_alt_az_deg(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            loc.target_alt_az_deg(f64::INFINITY),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            loc.local_radii_sep_deg(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
    });
}

#[test]
fn planet_occult_local_acquisition_serves_finite_words_over_fixture() {
    use crate::constants::{FLG_SWIEPH, JUPITER, VENUS};
    with_servings(|serving, star_serving| {
        let loc = SealedPlanetOccultLocalServing::try_new_planet(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH,
            VENUS,
            JUPITER,
            10.0,
            50.0,
            100.0,
        )
        .expect("local planet pair constructs");
        // Topocentric altitudes/azimuths are finite south-based
        // `azalt`-convention words.
        let (occ_alt, occ_az) = loc.occulting_alt_az_deg(1.0).unwrap();
        assert!(occ_alt.is_finite() && occ_az.is_finite());
        assert!((-90.0..=90.0).contains(&occ_alt));
        assert!((0.0..360.0).contains(&occ_az));
        let (tgt_alt, tgt_az) = loc.target_alt_az_deg(1.0).unwrap();
        assert!(tgt_alt.is_finite() && tgt_az.is_finite());
        // The topocentric radii/separation composition mirrors the
        // release maximum branch: both discs finite with a small
        // separation on the conjunction fixture.
        let (occ_r, tgt_r, sep) = loc.local_radii_sep_deg(1.0).unwrap();
        assert!(occ_r > 0.0 && tgt_r > 0.0);
        assert!((occ_r - 0.0116).abs() < 0.01);
        assert!((tgt_r - 0.0268).abs() < 0.01);
        assert!(sep < occ_r + tgt_r);
        // The sealed local wrapper runs the full walk: either the
        // fixture geometry is visible from the site (composed
        // event) or it is not (typed exhaustion) — never a panic
        // or a fabricated event.
        match planet_occult_when_loc_sealed(&loc, 1.0) {
            Ok((retflag, times, attr)) => {
                assert!(retflag != 0);
                assert!(times[0] > 0.0);
                assert_eq!(attr.len(), 20);
            }
            Err(Error::EclipseSearchExhausted { .. }) => {}
            Err(other) => panic!("unexpected local walk error: {other:?}"),
        }
    });
}

#[test]
fn planet_occult_servings_escalate_provider_misses_without_panic() {
    use crate::constants::{FLG_SWIEPH, JUPITER, SATURN, VENUS};
    with_servings(|serving, star_serving| {
        let glob = SealedPlanetOccultGlobalServing::try_new_planet(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH,
            VENUS,
            JUPITER,
        )
        .expect("planet pair constructs");
        // Past fixture coverage the miss escalates to the typed
        // range contract (never a panic, never a fabricated event).
        let err = glob
            .occulting_ra_dec_deg(2460000.5)
            .expect_err("out-of-range read must miss");
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
        // A classical planet absent from the asset misses the same
        // way inside coverage.
        let saturn = SealedPlanetOccultGlobalServing::try_new_planet(
            &serving,
            &star_serving,
            &zero_dt,
            FLG_SWIEPH,
            SATURN,
            JUPITER,
        )
        .expect("classical pair constructs without asset contact");
        let err = saturn
            .occulting_ra_dec_deg(1.0)
            .expect_err("unserved channel must miss");
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
        // The sealed search wrappers surface the miss as exhaustion
        // instead of an event when the walk steps off the asset.
        let err = planet_occult_when_glob_sealed(&glob, 2460000.5, 0)
            .expect_err("out-of-range walk must not yield an event");
        assert!(
            matches!(err, Error::EclipseSearchExhausted { .. }),
            "unexpected walk error: {err:?}"
        );
    });
}

/// Shared-core diagnostic (Venus-Jupiter residual): serving-objective words for
/// the 2065-11-22 Venus-Jupiter global maximum over the medium asset.
/// Emits the walk's exact objective (`planet_occult_separation_via`
/// over the sealed global serving) plus the underlying RA/Dec words
/// on a dense grid around the banked release maximum, so the
/// serving acquisition can be compared word-for-word against the
/// release and the C-ABI `calc_ut` path. Asset-gated: skips when
/// the machine-local asset is absent, never fails for missing data.
#[test]
fn planet_occult_vj2065_serving_objective_words() {
    use crate::constants::{FLG_SWIEPH, JUPITER, VENUS};
    let path = "/Volumes/ExFAT/libephemeris/leb/medium_core.leb2";
    let bytes = match std::fs::read(path) {
        Ok(bytes) => bytes,
        Err(_) => {
            eprintln!("planet_occult_vj2065_serving_objective_words: SKIP (missing {path})");
            return;
        }
    };
    let index = crate::leb::FileIndex::decode(&bytes).expect("medium LEB decodes");
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
    let glob = SealedPlanetOccultGlobalServing::try_new_planet(
        &serving,
        &star_serving,
        &deltat_at,
        FLG_SWIEPH,
        VENUS,
        JUPITER,
    )
    .expect("VJ pair constructs");
    // Banked release global maximum (O-B pin side).
    let center = f64::from_bits(0x4142e32e03ddba68);
    for k in -8..=8 {
        let t = center + k as f64 * 5e-8;
        let (ora, odec) = glob.occulting_ra_dec_deg(t).expect("occ ra/dec serves");
        let (tra, tdec) = glob.target_ra_dec_deg(t).expect("tgt ra/dec serves");
        let sep = planet_occult_separation_via(&glob, t).expect("sep serves");
        eprintln!(
            "VJOBJ t={:016x} ora={:016x} odec={:016x} tra={:016x} tdec={:016x} sep={:016x}",
            t.to_bits(),
            ora.to_bits(),
            odec.to_bits(),
            tra.to_bits(),
            tdec.to_bits(),
            sep.to_bits()
        );
    }
    // Full sealed walk ground truth on this build from the O-B start.
    let (retflag, tret) =
        planet_occult_when_glob_sealed(&glob, 2475590.5, 0).expect("VJ walk serves");
    eprintln!(
        "VJWALK ret={} tret=({:016x},{:016x},{:016x},{:016x},{:016x},{:016x})",
        retflag,
        tret[0].to_bits(),
        tret[1].to_bits(),
        tret[2].to_bits(),
        tret[3].to_bits(),
        tret[4].to_bits(),
        tret[5].to_bits(),
    );
    // Internal-consistency leg: the production `calc` engine over the
    // same asset bytes must serve bitwise-identical equatorial words
    // at the same TT (LOC-CLOSE discriminant shape). Any gap here is
    // a Rust serving-path defect, not release behavior. The center
    // serving words are fetched before the asset move below: `glob`
    // borrows the bytes through the composite/nutation closures.
    let center_tt = glob.jd_tt_for_ut(center).expect("tt serves");
    let (center_ora, center_odec) = glob
        .occulting_ra_dec_deg(center)
        .expect("center occ serves");
    let (center_tra, center_tdec) = glob.target_ra_dec_deg(center).expect("center tgt serves");
    // Frame-matrix leg: native Vondrak PN matrix + nutation words at
    // the event TT, for word comparison against the release
    // `_get_leb_frame_data` ( frame-rotation bisection). Runs
    // before the asset move below: the closure borrows the bytes.
    let (dpsi, deps) = nutation_at(center_tt).expect("nutation serves");
    let frames =
        crate::positions::frame_matrices_for_date(center_tt, dpsi, deps).expect("frames serve");
    let pn = frames.pn();
    eprintln!(
        "VJFRM dpsi={:016x} deps={:016x} pn0={:016x},{:016x},{:016x} pn1={:016x},{:016x},{:016x} pn2={:016x},{:016x},{:016x}",
        dpsi.to_bits(),
        deps.to_bits(),
        pn[0][0].to_bits(),
        pn[0][1].to_bits(),
        pn[0][2].to_bits(),
        pn[1][0].to_bits(),
        pn[1][1].to_bits(),
        pn[1][2].to_bits(),
        pn[2][0].to_bits(),
        pn[2][1].to_bits(),
        pn[2][2].to_bits(),
    );
    // Frame-bisection legs retain the actual mean-obliquity and nutation words.
    let eps_mean = crate::positions::vondrak_mean_obliquity_rad(center_tt).expect("eps serves");
    let nmat = crate::positions::nutation_matrix(eps_mean, dpsi, deps).expect("n serves");
    eprintln!(
        "VJNPM eps={:016x} n0={:016x},{:016x},{:016x} n1={:016x},{:016x},{:016x} n2={:016x},{:016x},{:016x}",
        eps_mean.to_bits(),
        nmat[0][0].to_bits(),
        nmat[0][1].to_bits(),
        nmat[0][2].to_bits(),
        nmat[1][0].to_bits(),
        nmat[1][1].to_bits(),
        nmat[1][2].to_bits(),
        nmat[2][0].to_bits(),
        nmat[2][1].to_bits(),
        nmat[2][2].to_bits(),
    );
    for tt in [center_tt, 2451545.0] {
        let pmat_j = crate::positions::vondrak_precession_matrix(tt, true).expect("p serves");
        eprintln!(
            "VJPRM t={:016x} p0={:016x},{:016x},{:016x} p1={:016x},{:016x},{:016x} p2={:016x},{:016x},{:016x}",
            tt.to_bits(),
            pmat_j[0][0].to_bits(),
            pmat_j[0][1].to_bits(),
            pmat_j[0][2].to_bits(),
            pmat_j[1][0].to_bits(),
            pmat_j[1][1].to_bits(),
            pmat_j[1][2].to_bits(),
            pmat_j[2][0].to_bits(),
            pmat_j[2][1].to_bits(),
            pmat_j[2][2].to_bits(),
        );
    }
    let leb = crate::leb::LebFile::from_bytes("mem://vj2065.leb2".to_string(), bytes)
        .expect("asset attaches");
    let mut ctx =
        crate::context::EphemerisContext::new().with_mode(crate::context::ProviderMode::SealedLeb);
    ctx.attach_leb_file(leb);
    let (calc_v, calc_j, ecl_v, ecl_j) = crate::positions::with_sealed_calc_engine(
        &ctx,
        crate::time::DeltaTModel::EspenakMeeus,
        false,
        |engine| {
            use crate::constants::{FLG_EQUATORIAL, FLG_SPEED, FLG_SWIEPH};
            let (venus, _) = crate::positions::calc(
                engine,
                center_tt,
                VENUS,
                FLG_SWIEPH | FLG_EQUATORIAL | FLG_SPEED,
            )?;
            let (jupiter, _) = crate::positions::calc(
                engine,
                center_tt,
                JUPITER,
                FLG_SWIEPH | FLG_EQUATORIAL | FLG_SPEED,
            )?;
            // Ecliptic leg for the bisection: lon/lat words
            // isolate the LEB-core acquisition from the equatorial
            // frame rotation.
            let (ecl_venus, _) =
                crate::positions::calc(engine, center_tt, VENUS, FLG_SWIEPH | FLG_SPEED)?;
            let (ecl_jupiter, _) =
                crate::positions::calc(engine, center_tt, JUPITER, FLG_SWIEPH | FLG_SPEED)?;
            Ok((venus, jupiter, ecl_venus, ecl_jupiter))
        },
    )
    .expect("shared calc engine serves");
    eprintln!(
        "VJCALC t={:016x} vra={:016x} vdec={:016x} jra={:016x} jdec={:016x}",
        center_tt.to_bits(),
        calc_v[0].to_bits(),
        calc_v[1].to_bits(),
        calc_j[0].to_bits(),
        calc_j[1].to_bits(),
    );
    eprintln!(
        "VJECL t={:016x} vlon={:016x} vlat={:016x} jlon={:016x} jlat={:016x}",
        center_tt.to_bits(),
        ecl_v[0].to_bits(),
        ecl_v[1].to_bits(),
        ecl_j[0].to_bits(),
        ecl_j[1].to_bits(),
    );
    // Is a Rust serving-path defect; agreement exonerates the
    assert_eq!(
        center_ora.to_bits(),
        calc_v[0].to_bits(),
        "occulting RA serving-vs-calc word"
    );
    assert_eq!(
        center_odec.to_bits(),
        calc_v[1].to_bits(),
        "occulting Dec serving-vs-calc word"
    );
    assert_eq!(
        center_tra.to_bits(),
        calc_j[0].to_bits(),
        "target RA serving-vs-calc word"
    );
    assert_eq!(
        center_tdec.to_bits(),
        calc_j[1].to_bits(),
        "target Dec serving-vs-calc word"
    );
    // Independent pinned values remain angular comparisons. Propagation
    // retains a fractional emission epoch, so bit identity is not physical.
    for (name, actual, expected) in [
        (
            "Venus longitude",
            ecl_v[0],
            f64::from_bits(0x406d1e584f333415),
        ),
        (
            "Jupiter longitude",
            ecl_j[0],
            f64::from_bits(0x406d1e5956f003b1),
        ),
        (
            "nutation longitude",
            dpsi.to_degrees(),
            f64::from_bits(0x3f12b3187d4005fa).to_degrees(),
        ),
        (
            "nutation obliquity",
            deps.to_degrees(),
            f64::from_bits(0x3eefcbf98bec8e5b).to_degrees(),
        ),
    ] {
        assert!(
            (actual - expected).abs() * 3600.0 <= 0.001,
            "{name}: {actual} vs pinned {expected}"
        );
    }
    assert_eq!(retflag, 16, "VJ2065 sealed global walk retflag");
}
